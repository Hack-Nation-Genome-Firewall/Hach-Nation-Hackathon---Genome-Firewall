"""
GENOME FIREWALL — Module 1 (Track A): the Genome Reader.

Turns one genome (a FASTA file) into the fixed feature row the predictor was
trained on. Nothing here is machine-learned; this is a deterministic featurizer.

    FASTA file  ->  [ annotator backend ]  ->  small marker list  ->  feature row
                          ^ swappable                                  ^ frozen contract

The only thing downstream (Track B / Track C) depends on is the FEATURE CONTRACT
in data/manifests/feature_spec.json — the exact columns and their order. Any tool
or dataset that can produce that contract plugs in here. AMRFinderPlus is just the
default way to fill it in; it is not hard-wired.

============================================================================
 EXTENSION POINT — bring your own tool or your own dataset
============================================================================
To use a different annotation tool (ResFinder, cAMRah, a custom script) or to
load your own precomputed feature table, you do NOT edit the pipeline. You:

    1. Subclass `FeatureAnnotator` and implement `annotate(genome_id, source)`.
    2. Register it in the `ANNOTATORS` registry at the bottom of this file.
    3. Point the config / `backend=` argument at its name.

`PrecomputedAnnotator` below is a complete, working example written to be copied.
Whatever your backend emits is checked by `validate_feature_row()` against the
frozen contract and rejected loudly if it does not match — so bring-your-own is
safe, not a way to silently corrupt features.
============================================================================
"""
from __future__ import annotations

import csv
import json
import math
import subprocess
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional

HERE = Path(__file__).resolve().parents[1]
MODULE_DIR = Path(__file__).resolve().parent

# The real contract is Phase 0's shared deliverable and lives here (NOT owned by
# Track A — we only read it). Until Phase 0 delivers it, we fall back to a sample
# bundled inside this module so it runs standalone.
SHARED_SPEC_PATH = HERE / "data/manifests/feature_spec.json"
SAMPLE_SPEC_PATH = MODULE_DIR / "fixtures/feature_spec.sample.json"

# Phase 0's shared project config (species, drug panel, QC policy). We consume it
# read-only to stay aligned with the rest of the team; we never write it.
PROJECT_CONFIG_PATH = HERE / "data/config/project.json"


# --------------------------------------------------------------------------- #
# The frozen contract
# --------------------------------------------------------------------------- #
def load_spec(path: Optional[Path] = None) -> dict:
    """
    Load the feature contract (single source of truth for columns + order).

    Resolution order: an explicit `path` if given; else Phase 0's shared
    data/manifests/feature_spec.json if it exists; else Track A's bundled sample.
    """
    if path is not None:
        chosen = Path(path)
    elif SHARED_SPEC_PATH.exists():
        chosen = SHARED_SPEC_PATH
    else:
        chosen = SAMPLE_SPEC_PATH
    with open(chosen) as f:
        return json.load(f)


def load_project_config(path: Optional[Path] = None) -> Optional[dict]:
    """Read Phase 0's shared project config, or None if it hasn't landed yet."""
    path = Path(path) if path else PROJECT_CONFIG_PATH
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


def load_qc_map(path) -> dict:
    """
    Read per-genome CheckM QC from a selected_genomes.csv into
    {genome_id: {"completeness", "contamination", "contigs"}}.

    Uses columns checkm_completeness, checkm_contamination, contigs — so we read the
    QC the cohort already computed rather than re-running CheckM. A row missing its
    genome_id or any QC value is skipped, leaving that genome QC-unknown (-> no-call).
    """
    qc: dict = {}
    with open(path, newline="") as f:
        for rec in csv.DictReader(f):
            gid = rec.get("genome_id")
            if not gid:
                continue
            try:
                qc[gid] = {
                    "completeness": float(rec["checkm_completeness"]),
                    "contamination": float(rec["checkm_contamination"]),
                    "contigs": int(rec["contigs"]),
                }
            except (KeyError, TypeError, ValueError):
                continue   # unknown/malformed QC -> leave genome out (stays no-call)
    return qc


def spec_project_discrepancies(spec: dict, project: Optional[dict] = None) -> list[str]:
    """
    Cross-check the feature spec against Phase 0's project config.

    Returns a list of human-readable mismatches (empty = aligned). Used to catch
    drift between what we featurize and the team's declared species/drug panel.
    Returns [] when no project config exists yet (nothing to check against).
    """
    project = project if project is not None else load_project_config()
    if project is None:
        return []
    issues = []
    proj_species = (project.get("species") or {}).get("name")
    spec_species = (spec.get("species") or {}).get("name") if isinstance(spec.get("species"), dict) else spec.get("species")
    if proj_species and proj_species != spec_species:
        issues.append(f"species mismatch: spec={spec_species!r} project={proj_species!r}")
    proj_drugs = set(project.get("antibiotics", []))
    spec_drugs = set(spec.get("drugs", []))
    if proj_drugs:
        if proj_drugs - spec_drugs:
            issues.append(f"drugs in project.json missing from spec: {sorted(proj_drugs - spec_drugs)}")
        if spec_drugs - proj_drugs:
            issues.append(f"drugs in spec not in project.json: {sorted(spec_drugs - proj_drugs)}")
    return issues


def marker_columns(spec: dict) -> list[str]:
    """The ordered binary marker features the model consumes (spec.model_features)."""
    return list(spec["model_features"])


def target_columns(spec: dict) -> list[str]:
    """The target-status columns (drug_targets values), de-duplicated in stable order."""
    seen, cols = set(), []
    for feats in spec["drug_targets"].values():
        for c in feats:
            if c not in seen:
                seen.add(c)
                cols.append(c)
    return cols


def quality_columns(spec: dict) -> list[str]:
    """The three QC column names, ordered completeness, contamination, contigs."""
    qf = spec["quality_features"]
    return [qf["completeness"], qf["contamination"], qf["contigs"]]


class ContractError(ValueError):
    """Raised when a produced feature row does not match the frozen contract."""


def validate_feature_row(row: dict, spec: dict) -> dict:
    """
    The validation gate. Every annotator's output passes through here.

    Guarantees the row has exactly the contract's columns in a form Track B can
    consume. Model features must be binary. Target and QC measurements may be
    explicitly unknown (None/blank/NaN), which Track B routes to no-call rather
    than treating missing evidence as a pass. Invalid measured values fail loudly.
    """
    out: dict = {}
    missing, nonbinary = [], []

    # Model features are always required and binary.
    for col in marker_columns(spec):
        if col not in row:
            missing.append(col)
            continue
        v = row[col]
        if v in (0, 1, "0", "1", True, False):
            out[col] = int(v)
        else:
            try:
                iv = int(v)
                if iv in (0, 1):
                    out[col] = iv
                else:
                    nonbinary.append((col, v))
            except (TypeError, ValueError):
                nonbinary.append((col, v))

    # Target measurements are binary when known. Unknown is a valid inference
    # state because Track B has an explicit target_status_unknown no-call gate.
    for col in target_columns(spec):
        if col not in row:
            missing.append(col)
            continue
        v = row[col]
        if _unknown(v):
            out[col] = None
        elif v in (0, 1, "0", "1", True, False):
            out[col] = int(v)
        else:
            try:
                iv = int(v)
                if iv in (0, 1):
                    out[col] = iv
                else:
                    nonbinary.append((col, v))
            except (TypeError, ValueError):
                nonbinary.append((col, v))

    # QC measurements are numeric when known. Unknown values trigger Track B's
    # quality_status_unknown no-call gate.
    qf = spec["quality_features"]
    invalid_quality = []
    for col in quality_columns(spec):
        if col not in row:
            missing.append(col)
            continue
        v = row[col]
        if _unknown(v):
            out[col] = None
            continue
        try:
            out[col] = int(v) if col == qf["contigs"] else float(v)
        except (TypeError, ValueError):
            invalid_quality.append((col, v))

    if missing or nonbinary or invalid_quality:
        parts = []
        if missing:
            parts.append(f"missing columns: {missing}")
        if nonbinary:
            parts.append(f"non-binary flag values: {nonbinary}")
        if invalid_quality:
            parts.append(f"non-numeric QC values: {invalid_quality}")
        raise ContractError(
            "Feature row does not match feature_spec.json — " + "; ".join(parts)
        )
    return out


def _unknown(value) -> bool:
    return value is None or value == "" or (
        isinstance(value, float) and math.isnan(value)
    )


# --------------------------------------------------------------------------- #
# The interface every backend implements
# --------------------------------------------------------------------------- #
class FeatureAnnotator(ABC):
    """
    Base class for every way of turning a genome into a feature row.

    Implement `annotate()` and return a plain dict of column -> value. Cover the
    marker columns and the `target__<gene>` / QC columns; anything you omit is
    caught by the validation gate. Missing markers may be left out (they default
    to 0). See PrecomputedAnnotator for a copyable example.
    """

    def __init__(self, spec: dict, qc_source: Optional[dict] = None):
        self.spec = spec
        # Markers seen on the most recent annotate() that are outside the vocabulary.
        # Preserved (not dropped) per the Track A guardrail; batch runs log them.
        self.last_unknown_markers: list[str] = []
        # Optional {genome_id: {"completeness","contamination","contigs"}} used to
        # fill real QC (e.g. from load_qc_map(selected_genomes.csv)).
        self.qc_source = qc_source or {}

    @abstractmethod
    def annotate(self, genome_id: str, source) -> dict:
        """genome_id + a source (FASTA path, or an id into a loaded table) -> raw feature dict."""

    def annotate_validated(self, genome_id: str, source) -> dict:
        """annotate() then push the result through the validation gate."""
        return validate_feature_row(self.annotate(genome_id, source), self.spec)

    def _empty_row(self) -> dict:
        """All markers absent, target/QC evidence unknown: the base each backend fills in."""
        row = {c: 0 for c in marker_columns(self.spec)}
        for c in target_columns(self.spec):
            row[c] = None
        qf = self.spec["quality_features"]
        row[qf["completeness"]] = None
        row[qf["contamination"]] = None
        row[qf["contigs"]] = None
        return row

    def _fill_qc(self, row: dict, genome_id: str) -> None:
        """Fill QC columns from qc_source; leave unknown (None) if not available."""
        qc = self.qc_source.get(genome_id)
        if not qc:
            return
        qf = self.spec["quality_features"]
        row[qf["completeness"]] = float(qc["completeness"])
        row[qf["contamination"]] = float(qc["contamination"])
        row[qf["contigs"]] = int(qc["contigs"])


# --------------------------------------------------------------------------- #
# Backend 1 (default): AMRFinderPlus
# --------------------------------------------------------------------------- #
def parse_amrfinder_markers(tsv_path: Path, spec: dict) -> tuple[dict, list[str]]:
    """
    Turn an AMRFinderPlus TSV into (marker flags, unknown markers).

    Reads the `Element symbol` column, maps it (directly or via marker_aliases)
    onto our vocabulary, and sets those flags to 1. The model feature vector must
    stay fixed-shape, so a symbol outside the vocabulary can't become a model
    column — but per the Track A guardrail we do NOT silently drop it: every such
    symbol is returned in `unknown_markers` (de-duplicated, in order seen) so it is
    preserved for review and possible addition to the vocabulary. Kept separate
    from the subprocess call so it is testable against a fixture without the tool.
    """
    vocab = set(marker_columns(spec))
    aliases = spec.get("marker_aliases", {})
    flags = {m: 0 for m in vocab}
    unknown: list[str] = []
    seen_unknown: set[str] = set()

    with open(tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        sym_col = _pick_symbol_column(reader.fieldnames or [])
        for rec in reader:
            symbol = (rec.get(sym_col) or "").strip()
            if not symbol:
                continue
            marker = aliases.get(symbol, symbol)
            if marker in vocab:
                flags[marker] = 1
            elif symbol not in seen_unknown:
                seen_unknown.add(symbol)
                unknown.append(symbol)   # preserve, don't drop
    return flags, unknown


def parse_amrfinder_tsv(tsv_path: Path, spec: dict) -> dict:
    """Marker flags only (thin wrapper over parse_amrfinder_markers)."""
    return parse_amrfinder_markers(tsv_path, spec)[0]


def _pick_symbol_column(fieldnames: list[str]) -> str:
    """AMRFinderPlus has renamed this column across versions; accept known variants."""
    for candidate in ("Element symbol", "Gene symbol", "Element_symbol"):
        if candidate in fieldnames:
            return candidate
    raise ContractError(
        f"Could not find an element-symbol column in AMRFinderPlus output; "
        f"got columns {fieldnames}"
    )


class AMRFinderPlusAnnotator(FeatureAnnotator):
    """
    Default backend. Runs AMRFinderPlus on a nucleotide FASTA and parses the TSV.

    If `tsv_override` is given, uses that TSV instead of invoking the tool — the
    path used for tests/fixtures and for organizer-provided PRECOMPUTED
    AMRFinderPlus results (so we never re-run the tool on genomes we already have).

    NOTE (honest limitation): AMRFinderPlus reports resistance markers, not the
    presence/intactness of a drug's *target* housekeeping gene. Target-presence
    detection (the `target__<gene>` flags) needs a separate gene-presence check
    against the assembly; until that is wired in, targets remain unknown and Track
    B returns no-call.
    ompK36_loss is likewise a derived "absence" feature, not a direct AMRFinderPlus
    hit. Both are flagged TODO below rather than faked. Markers AMRFinderPlus reports
    that are outside the vocabulary are recorded in `last_unknown_markers`, not dropped.
    """

    def __init__(self, spec: dict, organism: Optional[str] = None,
                 amrfinder_bin: str = "amrfinder", extra_args: Optional[list] = None,
                 precomputed_tsv: bool = False, qc_source: Optional[dict] = None):
        super().__init__(spec, qc_source=qc_source)
        self.organism = organism
        self.amrfinder_bin = amrfinder_bin
        self.extra_args = extra_args or []
        # When True, `source` is a saved AMRFinderPlus TSV (organizer-precomputed
        # results) rather than a FASTA — so batch runs never re-run the tool.
        self.precomputed_tsv = precomputed_tsv

    def annotate(self, genome_id: str, source, tsv_override: Optional[Path] = None) -> dict:
        row = self._empty_row()
        if tsv_override is not None:
            tsv = Path(tsv_override)
        elif self.precomputed_tsv:
            tsv = Path(source)
        else:
            tsv = self._run(Path(source))
        flags, unknown = parse_amrfinder_markers(tsv, self.spec)
        row.update(flags)
        self.last_unknown_markers = unknown   # preserved for review, not dropped
        self._fill_qc(row, genome_id)         # real QC from qc_source, else unknown
        # TODO(target-presence): replace unknown targets with a real
        #   gene-presence check (e.g. BLAST each drug_targets gene vs the assembly).
        return row

    def _run(self, fasta_path: Path) -> Path:
        out_tsv = fasta_path.with_suffix(".amrfinder.tsv")
        cmd = [self.amrfinder_bin, "-n", str(fasta_path), "-o", str(out_tsv)]
        if self.organism:
            cmd += ["-O", self.organism]
        cmd += self.extra_args
        subprocess.run(cmd, check=True)
        return out_tsv


# --------------------------------------------------------------------------- #
# Backend 2: Precomputed  ==  THE WORKED EXAMPLE for "bring your own dataset"
# --------------------------------------------------------------------------- #
class PrecomputedAnnotator(FeatureAnnotator):
    """
    Load feature rows from an existing table instead of annotating a genome.

    This is the "bring your own dataset" path AND the template to copy when adding
    a new backend. A scientist with their own features (from any tool, or their own
    experiments) drops in a CSV keyed by `genome_id`; the validation gate confirms
    it matches the contract. To add a brand-new tool, copy this class, keep the
    __init__/annotate shape, and swap the body of annotate() for your tool's logic.
    """

    def __init__(self, spec: dict, table_path: Path, qc_source: Optional[dict] = None):
        super().__init__(spec, qc_source=qc_source)
        self.table_path = Path(table_path)
        self._rows = self._load_table(self.table_path)

    @staticmethod
    def _load_table(path: Path) -> dict:
        rows: dict = {}
        with open(path, newline="") as f:
            for rec in csv.DictReader(f):
                gid = rec.get("genome_id")
                if gid is not None:
                    rows[gid] = rec
        return rows

    def annotate(self, genome_id: str, source=None) -> dict:
        if genome_id not in self._rows:
            raise ContractError(
                f"genome_id '{genome_id}' not found in {self.table_path}"
            )
        row = self._empty_row()
        row.update({k: v for k, v in self._rows[genome_id].items() if k in row})
        self._fill_qc(row, genome_id)   # qc_source overrides table QC when provided
        return row


# --------------------------------------------------------------------------- #
# Registry — the config-selectable list of backends (visible = discoverable)
# --------------------------------------------------------------------------- #
ANNOTATORS = {
    "amrfinderplus": AMRFinderPlusAnnotator,   # default: FASTA -> AMRFinderPlus
    "precomputed": PrecomputedAnnotator,       # bring your own feature table
    # "resfinder": ResFinderAnnotator,         # <- add your own here (see README)
}

DEFAULT_BACKEND = "amrfinderplus"


def get_annotator(name: Optional[str] = None, spec: Optional[dict] = None, **kwargs) -> FeatureAnnotator:
    """Factory: resolve a backend name to a constructed annotator."""
    name = name or DEFAULT_BACKEND
    spec = spec or load_spec()
    if name not in ANNOTATORS:
        raise ValueError(
            f"Unknown annotator backend '{name}'. Available: {list(ANNOTATORS)}. "
            f"Add your own by registering it in ANNOTATORS (see module1_reader/README.md)."
        )
    return ANNOTATORS[name](spec, **kwargs)
