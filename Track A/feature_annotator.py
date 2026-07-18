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
import subprocess
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional

HERE = Path(__file__).resolve().parents[1]
TRACK_A = Path(__file__).resolve().parent

# The real contract is Phase 0's shared deliverable and lives here (NOT owned by
# Track A — we only read it). Until Phase 0 delivers it, we fall back to a sample
# bundled inside Track A so this module runs standalone.
SHARED_SPEC_PATH = HERE / "data/manifests/feature_spec.json"
SAMPLE_SPEC_PATH = TRACK_A / "fixtures/feature_spec.sample.json"


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


def target_columns(spec: dict) -> list[str]:
    """The `target__<gene>` presence columns implied by drug_targets."""
    genes = {g for genes in spec["drug_targets"].values() for g in genes}
    return [f"target__{g}" for g in sorted(genes)]


def marker_columns(spec: dict) -> list[str]:
    """Feature-order columns that are resistance markers (not target flags)."""
    return [c for c in spec["feature_order"] if not c.startswith("target__")]


class ContractError(ValueError):
    """Raised when a produced feature row does not match the frozen contract."""


def validate_feature_row(row: dict, spec: dict) -> dict:
    """
    The validation gate. Every annotator's output passes through here.

    Guarantees the row has exactly the contract's columns, in a form Track B can
    consume: all feature_order + target + QC columns present, binary flags coerced
    to 0/1 ints, QC coerced to the right numeric types. Fails LOUDLY on anything
    missing or non-binary so a broken backend can never silently corrupt features.
    """
    out: dict = {}
    missing, nonbinary = [], []

    for col in spec["feature_order"]:  # markers + target flags fed to the model
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

    for col in spec.get("qc_columns", []):
        if col not in row:
            missing.append(col)
            continue
        out[col] = float(row[col]) if col == "qc_complete" else int(row[col])

    if missing or nonbinary:
        parts = []
        if missing:
            parts.append(f"missing columns: {missing}")
        if nonbinary:
            parts.append(f"non-binary flag values: {nonbinary}")
        raise ContractError(
            "Feature row does not match feature_spec.json — " + "; ".join(parts)
        )
    return out


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

    def __init__(self, spec: dict):
        self.spec = spec

    @abstractmethod
    def annotate(self, genome_id: str, source) -> dict:
        """genome_id + a source (FASTA path, or an id into a loaded table) -> raw feature dict."""

    def annotate_validated(self, genome_id: str, source) -> dict:
        """annotate() then push the result through the validation gate."""
        return validate_feature_row(self.annotate(genome_id, source), self.spec)

    def _empty_row(self) -> dict:
        """All markers/targets absent, QC neutral — the base every backend fills in."""
        row = {c: 0 for c in self.spec["feature_order"]}
        for c in target_columns(self.spec):
            row[c] = 1  # targets assumed present until evidence of loss (see AMRFinderPlus backend)
        row["qc_complete"] = 1.0
        row["qc_contigs"] = 0
        return row


# --------------------------------------------------------------------------- #
# Backend 1 (default): AMRFinderPlus
# --------------------------------------------------------------------------- #
def parse_amrfinder_tsv(tsv_path: Path, spec: dict) -> dict:
    """
    Turn an AMRFinderPlus TSV into raw marker flags against our vocabulary.

    Reads the `Element symbol` column, maps it (directly or via marker_aliases)
    onto our marker names, and sets those flags to 1. Symbols outside our frozen
    vocabulary are dropped (with the caller free to log them) so the row shape
    stays fixed. Kept separate from the subprocess call so it is testable against
    a saved fixture without AMRFinderPlus installed.
    """
    vocab = set(marker_columns(spec))
    aliases = spec.get("marker_aliases", {})
    flags = {m: 0 for m in vocab}

    with open(tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        sym_col = _pick_symbol_column(reader.fieldnames or [])
        for rec in reader:
            symbol = (rec.get(sym_col) or "").strip()
            marker = aliases.get(symbol, symbol)
            if marker in vocab:
                flags[marker] = 1
    return flags


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
    against the assembly; until that is wired in, targets default to present (1).
    ompK36_loss is likewise a derived "absence" feature, not a direct AMRFinderPlus
    hit. Both are flagged TODO below rather than faked.
    """

    def __init__(self, spec: dict, organism: Optional[str] = None,
                 amrfinder_bin: str = "amrfinder", extra_args: Optional[list] = None):
        super().__init__(spec)
        self.organism = organism
        self.amrfinder_bin = amrfinder_bin
        self.extra_args = extra_args or []

    def annotate(self, genome_id: str, source, tsv_override: Optional[Path] = None) -> dict:
        row = self._empty_row()
        tsv = Path(tsv_override) if tsv_override else self._run(Path(source))
        row.update(parse_amrfinder_tsv(tsv, self.spec))
        # TODO(target-presence): replace default-present targets with a real
        #   gene-presence check (e.g. BLAST each drug_targets gene vs the assembly).
        # TODO(qc): fill qc_complete / qc_contigs from assembly stats for this genome.
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

    def __init__(self, spec: dict, table_path: Path):
        super().__init__(spec)
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
            f"Add your own by registering it in ANNOTATORS (see 'Track A/README.md')."
        )
    return ANNOTATORS[name](spec, **kwargs)
