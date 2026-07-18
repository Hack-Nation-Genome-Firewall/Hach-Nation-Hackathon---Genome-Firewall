"""BV-BRC API client and pure cohort-building functions.

The API client records exact response bytes and request metadata. Cleaning and
selection are separate pure functions so scientific policy can be tested without
network access.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence
from urllib.parse import quote

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


AMR_REQUIRED_FIELDS = {
    "id",
    "genome_id",
    "antibiotic",
    "resistant_phenotype",
    "evidence",
    "taxon_id",
}
GENOME_METADATA_FIELDS = [
    "genome_id",
    "genome_name",
    "taxon_id",
    "contigs",
    "genome_length",
    "checkm_completeness",
    "checkm_contamination",
    "genome_quality_flags",
    "genome_status",
    "assembly_accession",
    "date_modified",
]
CONFLICT_COLUMNS = [
    "genome_id",
    "antibiotic",
    "phenotypes",
    "record_count",
    "record_ids",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_json(value: Any) -> str:
    body = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return sha256_bytes(body)


def load_config(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    required = {"schema_version", "source", "species", "antibiotics", "label_policy"}
    missing = sorted(required - config.keys())
    if missing:
        raise ValueError(f"Configuration is missing keys: {missing}")
    if config["source"].get("evidence") != "Laboratory Method":
        raise ValueError("The source evidence filter must be exactly 'Laboratory Method'")
    if not config["antibiotics"]:
        raise ValueError("At least one antibiotic must be configured")
    return config


def rql_eq(field: str, value: Any) -> str:
    return f"eq({field},{quote(str(value), safe='._-/')})"


def rql_in(field: str, values: Sequence[Any]) -> str:
    encoded = ",".join(quote(str(value), safe="._-/") for value in values)
    return f"in({field},({encoded}))"


@dataclass(frozen=True)
class RequestRecord:
    url: str
    retrieved_at: str
    status_code: int
    etag: str | None
    content_range: str | None
    sha256: str
    bytes: int
    records: int
    snapshot_path: str | None


class BvbrcClient:
    """Small read-only client for the public BV-BRC HTTPS API."""

    def __init__(
        self,
        api_base: str,
        *,
        timeout_seconds: float = 60.0,
        page_size: int = 1000,
        session: requests.Session | None = None,
    ) -> None:
        self.api_base = api_base.rstrip("/")
        self.timeout_seconds = timeout_seconds
        self.page_size = page_size
        self.session = session or requests.Session()
        retry = Retry(
            total=5,
            connect=5,
            read=5,
            status=5,
            backoff_factor=0.75,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset({"GET"}),
            respect_retry_after_header=True,
        )
        adapter = HTTPAdapter(max_retries=retry)
        self.session.mount("https://", adapter)
        self.session.headers.update(
            {
                "Accept": "application/json",
                "User-Agent": "GenomeFirewall/0.1 reproducible-research-pipeline",
            }
        )

    def build_url(self, resource: str, clauses: Sequence[str]) -> str:
        return f"{self.api_base}/{resource.strip('/')}/?{'&'.join(clauses)}"

    def get_json(self, resource: str, clauses: Sequence[str]) -> tuple[list[dict[str, Any]], RequestRecord]:
        url = self.build_url(resource, clauses)
        response = self.session.get(url, timeout=self.timeout_seconds)
        response.raise_for_status()
        try:
            records = response.json()
        except requests.JSONDecodeError as exc:
            raise ValueError(f"BV-BRC returned invalid JSON for {url}") from exc
        if not isinstance(records, list):
            raise ValueError(f"BV-BRC returned {type(records).__name__}, expected a list")
        request_record = RequestRecord(
            url=url,
            retrieved_at=utc_now(),
            status_code=response.status_code,
            etag=response.headers.get("ETag"),
            content_range=response.headers.get("Content-Range"),
            sha256=sha256_bytes(response.content),
            bytes=len(response.content),
            records=len(records),
            snapshot_path=None,
        )
        return records, request_record

    def fetch_pages(
        self,
        resource: str,
        filters: Sequence[str],
        *,
        select: Sequence[str] | None = None,
        sort: Sequence[str] = ("genome_id", "id"),
        snapshot_dir: Path | None = None,
        snapshot_prefix: str = "page",
        overwrite: bool = False,
    ) -> tuple[list[dict[str, Any]], list[RequestRecord]]:
        all_records: list[dict[str, Any]] = []
        requests_log: list[RequestRecord] = []
        offset = 0
        page_index = 0
        total: int | None = None
        if snapshot_dir:
            snapshot_dir.mkdir(parents=True, exist_ok=True)

        while total is None or offset < total:
            clauses = list(filters)
            if select:
                clauses.append(f"select({','.join(select)})")
            if sort:
                sort_fields = ",".join(f"%2B{field}" for field in sort)
                clauses.append(f"sort({sort_fields})")
            clauses.append(f"limit({self.page_size},{offset})")
            url = self.build_url(resource, clauses)
            response = self.session.get(url, timeout=self.timeout_seconds)
            response.raise_for_status()
            try:
                page = response.json()
            except requests.JSONDecodeError as exc:
                raise ValueError(f"BV-BRC returned invalid JSON for {url}") from exc
            if not isinstance(page, list):
                raise ValueError(f"BV-BRC returned {type(page).__name__}, expected a list")

            snapshot_path: Path | None = None
            if snapshot_dir:
                snapshot_path = snapshot_dir / f"{snapshot_prefix}_{page_index:05d}.json"
                if snapshot_path.exists() and not overwrite:
                    raise FileExistsError(
                        f"Snapshot already exists: {snapshot_path}. Use --overwrite intentionally."
                    )
                _write_bytes_atomic(snapshot_path, response.content)

            content_range = response.headers.get("Content-Range")
            total = _total_from_content_range(content_range)
            request_record = RequestRecord(
                url=url,
                retrieved_at=utc_now(),
                status_code=response.status_code,
                etag=response.headers.get("ETag"),
                content_range=content_range,
                sha256=sha256_bytes(response.content),
                bytes=len(response.content),
                records=len(page),
                snapshot_path=str(snapshot_path) if snapshot_path else None,
            )
            requests_log.append(request_record)
            all_records.extend(page)

            if not page or len(page) < self.page_size:
                break
            offset += len(page)
            page_index += 1

        return all_records, requests_log


def _total_from_content_range(value: str | None) -> int | None:
    if not value:
        return None
    match = re.search(r"/(\d+)\s*$", value)
    return int(match.group(1)) if match else None


def _write_bytes_atomic(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_bytes(content)
    os.replace(temporary, path)


def write_json_atomic(path: Path, value: Any) -> None:
    encoded = (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True) + "\n").encode("utf-8")
    _write_bytes_atomic(path, encoded)


def write_csv_atomic(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    frame.to_csv(temporary, index=False, lineterminator="\n")
    os.replace(temporary, path)


def clean_amr_records(
    records: Iterable[dict[str, Any]], config: dict[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Validate, filter, and collapse BV-BRC AMR records.

    Conflicting S/R records for one genome/drug pair are returned separately and
    never enter the modeling label frame.
    """

    expected_taxon = int(config["species"]["taxon_id"])
    expected_evidence = config["source"]["evidence"]
    configured_drugs = {drug.casefold(): drug for drug in config["antibiotics"]}
    included = set(config["label_policy"]["included"])
    excluded = set(config["label_policy"]["excluded"])
    stats = {
        "source_records": 0,
        "laboratory_records": 0,
        "wrong_evidence": 0,
        "wrong_taxon": 0,
        "unconfigured_antibiotic": 0,
        "excluded_phenotype": 0,
        "invalid_phenotype": 0,
        "accepted_records_before_collapse": 0,
        "duplicate_records_collapsed": 0,
        "conflicting_pairs_excluded": 0,
        "clean_pairs": 0,
    }
    accepted: list[dict[str, Any]] = []

    for raw in records:
        stats["source_records"] += 1
        missing = sorted(field for field in AMR_REQUIRED_FIELDS if field not in raw)
        if missing:
            raise ValueError(f"AMR record {raw.get('id', '<unknown>')} is missing fields: {missing}")
        if raw["evidence"] != expected_evidence:
            stats["wrong_evidence"] += 1
            continue
        stats["laboratory_records"] += 1
        if int(raw["taxon_id"]) != expected_taxon:
            stats["wrong_taxon"] += 1
            continue
        drug_key = str(raw["antibiotic"]).strip().casefold()
        if drug_key not in configured_drugs:
            stats["unconfigured_antibiotic"] += 1
            continue
        phenotype = str(raw["resistant_phenotype"]).strip()
        if phenotype in excluded:
            stats["excluded_phenotype"] += 1
            continue
        if phenotype not in included:
            stats["invalid_phenotype"] += 1
            continue
        normalized = dict(raw)
        normalized["genome_id"] = str(raw["genome_id"])
        normalized["antibiotic"] = configured_drugs[drug_key]
        normalized["resistant_phenotype"] = phenotype
        accepted.append(normalized)

    stats["accepted_records_before_collapse"] = len(accepted)
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for record in accepted:
        key = (record["genome_id"], record["antibiotic"])
        grouped.setdefault(key, []).append(record)

    clean_rows: list[dict[str, Any]] = []
    conflicts: list[dict[str, Any]] = []
    for (genome_id, antibiotic), group in sorted(grouped.items()):
        phenotypes = sorted({record["resistant_phenotype"] for record in group})
        record_ids = sorted({str(record["id"]) for record in group})
        if len(phenotypes) != 1:
            conflicts.append(
                {
                    "genome_id": genome_id,
                    "antibiotic": antibiotic,
                    "phenotypes": ";".join(phenotypes),
                    "record_count": len(group),
                    "record_ids": ";".join(record_ids),
                }
            )
            continue
        phenotype = phenotypes[0]
        clean_rows.append(
            {
                "genome_id": genome_id,
                "antibiotic": antibiotic,
                "phenotype": phenotype,
                "label": 1 if phenotype == "Resistant" else 0,
                "evidence": expected_evidence,
                "record_count": len(group),
                "record_ids": ";".join(record_ids),
                "measurement": _joined_values(group, "measurement"),
                "measurement_unit": _joined_values(group, "measurement_unit"),
                "laboratory_typing_method": _joined_values(group, "laboratory_typing_method"),
                "testing_standard": _joined_values(group, "testing_standard"),
                "testing_standard_year": _joined_values(group, "testing_standard_year"),
                "pmid": _joined_values(group, "pmid"),
                "source_date_modified": _joined_values(group, "date_modified"),
            }
        )
        stats["duplicate_records_collapsed"] += len(group) - 1

    stats["conflicting_pairs_excluded"] = len(conflicts)
    stats["clean_pairs"] = len(clean_rows)
    clean = pd.DataFrame(clean_rows)
    conflict_frame = pd.DataFrame(conflicts, columns=CONFLICT_COLUMNS)
    return clean, conflict_frame, stats


def _joined_values(records: Sequence[dict[str, Any]], field: str) -> str:
    values: set[str] = set()
    for record in records:
        value = record.get(field)
        if value is None or value == "":
            continue
        items = value if isinstance(value, list) else [value]
        for item in items:
            if item is not None and str(item).strip():
                values.add(str(item).strip())
    return ";".join(sorted(values))


def evaluate_genome_quality(
    metadata: pd.DataFrame, quality_policy: dict[str, Any]
) -> pd.DataFrame:
    frame = metadata.copy()
    required = [
        "checkm_completeness",
        "checkm_contamination",
        "contigs",
        "genome_length",
    ]
    for field in required:
        if field not in frame.columns:
            frame[field] = pd.NA

    evaluations: list[tuple[bool, str]] = []
    for row in frame.to_dict(orient="records"):
        reasons: list[str] = []
        missing = [field for field in required if _is_missing(row.get(field))]
        if missing and quality_policy.get("missing_quality_policy") == "exclude":
            reasons.append(f"missing:{','.join(missing)}")
        if not _is_missing(row.get("checkm_completeness")) and float(row["checkm_completeness"]) < float(
            quality_policy["minimum_checkm_completeness"]
        ):
            reasons.append("completeness_below_minimum")
        if not _is_missing(row.get("checkm_contamination")) and float(row["checkm_contamination"]) > float(
            quality_policy["maximum_checkm_contamination"]
        ):
            reasons.append("contamination_above_maximum")
        if not _is_missing(row.get("contigs")) and int(row["contigs"]) > int(quality_policy["maximum_contigs"]):
            reasons.append("contigs_above_maximum")
        if not _is_missing(row.get("genome_length")):
            length = int(row["genome_length"])
            if length < int(quality_policy["minimum_genome_length"]):
                reasons.append("genome_too_short")
            if length > int(quality_policy["maximum_genome_length"]):
                reasons.append("genome_too_long")
        evaluations.append((not reasons, ";".join(reasons)))

    frame["quality_pass"] = [value[0] for value in evaluations]
    frame["quality_reasons"] = [value[1] for value in evaluations]
    return frame


def _is_missing(value: Any) -> bool:
    return value is None or (isinstance(value, float) and math.isnan(value)) or value is pd.NA


def select_cohort(
    labels: pd.DataFrame,
    metadata: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    if labels.empty:
        raise ValueError("No clean laboratory labels are available")
    if metadata.empty:
        raise ValueError("No genome metadata is available")
    if metadata["genome_id"].astype(str).duplicated().any():
        raise ValueError("Genome metadata contains duplicate genome_id rows")

    quality = evaluate_genome_quality(metadata, config["quality_policy"])
    quality["genome_id"] = quality["genome_id"].astype(str)
    eligible_ids = set(quality.loc[quality["quality_pass"], "genome_id"])
    candidate_labels = labels[labels["genome_id"].astype(str).isin(eligible_ids)].copy()
    counts = candidate_labels.groupby("genome_id")["antibiotic"].nunique()
    minimum_drugs = int(config["selection_policy"]["minimum_drugs_per_genome"])
    eligible_panel_ids = set(counts[counts >= minimum_drugs].index.astype(str))

    seed = int(config["selection_policy"]["stable_seed"])
    ranked = quality[quality["genome_id"].isin(eligible_panel_ids)].copy()
    ranked["tested_drug_count"] = ranked["genome_id"].map(counts).astype(int)
    ranked["stable_rank"] = ranked["genome_id"].map(
        lambda genome_id: hashlib.sha256(f"{seed}:{genome_id}".encode("utf-8")).hexdigest()
    )
    ranked = ranked.sort_values(
        ["tested_drug_count", "stable_rank", "genome_id"],
        ascending=[False, True, True],
        kind="mergesort",
    )
    maximum = int(config["selection_policy"]["maximum_genomes"])
    selected_genomes = ranked.head(maximum).reset_index(drop=True)
    selected_ids = set(selected_genomes["genome_id"])
    selected_labels = candidate_labels[candidate_labels["genome_id"].isin(selected_ids)].copy()
    selected_labels = selected_labels.sort_values(["genome_id", "antibiotic"], kind="mergesort").reset_index(drop=True)

    class_counts = (
        selected_labels.groupby(["antibiotic", "phenotype"]).size().rename("count").reset_index()
    )
    for drug in config["antibiotics"]:
        drug_classes = set(class_counts.loc[class_counts["antibiotic"] == drug, "phenotype"])
        if drug_classes != {"Resistant", "Susceptible"}:
            raise ValueError(f"Selected cohort does not contain both classes for {drug}: {sorted(drug_classes)}")

    stats = {
        "metadata_rows": int(len(metadata)),
        "quality_pass_genomes": int(quality["quality_pass"].sum()),
        "quality_fail_genomes": int((~quality["quality_pass"]).sum()),
        "minimum_drugs_per_genome": minimum_drugs,
        "panel_eligible_genomes": int(len(ranked)),
        "selected_genomes": int(len(selected_genomes)),
        "selected_label_pairs": int(len(selected_labels)),
        "class_counts": class_counts.to_dict(orient="records"),
    }
    return selected_genomes, selected_labels, stats


def fetch_genome_metadata(
    client: BvbrcClient,
    genome_ids: Sequence[str],
    resource: str,
    snapshot_dir: Path,
    *,
    chunk_size: int = 100,
    overwrite: bool = False,
) -> tuple[list[dict[str, Any]], list[RequestRecord]]:
    records: list[dict[str, Any]] = []
    requests_log: list[RequestRecord] = []
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    ordered_ids = sorted(set(str(genome_id) for genome_id in genome_ids))
    for offset in range(0, len(ordered_ids), chunk_size):
        chunk = ordered_ids[offset : offset + chunk_size]
        clauses = [
            rql_in("genome_id", chunk),
            f"select({','.join(GENOME_METADATA_FIELDS)})",
            "sort(%2Bgenome_id)",
            f"limit({len(chunk)},0)",
        ]
        page, request_record = client.get_json(resource, clauses)
        path = snapshot_dir / f"metadata_{offset // chunk_size:05d}.json"
        if path.exists() and not overwrite:
            raise FileExistsError(f"Snapshot already exists: {path}. Use --overwrite intentionally.")
        write_json_atomic(path, page)
        request_record = RequestRecord(
            **{
                **asdict(request_record),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "snapshot_path": str(path),
            }
        )
        records.extend(page)
        requests_log.append(request_record)
    return records, requests_log


def current_git_commit(root: Path) -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def runtime_versions() -> dict[str, str]:
    return {
        "python": sys.version.split()[0],
        "pandas": pd.__version__,
        "requests": requests.__version__,
    }


def request_records_to_dict(records: Iterable[RequestRecord]) -> list[dict[str, Any]]:
    return [asdict(record) for record in records]


def path_for_record(path: Path, root: Path) -> str:
    """Return a portable repository-relative path when possible."""

    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())
