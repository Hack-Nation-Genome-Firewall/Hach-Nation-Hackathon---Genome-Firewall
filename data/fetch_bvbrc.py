"""Fetch laboratory AMR records and build a frozen BV-BRC cohort manifest."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import pandas as pd

from data.bvbrc import (
    BvbrcClient,
    clean_amr_records,
    current_git_commit,
    fetch_genome_metadata,
    load_config,
    path_for_record,
    request_records_to_dict,
    rql_eq,
    runtime_versions,
    select_cohort,
    sha256_file,
    sha256_json,
    utc_now,
    write_csv_atomic,
    write_json_atomic,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "data/config/project.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--raw-dir", type=Path, default=ROOT / "data/raw/bvbrc")
    parser.add_argument("--generated-dir", type=Path, default=ROOT / "data/generated/bvbrc")
    parser.add_argument("--manifest-dir", type=Path, default=ROOT / "data/manifests")
    parser.add_argument("--page-size", type=int, default=1000)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def run(args: argparse.Namespace) -> dict[str, Any]:
    started_at = utc_now()
    config_path = args.config.resolve()
    config = load_config(config_path)
    source = config["source"]
    client = BvbrcClient(source["api_base"], page_size=args.page_size)
    raw_dir = args.raw_dir.resolve()
    generated_dir = args.generated_dir.resolve()
    manifest_dir = args.manifest_dir.resolve()

    all_records: list[dict[str, Any]] = []
    request_log = []
    for drug in config["antibiotics"]:
        filters = [
            rql_eq("evidence", source["evidence"]),
            rql_eq("taxon_id", config["species"]["taxon_id"]),
            rql_eq("antibiotic", drug),
        ]
        records, requests_for_drug = client.fetch_pages(
            source["amr_resource"],
            filters,
            sort=("genome_id", "antibiotic", "id"),
            snapshot_dir=raw_dir / "amr" / drug,
            snapshot_prefix="amr",
            overwrite=args.overwrite,
        )
        all_records.extend(records)
        request_log.extend(requests_for_drug)
        print(f"fetched {drug}: {len(records)} laboratory records")

    clean_labels, conflicts, cleaning_stats = clean_amr_records(all_records, config)
    generated_dir.mkdir(parents=True, exist_ok=True)
    clean_path = generated_dir / "clean_labels_all.csv"
    conflict_path = manifest_dir / "label_conflicts.csv"
    write_csv_atomic(clean_path, clean_labels)
    write_csv_atomic(conflict_path, conflicts)

    genome_ids = clean_labels["genome_id"].astype(str).unique().tolist()
    metadata_records, metadata_requests = fetch_genome_metadata(
        client,
        genome_ids,
        source["genome_resource"],
        raw_dir / "genome_metadata",
        overwrite=args.overwrite,
    )
    request_log.extend(metadata_requests)
    metadata = pd.DataFrame(metadata_records)
    if metadata.empty:
        raise ValueError("BV-BRC returned no genome metadata")
    metadata["genome_id"] = metadata["genome_id"].astype(str)
    metadata_path = generated_dir / "genome_metadata_all.csv"
    write_csv_atomic(metadata_path, metadata.sort_values("genome_id", kind="mergesort"))

    selected_genomes, selected_labels, selection_stats = select_cohort(clean_labels, metadata, config)
    manifest_dir.mkdir(parents=True, exist_ok=True)
    genomes_path = manifest_dir / "selected_genomes.csv"
    labels_path = manifest_dir / "labels.csv"
    download_path = manifest_dir / "download_manifest.csv"
    write_csv_atomic(genomes_path, selected_genomes)
    write_csv_atomic(labels_path, selected_labels)

    download = selected_genomes[
        [
            "genome_id",
            "genome_name",
            "tested_drug_count",
            "genome_length",
            "contigs",
            "checkm_completeness",
            "checkm_contamination",
        ]
    ].copy()
    download["fasta_relative_path"] = download["genome_id"].map(lambda value: f"data/genomes/{value}.fna")
    download["fasta_sha256"] = ""
    write_csv_atomic(download_path, download)

    artifacts = {}
    for path in [clean_path, metadata_path, conflict_path, genomes_path, labels_path, download_path]:
        artifacts[path_for_record(path, ROOT)] = {
            "sha256": sha256_file(path),
            "bytes": path.stat().st_size,
        }

    request_records = request_records_to_dict(request_log)
    for record in request_records:
        if record["snapshot_path"]:
            record["snapshot_path"] = path_for_record(Path(record["snapshot_path"]), ROOT)

    provenance = {
        "schema_version": 1,
        "started_at": started_at,
        "completed_at": utc_now(),
        "status": config.get("status"),
        "config_path": path_for_record(config_path, ROOT),
        "config_sha256": sha256_json(config),
        "git_commit": current_git_commit(ROOT),
        "runtime": runtime_versions(),
        "source": source,
        "species": config["species"],
        "antibiotics": config["antibiotics"],
        "cleaning_stats": cleaning_stats,
        "selection_stats": selection_stats,
        "requests": request_records,
        "artifacts": artifacts,
    }
    provenance_path = manifest_dir / "provenance.json"
    write_json_atomic(provenance_path, provenance)
    print(
        f"selected {selection_stats['selected_genomes']} genomes and "
        f"{selection_stats['selected_label_pairs']} genome/drug labels"
    )
    print(f"wrote provenance: {provenance_path}")
    return provenance


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
