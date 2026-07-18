"""Download selected BV-BRC assemblies through the HTTPS API as FASTA files."""

from __future__ import annotations

import argparse
import concurrent.futures
import os
from pathlib import Path
from typing import Any

import pandas as pd

from data.bvbrc import BvbrcClient, load_config, rql_eq, sha256_file, write_csv_atomic


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "data/config/project.json")
    parser.add_argument("--manifest", type=Path, default=ROOT / "data/manifests/download_manifest.csv")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/genomes")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def format_fasta(records: list[dict[str, Any]]) -> bytes:
    if not records:
        raise ValueError("Genome has no sequence records")
    lines: list[str] = []
    for record in sorted(records, key=lambda value: str(value.get("sequence_id", ""))):
        sequence_id = str(record.get("sequence_id", "")).strip()
        sequence = str(record.get("sequence", "")).replace("\n", "").replace("\r", "").upper()
        if not sequence_id or not sequence:
            raise ValueError("Sequence record is missing sequence_id or sequence")
        invalid = set(sequence) - set("ACGTNRYKMSWBDHV")
        if invalid:
            raise ValueError(f"Sequence {sequence_id} contains invalid FASTA symbols: {sorted(invalid)}")
        description = str(record.get("description", "")).strip()
        header = f">{sequence_id}" + (f" {description}" if description else "")
        lines.append(header)
        lines.extend(sequence[offset : offset + 80] for offset in range(0, len(sequence), 80))
    return ("\n".join(lines) + "\n").encode("ascii")


def download_one(
    genome_id: str,
    output_dir: Path,
    api_base: str,
    sequence_resource: str,
    overwrite: bool,
) -> dict[str, Any]:
    output = output_dir / f"{genome_id}.fna"
    if output.exists() and not overwrite:
        return {
            "genome_id": genome_id,
            "status": "existing",
            "fasta_relative_path": str(output.relative_to(ROOT)),
            "fasta_sha256": sha256_file(output),
            "bytes": output.stat().st_size,
        }
    client = BvbrcClient(api_base, page_size=10000)
    clauses = [
        rql_eq("genome_id", genome_id),
        "select(genome_id,sequence_id,accession,description,sequence)",
        "sort(%2Bsequence_id)",
        "limit(10000,0)",
    ]
    records, _request = client.get_json(sequence_resource, clauses)
    content = format_fasta(records)
    output_dir.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    temporary.write_bytes(content)
    os.replace(temporary, output)
    return {
        "genome_id": genome_id,
        "status": "downloaded",
        "fasta_relative_path": str(output.relative_to(ROOT)),
        "fasta_sha256": sha256_file(output),
        "bytes": output.stat().st_size,
    }


def main() -> None:
    args = parse_args()
    config = load_config(args.config.resolve())
    manifest = pd.read_csv(args.manifest, dtype={"genome_id": str})
    genome_ids = manifest["genome_id"].drop_duplicates().tolist()
    if args.limit is not None:
        genome_ids = genome_ids[: args.limit]
    output_dir = args.output_dir.resolve()
    source = config["source"]
    results: list[dict[str, Any]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                download_one,
                genome_id,
                output_dir,
                source["api_base"],
                source["sequence_resource"],
                args.overwrite,
            ): genome_id
            for genome_id in genome_ids
        }
        for future in concurrent.futures.as_completed(futures):
            genome_id = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                raise RuntimeError(f"Failed to download genome {genome_id}") from exc
            results.append(result)
            print(f"{result['status']}: {genome_id}")
    checksums = pd.DataFrame(results).sort_values("genome_id", kind="mergesort")
    write_csv_atomic(args.manifest.parent / "fasta_checksums.csv", checksums)


if __name__ == "__main__":
    main()
