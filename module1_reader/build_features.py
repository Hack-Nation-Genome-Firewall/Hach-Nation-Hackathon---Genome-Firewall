"""
GENOME FIREWALL — Module 1 (Track A) public entry points.

Two ways in, both returning rows that already satisfy the frozen feature contract:

  run_genome_reader(fasta)      one genome -> validated feature row   (used live by the app)
  build_features_table(...)     many genomes -> features.csv          (used to build training data)

Both are backend-agnostic: pass backend="amrfinderplus" (default) to run the tool,
backend="precomputed" to load an existing table, or your own registered backend.
See feature_annotator.py for the extension point and module1_reader/README.md for
how to add a source.

CLI:
    # single genome via AMRFinderPlus
    python module1_reader/build_features.py --fasta genome.fasta --genome-id G1

    # single genome from a saved AMRFinderPlus TSV (no tool needed)
    python module1_reader/build_features.py --genome-id G1 \
        --backend amrfinderplus --tsv module1_reader/fixtures/sample_amrfinder.tsv

    # batch, bringing your own precomputed feature table
    python module1_reader/build_features.py --backend precomputed \
        --table my_features.csv --out data/manifests/features.csv
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Optional

from feature_annotator import (
    MODULE_DIR, get_annotator, load_spec, marker_columns, target_columns,
    quality_columns, ContractError,
)


def run_genome_reader(fasta: Optional[str] = None, *, genome_id: str = "query",
                      backend: str = "amrfinderplus", spec: Optional[dict] = None,
                      tsv_override: Optional[str] = None, **backend_kwargs) -> dict:
    """
    One genome -> one validated feature row (dict). The function the app calls at
    inference. Guarantees the same column order the models were trained on.
    """
    spec = spec or load_spec()
    annotator = get_annotator(backend, spec, **backend_kwargs)
    # AMRFinderPlus backend accepts a saved TSV so we can run without the tool.
    if tsv_override is not None and backend == "amrfinderplus":
        from feature_annotator import validate_feature_row
        raw = annotator.annotate(genome_id, fasta, tsv_override=Path(tsv_override))
        return validate_feature_row(raw, spec)
    return annotator.annotate_validated(genome_id, fasta)


def build_features_table(genomes: list[dict], *, backend: str = "amrfinderplus",
                         out_path: Optional[Path] = None, spec: Optional[dict] = None,
                         **backend_kwargs) -> Path:
    """
    Build features.csv from many genomes.

    `genomes` is a list of {"genome_id": ..., "source": <fasta path or table id>}.
    Writes genome_id + model_features + target columns + QC columns in contract order.
    """
    spec = spec or load_spec()
    annotator = get_annotator(backend, spec, **backend_kwargs)
    # Default output stays INSIDE this module. At integration, pass out_path=
    # "data/manifests/features.csv" so Track B/C consume it from the shared location.
    out_path = Path(out_path) if out_path else MODULE_DIR / "out/features.csv"

    header = (["genome_id"] + marker_columns(spec)
              + target_columns(spec) + quality_columns(spec))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        for g in genomes:
            row = annotator.annotate_validated(g["genome_id"], g.get("source"))
            row["genome_id"] = g["genome_id"]
            writer.writerow(row)
    return out_path


def _main() -> None:
    p = argparse.ArgumentParser(description="Genome Firewall — Module 1 feature builder")
    p.add_argument("--fasta", help="path to a nucleotide FASTA (amrfinderplus backend)")
    p.add_argument("--genome-id", default="query")
    p.add_argument("--backend", default="amrfinderplus")
    p.add_argument("--tsv", help="use a saved AMRFinderPlus TSV instead of running the tool")
    p.add_argument("--table", help="precomputed feature table (precomputed backend)")
    p.add_argument("--organism", help="AMRFinderPlus --organism, e.g. Klebsiella_pneumoniae")
    p.add_argument("--out", help="write a features.csv here instead of printing one row")
    args = p.parse_args()

    spec = load_spec()
    kwargs = {}
    if args.backend == "precomputed":
        if not args.table:
            p.error("--table is required for the precomputed backend")
        kwargs["table_path"] = args.table
    if args.backend == "amrfinderplus" and args.organism:
        kwargs["organism"] = args.organism

    try:
        row = run_genome_reader(
            args.fasta, genome_id=args.genome_id, backend=args.backend,
            spec=spec, tsv_override=args.tsv, **kwargs,
        )
    except ContractError as e:
        p.error(str(e))

    present = [c for c in marker_columns(spec) if row.get(c) == 1]
    tgt_absent = [c for c in target_columns(spec) if row.get(c) == 0]
    qf = spec["quality_features"]
    print(f"genome_id: {args.genome_id}  (backend={args.backend})")
    print(f"markers present ({len(present)}): {present}")
    print(f"targets absent/disrupted: {tgt_absent or 'none'}")
    print(f"QC: completeness={row.get(qf['completeness'])} "
          f"contamination={row.get(qf['contamination'])} contigs={row.get(qf['contigs'])}")


if __name__ == "__main__":
    _main()
