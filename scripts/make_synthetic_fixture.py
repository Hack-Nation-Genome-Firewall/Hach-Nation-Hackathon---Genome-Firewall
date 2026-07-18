"""Generate deterministic contract fixtures for Track B integration tests.

The output is deliberately synthetic and must never be presented as biological
or clinical model performance.
"""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from data.bvbrc import write_csv_atomic, write_json_atomic  # noqa: E402


OUTPUT = ROOT / "data/synthetic"
DRUGS = ["meropenem", "ciprofloxacin", "gentamicin", "ceftazidime"]
KNOWN_MARKERS = {drug: f"marker__known__{drug}" for drug in DRUGS}
MODEL_FEATURES = [*KNOWN_MARKERS.values(), "marker__lineage_a", "marker__lineage_b", "marker__background"]


def build() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    feature_rows = []
    label_rows = []
    split_rows = []
    for cluster_number in range(60):
        cluster_id = f"cluster_{cluster_number:03d}"
        if cluster_number < 40:
            split = "train"
        elif cluster_number < 50:
            split = "calibration"
        else:
            split = "test"
        cluster_markers = {
            marker: int((cluster_number + offset) % 5 in {0, 1})
            for offset, marker in enumerate(KNOWN_MARKERS.values())
        }
        cluster_markers["marker__lineage_a"] = int(cluster_number % 3 == 0)
        cluster_markers["marker__lineage_b"] = int(cluster_number % 4 == 0)

        for replicate in range(4):
            genome_id = f"synthetic_{cluster_number:03d}_{replicate}"
            features = {
                "genome_id": genome_id,
                **cluster_markers,
                "marker__background": int((cluster_number + replicate) % 7 == 0),
                **{f"target__{drug}": 1 for drug in DRUGS},
                "qc_completeness": 98.0,
                "qc_contamination": 1.0,
                "qc_contigs": 30,
            }
            if cluster_number == 50 and replicate == 0:
                features["target__meropenem"] = 0
            if cluster_number == 51 and replicate == 0:
                features["qc_completeness"] = 70.0
                features["qc_contigs"] = 900
            if cluster_number == 59:
                for marker in MODEL_FEATURES:
                    features[marker] = 1
            feature_rows.append(features)
            split_rows.append({"genome_id": genome_id, "cluster_id": cluster_id, "split": split})

            for offset, drug in enumerate(DRUGS):
                known = bool(cluster_markers[KNOWN_MARKERS[drug]])
                statistical = bool(cluster_markers["marker__lineage_a"] and (cluster_number + offset) % 7 == 0)
                resistant = known or statistical
                if (cluster_number * 7 + offset) % 23 == 0:
                    resistant = not resistant
                label_rows.append(
                    {
                        "genome_id": genome_id,
                        "antibiotic": drug,
                        "phenotype": "Resistant" if resistant else "Susceptible",
                        "evidence": "Synthetic fixture",
                    }
                )

    marker_evidence = {
        marker: {
            "type": "synthetic_marker",
            "drugs": [drug],
            "source": "Synthetic integration fixture; not biological evidence",
            "amrfinder_method": "SYNTHETIC",
        }
        for drug, marker in KNOWN_MARKERS.items()
    }
    spec = {
        "schema_version": 1,
        "synthetic": True,
        "status": "integration_fixture_not_biological_data",
        "species": {"name": "Klebsiella pneumoniae", "taxon_id": 573},
        "drugs": DRUGS,
        "model_features": MODEL_FEATURES,
        "marker_evidence": marker_evidence,
        "drug_targets": {drug: [f"target__{drug}"] for drug in DRUGS},
        "quality_features": {
            "completeness": "qc_completeness",
            "contamination": "qc_contamination",
            "contigs": "qc_contigs",
        },
        "quality_policy": {
            "minimum_completeness": 90.0,
            "maximum_contamination": 5.0,
            "maximum_contigs": 500,
        },
        "expected_label_evidence": "Synthetic fixture",
    }
    return (
        pd.DataFrame(feature_rows).sort_values("genome_id", kind="mergesort"),
        pd.DataFrame(label_rows).sort_values(["genome_id", "antibiotic"], kind="mergesort"),
        pd.DataFrame(split_rows).sort_values("genome_id", kind="mergesort"),
        spec,
    )


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    features, labels, splits, spec = build()
    write_csv_atomic(OUTPUT / "features.csv", features)
    write_csv_atomic(OUTPUT / "labels.csv", labels)
    write_csv_atomic(OUTPUT / "split_manifest.csv", splits)
    write_json_atomic(OUTPUT / "feature_spec.json", spec)
    print(f"wrote {len(features)} genomes and {len(labels)} labels to {OUTPUT}")


if __name__ == "__main__":
    main()
