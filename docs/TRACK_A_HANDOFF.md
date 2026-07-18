# Track A to Track B handoff

Track B is runnable against `data/synthetic/`, but real scientific training is
intentionally blocked until Track A provides AMRFinderPlus-derived features and
a homology-grouped split. Track B will not invent missing features or assume
missing targets are present.

## Required files

### `features.csv`

One row per genome with:

- `genome_id` matching `data/manifests/selected_genomes.csv`;
- every binary column in `feature_spec.json:model_features`;
- every target-status column named under `drug_targets`;
- the three QC columns mapped by `quality_features`.

Marker features must be binary and complete. Unknown/missing marker columns are
a schema error, not absence. A missing target or QC value produces no-call at
inference.

### `split_manifest.csv`

Exactly one row per feature genome:

```text
genome_id,cluster_id,split
573.x,cluster_001,train
```

Allowed splits are exactly `train`, `calibration`, and `test`. One homology
cluster may occur in only one split. Track B asserts this before fitting.

### `feature_spec.json`

Use `data/synthetic/feature_spec.json` as the schema example, then replace every
synthetic marker and target with reviewed real definitions. Required sections:

- supported species and drugs;
- ordered `model_features`;
- `marker_evidence` with marker type, applicable drugs, source, and
  AMRFinderPlus detection method;
- `drug_targets`, expressed as input feature names;
- QC feature mapping and thresholds;
- `expected_label_evidence` set to `Laboratory Method` for the real cohort.

The domain expert must approve target definitions and marker-to-drug mappings.
AMRFinderPlus does not by itself prove that every ordinary molecular target is
present and intact, so Track A must document the separate target-detection rule.

## Real training command

```bash
python -m module2_predictor.train \
  --features path/to/features.csv \
  --labels data/manifests/labels.csv \
  --splits path/to/split_manifest.csv \
  --spec path/to/feature_spec.json \
  --output models/kpneumoniae_v1.joblib \
  --model-version kpneumoniae-v1
```

Then evaluate using the same four input contracts:

```bash
python -m module2_predictor.evaluate \
  --features path/to/features.csv \
  --labels data/manifests/labels.csv \
  --splits path/to/split_manifest.csv \
  --spec path/to/feature_spec.json \
  --bundle models/kpneumoniae_v1.joblib
```

Do not use the organizer's hidden labels during development. The local `test`
split is a frozen grouped holdout; the organizer's hidden set remains external.
