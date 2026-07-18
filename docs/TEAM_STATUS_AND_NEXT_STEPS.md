# Genome Firewall: Team Status and Next Steps

**Status snapshot:** 18 July 2026 (Reykjavik)

**Repository:** https://github.com/liiandy/Hach-Nation-Hackathon---Genome-Firewall

**Current main commit:** `c4b9aa4`

## Executive summary

We have a reproducible BV-BRC laboratory-label cohort and a tested Track B
software pipeline. We do **not** yet have real AMRFinderPlus features, a
homology-grouped split, or scientifically meaningful model results.

The critical path is now Track A:

```text
real FASTA
  -> AMRFinderPlus and target/QC features
  -> homology clustering and grouped split
  -> real Track B training and calibration
  -> Streamlit decision report
```

Do not report the current synthetic metrics as biological performance.

## Organizer clarification

The challenge discussion confirmed:

- The organizers are not providing the fixed, cleaned dataset described in the
  appendix.
- Teams may use any dataset that fits the prompt.
- We must build and document our own cohort.

This makes data provenance, cleaning, grouped evaluation, and honest scope part
of the submission value. It also means cross-team headline accuracy comparisons
will not be directly fair because teams may use different cohorts.

## Provisional scientific scope

- Species: *Klebsiella pneumoniae* (`taxon_id=573`)
- Antibiotics:
  - meropenem
  - ciprofloxacin
  - gentamicin
  - ceftazidime
- Input boundary: assembled, quality-checked bacterial FASTA only
- Output: `likely_to_work`, `likely_to_fail`, or `no_call`
- Safety: research decision support only; every result requires standard
  laboratory confirmation

The domain expert still needs to formally approve the species, drugs, label
policy, target rules, marker mappings, QC thresholds, and homology threshold.

## What is complete

### 1. Reproducible BV-BRC acquisition

The pipeline:

- queries the public BV-BRC HTTPS API;
- accepts only `evidence == "Laboratory Method"`;
- includes categorical `Susceptible` and `Resistant` labels;
- audits and excludes `Intermediate` and phenotype-free/MIC-only records;
- collapses duplicate records that agree;
- excludes conflicting S/R pairs instead of silently resolving them;
- retrieves genome QC metadata;
- creates a deterministic cohort;
- records exact URLs, ETags, timestamps, hashes, filtering counts, runtime
  versions, and artifact checksums.

Frozen cohort summary:

| Item | Count |
|---|---:|
| Laboratory source records examined | 20,375 |
| Intermediate records excluded | 475 |
| Missing/uncategorized phenotype records excluded | 4,219 |
| Agreeing duplicate records collapsed | 93 |
| Conflicting S/R genome/drug pairs | 0 |
| Clean unique genome/drug pairs | 15,588 |
| Genomes with metadata | 4,417 |
| Genomes passing configured QC | 4,343 |
| Selected genomes | 3,000 |
| Selected genome/drug labels | 11,900 |

Selected class counts:

| Antibiotic | Resistant | Susceptible |
|---|---:|---:|
| Ceftazidime | 2,037 | 926 |
| Ciprofloxacin | 1,918 | 1,043 |
| Gentamicin | 1,049 | 1,945 |
| Meropenem | 807 | 2,175 |

Important files:

- `data/config/project.json`
- `data/manifests/CENSUS.md`
- `data/manifests/labels.csv`
- `data/manifests/selected_genomes.csv`
- `data/manifests/download_manifest.csv`
- `data/manifests/provenance.json`

Raw API pages and genome FASTA files are intentionally excluded from Git.

### 2. Track B software

Track B now supports:

- one regularized logistic-regression model per drug;
- strict feature and label contract validation;
- an assertion that no homology cluster crosses data splits;
- dedicated probability calibration using sigmoid or isotonic calibration;
- per-drug no-call threshold selection;
- deterministic target-status gating;
- assembly-QC no-calls;
- out-of-distribution checks;
- known-marker/model-conflict no-calls;
- conventional calibrated `P(fail)` and predicted-verdict probability;
- honest known-marker, statistical-only, and no-signal evidence categories;
- model version and feature-spec checksum in every prediction;
- overall and per-cluster held-out evaluation;
- Brier, balanced accuracy, R/S recall, F1, AUROC, PR-AUC, no-call rate,
  coverage, and accuracy on called cases.

### 3. Synthetic integration fixture

The deterministic fixture proves that the contracts and code run end to end:

```bash
python scripts/make_synthetic_fixture.py
python -m module2_predictor.train
python -m module2_predictor.predict
python -m module2_predictor.evaluate
```

It contains synthetic markers and labels. Its metrics are software test results,
not biological evidence.

### 4. Track C compatibility

The Streamlit app and Track C evaluation entry point now consume the current
Track B bundle and record schema. The app visibly labels synthetic mode and
retains the mandatory laboratory-confirmation warning.

Automated status:

```text
12 tests passed
Synthetic generate -> train -> predict -> evaluate passed
Streamlit application test executed without exceptions
```

## What is not complete

1. The 3,000 selected FASTA assemblies have not all been downloaded.
2. AMRFinderPlus has not been installed, version-pinned, and run on the cohort.
3. Real resistance marker features do not exist yet.
4. Molecular target presence/intactness rules are not finalized.
5. Mash/Sourmash homology clusters and grouped splits do not exist yet.
6. Track B has not been trained on real biological features.
7. There are no real held-out performance or calibration results.
8. The app does not yet process a real uploaded FASTA through Track A.
9. Attribution, final model card, demo narrative, and final scientific review
   remain incomplete.

## Next milestone

The next milestone is:

> Five real selected FASTAs successfully downloaded, annotated by a pinned
> AMRFinderPlus installation, and converted into feature rows matching the
> Track B contract.

Start with:

```bash
git switch main
git pull
python -m pytest
python -m data.download_genomes --limit 5
```

Verify the generated `data/manifests/fasta_checksums.csv`, inspect all five
FASTA files, and run the same AMRFinderPlus command on each.

Do not immediately download all 3,000 genomes. A typical assembly is roughly
5-6 MB, so the full cohort may require around 15-18 GB before annotations and
intermediate output. Measure download time, disk usage, and annotation runtime
on the five-genome pilot first. If necessary, create a documented cohort v2
capped at 1,000-1,500 genomes.

## Track A deliverables

Track A must provide:

```text
features.csv
feature_spec.json
split_manifest.csv
```

Detailed requirements are in `docs/TRACK_A_HANDOFF.md`.

Track A tasks:

1. Pin AMRFinderPlus binary and database versions.
2. Run FASTA -> AMRFinderPlus TSV reproducibly.
3. Construct binary resistance gene/mutation features.
4. Preserve unknown markers rather than silently dropping them.
5. Implement and document target-presence/intactness features separately.
6. Include QC features using the agreed units.
7. Sketch genomes using Mash or Sourmash.
8. Cluster related genomes using an expert-approved threshold.
9. Assign complete clusters to train, calibration, and test.
10. Prove no cluster crosses splits.

## Real Track B commands

After Track A delivers its contracts:

```bash
python -m module2_predictor.train \
  --features path/to/features.csv \
  --labels data/manifests/labels.csv \
  --splits path/to/split_manifest.csv \
  --spec path/to/feature_spec.json \
  --output models/kpneumoniae_v1.joblib \
  --model-version kpneumoniae-v1
```

```bash
python -m module2_predictor.evaluate \
  --features path/to/features.csv \
  --labels data/manifests/labels.csv \
  --splits path/to/split_manifest.csv \
  --spec path/to/feature_spec.json \
  --bundle models/kpneumoniae_v1.joblib
```

Thresholds may be tuned on the calibration split only. Do not inspect or tune
against organizer-hidden labels if a hidden evaluation is later provided.

## Suggested ownership

| Owner | Immediate responsibility |
|---|---|
| Domain expert | Approve scope, labels, marker/drug mappings, target rules, QC, and homology threshold |
| Track A owner | FASTA download, AMRFinderPlus, real features, Mash/Sourmash clustering |
| Track B owner | Real training, calibration, no-call tuning, grouped evaluation |
| Track C owner | Real upload workflow, report UI, plots, demo reliability |
| Shared | Safety review, attribution, model card, README, presentation, backup demo |

## Guardrails

- Use one species unless the team explicitly accepts the added scientific and
  evaluation complexity of multiple species.
- Use laboratory labels only.
- Never infer S/R from MIC without an expert-approved breakpoint policy.
- Never use a random row split.
- Never treat missing features as absence.
- Never treat missing target information as target present.
- Never present SHAP or model coefficients as biological causation.
- Never present synthetic metrics as real model performance.
- Never claim clinical readiness.
- Do not start a genomic-language-model experiment until the full baseline works
  on real held-out data.

## Definition of done

- Real assembled FASTA -> real annotation -> real features works reproducibly.
- Homology-grouped train/calibration/test manifests are frozen.
- Real calibrated per-drug models are evaluated on held-out clusters.
- No-call coverage and accuracy-on-called are reported.
- Known biology and statistical associations are displayed separately.
- The app processes real supported input and shows mandatory safety language.
- Data, tool, database, and literature attribution is complete.
- The expert signs off on every scientific claim shown to judges.
