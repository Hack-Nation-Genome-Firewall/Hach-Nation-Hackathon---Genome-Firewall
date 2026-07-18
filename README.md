# Genome Firewall

Research prototype for predicting antibiotic resistance from an assembled
bacterial genome. The project is strictly defensive decision support and is not
for clinical use. Every result must be confirmed with standard laboratory
testing.

## Current scope

The provisional data scope is **Klebsiella pneumoniae** (`taxon_id=573`) and:

- meropenem
- ciprofloxacin
- gentamicin
- ceftazidime

This scope is configuration, not a claim of validated coverage. The domain
expert must approve the final species, drugs, marker mappings, target rules,
quality thresholds, and label policy before submission.

## Repository status

Work is organized around three contracts:

1. `data/config/project.json` records the source query and scientific policy.
2. `data/manifests/feature_spec.json` fixes model feature names and evidence
   metadata.
3. Track B emits versioned JSON-compatible prediction records consumed by the
   app and evaluator.

Raw FASTA files, API response snapshots, and fitted models are generated
artifacts and are not committed. Frozen selection manifests, checksums,
provenance, code, and tests are committed so another team can reproduce the
same process.

## Development setup

Python 3.11 is the supported runtime.

```bash
python -m venv .venv
.venv/Scripts/activate       # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python -m pytest
```

Build a frozen BV-BRC laboratory-label cohort:

```bash
python -m data.fetch_bvbrc
```

This writes committed manifests under `data/manifests/` and keeps large raw
responses under ignored `data/raw/`. To download selected assembled genomes
through the BV-BRC HTTPS API:

```bash
python -m data.download_genomes --workers 4
```

Use `--limit 5` for a small end-to-end download check. Acquisition refuses to
overwrite existing raw snapshots unless `--overwrite` is explicitly supplied.

## Safety boundary

The system accepts an already assembled, quality-checked genome. It does not
collect samples, process raw sequencing reads, identify species, make treatment
decisions, or design or modify organisms.
