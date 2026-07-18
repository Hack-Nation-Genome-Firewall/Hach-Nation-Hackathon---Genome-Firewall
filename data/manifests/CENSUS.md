# BV-BRC cohort census

Generated on 2026-07-18 by commit
`e9a931417a0d300579cfcdc9e4a1253fccf5b2d1` using
`data/config/project.json`. The complete machine-readable record is
`provenance.json`.

## Scope

- Species: *Klebsiella pneumoniae* (`taxon_id=573`)
- Evidence: `Laboratory Method` only
- Included labels: `Susceptible`, `Resistant`
- Excluded labels: `Intermediate`, missing/uncategorized phenotype
- Selected cohort cap: 3,000 genomes
- Minimum observed drugs per selected genome: 2

## Filtering audit

| Stage | Count |
|---|---:|
| Laboratory source records retrieved | 20,375 |
| Intermediate records excluded | 475 |
| Missing/uncategorized phenotype records excluded | 4,219 |
| S/R records accepted before duplicate collapse | 15,681 |
| Duplicate agreeing records collapsed | 93 |
| Conflicting genome/drug pairs excluded | 0 |
| Clean unique genome/drug pairs | 15,588 |
| Genomes with metadata | 4,417 |
| Genomes passing configured QC | 4,343 |
| Genomes passing QC with at least two configured drug results | 4,186 |
| Deterministically selected genomes | 3,000 |
| Selected genome/drug labels | 11,900 |

No phenotype was inferred from an MIC value. Records without a categorical S/R
phenotype remain excluded pending an expert-approved breakpoint policy.

## Selected class counts

| Antibiotic | Resistant | Susceptible | Total |
|---|---:|---:|---:|
| Ceftazidime | 2,037 | 926 | 2,963 |
| Ciprofloxacin | 1,918 | 1,043 | 2,961 |
| Gentamicin | 1,049 | 1,945 | 2,994 |
| Meropenem | 807 | 2,175 | 2,982 |

These are labels in the selected cohort, not model results. No accuracy or
clinical-performance claim can be made from this census.

## Reproducibility artifacts

- `labels.csv`: frozen laboratory labels and source record provenance.
- `selected_genomes.csv`: frozen cohort and BV-BRC quality metadata.
- `download_manifest.csv`: genome IDs and intended local FASTA paths.
- `label_conflicts.csv`: auditable excluded conflicts; header-only in this run.
- `provenance.json`: exact query URLs, ETags, response hashes, filtering counts,
  runtime versions, configuration hash, and artifact hashes.

Raw API pages and FASTA assemblies are intentionally not committed. Their
identities are recorded by the provenance and checksum manifests.
