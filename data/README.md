# Data workflow

BV-BRC is the primary source. Labels come only from `genome_amr` records whose
`evidence` is exactly `Laboratory Method`. Computational phenotype predictions
are never accepted as labels.

## Committed

- `config/project.json`: query, label, quality, and selection policy.
- `manifests/*.csv`: the frozen selected genome/label rows used by the team.
- `manifests/provenance.json`: query URLs, timestamps, counts, and SHA-256
  checksums.
- Download and preparation code.

## Not committed

- `raw/`: exact API responses, retained locally and checksum-recorded.
- `genomes/`: assembled FASTA files.
- `generated/`: intermediate metadata and reports.

Generated content is excluded because genome collections are large. The frozen
manifest and checksums identify every input without putting gigabytes in Git.

## Scientific policy

- Include `Susceptible` and `Resistant` laboratory records.
- Exclude `Intermediate` from the binary baseline.
- Collapse repeated records that agree.
- Exclude genome/drug pairs with conflicting S/R records.
- Preserve MIC, testing-standard, method, publication, and source metadata.
- Apply genome-quality criteria before selection.
- Prefer genomes with results for more configured drugs, then use a seeded
  SHA-256 rank to cap the cohort deterministically.
- Split by a separately computed homology cluster; never randomly split rows.

The current policy is provisional until approved by the domain expert.
