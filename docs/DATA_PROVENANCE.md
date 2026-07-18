# Data provenance and reproducibility

The pipeline separates **source acquisition**, **scientific cleaning**, and
**dataset selection**. Each stage writes an artifact and a SHA-256 checksum.
Changing the source response, scientific policy, or selection code therefore
changes the recorded provenance.

## Source of truth

The machine-readable policy is `data/config/project.json`. A run records:

- retrieval time in UTC;
- exact BV-BRC API URLs and query parameters;
- HTTP ETag and content range where available;
- raw-response SHA-256 checksums;
- counts before and after every filter;
- selected genome IDs and phenotype labels;
- local FASTA SHA-256 checksums;
- configuration SHA-256 checksum;
- pipeline version.

## Label boundary

Only BV-BRC records with `evidence == "Laboratory Method"` may become labels.
Records with computational evidence are rejected even when they contain a
plausible phenotype. Intermediate results are retained in audit counts but not
used by the binary baseline. Conflicting S/R records for the same genome and
antibiotic are excluded rather than silently resolved.

## Reproducing a frozen dataset

The committed selection manifest is the immutable description of the modeling
cohort. A future API query may return additional or corrected records, so a new
retrieval is a new dataset version. It must not silently overwrite an existing
manifest used for reported results.

Genome FASTA files are regenerated from the selected genome IDs and verified
against the recorded checksums. They are not stored in Git.
