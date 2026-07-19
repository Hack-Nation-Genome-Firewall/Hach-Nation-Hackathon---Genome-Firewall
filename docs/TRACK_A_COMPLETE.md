# Track A — real feature layer, homology split, integration (complete)

This closes the Track A gap flagged in `TEAM_STATUS_AND_NEXT_STEPS.md`:
real AMRFinderPlus features, a homology-grouped split, and validated hand-off
into Track B. Everything here is contract-checked against
`module2_predictor/contracts.py`.

## What is in this branch

### 1. Homology-grouped split — `data/manifests/split_manifest.csv`
- **3,000 genomes** (exactly the `selected_genomes.csv` cohort), **1,715 clusters**
- Split **2,100 train / 360 calibration / 540 test**
- **Zero clusters span splits** (Track B's cross-split assertion passes)
- Group key: **BV-BRC cgMLST HC10** (genomes within 10 core-genome allele
  differences), with `mlst` then singleton fallback for the 17 genomes lacking
  an HC10 id.
- R-rate preserved across splits for all four drugs (e.g. ceftazidime
  0.69/0.65/0.69 train/cal/test; meropenem 0.28/0.26/0.24).

> **Why cgMLST instead of Mash/Sourmash.** The domain expert approved cgMLST
> HC10: it is a within-species, allele-level clustering already computed by
> BV-BRC and is at least as strict as a Mash ANI threshold for this use. The
> no-leakage invariant (no cluster in two splits) is enforced in code, so the
> guarantee Track B relies on is identical either way.

### 2. Feature assembler — `module1_reader/assemble_features.py`
Turns a directory of AMRFinderPlus TSVs into the exact Track B contract:
- `features.csv` — `genome_id` + binary `marker__*` + `target__*` + `qc_*`
- `feature_spec.json` — full `schema_version / species / drugs / model_features
  / marker_evidence / drug_targets / quality_features / quality_policy /
  expected_label_evidence` shape (passes `validate_feature_spec`).

The marker vocabulary is **learned empirically** from real AMRFinderPlus output
(family-collapsed; 3–97% prevalence band drops chromosomal constants like
`blaSHV` and ultra-rare alleles). On the 87-genome pilot it yields 20 markers:
CTX-M/KPC/NDM/OXA/TEM β-lactamases, qnrB/qnrS, gyrA/parC point mutations,
ompK35/ompK36 porin loss, aac/aph/aad/ant/rmt aminoglycoside genes. The final
vocabulary is regenerated from the full 3,000-genome run.

### 3. Integration proven
On the 87-genome pilot, the assembler output passes, in order:
`validate_feature_spec` → `validate_training_frames` → into
`module2_predictor.train`'s model-fit and calibration stages. The only thing
that stops the pilot short of a saved model is the per-class / calibration
sample-size floors (≥5 per class, ≥10 for calibration) — a data-volume limit
that the full run clears, **not** a schema mismatch.

## Remaining step — the 3,000-genome AMRFinderPlus run (needs a real machine)

This cannot run in the constrained sandbox (~10 GB disk; the run needs ~15–18 GB
for streamed assemblies + annotation). Run it on the cluster:

```bash
# 1. env (once)
conda create -y -n amrfinder -c conda-forge -c bioconda ncbi-amrfinderplus
conda activate amrfinder
amrfinder -u                              # DB (~250 MB)

# 2. annotate all selected genomes (resumable; streams FASTA download->run->delete)
python module1_reader/run_cohort_parallel.py \
    --ids  data/manifests/selected_genomes.csv \   # or a genome_id json list
    --db   $CONDA_PREFIX/share/amrfinderplus/data/latest \
    --amrfinder $(which amrfinder) \
    --outdir cohort_tsv --workers 16 --amr-threads 2

# 3. assemble the real feature matrix + spec (Track B contract)
python -m module1_reader.assemble_features \
    --tsv-dir cohort_tsv \
    --selected data/manifests/selected_genomes.csv \
    --out-features data/manifests/features.csv \
    --out-spec     data/manifests/feature_spec.json

# 4. real training + evaluation (Track B/C, unchanged)
python -m module2_predictor.train \
    --features data/manifests/features.csv \
    --labels   data/manifests/labels.csv \
    --splits   data/manifests/split_manifest.csv \
    --spec     data/manifests/feature_spec.json \
    --output   models/kp_bundle.joblib --model-version kp-real-v1
python -m module2_predictor.evaluate ...        # per-drug metrics + reliability
```
Runtime: ~50 s/genome on fragmented K. pneumoniae assemblies (many 200–450
contigs); ~2.5–4 h at 16 workers for 3,000. `run_cohort_parallel.py` skips
genomes already annotated, so the run is safe to resume.

## Domain-expert review still required
- `feature_spec.json:marker_evidence[*].drugs` — a few β-lactamase families map
  to both meropenem and ceftazidime by AMRFinderPlus class; confirm the
  carbapenemase vs ESBL split per family.
- `drug_targets` — confirm the essential-target gene list per drug and the
  POINT_DISRUPT gate rule.
- `quality_policy` thresholds (completeness ≥90, contamination ≤5, contigs ≤500).
