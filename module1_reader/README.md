# Module 1 — The Genome Reader (Track A)

Turns one genome (a FASTA file) into the fixed feature row the predictor was
trained on. Deterministic, not machine-learned.

```
FASTA file  ->  [ annotator backend ]  ->  marker list  ->  validated feature row
                     ^ swappable                              ^ frozen contract
```

The **only** thing Tracks B and C depend on is the feature contract — the exact
columns and their order. Any tool or dataset that can produce that contract plugs
in here — **AMRFinderPlus is the default, not a hard requirement.**

**Where the contract lives (ownership):** the real `feature_spec.json` is *Phase 0's*
shared deliverable (the domain expert's), and belongs at
`data/manifests/feature_spec.json` — it is **not** Track A's to own. Track A only
*reads* it. So this track doesn't create that shared file. Until Phase 0 delivers
it, Track A falls back to a bundled dev stand-in, `fixtures/feature_spec.sample.json`,
purely so this module runs standalone. `load_spec()` prefers the shared file when it
exists and uses the sample otherwise.

**Connection to the shared project config:** this module also reads (read-only)
Phase 0's `data/config/project.json` — the declared species, drug panel, and QC
policy. `spec_project_discrepancies()` cross-checks the spec against it, and a test
fails if our sample drifts from the team's declared species/drugs. The sample is kept
aligned with that config (currently 4 drugs: meropenem, ciprofloxacin, gentamicin,
ceftazidime). We never write that file — Phase 0 owns it.

## Run it

```bash
# single genome via AMRFinderPlus (needs the tool installed)
python module1_reader/build_features.py --fasta genome.fasta --genome-id G1 \
    --organism Klebsiella_pneumoniae

# single genome from a saved AMRFinderPlus TSV — no tool needed (uses the fixture)
python module1_reader/build_features.py --genome-id G1 \
    --backend amrfinderplus --tsv module1_reader/fixtures/sample_amrfinder.tsv

# bring your own precomputed feature table -> features.csv
# (default output stays inside module1_reader/out/; point --out at data/manifests/ only at integration)
python module1_reader/build_features.py --backend precomputed \
    --table my_features.csv --out module1_reader/out/features.csv
```

At inference the app calls one function:

```python
from module1_reader.build_features import run_genome_reader
row = run_genome_reader("genome.fasta", genome_id="G1")   # -> validated feature dict
```

## Tests

Fast, deterministic, no AMRFinderPlus needed (they use the bundled fixtures):

```bash
pip install pytest
pytest module1_reader/tests
```

They lock down the fragile parts: the exact feature-column order (the tripwire for
a silent contract break), TSV parsing + alias mapping, unknown-marker preservation,
the validation gate's rejections, the precomputed backend, the spec fallback, and
compatibility with Track B's own contract validators.

## Unknown markers are preserved, not dropped

The model vector is fixed-shape, so a marker AMRFinderPlus reports that isn't in the
vocabulary can't become a model column — but per the Track A guardrail we never
silently drop it. Each is recorded: surfaced via `run_genome_reader(...,
unknown_markers_out=[])`, printed by the CLI, and written to an `unknown_markers.csv`
sidecar next to `features.csv` in batch runs (`genome_id,unknown_marker`). That keeps
them available for the domain expert to review and possibly add to the vocabulary.

## Adding your own annotation source (the door is already open)

You do **not** edit the pipeline. Three steps:

1. **Subclass `FeatureAnnotator`** in `feature_annotator.py` and implement
   `annotate(genome_id, source) -> dict`. Return a dict of column → value; markers
   you don't set default to 0. Copy **`PrecomputedAnnotator`** — it's written to be
   the template.
2. **Register it** in the `ANNOTATORS` dict at the bottom of `feature_annotator.py`,
   e.g. `"resfinder": ResFinderAnnotator`.
3. **Select it** with `backend="resfinder"` (in code) or `--backend resfinder` (CLI).

Whatever your backend emits passes through **`validate_feature_row()`**, which
checks it against `feature_spec.json` and fails loudly if columns are missing or a
flag isn't 0/1. So bring-your-own is safe — a wrong table is rejected at the door,
never silently used.

Two "bring your own" cases, both already covered:
- **your own tool** (still start from FASTA) → new `FeatureAnnotator` subclass.
- **your own dataset / features** (skip annotation) → the `precomputed` backend.

## Homology-grouped split manifest (`split_manifest.py`)

Produces `split_manifest.csv` (`genome_id,cluster_id,split`) — required by Track B,
which asserts no homology cluster crosses splits. Two steps:

1. **cluster genomes by DNA similarity** → `cluster_id` (Mash/sourmash) — *pluggable*,
   `MashClusterer` is a documented stub; use precomputed clusters via `load_clusters_csv()`.
2. **assign whole clusters to train/calibration/test** — *done*: deterministic (seeded),
   greedy largest-cluster-first, whole clusters only, so no cluster can straddle splits.

```python
from split_manifest import build_split_manifest, load_clusters_csv
clusters = load_clusters_csv("clusters.csv")   # genome_id,cluster_id
build_split_manifest(clusters, "split_manifest.csv", seed=20260718,
                     feature_genome_ids=feature_ids)
```

A test builds a manifest and runs it through Track B's own
`validate_training_frames` (including the no-leak assertion), so the output is
provably contract-valid.

## Known TODOs (honest limitations)

- **Mash clustering** (step 1 above): needs Mash + the genome FASTAs (not in the repo).
  Until then, feed precomputed clusters. The threshold is an expert-approved knob.
- **Target-presence flags** (`target__<gene>`): AMRFinderPlus reports resistance
  markers, not whether a drug's *target* gene is present/intact. These default to
  present (1) until a real gene-presence check (e.g. BLAST the target gene vs the
  assembly) is wired into `AMRFinderPlusAnnotator.annotate()`.
- **`ompK36_loss`** is a derived *absence* feature, not a direct AMRFinderPlus hit.
- **QC columns** (`qc_completeness`, `qc_contamination`, `qc_contigs`): *done* — read
  straight from `data/manifests/selected_genomes.csv` (CheckM completeness/contamination
  + contig count) via `load_qc_map()` / the `--selected-genomes` flag / `qc_source=`.
  When a genome's QC is unknown, the columns are left **None**, not faked clean, so
  Track B routes it to no-call (guardrail: a missing QC value must not read as good).
  A test checks this against Track B's own `_quality_gate` (pass / fail / unknown).
- **Pin the DB**: record `amrfinder -V` into `feature_spec.json`
  (`annotation_tool_version`, `database_version`) at setup and never update mid-event.
