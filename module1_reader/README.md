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
a silent contract break), TSV parsing + alias mapping, unknown-marker dropping, the
validation gate's rejections, the precomputed backend, and the spec fallback.

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

## Known TODOs (honest limitations)

- **Target-presence flags** (`target__<gene>`): AMRFinderPlus reports resistance
  markers, not whether a drug's *target* gene is present/intact. These default to
  present (1) until a real gene-presence check (e.g. BLAST the target gene vs the
  assembly) is wired into `AMRFinderPlusAnnotator.annotate()`.
- **`ompK36_loss`** is a derived *absence* feature, not a direct AMRFinderPlus hit.
- **QC columns** (`qc_complete`, `qc_contigs`) need to be filled from assembly stats.
- **Pin the DB**: record `amrfinder -V` into `feature_spec.json`
  (`annotation_tool_version`, `database_version`) at setup and never update mid-event.
