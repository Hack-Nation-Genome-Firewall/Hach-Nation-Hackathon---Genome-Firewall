# Differentiators — temporal validation + stewardship layer

Two additions that raise the submission above a standard "predict resistance"
model, both grounded in real cohort data and neither duplicating the existing
Track B/C safety machinery (calibration, OOD abstention, mechanism evidence,
target/QC gates are already in `module2_predictor/predict.py`).

---

## 1. Prospective (temporal) validation — `data/manifests/split_manifest_temporal.csv`

A homology-grouped split proves the model does not memorise near-identical
strains. A **temporal** split proves something a deployment actually needs: that
a model trained on *today's* isolates still works on *next year's*.

- **Train:** isolates collected **2001–2014** (1,782 genomes)
- **Calibration:** a grouped 15% holdout of the same pre-2015 period (221)
- **Test:** isolates collected **2015–2018** (813 genomes)
- Collection years pulled from BV-BRC (`collection_year` / parsed
  `collection_date`; 94% of the cohort has a usable year).

### Why this matters — resistance genuinely drifts in this cohort
R-rate, pre-2015 vs 2015–2018:

| drug | past R-rate | future R-rate | drift |
|---|---|---|---|
| ceftazidime | 0.75 | 0.54 | **−21 pts** |
| ciprofloxacin | 0.71 | 0.50 | **−22 pts** |
| meropenem | 0.28 | 0.23 | −5 pts |
| gentamicin | 0.35 | 0.34 | −1 pt |

The ceftazidime/ciprofloxacin shifts are large — a random split would hide this
distribution change entirely and report optimistic metrics. Reporting the
temporal-split numbers alongside the grouped-split numbers is an honesty signal:
it shows the real-world generalisation gap instead of concealing it.

### Honest caveat (documented, not hidden)
This is a **pure temporal** cut. 129/813 test genomes (15%) share a cgMLST HC10
cluster with a pre-2015 isolate, so it is not simultaneously a strict homology
split. Two legitimate framings, both reported:
- *temporal* generalisation: this split (train on the past, test on the future);
- *phylogenetic* generalisation: the grouped `split_manifest.csv`.
A combined "future AND novel-lineage" split is stricter still but shrinks the
test set sharply (spanning clusters are large); the two-axis reporting is the
more informative choice for a submission.

Use it by pointing the existing Track B commands at this manifest:
```bash
python -m module2_predictor.train  --splits data/manifests/split_manifest_temporal.csv ...
python -m module2_predictor.evaluate --splits data/manifests/split_manifest_temporal.csv ...
```

---

## 2. Stewardship recommendation layer — `module2_predictor/stewardship.py`

Track B outputs a per-drug antibiogram. A clinician still has to choose an agent.
This layer makes that choice explicit **and stewardship-aware**:

- recommends the **narrowest-spectrum** agent predicted to work
  (gentamicin < ciprofloxacin < ceftazidime < meropenem), sparing last-resort
  carbapenems unless nothing narrower is predicted effective;
- carries the **mechanism + calibrated confidence** into the recommendation, so
  it is explainable;
- **safe by construction**: if every agent is `no_call` or `likely_to_fail`, it
  does not invent a choice — it returns `escalate_to_laboratory_ast`. A `no_call`
  is never overridden.

```python
from module2_predictor.stewardship import recommend_from_genome
out = recommend_from_genome(row, bundle, spec)
# out["stewardship"]["primary_choice"] -> "gentamicin"  (narrowest workable)
#   or                                  -> escalate_to_laboratory_ast
```
Tested in `tests/test_stewardship.py` (narrowest-agent selection; all-fail and
all-no_call both escalate; confidence floor only tightens; safety notice always
present).

### Why it helps the pitch
It connects the model to the actual clinical problem — **antibiotic
stewardship**. Narrow-spectrum-first prescribing is exactly what slows the
resistance this tool detects, so the recommendation logic is on-message with the
"defensive firewall" thesis, not a bolt-on.

> Decision support only. Every recommendation requires confirmation by standard
> laboratory antimicrobial susceptibility testing.

---

## Explicitly out of scope (by design)
- **Longitudinal evolution monitoring / combo-therapy escalation** — needs
  serial isolates per patient/lineage the cohort doesn't have, and the brief
  marks it out of scope. The temporal-drift result above *is* the evolution
  evidence; live surveillance belongs on the roadmap slide, not in the build.
- **De-novo target / sequence generation** — this is a defensive firewall; the
  rubric rewards calibrated abstention, not speculative generation. The
  on-theme version of "help beyond detection" is the stewardship layer above.
- **Multi-tool feature merging** — the K. pneumoniae benchmark shows combining
  annotators *lowers* ML performance; one well-curated tool is the correct choice
  (see the annotator note in the data discussion).
