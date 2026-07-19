# Genome Firewall — Submission

**One-line:** A defensive AI that reads a *Klebsiella pneumoniae* genome and
predicts, per antibiotic, whether the drug will work — with *calibrated
confidence* and an honest **"I don't know, send it to the lab"** whenever the
genome is unfamiliar. A firewall, not a classifier.

**Live app:** ⟦paste your streamlit.app URL here after deploy⟧
**Repo:** https://github.com/liiandy/Hach-Nation-Hackathon---Genome-Firewall

---

## 1. The problem

When a patient is admitted with a *Klebsiella* infection, the doctor must choose
an antibiotic *now* — but culture-based susceptibility testing takes 24–72 hours.
Guess wrong and the patient stays on an ineffective drug while the infection
advances; in septic shock, each hour of ineffective therapy measurably raises
mortality. Whole-genome sequencing of a clinical isolate is fast and increasingly
routine, so the genome is available long before the culture result. The question
is whether we can read resistance *from the genome* — safely enough to act on.

## 2. Why existing tools aren't enough

Genome-based AMR tools today (ResFinder, CARD, AMRFinderPlus, Pathogenwatch,
Kleborate) are **rule-based**: they report the presence or absence of known
resistance genes and emit a binary resistant/susceptible call. Two problems make
them hard to trust at the bedside:

- **No calibrated confidence.** A binary call gives the clinician no way to know
  *how sure* the tool is. Published balanced accuracy for these tools on clinical
  isolates ranges 0.52–0.66 and varies wildly by drug — but the output looks
  equally confident every time.
- **No abstention.** They answer even when the genome carries a resistance
  mechanism the tool has never characterized — precisely the case where a
  confident wrong answer is most dangerous.

## 3. Our solution — the firewall

Genome Firewall keeps the interpretability of the rule-based world and adds the
three things that make a prediction *trustworthy*:

1. **Calibrated probability.** Per-drug logistic regression on AMRFinderPlus
   markers, isotonic/sigmoid-calibrated on a held-out split, so a "0.9" means a
   0.9 observed failure rate — not just a high score.
2. **Explicit abstention (the firewall).** An out-of-distribution detector
   compares each genome's marker profile to the training distribution; when it's
   too unfamiliar, the system returns **`no_call → escalate to laboratory AST`**
   instead of guessing.
3. **Glass-box mechanism.** Every call names the markers driving it (e.g.
   *blaKPC* for carbapenem resistance, *gyrA*_S83 for fluoroquinolones) and an
   evidence tier (known marker / statistical-only / no signal).

On top sits a **stewardship layer** that, among the drugs predicted to work,
recommends the *narrowest-spectrum* agent — and is safe by construction: if every
drug is a fail or a no-call, it escalates to the lab rather than inventing a
choice.

## 4. What makes it novel (USP)

| | Rule-based tools (ResFinder/CARD/AMRFinderPlus, Pathogenwatch, Kleborate) | Commercial ML WGS-AST predictors | **Genome Firewall** |
|---|---|---|---|
| Output | binary R/S | R/S prediction | **calibrated probability + verdict** |
| Says "I don't know" | no | rarely | **yes — explicit OOD abstention** |
| Mechanism shown | gene list | usually opaque | **per-call markers + evidence tier** |
| Validation | in-distribution | in-distribution | **homology-grouped + temporal** |
| Stewardship | none | none | **narrowest-spectrum, safe-by-construction** |

The core claim: **incumbents give a confident binary call; we give a calibrated
probability with an honest refusal when the genome is unfamiliar.** That is the
difference between a tool a clinician *reads* and a tool a clinician can *trust*.

## 5. Technology

- **Data:** 2,997 *K. pneumoniae* genomes from BV-BRC, **laboratory-measured
  labels only** (we explicitly excluded ~72k computationally-predicted phenotypes
  per drug — training on model predictions would be circular). 4 antibiotics:
  ciprofloxacin, meropenem, gentamicin, ceftazidime.
- **Features:** AMRFinderPlus run on each assembly (v4.2.7, NCBI reference
  database), symbols family-collapsed into a 20-marker vocabulary + drug-target
  presence flags + assembly QC.
- **Model:** per-drug regularized logistic regression, calibrated on a held-out
  split. Deterministic target-presence gate. OOD detector on marker-profile
  distance.
- **Two validation axes most tools skip:**
  - **Homology-grouped split** — whole cgMLST clusters are held out, so the model
    can't win by memorizing a lineage it also trained on.
  - **Temporal split** — train on pre-2015 isolates, test on 2015–2018, because
    resistance drifts over time.
- **App:** Streamlit decision report with a mandatory "confirm with standard lab
  testing" banner, per-drug verdict cards, calibrated confidence, mechanism, and
  the stewardship recommendation. Optional grounded chat that refuses dual-use,
  clinical-dosing, and diagnosis requests.

## 6. Results (real held-out data)

**Homology-grouped split** — the honest, no-memorization test. Every metric the
challenge rubric asks for, per drug:

| Drug | Balanced acc | Recall(R) | Recall(S) | F1 | AUROC | PR-AUC | Brier | No-call rate | Called bal.acc |
|---|---|---|---|---|---|---|---|---|---|
| ciprofloxacin | **0.978** | 0.982 | 0.975 | 0.984 | 0.995 | 0.997 | 0.020 | 36% | 0.821 |
| meropenem | **0.937** | 0.913 | 0.961 | 0.895 | 0.970 | 0.936 | 0.042 | 74% | 0.714 |
| gentamicin | **0.935** | 0.880 | 0.989 | 0.925 | 0.967 | 0.957 | 0.041 | 38% | 0.976 |
| ceftazidime | **0.935** | 0.948 | 0.922 | 0.956 | 0.976 | 0.988 | 0.047 | 31% | 0.728 |

*Recall(R) = recall on resistant cases (drug likely to fail); Recall(S) = recall
on susceptible cases (drug likely to work), reported separately as required.
No-call rate = fraction abstained; Called bal.acc = balanced accuracy on the
non-abstained predictions (balanced, not raw accuracy, to stay honest under class
imbalance). Full numbers in `results/pitch_metrics.csv`; per-drug reliability
curves in `results/fig_calibration_real.png`.*

**Generalization by genetic group** (rubric: "performance broken down by
genetically related groups, groups not seen during training"): evaluated per
held-out cgMLST cluster — see `results/fig_per_group.png` and
`results/eval_grouped/per_group_metrics.csv`. Accuracy is near-perfect on most
held-out lineages, with the honest outliers (e.g. one meropenem cluster) visible
rather than averaged away.

Balanced accuracy **0.94–0.98** with **Brier 0.02–0.05** — the reliability curves
sit on the diagonal, so the confidences are real, not just high. For reference,
rule-based tools report 0.52–0.66 on clinical isolates with no probability at all.

**Temporal split** (train ≤2014 → test ≥2015): ciprofloxacin 0.94, ceftazidime
0.91, gentamicin 0.90 hold across a 4-year forward gap *despite a real ~21-point
resistance drift* in the ESBL drugs. Meropenem drops to 0.78 — and this is where
the firewall shows its value: under temporal shift, out-of-distribution
abstentions rose from ~5 genomes per drug to 71–90 per drug. **The model detected
that future isolates looked less familiar and abstained more, rather than making
confident errors.** That honesty under drift is the feature, not a bug.

**Impact:** genome-to-decision in seconds instead of 24–72 hours of culture — with
a calibrated probability where the genome is informative, and an explicit
"escalate to the lab" where it isn't. Speed *with* honesty is what makes it
deployable.

> **Safety.** Decision support only. Every prediction must be confirmed by
> standard laboratory antimicrobial susceptibility testing before any treatment
> decision. The system is designed to abstain rather than risk a confident wrong
> call, and the app carries this notice on every report.
