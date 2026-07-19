# Genome Firewall — Complete Project Recap

*A single reference for prepping a pitch deck or a separate presentation. Everything
about the project in one place: problem, solution, tech, results, team, and the
story behind the choice.*

---

## 1. The one-liner

Genome Firewall reads a *Klebsiella pneumoniae* patient's genome and predicts, per
antibiotic, whether the drug will work — with **calibrated confidence** and an
honest **"I don't know, send it to the lab"** whenever the genome is unfamiliar.
A firewall, not a classifier.

---

## 2. The problem (and why it's personal)

When a patient is admitted with a *Klebsiella* infection, the clinician must choose
an antibiotic immediately — but culture-based susceptibility testing takes 24–72
hours. In severe infection, that delay is measured in lives.

The numbers (all from peer-reviewed literature, citable):
- Carbapenem-resistant *K. pneumoniae* (CRKP) infections carry **~42% pooled
  mortality vs ~21%** for susceptible infections; **~54% in bloodstream
  infections**, **~49% in ICU-admitted patients**. In ICU cohorts with the most
  resistant strains, mortality has been reported **as high as ~70%**.
- The decisive factor is *speed to the right drug*: **30-day mortality was 77%
  when the initial antibiotic was inappropriate, versus 9% when it was
  appropriate.** Inappropriate initial therapy, septic shock, and ICU admission
  are independent risk factors for death.

That 77%-vs-9% gap is the entire thesis: **getting the first antibiotic right is
the difference between living and dying — and that is exactly the decision Genome
Firewall accelerates**, from days to seconds.

**Why Klebsiella, personally.** This problem is not abstract for our team. Our
domain lead chose *Klebsiella pneumoniae* because she lost a family member to a
multidrug-resistant *Klebsiella* infection acquired in the ICU — one of the ~70%
who did not survive. Building a tool that could have given that clinician a fast,
honest answer is the reason this project exists. *(Include or omit per comfort;
it is a genuine and powerful motivation for the Creativity/Presentation scoring.)*

---

## 3. Why existing tools aren't enough

Genome-based AMR tools today (ResFinder, CARD, AMRFinderPlus, Pathogenwatch,
Kleborate) are **rule-based**: they report presence/absence of known resistance
genes and emit a binary resistant/susceptible call. Two trust problems:
- **No calibrated confidence** — a binary call hides how sure the tool is;
  published balanced accuracy on clinical isolates ranges 0.52–0.66.
- **No abstention** — they answer even when the genome carries a mechanism they've
  never characterized, precisely when a confident wrong answer is most dangerous.

---

## 4. Our solution — the firewall

Keeps the interpretability of the rule-based world, adds the three things that make
a prediction trustworthy:

1. **Calibrated probability** — per-drug logistic regression on AMRFinderPlus
   markers, isotonically calibrated on a held-out split, so "0.9" means a 0.9
   observed failure rate.
2. **Explicit abstention (the firewall)** — an out-of-distribution detector returns
   `no_call → escalate to laboratory AST` when the genome is too unfamiliar.
3. **Glass-box mechanism** — every call names the markers driving it (e.g. blaNDM
   for carbapenems, gyrA for fluoroquinolones) with an evidence tier.

On top: a **stewardship layer** that recommends the narrowest-spectrum effective
agent, and is safe by construction — if every drug is a fail or a no-call, it
escalates rather than inventing a choice.

---

## 5. Technology / stack (name by name)

- **Language:** Python
- **ML:** scikit-learn (per-drug logistic regression + isotonic calibration),
  NumPy, pandas, joblib
- **Genomics:** AMRFinderPlus (NCBI, v4.2.7) for marker extraction; BV-BRC as the
  data source; cgMLST for homology-grouped splitting
- **App:** Streamlit (decision-report UI), Plotly + Matplotlib (figures)
- **Optional chat:** OpenAI API (grounded, guardrail-tested assistant)
- **Reports/testing:** fpdf2 (PDF export), pytest (contract + guardrail tests)
- **No JavaScript/React** — the app is pure Python/Streamlit.

**Tag list (for the submission form):** `Python` · `scikit-learn` · `Streamlit` ·
`AMRFinderPlus` · `BV-BRC` · `NumPy` · `pandas` · `Plotly` · `Matplotlib` ·
`Machine Learning` · `Bioinformatics` · `Antimicrobial Resistance`

---

## 6. Results — real held-out data (2,997 K. pneumoniae genomes)

BV-BRC, **laboratory-measured labels only** (~72k computational predictions per
drug excluded to avoid circular training). 4 drugs: ciprofloxacin, meropenem,
gentamicin, ceftazidime.

**Homology-grouped split (whole cgMLST clusters held out — no memorization):**

| Drug | Bal.Acc | Recall(R) | Recall(S) | F1 | AUROC | PR-AUC | Brier | No-call | Called bal.acc |
|---|---|---|---|---|---|---|---|---|---|
| ciprofloxacin | 0.978 | 0.982 | 0.975 | 0.984 | 0.995 | 0.997 | 0.020 | 36% | 0.821 |
| meropenem | 0.937 | 0.913 | 0.961 | 0.895 | 0.970 | 0.936 | 0.042 | 74% | 0.714 |
| gentamicin | 0.935 | 0.880 | 0.989 | 0.925 | 0.967 | 0.957 | 0.041 | 38% | 0.976 |
| ceftazidime | 0.935 | 0.948 | 0.922 | 0.956 | 0.976 | 0.988 | 0.047 | 31% | 0.728 |

- Balanced accuracy **0.94–0.98**, Brier **0.02–0.05** (calibration curves on the
  diagonal). Rule-based tools: 0.52–0.66 with no probability.

**Temporal split (train ≤2014 → test ≥2015, real ~21-pt resistance drift):**
ciprofloxacin 0.935, ceftazidime 0.914, gentamicin 0.903 hold; **meropenem drops
to 0.779** — and under drift the firewall's OOD abstentions rose from ~5/drug to
71–90/drug. The model detected unfamiliar future isolates and abstained more
rather than making confident errors. *Honest degradation is the feature.*

**Figures:** `results/fig_architecture.png`, `fig_calibration_real.png`,
`fig_grouped_vs_temporal.png`, `fig_per_group.png`. Data: `results/pitch_metrics.csv`,
`results/eval_grouped/`.

---

## 7. What makes it novel (USP)

- Incumbents give a **confident binary call**; we give a **calibrated probability
  with an honest refusal** when the genome is unfamiliar.
- We validate on **two axes most tools skip** — unseen lineages (grouped) *and*
  future isolates (temporal).
- We add **stewardship** (narrowest-spectrum, safe-by-construction) that no
  rule-based tool provides.

The difference between a tool a clinician *reads* and a tool a clinician can
*trust*.

---

## 8. The team — who did what

Four people: one domain/clinical lead + three ML engineers. Track structure:

| Track | Scope | Owner |
|---|---|---|
| **A — Genome reader, data & splits** | BV-BRC cohort curation, AMRFinderPlus feature assembly, cgMLST homology-grouped split, temporal split, contracts | **Nicole (domain lead)** — ⟦confirm/fill name⟧ |
| **B — Calibrated models & abstention** | Per-drug logistic regression, isotonic calibration, OOD detector, evaluation metrics | ⟦teammate 2 — fill name⟧ |
| **C — App & grounded chat** | Streamlit decision-report UI, FASTA upload wiring, per-drug verdict cards, grounded assistant | ⟦teammate 3 — fill name⟧ |
| **D — Clinical framing & safety (shared)** | Problem framing, stewardship logic, safety banners, documentation, submission | Shared across all four |

*(Fill in your three teammates' real names before the pitch. The domain lead —
you — owned Track A and the clinical framing; the personal Klebsiella motivation
is yours.)*

**Suggested one-line pitch-deck credits:**
> Four builders — one with a clinical/domain background who lived the problem, and
> three ML engineers — turning a genome into an antibiotic decision a hospital
> could trust tomorrow.

---

## 9. Live app

- **Deploy:** Streamlit Community Cloud (share.streamlit.io) → entry point
  `streamlit_app.py`, branch `main`. See `docs/DEPLOY.md`.
- **Live URL:** ⟦paste after deploy⟧
- **Demo genomes (verified):** `573.13252` (MDR — all four drugs likely_to_fail,
  blaNDM/blaCTX-M markers) and `573.46111` (clean — all four likely_to_work).

---

## 10. Safety statement

Decision support only. Every prediction must be confirmed by standard laboratory
antimicrobial susceptibility testing before any treatment decision. The system is
designed to abstain rather than risk a confident wrong call, and the app carries
this notice on every report.
