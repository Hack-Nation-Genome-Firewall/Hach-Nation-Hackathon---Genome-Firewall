# Genome Firewall — Project Description

Structured submission answers. Numbers in `⟦…⟧` are filled from
`results/pitch_metrics.csv` once the real 3,000-genome run completes.

---

## 1. Problem & Challenge
A patient is admitted with a *Klebsiella pneumoniae* infection. The clinician
must choose an effective antibiotic **now** — but standard antimicrobial
susceptibility testing (AST) requires culturing the isolate, which takes
**24–72 hours**. In sepsis, every hour to effective therapy raises mortality
(~7.6%/hour in septic shock; Kumar 2006). The clinician is forced to either wait
(dangerous) or prescribe broad-spectrum empirically (drives resistance, and may
still miss).

Whole-genome sequencing of the isolate is fast and increasingly routine, and the
genome already encodes most of the resistance information. The challenge: turn a
genome into a **trustworthy** per-drug susceptibility prediction — one a
clinician can actually rely on, which means it must *know when it doesn't know*.

## 2. Target Audience
- **Primary:** clinical microbiologists and infectious-disease clinicians in
  hospitals, at the point of a positive culture, needing faster susceptibility
  guidance than culture-based AST provides.
- **Secondary:** antimicrobial stewardship teams (narrowest-effective-drug
  selection) and public-health / surveillance labs tracking resistance spread.
- **Explicitly not:** industrial QC or contamination screening (a different
  problem — see USP vs Spore.bio below).

## 3. Solution & Core Features
**Input:** one file — the assembled genome (FASTA) of the patient's isolate.
Nothing else is required from the user. **Output:** a per-drug decision report.

The system:
1. **Reads the genome** (Module 1) — runs AMRFinderPlus, extracts a curated,
   empirically-validated set of resistance markers (acquired genes + resistance
   point mutations) and a genome-quality profile.
2. **Predicts per drug** (Module 2) — a calibrated model per antibiotic
   (meropenem, ciprofloxacin, gentamicin, ceftazidime) returns
   `likely_to_work / likely_to_fail / no_call`, each with a **calibrated
   probability** and the **mechanism** (which marker drove the call).
3. **Decides safely** (Module 3) — a Streamlit report + grounded chat assistant.
   Core safety features, all *by construction*:
   - **Calibrated confidence** — a 0.9 means ~90% in held-out reality (reliability
     curve, not a raw softmax).
   - **Honest abstention (`no_call`)** — on out-of-distribution genomes or
     thin evidence, it refuses to guess and routes to the lab.
   - **Deterministic target-presence gate** — if the drug's target gene is
     absent/disrupted, that is stated as fact, independent of the model.
   - **Stewardship layer** — recommends the *narrowest-spectrum* effective agent.
   - **Mandatory safety banner + grounded chat** — every screen says "confirm
     with laboratory AST"; the assistant refuses clinical/dosing advice and
     dual-use requests, and never invents numbers.

## 4. Unique Selling Proposition (USP) — landscape analysis
**The one-line USP: existing genome-AMR tools output a *call*; we output a
*calibrated call plus an honest "I don't know."* We are a firewall, not a
classifier.**

| Tool / company | Approach | What it outputs | The gap we fill |
|---|---|---|---|
| **ResFinder, CARD, AMRFinderPlus** | rule-based presence/absence of known genes | binary R/S (often by drug *class*) | no calibrated probability, no abstention; balanced accuracy 0.52–0.66 on clinical isolates, varies wildly by drug |
| **Pathogenwatch, Kleborate** | curated species-specific rules + QC report | R/S profile, some MICs | QC-gated but still a hard call; no calibrated uncertainty or OOD abstention |
| **Commercial ML WGS-AST predictors** (e.g. diagnostics startups training ML on whole genomes) | ML on whole genome | R/S prediction | closest peer class, typically proprietary; our differentiator is calibration + explicit no-call + mechanism transparency. *(Verify specific vendors before naming any in the video.)* |
| **Spore.bio / Spore.Labs** | light/spectroscopy, ML on optical signal | contamination presence/quantity (industrial QC); Spore.Labs now open-source AMR | different modality & setting (factory QC, not clinical genome→drug decision) |

Concretely, our four differentiators — none of which the rule-based incumbents
provide together:
1. **Calibrated uncertainty** (reliability curves per drug), so a probability is
   clinically meaningful.
2. **Honest abstention** — an out-of-distribution detector + evidence-tier logic
   that returns `no_call` instead of a confident wrong answer. The *safety* case,
   not the accuracy case, is the headline.
3. **Mechanism transparency** — every call names the marker; glass-box, not
   black-box.
4. **Two-axis validation rigor** — homology-grouped split (can't memorize
   lineages) *and* a temporal split (train on 2001–2014, test on 2015–2018),
   because resistance drifts (ceftazidime/ciprofloxacin R-rate fell ~21 pts
   across those eras in our own cohort). Most tools report neither.

## 5. Implementation & Technology
- **Data:** BV-BRC (PATRIC), **laboratory-measured phenotypes only** (~70k
  computational predictions per drug explicitly excluded). 3,000 QC-selected
  *K. pneumoniae* genomes, conflict-excluded.
- **Genome reader:** AMRFinderPlus 4.2.7 (DB 2026-05-15.1), run on the assembled
  FASTA; a custom symbol-mapping layer collapses allele families into a stable,
  empirically-derived marker vocabulary (20 markers + target-gene flags + QC).
- **Model:** one L2-regularized logistic regression per antibiotic (the brief's
  baseline, done well) with **isotonic calibration** on a held-out calibration
  split; a deterministic target-presence gate and an OOD (Hamming-distance)
  abstention layer wrap the model.
- **Splits:** cgMLST HC10 homology-grouped split (no cluster spans train/test) +
  a temporal split. Zero cluster leakage enforced in code.
- **Evaluation:** balanced accuracy, recall(R) and recall(S) separately, F1,
  AUROC, PR-AUC, Brier, reliability curve, and no-call rate — reported on both
  splits.
- **App:** Streamlit; grounded LLM report-assistant (gpt-4o-mini, optional,
  behind one swappable function; key gitignored). Runs fully without the key.
- **Compute:** the 3,000-genome annotation runs on an HPC cluster (SLURM), ~1s
  per genome parallelized; the model and app run on a laptop.
- **Reproducibility:** frozen feature-spec contract; every module validates
  against it; scripted end-to-end (`scripts/annotate_cluster.sh`,
  `scripts/make_pitch_results.py`).

## 6. Results & Impact
*Real held-out numbers from `results/pitch_metrics.csv` — 2,997 K. pneumoniae
genomes (BV-BRC, laboratory-measured labels only), AMRFinderPlus features.*

**Held-out performance — homology-grouped split** (whole cgMLST clusters held
out; the model cannot memorize a lineage):

| Drug | Balanced acc | Recall(R) | Recall(S) | AUROC | Brier | No-call rate |
|---|---|---|---|---|---|---|
| ciprofloxacin | **0.978** | 0.98 | 0.97 | 0.995 | 0.020 | 36% |
| meropenem | **0.937** | 0.91 | 0.96 | 0.970 | 0.042 | 74% |
| gentamicin | **0.935** | 0.88 | 0.99 | 0.967 | 0.041 | 38% |
| ceftazidime | **0.935** | 0.95 | 0.92 | 0.976 | 0.047 | 31% |

- **Calibration:** reliability curves hug the diagonal (Brier 0.02–0.05) — the
  confidences are trustworthy, not just high. Rule-based incumbents report
  balanced accuracy 0.52–0.66 on clinical isolates with no probability at all.
- **Temporal generalization:** trained on pre-2015 isolates, tested on
  2015–2018 despite a real ~21-point resistance drift in the ESBL drugs —
  balanced accuracy holds at ciprofloxacin 0.94, ceftazidime 0.91, gentamicin
  0.90; meropenem drops to 0.78 (honest degradation, see below). A test most
  tools never run.
- **Safety behavior under shift:** out-of-distribution abstentions rise from
  ~5 genomes per drug on the grouped split to 71–90 per drug on the temporal
  split — the firewall *detects* that future isolates look less familiar and
  abstains more, rather than guessing. Meropenem's temporal drop is absorbed by
  a rising no-call rate (recall(S)=1.00 on the calls it does make), not by
  confident errors. Every susceptible case yields a narrowest-spectrum
  stewardship recommendation.

**Impact:** genome-to-decision in seconds instead of 24–72 h of culture, with
calibrated confidence and an explicit "escalate to the lab" whenever the model
is uncertain — faster effective therapy where the genome is informative, and no
false reassurance where it isn't. That combination — speed *with* honesty — is
what makes it deployable, and what distinguishes it from every rule-based
incumbent.

> Decision support only. Every prediction must be confirmed by standard
> laboratory antimicrobial susceptibility testing before any treatment decision.
