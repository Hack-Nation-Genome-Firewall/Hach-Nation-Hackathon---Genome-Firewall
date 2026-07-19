# Genome Firewall — pitch narrative + video scripts + demo script

The thesis in one line: **this is a firewall, not a classifier.** A classifier
guesses on every input; a firewall is trusted because it *refuses to guess when
it shouldn't*. Everything below sells that distinction.

Fill the `⟦…⟧` slots from `results/pitch_metrics.csv` once the real run finishes.

> **Submission format (confirmed):** NO slide deck. Three ≤60-second videos
> (Demo / Tech / Team) + the structured written description in
> `docs/PROJECT_DESCRIPTION.md`. Judging = Creativity · Technical depth ·
> Presentation. The scripts below are written to those three videos.

---

## 0. The three 60-second video scripts

### Video A — DEMO (≤60s, UI/UX & product flow)
Screen-record the Streamlit app; voiceover:
- (0–10s) "A patient has a *Klebsiella* infection. Which antibiotic works?
  Culture takes 1–3 days. We answer from the genome in seconds."
- (10–35s) Upload/select a genome → the report appears. Point at one
  `likely_to_fail` with its mechanism (blaKPC), one `likely_to_work`, and the
  **calibrated confidence** bar. Say "each call shows *why*, and how sure it is."
- (35–50s) The **money shot**: the out-of-distribution genome → **no_call →
  escalate to lab.** "When it hasn't seen anything like this strain, it refuses
  to guess. That's the firewall."
- (50–60s) Ask the chat "why is meropenem a no_call?" → grounded plain-English
  answer. "And it explains itself — no jargon, no invented numbers."

### Video B — TECH (≤60s, stack / architecture / implementation)
Screen = architecture diagram + a scroll of the repo:
- (0–15s) "Input: one genome FASTA. AMRFinderPlus extracts resistance markers →
  a calibrated logistic-regression per drug → a decision layer."
- (15–35s) "Three things make it trustworthy: **isotonic calibration** so
  probabilities are real; an **out-of-distribution detector** that abstains; a
  **deterministic target-presence gate**. Labels are BV-BRC
  laboratory-measured only — no computational phenotypes."
- (35–55s) "We validate two ways most tools skip: a **homology-grouped** split so
  it can't memorize lineages, and a **temporal** split — train on the past,
  test on the future — because resistance drifts. Real held-out numbers:
  balanced accuracy ⟦⟧, Brier ⟦⟧, no-call rate ⟦⟧."
- (55–60s) "Every module validates against one frozen feature contract.
  Reproducible end to end."

### Video C — TEAM (≤60s, who built it)
- Each member: name, role, one sentence. Suggested roles:
  Track A (genome reader + data + splits) · Track B (calibrated models +
  abstention) · Track C (app + grounded chat) · domain/clinical framing & safety.
- Close on the shared line: "Four of us — [bio/domain] + three ML engineers —
  built a tool a hospital could use tomorrow."

---

## 1. The 90-second story (say this, in order)

1. **The problem, made visceral.** A patient is septic. The clinician must pick
   an antibiotic *now*. Every hour to effective therapy raises mortality
   (~7.6%/hr in septic shock, Kumar 2006). Wait for the lab culture → 24–48 h
   lost. Guess → you might pick a drug the bug already resists.
2. **What we built.** From the bacterial genome alone, in seconds, a per-drug
   verdict — *likely to work / likely to fail / no call* — for four antibiotics,
   each with a **calibrated** confidence and the **mechanism** behind it.
3. **Why it can be trusted — the firewall move.** It abstains. On a genome
   unlike anything it trained on, or where evidence is thin, it returns
   **no_call** and routes to the lab, instead of a confident wrong answer. Our
   headline metric isn't accuracy — it's *accuracy on the calls it chooses to
   make, plus how honestly it abstains on the rest.*
4. **The proof it generalizes.** Two independent stress tests: (a) a
   **homology-grouped** split — no near-identical strain is in both train and
   test, so it can't cheat by memorizing lineages; (b) a **temporal** split —
   train on 2001–2014 isolates, test on 2015–2018, because resistance *drifts*
   (ceftazidime & ciprofloxacin resistance fell ~21 pts across those eras in our
   own cohort). A model that survives both is one you'd actually deploy.
5. **The clinical payoff.** It doesn't stop at prediction — it recommends the
   **narrowest-spectrum** drug predicted to work, sparing last-resort
   carbapenems. That's antibiotic stewardship: the single best lever against the
   very resistance this tool detects.
6. **The honesty close.** Every screen carries "confirm with laboratory AST."
   We are decision *support*, and we designed the system to know its limits. That
   is exactly what makes it safe to put in front of a clinician.

---

## 2. Live demo (3 genomes, ~2 min) — rehearse this exact path

Run the Streamlit app (`TrackC/app.py`); upload/select three genomes chosen to
show the three behaviors:

| genome | what it shows | judge takeaway |
|---|---|---|
| a clearly resistant one (e.g. carbapenemase+) | `likely_to_fail` on meropenem, **mechanism shown** (blaKPC), high confidence | "it explains itself" |
| a susceptible one | `likely_to_work`, stewardship picks the **narrowest** agent | "it's clinically useful, not just a label" |
| an out-of-distribution / thin-evidence one | **no_call → escalate to lab** | "it refuses to guess — that's the firewall" |

Pick the three concrete genome IDs from the held-out test set after the run;
the third is the money shot — **make the app say "I don't know, go to the lab"
on stage.** Judges remember the model that admits uncertainty.

---

## 3. The metrics slide (numbers, no adjectives)

Pull from `results/pitch_metrics.csv`. Show a compact table, both splits:

| drug | split | balanced acc | recall(R) | recall(S) | AUROC | Brier | no-call % |
|---|---|---|---|---|---|---|---|
| ceftazidime | grouped | ⟦⟧ | ⟦⟧ | ⟦⟧ | ⟦⟧ | ⟦⟧ | ⟦⟧ |
| ceftazidime | temporal | ⟦⟧ | ⟦⟧ | … | | | |
| … four drugs × two splits … | | | | | | | |

Plus two figures: the **reliability/calibration curve** (predicted vs observed —
proves "calibrated" isn't a buzzword) and **grouped-vs-temporal balanced
accuracy** (`results/fig_grouped_vs_temporal.png`). If temporal is close to
grouped, say so — "it holds up on future isolates." If it drops, say *that* too
and show the no-call rate rising to absorb the uncertainty — that's the firewall
working, not a weakness.

---

## 4. Anticipated judge questions (have the answer ready)

- *"How do you know it's not memorizing?"* → grouped split, no shared clusters;
  we also report the temporal split, a harder test most teams skip.
- *"What about a resistance mechanism you've never seen?"* → out-of-distribution
  detector abstains (`no_call`, reason `out_of_distribution`) instead of
  guessing. Show the third demo genome.
- *"Why should a clinician trust a probability?"* → calibration curve; a 0.9 from
  us means ~90% in held-out reality, not just a high softmax.
- *"Isn't this just an AMR gene lookup?"* → the target-presence gate is the
  deterministic floor; the calibrated model + abstention + stewardship on top is
  what makes it a decision tool, not a BLAST hit.
- *"Labeled data source?"* → BV-BRC, **laboratory-measured phenotypes only**
  (we explicitly excluded ~70k computational predictions per drug), one species
  done rigorously.

---

## 5. Explicitly scoped OUT (say it before they ask — it reads as discipline)
- No de-novo target/sequence generation — this is a *defensive* tool.
- No longitudinal per-patient evolution tracking — needs serial isolates we
  don't have; the temporal-drift result *is* the evolution evidence, live
  surveillance is roadmap.
- No multi-tool feature merging — a K. pneumoniae benchmark (cited by our team)
  shows it *lowers* ML performance; one well-curated annotator is the right call.
- Single species (K. pneumoniae) by design; the pluggable annotator + frozen
  feature contract are what would let it generalize — shown as architecture, not
  overclaimed as done.

---

## 6. Roadmap slide (the vision, 20 seconds)
Kleborate as the K. pneumoniae-specialist annotator · multi-species via the same
contract · lineage surveillance dashboard for resistance drift · stewardship-aware
combination-therapy suggestions. One slide, spoken fast, then stop.
