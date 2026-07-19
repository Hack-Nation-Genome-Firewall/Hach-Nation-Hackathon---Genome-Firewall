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

### Video A — DEMO (~55s, UI/UX & product flow) — FINAL READ-ALOUD SCRIPT
*Screen: the live app throughout. Have genome 573.13252 (resistant) and 573.46111
(clean) ready in the Demo tab. Read the quoted text; italics are screen cues.*

> "When a patient is admitted with a Klebsiella infection, the doctor has to
> choose an antibiotic immediately — but the lab culture that tells you which drug
> will work takes one to three days. Genome Firewall answers in seconds, straight
> from the bacterium's genome.
>
> *(select genome 573.13252)* This is a real resistant isolate. Meropenem —
> likely to fail, and it shows you why: the blaNDM carbapenemase gene.
> Ciprofloxacin — likely to fail, driven by a gyrA mutation. Every call comes with
> a calibrated confidence, not just a yes or no.
>
> *(select genome 573.46111)* Now a clean isolate. All four antibiotics — likely
> to work, and it recommends the narrowest-spectrum option so we don't waste the
> last-line drugs.
>
> *(show a no_call verdict)* And when a genome looks unlike anything it was trained
> on, it refuses to guess — it returns no-call, escalate to the lab. Every report
> carries the same reminder: confirm with standard lab testing.
>
> Genome Firewall — antibiotic decisions in seconds, honest about what it doesn't
> know. Usable in a hospital tomorrow."

*(All verdicts verified against the real bundle: 573.13252 → all four
likely_to_fail, blaNDM + gyrA markers present; 573.46111 → all four
likely_to_work, no markers.)*

### Video B — TECH (~55s, stack / architecture / implementation)
Maps to the judges' 5 asks: (1) what it does, (2) architecture/APIs, (3) live
proof, (4) challenges+metrics, (5) clear visuals.
- **(0–8s) what it does — show: architecture diagram.** "Genome Firewall reads a
  *Klebsiella* patient's genome and predicts which antibiotics will work — with
  calibrated confidence, and an honest refusal when it isn't sure."
- **(8–22s) architecture & stack — show: architecture diagram, trace the flow.**
  "The pipeline: AMRFinderPlus extracts resistance markers from the genome; a
  per-drug logistic-regression model, isotonic-calibrated, produces a
  probability; a decision layer adds an out-of-distribution detector and a
  target-presence gate. Data is BV-BRC — laboratory-measured labels only."
- **(22–32s) live proof — show: ~5s screencast, Demo tab, genome 573.13252.**
  "Here it flags meropenem resistance on a real isolate and names the gene
  driving it — blaNDM — with its confidence." *(Verified: genome 573.13252 →
  meropenem likely_to_fail, p_fail 0.99, blaNDM among supporting markers. Use
  this exact genome so what's on screen matches the words.)*
- **(32–48s) highlights + metrics — show: calibration figure, then
  grouped-vs-temporal.** "Two hard problems we solved: calibration and honest
  generalization. On **unseen lineages**, balanced accuracy is **0.94 to 0.98**
  with Brier **0.02 to 0.05**, versus 0.52 to 0.66 for rule-based tools. On
  **future isolates** it's lower — balanced accuracy **0.78 to 0.94**, Brier up
  to **0.07** — and that's the point: under time drift the model abstains more
  instead of guessing wrong."
- **(48–55s) close — show: a no_call verdict in the app.** "A firewall that knows
  when it doesn't know."

> NOTE ON NUMBERS (do not misstate): every metric range is **split-specific**.
> Homology-grouped (unseen-lineage): balanced accuracy 0.94–0.98, Brier
> 0.02–0.05. Temporal (future-isolate): balanced accuracy 0.78 (meropenem) to
> 0.94 (ciprofloxacin), Brier 0.05–0.068. Never quote a grouped number as if it
> covers both splits — say which split each figure is from.

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
