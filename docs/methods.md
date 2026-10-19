# Scoring rules and composite definitions

## Statistical learning tasks

Two tasks, each given at both visits:

- **Letter SL** (visual): letter triplets GJA, MDH, FKC, LBE embedded in a continuous stream.
- **Syllable SL** (auditory): trisyllabic words da-ro-pi, go-la-tu, pa-bi-ku, ti-bu-do embedded in a continuous stream.

### 2AFC (two-alternative forced choice)

- The correct option on each trial is the one that is a familiarization word; answer keys are derived from the
  stimuli shown, not from hard-coded lists.
- **Letter:** 16 test trials. Trial 1 follows the practice screen and is treated as practice; trial 17 repeats trial 1.
  Only the first click per trial counts, and timeouts are excluded.
- **Syllable:** trial 7 is dropped for everyone (in the older task version neither option is a word), leaving 15 trials.
- Score: proportion correct (chance = .50).

### SICR (statistically induced chunking recall)

- Each trial presents six items; the child recalls them in order. Target trials are made of two familiarization
  words; foil trials use the same items in other orders.
- **Letter keys** are the letter sequences actually displayed. The response on trial *t* answers the display on trial
  *t* − 1, and the first display is a warm-up. One task version showed JGA in recall although the learned word is GJA;
  these runs are scored as displayed and flagged.
- **Syllable** responses come from transcriptions of audio recordings.
- Short responses are padded with "xx" in the arrangement that maximizes items in the correct position.
- Scores per trial: items in correct serial position (out of 6), key bigrams present (positions 1–2, 2–3, 4–5,
  5–6; out of 4), trigrams in position (out of 2), and items present anywhere. Each is averaged separately over
  target and foil trials.
- A child's SICR scores are left blank if they gave 2 or fewer responses of either trial type, or if every
  transcribed response is "xx".


### Linguistic SL composite

For each task at each session:

> task composite = mean of z(2AFC accuracy) and z(SICR bigram recall on target trials)

> linguistic SL composite = mean of the letter and syllable task composites

Each measure is z-scored within its session across the analysis sample. A composite averages the components a child
has. Models with the SL composite as a predictor adjust for SICR bigram recall on foil trials (the mean of the two
tasks' within-session z-scores).

## Assessment composites

Means of available standard scores (mean 100, SD 15); at least two are required. TILLS subtests (mean 10, SD 3) are
converted with 100 + 5 × (score − 10). Scores outside the valid range are excluded.

- **Language, visit 1:** TILLS Vocabulary Awareness, Phonemic Awareness and Nonword Repetition; NIH Toolbox Picture
  Vocabulary; TROG-2.
- **Language, visit 2:** NIH Toolbox Picture Vocabulary and TROG-2 (TILLS was not given at visit 2).
- **Reading:** WRMT-III Letter Identification, Word Identification and Word Attack; NIH Toolbox Oral Reading
  Recognition.

### Language/literacy disorder status

The TILLS Identification Core for ages 6;0–7;11 is the sum of the Vocabulary Awareness, Phonemic Awareness and
Nonword Repetition standard scores. A sum below **24** is classified as a language/literacy disorder (TILLS Technical
Manual, Tables 2.2–2.3).
