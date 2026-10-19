# plausl

Code for **Predicting Language and Literacy Growth in Children With ASD Using Statistical Learning (PLAUSL)**.
ClinicalTrials.gov [NCT06332144](https://clinicaltrials.gov/study/NCT06332144); funded by NIH R56DC020208.

First-grade children with and without autism completed syllable and letter statistical learning (SL) tasks, each
with a two-alternative forced-choice test (2AFC) and statistically induced chunking recall (SICR), and standardized
language and reading assessments at two visits about three months apart.

## Contents

| Path | What it does |
|---|---|
| `scripts/score_plausl_sl.py` | Scores the SL tasks: 2AFC accuracy and SICR recall for target and foil trials. |
| `scripts/merge_sl_into_clinical.py` | Computes the language and reading composites and the TILLS language/literacy disorder status. |
| `analysis/sap_aim1.R` | Aim 1 analyses: learning checks, reliability, concurrent and growth models, group comparisons, sensitivity analyses. |
| `docs/methods.md` | Scoring rules and composite definitions. |

## Task materials

| Task | Gorilla Open Materials |
|---|---|
| Syllable statistical learning (SSL) | https://app.gorilla.sc/openmaterials/1351289 |
| Phonological serial recall task (PSRT) | https://app.gorilla.sc/openmaterials/1352806 |

## Reproducing the results

Download the de-identified dataset from Zenodo: [doi:10.5281/zenodo.23462288](https://doi.org/10.5281/zenodo.23462288). From the folder holding `plausl_participants.csv`
and `plausl_sl_trials.csv`:

```bash
# 1. SL scores from the trial file
python3 /path/to/plausl/scripts/score_plausl_sl.py --from-trials plausl_sl_trials.csv rescored/
# 2. composites and TILLS status from the standard scores
python3 /path/to/plausl/scripts/merge_sl_into_clinical.py --from-shared plausl_participants.csv derived_scores.csv
# 3. Aim 1 analyses
PLAUSL_SHARED=. Rscript /path/to/plausl/analysis/sap_aim1.R > sap_aim1_output.txt
```

Steps 1 and 2 reproduce the scores in `plausl_participants.csv`. Step 3 reproduces the Aim 1 estimates; permutation
and bootstrap results can differ slightly between runs on differently ordered data.

Requirements: Python 3.8+ (standard library only); R 4.4 with readxl, readr, jsonlite, lme4, lmerTest, sandwich,
dplyr and tidyr.

## License

Code (`scripts/`, `analysis/`): [MIT](LICENSE). Documentation (`docs/`, this README):
[CC BY 4.0](LICENSE-CC-BY-4.0.txt).

Please cite this repository as described in [CITATION.cff](CITATION.cff).
