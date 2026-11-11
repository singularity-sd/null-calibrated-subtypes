# Are cognitive subtypes real? Null-calibrated consensus clustering for test batteries: a simulation study

Null-calibrated consensus clustering for cognitive test batteries: which procedures invent subtypes, from a Singularity v3 forecast.

**Paper:** [paper.pdf](paper.pdf) · **Page:** [singularitysd.com/research/subtypes.html](https://singularitysd.com/research/subtypes.html) · Preprint, not yet peer reviewed.

## Abstract

Subtyping studies in cognitive and clinical psychology cluster test scores and report profiles. Most cognitive batteries, however, share one general factor, and many scores are skewed or have ceiling effects. We ask how often four procedures report subtypes in simulated batteries (300 participants, 8 tests) that have none, and how often they find subtypes that are there: $k$-means with a silhouette threshold, latent profile analysis (diagonal Gaussian mixtures chosen by BIC), full-covariance Gaussian mixtures chosen by BIC, and consensus clustering calibrated against a correlation-preserving null — an approach developed in cancer genomics and, to our knowledge, not used in psychology. Latent profile analysis reported subtypes in 100% of datasets, with or without real subtypes. Full-covariance mixtures reported subtypes in 100% of skewed or ceiling-limited null datasets and in none of the datasets with real subtypes. $k$-means reported false subtypes in 45–55% of non-normal null datasets and found real ones in at most 20%. Calibrated consensus clustering was the only procedure whose rate differed between the two situations: 15–30% false discoveries (above its nominal 9%) and 30–95% detection, increasing with separation. The common practice of choosing the number of latent profiles by BIC cannot, in this setting, tell subtypes from a single continuous ability.

## How it was made

This study started from the Singularity forecast **T3v3-2026-0038**: a method predicted to enter a new field. It was then tested once, with the plan fixed before looking at the results.

## What is in this repository

| file | what it is |
|---|---|
| `experiment.py` | reproduces the study end to end |
| `fig_rates.png` | figure |
| `main.tex` | LaTeX source of the paper |
| `paper.pdf` | the paper |
| `results.json` | every number reported in the paper |

## Reproduce

```bash
pip install matplotlib numpy scikit-learn
python experiment.py
```

It writes `results.json` and the figures next to the script.

## Cite

```bibtex
@misc{singularity_null_calibrated_subtypes_2026,
  title  = {Are cognitive subtypes real? Null-calibrated consensus clustering for test batteries: a simulation study},
  author = {Singularity},
  year   = {2026},
  url    = {https://github.com/singularity-sd/null-calibrated-subtypes}
}
```

## License

Code: MIT (`LICENSE`). Paper, text and figures: CC BY 4.0 (`LICENSE-paper.md`). Third-party data keep their own licence where stated.

---

Part of [Singularity](https://singularitysd.com): forecasting where science goes next, then testing it.
