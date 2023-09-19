
#### Most recent note/work update

- September 19: Full pipeline for multiple (simple) models with
  summaries working. Next step to explore some coverage and then start
  expanding the complexity of the stan models

### Simulations to compare alternative methods of analyzing serology data

#### Project Vision/Objectives

- Show under coverage of two-step cluster + regression
- Illustrate potential of Bayesian clustering models that simultaneously
  consider covariates that affect grouping and values within those
  groups

#### Repository Structure and Reproducibility

- `data/` contains data (empirical and/or saved simulated)
- `R/` contains functions used in this analysis.
- `stan_models/` contains stan model definitions
- `reports/` contains literate code for R Markdown reports generated in
  the analysis
- `outputs/` contains compiled reports and figures.
- This project uses the `targets` package to create its analysis
  pipeline. The steps are defined in the `_targets.R` file and the
  workflow can be executed by running `targets::tar_make()`.

#### A coverage plot from some initial simulations

<img src="outputs/examp_sims.png" width="75%" />
