
#### Code

- To use this code open \_targets.R and work from there (code is
  extensively commented in-line).
- This project uses the targets package to create its analysis pipeline.
  The steps are defined in the \_targets.R file and the workflow can be
  executed by running targets::tar_make(). For details on the use of
  targets see <https://books.ropensci.org/targets/>.

#### Data

- All parameters to recreate the 500 parameter sets analyzed for this
  manuscript are given in the /data folder (see \_targets.R for using
  these to recreate the complete analyzed dataset).

#### Repository Structure and Reproducibility

- `data/` contains csv files of parameter values
- `R/` contains functions used in primary \_targets R pipeline (main
  analysis pipeline)
- `scr_R/` contains R scripts for data exploration and
  manuscript/supplement figure creation
- `stan_models/` contains stan model definitions
