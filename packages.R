######################## LOAD R PACKAGES #######################################

################################################################################
#
#' R packages needed to run any/most {targets} workflows
#
################################################################################

library(targets)
library(tarchetypes)
library(tidyverse)
library(here)
library(knitr)
library(rmarkdown)
library(future)

################################################################################
#
#' Additional R packages needed to run your specific workflow
#' 
#' * Delete or hash out lines of code for R packages not needed in your workflow
#' * Insert code here to load additional R packages that your workflow requires
#
################################################################################

library(dplyr)
library(tidyr)
library(magrittr)
library(MASS)
library(ggplot2)
library(GGally)
library(mclust)
library(factoextra)
library(rstan)
library(cmdstanr)
library(effects)
library(pomp)

## Needed for internal calls in stan | targets
conflicted::conflicts_prefer(rstan::extract)
conflicted::conflicts_prefer(stats::lag)
conflicted::conflicts_prefer(purrr::map)
conflicted::conflicts_prefer(future::run)
