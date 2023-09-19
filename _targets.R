################################################################################
#
# Project build script
#
################################################################################

# Load packages (in packages.R) and load project-specific functions in R folder
suppressPackageStartupMessages(source("packages.R"))
for (f in list.files(here::here("R"), full.names = TRUE)) source (f)

# Setup  ------------------------------------------------------------

## Setup for future package for parallelization. This is the max allowable number of workers
nworker <- 6
nthread <- 6
future::plan(list(tweak(future.callr::callr, workers = nworker), tweak(future.callr::callr, workers = nthread)))
future::plan(future.callr::callr, workers = nworker)

# Groups of targets ------------------------------------------------------------

## Group of targets associated with parameter choices that the user must define
setup_targets <- tar_plan(
  
    ## Establish what stan models to fit
    tar_target(models_to_fit,
      c(
   #    "cluster_regression_base.stan"
        "cluster_regression_with_beta.stan"
      , "cluster_regression_with_beta_theta.stan"
      )
    )
    
    ## Into list for future and purrr::pmap
  , tar_target(stan_models.l, {
      data.frame(model = models_to_fit) %>% mutate(index = seq(n())) %>% split_tibble(., "index")
    })
    
    ## Establish parameters for the simulations
  , tar_target(sim.params,
      establish_parameters(
        n_param_sets   = 3
      , n_sims_per_set = 2
      , n_samps        = 1000
      , mu_neg         = -2.75
      , sd_neg         = 1
      , mu_pos_delta   = c(1, 5)
      , sd_pos_delta   = 0.5
      , beta_base      = 0.2 ## 20% of adults are seropositive
      , beta_age_delta = 0.2 ## 40% of juveniles are seropositive
      , mu_theta_age   = 1
      , age_prop       = 0.5
     )
   )
  
)

## Simulate data from established parameters
simulation_targets <- tar_plan(
  
    ## Simulate data
    tar_target(sim.data,
      simulate_data(
        param_sets = sim.params
      )
    )
    
    ## Into list for future and purrr::pmap
  , tar_target(simulated_data.l, {
      sim.data %>% split_tibble(., c("param_set", "sim_num"))
    })
  
    ## Plot raw data
  , tar_target(data_plot,
      examine_data(
        simulated_data = sim.data
      )
    )
  
)

## fitting in a few ways
fitting_targets <- tar_plan(
  
    ## -- Rough 3SD cutoff above the mean of the negative control
     ## we dont have a negative control so use 3sd above the mean of the left group
    
    ## get group assignment based on 3Sd cutoff
    tar_target(three_sd.groups,
      group_via_3sd(
        simulated_data = sim.data
      )
    )
    
    ## then fit the regressions on these group assignments
  , tar_target(three_sd.groups.regression,
      fit_regression(
        groups         = three_sd.groups
      )
    )
  
    ## -- Two stage via mclust and then regression
  
    ## stage one: run mclust
  , tar_target(mculst.groups,
      group_via_mculst(
        simulated_data = sim.data
      )
    )
  
    ## stage two: run regression
  , tar_target(mculst.groups.regression, 
      fit_regression(
        groups         = mculst.groups
      )
    )
  
    ## -- One stage stan model
  
    ## Fir the stan models and return all the raw models in a list
  , tar_target(stan_fits.l, 
      fit_stan_models(
        simulated_data = simulated_data.l
      , param_sets     = sim.params
      , model_names    = stan_models.l
      )
    , pattern   = cross(simulated_data.l, stan_models.l)
    , iteration = "list"
   )
  
    ## Sort these output stan models into a sensible list object to combine with parameters for cleaning
  , tar_target(stan.fits,
      sort_stan_fits(
        stan_fits.l   = stan_fits.l
      , models_to_fit = models_to_fit
      )
    )
  
)

## Tidy up code returned by fitting targets
cleanup_targets <- tar_plan(
  
   ## Summarize regression fits for 3sd
   tar_target(three_sd.groups.regression.summary,
     sort_regression(
       fitted_regressions = three_sd.groups.regression
     , param_sets         = sim.params
     ) %>% mutate(
       model = paste("3sd -- ", model, sep = "")
     )
   )
  
   ## Summarize regression fits for mclust
 , tar_target(mculst.groups.regression.summary,
     sort_regression(
       fitted_regressions = mculst.groups.regression
     , param_sets         = sim.params
     ) %>% mutate(
       model = paste("mclust -- ", model, sep = "")
     )
   )
  
   ## Summarize stan fits and combine with parameter values
 , tar_target(stan.summary,
      summarize_stan_fits(
        model_fits     = stan.fits
      , param_sets     = sim.params
      , simulated_data = sim.data
      )
    )
 
   ## Clean up individual level group ID assignment predictions
 , tar_target(group_assignment,
     calculate_group_assignments(
        three_sd.g     = three_sd.groups
      , mclust.g       = mculst.groups
      , stan.g         = stan.summary$group_pred
     )
   )
 
   ## Get and group estimates on population level seropositivity
 , tar_target(pop_seropositivity,
     calculate_population_seropositivity(
        three_sd.g     = three_sd.groups
      , mclust.g       = mculst.groups
      , stan.g         = stan.summary$prop_seropos
     )
   )

)

## collate targets
collate_targets <- tar_plan(
  
  ## pull together all output for convenience
  tar_target(all.out,
     collate_outputs(
        pop_seropositivity = pop_seropositivity
      , group_assignment   = group_assignment
      , three_sd.sum       = three_sd.groups.regression.summary
      , mclust.sum         = mculst.groups.regression.summary
      , stan.sum           = stan.summary$coef
     )
   )
  
)


## plot output
plotting_targets <- tar_plan(
  
  ## Explore fits for the regression coefficients 
  tar_target(fit.plot,
    plot_summary(
      coef_ests   = all.out$coefficient_ests
    , param_sets  = sim.params
    , coverage    = all.out$coverage
    )
  )
  
  ## Explore individual-level group assignments
, tar_target(group_id.plot,
    plot_group_assignments(
      group_assignment   = group_assignment
    )
  )

  ## Plot population level seropositivity estimates
, tar_target(sero.plot,
    plot_pop_seropos(
      pop_seropositivity = pop_seropositivity
    )
  )

)


# List targets -----------------------------------------------------------------

list(
  setup_targets
, simulation_targets
, fitting_targets
, cleanup_targets
, collate_targets
, plotting_targets
  )
