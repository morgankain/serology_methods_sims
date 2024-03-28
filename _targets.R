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
# future::plan(list(tweak(future.callr::callr, workers = nworker), tweak(future.callr::callr, workers = nthread)))
future::plan(future.callr::callr, workers = nworker)

# Groups of targets ------------------------------------------------------------

## Group of targets associated with parameter choices that the user must define
setup_targets <- tar_plan(
  
    ## Control complexity level of data simulation.
    ## NOTE: setting data_complexity to 1 will ignore some parameters in "sim.params" and
     ## make it impossible to fit certain models (those with an _2 at the end of their names)
      ## 1: One categorical predictor affecting group identity
       ## : One categorical predictor affecting within-group values
      ## 2: One categorical and one continuous predictor affecting group identity
       ## : One categorical fixed and one categorical random effect affecting group identity
      tar_target(data_complexity, 2)
  
  , tar_target(models_to_compile,
      establish_models(       
        model_set = c(
         "publication_model_normal_2.stan"
       , "publication_model_lognormal_2.stan"
       , "publication_model_skew_normal_2.stan"
       , "publication_model_skew_normal_wf_2.stan"
       )
     , complexity = data_complexity
     )
   )
      
    ## Establish what stan models to fit
    ## NOTE: _X controls -the minimal- data complexity to fit this model. Make sure if a model is listed
     ## here with _X, data_complexity >= X
  , tar_target(models_to_fit.l,
      establish_models(       
        model_set = c(
         "publication_model_2.stan"
       , "publication_model_skew_2.stan"
       )
     , complexity = data_complexity
     ) %>% as.data.frame() %>% rename(model_base_names = 1) %>% split_tibble(., "model_base_names")
   )
  
    ## Compile stan models and stick them into a tibble that becomes a list for future and purrr::pmap
  , tar_target(stan_models,
      compile_stan_models(
        model_set = models_to_compile
      ) 
   )
   
    ## Establish parameters for the simulations
  , tar_target(sim.params,
      establish_parameters_for_pub(
        ## Complexity, which is used to make sure the correct parameters are listed here
        complexity       = data_complexity
        
        ## Simulation and sample size
      , n_param_sets     = 10
      , n_sims_per_set   = 1
      , n_samps          = c(100, 2000)
      
        ## Sample composition
         ## catx as a generic stand-in for some categorical difference in the sample
      , cat1f_prop       = 0.5  
      , cat2f_prop       = 0.5
      , cat1r_count      = 10
      , con1f_sd         = 2
    # , con2f_sd         = 2
    
        ## Group identity covariates (all on logit scale)
      , beta_base        = c(-5, 0)
      , beta_cat1f_delta = c(0.2, 2)
      , beta_cat2f_delta = c(0.2, 2)
      , beta_con1f_delta = c(0.1, 1)
      
      , mu_neg            = c(-7, -3)
      , sd_neg            = c(0.2, 1.5) 
      , mu_pos_delta      = c(0.1, 5)
      , sd_pos_delta      = c(1, 2)
    #  , theta_con2f_delta = c(0, 0.5)
    #  , theta_con1f_delta = c(0.1, 1)
    #  , theta_cat2f_mu   = 0 
    #  , theta_cat1r_sd   = 0.5 
  
      , logit_1          = 30000
      , logit_2          = -1
      , logit_3          = 1
     )
   )
  
)

## Simulate data from established parameters
simulation_targets <- tar_plan(
  
    ## Simulate data
    tar_target(sim.data,
     simulate_data_for_pub(
        param_sets = sim.params
      , complexity = data_complexity
      ## take the log of the data in every simulation -- will fit to both
      ) %>% mutate(log_mfi = log(mfi)) %>% 
        pivot_longer(., c(mfi, log_mfi), names_to = "log_mfi", values_to = "mfi")
    )
    
    ## Into list for future and purrr::pmap
  , tar_target(simulated_data.l, {
      sim.data %>% split_tibble(., c("param_set", "sim_num", "log_mfi"))
    })
  
    ## Skew of raw data
  , tar_target(sim.data.summaries,
      calc_sim_summaries(
        simulated_data = sim.data
      , param_sets     = sim.params
      )
    )
  
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
      , param_sets     = sim.params
      , groupings      = c("param_set", "sim_num", "log_mfi")
      )
    )
    
    ## then fit the regressions on these group assignments
  , tar_target(three_sd.groups.regression,
      fit_regression(
        groupings    = three_sd.groups
      , gam_formula  = {
        if (data_complexity == 1) {
          "assigned_group ~ cat1f"
        } else if (data_complexity == 2) {
          "assigned_group ~ cat1f + cat2f + con1f"
        } else {
          NULL
        }
      }
      , complexity   = data_complexity
      , groupings1   = c("param_set")
      , groupings2   = c("log_mfi", "sim_num", "sd_method")
      , method       = "sd"
      )
    )
  
    ## -- Two stage via mclust and then regression
  
    ## stage one: run mclust
  , tar_target(mclust.groups,
      group_via_mclust2(
        simulated_data = sim.data
      , groupings      = c("param_set", "sim_num", "log_mfi")
      )
    )

    ## stage two: run regression
  , tar_target(mclust.groups.regression, 
      fit_regression(
        groupings    = mclust.groups
      , gam_formula  = {
        if (data_complexity == 1) {
          "assigned_group ~ cat1f"
        } else if (data_complexity == 2) {
          "assigned_group ~ cat1f + cat2f + con1f"
        } else {
          NULL
        }
      }
      , complexity   = data_complexity
      , groupings1   = c("param_set")
      , groupings2   = c("log_mfi", "sim_num", "method")
      , method       = "mclust"
      )
    )
  
    ## -- One stage stan model
     ## NOTE: Stan pipeline here does not perfectly parallel other options in terms of fit then summarize because of memory
     ## considerations (cant load all the stan models at once and then summarize), so instead, this function fits the model
     ## and does the summary right away
  , tar_target(needed_stan_fits,
     establish_stan_combos(
       simulated_data  = simulated_data.l
     , models_to_fit   = models_to_fit.l
     )
  )
  
  , tar_target(stan.summary.l, 
     fit_stan_models_for_pub(
        simulated_data  = simulated_data.l
      , param_sets      = sim.params
      , compiled_models = stan_models
      , model_names     = models_to_fit.l
      , data_complexity = data_complexity
      , max_time        = 60*20
      )
    , pattern   = cross(simulated_data.l, models_to_fit.l)
    , iteration = "list"
    , error     = "null" 
   )
  
)

## Tidy up code returned by fitting targets
cleanup_targets <- tar_plan(
  
   ## Summarize regression fits for 3sd
   tar_target(three_sd.groups.regression.summary,
     sort_regression(
       fitted_regressions = three_sd.groups.regression
     , param_sets         = sim.params
     , complexity         = data_complexity
     , groupings1         = c("param_set")
     , groupings2         = c("log_mfi", "sim_num", "method")
     ) %>% mutate(
       model = paste("3sd -- ", model, sep = "")
     )
   )
  
   ## Summarize regression fits for mclust
 , tar_target(mclust.groups.regression.summary,
     sort_regression(
       fitted_regressions = mclust.groups.regression
     , param_sets         = sim.params
     , complexity         = data_complexity
     , groupings1         = c("param_set")
     , groupings2         = c("log_mfi", "sim_num", "method")
     ) %>% mutate(
       model = paste("mclust -- ", model, sep = "")
     )
   )
 
   ## Combine all of the stan summaries
 , tar_target(stan.summary,
      summarize_stan_summary(
        stan_summary   = stan.summary.l
      , param_sets     = sim.params
      )
    )
 
   ## Clean up individual level group ID assignment predictions
 , tar_target(group_assignment,
     calculate_group_assignments(
        three_sd.g     = three_sd.groups
      , mclust.g       = mclust.groups
      , stan.g         = stan.summary$group_pred
      , param_sets     = sim.params
     )
   )
 
   ## Get and group estimates on population level seropositivity
 , tar_target(pop_seropositivity,
     calculate_population_seropositivity(
        three_sd.g     = three_sd.groups
      , mclust.g       = mclust.groups
      , stan.g         = stan.summary$prop_seropos
      , param_sets     = sim.params
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
      , mclust.sum         = mclust.groups.regression.summary
      , stan.sum           = stan.summary$coef
      , coef_name_vec      = c(
          "beta_base"
        , "beta_cat1f_delta"
        , "beta_cat2f_delta"
        , "beta_con1f_delta"
      )
     )
   )
  
)

## plot output
plotting_targets <- tar_plan(
  
  ## Explore fits for the regression coefficients 
  tar_target(fit.plot,
    plot_summary(
      coef_ests     = all.out$coefficient_ests
    , param_sets    = sim.params
    , coverage      = all.out$coverage
    , coef_name_vec = c(
        "beta_base"
      , "beta_cat1f_delta"
      , "beta_cat2f_delta"
      , "beta_con1f_delta"
      , "theta_con2f_delta"
      )  
    )
  )
  
  ## Explore individual-level group assignments
, tar_target(group_id.plot,
    plot_group_assignment_summary(
      group_assignment   = group_assignment
    )
  )

, ## Plot individual-level group assignment probabilities
  tar_target(ind_group_prob.plot,
    plot_individual_group_prob(
      three_sd.g     = three_sd.groups
    , mclust.g       = mclust.groups
    , stan.g         = stan.summary$group_pred
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
