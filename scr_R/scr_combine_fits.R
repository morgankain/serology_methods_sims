sort_stan_fits <- function(stan_fits.l) {
  ## Need a list of lists
  ## outer list is of length n-models
  ## each of these has length = n_param_sets
  model_names <- lapply(stan_fits.l, FUN = function(x) {
    names(x[[1]])
  }) %>% unlist()
  
  fit_details <- lapply(model_names %>% as.list(), FUN = function(x) {
    strsplit(x, " -- ") %>% unlist() %>% t() %>% as.data.frame()
  }) %>% do.call("rbind", .) %>% 
    rename(model = V1, param_set = V2, sim_num = V3, log_mfi = V4) %>% as_tibble()
  
  resorted_fits <- fit_details %>% 
    mutate(
      fitted_model = lapply(stan_fits.l, FUN = function(x) {
        x[[1]]
      })
      , fit_details = lapply(stan_fits.l, FUN = function(x) {
        x[[2]]
      })
      , fit_success = lapply(stan_fits.l, FUN = function(x) {
        ifelse(nrow(x[[2]][[1]]) == 1, 0, 1)
      }) %>% unlist()
    )
  
  return(resorted_fits)
  
}

fits_to_load <- list.files("_targets/objects")
fits_to_load <- fits_to_load[grepl("stan_fits.l", fits_to_load)]

for (i in seq_along(fits_to_load)) {
  temp_name <- fits_to_load[i]
  tar_load(temp_name)
  sorted_out <- sort_stan_fits(stan_fits.l = get(temp_name) %>% list())
  
  temp_sort <- summarize_stan_fits_for_pub(
      model_fits     = sorted_out
    , param_sets     = sim.params
    , simulated_data = sim.data
    , complexity     = data_complexity
  )
  
  rm(list=grep("stan_fits.l_",ls(),value=TRUE,invert=FALSE))
  
  if (((i / 10) %% 1) == 0) { gc() }
  
  if (i == 1) {
    final_sort <- temp_sort
  } else {
    final_sort$coef         <- rbind(final_sort$coef, temp_sort$coef)
    final_sort$group_pred   <- rbind(final_sort$group_pred, temp_sort$group_pred)
    final_sort$prop_seropos <- rbind(final_sort$prop_seropos, temp_sort$prop_seropos)
    final_sort$fit_details  <- rbind(final_sort$fit_details, temp_sort$fit_details)
  }
  
  print(i)
  
}

stan.summary <- final_sort

saveRDS(list(
  sim.params       = sim.params
, sim.data         = sim.data
, all.out          = all.out
, group_assignment = group_assignment
, three_sd.groups  = three_sd.groups
, mculst.groups    = mculst.groups
, stan.summary     = stan.summary
, mculst.groups.regression.summary   = mculst.groups.regression.summary
, three_sd.groups.regression.summary = three_sd.groups.regression.summary
, pop_seropositivity                 = pop_seropositivity
), "test_fits_v3.Rds")

