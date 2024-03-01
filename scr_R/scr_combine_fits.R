fits_to_load <- list.files("_targets/objects")
fits_to_load <- fits_to_load[-c(1:9)]
fits_to_load <- fits_to_load[-c(506:510)]

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
  }
  
  print(i)
  
}

saveRDS(list(
  sim.params = sim.params
, sim.data   = sim.data
, all.out    = all.out
, group_assignment = group_assignment
, three_sd.groups = three_sd.groups
, mculst.groups = mculst.groups
, stan.summary = stan.summary
, mculst.groups.regression.summary = mculst.groups.regression.summary
, three_sd.groups.regression.summary = three_sd.groups.regression.summary
, pop_seropositivity = pop_seropositivity
), "test_fits.Rds")



