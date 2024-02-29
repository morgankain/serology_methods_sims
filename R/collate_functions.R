## collate_outputs
collate_outputs         <- function(
    pop_seropositivity, group_assignment
  , three_sd.sum, mclust.sum, stan.sum
  , coef_name_vec) {
  
  all.out <- rbind(
    stan.sum %>% dplyr::select(-n_samps) %>% relocate(log_mfi, .after = model) %>%
      relocate(sim_num, .before = param_set) %>%
      mutate(method = "Bayesian LCR", .after = param_set)
  , mclust.sum %>% relocate(model, .before = 1) %>%
      mutate(lwr_n = NA, .after = lwr) %>% mutate(upr_n = NA, .after = mid)
  , three_sd.sum %>% relocate(model, .before = 1) %>%
      mutate(lwr_n = NA, .after = lwr) %>% mutate(upr_n = NA, .after = mid)
  )
  
  coverage <- all.out %>% 
    filter(name %in% coef_name_vec) %>%
    group_by(param_set, model, method, log_mfi, name) %>%
    summarize(
      coverage = mean(cover, na.rm = T)
    )
  
  return(
    list(
      pop_seropositivity = pop_seropositivity
    , group_assignment   = group_assignment
    , coefficient_ests   = all.out
    , coverage           = coverage
    )
  )
  
}

