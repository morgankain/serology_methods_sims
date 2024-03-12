## Simulate data for all parameter sets
simulate_data           <- function(param_sets, complexity) {
  
param_sets %<>% split_tibble(., c("param_set", "sim_num"))
  
if (complexity == 1) {
  
lapply(param_sets, FUN = function(x) {
  
  mu_vec <- with(x, c(mu_neg, mu_pos))
  sd_vec <- with(x, c(sd_neg, sd_pos))
  
  data.frame(
    cat1f  = with(x, rbinom(n_samps, 1, cat1f_prop))
  , cat2f  = with(x, rbinom(n_samps, 1, cat2f_prop))
  ) %>% mutate(
    group = rbinom(n(), 1
                   , plogis(
                      x$beta_base + 
                      x$beta_cat1f_delta * cat1f
                     )
                   ) + 1
    , titer = rnorm(n()
                    , mu_vec[group] + 
                      (cat2f * x$theta_cat2f_mu * (group - 1))
                    , sd_vec[group]
    )
    , mfi   = logit2(x$logit_1, x$logit_2, x$logit_3, titer)
  ) %>% mutate(
      param_set = x$param_set
    , sim_num   = x$sim_num
    , .before   = 1
    )

}) %>% do.call("rbind", .)
  
} else if (complexity == 2) {
  
lapply(param_sets, FUN = function(x) {
  
  mu_vec   <- with(x, c(mu_neg, mu_pos))
  sd_vec   <- with(x, c(sd_neg, sd_pos))
  rand_dev <- rnorm(x$cat1r_count, 0, x$theta_cat1r_sd)

data.frame(
    cat1f  = with(x, rbinom(n_samps, 1, cat1f_prop))
  , cat2f  = with(x, rbinom(n_samps, 1, cat2f_prop))
  ) %>% mutate(
    cat1r  = with(x, sample(seq(cat1r_count), n(), replace = T))
  , con1f  = with(x, rnorm(n_samps, 0, con1f_sd))
  ) %>% 
  mutate( 
    cat1r_dev = rand_dev[cat1r]
  ) %>% 
  mutate(
    group = rbinom(n(), 1
                   , plogis(
                      x$beta_base + 
                      x$beta_cat1f_delta * cat1f + 
                      con1f * x$beta_con1f_delta
                     )
                   ) + 1
    , titer = rnorm(n()
                    , mu_vec[group] + 
                      (cat2f * x$theta_cat2f_mu * (group - 1)) +
                      cat1r_dev
                    , sd_vec[group]
    )
    , mfi   = logit2(x$logit_1, x$logit_2, x$logit_3, titer)
  ) %>% mutate(
      param_set = x$param_set
    , sim_num   = x$sim_num
    , .before   = 1
    )

}) %>% do.call("rbind", .)
  
} else {
  stop("Complexity option not yet supported")
}
  
}

## Modification of simulation function specifically for simulations for publication
simulate_data_for_pub   <- function(param_sets, complexity) {
  
  param_sets %<>% split_tibble(., c("param_set", "sim_num"))
  
    lapply(param_sets, FUN = function(x) {
      
      mu_vec   <- with(x, c(mu_neg, mu_pos))
      sd_vec   <- with(x, c(sd_neg, sd_pos))
      
      data.frame(
          cat1f  = with(x, rbinom(n_samps, 1, cat1f_prop))
        , cat2f  = with(x, rbinom(n_samps, 1, cat2f_prop))
      ) %>% mutate(
          con1f  = with(x, rnorm(n_samps, 0, con1f_sd))
      #  , con2f  = with(x, rnorm(n_samps, 0, con2f_sd))
      ) %>% 
        mutate(
          group = rbinom(n(), 1
                         , plogis(
                           x$beta_base + 
                             x$beta_cat1f_delta * cat1f +
                             x$beta_cat2f_delta * cat2f +
                             con1f * x$beta_con1f_delta
                         )
          ) + 1
          , titer = rnorm(n()
                          , mu_vec[group] # + (con2f * x$theta_con2f_delta * (group - 1))
                          , sd_vec[group]
          )
          , mfi   = logit2(x$logit_1, x$logit_2, x$logit_3, titer)
        ) %>% mutate(
          param_set = x$param_set
          , sim_num   = x$sim_num
          , .before   = 1
        )
      
    }) %>% do.call("rbind", .)
    
}

## Calculate skew of the simulated distribution
calc_sim_summaries      <- function(simulated_data, param_sets) {
  
  ## figure out some summaries of these distributions
  simulated_data.s <- simulated_data %>% 
    group_by(param_set, sim_num, log_mfi) %>%
    arrange(mfi) %>%
    mutate(
      g1mfi         = ifelse(group == 1, mfi, NA) 
    , max_neg_mfi   = max(g1mfi, na.rm = T)
    , quant_neg_mfi = quantile(g1mfi, 0.975, na.rm = T)
    ) %>% mutate(
      overlap_max   = ifelse(mfi < max_neg_mfi, 1, 0)  
    , overlap_quant = ifelse(mfi < quant_neg_mfi, 1, 0) 
    )
  
  prop_overlap <- simulated_data.s %>% 
    filter(group == 2) %>% 
    summarize(
      prop_overlap_max   = sum(overlap_max) / n()
    , prop_overlap_quant = sum(overlap_quant) / n()
    )
  
  mean_mfi <- simulated_data %>% 
    group_by(param_set, group, log_mfi) %>%
    summarize(mv = mean(mfi)) %>% 
    mutate(group = as.factor(group)) %>%
    ungroup(group) %>%
    pivot_wider(., names_from = "group", values_from = "mv") %>%
    rename(mean_mfi_1 = `1`, mean_mfi_2 = `2`) %>%
    mutate(diff_mean = mean_mfi_2 - mean_mfi_1)
  
  skew <- simulated_data %>% 
    group_by(param_set, sim_num, group, log_mfi) %>%
    summarize(
      titer_var  = var(titer)
    , titer_skew = skewness(titer)
    , mfi_var    = var(mfi)
    , mfi_skew   = skewness(mfi)
    )
  
  trunc_sev <- simulated_data %>% 
    group_by(param_set, sim_num, group, log_mfi) %>%
    mutate(median_titer = median(titer)) %>%
    mutate(abov_med_tit = ifelse(titer > median_titer, 1, 0)) %>%
    group_by(param_set, sim_num, group, log_mfi, abov_med_tit) %>%
    summarize(range_mfi = max(mfi) - min(mfi)) %>%
    pivot_wider(names_from = abov_med_tit, values_from = range_mfi) %>%
    mutate(prop_squish = `1` / `0`) %>%
    ungroup() %>% dplyr::select(-c(`0`, `1`))
  
  trunc_sev2 <- simulated_data %>% 
    group_by(param_set, sim_num, group, log_mfi) %>%
    mutate(
      mmode = modeest::mlv(mfi, method = "meanshif")
      , abv_m = ifelse(mfi > mmode, 1, 0)
    )
  
  trunc_sev2.s <- trunc_sev2 %>%
    summarize(
      p_abv_m = sum(abv_m) / n()
    ) 
  
  simulated_data.s <- left_join(
    prop_overlap
    , mean_mfi
    , by = c("param_set", "log_mfi")) %>%
    left_join(
      .
      , skew %>% ungroup() %>% 
        dplyr::select(-sim_num) %>% 
        pivot_wider(
          .
          , id_cols = c(param_set, log_mfi)
          , names_from = group
          , values_from = c(titer_var, titer_skew, mfi_var, mfi_skew))
      , by = c("param_set", "log_mfi")
    ) %>% left_join(
      .
      , trunc_sev %>%
        dplyr::select(-sim_num) %>% 
        pivot_wider(
          .
          , id_cols = c(param_set, log_mfi)
          , names_from = group
          , values_from = prop_squish) %>%
        rename(prop_squish_1 = `1`, prop_squish_2 = `2`)
      , by = c("param_set", "log_mfi")
    ) %>%
    left_join(
      .
      , trunc_sev2.s %>% ungroup() %>% 
        dplyr::select(-sim_num) %>% 
        pivot_wider(
          .
          , id_cols = c(param_set, log_mfi)
          , names_from = group
          , values_from = p_abv_m) %>%
        rename(p_abv_m_1 = `1`, p_abv_m_2 = `2`)
      , by = c("param_set", "log_mfi")
    )
  
  return(
    simulated_data.s %>% left_join(
      .
    , param_sets %>% dplyr::select(-sim_num)
    )
    )
  
}

## Quick and dirty plot of data -- could eventually become more complicated
examine_data            <- function(simulated_data) {
  
if (n_distinct(interaction(simulated_data$param_set, simulated_data$sim_num)) <= 30) {
  
 data_plot <- simulated_data %>% mutate(
  group = as.factor(group)
, int   = interaction(cat1f, group)) %>% {
  ggplot(., aes(x = mfi)) + 
    geom_histogram(aes(colour = int, fill = int)) + 
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    facet_grid(sim_num~param_set)
}
} else {
  data_plot <- c("Too many prameters to plot data. Will only plot if n_param_sets <= 30")
}
  
  return(data_plot)
  
}
