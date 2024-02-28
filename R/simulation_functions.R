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
        , con2f  = with(x, rnorm(n_samps, 0, con2f_sd))
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
                          , mu_vec[group] + 
                            (con2f * x$theta_con2f_delta * (group - 1))
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
calc_skew               <- function(sim.data) {
  sim.data %>% 
    group_by(param_set, sim_num, group, log_mfi) %>%
    summarize(
      titer_var  = var(titer)
    , titer_skew = skewness(titer)
    , mfi_var    = var(mfi)
    , mfi_skew   = skewness(mfi)
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


