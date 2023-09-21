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
  , mfi   = rnorm(n()
                  , mu_vec[group] + 
                    (cat2f * x$theta_cat2f_mu * (group - 1))
                  , sd_vec[group]
                  )
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
  ) %>% mutate(
    group = rbinom(n(), 1
                   , plogis(
                      x$beta_base + 
                      x$beta_cat1f_delta * cat1f + 
                      con1f * x$beta_con1f_delta
                     )
                   ) + 1
  , mfi   = rnorm(n()
                  , mu_vec[group] + 
                    (cat2f * x$theta_cat2f_mu * (group - 1)) + 
                    rand_dev[group]
                  , sd_vec[group]
                  )
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
