## Simulate data for all parameter sets
simulate_data           <- function(param_sets) {
  
param_sets %<>% split_tibble(., c("param_set", "sim_num"))
  
lapply(param_sets, FUN = function(x) {
  mu_vec     <- with(x, c(mu_neg, mu_pos))
  sd_vec     <- with(x, c(sd_neg, sd_pos))
data.frame(
    age   = with(x, rbinom(n_samps, 1, age_prop))
  ) %>% mutate(
    group = rbinom(n(), 1, x$beta_base + x$beta_age_delta * abs(age - 1)) + 1
  , mfi   = rnorm(n(), mu_vec[group] + (age * x$mu_theta_age * (group - 1)), sd_vec[group])
  ) %>% mutate(
      param_set = x$param_set
    , sim_num   = x$sim_num
    , .before   = 1
    )
}) %>% do.call("rbind", .)
  
}

## Quick and dirty plot of data -- could eventually become more complicated
examine_data            <- function(simulated_data) {
  
if (n_distinct(interaction(simulated_data$param_set, simulated_data$sim_num)) <= 30) {
  
 data_plot <- simulated_data %>% mutate(
  group = as.factor(group)
, int   = interaction(age, group)) %>% {
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

