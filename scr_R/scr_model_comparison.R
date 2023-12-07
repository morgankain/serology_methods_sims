simulated_data <- simulated_data[[1]]
model_name     <- model_names[[3]]
param_set      <- param_sets %>% filter(
  param_set == unique(simulated_data$param_set)
, sim_num   == unique(simulated_data$sim_num)
)

stan_fit <- stan(
  file    = "stan_models/cluster_regression_with_beta_1.stan"
, data    = list(
   N           = param_set$n_samps
 , y           = simulated_data$mfi
 , cat1f       = simulated_data$cat1f
 , cat1f_index = simulated_data$cat1f + 1
 )
, pars    = c("membership_l", "ind_sero", "log_beta", "beta_vec") 
, include = FALSE
, chains  = 4
, seed    = 483892929
, refresh = 2000
, cores   = 1
)

stan_fit <- stan(
  file    = "stan_models/cluster_regression_with_beta_theta_2.stan"
, data    = list(
   N           = param_set$n_samps
 , N_cat1r     = param_set$cat1r_count
 , y           = simulated_data$mfi
 , cat1f       = simulated_data$cat1f
 , cat2f       = simulated_data$cat2f
 , con1f       = simulated_data$con1f
 , cat1r       = simulated_data$cat1r
 )
, pars    = c("membership_l", "ind_sero", "log_beta", "beta_vec") 
, include = FALSE
, chains  = 4
, seed    = 483892929
, refresh = 2000
, cores   = 1
)

stan_fit2 <- stan(
  file    = "stan_models/cluster_regression_with_beta_theta_1.stan"
, data    = list(
   N           = param_set$n_samps
 , y           = simulated_data$mfi
 , cat1f       = simulated_data$cat1f
 , cat2f       = simulated_data$cat2f
 )
, pars    = c("membership_l", "ind_sero", "log_beta", "beta_vec"
            , "theta_cat1r_eps", ) 
, include = FALSE
, chains  = 4
, seed    = 483892929
, refresh = 2000
, cores   = 1
)

simulated_data %>% 
  group_by(group) %>% 
  summarize(n_entry = n())

e1 <- extract(stan_fit)
e2 <- extract(stan_fit2)

data.frame(
   mod1 = e1$beta_cat1f_delta
 , mod2 = e2$beta_cat1f_delta
 ) %>% mutate(
   index = seq(n())
 ) %>% pivot_longer(
   -index
 ) %>% {
   ggplot(., aes(x = value)) + 
     geom_histogram() +
     facet_wrap(~name, ncol = 1)
 }

param_set$mu_neg



