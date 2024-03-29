#simulated_data <- simulated_data[[16]]
simulated_data <- simulated_data[[2]]
#simulated_data <- simulated_data[[177]]
#simulated_data <- simulated_data[[180]]
model_name     <- model_names[[2]] %>% pull(model_base_names)
param_set      <- param_sets %>% filter(
    param_set == unique(simulated_data$param_set)
  , sim_num   == unique(simulated_data$sim_num)
)

simulated_data %>% 
  mutate(group = as.factor(group)) %>% {
  ggplot(., aes(x = mfi)) + 
    geom_density(aes(colour = group, fill = group), alpha = 0.2) + 
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    facet_wrap(~log_mfi, scales = "free", nrow = 3)
  }

simulated_data %>% group_by(group) %>% summarize(n(), mean(mfi))

if (model_name == "publication_model_2.stan") {
  if (unique(simulated_data$log_mfi) == "mfi") {
    model_name <- "publication_model_lognormal_2.stan"
  } else if (unique(simulated_data$log_mfi) == "log_mfi") {
    model_name <- "publication_model_normal_2.stan"
  } else {
    stop("Unknown model name / mfi / log mfi combo")
  }
} else if (model_name == "publication_model_skew_2.stan") {
  if (unique(simulated_data$log_mfi) == "mfi") {
    model_name <- "publication_model_skew_lognormal_2.stan"
  } else if (unique(simulated_data$log_mfi) == "log_mfi") {
    model_name <- "publication_model_skew_normal_2.stan"
  } else {
    stop("Unknown model name / mfi / log mfi combo")
  }
} else {
  stop("Unknown model name / mfi / log mfi combo")
}

model_name <- "publication_model_skew_normal_wf_2.stan" 
this_model  <- compiled_models %>% filter(model == model_name) %>% pull(compiled_model)
stan_priors <- build_stan_priors(
    simulated_data = simulated_data
  , skew_fit = ifelse(grepl("skew", model_name), T, F)
  , fit_attempt = 1
)

stan_fit    <- this_model[[1]]$sample(
  data    = list(
      N           = length(simulated_data$mfi)
    , N_cat1r     = param_set$cat1r_count
    , y           = simulated_data$mfi
    , cat1f       = simulated_data$cat1f
    , cat2f       = simulated_data$cat2f
    , con1f       = simulated_data$con1f
    , beta_base_prior_m   = stan_priors$priors %>% filter(param == "beta_base_prior_m") %>% pull(prior)
    , beta_base_prior_v   = stan_priors$priors %>% filter(param == "beta_base_prior_v") %>% pull(prior)
    , mu_base_prior_m     = stan_priors$priors %>% filter(param == "mu_base_prior_m") %>% pull(prior)
    , mu_diff_prior_m     = stan_priors$priors %>% filter(param == "mu_diff_prior_m") %>% pull(prior)
    , sigma_base_prior_m  = stan_priors$priors %>% filter(param == "sigma_base_prior_m") %>% pull(prior)
    , sigma_diff_prior_m  = stan_priors$priors %>% filter(param == "sigma_diff_prior_m") %>% pull(prior)
    , mu_base_prior_v     = stan_priors$priors %>% filter(param == "mu_base_prior_v") %>% pull(prior)
    , mu_pos_prior_v      = stan_priors$priors %>% filter(param == "mu_pos_prior_v") %>% pull(prior)
    , sigma_base_prior_v  = stan_priors$priors %>% filter(param == "sigma_base_prior_v") %>% pull(prior)
    , sigma_diff_prior_v  = stan_priors$priors %>% filter(param == "sigma_diff_prior_v") %>% pull(prior)
    ## Only used in skew model
    , skew_pos_prior_m  = min(simulated_data %>% filter(group == 2) %>% pull(mfi) %>% skewness(), -1)
    , skew_pos_prior_v  = 3
    , skew_neg_prior_m  = 0
    , skew_neg_prior_v  = 0.5
  )
  , init            = stan_priors$starting_conditions
  , chains          = 4
  , parallel_chains = 4
  , max_treedepth   = 12
  , iter_warmup     = 2000
  , iter_sampling   = 1000
  , adapt_delta     = 0.96
  , seed            = 483892929
  , refresh         = 500
)

stanfit     <- rstan::read_stan_csv(stan_fit$output_files())
samps       <- rstan::extract(stanfit)
stansummary <- summary(stanfit)
lwr_upr     <- quantile(samps$beta_base, c(0.025, 0.975))
(max(stansummary$summary[, 10], na.rm = T) > 2) | ((param_set$beta_base < lwr_upr[1]) | (param_set$beta_base > lwr_upr[2]))

samps$beta_base %>% quantile(., c(0.025, 0.50, 0.975))
param_set$beta_base
hist(samps$beta_base, breaks = 100)
hist(samps$mu[, 2], breaks = 100)
hist(samps$pop_sero, breaks = 100)
hist(samps$skew_pos, breaks = 100)

plot(
  samps$beta_base
, samps$mu[, 2]
)

plot(
  samps$mu[, 1]
, samps$mu[, 2]
)

plot(
    samps$beta_base
  , samps$beta_cat1f_delta
)

plot(
  samps$beta_base
, samps$skew_pos
)

simulated_data %>% 
  group_by(group) %>%
  summarize(mean(mfi), n())

plot(
    samps$beta_base
  , samps$mu[, 2]
)
abline(v = param_set$beta_base)
abline(h = 9.84)

beta_base <- samps$beta_base[samps$beta_base < 0]
mu1 <- samps$mu[, 1][samps$beta_base < 0]
mu2 <- samps$mu[, 2][samps$beta_base < 0]
sig1 <- samps$sigma[, 1][samps$beta_base < 0]
sig2 <- samps$sigma[, 2][samps$beta_base < 0]

skew2 <- samps$skew_pos[samps$beta_base < 0]

param_set$pop_sero

left_dist <- data.frame(
  samps$mu[, 1]
, samps$sigma_base
, samps$skew_neg
) %>% as.matrix() %>% 
  apply(., 1, FUN = function(x) {
    sn::rsn(1, x[1], x[2], x[3]) 
  }) %>% as.data.frame( ) %>% mutate(
    dist = "negative"
  )

right_dist <- data.frame(
    mu2
  , sig2
  , skew2
) %>% as.matrix() %>% 
  apply(., 1, FUN = function(x) {
    sn::rsn(1, x[1], x[2], x[3]) 
  }) %>% as.data.frame() %>% mutate(
    dist = "positive"
  )

left_dist <- data.frame(
    mu1
  , sig1
) %>% as.matrix() %>% 
  apply(., 1, FUN = function(x) {
    rnorm(1, x[1], x[2]) 
  }) %>% as.data.frame( ) %>% mutate(
    dist = "negative"
  )

left_dist <- data.frame(
    mu1
  , sig1
) %>% as.matrix() %>% 
  apply(., 1, FUN = function(x) {
    rlnorm(1, x[1], x[2]) 
  }) %>% as.data.frame( ) %>% mutate(
    dist = "negative"
  )

right_dist <- data.frame(
    mu2
  , sig2
) %>% as.matrix() %>% 
  apply(., 1, FUN = function(x) {
    rnorm(1, x[1], x[2]) 
  }) %>% as.data.frame( ) %>% mutate(
    dist = "positive"
  )

rbind(left_dist, right_dist) %>% rename(val = 1) %>% {
  ggplot(., aes(x = val)) + 
    geom_density(aes(colour = dist, fill = dist), alpha = 0.2) + 
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    scale_x_log10()
}

hist(, breaks = 100)
hist(sn::rsn(1, samps$mu[2, ], samps$sigma_base + samps$sigma_diff, samps$skew_pos), breaks = 100)



fits$sim.data %>%
  group_by(log_mfi, group, param_set) %>%
  summarize(mv = var(mfi)) %>%
  pivot_wider(values_from = mv, names_from = group) %>%
  mutate(pop_dif = `2` / `1`) %>% {
    ggplot(., aes(x = pop_dif)) + 
      geom_density() +
      facet_wrap(~log_mfi, scales = "free")
  }

fits$sim.data %>%
  group_by(log_mfi, group, param_set) %>%
  summarize(mv = sd(mfi)) %>%
  pivot_wider(values_from = mv, names_from = group) %>%
  mutate(pop_dif = `2` / `1`) %>% 
  ungroup() %>%
  group_by(log_mfi) %>%
  summarize(
    lwr   = quantile(pop_dif, 0.025)
  , lwr_n = quantile(pop_dif, 0.200)
  , mid   = quantile(pop_dif, 0.500)
  , upr_n = quantile(pop_dif, 0.800)
  , upr   = quantile(pop_dif, 0.975)
  )

fits$sim.data %>%
  group_by(log_mfi, group, param_set) %>%
  summarize(mv = sd(mfi)) %>%
  ungroup() %>%
  group_by(log_mfi, group) %>%
  summarize(
      lwr   = quantile(mv, 0.025)
    , lwr_n = quantile(mv, 0.200)
    , mid   = quantile(mv, 0.500)
    , upr_n = quantile(mv, 0.800)
    , upr   = quantile(mv, 0.975)
  )

fits$sim.data %>%
  group_by(log_mfi, group, param_set) %>%
  summarize(mm = mean(mfi), mv = var(mfi)) %>%
  mutate(group = as.factor(group)) %>% {
    ggplot(., aes(x = mm)) + 
      geom_density(aes(colour = group, fill = group), alpha = 0.2) +
      facet_wrap(~log_mfi, scales = "free")
  }

