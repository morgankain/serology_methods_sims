#### ---- Logistic functions - both directions ---- ####
logit   <- function(L, b, k, x) {
  L / (
    1 + b * exp(-k*x)
  )
}
i_logit <- function(L, b, k, y) {
  log((L - y) / (y * b)) / -k
}

#### ---- Function to fit cluster models on MFI and latent titer and compare fits ---- ####
compare_cluster_methods <- function(
    sera_data_path = "data/Rabbit_sera_filo_controls.csv"
  , n_samps, prop_each_group
  , mu_vec, sd_vec, base_pos, beta_cat1f_delta
  , stop_at_data
  ) {
  
## Load data
sera <- read.csv(sera_data_path) %>% 
  filter(Control == "EBOV") %>% 
  mutate(EBOV = as.numeric(EBOV)) %>%
  dplyr::select(Rep, Conc, Control, EBOV) %>% 
  ## Some modifications first to ease mle2
  mutate(
    EBOV = EBOV / 10000
  , Conc = log(Conc)
 # , Conc = as.factor(Conc) %>% as.numeric()
 # , Conc = scale(Conc) * 2
  )

## Fit logistic curve to calibration curve
m1 <- bbmle::mle2(
  EBOV ~ dlnorm(
    meanlog = log(logit(L, b, k, Conc))
  , sdlog   = exp(logsdlog)
  )
, method = "L-BFGS-B"
, lower  = c(L = 1e-2, b = 1e-2, k = 1e-2, logsdlog = -10)
, start  = list(L = 1, b = 1, k = 1, logsdlog = -3)
, data   = sera 
)

point_est <- summary(m1)@coef[, 1]

## Simulate titer data and back transform it to MFI with mle2 fit
titer.out <- data.frame(
 cat1f = rbinom(n_samps, 1, 0.5)
) %>% mutate(
  group = rbinom(n(), 1
                   , plogis(
                      base_pos + 
                      beta_cat1f_delta * cat1f
                   )
                   ) + 1
, titer = rnorm(n()
                  , mu_vec[group] 
                  , sd_vec[group]
                  )
) %>% mutate(
  MFI = logit(point_est[1], point_est[2], point_est[3], titer)
) %>% mutate(
  titer_inv_test = i_logit(point_est[1], point_est[2], point_est[3], MFI)
) %>% mutate(
  ind = seq(n()), .before = 1
)

if (stop_at_data) {
  return(titer.out)
}

titer.out %>% mutate(group = as.factor(group)) %>% {ggplot(., aes(x = MFI)) + geom_density(aes(fill = group))}
titer.out %>% mutate(group = as.factor(group)) %>% {ggplot(., aes(x = titer)) + geom_density(aes(fill = group))}
titer.out %>% {ggplot(., aes(titer, titer_inv_test)) + geom_point()}

## Take summaries of the simulated data and store parameter values
mfi_params <- titer.out %>% group_by(group) %>% summarize(m_mfi = mean(MFI), sd_mfi = sd(MFI))

sim_vals <- data.frame(
    type             = "true"
  , mu_base_titer    = mu_vec[1]
  , mu_pos_titer     = mu_vec[2]
  , sigma_base_titer = sd_vec[1]
  , sigma_pos_titer  = sd_vec[2]
  , beta_base  = plogis(base_pos)
  , beta_cat1f = plogis(base_pos + beta_cat1f_delta)
  , beta_cat1f_delta = beta_cat1f_delta
  , mu_base_MFI    = mfi_params$m_mfi[1]
  , mu_pos_MFI     = mfi_params$m_mfi[2]
  , sigma_base_MFI = mfi_params$sd_mfi[1]
  , sigma_pos_MFI  = mfi_params$sd_mfi[2] - mfi_params$sd_mfi[1]
) %>% pivot_longer(-type, values_to = "mid") %>%
  mutate(
    lwr = NA, upr = NA
  )

titer.out.just_groups <- titer.out %>% dplyr::select(ind, group) %>% 
  pivot_longer(-ind) %>%
  mutate(model = "true", type = "mid", .before = 1) %>% 
  dplyr::select(-name)

## Fit stan model to calibration data and estimate MFI values from the titer values as a check on the fit
stan_model <- cmdstanr::cmdstan_model("stan_models/logistic_curve.stan", pedantic = FALSE)

stan_fit <- stan_model$sample(
  data    = list(
    N        = nrow(titer.out)
  , MFI      = titer.out$MFI
  , N_cal    = nrow(sera)
  , MFI_cal  = sera$EBOV
  , conc_cal = sera$Conc[, 1]
  , min_plat = max(c(titer.out$MFI, sera$EBOV))
  )
, chains  = 4
, parallel_chains = 1
, seed    = 483892929
, refresh = 2000
)

stanfit <- rstan::read_stan_csv(stan_fit$output_files())
samps   <- rstan::extract(stanfit)

## values for the logistic curve for use later
logistic_curve_priors <- c(
    cal_L_prior_m  = fitdistr(samps$cal_L, "normal")[[1]][1]
  , cal_L_prior_sd = fitdistr(samps$cal_L, "normal")[[1]][2]
  , cal_b_prior_m  = fitdistr(samps$cal_b, "normal")[[1]][1]
  , cal_b_prior_sd = fitdistr(samps$cal_b, "normal")[[1]][2]
  , cal_k_prior_m  = fitdistr(samps$cal_k, "normal")[[1]][1]
  , cal_k_prior_sd = fitdistr(samps$cal_k, "normal")[[1]][2]
)

#### ---- First fit simple cluster model just on the MFI values ---- #### 
stan_model <- cmdstanr::cmdstan_model("stan_models/cluster_regression_with_beta_ln_1.stan", pedantic = FALSE)

stan_fit <- stan_model$sample(
  data    = list(
    N        = nrow(titer.out)
  , y        = titer.out$MFI
  , cat1f    = titer.out$cat1f
  )
, chains  = 4
, parallel_chains = 4
, seed    = 483892929
, refresh = 2000
)

stanfit <- rstan::read_stan_csv(stan_fit$output_files())
samps   <- rstan::extract(stanfit)

## Extract the relevant parameters and summarize
samps_out <- with(samps, data.frame(
    mu_base_MFI      = mu_base
  , mu_pos_MFI       = mu[, 2]
  , sigma_base_MFI   = sigma_base
  , sigma_pos_MFI    = sigma[, 2]
  , beta_base        = plogis(beta_base)
  , beta_cat1f       = plogis(beta_base + beta_cat1f_delta)
  , beta_cat1f_delta = beta_cat1f_delta
  ))

samps_out %<>% 
  mutate(samp = seq(n()), .before = 1) %>%
  pivot_longer(-samp) %>% group_by(name) %>%
  summarize(
    lwr   = quantile(value, 0.025) 
  , mid   = quantile(value, 0.500)
  , upr   = quantile(value, 0.975)
  ) %>% mutate(
    type = "cluster on transformed MFI"
  )

sim_vals %<>% rbind(., samps_out)

pred_pos <- samps$membership_p %>% 
  reshape2::melt() %>% 
  rename(clust = Var2, ind = Var3) %>%
  pivot_wider(values_from = value, names_from = clust) %>%
  group_by(ind) %>%
  summarize(
    lwr   = quantile(`1`, 0.025)
  , mid   = quantile(`1`, 0.500)
  , upr   = quantile(`1`, 0.975)
  ) %>% pivot_longer(-ind, names_to = "type") %>%
  relocate(type, .before = 1) %>% 
  mutate(model = "cluster on transformed MFI", .before = 1)

titer.out.just_groups %<>% rbind(., pred_pos)

pop_seropos <- (samps$pop_sero / nrow(titer.out)) %>% 
  quantile(c(0.025, 0.200, 0.500, 0.800, 0.975)) %>% 
  t() %>% as.data.frame()

names(pop_seropos) <- c("lwr", "lwr_n", "mid", "upr_n", "upr")

pop_seropos %<>% 
  mutate(
    stan_model = "cluster on transformed MFI"
  , true       = (titer.out %>% mutate(group = group - 1) %>% pull(group) %>% sum()) / nrow(titer.out)
  )

pop_seropos.f <- pop_seropos

#### ---- Explore stan cluster model with logistic curve with priors from previous model (no fitting in model) ---- ####
stan_model <- cmdstanr::cmdstan_model("stan_models/cluster_regression_with_beta_ln_1_titer_just_prior.stan", pedantic = FALSE)

stan_fit <- stan_model$sample(
  data    = list(
    N        = nrow(titer.out)
  , MFI      = titer.out$MFI
  , N_cal    = nrow(sera)
  , cat1f    = titer.out$cat1f
  , MFI_cal  = sera$EBOV
  , conc_cal = sera$Conc[, 1]
  , min_plat = max(c(titer.out$MFI, sera$EBOV))
  , cal_L_prior_m  = logistic_curve_priors[1]
  , cal_L_prior_sd = logistic_curve_priors[2]
  , cal_b_prior_m  = logistic_curve_priors[3]
  , cal_b_prior_sd = logistic_curve_priors[4]
  , cal_k_prior_m  = logistic_curve_priors[5]
  , cal_k_prior_sd = logistic_curve_priors[6]
  )
, chains  = 4
, parallel_chains = 4
, seed    = 483892927
, refresh = 2000
, adapt_delta   = 0.98
, max_treedepth = 15
)

stanfit <- rstan::read_stan_csv(stan_fit$output_files())
samps   <- rstan::extract(stanfit)

samps_out <- with(samps, data.frame(
    mu_base_titer      = mu_base
  , mu_pos_titer       = mu[, 2]
  , sigma_base_titer   = sigma_base
  , sigma_pos_titer    = sigma[, 2]
  , beta_base        = plogis(beta_base)
  , beta_cat1f       = plogis(beta_base + beta_cat1f_delta)
  , beta_cat1f_delta = beta_cat1f_delta
  ))

samps_out %<>% 
  mutate(samp = seq(n()), .before = 1) %>%
  pivot_longer(-samp) %>% group_by(name) %>%
  summarize(
    lwr   = quantile(value, 0.025) 
  , mid   = quantile(value, 0.500)
  , upr   = quantile(value, 0.975)
  ) %>% mutate(
    type = "cluster on latent titer -- just prior"
  )

sim_vals %<>% rbind(., samps_out)

pred_pos <- samps$membership_p %>% 
  reshape2::melt() %>% 
  rename(clust = Var2, ind = Var3) %>%
  pivot_wider(values_from = value, names_from = clust) %>%
  group_by(ind) %>%
  summarize(
    lwr   = quantile(`1`, 0.025)
  , mid   = quantile(`1`, 0.500)
  , upr   = quantile(`1`, 0.975)
  ) %>% pivot_longer(-ind, names_to = "type") %>%
  relocate(type, .before = 1) %>% 
  mutate(model = "cluster on latent titer -- just prior", .before = 1)

titer.out.just_groups %<>% rbind(., pred_pos)

pop_seropos <- (samps$pop_sero / nrow(titer.out)) %>% 
  quantile(c(0.025, 0.200, 0.500, 0.800, 0.975)) %>% 
  t() %>% as.data.frame()

names(pop_seropos) <- c("lwr", "lwr_n", "mid", "upr_n", "upr")

pop_seropos %<>% 
  mutate(
    stan_model = "cluster on latent titer -- just prior"
  , true       = (titer.out %>% mutate(group = group - 1) %>% pull(group) %>% sum()) / nrow(titer.out)
  )

pop_seropos.f <- rbind(pop_seropos.f, pop_seropos)

return(
  list(
    sim_vals
  , titer.out.just_groups
  , pop_seropos.f
  )
)
  
}

plot_cluster_method_comparison <- function(cluster_methods.out) {
  
  sim_vals              <- cluster_methods.out[[1]]
  titer.out.just_groups <- cluster_methods.out[[2]]
  ind_prob.for.gg       <- cluster_methods.out[[3]]
  
  gg.1 <- sim_vals %>% {
  ggplot(., aes(name, mid)) + 
    geom_point(aes(colour = type)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr, colour = type)) +
    scale_color_brewer(palette = "Dark2") +
    theme(axis.text.x = element_text(size = 10, angle = 300, hjust = 0))
}

titer.out.just_groups[titer.out.just_groups$model == "true", ]$value <- titer.out.just_groups[titer.out.just_groups$model == "true", ]$value - 1

ind_prob.for.gg <- titer.out.just_groups %>% 
  filter(model != "true") %>% 
  mutate(model = plyr::mapvalues(model, from = unique(model), to = c("MFI", "Titer"))) %>%
  mutate(m.t = interaction(model, type)) %>%
  dplyr::select(-c(model, type)) %>%
  pivot_wider(names_from = m.t, values_from = value) %>% 
  left_join(., titer.out.just_groups %>% filter(model == "true"))

gg.2 <- ind_prob.for.gg %>% mutate(value = as.factor(value)) %>% {
  ggplot(., aes(MFI.mid, Titer.mid, colour = value)) + 
  geom_point() +
  geom_errorbar(aes(ymin = Titer.lwr, ymax = Titer.upr)) +
  geom_errorbarh(aes(xmin = MFI.lwr, xmax = MFI.upr)) +
  scale_colour_brewer(palette = "Dark2", name = "True Group") +
  geom_abline(aes(intercept = 0, slope = 1)) +
      xlab("Estimated Probability Positive -- Model on MFI") +
      ylab("Estimated Probability Positive -- Model on Titer")
}

ind_prob.for.gg <- titer.out.just_groups %>% 
  filter(model != "true") %>% 
  mutate(model = plyr::mapvalues(model, from = unique(model), to = c("MFI", "Titer"))) %>%
  mutate(m.t = interaction(model, type)) %>%
  dplyr::select(-c(model, type)) %>%
  filter(grepl("mid", m.t)) %>%
 # pivot_wider(names_from = m.t, values_from = value) %>% 
  left_join(., titer.out.just_groups %>% filter(model == "true") %>% rename(true = value))

gg.3 <- ind_prob.for.gg %>% {
  ggplot(., aes(true, value)) + 
  geom_jitter(width = 0.1, height = 0.1, aes(colour = m.t)) +
  xlab("True Group") +
  ylab("Estimated Probability Positive -- Model on Titer") +
  scale_colour_brewer(palette = "Dark2", name = "Model") +
  geom_line(
    data = ind_prob.for.gg %>% group_by(m.t, true) %>% summarize(mm = mean(value))
  , aes(true, mm, colour = m.t)
  )
}

gg.4 <- pop_seropos.f %>% {
  ggplot(., aes(mid, stan_model)) + 
    geom_point() +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) +
    geom_vline(aes(xintercept = true)) +
    xlab("Population Seropositivity Estimate") +
    ylab("Model")
}

return(
  list(
    gg.1, gg.2, gg.3, gg.4
  )
)

}

cluster_methods.out <- compare_cluster_methods(
    sera_data_path   = "data/Rabbit_sera_filo_controls.csv"
  , n_samps          = 1000
  , prop_each_group  = 0.5
  , mu_vec           = c(-0.5, 1.5)
  , sd_vec           = c(0.25, 1)
  , base_pos         = -1
  , beta_cat1f_delta = 0.75
  , stop_at_data     = F
)

cluster_methods.out[[1]] %>% as.data.frame()

cluster_methods.out %>% mutate(group = as.factor(group)) %>% {ggplot(., aes(x = MFI)) + geom_density(aes(fill = group))}
cluster_methods.out %>% mutate(group = as.factor(group)) %>% {ggplot(., aes(x = titer)) + geom_density(aes(fill = group))}
cluster_methods.out %>% {ggplot(., aes(titer, titer_inv_test)) + geom_point()}

cluster_plots <- plot_cluster_method_comparison(cluster_methods.out)

cluster_plots[[2]]


