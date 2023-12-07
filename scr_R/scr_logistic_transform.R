## NOTE: !! The titer that you get out of the stan model (when you are fitting not simulating
 ## and therefore do not know "real" titer, will be on the same scale as the concentrations
  ## in the logistic curve that you send in. These concentrations can be scaled however you want,
   ## but this makes the titer that gets estimated entirely arbitrary except that they will match
    ## the values of the concentrations)


### Notes from exploration
 ## 1) Basically what happens when you try to use priors on the logistic and fit clusters on latent titer
 ## is that titer as a free parameter pulls the logistic such that the clusters of titer match the clusters
 ## in MFI as close as possible. This results in alright recovery of many parameters, but not underlying titer at all
  ## -- This results in lower probabilities assigned to true positives than the model on MFI
 ## 2) It doesn't seem to really matter the parameterization of the logistic. When bounds are placed on the
 ## parameters to really force the logistic into the structure you want from the priors, the model starts
 ## to get divergent transitions and the probabilities of assignment for individuals starts to get really
 ## wonky with huge errorbars
 ## 3) It is really hard to think about constraints of any type, because the whole point is to allow the model
 ## to determine what the distributions are. The model already does well at the lower mean titer. Hard because
 ## the titer is in the context of the calibration concentrations. Fine if they are arbitrary, but they are 
 ## still pegged together. How to specify where the negative group in the empirical data lies?

#### ---- Logistic functions - both directions ---- ####
logit    <- function(L, b, k, x) {
  L / (
    1 + b * exp(-k*x)
  )
}
logit2   <- function(L, b, k, x) {
  L / (
    1 + exp(-k*(x - b))
  )
}
i_logit  <- function(L, b, k, y) {
  log((L - y) / (y * b)) / -k
}
i_logit2 <- function(L, b, k, y) {
  (log((L - y) / y) / -k) + b
}

## First, a quick exploration of log vs non-log mfi values
test_conc <- seq(-7, 5, by = 0.01)
test_dat <- data.frame(
  conc = test_conc
, mfi = logit2(30000, -1, 1, test_conc)
) %>% mutate(
  l_conc = dplyr::lag(conc, 1)
, l_mfi  = dplyr::lag(mfi, 1)
) %>% mutate(
  d_conc = l_conc - conc
, d_mfi  = l_mfi - mfi
) %>% mutate(
  slo = d_mfi / d_conc
)

test_dat %>% {
  ggplot(., aes(conc, mfi)) + geom_point() +
    geom_abline(slope = 1, intercept = 10)
}

test_dat %>% {
  ggplot(., aes(conc, log(mfi))) + geom_point() +
    geom_abline(slope = 1, intercept = 10)
}

test_dat %>% {
  ggplot(., aes(conc %>% exp(), mfi)) + geom_point()
}

plot(((test_conc / max(test_conc)) + 1) / 2, logit2(1, -1, 1, test_conc))

logit2(3, -1, 1, 4)
i_logit2(3, -1, 1, logit2(3, -1, 1, 5))

ttt <- read.csv("data/Rabbit_sera_filo_controls.csv") %>% 
  dplyr::select(Rep, Conc, MFI, Virus) %>% 
  mutate(Rep = as.factor(Rep)) 
  
ttt %>% {
    ggplot(., aes(Conc, MFI)) + geom_line(aes(group = interaction(Rep, Virus), colour = Virus)) + 
      scale_x_log10() #+
     # geom_vline(xintercept = .5, linetype = "dashed", linewidth = 0.5)
}

ttt %>% filter(Conc < 5, Conc > 0.1) %>% {
  ggplot(., aes(Conc, MFI)) + geom_line(aes(group = interaction(Rep, Virus), colour = Virus)) + 
    scale_x_log10() + scale_y_log10()
}

ttt %>% filter(Conc < 5) %>% {
  ggplot(., aes(Conc, MFI)) + geom_line(aes(group = interaction(Rep, Virus), colour = Virus)) 
}

base_mean <- ttt %>% filter(Conc == 0.05) %>% 
  group_by(Virus) %>% summarize(mm = mean(MFI))

test_samps <- expand.grid(
  Virus  = unique(base_mean$Virus)
, values = seq(100)
) %>% left_join(
  base_mean
) %>% mutate(
  samp_val = rnorm(n(), mm, mm / 10)
)

test_samps %>% {
  ggplot(., aes(x = samp_val)) + 
    geom_density(aes(fill = Virus)) +
    xlab("Sample Positive MFI value") +
    ylab("Density")
} 

#### ---- Take a look at sera data ---- ####
sera <- read.csv("data/Rabbit_sera_filo_controls.csv") %>% 
  dplyr::select(!contains("X")) %>%
  rename(Control = Virus) %>%
  filter(Control == "EBOV") %>% 
  rename(EBOV = MFI) %>%
  mutate(EBOV = as.numeric(EBOV)) %>%
  dplyr::select(Rep, Conc, Control, EBOV) %>% 
  ## Some modifications first to ease mle2
  mutate(
    EBOV = EBOV / 10000
  , Conc = log(Conc)
  #, Conc = as.factor(Conc) %>% as.numeric()
  #, Conc = scale(Conc) * 2
  ) %>%
  mutate(
   test_fit_vals = (
     logit(3, 0.5, 3, Conc)
   )
  )

sera %>% {#mutate(Conc = log(Conc)) %>% {
  #mutate(Conc = Conc - min(Conc)) %>% {
  ggplot(., aes(Conc, EBOV)) + 
    geom_point(aes(group = Rep)) 
}

sera %>% {
  ggplot(., aes(Conc, EBOV)) + 
    geom_line(aes(group = Rep)) 
}

sera %>% {
  ggplot(., aes(EBOV, Conc)) + 
    geom_line(aes(group = Rep)) 
}

#### ---- Fit logistic curve with bbmle::mle2 ---- ####
m1 <- bbmle::mle2(
  EBOV ~ dlnorm(
    meanlog = log(logit(L, b, k, Conc))
  , sdlog   = exp(logsdlog)
  )
, method = "L-BFGS-B"
, lower  = c(L = 1e-2, b = 1e-2, k = 1e-2, logsdlog = -10)
, start  = list(L = 1, b = 1, k = 1, logsdlog = -3)
, data   = sera #%>% mutate(EBOV = EBOV / 10000)
)
m2 <- bbmle::mle2(
  EBOV ~ dlnorm(
    meanlog = log(logit2(L, b, k, Conc))
  , sdlog   = exp(logsdlog)
  )
, method = "L-BFGS-B"
, lower  = c(L = 1e-2, b = -5, k = 1e-2, logsdlog = -10)
, start  = list(L = 1, b = 1, k = 1, logsdlog = -3)
, data   = sera #%>% mutate(EBOV = EBOV / 10000)
)

point_est1 <- summary(m1)@coef[, 1]
point_est2 <- summary(m2)@coef[, 1]

#### ---- Examine fitted relationship ---- ####
sera %<>% mutate(
  mle2_fit_vals1 = logit(point_est1[1], point_est1[2], point_est1[3], Conc)
, mle2_fit_vals2 = logit2(point_est2[1], point_est2[2], point_est2[3], Conc)
)

sera %>% {
  ggplot(., aes(Conc, EBOV)) + 
    geom_line(aes(group = Rep)) +
    geom_line(aes(Conc, mle2_fit_vals1), colour = "firebrick3", linewidth = 2) +
    geom_line(aes(Conc, mle2_fit_vals2), colour = "dodgerblue2", linewidth = 1)
}

#### ---- Simulate some data for "equivalent" model in stan ---- ####
mu_vec   <- c(-1, 0.5)
sd_vec   <- c(0.25, 1)
base_pos <- -1
beta_cat1f_delta <- 0.75

titer.out <- data.frame(
 cat1f = rbinom(1000, 1, 0.5)
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

titer.out %>% mutate(group = as.factor(group)) %>% {ggplot(., aes(x = MFI)) + geom_density(aes(fill = group))}
titer.out %>% mutate(group = as.factor(group)) %>% {ggplot(., aes(x = titer)) + geom_density(aes(fill = group))}
titer.out %>% {ggplot(., aes(titer, titer_inv_test)) + geom_point()}

## transformed values to MFI
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

#### ---- Explore stan model of just logistic curve ---- ####
stan_model <- cmdstanr::cmdstan_model("stan_models/logistic_curve2.stan", pedantic = FALSE)

stan_fit <- stan_model$sample(
  data    = list(
    N        = nrow(titer.out)
  , MFI      = titer.out$MFI
  , N_cal    = nrow(sera)
  , MFI_cal  = sera$EBOV
  , conc_cal = sera$Conc
  , min_plat = max(c(titer.out$MFI, sera$EBOV))
  )
, chains  = 4
, parallel_chains = 1
, seed    = 483892929
, refresh = 2000
)

stanfit <- rstan::read_stan_csv(stan_fit$output_files())
samps   <- rstan::extract(stanfit)

## values for use later
logistic_curve_priors <- c(
    cal_L_prior_m  = fitdistr(samps$cal_L, "normal")[[1]][1]
  , cal_L_prior_sd = fitdistr(samps$cal_L, "normal")[[1]][2]
  , cal_b_prior_m  = fitdistr(samps$cal_b, "normal")[[1]][1]
  , cal_b_prior_sd = fitdistr(samps$cal_b, "normal")[[1]][2]
  , cal_k_prior_m  = fitdistr(samps$cal_k, "normal")[[1]][1]
  , cal_k_prior_sd = fitdistr(samps$cal_k, "normal")[[1]][2]
  , cal_var_m      = fitdistr(samps$cal_var, "normal")[[1]][2]
  , cal_var_sd     = fitdistr(samps$cal_var, "normal")[[1]][2]  
)

sera %<>% mutate(
  stan_pred_mfi   = samps$titer_vals %>% colMeans()
, stan_pred_mfi_v = apply(samps$titer_vals, 2, var)
)
titer.out %<>% mutate(
  stan_pred_titer = samps$titer %>% colMeans(na.rm = T)
, pred_titer_v    = apply(samps$titer, 2, var)
)

sera %>% {
  ggplot(., aes(Conc, EBOV)) + 
    geom_line(aes(group = Rep)) +
    geom_line(aes(Conc, stan_pred_mfi), colour = "firebrick3", linewidth = 2) +
    geom_line(aes(Conc, mle2_fit_vals2), colour = "dodgerblue4", linewidth = 2)
}

sera %>% {
  ggplot(., aes(EBOV, stan_pred_mfi_v)) + geom_point()
}

titer.out %>% 
  dplyr::select(group, titer, stan_pred_titer) %>%
  pivot_longer(-group) %>%
  mutate(group = as.factor(group)) %>% {
    ggplot(., aes(x = value)) + 
      geom_density(aes(fill = group)) +
      facet_wrap(~name, ncol = 1)
  }

titer.out %>% {
  ggplot(., aes(MFI, pred_titer_v)) + geom_point()
}

titer.out %>% mutate(group = as.factor(group)) %>% {ggplot(., aes(x = MFI)) + geom_density(aes(fill = group))}

titer.out %>% {
  ggplot(., aes(titer, stan_pred_titer)) + geom_point() +
    geom_abline(intercept = 0, slope = 1)
}

#### ---- First fit simple cluster model just on the MFI values ---- #### 
stan_model <- cmdstanr::cmdstan_model("stan_models/cluster_regression_with_beta_ln_1.stan", pedantic = FALSE)

stan_fit <- stan_model$sample(
  data    = list(
    N        = nrow(titer.out)
  , y        = titer.out$MFI
  , cat1f    = titer.out$cat1f
  )
, chains  = 4
, parallel_chains = 1
, seed    = 483892921
, refresh = 2000
)

stanfit <- rstan::read_stan_csv(stan_fit$output_files())
samps   <- rstan::extract(stanfit)

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
  , cat1f    = titer.out$cat1f
  , min_plat = max(c(titer.out$MFI, sera$EBOV))
  , cal_L_prior_m  = logistic_curve_priors[1]
  , cal_L_prior_sd = logistic_curve_priors[2]
  , cal_b_prior_m  = logistic_curve_priors[3]
  , cal_b_prior_sd = logistic_curve_priors[4]
  , cal_k_prior_m  = logistic_curve_priors[5]
  , cal_k_prior_sd = logistic_curve_priors[6]
  )
, chains  = 4
, parallel_chains = 1
, seed    = 4838929
, refresh = 2000
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

sera %<>% mutate(
 cluster_log_pred = logit2(mean(samps$cal_L), mean(samps$cal_b), mean(samps$cal_k), Conc)
)

sera %>% {
  ggplot(., aes(Conc, EBOV)) + 
    geom_line(aes(group = Rep)) +
    geom_line(aes(Conc, cluster_log_pred), colour = "firebrick3", linewidth = 2) +
    geom_line(aes(Conc, mle2_fit_vals2), colour = "dodgerblue2", linewidth = 2)
}

titer.out %<>% mutate(
    titer_pred_stan_prior = samps$titer %>% colMeans()
  )

titer.out %>% {
  ggplot(., aes(titer, titer_pred_stan_prior)) + geom_point()
}
titer.out %>% mutate(group = as.factor(group)) %>% {
  ggplot(., aes(x = MFI)) + geom_density(aes(fill = group))
}
titer.out %>% mutate(group = as.factor(group)) %>% {
  ggplot(., aes(x = titer_pred_stan_prior)) + geom_density(aes(fill = group))
}
titer.out %>% mutate(group = as.factor(group)) %>% {
  ggplot(., aes(x = titer)) + geom_density(aes(fill = group))
}

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

#### ---- See how we are doing up until this point ---- ####

sim_vals %>% {
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

ind_prob.for.gg %>% mutate(value = as.factor(value)) %>% {
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

ind_prob.for.gg %>% {
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

pop_seropos.f %>% {
  ggplot(., aes(mid, stan_model)) + 
    geom_point() +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) +
    geom_vline(aes(xintercept = true)) +
    xlab("Population Seropositivity Estimate") +
    ylab("Model")
}

#### ---- Explore stan cluster model with internal logistic curve fit ---- ####
stan_model <- cmdstanr::cmdstan_model("stan_models/cluster_regression_with_beta_ln_1_titer.stan", pedantic = FALSE)

stan_fit <- stan_model$sample(
  data    = list(
    N        = nrow(titer.out)
  , MFI      = titer.out$MFI 
  , N_cal    = nrow(sera)
  , cat1f    = titer.out$cat1f
  , MFI_cal  = sera$EBOV
  , conc_cal = sera$Conc
  , min_plat = max(c(titer.out$MFI, sera$EBOV))
  )
, chains  = 4
, parallel_chains = 4
, max_treedepth = 14
, adapt_delta = 0.95
, seed    = 483892927
, refresh = 2000
)

stanfit <- rstan::read_stan_csv(stan_fit$output_files())
samps   <- rstan::extract(stanfit)

sera %<>% mutate(
  stan_pred_mfi2   = samps$titer_vals %>% colMeans()
)
titer.out %<>% mutate(
  stan_pred_titer2 = samps$titer %>% colMeans(na.rm = T)
)

titer.out %>% dplyr::select(group, titer, stan_pred_titer2) %>%
  pivot_longer(-group) %>%
  mutate(group = as.factor(group)) %>% {
    ggplot(., aes(x = value)) + geom_density(aes(fill = group)) +
      facet_wrap(~name, ncol = 1)
  }

titer.out %>% mutate(group = as.factor(group)) %>% {ggplot(., aes(x = MFI)) + geom_density(aes(fill = group))}

titer.out %>% {
  ggplot(., aes(titer, stan_pred_titer2)) + geom_point()
}

titer.out %>% {
  ggplot(., aes(stan_pred_titer1, stan_pred_titer2)) + geom_point()
}

plot(samps$beta_cat1f_delta)

samps$beta_cat1f_delta %>% hist(breaks = 50)

mu_vec
sd_vec 


sera %<>% mutate(
  stan_pred_mfi   = samps$titer_vals %>% colMeans()
)
titer.out %<>% mutate(
  stan_pred_titer = samps$titer %>% colMeans(na.rm = T)
)

sera %>% {
  ggplot(., aes(Conc, EBOV / 10000)) + 
    geom_line(aes(group = Rep)) +
    geom_line(aes(Conc, stan_pred_mfi), colour = "firebrick3", linewidth = 2)
}

titer.out %>% dplyr::select(group, titer, stan_pred_titer) %>%
  pivot_longer(-group) %>%
  mutate(group = as.factor(group)) %>% {
    ggplot(., aes(x = value)) + geom_density(aes(fill = group)) +
      facet_wrap(~name, ncol = 1)
  }

titer.out %>% mutate(group = as.factor(group)) %>% {ggplot(., aes(x = MFI)) + geom_density(aes(fill = group))}

titer.out %>% {
  ggplot(., aes(titer, stan_pred_titer)) + geom_point()
}




