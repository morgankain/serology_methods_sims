tar_load(models_to_fit)
tar_load(sim.data)
tar_load(sim.params)
tar_load(stan_models.l)

this_set   <- 311

sim.data %>% 
  filter(param_set == this_set) %>% 
  mutate(
    group = as.factor(group)
    , int = interaction(cat1f, group)
  ) %>% {
    ggplot(., aes(x = mfi)) + 
      geom_density(aes(colour = group, fill = group), alpha = 0.2) + 
      scale_colour_brewer(palette = "Dark2") +
      scale_fill_brewer(palette = "Dark2") +
      facet_wrap(param_set~log_mfi, scales = "free", nrow = 3)
  }

stan_models.l <- stan_models.l[c(1,3,5)]
fit_list2     <- vector("list", length(stan_models.l))
p <- 1

prior_dat <- data.frame(
    max_mfi           = 30000 
  , mu_base_prior     = 10000
  , mu_diff_prior     = 10000
  , sigma_base_prior  = 1000 
  , sigma_diff_prior  = 1000
) %>% t() %>% 
  as.data.frame() %>%
  rename(mfi = V1) %>%
  mutate(param = rownames(.), .before = 1) %>%
  mutate(log_mfi = log(mfi))

for (i in 1:length(stan_models.l)) {
  for (j in 1) {
  
    if (j == 1) {
      fit_log     <- "log_mfi"
      prior_dat.t <- prior_dat %>% dplyr::select(param, log_mfi) %>% rename(prior = log_mfi)
    } else {
      fit_log <- "mfi"
      prior_dat.t <- prior_dat %>% dplyr::select(param, mfi) %>% rename(prior = mfi)
    }
    
    simulated_data.t <- sim.data %>% filter(param_set == this_set, log_mfi == fit_log)
    param_sets.t     <- sim.params %>% filter(param_set == this_set)
    model_names.t    <- stan_models.l[[i]]  
    
    this_name <- paste("stan_fit", fit_log, model_names.t$model, sep = ".") 
    
    tfit <- model_names.t$compiled_model[[1]]$sample(
      data    = list(
          N           = param_sets.t$n_samps %>% round()
        , N_cat1r     = param_sets.t$cat1r_count
        , y           = simulated_data.t$mfi
        , cat1f       = simulated_data.t$cat1f
        , cat2f       = simulated_data.t$cat2f
        , con1f       = simulated_data.t$con1f
        , max_mfi     = prior_dat.t %>% filter(param == "max_mfi") %>% pull(prior)
        
        , mu_base_prior     = prior_dat.t %>% filter(param == "mu_base_prior") %>% pull(prior)
        , mu_diff_prior     = prior_dat.t %>% filter(param == "mu_diff_prior") %>% pull(prior)
        , sigma_base_prior  = prior_dat.t %>% filter(param == "sigma_base_prior") %>% pull(prior)
        , sigma_diff_prior  = prior_dat.t %>% filter(param == "sigma_diff_prior") %>% pull(prior)
      )
      , chains  = 4
      , parallel_chains = 4
      , seed    = 483892929
      , refresh = 2000
    )
    
    fit_list2[[p]] <- tfit
    p <- p + 1
    
  }
}

for (i in 1:3) {
 
  stanfit <- rstan::read_stan_csv(fit_list2[[i]]$output_files())
  samps   <- rstan::extract(stanfit)
  
  samps_out <- with(samps, data.frame(
      beta_base         = beta_base
    , beta_cat1f_delta  = beta_cat1f_delta
    , beta_cat2f_delta  = beta_cat2f_delta
    , beta_con1f_delta  = beta_con1f_delta
  ))
  
#  mods <- data.frame(
#    mod_name = rep(lapply(stan_models.l, FUN = function(x)x$model) %>% unlist(), each = 2)
#  ) %>% mutate(
#    log_mfi = rep(c("log", "mfi"), 5)
#  )
  
  mods <- data.frame(
    mod_name = rep(lapply(stan_models.l, FUN = function(x)x$model) %>% unlist())
  ) %>% mutate(
    log_mfi = rep("log", 3)
  )
  
t.o <- param_sets.t %>% 
    dplyr::select(param_set, beta_base, beta_cat1f_delta, beta_cat2f_delta, beta_con1f_delta) %>% 
    pivot_longer(-param_set, values_to = "true") %>% left_join(
      .
      , samps_out %>% 
        mutate(samp = seq(n()), .before = 1) %>%
        pivot_longer(-samp) %>% 
        group_by(name) %>%
        summarize(
            lwr   = quantile(value, 0.025) 
          , lwr_n = quantile(value, 0.200)
          , mid   = quantile(value, 0.500)
          , upr_n = quantile(value, 0.800)
          , upr   = quantile(value, 0.975)
        ) 
    ) %>% mutate(
      method = paste(mods$mod_name[i], mods$log_mfi[i], sep = "."), .before = 1
    )
  
  if (i == 1) {
    t.f <- t.o
  } else {
    t.f <- rbind(t.f, t.o)
  }
   
}

t.f %>% {
  ggplot(., aes(mid, method)) +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, linewidth = 0.5) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0, linewidth = 1.0) +
    geom_point() +
    geom_vline(aes(xintercept = true)) +
    facet_wrap(~name) 
}

samps_out <- with(samps, data.frame(
    mu_base     = mu_base
  , mu_diff     = mu_diff
  , sigma_base  = sigma_base
  , sigma_diff  = sigma_diff
  , skew_neg    = skew_neg
  , skew_pos    = skew_pos
))

samps_out.s <- samps_out %>% 
  mutate(samp = seq(n()), .before = 1) %>%
  pivot_longer(-samp) %>% 
  group_by(name) %>%
  summarize(
      lwr   = quantile(value, 0.025) 
    , lwr_n = quantile(value, 0.200)
    , mid   = quantile(value, 0.500)
    , upr_n = quantile(value, 0.800)
    , upr   = quantile(value, 0.975)
  ) 

for (i in 1:nrow(samps_out)) {
  pred_pos_mfi <- c(
    pred_pos_mfi
  , rsn(
    10
  , samps_out[i, ]$mu_base + samps_out[i, ]$mu_diff
  , samps_out[i, ]$sigma_base + samps_out[i, ]$sigma_diff
  , samps_out[i, ]$skew_neg
  ))
}

real_est <- data.frame(
  real = sim.data %>% filter(param_set == 44, group == 2, log_mfi == "log_mfi") %>% pull(mfi)
, sim  = sample(pred_pos_mfi, length(sim.data %>% filter(param_set == 44, group == 2, log_mfi == "log_mfi") %>% pull(mfi)) * 6)
) %>% mutate(samp = seq(n()), .before = 1) %>%
  pivot_longer(-samp)

real_est %>% {
  ggplot(., aes(x = value)) + 
    geom_density(aes(colour = name, fill = name), alpha = 0.2) + 
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")
}



