## organize regressions for downstream analyses
sort_regression         <- function(fitted_regressions, param_sets, complexity) {

  param_sets.l <- param_sets %>% split_tibble(., "param_set")
  
regression.pred <- purrr::pmap(list(fitted_regressions, param_sets.l), .f = function(x, y) {
  
  param_sets_sims.l <- y %>% split_tibble(., "sim_num")
  
  regression_sims.pred <- purrr::pmap(list(x, param_sets_sims.l), .f = function(v, w) {
    
  if (complexity == 1) {
    
  true_vals <- w %>%
   dplyr::select(param_set, sim_num, beta_base, beta_cat1f_delta) %>% 
   mutate(
       beta_cat1f = plogis(beta_base + beta_cat1f_delta)
     , beta_base  = plogis(beta_base)
     ) %>%
   dplyr::select(-beta_cat1f_delta) %>%
   pivot_longer(-c(param_set, sim_num), values_to = "true")
 
    pred.out <- predictorEffect("cat1f", v[[1]]) %>% summary()
    pred.out <- with(pred.out, data.frame(
      lwr = lower
    , mid = effect
    , upr = upper
    )) %>% mutate(
      name = c("beta_base", "beta_cat1f")
    , .before = 1
    )
    
    out1 <- left_join(true_vals, pred.out, by = "name") %>% mutate(model = "no_variance")
    
    pred.out <- predictorEffect("cat1f", v[[2]]) %>% summary()
    pred.out <- with(pred.out, data.frame(
      lwr = lower
    , mid = effect
    , upr = upper
    )) %>% mutate(
      name = c("beta_base", "beta_cat1f")
    , .before = 1
    )
    
    out2 <- left_join(true_vals, pred.out, by = "name") %>% mutate(model = "positive_probability")
    
  } else if (complexity == 2) {
    
   true_vals <- w %>%
     dplyr::select(param_set, sim_num, beta_base, beta_cat1f_delta, beta_con1f_delta) %>% 
     pivot_longer(-c(param_set, sim_num), values_to = "true")
 
   ci_attempt <- try(
     {
     confint(v[[1]])
     }, silent = T
   )
   
   if (any(class(ci_attempt) == "try-error")) {
     pred.out <- data.frame(
        lwr = rep(NA, 3)
      , mid = rep(NA, 3)
      , upr = rep(NA, 3)
     ) %>% 
    mutate(
      name = c("beta_base", "beta_cat1f_delta", "beta_con1f_delta")
    , .before = 1
    )
   } else {
     pred.out <- ci_attempt %>% as.data.frame() %>%
      mutate(mid = coef(v[[1]])) %>% 
      rename(lwr = "2.5 %", upr = "97.5 %") %>%
      relocate(mid, .after = lwr) %>%
      mutate(
      name = c("beta_base", "beta_cat1f_delta", "beta_con1f_delta")
    , .before = 1
    )
   }
   
    out1 <- left_join(true_vals, pred.out, by = "name") %>% mutate(model = "no_variance")
    
   ci_attempt <- try(
     {
     confint(v[[2]])
     }, silent = T
   )

   if (any(class(ci_attempt) == "try-error")) {
     pred.out <- data.frame(
        lwr = rep(NA, 3)
      , mid = rep(NA, 3)
      , upr = rep(NA, 3)
     ) %>% 
    mutate(
      name = c("beta_base", "beta_cat1f_delta", "beta_con1f_delta")
    , .before = 1
    )
   } else {
    pred.out <- confint(v[[2]]) %>% as.data.frame() %>%
      mutate(mid = coef(v[[2]])) %>% 
      rename(lwr = "2.5 %", upr = "97.5 %") %>%
      relocate(mid, .after = lwr) %>%
      mutate(
      name = c("beta_base", "beta_cat1f_delta", "beta_con1f_delta")
    , .before = 1
    )
   }
    
    out2 <- left_join(true_vals, pred.out, by = "name") %>% mutate(model = "positive_probability")
  
  } else {
    stop("Model complexity not -yet- supported")
  }
    
    return(
      rbind(
        out1, out2
      )
    )
    
  }) %>% do.call("rbind", .)
    
  return(regression_sims.pred)
  
  }) %>% do.call("rbind", .)

return(regression.pred %>% mutate(
  cover  = ifelse(true > lwr & true < upr, 1, 0)
, CI_wid = upr - lwr
, m_diff = abs(mid - true)
  ))
  
}

## Summarize stan fits
summarize_stan_fits     <- function(model_fits, param_sets, simulated_data) {
  
for (i in 1:nrow(model_fits)) {
  
y <- param_sets %>% filter(
  param_set == model_fits$param_set[i]
, sim_num   == model_fits$sim_num[i]
)

x     <- model_fits[i, ]$fitted_model[[1]]
samps <- rstan::extract(x[[1]])

z <- simulated_data %>% filter(
  param_set == model_fits$param_set[i]
, sim_num   == model_fits$sim_num[i]
)

if (model_fits[i, ]$model == "cluster_regression_base_1.stan") {
  
samps_out <- with(samps, data.frame(
    mu_base    = mu_base
  , mu_pos     = mu[, 2]
  , sigma_base = sigma_base
  , sigma_pos  = sigma[, 2]
  , beta       = beta
  ))

true_vals <- y %>% 
  mutate(beta = (plogis(beta_base)*cat1f_prop + (plogis(beta_base + beta_cat1f_delta)*age_prop))) %>% 
  dplyr::select(param_set, sim_num, mu_neg, mu_pos, sd_neg, sd_pos, beta) %>% 
  pivot_longer(-c(param_set, sim_num), values_to = "true")

} else if (model_fits[i, ]$model == "cluster_regression_with_beta_1.stan") {
  
samps_out <- with(samps, data.frame(
    mu_base    = mu_base
  , mu_pos     = mu[, 2]
  , sigma_base = sigma_base
  , sigma_pos  = sigma[, 2]
  , beta_base  = beta_base
  , beta_cat1f_delta = beta_cat1f_delta
  ))

true_vals <- y %>% 
  dplyr::select(param_set, sim_num, mu_neg, mu_pos, sd_neg, sd_pos, beta_base, beta_cat1f_delta) %>% 
  pivot_longer(-c(param_set, sim_num), values_to = "true")
  
} else if (model_fits[i, ]$model == "cluster_regression_with_beta_theta_1.stan") {
  
samps_out <- with(samps, data.frame(
    mu_base    = mu_base
  , mu_pos     = mu[, 2]
  , sigma_base = sigma_base
  , sigma_pos  = sigma[, 2]
  , beta_base  = beta_base
  , beta_cat1f_delta = beta_cat1f_delta
  , theta_cat2f_mu   = theta_cat2f_mu
  ))

true_vals <- y %>% 
  dplyr::select(
    param_set, sim_num, mu_neg, mu_pos, sd_neg, sd_pos
  , beta_base, beta_cat1f_delta, theta_cat2f_mu
  ) %>% 
  pivot_longer(-c(param_set, sim_num), values_to = "true")
  
} else if (model_fits[i, ]$model == "cluster_regression_with_beta_theta_2.stan") { 

samps_out <- with(samps, data.frame(
    mu_base    = mu_base
  , mu_pos     = mu[, 2]
  , sigma_base = sigma_base
  , sigma_pos  = sigma[, 2]
  , beta_base  = beta_base
  , beta_cat1f_delta = beta_cat1f_delta
  , theta_cat2f_mu   = theta_cat2f_mu
  , theta_cat1r_sd   = theta_cat1r_sd
  , beta_con1f_delta = beta_con1f_delta
  ))

true_vals <- y %>% 
  dplyr::select(
    param_set, sim_num, mu_neg, mu_pos, sd_neg, sd_pos
  , beta_base, beta_cat1f_delta, theta_cat2f_mu
  , theta_cat1r_sd, beta_con1f_delta
  ) %>% 
  pivot_longer(-c(param_set, sim_num), values_to = "true")
    
} else {
  stop("Stan model name not known")
}

samps_out %<>% 
  mutate(samp = seq(n()), .before = 1) %>%
  pivot_longer(-samp) %>% group_by(name) %>%
  summarize(
    lwr   = quantile(value, 0.025) 
  , lwr_n = quantile(value, 0.200)
  , mid   = quantile(value, 0.500)
  , upr_n = quantile(value, 0.800)
  , upr   = quantile(value, 0.975)
  ) %>% mutate(
    name = plyr::mapvalues(
      name
    , from = c("mu_base", "sigma_base", "sigma_pos")
    , to = c("mu_neg", "sd_neg", "sd_pos"))
  ) 

out.clean.t <- left_join(true_vals, samps_out, by = "name") %>% 
  left_join(model_fits[i, ] %>% dplyr::select(model, param_set, sim_num) %>%
              mutate(param_set = as.numeric(param_set), sim_num = as.numeric(sim_num))
            , ., by = c("param_set", "sim_num")) %>%
  left_join(., y %>% dplyr::select(param_set, sim_num, n_samps), by = c("param_set", "sim_num"))

pred_pos <- samps$membership_p %>% 
  reshape2::melt() %>% 
  rename(clust = Var2, samp = Var3) %>%
  pivot_wider(values_from = value, names_from = clust) %>%
  group_by(samp) %>%
  summarize(
    lwr   = quantile(`1`, 0.025)
  , lwr_n = quantile(`1`, 0.200)
  , mid   = quantile(`1`, 0.500)
  , upr_n = quantile(`1`, 0.800)
  , upr   = quantile(`1`, 0.975)
  )

z %<>% mutate(samp = seq(n()), .before = 1) %>% 
  left_join(., pred_pos, by = "samp") %>% 
  mutate(stan_model = model_fits[i, ]$model, .before = 1)

pop_seropos <- (samps$pop_sero / y$n_samps) %>% 
  quantile(c(0.025, 0.200, 0.500, 0.800, 0.975)) %>% 
  t() %>% as.data.frame()

names(pop_seropos) <- c("lwr", "lwr_n", "mid", "upr_n", "upr")

pop_seropos %<>% 
  mutate(
    stan_model = model_fits[i, ]$model, .before = 1
  , param_set  = y$param_set
  , sim_num    = y$sim_num
  , true       = (z %>% mutate(group = group - 1) %>% pull(group) %>% sum()) / y$n_samps
  )

if (i == 1) {
  out.clean     <- out.clean.t
  z.f           <- z
  pop_seropos.f <- pop_seropos
} else {
  out.clean     <- rbind(out.clean, out.clean.t)
  z.f           <- rbind(z.f, z)
  pop_seropos.f <- rbind(pop_seropos.f, pop_seropos)
}

}
  
return(
  list(
  coef = out.clean %>% mutate(
      cover = ifelse(true > lwr & true < upr, 1, 0)
    , CI_wid = upr - lwr
    , m_diff = abs(mid - true)
  )
, group_pred   = z.f %>% mutate(group = group - 1)
, prop_seropos = pop_seropos.f
 )
)

}

## Summarize group assignment predictions
calculate_group_assignments <- function(three_sd.g, mclust.g, stan.g) {

  three_sd.g %<>%
   dplyr::select(-assigned_group) %>%
   group_by(model, param_set, sim_num, group) %>%
   summarize(
     prob = mean(V2)
   ) %>% mutate(
     quantile = "mid", .before = prob
   )
  
  mclust.g %<>%
   dplyr::select(-assigned_group) %>%
   group_by(model, param_set, sim_num, group) %>%
   summarize(
     prob = mean(V2)
   ) %>% mutate(
     quantile = "mid", .before = prob
   )
  
  stan.g %<>% 
   rename(model = stan_model) %>%
   dplyr::select(-c(samp, cat1f, cat2f, mfi)) %>%
   pivot_longer(-c(model, param_set, sim_num, group)
               , names_to = "quantile") %>%
   group_by(model, param_set, sim_num, group, quantile) %>%
   summarize(
     prob = mean(value)
   )  
  
  return(
    rbind(
      three_sd.g, mclust.g, stan.g
      )
    )
  
}

## Collate and summarize population level seropositivity estimates
calculate_population_seropositivity <- function(three_sd.g, mclust.g, stan.g) {
  
  three_sd.g %<>% 
    group_by(param_set, sim_num) %>% 
    summarize(
      true     = mean(group)
    , prop_pos = mean(assigned_group)
    ) %>% mutate(
      prop_pos_diff = prop_pos - true 
    ) %>% ungroup() %>% 
    mutate(model = "3sd", .before = 1) %>%
    mutate(quantile = "mid", .after = "true")
  
  mclust.g %<>% 
    group_by(param_set, sim_num) %>% 
    summarize(
      true     = mean(group)
    , prop_pos = sum(V2) / n()
    ) %>% mutate(
      prop_pos_diff = prop_pos - true 
    ) %>% ungroup() %>% 
    mutate(model = "mclust", .before = 1) %>%
    mutate(quantile = "mid", .after = "true")
  
  stan.g %<>% 
    pivot_longer(-c(stan_model, param_set, sim_num, true)
                 , names_to = "quantile", values_to = "prop_pos") %>%
    mutate(prop_pos_diff = prop_pos - true) %>%
    rename(model = stan_model)
    
  all.g <- rbind(three_sd.g, mclust.g, stan.g)
  
  return(all.g)
  
}

