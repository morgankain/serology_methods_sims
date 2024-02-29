## organize regressions for downstream analyses
sort_regression         <- function(fitted_regressions, param_sets, complexity, groupings1, groupings2) {

  param_sets.l <- param_sets %>% split_tibble(., groupings1)
  
regression.pred <- purrr::pmap(list(fitted_regressions, param_sets.l), .f = function(x, y) {
  
  param_sets_sims.l <- y %>% split_tibble(., "sim_num")
  
  regression_fits <- names(x) %>% matrix() %>% apply(., 1, FUN = function(xx) {
    xx %>% strsplit(., "[.]") %>% unlist() %>% t() %>% as.data.frame() 
  }) %>% 
    do.call("rbind", .) %>% 
    mutate(index = seq(n())) %>%
    rename(log_mfi = 1, sim_num = 2, method = 3) %>% 
    mutate(sim_num = as.numeric(sim_num)) %>%
    left_join(., y, by = "sim_num") %>% 
    split_tibble(., groupings2)
  
  regression_sims.pred <- purrr::pmap(list(x, regression_fits), .f = function(v, w) {
    
    true_vals <- w %>%
      dplyr::select(param_set, sim_num, beta_base, beta_cat1f_delta, beta_con1f_delta) %>% 
      pivot_longer(-c(param_set, sim_num), values_to = "true")
    
    for (vv in 1:length(v)) {
      
      t_name <- paste("out", vv, sep = "")
      
      ci_attempt <- try(
        {
          suppressMessages(confint(v[[vv]]))
        }, silent = T
      )
      
      if (any(class(ci_attempt) == "try-error")) {
        pred.out <- data.frame(
            lwr = rep(NA, 4)
          , mid = rep(NA, 4)
          , upr = rep(NA, 4)
        ) %>% 
          mutate(
            name = c("beta_base", "beta_cat1f_delta", "beta_cat2f_delta", "beta_con1f_delta")
            , .before = 1
          )
      } else {
        if (complexity == 1) {
          pred.out <- ci_attempt %>% as.data.frame() %>%
            mutate(mid = coef(v[[vv]])) %>% 
            rename(lwr = "2.5 %", upr = "97.5 %") %>%
            relocate(mid, .after = lwr) %>%
            mutate(
              name = c("beta_base", "beta_cat1f_delta")
              , .before = 1
            )
        } else {
          pred.out <- ci_attempt %>% as.data.frame() %>%
            mutate(mid = coef(v[[vv]])) %>% 
            rename(lwr = "2.5 %", upr = "97.5 %") %>%
            relocate(mid, .after = lwr) %>%
            mutate(
              name = c("beta_base", "beta_cat1f_delta", "beta_cat2f_delta", "beta_con1f_delta")
              , .before = 1
            )
        }
      }
      
      yes_var <- ifelse(vv == 1, "no_variance", "variance")
      
      assign(t_name, left_join(
        w %>% dplyr::select(log_mfi, sim_num, param_set, method)
        , true_vals
        , by = c("param_set", "sim_num")
      ) %>% left_join(
        .
        , pred.out
        , by = "name"
      ) %>% mutate(model = yes_var, .after = "method")
      )
      
    }
    
    if (length(v) == 1) {
      return(out1)
    } else {
      return(rbind(out1, out2))
    }
    
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
summarize_stan_fits     <- function(model_fits, param_sets, simulated_data, complexity) {
  
for (i in 1:nrow(model_fits)) {
  
y <- param_sets %>% filter(
  param_set == model_fits$param_set[i]
, sim_num   == model_fits$sim_num[i]
)

samps <- model_fits[i, ]$fitted_model[[1]][[1]]

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
  mutate(beta = (plogis(beta_base)*cat1f_prop + (plogis(beta_base + beta_cat1f_delta)*(1-cat1f_prop)))) %>% 
  dplyr::select(param_set, sim_num, mu_neg, mu_pos, sd_neg, sd_pos, beta) %>% 
  pivot_longer(-c(param_set, sim_num), values_to = "base_true")

z.p <- z %>% 
  group_by(param_set, sim_num) %>%
  filter(group == 2) %>%
  summarize(
    mu_pos = mean(mfi)
    , sd_pos = sd(mfi) 
    , beta   = n() / nrow(z)
  ) %>% pivot_longer(-c(param_set, sim_num), values_to = "true")

z.n <- z %>% 
  group_by(param_set, sim_num) %>%
  filter(group == 1) %>%
  summarize(
    mu_neg = mean(mfi)
    , sd_neg = sd(mfi) 
  ) %>% pivot_longer(-c(param_set, sim_num), values_to = "true")

true_vals %<>% left_join(., rbind(z.p, z.n)) %>% ungroup()

} else if (model_fits[i, ]$model == "cluster_regression_with_beta_1.stan") {
  
if (complexity == 1) {
  
samps_out <- with(samps, data.frame(
    mu_base    = mu_base
  , mu_pos     = mu[, 2]
  , sigma_base = sigma_base
  , sigma_pos  = sigma[, 2]
  , beta_base  = plogis(beta_base)
  , beta_cat1f = plogis(beta_base + beta_cat1f_delta)
  , beta_cat1f_delta = beta_cat1f_delta
  ))

true_vals <- y %>% 
  mutate(
    beta_cat1f = plogis(beta_base + beta_cat1f_delta)
  , beta_base  = plogis(beta_base)) %>% 
  dplyr::select(param_set, sim_num, mu_neg, mu_pos, sd_neg, sd_pos, beta_cat1f, beta_base, beta_cat1f_delta) %>% 
  pivot_longer(-c(param_set, sim_num), values_to = "true")
  
} else if (complexity == 2) {
  
samps_out <- with(samps, data.frame(
    mu_base    = mu_base
  , mu_pos     = mu[, 2]
  , sigma_base = sigma_base
  , sigma_pos  = sigma[, 2]
  , beta_base  = beta_base
  , beta_cat1f_delta = beta_cat1f_delta
  ))

true_vals <- y %>% 
  dplyr::select(
    param_set, sim_num, mu_neg, mu_pos, sd_neg, sd_pos
  , beta_base, beta_cat1f_delta
  ) %>% 
  pivot_longer(-c(param_set, sim_num), values_to = "true")  
  
} else {
  stop("Complexity not -yet- supported")
}
  
} else if (model_fits[i, ]$model == "cluster_regression_with_beta_theta_ln_1.stan") {
  
if (complexity == 1) {
  
samps_out <- with(samps, data.frame(
    mu_base    = mu_base
  , mu_pos     = mu[, 2]
  , sigma_base = sigma_base
  , sigma_pos  = sigma[, 2]
  , beta_base  = plogis(beta_base)
  , beta_cat1f = plogis(beta_base + beta_cat1f_delta)
  , beta_cat1f_delta = beta_cat1f_delta
  , theta_cat2f_mu   = theta_cat2f_mu
  ))
  
true_vals <- y %>% 
  mutate(
    beta_cat1f = plogis(beta_base + beta_cat1f_delta)
  , beta_base  = plogis(beta_base)
  ) %>% 
  dplyr::select(param_set, sim_num, mu_neg, mu_pos, sd_neg, sd_pos
              , beta_cat1f, beta_base, beta_cat1f_delta, theta_cat2f_mu
                ) %>% 
  pivot_longer(-c(param_set, sim_num), values_to = "true")
  
} else if (complexity == 2) {
  
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
  
} else {
  stop("Complexity not -yet- supported")
}

  
} else if (model_fits[i, ]$model == "cluster_regression_with_beta_theta_ln_2.stan") { 
 
samps_out <- with(samps, data.frame(
    mu_base      = mu_base
  , mu_pos       = mu_base + mu_diff
  , mu_pos_delta = mu_diff
  , sigma_base   = sigma_base
  , sigma_pos    = sigma[, 2]
  , beta_base    = beta_base
  , beta_cat1f_delta = beta_cat1f_delta
  , theta_cat2f_mu   = theta_cat2f_mu
  , theta_cat1r_sd   = theta_cat1r_sd
  , beta_con1f_delta = beta_con1f_delta
  ))

true_vals <- y %>% 
  dplyr::select(
    param_set, sim_num, mu_neg, mu_pos, mu_pos_delta, sd_neg, sd_pos
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

## Summarize stan fits for pub (reduced function for just the model of relevance for the pub)
summarize_stan_fits_for_pub <- function(model_fits, param_sets, simulated_data, complexity) {
  
  for (i in 1:nrow(model_fits)) {
    
    y <- param_sets %>% filter(
        param_set == model_fits$param_set[i]
      , sim_num   == model_fits$sim_num[i]
    )
    
    samps <- model_fits[i, ]$fitted_model[[1]][[1]]
    
    z <- simulated_data %>% filter(
        param_set == model_fits$param_set[i]
      , sim_num   == model_fits$sim_num[i]
      , log_mfi   == model_fits$log_mfi[i]
    )
    
    samps_out <- with(samps, data.frame(
        mu_base      = mu_base
      , mu_pos       = mu_base + mu_diff
      , mu_pos_delta = mu_diff
      , sigma_base   = sigma_base
      , sigma_pos    = sigma[, 2]
      , beta_base    = beta_base
      , beta_cat1f_delta  = beta_cat1f_delta
      , beta_cat2f_delta  = beta_cat2f_delta
      , beta_con1f_delta  = beta_con1f_delta
      , theta_con2f_delta = theta_con2f_delta
    ))
    
      true_vals <- y %>% 
        dplyr::select(
            param_set, sim_num, mu_neg, mu_pos, mu_pos_delta, sd_neg, sd_pos
          , beta_base, beta_cat1f_delta, beta_cat2f_delta, beta_con1f_delta
          , theta_con2f_delta
        ) %>% 
        pivot_longer(-c(param_set, sim_num), values_to = "true")
      
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
      left_join(model_fits[i, ] %>% dplyr::select(model, param_set, sim_num, log_mfi) %>%
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
        , log_mfi    = model_fits[i, ]$log_mfi
        , true       = (z %>% mutate(group = group - 1, group = ifelse(group > 0, 1, 0)) %>% pull(group) %>% sum()) / y$n_samps
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

## Further cleanup of stan fits (adding summaries and adding sim parameters)
summarize_stan_summary  <- function(stan_summary, param_sets) {
  
lapply(stan_summary, FUN = function(x) {
  x %>% left_join(., param_sets %>% dplyr::select(
    param_set, sim_num, n_samps, beta_base, mu_neg, sd_neg
  , mu_pos, sd_pos, mu_pos_delta, sd_pos_delta
  ))
})
  
}

## Summarize group assignment predictions
calculate_group_assignments <- function(three_sd.g, mclust.g, stan.g, param_sets) {

  three_sd.g %<>%
    dplyr::select(-assigned_group) %>%
    mutate(
        false_pos = ifelse(group == 0 & V2 == 1, 1, 0)
      , false_neg = ifelse(group == 1 & V1 == 1, 1, 0)
    ) %>%
    group_by(model, param_set, sim_num, sd_method, log_mfi, group) %>%
    summarize(
        prob        = mean(V2)
      , false_pos_p = length(which(false_pos == 1)) / n() 
      , false_neg_p = length(which(false_neg == 1)) / n()
    ) %>% mutate(
      misclass_error_p = ifelse(group == 0, false_pos_p, false_neg_p)
    ) %>% dplyr::select(-c(false_pos_p, false_neg_p)) %>% mutate(
      quantile = "mid", .before = prob
    ) %>% left_join(., param_sets %>% dplyr::select(
        param_set, sim_num, n_samps, beta_base, mu_neg, sd_neg
      , mu_pos, sd_pos, mu_pos_delta, sd_pos_delta
    ))
  
  mclust.g %<>%
    dplyr::select(-assigned_group) %>%
    mutate(
      misclass_error = ifelse(group == 0, V2_adj, V1)
    ) %>%
    group_by(model, param_set, sim_num, method, log_mfi, group) %>%
    summarize(
        prob             = mean(V2_adj)
      , misclass_error_p = mean(misclass_error, na.rm = T)
    ) %>% mutate(
      quantile = "mid", .before = prob
    ) %>% left_join(., param_sets %>% dplyr::select(
      param_set, sim_num, n_samps, beta_base, mu_neg, sd_neg
      , mu_pos, sd_pos, mu_pos_delta, sd_pos_delta
    ))
  
  stan.g %<>% 
    rename(model = stan_model) %>%
    dplyr::select(-c(samp, cat1f, cat2f, mfi)) %>%
    pivot_longer(-c(model, param_set, sim_num, log_mfi, group, titer)
                 , names_to = "quantile") %>%
    group_by(model, param_set, sim_num, group, quantile) %>%
    summarize(
      prob = mean(value)
    ) %>% mutate(misclass_error_p = ifelse(group == 0, prob, 1 - prob)) %>%
    left_join(., param_sets %>% dplyr::select(
        param_set, sim_num, n_samps, beta_base, mu_neg, sd_neg
      , mu_pos, sd_pos, mu_pos_delta, sd_pos_delta
    ))
  
  return(
    rbind(
      three_sd.g, mclust.g, stan.g
      )
    )
  
}

## Collate and summarize population level seropositivity estimates
calculate_population_seropositivity <- function(three_sd.g, mclust.g, stan.g, param_sets) {
  
  three_sd.g %<>% 
    group_by(param_set, sim_num, sd_method, log_mfi) %>% 
    summarize(
      true     = mean(group)
    , prop_pos = mean(assigned_group)
    ) %>% mutate(
      prop_pos_diff = prop_pos - true 
    ) %>% ungroup() %>% 
    mutate(model = "3sd", .before = 1) %>%
    mutate(quantile = "mid", .after = "true") %>%
    rename(method = sd_method)
  
  mclust.g %<>% 
    group_by(param_set, sim_num, method, log_mfi) %>% 
    summarize(
      true     = mean(group)
    , prop_pos = sum(V2_adj) / n()
    ) %>% mutate(
      prop_pos_diff = prop_pos - true 
    ) %>% ungroup() %>% 
    mutate(model = "mclust", .before = 1) %>%
    mutate(quantile = "mid", .after = "true") 
  
  stan.g %<>% 
    pivot_longer(-c(stan_model, param_set, sim_num, log_mfi, true)
                 , names_to = "quantile", values_to = "prop_pos") %>%
    mutate(prop_pos_diff = prop_pos - true) %>%
    rename(model = stan_model) %>% 
    mutate(method = "Bayesian LCR", .after = sim_num)
    
  all.g <- rbind(three_sd.g, mclust.g, stan.g) %>% left_join(
    .
  , param_sets %>% dplyr::select(
     param_set, sim_num, n_samps, beta_base, mu_neg, sd_neg
   , mu_pos, sd_pos, mu_pos_delta, sd_pos_delta
  ))
  
  return(all.g)
  
}
