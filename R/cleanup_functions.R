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
      dplyr::select(param_set, sim_num, beta_base
                    , beta_cat1f_delta, beta_cat2f_delta, beta_con1f_delta) %>% 
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
    
    ## Fitting details
    fit_details <- model_fits$fit_details[[1]] %>% mutate(
        chain = seq(n())
      , .before = 1
    )
    
    z <- simulated_data %>% filter(
        param_set == model_fits$param_set[i]
      , sim_num   == model_fits$sim_num[i]
      , log_mfi   == model_fits$log_mfi[i]
    )
    
    y <- param_sets %>% filter(
        param_set == model_fits$param_set[i]
      , sim_num   == model_fits$sim_num[i]
    )
    
    true_vals <- y %>% 
      dplyr::select(
          param_set, sim_num, mu_neg, mu_pos, mu_pos_delta, sd_neg, sd_pos
        , beta_base, beta_cat1f_delta, beta_cat2f_delta, beta_con1f_delta
      ) %>% 
      pivot_longer(-c(param_set, sim_num), values_to = "true")
    
    
    if (model_fits$fit_success[i] == 1) {
      
    samps <- model_fits[i, ]$fitted_model[[1]]
    
    samps_out <- with(samps, data.frame(
        mu_base      = mu[, 1]
      , mu_pos       = mu[, 2]
      , mu_pos_delta = mu[, 2] - mu[, 1]
      , sigma_base   = sigma[, 1]
      , sigma_pos    = sigma[, 2]
      , beta_base    = beta_base
      , beta_cat1f_delta  = beta_cat1f_delta
      , beta_cat2f_delta  = beta_cat2f_delta
      , beta_con1f_delta  = beta_con1f_delta
    ))
  
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
    
    } else {
      samps_out <- data.frame(name = true_vals$name, lwr = NA, lwr_n = NA, mid = NA, upr_n = NA, upr = NA)
    }
    
    out.clean.t <- left_join(true_vals, samps_out, by = "name") %>% 
      left_join(model_fits[i, ] %>% dplyr::select(model, param_set, sim_num, log_mfi) %>%
                  mutate(param_set = as.numeric(param_set), sim_num = as.numeric(sim_num))
                , ., by = c("param_set", "sim_num")) %>%
      left_join(., y %>% dplyr::select(param_set, sim_num, n_samps), by = c("param_set", "sim_num"))
    
    if (model_fits$fit_success[i] == 1) {
    
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
    
    } else {
      pred_pos <- data.frame(samp = seq(out.clean.t$n_samps[1]), lwr = NA, lwr_n = NA, mid = NA, upr_n = NA, upr = NA)
    }
    
    z %<>% mutate(samp = seq(n()), .before = 1) %>% 
      left_join(., pred_pos, by = "samp") %>% 
      mutate(stan_model = model_fits[i, ]$model, .before = 1)
    
    if (model_fits$fit_success[i] == 1) {
      
    pop_sero <- apply(samps$beta_vec, 1, FUN = function(x) {
      rbinom(length(x), 1, x) %>% sum()
    })
    
    pop_seropos <- (pop_sero / y$n_samps) %>% 
      quantile(c(0.025, 0.200, 0.500, 0.800, 0.975)) %>% 
      t() %>% as.data.frame()
    
    names(pop_seropos) <- c("lwr", "lwr_n", "mid", "upr_n", "upr")
    
    } else {
      pop_seropos <- data.frame(lwr = NA, lwr_n = NA, mid = NA, upr_n = NA, upr = NA)
    }
    
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
      fit_details.f <- fit_details
    } else {
      out.clean     <- rbind(out.clean, out.clean.t)
      z.f           <- rbind(z.f, z)
      pop_seropos.f <- rbind(pop_seropos.f, pop_seropos)
      fit_details.f <- rbind(fit_details.f, fit_details)
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
    , prop_seropos = pop_seropos.f %>% mutate(
      cover = ifelse(true > lwr & true < upr, 1, 0)
        , CI_wid = upr - lwr
        , m_diff = abs(mid - true)
      )
    , fit_details  = fit_details.f
    )
  )
  
}

## Further cleanup of stan fits (adding summaries and adding sim parameters)
summarize_stan_summary  <- function(stan_summary, param_sets) {

  coef         <- lapply(stan_summary, FUN = function(x) {x$coef}) %>% do.call("rbind", .)
  group_pred   <- lapply(stan_summary, FUN = function(x) {x$group_pred}) %>% do.call("rbind", .)
  prop_seropos <- lapply(stan_summary, FUN = function(x) {x$prop_seropos}) %>% do.call("rbind", .)
  fit_details  <- lapply(stan_summary, FUN = function(x) {x$fit_details}) %>% do.call("rbind", .)
  
return(
  list(
    coef         = coef
  , group_pred   = group_pred
  , prop_seropos = prop_seropos
  , fit_details  = fit_details
  )
)
  
}

## Summarize group assignment predictions
calculate_group_assignments <- function(three_sd.g, mclust.g, stan.g, param_sets) {

  three_sd.g %<>%
    dplyr::select(-assigned_group) %>%
    mutate(
        true_pos  = ifelse(group == 0 & V2 == 0, 1, 0)
      , true_neg  = ifelse(group == 1 & V2 == 1, 1, 0)
      , false_pos = ifelse(group == 0 & V2 == 1, 1, 0)
      , false_neg = ifelse(group == 1 & V1 == 1, 1, 0)
    ) %>%
    group_by(model, param_set, sim_num, sd_method, log_mfi, group) %>%
    summarize(
        prob        = mean(V2)
      , false_pos_p = length(which(false_pos == 1)) / n() 
      , false_neg_p = length(which(false_neg == 1)) / n()
      , true_pos_p  = length(which(true_pos == 1)) / n()
      , true_neg_p  = length(which(true_neg == 1)) / n()
    ) %>% mutate(
        misclass_error_p = false_pos_p + false_neg_p
      , correct_class_p  = 1 - misclass_error_p
    ) %>% mutate(
      quantile = "mid", .before = prob
    ) %>% left_join(., param_sets %>% dplyr::select(
        param_set, sim_num, n_samps, beta_base, mu_neg, sd_neg
      , mu_pos, sd_pos, mu_pos_delta, sd_pos_delta
    )) %>% rename(method = sd_method)
  
  mclust.g %<>%
    dplyr::select(-assigned_group) %>%
    mutate(
      misclass_error = ifelse(group == 0, V2_adj, V1)
    ) %>%
    group_by(model, param_set, sim_num, method, log_mfi, group) %>%
    summarize(
        prob             = mean(V2_adj)
      , misclass_error_p = mean(misclass_error, na.rm = T)
      , correct_class_p  = 1 - misclass_error_p
    ) %>% mutate(
        false_pos_p = NA
      , false_neg_p = NA
      , true_pos_p  = NA
      , true_neg_p  = NA
      , .before = "misclass_error_p"
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
    group_by(model, param_set, sim_num, log_mfi, group, quantile) %>%
    summarize(
      prob = mean(value)
    ) %>% mutate(
       misclass_error_p = ifelse(group == 0, prob, 1 - prob)
     , correct_class_p  = 1 - misclass_error_p
    ) %>% mutate(
        false_pos_p = NA
      , false_neg_p = NA
      , true_pos_p  = NA
      , true_neg_p  = NA
      , .before = "misclass_error_p"
    ) %>%
    left_join(., param_sets %>% dplyr::select(
        param_set, sim_num, n_samps, beta_base, mu_neg, sd_neg
      , mu_pos, sd_pos, mu_pos_delta, sd_pos_delta
    )) %>% mutate(
      method = "Bayesian LCR", .after = sim_num
    )
  
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
      , tot_n    = n()
      , tot_pos = sum(V2)
    ) %>% ungroup() %>%
    mutate(index = seq(n())) %>%
    group_by(index) %>%
    mutate(
        lwr = (prop.test(tot_pos, tot_n))$conf.int[1]
      , upr = (prop.test(tot_pos, tot_n))$conf.int[2]
    ) %>%
    rename(
      mid = prop_pos
    ) %>%
    ungroup() %>% 
    mutate(model = "3sd", .before = 1) %>%
    rename(method = sd_method) %>%
    pivot_longer(., c(lwr, mid, upr), names_to = "quantile", values_to = "prop_pos") %>%
    ungroup() %>%
    mutate(prop_pos_diff = prop_pos - true) %>%
    dplyr::select(-c(tot_n, tot_pos, index)) %>%
    relocate(quantile, .after = true)
  
  mclust.g %<>% 
    group_by(param_set, sim_num, method, log_mfi) %>% 
    summarize(
      true     = mean(group)
      , prop_pos = sum(V2_adj) / n()
      , tot_n    = n()
      , tot_pos = sum(assigned_group)
    ) %>% ungroup() %>%
    mutate(index = seq(n())) %>%
    group_by(index) %>%
    mutate(
      lwr = (prop.test(tot_pos, tot_n))$conf.int[1]
      , upr = (prop.test(tot_pos, tot_n))$conf.int[2]
    ) %>%
    rename(
      mid = prop_pos
    ) %>%
    ungroup() %>% 
    mutate(model = "mclust", .before = 1) %>%
    pivot_longer(., c(lwr, mid, upr), names_to = "quantile", values_to = "prop_pos") %>%
    ungroup() %>%
    mutate(prop_pos_diff = prop_pos - true) %>%
    dplyr::select(-c(tot_n, tot_pos, index)) %>%
    relocate(quantile, .after = true)
  
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
