################################################
## Functions for mixture model for MFI values ##
################################################

## Build parameter sets from input parameter ranges
establish_parameters    <- function(n_param_sets, ...) {
  
needed_params      <- c(
  "n_samps", "n_sims_per_set"
, "mu_neg", "sd_neg", "mu_pos_delta", "sd_pos_delta"
, "beta_base", "beta_age_delta"
, "mu_theta_age"
, "age_prop"
)
model_input_params <- list(...) 
  
if (any(needed_params %notin% names(model_input_params))) {
  stop(paste(
    "Need parameter[s]:"
  , paste(needed_params[which(needed_params %notin% names(model_input_params))], collapse = " -- ")
  ))
}
  
point_vals <- which(lapply(model_input_params, length) %>% unlist() == 1)

model_input_params.p <- model_input_params[point_vals]
model_input_params.r <- model_input_params[-point_vals]
  
if (length(model_input_params.r) != 0) {
  if (n_param_sets > 1) {
     param_input_vals <- pomp::sobol_design(
       lower = lapply(model_input_params.r, FUN = function(x) x[1]) %>% unlist()
     , upper = lapply(model_input_params.r, FUN = function(x) x[2]) %>% unlist()
     , nseq  = n_param_sets
    )
  } else {
    param_input_vals <- lapply(model_input_params.r, FUN = function(x) median(x)) %>% as.data.frame()
  }
  
  model_params <- cbind(param_input_vals, data.frame(model_input_params.p))
  
} else {
  model_params <- data.frame(model_input_params.p)
}

model_params %<>% 
  dplyr::mutate(
    mu_pos = mu_neg + mu_pos_delta
  , sd_pos = sd_neg + sd_pos_delta) %>% 
  dplyr::mutate(param_set = seq(n()), .before = 1) %>%
  group_by(param_set)

purrr::map_dfr(seq_len(model_input_params$n_sims_per_set), ~model_params) %>%
  mutate(sim_num = seq(n()), .after = "param_set")
  
}

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

## Determine groupings with 3sd
group_via_3sd           <- function(simulated_data, param_sets) {
  
 simulated_data.l <- simulated_data %>% split_tibble(
   .
 , c("param_set", "sim_num")
 )
 
 three_sd.g <- lapply(simulated_data.l, FUN = function(x) {

   ## Not particularly realistic here, but lets pretend we can use
    ## some of the group 1 mfi values as "negative controls" 
     ## (until I can come up with something better)
   sd_all <- x %>% summarize(sd_all = sd(mfi)) %>%
     pull(sd_all)
    
   mean_neg <- x %>% filter(group == 1) %>% 
     slice(sample(seq(n()), n() / 5)) %>% 
     summarize(mean_neg = mean(mfi)) %>%
     pull(mean_neg)
   
   x %>% mutate(
    assigned_group = ifelse(
      mfi > mean_neg + 3 * sd_all
    , 2
    , 1)
   )
  
 }
 ) %>% do.call("rbind", .)
 
 three_sd.g %>% mutate(
    group = group - 1
   , V1 = ifelse(assigned_group == 1, 1, 0)
   , V2 = 1 - V1
  ) %>% mutate(model = "3sd", .before = 1)
  
}

## Determine groupings with mclust
group_via_mculst        <- function(simulated_data) {
  
 simulated_data.l <- simulated_data %>% split_tibble(
   .
 , c("param_set", "sim_num")
 )
  
mclust.g <- lapply(simulated_data.l, FUN = function(x) {
   ## !! Note: probabilities conditional on group mean and variance (which does has uncertainty)
     ## -- calculate probabilities over mclust posterior (uncertainty in what the mean and
      ## variance are among clusters)
   clust.fit <- Mclust(x$mfi, G = 2, modelNames = "V", verbose = FALSE)
   x %>% cbind(., as.data.frame(clust.fit$z)) %>% 
  mutate(
   assigned_group = clust.fit$classification 
  )
 }) %>% do.call("rbind", .)
 
mclust.g %>% mutate(
    group = group - 1
  ) %>% mutate(model = "mclust", .before = 1)
  
}

## Second phase regression models
fit_regression          <- function(groups) {

groups.l     <- groups %>% split_tibble(., c("param_set"))

regression.pred <- lapply(groups.l, FUN = function(x) {

  groups.sims.l <- x %>% split_tibble(., "sim_num")
  
 regression_sims.pred <- lapply(groups.sims.l, FUN = function(y) {
    
    y %<>% 
      mutate(
        assigned_group = assigned_group - 1
      , age            = as.factor(age)
        )
    
    no_variance <- glm(
      assigned_group ~ age
    , family = "binomial"
    , data   = y
    )

    ## !! NOTE: Uncertain about this -- check. How?
     ## !! Beta?
    with_variance <- glm(
      V2 ~ age
    , family = "binomial"
    , data   = y
    )

    return(
     list(
       no_variance, with_variance
     )
    )
    
  })
  
  return(regression_sims.pred)
    
  })

return(regression.pred)
  
}
sort_regression         <- function(fitted_regressions, param_sets) {

  param_sets.l <- param_sets %>% split_tibble(., "param_set")
  
regression.pred <- purrr::pmap(list(fitted_regressions, param_sets.l), .f = function(x, y) {
  
  param_sets_sims.l <- y %>% split_tibble(., "sim_num")
  
  regression_sims.pred <- purrr::pmap(list(x, param_sets_sims.l), .f = function(v, w) {
    
  true_vals <- w %>%
   dplyr::select(param_set, sim_num, beta_base, beta_age_delta) %>% 
   mutate(beta_age = beta_base + beta_age_delta) %>%
   dplyr::select(-beta_age_delta) %>%
   pivot_longer(-c(param_set, sim_num), values_to = "true")
 
    pred.out <- predictorEffect("age", v[[1]]) %>% summary()
    pred.out <- with(pred.out, data.frame(
      lwr = lower
    , mid = effect
    , upr = upper
    )) %>% mutate(
      name = c("beta_age", "beta_base")
    , .before = 1
    )
    
    out1 <- left_join(true_vals, pred.out, by = "name") %>% mutate(model = "no_variance")
    
    pred.out <- predictorEffect("age", v[[2]]) %>% summary()
    pred.out <- with(pred.out, data.frame(
      lwr = lower
    , mid = effect
    , upr = upper
    )) %>% mutate(
      name = c("beta_age", "beta_base")
    , .before = 1
    )
    
    out2 <- left_join(true_vals, pred.out, by = "name") %>% mutate(model = "positive_probability")
    
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

## Fit all of the different possible stan models
fit_stan_models         <- function(simulated_data, param_sets, model_names) {
  
simulated_data <- simulated_data[[1]]
model_name     <- model_names[[1]]
param_set      <- param_sets %>% filter(
  param_set == unique(simulated_data$param_set)
, sim_num   == unique(simulated_data$sim_num)
)

if (model_name$model == "cluster_regression_base.stan") {
  
stan_fit <- stan(
  file    = "stan_models/cluster_regression_base.stan"
, data    = list(N = param_set$n_samps, y = simulated_data$mfi)
, pars    = c("membership_l", "ind_sero") 
, include = FALSE
, chains  = 4
, seed    = 483892929
, refresh = 2000
, cores   = 1
)

} else if (model_name$model == "cluster_regression_with_beta.stan") {
  
stan_fit <- stan(
  file    = "stan_models/cluster_regression_with_beta.stan"
, data    = list(
   N = param_set$n_samps, y = simulated_data$mfi
 , age = simulated_data$age, age_index = simulated_data$age + 1
 )
, pars    = c("membership_l", "ind_sero") 
, include = FALSE
, chains  = 4
, seed    = 483892929
, refresh = 2000
, cores   = 1
)
  
} else if (model_name$model == "cluster_regression_with_beta_theta.stan") {
  
stan_fit <- stan(
  file    = "stan_models/cluster_regression_with_beta_theta.stan"
, data    = list(
   N = param_set$n_samps, y = simulated_data$mfi
 , age = simulated_data$age, age_index = simulated_data$age + 1
 )
, pars    = c("membership_l", "ind_sero") 
, include = FALSE
, chains  = 4
, seed    = 483892929
, refresh = 2000
, cores   = 1
)
  
} else {
  return(NULL)
}

stan_fit        <- list(stan_fit)
names(stan_fit) <- paste(
  model_name$model
, unique(simulated_data$param_set)
, unique(simulated_data$sim_num)
, sep = " -- ")

return(stan_fit)
  
}

## ^^ Function to deal with what is returned from fit_stan_models_f to get into structure for
 ## how I built summarize_stan_fits.
sort_stan_fits          <- function(stan_fits.l, models_to_fit) {
  ## Need a list of lists
   ## outer list is of length n-models
    ## each of these has length = n_param_sets
 model_names <- lapply(stan_fits.l, FUN = function(x) {
    names(x)
  }) %>% unlist()
 
fit_details <- lapply(model_names %>% as.list(), FUN = function(x) {
   strsplit(x, " -- ") %>% unlist() %>% t() %>% as.data.frame()
  }) %>% do.call("rbind", .) %>% 
  rename(model = V1, param_set = V2, sim_num = V3) %>% as_tibble()

resorted_fits <- fit_details %>% mutate(fitted_model = stan_fits.l)
 
return(resorted_fits)
 
}

## Summarize fits
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

if (model_fits[i, ]$model == "cluster_regression_base.stan") {
  
samps_out <- with(samps, data.frame(
    mu_base    = mu_base
  , mu_pos     = mu[, 2]
  , sigma_base = sigma_base
  , sigma_pos  = sigma[, 2]
  , beta       = beta
  ))

true_vals <- y %>% 
  mutate(beta = 1 - (beta_base*age_prop + (beta_base + beta_age_delta)*age_prop)) %>% 
  dplyr::select(param_set, sim_num, mu_neg, mu_pos, sd_neg, sd_pos, beta) %>% 
  pivot_longer(-c(param_set, sim_num), values_to = "true")

} else if (model_fits[i, ]$model == "cluster_regression_with_beta.stan") {
  
samps_out <- with(samps, data.frame(
    mu_base    = mu_base
  , mu_pos     = mu[, 2]
  , sigma_base = sigma_base
  , sigma_pos  = sigma[, 2]
  , beta_base  = 1 - beta_vec[, 2]
  , beta_age   = 1 - beta_vec[, 1]
  , beta_age_delta = beta_age_delta
  ))

true_vals <- y %>% 
  mutate(beta_age = beta_base + beta_age_delta) %>% 
  dplyr::select(param_set, sim_num, mu_neg, mu_pos, sd_neg, sd_pos, beta_age, beta_base, beta_age_delta) %>% 
  pivot_longer(-c(param_set, sim_num), values_to = "true")
  
} else if (model_fits[i, ]$model == "cluster_regression_with_beta_theta.stan") {
  
samps_out <- with(samps, data.frame(
    mu_base    = mu_base
  , mu_pos     = mu[, 2]
  , sigma_base = sigma_base
  , sigma_pos  = sigma[, 2]
  , beta_base  = 1 - beta_vec[, 2]
  , beta_age   = 1 - beta_vec[, 1]
  , beta_age_delta = beta_age_delta
  , mu_theta_age   = mu_theta_age
  ))

true_vals <- y %>% 
  mutate(beta_age = beta_base + beta_age_delta) %>% 
  dplyr::select(param_set, sim_num, mu_neg, mu_pos, sd_neg, sd_pos, beta_age, beta_base, beta_age_delta, mu_theta_age) %>% 
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
              mutate(param_set = as.numeric(param_set), sim_num = as.numeric(sim_num)), ., by = c("param_set", "sim_num")) %>%
  left_join(., y %>% dplyr::select(param_set, sim_num, n_samps), by = c("param_set", "sim_num"))

pred_pos <- samps$membership_p %>% 
  reshape2::melt() %>% 
  rename(clust = Var2, samp = Var3) %>%
  pivot_wider(values_from = value, names_from = clust) %>%
  group_by(samp) %>%
  summarize(
    lwr   = quantile(`2`, 0.025)
  , lwr_n = quantile(`2`, 0.200)
  , mid   = quantile(`2`, 0.500)
  , upr_n = quantile(`2`, 0.800)
  , upr   = quantile(`2`, 0.975)
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
  
  stan.g %<>% rename(model = stan_model) %>%
    dplyr::select(-c(samp, age, mfi)) %>%
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
    mutate(
      group = group - 1
    , assigned_group = assigned_group - 1
    ) %>%
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
    mutate(
      group = group - 1
    , assigned_group = assigned_group - 1
    ) %>%
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

## collate_outputs
collate_outputs         <- function(
    pop_seropositivity, group_assignment
  , three_sd.sum, mclust.sum, stan.sum) {
  
  all.out <- rbind(
    stan.sum %>% dplyr::select(-n_samps)
  , mclust.sum %>% relocate(model, .before = 1) %>%
      mutate(lwr_n = NA, .after = lwr) %>% mutate(upr_n = NA, .after = mid)
  , three_sd.sum %>% relocate(model, .before = 1) %>%
      mutate(lwr_n = NA, .after = lwr) %>% mutate(upr_n = NA, .after = mid)
  )
  
  coverage <- all.out %>% filter(name %in% c("beta_base", "beta_age")) %>%
    group_by(param_set, model, name) %>%
    summarize(
      coverage = mean(cover)
    )
  
  return(
    list(
      pop_seropositivity = pop_seropositivity
    , group_assignment   = group_assignment
    , coefficient_ests   = all.out
    , coverage           = coverage
    )
  )
  
}

## Explore fits for the regression coefficients 
plot_summary            <- function(coef_ests, param_sets, coverage) {
  
  stan_all.gg <- coef_ests %>% filter(!grepl("3sd|mclust", model)) %>% {
    ggplot(., aes(mid, name)) + 
      geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0, linewidth = 1) +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, linewidth = 0.3) +
      geom_point(aes(true, name), colour = "firebrick3") +
      facet_wrap(~model)
  }
  
  all_out.gg <- coef_ests %>% filter(name %in% c("beta_base", "beta_age")) %>% {
    ggplot(., aes(mid, model)) + 
      geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0, linewidth = 1) +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, linewidth = 0.3) +
      geom_point(aes(true, model), colour = "firebrick3") +
      facet_wrap(~name)
  }
  
  all_out.theta <- coef_ests %>% left_join(., param_sets, by = c("param_set", "sim_num")) %>%
   filter(name %in% c("beta_base", "beta_age"))
  
  cov1.gg <- all_out.theta %>% {
    ggplot(., aes(mu_pos_delta, cover)) + 
      geom_jitter(height = 0.05) +
      facet_wrap(~model)
  }
  
  cov2.gg <- all_out.theta %>% 
  mutate(mu_pos_delta_r = plyr::round_any(mu_pos_delta, 0.5)) %>%
  group_by(model, name, mu_pos_delta_r) %>% 
  summarize(m_cover = mean(cover)) %>% {
    ggplot(., aes(mu_pos_delta_r, m_cover)) + 
      geom_point(size = 2) +
      geom_line() +
      facet_grid(name~model)
  }
  
  cov3.gg <- all_out.theta %>% {
    ggplot(., aes(mu_pos_delta, CI_wid)) + 
      geom_point() +
      facet_wrap(~model)
  }
  
  cov4.gg <- all_out.theta %>% {
    ggplot(., aes(mu_pos_delta, m_diff)) + 
      geom_point() +
      facet_wrap(~model)
  }
  
  cov5.gg <- coef_ests %>% filter(name %in% c("beta_base", "beta_age")) %>%
   left_join(., param_sets, by = c("param_set", "sim_num")) %>% 
   arrange(desc(mu_pos_delta)) %>%
   mutate(mu_pos_delta = factor(mu_pos_delta, levels = unique(mu_pos_delta))) %>% {
     ggplot(., aes(mid, mu_pos_delta)) + 
       geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0, linewidth = 1) +
       geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, linewidth = 0.3) +
       geom_vline(aes(xintercept = true), colour = "firebrick3", linewidth = 0.3) +
       facet_grid(model~name) +
       xlab("Estimate") +
       ylab("Difference in -mean- between positive/negative") +
       theme(axis.text.y = element_text(size = 9))
   }
  
  cov6.gg <- coverage %>%
    ungroup() %>%
    group_by(model, name) %>% 
    summarize(coverage = mean(coverage)) %>% {
    ggplot(., aes(coverage, model)) +
    geom_point() +
    facet_wrap(~name) 
  }

  return(
    list(
      stan_all.gg
    , all_out.gg
    , cov1.gg
    , cov2.gg
    , cov3.gg
    , cov4.gg
    , cov5.gg
    , cov6.gg
    )
  )
  
}

## Explore individual-level group assignments
plot_group_assignments  <- function(group_assignment) {
  
group_assignment %>% 
    ungroup() %>%
    group_by(model, group, quantile) %>%
    summarize(
      lwr   = quantile(prob, 0.025)
    , lwr_n = quantile(prob, 0.200)
    , mid   = quantile(prob, 0.500)
    , upr_n = quantile(prob, 0.800)
    , upr   = quantile(prob, 0.975)
    ) %>% filter(quantile == "mid") %>% 
    mutate(group = as.factor(group)) %>% {
      ggplot(., aes(mid, group)) +
        geom_errorbar(aes(xmin = lwr_n, xmax = upr_n), linewidth = 1, width = 0) +
        geom_errorbar(aes(xmin = lwr, xmax = upr), linewidth = 0.5, width = 0.2) +
        geom_point() +
        facet_wrap(~model, ncol = 1) +
        xlab("Probability of Group 1 Affiliation") +
        ylab("True Group Affiliation")
    }
  
}

## Plot population level seropositivity
plot_pop_seropos        <- function(pop_seropositivity) {
  
    pop_seropositivity %>% dplyr::select(-prop_pos_diff) %>%
    pivot_wider(c(model, param_set, sim_num, true)
                , names_from = "quantile"
                , values_from = "prop_pos") %>% 
    mutate(
      run = interaction(param_set, sim_num)
    , sim_num = as.factor(sim_num)
    ) %>% {
    ggplot(., aes(mid, sim_num)) + 
        geom_errorbar(aes(xmin = lwr_n, xmax = upr_n, colour = sim_num), linewidth = 1, width = 0) +
        geom_errorbar(aes(xmin = lwr, xmax = upr, colour = sim_num), linewidth = 0.5, width = 0.2) +
        geom_point(aes(colour = sim_num)) +
        geom_vline(aes(xintercept = true, colour = sim_num)) +
        scale_colour_brewer(palette = "Dark2") +
        facet_grid(model~param_set) +
        theme(axis.text.x = element_text(size = 10)) +
        xlab("Estimate") + ylab("Simulation")
    }
  
}
