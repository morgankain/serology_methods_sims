################################################
## Functions for mixture model for MFI values ##
################################################

## Build parameter sets from input parameter ranges
establish_parameters    <- function(n_param_sets, ...) {
  
needed_params      <- c(
  "n_samps", "lambda", "mu_neg", "sd_neg", "mu_pos_delta", "sd_pos_delta"
, "beta_age", "age_prop"
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
  
  model_input_params <- cbind(param_input_vals, data.frame(model_input_params.p))
  
} else {
  model_input_params <- data.frame(model_input_params.p)
}

model_input_params %>% 
  dplyr::mutate(
    mu_pos = mu_neg + mu_pos_delta
  , sd_pos = sd_neg + sd_pos_delta) %>% 
  dplyr::mutate(param_set = seq(n()), .before = 1)

}

## Simulate data for all parameter sets
simulate_data           <- function(param_sets) {
  
param_sets %<>% split_tibble(., "param_set")
  
lapply(param_sets, FUN = function(x) {
  mu_vec     <- with(x, c(mu_neg, mu_pos))
  sd_vec     <- with(x, c(sd_neg, sd_pos))
data.frame(
    age   = with(x, rbinom(n_samps, 1, age_prop))
  ) %>% mutate(
    group = rbinom(n(), 1, x$lambda + x$lambda_age * abs(age - 1)) + 1
  , mfi   = rnorm(n(), mu_vec[group] + (age * x$beta_age * (group - 1)), sd_vec[group])
  ) %>% mutate(param_set = x$param_set, .before = 1)
}) %>% do.call("rbind", .)
  
}

## Quick and dirty plot of data -- could eventually become more complicated
examine_data            <- function(simulated_data) {
  
if (n_distinct(simulated_data$param_set) <= 30) {
  
 data_plot <- simulated_data %>% mutate(
  group = as.factor(group)
, int   = interaction(age, group)) %>% {
  ggplot(., aes(x = mfi)) + 
    geom_histogram(aes(colour = int, fill = int)) + 
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    facet_wrap(~param_set)
}
} else {
  data_plot <- c("Too many prameters to plot data. Will only plot if n_param_sets <= 30")
}
  
  return(data_plot)
  
}

## Determine groupings with mclust
group_via_mculst        <- function(simulated_data, param_sets) {
  
 simulated_data.l <- simulated_data %>% split_tibble(., "param_set")
  
 lapply(simulated_data.l, FUN = function(x) {
   clust.fit <- Mclust(x$mfi, G = 2, modelNames = "V", verbose = FALSE)
   x %>% cbind(., as.data.frame(clust.fit$z)) %>% 
  mutate(
   assigned_group = clust.fit$classification 
  )
 }) %>% do.call("rbind", .)
  
}

## Second phase regression models
fit_regression          <- function(groups, param_sets) {

groups.l     <- groups %>% split_tibble(., "param_set")

regression.pred <- lapply(groups.l, FUN = function(x) {

    x %<>% 
      mutate(
        assigned_group = assigned_group - 1
      , age = as.factor(age)
        )
    
    no_variance <- glm(
      assigned_group ~ age
    , family = "binomial"
    , data = x
    )

    with_variance <- glm(
      V2 ~ age
    , family = "binomial"
    , data = x
    )

    return(
     list(
       no_variance, with_variance
     )
    )
    
  })

return(regression.pred)
  
}
sort_regression         <- function(fitted_regressions, param_sets) {

  param_sets.l <- param_sets %>% split_tibble(., "param_set")
  
regression.pred <- purrr::pmap(list(fitted_regressions, param_sets.l), .f = function(x, y) {
    
  true_vals <- y %>%
   dplyr::select(param_set, lambda, lambda_age) %>% 
   rename(theta_age2 = lambda) %>%
   mutate(theta_age1 = theta_age2 + lambda_age) %>%
   dplyr::select(-lambda_age) %>%
   pivot_longer(-param_set, values_to = "true")
 
    pred.out <- predictorEffect("age", x[[1]]) %>% summary()
    pred.out <- with(pred.out, data.frame(
      lwr = lower
    , mid = effect
    , upr = upper
    )) %>% mutate(
      name = c("theta_age1", "theta_age2")
    , .before = 1
    )
    
    out1 <- left_join(true_vals, pred.out, by = "name") %>% mutate(model = "no_variance")
    
    pred.out <- predictorEffect("age", x[[2]]) %>% summary()
    pred.out <- with(pred.out, data.frame(
      lwr = lower
    , mid = effect
    , upr = upper
    )) %>% mutate(
      name = c("theta_age1", "theta_age2")
    , .before = 1
    )
    
    out2 <- left_join(true_vals, pred.out, by = "name") %>% mutate(model = "positive_probability")
    
    return(
      rbind(
        out1, out2
      )
    )
    
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
param_set      <- param_sets %>% filter(param_set == unique(simulated_data$param_set))

if (model_name$model == "mixing_simple_diff.stan") {
  
stan_fit <- stan(
  file    = "stan_models/mixing_simple_diff.stan"
, data    = list(N = param_set$n_samps, y = simulated_data$mfi)
, chains  = 4
, seed    = 483892929
, refresh = 2000
, cores   = 1
)

} else if (model_name$model == "mixing_simple_diff_beta.stan") {
  
stan_fit <- stan(
  file    = "stan_models/mixing_simple_diff_beta.stan"
, data    = list(N = param_set$n_samps, y = simulated_data$mfi, age = simulated_data$age)
, chains  = 4
, seed    = 483892929
, refresh = 2000
, cores   = 1
)
  
} else if (model_name$model == "mixing_simple_diff_beta2.stan") {
  
stan_fit <- stan(
  file    = "stan_models/mixing_simple_diff_beta2.stan"
, data    = list(N = param_set$n_samps, y = simulated_data$mfi, age = simulated_data$age, age_index = simulated_data$age + 1)
, chains  = 4
, seed    = 483892929
, refresh = 2000
, cores   = 1
)
  
} else {
  return(NULL)
}

stan_fit        <- list(stan_fit)
names(stan_fit) <- paste(model_name$model, unique(simulated_data$param_set), sep = " -- ")

return(stan_fit)
  
}

## Summarize fits
summarize_stan_fits     <- function(model_fits, param_sets, stan_models) {
  
param_sets.l     <- split_tibble(param_sets, "param_set") 

for (i in seq_along(model_fits)) {
  
out.clean.t <- purrr::pmap(list(model_fits[[i]], param_sets.l), .f = function(x, y) {

samps     <- rstan::extract(x[[1]])

if (stan_models[i] == "mixing_simple_diff.stan") {
  
samps_out <- with(samps, data.frame(
    mu_base    = mu[, 1]
  , mu_pos     = mu[, 2]
  , sigma_base = sigma[, 1]
  , sigma_pos  = sigma[, 2]
  , theta      = theta
  ))

true_vals <- y %>% mutate(theta = 1 - (lambda*age_prop + (lambda + lambda_age)*age_prop)) %>% 
  dplyr::select(param_set, mu_neg, mu_pos, sd_neg, sd_pos, theta) %>% 
  pivot_longer(-param_set, values_to = "true")
  
} else if (stan_models[i] == "mixing_simple_diff_beta.stan") {
  
samps_out <- with(samps, data.frame(
    mu_base    = mu[, 1]
  , mu_pos     = mu[, 2]
  , sigma_base = sigma[, 1]
  , sigma_pos  = sigma[, 2]
  , theta      = theta
  , beta_age   = beta_age
  ))

true_vals <- y %>% mutate(theta = 1 - (lambda*age_prop + (lambda + lambda_age)*age_prop)) %>% 
  dplyr::select(param_set, mu_neg, mu_pos, sd_neg, sd_pos, theta, beta_age) %>% 
  pivot_longer(-param_set, values_to = "true")
  
} else if (stan_models[i] == "mixing_simple_diff_beta2.stan") {
  
samps_out <- with(samps, data.frame(
    mu_base    = mu[, 1]
  , mu_pos     = mu[, 2]
  , sigma_base = sigma[, 1]
  , sigma_pos  = sigma[, 2]
  , theta_age1 = 1 - theta[, 1]
  , theta_age2 = 1 - theta[, 2]
  , theta_age  = theta_age
  , beta_age   = beta_age
  ))

true_vals <- y %>%
  dplyr::select(param_set, mu_neg, mu_pos, sd_neg, sd_pos, beta_age, lambda, lambda_age) %>% 
  rename(theta_age2 = lambda) %>%
  mutate(theta_age1 = theta_age2 + lambda_age) %>%
  rename(theta_age = lambda_age) %>%
#  mutate(
#    theta_age1 = 1 - theta_age1
#  , theta_age2 = 1 - theta_age2
#  ) %>%
  pivot_longer(-param_set, values_to = "true")
  
} else {
  stop("Stan model name not known")
}

samps_out %<>% mutate(samp = seq(n()), .before = 1) %>%
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

left_join(true_vals, samps_out, by = "name")

}) 
out.clean.t %<>% do.call("rbind", .) %>% 
  left_join(param_sets %>% dplyr::select(param_set, n_samps), ., by = "param_set") %>% 
  mutate(model_name = stan_models[i], .after = param_set)

if (i == 1) {
  out.clean <- out.clean.t
} else {
  out.clean <- rbind(out.clean, out.clean.t)
}

}

return(out.clean %>% mutate(
  cover = ifelse(true > lwr & true < upr, 1, 0)
, CI_wid = upr - lwr
, m_diff = abs(mid - true)
  ))
  
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
 
 list_entries <- lapply(models_to_fit %>% as.list(), FUN = function(x) {
   which(grepl(x, model_names))
 })
 
 resorted_fits <- lapply(list_entries, FUN = function(x) {
   stan_fits.l[x]
 })
 
 return(resorted_fits)
 
}

## Explore coverage of the fits for the beta 
plot_summary            <- function(stan.sum, mclust.sum, param_sets) {
  
  stan_all.gg <- stan.sum %>% {
    ggplot(., aes(mid, name)) + 
      geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0, linewidth = 1) +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, linewidth = 0.3) +
      geom_point(aes(true, name), colour = "firebrick3") +
      facet_wrap(~model_name)
  }
  
  all.out <- rbind(
    stan.sum %>% dplyr::select(-n_samps)
  , mclust.sum %>% rename(model_name = model) %>% relocate(model_name, .after = "param_set") %>%
      mutate(lwr_n = NA, .after = lwr) %>% mutate(upr_n = NA, .after = mid)
  )
  
  all_out.gg <- all.out %>% filter(name %in% c("theta_age1", "theta_age2")) %>% {
    ggplot(., aes(mid, model_name)) + 
      geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0, linewidth = 1) +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, linewidth = 0.3) +
      geom_point(aes(true, model_name), colour = "firebrick3") +
      facet_wrap(~name)
  }
  
  all_out.theta <- all.out %>% left_join(., param_sets, by = "param_set") %>%
  filter(name %in% c("theta_age1", "theta_age2"))
  
  cov1.gg <- all_out.theta %>% {
    ggplot(., aes(mu_pos_delta, cover)) + 
      geom_jitter(height = 0.05) +
      facet_wrap(~model_name)
  }
  
  cov2.gg <- all_out.theta %>% 
  mutate(mu_pos_delta_r = plyr::round_any(mu_pos_delta, 0.5)) %>%
  group_by(model_name, name, mu_pos_delta_r) %>% 
  summarize(m_cover = mean(cover)) %>% {
    ggplot(., aes(mu_pos_delta_r, m_cover)) + 
      geom_point(size = 2) +
      geom_line() +
      facet_grid(name~model_name)
  }
  
  cov3.gg <- all_out.theta %>% {
    ggplot(., aes(mu_pos_delta, CI_wid)) + 
      geom_point() +
      facet_wrap(~model_name)
  }
  
  cov4.gg <- all_out.theta %>% {
    ggplot(., aes(mu_pos_delta, m_diff)) + 
      geom_point() +
      facet_wrap(~model_name)
  }
  
  cov5.gg <- all.out %>% filter(name %in% c("theta_age1", "theta_age2")) %>%
  left_join(., param_sets, by = "param_set") %>% 
  arrange(desc(mu_pos_delta)) %>%
  mutate(mu_pos_delta = factor(mu_pos_delta, levels = unique(mu_pos_delta))) %>% {
    ggplot(., aes(mid, mu_pos_delta)) + 
      geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0, linewidth = 1) +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, linewidth = 0.3) +
      geom_vline(aes(xintercept = true), colour = "firebrick3", linewidth = 0.3) +
      facet_grid(model_name~name) +
      xlab("Estimate") +
      ylab("Difference in -mean- between positive/negative") +
      theme(axis.text.y = element_text(size = 9))
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
    )
  )
  
}

