build_stan_priors <- function(simulated_data, skew_fit, fit_attempt) {

  if (!skew_fit) {
    
if (fit_attempt == 1) {
  prior_dat <- data.frame(
      max_mfi             = 30000 %>% log()
    , beta_base_prior_m   = -3
    , beta_base_prior_v   = 2
    , mu_base_prior_m     = 6
    , mu_base_prior_v     = 3
    , mu_diff_prior_m     = 2
    , mu_pos_prior_v      = 3
    , sigma_base_prior_m  = 1
    , sigma_base_prior_v  = 2
    , sigma_diff_prior_m  = 0
    , sigma_diff_prior_v  = 2
    , skew_pos_prior_v    = 3
  ) %>% t() %>% 
    as.data.frame() %>%
    rename(prior = V1) %>%
    mutate(param = rownames(.), .before = 1)
  } else {
  prior_dat <- data.frame(
      max_mfi             = 30000 %>% log()
    , beta_base_prior_m   = -3
    , beta_base_prior_v   = 0.3
    , mu_base_prior_m     = 6
    , mu_base_prior_v     = 0.5
    , mu_diff_prior_m     = 2
    , mu_pos_prior_v     = 0.3
    , sigma_base_prior_m  = 1
    , sigma_base_prior_v  = 2
    , sigma_diff_prior_m  = 0
    , sigma_diff_prior_v  = 2
    , skew_pos_prior_v    = 1
  ) %>% t() %>% 
    as.data.frame() %>%
    rename(prior = V1) %>%
    mutate(param = rownames(.), .before = 1)
  }
    
    p1 <- list(beta_base = rnorm(1, -5, 0.3), mu = c(rnorm(1, 4, 0.3), rnorm(1, 9, 0.3)))
    p2 <- list(beta_base = rnorm(1, -5, 0.3), mu = c(rnorm(1, 4, 0.3), rnorm(1, 9, 0.3)))
    p3 <- list(beta_base = rnorm(1, -5, 0.3), mu = c(rnorm(1, 4, 0.3), rnorm(1, 9, 0.3)))
    p4 <- list(beta_base = rnorm(1, -5, 0.3), mu = c(rnorm(1, 4, 0.3), rnorm(1, 9, 0.3)))
    
  } else {
    
    if (fit_attempt == 1) {
      prior_dat <- data.frame(
          max_mfi             = 30000 %>% log()
        , beta_base_prior_m   = -3
        , beta_base_prior_v   = 2
        , mu_diff_prior_m     = 3
        , mu_diff_prior_v     = 2
        , mu_pos_prior_m      = 9
        , mu_pos_prior_v      = 2
        , sigma_base_prior_m  = 1
        , sigma_base_prior_v  = 2
        , sigma_diff_prior_m  = 0
        , sigma_diff_prior_v  = 2
        , skew_pos_prior_v    = 3
      ) %>% t() %>% 
        as.data.frame() %>%
        rename(prior = V1) %>%
        mutate(param = rownames(.), .before = 1)
    } else {
      prior_dat <- data.frame(
          max_mfi             = 30000 %>% log()
        , beta_base_prior_m   = -3
        , beta_base_prior_v   = 0.3
        , mu_diff_prior_m     = 3
        , mu_diff_prior_v     = 0.5
        , mu_pos_prior_m      = 9
        , mu_pos_prior_v      = 1
        , sigma_base_prior_m  = 1
        , sigma_base_prior_v  = 2
        , sigma_diff_prior_m  = 0
        , sigma_diff_prior_v  = 2
        , skew_pos_prior_v    = 1
      ) %>% t() %>% 
        as.data.frame() %>%
        rename(prior = V1) %>%
        mutate(param = rownames(.), .before = 1)
    }
    
    p1 <- list(beta_base = rnorm(1, -5, 0.3), mu = c(rnorm(1, 4, 0.3), rnorm(1, 9, 0.3))
             , sigma     = c(rnorm(1, 1, 0.3), rnorm(1, 1, 0.3))
             , skew_pos  = rnorm(1, -1, 0.3)
             )
    p2 <- list(beta_base = rnorm(1, -5, 0.3), mu = c(rnorm(1, 4, 0.3), rnorm(1, 9, 0.3))
             , sigma     = c(rnorm(1, 1, 0.3), rnorm(1, 1, 0.3))
             , skew_pos  = rnorm(1, -1, 0.3)
    )
    p3 <- list(beta_base = rnorm(1, -5, 0.3), mu = c(rnorm(1, 4, 0.3), rnorm(1, 9, 0.3))
             , sigma     = c(rnorm(1, 1, 0.3), rnorm(1, 1, 0.3))
             , skew_pos  = rnorm(1, -1, 0.3)
    )
    p4 <- list(beta_base = rnorm(1, -5, 0.3), mu = c(rnorm(1, 4, 0.3), rnorm(1, 9, 0.3))
             , sigma     = c(rnorm(1, 1, 0.3), rnorm(1, 1, 0.3))
             , skew_pos  = rnorm(1, -1, 0.3)
    )
  }
  
  return(
    list(
      priors              = prior_dat
    , starting_conditions = list(p1, p2, p3, p4)
    )
  )
  
}

fit_a_stan_model  <- function(stan_priors, simulated_data, this_model, skew_fit, param_set, max_time) {
  
  if(!skew_fit) {
  
  stan_fit    <- try(R.utils::withTimeout(this_model[[1]]$sample(
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
    , parallel_chains = 1
    , max_treedepth   = 12
    , iter_warmup     = 2000
    , iter_sampling   = 1000
    , adapt_delta     = 0.96
    , seed            = 483892929
    , refresh         = 500
  ), timeout = max_time
  ), silent  = TRUE)
  
  } else {
    
  stan_fit    <- try(R.utils::withTimeout(this_model[[1]]$sample(
    data    = list(
        N           = length(simulated_data$mfi)
      , N_cat1r     = param_set$cat1r_count
      , y           = simulated_data$mfi
      , cat1f       = simulated_data$cat1f
      , cat2f       = simulated_data$cat2f
      , con1f       = simulated_data$con1f
      , beta_base_prior_m   = stan_priors$priors %>% filter(param == "beta_base_prior_m") %>% pull(prior)
      , beta_base_prior_v   = stan_priors$priors %>% filter(param == "beta_base_prior_v") %>% pull(prior)
      , mu_diff_prior_m     = stan_priors$priors %>% filter(param == "mu_diff_prior_m") %>% pull(prior)
      , mu_diff_prior_v     = stan_priors$priors %>% filter(param == "mu_diff_prior_v") %>% pull(prior)
      , mu_pos_prior_m      = stan_priors$priors %>% filter(param == "mu_pos_prior_m") %>% pull(prior)
      , mu_pos_prior_v      = stan_priors$priors %>% filter(param == "mu_pos_prior_v") %>% pull(prior)
      , sigma_base_prior_m  = stan_priors$priors %>% filter(param == "sigma_base_prior_m") %>% pull(prior)
      , sigma_diff_prior_m  = stan_priors$priors %>% filter(param == "sigma_diff_prior_m") %>% pull(prior)
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
    , parallel_chains = 1
    , max_treedepth   = 12
    , iter_warmup     = 2000
    , iter_sampling   = 1000
    , adapt_delta     = 0.96
    , seed            = 483892929
    , refresh         = 500
  ), timeout = max_time
  ), silent  = TRUE)
    
  }
  
  return(stan_fit)
  
}

skew_mean         <- function(loc, scal, shap) {
  a  <- shap / (sqrt(1 + shap^2))
  mm <- loc + scal * a * sqrt(2/base::pi)
  return(mm)
}
