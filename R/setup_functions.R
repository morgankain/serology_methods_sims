## Build parameter sets from input parameter ranges
establish_parameters        <- function(n_param_sets, complexity, ...) {
  
if (complexity == 1) {
needed_params      <- c(
  
  "n_sims_per_set", "n_samps"
  
, "cat1f_prop", "cat2f_prop"

, "beta_base", "beta_cat1f_delta"

, "mu_neg", "sd_neg", "mu_pos_delta", "sd_pos_delta"

, "theta_cat2f_mu"

)
} else if (complexity == 2) {
needed_params      <- c(
  
  "n_sims_per_set", "n_samps"
  
, "cat1f_prop", "cat2f_prop", "cat1r_count", "con1f_sd"

, "beta_base", "beta_cat1f_delta", "beta_con1f_delta"

, "mu_neg", "sd_neg", "mu_pos_delta", "sd_pos_delta"

, "theta_cat2f_mu", "theta_cat1r_sd"

)
} else {
  stop("Model complexity not yet supported")
}

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
  , sd_pos = sd_neg * sd_pos_delta) %>% 
  dplyr::mutate(param_set = seq(n()), .before = 1) %>%
  group_by(param_set)

purrr::map_dfr(seq_len(model_input_params$n_sims_per_set), ~model_params) %>%
  mutate(sim_num = seq(n()), .after = "param_set")
  
}

## Build parameter sets from input parameter ranges, specifically for publication
establish_parameters_for_pub <- function(n_param_sets, complexity, ...) {
 
    needed_params      <- c(
        "n_sims_per_set", "n_samps"
      , "cat1f_prop", "con1f_sd"
      , "beta_base", "beta_cat1f_delta", "beta_con1f_delta"
      , "mu_neg", "sd_neg", "mu_pos_delta", "sd_pos_delta"
      , "theta_con2f_delta"
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
      , sd_pos = sd_neg * sd_pos_delta) %>% 
    dplyr::mutate(param_set = seq(n()), .before = 1) %>%
    group_by(param_set)
  
  purrr::map_dfr(seq_len(model_input_params$n_sims_per_set), ~model_params) %>%
    mutate(sim_num = seq(n()), .after = "param_set")
  
}

## Set up vector of models and check if they can be fit | chosen complexity
establish_models        <- function(model_set, complexity) {
  
  min_complexities <- apply(model_set %>% matrix(), 1, FUN = function(x) {
    tc <- strsplit(x, "[.]")[[1]][1]
    tc <- strsplit(tc, "_")[[1]]
    tc[length(tc)]
  })
  
  if (any(min_complexities  > complexity)) {
    stop("All models must have complexity <= data_complexity (see number directly preceeding `.stan`)")
  } else {
    return(model_set)
  }

}

## Compile the chosen stan models
compile_stan_models     <- function(model_set) {
  
 compiled_stan_models <- model_set %>% as.list() %>% lapply(., FUN = function(x) {
   model_path <- paste0("stan_models/", x)
   stan_model <- cmdstanr::cmdstan_model(model_path, pedantic = FALSE)
   return(stan_model)
 })
 
 model_storage <- tibble(
   model = model_set
 ) %>% mutate(
    index = seq(n())
  , compiled_model = compiled_stan_models
 )
 
 return(model_storage)
    
}

