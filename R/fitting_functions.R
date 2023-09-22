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
   , assigned_group = assigned_group - 1
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
 
mclust.g %>% 
  mutate(
    group          = group - 1
  , assigned_group = assigned_group - 1
  ) %>% mutate(
    model = "mclust", .before = 1
  )
  
}

## Second phase regression models
fit_regression          <- function(groups, complexity) {

groups.l     <- groups %>% split_tibble(., c("param_set"))

regression.pred <- lapply(groups.l, FUN = function(x) {

  groups.sims.l <- x %>% split_tibble(., "sim_num")
  
 regression_sims.pred <- lapply(groups.sims.l, FUN = function(y) {
    
   if (complexity == 1) {
     
    y %<>% mutate(cat1f = as.factor(cat1f))
    
    no_variance <- glm(
      assigned_group ~ cat1f
    , family = "binomial"
    , data   = y
    )

    ## !! NOTE: Uncertain about this -- check. How?
     ## !! Beta?
    with_variance <- glm(
      V2 ~ cat1f
    , family = "binomial"
    , data   = y
    )
     
   } else if (complexity == 2) {
     
    y %<>% mutate(cat1f = as.factor(cat1f))
    
    no_variance <- glm(
      assigned_group ~ cat1f + con1f
    , family = "binomial"
    , data   = y
    )

    ## !! NOTE: Uncertain about this -- check. How?
     ## !! Beta?
    with_variance <- glm(
      V2 ~ cat1f + con1f
    , family = "binomial"
    , data   = y
    )
     
   } else {
     stop("Complexity not -yet- supported")
   }
   
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

## Fit all of the different possible stan models
fit_stan_models         <- function(simulated_data, param_sets, model_names) {
  
simulated_data <- simulated_data[[1]]
model_name     <- model_names[[1]]
param_set      <- param_sets %>% filter(
  param_set == unique(simulated_data$param_set)
, sim_num   == unique(simulated_data$sim_num)
)

if (model_name$model == "cluster_regression_base_1.stan") {
  
stan_fit <- stan(
  file    = "stan_models/cluster_regression_base_1.stan"
, data    = list(
  N = param_set$n_samps
, y = simulated_data$mfi
)
, pars    = c("membership_l", "ind_sero") 
, include = FALSE
, chains  = 4
, seed    = 483892929
, refresh = 2000
, cores   = 1
)

} else if (model_name$model == "cluster_regression_with_beta_1.stan") {
  
stan_fit <- stan(
  file    = "stan_models/cluster_regression_with_beta_1.stan"
, data    = list(
   N           = param_set$n_samps
 , y           = simulated_data$mfi
 , cat1f       = simulated_data$cat1f
 )
, pars    = c("membership_l", "ind_sero", "log_beta", "beta_vec") 
, include = FALSE
, chains  = 4
, seed    = 483892929
, refresh = 2000
, cores   = 1
)
  
} else if (model_name$model == "cluster_regression_with_beta_theta_ln_1.stan") {
  
stan_fit <- stan(
  file    = "stan_models/cluster_regression_with_beta_theta_ln_1.stan"
, data    = list(
   N           = param_set$n_samps
 , y           = simulated_data$mfi
 , cat1f       = simulated_data$cat1f
 , cat2f       = simulated_data$cat2f
 )
, pars    = c("membership_l", "ind_sero", "log_beta", "beta_vec") 
, include = FALSE
, chains  = 4
, seed    = 483892929
, refresh = 2000
, cores   = 1
)
  
} else if (model_name$model == "cluster_regression_with_beta_theta_ln_2.stan") {
  
stan_fit <- stan(
  file    = "stan_models/cluster_regression_with_beta_theta_ln_2.stan"
, data    = list(
   N           = param_set$n_samps
 , N_cat1r     = param_set$cat1r_count
 , y           = simulated_data$mfi
 , cat1f       = simulated_data$cat1f
 , cat2f       = simulated_data$cat2f
 , con1f       = simulated_data$con1f
 , cat1r       = simulated_data$cat1r
 )
, pars    = c("membership_l", "ind_sero", "log_beta", "beta_vec", "theta_cat1r_eps") 
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
