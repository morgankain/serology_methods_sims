## Determine groupings with 3sd
group_via_3sd           <- function(simulated_data, param_sets, groupings) {
  
 simulated_data.l <- simulated_data %>% split_tibble(., groupings)
 
 three_sd.g <- lapply(simulated_data.l, FUN = function(x) {

   ## Lets first assume we use a robust mean and sd estimator to downweight extreme values
   neg_set_a <- x %>% 
     summarize(
       mean_neg = dplR::tbrm(mfi)
     , sd_neg   = jointseg::estimateSd(mfi)
     ) %>%
     dplyr::select(mean_neg, sd_neg) %>% c() %>% unlist()
   
   ## As an alternative, lets try and model a true control (i.e., blanks or some lab colony or something like that),
    ## which we expect to have lower MFI than many wild animals. So instead, here I sort the MFI values and then sample
     ## 20 from the lowest 50% of MFI values
   ## NOTE: this will obviously have a huge impact on the standard deviation, which could result in a quite different
    ## cutoff and thus quite different designations for what is a seropositive and seronegative individual
   ## BUT: to a large degree this is the point of the whole paper. Either of these seem defensible enough and have
    ## large implications
   neg_set_b <- x %>% 
     arrange(mfi) %>%
     slice(1:(max(c(round((1/2)*n(), 1), 20)))) %>% 
     slice(sample(seq(n()), min(20, n()))) %>%
     summarize(
       mean_neg = mean(mfi)
       , sd_neg   = sd(mfi)
     ) %>%
     dplyr::select(mean_neg, sd_neg) %>% c() %>% unlist()
   
   x %>% mutate(
     assigned_group_robust = ifelse(
       mfi > neg_set_a["mean_neg"] + 3 * neg_set_a["sd_neg"]
       , 2
       , 1
     )
     , assigned_group_control = ifelse(
       mfi > neg_set_b["mean_neg"] + 3 * neg_set_b["sd_neg"]
       , 2
       , 1
     )
   ) %>% pivot_longer(., c(assigned_group_robust, assigned_group_control)
                     , names_to = "sd_method", values_to = "assigned_group")
  
 }
 ) %>% do.call("rbind", .)
 
 three_sd.g %>% mutate(
     group  = group - 1
   , V1     = ifelse(assigned_group == 1, 1, 0)
   , V2     = 1 - V1
   , assigned_group  = assigned_group - 1
  ) %>% 
   mutate(model = "3sd", .before = 1)
  
}

## Explore grouping by 3sd across alternative strategies for determining sd
group_via_3sd_alt       <- function(simulated_data, param_sets, groupings) {
  
  perc_sd <- c(seq(0, 1, by = 0.05), 2)
  
  simulated_data.l <- simulated_data %>% split_tibble(., groupings)
  
  three_sd.g <- lapply(simulated_data.l, FUN = function(x) {
    
   for (z in seq_along(perc_sd)) {
    
    ## As an alternative, lets try and model a true control (i.e., blanks or some lab colony or something like that),
    ## which we expect to have lower MFI than many wild animals. So instead, here I sort the MFI values and then sample
    ## 20 from the lowest 10% of MFI values
    ## NOTE: this will obviously have a huge impact on the standard deviation, which could result in a quite different
    ## cutoff and thus quite different designations for what is a seropositive and seronegative individual
    ## BUT: to a large degree this is the point of the whole paper. Either of these seem defensible enough and have
    ## large implications
     
     if (perc_sd[z] <= 1) {
    neg_set <- x %>% 
      arrange(mfi) %>%
      slice(1:(max(c(round(perc_sd[z]*n(), 1), 20)))) %>% 
      slice(sample(seq(n()), min(20, n()))) %>%
      summarize(
          mean_neg = mean(mfi)
        , sd_neg   = sd(mfi)
      ) %>%
      dplyr::select(mean_neg, sd_neg) %>% c() %>% unlist()
     } else {
    neg_set <- x %>%
      summarize(
          mean_neg = dplR::tbrm(mfi)
        , sd_neg   = jointseg::estimateSd(mfi)
      ) %>%
      dplyr::select(mean_neg, sd_neg) %>% c() %>% unlist()
     }
    
    x.t <- x %>% mutate(
      assigned_group_sample_mfi = ifelse(
        mfi > neg_set["mean_neg"] + 3 * neg_set["sd_neg"]
        , 2
        , 1
      )
    ) %>% pivot_longer(., assigned_group_sample_mfi
                       , names_to = "sd_method"
                       , values_to = "assigned_group") %>%
      mutate(perc_sd = perc_sd[z])
    
    if (z == 1) {
      x.f <- x.t
    } else {
      x.f <- rbind(x.f, x.t)
    }
    
   }
    
    x.f
    
  }
  ) %>% do.call("rbind", .)
  
  three_sd.g %>% mutate(
    group  = group - 1
    , V1     = ifelse(assigned_group == 1, 1, 0)
    , V2     = 1 - V1
    , assigned_group  = assigned_group - 1
  ) %>% 
    mutate(model = "3sd", .before = 1)
  
}

## Determine groupings with mclust
group_via_mculst        <- function(simulated_data, groupings) {
  
 simulated_data.l <- simulated_data %>% split_tibble(., groupings)
  
mclust.g <- lapply(simulated_data.l, FUN = function(x) {

  clust.fit_constrained   <- Mclust(x$mfi, G = 2, modelNames = "V", verbose = FALSE)
  clust.fit_unconstrained <- Mclust(x$mfi, modelNames = "V", verbose = FALSE)
  
  dif_groups <- clust.fit_unconstrained$G - clust.fit_constrained$G
  
  unconstrained_z <- clust.fit_unconstrained$z
  constrained_z   <- clust.fit_constrained$z
  
  unconstrained_class <- clust.fit_unconstrained$classification
  constrained_class   <- clust.fit_constrained$classification
  
  if (dif_groups > 0) {
    constrained_z <- cbind(constrained_z, matrix(data = NA, nrow = nrow(constrained_z), ncol = dif_groups))
  }
  if (dif_groups < 0) {
    unconstrained_z <- cbind(unconstrained_z, matrix(data = NA, nrow = nrow(unconstrained_z), ncol = abs(dif_groups)))
  }
  
  constrained_x   <- x %>% mutate(method = "constrained_mclust") %>% 
    cbind(., as.data.frame(constrained_z)) %>% 
    mutate(assigned_group = constrained_class)
  
  unconstrained_x <- x %>% mutate(method = "unconstrained_mclust") %>% 
    cbind(., as.data.frame(unconstrained_z)) %>% 
    mutate(assigned_group = unconstrained_class)
  
  rbind(constrained_x, unconstrained_x)
   
 })

max_v <- lapply(mclust.g, FUN = function(x) {
  x %>% dplyr::select(contains("V")) %>% ncol()
}) %>% unlist() %>% max()

mclust.g.c <- lapply(mclust.g, FUN = function(x) {
  t_num       <- x %>% dplyr::select(contains("V")) %>% ncol()
  if (t_num < max_v) {
    needed_cols <- paste("V", seq((t_num + 1), (t_num + (max_v - t_num))), sep = "")
    highest_v   <- paste("V", t_num, sep = "")
    for (i in seq_along(needed_cols)){
      x %<>% mutate(
        !!needed_cols[i] := NA
        , .after = highest_v
      )
      highest_v <- needed_cols[i]
    } 
  }
  x
}) %>% do.call("rbind", .) 

## Want to retain all of the raw probabilities for all the unconstrained groups to make the point about
 ## an unconstrained cluster model, but fundamentally the goal here is to compare non-exposure to exposure
  ## (we could think about this as a coarse / hacky way of getting "no exposure" vs "varying levels of exposure")
   ## so still want to collapse YES vs NO for the unconstrained cluster model for the purpose of comparing to the
    ## rest of the methods. While this could conceivably be done in a number of ways, we don't want the permutations
     ## to become --too-- large here, so here we simply compare the first to all other clusters
      ## (but again, retaining the raw cluster probabilities to make a different point later)
 
mclust.g.c %>% 
  mutate(
    V2_adj         = 1 - V1
  , group          = group - 1
  , assigned_group = assigned_group - 1
  , assigned_group = ifelse(assigned_group > 0, 1, 0)
  ) %>% mutate(
    model = "mclust", .before = 1
  )
  
}

## Second phase regression models
fit_regression          <- function(groupings, gam_formula, complexity, groupings1, groupings2, method) {

groups.l     <- groupings %>% split_tibble(., groupings1)

regression.pred <- lapply(groups.l, FUN = function(x) {

  groups.sims.l <- x %>% split_tibble(., groupings2)
  
 regression_sims.pred <- lapply(groups.sims.l, FUN = function(y) {
    
     y %<>% mutate(
         cat1f = as.factor(cat1f)
       , cat2f = as.factor(cat2f)
     )
   
     if (method == "mclust") {
      
       y %<>% mutate(ww = ifelse(V1 > V2_adj, V1, V2_adj))
       
       no_variance <- glm(
           formula = gam_formula %>% as.formula()
         , family = "binomial"
         , data   = y
       )
      
       with_variance <- glm(
           formula = gam_formula %>% as.formula()
         , family  = "binomial"
         , weights = ww
         , data    = y
       )
       
       return(list(no_variance, with_variance))
       
     } else if (method == "sd") {
       
       no_variance <- glm(
           formula = gam_formula %>% as.formula()
         , family = "binomial"
         , data   = y
       )
       
       return(list(no_variance))
       
     } else {
       stop("method not supported")
     }
     
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
  
stan_fit <- model_name$compiled_model[[1]]$sample(
  data    = list(
    N = param_set$n_samps
  , y = simulated_data$mfi
  )
, chains  = 4
, parallel_chains = 1
, seed    = 483892929
, refresh = 2000
)

stanfit <- rstan::read_stan_csv(stan_fit$output_files())
samps   <- rstan::extract(stanfit)
samps   <- samps[!grepl("membership_l|ind_sero|lp_", names(samps))]

} else if (model_name$model == "cluster_regression_with_beta_1.stan") {
  
stan_fit <- model_name$compiled_model[[1]]$sample(
  data    = list(
    N           = param_set$n_samps
  , y           = simulated_data$mfi
  , cat1f       = simulated_data$cat1f
  )
, chains  = 4
, parallel_chains = 1
, seed    = 483892929
, refresh = 2000
)

stanfit <- rstan::read_stan_csv(stan_fit$output_files())
samps   <- rstan::extract(stanfit)
samps   <- samps[!grepl("membership_l|ind_sero|log_beta|beta_vec", names(samps))]
  
} else if (model_name$model == "cluster_regression_with_beta_theta_ln_1.stan") {
  
stan_fit <- model_name$compiled_model[[1]]$sample(
  data    = list(
    N           = param_set$n_samps
  , y           = simulated_data$mfi
  , cat1f       = simulated_data$cat1f
  , cat2f       = simulated_data$cat2f
  )
, chains  = 4
, parallel_chains = 1
, seed    = 483892929
, refresh = 2000
)

stanfit <- rstan::read_stan_csv(stan_fit$output_files())
samps   <- rstan::extract(stanfit)
samps   <- samps[!grepl("membership_l|ind_sero|log_beta|beta_vec", names(samps))]
  
} else if (model_name$model == "cluster_regression_with_beta_theta_ln_2.stan") {

stan_fit <- model_name$compiled_model[[1]]$sample(
  data    = list(
    N           = param_set$n_samps
  , N_cat1r     = param_set$cat1r_count
  , y           = simulated_data$mfi
  , cat1f       = simulated_data$cat1f
  , cat2f       = simulated_data$cat2f
  , con1f       = simulated_data$con1f
  , cat1r       = simulated_data$cat1r
  )
, chains  = 4
, parallel_chains = 1
, seed    = 483892929
, refresh = 2000
)

stanfit <- rstan::read_stan_csv(stan_fit$output_files())
samps   <- rstan::extract(stanfit)
samps   <- samps[!grepl("membership_l|ind_sero|log_beta|beta_vec|theta_cat1r_eps", names(samps))]

} else {
  return(NULL)
}

stan_fit        <- list(samps)
names(stan_fit) <- paste(
  model_name$model
, unique(simulated_data$param_set)
, unique(simulated_data$sim_num)
, unique(simulated_data$log_mfi)
, sep = " -- ")

return(stan_fit)
  
}

## Above model fitting function simplified for the one model being fit for the pub
fit_stan_models_for_pub <- function(simulated_data, param_sets, compiled_models, model_names) {
  
  simulated_data <- simulated_data[[1]]
  model_name     <- model_names[[1]] %>% pull(model_base_names)
  param_set      <- param_sets %>% filter(
    param_set == unique(simulated_data$param_set)
  , sim_num   == unique(simulated_data$sim_num)
  )
  
  if (model_name == "publication_model_2.stan") {
    if (unique(simulated_data$log_mfi) == "mfi") {
      model_name <- "publication_model_lognormal_2.stan"
    } else if (unique(simulated_data$log_mfi) == "log_mfi") {
      model_name <- "publication_model_normal_2.stan"
    } else {
      stop("Unknown model name / mfi / log mfi combo")
    }
  }
  
  this_model <- compiled_models %>% filter(model == model_name) %>% pull(compiled_model)
  
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
  
  if (unique(simulated_data$log_mfi) == "log_mfi") {
    prior_dat %<>% dplyr::select(param, log_mfi) %>% rename(prior = log_mfi)
  } else if (unique(simulated_data$log_mfi) == "mfi") {
    prior_dat %<>% dplyr::select(param, mfi) %>% rename(prior = mfi)
  } else {
    stop("data transform not supported")
  }
  
  stan_fit <- this_model[[1]]$sample(
    data    = list(
        N           = length(simulated_data$mfi)
      , N_cat1r     = param_set$cat1r_count
      , y           = simulated_data$mfi
      , cat1f       = simulated_data$cat1f
      , cat2f       = simulated_data$cat2f
      , con1f       = simulated_data$con1f
      , max_mfi     = prior_dat %>% filter(param == "max_mfi") %>% pull(prior)
      
      , mu_base_prior     = prior_dat %>% filter(param == "mu_base_prior") %>% pull(prior)
      , mu_diff_prior     = prior_dat %>% filter(param == "mu_diff_prior") %>% pull(prior)
      , sigma_base_prior  = prior_dat %>% filter(param == "sigma_base_prior") %>% pull(prior)
      , sigma_diff_prior  = prior_dat %>% filter(param == "sigma_diff_prior") %>% pull(prior)
    )
    , chains  = 4
    , parallel_chains = 1
    , max_treedepth = 12
    , iter_warmup = 2000
    , iter_sampling = 1000
    , adapt_delta = 0.96
    , seed    = 483892929
    , refresh = 2000
  )

  stanfit <- rstan::read_stan_csv(stan_fit$output_files())
  samps   <- rstan::extract(stanfit)
  samps   <- samps[!grepl("membership_l|ind_sero|log_beta|beta_vec|theta_cat1r_eps", names(samps))]
  
  stan_fit        <- list(samps)
  names(stan_fit) <- paste(
      model_name
    , unique(simulated_data$param_set)
    , unique(simulated_data$sim_num)
    , unique(simulated_data$log_mfi)
    , sep = " -- ")
  
  return(stan_fit)
  
}

## ^^ Function to deal with what is returned from fit_stan_models_f to get into structure for
 ## how I built summarize_stan_fits.
sort_stan_fits          <- function(stan_fits.l) {
  ## Need a list of lists
   ## outer list is of length n-models
    ## each of these has length = n_param_sets
 model_names <- lapply(stan_fits.l, FUN = function(x) {
    names(x)
  }) %>% unlist()
 
fit_details <- lapply(model_names %>% as.list(), FUN = function(x) {
   strsplit(x, " -- ") %>% unlist() %>% t() %>% as.data.frame()
  }) %>% do.call("rbind", .) %>% 
  rename(model = V1, param_set = V2, sim_num = V3, log_mfi = V4) %>% as_tibble()

resorted_fits <- fit_details %>% mutate(fitted_model = stan_fits.l)
 
return(resorted_fits)
 
}
