gamma_param_solve <- function(a_base, b_base, new_mean, new_var) {
  test_vals <- expand.grid(
    a_diff = seq(0, 12, by = 0.01)
  , b_diff = seq(-4, 4, by = 0.01)
  ) %>% mutate(
    new_a  = a_base + a_diff
  , new_b  = b_base + b_diff
  ) %>% mutate(
    desired_mean = new_mean
  , desired_var  = new_var
  ) %>% mutate(
    calc_mean    = new_a * new_b
  , calc_var     = new_a * new_b^2
  ) %>% mutate(
    diff_mean    = desired_mean - calc_mean
  , diff_var     = desired_var  - calc_var
  ) %>% mutate(
    joint_min    = abs(diff_mean) + abs(diff_var)
  , joint_min_p  = joint_min / min(joint_min)
  )
  
  test_vals %>% arrange(joint_min_p) %>% slice(1) 
  
}
new_gamma_parms   <- gamma_param_solve(a_base = 1, b_base = 1, new_mean = 5, new_var = 4)

new_gamma_parms %<>% dplyr::select(new_a, new_b) %>% unlist()

rgamma(
    10000
  , shape = 1
  , scale = 1
  ) %>% var()

hist(
  rgamma(
    10000
  , shape = new_gamma_parms[1]
  , scale = new_gamma_parms[2]
  )
, breaks = 50
)

mu_vec <- matrix()

MASS::mvrnorm(
  n  = 100
, mu = matrix()
)


x <- param_sets[[1]]
x$mu_neg <- 1
x$sd_neg <- 1
x$mu_pos <- 5
x$sd_pos <- 1
x$theta_cat1r_sd <- 0.5
x$n_samps <- 5000

mu_vec   <- with(x, c(mu_neg, mu_pos))
sd_vec   <- with(x, c(sd_neg, sd_pos))
rand_dev <- rnorm(x$cat1r_count, 0, x$theta_cat1r_sd)

tt <- data.frame(
    cat1f  = with(x, rbinom(n_samps, 1, cat1f_prop))
  , cat2f  = with(x, rbinom(n_samps, 1, cat2f_prop))
  ) %>% mutate(
    cat1r  = with(x, sample(seq(cat1r_count), n(), replace = T))
  , con1f  = with(x, rnorm(n_samps, 0, con1f_sd))
  ) %>% 
  mutate( 
    cat1r_dev = rand_dev[cat1r]
  ) %>% 
  mutate(
    group = rbinom(n(), 1
                   , plogis(
                      x$beta_base + 
                      x$beta_cat1f_delta * cat1f + 
                      con1f * x$beta_con1f_delta
                     )
                   ) + 1
  , mfi   = rgamma(n()
                  , mu_vec[group] + 
                    (cat2f * x$theta_cat2f_mu * (group - 1)) + 
                    cat1r_dev
                  , sd_vec[group]
                  )
  ) %>% mutate(
      param_set = x$param_set
    , sim_num   = x$sim_num
    , .before   = 1
    )

tt %>% mutate(
  cat1r_dev = round(cat1r_dev, 2)
, cat1r_dev = as.factor(cat1r_dev)
, group     = as.factor(group)
) %>% {
  ggplot(., aes(x = mfi)) + 
    geom_density(aes(fill = group), alpha = 0.2) +
    facet_wrap(~cat1r_dev)
}
