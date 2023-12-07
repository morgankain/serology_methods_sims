complexity <- 1

param_set <- 1
sim_num   <- 1

n_samps <- 1000

logit_1 <- 30000
logit_2 <- -1
logit_3 <- 1

beta_base        <- -2

beta_cat1f_delta <- 0 #0.980
beta_con1f_delta <- 0.5

cat1r_count <- 10
theta_cat1r_sd <- 1

con1f_sd    <- 2

theta_cat2f_mu <- 2

cat1f_prop <- 0.5
cat2f_prop <- 0.5

mu_neg       <- -3
sd_neg       <- 1
mu_pos_delta <- 2
sd_pos_delta <- 1

mu_vec <- c(mu_neg, mu_neg + mu_pos_delta)
sd_vec <- c(sd_neg, sd_neg * sd_pos_delta)
  
if (complexity == 1) {

simulated_data <- data.frame(
    cat1f  = rbinom(n_samps, 1, cat1f_prop)
  , cat2f  = rbinom(n_samps, 1, cat2f_prop)
 ) %>% mutate(
    group = rbinom(n(), 1
                   , plogis(
                      beta_base + beta_cat1f_delta * cat1f
                   )
    ) + 1
    , titer = rnorm(n()
                    , mu_vec[group] + 
                      (cat2f * theta_cat2f_mu * (group - 1))
                    , sd_vec[group]
    )
    , mfi   = logit2(logit_1, logit_2, logit_3, titer)
  ) %>% mutate(
      param_set = param_set
    , sim_num   = sim_num
    , .before   = 1
  )

} else {
  
rand_dev <- rnorm(cat1r_count, 0, theta_cat1r_sd)
  
simulated_data <- data.frame(
      cat1f  = rbinom(n_samps, 1, cat1f_prop)
    , cat2f  = rbinom(n_samps, 1, cat2f_prop)
  ) %>% mutate(
      cat1r  = sample(seq(cat1r_count), n(), replace = T)
    , con1f  = rnorm(n_samps, 0, con1f_sd)
  ) %>% 
    mutate( 
      cat1r_dev = rand_dev[cat1r]
    ) %>% 
    mutate(
      group = rbinom(n(), 1
                     , plogis(
                        beta_base + beta_cat1f_delta * cat1f + con1f * beta_con1f_delta
                       )
      ) + 1
      , titer = rnorm(n()
                      , mu_vec[group] + 
                        (cat2f * theta_cat2f_mu * (group - 1)) +
                        cat1r_dev
                      , sd_vec[group]
      )
      , mfi   = logit2(logit_1, logit_2, logit_3, titer)
    ) %>% mutate(
        param_set = param_set
      , sim_num   = sim_num
      , .before   = 1
    )
}

simulated_data %>% group_by(cat1f, group) %>% summarize(n())
#simulated_data %>% {
#  ggplot(., aes(con1f, group)) + geom_jitter()
#}
simulated_data %>% filter(group == 2) %>% {
  ggplot(., aes(cat2f, mfi)) + geom_jitter()
}

simulated_data %>% mutate(
    group = as.factor(group)
  , int   = interaction(cat1f, group)) %>% {
    ggplot(., aes(x = mfi)) + 
      geom_histogram(aes(colour = int, fill = int), alpha = 0.3, bins = 50) + 
      scale_colour_brewer(palette = "Dark2") +
      scale_fill_brewer(palette = "Dark2") +
      facet_wrap(~cat1f)
  }

simulated_data %>% mutate(
    group = as.factor(group)
  , int   = interaction(cat1f, group)) %>% 
  mutate(
    mfi = log(mfi)
  ) %>% {
    ggplot(., aes(x = mfi)) + 
      geom_histogram(aes(colour = int, fill = int), alpha = 0.3, bins = 50) + 
      scale_colour_brewer(palette = "Dark2") +
      scale_fill_brewer(palette = "Dark2") +
      facet_grid(sim_num~param_set)
  }

