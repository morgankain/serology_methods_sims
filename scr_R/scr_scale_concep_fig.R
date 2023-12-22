param_set <- 1
sim_num   <- 1
n_samps   <- 10000
logit_1   <- 30000
logit_2   <- 0
logit_3   <- 1.75
beta_base <- -2

## 1
#mu_neg       <- -1
#sd_neg       <- 0.5
#mu_pos_delta <- 1.5
#sd_pos_delta <- 0.5
 
## 2
#mu_neg       <- -2
#sd_neg       <- 0.5
#mu_pos_delta <- 2.5
#sd_pos_delta <- 1.5

## 3
mu_neg       <- 0.1
sd_neg       <- 0.65
mu_pos_delta <- 1.4
sd_pos_delta <- 0.3

beta_cat1f_delta <- 0 
beta_con1f_delta <- 0
cat1r_count      <- 10
theta_cat1r_sd   <- 1
con1f_sd         <- 2
theta_cat2f_mu   <- 0
cat1f_prop       <- 0.5
cat2f_prop       <- 0.5

mu_vec <- c(mu_neg, mu_neg + mu_pos_delta)
sd_vec <- c(sd_neg, sd_neg * sd_pos_delta)

simulated_data <- data.frame(
    cat1f  = rbinom(n_samps, 1, cat1f_prop)
  , cat2f  = rbinom(n_samps, 1, cat2f_prop)
) %>% mutate(
  group = rbinom(n(), 1
                 , plogis(
                   beta_base + beta_cat1f_delta * cat1f
                 )
  ) + 1
  , titer = rlnorm(n()
                  , mu_vec[group] + 
                    (cat2f * theta_cat2f_mu * (group - 1))
                  , sd_vec[group]
  ) -3
  , mfi   = logit2(logit_1, logit_2, logit_3, titer)
  , log_mfi = log10(mfi)
) %>% mutate(
  param_set = param_set
  , sim_num   = sim_num
  , .before   = 1
)

simulated_data.3 <- simulated_data

simulated_data <- rbind(
  simulated_data.1
, simulated_data.2
, simulated_data.3
)

sig_curve <- data.frame(
  titer = seq(-4, 4, by = 0.1)
, MFI   = logit2(logit_1, logit_2, logit_3, seq(-4, 4, by = 0.1))
)

gg.1 <- simulated_data %>% 
  mutate(group = as.factor(group)) %>% {
  ggplot(., aes(x = titer)) + 
    #geom_density(aes(colour = group, fill = group), alpha = 0.3) +
    geom_histogram(aes(y = (..count..)/sum(..count..), colour = group, fill = group), alpha = 0.3, bins = 50, position = "identity") +
    scale_colour_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
    scale_fill_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
    scale_x_continuous(lim = c(-4, 4)) +
    xlab("Titer") + ylab("Frequency")
}

gg.2 <- sig_curve %>% {
  ggplot(., aes(titer, MFI)) + 
    geom_line() +
    theme(plot.margin = margin(0.3, 3, 0.3, 0.3, "cm")) +
    xlab("Titer") + ylab("MFI")
}

gg.3 <- simulated_data %>% mutate(group = as.factor(group)) %>% {
  ggplot(., aes(x = mfi)) + 
    geom_histogram(aes(y = (..count..)/sum(..count..), colour = group, fill = group), alpha = 0.3, bins = 50, position = "identity") +
    scale_colour_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
    scale_fill_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
    xlab("MFI") + ylab("Frequency")
}

gg.4 <- simulated_data %>% mutate(group = as.factor(group)) %>% {
  ggplot(., aes(x = log_mfi)) + 
    geom_histogram(aes(y = (..count..)/sum(..count..), colour = group, fill = group), alpha = 0.3, bins = 50, position = "identity") +
    scale_colour_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
    scale_fill_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
    xlab("log10(MFI)") + ylab("Frequency")
}

gridExtra::grid.arrange(gg.1, gg.2, gg.3, gg.4, ncol = 1)

gg.1 <- simulated_data %>% 
  mutate(
    group = as.factor(group)
  , int = interaction(group, param_set)) %>% {
    ggplot(., aes(x = titer)) + 
      geom_histogram(aes(y = (..count..)/sum(..count..), colour = group, fill = group), alpha = 0.3, bins = 50, position = "identity") +
      scale_colour_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
      scale_fill_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
      scale_x_continuous(lim = c(-4, 4)) +
      xlab("Titer") + ylab("Frequency") + 
      facet_wrap(~param_set) +
      theme(strip.text.x = element_blank())
  }

sig_curve <- rbind(
  sig_curve %>% mutate(param_set = 1)
, sig_curve %>% mutate(param_set = 2)
, sig_curve %>% mutate(param_set = 3)
)

gg.2 <- sig_curve %>% {
  ggplot(., aes(titer, MFI)) + 
    geom_line() +
    theme(plot.margin = margin(0.3, 3, 0.3, 0.3, "cm")
        , strip.text.x = element_blank()) +
    xlab("Titer") + ylab("MFI") +
    facet_wrap(~param_set) +
    scale_y_continuous(
      breaks = c(0, 10000, 20000, 30000)
    , labels = c("0", "10k", "20k", "30k")
    ) 
}

gg.3 <- simulated_data %>% mutate(group = as.factor(group)) %>% {
  ggplot(., aes(x = mfi)) + 
    geom_histogram(aes(y = (..count..)/sum(..count..), colour = group, fill = group), alpha = 0.3, bins = 50, position = "identity") +
    scale_colour_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
    scale_fill_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
    xlab("MFI") + ylab("Frequency") + 
    scale_y_continuous(breaks = c(0, 0.03, 0.06, 0.09, 0.12)) +
    scale_x_continuous(
      breaks = c(0, 10000, 20000, 30000)
      , labels = c("0", "10k", "20k", "30k")
    ) +
    theme(strip.text.x = element_blank()) +
    facet_wrap(~param_set)
}

gg.4 <- simulated_data %>% mutate(group = as.factor(group)) %>% {
  ggplot(., aes(x = log_mfi)) + 
    geom_histogram(aes(y = (..count..)/sum(..count..), colour = group, fill = group), alpha = 0.3, bins = 50, position = "identity") +
    scale_colour_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
    scale_fill_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
    xlab("log10(MFI)") + ylab("Frequency") + 
    theme(strip.text.x = element_blank()) +
    facet_wrap(~param_set)
}

gridExtra::grid.arrange(gg.1, gg.2, gg.3, gg.4, ncol = 1)


simulated_data %>% group_by(
  param_set, group
) %>% summarize(
  skew = skewness(mfi)
)


