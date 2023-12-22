param_set <- 2
sim_num   <- 1
n_samps   <- 10000
logit_1   <- 30000
logit_2   <- 0
logit_3   <- 1.75
beta_base <- -2

## 1
mu_neg       <- -2
sd_neg       <- 0.5
mu_pos_delta <- 1.5
sd_pos_delta <- 0.5

## 2
mu_neg       <- -3
sd_neg       <- 0.5
mu_pos_delta <- 2
sd_pos_delta <- 0.5

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
    , titer = rnorm(n()
                    , mu_vec[group] + 
                      (cat2f * theta_cat2f_mu * (group - 1))
                    , sd_vec[group]
    )
    , mfi   = logit2(logit_1, logit_2, logit_3, titer)
    , log10_mfi = log10(mfi)
    , log_mfi   = log(mfi)
  ) %>% mutate(
    param_set = param_set
    , sim_num   = sim_num
    , .before   = 1
  )
 
sd_thresh <- simulated_data %>% filter(
  group == 1
) %>% slice(
  sample(seq(n()), 3000)
) %>% dplyr::select(mfi, log_mfi, log10_mfi) %>%
  as.matrix() %>% apply(., 2, FUN = function(x) mean(x) + 3 * sd(x))

sd_thresh2 <- sd_thresh

simulated_data %<>% mutate(
  mfi_thresh       = ifelse(mfi > sd_thresh[1], 2, 1)
, mfi_log_thresh   = ifelse(log_mfi > sd_thresh[2], 2, 1)
, mfi_log10_thresh = ifelse(log10_mfi > sd_thresh[3], 2, 1)
)

simulated_data %<>% mutate(
  mfi_compare = ifelse(
  mfi_thresh == group, 0
, ifelse(mfi_thresh > group, 1, -1)
  )
, mfi_log_compare = ifelse(
  mfi_log_thresh == group, 0
, ifelse(mfi_log_thresh > group, 1, -1)
)
, mfi_log10_compare = ifelse(
  mfi_log10_thresh == group, 0
, ifelse(mfi_log10_thresh > group, 1, -1)
)
) 

sd_thresh <- rbind(
  sd_thresh1
  , sd_thresh2
) %>% as.data.frame() %>% mutate(
  param_set = c(1, 2)
)

simulated_data.2 <- simulated_data

simulated_data <- rbind(
  simulated_data.1, simulated_data.2
)

gg.1 <- simulated_data %>% 
  mutate(group = as.factor(group)) %>% {
    ggplot(., aes(x = mfi)) + 
      geom_histogram(aes(y = (..count..)/sum(..count..), colour = group, fill = group), alpha = 0.3, bins = 50, position = "identity") +
      scale_colour_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
      scale_fill_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
      xlab("MFI") + ylab("Frequency") +
      theme(
        strip.text.x = element_blank()
      , plot.margin = margin(0.3, 0.1, 0.3, 0.5, "cm")
      ) +
      geom_vline(data = sd_thresh, aes(xintercept = mfi), linetype = "dashed") +
      facet_wrap(~param_set) 
  }

gg.2 <- simulated_data %>% 
  mutate(group = as.factor(group)) %>% {
    ggplot(., aes(x = log_mfi)) + 
      geom_histogram(aes(y = (..count..)/sum(..count..), colour = group, fill = group), alpha = 0.3, bins = 50, position = "identity") +
      scale_colour_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
      scale_fill_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
      xlab("log(MFI)") + ylab("Frequency") +
      theme(strip.text.x = element_blank()) +
      geom_vline(data = sd_thresh, aes(xintercept = log_mfi), linetype = "dashed") +
      facet_wrap(~param_set) 
  }

gg.3 <- simulated_data %>% 
  mutate(group = as.factor(group)) %>% {
    ggplot(., aes(x = log10_mfi)) + 
      geom_histogram(aes(y = (..count..)/sum(..count..), colour = group, fill = group), alpha = 0.3, bins = 50, position = "identity") +
      scale_colour_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
      scale_fill_brewer(palette = "Dark2", name = "Serostatus", labels = c("Negative", "Positive")) +
      xlab("log10(MFI)") + ylab("Frequency") +
      theme(strip.text.x = element_blank()) +
      geom_vline(data = sd_thresh, aes(xintercept = log10_mfi), linetype = "dashed") +
      facet_wrap(~param_set) 
  }

gridExtra::grid.arrange(gg.1, gg.2, gg.3, ncol = 1)

simulated_data %>% 
  group_by(param_set, group, mfi_log_compare) %>% 
  summarize((estimate = n() / 10000) %>% round(3))

simulated_data %>% filter(
  mfi_compare == -1, group == 2
)

