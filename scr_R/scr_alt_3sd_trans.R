####
## Explore many alternative strategies for defining the negative control distribution
## for the 3sd method of determining positives and negatives
####

## Use alt function to test for groupings 
all_3sd_groupings <- group_via_3sd_alt(
  simulated_data = fits$sim.data
  , param_sets   = fits$sim.params
  , groupings    = c("param_set", "sim_num", "log_mfi")
)

all_3sd_groupings.s <- all_3sd_groupings %>%
  dplyr::select(-assigned_group) %>%
  mutate(
      true_pos  = ifelse(group == 0 & V2 == 0, 1, 0)
    , true_neg  = ifelse(group == 1 & V2 == 1, 1, 0)
    , false_pos = ifelse(group == 0 & V2 == 1, 1, 0)
    , false_neg = ifelse(group == 1 & V1 == 1, 1, 0)
  ) %>%
  group_by(model, param_set, sim_num, sd_method, log_mfi, perc_sd, group) %>%
  summarize(
      prob        = mean(V2)
    , false_pos_p = length(which(false_pos == 1)) / n() 
    , false_neg_p = length(which(false_neg == 1)) / n()
    , true_pos_p  = length(which(true_pos == 1)) / n()
    , true_neg_p  = length(which(true_neg == 1)) / n()
  ) %>% mutate(
      misclass_error_p = false_pos_p + false_neg_p
    , correct_class_p  = 1 - misclass_error_p
  ) %>% mutate(
    quantile = "mid", .before = prob
  ) %>% left_join(., fits$sim.params %>% dplyr::select(
    param_set, sim_num, n_samps, beta_base, mu_neg, sd_neg
    , mu_pos, sd_pos, mu_pos_delta, sd_pos_delta
  ))

all_3sd_groupings.p <- all_3sd_groupings %>% 
  group_by(param_set, sim_num, sd_method, perc_sd, log_mfi) %>% 
  summarize(
      true     = mean(group)
    , prop_pos = mean(assigned_group)
  ) %>% mutate(
    prop_pos_diff = prop_pos - true 
  ) %>% ungroup() %>% 
  mutate(model = "3sd", .before = 1) %>%
  mutate(quantile = "mid", .after = "true") %>%
  rename(method = sd_method) %>% 
  left_join(
    .
  , fits$sim.params %>% dplyr::select(
    param_set, sim_num, n_samps, beta_base, mu_neg, sd_neg
    , mu_pos, sd_pos, mu_pos_delta, sd_pos_delta
  )
  )

all_3sd_groupings.s %>% 
  filter(param_set == 1) %>%
  filter(group == 0) %>% {
  ggplot(., aes(perc_sd, correct_class_p)) + 
      geom_line(aes(colour = log_mfi))
  }

all_3sd_groupings.p %>% filter(
  perc_sd %in% seq(0, 1, by = 0.1)
) %>% {
    ggplot(., aes(true, prop_pos)) + 
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    facet_grid(perc_sd~log_mfi) +
    theme(axis.text.y = element_text(size = 8)) +
    xlab("True") +
    ylab("Estimated")
}

all_3sd_groupings.s %>% 
  left_join(., all_3sd_groupings.p %>% 
              dplyr::select(param_set, perc_sd, log_mfi, true)
  ) %>%
filter(
  perc_sd %in% seq(0, 1, by = 0.1)
) %>% filter(
  group == 0
) %>% {
  ggplot(., aes(true, correct_class_p)) + 
    geom_point() +
    facet_grid(perc_sd~log_mfi) +
    theme(axis.text.y = element_text(size = 8)) +
    xlab("True") +
    ylab("Estimated")
}

all_3sd_groupings.p %>% 
  mutate(true = plyr::round_any(true, 0.1)) %>%
  group_by(perc_sd, log_mfi, true) %>%
  summarize(
    m_error = mean(prop_pos_diff)
  ) %>% mutate(
    perc_sd = as.factor(perc_sd)
  ) %>% filter(
    !(perc_sd == 1 & true == 0.8)
  ) %>% {
    ggplot(., aes(true, m_error)) +
      geom_line(aes(colour = perc_sd)) +
      geom_hline(yintercept = 0) +
      xlab("True Proportion Seropositive (rounded to nearest 0.1)") +
      ylab("Bias in Serostatus Estimate") +
      facet_wrap(~log_mfi) 
  }

all_3sd_groupings.s %>% 
  left_join(., all_3sd_groupings.p %>% 
              dplyr::select(param_set, perc_sd, log_mfi, true)
            ) %>%
  mutate(true = plyr::round_any(true, 0.1)) %>%
  group_by(perc_sd, log_mfi, true, group) %>%
  summarize(
    m_error = mean(correct_class_p)
  ) %>% mutate(
    perc_sd = as.factor(perc_sd)
  ) %>% filter(
    !(perc_sd == 1 & true == 0.8)
  ) %>% {
    ggplot(., aes(true, m_error)) +
      geom_line(aes(colour = perc_sd)) +
      xlab("True Proportion Seropositive (rounded to nearest 0.1)") +
      ylab("Probability of Correct Individual Assignment") +
      facet_wrap(~group+log_mfi) 
  }
 