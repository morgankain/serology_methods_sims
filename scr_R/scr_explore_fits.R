fits.out <- readRDS("test_fits_v4.Rds")
all.out <- fits.out[[1]]
sim.data.summaries <- fits.out[[2]]
sim.params <- fits.out[[3]]
sim.data <- fits.out[[4]]
stan_models <- fits.out[[5]]
stan.summary <- fits.out[[6]]
mculst.groups <- fits.out[[7]]
three_sd.groups <- fits.out[[8]]

all.out$group_assignment <- all.out$group_assignment %>% 
  left_join(., all.out$pop_seropositivity %>% 
              dplyr::select(param_set, true)  %>% 
              group_by(param_set) %>% slice(1)
  )

sim.data.summaries.j <- sim.data.summaries %>%
  dplyr::select(
    param_set, sim_num, log_mfi
  , prop_overlap_quant, diff_mean, mfi_skew_2
  , prop_squish_2, prop_pos
  ) %>% rename(prop_pos_true = prop_pos)

## Add summary info to each summary
all.out$group_assignment <- all.out$group_assignment %>% 
  left_join(., sim.data.summaries.j)
all.out$pop_seropositivity <- all.out$pop_seropositivity %>% 
  left_join(., sim.data.summaries.j)
all.out$coefficient_ests <- all.out$coefficient_ests %>% 
  left_join(., sim.data.summaries.j)
all.out$coverage <- all.out$coverage %>% 
  left_join(., sim.data.summaries.j)

## Explore variation in summaries of the distributions
sim.data.summaries %>% ungroup() %>% 
  filter(log_mfi == "log_mfi") %>% dplyr::select(
  param_set, prop_pos, diff_mean, prop_overlap_quant
, mfi_skew_2, prop_squish_2, p_abv_m_2
  ) %>% {
  ggpairs(.) +
    theme(
      axis.text.x = element_text(size = 10, hjust = 0, angle = 300)
    , axis.text.y = element_text(size = 10)
    , strip.text.x = element_text(size = 10)
    )
  }

sim.data %>% 
  filter(param_set == 190) %>%
  mutate(
    group = as.factor(group)
    , int = interaction(cat1f, group)
  ) %>% {
    ggplot(., aes(x = mfi)) + 
      geom_density(aes(colour = group, fill = group), alpha = 0.2) + 
      scale_colour_brewer(palette = "Dark2") +
      scale_fill_brewer(palette = "Dark2") +
      facet_wrap(~log_mfi, scales = "free", nrow = 3)
  }

## What are the x-axes of interest?
 ## mfi_skew_2 -- skew in right distribution
 ## prop_pos -- proportion positive
 ## prop_overlap_quant -- proportion of positive distribution below upper 95% quantile of negative

####
## Plot for an individual simulation to look at distributions
####

## 68, 26, 170, 160, 130

sim.data %>% 
  filter(param_set == 26, log_mfi == "log_mfi") %>%
  pull(mfi) %>%
  dplR::tbrm()

sim.data %>% 
  filter(param_set == 26, log_mfi == "log_mfi") %>%
  pull(mfi) %>%
  jointseg::estimateSd()
  
sim.data %>% 
  filter(param_set == 126) %>%
  mutate(
      group = as.factor(group)
    , int = interaction(cat1f, group)
  ) %>% {
    ggplot(., aes(x = mfi)) + 
      geom_density(aes(colour = group, fill = group), alpha = 0.2) + 
      scale_colour_brewer(palette = "Dark2") +
      scale_fill_brewer(palette = "Dark2") +
      facet_wrap(~log_mfi, scales = "free", nrow = 3)
  }

####
## Individual-level probability of assignment 
####

all.out$group_assignment %>% 
  filter(model == "3sd") %>% 
  mutate(group = as.factor(group)) %>% {
    ggplot(., aes(true, correct_class_p)) + 
      geom_point(aes(colour = group)) +
      geom_abline(intercept = 0, slope = 1) +
      facet_grid(log_mfi~method)
  }

fits$all.out$group_assignment %>% 
  filter(model == "3sd") %>% 
  mutate(group = as.factor(group)) %>%
  

####
## Population-level probability of assignment 
####

names(fits$all.out$group_assignment)

## beta_base
## mu_neg
## sd_neg
## mu_pos
## sd_pos
## mu_pos_delta
## sd_pos_delta

all.out$pop_seropositivity %>%
  filter(model == "3sd", method == "assigned_group_robust") %>%
  arrange(prop_pos_diff)

all.out$pop_seropositivity %>% 
  filter(model == "3sd") %>% {
    ggplot(., aes(true, prop_pos)) + 
      geom_point() +
      geom_abline(intercept = 0, slope = 1) +
      facet_grid(log_mfi~method)
  }

## mfi_skew_2
## prop_pos_true
## prop_overlap_quant

all.out$pop_seropositivity %>% 
  filter(model == "3sd") %>% {
    ggplot(., aes(prop_pos_true, prop_pos_diff)) + 
      geom_point() +
      geom_hline(yintercept = 0) +
      facet_grid(log_mfi~method)
  }

all.out$pop_seropositivity %>% 
  filter(model == "3sd") %>%
  dplyr::select(log_mfi, param_set, method, true, prop_pos_diff) %>% 
  pivot_wider(names_from = log_mfi, values_from = prop_pos_diff) %>% {
    ggplot(., aes(mfi, log_mfi)) + 
      geom_point(aes(colour = method)) +
      geom_point(
        data = all.out$pop_seropositivity %>% 
          filter(model == "3sd") %>%
          dplyr::select(log_mfi, param_set, method, true, prop_pos_diff) %>% 
          pivot_wider(names_from = log_mfi, values_from = prop_pos_diff) %>%
          filter(param_set %in% c(126, 190))
      , aes(mfi, log_mfi), size = 5, shape = c(0, 0, 2, 2)
      ) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_hline(yintercept = 0, linetype = "dotted") +
      xlab("MFI") +
      ylab("Log(MFI)") +
      scale_colour_brewer(palette = "Dark2", name = "Method")
  }

all.out$pop_seropositivity %>% 
  filter(model == "3sd") %>%
  filter(method == "assigned_group_control") %>% 
  dplyr::select(log_mfi, param_set, sd_neg, true, prop_pos) %>% 
  pivot_wider(names_from = log_mfi, values_from = prop_pos) %>% {
    ggplot(., aes(sd_neg, log_mfi)) + 
      geom_point(colour = "dodgerblue3") +
      geom_point(aes(sd_neg, mfi), colour = "firebrick3") +
      geom_errorbar(aes(ymin = mfi, ymax = log_mfi), width = 0, alpha = 0.5)
  }

test_for_gg <- all.out$pop_seropositivity %>% 
  filter(model == "3sd") %>%
  filter(method == "assigned_group_control") %>%
  mutate(
    gg = ifelse(
      prop_pos_diff > -0.001 & prop_pos_diff < 0.001
    , 1, 0)
  ) %>% 
  filter(gg == 1) %>% 
  group_by(param_set) %>%
  mutate(n_entry = n()) %>% 
  ungroup()

better_log_option_A <- test_for_gg %>% 
  filter(n_entry == 1) %>%
  filter(log_mfi == "log_mfi")

better_log_option_B <- all.out$pop_seropositivity %>% 
  filter(model == "3sd") %>%
  filter(method == "assigned_group_control") %>% 
  dplyr::select(log_mfi, param_set, true, prop_pos) %>% 
  pivot_wider(names_from = log_mfi, values_from = prop_pos) %>%
  mutate(
    log_mfi = abs(log_mfi - true)
  , mfi     = abs(mfi - true)
  ) %>% mutate(diff = log_mfi - mfi) %>%
  arrange(desc(diff)) %>%
  left_join(
    .
  , all.out$pop_seropositivity %>% 
    filter(model == "3sd") %>%
    filter(method == "assigned_group_control") %>%
    filter(log_mfi == "mfi") %>%
    dplyr::select(
      param_set, true, n_samps, beta_base, mu_neg, sd_neg, mu_pos
    , sd_pos, mu_pos_delta, sd_pos_delta
    )
  )
  
better_log_option_B %>% {
  ggplot(., aes(mu_neg, diff)) + 
    geom_point() +
    geom_hline(yintercept = 0)
}

better_log_option_B %>% 
  mutate(gg = ifelse(log_mfi < 0.002, 1, 0) %>% as.factor()) %>% 
  {
  ggplot(., aes(sd_neg, mu_pos_delta)) + 
    geom_point(aes(colour = log_mfi, size = gg, shape = gg)) +
    scale_size_manual(values = c(1, 2)) 
  }

all.out$pop_seropositivity %>% 
  filter(model == "3sd") %>%
  filter(method == "assigned_group_control") %>%
  mutate(
    gg = ifelse(
      prop_pos_diff > -0.001 & prop_pos_diff < 0.001
      , 1, 0)
  ) %>% {
     ggplot(., aes(sd_neg, gg)) + 
        geom_jitter() +
        facet_wrap(~log_mfi, nrow = 1)
  }

## mfi_skew_2
## prop_pos_true
## prop_overlap_quant

all.out$pop_seropositivity %>% 
  filter(model == "3sd") %>%
  filter(method == "assigned_group_robust") %>%
  filter(prop_pos_true < 0.1) %>% {
    ggplot(., aes(prop_overlap_quant, prop_pos_diff)) + 
                  geom_point(aes(colour = log_mfi), size = 3) +
      scale_colour_brewer(palette = "Dark2")
              } 

#### ------

three_sd.groups %>% 
  filter(param_set == 1) %>% {
    ggplot(., aes(group, assigned_group)) + 
      geom_jitter(width = 0.2, height = 0.2) +
      scale_x_continuous(breaks = c(0, 1)) +
      scale_y_continuous(breaks = c(0, 1)) +
      facet_wrap(~sd_method+log_mfi, scales = "free_x") +
      xlab("True Group") +
      ylab("Assigned Group")
  }

## 126, 190

three_sd.groups %>% 
 filter(param_set == 190) %>% 
  mutate(group = as.factor(group)) %>% {
    ggplot(., aes(mfi, assigned_group)) + 
      geom_jitter(
        aes(colour = group)  
      , width = 0.2, height = 0.2) +
      scale_y_continuous(breaks = c(0, 1)) +
      facet_grid(sd_method~log_mfi, scales = "free_x") +
      ylab("True Group") +
      xlab("MFI") +
      scale_colour_manual(
        values = c("dodgerblue3", "firebrick3")
      )
  }

fits$all.out$group_assignment %>% 
  filter(model == "3sd") %>%
  filter(method == "assigned_group_sample_mfi") %>% {
    ggplot(., aes(true, prop_pos_diff)) + geom_point() +
      geom_hline(yintercept = 0)
  }

fits$all.out$group_assignment %>% 
  filter(model == "3sd") %>%
  filter(method == "assigned_group_sample_mfi") %>% {
    ggplot(., aes(sd_pos_delta, prop_pos_diff)) + geom_point() +
      geom_hline(yintercept = 0)
  }

fits$all.out$group_assignment %>% 
  filter(model == "3sd") %>%
  filter(method == "assigned_group_sample_mfi") %>%
  dplyr::select()

fits$all.out$group_assignment %>% 
  filter(model == "3sd") %>%
  filter(method == "assigned_group_control_mfi") %>% {
    ggplot(., aes(true, prop_pos)) + geom_point() +
      geom_abline(intercept = 0, slope = 1)
  }

fits$all.out$group_assignment %>% 
  filter(model == "3sd") %>%
  filter(method == "assigned_group_control_mfi") %>% {
    ggplot(., aes(true, prop_pos_diff)) + geom_point() +
      geom_hline(yintercept = 0)
  }

#### ------

## mfi_skew_2
## prop_pos_true
## prop_overlap_quant

all.out$pop_seropositivity %>% 
  filter(model == "mclust") %>%
  filter(method == "constrained_mclust") %>% {
    ggplot(., aes(true, prop_pos)) + geom_point() +
      geom_abline(intercept = 0, slope = 1)
  }

# "constrained_mclust"
# "unconstrained_mclust"
# "unconstrained_reduced_mclust"

all.out$pop_seropositivity %>% 
  filter(model == "mclust") %>%
  filter(method == "constrained_mclust") %>% {
    ggplot(., aes(true, prop_pos)) + 
      geom_point() +
      geom_abline(intercept = 0, slope = 1) +
      facet_grid(log_mfi~method)
  }

all.out$pop_seropositivity %>% 
  filter(model == "mclust") %>%
  filter(method == "constrained_mclust") %>% {
    ggplot(., aes(true, prop_pos_diff)) + 
      geom_point(aes(colour = log_mfi)) +
      scale_colour_brewer(palette = "Dark2") +
      geom_hline(yintercept = 0)
  }

## mfi_skew_2
## prop_pos_true
## prop_overlap_quant

all.out$pop_seropositivity %>% 
  filter(model == "mclust") %>%
  filter(method == "constrained_mclust") %>% {
    ggplot(., aes(prop_overlap_quant, prop_pos_diff)) + 
      geom_point(aes(colour = log_mfi)) +
      scale_colour_brewer(palette = "Dark2") +
      geom_hline(yintercept = 0)
  }

all.out$pop_seropositivity %>% 
  filter(model == "mclust") %>% {
    ggplot(., aes(prop_pos_true, prop_pos)) + 
      geom_point(aes(colour = prop_overlap_quant)) +
      geom_abline(intercept = 0, slope = 1) +
      scale_colour_viridis() +
      facet_grid(log_mfi~method)
  }

all.out$pop_seropositivity %>% 
  filter(model == "mclust") %>% {
    ggplot(., aes(prop_overlap_quant, prop_pos_diff)) + 
      geom_point() +
      geom_hline(yintercept = 0) +
      scale_colour_viridis() +
      facet_grid(log_mfi~method)
  }

all.out$pop_seropositivity %>% 
  filter(model == "mclust") %>%
  filter(method %in% c(
    "constrained_mclust", "unconstrained_reduced_mclust"
  )) %>%
  dplyr::select(log_mfi, param_set, method, true, prop_pos_diff) %>% 
  pivot_wider(names_from = method, values_from = prop_pos_diff) %>% {
    ggplot(., aes(constrained_mclust, unconstrained_reduced_mclust)) + 
      geom_point(aes(colour = log_mfi)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_hline(yintercept = 0, linetype = "dotted") +
      xlab("Constrained Mclust") +
      ylab("Unconstrained Reduced Mclust") +
      scale_colour_brewer(palette = "Dark2", name = "Method")
  }

all.out$pop_seropositivity %>% 
  filter(model == "mclust") %>%
  filter(method == "constrained_mclust") %>% 
  dplyr::select(log_mfi, param_set, true, prop_pos_diff) %>% 
  pivot_wider(names_from = log_mfi, values_from = prop_pos_diff) %>% {
    ggplot(., aes(mfi, log_mfi)) + 
      geom_point() +
      geom_abline(intercept = 0, slope = 1)
  }

all.out$pop_seropositivity %>% 
  filter(model == "mclust") %>%
  dplyr::select(log_mfi, param_set, method, true, prop_pos_diff) %>% 
  pivot_wider(names_from = log_mfi, values_from = prop_pos_diff) %>% {
    ggplot(., aes(mfi, log_mfi)) + 
      geom_point(aes(colour = method)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_hline(yintercept = 0, linetype = "dotted") +
      xlab("MFI") +
      ylab("Log(MFI)") +
      scale_colour_brewer(palette = "Dark2", name = "Method")
  }

all.out$pop_seropositivity %>% 
  filter(model == "mclust") %>%
  filter(method == "constrained_mclust") %>% 
  dplyr::select(log_mfi, param_set, sd_neg, true, prop_pos) %>% 
  pivot_wider(names_from = log_mfi, values_from = prop_pos) %>% {
    ggplot(., aes(sd_neg, log_mfi)) + 
      geom_point(colour = "dodgerblue3") +
      geom_point(aes(sd_neg, mfi), colour = "firebrick3") +
      geom_errorbar(aes(ymin = mfi, ymax = log_mfi), width = 0, alpha = 0.5)
  }

all.out$pop_seropositivity %>% 
  filter(model == "mclust") %>%
  filter(method == "unconstrained_mclust") %>% {
    ggplot(., aes(true, prop_pos)) + geom_point() +
      geom_abline(intercept = 0, slope = 1)
  }

all.out$pop_seropositivity %>% 
  filter(model == "mclust") %>%
  filter(method == "unconstrained_mclust") %>% {
    ggplot(., aes(true, prop_pos_diff)) + 
      geom_point() +
      geom_hline(yintercept = 0)
  }

#### ------

stan.fit_stats <- stan.summary$fit_details %>% 
  left_join(., sim.data.summaries.j) %>% 
  group_by(model_name, param_set, log_mfi) %>%
  summarize(
    prop_overlap_quant = mean(prop_overlap_quant)
    , diff_mean = mean(diff_mean)
    , mfi_skew_2 = mean(mfi_skew_2)
    , prop_squish_2 = mean(prop_squish_2)
    , prop_pos_true = mean(prop_pos_true)
    , m_R = mean(max_Rhat)
    , m_D = mean(divergent_transitions)
    , m_T = mean(time_to_fit)
  )

## mfi_skew_2
## prop_pos_true
## prop_overlap_quant

## publication_model_normal_2.stan
## publication_model_lognormal_2.stan
## publication_model_skewnormal_2.stan

stan.summary$fit_details %>% 
  filter(param_set == 133)

all.out$pop_seropositivity %>%
  filter(method == "Bayesian LCR") %>%
  filter(param_set == 133)

all.out$pop_seropositivity %>%
  filter(method == "Bayesian LCR") %>% 
  filter(quantile == "mid") %>% {
    ggplot(., aes(true, prop_pos)) + 
      geom_point(aes(colour = model)) +
      geom_abline(intercept = 0, slope = 1) +
      facet_wrap(~log_mfi, ncol = 1)
  }

all.out$pop_seropositivity %>%
  filter(method == "Bayesian LCR") %>% 
  filter(quantile == "mid") %>% {
    ggplot(., aes(true, prop_pos_diff)) + 
      geom_point(aes(colour = model)) +
      geom_hline(yintercept = 0) +
      facet_wrap(~log_mfi, ncol = 1)
  }

all.out$pop_seropositivity %>%
  filter(method == "Bayesian LCR") %>% 
  filter(quantile == "mid") %>% 
  left_join(
  .
, stan.fit_stats %>% dplyr::select(model_name, param_set, log_mfi, m_R, m_D, m_T) %>% rename(model = model_name)
  ) %>% 
  mutate(m_R = ifelse(is.infinite(m_R), 1000, m_R)) %>% {
    ggplot(., aes(m_D, prop_pos_diff)) + 
      geom_point(aes(colour = model)) +
      geom_hline(yintercept = 0) +
      facet_wrap(~log_mfi, ncol = 1) +
      scale_x_log10()
}

####
## Coefficient estimates from regression model
####

all.out$coefficient_ests %>% 
  filter(method == "Bayesian LCR", log_mfi == "log_mfi", name == "beta_base", !is.na(m_diff)) %>% 
  filter(param_set %notin% failed_fits) %>%
  group_by(model) %>% 
  summarize(
    m_diff = mean(m_diff)
  , cov    = sum(cover)/n()
    )

coeff <- all.out$coefficient_ests %>% 
  mutate(m_diff = ifelse(mid < true, m_diff * -1, m_diff)) %>%
  left_join(
  .
, sim.params
) %>% left_join(
  .
, sim.data.summaries.j
)

coeff %<>% left_join(
  .
, stan.fit_stats %>% 
  dplyr::select(model_name, param_set, log_mfi, m_R, m_D, m_T) %>%
  rename(model = model_name)
) %>% mutate(
  ok_fit = ifelse((is.na(m_R) | m_R < 1.10), 1, 0)
)

# beta_base
# beta_cat1f_delta
# beta_cat2f_delta
# beta_con1f_delta

coef.cover <- coeff %>% 
  mutate(cover = ifelse(is.na(cover), 0, cover)) %>%
  group_by(model, method, log_mfi, name, ok_fit) %>%
  filter(name %notin% c(
    "mu_neg", "mu_pos", "mu_pos_delta"
  , "sd_neg", "sd_pos", "theta_con2f_delta"
  )) %>% 
  summarize(perc_cov = sum(cover) / n()) %>% 
  ungroup() %>%
  arrange(desc(perc_cov)) %>%
  mutate(axx = interaction(model, method, log_mfi))

coef.cover.order <- coeff %>% 
  mutate(cover = ifelse(is.na(cover), 0, cover)) %>%
  group_by(model, method, log_mfi) %>%
  filter(name %notin% c(
    "mu_neg", "mu_pos", "mu_pos_delta"
    , "sd_neg", "sd_pos", "theta_con2f_delta"
  )) %>% 
  summarize(perc_cov = sum(cover) / n()) %>% 
  ungroup() %>%
  arrange(desc(perc_cov)) %>%
  mutate(axx = interaction(model, method, log_mfi))
  
coef.cover %<>% mutate(
  axx = factor(axx, levels = coef.cover.order$axx)
)

coef.cover %>% {
  ggplot(., aes(axx, perc_cov)) +
      geom_point(aes(
        colour = interaction(model, method) 
      , shape = log_mfi
      )) +
      facet_wrap(~name, nrow = 1) +
    theme(
      axis.text.x = element_blank()
    , axis.ticks.length.x = unit(0, "cm")
    , strip.text.x = element_text(size = 10)
      ) +
    xlab("") +
    ylab("Coverage") +
    facet_wrap(~ok_fit)
}

which_param <- "beta_con1f_delta"
which_model <- "3sd -- no_variance"

which_param <- "beta_cat1f_delta"
which_model <- "mclust -- variance"

which_param <- "beta_cat1f_delta"
which_model <- "publication_model_2.stan"

coeff %>% 
  filter(model == which_model) %>% 
  filter(name == which_param) %>% 
  filter(mid > -10, mid < 10) %>% 
  mutate(cover = as.factor(cover)) %>% {
    ggplot(., aes(true, mid)) + 
      geom_point(aes(colour = cover)) +
      geom_abline(intercept = 0, slope = 1) +
      facet_wrap(~log_mfi+method, scale = "free_y") +
      scale_colour_manual(values = c("firebrick3", "dodgerblue4"))
  }

which_param <- c(
  "beta_base"
, "beta_cat1f_delta"
)

coef.cover.a <- coeff %>% 
  filter(model == which_model) %>% 
  filter(name %in% which_param) %>% 
  filter(mid > -10, mid < 10) %>% 
  mutate(cover = as.factor(cover)) %>%
  dplyr::select(
    log_mfi, param_set, name, m_diff
  ) %>% pivot_wider(
    names_from = name
  , values_from = m_diff
  ) 

coef.cover.b <- coeff %>% 
  filter(model == which_model) %>% 
  filter(name %in% which_param) %>% 
  filter(mid > -10, mid < 10) %>% 
  mutate(cover = as.factor(cover)) %>%
  dplyr::select(
    log_mfi, param_set, name, cover
  ) %>% pivot_wider(
    names_from = name
    , values_from = cover
  ) %>% rename(
    beta_base_c = beta_base
  , beta_cat1f_delta_c = beta_cat1f_delta
  )

coef.cover.a %<>% left_join(
  .
, coef.cover.b
)

coef.cover.a %>% 
  mutate(cover_type = interaction(beta_base_c, beta_cat1f_delta_c)) %>% {
    ggplot(., aes(beta_base, beta_cat1f_delta)) + 
      geom_point(aes(colour = cover_type)) +
      facet_wrap(~log_mfi, scale = "free_y") +
      scale_colour_brewer(palette = "Dark2") +
   #   scale_x_log10() +
   #   scale_y_log10() +
      scale_y_continuous(trans = "pseudo_log") +
      scale_x_continuous(trans = "pseudo_log") +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      xlab("Beta Base: Estimate - True") +
      ylab("Cat1: Estimate - True")
  }

## n_samps
## beta_base
## mu_neg
## sd_neg
## mu_pos
## sd_pos
## mu_pos_delta
## sd_pos_delta

which_param <- "beta_base"

coeff %>% 
  filter(model == which_model) %>% 
  filter(name == which_param) %>% 
  filter(mid > -10, mid < 10) %>% {
    ggplot(., aes(sd_neg, m_diff)) + 
      geom_point() +
      facet_wrap(~log_mfi) +
      geom_hline(yintercept = 0)
  }

which_param <- "beta_con1f_delta"

coeff %>% 
  filter(model == which_model) %>% 
  filter(name == which_param) %>%
  filter(mid > -10, mid < 10) %>%
  group_by(log_mfi) %>%
  summarize(
    m_cov = mean(cover)
  , m_wid = mean(CI_wid)
  , m_bia = mean(m_diff)
  )

####
## Illustrative single parameter set figures for the publication
####

library(scales)

## Potential parameter sets:
 ## 2, 4, 12, 20, 21, 26, 28, 
 ## 30, 36, 

three_sd.groups %>% 
  filter(param_set == 4, log_mfi == "log_mfi") %>% 
  mutate(group = as.factor(group)) %>% {
    ggplot(., aes(x = mfi)) + 
      geom_density(aes(fill = group, colour = group), alpha = 0.4) +
      scale_fill_manual(
        values = c("dodgerblue3", "firebrick3")
      ) +
      scale_colour_manual(
        values = c("dodgerblue3", "firebrick3")
      ) +
      scale_x_continuous(
        breaks = seq(5, 10, 1)
      , labels = round(exp(seq(5, 10, 1)), 0)#label_math(e^.x)
      ) 
  }

####
## Stan convergence
####

stan.summary$fit_details %>% filter(param_set == 133)
sim.data.summaries %>% filter(param_set == 133) %>% as.data.frame()

# prop_overlap_quant diff_mean mfi_skew_2 prop_squish_2 prop_pos_true

stan.summary$fit_details %>% 
  left_join(., sim.data.summaries.j) %>% 
  mutate(max_Rhat = ifelse(is.infinite(max_Rhat), 3000, max_Rhat)) %>% {
    ggplot(., aes(prop_pos_true, max_Rhat)) +
      geom_jitter(aes(colour = model_name), height = 0.1) +
      scale_colour_brewer(palette = "Dark2") +
      scale_y_log10() +
      facet_wrap(~model_name, nrow = 1)
  }

inf_rhats <- stan.summary$fit_details %>% 
  left_join(., sim.data.summaries.j) %>% 
  arrange(desc(max_Rhat)) %>%
  filter(is.infinite(max_Rhat)) %>%
  filter(model_name == "publication_model_normal_2.stan") %>% 
  filter(chain == 1) %>% pull(param_set)

## No fit
failed_fits <- stan.summary$fit_details %>% 
  left_join(., sim.data.summaries.j) %>% 
  filter(log_mfi == "log_mfi") %>%
  group_by(model_name, param_set, log_mfi) %>%
  summarize(
      prop_overlap_quant = mean(prop_overlap_quant)
    , diff_mean = mean(diff_mean)
    , mfi_skew_2 = mean(mfi_skew_2)
    , prop_squish_2 = mean(prop_squish_2)
    , prop_pos_true = mean(prop_pos_true)
    , m_R = mean(max_Rhat, na.rm = T)
    , m_D = mean(divergent_transitions, na.rm = T)
    , m_T = mean(time_to_fit, na.rm = T)
  ) %>% ungroup() %>%
  mutate(m_R = ifelse(is.infinite(m_R), 1000, m_R)) %>%
  group_by(param_set) %>%
  filter(
    m_R == min(m_R, na.rm = T)
  ) %>% arrange(desc(m_R)) %>% 
  as.data.frame() %>% 
  group_by(param_set) %>%
  slice(1) %>%
  arrange(desc(m_R)) %>% 
  as.data.frame() %>%
  filter(m_R > 5) %>%
  pull(param_set)

inf_rhats
all_fails <- failed_fits[failed_fits %in% inf_rhats]

## mfi_skew_2
## prop_pos_true
## prop_overlap_quant

sim.data.summaries %>% 
  filter(log_mfi == "log_mfi") %>%
  mutate(failed_fit = ifelse(param_set %in% all_fails, 1, 0)) %>% {
    ggplot(., aes(prop_overlap_quant, failed_fit)) + geom_point()
  }

sim.data %>% filter(
  param_set %in% inf_rhats
, log_mfi == "log_mfi"
) %>% group_by(param_set) %>%
  summarize(n_pos = length(which(group == 2))/n()) %>% 
  arrange(desc(n_pos))

stan.summary$fit_details %>% 
  left_join(., sim.data.summaries.j) %>% 
  filter(model_name == "publication_model_skewnormal_2.stan") %>% {
    ggplot(., aes(prop_pos_true, divergent_transitions)) +
      geom_jitter(aes(colour = model_name), height = 1) +
      scale_colour_brewer(palette = "Dark2") +
      scale_y_log10() +
      facet_wrap(~log_mfi, nrow = 1)
  }

stan.summary$fit_details %>% 
  left_join(., sim.data.summaries.j) %>% 
  filter(model_name == "publication_model_normal_2.stan") %>%
  group_by(model_name, param_set, log_mfi) %>%
  summarize(
      prop_overlap_quant = mean(prop_overlap_quant)
    , diff_mean = mean(diff_mean)
    , mfi_skew_2 = mean(mfi_skew_2)
    , prop_squish_2 = mean(prop_squish_2)
    , prop_pos_true = mean(prop_pos_true)
    , m_R = mean(max_Rhat)
    , m_D = mean(divergent_transitions)
    , m_T = mean(time_to_fit)
  ) %>% 
  ungroup() %>% 
  mutate(m_R = ifelse(is.infinite(m_R), 1000, m_R)) %>% {
    ggplot(., aes(mfi_skew_2, m_R)) +
      geom_jitter(aes(colour = model_name), height = 1) +
      scale_colour_brewer(palette = "Dark2") +
      scale_y_log10() +
      facet_wrap(~log_mfi, nrow = 1)
  }

#### ------

fits$stan.summary$fit_details %>% 
  group_by(model_name, param_set) %>%
  summarize(
    m_R = mean(max_Rhat)
  , m_D = mean(divergent_transitions)
  , m_T = mean(time_to_fit)
  ) %>% ungroup() %>%
  filter(model_name %in% c(
    "publication_model_skewnormal_2.stan"
  , "publication_model_skewnormal_N_2.stan"
  )) %>% arrange(desc(m_D)) %>%
  mutate(param_set = factor(param_set, levels = unique(param_set))) %>% {
    ggplot(., aes(param_set, m_D)) + 
      geom_point(aes(colour = model_name)) +
      scale_y_log10()
  }

fits$stan.summary$fit_details %>% 
  filter(log_mfi == "log_mfi") %>%
  group_by(model_name, param_set) %>%
  summarize(
    m_R = mean(max_Rhat)
    , m_D = mean(divergent_transitions)
    , m_T = mean(time_to_fit)
  ) %>% ungroup() %>%
  filter(model_name %in% c(
    "publication_model_skewnormal_2.stan"
    , "publication_model_skewnormal_N_2.stan"
  )) %>% arrange(desc(m_D)) %>%
  mutate(param_set = factor(param_set, levels = unique(param_set))) %>% 
  dplyr::select(model_name, param_set, m_D) %>%
  pivot_wider(values_from = m_D, names_from = model_name) %>% 
  mutate(dif_D = publication_model_skewnormal_N_2.stan - publication_model_skewnormal_2.stan) %>%
  arrange(desc(dif_D))

fits$sim.data %>% 
  filter(param_set == 166) %>% 
  mutate(
    group = as.factor(group)
    , int = interaction(cat1f, group)
  ) %>% {
    ggplot(., aes(x = mfi)) + 
      geom_density(aes(colour = group, fill = group), alpha = 0.2) + 
      scale_colour_brewer(palette = "Dark2") +
      scale_fill_brewer(palette = "Dark2") +
      facet_wrap(param_set~log_mfi, scales = "free", nrow = 3)
  }
