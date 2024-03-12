fits <- readRDS("test_fits_adj.Rds")
names(fits)

fits$all.out$group_assignment <- calculate_group_assignments(
  fits$three_sd.groups
, fits$mculst.groups
, fits$stan.summary$group_pred
, fits$sim.params
) 

fits$all.out$group_assignment <- fits$all.out$group_assignment %>% 
  left_join(., fits$all.out$pop_seropositivity %>% 
              dplyr::select(param_set, true)  %>% 
              group_by(param_set) %>% slice(1)
  )

fits$all.out$pop_seropositivity <- calculate_population_seropositivity(
  fits$three_sd.groups
, fits$mculst.groups
, fits$stan.summary$prop_seropos
, fits$sim.params
)


####
## Individual-level probability of assignment 
####

fits$all.out$group_assignment %>% 
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

fits$all.out$group_assignment %>% 
  filter(model == "3sd") %>% {
    ggplot(., aes(true, prop_pos)) + 
      geom_point() +
      geom_abline(intercept = 0, slope = 1) +
      facet_grid(log_mfi~method)
  }

fits$all.out$group_assignment %>% 
  filter(model == "3sd") %>% {
    ggplot(., aes(true, prop_pos_diff)) + 
      geom_point() +
      geom_hline(yintercept = 0) +
      facet_grid(log_mfi~method)
  }

fits$all.out$group_assignment %>% 
  filter(model == "3sd") %>%
  filter(method == "assigned_group_sample_mfi") %>% 
  dplyr::select(log_mfi, param_set, true, prop_pos_diff) %>% 
  pivot_wider(names_from = log_mfi, values_from = prop_pos_diff) %>% {
    ggplot(., aes(mfi, log_mfi)) + 
      geom_point() +
      geom_abline(intercept = 0, slope = 1)
  }

fits$all.out$group_assignment %>% 
  filter(model == "3sd") %>%
  filter(method == "assigned_group_sample_mfi") %>% 
  dplyr::select(log_mfi, param_set, sd_neg, true, prop_pos) %>% 
  pivot_wider(names_from = log_mfi, values_from = prop_pos) %>% {
    ggplot(., aes(sd_neg, log_mfi)) + 
      geom_point(colour = "dodgerblue3") +
      geom_point(aes(sd_neg, mfi), colour = "firebrick3") +
      geom_errorbar(aes(ymin = mfi, ymax = log_mfi), width = 0, alpha = 0.5)
  }

test_for_gg <- fits$all.out$group_assignment %>% 
  filter(model == "3sd") %>%
  filter(method == "assigned_group_sample_mfi") %>%
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

better_log_option_B <- fits$all.out$group_assignment %>% 
  filter(model == "3sd") %>%
  filter(method == "assigned_group_sample_mfi") %>% 
  dplyr::select(log_mfi, param_set, true, prop_pos) %>% 
  pivot_wider(names_from = log_mfi, values_from = prop_pos) %>%
  mutate(
    log_mfi = abs(log_mfi - true)
  , mfi     = abs(mfi - true)
  ) %>% mutate(diff = log_mfi - mfi) %>%
  arrange(desc(diff)) %>%
  left_join(
    .
  , fits$all.out$group_assignment %>% 
    filter(model == "3sd") %>%
    filter(method == "assigned_group_sample_mfi") %>%
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

fits$all.out$group_assignment %>% 
  filter(model == "3sd") %>%
  filter(method == "assigned_group_sample_mfi") %>%
  mutate(
    gg = ifelse(
      prop_pos_diff > -0.001 & prop_pos_diff < 0.001
      , 1, 0)
  ) %>% {
     ggplot(., aes(sd_neg, gg)) + 
        geom_jitter() +
        facet_wrap(~log_mfi, nrow = 1)
  }

fits$all.out$group_assignment %>% 
  filter(model == "3sd") %>%
  filter(method == "assigned_group_sample_mfi") %>% 
  dplyr::select(log_mfi, n_samps, beta_base, mu_neg, sd_neg
              , mu_pos, sd_pos, mu_pos_delta, sd_pos_delta
              , prop_pos_diff, prop_pos) %>% {
    ggplot(., aes(sd_neg, prop_pos)) + 
                  geom_point() +
                  facet_wrap(~log_mfi, nrow = 1)
  } 

fits$three_sd.groups %>% 
  filter(param_set == 1) %>% {
    ggplot(., aes(group, assigned_group)) + 
      geom_jitter(width = 0.2, height = 0.2) +
      scale_x_continuous(breaks = c(0, 1)) +
      scale_y_continuous(breaks = c(0, 1)) +
      facet_wrap(~sd_method+log_mfi, scales = "free_x") +
      xlab("True Group") +
      ylab("Assigned Group")
  }

fits$three_sd.groups %>% 
 filter(param_set == 1) %>% 
  mutate(group = as.factor(group)) %>% {
    ggplot(., aes(mfi, assigned_group)) + 
      geom_jitter(
        aes(colour = group)  
      , width = 0.2, height = 0.2) +
      scale_y_continuous(breaks = c(0, 1)) +
      facet_wrap(~sd_method+log_mfi, scales = "free_x") +
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

fits$all.out$group_assignment %>% 
  filter(model == "mclust") %>%
  filter(method == "constrained_mclust") %>% {
    ggplot(., aes(true, prop_pos)) + geom_point() +
      geom_abline(intercept = 0, slope = 1)
  }

fits$all.out$group_assignment %>% 
  filter(model == "mclust") %>%
  filter(method == "constrained_mclust") %>% {
    ggplot(., aes(true, prop_pos)) + 
      geom_point() +
      geom_abline(intercept = 0, slope = 1) +
      facet_grid(log_mfi~method)
  }

fits$all.out$group_assignment %>% 
  filter(model == "mclust") %>%
  filter(method == "constrained_mclust") %>% {
    ggplot(., aes(true, prop_pos_diff)) + geom_point() +
      geom_hline(yintercept = 0)
  }

fits$all.out$group_assignment %>% 
  filter(model == "mclust") %>%
  filter(method == "constrained_mclust") %>% 
  dplyr::select(log_mfi, param_set, true, prop_pos_diff) %>% 
  pivot_wider(names_from = log_mfi, values_from = prop_pos_diff) %>% {
    ggplot(., aes(mfi, log_mfi)) + 
      geom_point() +
      geom_abline(intercept = 0, slope = 1)
  }

fits$all.out$group_assignment %>% 
  filter(model == "mclust") %>%
  filter(method == "constrained_mclust") %>% 
  dplyr::select(log_mfi, param_set, sd_neg, true, prop_pos) %>% 
  pivot_wider(names_from = log_mfi, values_from = prop_pos) %>% {
    ggplot(., aes(sd_neg, log_mfi)) + 
      geom_point(colour = "dodgerblue3") +
      geom_point(aes(sd_neg, mfi), colour = "firebrick3") +
      geom_errorbar(aes(ymin = mfi, ymax = log_mfi), width = 0, alpha = 0.5)
  }

fits$all.out$group_assignment %>% 
  filter(model == "mclust") %>%
  filter(method == "unconstrained_mclust") %>% {
    ggplot(., aes(true, prop_pos)) + geom_point() +
      geom_abline(intercept = 0, slope = 1)
  }

fits$all.out$group_assignment %>% 
  filter(model == "mclust") %>%
  filter(method == "unconstrained_mclust") %>% {
    ggplot(., aes(true, prop_pos_diff)) + 
      geom_point() +
      geom_hline(yintercept = 0)
  }

fits$all.out$group_assignment %>% 
  filter(model == "publication_model_2.stan") %>%
  filter(method == "Bayesian LCR") %>% 
  filter(quantile == "mid") %>% {
    ggplot(., aes(true, prop_pos)) + geom_point() +
      geom_abline(intercept = 0, slope = 1) +
      facet_wrap(~log_mfi, ncol = 1)
  }

fits$all.out$group_assignment %>% 
  filter(model == "publication_model_2.stan") %>%
  filter(method == "Bayesian LCR") %>% 
  filter(quantile == "mid") %>% {
    ggplot(., aes(true, prop_pos_diff)) + geom_point() +
      geom_hline(yintercept = 0)
  }



####
## Population-level seropositivity 
####

fits$all.out$pop_seropositivity %>% 
  filter(model == "3sd") %>% {
    ggplot(., aes(true, prop_pos)) + 
      geom_point() +
      geom_abline(intercept = 0, slope = 1) +
      facet_grid(log_mfi~method)
  }

fits$all.out$pop_seropositivity %>% 
  filter(model == "mclust") %>% {
    ggplot(., aes(true, prop_pos)) + 
      geom_point() +
      geom_abline(intercept = 0, slope = 1) +
      facet_grid(log_mfi~method)
  }

fits$all.out$pop_seropositivity %>% 
  filter(model == "publication_model_2.stan") %>% {
    ggplot(., aes(true, prop_pos)) + 
      geom_point() +
      geom_abline(intercept = 0, slope = 1) +
      facet_grid(log_mfi~method)
  }


####
## Coefficient estimates from regression model
####

coeff <- fits$all.out$coefficient_ests %>% 
  mutate(m_diff = ifelse(mid < true, m_diff * -1, m_diff)) %>%
  left_join(
  .
, fits$sim.params
)

# beta_base
# beta_cat1f_delta
# beta_cat2f_delta
# beta_con1f_delta

coef.cover <- coeff %>% 
  mutate(cover = ifelse(is.na(cover), 0, cover)) %>%
  group_by(model, method, log_mfi, name) %>%
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
    ylab("Coverage")
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

fits$three_sd.groups %>% 
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
      ) +
      geom_point(
        aes()
      )
  }

plot_individual_group_prob(
    three_sd.g = fits$three_sd.groups
  , mclust.g   = fits$mculst.groups
  , stan.g     = fits$stan.summary$group_pred
  , which_fits = c(4, 26, 36)
)

plot_group_assignment_summary(
  group_assignment   = fits$group_assignment
)

