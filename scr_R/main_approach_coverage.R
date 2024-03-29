all_coef_ests <- all.out$coefficient_ests %>% 
  mutate(cover = ifelse(is.na(cover), 0, cover)) %>%
  mutate(mod.meth = interaction(model, method)) %>%
  mutate(
    main_approach = plyr::mapvalues(
      method
      , from = unique(method)
      , to   = c("Bayesian LCR", "Mclust", "Mclust", "Mclust", "3sd", "3sd")
    )
    , .before = model
  ) %>%
  mutate(
    sub_approach = plyr::mapvalues(
      mod.meth
      , from = unique(mod.meth)
      , to = c(
          "normal-normal"
        , "normal-skewnormal"
        , "lognormal-lognormal"
        , "mclust (2 group constrained) + unweighted regression"
        , "mclust (2 group constrained) + probability weighted regression"
        , "mclust (Unconstrained) + unweighted regression"
        , "mclust (Unconstrained) + probability weighted regression"
        , "mclust (Unconstrained reduced) + unweighted regression"
        , "mclust (Unconstrained reduced) + probability weighted regression"
        , "3sd (Approximated control)"
        , "3sd (Robust Mean and SD)"
      ))
  )

param_coverage_by_model <- all_coef_ests %>%
  filter(name %in% c(
    "beta_base", "beta_cat1f_delta", "beta_cat2f_delta", "beta_con1f_delta"
  )) %>% group_by(
    param_set, log_mfi, main_approach, sub_approach
  ) %>% 
  summarize(
    cover  = sum(cover, na.rm = T)
  ) %>% ungroup() %>% 
  group_by(param_set, main_approach) %>% 
  filter(cover == max(cover))

ci_wid_by_model <- all_coef_ests %>% 
  filter(name %in% c(
    "beta_base", "beta_cat1f_delta", "beta_cat2f_delta", "beta_con1f_delta"
  )) %>% 
  dplyr::select(main_approach, sub_approach, log_mfi, name, param_set, CI_wid) %>%
  group_by(main_approach, sub_approach, log_mfi, param_set) %>%
  summarize(ci_wid_mean = mean(CI_wid, na.rm = T)) %>%
  ungroup()

the_best_fits <- param_coverage_by_model %>% 
  left_join(., ci_wid_by_model) %>%
  filter(ci_wid_mean == min(ci_wid_mean)) %>%
  ungroup() %>%
  mutate(this_fit = 1)

best_fits <- all_coef_ests %>% 
  left_join(., the_best_fits %>% rename(tot_cover = cover)) %>%
  filter(this_fit == 1)

best_fits %<>% 
  mutate(m_diff = ifelse(mid < true, m_diff * -1, m_diff)) %>%
  left_join(
    .
    , sim.params
  ) %>% left_join(
    .
    , sim.data.summaries.j
  ) %>% left_join(
    .
    , stan.fit_stats %>% 
      dplyr::select(model_name, param_set, log_mfi, m_R, m_D, m_T) %>%
      rename(model = model_name)
  ) %>% mutate(
    ok_fit = ifelse((is.na(m_R) | m_R < 1.10), 1, 0)
  )

coef.cover <- best_fits %>% 
  mutate(cover = ifelse(is.na(cover), 0, cover)) %>%
  group_by(main_approach, name) %>%
  filter(name %notin% c(
      "mu_neg", "mu_pos", "mu_pos_delta"
    , "sd_neg", "sd_pos", "theta_con2f_delta"
  )) %>% 
  summarize(
      perc_cov = sum(cover) / n()
    , m_bias   = mean(m_diff, na.rm = T)
  ) %>% 
  ungroup() %>%
  arrange(perc_cov) %>% 
  mutate(
    main_approach = factor(main_approach, levels = unique(main_approach))
  )

coef.cover %>% mutate(name = plyr::mapvalues(
  name
  , from = c("beta_base", "beta_cat1f_delta", "beta_cat2f_delta", "beta_con1f_delta")
  , to   = c("Baseline (Intercept)", "Categorical Covariate 1"
             , "Categorical Covariate 2", "Continuous Covariate 1")
)) %>% {
    ggplot(., aes(main_approach, perc_cov)) +
      geom_point(aes(
          colour = main_approach
        , shape = name
      ), size = 5) +
      scale_shape_manual(name = "Regression
Parameter", values = c(16, 6, 2, 12)) +
      scale_colour_manual(values = c(
          "#1b9e77"
        , "#7570b3"
        , "#e6ab02"
      ), name = "Approach") +
      theme(
          axis.text.x = element_blank()
        , axis.ticks.length.x = unit(0, "cm")
      ) +
      xlab("") +
      ylab("Coverage")
  }

best_fits %>% 
  group_by(main_approach, param_set, name) %>% 
  slice(1) %>%
  ungroup() %>%
  filter(name %notin% c(
    "mu_neg", "mu_pos", "mu_pos_delta"
    , "sd_neg", "sd_pos", "theta_con2f_delta"
  )) %>%
  mutate(name = plyr::mapvalues(
  name
  , from = c("beta_base", "beta_cat1f_delta", "beta_cat2f_delta", "beta_con1f_delta")
  , to   = c("Baseline (Intercept)", "Categorical Covariate 1"
             , "Categorical Covariate 2", "Continuous Covariate 1")
)) %>% mutate(main_approach = factor(
  main_approach
  , levels = c("3sd", "Mclust", "Bayesian LCR"))
  ) %>% {
  ggplot(., aes(main_approach, m_diff)) +
    geom_violin(aes(
        colour = main_approach, fill = main_approach
    )) +
    scale_shape_manual(name = "Regression
Parameter", values = c(16, 6, 2, 12)) +
    scale_colour_manual(values = c(
        "#1b9e77"
      , "#7570b3"
      , "#e6ab02"
    ), name = "Approach") +
    scale_fill_manual(values = c(
        "#1b9e77"
      , "#7570b3"
      , "#e6ab02"
    ), name = "Approach") +
    theme(
        axis.text.x = element_blank()
      , axis.ticks.length.x = unit(0, "cm")
      , strip.text.x = element_text(size = 12)
    ) +
    xlab("") +
    ylab("Bias") +
    facet_wrap(~name) +
    scale_y_continuous(
      trans = "pseudo_log"
    , breaks = c(50, 20, 5, 0, -5, -20, -50, -100, -250)) +
    geom_hline(yintercept = 0, linetype = "dashed")
}

####################################################################################

all_coef_ests <- all.out$coefficient_ests %>% 
  mutate(m_diff = ifelse(mid < true, m_diff * -1, m_diff)) %>%
  mutate(cover = ifelse(is.na(cover), 0, cover)) %>%
  mutate(mod.meth = interaction(model, method)) %>%
  filter(name %notin% c(
    "mu_neg", "mu_pos", "mu_pos_delta"
    , "sd_neg", "sd_pos", "theta_con2f_delta"
  )) %>%
  mutate(
    main_approach = plyr::mapvalues(
      method
      , from = unique(method)
      , to   = c("Bayesian LCR", "Mclust", "Mclust", "Mclust", "3sd", "3sd")
    )
    , .before = model
  ) %>%
  mutate(
    sub_approach = plyr::mapvalues(
      mod.meth
      , from = unique(mod.meth)
      , to = c(
          "normal-normal"
        , "normal-skewnormal"
        , "lognormal-lognormal"
        , "(2 group constrained) + unweighted regression"
        , "(2 group constrained) + probability weighted regression"
        , "(Unconstrained) + unweighted regression"
        , "(Unconstrained) + probability weighted regression"
        , "(Unconstrained reduced) + unweighted regression"
        , "(Unconstrained reduced) + probability weighted regression"
        , "(Approximated control)"
        , "(Robust Mean and SD)"
      ))
  )

coef.cover <- all_coef_ests %>% 
  mutate(cover = ifelse(is.na(cover), 0, cover)) %>%
  group_by(main_approach, sub_approach, name) %>%
  summarize(
      perc_cov = sum(cover) / n()
    , m_bias   = mean(m_diff, na.rm = T)
  ) %>% 
  ungroup() %>%
  arrange(perc_cov) %>% 
  mutate(m.s = interaction(main_approach, sub_approach, sep = " -- ")) %>%
  mutate(m.s = factor(m.s, levels = unique(m.s)))

coef.cover %>%  
  mutate(name = plyr::mapvalues(
    name
    , from = c("beta_base", "beta_cat1f_delta", "beta_cat2f_delta", "beta_con1f_delta")
    , to   = c("Baseline (Intercept)", "Categorical Covariate 1"
               , "Categorical Covariate 2", "Continuous Covariate 1")
  )) %>% {
    ggplot(., aes(m.s, perc_cov)) +
      geom_point(aes(colour = m.s, shape = name), size = 3) +
      scale_colour_discrete(name = "Approach") +
      scale_shape_manual(name = "Regression
Parameter", values = c(16, 6, 2, 12)) +
      theme(
        axis.text.x = element_blank()
      , axis.ticks.length.x = unit(0, "cm")
      , strip.text.x = element_text(size = 10)
      ) +
      xlab("") +
      ylab("Coverage") +
      scale_y_continuous(lim = c(0, 1)) #+
     # scale_y_continuous(breaks = c(2, 1, 0, -2, -4, -8)) #+
    #  geom_hline(yintercept = 0, linetype = "dashed")
  }
