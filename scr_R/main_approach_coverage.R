stan.fit_stats <- stan.summary$fit_details %>% 
  left_join(., sim.data.summaries.j) %>% 
  group_by(model_name, param_set, log_mfi) %>%
  summarize(
      prop_overlap_quant = mean(prop_overlap_quant)
    , diff_mean = mean(diff_mean)
    , mfi_skew_2 = mean(mfi_skew_2)
    , prop_squish_2 = mean(prop_squish_2)
    , prop_pos_true = mean(prop_pos_true)
#    , n_samps = mean(n_samps)
    , m_R = mean(max_Rhat)
    , m_D = mean(divergent_transitions)
    , m_T = mean(time_to_fit)
  )

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
  ungroup()# %>% 
 # left_join(., sim.data.summaries.j %>% dplyr::select(param_set, n_samps) %>% group_by(param_set) %>% slice(1) %>% ungroup())

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

gg.1 <- coef.cover %>% 
  mutate(main_approach = plyr::mapvalues(main_approach, from = "Mclust", to = "GMM")) %>%
  mutate(main_approach = factor(
    main_approach
    , levels = c("3sd", "Bayesian LCR", "GMM"))
  ) %>%
  mutate(name = plyr::mapvalues(
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
      ylab("Coverage") +
    geom_hline(yintercept = 0.95, linetype = "dashed")
  }

gg.2 <- best_fits %>% 
  mutate(main_approach = plyr::mapvalues(main_approach, from = "Mclust", to = "GMM")) %>%
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
  , levels = c("3sd", "Bayesian LCR", "GMM"))
  ) %>%
  filter(
    m_diff > -20
  ) %>% {
  ggplot(., aes(main_approach, m_diff)) +
    #geom_violin(aes(colour = main_approach, fill = main_approach)) +
    geom_boxplot(aes(fill = main_approach), alpha = 0.4, linewidth = 0.3) +
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
      , axis.text.y = element_text(size = 10)
    ) +
    xlab("") +
    ylab("Bias") +
    facet_wrap(~name) +
    scale_y_continuous(
      trans = "pseudo_log"
   # , breaks = c(50, 20, 5, 0, -5, -20, -50, -100, -250)
     , breaks = c(50, 20, 5, 0, -5, -20)
   ) +
    geom_hline(yintercept = 0, linetype = "dashed")
  }

gridExtra::grid.arrange(gg.1, gg.2, ncol = 1)

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
        , "(Unconstrained; sum of 2-n) + unweighted regression"
        , "(Unconstrained; sum of 2-n) + probability weighted regression"
        , "(Unconstrained; BIC collapsed) + unweighted regression"
        , "(Unconstrained; BIC collapsed) + probability weighted regression"
        , "(Approximated control)"
        , "(Robust Mean and SD)"
      ))
  )

coef.cover <- all_coef_ests %>% 
  left_join(., best_fits %>% dplyr::select(main_approach, sub_approach, log_mfi, param_set, ok_fit, name)) %>%
  mutate(ok_fit = ifelse(is.na(ok_fit), 1, 0)) %>%
 # filter(ok_fit == 1) %>%
  mutate(cover = ifelse(is.na(cover), 0, cover)) %>%
  group_by(main_approach, sub_approach, name, param_set) %>%
  summarize(cover = sum(cover)) %>%
  mutate(cover = ifelse(cover > 0, 1, 0)) %>%
  ungroup(param_set) %>%
  summarize(perc_cov = sum(cover) / n()) %>% 
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
  )) %>%
  mutate(main_approach = plyr::mapvalues(main_approach, from = "Mclust", to = "GMM")) %>%
  mutate(main_approach = factor(
    main_approach
    , levels = c("3sd", "Bayesian LCR", "GMM"))
  ) %>% mutate(
    
    m.s  = plyr::mapvalues(
      m.s
      , from = c(
          "Mclust -- (Unconstrained; sum of 2-n) + unweighted regression"
        , "Mclust -- (Unconstrained; sum of 2-n) + probability weighted regression"
        , "Mclust -- (2 group constrained) + unweighted regression"
        , "Mclust -- (2 group constrained) + probability weighted regression"
        , "Mclust -- (Unconstrained; BIC collapsed) + unweighted regression"
        , "Mclust -- (Unconstrained; BIC collapsed) + probability weighted regression"
      ) 
      , to   = c(
          "GMM -- (Unconstrained; sum of 2-n) + unweighted regression"
        , "GMM -- (Unconstrained; sum of 2-n) + probability weighted regression"
        , "GMM -- (2 group constrained) + unweighted regression"
        , "GMM -- (2 group constrained) + probability weighted regression"
        , "GMM -- (Unconstrained; BIC collapsed) + unweighted regression"
        , "GMM -- (Unconstrained; BIC collapsed) + probability weighted regression"
      )
    )
    
  ) %>% {
    ggplot(., aes(perc_cov, m.s)) +
      geom_point(aes(colour = main_approach, shape = name), size = 3) +
      scale_colour_manual(values = c(
          "#1b9e77"
        , "#7570b3"
        , "#e6ab02"
      ), name = "Approach") +
      
      scale_shape_manual(name = "Regression
Parameter", values = c(16, 6, 2, 12)) +
      theme(
          strip.text.x = element_text(size = 11)
        , plot.margin = unit(c(0.2,0.2,0.2,0.05), "cm")
        , axis.text.x = element_text(size = 12)
        , axis.text.y = element_text(size = 12)
      ) +
      ylab("Model") +
      xlab("Coverage") +
      scale_x_continuous(lim = c(0, 1)) +
      geom_vline(xintercept = 0.95, linetype = "dashed")
    # scale_y_continuous(breaks = c(2, 1, 0, -2, -4, -8)) #+
    #  geom_hline(yintercept = 0, linetype = "dashed")
  }

