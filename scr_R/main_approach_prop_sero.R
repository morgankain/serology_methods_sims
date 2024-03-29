all_coef_ests <- all.out$pop_seropositivity %>% 
  mutate(mod.meth = interaction(model, method)) %>%
  mutate(
    main_approach = plyr::mapvalues(
      method
      , from = unique(method)
      , to   = c("3sd", "3sd", "Mclust", "Mclust", "Mclust", "Bayesian LCR")
    )
    , .before = model
  ) %>%
  mutate(
    sub_approach = plyr::mapvalues(
      mod.meth
      , from = unique(mod.meth)
      , to = c(
          "3sd (Approximated control)"
        , "3sd (Robust Mean and SD)"
        , "mclust (2 group constrained)"
        , "mclust (Unconstrained)"
        , "mclust (Unconstrained reduced)"
        , "normal-normal"
        , "normal-skewnormal"
        , "lognormal-lognormal"
      ))
  ) %>% 
  mutate(m.s = interaction(main_approach, sub_approach, sep = " -- ")) %>%
  mutate(m.s = factor(m.s, levels = unique(m.s)))

coef.cover <- all_coef_ests %>% 
  group_by(main_approach, sub_approach) %>%
  summarize(m_bias = mean(prop_pos_diff, na.rm = T)) %>% 
  ungroup() %>%
  arrange(m_bias) %>% 
  mutate(m.s = interaction(main_approach, sub_approach, sep = " -- ")) %>%
  mutate(m.s = factor(m.s, levels = unique(m.s)))

all_coef_ests %>% mutate(m.s = factor(m.s, levels = coef.cover$m.s)) %>% 
  rename(Approach = m.s) %>% {
    ggplot(., aes(Approach, prop_pos_diff)) +
      geom_violin(aes(colour = Approach, fill = Approach), alpha = 0.5) +
      geom_point(data = coef.cover %>% rename(Approach = m.s, prop_pos_diff = m_bias), aes(colour = Approach), size = 3) +
      theme(
          axis.text.x = element_blank()
        , axis.ticks.length.x = unit(0, "cm")
        , strip.text.x = element_text(size = 12)
      ) +
      xlab("") +
      ylab("Bias") +
      geom_hline(yintercept = 0, linetype = "dashed")
  }

all_coef_ests %>% 
  filter(main_approach == "Bayesian LCR", quantile == "cover") %>% {
    ggplot(., aes(mfi_skew_2, prop_pos)) + 
      geom_jitter(aes(colour = sub_approach))
  }

all_coef_ests %>% 
  filter(main_approach == "Bayesian LCR") %>% 
  left_join(., the_best_fits %>% filter(main_approach == "Bayesian LCR")) %>%
  filter(this_fit == 1, quantile == "cover") %>%
  group_by(main_approach) %>%
  summarize(perc_cov = sum(prop_pos) / n()) %>% 
  ungroup() %>%
  arrange(perc_cov) %>% 
  mutate(
    main_approach = factor(main_approach, levels = unique(main_approach))
  )


