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
        , "mclust (Unconstrained; sum of 2-n)"
        , "mclust (Unconstrained; BIC collapsed)"
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

gg.2 <- all_coef_ests %>% mutate(m.s = factor(m.s, levels = coef.cover$m.s)) %>% 
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

all_coef_ests %<>% 
  dplyr::select(-prop_pos_diff) %>%
  pivot_wider(., values_from = prop_pos, names_from = quantile) %>%
  mutate(cover = ifelse(lwr < true & upr > true, 1, 0)) %>%
  mutate(
    prop_pos_diff = mid - true
  , CI_wid = upr - lwr
  )
  
param_coverage_by_model <- all_coef_ests %>%
  group_by(param_set, log_mfi, main_approach, sub_approach) %>% 
  summarize(cover = sum(cover, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(param_set, main_approach) %>% 
  filter(cover == max(cover))
  
ci_wid_by_model <- all_coef_ests %>% 
  dplyr::select(main_approach, sub_approach, log_mfi, param_set, CI_wid) %>%
  group_by(main_approach, sub_approach, log_mfi, param_set) %>%
  summarize(ci_wid_mean = mean(CI_wid, na.rm = T)) %>%
  ungroup()

param_coverage_by_model.core <- param_coverage_by_model %>%
  group_by(main_approach, param_set) %>%
  summarize(cover = sum(cover)) %>%
  ungroup() %>%
  group_by(main_approach) %>%
  summarize(cover = length(which(cover > 0)) / n())

param_coverage_by_model.all <- param_coverage_by_model %>%
  group_by(main_approach, sub_approach, param_set) %>%
  summarize(cover = sum(cover)) %>%
  ungroup() %>%
  group_by(main_approach, sub_approach) %>%
  summarize(cover = length(which(cover > 0)) / n())

param_coverage_by_model.all.adj <- param_coverage_by_model.all %>% 
  mutate(sub_approach = as.character(sub_approach)) %>%
  mutate(m.s = interaction(main_approach, sub_approach, sep = " -- ")) %>%
  mutate(m.s = as.character(m.s)) %>%
  arrange(cover) 
  
levvs <- param_coverage_by_model.all.adj %>% pull(m.s)
  
gg.1 <- param_coverage_by_model.all.adj %>%
  mutate(m.s = factor(m.s, levels = levvs)) %>%
  rename(Approach = m.s) %>% {
    ggplot(., aes(Approach, cover)) +
      geom_point(aes(colour = Approach), size = 3) +
     # geom_point(data = param_coverage_by_model.core, aes(main_approach, cover)) + 
      theme(
          axis.text.x = element_blank()
        , axis.ticks.length.x = unit(0, "cm")
        , strip.text.x = element_text(size = 12)
      ) +
      xlab("") +
      ylab("Coverage")
  }

gridExtra::grid.arrange(gg.1, gg.2, ncol = 1)
