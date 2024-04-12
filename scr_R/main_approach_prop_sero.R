failed_fits <- stan.summary$fit_details %>% 
  left_join(., sim.data.summaries.j) %>%
  group_by(model_name, param_set, log_mfi) %>%
  summarize(
      prop_overlap_quant = mean(prop_overlap_quant)
    , diff_mean = mean(diff_mean)
    , mfi_skew_2 = mean(mfi_skew_2)
    , prop_squish_2 = mean(prop_squish_2)
    , prop_pos_true = mean(prop_pos_true)
   # , n_samps = mean(n_samps)
    , m_R = mean(max_Rhat, na.rm = T)
    , m_D = mean(divergent_transitions, na.rm = T)
    , m_T = mean(time_to_fit, na.rm = T)
  ) %>% ungroup() %>%
  mutate(m_R = ifelse(is.infinite(m_R), 1000, m_R)) %>%
  filter(m_R > 5) %>% 
  mutate(failed = 1)

all_coef_ests <- all.out$pop_seropositivity %>% 
  mutate(prop_pos_diff = prop_pos/true) %>%
  left_join(., failed_fits %>% rename(model = model_name) %>% dplyr::select(model, param_set, log_mfi, failed)) %>% 
  mutate(failed = ifelse(is.na(failed), 0, 1)) %>%
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
  mutate(m.s = factor(m.s, levels = unique(m.s))) #%>%
  #filter(failed == 0)

coef.cover <- all_coef_ests %>% 
  group_by(main_approach, sub_approach) %>%
  summarize(m_bias = mean(prop_pos_diff, na.rm = T)) %>% 
  ungroup() %>%
  arrange(m_bias) %>% 
  mutate(m.s = interaction(main_approach, sub_approach, sep = " -- ")) %>%
  mutate(m.s = factor(m.s, levels = unique(m.s)))

all_coefs_for_gg.2 <- all_coef_ests 

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
        , plot.margin = unit(c(0.2,0.2,0.2,0.8), "cm")
      ) +
      xlab("") +
      ylab("Coverage")
  }

gg.2 <- all_coefs_for_gg.2 %>% 
  mutate(m.s = factor(m.s, levels = levvs)) %>% 
  rename(Approach = m.s) %>% {
    ggplot(., aes(Approach, prop_pos_diff)) +
      geom_violin(aes(colour = Approach, fill = Approach), alpha = 0.5) +
      geom_point(data = coef.cover %>% rename(Approach = m.s, prop_pos_diff = m_bias), aes(colour = Approach), size = 3) +
      theme(
          axis.text.x = element_blank()
        , axis.ticks.length.x = unit(0, "cm")
        , strip.text.x = element_text(size = 12)
        , plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
      ) +
      scale_y_log10(
        breaks = c(.001, .1, 1, 10, 50)
      , labels = c("0.001", "0.01", "1.00", "10", "50")) +
      xlab("") +
      ylab("Estimated / True Seropositivity") +
      geom_hline(yintercept = 1, linetype = "dashed")
  }

gridExtra::grid.arrange(gg.1, gg.2, ncol = 1)

param_coverage_by_model.all.adj %<>% 
  ungroup() %>% mutate(
  m.s = plyr::mapvalues(m.s, from = unique(m.s), to = c(
    "3sd (Robust Mean and SD)"                
  , "3sd (Approximated control)"              
  , "Mclust (Unconstrained; sum of 2-n)"   
  , "Mclust (Unconstrained; BIC collapsed)"
  , "Mclust (2 group constrained)"         
  , "Bayesian LCR (lognormal-lognormal)"            
  , "Bayesian LCR (normal-normal)"                  
  , "Bayesian LCR (normal-skewnormal)"   
  ))
) 

levvs <- param_coverage_by_model.all.adj$m.s

param_coverage_by_model.all.adj %<>% mutate(
  m.s = factor(m.s, levels = levvs)
)

all_coefs_for_gg.2 %<>% ungroup() %>% mutate(
  m.s = plyr::mapvalues(m.s, from = unique(m.s), to = c(
      "3sd (Approximated control)"
    , "3sd (Robust Mean and SD)"  
    , "Mclust (2 group constrained)"    
    , "Mclust (Unconstrained; sum of 2-n)"   
    , "Mclust (Unconstrained; BIC collapsed)"
    , "Bayesian LCR (normal-normal)" 
    , "Bayesian LCR (normal-skewnormal)" 
    , "Bayesian LCR (lognormal-lognormal)"            
  ))
) %>% mutate(
  m.s = factor(m.s, levels = unique(param_coverage_by_model.all.adj$m.s))
)

gg.1 <- param_coverage_by_model.all.adj %>%
  rename(Approach = m.s) %>% {
    ggplot(., aes(cover, Approach)) +
      geom_point(aes(colour = main_approach), size = 3) +
      scale_colour_manual(values = c(
          "#1b9e77"
        , "#7570b3"
        , "#e6ab02"
      ), name = "Approach") +
      # geom_point(data = param_coverage_by_model.core, aes(main_approach, cover)) + 
      theme(
          strip.text.x = element_text(size = 12)
        , plot.margin = unit(c(0.2,0.05,0.2,0.8), "cm")
        , axis.text.y = element_text(size = 10)
        , axis.text.x = element_text(size = 12)
        , axis.ticks.length.y = unit(0, "cm")
        , legend.position = "none"
      ) +
      scale_x_continuous(
        breaks = c(0.1, 0.35, 0.60, 0.85, 0.95)
      , labels = c("10%", "35%", "60%", "85%", "95%")) +
      xlab("Coverage") +
      ylab("Model") +
      geom_vline(xintercept = 0.95, linetype = "dashed")
  }

gg.2 <- all_coefs_for_gg.2 %>% 
  filter(quantile == "mid") %>%
  mutate(m.s = factor(m.s, levels = levvs)) %>% 
  rename(Approach = m.s) %>% {
    ggplot(., aes(prop_pos_diff, Approach)) +
      geom_boxplot(aes(fill = main_approach), alpha = 0.4) +
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
          strip.text.x = element_text(size = 11)
        , plot.margin = unit(c(0.2,0.2,0.2,0.05), "cm")
        , axis.text.y = element_blank()
        , axis.text.x = element_text(size = 12)
      ) +
      scale_x_log10(
          breaks = c(.001, 0.01, 0.1, 1, 10, 50)
        , labels = c("0.001", "0.01", "0.1", "1.00", "10", "50")) +
      ylab("") +
      xlab("Bias") +
      geom_vline(xintercept = 1, linetype = "dashed")
  }

gridExtra::grid.arrange(gg.1, gg.2, nrow = 1)
