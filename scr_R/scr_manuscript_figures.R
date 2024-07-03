#####
## Bit of extra cleanup
#####

tar_load(all.out); tar_load(sim.data.summaries)
tar_load(stan.summary); tar_load(mclust.groups)
tar_load(three_sd.groups); tar_load(sim.params)
tar_load(sim.data)

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


##### 
## Figure 1: 3sd approaches, log and linear
#####

all.out$pop_seropositivity %>% 
  mutate(prop_pos_diff = prop_pos / true) %>%
  mutate(prop_pos_diff = ifelse(prop_pos_diff == 0, 0.0001, prop_pos_diff)) %>% 
  filter(model == "3sd") %>%
  filter(quantile == "mid") %>%
  dplyr::select(log_mfi, param_set, method, true, prop_pos_diff) %>% 
  pivot_wider(names_from = log_mfi, values_from = prop_pos_diff) %>% 
  mutate(method = plyr::mapvalues(method, from = c(
    "assigned_group_control", "assigned_group_robust"
  ), to = c(
    "Approximated
Control", "Robust Mean
and SD"
  ))) %>% {
    ggplot(., aes(mfi, log_mfi)) + 
      geom_point(aes(colour = method)) +
      geom_point(
        data = all.out$pop_seropositivity %>% 
          mutate(prop_pos_diff = prop_pos / true) %>%
          mutate(prop_pos_diff = ifelse(prop_pos_diff == 0, 0.0001, prop_pos_diff)) %>% 
          filter(model == "3sd") %>%
          filter(quantile == "mid") %>%
          dplyr::select(log_mfi, param_set, method, true, prop_pos_diff) %>% 
          pivot_wider(names_from = log_mfi, values_from = prop_pos_diff) %>%
          mutate(method = plyr::mapvalues(method, from = c(
            "assigned_group_control", "assigned_group_robust"
          ), to = c(
            "Approximated
Control", "Robust Mean
and SD"
          ))) %>%
          filter(param_set %in% c(445, 412))
        , aes(mfi, log_mfi, colour = method), size = 12, shape = c(15, 15, 17, 17)
        , alpha = 0.3
      ) +
      scale_y_log10() +
      scale_x_log10() +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      geom_vline(xintercept = 1, linetype = "dotted") +
      geom_hline(yintercept = 1, linetype = "dotted") +
      xlab("Estimated / True Seropositivity (MFI)") +
      ylab("Estimated / True Seropositivity (Log MFI)") +
      scale_colour_brewer(palette = "Dark2", name = "Method") +
      theme(
        legend.key.size = unit(0.8, "cm")
      , legend.position = c(0.2, 0.85)) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  }

## 445, 412, 304, 345

three_sd.groups %>% 
  filter(param_set %in% c(445, 412)) %>% 
  mutate(param_set = as.factor(param_set)) %>%
  mutate(sd_method = plyr::mapvalues(sd_method, from = c(
    "assigned_group_control", "assigned_group_robust"
  ), to = c(
    "Approximated
Control", "Robust Mean
and SD"
  ))) %>%
  mutate(group = plyr::mapvalues(group, from = c(0, 1), to = c("Negative", "Positive"))) %>%
  mutate(group = as.factor(group)) %>% {
    ggplot(., aes(mfi, assigned_group)) + 
      geom_jitter(
         aes(colour = group, shape = sd_method)  
       , width = 0.2, height = 0.2) +
      facet_grid(param_set~log_mfi, scales = "free_x") +
      ylab("Assigned Serostatus") +
      xlab("MFI") +
      scale_y_continuous(breaks = c(0, 1), labels = c("Negative", "Positive")) +
      scale_colour_manual(
        values = c("dodgerblue3", "firebrick3")
      , name = "True
Serostatus"
      ) +
      scale_shape_manual(
        values = c(0, 2)
      , name = "Method"
      ) +
      theme(
        axis.text.y = element_text(size = 12)
      , axis.text.x = element_text(size = 10)
      , strip.text.y = element_blank()
      , strip.text.x = element_blank()
      , legend.key.size = unit(0.25, "cm")
      , legend.position = c(0.85, 0.25)
      , legend.text = element_text(size = 8))
  }

these_sets <- sample(seq(500), 32)

three_sd.groups %>%
  filter(log_mfi == "log_mfi") %>% 
  filter(
     # param_set %in% c(445, 412)
     param_set %in% these_sets
    ) %>% 
  mutate(
      group = as.factor(group)
    , int = interaction(cat1f, group)
  ) %>% mutate(group = plyr::mapvalues(group, from = c(0, 1), to = c("Negative", "Positive"))) %>% {
    ggplot(., aes(x = mfi)) + 
      geom_density(aes(colour = group, fill = group), alpha = 0.2) + 
      scale_colour_manual(
        values = c("dodgerblue3", "firebrick3")
        , name = "True
Serostatus"
      ) +
      scale_fill_manual(
        values = c("dodgerblue3", "firebrick3")
        , name = "True
Serostatus"
      ) +
      facet_wrap(param_set~log_mfi, scales = "free", nrow = 8) +
      ylab("Density") +
      xlab("Log MFI") +
      theme(
        strip.text.x = element_blank()
      , axis.text.y = element_text(size = 10)
      , axis.text.x = element_text(size = 10)
      , legend.key.size = unit(0.5, "cm")
    #  , legend.position = c(0.85, 0.25)
      , legend.text = element_text(size = 8)
        )
  }

##### 
## Figure 2 and Figure X: x-axis skew of distribution, y-axis bias in seropositivity estimates
## log MFI, mclust, N-N, N-SN
## biase and coverage of stan ceof estimates across skew
#####

all.out$pop_seropositivity %>% 
  mutate(prop_pos_diff = prop_pos / true) %>%
  filter(model %in% c(
    "publication_model_normal_2.stan"
  , "publication_model_skew_normal_wf_2.stan"
  , "mclust"
  )
, method %in% c(
  "Bayesian LCR"
, "unconstrained_reduced_mclust"
)) %>% mutate(mod.meth = interaction(model, method)) %>%
  mutate(Approach = plyr::mapvalues(
    mod.meth
    , from = c(
      "mclust.unconstrained_reduced_mclust"
    , "publication_model_normal_2.stan.Bayesian LCR"
    , "publication_model_skew_normal_wf_2.stan.Bayesian LCR"
    )
    , to = c(
        "Mclust -- (Unconstrained -- BIC collapsed)"
      , "Bayesian LCR -- normal-normal"
      , "Bayesian LCR -- normal-skewnormal"
    ))
    ) %>% #left_join(
    #  .
   # , sim.data.summaries %>% dplyr::select(
   #   param_set, log_mfi
   # , prop_overlap_quant, mfi_skew_2
   # , diff_mean
  #  )) %>% 
  filter(
      quantile == "mid"
    , log_mfi == "log_mfi"
    ) %>% mutate(
      mfi_skew_2 = plyr::round_any(mfi_skew_2, 0.1)
    ) %>% group_by(
      mfi_skew_2, Approach
    ) %>% summarize(
      prop_pos_diff = mean(prop_pos_diff, na.rm = T)
    ) %>% {
    ggplot(., aes(mfi_skew_2, prop_pos_diff)) +
        geom_jitter(aes(colour = Approach), size = 3) +
        scale_colour_brewer(palette = "Dark2") +
        xlab("Seropositive Distribution Skew") +
        ylab("Estimated / True Seropositivity") +
        geom_hline(yintercept = 0, linetype = "dashed") +
      #  geom_smooth(aes(colour = Approach), se = F) +
        scale_y_log10() +
        geom_hline(yintercept = 1, linetype = "dashed")
    }

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

all_coef_ests %>% 
  mutate(m_diff = ifelse(mid < true, m_diff * -1, m_diff)) %>% 
  left_join(., sim.data.summaries %>% dplyr::select(
    param_set, log_mfi
    , prop_overlap_quant, mfi_skew_2
    , diff_mean
  )) %>%
  filter(name %notin% c(
      "mu_neg", "mu_pos", "mu_pos_delta"
    , "sd_neg", "sd_pos", "theta_con2f_delta"
  )) %>% 
  filter(main_approach %in% c("Bayesian LCR"), log_mfi == "log_mfi") %>% 
  mutate(
    mfi_skew_2 = plyr::round_any(mfi_skew_2, 0.05)
  ) %>% 
  rename(Approach = sub_approach) %>%
  mutate(Approach = plyr::mapvalues(Approach, from = c(
    "normal-normal", "normal-skewnormal"
  ), to = c(
    "Bayesian LCR (normal-normal)", "Bayesian LCR (normal-skewnormal)"
  ))) %>%
  group_by(
    mfi_skew_2, Approach, name
  ) %>% summarize(
    m_diff = mean(m_diff, na.rm = T)
  ) %>% mutate(name = plyr::mapvalues(
    name
    , from = c("beta_base", "beta_cat1f_delta", "beta_cat2f_delta", "beta_con1f_delta")
    , to   = c("Baseline (Intercept)", "Categorical Covariate 1"
               , "Categorical Covariate 2", "Continuous Covariate 1")
  )) %>% {
    ggplot(., aes(mfi_skew_2, m_diff)) + 
      geom_point(aes(colour = Approach)) +
      scale_colour_brewer(palette = "Dark2") +
      geom_smooth(aes(colour = Approach), se = F) +
      facet_wrap(~name, scales = "free") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      xlab("Seropositive Skew") +
      ylab("Estimated - True Coefficient Estimate")
  }

all_coef_ests %>% 
  mutate(m_diff = ifelse(mid < true, m_diff * -1, m_diff)) %>% 
  left_join(., sim.data.summaries %>% dplyr::select(
    param_set, log_mfi
    , prop_overlap_quant, mfi_skew_2
    , diff_mean
  )) %>%
  filter(name %notin% c(
     "mu_neg", "mu_pos", "mu_pos_delta"
    , "sd_neg", "sd_pos", "theta_con2f_delta"
  )) %>% 
  filter(main_approach %in% c("Bayesian LCR"), log_mfi == "log_mfi") %>% 
  mutate(
    mfi_skew_2 = plyr::round_any(mfi_skew_2, 0.2)
  ) %>% 
  rename(Approach = sub_approach) %>%
  mutate(Approach = plyr::mapvalues(Approach, from = c(
    "normal-normal", "normal-skewnormal"
  ), to = c(
    "Bayesian LCR (normal-normal)", "Bayesian LCR (normal-skewnormal)"
  ))) %>%
  group_by(
    mfi_skew_2, Approach, name
  ) %>% summarize(
    cover = mean(cover, na.rm = T)
  ) %>% mutate(name = plyr::mapvalues(
    name
    , from = c("beta_base", "beta_cat1f_delta", "beta_cat2f_delta", "beta_con1f_delta")
    , to   = c("Baseline (Intercept)", "Categorical Covariate 1"
              , "Categorical Covariate 2", "Continuous Covariate 1")
  )) %>% {
    ggplot(., aes(mfi_skew_2, cover)) + 
      geom_point(aes(colour = Approach)) +
      scale_colour_brewer(palette = "Dark2") +
      geom_smooth(aes(colour = Approach), se = F) +
      facet_wrap(~name, scales = "free") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      xlab("Seropositive Skew") +
      ylab("Estimated - True Coefficient Estimate")
  }

##### 
## Figure X: Narrowing down to a small subset of parameter space for 3sd vs stan
#####

all_compare1 <- all.out$pop_seropositivity %>% 
  mutate(prop_pos_diff = prop_pos / true) %>%
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
  mutate(m.s = factor(m.s, levels = unique(m.s))) %>% 
  left_join(., sim.data.summaries %>% dplyr::select(
      param_set, log_mfi
    , prop_overlap_quant, mfi_skew_2
    , diff_mean, prop_pos
  ) %>% rename(prop_pos_true = prop_pos)) %>% 
  dplyr::select(-prop_pos_diff) %>%
  pivot_wider(., values_from = prop_pos, names_from = quantile) %>%
  mutate(cover = ifelse(lwr < true & upr > true, 1, 0)) %>%
  mutate(
      prop_pos_diff = mid / true
    , CI_wid = upr - lwr
  ) %>% 
  filter(prop_pos_true < 0.05, prop_overlap_quant < 0.2) %>%
  group_by(mod.meth, main_approach, log_mfi) %>%
  summarize(
      tot_cov = mean(cover, na.rm = T)
    , m_diff  = mean(prop_pos_diff, na.rm = T)
  )

all_compare2 <- all_coef_ests %>% 
  mutate(m_diff = ifelse(mid < true, m_diff * -1, m_diff)) %>% 
  left_join(., sim.data.summaries %>% dplyr::select(
    param_set, log_mfi
    , prop_overlap_quant, mfi_skew_2
    , diff_mean, prop_pos
  )) %>%
  filter(name %notin% c(
    "mu_neg", "mu_pos", "mu_pos_delta"
    , "sd_neg", "sd_pos", "theta_con2f_delta"
  )) %>% filter(prop_pos < 0.05, prop_overlap_quant < 0.2) %>% 
  group_by(mod.meth, log_mfi, main_approach, name) %>%
  summarize(
     tot_cov = mean(cover, na.rm = T)
    , m_diff  = mean(m_diff, na.rm = T)
  )

all_compare1 %<>% 
  droplevels() %>% 
  ungroup() %>%
  mutate(
    Approach = plyr::mapvalues(mod.meth, from = unique(mod.meth), to = c(
      "3sd (Approximated control)"                          
    , "3sd (Robust Mean and SD)"                           
    , "Bayesian LCR (lognormal-lognormal)"     
    , "Bayesian LCR (normal-normal)"        
    , "Bayesian LCR (normal-skewnormal)"
    , "Mclust (2 group constrained)"                           
    , "Mclust (Unconstrained; sum of 2-n)"                         
    , "Mclust (Unconstrained; BIC collapsed)"  
  ))
) %>% mutate(log_mfi = plyr::mapvalues(log_mfi, from = c(
  "log_mfi", "mfi"
  ), to = c("Log MFI", "Linear MFI"))
)

all_compare2 %<>% 
  droplevels() %>%
  ungroup() %>% mutate(
    Approach = plyr::mapvalues(mod.meth, from = unique(mod.meth), to = c(
      "3sd (Approximated control; Unweighted GLM)"                          
    , "3sd (Robust Mean and SD; Unweighted GLM)"                           
    , "Bayesian LCR (lognormal-lognormal)"     
    , "Bayesian LCR (normal-normal)"        
    , "Bayesian LCR (normal-skewnormal)"
    , "Mclust (2 group constrained -- Unweighted GLM)"  
    , "Mclust (2 group constrained -- Weighted GLM)"
    , "Mclust (Unconstrained (sum of 2-n) -- Unweighted GLM)"                         
    , "Mclust (Unconstrained (sum of 2-n) -- Weighted GLM)"  
    , "Mclust (Unconstrained (BIC collapsed) -- Unweighted GLM)"                    
    , "Mclust (Unconstrained (BIC collapsed) -- Weighted GLM)"
  ))) %>% mutate(log_mfi = plyr::mapvalues(log_mfi, from = c(
  "log_mfi", "mfi"
  ), to = c("Log MFI", "Linear MFI"))
) %>% mutate(name = plyr::mapvalues(
    name
  , from = c("beta_base", "beta_cat1f_delta", "beta_cat2f_delta", "beta_con1f_delta")
  , to   = c("Baseline (Intercept)", "Categorical Covariate 1"
             , "Categorical Covariate 2", "Continuous Covariate 1")
))

gg.1 <- all_compare1 %>% {
  ggplot(., aes(Approach, tot_cov)) + 
    geom_point(aes(colour = Approach, shape = log_mfi), size = 3) +
    xlab("") +
    ylab("Coverage") +
    theme(axis.text.x = element_blank()) +
    scale_shape_discrete(name = "MFI Scale")
}

gg.2 <- all_compare1 %>% {
  ggplot(., aes(Approach, m_diff)) + 
    geom_point(aes(colour = Approach, shape = log_mfi), size = 3) +
    xlab("") +
    ylab("Estimated / True Coefficient Estimate") +
    theme(axis.text.x = element_blank()) +
    scale_y_log10(breaks = c(0.7, 1, 3, 10)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_shape_discrete(name = "MFI Scale")
}

gridExtra::grid.arrange(gg.1, gg.2, ncol = 1)

gg.1 <- all_compare1 %>% 
  mutate(Approach = factor(Approach, levels = c(
    unique(param_coverage_by_model.all.adj$m.s)
  ))) %>% {
  ggplot(., aes(tot_cov, Approach)) + 
    geom_point(aes(colour = main_approach, shape = log_mfi), size = 3) +
    ylab("Model") +
    xlab("Coverage") +
    scale_colour_manual(values = c(
        "#1b9e77"
      , "#7570b3"
      , "#e6ab02"
    ), name = "Approach") +
    theme(
        strip.text.x = element_text(size = 11)
      , plot.margin = unit(c(0.2,0.2,0.2,0.05), "cm")
      , axis.text.x = element_text(size = 12)
      , axis.text.y = element_text(size = 10)
      , legend.position = "none"
    ) +
    scale_x_continuous(
        breaks = c(0.1, 0.35, 0.60, 0.85, 0.95)
      , labels = c("10%", "35%", "60%", "85%", "95%")
      , lim = c(0, 1)) +
    scale_shape_discrete(name = "MFI Scale") +
    geom_vline(xintercept = 0.95, linetype = "dashed")
}

gg.2 <- all_compare1 %>% 
  mutate(Approach = factor(Approach, levels = c(
    unique(param_coverage_by_model.all.adj$m.s)
  ))) %>% {
  ggplot(., aes(m_diff, Approach)) + 
    geom_point(aes(colour = main_approach, shape = log_mfi), size = 3) +
    ylab("") +
    xlab("Bias") +
      scale_colour_manual(values = c(
          "#1b9e77"
        , "#7570b3"
        , "#e6ab02"
      ), name = "Approach") +
      theme(
          strip.text.x = element_text(size = 11)
        , plot.margin  = unit(c(0.2,0.2,0.2,0.05), "cm")
        , axis.text.x  = element_text(size = 12)
        , axis.text.y  = element_blank()
        , axis.ticks.length.y = unit(0, "cm")
      ) +
    scale_x_log10(breaks = c(0.7, 1, 3, 10)) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_shape_discrete(name = "MFI Scale")
  }

gg.2 <- all_coefs_for_gg.2 %>% 
  mutate(m.s = factor(m.s, levels = levvs)) %>% 
  rename(Approach = m.s) %>%
  filter(quantile == "mid") %>%
  filter(true < 0.05, prop_overlap_quant < 0.2) %>% {
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
        breaks = c(0.5, 1, 2.5, 5, 10, 20)
        , labels = c("0.5", "1.00", "2.5", "5", "10", "20")) +
      ylab("") +
      xlab("Bias") +
      geom_vline(xintercept = 1, linetype = "dashed")
  }

gridExtra::grid.arrange(gg.1, gg.2, ncol = 2)

gg.1 <- all_compare2 %>% {
    ggplot(., aes(Approach, tot_cov)) + 
      geom_point(aes(colour = Approach, shape = log_mfi), size = 3) +
      facet_wrap(~name, nrow = 1) +
      scale_y_continuous(
        breaks = c(0, 0.25, 0.5, 0.75, 1)
      , labels = c("0.000", "0.250", "0.500", "0.750", "1.000")
      ) +
      xlab("") +
      ylab("Coverage") +
      theme(
        axis.text.x = element_blank()
      , strip.text.x = element_text(size = 10)
      ) +
      scale_shape_discrete(name = "MFI Scale")
}

to_vals <- c(
    "3sd -- (Approximated control)"                                             
  , "3sd -- (Robust Mean and SD)"
  , "Bayesian LCR -- lognormal-lognormal"
  , "Bayesian LCR -- normal-normal"
  , "Bayesian LCR -- normal-skewnormal"
  , "Mclust -- (2 group constrained) + unweighted regression"                   
  , "Mclust -- (2 group constrained) + probability weighted regression"
  , "Mclust -- (Unconstrained; sum of 2-n) + unweighted regression"             
  , "Mclust -- (Unconstrained; sum of 2-n) + probability weighted regression"   
  , "Mclust -- (Unconstrained; BIC collapsed) + unweighted regression"          
  , "Mclust -- (Unconstrained; BIC collapsed) + probability weighted regression"
)

gg.1 <- all_compare2 %>% 
  mutate(
  Approach = plyr::mapvalues(Approach, from = unique(Approach), to = to_vals)
) %>% mutate(
  Approach = factor(Approach, levels = c(
    coef.cover$m.s %>% levels()
  ))
) %>% {
    ggplot(., aes(Approach, tot_cov)) + 
      geom_point(aes(colour = main_approach, shape = log_mfi), size = 3) +
      facet_wrap(~name, nrow = 1) +
      xlab("") +
      ylab("Coverage") +
      scale_colour_manual(values = c(
          "#1b9e77"
        , "#7570b3"
        , "#e6ab02"
      ), name = "Approach") +
      theme(
          strip.text.x = element_text(size = 11)
        , plot.margin = unit(c(0.2,0.2,0.0,0.05), "cm")
        , axis.text.x = element_blank()
        , axis.ticks.x = element_blank()
        , axis.text.y = element_text(size = 10)
        , legend.position = "none"
      ) +
      scale_y_continuous(
        breaks = c(0.1, 0.35, 0.60, 0.85, 0.95)
        , labels = c("10%", "35%", "60%", "85%", "95%")
        , lim = c(0, 1)) +
      scale_shape_discrete(name = "MFI Scale") +
      geom_hline(yintercept = 0.95, linetype = "dashed")
  }

gg.2 <- all_compare2 %>% {
  ggplot(., aes(Approach, m_diff)) + 
    geom_point(aes(colour = Approach, shape = log_mfi), size = 3) +
    facet_wrap(~name, nrow = 1) +
    scale_y_continuous(breaks = c(-1.000, -0.500, 0.500, 1.000, 2.500, 5.000)) +
    xlab("") +
    ylab("Estimated - True Coefficient Estimate") +
    theme(
      axis.text.x = element_blank()
    , strip.text.x = element_blank()
    , legend.position = "none"
    , plot.margin = unit(c(0,9.9,1,0.59), "cm")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_shape_discrete(name = "MFI Scale")
}

parm_sets <- all_coefs_for_gg.2 %>% 
  mutate(m.s = factor(m.s, levels = levvs)) %>% 
  rename(Approach = m.s) %>%
  filter(quantile == "mid") %>%
  filter(true < 0.05, prop_overlap_quant < 0.2) %>%
  pull(param_set) %>%
  unique()

these_to_vals <- to_vals[c(
  4, 5, 3, 6, 7
, 8, 9, 10, 11
, 1, 2
)]

gg.2 <- all_coef_ests %>%
  filter(param_set %in% parm_sets) %>%
  mutate(
    Approach = plyr::mapvalues(mod.meth, from = unique(mod.meth), to = these_to_vals)
  ) %>% mutate(
    Approach = factor(Approach, levels = c(
      coef.cover$m.s %>% levels()
    ))
  ) %>% filter(
    m_diff > -5, m_diff < 5
  ) %>% {
  ggplot(., aes(Approach, m_diff)) +
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
          strip.text.x = element_blank()
        , plot.margin = unit(c(0.0,0.2,0.2,0.15), "cm")
       # , axis.text.x = element_text(size = 8, angle = 330, hjust = 0, vjust = 0.5)
        , axis.text.x = element_blank()
        , axis.text.y = element_text(size = 10)
        , legend.position = "none"
      ) +
    xlab("Model") +
    ylab("Bias") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~name, nrow = 1)
}

gridExtra::grid.arrange(gg.1, gg.2, ncol = 1)

#####
## Figure SX: Impact of sample size
#####

## See main_approach_prop_sero for required code

all_coefs_for_gg.2 %>% 
  mutate(m.s = factor(m.s, levels = levvs)) %>% 
  rename(Approach = m.s) %>% 
  filter(quantile == "mid") %>% {
    ggplot(., aes(n_samps, prop_pos_diff)) +
      geom_point(aes(colour = Approach), alpha = 0.5) +
      theme(
        strip.text.x = element_text(size = 12)
      , plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
      ) +
      scale_y_log10(
        breaks = c(.001, .1, 1, 10, 50)
      , labels = c("0.001", "0.01", "1.00", "10", "50")) +
      xlab("Sample Size") +
      ylab("Estimated / True Seropositivity") +
      geom_hline(yintercept = 1, linetype = "dashed") +
      facet_wrap(~Approach, ncol = 2)
  }

## See main_approach_coverage for required code

ci_wid_by_model %>% 
  mutate(m.s = interaction(main_approach, sub_approach, sep = " -- ")) %>%
  mutate(m.s = plyr::mapvalues(m.s, from = c(
    "Mclust -- mclust (2 group constrained) + unweighted regression"
  , "Mclust -- mclust (2 group constrained) + probability weighted regression"
  , "Mclust -- mclust (Unconstrained) + unweighted regression"
  , "Mclust -- mclust (Unconstrained) + probability weighted regression"
  , "Mclust -- mclust (Unconstrained reduced) + unweighted regression"
  , "Mclust -- mclust (Unconstrained reduced) + probability weighted regression"
  ), to = c(
    "Mclust -- mclust (2 group constrained) + 
unweighted regression"
  , "Mclust -- mclust (2 group constrained) + 
probability weighted regression"
  , "Mclust -- mclust (Unconstrained) + 
unweighted regression"
  , "Mclust -- mclust (Unconstrained) + 
probability weighted regression"
  , "Mclust -- mclust (Unconstrained reduced) + 
unweighted regression"
  , "Mclust -- mclust (Unconstrained reduced) + 
probability weighted regression"
  ))) %>%
  mutate(log_mfi = plyr::mapvalues(log_mfi, from = c(
    "log_mfi", "mfi"
  ), to = c("Log MFI", "Linear MFI"))) %>%
  mutate(m.s = factor(m.s, levels = unique(m.s))) %>%
  rename(
    Approach = m.s
  , `Main Approach` = main_approach
  , `MFI Scale` = log_mfi) %>% 
  filter(ci_wid_mean < 200) %>% {
    ggplot(., aes(n_samps, ci_wid_mean)) +
      geom_point(aes(colour = `Main Approach`, shape = `MFI Scale`), alpha = 0.5) +
      theme(
        strip.text.x = element_text(size = 9)
      , plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
      ) +
      scale_y_log10(
        breaks = c(0.4, 2, 10, 50, 200)
      , labels = c("0.4", "2.00", "10", "50", "200")) +
      xlab("Sample Size") +
      ylab("CI Width") +
      facet_wrap(~Approach, ncol = 3)
  }

sim.data.summaries.j <- sim.data.summaries %>%
  dplyr::select(
    param_set, sim_num, log_mfi
    , prop_overlap_quant, diff_mean, mfi_skew_2
    , prop_squish_2, prop_pos, n_samps
  ) %>% rename(prop_pos_true = prop_pos) %>%
  mutate(n_samps = round(n_samps))

param_coverage_by_model %>% 
  left_join(., left_join(., sim.data.summaries.j %>% dplyr::select(param_set, n_samps) %>% 
                           group_by(param_set) %>% slice(1) %>% ungroup())) %>%
  mutate(n_samps = plyr::round_any(n_samps, 100)) %>%
  group_by(n_samps, log_mfi, main_approach, sub_approach) %>%
  summarize(m_cov = mean(cover)) %>% 
  ungroup() %>%
  mutate(m.s = interaction(main_approach, sub_approach, sep = " -- ")) %>%
  mutate(m.s = plyr::mapvalues(m.s, from = c(
    "Mclust -- mclust (2 group constrained) + unweighted regression"
    , "Mclust -- mclust (2 group constrained) + probability weighted regression"
    , "Mclust -- mclust (Unconstrained) + unweighted regression"
    , "Mclust -- mclust (Unconstrained) + probability weighted regression"
    , "Mclust -- mclust (Unconstrained reduced) + unweighted regression"
    , "Mclust -- mclust (Unconstrained reduced) + probability weighted regression"
  ), to = c(
    "Mclust -- mclust (2 group constrained) + 
unweighted regression"
    , "Mclust -- mclust (2 group constrained) + 
probability weighted regression"
    , "Mclust -- mclust (Unconstrained) + 
unweighted regression"
    , "Mclust -- mclust (Unconstrained) + 
probability weighted regression"
    , "Mclust -- mclust (Unconstrained reduced) + 
unweighted regression"
    , "Mclust -- mclust (Unconstrained reduced) + 
probability weighted regression"
  ))) %>%
  mutate(log_mfi = plyr::mapvalues(log_mfi, from = c(
    "log_mfi", "mfi"
  ), to = c("Log MFI", "Linear MFI"))) %>%
  mutate(m.s = factor(m.s, levels = unique(m.s))) %>%
  rename(
      Approach = m.s
    , `Main Approach` = main_approach
    , `MFI Scale` = log_mfi) %>%
   mutate(m_cov = m_cov / 4) %>% {
      ggplot(., aes(n_samps, m_cov)) +
        geom_point(aes(colour = `Main Approach`, shape = `MFI Scale`)) +
        theme(
           strip.text.x = element_text(size = 9)
         , plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
        ) +
        xlab("Sample Size") +
        ylab("Coverage") +
       geom_hline(yintercept = 0.95, linetype = "dashed") +
        facet_wrap(~Approach, ncol = 3)
  }

##### 
## Supplement Figure 3: Bias in pop estimates across two distributional summaries
#####

all.out$pop_seropositivity %>% 
  filter(model == "3sd") %>% 
  mutate(method = plyr::mapvalues(method, from = c(
    "assigned_group_control", "assigned_group_robust"
  ), to = c(
    "Approximated
Control", "Robust Mean
and SD"
  ))) %>% 
  mutate(log_mfi = plyr::mapvalues(log_mfi, from = c(
    "log_mfi", "mfi"
  ), to = c("Log MFI", "Linear MFI"))) %>% {
    ggplot(., aes(true, prop_pos)) + 
      geom_point(aes(colour = method)) +
      geom_abline(intercept = 0, slope = 1) +
      facet_wrap(~log_mfi) +
      scale_colour_brewer(palette = "Dark2", name = "Method") +
      xlab("True Proportion Positive") +
      ylab("Estimated Proportion Positive")
  }

all.out$pop_seropositivity %>% 
  mutate(prop_pos_diff = prop_pos / true) %>%
  left_join(., sim.data.summaries %>% dplyr::select(
    param_set, log_mfi
  , prop_overlap_quant, mfi_skew_2
  , diff_mean
  )) %>% 
  filter(model == "3sd") %>%
  mutate(method = plyr::mapvalues(method, from = c(
    "assigned_group_control", "assigned_group_robust"
  ), to = c(
    "Approximated
Control", "Robust Mean
and SD"
  ))) %>% 
  mutate(log_mfi = plyr::mapvalues(log_mfi, from = c(
    "log_mfi", "mfi"
  ), to = c("Log MFI", "Linear MFI"))) %>% 
  dplyr::select(
    method, log_mfi, prop_pos_diff
  , prop_overlap_quant
  , mfi_skew_2
  , diff_mean
  ) %>%
  pivot_longer(-c(method, log_mfi, prop_pos_diff)) %>% 
  mutate(name = plyr::mapvalues(
    name
  , from = c("prop_overlap_quant", "mfi_skew_2", "diff_mean")
  , to   = c("Distribution Overlap", "Seropositive Skew", "Difference in MFI Means")
  )) %>% {
  ggplot(., aes(value, prop_pos_diff)) +
      geom_point(aes(colour = method, shape = log_mfi)) +
      scale_colour_brewer(palette = "Dark2", name = "Method") +
      facet_wrap(~name, scales = "free_x") +
      scale_shape_discrete(name = "MFI Scale") +
      xlab("Value") +
      ylab("Estimate / True Proportion Positive") +
      theme(axis.text.x = element_text(size = 12, angle = 300, hjust = 0)) +
      scale_y_log10() +
      geom_hline(yintercept = 1, linetype = "dashed")
}

##### 
## Supplement Figure X: More linkages between individual and population level seropositivity estimates
#####

source("scr_figure_plot.R")

all.out$pop_seropositivity %>% 
  filter(model == "3sd") %>%
  filter(true < 0.1) %>% 
  mutate(log_mfi = plyr::mapvalues(log_mfi, from = c(
    "log_mfi", "mfi"
  ), to = c("Log MFI", "Raw MFI"))) %>%
  mutate(method = plyr::mapvalues(method, from = c(
    "assigned_group_control", "assigned_group_robust"
  ), to = c(
    "Approximated Control", "Robust Mean and SD"
  ))) %>% {
    ggplot(., aes(prop_overlap_quant, prop_pos_diff)) + 
      geom_point(aes(colour = log_mfi, shape = method), size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_colour_brewer(palette = "Dark2", name = "MFI") +
      scale_shape_discrete(name = "Method") +
      xlab("Distributional Overlap") +
      ylab("Estimated - True Proportion Positive")
  } 

#####
## Supplemental Figure X and Y: comparing mclust approaches on pop seropositivity and coef estimates
#####

all.out$pop_seropositivity %>% 
  left_join(., sim.data.summaries %>% dplyr::select(
    param_set, log_mfi
    , prop_overlap_quant, mfi_skew_2
    , diff_mean
  )) %>%
  filter(model == "mclust") %>%
  mutate(method = plyr::mapvalues(method, from = c(
      "constrained_mclust"
    , "unconstrained_mclust"
    , "unconstrained_reduced_mclust"
  ), to = c(
      "2 group constrained"
    , "Unconstrained; sum of 2-n"
    , "Unconstrained; BIC collapsed"
  ))
  , log_mfi = plyr::mapvalues(log_mfi, from = c(
    "log_mfi", "mfi"
  ), to = c("Log MFI", "Linear MFI"))) %>% {
    ggplot(., aes(true, prop_pos)) + 
      geom_point(aes(colour = prop_overlap_quant)) +
      geom_abline(intercept = 0, slope = 1) +
      scale_colour_viridis(name = "Distributional
Overlap") +
      facet_grid(log_mfi~method) +
      theme(
        strip.text.x = element_text(size = 12)
      , strip.text.y = element_text(size = 14)
      ) +
      ylab("Estimated Proportion Positive") +
      xlab("True Proportion Positive")
  }

all.out$coefficient_ests %>% 
  mutate(m_diff = ifelse(mid < true, m_diff * -1, m_diff)) %>% 
  filter(mid > -10, mid < 3) %>%
  left_join(., sim.data.summaries %>% dplyr::select(
    param_set, log_mfi
    , prop_overlap_quant, mfi_skew_2
    , diff_mean
  )) %>%
  filter(
    model == "mclust -- variance"
  , name  == "beta_base"
  ) %>% mutate(method = plyr::mapvalues(method, from = c(
    "constrained_mclust"
    , "unconstrained_mclust"
    , "unconstrained_reduced_mclust"
  ), to = c(
    "2 group constrained"
    , "Unconstrained; sum of 2-n"
    , "Unconstrained; BIC collapsed"
  ))
  , log_mfi = plyr::mapvalues(log_mfi, from = c(
    "log_mfi", "mfi"
  ), to = c("Log MFI", "Linear MFI"))) %>% {
    ggplot(., aes(true, mid)) + 
      geom_point(aes(colour = prop_overlap_quant)) +
      geom_abline(intercept = 0, slope = 1) +
      scale_colour_viridis(name = "Distributional
Overlap") +
      facet_grid(log_mfi~method) +
      theme(
        strip.text.x = element_text(size = 12)
        , strip.text.y = element_text(size = 14)
      ) +
      ylab("Estimated Categorical Covariate Effect Size") +
      xlab("True Categorical Covariate Effect Size")
  }


#####
## Figure SX: Individual mclust probabilities
#####

this_param_set <- c(31, 198, 237, 298, 202)

mclust.groups %>% 
  filter(param_set %in% this_param_set) %>% 
  dplyr::select(group, mfi, method, param_set, log_mfi, V1, V2, V2_adj) %>%
  mutate(Approach = plyr::mapvalues(method, from = c(
    "constrained_mclust"
    , "unconstrained_mclust"
    , "unconstrained_reduced_mclust"
  ), to = c(
    "2 group constrained"
    , "Unconstrained; sum of 2-n"
    , "Unconstrained; BIC collapsed"
  ))
  , log_mfi = plyr::mapvalues(log_mfi, from = c(
    "log_mfi", "mfi"
  ), to = c("Log MFI", "Linear MFI"))) %>% {
    ggplot(., aes(mfi, V2_adj)) + 
      geom_jitter(aes(colour = Approach, shape = Approach)) +
      scale_colour_brewer(palette = "Dark2") +
      facet_grid(param_set~log_mfi, scales = "free_x") +
      scale_x_log10() +
      xlab("MFI") +
      ylab("Positive Probability") +
      theme(
        axis.text.y = element_text(size = 10)
        , strip.text.y = element_blank())
  }

#####
## Figure SX: various alternative methods of determining an optimal empirical MFI mean and se
#####

source("scr_alt_3sd_trans.R")

#####
## Figure SX: Stan fitting
#####

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

failed_fits <- stan.summary$fit_details %>% 
  left_join(., sim.data.summaries.j) %>%
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



sim.data.summaries.j %>% filter(
  param_set %in% failed_fits
)

three_sd.groups %>% 
  filter(param_set %in% failed_fits) %>% 
  mutate(
    group = as.factor(group)
    , int = interaction(cat1f, group)
  ) %>% mutate(group = plyr::mapvalues(group, from = c(0, 1), to = c("Negative", "Positive"))) %>% {
    ggplot(., aes(x = mfi)) + 
      geom_density(aes(colour = group, fill = group), alpha = 0.2) + 
      scale_colour_manual(
        values = c("dodgerblue3", "firebrick3")
        , name = "True
Serostatus"
      ) +
      scale_fill_manual(
        values = c("dodgerblue3", "firebrick3")
        , name = "True
Serostatus"
      ) +
      facet_wrap(param_set~log_mfi, scales = "free", ncol = 2) +
      ylab("Density") +
      xlab("MFI") +
      theme(
        strip.text.x = element_blank()
        , axis.text.y = element_text(size = 10)
        , axis.text.x = element_text(size = 10)
        , legend.key.size = unit(0.5, "cm")
        , legend.text = element_text(size = 8)
      )
  }

failed_fits <- stan.summary$fit_details %>% 
  left_join(., sim.data.summaries.j) %>%
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
  filter(m_R > 5) %>% 
  mutate(failed = 1)

stan.fit_stats %>% 
  left_join(., sim.data.summaries.j) %>%
  left_join(., failed_fits %>% dplyr::select(model_name, param_set, log_mfi)) %>%
  mutate(failed = ifelse(is.na(failed), 0, 1)) %>% 
  group_by(model_name, log_mfi, failed) %>%
  summarize(
    n_failed = n()
  , prop_overlap_quant = mean(prop_overlap_quant)
  , diff_mean = mean(diff_mean)
  , mfi_skew_2 = mean(mfi_skew_2) 
  , prop_squish_2 = mean(prop_squish_2)
  , prop_pos_true = mean(prop_pos_true)
  )

stan.fit_stats %>% 
  mutate(m_R = ifelse(m_R > 5, 5, m_R)) %>%
  mutate() %>% {
  ggplot(., aes(mfi_skew_2, m_R)) +
    geom_point(aes(colour = model_name, shape = log_mfi)) +
    scale_y_log10()
}


##### mclust method and log vs not log on pop seropositivity ----- 

library(viridis)

## mfi_skew_2
## prop_pos_true
## prop_overlap_quant

all.out$pop_seropositivity %>% 
  filter(model == "mclust") %>% 
  mutate(log_mfi = plyr::mapvalues(log_mfi, from = c(
    "log_mfi", "mfi"
  ), to = c("Log MFI", "Raw MFI"))) %>%
  mutate(method = plyr::mapvalues(method, from = c(
    "constrained_mclust", "unconstrained_mclust", "unconstrained_reduced_mclust"
  ), to = c(
    "2 groups fitted", "2-n group probabilities summed", "n groups collapsed to 2 groups"
  ))) %>% filter(method != "2-n group probabilities summed") %>%
  {
    ggplot(., aes(prop_pos_true, prop_pos)) + 
      geom_point(aes(colour = prop_overlap_quant)) +
      geom_abline(intercept = 0, slope = 1) +
      scale_colour_viridis(name = "Distributional
Overlap") +
      facet_grid(log_mfi~method) +
      xlab("True Proportion Positive") +
      ylab("Estimated Proportion Positive") +
      theme(
        strip.text.x = element_text(size = 12)
      , axis.text.x = element_text(size = 10)
      , axis.text.y = element_text(size = 10)
      )
  }

##### LCR method on pop seropositivity ----- 

all.out$pop_seropositivity %>% 


##### Comparing all methods on pop seropositivity ----- 

unique(all.out$pop_seropositivity$method)

all.out$pop_seropositivity %>% 
  filter(method %in% c("assigned_group_robust", "unconstrained_reduced_mclust")) %>%
  dplyr::select(
    param_set, log_mfi, method, true, prop_pos_diff
  , prop_pos_true, prop_overlap_quant, mfi_skew_2
    ) %>% 
  mutate(log_mfi = plyr::mapvalues(log_mfi, from = c(
    "log_mfi", "mfi"
  ), to = c("Log MFI", "Raw MFI"))) %>%
  pivot_wider(names_from = method, values_from = prop_pos_diff) %>% {
    ggplot(., aes(assigned_group_robust, unconstrained_reduced_mclust)) + 
      geom_point(aes(colour =  prop_overlap_quant)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_hline(yintercept = 0, linetype = "dotted") +
      xlab("Bias in Population Seropositivity -- 3SD") +
      ylab("Bias in Population Seropositivity -- Mclust") +
      scale_colour_viridis(name = "Distributional
Overlap") +
      theme(axis.text.x = element_text(hjust = 0, angle = 300)) +
      facet_wrap(~log_mfi)
  }

source("main_approach_pop_sero.R")

##### 3sd method and log vs not log on coef estimates ----- 

which_param <- "beta_base"
which_model <- "3sd -- no_variance"

all.out$coefficient_ests %>% 
  filter(model == which_model) %>% 
  filter(name == which_param) %>% 
  arrange(true) %>% 
  mutate(param_set = factor(param_set, levels = unique(param_set))) %>%
  filter(CI_wid < 10, mid > -20) %>%
  mutate(cover = as.factor(cover)) %>% {
    ggplot(., aes(true, param_set)) + 
      geom_errorbarh(aes(xmin = lwr, xmax = upr)) +
      geom_point(aes(colour = cover)) +
      facet_wrap(~log_mfi+method, scale = "free_y") +
      scale_colour_manual(values = c("firebrick3", "dodgerblue4")) +
      theme(axis.text.y = element_blank())
  }

##### mclust method and log vs not log on coef estimates ----- 



##### LCR method and log vs not log on coef estimates ----- 

stan.summary$prop_seropos %>% 
  left_join(., sim.data.summaries) %>%
  mutate(m_diff = ifelse(mid < true, m_diff * -1, m_diff)) %>% 
  mutate(cover = as.factor(cover)) %>% {
    ggplot(., aes(mfi_skew_2, m_diff)) +
      geom_point(aes(colour = stan_model, shape = cover)) +
      facet_wrap(~log_mfi, scales = "free_x") +
      geom_hline(yintercept = 0)
  }

stan.summary$prop_seropos %>% 
  left_join(., sim.data.summaries) %>%
  mutate(mfi_skew_2 = plyr::round_any(mfi_skew_2, 0.2)) %>%
  group_by(mfi_skew_2, log_mfi, stan_model) %>%
  summarize(cover = mean(cover)) %>% {
    ggplot(., aes(mfi_skew_2, cover)) +
      geom_point(aes(colour = stan_model), size = 3) +
      facet_wrap(~log_mfi, scales = "free_x")
  }

##### Comparing all methods on coefficient estimates ----- 

source("main_approach_coverage.R")

############################################################
############################################################
############################################################
############################################################

##### Parameter space exploration ----- 

sim.data.summaries %>% 
  ungroup() %>% 
  filter(log_mfi == "log_mfi") %>% dplyr::select(
      prop_pos, diff_mean, prop_overlap_quant
    , mfi_skew_2, prop_squish_2, p_abv_m_2
  ) %>% 
  rename(
      `Proportion 
Seropositive` = prop_pos
    , `Difference 
in MFI Means` = diff_mean
    , `Distribution 
Overlap` = prop_overlap_quant
    , `Seropositive 
Skew` = mfi_skew_2
    , `Quantile 
Compression` = prop_squish_2
    , `Proportion 
Above Median` = p_abv_m_2
  ) %>% {
    ggpairs(
      .
    , lower = list(continuous = wrap("points", alpha = 0.3), combo = wrap("dot_no_facet", alpha = 0.4))
    ) +
      theme(
          axis.text.x = element_text(size = 10, hjust = 0, angle = 300)
        , axis.text.y = element_text(size = 10)
        , strip.text.x = element_text(size = 10)
        , strip.text.y = element_text(size = 10)
      ) 
  }

##### LCR method on fitting success ----- 

sim.data %>% 
  filter(param_set %in% inf_rhats[1:20]) %>%
  filter(log_mfi == "log_mfi") %>%
  mutate(group = as.factor(group)) %>% {
    ggplot(., aes(x = mfi)) + 
      geom_density(aes(colour = group, fill = group), alpha = 0.2) + 
      scale_colour_brewer(palette = "Dark2") +
      scale_fill_brewer(palette = "Dark2") +
      facet_wrap(~param_set, scales = "free", nrow = 3)
  }

##### Mclust collapsing method on probability space ----- 

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



