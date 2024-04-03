##### 
## Figure 1: 3sd approaches, log and linear
#####

all.out$pop_seropositivity %>% 
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
          filter(model == "3sd") %>%
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
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_hline(yintercept = 0, linetype = "dotted") +
      xlab("Estimated - True Seropositivity (MFI)") +
      ylab("Estimated - True Seropositivity (Log MFI)") +
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

three_sd.groups %>% 
  filter(param_set %in% c(445, 412)) %>% 
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
      facet_wrap(param_set~log_mfi, scales = "free", nrow = 3) +
      ylab("Density") +
      xlab("MFI") +
      theme(
        strip.text.x = element_blank()
      , axis.text.y = element_text(size = 10)
      , axis.text.x = element_text(size = 10)
      , legend.key.size = unit(0.5, "cm")
      , legend.position = c(0.85, 0.25)
      , legend.text = element_text(size = 8)
        )
  }

##### 
## Figure 2 and Figure X: x-axis skew of distribution, y-axis bias in seropositivity estimates
## log MFI, mclust, N-N, N-SN
## log MFI, all mclust methods
#####

all.out$pop_seropositivity %>% 
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
    ) %>% left_join(
      .
    , sim.data.summaries %>% dplyr::select(
      param_set, log_mfi
    , prop_overlap_quant, mfi_skew_2
    , diff_mean
    )) %>% filter(
      quantile == "mid"
    , log_mfi == "log_mfi"
    ) %>% mutate(
      mfi_skew_2 = plyr::round_any(mfi_skew_2, 0.05)
    ) %>% group_by(
      mfi_skew_2, Approach
    ) %>% summarize(
      prop_pos_diff = mean(prop_pos_diff, na.rm = T)
    ) %>% {
    ggplot(., aes(mfi_skew_2, prop_pos_diff)) +
        geom_jitter(aes(colour = Approach), size = 3) +
        scale_colour_brewer(palette = "Dark2") +
        xlab("Seropositive Distribution Skew") +
        ylab("Estimated - True Seropositivity") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_smooth(aes(colour = Approach), se = F) 
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
  ))%>% {
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
      prop_pos_diff = mid - true
    , CI_wid = upr - lwr
  ) %>% filter(prop_pos_true < 0.05, prop_overlap_quant < 0.1) %>%
  group_by(mod.meth, log_mfi) %>%
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
  )) %>% filter(prop_pos < 0.05, prop_overlap_quant < 0.1) %>% 
  group_by(mod.meth, log_mfi, name) %>%
  summarize(
     tot_cov = mean(cover, na.rm = T)
    , m_diff  = mean(m_diff, na.rm = T)
  )

all_compare1 %>% {
  ggplot(., aes(mod.meth, tot_cov)) + 
    geom_point(aes(colour = mod.meth), size = 3) +
    facet_wrap(~log_mfi) +
    xlab("") +
    ylab("Estimated - True Coefficient Estimate") +
    theme(axis.text.x = element_blank())
}

all_compare1 %>% {
  ggplot(., aes(mod.meth, m_diff)) + 
    geom_point(aes(colour = mod.meth), size = 3) +
    facet_wrap(~log_mfi) +
    xlab("") +
    ylab("Estimated - True Coefficient Estimate") +
    theme(axis.text.x = element_blank())
}

all_compare2 %>% {
    ggplot(., aes(mod.meth, tot_cov)) + 
      geom_point(aes(colour = mod.meth), size = 3) +
      facet_grid(log_mfi~name) +
      xlab("") +
      ylab("Estimated - True Coefficient Estimate") +
      theme(axis.text.x = element_blank())
}

all_compare2 %>% {
  ggplot(., aes(mod.meth, m_diff)) + 
    geom_point(aes(colour = mod.meth), size = 3) +
    facet_grid(log_mfi~name) +
    xlab("") +
    ylab("Estimated - True Coefficient Estimate") +
    theme(axis.text.x = element_blank())
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
      ylab("Estimate - True Proportion Positive") +
      theme(
        axis.text.x = element_text(size = 12, angle = 300, hjust = 0)
      )
}

##### 
## Supplement Figure 4: More linkages between individual and population level seropositivity estimates
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
  filter(log_mfi == "mfi") %>% dplyr::select(
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
      ) +
      scale_x_log10() +
      scale_y_log10()
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



