##### 3sd method and log vs not log on pop seropositivity ----- 

all.out$pop_seropositivity %>% 
  filter(model == "3sd") %>%
  dplyr::select(log_mfi, param_set, method, true, prop_pos_diff) %>% 
  pivot_wider(names_from = log_mfi, values_from = prop_pos_diff) %>% 
  mutate(method = plyr::mapvalues(method, from = c(
    "assigned_group_control", "assigned_group_robust"
  ), to = c(
    "Approximated Control", "Robust Mean and SD"
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
            "Approximated Control", "Robust Mean and SD"
          ))) %>%
          filter(param_set %in% c(126, 190))
        , aes(mfi, log_mfi, colour = method), size = 5, shape = c(0, 0, 2, 2)
      ) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_hline(yintercept = 0, linetype = "dotted") +
      xlab("MFI") +
      ylab("Log(MFI)") +
      scale_colour_brewer(palette = "Dark2", name = "Method")
  }

## 126, 190

three_sd.groups %>% 
  filter(param_set %in% c(126, 190)) %>% 
  mutate(param_set = as.factor(param_set)) %>%
  mutate(sd_method = plyr::mapvalues(sd_method, from = c(
    "assigned_group_control", "assigned_group_robust"
  ), to = c(
    "Approximated Control", "Robust Mean and SD"
  ))) %>%
  mutate(group = plyr::mapvalues(group, from = c(0, 1), to = c("Negative", "Positive"))) %>%
  mutate(group = as.factor(group)) %>% {
    ggplot(., aes(mfi, assigned_group)) + 
      geom_jitter(
        aes(colour = group, shape = param_set)  
        , width = 0.2, height = 0.2) +
      facet_grid(sd_method~log_mfi, scales = "free_x") +
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
      , labels = c("", "")
      , name = "Parameter
Set"
      ) +
      theme(
        axis.text.y = element_text(size = 12)
      , strip.text.x = element_blank())
  }

three_sd.groups %>% 
  filter(param_set %in% c(126, 190)) %>% 
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
      facet_wrap(log_mfi~param_set, scales = "free", nrow = 3) +
      ylab("Density") +
      xlab("MFI") +
      theme(strip.text.x = element_blank())
  }

all.out$pop_seropositivity %>% 
  filter(model == "3sd") %>%
  filter(prop_pos_true < 0.1) %>% 
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



