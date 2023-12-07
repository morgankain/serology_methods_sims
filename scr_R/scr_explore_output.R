####################################################################################
## Exploring output, but also starting to construct sample figures for manuscript ##
####################################################################################

########################
## Group assignment: Misclassification error
## ------ tar_make(group_assignment); tar_load(group_assignment)
########################

####
## Exploring 3SD 
####

ggpairs(
  group_assignment %>% filter(model == "3sd") %>% ungroup() %>% dplyr::select(
    group, misclass_error_p, beta_base, mu_neg, sd_neg, mu_pos_delta, sd_pos_delta
  ) %>% mutate(group = as.factor(group)) %>% filter(misclass_error_p > 0.9)
  , ggplot2::aes(colour = group)
)

ggpairs(
  group_assignment %>% filter(model == "3sd") %>% ungroup() %>% dplyr::select(
    group, misclass_error_p, beta_base, mu_neg, sd_neg, mu_pos_delta, sd_pos_delta
  ) %>% mutate(group = as.factor(group)) %>% filter(group == 1)
  , ggplot2::aes(colour = group)
)

group_assignment %>% mutate(group = as.factor(group)) %>% 
  filter(model == "3sd") %>% {
    ggplot(., aes(mu_pos_delta, misclass_error_p)) + 
      geom_jitter(aes(colour = group, size = sd_neg)) +
      scale_colour_brewer(palette = "Dark2")
  }

group_assignment %>% mutate(group = as.factor(group)) %>% 
  filter(model == "3sd") %>% {
    ggplot(., aes(sd_neg, misclass_error_p)) + 
      geom_point(aes(colour = group, size = mu_pos_delta)) +
      scale_colour_brewer(palette = "Dark2")
  }

group_assignment %>% mutate(group = as.factor(group)) %>% 
  filter(model == "3sd") %>%
  arrange(desc(misclass_error_p))

head(group_assignment)

group_assignment %>% mutate(group = as.factor(group)) %>% 
  filter(model == "3sd") %>%
  filter(misclass_error_p > 0.8) %>% as.data.frame()


####
## Exploring mclust
####

group_assignment %>% 
  mutate(group = as.factor(group)) %>% 
  filter(model == "mclust") %>%
  group_by(group) %>% summarize(
    mean(misclass_error_p)
  )
  
ggpairs(
  group_assignment %>% filter(model == "mclust") %>% ungroup() %>% dplyr::select(
    group, misclass_error_p, beta_base, mu_neg, sd_neg, mu_pos_delta, sd_pos_delta
  ) %>% mutate(group = as.factor(group)) 
  , ggplot2::aes(colour = group)
)

group_assignment %>% mutate(group = as.factor(group)) %>% 
  filter(model == "mclust") %>% {
    ggplot(., aes(mu_pos_delta, misclass_error_p)) + 
      geom_jitter(aes(colour = group, size = sd_neg)) +
      scale_colour_brewer(palette = "Dark2")
  }

group_assignment %>% mutate(group = as.factor(group)) %>% 
  filter(model == "mclust") %>% {
    ggplot(., aes(sd_neg, misclass_error_p)) + 
      geom_point(aes(colour = group, size = mu_pos_delta)) +
      scale_colour_brewer(palette = "Dark2")
  }


####
## Exploring stan
####

group_assignment %>% 
  mutate(group = as.factor(group)) %>% 
  filter(model == "cluster_regression_base_1.stan") %>%
  group_by(group) %>% summarize(
    mean(misclass_error_p)
  )

ggpairs(
  group_assignment %>% filter(
    model == "cluster_regression_base_1.stan"
  , quantile == "mid") %>% ungroup() %>% dplyr::select(
    group, misclass_error_p, beta_base, mu_neg, sd_neg, mu_pos_delta, sd_pos_delta
  ) %>% mutate(group = as.factor(group)) 
  , ggplot2::aes(colour = group)
)

####
## Among model fit types
####

group_assignment %>% filter(model != "cluster_regression_base_1.stan") %>% 
  mutate(group = as.factor(group)) %>%
  dplyr::select(-c(sim_num, quantile, prob, n_samps)) %>% pivot_wider(
    values_from = misclass_error_p, names_from = model
  ) %>% mutate(
    group = plyr::mapvalues(
      group
    , from = c(1, 0)
    , to   = c("True Positive
Assigned Negative", "True Negative
Assigned Positive")
      )
  ) %>% {
   ggplot(., aes(`3sd`, mclust)) + geom_point(
     aes(colour = group)
   ) + 
      scale_colour_brewer(palette = "Dark2", name = "Seropositivity
Assessment") +
      geom_abline(intercept = 0, slope = 1) +
      xlab("Method: 3sd -- Proportion of Individuals") +
      ylab("Method: Clustering (mclust) -- Average Probability") +
      theme(
         legend.key.size = unit(.85, "cm")
      )
  }


########################
## Population-level seropositivity
## ------ tar_make(pop_seropositivity); tar_load(pop_seropositivity)
########################

ggpairs(
  pop_seropositivity %>% filter(model == "mclust") %>% ungroup() %>% dplyr::select(
    prop_pos_diff, beta_base, mu_neg, sd_neg, mu_pos_delta, sd_pos_delta
  ) 
)

ggpairs(
  pop_seropositivity %>% filter(model == "cluster_regression_base_1.stan") %>% filter(quantile == "mid") %>% ungroup() %>% dplyr::select(
    prop_pos_diff, beta_base, mu_neg, sd_neg, mu_pos_delta, sd_pos_delta
  ) 
)

three_vs_stan <- pop_seropositivity %>% filter(quantile == "mid", model != "mclust") %>% 
  dplyr::select(-prop_pos) %>% 
  mutate(
      beta_base = plogis(beta_base)
    , true = round(true, 3)
    ) %>%
  pivot_wider(names_from = model, values_from = prop_pos_diff) %>% 
  mutate(difff = cluster_regression_base_1.stan - `3sd`)
  
three_vs_stan %>% {
    ggplot(., aes(`3sd`, cluster_regression_base_1.stan)) + 
      geom_point(aes(size = sd_neg, colour = beta_base)) +
      scale_color_continuous(name = "True 
proportion 
positive") +
      scale_size_continuous(name = "SD of 
negative 
distribution") +
      xlab("Estimate - True; Method: 3sd") + 
      ylab("Estimate - True; Method: Bayesian normal, normal") +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      annotate(geom = "text", x = -0.42, y = 0.58, label = "A", fontface = "bold", size = 7) +
      annotate(geom = "text", x = -0.42, y = -0.13, label = "C", fontface = "bold", size = 7) +
      annotate(geom = "text", x = 0.1, y = -0.13, label = "D", fontface = "bold", size = 7) +
      annotate(geom = "text", x = 0.1, y = 0.58, label = "B", fontface = "bold", size = 7)
  }

three_vs_stan %>% {
  ggplot(., aes(mu_pos_delta, difff)) + geom_point()
}

stan.pop_sero <- pop_seropositivity %>% 
  filter(model == "cluster_regression_base_1.stan") %>%
  dplyr::select(-prop_pos_diff) %>% pivot_wider(
  names_from = quantile, values_from = prop_pos
) %>% mutate(
  cover = ifelse(lwr/1.02 < true & upr*1.02 > true, 1, 0)
, bias  = ifelse(cover == 0, ifelse(upr*1.02 < true, -1, 1), 0)
) %>% mutate(
  bias_ci = ifelse(bias == -1, upr*1.02 - true, lwr/1.02 - true)
) %>% mutate(
  bias_ci = ifelse(bias == 0, 0, bias_ci)
)

ggpairs(
  stan.pop_sero %>% filter(model == "cluster_regression_base_1.stan") %>% ungroup() %>% dplyr::select(
     lwr_bias, beta_base, mu_neg, sd_neg, mu_pos_delta, sd_pos_delta
  ) 
)

gg.1 <- stan.pop_sero %>% mutate(
    cover = as.factor(cover)
  , bias = as.factor(bias)
  ) %>% mutate(
    bias = plyr::mapvalues(
      bias
    , from = c(-1, 0, 1)
    , to   = c("Under-estimated", "Covered", "Over-estimated")
    )
  , beta_base = plogis(beta_base)
  ) %>% {
  ggplot(., aes(beta_base, bias_ci)) + 
    geom_point(aes(colour = bias, shape = cover, size = sd_neg)) +
    scale_shape_manual(values = c(16, 18), name = "Cover") +
    scale_color_brewer(palette = "Dark2", name = "Bias") +
    scale_size_continuous(name = "SD of 
negative 
distribution") +
    xlab("True proportion positive") + 
    ylab("Bias in proportion estimate from 95% CI")
}

gg.2 <- stan.pop_sero %>% mutate(
  cover = as.factor(cover)
  , bias = as.factor(bias)
) %>% mutate(
  bias = plyr::mapvalues(
    bias
    , from = c(-1, 0, 1)
    , to   = c("Under-estimated", "Covered", "Over-estimated")
  )
  , beta_base = plogis(beta_base)
) %>% {
  ggplot(., aes(mu_pos_delta, bias_ci)) + 
    geom_point(aes(colour = bias, shape = cover, size = sd_neg)) +
    scale_shape_manual(values = c(16, 18), name = "Cover") +
    scale_color_brewer(palette = "Dark2", name = "Bias") +
    scale_size_continuous(name = "SD of 
negative 
distribution") +
    xlab("Difference in mean between negative and positive group") + 
    ylab("Bias in proportion estimate from 95% CI")
}

gg.3 <- stan.pop_sero %>% mutate(
  cover = as.factor(cover)
  , bias = as.factor(bias)
) %>% mutate(
  bias = plyr::mapvalues(
    bias
    , from = c(-1, 0, 1)
    , to   = c("Under-estimated", "Covered", "Over-estimated")
  )
  , beta_base = plogis(beta_base)
) %>% {
  ggplot(., aes(sd_neg, bias_ci)) + 
    geom_point(aes(colour = bias, shape = cover, size = beta_base)) +
    scale_shape_manual(values = c(16, 18), name = "Cover") +
    scale_color_brewer(palette = "Dark2", name = "Bias") +
    scale_size_continuous(name = "True 
proportion 
positive") +
    xlab("Standard deviation (in titer) in negative distribution") + 
    ylab("Bias in proportion estimate from 95% CI")
}

gridExtra::grid.arrange(
  gg.1, gg.2, gg.3, ncol = 1
)

pop_seropositivity %>%
  filter(model == "3sd") %>% {
    ggplot(., aes(mu_pos_delta, prop_pos_diff)) + 
      geom_jitter(aes(size = sd_neg))
  }

pop_seropositivity %>% filter(model == "3sd") %>% {
  ggplot(., aes()) 
}



