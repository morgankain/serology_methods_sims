dat_for_dist <- mclust.g %>% mutate(
  group   = as.factor(group)
  , sim_num = as.factor(sim_num)
) %>% filter(
  method == "constrained_mclust"
  , log_mfi == "mfi"
)

dat_for_points_sd <- three_sd.g %>% 
  mutate(
    gp      = group - V2
    , sim_num = as.factor(sim_num)
  ) %>%
  filter(
    log_mfi == "mfi"
    , sd_method == "assigned_group_sample_mfi"
  ) 

sd_cutoffs_sd <- data.frame(
  param_set = rand_ps
  , lwr = dat_for_points_sd %>% 
    group_by(param_set) %>% 
    filter(V1 == 1) %>%
    filter(mfi == max(mfi)) %>%
    pull(mfi)
  , upr = dat_for_points_sd %>% 
    group_by(param_set) %>% 
    filter(V2 == 1) %>%
    filter(mfi == min(mfi)) %>%
    pull(mfi)
) %>% 
  mutate(cutoff = (lwr + upr) / 2)

dat_for_points_sd %<>% 
  filter(gp != 0) %>%
  mutate(
    gp = ifelse(
      gp > 0
    , gp / 1000
    , gp / 1000 #-0.0001
    ) 
  ) %>%
  mutate(group = as.factor(group))
  
gg.1 <- dat_for_dist %>% {
   ggplot(., aes(x = mfi
                 #, y = after_stat(count / sum(count))
          )) + 
    geom_vline(
        data = sd_cutoffs_sd
      , aes(xintercept = cutoff)
      , linetype = "dashed"
    ) +
#    geom_violinh(
#      data = dat_for_points_sd
#      , aes(x = mfi, y = gp
#            , group = gp
#      ) 
#    ) +
    geom_jitter(data = dat_for_points_sd 
              , aes(y = gp, colour = group)
              , height = 0.0001
    ) +
    geom_density(aes(fill = group, colour = group
                  , group = interaction(sim_num, group))
                 , alpha = 0.3
                   ) +  
      scale_fill_brewer(
        palette = "Dark2"
        , name = "True
Serostatus"
      , labels = c(
        "Negative"
      , "Positive"
      )
      ) +
      scale_colour_brewer(
        palette = "Dark2"
        , name = "True
Serostatus"
        , labels = c(
          "Negative"
        , "Positive"
        )
        ) +
    facet_wrap(~param_set) +
      xlab("MFI") +
      theme(
          axis.text.x = element_text(
            size = 10
          , angle = 300
          , hjust = 0)
        , axis.text.y = element_text(size = 10)
        , axis.title.y = element_text(size = 12)
        , panel.spacing.x = unit(1, "lines")
      ) +
    scale_y_continuous(
      name = "P(False Positive)     |  P(False Negative)"
      , breaks = c(
        -0.001, 0, 0.001
      )
      , labels = c(
        "1"
, ""
, "1"
      )
    , sec.axis = dup_axis(
      name = "Proportion of Population"
      , breaks = c(
        -.001, -0.0005, 0, 0.0005, 0.001
      )
      , labels = c("", "", "0", "0.05", "0.1")
    )
    )
}

dat_for_points_mclust <- mclust.g %>% 
  mutate(
      gp      = group - V2
    , sim_num = as.factor(sim_num)
    ) %>%
  filter(
      log_mfi == "mfi"
    , method == "constrained_mclust"
  ) 

dat_for_points_mclust %<>% 
  filter(
    !(abs(gp) < 0.02)
  ) %>%
  mutate(gp = gp / 1000) %>%
  mutate(group = as.factor(group))

gg.2 <- dat_for_dist %>% {
  ggplot(., aes(x = mfi
                #, y = after_stat(count / sum(count))
  )) + 
    geom_vline(
        data = sd_cutoffs_sd
      , aes(xintercept = cutoff)
      , linetype = "dashed"
    ) +
    geom_point(data = dat_for_points_mclust
                , aes(y = gp, colour = group)
    ) +
    geom_density(aes(fill = group, colour = group
                     , group = interaction(sim_num, group))
                 , alpha = 0.3
    ) +  
    scale_fill_brewer(
      palette = "Dark2"
      , name = "True
Serostatus"
      , labels = c(
        "Negative"
        , "Positive"
      )
    ) +
    scale_colour_brewer(
      palette = "Dark2"
      , name = "True
Serostatus"
      , labels = c(
        "Negative"
        , "Positive"
      )
    ) +
    facet_wrap(~param_set) +
    xlab("MFI") +
    theme(
      axis.text.x = element_text(
          size = 10
        , angle = 300
        , hjust = 0)
      , axis.text.y = element_text(size = 10)
      , axis.title.y = element_text(size = 12)
      , panel.spacing.x = unit(1, "lines")
    ) +
    scale_y_continuous(
      name = "P(False Positive)     |  P(False Negative)"
      , breaks = c(
        -.001, -0.0005, 0, 0.0005, 0.001
      )
      , labels = c(
        "1", "0.5", "0", "0.5", "1"
      )
, sec.axis = dup_axis(
  name = "Proportion of Population"
  , breaks = c(
    -.001, -0.0005, 0, 0.0005, 0.001
  )
  , labels = c("", "", "0", "0.05", "0.1")
)
    )
}

dat_for_points_stan <- stan.g %>% 
  mutate(
    V1 = 1 - mid
  , V2 = mid
  , gp = group - V2
  , sim_num = as.factor(sim_num)
  ) %>%
  filter(
    log_mfi == "mfi"
  ) 

dat_for_points_stan %<>% 
  filter(
    !(abs(gp) < 0.01)
  ) %>%
  mutate(gp = gp / 1000) %>%
  mutate(group = as.factor(group))

gg.3 <- dat_for_dist %>% {
  ggplot(., aes(x = mfi
                #, y = after_stat(count / sum(count))
  )) + 
    geom_vline(
      data = sd_cutoffs_sd
      , aes(xintercept = cutoff)
      , linetype = "dashed"
    ) +
    geom_point(data = dat_for_points_stan 
               , aes(y = gp, colour = group)
    ) +
    geom_density(aes(fill = group, colour = group
                     , group = interaction(sim_num, group))
                 , alpha = 0.3
    ) +  
    scale_fill_brewer(
      palette = "Dark2"
      , name = "True
Serostatus"
      , labels = c(
        "Negative"
        , "Positive"
      )
    ) +
    scale_colour_brewer(
      palette = "Dark2"
      , name = "True
Serostatus"
      , labels = c(
        "Negative"
        , "Positive"
      )
    ) +
    facet_wrap(~param_set) +
    xlab("MFI") +
    theme(
      axis.text.x = element_text(
        size = 10
        , angle = 300
        , hjust = 0)
      , axis.text.y = element_text(size = 10)
      , axis.title.y = element_text(size = 12)
      , panel.spacing.x = unit(1, "lines")
    ) +
    scale_y_continuous(
      name = "P(False Positive)     |  P(False Negative)"
      , breaks = c(
        -.001, -0.0005, 0, 0.0005, 0.001
      )
      , labels = c(
        "1", "0.5", "0", "0.5", "1"
      )
      , sec.axis = dup_axis(
        name = "Proportion of Population"
        , breaks = c(
          -.001, -0.0005, 0, 0.0005, 0.001
        )
        , labels = c("", "", "0", "0.05", "0.1")
      )
    )
}

gridExtra::grid.arrange(
  gg.1 + theme(
    axis.text.x = element_blank()
  , axis.title.x = element_blank()
  , axis.ticks.length.x = unit(0, "cm")
  , strip.text.x = element_blank()
  , plot.margin = unit(c(0.1,0.1,0.1,0.4), "cm")
  , axis.title.y = element_text(
    size = 10, vjust = 5 
    )
  )
, gg.2 + theme(
  axis.text.x = element_blank()
  , axis.title.x = element_blank()
  , axis.ticks.length.x = unit(0, "cm")
  , strip.text.x = element_blank()
  , plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm") 
  , axis.title.y = element_text(size = 10)
  )
, gg.3 + theme(
    strip.text.x = element_blank()
 # , axis.text.x = element_blank()
#  , axis.title.x = element_blank()
  , plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
  , axis.title.y = element_text(size = 10)
)
, ncol = 1
)

mclust_bias_summary <- dat_for_points_mclust %>% 
  group_by(group, param_set) %>%
  summarize(
    mean_bias = mean(gp)
  ) %>% mutate(
    method = "mclust"
  , .before = 1
  )
stan_bias_summary <- dat_for_points_stan %>% 
  group_by(group, param_set) %>%
  summarize(
    mean_bias = mean(gp)
  ) %>% mutate(
      method = "stan"
    , .before = 1
  )
bias_summary <- rbind(
  mclust_bias_summary
, stan_bias_summary
)

bias_summary %>% mutate(
  mean_bias = abs(mean_bias)
) %>% 
  pivot_wider(
    names_from = method
  , values_from = mean_bias
  ) %>%
  mutate(
    param_set = as.factor(param_set)
  , group = as.factor(group)
  ) %>% {
    ggplot(., aes(mclust, stan)) +
      geom_point(aes(
        colour = param_set
      , shape  = group
      ), size = 3) +
      geom_abline(
        intercept = 0, slope = 1
      )
  }

