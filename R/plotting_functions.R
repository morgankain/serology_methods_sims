## Explore fits for the regression coefficients 
plot_summary                   <- function(coef_ests, param_sets, coverage) {
  
  stan_all.gg <- coef_ests %>% filter(!grepl("3sd|mclust", model)) %>% {
    ggplot(., aes(mid, name)) + 
      geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0, linewidth = 1) +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, linewidth = 0.3) +
      geom_point(aes(true, name), colour = "firebrick3") +
      facet_wrap(~model)
  }
  
  all_out.gg <- coef_ests %>% filter(name %in% c("beta_base", "beta_age")) %>% {
    ggplot(., aes(mid, model)) + 
      geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0, linewidth = 1) +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, linewidth = 0.3) +
      geom_point(aes(true, model), colour = "firebrick3") +
      facet_wrap(~name)
  }
  
  all_out.theta <- coef_ests %>% left_join(., param_sets, by = c("param_set", "sim_num")) %>%
   filter(name %in% c("beta_base", "beta_age"))
  
  cov1.gg <- all_out.theta %>% {
    ggplot(., aes(mu_pos_delta, cover)) + 
      geom_jitter(height = 0.05) +
      facet_wrap(~model)
  }
  
  cov2.gg <- all_out.theta %>% 
  mutate(mu_pos_delta_r = plyr::round_any(mu_pos_delta, 0.5)) %>%
  group_by(model, name, mu_pos_delta_r) %>% 
  summarize(m_cover = mean(cover)) %>% {
    ggplot(., aes(mu_pos_delta_r, m_cover)) + 
      geom_point(size = 2) +
      geom_line() +
      facet_grid(name~model)
  }
  
  cov3.gg <- all_out.theta %>% {
    ggplot(., aes(mu_pos_delta, CI_wid)) + 
      geom_point() +
      facet_wrap(~model)
  }
  
  cov4.gg <- all_out.theta %>% {
    ggplot(., aes(mu_pos_delta, m_diff)) + 
      geom_point() +
      facet_wrap(~model)
  }
  
  cov5.gg <- coef_ests %>% filter(name %in% c("beta_base", "beta_age")) %>%
   left_join(., param_sets, by = c("param_set", "sim_num")) %>% 
   arrange(desc(mu_pos_delta)) %>%
   mutate(mu_pos_delta = factor(mu_pos_delta, levels = unique(mu_pos_delta))) %>% {
     ggplot(., aes(mid, mu_pos_delta)) + 
       geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0, linewidth = 1) +
       geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, linewidth = 0.3) +
       geom_vline(aes(xintercept = true), colour = "firebrick3", linewidth = 0.3) +
       facet_grid(model~name) +
       xlab("Estimate") +
       ylab("Difference in -mean- between positive/negative") +
       theme(axis.text.y = element_text(size = 9))
   }
  
  cov6.gg <- coverage %>%
    ungroup() %>%
    group_by(model, name) %>% 
    summarize(coverage = mean(coverage)) %>% {
    ggplot(., aes(coverage, model)) +
    geom_point() +
    facet_wrap(~name) 
  }

  return(
    list(
      stan_all.gg
    , all_out.gg
    , cov1.gg
    , cov2.gg
    , cov3.gg
    , cov4.gg
    , cov5.gg
    , cov6.gg
    )
  )
  
}

## Explore individual-level group assignments
plot_group_assignment_summary  <- function(group_assignment) {
  
group_assignment %>% 
    ungroup() %>%
    group_by(model, group, quantile) %>%
    summarize(
      lwr   = quantile(prob, 0.025)
    , lwr_n = quantile(prob, 0.200)
    , mid   = quantile(prob, 0.500)
    , upr_n = quantile(prob, 0.800)
    , upr   = quantile(prob, 0.975)
    ) %>% filter(quantile == "mid") %>% 
    mutate(group = as.factor(group)) %>% {
      ggplot(., aes(mid, group)) +
        geom_errorbar(aes(xmin = lwr_n, xmax = upr_n), linewidth = 1, width = 0) +
        geom_errorbar(aes(xmin = lwr, xmax = upr), linewidth = 0.5, width = 0.2) +
        geom_point() +
        facet_wrap(~model, ncol = 1) +
        xlab("Probability of Group 1 Affiliation") +
        ylab("True Group Affiliation")
    }
  
}

## Plot population level seropositivity
plot_pop_seropos               <- function(pop_seropositivity) {
  
    pop_seropositivity %>% dplyr::select(-prop_pos_diff) %>%
    pivot_wider(c(model, param_set, sim_num, true)
                , names_from = "quantile"
                , values_from = "prop_pos") %>% 
    mutate(
      run = interaction(param_set, sim_num)
    , sim_num = as.factor(sim_num)
    ) %>% {
    ggplot(., aes(mid, sim_num)) + 
        geom_errorbar(aes(xmin = lwr_n, xmax = upr_n, colour = sim_num), linewidth = 1, width = 0) +
        geom_errorbar(aes(xmin = lwr, xmax = upr, colour = sim_num), linewidth = 0.5, width = 0.2) +
        geom_point(aes(colour = sim_num)) +
        geom_vline(aes(xintercept = true, colour = sim_num)) +
        scale_colour_brewer(palette = "Dark2") +
        facet_grid(model~param_set) +
        theme(axis.text.x = element_text(size = 10)) +
        xlab("Estimate") + ylab("Simulation")
    }
  
}

## Plot individual-level group assignment probabilities
plot_individual_group_prob     <- function(three_sd.g, mclust.g, stan.g, num_ps = 5, num_sn = 5) {
  
if (n_distinct(three_sd.g$param_set) > num_ps) {
  rand_ps    <- sample(seq(n_distinct(three_sd.g$param_set)), num_ps)
  three_sd.g %<>% filter(param_set %in% rand_ps)
  mclust.g   %<>% filter(param_set %in% rand_ps)
  stan.g     %<>% filter(param_set %in% rand_ps)
}
  
if (n_distinct(three_sd.g$sim_num) > num_sn) {
  rand_sn    <- sample(seq(n_distinct(three_sd.g$sim_num)), num_sn)
  three_sd.g %<>% filter(sim_num %in% rand_sn)
  mclust.g   %<>% filter(sim_num %in% rand_sn)
  stan.g     %<>% filter(sim_num %in% rand_sn)
}
  
three_sd.g %<>% 
  mutate(model = "three_sd", .before = 1) %>%
  dplyr::select(-assigned_group)
  
mclust.g %<>% 
  mutate(model = "mclust", .before = 1) %>%
  dplyr::select(-assigned_group)
  
gg.1 <- mclust.g %>% mutate(
    group   = as.factor(group)
  , sim_num = as.factor(sim_num)) %>% {
     ggplot(., aes(mfi)) + 
        geom_density(aes(fill = group, colour = group
                         , group = interaction(sim_num, group)), alpha = 0.3) +
        facet_wrap(~param_set) +
        scale_fill_brewer(palette = "Dark2") +
        scale_colour_brewer(palette = "Dark2") +
     geom_jitter(data = three_sd.g %>% 
    mutate(
      gp    = group - V2
    , sim_num = as.factor(sim_num)
    ) %>% dplyr::select(-group)
    , aes(y = gp, shape = model)
    , height = 0.05) +
      ylab("True Group Assignment - Prob(pos)") +
      theme(
        axis.text.x = element_text(size = 9)
      , axis.text.y = element_text(size = 9)
      , axis.title.y = element_text(size = 10)
      )
 }
  
gg.2 <- mclust.g %>% mutate(
    group   = as.factor(group)
  , sim_num = as.factor(sim_num)) %>% {
     ggplot(., aes(mfi)) + 
        geom_density(aes(fill = group, colour = group
                         , group = interaction(sim_num, group)), alpha = 0.3) +
        facet_wrap(~param_set) +
        scale_fill_brewer(palette = "Dark2") +
        scale_colour_brewer(palette = "Dark2") +
     geom_point(data = mclust.g %>% 
    mutate(
      gp    = group - V2
    , sim_num = as.factor(sim_num)
    ), aes(y = gp, shape = model)) +
      ylab("True Group Assignment - Prob(pos)") +
      theme(
        axis.text.x = element_text(size = 9)
      , axis.text.y = element_text(size = 9)
      , axis.title.y = element_text(size = 10)
      )
 }
 
gg.0 <- stan.g %>% mutate(age = as.factor(age)) %>% {
  ggplot(., aes(mfi, mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = stan_model
                    , colour = stan_model, linetype = age
                    , group = interaction(sim_num, stan_model, age)), alpha = 0.2
                , linewidth = 0) +
    geom_line(aes(colour = stan_model, linetype = age
                , group = interaction(sim_num, stan_model, age))) +
    facet_wrap(~param_set)
  }
  
stan.g %<>% rename(model = stan_model) %>%
  mutate(V1 = 1 - mid, V2 = mid) %>%
  dplyr::select(-c(lwr, lwr_n, mid, upr_n, upr, samp))
  
gg.3 <- mclust.g %>% mutate(
    group   = as.factor(group)
  , sim_num = as.factor(sim_num)) %>% {
     ggplot(., aes(mfi)) + 
        geom_density(aes(fill = group, colour = group
                         , group = interaction(sim_num, group)), alpha = 0.3) +
        facet_wrap(~param_set) +
        scale_fill_brewer(palette = "Dark2") +
        scale_colour_brewer(palette = "Dark2") +
     geom_point(data = stan.g %>% 
    mutate(
      gp    = group - V2
    , sim_num = as.factor(sim_num)
    , model = plyr::mapvalues(model, from = c("cluster_regression_with_beta_theta.stan"
                                              , "cluster_regression_with_beta.stan")
                              , to = c("beta
theta", "beta"))
    ), aes(y = gp, shape = model)) +
      ylab("True Group Assignment - Prob(pos)") +
      theme(
        axis.text.x = element_text(size = 9)
      , axis.text.y = element_text(size = 9)
      , axis.title.y = element_text(size = 10)
      )
  }
 
return(
  list(
      gridExtra::grid.arrange(gg.1, gg.2, gg.3, ncol = 1)
    , gg.0
  )
)
  
}
