## Explore fits for the regression coefficients 
plot_summary            <- function(coef_ests, param_sets, coverage) {
  
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
plot_group_assignments  <- function(group_assignment) {
  
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
plot_pop_seropos        <- function(pop_seropositivity) {
  
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

