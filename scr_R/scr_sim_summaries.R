#tar_make(sim.data); tar_load(sim.data)
this_set <- sample(seq(500), 3)

# sn::rsn(10000, 0, 1, -5, 0) %>% hist(breaks = 100)

sim.data.summaries %>% {
  ggplot(., aes(sd_pos_delta, prop_squish_2)) +
    geom_point() +
    facet_wrap(~log_mfi, scales = "free")
}

## figure out some summaries of these distributions
sim.data.s <- sim.data %>% 
  group_by(param_set, log_mfi) %>%
  arrange(mfi) %>%
  mutate(
    g1mfi         = ifelse(group == 1, mfi, NA) 
  , max_neg_mfi   = max(g1mfi, na.rm = T)
  , quant_neg_mfi = quantile(g1mfi, 0.975, na.rm = T)
  ) %>% mutate(
    overlap_max   = ifelse(mfi < max_neg_mfi, 1, 0)  
  , overlap_quant = ifelse(mfi < quant_neg_mfi, 1, 0) 
  )

prop_overlap <- sim.data.s %>% 
  filter(group == 2) %>% 
  summarize(
    prop_overlap_max   = sum(overlap_max) / n()
  , prop_overlap_quant = sum(overlap_quant) / n()
  )

mean_mfi <- sim.data %>% 
  group_by(param_set, group, log_mfi) %>%
  summarize(mv = mean(mfi)) %>% 
  mutate(group = as.factor(group)) %>%
  ungroup(group) %>%
  pivot_wider(., names_from = "group", values_from = "mv") %>%
  rename(mean_mfi_1 = `1`, mean_mfi_2 = `2`) %>%
  mutate(diff_mean = mean_mfi_2 - mean_mfi_1)

skew <- sim.data %>% 
  group_by(param_set, sim_num, group, log_mfi) %>%
  summarize(
      titer_var  = var(titer)
    , titer_skew = skewness(titer)
    , mfi_var    = var(mfi)
    , mfi_skew   = skewness(mfi)
  )

trunc_sev <- sim.data %>% 
  group_by(param_set, sim_num, group, log_mfi) %>%
  mutate(median_titer = median(titer)) %>%
  mutate(abov_med_tit = ifelse(titer > median_titer, 1, 0)) %>%
  group_by(param_set, sim_num, group, log_mfi, abov_med_tit) %>%
  summarize(range_mfi = max(mfi) - min(mfi)) %>%
  pivot_wider(names_from = abov_med_tit, values_from = range_mfi) %>%
  mutate(prop_squish = `1` / `0`) %>%
  ungroup() %>% dplyr::select(-c(`0`, `1`))

trunc_sev2 <- sim.data %>% 
  group_by(param_set, sim_num, group, log_mfi) %>%
  mutate(
    mmode = modeest::mlv(mfi, method = "meanshif")
  , abv_m = ifelse(mfi > mmode, 1, 0)
  )

trunc_sev2.s <- trunc_sev2 %>%
  summarize(
    p_abv_m = sum(abv_m) / n()
  ) 

trunc_sev2.s %>% ungroup() %>% {
  ggplot(., aes(x = p_abv_m)) +
    geom_density() + 
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    facet_wrap(~log_mfi, scales = "free", nrow = 3) 
}

#trunc_sev2.s %>% filter(p_abv_m > 0.75) %>% arrange(desc(p_abv_m))

#trunc_sev %>% filter(param_set == 485) %>%
#  mutate(group = as.factor(group)) %>% {
#  ggplot(., aes(x = mfi)) + 
#    geom_density(aes(colour = group, fill = group), alpha = 0.2) + 
#    scale_colour_brewer(palette = "Dark2") +
#    scale_fill_brewer(palette = "Dark2") +
#    facet_wrap(param_set~log_mfi, scales = "free", nrow = 3) +
#      geom_vline(aes(xintercept = mmode))
#}

sim.data.s <- left_join(
  prop_overlap
, mean_mfi
, by = c("param_set", "log_mfi")) %>%
  left_join(
    .
  , skew %>% ungroup() %>% 
    dplyr::select(-sim_num) %>% 
    pivot_wider(
      .
    , id_cols = c(param_set, log_mfi)
    , names_from = group
    , values_from = c(titer_var, titer_skew, mfi_var, mfi_skew))
    , by = c("param_set", "log_mfi")
  ) %>% left_join(
    .
    , trunc_sev %>%
      dplyr::select(-sim_num) %>% 
      pivot_wider(
        .
        , id_cols = c(param_set, log_mfi)
        , names_from = group
        , values_from = prop_squish) %>%
      rename(prop_squish_1 = `1`, prop_squish_2 = `2`)
    , by = c("param_set", "log_mfi")
  ) %>%
  left_join(
    .
  , trunc_sev2.s %>% ungroup() %>% 
    dplyr::select(-sim_num) %>% 
    pivot_wider(
      .
      , id_cols = c(param_set, log_mfi)
      , names_from = group
      , values_from = p_abv_m) %>%
    rename(p_abv_m_1 = `1`, p_abv_m_2 = `2`)
  , by = c("param_set", "log_mfi")
  )

sim.data.s %>% {
  ggplot(., aes(mfi_skew_2, p_abv_m_2)) +
    geom_point(aes(colour = log_mfi))
}

sim.data %>% 
  filter(param_set == 419) %>%
  filter(group == 2) %>% 
  filter(log_mfi == "mfi") %>% 
  pull(mfi) 

prop_overlap %>% arrange(desc(prop_overlap_quant))
prop_overlap %>% arrange(prop_overlap_quant)

prop_overlap %>% filter(log_mfi == "log_mfi") %>% {
  ggplot(., aes(x = prop_overlap_quant)) + 
    geom_histogram(bins = 50)
}

sim.data %>% 
  filter(param_set == 6) %>% 
  mutate(
    group = as.factor(group)
    , int = interaction(cat1f, group)
  ) %>% {
    ggplot(., aes(x = mfi)) + 
      geom_density(aes(colour = group, fill = group), alpha = 0.2) + 
      scale_colour_brewer(palette = "Dark2") +
      scale_fill_brewer(palette = "Dark2") +
      facet_wrap(param_set~log_mfi, scales = "free", nrow = 3)
  }

sim.params %>% filter(param_set == 18) %>% as.data.frame()

mean_mfi %>% {
    ggplot(., aes(x = mv)) +
      geom_density(aes(colour = group, fill = group), alpha = 0.2) + 
      scale_colour_brewer(palette = "Dark2") +
      scale_fill_brewer(palette = "Dark2")
  }

tt2 <- sim.data %>% 
  mutate(group = group - 1) %>%
  filter(log_mfi == "mfi") %>%
  group_by(param_set) %>%
  summarize(
    n_pos = sum(group)
  , ss    = n()
  ) %>% mutate(
    prop_p = n_pos / ss
  ) %>% {
    ggplot(., aes(x = prop_p)) + 
      geom_histogram(bins = 50)
  }
