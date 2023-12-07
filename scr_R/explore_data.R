
sera <- read.csv("data/Rabbit_sera_filo_controls.csv")

sigmoid <- function(k, x) {
  k / (1 + exp(-x))
}

samp_x <- seq(1, 100, by = 0.1)
samp_x_s <- scale(samp_x, scale = F)[, 1]

plot(data.frame(
  samp_x_s, sigmoid(100, samp_x_s))
)

sera %>%
  mutate_at(-c(1:3), as.numeric) %>% 
  pivot_longer(-c(Rep, Conc, Control)) %>% 
  filter(Control == name) %>% {
 # mutate(Conc = scale(Conc, scale = F)[, 1]) %>% 
    ggplot(., aes(Conc, value)) + 
      geom_line(aes(colour = name, group = interaction(name, Rep))) +
      scale_x_log10() +
     # scale_y_log10() +
      facet_wrap(~Control)
  }

sera %>%
  mutate_at(-c(1:3), as.numeric) %>% 
  pivot_longer(-c(Rep, Conc, Control)) %>% 
  filter(Conc < 1) %>% {
    ggplot(., aes(Conc, value)) + 
      geom_line(aes(colour = name, group = interaction(name, Rep))) +
     # scale_x_log10() +
     # scale_y_log10() +
      facet_wrap(~Control)
  }

bioplex_logistic <- function(x, a, b, c, d, g) {
  
  d + (
   (a - d) /
   (1 + (x/c)^b)^g   
      )
  
  # x is the concentration
  # y is the response
  # a is the estimated response at infinite concentration
  # b is the slope of the tangent at midpoint
  # c is the midrange concentration or midpoint
  # d is the estimated response at zero concentration
  # g is an asymmetry factor
   
}

liberia <- read.csv("data/sero_examp.csv")

liberia %>% filter(sample_number > 0) %>% 
  mutate(Plate = as.factor(Plate)) %>% {
  ggplot(., aes(x = EBOV)) + 
    geom_density(aes(fill = Plate)) +
    scale_fill_brewer(palette = "Dark2") +
      #   scale_x_log10() +
      facet_wrap(~Plate, scales = "free")
  }

liberia %>% filter(sample_number > 0) %>% 
  mutate(Plate = as.factor(Plate)) %>% {
 # mutate(EBOV = ifelse(EBOV < 1, 1, EBOV)) %>% {
  ggplot(., aes(x = EBOV)) + 
    geom_density(aes(fill = Plate, colour = NULL), alpha = 0.2) +
   # scale_x_log10() +
    scale_fill_brewer(palette = "Dark2")
  }

liberia %>% filter(sample_number > 0) %>% 
  mutate(Plate = as.factor(Plate)) %>% 
  mutate(EBOV = ifelse(EBOV < 1, 1, EBOV)) %>%
  mutate(EBOV = log(EBOV)) %>% {
  ggplot(., aes(x = EBOV)) + 
    geom_density(aes(fill = Plate, colour = NULL), alpha = 0.2) +
    scale_fill_brewer(palette = "Dark2")
  }

liberia %>% filter(sample_number > 0) %>% 
  mutate(Plate = as.factor(Plate)) %>% 
  mutate(EBOV = ifelse(EBOV < 1, 1, EBOV)) %>% 
  group_by(Plate) %>%
  mutate(EBOV = scale(EBOV)[, 1]) %>% {
  ggplot(., aes(x = EBOV)) + 
    geom_density(aes(fill = Plate, colour = NULL), alpha = 0.2) +
    scale_fill_brewer(palette = "Dark2")
  }

dat_for_stan <- liberia %>% filter(sample_number > 0) %>% 
  mutate(Plate = as.factor(Plate)) %>% 
  mutate(EBOV = ifelse(EBOV < 1, 1, EBOV)) %>% 
  group_by(Plate) %>%
  mutate(EBOV = scale(EBOV)[, 1])

stan_fit <- stan(
  file    = "stan_models/cluster_regression_base_plate_random_1.stan"
, data    = list(
   N           = nrow(dat_for_stan)
 , N_cat1r     = n_distinct(dat_for_stan$Plate)
 , y           = dat_for_stan$EBOV
 , cat1r       = dat_for_stan$Plate
 )
, pars    = c("membership_l", "ind_sero", "log_beta") 
, include = FALSE
, chains  = 4
, seed    = 483892929
, refresh = 2000
, cores   = 4
)

stan_fit <- stan(
  file    = "stan_models/cluster_regression_base_1_test.stan"
, data    = list(
   N = nrow(dat_for_stan)
 , y = dat_for_stan$EBOV
)
, pars    = c("membership_l", "ind_sero") 
, include = FALSE
, chains  = 4
, seed    = 483892929
, refresh = 2000
, cores   = 4
)

samps <- extract(stan_fit)

samps$beta %>% hist(breaks = 100)

samps$membership_p[, 1, 3] %>% hist(breaks = 100)

dat_for_stan.p <- dat_for_stan %>% ungroup() %>%
  mutate(
    lwr_prob = (samps$membership_p %>% apply(., 2:3, FUN = function(x) quantile(x, 0.025)))[1, ] %>% c()
  ) %>%
  mutate(
    pos = ifelse(lwr_prob > 0.95, 1, 0)
  )
  

dat_for_stan.p %>% summarize(
  sp = sum(assigned_positive)
, mp = sum(pos)
)

dat_for_stan.p %>% filter(pos == 1)



