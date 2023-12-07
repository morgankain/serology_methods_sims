library(ggridges)

lwr <- replicate(3, rnorm(100, 3000, 2000))
lwr <- ifelse(lwr < 0, lwr + 1, lwr)
upr <- replicate(3, rnorm(5, 12000, 2000))
out <- rbind(lwr, upr) %>% reshape2::melt(.)

gg.1 <- out %>% {
  ggplot(., aes(
      x     = value
    , y     = Var2 %>% as.factor()
  )) +
geom_density_ridges2(
    jittered_points = TRUE
  , point_alpha = 0.6
  , alpha = 0.5
  , point_size = 1.5
  , size = .75
  , quantile_lines = TRUE
  , quantiles = 3
) + xlab("MFI") +
    ylab("Simulation")
}

gg.2 <- out %>% {
  ggplot(., aes(
      x     = value %>% log()
    , y     = Var2 %>% as.factor()
  )) +
    geom_density_ridges2(
      jittered_points = TRUE
      , point_alpha = 0.6
      , alpha = 0.5
      , point_size = 1.5
      , size = .75
      , quantile_lines = TRUE
      , quantiles = 3
    ) + xlab("MFI") +
    ylab("Simulation")
}

gridExtra::grid.arrange(gg.1, gg.2, ncol = 2)
  
test_conc <- data.frame(
  n_anti = exp(seq(1, 14, by = 0.1))
, flore   = 0
  )

n_beads <- 100

for (i in 1:nrow(test_conc)) {
  
  bead_attachment <- data.frame(
    bead = seq(n_beads)
  , anti = 0
  )
  
  bead_attach <- data.frame(bead_attach = sample(bead_attachment$bead, test_conc$n_anti[i] %>% round(), replace = T))
  bead_attach %<>% group_by(bead_attach) %>% summarize(n_attach = n()) %>% arrange(desc(n_attach))
  
  bead_attachment[bead_attach$bead_attach, ]$anti <- bead_attach$n_attach
  
  bead_attachment %<>% mutate(log_anti = ifelse(anti > 0, log(anti), 0))
  
  test_conc$flore[i] <- mean(bead_attachment$log_anti)
  
}

test_conc %>% {
  ggplot(., aes(n_anti, flore)) + geom_point() +
    scale_x_log10() + scale_y_log10()
}




