## needed packages
needed_packages <- c(
  "tidyverse"
, "magrittr"
, "ggplot2"
, "ggridges"
)

lapply(needed_packages, require, character.only = TRUE) %>% unlist()

## ggplot theme
theme_set(theme_bw())
suppressWarnings(
  theme_update(
    axis.text.x      = element_text(size = 10)
  , axis.text.y      = element_text(size = 10)
  , axis.title.x     = element_text(size = 16)
  , axis.title.y     = element_text(size = 16)
  , legend.title     = element_text(size = 12)
  , panel.grid.major = element_blank()
  , panel.grid.minor = element_blank()
  , strip.background = element_blank()
  , panel.margin     = unit(0, "lines")
  , legend.key.size  = unit(.55, "cm")
  , legend.key       = element_rect(fill = "white")
  , panel.margin.y   = unit(0.5, "lines")
  , panel.border     = element_rect(colour = "black", fill = NA, size = 1)
  , strip.text.x     = element_text(size = 16, colour = "black", face = "bold"))
)

## plotting function
mfi_plot <- function(
    ## What viruses get grouped and plotted (default being the three Nipah, but can define w/e you want)
  viruses    = unique(small_bat$Virus)[str_detect(unique(small_bat$Virus), "^Ni")]
    ## Define the "y-axis" what to compare vertically (location or species being sensible options I think)
, group_BY1  = "location"
    ## Define color, i.e., another grouping variable (location or species being sensible options I think)
, group_BY2  = "PlateNo"
, colour_BY  = "species"
    ## May or may not want to facet. Could help to form a grid of viruses (not implemented yet)
, facet_BY   = NA
    ## Plot log MFI on x or raw
, logMFI     = TRUE
    ## Include the R. viruses? (not sure what they are, but sample size is rather low)
, include.R  = FALSE
) {
  
  adj_data <- small_bat %>% filter(
    Virus %in% viruses
   ) %>% mutate(full_groups = 
                  {
                    if(!is.na(group_BY2)) {
                     interaction(get(group_BY1), get(group_BY2), get(colour_BY))
                    } else {
                     interaction(get(group_BY1), get(colour_BY)) 
                    }
                  }
                )
  
  if (!include.R) {
    adj_data %<>% filter(!str_detect(Virus, "^R."))
  }
    
 gg1 <- adj_data %>% {
   ggplot(., aes(
     x     = MFI
   , y     = get(group_BY1)
   , fill  = get(colour_BY)
   , group = full_groups)
   ) +
    geom_density_ridges2(
      jittered_points = TRUE
    , point_alpha = 0.6
    , alpha = 0.5
    , point_size = 1.5
    , size = .75
    , quantile_lines = TRUE
    , quantiles = 3
    ) +
    scale_fill_brewer(
      palette = "Dark2"
    , name = colour_BY
    ) +
    ylab(group_BY1) +
    ggtitle(paste(viruses, collapse = "-")) +
    theme(
      axis.text.y = element_text(size = 10)
    , axis.text.x = element_text(size = 10)
    )
 }
 
   if (!is.na(facet_BY)) {
    gg1 <- gg1 + facet_wrap(~get(facet_BY))
   }
 
   if (logMFI) {
     gg1 <- gg1 + scale_x_log10(breaks = c(1E0, 1E1, 1E2, 1E3, 1E4, 3E4))
   }
 
  gg1
 
}

## Load and some quick reorganization
small_bat <- read.csv("../data/small_bat.csv", comment.char="#") %>% 
  filter(!is.na(SLNO)) %>% 
  dplyr::select(-SLNO) %>% 
  pivot_longer(
      .
    , -c("merge_identity", "date_sampled", "species", "location", "Sex", "Age", "BCS", "PlateNo")
    , names_to = "Virus"
    , values_to = "MFI"
  ) %>% 
  mutate(
    location = plyr::mapvalues(
      location
      ## Clean up some names
    , from = c(
        " Sadar,Faridpur", "Mogalhat,Lalmonirhat", "Sreemangal", "  sadar, Khagrachari", "Gangachara, Rangpur"      
      , "Kakoni,Lalmonirhat", "Customs house,Lalmonirhat", " kaijuri, Rajbari", "Piljong,Bagerhat", "Abhaynagar, Jessore"      
      , "Rajbari", " kanaipur, Faripur")
    , to = c(
        "Sadar, Faridpur", "Mogalhat, Lalmonirhat", "Sreemangal", "Sadar, Khagrachari", "Gangachara, Rangpur"      
      , "Kakoni, Lalmonirhat", "Customs house, Lalmonirhat", "Kaijuri, Rajbari", "Piljong, Bagerhat", "Abhaynagar, Jessore"      
      , "Rajbari", "Kanaipur, Faripur")
    )
    ## Join redundant male sex labels
   , Sex = plyr::mapvalues(Sex, from = c("male ", "Male "), to = c("Male", "Male"))
  ) %>% 
  ## For consistency
  rename(sex = Sex, age = Age, bcs = BCS)

mfi_plot(
  viruses   = unique(small_bat$Virus)
, group_BY1  = "location"
, group_BY2  = NA
, facet_BY  = "Virus"
, logMFI    = TRUE
, include.R = FALSE
)

mfi_plot(
  viruses    = "MenVN"
, group_BY1  = "location"
, group_BY2  = "PlateNo"
, colour_BY  = "species"
, facet_BY   = NA
, logMFI     = FALSE
, include.R  = FALSE
)

#########################
### ---- Fitting ---- ###

## Potential covariates:
 ## Location
 ## Species
 ## Sex
 ## Age
 ## BCS (body condition?)
 ## Plate number

small_bat.t <- small_bat %>% 
  ## For now, just filos
  #filter(Virus %in% c("MarVGP", "EboVGP")) %>%
  #filter(Virus %in% "MarVGP") %>%
  filter(Virus %in% "EboVGP") %>%
  mutate(
    location_n = as.numeric(as.factor(location))
  , spec_n     = as.numeric(as.factor(species))
  ) %>% 
  mutate(
    entry = seq(n()), .before = 1
  ) %>% mutate(
    MFI = log10(MFI)
  ) %>% droplevels()

cov_mat <- model.matrix(~species+sex, small_bat.t)[, ]

conflicted::conflicts_prefer(rstan::lookup)

stan_fit   <- stan(
  file  = "cluster_regression.stan"
 , data = list(
      N       = nrow(small_bat.t)
    , N_plate = max(small_bat.t$PlateNo)
    , plate   = small_bat.t$PlateNo
    , N_loc   = max(small_bat.t$location_n)
    , loc     = small_bat.t$location_n
    , N_spec  = max(small_bat.t$spec_n)
    , spec    = small_bat.t$spec_n
    , y       = small_bat.t$MFI
    , mm      = cov_mat
    , N_cov   = ncol(cov_mat)
  )
 , iter    = 2000#6000            
 , warmup  = 500#2000
 , thin    = 1
 , chains  = 3
 , cores   = 3
 , seed    = 10001
 , refresh = 200
 , control = list(adapt_delta = 0.92, max_treedepth = 13)
)

## Extract samples
stan.fit.samples <- stan_fit %>% extract()

membership_p.mar  <- stan.fit.samples$membership_p
membership_p.ebo  <- stan.fit.samples$membership_p

small_bat.t.mar   <- small_bat.t
small_bat.t.ebo   <- small_bat.t

## Convert pop positive count to percentage
#tot_n <- small_bat.t %>% group_by(location_n, spec_n) %>% summarize(n_entry = n()) %>% ungroup()
#for (i in 1:nrow(tot_n)) {
#  stan.fit.samples$pop_sero[, tot_n$location_n[i], tot_n$spec_n[i]] <- stan.fit.samples$pop_sero[, tot_n$location_n[i], tot_n$spec_n[i]] / tot_n$n_entry[i]
#}

pop_sero.mar <- membership_p.mar %>% 
  reshape2::melt() %>% 
  filter(Var2 == 1) %>%
  dplyr::select(-Var2) %>%
  mutate(assign_pos = ifelse(value > 0.975, 1, 0)) %>% 
  rename(entry = Var3) %>%
  left_join(., small_bat.t, by = "entry") %>%
  group_by(species, location, iterations) %>%
  summarize(pos_perc = sum(assign_pos) / n()) %>%
  ungroup() %>%
  group_by(species, location) %>%
  summarize(
      prob_inf_lwr   = quantile(pos_perc, 0.025)
    , prob_inf_lwr_n = quantile(pos_perc, 0.200)
    , prob_inf_mid   = quantile(pos_perc, 0.500)
    , prob_inf_upr_n = quantile(pos_perc, 0.800)
    , prob_inf_upr   = quantile(pos_perc, 0.975)
  ) %>% mutate(
   Virus = "MarVGP"
  # Virus = "EboVGP"
   , .before = 1
  )

pop_sero.b <- rbind(
  pop_sero.mar
, pop_sero.ebo
)

prob_inf <- stan.fit.samples$membership_p %>% 
  reshape2::melt() %>% 
  filter(Var2 == 1) %>%
  group_by(Var3) %>%
  summarize(
    prob_inf_lwr   = quantile(value, 0.025)
  , prob_inf_lwr_n = quantile(value, 0.200)
  , prob_inf_mid   = quantile(value, 0.500)
  , prob_inf_upr_n = quantile(value, 0.800)
  , prob_inf_upr   = quantile(value, 0.975)
  ) %>%
  rename(entry = Var3)

## Summarize stan estimates and add to data
small_bat.t %<>% 
  left_join(., prob_inf, by = "entry") %>% 
  mutate(
    assigned_positive = ifelse(prob_inf_lwr > 0.975, 1, 0) %>% as.factor()
  , .after = MFI
  ) %>% dplyr::select(
    -contains("prob")
  )

small_bat.t.MAR <- small_bat.t
small_bat.t.EBO <- small_bat.t

small_bat.t.b <- rbind(
  small_bat.t.MAR
, small_bat.t.EBO
)

small_bat.t.b %>% {
  ggplot(., aes(MFI, assigned_positive)) + geom_jitter(aes(colour = Virus))
}

## Custom plot (can update sometime soonish to use function)
scale_point_size_continuous <- function(name = ggplot2::waiver(), breaks = ggplot2::waiver(), labels = ggplot2::waiver(),
                                        limits = NULL, range = c(1, 6),
                                        trans = "identity", guide = "legend", aesthetics = "point_size") {
  ggplot2::continuous_scale(aesthetics, "area", scales::area_pal(range), name = name,
                            breaks = breaks, labels = labels, limits = limits, trans = trans,
                            guide = guide)
}
scale_point_shape <- function(..., solid = TRUE, aesthetics = "point_shape") {
  discrete_scale(aesthetics, "shape_d", scales::shape_pal(solid), ...)
}

small_bat.t.b[small_bat.t.b$entry == 350, ]$assigned_positive <- as.factor(c(0, 0))

gg1 <- small_bat.t.b %>% 
  mutate(ps = as.numeric(assigned_positive)) %>%
  rename(Seropostatus= assigned_positive) %>%
  mutate(Seropostatus = plyr::mapvalues(
    Seropostatus
  , from = c(0, 1)
  , to   = c("Negative", "Positive")
  )) %>% {
  ggplot(., aes(
      x            = MFI
    , y            = location
    , fill         = species
    , group        = interaction(location, species)
    , point_shape  = Seropostatus
    , point_size   = ps
    )
  ) +
    geom_density_ridges(
        jittered_points = TRUE
      , point_alpha = 0.6
      , alpha = 0.5
    #  , point_size = 1.75
      , size = .75
      , quantile_lines = FALSE
      , quantiles = 3
    ) +
    scale_fill_brewer(
        palette = "Dark2"
      , name = "Species"
    ) +
    scale_point_size_continuous(
      range = c(1, 2)
    , guide = "none"
    ) +
    ylab("Location") +
    xlab("Log10(MFI)") +
    theme(
        axis.text.y = element_text(size = 10)
      , axis.text.x = element_text(size = 10)
    ) + facet_wrap(~Virus, scales = "free") +
    scale_x_log10()
}

gg1

small_bat.t.b %>%
  group_by(location, species) %>%
  summarize(
    tot_pos = (as.numeric(assigned_positive) - 1) %>% sum()
  , tot_n   = n()
  ) %>% write.csv(
    "positives.csv"
  )

pop_sero.b %>%
  filter(prob_inf_lwr > 0) %>%
  mutate(y_ax = interaction(species, location, sep = " -- ")) %>% {
  ggplot(., aes(prob_inf_mid, location)) +
    geom_errorbarh(aes(xmin = prob_inf_lwr, xmax = prob_inf_upr, colour = species)
                   , height = 0.3, linewidth = 0.75
                   , position = position_dodge(0.5)) +
    geom_errorbarh(aes(xmin = prob_inf_lwr_n, xmax = prob_inf_upr_n, colour = species)
                   , height = 0, linewidth = 1.5
                   , position = position_dodge(0.5)) +
    geom_point(aes(colour = species), size = 2
               , position = position_dodge(0.5)) +
    scale_colour_brewer(palette = "Dark2") +
    xlab("Population Seropositivity") +
    ylab("Location") +
    theme(panel.spacing = unit(1, "lines")) +
    facet_wrap(~Virus)
}


