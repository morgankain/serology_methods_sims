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
    axis.text.x = element_text(size = 10)
  , axis.text.y = element_text(size = 10)
  , axis.title.x = element_text(size = 16)
  , axis.title.y = element_text(size = 16)
  , legend.title = element_text(size = 12)
  , panel.grid.major = element_blank()
  , panel.grid.minor = element_blank()
  , strip.background = element_blank()
  , panel.margin = unit(0, "lines")
  , legend.key.size = unit(.55, "cm")
  , legend.key = element_rect(fill = "white")
  , panel.margin.y = unit(0.5, "lines")
  , panel.border = element_rect(colour = "black", fill = NA, size = 1)
  , strip.text.x = element_text(size = 16, colour = "black", face = "bold"))
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
small_bat <- read.csv("/data/small_bat.csv", comment.char="#") %>% filter(!is.na(SLNO)) %>% dplyr::select(-SLNO) %>% pivot_longer(
  .
, -c("merge_identity", "date_sampled", "species", "location", "Sex", "Age", "BCS", "PlateNo")
, names_to = "Virus"
, values_to = "MFI"
)

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
, logMFI     = TRUE
, include.R  = FALSE
)
