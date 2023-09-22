#### 
## Project Notes
####

### ----- To Do starting Sep 25
 ## Explore sims from over the weekend
 ## Work on better summaries across parameter sets
 ## Strategies to fit to real data better

 ## Start on more complicated regressions for beta
  ## -- Framework that builds up to random / time varying in lambda
    ## -- more complicated random effects
    ## -- time varying predictors??
    ## -- nonlinear relationship between titer and MFI

### ----- Things to think about
  ## -- Transform to log titer first or after logistic[?]
   ## -- !! NOTE: shift to actually simulate titer then transform to MFI
    ## REVISIT: sigmoidal on raw vs log scale?
      ## logistic transform built into bioplex software
        ## https://www.bio-rad.com/webroot/web/pdf/lsr/literature/10022815.pdf

### ----- Initial sims
 ## 1) Simple
 ## 2) Consider low positive numbers
 ## 3) Much more complicated regression on betas

### ----- Focus
 ## 1) coverage vs bias for beta_age
 ## 2) feasibility of using the "second best" options: When can you collapse and not get hosed?

### ----- Conceptualization
 ## 1) Conceptual model is that there are two clusters and not getting two clusters
  ## is due to non-linearity and other calibration stuff

