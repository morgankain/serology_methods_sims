#### 
## Project Notes
####

### ----- To Do starting Sep 20

 ## 4) Start on more complicated regressions for theta (which will be renamed beta)
  ## -- e.g., a few covariates, some continuous
  ## -- machine variance[?]
  ## -- plate random effects
  ## -- treating variation in MFI values in positive group as a random effect (e.g., by species)
  ## -- consider low lambda
  ## -- nonlinear relationship between titer and MFI
  ## --> Framework that builds up to random / time varying in lambda
   ## -- Transform to log titer first or after logistic[?]
    ## -- !! NOTE: shift to actually simulate titer then transform to MFI
      ## REVISIT: sigmoidal on raw vs log scale?
        ## logistic transform built into bioplex software
          ## https://www.bio-rad.com/webroot/web/pdf/lsr/literature/10022815.pdf

### ----- Focus
 ## 1) coverage vs bias for beta_age
 ## 2) feasibility of using the "second best" options: When can you collapse and not get hosed?

### ----- Conceptualization
 ## 1) Conceptual model is that there are two clusters and not getting two clusters
  ## is due to non-linearity and other calibration stuff

