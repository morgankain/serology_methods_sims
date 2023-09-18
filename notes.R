#### 
## Project Notes
####

### ----- To Do -short term- (as of Sep 18, based on meeting with Noam)
 ## 1) Rename all of my covariates 
 ## 2) Add functionality to repeat each sim n times and calculate bias / coverage across the n sims
 ## 3) Run simple sim and report bias

### ----- To Do -longer term- 

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
 ## 1) coverage vs bias for what I currently call lambda_age
 ## 2) feasibility of using the "second best" options: When can you collapse and not get hosed?

### ----- Conceptualization
 ## 1) Conceptual model is that there are two clusters and not getting two clusters
  ## is due to non-linearity and other calibration stuff



