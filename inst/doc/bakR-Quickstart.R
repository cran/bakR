## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dim = c(6, 4)
)

## ----setup--------------------------------------------------------------------
library(bakR)
set.seed(123)

## -----------------------------------------------------------------------------
# Load data
data("cB_small")
data("metadf")

# Create bakRData object
bakRData <- bakRData(cB_small, metadf)


## -----------------------------------------------------------------------------
# Simulate a nucleotide recoding dataset
sim_data <- Simulate_bakRData(500)
  # This will simulate 500 features, 2 experimental conditions
  # and 3 replicates for each experimental condition
  # See ?Simulate_bakRData for details regarding tunable parameters

# Extract simulated bakRData object
bakRData <- sim_data$bakRData

# Extract simualted ground truths
sim_truth <- sim_data$sim_list

# Run the efficient model
Fit <- bakRFit(bakRData)


## ---- eval = FALSE------------------------------------------------------------
#  # Run Hybrid model (This might take several minutes to run)
#  Fit <- bakRFit(Fit, HybridFit = TRUE)
#  
#  # Run Full model (This might take ~10-30 minutes to run)
#  Fit <- bakRFit(Fit, StanFit = TRUE)
#  

## ---- fig.align='center'------------------------------------------------------
## MA Plot with Fast Fit
bakR::plotMA(Fit, Model = "MLE")


## ---- fig.align='center'------------------------------------------------------
## Volcano Plot with Fast Fit; significance assessed relative to an FDR control of 0.05
plotVolcano(Fit$Fast_Fit)


## ---- fig.align='center'------------------------------------------------------
## 2D PCA plot with replicate fraction news
  # The equivalent function prior to version 1.0.0 is FnPCA, now deprecated in 
  # favor of FnPCA2.
FnPCA2(Fit, Model = "MLE")


