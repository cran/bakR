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
# Load GRAND-SLAM table and metadf
data("GS_table")
data("metadf")

# Create fns data frame from GRAND-SLAM output
fns <- GSprocessing(GS_table)


# Create bakRFnData object
bfndo <- bakRFnData(fns, metadf)


## -----------------------------------------------------------------------------
# Run the efficient model
Fit <- bakRFit(bfndo)


## ---- eval = FALSE------------------------------------------------------------
#  # Run Hybrid model (This might take several minutes to run)
#  Fit <- bakRFit(Fit, HybridFit = TRUE)

## ---- fig.align='center'------------------------------------------------------
## MA Plot with Fast Fit
bakR::plotMA(Fit, Model = "MLE")


## ---- fig.align='center'------------------------------------------------------
## Volcano Plot with Fast Fit; significance assessed relative to an FDR control of 0.05
plotVolcano(Fit$Fast_Fit)


## ---- fig.align='center'------------------------------------------------------
## 2D PCA plot with replicate fraction news
FnPCA2(Fit, Model = "MLE")


