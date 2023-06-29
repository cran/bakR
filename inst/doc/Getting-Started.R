## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dim = c(6, 4)
)

## ----setup--------------------------------------------------------------------
library(bakR)
set.seed(123)

## ---- results = "hide"--------------------------------------------------------
# Load small example cB
data("cB_small")

# Observe contents of cB
head(cB_small)


## ---- echo = FALSE, warning = FALSE-------------------------------------------

knitr::kable(head(cB_small), "pipe")


## ---- results = "hide"--------------------------------------------------------
# Load metadf data frame; will be loaded as metadf in global environment
data("metadf")

# Print the transpose of metadf
# Rows will be the columns and columns will be the rows
print(t(metadf))


## ---- echo = FALSE, warning = FALSE-------------------------------------------
knitr::kable(t(metadf))


## -----------------------------------------------------------------------------
# metadf row names
print(rownames(metadf))


# cB sample names
print(unique(cB_small$sample))


## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
# Run efficient model with known mutation rates
# Pass the Fit object rather than the bakRData object and set FastRerun to TRUE
Fit <- bakRFit(Fit,
                     FastRerun = TRUE,
                     pnew = rep(0.05, times = 6), 
                     pold = 0.001)

## ---- eval = FALSE------------------------------------------------------------
#  # Set StanRateEst to TRUE to use Stan to estimate rates
#  # low_reads and high_reads defines the read count cutoffs used to select features
#    # default = between 1000 and 5000 reads
#  # RateEst_size determines the number of features to use (default = 30)
#  Fit <- bakRFit(Fit,
#                       FastRerun = TRUE,
#                       StanRateEst = TRUE)
#  
#  

## ---- fig.align='center'------------------------------------------------------

# Features that made it past filtering
XFs <- unique(Fit$Fast_Fit$Effects_df$XF)

# Simulated logit(fraction news) from features making it past filtering
true_fn <- sim_truth$Fn_rep_sim$Logit_fn[sim_truth$Fn_rep_sim$Feature_ID %in% XFs]

# Estimated logit(fraction news)
est_fn <- Fit$Fast_Fit$Fn_Estimates$logit_fn

# Compare estimate to truth
plot(true_fn, est_fn, xlab = "True logit(fn)", ylab = "Estimated logit(fn)")
abline(0, 1, col = "red")


## ---- eval = FALSE------------------------------------------------------------
#  # Load options that will make running models more efficient
#  rstan::rstan_options(auto_write = TRUE)
#  options(mc.cores = parallel::detectCores())
#  
#  
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


