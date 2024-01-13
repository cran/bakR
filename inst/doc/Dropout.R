## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dim = c(6, 4)
)

## ----setup--------------------------------------------------------------------
library(bakR)
set.seed(123)

## -----------------------------------------------------------------------------
# Simulate a nucleotide recoding dataset
sim_data <- Simulate_relative_bakRData(1000, depth = 1000000, nreps = 2,
                                       p_do = 0.4)
  # This will simulate 500 features, 500,000 reads, 2 experimental conditions
  # and 2 replicates for each experimental condition.
  # 40% dropout is simulated.
  # See ?Simulate_relative_bakRData for details regarding tunable parameters

# Run the efficient model
Fit <- bakRFit(sim_data$bakRData)


## -----------------------------------------------------------------------------
# Correct dropout-induced biases
Fit_c <- CorrectDropout(Fit)
  # You can also overwite the existing bakRFit object.
  # I am creating a separate bakRFit object to make comparisons later in this vignette.


## ----fig.align='center'-------------------------------------------------------
# Correct dropout-induced biases
Vis_DO <- VisualizeDropout(Fit)

# Visualize dropout for 1st replicate of reference condition
Vis_DO$ExpID_1_Rep_1


## ----fig.align='center'-------------------------------------------------------

# Extract simualted ground truths
sim_truth <- sim_data$sim_list

# Features that made it past filtering
XFs <- unique(Fit$Fast_Fit$Effects_df$XF)

# Simulated logit(fraction news) from features making it past filtering
true_fn <- sim_truth$Fn_rep_sim$Logit_fn[sim_truth$Fn_rep_sim$Feature_ID %in% XFs]

# Estimated logit(fraction news)
est_fn <- Fit$Fast_Fit$Fn_Estimates$logit_fn

# Compare estimate to truth
plot(true_fn, est_fn, xlab = "True logit(fn)", ylab = "Estimated logit(fn)")
abline(0, 1, col = "red")


## ----fig.align='center'-------------------------------------------------------

# Features that made it past filtering
XFs <- unique(Fit_c$Fast_Fit$Effects_df$XF)

# Simulated logit(fraction news) from features making it past filtering
true_fn <- sim_truth$Fn_rep_sim$Logit_fn[sim_truth$Fn_rep_sim$Feature_ID %in% XFs]

# Estimated logit(fraction news)
est_fn <- Fit_c$Fast_Fit$Fn_Estimates$logit_fn

# Compare estimate to truth
plot(true_fn, est_fn, xlab = "True logit(fn)", ylab = "Estimated logit(fn)")
abline(0, 1, col = "red")


