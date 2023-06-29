## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dim = c(6, 4)
)

## ----setup--------------------------------------------------------------------
library(bakR)
library(ggplot2)
library(corrplot)
set.seed(123)

## -----------------------------------------------------------------------------
# Simulate a nucleotide recoding dataset
sim_data <- Simulate_relative_bakRData(1000, depth = 1000000,
                                       nreps = 2)
  # This will simulate 1000 features, 1000000 reads, 2 experimental conditions,
  # and 2 replicates for each experimental condition
  # See ?Simulate_relative_bakRData for details regarding tunable parameters

# Run the efficient model
Fit <- bakRFit(sim_data$bakRData)

# Quality control checks
QC <- QC_checks(Fit)


## ---- fig.align='center'------------------------------------------------------
# Barplots of raw mutation rates
QC$raw_mutrates

## ---- fig.align='center'------------------------------------------------------
# Barplots of inferred mutation rates
QC$conversion_rates

## ---- fig.align='center'------------------------------------------------------
# Barplots of inferred mutation rates
ggplot(Fit$Fast_Fit$Fn_Estimates, aes(x = logit_fn, color = as.factor(sample))) + 
  geom_density() + 
  theme_classic() +
  scale_color_viridis_d() + 
  xlab("logit(fn) estimates") + 
  ylab("density")

## ---- fig.align='center'------------------------------------------------------
# Barplots of inferred mutation rates
  # Numerical indices are ordered as they appear in QC_checks() output message
  # So this is for replicate 1 and 2 of experimental ID 1
QC$correlation_plots[[1]]

## ---- fig.align='center'------------------------------------------------------
# Using function from corrplot package
corrplot.mixed(QC$correlation_matrix, 
               upper = "square", lower = "number", 
               addgrid.col = "black", tl.col = "black")

## -----------------------------------------------------------------------------
# Seed for reproducibility
set.seed(321)

# Simulate a nucleotide recoding dataset
sim_data <- Simulate_bakRData(1000, nreps = 2, fn_mean = -4)
  # This will simulate 1000 features, 2 experimental conditions,
  # and 2 replicates for each experimental condition
  # The average logit(fn) will be -4, which corresponds to an average fn of just under 0.02.

# Run the efficient model
  # Check the pnew estimates, which should all be around 0.05
Fit <- bakRFit(sim_data$bakRData)


## -----------------------------------------------------------------------------
# Run QC_checks and read messages
QC <- QC_checks(Fit)


## -----------------------------------------------------------------------------
# Rerun with Stan-based pnew estimation
  # This will take a couple minutes to run
Fit_s <- bakRFit(Fit, FastRerun = TRUE, StanRateEst = TRUE)


## ---- fig.align='center'------------------------------------------------------

# Simulated ground truth
sim_truth <- sim_data$sim_list

# Features that made it past filtering
XFs <- unique(Fit$Fast_Fit$Effects_df$XF)

# Simulated logit(fraction news) from features making it past filtering
true_fn <- sim_truth$Fn_rep_sim$Logit_fn[sim_truth$Fn_rep_sim$Feature_ID %in% XFs &
                                           sim_truth$Fn_rep_sim$Exp_ID == 2]

# Estimated logit(fraction news)
est_fn <- Fit$Fast_Fit$Fn_Estimates$logit_fn[Fit$Fast_Fit$Fn_Estimates$Exp_ID == 2]

# Compare estimate to truth
plot(true_fn, est_fn, xlab = "True logit(fn)", ylab = "Estimated logit(fn)",
     main = "Default pnew estimation",
     xlim = c(-8, 6),
     ylim = c(-8, 6))
abline(0, 1, col = "red")

# Estimated logit(fraction news)
est_fn <- Fit_s$Fast_Fit$Fn_Estimates$logit_fn[Fit_s$Fast_Fit$Fn_Estimates$Exp_ID == 2]

# Compare estimate to truth
plot(true_fn, est_fn, xlab = "True logit(fn)", ylab = "Estimated logit(fn)",
     main = "Alternative pnew estimation",
     xlim = c(-8, 6),
     ylim = c(-8, 6))
abline(0, 1, col = "red")



## -----------------------------------------------------------------------------
# Rerun with Stan-based pnew estimation
  # This will take a couple minutes to run
Fit_u <- bakRFit(Fit, FastRerun = TRUE, pnew = 0.05)


## ---- fig.align='center'------------------------------------------------------

# Estimated logit(fraction news)
est_fn <- Fit_u$Fast_Fit$Fn_Estimates$logit_fn[Fit_u$Fast_Fit$Fn_Estimates$Exp_ID == 2]

# Compare estimate to truth
plot(true_fn, est_fn, xlab = "True logit(fn)", ylab = "Estimated logit(fn)",
     main = "User provided pnew",
     xlim = c(-8, 6),
     ylim = c(-8, 6))
abline(0, 1, col = "red")



