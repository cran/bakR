## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dim = c(6, 4)
)

## ----setup--------------------------------------------------------------------
library(bakR)

# Packages that are NOT automatically installed when bakR is installed
library(DESeq2)

# Packages which are installed when bakR is installed
library(dplyr) 
library(magrittr) 
library(ggplot2) 
library(stats)

# Set the seed for reproducibility
set.seed(123)

## -----------------------------------------------------------------------------
# Simulate a nucleotide recoding dataset
sim_data <- Simulate_bakRData(1000,
                         num_kd_DE = c(0, 0),
                         num_ks_DE = c(0, 200))
  # This will simulate 500 features, 2 experimental conditions
  # and 3 replicates for each experimental condition
  # See ?Simulate_bakRData for details regarding tunable parameters

# Extract simulated bakRData object
bakRData <- sim_data$bakRData

# Extract simualted ground truths
sim_truth <- sim_data$sim_list

## Run the efficient model
# We'll tell it what the mutation rates are just for efficiency's sake
Fit <- bakRFit(bakRData, pnew = rep(0.05, times = 6), pold = 0.001)



## ---- results = 'hide'--------------------------------------------------------

# Get the count matrix from bakR
Counts <- Fit$Data_lists$Count_Matrix

# Experimental conditions for each sample
# There are 6 s4U treated samples (3 replicates of each condition)
# In addition, there are 2 -s4U control samples (1 for each condition)

## s4U conditions
# 1st three samples are reference (ref) samples
# Next three samples are experimental (exp) samples
conditions_s4U <- as.factor(rep(c("ref", "exp"), each = 3))

## -s4U control conditions
# 1st sample is reference, next is experimental
conditions_ctl <- as.factor(c("ref", "exp"))

# Combined s4U and -s4U control conditions
conditions <- c(conditions_s4U, conditions_ctl)

# Make the colData input for DESeq2
colData <- data.frame(conditions = conditions)
rownames(colData) <- colnames(Counts)

# Take a look at the colData object
print(t(colData))


## ---- echo = FALSE, warning = FALSE-------------------------------------------
knitr::kable(t(colData))


## -----------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = Counts,
                              colData = colData,
                              design = ~conditions)


## -----------------------------------------------------------------------------
ddso <- DESeq(dds)

# Extract results of experimental vs. reference comparison
reso <- results(ddso, contrast = c("conditions", "exp", "ref"))

# Look at the column names of the reso object
colnames(as.data.frame(reso))
                

## -----------------------------------------------------------------------------
ksyn_df <- data.frame(L2FC = reso$log2FoldChange + Fit$Fast_Fit$Effects_df$L2FC_kdeg,
                      Gene = Fit$Fast_Fit$Effects_df$XF)
                

## -----------------------------------------------------------------------------
ksyn_df$se <- sqrt(reso$lfcSE^2 + (Fit$Fast_Fit$Effects_df$se*log2(exp(1)))^2 ) 
                

## -----------------------------------------------------------------------------
# Calculate p-value using asymptotic Wald test
ksyn_df <- ksyn_df %>%
  mutate(pval = 2*pnorm(-abs(L2FC/se)),
         padj = p.adjust(pval, method = "BH"))
                

## ---- fig.align='center'------------------------------------------------------

# Add conclusion at 0.01 FDR control
ksyn_df <- ksyn_df %>%
  mutate(conclusion = as.factor(ifelse(padj < 0.01, 
                             ifelse(L2FC < 0, "Decreased txn", "Increased txn"),
                             "Not sig.")))

# Make volcano plot
ksyn_volc <- ggplot(ksyn_df, aes(x = L2FC, y = -log10(padj), color = conclusion)) + 
  theme_classic() + 
  geom_point(size = 1) + 
  xlab("L2FC(ksyn)") + 
  ylab("-log10(padj)") + 
  scale_color_manual(values = c("blue", "orange", "gray"))

# Observe volcano plot
ksyn_volc


## ---- fig.align='center'------------------------------------------------------

# Add conclusion at 0.01 FDR control
ksyn_df <- ksyn_df %>%
  mutate(result = as.factor(ifelse(padj < 0.01, 
                             ifelse(as.integer(Gene) <= 800, "FP", "TP"),
                             ifelse(as.integer(Gene) <= 800, "TN", "FN"))))

# Make volcano plot
ksyn_results <- ggplot(ksyn_df, aes(x = L2FC, y = -log10(padj), color = result)) + 
  theme_classic() + 
  geom_point(size = 1) + 
  xlab("L2FC(ksyn)") + 
  ylab("-log10(padj)") + 
  scale_color_manual(values = c("black", "gray", "forestgreen", "blue"))

# Observe volcano plot
ksyn_results


