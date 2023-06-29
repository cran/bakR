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
library(pheatmap)

# Packages which are installed when bakR is installed
library(dplyr) 
library(magrittr) 
library(ggplot2) 
library(stats)

# Set the seed for reproducibility
set.seed(123)

## -----------------------------------------------------------------------------

### Only differential synthesis

# Simulate a nucleotide recoding dataset
sim_data <- Simulate_relative_bakRData(1000, 1000000,
                         num_ks_DE = c(0, 300),
                         num_kd_DE = c(0, 0))
  # This will simulate 1000 features, 1 million reads, 2 experimental conditions
  # and 3 replicates for each experimental condition. 300 features will be differentially
  # transcribed, and there will be no differential stability
  # See ?Simulate_bakRData for details regarding tunable parameters

# Extract simulated bakRData object
bakRData <- sim_data$bakRData

## Run the efficient model
Fit_s <- bakRFit(bakRData)


### Only differential stability

# Simulate a nucleotide recoding dataset
sim_data <- Simulate_relative_bakRData(1000, 1000000,
                         num_ks_DE = c(0, 0),
                         num_kd_DE = c(0, 300))
  # This will simulate 1000 features, 1 million reads, 2 experimental conditions
  # and 3 replicates for each experimental condition. 300 features will be differentially
  # transcribed, and there will be no differential stability
  # See ?Simulate_bakRData for details regarding tunable parameters

# Extract simulated bakRData object
bakRData <- sim_data$bakRData

## Run the efficient model
Fit_d <- bakRFit(bakRData)


## -----------------------------------------------------------------------------

### Only differential synthesis dataset

# Get the count matrix from bakR
Counts <- Fit_s$Data_lists$Count_Matrix

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

# Make DESeq2 data object
dds_s <- DESeqDataSetFromMatrix(countData = Counts,
                              colData = colData,
                              design = ~conditions)

# Fit DESeq2 model
ddso_s <- DESeq(dds_s)

# Extract results of experimental vs. reference comparison
reso_s <- results(ddso_s, contrast = c("conditions", "exp", "ref"))


### Only differential stability dataset

# Get the other count matrix from bakR
Counts <- Fit_d$Data_lists$Count_Matrix

# Make DESeq2 data object
dds_d <- DESeqDataSetFromMatrix(countData = Counts,
                              colData = colData,
                              design = ~conditions)

# Fit DESeq2 model
ddso_d <- DESeq(dds_d)

# Extract results of experimental vs. reference comparison
reso_d <- results(ddso_d, contrast = c("conditions", "exp", "ref"))



## ---- results = 'hide'--------------------------------------------------------

### Only differential synthesis: 

# Convert to data frame
reso_s <- as.data.frame(reso_s)

# Make data frame
DE_df_s <- data.frame(XF = row.names(reso_s),
                    L2FC_RNA = reso_s$log2FoldChange,
                    DE_score = reso_s$stat,
                    DE_se = reso_s$lfcSE,
                    DE_pval = reso_s$pval,
                    DE_padj = reso_s$padj)


### Only differential degradation:

# Convert to data frame
reso_d <- as.data.frame(reso_d)

# Make data frame
DE_df_d <- data.frame(XF = row.names(reso_d),
                    L2FC_RNA = reso_d$log2FoldChange,
                    DE_score = reso_d$stat,
                    DE_se = reso_d$lfcSE,
                    DE_pval = reso_d$pval,
                    DE_padj = reso_d$padj)



## ---- results = 'hide'--------------------------------------------------------

# Decreasing sims parameter to speed up; wouldn't normally suggest this if you
# want higher precision mechanism p-values, discussed later
Mechs_s <- DissectMechanism(Fit_s, DE_df_s,
                          sims = 1000000)

Mechs_d <- DissectMechanism(Fit_d, DE_df_d,
                          sims = 1000000)


## ---- fig.align='center'------------------------------------------------------

# Nice red to blue color gradient
  # Feel free to use any coloring your heart desires
col <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", 
         "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", 
         "#D6604D", "#B2182B", "#67001F")

# Heatmap for differential synthesis only dataset
pheatmap(Mechs_s$Heatmap_df, cluster_cols = FALSE, show_rownames = FALSE, color = col)

# Heatmap for differential degradation only dataset
pheatmap(Mechs_d$Heatmap_df, cluster_cols = FALSE, show_rownames = FALSE, color = col)


## ---- fig.align='center'------------------------------------------------------

# Scatter plot of L2FC(kdeg) vs. L2FC(RNA) colored by mechanism test stat
  # Gotta transform the mech_stat because it spans many orders of magnitude
ggplot(Mechs_s$Mechanism_df, aes(x = L2FC_kdeg, y = L2FC_RNA, color = log10(abs(mech_stat) + 1)*sign(mech_stat))) +
  geom_point() + 
  theme_classic() + 
  scale_color_viridis_c() + 
  xlab("L2FC(kdeg) from bakR") + 
  ylab("L2FC(RNA) from DESeq2") +
  labs(color = "Mechanism")


## ---- fig.align='center'------------------------------------------------------

# Scatter plot of L2FC(kdeg) vs. L2FC(RNA) colored by mechanism test stat
  # Gotta transform the mech_stat because it spans many orders of magnitude
ggplot(Mechs_d$Mechanism_df, aes(x = L2FC_kdeg, y = L2FC_RNA, color = log10(abs(mech_stat) + 1)*sign(mech_stat))) +
  geom_point() + 
  theme_classic() + 
  scale_color_viridis_c() + 
  xlab("L2FC(kdeg) from bakR") + 
  ylab("L2FC(RNA) from DESeq2") +
  labs(color = "Mechanism")


## ---- fig.align='center'------------------------------------------------------

## Grid to plot on

# Grid size in each axis
n <- 300

# DE z-scores
DE_score <- rep(seq(from = -10, to = 10, length.out = n), times = n)

# bakR z-scores
bakR_score <- rep(seq(from = -10, to = 10, length.out = n), each = n)

grid_df <- tibble(DE = DE_score,
                  bakR = bakR_score)


# Mechanism score
grid_df <- grid_df %>%
  mutate(Mech = ifelse(DE > 0, 
                       (bakR + 2)*DE,
                       (bakR - 2)*DE))

# Calculate p-value
grid_df <- grid_df %>%
  mutate(Mech_pval = stats::df(abs(Mech), df1 = 2, df2 = 2, ncp = 2))

# Visualize test statistic
ggplot(grid_df, aes(x = bakR, y = DE, z = 0.03*Mech)) + 
  geom_contour_filled() + 
  theme_classic() + 
  xlab("bakR score") + 
  ylab("DE score") 
  

