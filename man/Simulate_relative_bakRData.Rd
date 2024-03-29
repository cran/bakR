% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simulations.R
\name{Simulate_relative_bakRData}
\alias{Simulate_relative_bakRData}
\title{Simulating nucleotide recoding data with relative count data}
\usage{
Simulate_relative_bakRData(
  ngene,
  depth,
  num_conds = 2L,
  nreps = 3L,
  eff_sd = 0.75,
  eff_mean = 0,
  kdlog_mean = -1.8,
  kdlog_sd = 0.65,
  kslog_mean = 1,
  kslog_sd = 0.65,
  tl = 2,
  p_new = 0.05,
  p_old = 0.001,
  read_lengths = 200L,
  p_do = 0,
  noise_deg_a = -0.3,
  noise_deg_b = -1.5,
  noise_synth = 0.1,
  sd_rep = 0.05,
  low_L2FC_ks = -1,
  high_L2FC_ks = 1,
  num_kd_DE = c(0L, as.integer(rep(round(as.integer(ngene)/2), times =
    as.integer(num_conds) - 1))),
  num_ks_DE = rep(0L, times = as.integer(num_conds)),
  sim_read_counts = TRUE,
  a1 = 5,
  a0 = 0.01,
  nreads = 50L,
  alpha = 25,
  beta = 75,
  STL = FALSE,
  STL_len = 40,
  lprob_U_sd = 0,
  lp_sd = 0
)
}
\arguments{
\item{ngene}{Number of genes to simulate data for}

\item{depth}{Total number of reads to simulate}

\item{num_conds}{Number of experimental conditions (including the reference condition) to simulate}

\item{nreps}{Number of replicates to simulate}

\item{eff_sd}{Effect size; more specifically, the standard deviation of the normal distribution from which non-zero
changes in logit(fraction new) are pulled from.}

\item{eff_mean}{Effect size mean; mean of normal distribution from which non-zero changes in logit(fraction new) are pulled from.
Note, setting this to 0 does not mean that some of the significant effect sizes will be 0, as any exact integer is impossible
to draw from a continuous random number generator. Setting this to 0 just means that there is symmetric stabilization and destabilization}

\item{kdlog_mean}{Degradation rate constants will be drawn from lognormal distribution with this logmean}

\item{kdlog_sd}{Degradation rate constants will be drawn from lognormal distribution with this logsd}

\item{kslog_mean}{Synthesis rate constants will be drawn from a lognormal distribution with this mean}

\item{kslog_sd}{Synthesis rate constants will be drawn from a lognormal distribution with this logsd}

\item{tl}{metabolic label feed time}

\item{p_new}{metabolic label (e.g., s4U) induced mutation rate. Can be a vector of length num_conds}

\item{p_old}{background mutation rate}

\item{read_lengths}{Total read length for each sequencing read (e.g., PE100 reads correspond to read_lengths = 200)}

\item{p_do}{Rate at which metabolic label containing reads are lost due to dropout; must be between 0 and 1}

\item{noise_deg_a}{Slope of trend relating log10(standardized read counts) to log(replicate variability)}

\item{noise_deg_b}{Intercept of trend relating log10(standardized read counts) to log(replicate variability)}

\item{noise_synth}{Homoskedastic variability of L2FC(ksyn)}

\item{sd_rep}{Variance of lognormal distribution from which replicate variability is drawn}

\item{low_L2FC_ks}{Most negative L2FC(ksyn) that can be simulated}

\item{high_L2FC_ks}{Most positive L2FC(ksyn) that can be simulated}

\item{num_kd_DE}{Vector where each element represents the number of genes that show a significant change in stability relative
to the reference. 1st entry must be 0 by definition (since relative to the reference the reference sample is unchanged)}

\item{num_ks_DE}{Same as num_kd_DE but for significant changes in synthesis rates.}

\item{sim_read_counts}{Logical; if TRUE, read counts are simulated as coming from a heterodisperse negative binomial distribution}

\item{a1}{Heterodispersion 1/reads dependence parameter}

\item{a0}{High read depth limit of negative binomial dispersion parameter}

\item{nreads}{Number of reads simulated if sim_read_counts is FALSE}

\item{alpha}{shape1 parameter of the beta distribution from which U-contents (probability that a nucleotide in a read from a transcript is a U) are
drawn for each gene.}

\item{beta}{shape2 parameter of the beta distribution from which U-contents (probability that a nucleotide in a read from a transcript is a U) are
drawn for each gene.}

\item{STL}{logical; if TRUE, simulation is of STL-seq rather than a standard TL-seq experiment. The two big changes are that a short read length is required
(< 60 nt) and that every read for a particular feature will have the same number of Us. Only one read length is simulated for simplicity.}

\item{STL_len}{Average length of simulated STL-seq length. Since Pol II typically pauses about 20-60 bases
from the promoter, this should be around 40}

\item{lprob_U_sd}{Standard deviation of the logit(probability nt is a U) for each sequencing read. The number of Us in a
sequencing read are drawn from a binomial distribution with prob drawn from a logit-Normal distribution with this logit-sd.}

\item{lp_sd}{Standard deviation of logit(probability a U is mutated) for each U. The number of mutations in a given read is the sum of
nU Bernoulli random variables, where nU is the number of Us, and p is drawn from a logit-normal distribution with lp_sd standard deviation
on logit scale.}
}
\value{
A list containing a simulated \code{bakRData} object as well as a list of simulated kinetic parameters of interest.
The contents of the latter list are:
\itemize{
\item Effect_sim; Dataframe meant to mimic formatting of Effect_df that are part of \code{bakRFit(StanFit = TRUE)}, \code{bakRFit(HybridFit = TRUE)} and \code{bakRFit(bakRData object)} output.
\item Fn_mean_sim; Dataframe meant to mimic formatting of Regularized_ests that is part of \code{bakRFit(bakRData object)} output. Contains information
about the true fraction new simulated in each condition (the mean of the normal distribution from which replicate fraction news are simulated)
\item Fn_rep_sim; Dataframe meant to mimic formatting of Fn_Estimates that is part of \code{bakRFit(bakRData object)} output. Contains information
about the fraction new simulated for each feature in each replicate of each condition.
\item L2FC_ks_mean; The true L2FC(ksyn) for each feature in each experimental condition. The i-th column corresponds to the L2FC(ksyn) when comparing
the i-th condition to the reference condition (defined as the 1st condition) so the 1st column is always all 0s
\item RNA_conc; The average number of normalized read counts expected for each feature in each sample.
}
}
\description{
\code{Simulate_relative_bakRData} simulates a \code{bakRData} object. It's output also includes the simulated
values of all kinetic parameters of interest.
}
\details{
The main difference between \code{Simulate_relative_bakRData}
and \code{Simulate_bakRData} is that the former requires both the number of
genes (\code{ngene}) and the total number of reads (\code{depth}) has to be set.
In the latter, only the number of genes is set, and the number of reads for each
gene is simulated so that no matter how many genes are simulated, the number of
reads given default parameters is reflective of what is seen in 20,000,000 read
human RNA-seq libraries. The benefit of \code{Simulate_relative_bakRData} is that it is
easier to test the impact of depth on model performance. This can theoretically
be done by changing the synthesis rate constant parameters in \code{Simulate_bakRData},
but the relationship between these parameters and sequencing depth is unintuitive. The
benefit of \code{Simulate_bakRData} is that fewer genes can be simulated
while still yielding reasonable per-gene coverage without figuring out what the
total depth in the small gene subset should be. This is nice for testing bakR and
other analysis tools on small datasets. \code{Simulate_relative_bakRData} is a more
realistic simulation that better accounts for the relative nature of RNA-seq read
counts (i.e., expected number of reads from a given feature is related to proportion of RNA molecules
coming from that feature).

Another difference between \code{Simulate_relative_bakRData} and \code{Simulate_bakRData}
is that \code{Simulate_relative_bakRData} uses the label time and simulated degradation
rate constants to infer the fraction new, whereas \code{Simulate_bakRData} uses simulated
fraction news and the label time to infer the degradation rate constants. Thus,
\code{Simulate_relative_bakRData} is preferable for assessing the impact of label
time on model performance (since it will have a realistic impact on the fraction new,
and the distribution of fraction news has a major impact on model performance).
Similarly, \code{Simulate_bakRData} is preferable for directly assessing the impact of
fraction news on model performance, without having to think about how both the label
time and simulated degradation rate constant distribution.

If investigating dropout, only \code{Simulate_relative_bakRData} should be used, as the
accurate simulation of read counts as being a function of the relative abundance of
each RNA feature is crucial to accurately simulate dropout.

Function to simulate a \code{bakRData} object according to a realistic generative model
}
\examples{
\donttest{
# 2 replicate, 2 experimental condition, 1000 gene simulation
sim_2reps <- Simulate_relative_bakRData(ngene = 1000, depth = 100000,
                               nreps = 2)

# 3 replicate, 2 experimental condition, 1000 gene simulation
# with 100 instances of differential degradation kinetics
sim_3reps <- Simulate_relative_bakRData(ngene = 1000, depth = 100000,
                                        num_kd_DE = c(0, 100))

# 2 replicates, 3 experimental condition, 1000 gene simulation
# with 100 instances of differential degradation kinetics in the 1st
# condition and no instances of differential degradation kinetics in the
# 2nd condition
sim_3es <- Simulate_relative_bakRData(ngene = 1000, depth = 100000,
                             nreps = 2,
                             num_conds = 3,
                             num_kd_DE = c(0, 100, 0))

}
}
