% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualizations.R
\name{VisualizeDropout}
\alias{VisualizeDropout}
\title{Visualize dropout}
\usage{
VisualizeDropout(obj, keep_data = FALSE, no_message = FALSE)
}
\arguments{
\item{obj}{bakRFit or bakRFnFit object}

\item{keep_data}{Logical; if TRUE, will return data used to make plots along with
the plots themselves}

\item{no_message}{Logical; if TRUE, will not output message regarding estimated
rates of dropout in each sample}
}
\value{
If keep_data is FALSE, then a list of \code{ggplot} objects are returned, one
for each +s4U sample. The plots show the relationship between a feature's fraction new
and the difference between its +s4U and -s4U read coverage. Nonlinear-least squares fit
is plotted on top of points as a blue line. If keep_data is TRUE, then the data used
to make the plots is returned in addition to the list of plots.
}
\description{
\code{VisualizeDropout} fits dropout model with \code{QuantifyDropout},
reports the fit results, and then generates a ggplot object showing the
data used to infer the fit as well as the fitted nonlinear trend.
}
\examples{
\donttest{
# Simulate data for 500 genes and 2 replicates with 40\% dropout
sim <- Simulate_relative_bakRData(500, 100000, nreps = 2, p_do = 0.4)

# Fit data with fast implementation
Fit <- bakRFit(sim$bakRData)

# Quantify dropout
DO_plots <- VisualizeDropout(Fit)

}
}
