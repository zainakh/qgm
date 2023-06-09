% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot.quantile.dag}
\alias{plot.quantile.dag}
\title{Calculate the graph of variable relationships at a particular quantile level.}
\usage{
\method{plot}{quantile.dag}(data, tau, weights = "marginal", verbose = FALSE, split = FALSE)
}
\arguments{
\item{data}{A data frame of data with unique column/rownames}

\item{tau}{A particular quantile level (0 to 1, not inclusive)}

\item{weights}{Weights will either be "marginal" (marginal association)
or "conditional" (association conditioned on everything) - default is
"marginal" (otherwise they will not be used)}

\item{split}{If you should split the data for more efficient estimates (fits twice as many models)}
}
\value{
Nothing
}
\description{
It assumes the input dataframe has at least two columns (columns are variables) and that
the quantile level is between 0 and 1 (non-inclusive). Additionally it assumes all
column names are unique.
}
\details{
This function calls the pc algorithm with the quantile association hypothesis test.
Afterwards, it constructs the marginal associations of all variables and uses their
marginal associations as edge weights in the returned graph.

An igraph plot is plotted with the marginal associations at the quantile level
prescribed.
}
\examples{
plot.quantile.dag(df, tau=0.5)
}
