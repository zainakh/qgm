% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantile_association.R
\name{split.quantile.ztest}
\alias{split.quantile.ztest}
\title{Performs a hypothesis test of two variables given a conditioning set and returns
a p-value on if they are independent (p > 0.05) or not (p <= 0.05).}
\usage{
\method{split}{quantile.ztest}(x, y, S, suffStat)
}
\arguments{
\item{x}{Index of a column}

\item{y}{Index of a column (not equal to x)}

\item{S}{Conditioning set to be used to determine if conditional independence exists (can be empty set)}

\item{suffStat}{The dataframe of data (there is no sufficient statistic for this calculation)}
}
\value{
pvalue corresponding to if the quantile level
}
\description{
In contrast from the original quantile.ztest
method, this method splits up data into halves and calculates the statistic
using a model trained on either the first or second half and predicts the
conditional quantile on the opposite half.
}
\details{
This test should have additional robustness and lower variance compared
to the original quantile.ztest method, but fits twice as many models.
}
