#' Adds jitter to the input dataframe to break ties.
#'
#' Goes column by column in the dataframe and adds jitter from the distribution
#'
#' Unif(-a, a)
#'
#' where a = factor * d /5 (and d is the minimum distance between any two points in a vector).
#'
#'
#' @param data Input dataframe of data
#' @param factor A float to scale the jittering by (if greater than 2.5, there is a risk of quantile crossings)
#' @return Jittered version of your input dataframe
jitter.columns <- function(data, factor=0.1) {
  num_cols <- length((colnames(data)))

  for(x in 1:num_cols) {
    data[,x] <- jitter(data[,x], factor=factor)
  }
  return(data)
}


#' Performs a hypothesis test of two variables given a conditioning set and returns
#' a p-value on if they are independent (p > 0.05) or not (p <= 0.05).
#'
#' Uses a hypothesis test that considers the proportion of times two variables
#' have concordant (conditional) quantile regression residuals. Counts the
#' amount of times the residual pairs are concordant (jointly above or below their
#' respective conditional quantile regression lines) and normalizes that concordant
#' proportion. Normalization is based on what the proportions would be
#' under the null hypothesis. From the normalized variable, p values of the
#' relationship between columns x and y can be returned.
#'
#'
#' @param x Index of a column
#' @param y Index of a column (not equal to x)
#' @param S Conditioning set to be used to determine if conditional independence exists (can be empty set)
#' @param suffStat The dataframe of data (there is no sufficient statistic for this calculation)
#' @return pvalue corresponding to if the quantile level
quantile.ztest <- function (x, y, S, suffStat) {
  tau <- readRDS("tau.rds")
  n <- length(suffStat[,x])

  var1 <- suffStat[,x]
  var2 <- suffStat[,y]

  if(length(S) == 0){
    q1 <- quantreg::rq(var1 ~ 1, tau=tau)
    q2 <- quantreg::rq(var2 ~ 1, tau=tau)

    fit.q1 <- fitted(q1)
    fit.q2 <- fitted(q2)
  }
  else {
    Z <- sapply(suffStat[,S], as.numeric)

    q1 <- quantreg::rq(var1 ~ Z, tau=tau)
    q2 <- quantreg::rq(var2 ~ Z, tau=tau)

    fit.q1 <- fitted(q1)
    fit.q2 <- fitted(q2)
  }

  ptilde <- sum((var1 < fit.q1) & (var2 < fit.q2)) / n
  phat <- sum((var1 > fit.q1) & (var2 > fit.q2)) / n

  zstat_below <- (ptilde - tau^2) / sqrt( tau^2 * (1-tau)^2 / n)
  zstat_above <- (phat - (1 - tau)^2) / sqrt( tau^2 * (1-tau)^2 / n)

  p_val_below <- 2 * pnorm(abs(zstat_below), lower.tail = FALSE)
  p_val_above <- 2 * pnorm(abs(zstat_above), lower.tail = FALSE)

  if (tau < 0.5) {
    return(p_val_below)
  }
  else {
    return(p_val_above)
  }
}

#' Returns the  quantile association for all pairs of variables in a dataframe.
#'
#' Runs a double for loop across all columns and calculates the quantile association
#' statistic for all possible pairs. Depending on if the quantile level is above or
#' below the median, the statistic will consider if the residuals are jointly above or
#' jointly below their respective marginal regression lines.
#'
#' This test has two types where it considers marginal associations of variables or
#' considers the conditional case where each quantile regression is conditional w.r.t
#' all other variables in the dataframe.
#'
#' @param data Input dataframe of data
#' @param tau A particular quantile level (0 to 1, not inclusive)
#' @param type If weights of the regression should be marginal or conditional on all other variables
#' @return An n by n matrix (where n is the number of columns in data) that contains the marginal relationships of each pair of columns
pairwise.test <- function(data, tau, weights="marginal") {
  num_cols <- length((colnames(data)))

  zstat_above <- matrix(0, nrow=num_cols, ncol=num_cols)
  zstat_below <- matrix(0, nrow=num_cols, ncol=num_cols)

  for(x in 1:num_cols) {
    for(y in 1:num_cols) {

      n <- length(data[,x])

      var1 <- data[,x]
      var2 <- data[,y]

      if(weights == "marginal") {
        q1 <- quantreg::rq(var1 ~ 1, tau=tau)
        q2 <- quantreg::rq(var2 ~ 1, tau=tau)
      }
      else {
        col1 <- colnames(data)[x]
        col2 <- colnames(data)[y]
        q1 <- quantreg::rq(as.formula(paste0(col1, " ~ .")), tau=tau, data=data)
        q2 <- quantreg::rq(as.formula(paste0(col2, " ~ .")), tau=tau, data=data)
      }

      fit.q1 <- fitted(q1)
      fit.q2 <- fitted(q2)

      ptilde <- sum((var1 < fit.q1) & (var2 < fit.q2)) / n
      phat <- sum((var1 > fit.q1) & (var2 > fit.q2)) / n

      zstat_below[x, y] <- (ptilde - tau^2) / sqrt( tau^2 * (1-tau)^2 / n)
      zstat_above[x, y] <- (phat - (1 - tau)^2) / sqrt( tau^2 * (1-tau)^2 / n)

    }
  }

  if(tau >= 0.5){
    return( zstat_above )
  }
  return( zstat_below )
}
