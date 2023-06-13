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

#' Performs a hypothesis test of two variables given a conditioning set and returns
#' a p-value on if they are independent (p > 0.05) or not (p <= 0.05).
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
#' @param split If you should split the data for more efficient estimates (fits twice as many models)
#' @return An n by n matrix (where n is the number of columns in data) that contains
#' the marginal relationships of each pair of columns
pairwise.test <- function(data, tau, weights="marginal", split=FALSE) {
  num_cols <- length((colnames(data)))

  zstat_above <- matrix(0, nrow=num_cols, ncol=num_cols)
  zstat_below <- matrix(0, nrow=num_cols, ncol=num_cols)

  for(x in 1:num_cols) {
    for(y in 1:num_cols) {

      n <- length(data[,x])

      if(split) { # Split data
        data.first <- data[1:(n%/%2), ]
        data.second <- data[(1 + (n%/%2)):n, ]

        first.half.len <- length(data.first[,x])
        second.half.len <- length(data.second[,x])

        var1.first <- data.first[, x]
        var2.first <- data.first[, y]
        var1.second <- data.second[, x]
        var2.second <- data.second[, y]

        if(weights == "marginal") {
          q1.first <- quantreg::rq(as.formula(paste(colnames(data)[x], "~ 1")), tau = tau, data=data.first)
          q2.first <- quantreg::rq(as.formula(paste(colnames(data)[y], "~ 1")), tau = tau, data=data.first)
          q1.second <- quantreg::rq(as.formula(paste(colnames(data)[x], "~ 1")), tau = tau, data=data.second)
          q2.second <- quantreg::rq(as.formula(paste(colnames(data)[y], "~ 1")), tau = tau, data=data.second)
        }
        else {
          q1.first <- quantreg::rq(as.formula(paste(colnames(data)[x], "~ .")), tau = tau, data=data.first)
          q2.first <- quantreg::rq(as.formula(paste(colnames(data)[y], "~ .")), tau = tau, data=data.first)
          q1.second <- quantreg::rq(as.formula(paste(colnames(data)[x], "~ .")), tau = tau, data=data.second)
          q2.second <- quantreg::rq(as.formula(paste(colnames(data)[y], "~ .")), tau = tau, data=data.second)
        }

        # Even if the weights are marginal, this prediction will still work even though we give newdata when its a marginal formula
        pred.q1.first <- predict(q1.first, newdata=data.second[, -c(x), drop=FALSE])
        pred.q2.first <- predict(q2.first, newdata=data.second[, -c(y), drop=FALSE])
        pred.q1.second <- predict(q1.second, newdata=data.first[, -c(x), drop=FALSE])
        pred.q2.second <- predict(q2.second, newdata=data.first[, -c(y), drop=FALSE])

        ptilde_b.first <- sum((var1.first < pred.q1.second) & (var2.first < pred.q2.second)) * (1 / first.half.len)
        phat_a.first <- sum((var1.first > pred.q1.second) & (var2.first > pred.q2.second)) * (1 / first.half.len)
        ptilde_b.second <- sum((var1.second < pred.q1.first) & (var2.second < pred.q2.first)) * (1 / second.half.len)
        phat_a.second <- sum((var1.second > pred.q1.first) & (var2.second > pred.q2.first)) * (1 / second.half.len)

        zstat_b.first <- (ptilde_b.first - tau^2) / sqrt( tau^2 * (1-tau)^2 / (first.half.len) )
        zstat_a.first <- (phat_a.first - (1 - tau)^2) / sqrt( tau^2 * (1-tau)^2 / (first.half.len) )
        zstat_b.second <- (ptilde_b.second - tau^2) / sqrt( tau^2 * (1-tau)^2 / (second.half.len) )
        zstat_a.second <- (phat_a.second - (1 - tau)^2) / sqrt( tau^2 * (1-tau)^2 / (second.half.len) )

        zstat_above[x, y] <- ((zstat_a.first + zstat_a.second) / sqrt(2))
        zstat_below[x, y] <- ((zstat_b.first + zstat_b.second) / sqrt(2))
      }
      else { # Don't split data
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
  }

  if(tau >= 0.5){
    return( zstat_above )
  }
  return( zstat_below )
}


#' Performs a hypothesis test of two variables given a conditioning set and returns
#' a p-value on if they are independent (p > 0.05) or not (p <= 0.05).
#'
#' In contrast from the original quantile.ztest
#' method, this method splits up data into halves and calculates the statistic
#' using a model trained on either the first or second half and predicts the
#' conditional quantile on the opposite half.
#'
#' This test should have additional robustness and lower variance compared
#' to the original quantile.ztest method, but fits twice as many models.
#'
#' @param x Index of a column
#' @param y Index of a column (not equal to x)
#' @param S Conditioning set to be used to determine if conditional independence exists (can be empty set)
#' @param suffStat The dataframe of data (there is no sufficient statistic for this calculation)
#' @return pvalue corresponding to if the quantile level
split.quantile.ztest <- function(x, y, S, suffStat) {
  tau <- readRDS("tau.rds")
  n <- length(data[,x])

  data.first <- data[1:(n%/%2), ]
  data.second <- data[(1 + (n%/%2)):n, ]

  first.half.len <- length(data.first[,x])
  second.half.len <- length(data.second[,x])

  var1.first <- data.first[, x]
  var2.first <- data.first[, y]
  var1.second <- data.second[, x]
  var2.second <- data.second[, y]

  q1.first <- quantreg::rq(
    as.formula(paste(colnames(data)[x], "~",
                     paste(colnames(data)[S], collapse = "+"),
                     sep = ""
    )),
    tau = tau,
    data=data.first
  )
  q2.first <- quantreg::rq(
    as.formula(paste(colnames(data)[y], "~",
                     paste(colnames(data)[S], collapse = "+"),
                     sep = ""
    )),
    tau = tau,
    data=data.first
  )
  q1.second <- quantreg::rq(
    as.formula(paste(colnames(data)[x], "~",
                     paste(colnames(data)[S], collapse = "+"),
                     sep = ""
    )),
    tau = tau,
    data=data.second
  )
  q2.second <- quantreg::rq(
    as.formula(paste(colnames(data)[y], "~",
                     paste(colnames(data)[S], collapse = "+"),
                     sep = ""
    )),
    tau = tau,
    data=data.second
  )

  pred.q1.first <- predict(q1.first, newdata=data.second[, S, drop=FALSE])
  pred.q2.first <- predict(q2.first, newdata=data.second[, S, drop=FALSE])
  pred.q1.second <- predict(q1.second, newdata=data.first[, S, drop=FALSE])
  pred.q2.second <- predict(q2.second, newdata=data.first[, S, drop=FALSE])

  ptilde_b.first <- sum((var1.first < pred.q1.second) & (var2.first < pred.q2.second)) * (1 / first.half.len)
  phat_a.first <- sum((var1.first > pred.q1.second) & (var2.first > pred.q2.second)) * (1 / first.half.len)
  ptilde_b.second <- sum((var1.second < pred.q1.first) & (var2.second < pred.q2.first)) * (1 / second.half.len)
  phat_a.second <- sum((var1.second > pred.q1.first) & (var2.second > pred.q2.first)) * (1 / second.half.len)

  zstat_b.first <- (ptilde_b.first - tau^2) / sqrt( tau^2 * (1-tau)^2 / (first.half.len) )
  zstat_a.first <- (phat_a.first - (1 - tau)^2) / sqrt( tau^2 * (1-tau)^2 / (first.half.len) )
  zstat_b.second <- (ptilde_b.second - tau^2) / sqrt( tau^2 * (1-tau)^2 / (second.half.len) )
  zstat_a.second <- (phat_a.second - (1 - tau)^2) / sqrt( tau^2 * (1-tau)^2 / (second.half.len) )

  zstat.above <- ((zstat_a.first + zstat_a.second) / sqrt(2))
  zstat.below <- ((zstat_b.first + zstat_b.second) / sqrt(2))

  p_val_below <- 2 * pnorm(abs(zstat.below), lower.tail = FALSE)
  p_val_above <- 2 * pnorm(abs(zstat.above), lower.tail = FALSE)

  if (tau < 0.5) {
    return(p_val_below)
  }
  else {
    return(p_val_above)
  }
}
