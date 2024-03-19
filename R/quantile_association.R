#' Adds jitter to the input dataframe to break ties.
#'
#' Goes column by column in the dataframe and adds jitter from the distribution
#'
#' Unif(-a, a)
#'
#' where a = factor * d /5 (and d is the minimum distance between any two points in a vector).
#'
#' @param data Input data frame of data
#' @param factor A float to scale the jittering by (if greater than 2.5, there is a risk of quantile crossings)
#' @return Jittered version of your input data frame
jitter.columns <- function(data, factor=0.1) {
  num_cols <- length((colnames(data)))

  for(x in 1:num_cols) {
    data[,x] <- jitter(data[,x], factor=factor)
  }
  return(data)
}

#' Calculates one of three similarity metrics for use in evaluating graph prediction.
#'
#' @param actual Actual adjacency matrix as a vector (binary values only)
#' @param predicted Predicted adjacency matrix as a vector (binary values only)
#' @param precision Precision between actual/predicted
#' @param recall Recall of edge presence
#' @param hamming Hamming distance between the two
#' @return One of the three similarity metrics above
adjacency.similarity <- function(actual, predicted, precision=TRUE, recall=FALSE, hamming=FALSE) {
  if(length(actual) != length(predicted)) {
    return("Vectors need to be of same length")
  }
  if(precision) {
    if(sum(predicted) == 0) {
      return('Divide by 0 error')
    }
    res <- sum(actual & predicted) / sum(predicted)
  }

  if(recall) {
    if(sum(actual) == 0) {
      return('Divide by 0 error')
    }
    res <- sum(actual & predicted) / sum(actual)
  }
  if(hamming) {
    res <- mean(actual != predicted)
  }

  return(res)
}

#' Calculates skeleton of QGM graph based on different independence tests.
#'
#' @param data Input data frame of data
#' @param tau Quantile level of interest (from 0 to 1, non-inclusive)
#' @param quacc Use QuACC hypothesis test
#' @param correl Use Gaussian CI test
#' @param verbose Print verbose independence tests results
#' @param adj_vector Return vector form of adjacency matrix
#' @return Graph or adjacency vector of underlying relationships
calculate.skeleton <- function(data, tau, quacc=TRUE, correl=FALSE, verbose=FALSE, adj_vector=FALSE) {
  saveRDS(tau, "tau.rds")
  if(quacc) {
    pc_graph <- pcalg::skeleton(data, indepTest = linear.quacc, labels = colnames(data), alpha = 0.05, verbose = verbose, NAdelete=FALSE)
  }
  else {
    if(correl) {
      suffStat <- list(C = cor(data), n = nrow(data))
      pc_graph <- pcalg::skeleton(suffStat, indepTest = pcalg::gaussCItest, labels = colnames(data), alpha = 0.05, verbose = verbose)
    }
    else{
      pc_graph <- pcalg::skeleton(data, indepTest = orig.quantile.ztest, labels = colnames(data), alpha = 0.05, verbose = verbose)
    }
  }
  unlink("tau.rds")
  if(adj_vector) {
    adj.mat <- as(pc_graph, "amat")
    return(as.vector(adj.mat))
  }
  return(pc_graph)
}

#' Calculate density of a particular variable using the Hendricks-Koenker sandwich.
#'
#' @param rq.object Quantile regression object (currently supports linear quantile regression)
#' @param x Not the index, but the actual values of column x with a padded intercept
#' @param y Not the index, but the actual values of column y with a padded intercept
#' @param hs Which way to calculate the bandwidth of quantiles
#' @return Jittered version of your input data frame
koenker.sandwich <- function(rq.object, x, y, hs=TRUE) {
  # Get constants
  eps <- .Machine$double.eps^(1/2)
  tau <- rq.object$tau
  n <- length(y)

  # Check for valid h
  h <- quantreg::bandwidth.rq(tau, n, hs = hs)
  while((tau - h < 0) || (tau + h > 1)) h <- h/2

  # Calculate Hendricks-Koenker sandwich
  bhi <- quantreg::rq.fit(x, y, tau = tau + h, method = rq.object$method)$coef
  blo <- quantreg::rq.fit(x, y, tau = tau - h, method = rq.object$method)$coef

  dyhat <- as.matrix(x) %*% (bhi - blo)
  f <- pmax(0, (2 * h)/(dyhat - eps))
  return(f)
}


#' QuACC when normalized between -1 and 1 using linear regression estimators.
#'
#' @param x Index of a column
#' @param y Index of a column (not equal to x)
#' @param S Conditioning set to be used to determine if conditional independence exists (can be empty set)
#' @param suffStat The dataframe of data (there is no sufficient statistic for this calculation)
#' @return pvalue corresponding to if the quantile level
linear.quacc.rho <- function(x, y, S, suffStat) {

  tau <- readRDS("tau.rds")
  data <- suffStat

  n <- length(data[,1])
  data.train <- data[1:(n%/%2), ]
  data.test <- data[(1 + (n%/%2)):n, ]
  n.test <- length(data.test[,1])
  var1.test <- data.test[, x]
  var2.test <- data.test[, y]

  if(length(S) == 0){
    q1 <- quantreg::rq(as.formula(paste(colnames(data)[x], "~1")),
                       data.train,
                       tau=tau)

    q2 <- quantreg::rq(as.formula(paste(colnames(data)[y], "~1")),
                       data.train,
                       tau=tau)
  }
  else {
    q1 <- quantreg::rq(as.formula(paste(colnames(data)[x], "~",
                                        paste(colnames(data)[S], collapse = "+"), sep = "")),
                       data.train,
                       tau=tau)

    q2 <- quantreg::rq(as.formula(paste(colnames(data)[y], "~",
                                        paste(colnames(data)[S], collapse = "+"), sep = "")),
                       data.train,
                       tau=tau)
  }

  fit.q1 <- predict(q1, newdata=data.test[, S, drop=FALSE])
  fit.q2 <- predict(q2, newdata=data.test[, S, drop=FALSE])


  # Calculate QuACC and normalize
  if(tau <= 0.5) {
    c.below <- sum((var1.test < fit.q1) & (var2.test < fit.q2)) / n.test
    quacc <- c.below - tau^2

    if(quacc > 0) {
      quacc <- quacc / (tau - tau^2)
    }
    else{
      quacc <- quacc / tau^2
    }
  }
  else{
    c.above <- sum((var1.test > fit.q1) & (var2.test > fit.q2)) / n.test
    quacc <- c.above - (1 - tau)^2
    if(quacc > 0) {
      quacc <- quacc / ((1 - tau) - (1 - tau)^2)
    }
    else{
      quacc <- quacc / ((1 - tau)^2)
    }
  }

  return(quacc)
}


#' Tests QuACC given linear quantile regression estimators.
#'
#' @param x Index of a column
#' @param y Index of a column (not equal to x)
#' @param S Conditioning set to be used to determine if conditional independence exists (can be empty set)
#' @param suffStat The dataframe of data (there is no sufficient statistic for this calculation)
#' @return pvalue corresponding to if the quantile level
linear.quacc <- function(x, y, S, suffStat) {

  koenker.sandwich <- function(rq.object, x, y, hs=TRUE, filter=FALSE, filt.indices=c(0)) {
    # Get constants
    eps <- .Machine$double.eps^(1/2)
    tau <- rq.object$tau
    n <- length(y)

    # Check for valid h
    h <- quantreg::bandwidth.rq(tau, n, hs = hs)
    while((tau - h < 0) || (tau + h > 1)) h <- h/2

    # Calculate Hendricks-Koenker sandwich
    if(filter){ # Filter for y, the regressor

      x.filt <- x[filt.indices, , drop = FALSE]
      y.filt <- y[filt.indices]

      bhi <- quantreg::rq.fit(x.filt, y.filt, tau = tau + h, method = rq.object$method)$coef
      blo <- quantreg::rq.fit(x.filt, y.filt, tau = tau - h, method = rq.object$method)$coef
    }
    else{ # Assume independence in density of var y and another variable that dictates filt.indices
      bhi <- quantreg::rq.fit(x, y, tau = tau + h, method = rq.object$method)$coef
      blo <- quantreg::rq.fit(x, y, tau = tau - h, method = rq.object$method)$coef
    }

    dyhat <- as.matrix(x) %*% (bhi - blo)
    f <- pmax(0, (2 * h)/(dyhat - eps))
    return(f)
  }

  tau <- readRDS("tau.rds")
  data <- suffStat

  complete.columns <- c(x, y, S)
  complete_cases_indices <- complete.cases(data[, complete.columns])
  data <- data[complete_cases_indices, ]

  n <- length(data[,1])
  data.train <- data[1:(n%/%2), ] # Split in half train, half test
  data.test <- data[(1 + (n%/%2)):n, ]
  n.test <- length(data.test[,1])
  var1.test <- data.test[, x]
  var2.test <- data.test[, y]

  if(length(S) == 0){
    q1 <- quantreg::rq(as.formula(paste(colnames(data)[x], "~1")),
                       data.train,
                       tau=tau)

    q2 <- quantreg::rq(as.formula(paste(colnames(data)[y], "~1")),
                       data.train,
                       tau=tau)
  }
  else {
    q1 <- quantreg::rq(as.formula(paste(colnames(data)[x], "~",
                                        paste(colnames(data)[S], collapse = "+"), sep = "")),
                       data.train,
                       tau=tau)

    q2 <- quantreg::rq(as.formula(paste(colnames(data)[y], "~",
                                        paste(colnames(data)[S], collapse = "+"), sep = "")),
                       data.train,
                       tau=tau)
  }

  fit.q1 <- predict(q1, newdata=data.test[, S, drop=FALSE])
  fit.q2 <- predict(q2, newdata=data.test[, S, drop=FALSE])

  # CDF estimations
  F.ecdf <- ecdf(var1.test)
  G.ecdf <- ecdf(var2.test)

  # Calculate QuACC and normalize
  if(tau < 0.5) {
    c.below <- sum((var1.test < fit.q1) & (var2.test < fit.q2)) / n.test
    quacc <- c.below - tau^2

    var.term <- tau^2 * (1 - tau^2)

    filt.indices.var1 <- which(var2.test < fit.q2)
    filt.indices.var2 <- which(var1.test < fit.q1)

    C <- as.matrix(F.ecdf(fit.q1))
    D <- as.matrix(G.ecdf(fit.q2))
  }
  else{
    c.above <- sum((var1.test > fit.q1) & (var2.test > fit.q2)) / n.test
    quacc <- c.above - (1 - tau)^2

    var.term <- (1 - tau)^2 * (1 - (1 - tau)^2)

    filt.indices.var1 <- which(var2.test > fit.q2)
    filt.indices.var2 <- which(var1.test > fit.q1)

    C <- 1 - as.matrix(F.ecdf(fit.q1))
    D <- 1 - as.matrix(G.ecdf(fit.q2))
  }

  s1 <- summary(q1, cov=TRUE, se='nid')
  s2 <- summary(q2, cov=TRUE, se='nid')
  padded.z <- as.matrix(cbind(rep(1, n.test), data.test[, S, drop=FALSE]))

  # Density estimations
  A <- as.matrix(diag(koenker.sandwich(q1, x=padded.z, y=var1.test, filter=FALSE, filt.indices=filt.indices.var1)))
  B <- as.matrix(diag(koenker.sandwich(q2, x=padded.z, y=var2.test, filter=FALSE, filt.indices=filt.indices.var2)))

  # Compute variance terms in QuACC
  kappa.var1 <- (1 / n.test) * t(padded.z) %*% A %*% C
  kappa.var2 <- (1 / n.test) * t(padded.z) %*% B %*% D

  # Compute QuACC
  sigma1 <- tau * (1 - tau) * s1$Hinv * s1$J * s1$Hinv
  sigma2 <- tau * (1 - tau) * s2$Hinv * s2$J * s2$Hinv

  quacc <- quacc / sqrt( (var.term / n.test) + (t(kappa.var1) %*% sigma1 %*% kappa.var1)[1] + (t(kappa.var2) %*% sigma2 %*% kappa.var2)[1] )

  # Calculate p-value of QuACC
  p_val <- 2 * pnorm(abs(quacc), lower.tail = FALSE)
  return(p_val)
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
orig.quantile.ztest <- function (x, y, S, suffStat) {
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
#' @param quacc If you should use the linear QuACC or use the standard non train test split statistic
#' @return An n by n matrix (where n is the number of columns in data) that contains
#' the marginal relationships of each pair of columns
pairwise.test <- function(data, tau, weights="marginal", quacc=TRUE) {
  num_cols <- length((colnames(data)))

  zstat_above <- matrix(0, nrow=num_cols, ncol=num_cols)
  zstat_below <- matrix(0, nrow=num_cols, ncol=num_cols)

  for(x in 1:num_cols) {
    for(y in 1:x) {

      n <- length(data[,x])

      if(quacc) {
        data.train <- data[1:(n%/%2), ] # Split in half train, half test
        data.test <- data[(1 + (n%/%2)):n, ]
        n.test <- length(data.test[,1])
        var1.test <- data.test[, x]
        var2.test <- data.test[, y]
        col1 <- colnames(data)[x]
        col2 <- colnames(data)[y]

        if(weights == "marginal"){
          q1 <- quantreg::rq(as.formula(paste0(col1, " ~ 1")), tau=tau, data=data)
          q2 <- quantreg::rq(as.formula(paste0(col2, " ~ 1")), tau=tau, data=data)

          S <- c()
          fit.q1 <- predict(q1, newdata=data.test[, S, drop=FALSE])
          fit.q2 <- predict(q2, newdata=data.test[, S, drop=FALSE])

          padded.z <- as.matrix(rep(1, n.test))
        }
        else {
          S <- 1:num_cols
          S <- S[-c(x, y)]

          q1 <- quantreg::rq(as.formula(paste(colnames(data)[x], "~",
                                              paste(colnames(data)[S], collapse = "+"), sep = "")),
                             data.train,
                             tau=tau)

          q2 <- quantreg::rq(as.formula(paste(colnames(data)[y], "~",
                                              paste(colnames(data)[S], collapse = "+"), sep = "")),
                             data.train,
                             tau=tau)

          fit.q1 <- predict(q1, newdata=data.test[, S, drop=FALSE])
          fit.q2 <- predict(q2, newdata=data.test[, S, drop=FALSE])

          padded.z <- as.matrix(cbind(rep(1, n.test), data.test[, S, drop=FALSE]))
        }

        # Calculate QuACC and normalize
        c.below <- sum((var1.test < fit.q1) & (var2.test < fit.q2)) / n.test
        quacc.low <- c.below - tau^2
        c.above <- sum((var1.test > fit.q1) & (var2.test > fit.q2)) / n.test
        quacc.high <- c.above - (1 - tau)^2

        s1 <- summary(q1, cov=TRUE, se='nid')
        s2 <- summary(q2, cov=TRUE, se='nid')

        # Density estimations
        A <- as.matrix(diag(koenker.sandwich(q1, x=padded.z, y=var1.test))) # Assume independence between 2 for now
        B <- as.matrix(diag(koenker.sandwich(q2, x=padded.z, y=var2.test))) # Assume independence between 2 for now

        # CDF estimations
        F.ecdf <- ecdf(var1.test)
        G.ecdf <- ecdf(var2.test)
        C <- as.matrix(F.ecdf(fit.q1))
        D <- as.matrix(G.ecdf(fit.q2))

        kappa.var1 <- (1 / n.test) * t(padded.z) %*% A %*% C
        kappa.var2 <- (1 / n.test) * t(padded.z) %*% B %*% D

        sigma1 <- tau * (1 - tau) * s1$Hinv * s1$J * s1$Hinv
        sigma2 <- tau * (1 - tau) * s2$Hinv * s2$J * s2$Hinv

        zstat_below[x, y] <- quacc.low / sqrt( ((tau^2 * (1 - tau^2)) / n.test + (t(kappa.var1) %*% sigma1 %*% kappa.var1)[1] + (t(kappa.var2) %*% sigma2 %*% kappa.var2)[1]) )
        zstat_above[x, y] <- quacc.high / sqrt( ((tau^2 * (1 - tau^2)) / n.test + (t(kappa.var1) %*% sigma1 %*% kappa.var1)[1] + (t(kappa.var2) %*% sigma2 %*% kappa.var2)[1]) )

      }
      else { # Don't use QuACC
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

  for(x in 1:num_cols) {
    for(y in 1:x) {
      zstat_below[y, x] <- zstat_below[x, y]
      zstat_above[y, x] <- zstat_above[x, y]
    }
  }

  if(tau >= 0.5){
    return( zstat_above )
  }
  return( zstat_below )
}

