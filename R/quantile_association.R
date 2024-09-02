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
adjacency.similarity <- function(actual, predicted, precision=TRUE, recall=FALSE, hamming=FALSE, undirected=TRUE) {
  if(length(actual) != length(predicted)) {
    return("Vectors need to be of same length")
  }

  if(undirected) { # Don't double count edges
    actual.mat <- matrix(actual, nrow = sqrt(length(actual)), ncol = sqrt(length(actual)))
    actual <- actual.mat[upper.tri(actual.mat, diag = FALSE)]

    pred.mat <- matrix(predicted, nrow = sqrt(length(predicted)), ncol = sqrt(length(predicted)))
    predicted <- pred.mat[upper.tri(pred.mat, diag = FALSE)]
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
calculate.skeleton <- function(data, tau, m.max=Inf, quacc=TRUE, linear=TRUE, correl=FALSE, verbose=FALSE, adj_vector=FALSE, alpha=0.05) {
  saveRDS(tau, "tau.rds")
  if(quacc) {
    if(linear) {
      pc_graph <- pcalg::skeleton(data, indepTest = linear.quacc, labels = colnames(data), alpha = alpha, verbose = verbose, NAdelete=FALSE, m.max=m.max)
    }
    else{
      pc_graph <- pcalg::skeleton(data, indepTest = rf.quacc, labels = colnames(data), alpha = alpha, verbose = verbose, NAdelete=FALSE, m.max=m.max)
    }
  }
  else {
    if(correl) {
      suffStat <- list(C = cor(data), n = nrow(data))
      pc_graph <- pcalg::skeleton(suffStat, indepTest = pcalg::gaussCItest, labels = colnames(data), alpha = alpha, verbose = verbose, m.max=m.max)
    }
    else{
      pc_graph <- pcalg::skeleton(data, indepTest = orig.quantile.ztest, labels = colnames(data), alpha = alpha, verbose = verbose, m.max=m.max)
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
#' @param filter If we want to remove/filter out certain datapoints
#' @param filt.indices Which data point indices we want to keep if we choose to filter
#' @return Jittered version of your input data frame
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


#' QuACC when normalized between -1 and 1 using linear regression estimators.
#'
#' @param x Index of a column
#' @param y Index of a column (not equal to x)
#' @param S Conditioning set to be used to determine if conditional independence exists (can be empty set)
#' @param suffStat The dataframe of data (there is no sufficient statistic for this calculation)
#' @return pvalue corresponding to if the quantile level
linear.quacc.rho <- function(x, y, S, suffStat) {

  quacc.rho.singular <- function(x, y, S, tau, data, train.indices=1:(n%/%2)) {

    data.train <- data[train.indices, ]
    data.test <- data[-train.indices, ]
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

  tau <- readRDS("tau.rds")
  k <- 5
  quacc.vals <- rep(0, k)

  data <- suffStat
  complete.columns <- c(x, y, S)
  complete_cases_indices <- complete.cases(data[, complete.columns])
  data <- data[complete_cases_indices, ]

  set.seed(123)
  data <- data[sample(nrow(data)),] # Randomly shuffle with fixed seed
  set.seed(NULL)

  n <- length(data[,1])
  folds <- cut(seq(1, nrow(data)), breaks=k, labels=FALSE)

  for(i in 1:k) {
    fold.indices <- which(folds==i, arr.ind=TRUE)
    quacc.vals[i] <- quacc.rho.singular(x, y, S, tau, data, train.indices=fold.indices)
  }
  quacc <- mean(quacc.vals)
  #quacc <- median(quacc.vals)
  return(quacc)

}

#' General function for getting QuACC of dataset at particular tau level.
#' Not to be used with PC algorithm.
#'
#' @param x Index of a column
#' @param y Index of a column (not equal to x)
#' @param S Conditioning set to be used to determine if conditional independence exists (can be empty set)
#' @param data The dataframe of data (there is no sufficient statistic for this calculation)
#' @param tau The quantile level of interest
#' @return pvalue corresponding to if the quantile level
general.linear.quacc.rho <- function(x, y, S, data, tau, train.indices) {

  n <- length(data[,1])
  data.train <- data[train.indices, ]
  data.test <- data[-train.indices, ]
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



#' General function for getting QuACC of dataset at particular tau level.
#' Not to be used with PC algorithm.
#'
#' @param x Index of a column
#' @param y Index of a column (not equal to x)
#' @param S Conditioning set to be used to determine if conditional independence exists (can be empty set)
#' @param data The dataframe of data (there is no sufficient statistic for this calculation)
#' @param tau The quantile level of interest
#' @return pvalue corresponding to if the quantile level
general.linear.quacc <- function(x, y, S, data, tau, train.indices) {

  n <- length(data[,1])
  data.train <- data[train.indices, ]
  data.test <- data[-train.indices, ]
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

    var.term <- tau^2 * (1 - tau)^2

    filt.indices.var1 <- which(var2.test < fit.q2)
    filt.indices.var2 <- which(var1.test < fit.q1)

    C <- as.matrix(F.ecdf(fit.q1))
    D <- as.matrix(G.ecdf(fit.q2))
  }
  else{
    c.above <- sum((var1.test > fit.q1) & (var2.test > fit.q2)) / n.test
    quacc <- c.above - (1 - tau)^2

    var.term <- (1 - tau)^2 * (1 - (1 - tau))^2

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

  quacc.var <- ( (var.term) + (t(kappa.var1) %*% sigma1 %*% kappa.var1)[1] + (t(kappa.var2) %*% sigma2 %*% kappa.var2)[1] ) / (n.test)
  return(c(quacc, quacc.var))
}


#' Tests QuACC given linear quantile regression estimators is S is empty otherwise uses honest random forests
#' and then returns a p-value.
#'
#' @param x Index of a column
#' @param y Index of a column (not equal to x)
#' @param S Conditioning set to be used to determine if conditional independence exists (can be empty set)
#' @param suffStat The dataframe of data (there is no sufficient statistic for this calculation)
#' @return pvalue corresponding to if the quantile level
rf.quacc <- function(x, y, S, suffStat) {

  koenker.sandwich <- function(tau, x, y, hs=TRUE, filter=FALSE, filt.indices=c(0)) {
    # Get constants
    eps <- .Machine$double.eps^(1/2)
    tau <- tau
    n <- length(y)

    # Check for valid h
    h <- quantreg::bandwidth.rq(tau, n, hs = hs)
    while((tau - h < 0) || (tau + h > 1)) h <- h/2

    # Calculate Hendricks-Koenker sandwich
    if(filter){ # Filter for y, the regressor

      x.filt <- x[filt.indices, , drop = FALSE]
      y.filt <- y[filt.indices]

      bhi <- quantreg::rq.fit(x.filt, y.filt, tau = tau + h)$coef
      blo <- quantreg::rq.fit(x.filt, y.filt, tau = tau - h)$coef
    }
    else{ # Assume independence in density of var y and another variable that dictates filt.indices
      bhi <- quantreg::rq.fit(x, y, tau = tau + h)$coef
      blo <- quantreg::rq.fit(x, y, tau = tau - h)$coef
    }

    dyhat <- as.matrix(x) %*% (bhi - blo)
    f <- pmax(0, (2 * h)/(dyhat - eps))
    return(f)
  }

  singular.quacc <- function(x, y, S, tau, data, train.indices) {
    data.train <- data[train.indices, ]
    data.test <- data[-train.indices, ]
    n.train <- length(data.train[,1])
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

      fit.q1 <- predict(q1, newdata=data.test[, S, drop=FALSE])
      fit.q2 <- predict(q2, newdata=data.test[, S, drop=FALSE])
    }
    else {
      q1 <- grf::quantile_forest(X=data.train[, S, drop=FALSE],
                                 Y=data.train[, x],
                                 honesty = TRUE,
                                 quantiles = c(tau), min.node.size=50)

      q2 <- grf::quantile_forest(X=data.train[, S, drop=FALSE],
                                 Y=data.train[, y],
                                 honesty = TRUE,
                                 quantiles = c(tau), min.node.size=50)

      fit.q1 <- predict(q1, data.test[, S, drop=FALSE], quantiles = c(tau))$predictions
      fit.q2 <- predict(q2, data.test[, S, drop=FALSE], quantiles = c(tau))$predictions
    }

    # CDF estimations
    F.ecdf <- ecdf(var1.test)
    G.ecdf <- ecdf(var2.test)

    # Calculate QuACC and normalize
    if(tau < 0.5) {
      c.below <- sum((var1.test < fit.q1) & (var2.test < fit.q2)) / n.test
      quacc <- c.below - tau^2

      var.term <- tau^2 * (1 - tau)^2

      filt.indices.var1 <- which(var2.test < fit.q2)
      filt.indices.var2 <- which(var1.test < fit.q1)

      C <- as.matrix(F.ecdf(fit.q1))
      D <- as.matrix(G.ecdf(fit.q2))
    }
    else{
      c.above <- sum((var1.test > fit.q1) & (var2.test > fit.q2)) / n.test
      quacc <- c.above - (1 - tau)^2

      var.term <- (1 - tau)^2 * (1 - (1 - tau))^2

      filt.indices.var1 <- which(var2.test > fit.q2)
      filt.indices.var2 <- which(var1.test > fit.q1)

      C <- 1 - as.matrix(F.ecdf(fit.q1))
      D <- 1 - as.matrix(G.ecdf(fit.q2))
    }
    padded.z <- as.matrix(cbind(rep(1, n.test), data.test[, S, drop=FALSE]))

    # Density estimations
    A <- as.matrix( diag( koenker.sandwich(tau, x=padded.z, y=var1.test, filter=FALSE, filt.indices=filt.indices.var1) ) )
    B <- as.matrix( diag( koenker.sandwich(tau, x=padded.z, y=var2.test, filter=FALSE, filt.indices=filt.indices.var2) ) )
    A.vec <- as.matrix(diag(A))
    B.vec <- as.matrix(diag(B))

    # Compute kappa terms in QuACC
    kappa.var1 <- (1 / n.test) * t(A.vec) %*% C
    kappa.var2 <- (1 / n.test) * t(B.vec) %*% D

    # Compute terms in variance of RF estimator
    w1 <- apply(grf::get_forest_weights(q1), 2, sum) # Sum over the training
    w1.rq <- rq(data.train[, x] ~ as.matrix(data.train[, S, drop=FALSE]), weights=w1)
    resid1 <- residuals(w1.rq)

    w2 <- apply(grf::get_forest_weights(q2), 2, sum)
    w2.rq <- quantreg::rq(data.train[, y] ~ as.matrix(data.train[, S, drop=FALSE]), weights=w2)
    resid2 <- residuals(w2.rq)

    check_function <- function(u, tau) {
      tau * (u >= 0) - (1 - tau) * (u < 0)
    }

    loss.H1 <- check_function(resid1, tau)
    loss.H2 <- check_function(resid2, tau)

    H1 <- var(loss.H1)
    H2 <- var(loss.H2)

    # Compute QuACC
    s <- 0.5 * n.train # Default sample.fraction parameter in grf
    sigma1 <- (n.train / s) * (H1 / (t(A.vec) %*% A.vec))
    sigma2 <- (n.train / s) * (H2 / (t(B.vec) %*% B.vec))

    quacc.var <- ((var.term) + (t(kappa.var1) %*% sigma1 %*% kappa.var1)[1] + (t(kappa.var2) %*% sigma2 %*% kappa.var2)[1]) / (n.test)
    return(c(quacc, quacc.var))
  }

  tau <- readRDS("tau.rds")
  k <- 5
  quacc.vals <- rep(0, k)
  quacc.vars <- rep(0, k)

  data <- suffStat
  complete.columns <- c(x, y, S)
  complete_cases_indices <- complete.cases(data[, complete.columns])
  data <- data[complete_cases_indices, ]

  set.seed(123)
  data <- data[sample(nrow(data)),] # Randomly shuffle with fixed seed
  set.seed(NULL)

  n <- length(data[,1])
  folds <- cut(seq(1, nrow(data)), breaks=k, labels=FALSE)

  for(i in 1:k) {
    fold.indices <- which(folds!=i, arr.ind=TRUE) # Train on all but kth fold, evaluate on fold
    fold.res <- singular.quacc(x, y, S, tau, data, train.indices=fold.indices)
    quacc.vals[i] <- fold.res[1]
    quacc.vars[i] <- fold.res[2]
  }
  quacc <- sum(quacc.vals) / sqrt( sum(quacc.vars) )

  # Calculate p-value of QuACC
  p_val <- 2 * pnorm(abs(quacc), lower.tail = FALSE)
  return(quacc)
}



#' Tests QuACC given linear quantile regression estimators and returns a p-value.
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

    dyhat <- x %*% (bhi - blo)
    f <- pmax(0, (2 * h)/(dyhat - eps))
    return(f)
  }

  singular.quacc <- function(x, y, S, tau, data, train.indices) {
    data.train <- data[train.indices, ]
    data.test <- data[-train.indices, ]
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
      c.below <- mean( (var1.test < fit.q1) & (var2.test < fit.q2) )
      quacc <- c.below - tau^2

      var.term <- tau^2 * (1 - tau)^2

      filt.indices.var1 <- which(var2.test < fit.q2)
      filt.indices.var2 <- which(var1.test < fit.q1)

      C <- as.matrix(F.ecdf(fit.q1))
      D <- as.matrix(G.ecdf(fit.q2))
    }
    else{
      c.above <- mean( (var1.test > fit.q1) & (var2.test > fit.q2) )
      quacc <- c.above - (1 - tau)^2

      var.term <- (1 - tau)^2 * (1 - (1 - tau))^2

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

    quacc.var <- ((var.term) + (t(kappa.var1) %*% sigma1 %*% kappa.var1)[1] + (t(kappa.var2) %*% sigma2 %*% kappa.var2)[1]) / (n.test)
    return(c(quacc, quacc.var))
  }

  tau <- readRDS("tau.rds")
  k <- 5
  quacc.vals <- rep(0, k)
  quacc.vars <- rep(0, k)

  data <- suffStat
  complete.columns <- c(x, y, S)
  complete_cases_indices <- complete.cases(data[, complete.columns])
  data <- data[complete_cases_indices, ]

  set.seed(123)
  data <- data[sample(nrow(data)),] # Randomly shuffle with fixed seed
  set.seed(NULL)

  n <- length(data[,1])
  folds <- cut(seq(1, nrow(data)), breaks=k, labels=FALSE)

  for(i in 1:k) {
    fold.indices <- which(folds!=i, arr.ind=TRUE) # Train on all but kth fold, evaluate on fold
    fold.res <- singular.quacc(x, y, S, tau, data, train.indices=fold.indices)
    quacc.vals[i] <- fold.res[1]
    quacc.vars[i] <- fold.res[2]
  }
  quacc <- sum(quacc.vals) / sqrt( sum(quacc.vars) )
  #quacc.var.estimate <- median(quacc.vars + (quacc.vals - mean(quacc.vals))^2)
  #quacc <- median(quacc.vals) / sqrt( quacc.var.estimate )

  # Calculate p-value of QuACC
  p_val <- 2 * pnorm(abs(quacc), lower.tail = FALSE)
  return(p_val)
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
#' @param df Input dataframe of data
#' @param tau A particular quantile level (0 to 1, not inclusive)
#' @param type If weights of the regression should be marginal or conditional on all other variables
#' @param quacc If you should use the linear QuACC or use the standard non train test split statistic
#' @return An n by n matrix (where n is the number of columns in data) that contains
#' the marginal relationships of each pair of columns
pairwise.test <- function(data, tau, weights="marginal", quacc=TRUE, rho=FALSE) {
  num_cols <- length((colnames(data)))
  quacc.table <- matrix(0, nrow=num_cols, ncol=num_cols)

  for(x in 1:num_cols) {
    for(y in 1:x) {
      if(y != x) {
        if(quacc) {

          data.subset <- data
          if(weights == "marginal"){
            S <- c()
            complete.columns <- c(x, y)
            complete_cases_indices <- complete.cases(data.subset[, complete.columns])
            data.subset <- data.subset[complete_cases_indices, ]
          }
          else {
            S <- 1:num_cols
            S <- S[-c(x, y)]
            complete.columns <- c(x, y, S)
            complete_cases_indices <- complete.cases(data.subset[, complete.columns])
            data.subset <- data.subset[complete_cases_indices, ]
          }

          k <- 5
          quacc.vals <- rep(0, k)
          quacc.vars <- rep(0, k)

          complete.columns <- c(x, y, S)
          complete_cases_indices <- complete.cases(data.subset[, complete.columns])
          data.subset <- data.subset[complete_cases_indices, ]

          set.seed(123)
          data.subset <- data.subset[sample(nrow(data.subset)),] # Randomly shuffle with fixed seed
          set.seed(NULL)

          n <- length(data.subset[,1])
          folds <- cut(seq(1, nrow(data.subset)), breaks=k, labels=FALSE)

          if(rho) {
            for(i in 1:k) {
              fold.indices <- which(folds!=i, arr.ind=TRUE)
              quacc.vals[i] <- general.linear.quacc.rho(x=x, y=y, S=S, tau=tau, data=data.subset, train.indices=fold.indices)
            }
            quacc.table[x, y] <- mean(quacc.vals)
          }
          else{
            for(i in 1:k) {
              fold.indices <- which(folds!=i, arr.ind=TRUE)
              fold.res <- general.linear.quacc(x=x, y=y, S=S, tau=tau, data=data.subset, train.indices=fold.indices)
              quacc.vals[i] <- fold.res[1]
              quacc.vars[i] <- fold.res[2]
            }
            quacc.table[x, y] <- sum(quacc.vals) / sqrt( sum(quacc.vars) )
            #quacc.var.estimate <- median(quacc.vars + (quacc.vals - mean(quacc.vals))^2)
            #quacc.table[x, y] <- median(quacc.vals) / sqrt( quacc.var.estimate )
          }

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

          if(tau < 0.5) {
            quacc.table[x, y] <- (ptilde - tau^2) / sqrt( tau^2 * (1-tau)^2 / n)
          }
          else{
            quacc.table[x, y] <- (phat - (1 - tau)^2) / sqrt( tau^2 * (1-tau)^2 / n)
          }
        }
      }

    }
  }

  # Symmetrize the matrix
  for(x in 1:num_cols) {
    for(y in 1:x) {
      quacc.table[y, x] <- quacc.table[x, y]
    }
  }
  return( quacc.table )
}

