#' Calculate Objective Function Value
#'
#' @param beta A numeric vector of coefficients (length p).
#' @param Pi A numeric matrix (n x T) for the low-rank component.
#' @param y A numeric matrix (n x T) of observations.
#' @param x A numeric array (n x T x p) of predictors.
#' @param u_quantile The quantile level (scalar).
#' @return The scalar value of the objective function part.
objective_function_r <- function(beta, Pi, y, x, u_quantile) {
  n <- dim(y)[1]
  T_dim <- dim(y)[2]
  # p <- dim(x)[3] # p can be inferred from length(beta)

  # Calculate X %*% beta term: sum_k x[i,j,k] * beta[k]
  # This should result in an n x T matrix.
  # x is n x T x p. beta is p.
  # More efficient calculation:
  X_beta_prod <- apply(sweep(x, MARGIN = 3, STATS = as.vector(beta), FUN = "*"), MARGIN = c(1, 2), FUN = sum)

  er <- y - X_beta_prod - Pi
  rho_er <- u_quantile * er + pmax(-er, 0) # pmax is element-wise max
  val <- sum(rho_er) # sum over all elements of the matrix

  return(val)
}

#' Nuclear Norm Calculation
#'
#' @param mat A numeric matrix.
#' @return The nuclear norm (sum of singular values).
nuclear_norm_r <- function(mat) {
  if (any(is.na(mat)) || any(is.infinite(mat))) {
    return(Inf) # Or handle error appropriately
  }
  svd_decomp <- svd(mat)
  return(sum(svd_decomp$d))
}


#' Principal Component Analysis using SVD
#'
#' @param data An M x N matrix of input data (M features, N samples).
#' @return A list containing:
#'   \item{Newdata}{Data projected onto PCA space (PCs are rows).}
#'   \item{PCASpace}{Eigenvectors (principal components, PCs are columns).}
#'   \item{EigValues}{Eigenvalues.}
pcasvd_r <- function(data) {
  r <- nrow(data) # Number of features (M)
  c <- ncol(data) # Number of samples (N)

  # Compute the mean of each feature (row means)
  m <- rowMeans(data)

  # Subtract the mean from each feature (center the data)
  # d <- data - matrix(rep(m, c), nrow = r, ncol = c, byrow = FALSE)
  d <- sweep(data, 1, m, "-") # More idiomatic R

  # Construct the matrix Z
  # Note: In MATLAB d' is c x r. Z is c x r.
  # d is r x c. t(d) is c x r.
  Z <- t(d) / sqrt(c - 1)

  # Calculate SVD of Z ( Z = L S R')
  # R's svd(X) returns list with u, d, v such that X = u diag(d) v'
  # So, L = svd_Z$u, S_vec = svd_Z$d, R = svd_Z$v
  svd_Z <- svd(Z)

  # Eigenvalues of the covariance matrix (d %*% t(d))/(c-1) are svd_Z$d^2
  EigValues <- svd_Z$d^2

  # PCASpace (eigenvectors of covariance matrix) are columns of R (svd_Z$v)
  # However, the MATLAB code assigns PCASpace = R. Let's check dimensions.
  # Z is c x r. svd(Z)$v is r x r. These are the principal directions (eigenvectors of Z'Z).
  # Z'Z = (d/sqrt(c-1)) %*% (t(d)/sqrt(c-1)) = (d %*% t(d))/(c-1), which is r x r covariance matrix.
  # So, svd_Z$v are the eigenvectors of the covariance matrix.
  PCASpace <- svd_Z$v # Eigenvectors are columns

  # Project original centered data onto PCA space
  # Newdata = PCASpace' * d
  # PCASpace is r x r (eigenvectors). t(PCASpace) is r x r. d is r x c.
  # Newdata should be r x c (principal components as rows)
  Newdata <- t(PCASpace) %*% d # Each row is a principal component scores

  return(list(Newdata = Newdata, PCASpace = PCASpace, EigValues = EigValues))
}

#' ADMM for L1 and Nuclear Norm Regularized Quantile Regression
#'
#' @param y Numeric matrix (n x T) of responses.
#' @param x Numeric array (n x T x p) of predictors.
#' @param u_quantile Scalar, the quantile level (e.g., 0.5 for median regression).
#' @param lamb1_orig Scalar, regularization parameter for L1 norm of beta.
#' @param lamb2_orig Scalar, regularization parameter for nuclear norm of Pi.
#' @param tol Scalar, convergence tolerance.
#' @param Niter Integer, maximum number of iterations.
#' @param initial Integer, 0 for cold start, non-zero for warm start (uses provided initial values).
#' @param beta_tilde_init Optional initial value for beta_tilde (1 x p matrix or vector of length p).
#' @param z_beta_init Optional initial value for z_beta (1 x p matrix or vector of length p).
#' @param u_beta_init Optional initial value for u_beta (1 x p matrix or vector of length p).
#' @param pi_tilde_init Optional initial value for pi_tilde (n x T matrix).
#' @param z_pi1_init Optional initial value for z_pi1 (n x T matrix).
#' @param u_pi1_init Optional initial value for u_pi1 (n x T matrix).
#' @param v_init Optional initial value for v (n x T matrix).
#' @param u_v_init Optional initial value for u_v (n x T matrix).
#' @param w_init Optional initial value for w (n x T matrix).
#' @param u_w_init Optional initial value for u_w (n x T matrix).
#' @return A list containing the estimated parameters and intermediate variables.
admm_l1_nuclear_norm_qr_v7_r <- function(
    y, x, u_quantile, lamb1_orig, lamb2_orig, tol, Niter, initial,
    beta_tilde_init = NULL, z_beta_init = NULL, u_beta_init = NULL,
    pi_tilde_init = NULL, z_pi1_init = NULL, u_pi1_init = NULL,
    v_init = NULL, u_v_init = NULL, w_init = NULL, u_w_init = NULL
) {

  n <- dim(y)[1]
  T_dim <- dim(y)[2]
  p <- dim(x)[3]

  lamb1_scaled <- n * T_dim * lamb1_orig
  lamb2_scaled <- n * T_dim * lamb2_orig
  eta <- 1.0 # ADMM penalty parameter

  # Initialization
  if (initial == 0) {
    beta_tilde <- matrix(0, nrow = 1, ncol = p) # Treat as a 1xp row matrix
    z_beta <- matrix(0, nrow = 1, ncol = p)
    u_beta <- matrix(0, nrow = 1, ncol = p)

    pi_tilde <- matrix(0, nrow = n, ncol = T_dim)
    z_pi1 <- matrix(0, nrow = n, ncol = T_dim)
    u_pi1 <- matrix(0, nrow = n, ncol = T_dim)

    v <- matrix(0, nrow = n, ncol = T_dim)
    u_v <- matrix(0, nrow = n, ncol = T_dim)

    w <- matrix(0, nrow = n, ncol = T_dim)
    u_w <- matrix(0, nrow = n, ncol = T_dim)
  } else {
    # Use provided initial values (add checks for NULL if necessary)
    if(is.null(beta_tilde_init) || is.null(z_beta_init) || is.null(u_beta_init) ||
       is.null(pi_tilde_init) || is.null(z_pi1_init) || is.null(u_pi1_init) ||
       is.null(v_init) || is.null(u_v_init) || is.null(w_init) || is.null(u_w_init)) {
      stop("Initial values must be provided for warm start (initial != 0)")
    }
    beta_tilde <- beta_tilde_init
    z_beta <- z_beta_init
    u_beta <- u_beta_init
    pi_tilde <- pi_tilde_init
    z_pi1 <- z_pi1_init
    u_pi1 <- u_pi1_init
    v <- v_init
    u_v <- u_v_init
    w <- w_init
    u_w <- u_w_init
  }

  # Precompute A_tilde2 = X'X (flattened)
  x_reshaped <- matrix(x, nrow = n * T_dim, ncol = p)
  A_tilde2 <- t(x_reshaped) %*% x_reshaped

  obj_val_log <- numeric(Niter)

  calculate_X_beta_prod <- function(current_x, current_beta_vec, n_rows, t_cols) {
    apply(sweep(current_x, MARGIN = 3, STATS = current_beta_vec, FUN = "*"),
          MARGIN = c(1, 2), FUN = sum)
  }

  for (iter in 1:Niter) {
    val_prev <- objective_function_r(as.vector(beta_tilde), pi_tilde, y, x, u_quantile) +
      lamb1_scaled * sum(abs(as.vector(beta_tilde))) +
      lamb2_scaled * nuclear_norm_r(pi_tilde)
    obj_val_log[iter] <- val_prev

    v_prev <- v
    beta_tilde_prev <- beta_tilde
    pi_tilde_prev <- pi_tilde
    z_beta_prev <- z_beta
    z_pi1_prev <- z_pi1
    w_prev <- w

    aux_v <- w - u_v
    v <- pmax(0, aux_v - (u_quantile / eta) * matrix(1, n, T_dim)) +
      pmin(0, aux_v + ((1 - u_quantile) / eta) * matrix(1, n, T_dim))

    A_mat_beta <- A_tilde2 + diag(rep(1, p))
    A_tilde_for_beta <- w + z_pi1 + u_w - y

    b_vec_vals <- numeric(p)
    for (k_ in 1:p) {
      # CORRECTED LINE: Removed drop=FALSE
      b_vec_vals[k_] <- sum(x[,,k_] * A_tilde_for_beta)
    }

    b_vec_rhs <- -b_vec_vals + as.vector(z_beta) + as.vector(u_beta)
    beta_tilde <- matrix(solve(A_mat_beta, as.vector(b_vec_rhs)), nrow = 1)

    Q_pi <- z_pi1 + u_pi1
    tau_pi <- lamb2_scaled / eta

    svd_Q <- svd(Q_pi)
    SQ_diag_shrunk <- pmax(svd_Q$d - tau_pi, 0)

    if (length(SQ_diag_shrunk) == 1 && nrow(svd_Q$u) == 1 && ncol(svd_Q$v) == 1) { # Handle 1x1 Q_pi case
      pi_tilde <- svd_Q$u * SQ_diag_shrunk * t(svd_Q$v)
    } else if (length(SQ_diag_shrunk) == 1) {
      pi_tilde <- svd_Q$u %*% matrix(SQ_diag_shrunk, 1, 1) %*% t(svd_Q$v)
    } else {
      pi_tilde <- svd_Q$u %*% diag(SQ_diag_shrunk) %*% t(svd_Q$v)
    }

    tau_beta <- lamb1_scaled / eta
    q_beta <- as.vector(beta_tilde) - as.vector(u_beta)
    z_beta_vec <- sign(q_beta) * pmax(abs(q_beta) - tau_beta, 0)
    z_beta <- matrix(z_beta_vec, nrow = 1)

    X_beta_prod_term <- calculate_X_beta_prod(x, as.vector(beta_tilde), n, T_dim)

    A_tilde_for_z_w <- -y + X_beta_prod_term + u_w
    B_tilde_for_z_w <- -v - u_v
    C_tilde_for_z_w <- -pi_tilde + u_pi1

    z_pi1 <- (-A_tilde_for_z_w - 2 * C_tilde_for_z_w + B_tilde_for_z_w) / 3.0
    w <- -A_tilde_for_z_w - C_tilde_for_z_w - 2 * z_pi1

    r_pi1 <- -pi_tilde + z_pi1
    u_pi1 <- u_pi1 + r_pi1

    r_beta_vec <- -as.vector(beta_tilde) + as.vector(z_beta)
    u_beta <- matrix(as.vector(u_beta) + r_beta_vec, nrow = 1)

    r_v <- v - w
    u_v <- u_v + r_v

    X_beta_prod_term_for_rw <- calculate_X_beta_prod(x, as.vector(beta_tilde), n, T_dim)
    r_w <- w - y + X_beta_prod_term_for_rw + z_pi1
    u_w <- u_w + r_w

    norm_r_pi1_sq <- sum(r_pi1^2)
    norm_r_beta_sq <- sum(r_beta_vec^2)
    norm_r_v_sq <- sum(r_v^2)
    norm_r_w_sq <- sum(r_w^2)
    primal_residual_norm <- sqrt(norm_r_pi1_sq + norm_r_beta_sq + norm_r_v_sq + norm_r_w_sq)

    chg_beta_sq <- sum((as.vector(beta_tilde) - as.vector(beta_tilde_prev))^2)
    chg_pi_sq <- sum((pi_tilde - pi_tilde_prev)^2)
    chg_v_sq <- sum((v - v_prev)^2)
    chg_z_pi1_sq <- sum((z_pi1 - z_pi1_prev)^2)
    chg_w_sq <- sum((w - w_prev)^2)
    chg_z_beta_sq <- sum((as.vector(z_beta) - as.vector(z_beta_prev))^2)
    variables_change_norm <- sqrt(chg_beta_sq + chg_pi_sq + chg_v_sq + chg_z_pi1_sq + chg_w_sq + chg_z_beta_sq)

    val_current <- objective_function_r(as.vector(beta_tilde), pi_tilde, y, x, u_quantile) +
      lamb1_scaled * sum(abs(as.vector(beta_tilde))) +
      lamb2_scaled * nuclear_norm_r(pi_tilde)
    obj_change <- abs(val_current - val_prev)

    if (max(c(variables_change_norm, primal_residual_norm, obj_change), na.rm = TRUE) < tol) { # Added na.rm
      #print(paste("Converged at iteration:", iter))
      break
    }
  }

  if (iter == Niter && max(c(variables_change_norm, primal_residual_norm, obj_change), na.rm = TRUE) >= tol) { # Added na.rm
    print("Reached maximum number of iterations without convergence.")
  }

  beta_hat <- z_beta
  Pi_hat <- pi_tilde

  return(beta_hat)
}


#' @param data: Dataframe containing the nodes of the graph.
#' @param tau: The quantile level for quantile regression (e.g., 0.5 for the median).
#' @param lambda1: The L1 penalty parameter for the first Lasso quantile regression.
#' @param lambda2: The second penalty parameter for the QGM method.
#' @param coef_threshold: A threshold below which estimated coefficients are considered zero.
#' @return A list containing the true directed adjacency matrix and the estimated undirected adjacency matrix.
belloni_adj_matrix <- function(data, tau, u_q, lambda1 = 60, lambda2 = 10, tolerance = 1e-2, iterations = 4000, coef_threshold = 1e-5) {
  # This matrix will store estimated coefficients where entry (k, j) is the
  # QGM coefficient of Y_j in the regression of Y_k using Y_{\setminus k}.
  n <- dim(data)[1]
  t <- dim(data)[2]
  p <- dim(data)[3]
  estimated_beta_matrix <- matrix(0, nrow = p, ncol = p)

  for (k in 1:p) {
    y_k <- matrix(data[, , k, drop = FALSE], nrow = n, ncol = t)
    X_k <- data[, , -k, drop = FALSE] # Regress on ALL other variables
    # Fit ADMM method
    coeffs <- admm_l1_nuclear_norm_qr_v7_r(
      y = y_k,
      x = X_k,
      u_quantile = u_q,
      lamb1_orig = lambda1,
      lamb2_orig = lambda2,
      tol = tolerance,
      Niter = iterations,
      initial = 0
    )
    original_indices_of_predictors <- (1:p)[-k]

    estimated_beta_matrix[k, original_indices_of_predictors] <- coeffs
  }

  for (j in 1:p) {
    for (k in 1:p) {
      if (j != k) {
        # Check if coefficient of Yj in regression of Yk is non-zero (above threshold)
        # AND coefficient of Yk in regression of Yj is non-zero (above threshold)
        # estimated_beta_matrix[k, j] is the coefficient of Y_j when Y_k is the response
        # estimated_beta_matrix[j, k] is the coefficient of Y_k when Y_j is the response
        if (abs(estimated_beta_matrix[k, j]) < coef_threshold ||
            abs(estimated_beta_matrix[j, k]) < coef_threshold) {
          estimated_beta_matrix[k, j] <- 0
          estimated_beta_matrix[j, k] <- 0 # Adjacency matrix is symmetric for undirected graph
        }
      }
    }
  }
  return(estimated_beta_matrix)
}


#' @param adj_matrix: The square adjacency matrix.
#' @param node_names: Optional vector of names for the nodes.
#' @param ...: Additional arguments to pass to qgm.igraph.plot.
plot_adj_matrix_qgm <- function(adj_matrix, node_names = NULL, ...) {
  p <- nrow(adj_matrix)

  if (is.null(node_names)) {
    node_names <- paste0("V", 1:p)
  }

  edge_list <- list()
  for (i in 1:p) {
    neighbors <- which(adj_matrix[i, ] != 0)
    if (length(neighbors) > 0) {
      edge_list[[node_names[i]]] <- list(edges = node_names[neighbors])
    } else {
      edge_list[[node_names[i]]] <- list(edges = character(0))
    }
  }

  g_nel <- new("graphNEL", nodes = node_names, edgeL = edge_list, edgemode = "undirected")

  qgm.igraph.plot(g_nel, ...)
}
