# Perform PCA and estimate covariance on CLR-transformed compositional data
pca_and_covariance_on_clr <- function(clr_proportions) {
  pca_result  <- prcomp(na.omit(clr_proportions), center = TRUE, scale. = FALSE)
  K_estimated <- cov(clr_proportions)

  return(list(pca = pca_result, K_estimated = K_estimated))
}

# Gradient of Log-Posterior for Compositional PCA Scores
gradient_log_posterior_compositional <- function(scores, x_data_i, pca) {
  q   <- length(scores)
  Phi <- pca$rotation[, 1:q, drop = FALSE]      # N × q
  rho <- pca$center + Phi %*% scores            # N × 1

  exp_rho <- exp(rho)
  pi      <- exp_rho / sum(exp_rho)              # N × 1
  s_i     <- sum(x_data_i)

  # Gradient of log-likelihood: Φᵀ (x − s π)
  term    <- x_data_i - s_i * pi                 # N × 1
  grad_ll <- as.numeric(t(Phi) %*% term)         # q × 1

  # Gradient of log-prior: −Λ z  (Gaussian prior)
  sigma2     <- pca$sdev[1:q]^2
  grad_prior <- -scores / sigma2                 # q × 1

  return(grad_ll + grad_prior)                   # q × 1
}

# Log-Posterior Density for Compositional PCA Scores
log_posterior_density_compositional <- function(scores, x_data_i, pca) {
  # clr-transformed density based on the PCA rotation and scores
  rho_i <- pca$center + pca$rotation %*% scores

  # clr-inverse transformation for compositional data (use exp and row-wise sum)
  pi_i <- exp(rho_i) / sum(exp(rho_i))

  # Multinomial log-likelihood for observed data
  log_likelihood <- sum(x_data_i * log(pi_i))

  # Gaussian prior penalty on scores
  penalty <- - sum(0.5 * scores^2 / (pca$sdev^2))

  # Combine log-likelihood and penalty
  return(log_likelihood + penalty)
}

# FUNCTIONAL COMPOSITIONAL PCA METHOD - METHOD 3
#
# This function iteratively refines a PCA decomposition and covariance estimate
# for compositional data using simulated samples and importance weighting.
# It updates the PCA center, rotation, and singular values based on weighted proposals.
functional_compositional_pca <- function(pca, simulated_df, dim_reduction, r, max_iter, lambda, eps) {
  which_reduced <- rev(cumsum(rev(pca$sdev^2))/sum(pca$sdev^2) > dim_reduction)
  which_reduced <- which_reduced | c(TRUE, TRUE, rep(FALSE, length(which_reduced) - 2))

  pca$sdev     <- pca$sdev[which_reduced]
  pca$rotation <- pca$rotation[, which_reduced, drop = FALSE]

  proposal_scores <- vector("list", length(simulated_df))
  weights         <- vector("list", length(simulated_df))

  for (k in 1:max_iter) {
    # E-Step
    for (i in 1:nrow(simulated_df)) {
      optim_result <- optim(
        rep(0, length = length(pca$sdev)), 
        log_posterior_density_compositional, 
        gr       = gradient_log_posterior_compositional,
        x_data_i = unlist(round(simulated_df[i, ])), 
        pca      = pca,
        control  = list(fnscale = -1), 
        method   = "BFGS"
      )

      scores_median <- as.vector(optim_result$par)

      proposal_scores[[i]] <- sapply(1:(r * k), function(t) {
        matrix(rnorm(length(scores_median), mean = scores_median, sd = lambda * pca$sdev))
      })

      log_weights <- apply(proposal_scores[[i]], 2, function(scores) {
        log_posterior_density_compositional(scores, x_data_i = unlist(round(simulated_df[i, ])), pca = pca) -
          sum(dnorm(scores, mean = scores_median, sd = lambda * pca$sdev, log = TRUE))
      })

      log_weights  <- log_weights - mean(log_weights, na.rm = TRUE)
      weights[[i]] <- exp(log_weights) / sum(exp(log_weights))
    }

    # M-Step
    mu_scores <- rowMeans(sapply(seq_along(weights), function(i) {
      proposal_scores[[i]] %*% weights[[i]]
    }))

    # update pca and other parameters
    pca_old    <- pca
    pca$center <- pca$center + pca$rotation %*% mu_scores
    pca$center <- pca$center - mean(pca$center)

    Sigma <- Reduce("+", lapply(seq_along(weights), function(i) {
      Reduce("+", lapply(1:(r * k), function(t) {
        weights[[i]][t] * (proposal_scores[[i]][, t] - mu_scores) %*%
          t((proposal_scores[[i]][, t] - mu_scores))
      }))
    })) / length(weights)

    eigen_decomp <- eigen(Sigma)
    pca$sdev     <- sqrt(eigen_decomp$values)
    pca$rotation <- pca$rotation %*% eigen_decomp$vectors
    pca$rotation <- sweep(pca$rotation, 2, colMeans(pca$rotation), "-")
    pca$x        <- scale(clr_proportions, center = pca$center, scale = FALSE) %*% pca$rotation

    K_old <- Reduce("+", lapply(seq_along(pca_old$sdev), function(k) {
      pca_old$rotation[, k] %*% t(pca_old$rotation[, k]) * (pca_old$sdev[k]^2)
    }))

    K_new <- Reduce("+", lapply(seq_along(pca$sdev), function(k) {
      pca$rotation[, k] %*% t(pca$rotation[, k]) * (pca$sdev[k]^2)
    }))

    # Check convergence using single tolerance for both center and covariance changes
    center_change <- sqrt(sum((pca_old$center - pca$center)^2))
    cov_change    <- sqrt(sum((K_old - K_new)^2))
    converged     <- (center_change < eps) && (cov_change < eps)
    
    # Early stopping if converged
    if (converged) {
      break
    }
  }

  return(list(
    pca = pca,
    K_estimated = K_new
  ))
}

