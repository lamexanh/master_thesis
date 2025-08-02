# Simulate compositional count data based on real compositional probabilities
simulate_compositional_counts <- function(real_compositions) {
  num_states            <- nrow(real_compositions)
  num_samples_per_state <- 1 
  log_lambda            <- runif(num_states * num_samples_per_state, min = 1, max = 3)
  lambda                <- 10^log_lambda
  total_counts          <- sapply(lambda, function(l) rpois(1, l))
  total_counts          <- matrix(total_counts, nrow = num_samples_per_state, ncol = num_states)
  num_categories        <- ncol(real_compositions)
  state_names           <- rownames(real_compositions)
  party_names           <- colnames(real_compositions)  
  
  simulated_counts <- do.call(rbind, lapply(1:num_states, function(i) {
    state_prob <- as.numeric(real_compositions[i, ])
    state_data <- t(sapply(total_counts[, i], function(n_ij) rmultinom(1, size = n_ij, prob = state_prob)))
    state_data
  }))
  
  simulated_data           <- as.data.frame(simulated_counts)
  colnames(simulated_data) <- party_names
  rownames(simulated_data) <- state_names
  
  return(simulated_data)
}

# Impute zeros in compositional count data using a specified method
impute_compositional_zeros <- function(count_data, impute_method) {
  if (impute_method %in% c("CZM", "GBM", "pseudo_counts_0.5")) {
    # impute_method is valid; nothing to change
  } else {
    stop("Invalid impute_method. Use 'CZM', 'GBM', or 'pseudo_counts_0.5'.")
  }
  library(zCompositions)

  if (all(count_data != 0)) {
    return(count_data)
  }

  if (impute_method == "pseudo_counts_0.5") {
    count_data <- count_data + 0.5
  }
  else if (impute_method %in% c("CZM", "GBM")) {
    count_data <- cmultRepl(count_data, label = 0, method = impute_method, output = "p-counts", z.warning = 0.8, z.delete = FALSE)
  }

  return(count_data)
}

# Convert imputed compositional count data to compositional probabilities.
calculate_compositional_proportions <- function(imputed_counts) {
  row_sums    <- rowSums(imputed_counts)
  proportions <- imputed_counts / row_sums

  return(proportions)
}

# Apply centered log-ratio (CLR) transformation to compositional proportions.
clr_transform <- function(proportions) {
  log_proportions      <- log(proportions)
  mean_log_proportions <- rowMeans(log_proportions)
  clr_values           <- log_proportions - mean_log_proportions

  return(clr_values)
}
