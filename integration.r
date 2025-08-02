required_packages <- c("zCompositions", "ggplot2", "vegan", "reshape2", "gridExtra", "ggrepel")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

source("transformations.R")
source("comparisons.R")
source("methods.R") 
source("visualizations.R")

REAL_DATA             <- read.csv(file = "election_2025.csv", header = TRUE)[,2:8]
STATE_NAMES           <- read.csv(file = "election_2025.csv", header = TRUE)[,1]
row.names(REAL_DATA)  <- STATE_NAMES
REAL_DATA_PROPORTIONS <- calculate_compositional_proportions(REAL_DATA)
REAL_CLR_PROPORTIONS  <- clr_transform(REAL_DATA_PROPORTIONS)

NUMBER_OF_SAMPLES = 100
IMPUTE_METHODS    = c("CZM", "pseudo_counts_0.5", "GBM")  # List of imputation methods to compare

DIM_REDUCTION = 0.01
R             = 10
MAX_ITERATION = 200
LAMBDA        = 1
EPS           = 0.01

# ============================================================================
# STEP 1: GENERATE SAMPLES ONCE (same for all imputation methods)
# ============================================================================

set.seed(123)

samples_data <- list()

for (i in 1:NUMBER_OF_SAMPLES) {
  cat("Generating sample", i, "of", NUMBER_OF_SAMPLES, "\n")
  simulated_data    <- simulate_compositional_counts(REAL_DATA_PROPORTIONS)
  samples_data[[i]] <- simulated_data
}

# Save the base samples for reproducibility
saveRDS(samples_data, "base_samples_for_comparing_imputation.rds")

# ============================================================================
# STEP 2: APPLY EACH IMPUTATION METHOD TO THE SAME SAMPLES
# ============================================================================

# Store results for all methods
all_methods_results <- list()

for (method in IMPUTE_METHODS) {
  cat("Processing imputation method:", method, "\n")  # Add this line

  # Results for this specific method
  method_results <- list()

  # Apply this method to all samples
  for (i in 1:NUMBER_OF_SAMPLES) {
    cat("  Sample", i, "with method", method, "\n")  # Add this line

    # Use the same base sample
    simulated_data <- samples_data[[i]]

    # Apply the current imputation method
    imputed_data <- impute_compositional_zeros(simulated_data, method)

    # Transform to compositional and CLR
    compositional_proportions <- calculate_compositional_proportions(imputed_data)
    clr_proportions           <- clr_transform(compositional_proportions)

    # Get oracle output (same for all methods)
    oracle_output <- pca_and_covariance_on_clr(REAL_CLR_PROPORTIONS)

    # Get stardard PCA output
    vanilla_output <- pca_and_covariance_on_clr(clr_proportions)

    # Get Functional Compositional PCA output
    functional_compositional_pca_output <- functional_compositional_pca(
      vanilla_output$pca, 
      simulated_data,
      DIM_REDUCTION,
      R,
      MAX_ITERATION,
      LAMBDA,
      EPS
    )

    # Compare methods
    comparison <- compare(functional_compositional_pca_output, vanilla_output, oracle_output)

    # Store results for this sample and method
    method_results <- append(method_results, list(
      list(
        sample_number                       = i,
        impute_method                       = method,
        simulated_data                      = simulated_data,
        vanilla_output                      = vanilla_output,
        functional_compositional_pca_output = functional_compositional_pca_output,
        comparison                          = comparison
      )
    ))
  }

  # Save results for this method
  method_filename <- paste0("all_results_", method, ".rds")
  saveRDS(method_results, method_filename)
  cat("Results for", method, "saved as:", method_filename, "\n")  

  # Store in overall results
  all_methods_results[[method]] <- method_results

  # Save comparison CSV for this method
  csv_filename <- paste0("comparison_results_", method, ".csv")
  write.csv(
    do.call(rbind, lapply(method_results, function(x) {
      data.frame(
        sample_number                       = x$sample_number,
        impute_method                       = x$impute_method,
        mean_estimation_error_FCPCA         = x$comparison$mean_estimation_error_FCPCA, # critical_value_1
        mean_estimation_error_vanilla       = x$comparison$mean_estimation_error_vanilla, # critical_value_2
        covariance_estimation_error_FCPCA   = x$comparison$covariance_estimation_error_FCPCA, # critical_value_3
        covariance_estimation_error_vanilla = x$comparison$covariance_estimation_error_vanilla, # critical_value_4
        variance_explained_FCPCA            = x$comparison$variance_explained_FCPCA, # critical_value_5
        variance_explained_vanilla          = x$comparison$variance_explained_vanilla, # critical_value_6
        procrustes_distance_FCPCA           = x$comparison$procrustes_distance_FCPCA, # critical_value_7
        procrustes_distance_vanilla         = x$comparison$procrustes_distance_vanilla # critical_value_8
      )
    })),
    file = csv_filename,
    row.names = FALSE
  )
}

# ============================================================================
# STEP 3: CREATE COMBINED COMPARISON RESULTS
# ============================================================================

# Combine all methods into one CSV for easy comparison
combined_results <- do.call(rbind, lapply(IMPUTE_METHODS, function(method) {
  method_data <- all_methods_results[[method]]
  do.call(rbind, lapply(method_data, function(x) {
    data.frame(
      sample_number                       = x$sample_number,
      impute_method                       = x$impute_method,
      mean_estimation_error_FCPCA         = x$comparison$mean_estimation_error_FCPCA, # critical_value_1
      mean_estimation_error_vanilla       = x$comparison$mean_estimation_error_vanilla, # critical_value_2
      covariance_estimation_error_FCPCA   = x$comparison$covariance_estimation_error_FCPCA, # critical_value_3
      covariance_estimation_error_vanilla = x$comparison$covariance_estimation_error_vanilla, # critical_value_4
      variance_explained_FCPCA            = x$comparison$variance_explained_FCPCA, # critical_value_5
      variance_explained_vanilla          = x$comparison$variance_explained_vanilla, # critical_value_6
      procrustes_distance_FCPCA           = x$comparison$procrustes_distance_FCPCA, # critical_value_7
      procrustes_distance_vanilla         = x$comparison$procrustes_distance_vanilla # critical_value_8
    )
  }))
}))

write.csv(combined_results, "combined_comparison_all_methods.csv", row.names = FALSE)

# Save all methods results
saveRDS(all_methods_results, "all_methods_results.rds")

#============================================================================
# COMPARISON VISUALIZATION EVALUATION METRICS 
# TO COMPARE Functional Compositional PCA METHOD VS VANILLA METHOD, USING ORACLE METHOD AS BENCHMARK
# GROUP BY DIFFERENT IMPUTE-METHOD
# ============================================================================

create_performance_comparison_boxplots <- function() {
  library(ggplot2)
  # Read the combined results
  combined_results <- read.csv("combined_comparison_all_methods.csv")

  # 1. Critical Values 1 & 2: Mean estimation error 
  mean_distance_data <- data.frame(
    Comparison_Metric = rep(c("Functional Compositional PCA - Oracle", "Vanilla Method - Oracle"), each = nrow(combined_results)),
    Impute_Method     = rep(combined_results$impute_method, 2),
    Distance          = c(combined_results$mean_estimation_error_FCPCA, combined_results$mean_estimation_error_vanilla),
    Sample            = rep(combined_results$sample_number, 2)
  )

  p1 <- ggplot(mean_distance_data, aes(x = Impute_Method, y = Distance, fill = Comparison_Metric)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
    labs(
      title    = "Mean Estimation Error",
      subtitle = "Distance between estimated means vs true mean (Lower = Better)",
      x        = "Imputation Method", 
      y        = "Distance",
      fill     = "Comparison Metric"
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom", 
      plot.title      = element_text(size = 14, face = "bold"),
      plot.subtitle   = element_text(size = 12)
    )

  # 2. Critical Values 3 & 4: Covariance estimation error
  k_distance_data <- data.frame(
    Comparison_Metric = rep(c("Functional Compositional PCA - Oracle", "Vanilla - Oracle"), each = nrow(combined_results)),
    Impute_Method     = rep(combined_results$impute_method, 2),
    Distance          = c(combined_results$covariance_estimation_error_FCPCA, combined_results$covariance_estimation_error_vanilla),
    Sample            = rep(combined_results$sample_number, 2)
  )

  p2 <- ggplot(k_distance_data, aes(x = Impute_Method, y = Distance, fill = Comparison_Metric)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
    labs(
      title    = "Covariance Matrix Estimation Error",
      subtitle = "Distance between estimated covariances vs true covariance (Lower = Better)",
      x        = "Imputation Method", 
      y        = "Distance",
      fill     = "Comparison Metric"
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom", 
      plot.title      = element_text(size = 14, face = "bold"),
      plot.subtitle   = element_text(size = 12)
    )

  # 3. Critical Values 5 & 6: Compare variation explained by the first two PCs
  variation_data <- data.frame(
    Comparison_Metric           = rep(c("Functional Compositional PCA", "Vanilla PCA"), each = nrow(combined_results)),
    Impute_Method               = rep(combined_results$impute_method, 2),
    PC1_PC2_Cumulative_Variance = c(combined_results$variance_explained_FCPCA, combined_results$variance_explained_vanilla),
    Sample                      = rep(combined_results$sample_number, 2)
  )
  
  p3 <- ggplot(variation_data, aes(x = Impute_Method, y = PC1_PC2_Cumulative_Variance, fill = Comparison_Metric)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(title = "Explained Variance (Higher = Better)",
         subtitle = "Cumulative variation explained by first 2 PCs",
         x = "Imputation Method", 
         y = "PC1 & PC2 Cumulative Variance",
         fill = "Comparison Metric") +
    theme_classic() +
    theme(legend.position = "bottom", 
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12))
  
  # 4. Critical Values 7 & 8: Procrustes Distance Comparison
  procrustes_data <- data.frame(
    Comparison_Metric = rep(c("Functional Compositional PCA vs Oracle", "Vanilla vs Oracle"), each = nrow(combined_results)),
    Impute_Method = rep(combined_results$impute_method, 2),
    Procrustes_Distance = c(combined_results$procrustes_distance_FCPCA, combined_results$procrustes_distance_vanilla),
    Sample = rep(combined_results$sample_number, 2)
  )
  
  p4 <- ggplot(procrustes_data, aes(x = Impute_Method, y = Procrustes_Distance, fill = Comparison_Metric)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
    labs(
      title    = "Structural Similarity Between PCA Methods",
      subtitle = "Sum of squared distances between scores from Oracle and from other PCA being aligned to Oracle space (Lower = Better)",
      x        = "Imputation Method", 
      y        = "Sum of Squares of Procrustes Distances",
      fill     = "Comparison Metric"
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom", 
      plot.title      = element_text(size = 14, face = "bold"),
      plot.subtitle   = element_text(size = 12)
    )

  # Combine all plots
  library(gridExtra)
  combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
  
  # Save the combined plot
  ggsave("performance_comparison_boxplots_all_methods.png", combined_plot, width = 16, height = 12, dpi = 300)
  cat("Performance comparison boxplots saved as: performance_comparison_boxplots_all_methods.png\n")
  
  return(list(mean_distance = p1, k_distance = p2, variation = p3, procrustes = p4))
}

# If you assign the function output to a variable:
plots <- create_performance_comparison_boxplots()
saveRDS(plots, "performance_comparison_plots.rds")

# You can access individual plots:
plots$mean_distance    # Critical values 1 & 2 plot
plots$k_distance       # Critical values 3 & 4 plot  
plots$variation        # Critical values 5 & 6 plot
plots$procrustes       # Critical values 7 & 8 plot

# Save individual plots
ggsave("mean_distance_comparison.png", plots$mean_distance, width = 12, height = 8, dpi = 300)
ggsave("covariance_distance_comparison.png", plots$k_distance, width = 12, height = 8, dpi = 300)
ggsave("variance_explained_comparison.png", plots$variation, width = 12, height = 8, dpi = 300)
ggsave("procrustes_distance_comparison.png", plots$procrustes, width = 12, height = 8, dpi = 300)

# ============================================================================
# PLOT SCORES IN TWO FIRST PRINCIPAL COMPONENTS FOR Oracle, VANILLA, Functional Compositional PCA 
# (given 1 impute-method and 1 sample)
# ============================================================================

state_labels_plot <- create_plot_with_state_labels('all_results_pseudo_counts_0.5.rds', 1, REAL_CLR_PROPORTIONS)
ggsave("state_labels_comparison_sample1_pseudo_counts_0.5.png", state_labels_plot, width = 12, height = 8, dpi = 300)

# ============================================================================
# BIPLOT: PLOT ANY PCA METHOD WITH PARTY CONTRIBUTIONS 
# ============================================================================

oracle_pca          <- pca_and_covariance_on_clr(REAL_CLR_PROPORTIONS)$pca
oracle_scores       <- oracle_pca$x[, 1:2]
cluster_assignments <- kmeans(as.data.frame(oracle_scores), centers = 3, nstart = 20)$cluster

oracle_biplot <- create_pca_biplot(oracle_pca, "Oracle PCA", cluster_assignments = cluster_assignments)
cat("\n=== PLOTTING ORACLE PCA ===\n")
ggsave("oracle_biplot.png", oracle_biplot$plot, width = 10, height = 8, dpi = 300)


# Example 2: Plot Vanilla Method with Procrustes alignment (for a specific sample, make sure input suitable with NUMBER_OF_SAMPLES above)
cat("\n=== PLOTTING VANILLA METHOD (Sample 1, pseudo_counts_0.5) ===\n")
vanilla_pca    <- all_methods_results[["pseudo_counts_0.5"]][[1]]$vanilla_output$pca  # Sample 1, pseudo_counts_0.5 imputation
vanilla_biplot <- create_pca_biplot(vanilla_pca, "Vanilla Method", oracle_pca = oracle_pca, use_procrustes = TRUE, cluster_assignments = cluster_assignments)
ggsave("vanilla_biplot_sample1_pseudo_counts_0.5.png", vanilla_biplot$plot, width = 10, height = 8, dpi = 300)

# Example 3: Plot Functional Compositional PCA method with Procrustes alignment (for a specific sample, make sure input suitable with NUMBER_OF_SAMPLES above)
cat("\n=== PLOTTING Functional Compositional PCA method (Sample 1, pseudo_counts_0.5) ===\n")
functional_compositional_pca        <- all_methods_results[["pseudo_counts_0.5"]][[1]]$functional_compositional_pca_output$pca  # Sample 1, pseudo_counts_0.5 imputation
functional_compositional_pca_biplot <- create_pca_biplot(functional_compositional_pca, "Functional Compositional PCA Method", oracle_pca = oracle_pca, use_procrustes = TRUE, cluster_assignments = cluster_assignments)

ggsave("functional_compositional_PCA_biplot_sample1_pseudo_counts_0.5.png", functional_compositional_pca_biplot$plot, width = 10, height = 8, dpi = 300)
