# ============================================================================
# BIPLOT: PLOT 1 PCA METHOD WITH PARTY CONTRIBUTIONS
# Coordinate systems are from PC1 and PC2 of oracle PCA.
# Other methods was aligned using Procrustes transformation.
# ============================================================================

create_pca_biplot <- function(
  pca_object, method_name = "PCA", 
  oracle_pca              = NULL, 
  use_procrustes          = TRUE, 
  cluster_assignments # result from kmeans clustering
) {
  library(ggplot2)
  library(ggrepel)
  library(vegan)  # For procrustes transformation
  
  # Determine if we need Procrustes transformation
  is_oracle_method <- grepl("oracle", tolower(method_name)) || is.null(oracle_pca)
  
  # Get the PCA scores
  if (is_oracle_method || !use_procrustes) {
    # Use oracle coordinates
    pca_scores_matrix <- pca_object$x[, 1:2]
    coord_note        <- ""
  } else {
    proc_result       <- procrustes(X = oracle_pca$x[, 1:2], Y = pca_object$x[, 1:2], scale = TRUE)
    pca_scores_matrix <- proc_result$Yrot
    coord_note        <- " (State positions aligned to Oracle space)"
  }
  
  # Create data frame with state names and PC1-PC2 coordinates for plotting
  pc_scores <- data.frame(
    State   = rownames(pca_scores_matrix),
    PC1     = pca_scores_matrix[, 1],
    PC2     = pca_scores_matrix[, 2],
    Cluster = as.factor(cluster_assignments)
  )
  
  # Define consistent cluster colors
  cluster_colors <- c(
    "1" = "#ff7700", 
    "2" = "#f303f7", 
    "3" = "#0dcd40da", 
    "4" = "#e3b11a", 
    "5" = "#1F78B4", 
    "6" = "#06e0fd"
  )
  # Adjust this if you have more than 6 clusters
  
  # Extract loadings (party contributions) - first two components in native spaces
  loading_scale  <- 3  # Adjust this to make arrows more visible
  party_loadings <- data.frame(
    Party       = rownames(pca_object$rotation),
    PC1_loading = pca_object$rotation[, 1] * loading_scale,
    PC2_loading = pca_object$rotation[, 2] * loading_scale
  )
  
  # Calculate variance explained by PC1 and PC2
  variance_explained <- summary(pca_object)$importance["Proportion of Variance", 1:2]
  pc1_var            <- round(variance_explained[1] * 100, 1)
  pc2_var            <- round(variance_explained[2] * 100, 1)
  
  # Create biplot
  biplot <- ggplot() +
    # Plot state scores as points colored by cluster
    geom_point(
      data = pc_scores, 
      aes(x = PC1, y = PC2, color = Cluster), 
      size = 3, 
      alpha = 0.7
    ) +
    
    # Add state labels colored by cluster
    geom_text_repel(
      data = pc_scores, 
      aes(x = PC1, y = PC2, label = State, color = Cluster),
      size = 3, 
      alpha = 0.8,
      box.padding = 0.3, 
      point.padding = 0.3,
      max.overlaps = 20, 
      seed = 123
    ) +
    
    # Add cluster color scale
    scale_color_manual(values = cluster_colors, name = "Cluster") +
    
    # Add party loading arrows
    geom_segment(
      data = party_loadings, 
      aes(x = 0, y = 0, xend = PC1_loading, yend = PC2_loading),
      arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
      color = "red", 
      linewidth = 1, 
      alpha = 0.8
    ) +
    
    # Add party labels at arrow tips
    geom_text_repel(
      data = party_loadings, 
      aes(x = PC1_loading, y = PC2_loading, label = Party),
      size = 3.5, 
      color = "red", 
      fontface = "bold",
      box.padding = 0.5, 
      point.padding = 0.3,
      max.overlaps = 20, 
      seed = 456
    ) +
    
    # Add axes through origin
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    
    # Labels and theme
    labs(
      title    = paste(method_name, " - Biplot for German federal election data in 2025"),
      subtitle = paste("Points = States (colored by true cluster), Red arrows = Parties", coord_note),
      x        = paste0("PC1 (", pc1_var, "% of variance)"),
      y        = paste0("PC2 (", pc2_var, "% of variance)")
    ) +

    theme_classic() +

    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "right",
      plot.caption = element_text(size = 10, face = "italic"),
      panel.grid.minor = element_blank()
    )
  
    return(list(
      plot = biplot,
      pca_object = pca_object
    )
  )
}

# ============================================================================
# PLOT SCORES IN TWO FIRST PRINCIPAL COMPONENTS FOR ORACLE, VANILLA, Functional Compositional PCA 
# (given 1 impute-method and 1 sample)
# ============================================================================

create_plot_with_state_labels <- function(results_file, sample_number, real_clr_proportions) {
  library(ggrepel)
  library(ggplot2)
  library(vegan)  # For procrustes function

  # Load saved results
  results       <- readRDS(results_file)
  impute_method <- gsub("all_results_", "", gsub(".rds", "", basename(results_file)))
  
  # Get the PCA results
  oracle_pca                   <- pca_and_covariance_on_clr(real_clr_proportions)$pca
  vanilla_pca                  <- results[[sample_number]]$vanilla_output$pca
  functional_compositional_pca <- results[[sample_number]]$functional_compositional_pca_output$pca
  
  # Apply Procrustes transformation
  vanilla_proc                      <- procrustes(X = oracle_pca$x[, 1:2], Y = vanilla_pca$x[, 1:2], scale = TRUE)
  functional_compositional_pca_proc <- procrustes(X = oracle_pca$x[, 1:2], Y = functional_compositional_pca$x[, 1:2], scale = TRUE)
  
  # Create coordinates with state labels
  num_states <- nrow(oracle_pca$x)
  all_coords <- data.frame(
    PC1 = c(oracle_pca$x[, 1], vanilla_proc$Yrot[, 1], functional_compositional_pca_proc$Yrot[, 1]),
    PC2 = c(oracle_pca$x[, 2], vanilla_proc$Yrot[, 2], functional_compositional_pca_proc$Yrot[, 2]),
    Method = rep(c("Oracle (Target)", "Vanilla PCA", "Functional Compositional PCA"), each = num_states),
    State = rep(rownames(oracle_pca$x), 3)
  )
  
  # Create plot with state labels
  comparison_plot <- ggplot(all_coords, aes(x = PC1, y = PC2, color = Method, shape = Method)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text_repel(
      aes(label = State), 
      size = 3, 
      alpha = 0.9,
      box.padding = 0.3,
      point.padding = 0.3,
      max.overlaps = 20,
      seed = 123
    ) +
    scale_color_manual(values = c("Oracle (Target)" = "#E31A1C", "Vanilla PCA" = "#0509ea", "Functional Compositional PCA" = "#1add0c")) +
    scale_shape_manual(values = c("Oracle (Target)" = 16, "Vanilla PCA" = 17, "Functional Compositional PCA" = 15)) +
    labs(
      title    = paste("State Positions in Oracle PC1-PC2 Space (Procrustes-aligned)"),
      subtitle = paste("Imputation Method:", impute_method, "| Sample", sample_number, "| Red = Oracle (Target), Blue = Vanilla PCA, Green = Functional Compositional PCA"),
      x        = "First Principal Component (PC1)",
      y        = "Second Principal Component (PC2)",
      color    = "Method",
      shape    = "Method"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  return(comparison_plot)
}
