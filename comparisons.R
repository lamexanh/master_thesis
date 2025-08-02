cummulative_variations  <- function(pca) {
  result <- summary(pca)$importance["Cumulative Proportion", 2]

  return(result)
}

procrustes_analysis<- function(pca, pca_base) {
  scores_full    <- pca_base$x[, 1:2]
  scores_sparse  <- pca$x[, 1:2]
  proc           <- procrustes(X = scores_full, Y = scores_sparse, scale = TRUE)
  sum_of_squares <- proc$ss 

  return(sum_of_squares)
}

compare <- function(output_1, output_2, output_base) {
  # output_1: FCPCA output
  # output_2: Vanilla PCA output
  # output_base: Oracal PCA output (benchmark)
  # Evaluation 1: compare estimated parameters vs true parameters (mean and covariance)
  mean_estimation_error_FCPCA         <- sqrt(sum((output_1$pca$center - output_base$pca$center)^2)) 
  mean_estimation_error_vanilla       <- sqrt(sum((output_2$pca$center - output_base$pca$center)^2)) 
  covariance_estimation_error_FCPCA   <- sqrt(sum((output_1$K_estimated - output_base$K_estimated)^2))
  covariance_estimation_error_vanilla <- sqrt(sum((output_2$K_estimated - output_base$K_estimated)^2)) 

  # Evaluation 2: compare variation explained by first two PCs (on the same simulated data)
  variance_explained_FCPCA   <- cummulative_variations(output_1$pca) 
  variance_explained_vanilla <- cummulative_variations(output_2$pca) 

  # Evaluation 3: assess the structural similarity between the PCA, using Procrustes sum of squares 
  library(vegan)  # Ensure vegan is loaded for procrustes function
  procrustes_distance_FCPCA   <- procrustes_analysis(output_1$pca, output_base$pca) 
  procrustes_distance_vanilla <- procrustes_analysis(output_2$pca, output_base$pca) 
  
  return(list(
    mean_estimation_error_FCPCA = mean_estimation_error_FCPCA, # critical_value_1
    mean_estimation_error_vanilla = mean_estimation_error_vanilla, # critical_value_2
    covariance_estimation_error_FCPCA = covariance_estimation_error_FCPCA, # critical_value_3
    covariance_estimation_error_vanilla = covariance_estimation_error_vanilla, # critical_value_4
    variance_explained_FCPCA = variance_explained_FCPCA, # critical_value_5
    variance_explained_vanilla = variance_explained_vanilla, # critical_value_6
    procrustes_distance_FCPCA = procrustes_distance_FCPCA, # critical_value_7
    procrustes_distance_vanilla = procrustes_distance_vanilla # critical_value_8
  ))
}

