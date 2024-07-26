clara_fun <- function(dataset, n, nsamples, samplesize, ncores, ncluster) {
  
  cost_medoids_mat_sum <- NULL
  medoids <- list()
  
  for (i in 1:nsamples) {
    
    medoid_indices <- "Cluster_Medoids_error"
    while (is.character(medoid_indices)) {
      binarydata_sample <- dataset[sample(1:n, samplesize, replace = FALSE), ]
      medoid_indices <- tryCatch(Cluster_Medoids(binarydata_sample, ncluster, distance_metric = "hamming", threads = ncores, seed = i)$medoid_indices,
                                 error = function(x) {
                                   message("trycatch hat gegriffen clara")
                                   medoid_indices <- "Cluster_Medoids_error"
                                   return(medoid_indices)
                                 })
    }
    
    medoids[[i]] <- binarydata_sample[medoid_indices,]
    cost_medoids_mat <- rdist:::hamming_cdist(dataset, medoids[[i]])
    cost_medoids_mat_sum[i] <- sum(rowMins(cost_medoids_mat))
    
    if (i == 1) {
      cluster_clara <- max.col(-cost_medoids_mat, ties="first")
    } else {
      if (cost_medoids_mat_sum[i] < cost_medoids_mat_sum[i - 1]) {
        cluster_clara <- max.col(-cost_medoids_mat, ties="first")
      }
    }
    rm(cost_medoids_mat)
  }
  
  return(cluster_clara)
}