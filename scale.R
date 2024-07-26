delta_add_scale <- function(cf_i, x, cf, N_all) {
  occurence_new <- cf[[cf_i]]$occurence + x
  S_new <- cf[[cf_i]]$S + sum(x)
  N_new <- cf[[cf_i]]$N + 1
  WCD_new <- sum(occurence_new^2) / (S_new * N_new)
  EWCD <- sum(sapply(cf[-cf_i], function(y, N_all) {y$N / N_all * y$WCD}, N_all)) + ((N_new / N_all) * WCD_new)
  
  return(list(occurence = occurence_new, S = S_new, N = N_new, WCD = WCD_new, EWCD = EWCD))
}


scale_cluster <- function(dataset, nsample, ncluster, ncores, seed, maxiteration, maxmoved, chunk_n) {
  
  n_data <- nrow(dataset)
  
  # set.seed(seed)
  
  # draw sample for initial clustering via pam
  dataset_sample_index <- sample(1:n_data, nsample)
  dataset_sample <- dataset[dataset_sample_index,]
  
  dataset_sample_dist <- rdist(dataset_sample, metric = "hamming")
  sample_cluster <- pam(dataset_sample, k = ncluster, cluster.only = TRUE)
  # sample_cluster <- Cluster_Medoids(dataset_sample, ncluster, distance_metric = "jaccard_coefficient", threads = ncores)$clusters
  
  # initial clusters
  clusters <- NULL
  # clusters[dataset_sample_index] <- sample_cluster
  clusterfeatures <- list()
  for (i in 1:ncluster) {
    clusterfeatures[[i]] <- list(occurence = colSums(dataset_sample[sample_cluster == i,, drop = FALSE]),
                                 S = sum(dataset_sample[sample_cluster == i,, drop = FALSE]),
                                 N = sum(sample_cluster == i))
    
    clusterfeatures[[i]]$WCD <- sum(clusterfeatures[[i]]$occurence^2) / (as.numeric(clusterfeatures[[i]]$S) * as.numeric(clusterfeatures[[i]]$N))
  }
  
  
  
  # generate random sequence for processing observations
  rnd_orig <- sample(1:nrow(dataset), nrow(dataset))
  # minus the observations in sample_cluster (they are already in a cluster)
  rnd <- setdiff(rnd_orig, dataset_sample_index)
  
  # dataset <- dataset[-dataset_sample_index,]
  
  # create chunks for processing in parallel. 
  rndseq <- seq(1, length(rnd), ncores * chunk_n)
  
  # phase 1
  # anfang <- Sys.time()
  for (i in 1:length(rndseq)) {
    
    if (i != length(rndseq)) {
      dataiter <- isplitRows(dataset[rnd[rndseq[i]:(rndseq[i + 1]- 1)],], chunks = ncores)
    } else {
      dataiter <- isplitRows(dataset[rnd[rndseq[i]:length(rnd)],], chunks = ncores)
    }
    
    # chunk_n observations get processed at a time
    cluster_result <- foreach(bdata = dataiter, .combine = c, .export = c("clusterfeatures")) %dopar% {
      clusterwhich <- NULL
      for (j in 1:nrow(bdata)) {
        EWCD_result <- lapply(1:length(clusterfeatures), delta_add_scale, bdata[j,], clusterfeatures, n_data)
        EWCD <- sapply(EWCD_result, function(x) x$EWCD)
        # to which cluster will observation be added
        clusterwhich[j] <- which.max(EWCD)
        
        # update clusterfeatures
        clusterfeatures[[clusterwhich[j]]] <- list(occurence = EWCD_result[[clusterwhich[j]]]$occurence,
                                                   S = EWCD_result[[clusterwhich[j]]]$S,
                                                   N = EWCD_result[[clusterwhich[j]]]$N, 
                                                   WCD = EWCD_result[[clusterwhich[j]]]$WCD)
      }
      clusterwhich
      
    }
    
    
    if (i != length(rndseq)) {
      valid_ids <- rnd[rndseq[i]:(rndseq[i + 1] - 1)]
    } else {
      valid_ids <- rnd[rndseq[i]:(length(rnd))]
    }
    
    
    clusters[valid_ids] <- cluster_result
    
    # get results for resulting clusterfeatures for clusters, where the observation has been added
    # cf_results <- mapply(function(x, clusterwhich) x[[clusterwhich]], lapply(v, function(x) x$EWCD_result), resulting_cluster_ids)
    
    # update clusterfeatures with processed cases
    for (j in unique(cluster_result)) {
      clusterfeatures[[j]] <- list(occurence = clusterfeatures[[j]]$occurence + colSums(dataset[valid_ids,][which(cluster_result == j),]),
                                   S = clusterfeatures[[j]]$S + sum(dataset[valid_ids,][which(cluster_result == j),]),
                                   N = clusterfeatures[[j]]$N + nrow(dataset[valid_ids,][which(cluster_result == j),]))
      clusterfeatures[[j]]$WCD <- sum(clusterfeatures[[j]]$occurence^2) / (as.numeric(clusterfeatures[[j]]$S) * as.numeric(clusterfeatures[[j]]$N))
    }
    
  }
  # ende <- Sys.time()
  
  
  clusters[dataset_sample_index] <- sample_cluster
  
  
  # phase 2
  
  maxiter <- 0
  moved <- Inf
  while (maxiter <= maxiteration & moved > maxmoved) {
    maxiter <- maxiter + 1
    moved <- 0
    
    
    for (i in 1:length(rndseq)) {
      
      if (i != length(rndseq)) {
        dataiter <- isplitRows(dataset[rnd[rndseq[i]:(rndseq[i + 1]- 1)],], chunks = ncores)
      } else {
        dataiter <- isplitRows(dataset[rnd[rndseq[i]:length(rnd)],], chunks = ncores)
      }
      
      
      if (i != length(rndseq)) {
        cluster_temp <- clusters[rnd[rndseq[i]:(rndseq[i + 1]- 1)]]
      } else {
        cluster_temp <- clusters[rnd[rndseq[i]:length(rnd)]]
      }
      
      
      
      # subtract observation from cluster
      for (j in unique(cluster_temp)) {
        if (i != length(rndseq)) {
          clusterfeatures[[j]] <- list(occurence = clusterfeatures[[j]]$occurence - colSums(dataset[rnd[rndseq[i]:(rndseq[i + 1]- 1)][which(clusters[rnd[rndseq[i]:(rndseq[i + 1]- 1)]] == j)],]),
                                       S = clusterfeatures[[j]]$S - sum(dataset[rnd[rndseq[i]:(rndseq[i + 1]- 1)][which(clusters[rnd[rndseq[i]:(rndseq[i + 1]- 1)]] == j)],]),
                                       N = clusterfeatures[[j]]$N - length(which(clusters[rnd[rndseq[i]:(rndseq[i + 1]- 1)]] == j)))
          clusterfeatures[[j]]$WCD <- sum(clusterfeatures[[j]]$occurence^2) / (as.numeric(clusterfeatures[[j]]$S) * as.numeric(clusterfeatures[[j]]$N))
        } else {
          clusterfeatures[[j]] <- list(occurence = clusterfeatures[[j]]$occurence - colSums(dataset[rnd[rndseq[i]:length(rnd)][which(clusters[rnd[rndseq[i]:length(rnd)]] == j)],]),
                                       S = clusterfeatures[[j]]$S - sum(dataset[rnd[rndseq[i]:length(rnd)][which(clusters[rnd[rndseq[i]:length(rnd)]] == j)],]),
                                       N = clusterfeatures[[j]]$N - length(which(clusters[rnd[rndseq[i]:length(rnd)]] == j)))
          clusterfeatures[[j]]$WCD <- sum(clusterfeatures[[j]]$occurence^2) / (as.numeric(clusterfeatures[[j]]$S) * as.numeric(clusterfeatures[[j]]$N))
        }
      }
      
      # chunk_n observations get processed at a time
      cluster_result <- foreach(bdata = dataiter, .combine = c, .export = c("clusterfeatures")) %dopar% {
        clusterwhich <- NULL
        for (j in 1:nrow(bdata)) {
          
          
          EWCD_result <- lapply(1:length(clusterfeatures), delta_add_scale, bdata[j,], clusterfeatures, n_data)
          EWCD <- sapply(EWCD_result, function(x) x$EWCD)
          # to which cluster will observation be added
          clusterwhich[j] <- which.max(EWCD)
          
          # update clusterfeatures
          clusterfeatures[[clusterwhich[j]]] <- list(occurence = EWCD_result[[clusterwhich[j]]]$occurence,
                                                     S = EWCD_result[[clusterwhich[j]]]$S,
                                                     N = EWCD_result[[clusterwhich[j]]]$N, 
                                                     WCD = EWCD_result[[clusterwhich[j]]]$WCD)
        }
        clusterwhich
        
      }
      
      
      if (any(cluster_temp != cluster_result)) {
        moved <- moved + sum(cluster_temp != cluster_result)
      }
      
      
      if (i != length(rndseq)) {
        valid_ids <- rnd[rndseq[i]:(rndseq[i + 1] - 1)]
      } else {
        valid_ids <- rnd[rndseq[i]:(length(rnd))]
      }
      
      
      clusters[valid_ids] <- cluster_result
      
      
      # update clusterfeatures with processed cases
      for (j in unique(cluster_result)) {
        clusterfeatures[[j]] <- list(occurence = clusterfeatures[[j]]$occurence + colSums(dataset[valid_ids,][which(cluster_result == j),]),
                                     S = clusterfeatures[[j]]$S + sum(dataset[valid_ids,][which(cluster_result == j),]),
                                     N = clusterfeatures[[j]]$N + nrow(dataset[valid_ids,][which(cluster_result == j),]))
        clusterfeatures[[j]]$WCD <- sum(clusterfeatures[[j]]$occurence^2) / (as.numeric(clusterfeatures[[j]]$S) * as.numeric(clusterfeatures[[j]]$N))
      }
      
    }
    
  }
  
  return(clusters)
}