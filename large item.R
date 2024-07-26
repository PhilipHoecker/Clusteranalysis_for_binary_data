cost <- function(cf_i, x, minsupport, cf, w, cf_length) {
  support_new <- cf[[cf_i]]$support + x
  N_new <- cf[[cf_i]]$N + 1
  minsup_new <- N_new * minsupport
  large_new <- which(support_new >= minsup_new)
  small_new <- which(support_new < minsup_new & support_new > 0)
  
  # if there is only one cluster
  if (cf_length == 1) {
    
    intra <- length(small_new)
    # since there is only one cluster sum(large) - union(large) is always 0
    inter <- 0
    cost_return <- w * intra + inter
    
  } else {
    
    small_cfs <- unique(unlist(lapply(cf[-cf_i], function(y) y$small)))
    intra <- length(union(small_new, small_cfs))
    large_cfs <- unique(unlist(lapply(cf[-cf_i], function(y) y$large)))
    inter <- sum(sapply(cf[-cf_i], function(x) length(x$large)) + length(large_new)) - length(union(large_new, large_cfs))
    cost_return <- w * intra + inter
  }
  
  return(list(support = support_new, N = N_new, minsup = minsup_new, large = large_new, small = small_new, cost_return = cost_return))
}

cost_single <- function(x, minsupport, cf, w, cf_length) {
  support <- x
  N <- 1
  minsup <- minsupport * 1
  large <- which(x == 1)
  small <- integer(0)
  
  if (cf_length == 1) {
    intra <- length(union(cf[[1]]$small, small))
    inter <- (length(cf[[1]]$large) + length(large)) - length(union(cf[[1]]$large, large))
    cost_return <- w * intra + inter
    
  } else {
    small_cfs <- unique(unlist(lapply(cf, function(y) y$small)))
    intra <- length(union(small, small_cfs))
    large_cfs <- unique(unlist(lapply(cf, function(y) y$large)))
    inter <- sum(sapply(cf, function(x) length(x$large)) + length(large)) - length(union(large, large_cfs))
    cost_return <- w * intra + inter
  }
  
  return(list(support = support, N = N, minsup = minsup, large = large, small = small, cost_return = cost_return))
}

cost_cluster <- function(cf_i, x, minsupport, cf, w, cf_length) {
  support_new <- cf[[cf_i]]$support + x$support
  N_new <- cf[[cf_i]]$N + x$N
  minsup_new <- N_new * minsupport
  large_new <- which(support_new >= minsup_new)
  small_new <- which(support_new < minsup_new & support_new > 0)
  
  # if there is only one cluster
  if (cf_length == 1) {
    
    intra <- length(small_new)
    # since there is only one cluster sum(large) - union(large) is always 0
    inter <- 0
    cost_return <- w * intra + inter
    
  } else {
    
    small_cfs <- unique(unlist(lapply(cf[-cf_i], function(y) y$small)))
    intra <- length(union(small_new, small_cfs))
    large_cfs <- unique(unlist(lapply(cf[-cf_i], function(y) y$large)))
    inter <- sum(sapply(cf[-cf_i], function(x) length(x$large)) + length(large_new)) - length(union(large_new, large_cfs))
    cost_return <- w * intra + inter
  }
  
  return(list(support = support_new, N = N_new, minsup = minsup_new, large = large_new, small = small_new, cost_return = cost_return))
}



cost_single_cluster <- function(x, minsupport, cf, w, cf_length) {
  support <- x$support
  N <- x$N
  minsup <- minsupport * x$N
  large <- which(support >= minsup)
  small <- which(support < minsup & support > 0)
  
  if (cf_length == 1) {
    intra <- length(union(cf[[1]]$small, small))
    inter <- (length(cf[[1]]$large) + length(large)) - length(union(cf[[1]]$large, large))
    cost_return <- w * intra + inter
    
  } else {
    small_cfs <- unique(unlist(lapply(cf, function(y) y$small)))
    intra <- length(union(small, small_cfs))
    large_cfs <- unique(unlist(lapply(cf, function(y) y$large)))
    inter <- sum(sapply(cf, function(x) length(x$large)) + length(large)) - length(union(large, large_cfs))
    cost_return <- w * intra + inter
  }
  
  return(list(support = support, N = N, minsup = minsup, large = large, small = small, cost_return = cost_return))
}




largeitem_par <- function(dataset, w, minsupport, maxiteration, maxmoved, ncores) {
  
  
  
  # binaryd <- isplitRows(dataset, chunks = ncores)
  
  
  results_largeitem <- foreach(dataset_iter = isplitRows(dataset, chunks = ncores), .combine = append) %dopar% {
    # clusterids
    clusters <- NULL
    clusters[1] <- 1  
    maxcluster <- 1
    
    # first observation constitutes the first cluster
    clusterfeatures <- list()
    clusterfeatures[[1]] <- list(support = dataset_iter[1,],
                                 N = 1,
                                 minsup = minsupport * 1,
                                 large = which(dataset_iter[1,] == 1),
                                 small = integer(0))
    
    for (i in 2:nrow(dataset_iter)) {
      
      # calculate cost for every cluster if observation is added
      cost_result <- lapply(1:length(clusterfeatures), cost, dataset_iter[i,], minsupport, clusterfeatures, w, cf_length = length(clusterfeatures))
      costs <- sapply(cost_result, function(x) x$cost_return)
      # cost if observations forms new cluster 
      cost_single_result <- cost_single(dataset_iter[i,], minsupport, clusterfeatures, w, length(clusterfeatures))
      costs <- append(costs, cost_single_result$cost_return)
      
      
      clusterwhich <- which.min(costs)
      clusters[i] <- clusterwhich
      
      # if observation gets added to existing cluster, update the clusterfeatures
      if (clusterwhich <= maxcluster) {
        clusterfeatures[[clusterwhich]] <- list(support = cost_result[[clusterwhich]]$support,
                                                N = cost_result[[clusterwhich]]$N,
                                                minsup = cost_result[[clusterwhich]]$minsup, 
                                                large = cost_result[[clusterwhich]]$large,
                                                small = cost_result[[clusterwhich]]$small)
        
        # if observation forms new cluster add new entry to clusterfeatures
      } else {
        clusterfeatures[[clusterwhich]] <- list(support = cost_single_result$support,
                                                N = cost_single_result$N,
                                                minsup = cost_single_result$minsup, 
                                                large = cost_single_result$large,
                                                small = cost_single_result$small)
        
        maxcluster <- clusterwhich
        
      }
    }
    
    maxiter <- 0
    moved <- Inf
    while (maxiter <= maxiteration & moved > maxmoved) {
      maxiter <- maxiter + 1
      moved <- 0
      
      for (i in 1:nrow(dataset_iter)) {
        
        # clus_temp for checking if moved
        clus_temp <- clusters[i]
        
        # subract observation from cluster
        clusterfeatures[[clusters[i]]] <- list(support = clusterfeatures[[clusters[i]]]$support - dataset_iter[i,],
                                               N = clusterfeatures[[clusters[i]]]$N - 1,
                                               minsup = (clusterfeatures[[clusters[i]]]$N - 1) * minsupport, 
                                               large = which(clusterfeatures[[clusters[i]]]$support - dataset_iter[i,] >= (clusterfeatures[[clusters[i]]]$N - 1) * minsupport),
                                               small = which(clusterfeatures[[clusters[i]]]$support - dataset_iter[i,] < (clusterfeatures[[clusters[i]]]$N - 1) * minsupport & 
                                                               clusterfeatures[[clusters[i]]]$support - dataset_iter[i,] > 0))
        
        # calculate cost_result for every clustering if observation is added to any cluster
        cost_result <- lapply(1:length(clusterfeatures), cost, dataset_iter[i,], minsupport, clusterfeatures, w, cf_length = length(clusterfeatures))
        costs <- sapply(cost_result, function(x) x$cost_return)
        
        clusterwhich <- which.min(costs)
        clusters[i] <- clusterwhich
        
        clusterfeatures[[clusterwhich]] <- list(support = cost_result[[clusterwhich]]$support,
                                                N = cost_result[[clusterwhich]]$N,
                                                minsup = cost_result[[clusterwhich]]$minsup, 
                                                large = cost_result[[clusterwhich]]$large,
                                                small = cost_result[[clusterwhich]]$small)
        
        if (clus_temp != clusterwhich) {
          moved <- moved + 1
        }
        
      }
    }
    list(clusterfeatures = clusterfeatures, clusters = clusters)
  }
  
  clusterfeatures <- results_largeitem[seq(1, (2 * ncores) - 1, 2)]
  
  clusters <- list()
  # make clusterlabels unique
  h <- 0
  for (i in 1:length(results_largeitem[seq(2, 2 * ncores, 2)])) {
    results_largeitem[seq(2, 2 * ncores, 2)][[i]] <- results_largeitem[seq(2, 2 * ncores, 2)][[i]] + h
    h <- max(results_largeitem[seq(2, 2 * ncores, 2)][[i]])
    
  }
  
  clusters <- as.integer(unlist(results_largeitem[seq(2, 2 * ncores, 2)]))
  
  clusterf <- list()
  h <- 0
  for (i in 1:length(clusterfeatures)) {
    
    clusterf[(h + 1):(h + length(clusterfeatures[[i]]))] <- clusterfeatures[[i]][1:(length(clusterfeatures[[i]]))]
    h <- h + length(clusterfeatures[[i]]) 
    
  }
  
  clusterfeatures2 <- list()
  clusterfeatures2[[1]] <- clusterf[[1]]
  maxcluster <- 1
  
  
  for (i in 2:length(clusterf)) {
    # calculate cost for every cluster if cluster is added
    cost_result <- lapply(1:length(clusterfeatures2), cost_cluster, clusterf[[i]], minsupport, clusterfeatures2, w, cf_length = length(clusterfeatures2))
    costs <- sapply(cost_result, function(x) x$cost_return)
    # cost if cluster forms new cluster 
    cost_single_result <- cost_single_cluster(clusterf[[i]], minsupport, clusterfeatures2, w, length(clusterfeatures2))
    costs <- append(costs, cost_single_result$cost_return)
    
    
    clusterwhich <- which.min(costs)
    
    
    # if cluster gets added to existing cluster, update the clusterfeatures
    if (clusterwhich <= maxcluster) {
      clusterfeatures2[[clusterwhich]] <- list(support = cost_result[[clusterwhich]]$support,
                                               N = cost_result[[clusterwhich]]$N,
                                               minsup = cost_result[[clusterwhich]]$minsup, 
                                               large = cost_result[[clusterwhich]]$large,
                                               small = cost_result[[clusterwhich]]$small)
      
      clusters[clusters == i] <- unique(clusters)[clusterwhich]
      
      # if cluster forms new cluster add new entry to clusterfeatures
    } else {
      clusterfeatures2[[clusterwhich]] <- list(support = cost_single_result$support,
                                               N = cost_single_result$N,
                                               minsup = cost_single_result$minsup, 
                                               large = cost_single_result$large,
                                               small = cost_single_result$small)
      
      maxcluster <- clusterwhich
      
    }
  }
  
  return(clusters)
}