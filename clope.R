deltaadd <- function(cf, x, r) {
  Snew <- cf$S + sum(x)
  Worignew <- cf$Worig + x
  Wnew <- length(which(Worignew != 0))
  Nnew <- cf$N + 1
  delta <- (Snew * Nnew / Wnew^r) - (cf$S * cf$N / cf$W^r)
  return(list(Snew = Snew, Wnew = Wnew, Worignew = Worignew, Nnew = Nnew, delta = delta))
}

deltaaddcluster <- function(cf, x) {
  # cf <- cf[[1]]
  Snew <- cf$S + x$S
  Worignew <- cf$Worig + x$Worig
  Wnew <- length(which(Worignew != 0))
  Nnew <- cf$N + x$N
  delta <- (Snew * Nnew / Wnew^r) - (cf$S * cf$N / cf$W^r)
  return(list(Snew = Snew, Wnew = Wnew, Worignew = Worignew, Nnew = Nnew, delta = delta))
}

r <- 2

clope_par <- function(dataset, r, ncores, maxiteration, maxmoved, clusterstop) {
  
  # binaryd <- isplitRows(dataset, chunks = ncores)
  
  results_clope <- foreach(dataset_iter = isplitRows(dataset, chunks = ncores), .combine = append) %dopar% {
    # clusterids
    clusters <- NULL
    clusters[1] <- 1  
    maxcluster <- 1
    
    # first observation constitutes the first cluster
    clusterfeatures <- list()
    clusterfeatures[[1]] <- list(S = sum(dataset_iter[1,]), 
                                 W = sum(dataset_iter[1,]),
                                 Worig = dataset_iter[1,],
                                 N = 1)
    
    for (i in 2:nrow(dataset_iter)) {
      # print(i)
      # results for deltaadd
      deltaadd_results <- lapply(clusterfeatures, deltaadd, dataset_iter[i,], r)
      deltas <- sapply(deltaadd_results, function(x) x$delta)
      # results for deltaadd in case of new cluster
      if (length(deltas) < clusterstop) {
        deltanewcluster <- (sum(dataset_iter[i,]) / (sum(dataset_iter[i,]))^r)
        # to which cluster will the observation be added
        clusterwhich <- which.max(c(deltas, deltanewcluster))
      } else {
        # to which cluster will the observation be added
        clusterwhich <- which.max(deltanewcluster)
      }
      
      # in case of new cluster, add new cluster to clusterfeatures and clusters
      if (clusterwhich > maxcluster) {
        clusterfeatures[[clusterwhich]] <- list(S = sum(dataset_iter[i,]), 
                                                W = sum(dataset_iter[i,]),
                                                Worig = dataset_iter[i,],
                                                N = 1)
        clusters[i] <- clusterwhich
        maxcluster <- clusterwhich
      } else {
        # allocate observation to cluster
        clusters[i] <- clusterwhich
        # update clusterfeatures
        clusterfeatures[[clusterwhich]] <- list(S = deltaadd_results[[clusterwhich]]$Snew,
                                                W = deltaadd_results[[clusterwhich]]$Wnew,
                                                Worig = deltaadd_results[[clusterwhich]]$Worignew, 
                                                N = deltaadd_results[[clusterwhich]]$Nnew)
      }
    }
    
    maxiter <- 0
    moved <- Inf
    while (maxiter <= maxiteration & moved > maxmoved) {
      maxiter <- maxiter + 1
      moved <- 0
      # print(maxiter)
      
      for (i in 1:nrow(dataset_iter)) {
        
        # clus_temp for checking if moved
        clus_temp <- clusters[i]
        
        # subtract case i from clusterfeature
        clusterfeatures[[clusters[i]]]$S <- clusterfeatures[[clusters[i]]]$S - sum(dataset_iter[i,])
        clusterfeatures[[clusters[i]]]$Worig <- clusterfeatures[[clusters[i]]]$Worig - dataset_iter[i,]
        clusterfeatures[[clusters[i]]]$W <- length(which(clusterfeatures[[clusters[i]]]$Worig != 0))
        clusterfeatures[[clusters[i]]]$N <- clusterfeatures[[clusters[i]]]$N - 1
        
        # if cf is empty now -> delete
        if (clusterfeatures[[clusters[i]]]$N == 0) {
          clusterfeatures[[clusters[i]]] <- NULL
          # subtract 1 from clusters higher than clusters[i] because the cluster was deleted
          clusters[clusters > clusters[i]] <- clusters[clusters > clusters[i]] - 1
        }
        
        
        # results for deltaadd
        deltaadd_results <- lapply(clusterfeatures, deltaadd, dataset_iter[i,], r)
        deltas <- sapply(deltaadd_results, function(x) x$delta)
        # results for deltaadd in case of new cluster
        if (length(deltas) < clusterstop) {
          deltanewcluster <- (sum(dataset_iter[i,]) / (sum(dataset_iter[i,]))^r)
          # to which cluster will the observation be added
          clusterwhich <- which.max(c(deltas, deltanewcluster))
        } else {
          # to which cluster will the observation be added
          clusterwhich <- which.max(deltanewcluster)
        }
        
        # check if moved
        if (clusterwhich != clus_temp) {
          moved <- moved + 1
        }
        
        # in case of new cluster, add new cluster to clusterfeatures and clusters
        if (clusterwhich > length(clusterfeatures)) {
          clusterfeatures[[clusterwhich]] <- list(S = sum(dataset_iter[i,]), 
                                                  W = sum(dataset_iter[i,]),
                                                  Worig = dataset_iter[i,],
                                                  N = 1)
          clusters[i] <- clusterwhich
          maxcluster <- clusterwhich
        } else {
          # allocate observation to cluster
          clusters[i] <- clusterwhich
          # update clusterfeatures
          clusterfeatures[[clusterwhich]] <- list(S = deltaadd_results[[clusterwhich]]$Snew,
                                                  W = deltaadd_results[[clusterwhich]]$Wnew,
                                                  Worig = deltaadd_results[[clusterwhich]]$Worignew, 
                                                  N = deltaadd_results[[clusterwhich]]$Nnew)
        }
      }
      # print(moved)
    }
    list(clusterfeatures = clusterfeatures, clusters = clusters)
  }
  
  clusterfeatures <- results_clope[seq(1, (2 * ncores) - 1, 2)]
  # clusters <- as.integer(unlist(results_clope[seq(2, 24, 2)]))
  
  clusters <- list()
  # make clusterlabels unique
  h <- 0
  for (i in 1:length(results_clope[seq(2, 2 * ncores, 2)])) {
    results_clope[seq(2, 2 * ncores, 2)][[i]] <- results_clope[seq(2, 2 * ncores, 2)][[i]] + h
    h <- max(results_clope[seq(2, 2 * ncores, 2)][[i]])
    
  }
  
  clusters <- as.integer(unlist(results_clope[seq(2, 2 * ncores, 2)]))
  
  
  
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
    # for (i in 2:14) {
    
    deltaadd_results <- lapply(clusterfeatures2, deltaaddcluster, clusterf[[i]])
    deltas <- sapply(deltaadd_results, function(x) x$delta)
    
    if (length(deltas) < clusterstop) {
      deltanewcluster <- (clusterf[[i]]$S * clusterf[[i]]$N) / ((clusterf[[i]]$W)^r)
      # to which cluster will the observation be added
      clusterwhich <- which.max(c(deltas, deltanewcluster))
    } else {
      # to which cluster will the observation be added
      clusterwhich <- which.max(deltanewcluster)
    }
    
    # in case of new cluster, add new cluster to clusterfeatures2
    if (clusterwhich > maxcluster) {
      clusterfeatures2[[clusterwhich]] <- clusterf[[i]]
      
      # clusters[i] <- clusterwhich
      maxcluster <- clusterwhich
    } else {
      # allocate observation to cluster
      clusters[clusters == i] <- unique(clusters)[clusterwhich]
      # update clusterfeatures
      clusterfeatures2[[clusterwhich]] <- list(S = deltaadd_results[[clusterwhich]]$Snew,
                                               W = deltaadd_results[[clusterwhich]]$Wnew,
                                               Worig = deltaadd_results[[clusterwhich]]$Worignew, 
                                               N = deltaadd_results[[clusterwhich]]$Nnew)
    }
    
  }
  
  
  return(clusters)
}