sorted_neighborhood <- function(dataset, windowsize, threshold = mean(pdist(dataset[sample(1:nrow(dataset), 100),], metric = "hamming")),
                                runs, allrandom = TRUE, ncores) {
  
  rownames(dataset) <- 1:nrow(dataset)
  clusids <- list()
  
  run <- 1
  while (run <= runs) {
    if (run == 1 && !allrandom) {
      dataset <- dataset[order(rowSums(dataset)),]
    } else {
      dataset <- dataset[sample(1:nrow(dataset), nrow(dataset)),]
    }
    
    # iterdata <- isplitRows(dataset, chunks = ncores)
    
    clusids_par <- foreach(dataset_iter = isplitRows(dataset, chunks = ncores), .combine = append) %dopar% {
      clusters <- 1
      dist_i <- rdist:::hamming_cdist(dataset_iter[2,, drop = FALSE], dataset_iter[1,, drop = FALSE])
      which_cases <- which(dist_i < threshold)
      clusters[2] <- ifelse(!length(which_cases), 2, clusters[which_cases[1]])
      
      for (i in 3:windowsize) {
        dist_i <- rdist:::hamming_cdist(dataset_iter[i,, drop = FALSE], dataset_iter[1:(i - 1),])
        which_case_min <- which.min(dist_i)
        tresh <- dist_i[which_case_min] < threshold
        clusters[i] <- ifelse(tresh, clusters[which_case_min], i)
      }
      
      for (i in (windowsize + 1):nrow(dataset_iter)) {
        dist_i <- rdist:::hamming_cdist(dataset_iter[i,, drop = FALSE], dataset_iter[(i - windowsize):(i - 1),])
        which_case_min <- which.min(dist_i)
        tresh <- dist_i[which_case_min] < threshold
        clusters[i] <- ifelse(tresh, clusters[which_case_min + (i - windowsize - 1)], i)
      }
      
      clusterids <- list()
      h <- 0
      for (i in unique(clusters)) {
        h <- h + 1
        clusterids[[h]] <- rownames(dataset_iter[clusters == i,, drop = FALSE])
      }
      clusterids
    }
    
    
    
    clusids <- append(clusids, clusids_par)
    run <- run + 1
  }
  
  # transitive hull ---------------------------------------------------------
  
  # remove singletons (not needed for transitive hull)
  clusids_without_singletons <- clusids[sapply(clusids, length) > 1]
  
  transmat <- matrix(nrow = length(clusids_without_singletons), ncol = length(clusids_without_singletons))
  
  transmat <- foreach (i = 1:length(clusids_without_singletons), .combine = rbind, .export = "transmat") %dopar% {
    for (j in i:length(clusids_without_singletons)) {
      if (i != j && any(clusids_without_singletons[[i]] %in% clusids_without_singletons[[j]])) {
        transmat[i, j] <- 1
      }
    }
    transmat[i,]
  }
  
  
  transmat <- as.matrix(forceSymmetric(transmat))
  
  
  transhull <- allShortestPaths(transmat)$length
  
  # union of the ids that belong together according to transitive hull
  unionids <- list()
  donelist <- NULL
  h <- 0
  for (i in 1:ncol(transhull)) {
    if (!i %in% donelist) {
      h <- h + 1
      unionids[[h]] <- unique(unlist(clusids_without_singletons[which(transhull[, i] > 0)]))
      donelist <- append(donelist, which(transhull[, i] > 0))
    }
  }
  
  # add the singletons
  unionids <- append(unionids, clusids[sapply(clusids, length) == 1])
  
  
  for (i in 1:length(unionids)) {
    for (j in i:length(unionids)) {
      if (i != j && length(unionids[[j]]) == 1) {
        if (unionids[[j]] %in% unionids[[i]]) {
          unionids[[j]] <- "0"
        }
      } 
    }
  }
  
  
  unionids[which(sapply(unionids, function(x) x) == "0")] <- NULL
  
  unionids <- lapply(unionids, as.integer)
  clusters <- NULL
  for (i in 1:length(unionids)) {
    clusters[unionids[[i]]] <- i
  }
  
  return(clusters)
}


# cross

cross_sorted <- function(dataset, ncluster, orig, ncores) {
  
  n <- nrow(dataset)
  # bound <- round(n / ncluster)
  # bounds <- seq(0, n, bound)
  
  datasample <- matrix(nrow = 2000, ncol = ncol(dataset))
  boundsample <- round(2000 / ncluster)
  # boundssample <- seq(0, 2000, boundsample)
  # boundssample[length(boundssample)] <- 2000
  boundssample <- round(seq(0, 2000, length.out = ncluster + 1))
  
  for (i in 1:ncluster) {
    # which_sample <- sample(1:round(n / ncluster), round(2000 / ncluster))
    which_sample <- sample(1:floor(n / ncluster), boundssample[i + 1] - boundssample[i])
    datasample[(boundssample[i] + 1):boundssample[i + 1],] <- dataset[which(orig == i),][which_sample,]
  }
  
  shuf_index <- sample(1:nrow(datasample), nrow(datasample))
  datasample <- datasample[shuf_index,]
  
  orig_sample <- integer()
  for (i in 1:ncluster) {
    orig_sample <- c(orig_sample, rep(i, round(nrow(datasample) / ncluster)))
  }
  
  orig_sample <- orig_sample[shuf_index]
  
  thresholdlist <- seq(0.3, 0.8, 0.1)
  
  randlist <- NULL
  
  for (i in 1:length(thresholdlist)) {
    sorted_neighborhood_result <- sorted_neighborhood(datasample, windowsize = 30, thresholdlist[i], runs = 3, allrandom = TRUE, ncores)
    randlist[i] <- rand(sorted_neighborhood_result, orig_sample, nrow(datasample))
  }
  
  threshold <- thresholdlist[which.max(randlist)]
  return(threshold)
}