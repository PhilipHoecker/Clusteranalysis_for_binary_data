canopy_cluster <- function(dataset, t1, t2) {
  int2 <- sample(1:nrow(dataset), nrow(dataset))
  canopys <- list()
  centers <- NULL
  
  h <- 1
  
  # if legnth(int2) == 0 every observation is in canopy
  while (length(int2) > 0) {
    # i is random observation
    i <- int2[1]
    # dist from center to all eligible observations
    dmat <- rdist:::hamming_cdist(dataset[i,, drop = FALSE], dataset[int2,, drop = FALSE])
    colnames(dmat) <- int2
    # canopys is list with canopys as list elements, which conist of the indezes of the observation inside 
    canopys[[h]] <- as.integer(colnames(dmat)[which(dmat < t1)])
    h <- h + 1
    # which observation must be removed
    int2 <- int2[-which(int2 %in% as.integer(colnames(dmat)[which(dmat < t2)]))]
  }
  
  
  clusters <- NULL
  for (i in 1:length(canopys)) {
    clusters[canopys[[i]]] <- i
  }
  return(clusters)
}


canopy_cluster_par <- function(dataset, t1, t2, ncluster, ncores) {
  
  # ii <- isplitRows(dataset, chunks = ncores)
  
  cano <- foreach(dataset_iter = isplitRows(dataset, chunks = ncores)) %dorng% {
    
    int2 <- sample(1:nrow(dataset_iter), nrow(dataset_iter))
    canopys <- list()
    centers <- NULL
    h <- 1
    while (length(int2) > 0) {
      # i is random observation
      i <- int2[1]
      # dist from center to all eligible observations
      if (h == 1) {
        dmat <- rdist:::hamming_cdist(dataset_iter[i,, drop = FALSE], dataset_iter)
        colnames(dmat) <- 1:length(dmat)
      } else {
        dmat <- rdist:::hamming_cdist(dataset_iter[i,, drop = FALSE], dataset_iter[int2,, drop = FALSE])
        colnames(dmat) <- int2
      }
      # canopys is list with canopys as list elements, which conist of the indezes of the observation inside
      canopys[[h]] <- as.integer(colnames(dmat)[which(dmat < t1)])
      names(canopys)[h] <- as.character(i)
      h <- h + 1
      int2 <- int2[-which(int2 %in% as.integer(colnames(dmat)[which(dmat < t2)]))]
    }
    canopys
  }
  
  chunksizes <- sapply(cano, function(x) {
    y <- NULL
    for (i in 1:length(x)) {
      y[i] <- length(x[[i]])
    }
    return(sum(y))
  })
  
  
  chunksizes_cum <- cumsum(chunksizes)
  chunksizes_cum <- chunksizes_cum[-length(chunksizes_cum)]
  
  for (i in 2:length(cano)) {
    cano[[i]] <- rapply(cano[[i]], function(x, y) x + y, how = "replace", y = chunksizes_cum[i - 1])
  }
  
  centers <- as.integer(unlist(lapply(cano, names)))
  
  clustered_canopies <- kmeans(dataset[centers,], ncluster)$cluster
  
  
  clusters <- integer(nrow(dataset))
  h <- 1
  for (i in 1:length(cano)) {
    if (length(cano[[i]]) == 1) {
      clusters[unlist(cano[[i]])] <- clustered_canopies[h]
      h <- h + 1
    } else {
      for (j in 1:length(cano[[i]])) {
        clusters[cano[[i]][[j]]] <- clustered_canopies[h]
        h <- h + 1
      }
    }
    
  }
  return(clusters)
}


# cross

cross_canopy <- function(dataset, ncluster, orig) {
  n <- nrow(dataset)
  
  datasample <- matrix(nrow = 2000, ncol = ncol(dataset))
  boundsample <- round(2000 / ncluster)
  boundssample <- round(seq(0, 2000, length.out = ncluster + 1))
  
  for (i in 1:ncluster) {
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
  
  thresholdlist <- seq(0.3, 0.6, 0.05)
  
  randlist <- NULL
  
  for (i in 1:length(thresholdlist)) {
    canopy_result <- canopy_cluster(datasample, thresholdlist[i], thresholdlist[i])
    randlist[i] <- rand(canopy_result, orig_sample, nrow(datasample))
  }
  
  threshold <- thresholdlist[which.max(randlist)]
  return(threshold)
}