g_mod_bin_data2 <- function(dataset, ncluster) {
  
  if (!is.matrix(dataset)) {
    dataset <- as.matrix(dataset)
  }
  
  
  A <- sparseMatrix(1:nrow(dataset), sample(1:ncluster, nrow(dataset), replace = TRUE))
  
  nk <- colSums(A)
  
  
  Y <- matrix(nrow = ncluster, ncol = ncol(dataset))
  for (i in seq_len(ncluster)) {
    dd <- matrix(A[, i], nrow = nrow(dataset), ncol = ncol(dataset))
    Y[i,] <- (1/nk[i]) * colSums(dd * dataset)
  }
  
  
  B <- matrix(nrow = ncluster, ncol = ncol(dataset))
  # B[Y > 0.5] <- 1
  # B[Y <= 0.5] <- 0
  B[Y > 1 / ncluster] <- 1L
  B[Y <= 1 / ncluster] <- 0L
  
  # objective criterion
  
  Oab1 <- sum(sqrt((dataset - A %*% B)^2))
  
  Oab0 <- Inf
  
  while (Oab1 < Oab0) {
    
    Oab0 <- Oab1
    
    
    Abkj <- matrix(nrow = nrow(dataset), ncol = ncluster)
    for (j in 1:ncluster) {
      BJ <- matrix(rep(B[j,], each = nrow(dataset)), nrow = nrow(dataset))
      Abkj[, j] <- rowSums((dataset - BJ)^2)
    }
    
    whichcluster <- apply(Abkj, 1, which.min)
    A <- matrix(0L, nrow = nrow(dataset), ncol = ncluster)
    for (i in 1:ncluster) {
      A[which(whichcluster == i), i] <- 1L
    }
    
    # compute B
    nk <- colSums(A)
    if (any(nk == 0)) {
      nk[which(nk == 0)] <- 1L
    }
    
    Y <- matrix(nrow = ncluster, ncol = ncol(dataset))
    for (i in seq_len(ncluster)) {
      dd <- matrix(A[, i], nrow = length(A[, i]), ncol = ncol(dataset))
      Y[i,] <- (1/nk[i]) * colSums(dd * dataset)
    }
    
    B <- matrix(nrow = ncluster, ncol = ncol(dataset))
    B[Y > 1 / ncluster] <- 1L
    B[Y <= 1 / ncluster] <- 0L
    
    # objective criterion. not too sure. rolling with it for now
    
    Oab1 <- sum(sqrt((dataset - A %*% B)^2))
    
  }
  
  clusters <- apply(A, 1, which.max)
  clusters
}


g_mod_bin_data2_par <- function(dataset, ncluster, ncores) {
  if (!is.matrix(dataset)) {
    dataset <- as.matrix(dataset)
  }
  # ii <- isplitRows(dataset, chunks = ncores)
  
  results <- foreach(dat = isplitRows(dataset, chunks = ncores)) %dopar% {
    g_mod_bin_data2(dat, ncluster)
  }
  
  # make clusterids unique
  maxclusid <- max(unique(results[[1]]))
  for (i in 2:length(results)) {
    results[[i]] <- results[[i]] + maxclusid
    maxclusid <- max(unique(results[[i]]))
  }
  
  results <- unlist(results)
  
  
  medoids <- matrix(nrow = max(results), ncol = ncol(dataset))
  for (i in unique(results)) {
    clustersample <- which(results == i)
    clustersample <- dataset[clustersample,, drop = FALSE][sample(1:length(clustersample), ifelse(length(clustersample) < 1000, length(clustersample), 1000)),, drop = FALSE]
    cs_dist <- pdist(clustersample, metric = "hamming")
    medoids[i,] <- clustersample[which.min(rowSums2(cs_dist)),]
  }
  
  medoids_diss <- pdist(medoids, metric = "hamming")
  result_medoids <- cluster::pam(medoids_diss, ncluster, cluster.only = TRUE)
  
  # result_medoids <- Cluster_Medoids(medoids, ncluster, distance_metric = "hamming", threads = 12L)$clusters
  
  for (i in seq_len(ncluster)) {
    results[which(results %in% which(result_medoids == i))] <- i
  }
  results
}


g_mod_bin_data1 <- function(dataset, ncluster, nclusterf) {
  
  if (!is.matrix(dataset)) {
    dataset <- as.matrix(dataset)
  }
  
  # initialize A
  A <- matrix(0L, nrow = nrow(dataset), ncol = ncluster)
  clusterseq <- floor(seq(0, nrow(dataset), length.out = ncluster + 1))
  for (i in 1:ncluster) {
    A[(clusterseq[i] + 1):clusterseq[i + 1], i] <- 1L
  }
  
  
  # initialize B
  B <- matrix(0L, nrow = nclusterf, ncol = ncol(dataset))
  clusterseq <- floor(seq(0, ncol(dataset), length.out = nclusterf + 1))
  for (i in 1:nclusterf) {
    B[i, (clusterseq[i] + 1):clusterseq[i + 1]] <- 1L
  }
  
  
  
  # compute X using eq 5
  
  pk <- colSums2(A)
  qc <- rowSums2(B)
  
  # if there are zeros, they will be replaced with a low number (cant devide by 0 in next step)
  if (any(pk == 0)) {
    pk[pk == 0] <- 0.00000001
  }
  
  if (any(qc == 0)) {
    qc[qc == 0] <- 0.00000001
  }
  
  X <- matrix(nrow = ncluster, ncol = nclusterf)
  for (i in seq_len(ncluster)) {
    for (j in seq_len(nclusterf)) {
      X[i, j] <- (1 / (pk[i] * qc[j])) * sum(dataset[which(A[, i] == 1), which(B[j,] == 1)])
    }
  }
  
  
  What <- A %*% X %*% B
  Oaxb1 <- norm(dataset - What, type = "F")
  
  Oaxb0 <- Inf
  
  
  while (Oaxb1 < Oaxb0) {
    
    Oaxb0 <- Oaxb1
    
    # update A using eq 6
    Aik <- matrix(nrow = nrow(dataset), ncol = ncluster)
    for (i in 1:ncluster) {
      for (j in 1:nclusterf) {
        Aik[, i] <- rowSums2((dataset[, which(B[j, ] == 1), drop = FALSE] - X[i, j])^2)
      }
    }
    
    whichcluster <- max.col(-Aik, ties = "first")
    A <- matrix(0L, nrow = nrow(dataset), ncol = ncluster)
    for (i in 1:ncluster) {
      A[which(whichcluster == i), i] <- 1L
    }
    
    # update B using eq 6
    Bjc <- matrix(nrow = nclusterf, ncol = ncol(dataset))
    for (i in 1:ncluster) {
      for (j in 1:nclusterf) {
        Bjc[j,] <- colSums2((dataset[which(A[, i] == 1), , drop = FALSE] - X[i, j])^2)
      }
    }
    
    # whichcluster <- apply(Bjc, 2, which.min)
    whichcluster <- max.col(t(-Bjc), ties = "first")
    B <- matrix(0L, nrow = nclusterf, ncol = ncol(dataset))
    for (j in 1:nclusterf) {
      B[j, which(whichcluster == j)] <- 1L
    }
    
    # compute X using eq 5
    pk <- colSums2(A)
    qc <- rowSums2(B)
    
    # if there are zeros, they will be replaced with a low number (cant devide by 0 in next step)
    if (any(pk == 0)) {
      pk[pk == 0] <- 0.00000001
    }
    
    if (any(qc == 0)) {
      qc[qc == 0] <- 0.00000001
    }
    
    X <- matrix(nrow = ncluster, ncol = nclusterf)
    for (i in seq_len(ncluster)) {
      for (j in seq_len(nclusterf)) {
        X[i, j] <- (1 / (pk[i] * qc[j])) * sum(dataset[which(A[, i] == 1), which(B[j,] == 1)])
      }
    }
    
    What <- A %*% X %*% B
    Oaxb1 <- norm(dataset - What, type = "F")
    
  }
  
  clusters <- max.col(A, ties = "first")
  clusters
}

g_mod_bin_data1_par <- function(dataset, ncluster, nclusterf, ncores) {
  if (!is.matrix(dataset)) {
    dataset <- as.matrix(dataset)
  }
  # ii <- isplitRows(dataset, chunks = ncores)
  
  results <- foreach(dat = isplitRows(dataset, chunks = ncores)) %dopar% {
    g_mod_bin_data1(dat, ncluster, nclusterf)
  }
  
  # make clusterids unique
  maxclusid <- max(unique(results[[1]]))
  for (i in 2:length(results)) {
    results[[i]] <- results[[i]] + maxclusid
    maxclusid <- max(unique(results[[i]]))
  }
  
  results <- unlist(results)
  
  
  medoids <- matrix(nrow = max(results), ncol = ncol(dataset))
  for (i in unique(results)) {
    clustersample <- which(results == i)
    clustersample <- dataset[clustersample,, drop = FALSE][sample(1:length(clustersample), ifelse(length(clustersample) < 1000, length(clustersample), 1000)),, drop = FALSE]
    cs_dist <- pdist(clustersample, metric = "hamming")
    medoids[i,] <- clustersample[which.min(rowSums2(cs_dist)),]
  }
  
  medoids_diss <- pdist(medoids, metric = "hamming")
  result_medoids <- cluster::pam(medoids_diss, ncluster, cluster.only = TRUE)
  
  
  for (i in seq_len(ncluster)) {
    results[which(results %in% which(result_medoids == i))] <- i
  }
  results
}