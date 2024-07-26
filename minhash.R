minhash <- function(dataset, n_hash, ncores) {
  
  n <- nrow(dataset)
  ncols <- ncol(dataset)
  # hashmat <- foreach(i = 1:n_hash, .combine = cbind) %dorng% {
  #   sample(1:ncols, ncols, replace = FALSE)
  # }
  
  pri <- generate_primes(min = ncols + 1, ncols + 1000)
  p <- pri[2]
  hashmat <- foreach(i = 1:n_hash, .combine = cbind) %dorng% {
    a <- sample(seq(1, p - 1, 2), 1)
    b1 <- sample(0:(p - 1), 1)
    m <- ncols + 1
    ((a*(1:ncols)+b1) %% p) %% m
  }
  
  
  # positions of 1s in every observation
  whichhashmat <- apply(dataset, 1, function(x) which(x == 1))
  
  
  ii <- isplitVector(whichhashmat, chunks = ncores)
  
  
  sigmat <- foreach(i = ii, .combine = cbind, .packages = "matrixStats" ) %dopar% {
    sapply(i, function(x) colMins(hashmat[x,, drop = FALSE]))
  }
  
  colnames(sigmat) <- 1:ncol(sigmat)
  ii <- isplitCols(sigmat, chunks = ncores)
  
  buckets <- foreach(dat = ii) %dopar% {
    data.table(h = apply(dat, 2, fastdigest), id = as.integer(colnames(dat)))
  }
  
  buckets <- rbindlist(buckets)
  buckets[, clusters := .GRP, by = h]
  return(buckets$clusters)
}



# cross minhash

cross_minhash <- function(dataset, ncluster, orig, ncores) {
  
  n <- nrow(dataset)
  # bound <- round(n / ncluster)
  # bounds <- seq(0, n, bound)
  
  datasample <- matrix(nrow = 2000, ncol = ncol(dataset))
  boundsample <- round(2000 / ncluster)
  # boundssample <- seq(0, 2000, boundsample)
  boundssample <- round(seq(0, 2000, length.out = ncluster + 1))
  for (i in 1:ncluster) {
    datasample[(boundssample[i] + 1):boundssample[i + 1],] <- dataset[which(orig == i),][sample(1:length(which(orig == i)), length((boundssample[i] + 1):boundssample[i + 1])),]
  }
  
  shuf_index <- sample(1:nrow(datasample), nrow(datasample))
  datasample <- datasample[shuf_index,]
  
  orig_sample <- integer()
  for (i in 1:ncluster) {
    orig_sample <- c(orig_sample, rep(i, round(nrow(datasample) / ncluster)))
  }
  
  orig_sample <- orig_sample[shuf_index]
  
  n_hash_list <- seq(2, 50, 2)
  bestlist <- NULL
  h <- 0
  for (i in n_hash_list) {
    h <- h + 1
    minhash_result <- minhash(datasample, i, ncores)
    bestlist[h] <- rec_pre_f1(orig_sample, minhash_result, nrow(datasample))$f1
  }
  
  whichbest <- n_hash_list[which.max(bestlist)]
  
  return(whichbest)
}