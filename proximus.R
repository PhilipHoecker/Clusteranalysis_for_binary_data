prox <- function(dataset, maxradius, ncores, ncluster) {
  
  dataset <- as.data.table(dataset)
  for (i in 1:ncol(dataset)) {
    set(dataset, j = i, value = as.logical(dataset[[i]]))
  }
  dataset <- as.matrix(dataset)
  
  # ii <- isplitRows(dataset, chunks = ncores)
  
  par_result <- foreach(dat = isplitRows(dataset, chunks = ncores)) %dorng% {
    pr <- cba::proximus(dat, max.radius = maxradius)
    fits <- as.integer(fitted(pr)$pl)
    patterns <- lapply(pr$a, function(x) x$y)
    list(pr = fits, patterns = patterns)
  }
  print(paste("prox par done", Sys.time()))
  
  patterns <- lapply(par_result, function(x) x$patterns)
  patterns <- unlist(patterns, recursive = FALSE)
  
  patternmat <- matrix(0, nrow = length(patterns), ncol = ncol(dataset))
  
  for (i in 1:length(patterns)) {
    patternmat[i, patterns[[i]]] <- 1
  }
  
  par_result <- lapply(par_result, function(x) x$pr)
  
  # make clusterids unique
  maxclusid <- max(unique(par_result[[1]]))
  for (i in 2:length(par_result)) {
    par_result[[i]] <- par_result[[i]] + maxclusid
    maxclusid <- max(unique(par_result[[i]]))
  }
  
  clusters <- unlist(par_result)
  
  result_cluster <- "more cluster centers than distinct data points"
  clustercenters <- ncluster
  print(paste("prox while", Sys.time()))
  while (is.character(result_cluster)) {
    result_cluster <- tryCatch(kmeans(patternmat, clustercenters)$cluster,
                               error = function(x) return("more cluster centers than distinct data points"))
    clustercenters <- clustercenters - 1
  }
  print(paste("prox while done", Sys.time()))
  rm(clustercenters)
  
  
  
  clusters <- as.data.table(clusters)
  
  print(paste("prox clusters", Sys.time()))
  clusters[, clusters := result_cluster[clusters[1]], by = clusters]
  print(paste("prox clusters done", Sys.time()))
  
  return(clusters)
}

cross_prox <- function(dataset, ncluster, orig) {
  
  n <- nrow(dataset)
  
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
  
  datasample <- as.data.table(datasample)
  for (i in 1:ncol(datasample)) {
    set(datasample, j = i, value = as.logical(datasample[[i]]))
  }
  datasample <- as.matrix(datasample)
  
  maxradiusmean <- round(mean(rowSums(datasample)))
  
  maxradiusmean10 <- maxradiusmean / 10
  
  maxradiuslist <- c(seq(maxradiusmean - 5 * maxradiusmean10, maxradiusmean, maxradiusmean10), 
                     seq(maxradiusmean + maxradiusmean10, maxradiusmean + 5 * maxradiusmean10, maxradiusmean10))
  
  randlist <- NULL
  
  for (i in 1:length(maxradiuslist)) {
    prox_result <- cba::proximus(datasample, max.radius = maxradiuslist[i])
    prox_result <- as.integer(fitted(prox_result)$pl)
    randlist[i] <- rand(prox_result, orig_sample, nrow(datasample))
  }
  
  maxradius <- maxradiuslist[which.max(randlist)]
  return(maxradius)
}