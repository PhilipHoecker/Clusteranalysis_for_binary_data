find_balanced_variable <- function(dataset) {
  paste0("V", as.character(which.min(abs(dataset[, lapply(.SD, sum)] - dataset[, .N] / 2))))
}

mbt <- function(dataset, threshold, nclusters) {
  size <- dataset[, .N]
  partitions <- data.table()
  nrow_partition <- threshold + 1
  
  central_variable <- find_balanced_variable(select(dataset, starts_with("V")))
  dataset[, partition := fcase(get(central_variable) == 0, 0,
                               get(central_variable) == 1, 1)]
  
  stopit <- FALSE
  nrow_partition2 <- 1
  while (any(nrow_partition >= threshold)) {
    n_prior <- dataset[, .N]
    dataset <- split(dataset, by = "partition")
    ii <- iter(dataset)
    dataset <- foreach(dataset = ii) %do% {
      central_variable <- find_balanced_variable(select(dataset, starts_with("V")))
      dataset[, partition := fcase(get(central_variable) == 0, sample(1:1e10, 1),
                                   get(central_variable) == 1, sample(1:1e10, 1))]
    }
    
    dataset <- rbindlist(dataset)
    nrow_partition <- dataset[, .N, by = partition]
    
    # if central variable is constant, it will be an infinite loop, because it doesnt get split nrow_partition[, .N] will be 1
    # this will stop it
    if (nrow_partition[, .N] == 1) {
      stopit <- TRUE
    }
    
    if (any(nrow_partition$N <= threshold)) {
      partitions <- rbindlist(list(partitions, dataset[partition %in% nrow_partition[N <= threshold, partition]]))
      dataset <- dataset[partition %in% nrow_partition[N > threshold, partition]]
    }
    nrow_partition <- dataset[, .N, by = partition]$N
    
    # if central variable is constant, it will be an infinite loop, because it doesnt get split. if more than one partition does not get split anymore
    # this will stop it. the remaining partitions will then be clusters
    if (identical(nrow_partition, nrow_partition2)) {
      stopit <- TRUE
    }
    
    # needed to check if it does not split
    nrow_partition2 <- nrow_partition
    
    
    if (stopit) {
      nrow_partition <- threshold - 1
    }
    
  }
  
  
  # in case of central variable is constant, the remaining data will not get split and regarded as one cluster
  # needs to be rbinded to partitions
  if (partitions[, .N] < size) {
    dataset[, partition := sample(1:1e10, 1), by = partition]
    partitions <- rbindlist(list(partitions, dataset))
  }
  
  # set.seed(1)
  medoids <- partitions[, .SD[sample(1:.N, 1)], by = partition]
  result <- kmeans(medoids, nclusters)$cluster
  
  clustertable <- data.table(medoids = medoids$partition, clusters = result)
  
  partitions <- merge(partitions, clustertable, by.x = "partition", by.y = "medoids")
  
  return(partitions$clusters)
}

# parallel

mbt_par <- function(dataset, threshold, nclusters) {
  
  # dataiter <- isplitRows(dataset, chunks = ncores)
  
  ergebnis <- foreach(dataset_iter = isplitRows(dataset, chunks = ncores)) %dorng% {
    size <- dataset_iter[, .N]
    partitions <- data.table()
    nrow_partition <- threshold + 1
    
    central_variable <- find_balanced_variable(select(dataset_iter, starts_with("V")))
    dataset_iter[, partition := fcase(get(central_variable) == 0, 0,
                                      get(central_variable) == 1, 1)]
    
    stopit <- FALSE
    nrow_partition2 <- 1
    while (any(nrow_partition >= threshold)) {
      n_prior <- dataset_iter[, .N]
      dataset_iter <- split(dataset_iter, by = "partition")
      ii <- iter(dataset_iter)
      dataset_iter <- foreach(dataset_iter2 = ii) %do% {
        central_variable <- find_balanced_variable(select(dataset_iter2, starts_with("V")))
        dataset_iter2[, partition := fcase(get(central_variable) == 0, sample(1:1e10, 1),
                                           get(central_variable) == 1, sample(1:1e10, 1))]
      }
      
      dataset_iter <- rbindlist(dataset_iter)
      nrow_partition <- dataset_iter[, .N, by = partition]
      
      # if central variable is constant, it will be an infinite loop, because it doesnt get split nrow_partition[, .N] will be 1
      # this will stop it
      if (nrow_partition[, .N] == 1) {
        stopit <- TRUE
      }
      
      if (any(nrow_partition$N <= threshold)) {
        partitions <- rbindlist(list(partitions, dataset_iter[partition %in% nrow_partition[N <= threshold, partition]]))
        dataset_iter <- dataset_iter[partition %in% nrow_partition[N > threshold, partition]]
      }
      nrow_partition <- dataset_iter[, .N, by = partition]$N
      
      # if central variable is constant, it will be an infinite loop, because it doesnt get split. if more than one partition does not get split anymore
      # this will stop it. the remaining partitions will then be clusters
      if (identical(nrow_partition, nrow_partition2)) {
        stopit <- TRUE
      }
      
      # needed to check if it does not split
      nrow_partition2 <- nrow_partition
      
      
      if (stopit) {
        nrow_partition <- threshold - 1
      }
      
    }
    
    
    # in case of central variable is constant, the remaining data will not get split and regarded as one cluster
    # needs to be rbinded to partitions
    if (partitions[, .N] < size) {
      dataset_iter[, partition := sample(1:1e10, 1), by = partition]
      partitions <- rbindlist(list(partitions, dataset_iter))
    }
    
    # set.seed(1)
    medoids <- partitions[, .SD[sample(1:.N, 1)], by = partition]

    n_partitions <- length(unique(partitions$partition))
    if (n_partitions > nclusters) {
      result <- kmeans(medoids[, -c("partition")], nclusters)$cluster
      clustertable <- data.table(medoids = medoids$partition, clusters = result)
    } else {
      clustertable <- data.table(medoids = medoids$partition, clusters = 1:n_partitions)
    }
    
    
    partitions <- merge(partitions, clustertable, by.x = "partition", by.y = "medoids")
    partitions$clusters
  }
  
  
  # cluster results -> draw samples in every returned cluster from the cores -> medoids -> cluster
  
  # make clusterlabels unique
  h <- 0
  for (i in 1:length(ergebnis)) {
    ergebnis[[i]] <- ergebnis[[i]] + h
    h <- max(ergebnis[[i]])
  }
  
  ergebnis <- unlist(ergebnis)
  
  # draw samples
  sample_list <- list()
  for (i in 1:max(ergebnis)) {
    which_i <- which(ergebnis == i)
    sampleindex <- sample(1:length(which_i), ifelse(length(which_i) < 2000, length(which_i), 2000), replace = FALSE)
    sample_list[[i]] <- dataset[which_i,][sampleindex,]
  }
  
  # get medoids
  medoids <- matrix(NA_integer_, nrow = length(sample_list), ncol = ncol(dataset))
  for (i in 1:length(sample_list)) {
    medoidmat <- pdist(sample_list[[i]], metric = "hamming")
    medoids[i,] <- as.integer(sample_list[[i]][which.min(rowSums2(medoidmat)), ])
  }
  
  # cluster medoids
  medoid_result <- kmeans(medoids, nclusters)$cluster
  
  
  # assign every observation in cluster the clusterlabel of the clustered medoids
  for (i in 1:nrow(medoids)) {
    ergebnis[which(ergebnis == i)] <- medoid_result[i]
  }
  
  return(ergebnis)
}