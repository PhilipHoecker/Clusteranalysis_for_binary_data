clarans_fun <- function(dataset, numlocal, ncluster, maxneighbor, ncores) {
  numlocal_list <- list()
  
  for (i in 1:numlocal) {
    # print(i)
    current <- dataset[sample(1:n, ncluster, replace = FALSE),, drop = FALSE]
    # cost_current <- sum(apply(binarydata, 1, cost_medoid, current))
    cost_current_mat <- rdist:::hamming_cdist(dataset, current)
    # cost_current <- sum(pmin.int(cost_current_mat[, 1], cost_current_mat[, 2]))
    cost_current <- sum(rowMins(cost_current_mat))
    rm(cost_current_mat)
    cost_differential <- 1
    j <- 1
    while (j < maxneighbor + ncores) {
      # print(j)
      while (cost_differential > 0) {
        cost_list <- foreach(h = 1:ncores, .combine = c, .packages = c("rdist", "matrixStats")) %dorng% {
          neighbor <- current
          neighbor[sample(1:ncluster, 1),] <- dataset[sample(1:n, 1),]
          # cost_neighbor <- sum(apply(binarydata, 1, cost_medoid, neighbor))
          cost_neighbor_mat <- rdist:::hamming_cdist(dataset, neighbor)
          # cost_neighbor <- sum(pmin.int(cost_neighbor_mat[, 1], cost_neighbor_mat[, 2]))
          cost_neighbor <- sum(rowMins(cost_neighbor_mat))
          rm(cost_neighbor_mat)
          list(neighbor, cost_neighbor)
        }
        new_best_medoid <- cost_list[[which.min(unlist(cost_list[seq(2, length(cost_list), 2)], use.names = FALSE)) * 2 - 1]]
        new_best_cost <- cost_list[[which.min(unlist(cost_list[seq(2, length(cost_list), 2)], use.names = FALSE)) * 2]]
        cost_differential <- cost_current - new_best_cost
        if (cost_differential > 0) {
          current <- new_best_medoid
          cost_current <- new_best_cost
          j <- 1
        } 
      } 
      if (cost_differential <= 0) {
        j <- j + ncores
      }
    }
    numlocal_list <- c(numlocal_list, list(cost_current, current))
  } # i
  
  medoids <- numlocal_list[[which.min(unlist(numlocal_list[seq(1, length(numlocal_list), 2)], use.names = FALSE)) * 2]]
  cost_medoids_mat <- rdist:::hamming_cdist(dataset, medoids)
  cluster_clarans <- max.col(-cost_medoids_mat, ties="first")
  
  cost_result <- numlocal_list[[which.min(unlist(numlocal_list[seq(1, length(numlocal_list), 2)], use.names = FALSE)) * 2 - 1]]
  return(list(medoids, cluster_clarans, cost_result))
}