find_central_variable3 <- function(dataset_fcv, ohne_split = NULL, ncores, ncols) {
  
  vars <- names(dataset_fcv)
  
  dataset_fcv <- as.matrix(dataset_fcv)
  colnames(dataset_fcv) <- NULL
  
  ohne <- ohne_split
  if (!is.null(ohne)) {
    ohne <- as.numeric(gsub("V", "", ohne))
    ohne_split <- as.numeric(gsub("V", "", ohne_split))
  }
  
  
  cluster_n <- nrow(dataset_fcv)
  if (is.null(ohne)) {
    nvars <- ncols
    a <- colSums2(dataset_fcv)
  } else {
    nvars <- length((1:ncols)[-ohne])
    a <- colSums2(dataset_fcv[, (1:ncols)[-ohne]])
  }
  
  I_left <- nvars * cluster_n * log(cluster_n)
  I_right <- sum(a * log(a) + (cluster_n - a) * log(cluster_n - a))
  I <- I_left - I_right
  
  if (is.null(ohne)) {
    ii <- iter(1:ncols)
  } else {
    ii <- iter((1:ncols)[-ohne])
  }
  
  
  Is <- foreach(i = ii, .combine = c, .packages = "data.table") %dopar% {
    ohne_loop <- i
    ohne <- c(ohne_split, ohne_loop)
    
    
    which_n_ix <- which(dataset_fcv[, i] == 0)
    cluster_n_Ix <- length(which_n_ix)
    
    nvars <- ncol(dataset_fcv[, -ohne])
    I_left_Ix <- nvars * cluster_n_Ix * log(cluster_n_Ix)
    
    which_cols <- (1:ncols)[-ohne]
    
    a_Ix <- colSums2(dataset_fcv[which_n_ix, which_cols])
    
    if (any(a_Ix == 0)) {
      a_Ix[which(a_Ix == 0)] <- 1
    }
    I_right_Ix <- sum(a_Ix * log(a_Ix) + (cluster_n_Ix - a_Ix) * log(cluster_n_Ix - a_Ix))
    I_Ix <- I_left_Ix - I_right_Ix
    
    
    
    
    which_n_iy <- which(dataset_fcv[, i] == 1)
    cluster_n_Iy <- length(which_n_iy)
    
    I_left_Iy <- nvars * cluster_n_Iy * log(cluster_n_Iy)
    
    a_Iy <- colSums2(dataset_fcv[which_n_iy, which_cols])
    
    if (any(a_Iy == 0)) {
      a_Iy[which(a_Iy == 0)] <- 1
    }
    I_right_Iy <- sum(a_Iy * log(a_Iy) + (cluster_n_Iy - a_Iy) * log(cluster_n_Iy - a_Iy))
    I_Iy <- I_left_Iy - I_right_Iy
    
    I - I_Ix - I_Iy
  }
  
  if (is.null(ohne)) {
    names(Is) <- (1:ncols)
  } else {
    names(Is) <- (1:ncols)[-ohne]
  }
  
  return(vars[as.numeric(names(Is)[which.max(Is)])])
}

mona <- function(dataset, n_splits, nvars, ncores, ncols) {
  ohne_split <- NULL
  central_variable <- NULL
  for (i in 1:n_splits) {
    splitvars <- append(ohne_split, central_variable)
    which_splitvars <- NULL
    x <- i
    while (floor(x / 2) >= 1) {
      which_splitvars <- append(which_splitvars, floor(x / 2))
      x <- floor(x / 2)
    }
    
    if (is.null(splitvars)) {
      ohne_split <- NULL
    } else {
      ohne_split <- splitvars[which_splitvars]
    }
    
    
    dataset_fcv <- select(dataset, starts_with("V"))
    
    vars <- names(dataset_fcv)
    
    
    ohne <- ohne_split
    if (!is.null(ohne)) {
      ohne <- as.numeric(gsub("V", "", ohne))
      ohne_split <- as.numeric(gsub("V", "", ohne_split))
    }
    
    
    cluster_n <- dataset_fcv[, .N]
    if (is.null(ohne)) {
      nvars <- ncols
      a <- dataset_fcv[, lapply(.SD, sum)]
    } else {
      nvars <- length((1:ncols)[-ohne])
      a <- dataset_fcv[, lapply(.SD, sum), .SDcols = names(dataset_fcv)[-ohne]]
    }
    
    I_left <- nvars * as.numeric(cluster_n) * log(cluster_n)
    I_right <- sum(a * log(a) + as.numeric((cluster_n - a)) * log(cluster_n - a))
    I <- I_left - I_right
    
    
    Is <- foreach(i = if (is.null(ohne)) iter(1:ncols) else iter((1:ncols)[-ohne]), .combine = c) %dopar% {
      ohne_loop <- i
      ohne1 <- c(ohne_split, ohne_loop)
      
      
      which_n_ix <- dataset_fcv[, which(get(paste0("V", i)) == 0)]
      cluster_n_Ix <- length(which_n_ix)
      
      nvars1 <- dataset_fcv[, ncol(.SD), .SDcols = names(dataset_fcv)[-ohne1]]
      I_left_Ix <- as.numeric(nvars1) * as.numeric(cluster_n_Ix) * as.numeric(log(cluster_n_Ix))
      
      which_cols <- (1:ncols)[-ohne1]
      
      a_Ix <- dataset_fcv[which_n_ix, lapply(.SD, sum), .SDcols = names(dataset_fcv)[which_cols]]
      
      if (any(a_Ix == 0)) {
        a_Ix[which(a_Ix == 0)] <- 1
      }
      I_right_Ix <- sum(a_Ix * log(a_Ix) + as.numeric((cluster_n_Ix - a_Ix)) * log(cluster_n_Ix - a_Ix))
      I_Ix <- I_left_Ix - I_right_Ix
      
      
      
      
      which_n_iy <- dataset_fcv[, which(get(paste0("V", i)) == 1)]
      cluster_n_Iy <- length(which_n_iy)
      
      I_left_Iy <- nvars1 * as.numeric(cluster_n_Iy) * log(cluster_n_Iy)
      
      a_all <- dataset_fcv[, lapply(.SD, sum), .SDcols = names(dataset_fcv)[which_cols]]
      a_Iy <- a_all - a_Ix
      
      if (any(a_Iy == 0)) {
        a_Iy[which(a_Iy == 0)] <- 1
      }
      I_right_Iy <- sum(a_Iy * as.numeric(log(a_Iy)) + as.numeric((cluster_n_Iy - a_Iy)) * log(cluster_n_Iy - a_Iy))
      I_Iy <- I_left_Iy - I_right_Iy
      
      I - I_Ix - I_Iy
    }
    
    if (is.null(ohne)) {
      names(Is) <- (1:ncols)
    } else {
      names(Is) <- (1:ncols)[-ohne]
    }
    
    central_variable <- vars[as.numeric(names(Is)[which.max(Is)])]
    
    if (i == 1) {
      dataset[, paste0("split", i) := fcase(get(central_variable) == 0, 0,
                                            get(central_variable) == 1, 1)]
    } else {
      dataset[get(paste0("split", floor(i / 2))) == ifelse(i %% 2 == 1, 1, 0), paste0("split", i) := fcase(get(central_variable) == 0, 0,
                                                                                                           get(central_variable) == 1, 1)]
    }
    
  }
  
  splits_for_clusters <- (n_splits + 1) / 2
  splits_for_clusters <- names(dataset)[(length(names(dataset)) - splits_for_clusters + 1):length(names(dataset))]
  dataset[, clusters := .GRP, by = c(splits_for_clusters)]
  return(dataset$clusters)
}