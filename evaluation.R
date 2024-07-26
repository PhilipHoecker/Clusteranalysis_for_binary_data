rec_pre_f1 <- function(orig, result, n) {
  
  n2 <- (n * (n - 1)) / 2
  m <- cbind(orig, result)
  rt <- bigtabulate(m, c(1,2))
  # rt <- table(orig, result)
  
  tp <- sum(apply(rt, 1, choose, 2))
  fptp <- sum(sapply(colSums(rt), choose, 2))
  fp <- fptp - tp
  
  
  fn <- sum(apply(rt, 1, function(x) {
    xx <- x
    sumx <- numeric(sum(xx > 1))
    i <- 0
    while(max(xx) > 1) {
      i <- i + 1
      wm <- which.max(xx)
      sumx[i] <- as.numeric(xx[wm]) * as.numeric(sum(xx[-wm]))
      xx <- xx[-wm]
    }
    return(sum(sumx))
    
  }))
  
  recall <- tp / (tp + fn)
  precision <- tp / (tp + fp)
  f1 <- 2 * ((precision * recall) / (precision + recall))
  
  return(data.frame(recall = recall, precision = precision, f1 = f1))
}


Entropy <- function(orig, result, n) {
  result <- frank(result, ties.method = "dense")
  maxresult <- max(result)
  maxorig <- max(orig)
  
  hw <- numeric(maxresult)
  for (i in 1:maxresult) {
    hw_temp <- numeric(maxorig)
    for (j in 1:maxorig) {
      which_i <- which(result == i)
      wc_nw <- length(which(orig[which_i] == j)) / length(which_i)
      if (wc_nw != 0) {
        hw_temp[j] <- wc_nw * log2(wc_nw)
      }
    }
    
    hw[i] <- -sum(hw_temp)
  }
  
  h_omega <- numeric(maxresult)
  for (i in 1:maxresult) {
    h_omega[i] <- hw[i] * (length(result[which(result == i)]) / n)
  }
  sum(h_omega)
}


purity <- function(orig, cluster_results, n, ncluster) {
  csize <- length(orig) / ncluster
  steps <- c(0, seq(csize, n, csize))
  
  clusterlabels <- unique(cluster_results)
  
  Si <- numeric(ncluster)
  for (i in 1:ncluster) {
    maxnij <- max(table(cluster_results[(steps[i] + 1):steps[i + 1]]))
    csize_i <- sum(cluster_results == clusterlabels[i], na.rm = TRUE)
    Si[i] <- (maxnij / csize_i) * (csize_i / n)
  }
  
  sum(Si, na.rm = TRUE)
}

rand <- function(pred, orig, n) {
  m <- cbind(pred, orig)
  tab <- bigtabulate(m, c(1,2))
  ra <- 1 + (sum(tab^2) - 0.5 * (sum(rowSums(tab)^2) + sum(colSums(tab)^2))) / choose(n, 2)
  return(ra)
}