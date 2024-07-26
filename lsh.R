rec_pre_f1_lsh <- function(result, ncluster, n, orig) {
  
  # tp <- lsh_result[(id1 < 501 & id2 < 501) | (id1 > 500 & id2 > 500), .N]
  
  tp_lsh <- function(result, ncluster, n, orig) {
    # bound <- round(n / ncluster)
    # bounds <- seq(0, n, bound)
    tp <- integer(ncluster)
    for (i in 1:ncluster) {
      # tp[i] <-  result[(id1 %in% (bounds[i] + 1):bounds[i + 1] & id2 %in% (bounds[i] + 1):bounds[i + 1]), .N]
      whichlabel <- which(orig == i)
      tp[i] <-  result[(id1 %in% whichlabel & id2 %in% whichlabel), .N]
    }
    sum(tp)
  } 
  
  tp <- tp_lsh(result, ncluster, n, orig)
  
  # fp <- lsh_result[(id1 < 501 & id2 > 500) | (id1 > 500 & id2 < 501), .N]
  
  fp_lsh <- function(result, ncluster, n, orig) {
    # bound <- round(n / ncluster)
    # bounds <- seq(0, n, bound)
    fp <- integer(ncluster)
    for (i in 1:ncluster) {
      whichlabel <- which(orig == i)
      # fp[i] <-  result[(id1 %in% (bounds[i] + 1):bounds[i + 1] & !id2 %in% (bounds[i] + 1):bounds[i + 1]), .N]
      fp[i] <-  result[(id1 %in% whichlabel & !id2 %in% whichlabel), .N]
    }
    sum(fp)
  }
  
  fp <- fp_lsh(result, ncluster, n, orig)
  
  # cluster coef brauch ich eigentlich nicht wenn ich precision habe
  # cluster_coef <- tp / fp
  
  
  # tn <- alle moeglichen paare aus unterschiedlichen clustern - fp
  # clustergroeße^anzahl cluster
  # tn <- ((binarydata[, .N] / ncluster)^ncluster) - fp
  tn <- ((n / ncluster)^ncluster) - fp
  
  # fn <- alle möglichen Paare aus selben cluster - tn
  # alle möglichen paare aus einem cluster -> choose(ncluster, 2) * anzahl cluster
  
  # fn <- (choose(binarydata[, .N] / ncluster, 2) * ncluster) - tp
  fn <- (choose(n / ncluster, 2) * ncluster) - tp
  
  
  
  
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  f1 <- 2 * ((precision * recall) / (precision + recall))
  # rand index
  # rand <- (tp + tn) / choose(binarydata[, .N], 2)
  rand <- (tp + tn) / choose(n, 2)
  
  return(data.frame(recall = recall, precision = precision, f1 = f1, rand = rand))
}

pairsfun <- function(d) {
  a <- comboGeneral(d, 2)
  setDT(list(id1 = a[, 1], id2 = a[, 2]))
}

# lsh
lsh_fun <- function(dataset, n_hash, b, ncores) {
  
  
  # minhash
  
  n <- nrow(dataset)
  ncols <- ncol(dataset)
  n_hash <- n_hash
  
  
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
  
  # the minimum of the hashvalues in the hashmatcolumns form the ultimate hashvalue for hashfunction r -> sigmat 
  ii <- isplitVector(whichhashmat, chunks = ncores)
  
  sigmat <- foreach(i = ii, .combine = cbind, .packages = "matrixStats" ) %dopar% {
    sapply(i, function(x) colMins(hashmat[x,, drop = FALSE]))
  }
  
  ii <- isplitRows(sigmat, chunks = b)
  
  buckets <- foreach(dat = ii) %dopar% {
    dat2 <- data.table(h = apply(dat, 2, fastdigest), id = 1:ncol(dat))
    dat2
  }
  
  return(buckets)
}


# cross validation 

cross_lsh <- function(dataset, ncluster, orig, ncores) {
  
  n <- nrow(dataset)
  # bound <- round(n / ncluster)
  # bounds <- seq(0, n, bound)
  
  datasample <- matrix(nrow = 2000, ncol = ncol(dataset))
  boundsample <- round(2000 / ncluster)
  # boundssample <- seq(0, 2000, boundsample)
  boundssample <- round(seq(0, 2000, length.out = ncluster + 1))
  for (i in 1:ncluster) {
    # datasample[(boundssample[i] + 1):boundssample[i + 1],] <- dataset[which(orig == i),][sample(1:round(n / ncluster), round(2000 / ncluster)),]
    datasample[(boundssample[i] + 1):boundssample[i + 1],] <- dataset[which(orig == i),][sample(1:length(which(orig == i)), length((boundssample[i] + 1):boundssample[i + 1])),]
  }
  
  shuf_index <- sample(1:nrow(datasample), nrow(datasample))
  datasample <- datasample[shuf_index,]
  
  orig_sample <- integer()
  for (i in 1:ncluster) {
    orig_sample <- c(orig_sample, rep(i, round(nrow(datasample) / ncluster)))
  }
  
  orig_sample <- orig_sample[shuf_index]
  
  # bandslist <- 2:15
  # rowslist <- 2:15
  # bestlist <- NULL
  # h <- 0
  # for (i in bandslist) {
  #   print(i)
  #   for (j in rowslist) {
  #     h <- h + 1
  #     lsh_result <- get_similar_pairs(datasample, bands_number = i, rows_per_band = j, distance = "jaccard")
  #     bestlist[h] <- rec_pre_f1_lsh(lsh_result, ncluster, nrow(datasample), orig_sample)$f1
  #     names(bestlist)[h] <- paste0(i, j)
  #   }
  # }
  
  bandslist <- 2:15
  n_hash_list <- seq(4, 50, 2)
  bestlist <- NULL
  h <- 0
  for (i in bandslist) {
    for (j in n_hash_list) {
      if ((j / i) %% 2 == 0) {
        # print(h)
        h <- h + 1
        lsh_result <- lsh_fun(datasample, n_hash = j, b = i, ncores = ncores)
        
        lsh_result <- lapply(lsh_result, function(x) {
          whichdup <- which(duplicated(x$h) | duplicated(x$h, fromLast=TRUE))
          if (length(whichdup) > 0) {
            x[whichdup, id, by = h][, pairsfun(id), by = h][, .(id1, id2)]
          } else {
            data.table()
          }
        })
        
        lsh_result <- rbindlist(lsh_result)
        
        # deduplicate in case a pair is in more than one bucket
        lsh_result <- unique(lsh_result)
        
        if (lsh_result[,.N] > 0) {
          bestlist[h] <- rec_pre_f1_lsh(lsh_result, ncluster, nrow(datasample), orig_sample)$f1
          names(bestlist)[h] <- paste0(i, " ", j)
        }
        
      }
      
    }
  }
  
  # whichbest <- names(bestlist[which.max(bestlist)])
  whichbest <- as.integer(strsplit(names(bestlist[which.max(bestlist)]), " ")[[1]])
  # bands_number <- as.integer(substr(whichbest, 1, 1))
  bands_number <- whichbest[1]
  # n_hash <- as.integer(substr(whichbest, 2, 2))
  n_hash <- whichbest[2]
  return(data.frame(bands_number = bands_number, n_hash = n_hash))
}


lsh_eval_fun <- function(buckets, n_cluster, n_data_sample, orig) {
  
  # clusterseq <- seq(0, nrow(dataset), length.out = n_cluster + 1)
  clusterseq <- round(seq(0, buckets[[1]][, .N], length.out = n_cluster + 1))
  clusterseq[length(clusterseq)] <- buckets[[1]][, .N]
  
  data_sample <- NULL
  for (i in 1:n_cluster) {
    data_sample <- c(data_sample, sample((clusterseq[i] + 1):clusterseq[i + 1], n_data_sample / n_cluster, replace = FALSE))
  }
  
  # find sampled observations in the buckets
  buckets2 <- lapply(buckets, function(x, data_sample) x[id %in% data_sample], data_sample)
  
  # form candidate pairs in the buckets of the sampled observations
  buckets2 <- lapply(buckets2, function(x) {
    whichdup <- which(duplicated(x$h) | duplicated(x$h, fromLast=TRUE))
    if (length(whichdup) > 0) {
      x[whichdup, id, by = h][, pairsfun(id), by = h][, .(id1, id2)]
    } else {
      data.table()
    }
  })
  
  buckets2 <- rbindlist(buckets2)
  
  # deduplicate in case a pair is in more than one bucket
  buckets2 <- unique(buckets2)
  
  rec_pre_f1_lsh(buckets2, n_cluster, n_data_sample, orig)
  
}