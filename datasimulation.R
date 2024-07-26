corrbinary <- function(n, means, sigma, sigma_correction = TRUE, hashs = NULL, cores) {
  
  ms <- means
  nvars <- length(means)
  
  nearPSD <- function(c) {
    
    n <- dim(c)[1]
    e <- eigen(c, sym = TRUE)
    val <- e$values * (e$values > 0)
    vec <- e$vectors
    T <- sqrt(diag(as.vector(1/(vec^2 %*% val)), n, n))
    B <- T %*% vec %*% diag(as.vector(sqrt(val)), n, n)
    out <- B %*% t(B)
    
    return(out)
  }
  
  
  sigma.star_fun <- function(ms, sigma, cores) {
    no.bin <- length(ms)
    p <- ms
    
    
    # cl <- makeCluster(cores)
    # registerDoSNOW(cl)
    
    result <- foreach(i = 1:no.bin, .packages = c("psych", "foreach"), .combine = rbind) %dopar%  {
      foreach(j = 1:no.bin, .combine = c) %do% {
        suppressMessages(phi2tetra(sigma[i, j], c(p[i], p[j])))
      }
    }
    
    # stopCluster(cl)
    
    if (!all(diag(result) == 1)) {
      diag(result) <- 1
    }
    
    if (!is.symmetric.matrix(result)) {
      result <- as.matrix(forceSymmetric(result, uplo = "L"))
    }
    
    
    if (!is.positive.semi.definite(result)) {
      result <- nearPSD(result)
    }
    
    result <- (result + t(result)) / 2
    return(result)
  }
  
  
  
  if (sigma_correction) {
    
    # Korrelationen die die Grenzen ueber- bzw. unterschreiten,
    # werden auf den maximal bzw. minimal zulaessigen Wert gesetzt.
    
    # setup zulaessige Korrelationen
    lows <- matrix(nrow = nvars, ncol = nvars)
    highs <- matrix(nrow = nvars, ncol = nvars)
    p <- ms
    q <- 1 - p
    for (i in 1:length(p)) {
      for (j in 1:length(p)) {
        if (i == j) {
          lows[i, j] <- 1L
          highs[i, j] <- 1L
        } else {
          lows[i, j] <- max(-sqrt((p[i] * p[j])/(q[i] * q[j])), -sqrt((q[i] * q[j])/(p[i] * p[j])))
          highs[i, j] <- min(sqrt((p[i] * q[j])/(q[i] * p[j])), sqrt((q[i] * p[j])/(p[i] * q[j])))
          sigma[i, j] <- ifelse(sigma[i, j] < lows[i, j], lows[i, j], ifelse(sigma[i, j] > highs[i, j], highs[i, j], sigma[i, j]))
        }
      }
    }
  } else {
    toolowstop <- 0
    toohighstop <- 0
    toolow <- data.table()
    toohigh <- data.table()
    
    # setup zulaessige Korrelationen
    lows <- matrix(nrow = nvars, ncol = nvars)
    highs <- matrix(nrow = nvars, ncol = nvars)
    p <- ms
    q <- 1 - p
    for (i in 1:length(p)) {
      for (j in 1:length(p)) {
        if (i == j) {
          lows[i, j] <- 1L
          highs[i, j] <- 1L
        } else {
          lows[i, j] <- max(-sqrt((p[i] * p[j])/(q[i] * q[j])), -sqrt((q[i] * q[j])/(p[i] * p[j])))
          highs[i, j] <- min(sqrt((p[i] * q[j])/(q[i] * p[j])), sqrt((q[i] * p[j])/(p[i] * q[j])))
          if (sigma[i, j] < lows[i, j]) {
            toolow <- rbindlist(list(toolow, data.table(row = i, column = j, lower.bound = lows[i, j])))
            toolowstop <- 1L
          }
          if (sigma[i, j] > highs[i, j])  {
            toohigh <- rbindlist(list(toohigh, data.table(row = i, column = j, upper.bound = highs[i, j])))
            toohighstop <- 1L
          }
        }
      }
    }
    if (toolowstop == 1 & toohighstop == 0) {
      setDF(toolow)
      print(toolow)
      stop("There were correlations outside the acceptable range. Check output for positions and bounds that were exceeded.")
    }
    if (toohighstop == 1 & toolowstop == 0) {
      setDF(toohigh)
      print(toohigh)
      stop("There were correlations outside the acceptable range. Check output for positions and bounds that were exceeded.")
    }
    if (toohighstop == 1 & toolowstop == 1) {
      setDF(toohigh)
      setDF(toolow)
      print(toohigh)
      print(toolow)
      stop("There were correlations outside the acceptable range. Check output for positions and bounds that were exceeded.")
    }
  }
  
  # sollte die Korrelationsmatrix nicht positiv semi-definite sein,
  # wird sie dazu gemacht.
  
  if (!is.positive.semi.definite(sigma)) {
    sigma <- nearPSD(sigma)
  }
  
  
  sigma_star <- sigma.star_fun(ms, sigma, cores)
  
  
  if (!is.positive.semi.definite(sigma_star)) {
    sigma_star <- nearPSD(sigma_star)
  }
  
  if (!is.positive.definite(sigma_star)) {
    sigma_star <- as.matrix(nearPD(sigma_star, corr = TRUE, keepDiag = TRUE)$mat)
  }
  
  
  data <- as.data.table(rmvn(n, mu = rep(0, nvars), sigma = sigma_star, ncores = detectCores() - 10))
  
  
  for (i in 1:nvars) {
    data[get(paste0("V", i)) < qnorm(1 - p[i]), c(paste0("V", i)) := 0L]
    data[get(paste0("V", i)) > qnorm(1 - p[i]), c(paste0("V", i)) := 1L]
  }
  
  
  for (j in 1:nvars) {
    set(data, j = j, value = as.integer(data[[j]]))
  }
  
  
  if (!is.null(hashs)) {
    hashs2 <- data[, fastdigest(.SD), by = seq_len(n)]$V1
    substitute <- which(hashs2 %in% hashs)
    while (length(substitute > 0)) {
      data_substitute <- as.data.table(rmvn(length(substitute), mu = rep(0, nvars), sigma = sigma_star, ncores = detectCores() - 10))
      
      for (i in 1:nvars) {
        data_substitute[get(paste0("V", i)) < qnorm(1 - p[i]), c(paste0("V", i)) := 0L]
        data_substitute[get(paste0("V", i)) > qnorm(1 - p[i]), c(paste0("V", i)) := 1L]
      }
      
      for (j in 1:nvars) {
        set(data_substitute, j = j, value = as.integer(data_substitute[[j]]))
      }
      
      # data[substitute] <- data_substitute
      data[substitute, names(data) := data_substitute]
      hashs2[substitute] <- data[substitute, fastdigest(.SD), by = seq_len(data[substitute, .N])]$V1
      substitute <- which(hashs2 %in% hashs)
    }
  }
  
  return(data)
  
}

simcor <- function(k = 6, size = c(10,5,8,2,15,50), rho = c(0.7,0.7,0.5,0.9,0.85,0.4),
                   delta = 0.39, epsilon = 0.99 - max(rho), eidim = 2){
  ndim <- sum(size)
  bigcor <- matrix(rep(delta, ndim * ndim), ncol = ndim)
  for (i in 1:k) {
    cor <- matrix(rep(rho[i], size[i] * size[i]), ncol = size[i])
    if (i == 1) {bigcor[1:size[1], 1:size[1]] <- cor}
    if (i != 1) {bigcor[(sum(size[1:(i - 1)]) + 1):sum(size[1:i]),
                        (sum(size[1:(i - 1)]) + 1):sum(size[1:i])] <- cor}
  }
  diag(bigcor) <- 1 - epsilon
  eivect <- c()
  for (i in 1:ndim) {
    ei <- runif(eidim, -1, 1)
    eivect <- cbind(eivect, sqrt(epsilon) * ei/sqrt(sum(ei^2)))
  }
  bigE <- t(eivect) %*% eivect
  cor.nz <- bigcor + bigE
  cor.nz
}