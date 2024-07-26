levelupdate <- function(x, rootid = "1") {
  if (x$id != rootid) {
    x$level <- x$level + 1
  }
  return(x)
}



bubble <- function(dataset, ncores, threshold, L) {
  
  nvars <- ncol(dataset)
  
  cftree <- list()
  
  
  nearest_history <- NULL
  clustroid_change <- FALSE
  split_node <- FALSE
  
  h <- 0
  id <- 0
  R <- 100
  p <- 10
  plist <- vector(mode = "list", length = L)
  
  # binaryd <- isplitRows(dataset, chunks = ncores)
  
  main_cftree <- foreach(binarydata = isplitRows(dataset, chunks = ncores), .packages = c("rlist", "rdist", "stringr", "data.table"),
                         .export = c("cftree", "clusterlabels", "clustroid_change", "split_node", "threshold",
                                     "h", "id", "R", "p", "L", "plist")) %dorng% {
                                       
                                       clusterlabels <- data.table(label = integer(nrow(binarydata)))
                                       for (i in 1:nrow(binarydata)) {
                                         # if a node was split, the observation must be inserted into the tree again
                                         run_case_again <- TRUE
                                         while (run_case_again) {
                                           run_case_again <- FALSE
                                           if (length(cftree) == 0) {
                                             h <- h + 1
                                             rowsum_clustroid <- 0
                                             id <- id + 1
                                             idparent <- "1"
                                             cftree[[h]] <- list(level = 2, leaf = TRUE, id = paste0(2, idparent, 1), L = 1, cfs = list(n = 1, clustroid = matrix(binarydata[i,], ncol = nvars),
                                                                                                                                        rowsum_clustroid = rowsum_clustroid,
                                                                                                                                        p_near = plist, rowsums_pnear = plist, 
                                                                                                                                        radius = NULL))
                                             
                                             h <- h + 1
                                             id <- id + 1
                                             
                                             cftree[[h]] <- list(level = 1, leaf = FALSE, id = idparent, L = 1, cfs = list(clustroid = t(binarydata[i,]), rowsum_clustroid = rowsum_clustroid),
                                                                 child = paste0(2, idparent, 1))
                                             idlist <- c(cftree[[2]]$id, cftree[[1]]$id)
                                             clusterlabels[i] <- as.integer(paste0(cftree[[1]]$id, 1))
                                           } else {
                                             leaf_reached <- FALSE
                                             root <- TRUE
                                             while (!leaf_reached) {
                                               if (root) {
                                                 clustroids_compare <- cftree[[list.which(cftree, id == "1")]]$cfs$clustroid
                                                 nextid <- cftree[[list.which(cftree, id == "1")]]$id
                                                 root <- FALSE
                                                 nearest_history <- NULL
                                                 parentid <- NULL
                                               } else {
                                                 clustroids_compare <- cftree[[list.which(cftree, id == nextid)]]$cfs$clustroid
                                               }
                                               
                                               
                                               leaf <- cftree[[list.which(cftree, id == nextid)]]$leaf
                                               
                                               dist_mat_clustroid <- rdist:::hamming_cdist(binarydata[i,, drop = FALSE], clustroids_compare)
                                               
                                               nearest <- which.min(dist_mat_clustroid)
                                               
                                               # save nearest in non-leasf for updating clustroids in parent-nodes
                                               nearest_history <- c(nearest_history, nearest)
                                               
                                               dist_to_clust <- min(dist_mat_clustroid)
                                               
                                               if (!leaf) {
                                                 parentid <- c(parentid, cftree[[list.which(cftree, id == nextid)]]$id)
                                                 nextid <- cftree[[list.which(cftree, id == nextid)]]$child[nearest]
                                               }
                                               
                                               # maintain leafnote
                                               if (leaf) {
                                                 leaf_reached <- TRUE
                                                 current_node <- cftree[[list.which(cftree, id == nextid)]]
                                                 
                                                 
                                                 if (dist_to_clust <= threshold) {
                                                   
                                                   # give clusterlabels
                                                   # clusterlabels[i] <- as.integer(paste0(current_node$id, nearest))
                                                   clusterlabels[i, label := as.integer(paste0(current_node$id, nearest))]
                                                   
                                                   # calculate rowsum_o
                                                   if (current_node$cfs$n[nearest] <= p) {
                                                     R_cases <- rbind(current_node$cfs$p_near[[nearest]], current_node$cfs$clustroid[nearest,])
                                                     
                                                     rowsum_o <- sum(rdist:::hamming_cdist(binarydata[i,, drop = FALSE], R_cases))
                                                     
                                                     current_node$cfs$radius[nearest] <- sum(rdist:::hamming_cdist(R_cases, current_node$cfs$clustroid[nearest,, drop = FALSE])) #/ current_node$cfs$n[nearest]
                                                     
                                                   } else {
                                                     # radius
                                                     clustroid <- current_node$cfs$clustroid[nearest,, drop = FALSE]
                                                     
                                                     
                                                     
                                                     # approximation rowsum_o
                                                     n_cl <- current_node$cfs$n[nearest]
                                                     dist_obs_clustroid <- rdist:::hamming_cdist(binarydata[i,, drop = FALSE], clustroid)
                                                     
                                                     radius <- current_node$cfs$radius[nearest] / n_cl
                                                     # zum quadrat kann weg, weil nur positiv
                                                     rowsum_o <- n_cl * radius + n_cl * dist_obs_clustroid
                                                     
                                                     current_node$cfs$radius[nearest] <- current_node$cfs$radius[nearest] + dist_obs_clustroid
                                                   }
                                                   
                                                   # update clustroid_rowsum and rowsum_pnear
                                                   # rowsum clustroid needs to be updated, as does rowsums_pnear
                                                   if (is.null(current_node$cfs$p_near[[nearest]])) {
                                                     
                                                     # update rowsum_clustroid
                                                     current_node$cfs$rowsum_clustroid[nearest] <- rdist:::hamming_cdist(current_node$cfs$clustroid[nearest,, drop = FALSE], binarydata[i,, drop = FALSE])
                                                     
                                                   } else { # p_near is not NULL
                                                     # update clustroid_rowsum and rowsums_pnear
                                                     
                                                     # update clustroid_rowsum by adding distance from clustroid to binarydata[i,] to clustroid_rowsum. Same with pfar and pnear.
                                                     current_node$cfs$rowsum_clustroid[nearest] <- current_node$cfs$rowsum_clustroid[nearest] +
                                                       rdist:::hamming_cdist(current_node$cfs$clustroid[nearest,, drop = FALSE], binarydata[i,, drop = FALSE])
                                                     
                                                     
                                                     # update rowsums_pnear
                                                     current_node$cfs$rowsums_pnear[[nearest]] <- rdist:::hamming_cdist(binarydata[i,, drop = FALSE], current_node$cfs$p_near[[nearest]]) + current_node$cfs$rowsums_pnear[[nearest]]
                                                   }
                                                   
                                                   
                                                   
                                                   
                                                   
                                                   clustroid <- current_node$cfs$clustroid[nearest,]
                                                   clustroid_rowsum <- current_node$cfs$rowsum_clustroid[nearest]
                                                   
                                                   # check if any of pnear or O have lower rowsum than clustroid
                                                   # it can happen, that the cluster consists of only the clustroid -> nearest > length(current_node$cfs$rowsums_pnear)
                                                   if (ifelse(nearest > length(current_node$cfs$rowsums_pnear), FALSE, 
                                                              ifelse(any(append(current_node$cfs$rowsums_pnear[[nearest]], rowsum_o) < clustroid_rowsum), TRUE, FALSE))) {
                                                     # indicator for clustroid change for updating clustroids of parent-nodes
                                                     clustroid_change <- TRUE
                                                     
                                                     
                                                     # old clustroid to pnearest
                                                     if (nrow(current_node$cfs$p_near[[nearest]]) < p) {
                                                       # position in pnearest to fill with old clustroid
                                                       # tofill_p <- length(current_node$cfs$p_near[[nearest]]) + 1
                                                       # current_node$cfs$p_near[[nearest]][tofill_p,] <- clustroid
                                                       current_node$cfs$p_near[[nearest]] <- rbind(current_node$cfs$p_near[[nearest]], clustroid)
                                                       # current_node$cfs$rowsums_pnear[[nearest]][tofill_p] <- clustroid_rowsum
                                                       current_node$cfs$rowsums_pnear[[nearest]] <- append(current_node$cfs$rowsums_pnear[[nearest]], clustroid_rowsum)
                                                     } else {
                                                       # if pnearest is already saturated
                                                       rowsums_pnear <- current_node$cfs$rowsums_pnear[[nearest]]
                                                       current_node$cfs$p_near[[nearest]][which.max(rowsums_pnear),] <- clustroid
                                                       current_node$cfs$rowsums_pnear[[nearest]][which.max(rowsums_pnear)] <- clustroid_rowsum
                                                     }
                                                     
                                                     # check which of pnear and O has lowest rowsum
                                                     which_new_clustroid <- which.min(append(rowsum_o, current_node$cfs$rowsums_pnear[[nearest]]))
                                                     
                                                     if (which_new_clustroid == 1) {
                                                       current_node$cfs$clustroid[nearest,] <- binarydata[i,]
                                                       current_node$cfs$rowsum_clustroid[nearest] <- rowsum_o
                                                     } else {
                                                       # -1 because the first one is rowsum_o
                                                       current_node$cfs$clustroid[nearest,] <- current_node$cfs$p_near[[nearest]][which_new_clustroid - 1,]
                                                       current_node$cfs$rowsum_clustroid[nearest] <- current_node$cfs$rowsums_pnear[[nearest]][which_new_clustroid - 1]
                                                     }
                                                     
                                                     
                                                     
                                                     # UPDATE PARENTNODE
                                                     
                                                   } else { # rowsum_o > clustroid_rowsum
                                                     
                                                     # pnear
                                                     if (is.null(current_node$cfs$p_near[[nearest]])) {
                                                       
                                                       # it can happen, that the cluster consists of only the clustroid -> nearest > length(current_node$cfs$rowsums_pnear)
                                                       if (nearest > length(current_node$cfs$rowsums_pnear)) {
                                                         current_node$cfs$rowsums_pnear <- append(current_node$cfs$rowsums_pnear, rowsum_o)
                                                       } else {
                                                         current_node$cfs$rowsums_pnear[[nearest]] <- rowsum_o
                                                       }
                                                       current_node$cfs$p_near[[nearest]] <- binarydata[i,, drop = FALSE]
                                                       
                                                     } else {
                                                       # if pnear not saturated. add binarydata[i,] to pnear
                                                       if (nrow(current_node$cfs$p_near[[nearest]]) < p) {
                                                         current_node$cfs$p_near[[nearest]] <- rbind(current_node$cfs$p_near[[nearest]], binarydata[i,])
                                                         current_node$cfs$rowsums_pnear[[nearest]] <- c(current_node$cfs$rowsums_pnear[[nearest]], rowsum_o)
                                                       }
                                                       
                                                       # if pnear already saturated, substitute pnearest case with highest rowsum with binarydata[i,]
                                                       if (rowsum_o < max(current_node$cfs$rowsums_pnear[[nearest]])) {
                                                         which_to_fill <- which.max(current_node$cfs$rowsums_pnear[[nearest]])
                                                         current_node$cfs$p_near[[nearest]][which_to_fill,] <- binarydata[i,]
                                                         current_node$cfs$rowsums_pnear[[nearest]][which_to_fill] <- rowsum_o
                                                       }
                                                     }
                                                   }
                                                   # n + 1
                                                   current_node$cfs$n[nearest] <- current_node$cfs$n[nearest] + 1
                                                 } else { # dist_to_clust > threshold
                                                   
                                                   
                                                   if (current_node$L < L) {
                                                     current_L <- current_node$L + 1
                                                     current_node$cfs$n[current_L] <- 1
                                                     current_node$cfs$clustroid <- rbind(current_node$cfs$clustroid, binarydata[i,])
                                                     current_node$cfs$rowsum_clustroid[current_L] <- 0
                                                     # current_node$cfs$R_cases[[current_L]] <- binarydata[i,]
                                                     current_node$L <- current_L
                                                     # add clustroid to clustroids parentnode
                                                     clustroid_change <- TRUE
                                                     
                                                     # give clusterlabels
                                                     # clusterlabels[i] <- as.integer(paste0(current_node$id, current_L))
                                                     clusterlabels[i, label := as.integer(paste0(current_node$id, current_L))]
                                                     
                                                   } else {
                                                     split_node <- TRUE
                                                     # indicator to run binarydata[i,] again through the tree, because of splitting.
                                                     run_case_again <- TRUE
                                                   }
                                                   
                                                 } # threshold
                                                 
                                                 # update the actual node with current_node
                                                 cftree[[list.which(cftree, id == nextid)]] <- current_node
                                                 # get current id for splitting later
                                                 current_id <- cftree[[list.which(cftree, id == nextid)]]$id
                                                 # root to true to enter root node on next iteration
                                                 root <- TRUE
                                                 
                                                 # update parentnode
                                                 if (clustroid_change) {
                                                   clustroid_change <- FALSE
                                                   # which clustroid of current_node is clustroid of parentnode
                                                   which_clustroid <- which.min(rowSums(rdist:::hamming_pdist(current_node$cfs$clustroid)))
                                                   for (pi in length(parentid):1) {
                                                     if (pi == length(nearest_history) - 1) {
                                                       # if leaf_node (length(nearest_history) - 1) change clustroid one node higher with updated clustroid
                                                       cftree[[list.which(cftree, id == parentid[pi])]]$cfs$clustroid[nearest_history[pi],] <- current_node$cfs$clustroid[which_clustroid,]
                                                     } else {
                                                       # if interior_node update clustroid of childnode with the clustroid of the clustroids of the child node
                                                       which_clustroid <- which.min(rowSums(rdist:::hamming_pdist(cftree[[list.which(cftree, id == parentid[pi + 1])]]$cfs$clustroid)))
                                                       cftree[[list.which(cftree, id == parentid[pi])]]$cfs$clustroid[nearest_history[pi],] <- cftree[[list.which(cftree, id == parentid[pi + 1])]]$cfs$clustroid[which_clustroid,]
                                                     }
                                                   }
                                                 }
                                                 
                                                 # split nodes
                                                 if (split_node) {
                                                   split_node <- FALSE
                                                   ids <- c(parentid, current_id)
                                                   for (spl in rev(ids)) {
                                                     if (cftree[[list.which(cftree, id == spl)]]$L == L) {
                                                       # find clustroids that are furthes apart
                                                       index_furthes <- which(rdist:::hamming_pdist(cftree[[list.which(cftree, id == spl)]]$cfs$clustroid) == max(
                                                         rdist:::hamming_pdist(cftree[[list.which(cftree, id == spl)]]$cfs$clustroid)), arr.ind = T)[1]
                                                       # find L/2 closest clustroid to clustroid that is furthest out
                                                       dist_clus_to_clus <- rdist:::hamming_cdist(cftree[[list.which(cftree, id == spl)]]$cfs$clustroid[index_furthes,, drop = FALSE],
                                                                                                  cftree[[list.which(cftree, id == spl)]]$cfs$clustroid)
                                                       # positions of L/2 closest
                                                       closest <- order(dist_clus_to_clus)[1:(L/2)]
                                                       # positions of L/2 farthest
                                                       farthest <- order(dist_clus_to_clus)[(L/2 + 1):length(dist_clus_to_clus)]
                                                       
                                                       # new node
                                                       # if leaf_node cfs need to be allocated
                                                       if (cftree[[list.which(cftree, id == spl)]]$leaf) {
                                                         
                                                         cfs_new_node <- list(n = cftree[[list.which(cftree, id == spl)]]$cfs$n[closest],
                                                                              clustroid = cftree[[list.which(cftree, id == spl)]]$cfs$clustroid[closest,],
                                                                              rowsum_clustroid = cftree[[list.which(cftree, id == spl)]]$cfs$rowsum_clustroid[closest],
                                                                              p_near = plist,
                                                                              rowsums_pnear = cftree[[list.which(cftree, id == spl)]]$cfs$rowsums_pnear[closest],
                                                                              radius = cftree[[list.which(cftree, id == spl)]]$cfs$radius[closest])
                                                         
                                                         cfs_new_node$p_near[1:length(closest)] <-  cftree[[list.which(cftree, id == spl)]]$cfs$p_near[closest]
                                                         
                                                         
                                                         
                                                         
                                                         cfs_old_node <- list(n = cftree[[list.which(cftree, id == spl)]]$cfs$n[farthest],
                                                                              clustroid = cftree[[list.which(cftree, id == spl)]]$cfs$clustroid[farthest,],
                                                                              rowsum_clustroid = cftree[[list.which(cftree, id == spl)]]$cfs$rowsum_clustroid[farthest],
                                                                              p_near = plist,
                                                                              rowsums_pnear = cftree[[list.which(cftree, id == spl)]]$cfs$rowsums_pnear[farthest],
                                                                              radius = cftree[[list.which(cftree, id == spl)]]$cfs$radius[farthest])
                                                         
                                                         cfs_old_node$p_near[1:length(farthest)] <-  cftree[[list.which(cftree, id == spl)]]$cfs$p_near[farthest]
                                                         
                                                         
                                                         new_node <- list(level = cftree[[list.which(cftree, id == spl)]]$level,
                                                                          leaf = cftree[[list.which(cftree, id == spl)]]$leaf,
                                                                          id = as.character(as.numeric(idlist[length(idlist)]) + 1),
                                                                          L = L / 2,
                                                                          cfs = cfs_new_node)
                                                         
                                                         old_node <- list(level = cftree[[list.which(cftree, id == spl)]]$level,
                                                                          leaf = cftree[[list.which(cftree, id == spl)]]$leaf,
                                                                          id = cftree[[list.which(cftree, id == spl)]]$id,
                                                                          L = L / 2,
                                                                          cfs = cfs_old_node)
                                                         
                                                         # children in parent_node need to be updated. gets last id + 1
                                                         cftree[[list.which(cftree, id == parentid[length(parentid)])]]$child <- append(cftree[[list.which(cftree, id == parentid[length(parentid)])]]$child,
                                                                                                                                        as.character(as.numeric(idlist[length(idlist)]) + 1),
                                                                                                                                        after = which(cftree[[list.which(cftree, id == parentid[length(parentid)])]]$child == old_node$id))
                                                         
                                                         # add new id to idlist
                                                         idlist <- append(idlist, as.character(as.numeric(idlist[length(idlist)]) + 1))
                                                         
                                                         ## update clustroid in parentnode
                                                         # old node first
                                                         # which to update
                                                         which_clustroid_to_update <- which(cftree[[list.which(cftree, id == parentid[length(parentid)])]]$child == cftree[[list.which(cftree, id == spl)]]$id)
                                                         replacement_clustroid <- which.min(rowSums(rdist:::hamming_pdist(old_node$cfs$clustroid)))
                                                         cftree[[list.which(cftree, id == parentid[length(parentid)])]]$cfs$clustroid[which_clustroid_to_update,] <- old_node$cfs$clustroid[replacement_clustroid,]
                                                         
                                                         # increase L in parentnode
                                                         cftree[[list.which(cftree, id == parentid[length(parentid)])]]$L <- cftree[[list.which(cftree, id == parentid[length(parentid)])]]$L + 1
                                                         
                                                         # new node
                                                         # clustroid of new node gets added
                                                         replacement_clustroid <- which.min(rowSums(rdist:::hamming_pdist(new_node$cfs$clustroid)))
                                                         cftree[[list.which(cftree, id == parentid[length(parentid)])]]$cfs$clustroid <- rbind(cftree[[list.which(cftree, id == parentid[length(parentid)])]]$cfs$clustroid,
                                                                                                                                               new_node$cfs$clustroid[replacement_clustroid,])
                                                         
                                                         # old node gets replaced. new node gets appended
                                                         cftree[[list.which(cftree, id == spl)]] <- old_node
                                                         cftree <- list.append(cftree, new_node)
                                                         
                                                         # update clusterlabels of observations in new node
                                                         # which_to_change <- which(clusterlabels %in% as.integer(paste0(old_node$id, closest)))
                                                         which_to_change <- which(clusterlabels$label %in% as.integer(paste0(old_node$id, closest)))
                                                         # new_labels <- as.character(clusterlabels[which_to_change])
                                                         new_labels <- as.character(clusterlabels$label[which_to_change])
                                                         # change id part of clusterlabel to new nodeid
                                                         # substr(new_labels, 1, nchar(new_node$id)) <- new_node$id
                                                         new_labels <- str_replace(new_labels, paste0("\\d{", nchar(new_node$id), "}"), new_node$id)
                                                         # change last digit to clusterposition in the new node
                                                         labelchange <- 1:length(new_labels)
                                                         for (ld in as.character(closest)) {
                                                           labelchange2 <- which(str_detect(new_labels[labelchange], paste0(new_node$id, ld, "(?!.)")))
                                                           new_labels[labelchange] <- str_replace(new_labels[labelchange], paste0(new_node$id, ld, "(?!.)"),
                                                                                                  paste0(new_node$id, as.character(which(as.character(closest) == ld))))
                                                           labelchange <- labelchange[-labelchange2]
                                                         }
                                                         
                                                         
                                                         clusterlabels[which_to_change, label := as.integer(new_labels)]
                                                         
                                                         
                                                         # update clusterlabels of observations in old node
                                                         which_to_change <- which(clusterlabels$label %in% as.integer(paste0(old_node$id, farthest)))
                                                         new_labels <- as.character(clusterlabels$label[which_to_change])
                                                         # change last digit to clusterposition in the old node
                                                         labelchange <- 1:length(new_labels)
                                                         for (ld in as.character(farthest)) {
                                                           labelchange2 <- which(str_detect(new_labels[labelchange], paste0(old_node$id, ld, "(?!.)")))
                                                           new_labels[labelchange] <- str_replace(new_labels[labelchange], paste0(old_node$id, ld, "(?!.)"),
                                                                                                  paste0(old_node$id, as.character(which(as.character(farthest) == ld))))
                                                           labelchange <- labelchange[-labelchange2]
                                                         }
                                                         
                                                         clusterlabels[which_to_change, label := as.integer(new_labels)]
                                                         
                                                         
                                                       } else { # non leaf
                                                         
                                                         # rootnode needs to be split
                                                         if (cftree[[list.which(cftree, id == spl)]]$id == parentid[1]) {
                                                           if (cftree[[list.which(cftree, id == spl)]]$L == L) {
                                                             
                                                             # only the clustroids are needed in the cfs
                                                             cfs_new_node <- list(clustroid = cftree[[list.which(cftree, id == spl)]]$cfs$clustroid[closest,])
                                                             cfs_old_node <- list(clustroid = cftree[[list.which(cftree, id == spl)]]$cfs$clustroid[farthest,])
                                                             
                                                             # contrary to leaf nodes, the children need to be split
                                                             new_node <- list(level = cftree[[list.which(cftree, id == spl)]]$level,
                                                                              leaf = cftree[[list.which(cftree, id == spl)]]$leaf,
                                                                              id = as.character(as.numeric(idlist[length(idlist)]) + 1),
                                                                              L = L / 2,
                                                                              cfs = cfs_new_node,
                                                                              child = cftree[[list.which(cftree, id == spl)]]$child[closest])
                                                             
                                                             old_node <- list(level = cftree[[list.which(cftree, id == spl)]]$level,
                                                                              leaf = cftree[[list.which(cftree, id == spl)]]$leaf,
                                                                              id = as.character(as.numeric(idlist[length(idlist)]) + 2),
                                                                              L = L / 2,
                                                                              cfs = cfs_old_node,
                                                                              child = cftree[[list.which(cftree, id == spl)]]$child[farthest])
                                                             
                                                             
                                                             
                                                             clustroid_old <- which.min(rowSums(rdist:::hamming_pdist(old_node$cfs$clustroid)))
                                                             clustroid_new <- which.min(rowSums(rdist:::hamming_pdist(new_node$cfs$clustroid)))
                                                             clustroids_new_root <- list(clustroid = rbind(new_node$cfs$clustroid[clustroid_new,], old_node$cfs$clustroid[clustroid_old,]))
                                                             
                                                             # initialize new rootnode
                                                             new_root <- list(level = 1,
                                                                              leaf = FALSE,
                                                                              id = "1",
                                                                              L = 2,
                                                                              cfs = clustroids_new_root,
                                                                              child = c(as.character(as.numeric(idlist[length(idlist)]) + 1), as.character(as.numeric(idlist[length(idlist)]) + 2)))
                                                             
                                                             # add new ids to idlist
                                                             idlist <- append(idlist, c(new_node$id, old_node$id))
                                                             
                                                             
                                                             # old node gets replaced. new node gets appended
                                                             cftree[[list.which(cftree, id == spl)]] <- old_node
                                                             cftree <- list.append(cftree, new_node)
                                                             
                                                             # add new root
                                                             cftree <- list.append(cftree, new_root)
                                                             
                                                             # levels need to be increased by 1, as well as the first character in id
                                                             cftree <- lapply(cftree, levelupdate)
                                                           }
                                                         } else { # not root or leaf
                                                           
                                                           # only the clustroids are needed in the cfs
                                                           cfs_new_node <- list(clustroid = cftree[[list.which(cftree, id == spl)]]$cfs$clustroid[closest,])
                                                           cfs_old_node <- list(clustroid = cftree[[list.which(cftree, id == spl)]]$cfs$clustroid[farthest,])
                                                           
                                                           # contrary to leaf nodes, the children need to be split
                                                           new_node <- list(level = cftree[[list.which(cftree, id == spl)]]$level,
                                                                            leaf = cftree[[list.which(cftree, id == spl)]]$leaf,
                                                                            id = as.character(as.numeric(idlist[length(idlist)]) + 1),
                                                                            L = L / 2,
                                                                            cfs = cfs_new_node,
                                                                            child = cftree[[list.which(cftree, id == spl)]]$child[closest])
                                                           
                                                           old_node <- list(level = cftree[[list.which(cftree, id == spl)]]$level,
                                                                            leaf = cftree[[list.which(cftree, id == spl)]]$leaf,
                                                                            id = cftree[[list.which(cftree, id == spl)]]$id,
                                                                            L = L / 2,
                                                                            cfs = cfs_old_node,
                                                                            child = cftree[[list.which(cftree, id == spl)]]$child[farthest])
                                                           
                                                           # children in parent_node need to be updated. gets last id + 1
                                                           cftree[[list.which(cftree, id == rev(ids)[which(rev(ids) == spl) + 1])]]$child <- append(cftree[[list.which(cftree, id == rev(ids)[which(rev(ids) == spl) + 1])]]$child,
                                                                                                                                                    as.character(as.numeric(idlist[length(idlist)]) + 1),
                                                                                                                                                    after = which(cftree[[list.which(cftree, id == rev(ids)[which(rev(ids) == spl) + 1])]]$child == old_node$id))
                                                           
                                                           # add new id to idlist
                                                           idlist <- append(idlist, as.character(as.numeric(idlist[length(idlist)]) + 1))
                                                           
                                                           # increase L in parentnode
                                                           cftree[[list.which(cftree, id == rev(ids)[which(rev(ids) == spl) + 1])]]$L <- cftree[[list.which(cftree, id == rev(ids)[which(rev(ids) == spl) + 1])]]$L + 1
                                                           
                                                           ## update clustroid in parentnode
                                                           # old node first
                                                           # which to update
                                                           which_clustroid_to_update <- which(cftree[[list.which(cftree, id == parentid[which(parentid == spl) - 1])]]$child == cftree[[list.which(cftree, id == spl)]]$id)
                                                           replacement_clustroid <- which.min(rowSums(rdist:::hamming_pdist(old_node$cfs$clustroid)))
                                                           cftree[[list.which(cftree, id == parentid[which(parentid == spl) - 1])]]$cfs$clustroid[which_clustroid_to_update,] <- old_node$cfs$clustroid[replacement_clustroid,]
                                                           
                                                           
                                                           # new node
                                                           # clustroid of new node gets added
                                                           new_clustroid <- which.min(rowSums(rdist:::hamming_pdist(new_node$cfs$clustroid)))
                                                           cftree[[list.which(cftree, id == parentid[which(parentid == spl) - 1])]]$cfs$clustroid <- rbind(cftree[[list.which(cftree, id == parentid[which(parentid == spl) - 1])]]$cfs$clustroid,
                                                                                                                                                           new_node$cfs$clustroid[new_clustroid,])
                                                           
                                                           
                                                           # old node gets replaced. new node gets appended
                                                           cftree[[list.which(cftree, id == spl)]] <- old_node
                                                           cftree <- list.append(cftree, new_node)
                                                           
                                                         }
                                                       }
                                                       
                                                     } # node$L == L
                                                     
                                                   } # spl
                                                 } # split_node
                                                 
                                               } # leaf
                                             } # while (!leaf_reached)
                                           } # else beginning
                                         } # while (run_case_again)
                                       } # i
                                       clusterlabels[, label := as.numeric(paste0(Sys.getpid(), label))]
                                       list(cftree, clusterlabels, Sys.getpid())
                                     } # foreach
  
  # extract and combine clusterlabels
  clusterlabels_combined <- rbindlist(lapply(main_cftree, function(x) x[[2]]))
  
  
  # extract and combine clustroids of leafnodes
  extractclustroids <- function(x) {
    xlist <- x[[1]][list.which(x[[1]], leaf == TRUE)]
    lapply(xlist, function(y) y$cfs$clustroid)
  }
  
  clustroids_combined <- lapply(main_cftree, extractclustroids)
  
  # give rownames to the clustroids for identification for final clustering
  
  extractids <- function(x) {
    xlist <- x[[1]][list.which(x[[1]], leaf == TRUE)]
    sapply(xlist, function(y) y$id)
  }
  
  ids <- lapply(main_cftree, extractids)
  
  # extract sessionids according to ncores for naming
  extractsessionids <- function(x) {
    x[[3]]
  }
  
  sessionids <- sapply(main_cftree, extractsessionids)
  
  # name listelements (cftrees) according to ncores
  names(clustroids_combined) <- sessionids
  
  
  nameclustroids <- function(x, id) {
    xlist <- x[1]
    for (i in 1:length(xlist[[1]])) {
      rownames(xlist[[1]][[i]]) <- paste0(names(xlist), id[i], seq(1, nrow(xlist[[1]][[i]])))
    }
    return(xlist)
  }
  
  
  for (i in 1:length(clustroids_combined)) {
    clustroids_combined[[i]] <- nameclustroids(clustroids_combined[i], ids[[i]])
  }
  
  clustroids_final <- NULL
  for (i in 1:length(clustroids_combined)) {
    for (j in 1:length(clustroids_combined[[i]][[1]])) {
      clustroids_final <- rbind(clustroids_final, clustroids_combined[[i]][[1]][[j]])
    }
  }
  
  return(list(clustroids = clustroids_final, clusterlabels = clusterlabels_combined))
}

# cross validation for threshold

cross_bubble <- function(dataset, ncluster, orig) {
  
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
  
  L_list <- c(10, 20, 50, 100)
  thresholdlist <- seq(0.5, 0.8, 0.1)
  flist <- NULL
  h <- 0
  for (thresh in 1:length(thresholdlist)) {
    for (L in 1:length(L_list)) {
      h <- h + 1
      bubble_result <- bubble(datasample, 1, thresholdlist[thresh], L_list[L])
      
      if (nrow(bubble_result$clustroids) > ncluster) {
        clustroids_dist <- rdist(bubble_result$clustroids, metric = "hamming")
        # result_pam <- Cluster_Medoids(bubble_result$clustroids, clusters = ncluster, distance_metric = "hamming", threads = ncores)
        result_pam <- pam(clustroids_dist, k = ncluster, cluster.only = TRUE)
        for (i in 1:ncluster) {
          set(bubble_result$clusterlabels, which(bubble_result$clusterlabels$label %in% rownames(bubble_result$clustroids[which(result_pam == i),])), j = 1L, value = i)
        }
      } else {
        for (i in 1:ncluster) {
          bubble_result$clusterlabels[label == as.integer(rownames(bubble_result$clustroids))[i], label := i]
        }
      }
      
      flist[h] <- rec_pre_f1(orig_sample, bubble_result$clusterlabels$label, n)$f1
      names(flist)[h] <- paste(thresholdlist[thresh], L_list[L])
    }
    
  }
  
  # threshold <- thresholdlist[which.max(flist)]
  threshold <- as.numeric(unlist(str_split(names(flist[which.max(flist)]), " "))[1])
  L <- as.numeric(unlist(str_split(names(flist[which.max(flist)]), " "))[2])
  return(data.table(threshold = threshold, L = L))
}