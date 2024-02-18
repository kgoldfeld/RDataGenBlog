
replication <- function(n_arms, n_eds) {
  
  randomize <- function(n_per_arm, n_arms, dd.hs) {
    
    dd.hs <- dd.hs[order(-N)]
    n_eds <- copy(dd.hs$N)
    
    dd.hs[, ed_index := .I ]
    setkey(dd.hs, hs_index)
    
    A <- 1:n_arms
    mat.assign <- matrix(0, nrow = length(n_eds), ncol = n_arms)
    
    for (i in seq_along(n_eds)) {
      
      # keeping only arms that haven't reached threshold
      
      B <- A[apply(mat.assign, 2, sum) < n_per_arm] 
      
      if (length(B) == 1) {       # only one possible arm remains
        
        mat.assign[i, B] <- 1
        
      } else {                    # multiple arms remain
        
        if (length(B) < n_eds[i]) {
          return(NULL) # added to prevent sampling error
        }
        a <- sample(B, n_eds[i], replace = FALSE)
        mat.assign[i, a] <- 1
      }
    }
    
    return(mat.assign[dd.hs$ed_index,])
  }
    
  dd.hs <- data.table(hsid = 1:length(n_eds), N = n_eds )
  dd.ed <- genCluster(dd.hs, "hsid", "N", "id")
  
  n_groups <- ceiling(n_eds / n_arms)
  
  dd.split <- data.table()
  
  for (i in seq_along(n_groups)) {
    eds <- dd.ed[hsid == i, id]
    n_to_samp <- min(8, dd.ed[hsid == i, .N])
    for (j in 1:n_groups[i]) {
      
      if (length(eds) > 1) sample_eds <- sample(eds, n_to_samp)
      else sample_eds <- eds
      
      dd.split <- rbind(dd.split, data.table(hs = i, group = j, ed = sample_eds) ) 
      eds <- eds[!(eds %in% sample_eds)]
      n_to_samp <- min(8, length(eds))
    }
  }
  
  dd.split[, ed_index := 1:.N, keyby = .(hs, group)]
  dd.hs <- dd.split[, .N, keyby = .(hs, group)][, hs_index := .I][]
  dd.split <- merge(dd.split, dd.hs, by = c("hs", "group"))
  
  n_per_arm <- dd.hs[, sum(N)] / n_arms
  
  mat.assign <- randomize(n_per_arm, n_arms, dd.hs)
  if (is.null(mat.assign)) return(data.table())
  
  rownames(mat.assign) <- 1:nrow(mat.assign)
  colnames(mat.assign) <- 1:ncol(mat.assign)
  
  dd.assign <- data.table(as.table(mat.assign)) # vectorize
  dd.assign <- dd.assign[N == 1]
  
  dd.assign <- dd.assign[, .(hs_index = V1, arm = V2)]
  dd.assign[,`:=`(hs_index = as.numeric(hs_index), arm = as.numeric(arm))]
  
  setkey(dd.assign, hs_index)
  dd.assign[, ed_index := 1:.N, keyby = hs_index]
  
  dd.assign <- merge(dd.assign, dd.split, by = c("hs_index", "ed_index"))
  
  dd.assign <- dd.assign[, .(hs, group, ed, arm)]
  setkey(dd.assign, ed)
  
  return(dd.assign)
}

n_arms <- 8
# n_eds = c(4, 7, 5, 2, 9, 6, 3, 3, 1) 
n_eds = c(12, 1, 11, 5, 4, 7) 

replication(n_arms, n_eds)

res <- parallel::mclapply(1:56000, function(x) replication(n_arms, n_eds))
resOK <- res[unlist(lapply(res, function(x) length(x) > 0))]

resbysite <- parallel::mclapply(1:sum(n_eds), function(x) unlist(lapply(resOK, function(a) a[ed == x, arm])))
probs <- lapply(resbysite, function(x) round(prop.table(table(x)), 3))
probs <- matrix(unlist(probs), nrow = length(probs), byrow = TRUE)

getrate <- function(h) {
  
  getdist <- function(x, h) {
    tabulate(x[hs == h][, arm], nbins = 8)
  }
  
  y <- lapply(resOK, function(x) getdist(x, h))
  round(apply(do.call(rbind, y), 2, mean), 3)
}

rates <- do.call(rbind, lapply(1:length(n_eds), function(x) getrate(x)))
save(probs, rates, file = "data/probs.Rdata")
