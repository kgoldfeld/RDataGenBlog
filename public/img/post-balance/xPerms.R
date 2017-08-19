xPerms <- function(nsites, nRx, N = NULL) {
  
  if (is.null(N)) {
    x <- t(combn(nsites, nRx))  
  } else {
    x <- cppPerm(nsites, nRx, N)
    x <- unique(x)
  }
  
  xmat <- convert01(x, nsites)
  
  dmat <- data.table(xmat)
  
  cc <- "c(varlist)"
  y <- names(dmat)
  varlist = paste( y, collapse=",")
  cc <- gsub("varlist", varlist, cc)
  
  setattr(dmat, "varlist" , cc)
  dmat[, id := .I]
  
  return(dmat)
  
}