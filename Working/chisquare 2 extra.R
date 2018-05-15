### Addendum

```{r}
drawUrn <- function(urn, n) {
  
  dtUrn <- data.table(id = 1:length(urn), urn)
  ids <- dtUrn[, sample(id, n, replace = FALSE)]
  selected <- dtUrn[ids, urn]
  urnafter <- dtUrn[-ids, urn]
  
  return(list(selected = selected, urnafter = urnafter))
  
}

genCtable <- function(row.margins, col.margins) {
  
  nR <- length(row.margins)
  nC <- length(col.margins)
  
  urn <- rep(c(1:nC), times = col.margins)
  
  d <- list()
  
  for (i in 1:nR) {
    
    draw <- drawUrn(urn, row.margins[i])
    d[[i]] <- draw$selected
    urn <- draw$urnafter
    
  }
  
  return(matrix(unlist(lapply(d, function(x) tabulate(x, nbins = nC))), nR, nC, byrow = T))
  
}

xm <- genCtable(row, col)
addmargins(xm)

```


```{r, eval=FALSE, echo=FALSE}

# Conditioning on column totals

prob <- row/N

condCol <- lapply(seq_len(length(col)), 
                  function(i) t(rmultinom(10000, size = col[i], prob=prob) ))

condCm <- lapply(seq_len(10000), 
                 
                 function(i) {
                   t(do.call(rbind, lapply(condCol, function(x) x[i,])))
                 }
)

addmargins(condCm[[1]])
addmargins(condCm[[2]])

sumC <- avgMatrix(condCm, sLabel = "C")

X2 <- sapply(condCm, function(x) estX2(x, expected))


round(sumC$sampAvg, 0)
round(sumC$sampVar, 0)

trueChisq <- rchisq(10000, 9)
# Comparing means
round(c( mean(X2), mean(trueChisq)), 1)

# Comparing variance
round(c( var(X2), var(trueChisq)), 1)
```