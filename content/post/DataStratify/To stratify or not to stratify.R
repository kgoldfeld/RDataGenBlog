
setDefs <- function(pX, precX, varRE, baseLO, effLOR, xLOR) {
  
  defc <- defData(varname = "p", formula = pX, variance = precX, 
                  dist = "beta", id = "site")
  defc <- defData(defc, varname = "a", formula = 0, variance = varRE)
  
  form <- genFormula(c(baseLO, effLOR, xLOR, 1), vars = c("rx", "x", "a"))
  
  defi1 <- defDataAdd(varname = "x", formula = "p", dist = "binary")
  defi2 <- defDataAdd(varname = "y", formula = form, dist = "binary", link = "logit")
  
  return(list(defc = defc, defi1 = defi1, defi2 = defi2))
  
}

genEsts <- function(strata, nclust, clustsize, pX, precX, 
                    varRE, baseLO, effLOR, xLOR) {
  
  defs <- setDefs(pX, precX, varRE, baseLO, effLOR, xLOR)
  
  dc <- genData(nclust, defs$defc)
  
  dx <- genCluster(dc, "site", clustsize , "id")
  dx <- addColumns(defs$defi1, dx)
  
  dx <- try(trtAssign(dx, strata = strata, grpName = "rx"), silent = TRUE)
  
  if ( (class(dx)[1]) == "try-error") {
    return(NULL)
  }
  
  dx <- addColumns(defs$defi2, dx)
  
  glmfit <- glm(y~factor(site) + rx + x, data = dx, family = "binomial")
  
  estrx <- t(coef(summary(glmfit))["rx", ])
  
  return(data.table(estrx))
}

forFunction <- function(strata, nclust, clustsize, pX, precX, 
                         varRE, baseLO, effLOR, xLOR) {
  
  res <- rbindlist(mclapply(1:2500, function(x) 
    genEsts(strata, nclust, clustsize, pX, precX, varRE, baseLO, effLOR, xLOR)))
  
  data.table(strata = length(strata), nclust, clustsize, pX, precX, 
             varRE, baseLO, effLOR, xLOR, 
             est = res[, mean(Estimate)],
             se.obs = res[, sd(Estimate)],
             se.est = res[, mean(`Std. Error`)],
             pval = res[, mean(`Pr(>|z|)` < 0.05)]
             )
}

####

library(parallel)

strata <- list("site", c("site", "x"))
nclust <- 8
clustsize <- c(30, 40, 50, 60, 70, 80)
pX <- 0.35
precX <- 30
varRE <- .5
baseLO <- c(-1.5, -1.25, -1.0, -0.5)
effLOR <- seq(0.5, 0.8, by = .05)
xLOR <- c(.75)

dparam <- data.table(expand.grid(strata, nclust, clustsize, pX, precX, 
              varRE, baseLO, effLOR, xLOR))

rm(strata, nclust, clustsize, pX, precX, varRE, baseLO, effLOR, xLOR)

setnames(dparam, c("strata","nclust", "clustsize", "pX", "precX", 
                   "varRE", "baseLO", "effLOR", "xLOR"))

####

resStrata <- mclapply(1:nrow(dparam), function(x) with(dparam[x,],  
    forFunction(strata[[1]], nclust, clustsize, pX, precX, varRE, baseLO, effLOR, xLOR)))

resStrata <- rbindlist(resStrata)
save(resStrata, file = "Working/resStrata.Rdata")

#### Plot

res1 <- resStrata[strata == 1, .(pval1 = pval, se.obs1 = se.obs, se.est1 = se.est,
                                 est1 = est, clustsize)]
res2 <- resStrata[strata == 2, .(pval2 = pval, se.obs2 = se.obs, se.est2 = se.est,
                                 est2 = est)]

dp <- data.table(cbind(res1, res2))

ggplot(data = dp[clustsize > 20], aes(x = pval1, y = pval2)) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_point(color = "grey20", fill = "grey20", size = .9, alpha = .5) +
  theme(panel.grid = element_blank()) +
  facet_wrap(clustsize ~ .)

ggplot(data = dp[clustsize > 20], aes(x = se.obs1, y = se.obs2)) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_point(color = "grey20", fill = "grey20", size = .9, alpha = .5) +
  theme(panel.grid = element_blank()) +
  facet_wrap(clustsize ~ .)

ggplot(data = dp[clustsize > 20], aes(x = se.est1, y = se.est2)) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_point(color = "grey20", fill = "grey20", size = .9, alpha = .5) +
  theme(panel.grid = element_blank()) +
  facet_wrap(clustsize ~ .)

ggplot(data = dp[clustsize > 20], aes(x = est1, y = est2)) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_point(color = "grey20", fill = "grey20", size = .9, alpha = .5) +
  theme(panel.grid = element_blank()) +
  facet_wrap(clustsize ~ .)



dp[, mean(pval2 > pval1), keyby = clustsize]
