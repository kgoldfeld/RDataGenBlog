library(lme4)
library(parallel)

Covar <- function(dx, clust, period1, period2, x_0, x_1) {
  
  v0 <- dx[ctemp == clust & period == period1, Y - x_0]
  v1 <- dx[ctemp == clust & period == period2, Y - x_1]
  sum(v0 %*% t(v1))
  
}

calcBP <- function(dx, period1, period2) {
  
  # dx <- copy(d2)
  
  # create cluster numbers starting from 1
  
  tt <- dx[, .N, keyby = cluster]
  nclust <- nrow(tt)
  dx[, ctemp := rep(1:nclust, times = tt$N)]
  
  dx <- dx[period %in% c(period1, period2)]
  
  ## Grand means
  
  dg <- dx[, .(m=.N, mu = mean(Y)), keyby = .(ctemp, period)]
  dg <- dcast(dg, formula = ctemp ~ period, value.var = c("m","mu"))
  
  setnames(dg, c("ctemp", "m_0", "m_1", "mu_0", "mu_1"))
  
  x_0 <- dg[, sum(m_0 * m_1 * mu_0)/sum(m_0 * m_1)]
  x_1 <- dg[, sum(m_0 * m_1 * mu_1)/sum(m_0 * m_1)]
  
  ## Variance (denominator)
  
  dss_0 <- dx[period == period1, .(ss_0 = sum((Y - x_0)^2)), 
              keyby = ctemp]
  dss_0[, m_1 := dg[, m_1]]
  v_0 <- dss_0[, sum(m_1 * ss_0)]
  
  dss_1 <- dx[period == period2, .(ss_1 = sum((Y - x_1)^2)), 
              keyby = ctemp]
  dss_1[, m_0 := dg[, m_0]]
  v_1 <- dss_1[, sum(m_0 * ss_1)]
  
  ## Covariance
  
  v0v1 <- sapply(1:nclust, 
                 function(x) Covar(dx, x, period1, period2, x_0, x_1))
  
  bp.icc <- sum(v0v1)/sqrt(v_0 * v_1)
  bp.icc
  
}

btwnPerICC <- function(dd, period1, period2, byWave = FALSE) {
  
  if (byWave) {
    waves <- dd[, unique(startTrt)]
    bpICCs <- sapply(waves, function(x)  
      calcBP(dd[startTrt==x], period1, period2))
    return(mean(bpICCs))
  } else {
    calcBP(dd, period1, period2)
  }
}

withinPerICC <- function(dx) {
  
  lmerfit <- lmer(Y~rx + (1|cluster), data = dx)
  vars <- as.data.table(VarCorr(lmerfit))[, vcov]
  vars[1]/sum(vars)
  
}

genPairs <- function(n) {
  x <- combn(x = c(1:n-1), 2)
  lapply(seq_len(ncol(x)), function(i) x[,i])
}

iccs <- function(dd, byWave = FALSE) {
  
  system("echo x")
  
  nperiods <- dd[, length(unique(period))]
  
  bperiods <- genPairs(nperiods)
  names <- 
    unlist(lapply(bperiods, function(x) paste0("bp", x[1], x[2])))
  
  bp.icc <- sapply(bperiods, 
                   function(x) btwnPerICC(dd, x[1], x[2], byWave))
  
  bdd.per <- lapply(1:nperiods - 1, function(x) dd[period == x])
  
  wp.icc <- lapply(bdd.per, 
                   function(x) withinPerICC(x))
  wp.icc <- unlist(wp.icc)
  nameswp <- sapply(1:nperiods - 1, function(x) paste0("wp", x))
  
  do <- data.table(t(c(bp.icc, wp.icc)))
  setnames(do, c(names, nameswp))
  
  return(do[])
  
}

## Constant

defc <- defData(varname = "ceffect", formula = 0, variance = 0.15, 
                dist = "normal", id = "cluster")
defc <- defData(defc, "m", formula = 10, dist = "nonrandom")

defa <- defDataAdd(varname = "Y", 
                   formula = "0 + 0.10  * period + 1 * rx +  ceffect", 
                   variance = 2, dist = "normal")

genDD <- function(defc, defa, nclust, nperiods, waves, len, start) {
  
  dc <- genData(nclust, defc)
  dp <- addPeriods(dc, nperiods, "cluster")
  dp <- trtStepWedge(dp, "cluster", nWaves = waves, 
                     lenWaves = len, startPer = start)
  dd <- genCluster(dp, cLevelVar = "timeID", "m", "id")
  dd <- addColumns(defa, dd)
  return(dd[])
  
}

icc <- mclapply(1:200, 
                function(x) iccs(genDD(defc, defa, 100, 7, 4, 1, 2), byWave = T),
                mc.cores = 4
)

observed <- sapply(rbindlist(icc), function(x) mean(x))
expected <- 0.15/(0.15 + 2)

pairs1 <- data.table(matrix(unlist(genPairs(7)), ncol = 2, byrow = TRUE))
pairs2 <- data.table(matrix(rep((1:7) - 1, each = 2), 
                            ncol = 2, byrow = TRUE))

results <- cbind(rbind(pairs1, pairs2), exp=expected, obs=observed)
setnames(results, c("V1", "V2"), c("Per1", "Per2"))

pres1 <- melt(results, id.vars = c("Per1", "Per2"))
pres1[Per1 == Per2, type := "within"]
pres1[Per1 != Per2, type := "between"]
pres1[variable == "exp", type := "expected"]

p1 <- ggplot(data = pres1, aes(x = variable, y = value, fill = type)) +
  geom_bar(stat = "identity") +
  facet_grid(Per1 ~ Per2) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank()) +
  scale_y_continuous(limits = c(0, 0.21), breaks = c(.05, .15)) +
  scale_fill_manual(values = c("#82a473", "grey80", "#9573a4"), 
                    breaks = c("between", "within"))

## Exchangeable

defc <- defData(varname = "mu", formula = 0, 
                dist = "nonrandom", id = "cluster")
defc <- defData(defc, "s2", formula = 0.15, dist = "nonrandom")

defa <- defDataAdd(varname = "Y", 
                   formula = "0 + 0.10  * period + 1 * rx + cteffect", 
                   variance = 2, dist = "normal")

genDD <- function(defc, defa, nclust, nperiods, waves, len, start) {
  dc <- genData(nclust, defc)
  dp <- addPeriods(dc, nperiods, "cluster")
  dp <- trtStepWedge(dp, "cluster", nWaves = waves, lenWaves = len, 
                     startPer = start)
  dp <- addCorGen(dtOld = dp, nvars = nperiods, idvar = "cluster", 
                  rho = 0.6, corstr = "cs", dist = "normal", 
                  param1 = "mu", param2 = "s2", cnames = "cteffect")
  
  dd <- genCluster(dp, cLevelVar = "timeID", numIndsVar = 10, 
                   level1ID = "id")
  dd <- addColumns(defa, dd)
  dd[]
}

icc <- mclapply(1:200, 
                function(x) iccs(genDD(defc, defa, 100, 7, 4, 1, 2), byWave = T),
                mc.cores = 4
)

observed <- sapply(rbindlist(icc), function(x) mean(x))

bp <- rep((0.15/(0.15 + 2)) * 0.6, length(genPairs(7)))
wp <- rep(0.15/(0.15 + 2), 7)
expected <- round(c(bp, wp), 4)

pairs1 <- data.table(matrix(unlist(genPairs(7)), ncol = 2, byrow = TRUE))
pairs2 <- data.table(matrix(rep((1:7) - 1, each = 2), 
                            ncol = 2, byrow = TRUE))

results <- cbind(rbind(pairs1, pairs2), exp=expected, obs=observed)
setnames(results, c("V1", "V2"), c("Per1", "Per2"))

pres2 <- melt(results, id.vars = c("Per1", "Per2"))
pres2[Per1 == Per2, type := "within"]
pres2[Per1 != Per2, type := "between"]
pres2[variable == "exp", type := "expected"]

p2 <- ggplot(data = pres2, aes(x = variable, y = value, fill = type)) +
  geom_bar(stat = "identity") +
  facet_grid(Per1 ~ Per2) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank()) +
  scale_y_continuous(limits = c(0, 0.21), breaks = c(.05, .15)) +
  scale_fill_manual(values = c("#82a473", "grey80", "#9573a4"), 
                    breaks = c("between", "within"))

### AR1

defc <- defData(varname = "mu", formula = 0, 
                dist = "nonrandom", id = "cluster")
defc <- defData(defc, "s2", formula = 0.15, dist = "nonrandom")

defa <- defDataAdd(varname = "Y", 
                   formula = "0 + 0.10  * period + 1 * rx + cteffect", 
                   variance = 2, dist = "normal")

genDD <- function(defc, defa, nclust, nperiods, waves, len, start) {
  dc <- genData(nclust, defc)
  dp <- addPeriods(dc, nperiods, "cluster")
  dp <- trtStepWedge(dp, "cluster", nWaves = waves, lenWaves = len, 
                     startPer = start)
  dp <- addCorGen(dtOld = dp, nvars = nperiods, idvar = "cluster", 
                  rho = 0.6, corstr = "ar1", dist = "normal", 
                  param1 = "mu", param2 = "s2", cnames = "cteffect")
  
  dd <- genCluster(dp, cLevelVar = "timeID", numIndsVar = 10, 
                   level1ID = "id")
  dd <- addColumns(defa, dd)
  dd[]
}

icc <- mclapply(1:200, 
                function(x) iccs(genDD(defc, defa, 100, 7, 4, 1, 2), byWave = T),
                mc.cores = 4
)

observed <- sapply(rbindlist(icc), function(x) mean(x))

ts <- combn(0:6, m = 2)
t1 <- ts[1,]
t2 <- ts[2,]

bp <- (0.15/(0.15 + 2)) * (.6 ^ (t2 - t1))
wp <- rep(0.15/(0.15 + 2), 7)
expected <- round(c(bp, wp), 4)

pairs1 <- data.table(matrix(unlist(genPairs(7)), ncol = 2, byrow = TRUE))
pairs2 <- data.table(matrix(rep((1:7) - 1, each = 2), 
                            ncol = 2, byrow = TRUE))

results <- cbind(rbind(pairs1, pairs2), exp=expected, obs=observed)
setnames(results, c("V1", "V2"), c("Per1", "Per2"))

pres3 <- melt(results, id.vars = c("Per1", "Per2"))
pres3[Per1 == Per2, type := "within"]
pres3[Per1 != Per2, type := "between"]
pres3[variable == "exp", type := "expected"]

p3 <- ggplot(data = pres3, aes(x = variable, y = value, fill = type)) +
  geom_bar(stat = "identity") +
  facet_grid(Per1 ~ Per2) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank()) +
  scale_y_continuous(limits = c(0, 0.21), breaks = c(.05, .15)) +
  scale_fill_manual(values = c("#82a473", "grey80", "#9573a4"), 
                    breaks = c("between", "within"))

### Random slope

defc <- defData(varname = "ceffect", formula = 0, variance = 0.15, 
                dist = "normal", id = "cluster")
defc <- defData(defc, "cteffect", formula = 0, variance = 0.01, 
                dist = "normal")

defa <- defDataAdd(varname = "Y", 
                   formula = "0 + ceffect + 0.10  * period + cteffect * period + 1 * rx", 
                   variance = 2, dist = "normal")

genDD <- function(defc, defa, nclust, nperiods, waves, len, start) {
  dc <- genData(nclust, defc)
  dp <- addPeriods(dc, nperiods, "cluster")
  dp <- trtStepWedge(dp, "cluster", nWaves = waves, 
                     lenWaves = len, startPer = start)
  
  dd <- genCluster(dp, cLevelVar = "timeID", numIndsVar = 10, 
                   level1ID = "id")
  dd <- addColumns(defa, dd)
  dd[]
}

icc <- mclapply(1:200, 
                function(x) iccs(genDD(defc, defa, 100, 7, 4, 1, 2), byWave = T),
                mc.cores = 4
)

observed <- sapply(rbindlist(icc), function(x) mean(x))

ts <- combn(0:6, m = 2)
t1 <- ts[1,]
t2 <- ts[2,]
tw <- c(0:6)

bp <- (0.15 + t1*t2*0.01)/sqrt((0.15 + t1^2*0.01 + 2)*(0.15 + t2^2*0.01 + 2))
wp <- (0.15 + tw^2*0.01)/(0.15 + tw^2*0.01 + 2)
expected <- round(c(bp, wp), 4)

pairs1 <- data.table(matrix(unlist(genPairs(7)), ncol = 2, byrow = TRUE))
pairs2 <- data.table(matrix(rep((1:7) - 1, each = 2), 
                            ncol = 2, byrow = TRUE))

results <- cbind(rbind(pairs1, pairs2), exp=expected, obs=observed)
setnames(results, c("V1", "V2"), c("Per1", "Per2"))

pres4 <- melt(results, id.vars = c("Per1", "Per2"))
pres4[Per1 == Per2, type := "within"]
pres4[Per1 != Per2, type := "between"]
pres4[variable == "exp", type := "expected"]

p4 <- ggplot(data = pres4, aes(x = variable, y = value, fill = type)) +
  geom_bar(stat = "identity") +
  facet_grid(Per1 ~ Per2) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank()) +
  scale_y_continuous(limits = c(0, 0.21), breaks = c(.05, .15)) +
  scale_fill_manual(values = c("#82a473", "grey80", "#9573a4"), 
                    breaks = c("between", "within"))

save(pres1, pres2, pres3, pres4, 
  file = "content/post/DataICC/results.rdata")

ggsave("static/img/post-iccvary/p1.png", p1, height = 4.5, width = 8)
ggsave("static/img/post-iccvary/p2.png", p2, height = 4.5, width = 8)
ggsave("static/img/post-iccvary/p3.png", p3, height = 4.5, width = 8)
ggsave("static/img/post-iccvary/p4.png", p4, height = 4.5, width = 8)



