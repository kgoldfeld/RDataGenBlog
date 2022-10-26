setwd("~/git R projects/rdatagen/content/post/2022-10-11-modeling-secular-trend-in-CRT-using-gam")

library(simstudy)
library(ggplot2)
library(cowplot)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(parallel)
library(pbapply)

def <- defData(varname = "b0", formula = 0, variance = 6)
def <- defData(def, varname = "A", formula = "1;1", dist = "trtAssign")
def <- defData(def, varname = "mu", formula = 0, dist = "nonrandom")
def <- defData(def, varname = "s2", formula = 16, dist = "nonrandom")

defOut <- defDataAdd(varname = "y", 
                     formula = "100 + b0 + b1k - 0.1 * k^2 + 5*A", 
                     variance = 9)

s_generate <- function() {
  
  dd <- genData(48, def, id = "site")
  dd <- addPeriods(dd, 20, "site", perName = "k")
  
  dd <- addCorGen(dtOld = dd, idvar = "site", nvars = 20, 
                  rho = .7, corstr = "ar1",
                  dist = "normal", param1 = "mu", param2 = "s2", cnames = "b1k")
  
  dd <- genCluster(dd, "timeID", numIndsVar = 30, level1ID = "id")
  dd <- addColumns(defOut, dd)
  
  dd[, normk := (k - min(k))/(max(k) - min(k))]
  dd[, site := as.factor(site)]
  
  dd[]
}

replicate <- function(){
  
  dd <- s_generate()
  
  linear <- lmer(y ~ A + k + ( 1  | site) , data = dd)
  est.lin <- summary(linear)$coefficients["A", c("Estimate", "Std. Error")]
  
  fix_cs <- lmer(y ~ A + bs(normk) + ( 1  | site) , data = dd)
  est.fcs <- summary(fix_cs)$coefficients["A", c("Estimate", "Std. Error")]
  
  ran_cs <- lmer(y ~ A + ( bs(normk) | site) , data = dd)
  est.rcs <- summary(ran_cs)$coefficients["A", c("Estimate", "Std. Error")]
  
  gam <- gamm(y ~ A + s(k, site, bs = "fs", k = 5), data = dd, method="REML")
  est.gam <- cbind(summary(gam$gam)$p.coeff, summary(gam$gam)$se)[2,]
  
  dres <- data.table(t(est.lin), t(est.fcs), t(est.rcs), t(est.gam))
  setnames(dres, c("d.lin", "se.lin","d.fcs","se.fcs","d.rcs","se.rcs", "d.gam", "se.gam"))
  
  dres[]
}

res <- rbindlist(pblapply(1:1000, function(x) replicate()))
save(res, file = "replication_code/res.rdata")
