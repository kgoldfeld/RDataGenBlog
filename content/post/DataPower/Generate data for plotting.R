library(parallel)
library(lme4)

defRE <- function(icc, dist = "binary", varW = NULL) {
  
  setVar <- iccRE(ICC = icc, dist = dist, varWithin = varW)
  def <- defData(varname = "a", formula = 0, variance = setVar, id = "cluster")
  
  return(def)
}

defBinOut <- function(p1, pctdelta) {
  
  p2 <- (1 - pctdelta) * p1
  
  int <- round(log( p1/(1-p1) ), 4)
  effect <- round(log( (p2/(1-p2)) / (p1/(1-p1) )), 4)
  formula <- genFormula( c(int, effect, 1), c("rx","a") )
  
  def <- defDataAdd(varname = "y", formula = formula, dist = "binary", 
                    link = "logit")
  return(def)
}

genDataSet <- function(nclust, clustsize, re.def, out.def) {
  
  dClust <- genData(nclust, re.def)
  dClust <- trtAssign(dClust, grpName = "rx")
  
  dPat <- genCluster(dtClust = dClust, cLevelVar = "cluster", 
                     numIndsVar = clustsize, level1ID = "id")
  dPat <- addColumns(out.def, dPat)
  
  return(dPat)
}

genBinEsts <- function(nclust, clustsize, re.def, out.def, fast = FALSE) {
  
  dP <- genDataSet(nclust, clustsize, re.def, out.def)
  
  mod.re <- glmer(y ~ rx + (1|cluster), data = dP, family = binomial,
                  control = glmerControl( optimizer = "bobyqa", calc.derivs = !(fast) ))
  
  convStatus <- as.numeric(length(summary(mod.re)$optinfo$conv$lme4))
  
  res <- data.table(convStatus, re = VarCorr(mod.re)$cluster,
                    t(coef(summary(mod.re))["rx",]))
  
  return(res)
}

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

nIters <- 2500
results <- NULL

ICC <- seq(0.025, 0.100, 0.025)
SS <- seq(800, 1600, 200)
nCLUST <- 40
ctlPROB <- c(0.10, 0.40)
pctDELTA <- 0.30

for (icc in ICC) {
  for (ss in SS) {
    for (nclust in nCLUST) {
      for (p1 in ctlPROB) {
        for (pdelta in pctDELTA) {
          
          clustsize <- ss %/% nclust
          p2 <- p1 * (1 - pdelta)
          
          defa <- defRE(icc)
          defy <- defBinOut(p1, pdelta)
          
          res <- rbindlist(mclapply(1:nIters, 
                        function(x) genBinEsts(nclust, clustsize, defa, defy)))
          
          dres <- data.table(icc, ss, nclust, clustsize, p1, p2, pdelta,
            converged = res[, mean(convStatus == 0)],
            # p.conv = res[convStatus == 0, mean(`Pr(>|z|)` < 0.05)],
            # p.all = res[convStatus != 2, mean(`Pr(>|z|)` < 0.05)],
            # avg.var.re.c = res[convStatus == 0, mean(`re.(Intercept)`)],
            # avg.est = res[convStatus != 2, mean(Estimate)], 
            # avg.est.c = res[convStatus == 0, mean(Estimate)], 
            obs.se.c = res[convStatus == 0 & `Std. Error` < 4, sd(Estimate)], 
            est.se.c = res[convStatus == 0 & `Std. Error` < 4, mean(`Std. Error`)], 
            sd.se.c = res[convStatus == 0 & `Std. Error` < 4, sd(`Std. Error`)],
            n.se = res[convStatus == 0 & `Std. Error` < 4, .N],
            n.bad = res[convStatus == 0 & `Std. Error` >= 4, .N]
          )
          
          print(c(icc, ss, nclust, clustsize, p1, p2, pdelta))
          results <- rbind(results, dres)
          
        }
      }
    }
  }
}

bin.clust.large <- copy(results)

save(bin.clust.large, file = "content/post/DataPower/binclustlarge.Rdata")

### Plot

library(paletteer)

ggplot(data = bin.clust.10, aes(x = pdelta, y = p.all)) +
  geom_line(aes(group = icc, color = factor(icc))) +
  facet_wrap(ss~., nrow = 2) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 11)) +
  scale_color_paletteer_d("jcolors::pal7", name = "ICC") +
  ylab("power") +
  xlab("effect size: percentage reduction") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1))

bin.clust.large[, `:=`(ymin = est.se.c - sd.se.c, ymax =  est.se.c + sd.se.c)]

ggplot(data = bin.clust.large, aes(x = ss, y = est.se.c)) +
  geom_line(aes(x = ss, y = obs.se.c, color = factor(icc))) +
  geom_errorbar(aes(x = ss, ymin = ymin, ymax = ymax, color = factor(icc)), 
                width = 0, alpha = .25) +
  geom_point(aes(color = factor(icc))) +
  facet_grid(icc~p1) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 11),
        legend.position = "none") +
 # scale_x_continuous(name = "sample size", limits = c(200, 480), 
#                     breaks = seq(200, 480, 80)) +
  scale_y_continuous(name = "standard error") +
  scale_color_paletteer_d("jcolors::pal7") +
  ggtitle("Estimated and observed standard errors by ICC and sample size (control proportion = 10% and percentage reduction = 40%)")
