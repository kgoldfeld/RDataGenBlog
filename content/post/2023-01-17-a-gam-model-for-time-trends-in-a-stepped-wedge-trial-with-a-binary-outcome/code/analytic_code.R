library(simstudy)
library(ggplot2)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(gratia)
library(mgcViz)
library(DHARMa)

####

s_define <- function() {
  
  def <- defData(varname = "a", formula = 0, variance = .6)
  def <- defData(def, varname = "mu_b", formula = 0, dist = "nonrandom")
  def <- defData(def, varname = "s2_b", formula = .10, dist = "nonrandom")
  
  defOut <- defDataAdd(varname = "y", formula = "-1.5 + a + b + 0 * A", dist = "binary", link="logit")
  
  return(list(def = def, defOut = defOut)) 
}

s_generate <- function(list_of_defs) {
  
  list2env(list_of_defs, envir = environment())
  
  #--- add data generation code ---#
  
  ds <- genData(24, def, id = "site")
  ds <- addPeriods(ds, 25, "site", perName = "k")
  ds <- addCorGen(dtOld = ds, idvar = "site", 
                  rho = 0.95, corstr = "ar1",
                  dist = "normal", param1 = "mu_b", param2 = "s2_b", cnames = "b")
  ds <- trtStepWedge(ds, "site", nWaves = 24, lenWaves = 1, startPer = 1, 
                     grpName = "A", perName = "k")
  ds$site <- as.factor(ds$site)
  
  dd <- genCluster(ds, "timeID", numIndsVar = 100, level1ID = "id")
  dd <- addColumns(defOut, dd)
  dd[, normk := (k - min(k))/(max(k) - min(k))]
  
  return(dd) #  generated_data is a data.table
}

s_model <- function(dd) {
  
  # fitgam <- gamm(
  #   y ~ A + s(k, k=4) + s(k, site, bs = "fs", k=4),
  #   data = dd,
  #   family = "binomial"
  # )
  
  # res_fitgam <-  coef(summary(fitgam$lme))[2, c(1,2)]
  # model_results <- data.table(t(res_fitgam))
  
  fitbam <- bam(
    y ~ A + s(k, k=4) + s(k, site, bs = "fs", k=4), 
    data = dd, 
    method = "fREML",
    family = "binomial"
  )
  
  model_results <-  data.table(summary(fitbam)$p.coeff["A"], 
                               summary(fitbam)$se["A"],
                               summary(fitbam)$p.pv["A"])
  
  setnames(model_results, c("est.gamm", "se.gamm", "p.gamm"))
  
  return(model_results) # model_results is a data.table
}

s_single_rep <- function(list_of_defs) {
  
  generated_data <- s_generate(list_of_defs)
  model_results <- s_model(generated_data)
  
  return(model_results)
}

s_replicate <- function(nsim) {
  
  list_of_defs <- s_define()
  
  model_results <- rbindlist(
    pbapply::pblapply(
      X = 1 : nsim, 
      FUN = function(x) s_single_rep(list_of_defs), 
      cl = 4)
  )
  
  #--- add summary statistics code ---#
  
  return(model_results) # summary_stats is a data.table
}

#### Replication experiment

dres <- s_replicate(500)
dres[, .(est = mean(est.gamm), obs.se = sd(est.gamm), est.se = mean(se.gamm)) ]

dres[est.gamm > 0.44 & est.gamm < .84, .(est = mean(est.gamm), obs.se = sd(est.gamm), est.se = mean(se.gamm)) ]
dres[, quantile(est.gamm, p = c(.025, .975))]

#### Single data set

dd <- s_generate(s_define())

fitgam.A <- gamm(
  y ~ A + s(k, k = 4) + s(k, site, k = 4, bs = "fs"), 
  data = dd, 
  family = "binomial"
)

fitgam.noA <- gamm(
  y ~ s(k, k = 4) +  s(k, site, k = 4, bs = "fs") , 
  data = dd, 
  family = "binomial"
)

fitgam.maincurve <- gamm(
  y ~ A + s(k, k = 4)  , 
  data = dd, 
  family = "binomial"
)


summary(fitgam.A$gam)
summary(fitgam.A$lme)

gam.check(fitgam.A$gam)

# bam

fitbam.A <- bam(
  y ~ A + s(k) + s(k, site, bs = "fs"), 
  data = dd, 
  method = "fREML",
  family = "binomial"
)

fitbam.noA <- bam(
  y ~ s(k) + s(k, site, bs = "fs"), 
  data = dd, 
  method = "fREML",
  family = "binomial"
)

fitbam.1curve <- bam(
  y ~ A + s(k, k = 4)  , 
  data = dd, 
  method = "fREML",
  family = "binomial"
)


summary(fitbam.A)
gam.check(fitbam.A)

### Plot predicterd vs actual

dp <- dd[, .(.N, p = mean(y)), keyby = .(site, k, A)]

ddpred <- dd[, .SD[1,] , keyby = .(site, k)][, .(site, k, A)]
dp$ppred <- predict(fitgam.A$gam, ddpred, type = "response")

# dp$ppred <- predict(fitbam.A, ddpred, type = "response")

ggplot(data = dp, aes(x = k, y = ppred)) +
  geom_point(aes(x = k, y = p), size = .5, color = "grey75") +
  geom_line(aes(group = A, color = factor(A, labels = c("Control", "Intervention")))) +
  scale_color_manual(values = c("#d07b7c", "#7ba7d0")) +
  facet_wrap(~site, ncol = 8) +
  facet_wrap(~site, ncol = 8) +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 7),
        strip.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  ) 

### plot site specific curves

m <- fitgam.A$gam
# m <- fitbam.A
draw(m)

# customize curve

sm <- smooth_estimates(m)
data.table(sm)[smooth == "s(k,site)" & site == 2]
data.table(sm)[smooth == "s(k)"]
smooths(m)

sm <- add_confint(sm)

ggplot(data = data.table(sm)[smooth == "s(k)"], aes(y = est, x = k)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              alpha = 0.2, fill = "forestgreen") +
  geom_line(colour = "forestgreen", size = 1.5) +
  labs(y = "Partial effect",
       title = expression("Partial effect of" ~ f(x[2])),
       x = expression(x[2]))


### Compere different models

itsadug::compareML(fitbam.A, fitbam.1curve)
itsadug::compareML(fitbam.A, fitbam.noA)
itsadug::compareML(fitbam.1curve, fitbam.noA)

### Plot observed vs. model

sim <- simulate.gam(fitgam.A$gam, nsim = 1000)
# sim <- simulate.gam(fitgam.noA$gam, nsim = 1000)
# sim <- simulate.gam(fitgam.maincurve$gam, nsim = 1000)

ls <- split(sim, rep(1:ncol(sim), each = nrow(sim)))

dq <- lapply(ls, 
             function(x) cbind(dd, sim = x)[, .(obs = mean(y), sim = mean(sim)), keyby = .(site, k)]
)

dl <- rbindlist(dq, idcol = ".id")

df <- dl[, .(obs = mean(obs), min = quantile(sim, p = 0.025), 
             max = quantile(sim, 0.975)), keyby = .(site, k)]

ggplot(data = df, aes(x= k, y = obs)) +
  geom_ribbon(aes(x = k, ymin = min, ymax = max),
              alpha = 0.2, fill = "forestgreen") +
  geom_point(color = "forestgreen", size = 1) +
  
  facet_wrap( ~ site, ncol = 6) +
  theme(panel.grid = element_blank())


#### Evaluate goodness of fit

simResp <- matrix(dl$sim, nrow = 600)
obsResp <- dq[[1]]$obs

DHARMaRes = createDHARMa(simulatedResponse = simResp, observedResponse = obsResp, 
                         integerResponse = F)

plot(DHARMaRes, quantreg = T, testDispersion = F)

plotQQunif(DHARMaRes, testDispersion = F, testOutliers = F)
plotResiduals(DHARMaRes, quantreg = T)

#### Different curve by site/treatment

fitgam.test <- gamm(
  y ~ A + s(A, k, site, bs = "fs"), # + s(k, site, k = 4, bs = "fs"), 
  data = dd, 
  family = "binomial"
)

fitbam.test <- bam(
  y ~ A + s(A, k, site, bs = "fs"), # + s(k, site, k = 4, bs = "fs"), 
  data = dd, 
  method = "REML",
  family = "binomial"
)

summary(fitgam.test$gam)
summary(fitbam.test)

itsadug::compareML(fitbam.A, fitbam.test)

summary(fitbam.A)
summary(fitbam.test)





