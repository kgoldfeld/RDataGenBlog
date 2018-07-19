library(lme4)
library(geepack)
library(broom)

iter <- 500
i <- 0

resWithin <- list()
resByClust <- list()

for (icc in c(.05, .20)) {
  
  cvar <- iccRE(icc, dist = "binary")
  
  # define data
  
  d <- defData(varname = "a", formula = 0, variance = cvar, 
               dist = "normal", id = "cid")
  d <- defData(d, varname = "nper", formula = 100, dist = "nonrandom")
  
  da <- defDataAdd(varname = "y", formula = "-1 + .4*rx + a", dist="binary", 
                   link = "logit")
  
  for (nCluster in c(30, 100)) {
    
    ### Randomization within cluster
    
    for (j in 1:iter) {
      
      i <- i + 1
      
      dc <- genData(nCluster, d)
      
      di <- genCluster(dc, "cid", "nper", "id")
      di <- trtAssign(di, strata = "cid", grpName = "rx")
      di <- addColumns(da, di)
      
      glmerfit <- glmer(y~rx + (1|cid), data = di, family = binomial)
      estC <- data.table(tidy(glmerfit))[term == "rx", estimate]
      
      obspc <- di[, .(estp = mean(y)), keyby = .(cid, rx)]
      obspc[estp == 0, estp := 0.005]
      obspc[estp == 1, estp := 0.995]
      obspc[, lo := log(odds(estp))]
      rawC <- obspc[rx == 1, mean(lo)] - obspc[rx == 0, mean(lo)] 
      
      theoryC <- .4
      
      glmfit <- glm(y ~ rx, data = di, family = binomial)
      estM <- data.table(tidy(glmfit))[term == "rx", estimate]
      
      rawM <- log(odds(di[rx==1, mean(y)])/odds(di[rx==0, mean(y)]))
      
      resWithin[[i]] <- data.table(i, random = "within", icc, nCluster, 
                                   estC, rawC, theoryC, estM, rawM)
      
    ### Randomization by cluster
    
      dc <- genData(nCluster, d)
      dc <- trtAssign(dc, grpName = "rx")
      
      di <- genCluster(dc, "cid", "nper", "id")
      di <- addColumns(da, di)
      
      glmerfit <- glmer(y~rx + (1|cid), data = di, family = binomial)
      estC <- data.table(tidy(glmerfit))[term == "rx", estimate]
      
      obspc <- di[, .(estp = mean(y)), keyby = .(cid, rx)]
      obspc[estp == 0, estp := 0.005]
      obspc[estp == 1, estp := 0.995]
      obspc[, lo := log(odds(estp))]
      rawC <- obspc[rx == 1, mean(lo)] - obspc[rx == 0, mean(lo)] 
      
      theoryC <- .4 + dc[rx == 1, mean(a)] - dc[rx == 0, mean(a)] 
      
      glmfit <- glm(y ~ rx, data = di, family = binomial)
      estM <- data.table(tidy(glmfit))[term == "rx", estimate]
      
      rawM <- log(odds(di[rx==1, mean(y)])/odds(di[rx==0, mean(y)]))
      
      resByClust[[i]] <- data.table(i, random = "bycluster", icc, nCluster, 
                                    estC, rawC, theoryC, estM, rawM)
    }
  }
}

resWithinB <- rbindlist(resWithin)
resByClustB <- rbindlist(resByClust)
resB <- rbindlist(list(resWithinB, resByClustB))

dm <- melt(resB, id.vars = c("i", "random", "icc", "nCluster"))
dm <- dm[variable != "theoryC"]

dm[variable %in% c("estC", "rawC", "theoryC"), type := "Conditional"]
dm[!(variable %in% c("estC", "rawC", "theoryC")), type := "Marginal"]

dm[, type := factor(type, levels = c("Conditional", "Marginal"))]

dm[variable %in% c("estC", "estM"), est := "Model"]
dm[variable %in% c("rawC", "rawM"), est := "Calculation"]

dm[, est := factor(est, levels = c("Model", "Calculation"))]

# Model only (not calculation)

dM <- dm[est == "Model" & nCluster == 100]
dM <- dm[est == "Model" & nCluster == 30]

dM[, params := paste0(type, "\n", random)]
dM[, icc := paste0("icc = ", format(icc, digits = 2))]

davg <- dM[, .(avg=round(mean(value),2), sd=round(sd(value), 2)), 
           keyby = .(icc, params)]

davg[, text := paste0("avg (sd) = ", format(avg, digits = 2), 
                       " (", 
                       format(sd, digits=2), 
                       ")")]

ggplot(data = dM, aes(x=value)) +
  geom_histogram(binwidth = 0.04, aes(fill = type)) +
  # geom_vline(aes(xintercept = avg), data= davg, color = "grey92") +
  geom_text(x = 1, y = 160, aes(label = text), data= davg, size = 2.5) +
  scale_x_continuous(limits = c(-.5, 1.5), 
                     breaks = seq(-.4, 1.2, by = .4),
                     name = "estimated effect") +
  scale_y_continuous( limits = c(0, 200)) +
  scale_fill_manual(values = c("#8594e6","#e6d785")) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey92"),
        legend.position = "none") +
  facet_grid(icc ~ params) 
  



# Random and Model Type

dM <- dm[est == "Model" & icc == 0.20 & nCluster == 30]
dM[, random := factor(random, levels = c("within", "bycluster"),
                     labels = c("Randomize within cluster", "Randomize by cluster"))]

davg <- dM[, .(avg=round(mean(value),2), sd=round(sd(value), 2)), 
           keyby = .(type, random)]

davg[, text := paste0("avg (sd) = ", format(avg, digits = 2), 
                      " (", 
                      format(sd, digits=2), 
                      ")")]

pRT <- ggplot(data = dM, aes(x=value)) +
  geom_histogram(binwidth = 0.04, aes(fill = type)) +
  # geom_vline(aes(xintercept = avg), data= davg, color = "grey92") +
  geom_text(x = 1, y = 160, aes(label = text), data= davg, size = 2) +
  scale_x_continuous(limits = c(-.5, 1.5), 
                     breaks = seq(-.4, 1.2, by = .4),
                     name = "estimated effect") +
  scale_y_continuous( limits = c(0, 200)) +
  scale_fill_manual(values = c("#8594e6","#e6d785")) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey92"),
        legend.position = "none",
        text = element_text(size = 8)) +
  facet_grid(random ~ type) 

# ICC & Model Type

dM <- dm[est == "Model" & random == "bycluster" & nCluster == 30]
dM[, icc := paste0("icc = ", format(icc, digits = 2))]

davg <- dM[, .(avg=round(mean(value),2), sd=round(sd(value), 2)), 
           keyby = .(type, icc)]

davg[, text := paste0("avg (sd) = ", format(avg, digits = 2), 
                      " (", 
                      format(sd, digits=2), 
                      ")")]

pIT <- ggplot(data = dM, aes(x=value)) +
  geom_histogram(binwidth = 0.04, aes(fill = type)) +
  geom_text(x = 1, y = 160, aes(label = text), data= davg, size = 2) +
  scale_x_continuous(limits = c(-.5, 1.5), 
                     breaks = seq(-.4, 1.2, by = .4),
                     name = "estimated effect") +
  scale_y_continuous( limits = c(0, 200)) +
  scale_fill_manual(values = c("#8594e6","#e6d785")) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey92"),
        legend.position = "none",
        text = element_text(size = 8)) +
  facet_grid(icc ~ type) 

# Number clusters & Model Type

dM <- dm[est == "Model" & random == "bycluster" & icc == 0.20]
dM[, nCluster := paste0("n = ", nCluster)]

davg <- dM[, .(avg=round(mean(value),2), sd=round(sd(value), 2)), 
           keyby = .(type, nCluster)]

davg[, text := paste0("avg (sd) = ", format(avg, digits = 2), 
                      " (", 
                      format(sd, digits=2), 
                      ")")]

pNT <- ggplot(data = dM, aes(x=value)) +
  geom_histogram(binwidth = 0.04, aes(fill = type)) +
  geom_text(x = 1, y = 160, aes(label = text), data= davg, size = 2) +
  scale_x_continuous(limits = c(-.5, 1.5), 
                     breaks = seq(-.4, 1.2, by = .4),
                     name = "estimated effect") +
  scale_y_continuous( limits = c(0, 200)) +
  scale_fill_manual(values = c("#8594e6","#e6d785")) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey92"),
        legend.position = "none",
        text = element_text(size = 8)) +
  facet_grid(nCluster ~ type) 


ggsave("static/img/post-condmarg/pNT.png", pNT, height = 3, width = 5)
ggsave("static/img/post-condmarg/pRT.png", pRT, height = 3, width = 5)
ggsave("static/img/post-condmarg/pIT.png", pIT, height = 3, width = 5)


