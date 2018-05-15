library(broom)

nsites <- 6
iter <- 2000

vars <-seq(0, 0.42, by = .06)^2
probreject <- data.table()

for (cvar in vars) {
  
  defC <- defData(varname = "siteFix", 
                  formula = 0, variance = cvar, dist = "normal", id = "cID")
  
  defC <- defData(defC, varname = "nsite", formula = 50, dist = "nonrandom")
  
  defI <- defDataAdd(varname = "y", formula = "siteFix", variance = 1 )
  
  estp <- vector(mode = "numeric", length = iter)
  estpre <- vector(mode = "numeric", length = iter)
  
  
  for (i in 1:iter) {
    
    dtC <- genData(nsites, defC) # random effect?
    
    dtCrx <- trtAssign(dtC)
    
    dtI <- genCluster(dtClust = dtCrx, cLevelVar = "cID", 
                      numIndsVar = "nsite", level1ID = "id")
    dtI <- addColumns(defI, dtI)
    
    lmfit <- lm(y~ trtGrp + factor(cID), data = dtI)
    # lmerfit <- lme4::lmer(y ~ trtGrp + (1|cID), data = dtI)
    estp[i] <- tidy(lmfit)[2, 5]
    # estpre[i] <- tidy(lmerfit)[2, 4]
    
  }
  
  probreject <- rbind(probreject, 
                      data.table(cvar, 
                                 preject = mean(estp <= .05),
                                 prejectre = mean(abs(estpre) >= 1.96)))
  
}

probc <- probreject



