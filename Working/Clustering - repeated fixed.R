library(broom)

nsites <- 24
iter <- 100

ICCs <- c(0.000, 0.025, 0.050, 0.075, 0.100)
vars <- ICCs/(1-ICCs)

probreject <- data.table()

for (cvar in vars) {
  
  defC <- defData(varname = "siteFix", 
                  formula = 0, variance = cvar, dist = "normal", id = "cID")
  
  defC <- defData(defC, varname = "nsite", formula = 50, dist = "nonrandom")
  
  defI <- defDataAdd(varname = "y", formula = "siteFix", variance = 1 )
  
  for (j in 1:100) {
    
    estp <- vector(mode = "numeric", length = iter)
    estpre <- vector(mode = "numeric", length = iter)
    
    dtC <- genData(nsites, defC) # random effect?
    
    for (i in 1:iter) {
      
      dtCrx <- trtAssign(dtC)
      
      dtI <- genCluster(dtClust = dtCrx, cLevelVar = "cID", 
                        numIndsVar = "nsite", level1ID = "id")
      dtI <- addColumns(defI, dtI)
      
      lmfit <- lm(y ~ trtGrp + factor(cID), data = dtI)
      lmerfit <- lme4::lmer(y ~ trtGrp + (1|cID), data = dtI)
      estp[i] <- tidy(lmfit)[2, 5]
      estpre[i] <- tidy(lmerfit)[2, 4]
      
    }
    
    probreject <- rbind(probreject, 
                        data.table(cvar, 
                                   j,
                                   preject = mean(estp <= .05),
                                   prejectre = mean(abs(estpre) >= 1.96)))
    
  }
    
}
  
