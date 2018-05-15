library(broom)

nsites <- 6
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

    dtC <- genData(nsites, defC) # random effect?
    
    for (i in 1:iter) {
      
      dtI <- genCluster(dtClust = dtC, cLevelVar = "cID", 
                        numIndsVar = "nsite", level1ID = "id")
      dtI <- addColumns(defI, dtI)
      dtI <- trtAssign(dtI, strata = "cID")
      
      lmfit <- lm(y ~ trtGrp + factor(cID), data = dtI)
      estp[i] <- tidy(lmfit)[2, 5]
    }
    
    probreject <- rbind(probreject, 
                        data.table(cvar, 
                                   j,
                                   preject = mean(estp <= .05)))
    
  }
    
}
  
