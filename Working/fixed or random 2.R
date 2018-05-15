library(broom)

iter <- 2500

probreject <- data.table()

for (cvar in vars) {
  
  def
  
  defC <- defData(varname = "siteFix", 
                  formula = 0, variance = cvar, dist = "normal", id = "cID")
  
  defC <- defData(defC, varname = "nsite", formula = 50, dist = "nonrandom")
  
  defI <- defDataAdd(varname = "y", formula = "siteFix", variance = 1 )
  
  dtC <- genData(6,defC) # fixed effect?
  
  estp <- vector(mode = "numeric", length = iter)
  
  for (i in 1:iter) {
    
    dtCrx <- trtAssign(dtC)
    
    dtI <- genCluster(dtClust = dtCrx, cLevelVar = "cID", 
                      numIndsVar = "nsite", level1ID = "id")
    dtI <- addColumns(defI, dtI)
    
    lmfit <- lm(y~ trtGrp + factor(cID) - 1, data = dtI)
    estp[i] <- tidy(lmfit)[1,5]
    
  }
  
  probreject <- rbind(probreject, data.table(cvar, preject = mean(estp <= .05)))
  
}




