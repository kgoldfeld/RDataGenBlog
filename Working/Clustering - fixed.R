library(broom)

nsites <- 100
iter <- 2000

probreject <- data.table()

for (j in seq(0, 0.7, by = .1)) {
  
  fixed <- seq(-j, j, length = nsites)
  
  conditionCol <- paste("cID ==", 1:nsites)
  
  dC <- data.table(condition = conditionCol, formula = fixed, 
                   variance = 0, dist="nonrandom", link="identity")
  
  defC <- defData(varname = "nsite", formula = 50, dist = "nonrandom", id = "cID")
  
  defI <- defDataAdd(varname = "y", formula = "siteFix", variance = 1 )
  
  dtC <- genData(nsites,defC) # fixed effect?
  dtC <- addCondition(dC, dtC, "siteFix")
  
  estp <- vector(mode = "numeric", length = iter)
  estpre <- vector(mode = "numeric", length = iter)
  
  for (i in 1:iter) {
    
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
                      data.table(j, 
                                 preject = mean(estp <= .05),
                                 prejectre = mean(abs(estpre) >= 1.96)))
}

probf <- probreject


