library(lme4)
library(gee)
library(gridExtra)
library(simstudy)

genFunc <- function(nClusters, effVar) {
  
  # Define the data definitions for cluster and individual level data
  
  def1 <- defData(varname = "clustEff", formula = 0, variance = effVar, id = "cID")
  def1 <- defData(def1, varname = "nInd", formula = 100, dist = "noZeroPoisson")
  
  def2 <- defDataAdd(varname = "age", formula = 0, variance = 2)
  def2 <- defDataAdd(def2, varname = "Y", formula = "-4 + clustEff + 2*grp + 2*age", 
                     dist = "binary", link = "logit")
  
  
  # Generate cluster level data
  
  dtC <- genData(nClusters, def1)
  dtC <- trtAssign(dtC, grpName = "grp")
  
  # Generate individual level data
  
  dt <- genCluster(dtClust = dtC, cLevelVar = "cID", numIndsVar = "nInd", 
                   level1ID = "id")
  dt <- addColumns(def2, dt)
  
  return(dt)
  
}

glmerFit1 <- glmer(Y ~ grp + age + (1 | cID), data = dt, family="binomial")

geeFit1 <- gee(Y ~ grp + age, family = binomial, data = dt, 
               corstr = "exchangeable", id = dt$cID)

dtnew <- expand.grid(cID = unique(dt1$cID), age=seq(-4, 4, by =.1))
dtnew$grp <- 1


dtnewgee <- data.table(grp = 1, age = seq(-4, 4, by = .1))
glmFit1 <- glm(Y ~ grp + age, family = binomial, data = dt)

dtnew$dtpred <- predict(glmerFit1, dtnew, type = "response")
dtnewgee$pred <- predict(object = glmFit1, newdata = dtnewgee, type = "response")
dtnewgee$glmer <- predict(glmerFit1, dtnewgee[,c(1,2)], re.form = NA, type="response")

p1 <- ggplot(aes(x = age, y = dtpred), data=dtnew) + 
  geom_line(color="grey", aes(group = cID)) +
  geom_line(data=dtnewgee, aes(x = age, y = pred), color = "red", size = 1) +
  geom_line(data=dtnewgee, aes(x = age, y = glmer), color = "black", size = 1) +
  ggtitle("Intervention group") +
  xlab("Age") +
  ylab("Probability") 

    geom_vline(lty=3, xintercept = .8) +
  my_theme()
    

