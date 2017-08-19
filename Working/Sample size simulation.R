library(lme4)
library(multcomp)

# Definitions

genImpact <- function(nsites, spp, ips = 500) {
  
  # spp - sites per period
  # ips - individuals per site
  
  # starting time depends on number of sites and number of sites starting per period
  
  starts <- "rep(c(2 : (end+1)), each = times)"
  starts <- gsub("end", nsites/spp, starts )
  starts <- gsub("times", spp, starts )
  
  # define site level data using simstudy function
  
  siteDef <- defData(varname = "bj", dist = "normal", formula = 0, variance = .01, id="site")
  siteDef <- defData(siteDef, varname = "sj", dist = "nonrandom", formula = starts)
  siteDef <- defData(siteDef, varname = "ips", dist = "nonrandom", formula = ips)
  
  # define individual level data - might not need this for your study
  
  indDef <- defDataAdd(varname = "bi", dist = "normal", formula = 0, variance = .01) # .01
  
  # define outcome
  
  trtDef <- defDataAdd(varname = "Ijt" , formula = "as.numeric(period >= sj)", dist = "nonrandom")
  trtDef <- defDataAdd(trtDef,
                       varname = "Yijt" ,
                       formula = "0.85 + .005*period + 0.675*Ijt + 0.001*Ijt*(period-sj) + bi + bj",
                       dist = "binary",
                       link = "logit"
  )
  
  # Generate data
  
  dtSite <- genData(nsites, siteDef)
  dtSite <- genCluster(dtSite, cLevelVar = "site", numIndsVar = "ips", level1ID = "id")
  dtSite <- addColumns(indDef, dtSite)
  
  dtSiteTm <- addPeriods(dtSite, nPeriods = 12, idvars = "id")
  dtSiteTm <- addColumns(trtDef, dtSiteTm)
  
  return(dtSiteTm)
  
}

#### Exammple data

set.seed(435)

dtAll <- genImpact(12, 3, 50)


# Fit model

glmfit <- glmer(data = dtAll, 
                Yijt ~ period + Ijt + I(Ijt*(period - sj)) + (1|id) + (1|site), 
                family="binomial" )
summary(glmfit)

# Plot

dt <- dtAll[, .(Y = mean(Yijt)), keyby = .(site, period, Ijt, sj)] # summary by site

ggplot(data = dt, aes(x=period, y=Y, group=site)) +
  geom_hline(yintercept = c(.7, .82),  color = "grey99") +
  geom_line(aes(color=factor(site))) +
  geom_point(data = dt[sj == period], color="grey50") +
  theme_ksg("grey90") +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  ylab("Proportion controlled") +
  scale_x_continuous(breaks = seq(0, 10, by = 2), 
                     labels = c("Baseline", paste("Year", c(1:5)))) +
  scale_y_continuous(limits = c(.5, 1), 
                     breaks = c(.5, .6, .7, .8, .9, 1)) +
  ggtitle("Scenario 1: immediate effect")


#### Iterate for power calculations

result <- NULL

i=1

while (i < 250) {
  
  dtSite <- genImpact(nsites = 15, spp = 3, ips = 30)
  
  glmfit <- tryCatch(glmer(data = dtSite, 
                  Yijt ~ period + Ijt + I(Ijt*(period - sj)) + (1|id) + (1|site), 
                  family="binomial" ),
                  warning = function(w) {
                    "warning"
                  }
  )
  
  if (! is.character(glmfit)) {
    
    K <- rbind( "Effect after 1 year" = c(0, 0, 1, 3))
    pvalue <- summary(glht(glmfit, K))$test$pvalues
    
    # pvalue <- coef(summary(glmfit))["Ijt", "Pr(>|z|)"] # 1st period effect only
    
    result <- c(result, pvalue)
    i <- i + 1
  }
  
}

mean(result < .05)


