# define site level data using simstudy function

nsites = 40
spp = 10
starts <- paste0("rep(c(2 : ", (nsites/spp)+1, "), each = ", spp, ")")

siteDef <- defData(varname = "bj", dist = "normal", formula = 0, 
                   variance = .01, id="site")
siteDef <- defData(siteDef, varname = "sj", dist = "nonrandom", 
                   formula = starts)
siteDef <- defData(siteDef, varname = "ips", dist = "nonrandom", 
                   formula = 100)

indDef <- defDataAdd(varname = "bi", dist = "normal", formula = 0, 
                     variance = .01) # .01

trtDef <- defDataAdd(varname = "Ijt" , 
                     formula = "as.numeric(period >= sj)", 
                     dist = "nonrandom")

trtDef <- defDataAdd(trtDef, varname = "Yijt", 
                     formula = "0.85 + .005*period + 0.75*Ijt + 0.001*Ijt*(period-sj) + bi + bj", 
                     dist = "binary", link = "logit")

siteDef
indDef
trtDef

# Generate data

set.seed(789)

dtSite <- genData(nsites, siteDef)
dtSite <- genCluster(dtSite, cLevelVar = "site", numIndsVar = "ips", level1ID = "id")
dtSite <- addColumns(indDef, dtSite)

dtSiteTm <- addPeriods(dtSite, nPeriods = 7, idvars = "id")
dtSiteTm <- addColumns(trtDef, dtSiteTm)

# Plot data

dt <- dtSiteTm[, .(Y = mean(Yijt)), keyby = .(site, period, Ijt, sj)] # summary by site

ggplot(data = dt, aes(x=period, y=Y, group=site)) +
  geom_hline(yintercept = c(.7, .83),  color = "grey99") +
  geom_line(aes(color=factor(site))) +
  geom_point(data = dt[sj == period], color="grey50") +
  theme(panel.background = element_rect(fill = "grey90"), 
        panel.grid = element_blank(), 
        plot.title = element_text(size = 10, hjust = 0), 
        panel.border = element_rect(fill = NA, colour = "gray90"),
        legend.position = "none",
        axis.title.x = element_blank()
  ) +
  ylab("Proportion controlled") +
  scale_x_continuous(breaks = seq(0, 10, by = 2), 
                     labels = c("Baseline", paste("Year", c(1:5)))) +
  scale_y_continuous(limits = c(.5, 1), 
                     breaks = c(.5, .6, .7, .8, .9, 1)) +
  ggtitle("Stepped-wedge design with immediate effect") +
  facet_grid(sj~.)

# Not to be shown

nsites = 20
spp = 5
starts <- paste0("rep(c(2 : ", (nsites/spp)+1, "), each = ", spp, ")")

siteDef <- defData(varname = "bj", dist = "normal", formula = 0, 
                   variance = .01, id="site")
siteDef <- defData(siteDef, varname = "sj", dist = "nonrandom", 
                   formula = starts)
siteDef <- defData(siteDef, varname = "ips", dist = "nonrandom", 
                   formula = 50)

indDef <- defDataAdd(varname = "bi", dist = "normal", formula = 0, 
                     variance = .01) # .01

trtDef <- defDataAdd(varname = "Ijt" , 
                     formula = "as.numeric(period >= sj)", 
                     dist = "nonrandom")

trtDef <- defDataAdd(trtDef, varname = "Yijt", 
                     formula = "0.85 + .005*period + 0.3*Ijt + 0.001*Ijt*(period-sj) + bi + bj", 
                     dist = "binary", link = "logit")

dtSite <- genData(nsites, siteDef)
dtSite <- genCluster(dtSite, cLevelVar = "site", numIndsVar = "ips", level1ID = "id")
dtSite <- addColumns(indDef, dtSite)

dtSiteTm <- addPeriods(dtSite, nPeriods = 7, idvars = "id")
dtSiteTm <- addColumns(trtDef, dtSiteTm)

# Fit model - show this
library(lme4)

glmfit <- glmer(data = dtSiteTm, 
                Yijt ~ period + Ijt + I(Ijt*(period - sj)) + (1|id) + (1|site), 
                family="binomial" )
summary(glmfit)

gData <- function() {
  
  dtSite <- genData(nsites, siteDef)
  dtSite <- genCluster(dtSite, cLevelVar = "site", numIndsVar = "ips", level1ID = "id")
  dtSite <- addColumns(indDef, dtSite)
  
  dtSiteTm <- addPeriods(dtSite, nPeriods = 7, idvars = "id")
  dtSiteTm <- addColumns(trtDef, dtSiteTm)
  
  return(dtSiteTm)
}

result <- NULL

i=1

while (i < 100) {
  
  dtSite <- gData()
  
  glmfit <- tryCatch(glmer(data = dtSite, 
      Yijt ~ period + Ijt + I(Ijt*(period - sj)) + (1|id) + (1|site), 
      family="binomial" ),
    warning = function(w) {"warning"}
  )
  
  if (! is.character(glmfit)) {
    
    pvalue <- coef(summary(glmfit))["Ijt", "Pr(>|z|)"]
    result <- c(result, pvalue)
    i <- i + 1
  }
  
}

mean(result < .05) # 0.91

