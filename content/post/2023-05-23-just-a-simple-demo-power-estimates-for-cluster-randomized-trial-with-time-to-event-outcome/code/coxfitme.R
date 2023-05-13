library(simstudy)
library(data.table)
library(coxme)

points <- list(c(30, 0.90), c(365, .50))
r <- survGetParams(points)

defc <- defData(varname = "b", formula = 0, variance = 0.05)
defc <- defData(defc, varname = "rx", formula = "1;1", dist = "trtAssign")

defa <- defDataAdd(varname = "start_day", formula = "1;182", dist = "uniformInt")
defa <- defDataAdd(defa, varname = "censor", 
                   formula = "365 - start_day ", dist = "nonrandom")

defs <- defSurv(varname = "ttc", formula = "r[1] + 0.4 * rx + b", shape = r[2])

dc <- genData(1000, defc, id = "site")
dd <- genCluster(dc, "site", numIndsVar = 500, "id")
dd <- addColumns(defa, dd)
dd <- genSurv(dd, defs, digits = 0)
dd <- addCompRisk(dd, events = c("ttc", "censor"), 
                  timeName = "time", censorName = "censor", keepEvents = TRUE)

fit_coxme <-summary(coxme(Surv(time, event) ~ rx + (1 | site), data = dd))

save(fit_coxme, file = "data/coxfit.rda")
