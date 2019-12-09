library(geepack)

iter <- function(n, np, defM) {
  
  dx <- genData(n)
  dx <- trtAssign(dx, grpName = "rx")
  dx <- addPeriods(dx, nPeriods = np)
  
  def <- defDataAdd(varname = "p", formula = "-2 + 0.2*period + .5*rx", 
                    dist = "nonrandom", link = "logit")
  
  dx <- addColumns(def, dx)
  dx <- addCorGen(dtOld = dx, idvar = "id", nvars = np, rho = .6, 
                  corstr = "ar1", dist = "binary", param1 = "p", 
                  method = "ep", formSpec = "-2 + 0.2*period + .5*rx",
                  cnames = "y")
  
  dm <- genMiss(dx, defM, "id", repeated = TRUE, periodvar = "period")
  nmiss <- length(dm[y==1, unique(id)])
  tmiss <- dm[, sum(y)]
  
  dObs <- genObs(dx, dm, idvars = "id")
  
  fit.f <- geeglm(y ~ period + rx, id = id, family = binomial, 
         data = dx, corstr = "ar1")
  
  fit.m <- geeglm(y ~ period + rx, id = id, family = binomial, 
         data = dObs, corstr = "ar1")
  
  fit.l <- glm(y ~ rx, data = dObs[period == (np - 1)], family = binomial)
  
  return(data.table(full = coef(fit.f)["rx"], 
                    miss = coef(fit.m)["rx"],
                    last = coef(fit.l)["rx"])
         )
}

## defM

MCAR <- defMiss(varname = "y", formula = "-2.6",
                logit.link = TRUE, monotonic = TRUE
)

MAR <- defMiss(varname = "y", 
               formula = "-2.9 + 0.2*period - 2*rx*LAG(y)",
               logit.link = TRUE, monotonic = TRUE
)

NMAR <- defMiss(varname = "y", 
                formula = "-2.9 + 0.2*period - 2*rx*y",
                logit.link = TRUE, monotonic = TRUE
)

##

iter(200, 5, MCAR)

##



library(parallel)

niter <- 2500

resMCAR <- rbindlist(mclapply(1:niter, function(x) iter(200, 5, MCAR)))
resMAR <- rbindlist(mclapply(1:niter, function(x) iter(200, 5, MAR)))
resNMAR <- rbindlist(mclapply(1:niter, function(x) iter(200, 5, NMAR)))

save(resMCAR, resMAR, resNMAR, 
     file = "content/post/DataRepeated/repeated.rdata" )

mrx <- function(res, n) {
  
  dm <- melt(res, measure.vars = c("full","miss","last")) 
  
  dm[, .(n,
         l95 = quantile(value, 0.025), 
         mu = mean(value), 
         u95 = quantile(value, 0.975)), 
     keyby = variable]
  
}

mylist <- list(resMCAR = resMCAR, resMAR = resMAR, resNMAR = resNMAR)

dsum <- lapply(seq_along(mylist), 
  function(y, n, i) mrx(y[[i]], n[[i]]), y = mylist, n=names(mylist))
dsum <- rbindlist(dsum)
dsum[ , variable := factor(variable, levels = c("last", "miss", "full"),
        labels = c("last period only", "observed", "actual"))]
dsum[, n := factor(n, levels = c("resMCAR", "resMAR", "resNMAR"), 
                   labels = c("MCAR", "MAR", "NMAR"))]

ggplot(data = dsum, aes(x=variable, y = mu)) +
  geom_segment(x = -Inf, xend = Inf, y = .5, yend = .5, 
               size = 1, color = "white") +
  geom_segment(aes(y = u95, yend = l95, x = variable, xend = variable),
               color = "#1a5dff") +
  geom_point(color = "#002d99") +
  coord_flip() +
  facet_grid(n ~ .) +
  theme(panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.title.y  = element_blank(),
        axis.ticks = element_blank()) +
  scale_y_continuous(limits  = c(-.75, 1.6 ), 
                     breaks = c(-.5, 0, .5, 1, 1.5),
                     name = "\nestimated treatment effect")





       