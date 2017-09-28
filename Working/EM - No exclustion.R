library(Weighted.Desc.Stat)
library(ivpack)
library(simstudy) # This is a pacakge I wrote to gnerate simulated data ...

# E-step

estep <- function(params, y, z, m) {
  
  piC <- 0
  piN <- 0
  piA <- 0
  
  if (z == 0 & m == 0) {
    
    gC0 <- dnorm((y - params$mC0)/params$sC0) / params$sC0
    gN0 <- dnorm((y - params$mN0)/params$sN0) / params$sN0
    
    piC <- params$pC * gC0 / ( params$pC * gC0 + params$pN * gN0)
    piN <- 1- piC
    
  }
  
  if (z == 0 & m == 1) {
    piA <- 1
  }
  
  if (z == 1 & m == 0) {
    piN <- 1
  }
  
  if (z == 1 & m == 1) {
    
    gC1 <- dnorm((y - params$mC1)/params$sC1) / params$sC1
    gA1 <- dnorm((y - params$mA1)/params$sA1) / params$sA1
    
    piC <- params$pC * gC1 / ( params$pC * gC1 + params$pA * gA1)
    piA <- 1 - piC
  }
  
  return(list(piC = piC, piN = piN, piA = piA))
  
}

mstep <- function(params, dx) {
  
  params$mN0 <- dx[z == 0 & m == 0, w.mean(y, piN)] # never-taker
  params$sN0 <- dx[z == 0 & m == 0, sqrt(w.var(y, piN))] # never-taker
  
  params$mN1 <- dx[z == 1 & m == 0, w.mean(y, piN)] # never-taker
  params$sN1 <- dx[z == 1 & m == 0, sqrt(w.var(y, piN))] # never-taker
  
  params$mA0 <- dx[z == 0 & m == 1, w.mean(y, piA)]# always-taker
  params$sA0 <- dx[z == 0 & m == 1, sqrt(w.var(y, piA))] # always-taker
  
  params$mA1 <- dx[z == 1 & m == 1, w.mean(y, piA)]# always-taker
  params$sA1 <- dx[z == 1 & m == 1, sqrt(w.var(y, piA))] # always-taker
  
  params$mC0 <- dx[z == 0 & m == 0, w.mean(y, piC)] # complier, z=0
  params$sC0 <- dx[z == 0 & m == 0, sqrt(w.var(y, piC))] # complier, z=0
  
  params$mC1 <- dx[z == 1 & m == 1, w.mean(y, piC)] # complier, z=1
  params$sC1 <- dx[z == 1 & m == 1, sqrt(w.var(y, piC))] # complier, z=1
  
  nC <- dx[, sum(piC)]
  nN <- dx[, sum(piN)]
  nA <- dx[, sum(piA)]
  
  params$pC <- (nC / sum(nC, nN, nA))
  params$pN <- (nN / sum(nC, nN, nA))
  params$pA <- (nA / sum(nC, nN, nA))
  
  return(params)
}

like.i <- function(params, y, z, m) {
  
  if (z == 0 & m == 0) {
    l <- params$pC * dnorm(x = y, mean = params$mC0, sd = params$sC0) +
      params$pN * dnorm(x = y, mean = params$mN0, sd = params$sN0)
  }
  
  if (z == 0 & m == 1) {
    l <- params$pA * dnorm(x = y, mean = params$mA0, sd = params$sA0)
  }
  
  if (z == 1 & m == 0) {
    l <- params$pN * dnorm(x = y, mean = params$mN1, sd = params$sN1)
  }
  
  if (z == 1 & m == 1) {
    l <- params$pC * dnorm(x = y, mean = params$mC1, sd = params$sC1) +
      params$pA * dnorm(x = y, mean = params$mA1, sd = params$sA1)
  }
  
  return(l)
}

loglike <- function(dt, params){
  
  dl <- dt[, .(l.i = like.i(params, y, z, m)), keyby = id]
  return(dl[, sum(log(l.i))])
  
}

initparams <- function() {
  
  params = list(pC = 1/3, pN = 1/3, pA = 1/3, 
                mC0 = rnorm(1,0,.1), sC0 = 0.2,
                mC1 = rnorm(1,0,.1), sC1 = 0.2, 
                mN0 = rnorm(1,0,.1), sN0 = 0.2,
                mN1 = rnorm(1,0,.1), sN1 = 0.2,
                mA0 = rnorm(1,0,.1), sA0 = 0.2,
                mA1 = rnorm(1,0,.1), sA1 = 0.2)
  
  return(params)
}

### Define data distributions

# Status :

# 1 = A(lways taker)
# 2 = N(ever taker)
# 3 = C(omplier)

def <- defDataAdd(varname = "Status", 
                  formula = "0.25; 0.40; 0.35", dist = "categorical")

# potential outcomes (PO) for intervention depends on group status

def <- defDataAdd(def, varname = "M0", 
                  formula = "(Status == 1) * 1", dist = "nonrandom")
def <- defDataAdd(def, varname = "M1", 
                  formula = "(Status != 2) * 1", dist = "nonrandom")

# observed intervention status based on randomization and PO

def <- defDataAdd(def, varname = "m", 
                  formula = "(z==0) * M0 + (z==1) * M1", dist = "nonrandom")

# potential outcome for Y (depends group status - A, N, or C)

defY0 <- defCondition(condition = "Status == 1",
                      formula = 0.3, variance = .20, dist = "normal")
defY0 <- defCondition(defY0, condition = "Status == 2",
                      formula = 0.0, variance = .36, dist = "normal")
defY0 <- defCondition(defY0, condition = "Status == 3",
                      formula = 0.1, variance = .16, dist = "normal")

defY1 <- defCondition(condition = "Status == 1",
                      formula = 0.7, variance = .25, dist = "normal")
defY1 <- defCondition(defY1, condition = "Status == 2",
                      formula = 0.2, variance = .40, dist = "normal")
defY1 <- defCondition(defY1, condition = "Status == 3",
                      formula = 0.9, variance = .49, dist = "normal")

# potential outcome for Y (depends group status - A, N, or C)

# defY0 <- defCondition(condition = "Status == 1",
#                       formula = 0.5, variance = .25, dist = "normal")
# defY0 <- defCondition(defY0, condition = "Status == 2",
#                       formula = 0.0, variance = .36, dist = "normal")
# defY0 <- defCondition(defY0, condition = "Status == 3",
#                       formula = 0.1, variance = .16, dist = "normal")
# 
# defY1 <- defCondition(condition = "Status == 1",
#                       formula = 0.5, variance = .25, dist = "normal")
# defY1 <- defCondition(defY1, condition = "Status == 2",
#                       formula = 0.0, variance = .36, dist = "normal")
# defY1 <- defCondition(defY1, condition = "Status == 3",
#                       formula = 0.9, variance = .49, dist = "normal")

# observed outcome function of actual treatment

defy <- defDataAdd(varname = "y", 
                   formula = "(z == 0) * Y0 + (z == 1) * Y1", dist = "nonrandom")


createDT <- function(n, def, defY0, defY1, defy) {
  
  dt <- genData(n)
  dt <- trtAssign(dt, n=2, grpName = "z")
  dt <- addColumns(def, dt)
  
  genFactor(dt, "Status", labels = c("Always-taker","Never-taker", "Complier"), prefix = "A")
  
  dt <- addCondition(defY0, dt, "Y0")
  dt <- addCondition(defY1, dt, "Y1")
  dt <- addColumns(defy, dt)
  
}

# Simulation 

set.seed(123)

estCE <- data.table()
paramsEst <- data.table()

for (i in 1:500) {
  
  dt <- createDT(600, def, defY0, defY1, defy)
  
  params <- initparams()
  
  steps <- data.table()
  prev.loglike <- -Inf
  continue <- TRUE
  
  while (continue) {
    
    dtPIs <- dt[, estep(params, y, z, m), keyby = id]
    dx <- dt[dtPIs]
    
    params <- mstep(params, dx)
    current.loglike <- loglike(dt, params)
    
    steps <- rbind(steps,
                   data.table(AACE = params$mA1 - params$mA0,
                              NACE = params$mN1 - params$mN0,
                              CACE = params$mC1 - params$mC0,
                              current.loglike)
    )
    
    diff <- current.loglike - prev.loglike
    prev.loglike <- current.loglike
    
    if ( diff < 1.00e-05 ) continue = FALSE
    
  }
  
  ivmodel <- ivreg(formula = y ~ m | z, data = dt, x = TRUE)
  
  ITT <- dt[z==1, mean(y)] - dt[z==0, mean(y)] #ITT
  ITTm <- dt[z==1, mean(m)] - dt[z==0, mean(m)] # strength of instrument
  
  iterCE <- data.table(truthA = dt[AStatus == "Always-taker", mean(Y1 - Y0)],
                        truthN = dt[AStatus == "Never-taker", mean(Y1 - Y0)],
                        truthC = dt[AStatus == "Complier", mean(Y1 - Y0)],
                        MM = ITT / ITTm, 
                        IV = coef(ivmodel)[2],
                        PS.a = steps[nrow(steps), AACE],
                        PS.n = steps[nrow(steps), NACE],
                        PS.c = steps[nrow(steps), CACE])
  
  
  estCE <- rbind(estCE, iterCE)
  paramsEst <- rbind(paramsEst, params)
  
}

estCE[, mean((PS.c - 0.8)^2)]
estCE[, mean((IV - .8)^2)]


estCE[, mean(PS.c - .8)]^2 + estCE[, var(PS.c)]
estCE[, mean(IV - .8)]^2 + estCE[, var(IV)]

estCE[, mean(PS.c - .8)]
estCE[, mean(IV - .8)]

estCE[, sd(IV)]
estCE[, sd(PS.c)]
estCE[, sd(truthC)]

dp <-rbind(estCE[,.(method = "EM", trueCACE = truthC, est = PS.c)],
           estCE[,.(method = "IV", trueCACE = truthC, est = IV)])

p2 <- ggplot(data=dp, aes(y = est, x = trueCACE)) +
  geom_point(aes(color=method, alpha = method), size = 1, alpha = .5) +
  geom_hline(yintercept = .8, lty = 3, color = "grey40") +
  geom_vline(xintercept = .8, lty = 3, color = "grey40") +
  scale_color_manual(values = c("#974b0a", "#0a5697")) +
  scale_x_continuous(limits=c(0.4, 1.2), breaks = seq(.4, 1.2, by=.4), name = "Actual CACE of sample") +
  scale_y_continuous(limits=c(0, 2), breaks = seq(0,2, by = .4),name = "Estimated CACE") +
  theme(legend.title = element_blank(),
        legend.position = c(.9, .85),
        panel.grid = element_blank()
  ) +
  ggtitle("Exclusion restriction violated")

png("No_exclusion_restriction.png", width = 400, height=575)
p2
dev.off()


