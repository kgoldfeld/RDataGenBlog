library(simstudy)
library(mgcv)
library(gamm4)

# test 1

def <- defData(varname = "a", formula = 0, variance = 0.6)
def <- defData(def, varname = "mu_b", formula = 0, dist = "nonrandom")
def <- defData(def, varname = "s2_b", formula = 0.1, dist = "nonrandom")

defOut <- defDataAdd(varname = "x", formula = .4, dist = "binary")
defOut <- defDataAdd(defOut, varname = "y", 
                     formula = "-1.5 + a + b + 0.4 * x + 0.65 * A", 
                     dist = "binary", link="logit"
)

ds <- genData(30, def, id = "site")
ds <- addPeriods(ds, 250, "site", perName = "k")
ds <- addCorGen(
  dtOld = ds, idvar = "site", 
  rho = 0.99, corstr = "ar1",
  dist = "normal", param1 = "mu_b", param2 = "s2_b", cnames = "b"
)

ds <- trtStepWedge(ds, "site", nWaves = 30, 
                   lenWaves = 3, startPer = 100, 
                   grpName = "A", perName = "k"
)

ds$site <- as.factor(ds$site)

dd <- genCluster(ds, "timeID", numIndsVar = 80, level1ID = "id")
dd <- addColumns(defOut, dd)

system.time(fit.A <- bam(
  y ~ x + A + s(k) + s(k, site, bs = "fs"), 
  data = dd, 
  method = "fREML",
  family = "binomial")
)

summary(fit.A)

# test 2

d <- defData(varname = "b", formula = 0, variance = .15)
d <- defData(d, varname = "n_meas", formula = 1, dist = "noZeroPoisson")
d <- defData(d, varname = "x", formula = 0.4, dist = "binary")

d2 <- defDataAdd(varname = "y", formula = "-1.5 + 0.5*x + b", dist = "binary", link="logit")

dd <- genData(10000, d)
dd <- genCluster(dd,cLevelVar = "id", numIndsVar = "n_meas", level1ID = "measure")
dd <- addColumns(d2, dd)

dd$id <- factor(dd$id)

system.time(fit <- bam(
  y ~ x + s(id, bs = "re"), 
  data = dd, 
  method = "fREML",
  family = "binomial")
)

summary(fit)

fit <- lme4::glmer(y ~ x + (1 | id), data = dd, family = binomial)
summary(fit)

system.time(fit <- gamm4::gamm4(
  y ~ x + s(id, bs = "re"), 
  data = dd, 
  method = "fREML",
  family = "binomial")
)

# test 3

def <- defData(varname = "a", formula = 0, variance = 0.6)
def <- defData(def, varname = "A", formula = "1;1", dist = "trtAssign")

defOut <- defDataAdd(varname = "y", 
                     formula = "-1.5 + a + 0.65 * A", 
                     dist = "binary", link="logit"
)

ds <- genData(30, def, id = "site")
dd <- genCluster(ds, "site", 100000, "id")
dd <- addColumns(defOut, dd)

dd$site <- factor(dd$site)

system.time(fit <- bam(
  y ~ A + s(site, bs = "re"), 
  data = dd, 
  method = "fREML",
  family = "binomial")
)

summary(fit)

fit <- lme4::glmer(y ~ A + (1 | site), data = dd, family = binomial)
summary(fit)
