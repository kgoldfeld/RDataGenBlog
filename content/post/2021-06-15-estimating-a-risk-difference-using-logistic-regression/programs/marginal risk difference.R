library(simstudy)
library(data.table)
library(parallel)
library(targeted)

# def <- defData(varname = "x1", formula = 0, variance = 1)
# def <- defData(def, "x2", formula = 0.3, dist="binary")
# def <- defData(def, "rx", formula = .5, dist = "binary")
# def <- defData(def, "y", formula = "-0.5 + 1*rx + 1*x1 - 1*x2", 
#                dist = "binary", link = "logit")

def <- defDataAdd(varname = "x1", formula = "..mu_x", variance = 2, dist = "beta")
def <- defDataAdd(def, varname = "x", formula = "x1 - 1", dist = "nonrandom")
def <- defDataAdd(def, varname = "z1", formula = "..mu_z", variance = 9, dist = "beta")
def <- defDataAdd(def, varname = "z", formula = "z1 - 1", dist = "nonrandom")
def <- defDataAdd(def, "y", formula = "-2 + 1.2 * rx + .5 * x + .5 * z",
              dist = "binary", link="logit")

mu_x <- .5
mu_z <- .8

mu_x <- .1
mu_z <- .2

dd <- genData(300)
dd <- trtAssign(dd, grpName = "rx")
dd <- addColumns(def, dd)


hist(dd$x)
hist(dd$z)

dd[, mean(y), keyby =.(rx)]
dd[rx == 1, mean(y)] - dd[rx == 0, mean(y)]
log(dd[rx == 1, mean(y)] / dd[rx == 0, mean(y)])

glmfit <- glm(y ~ rx + x + z, data = dd, family = "binomial")
summary(glmfit)

newdata <- dd[, .(rx=1, x, z)]
p1 <- mean(predict(glmfit, newdata, type = "response"))

newdata <- dd[, .(rx=0, x, z)]
p0 <- mean(predict(glmfit, newdata, type = "response"))

p1 - p0
log(p1/p0)

bootdif <- function(dd) {
  
  db <- dd[, .(id = sample(id, replace = TRUE)), 
              keyby = rx]
  
  db <- merge(db[, id, rx], dd, by = c("id", "rx"))
  
  bootfit <- glm(y ~ rx + x + z, data = db, family = "binomial")
  
  newdata <- db[, .(rx=1, x,  z)]
  p1 <- mean(predict(bootfit, newdata, type = "response"))
  newdata <- db[, .(rx=0, x, z)]
  p0 <- mean(predict(bootfit, newdata, type = "response"))
  
  data.table(rd = p1 - p0, rr = p1/p0)
}

bootest <- rbindlist(mclapply(1:999, function(x) bootdif(dd), mc.cores = 4))
bootest[, quantile(rd, c(0.025, 0.975))]
bootest[, log(quantile(rr, c(0.025, 0.975)))]

x1 <- model.matrix(~x1 + x2, dd)
x2 <- model.matrix(~x1 + x2, dd)

fit1 <- with(dd, riskreg_mle(y, rx, x1, x2, type="rd"))
fit1

fit2 <- with(dd, riskreg_mle(y, rx, x1, x2, type="rr"))
fit2



dd <- genData(10000, d)



