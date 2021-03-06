rbinom(3,c(100, 100, 100) ,c(.5, .4, .3))

genBinomial <- function(p, n) {
  rbinom(n= 1, size = n, p = p)
}

# genCat - specify label names as 2nd parameter (variance)

mapData <- function(dt, map, asFactor = FALSE, replace = FALSE) {
  
  # need to check if map column exists,
  # make sure all possible maps are specified?
  
  dt2 <- copy(dt)
  origkey <- key(dt)
  
  newcol <- names(map)[2]
  
  varCat <- names(map)[1]
  
  setkeyv(dt2, varCat)
  setkeyv(map, varCat)
  
  dt2 <- map[dt2]
  setkeyv(dt2, origkey)
  
  setcolorder(dt2, c(names(dt), newcol))
  
  if (asFactor) {
    dt2[, (newcol) := as.factor(get(newcol))]
  }
  
  if (replace) {
    dt2[, (varCat) := get(newcol)]
    dt2[, (newcol) := NULL]
  }
  
  return(dt2)
}

def <- defData(varname = "p", formula = .5, dist = "nonrandom")
def <- defData(def,
               varname = "size", 
               formula = "0.35;0.20;0.25", 
               dist="categorical")

map <- data.table(size = c(1,2,3,4), 
                  fsize = c(100, 200, 300, 400))

# generate data

dt <- genData(1000, def)
dt <- mapData(dt, map, asFactor = FALSE, replace = TRUE)

dt[, xbin := genBinomial(p, size), keyby = id]
dt[, xpois := rpois(1, p*size), keyby = id]

dt[, .(avgbin = mean(xbin), varbin = var(xbin),
       avgpois = mean(xpois), varpois = var(xpois)), keyby = size]

dm <- melt(dt, id.vars = c("id", "size"), measure.vars = c("xbin", "xpois"))

ggplot(data = dm, aes(x = size, y = value)) +
  geom_jitter(size = 0.5, width = 25) +
  ylim(-5,250) +
  theme(panel.grid.minor = element_blank()) +
  facet_grid(~variable)

# true model is binomial

binfit <- glm(cbind(xbin, size-xbin) ~ 1, family = binomial, data = dt)
poisfit <- glm(xbin ~ 1, family = poisson, offset = log(size), data = dt)

dt[, simbin := simulate(binfit, nsim = 1 )[[1]][,1]]
dt[, simpois := simulate(poisfit, nsim = 1)]

dmsim <- melt(dt, id.vars = c("id", "size"), 
           measure.vars = c("xbin", "simbin", "simpois"))


ggplot(data = dmsim, aes(x = size, y = value)) +
  geom_jitter(size = 0.5, width = 25) +
  ylim(-5,250) +
  theme(panel.grid.minor = element_blank()) +
  facet_grid(~variable)

# true model is Poisson

binfit <- glm(cbind(xpois, size-xpois) ~ 1, family = binomial, data = dt)
poisfit <- glm(xpois ~ 1, family = poisson, offset = log(size), data = dt)

dt[, simbin := simulate(binfit, nsim = 1 )[[1]][,1]]
dt[, simpois := simulate(poisfit, nsim = 1)]

dmsim <- melt(dt, id.vars = c("id", "size"), 
              measure.vars = c("xpois", "simbin", "simpois"))


ggplot(data = dmsim, aes(x = size, y = value)) +
  geom_jitter(size = 0.5, width = 25) +
  ylim(-5,250) +
  theme(panel.grid.minor = element_blank()) +
  facet_grid(~variable)

