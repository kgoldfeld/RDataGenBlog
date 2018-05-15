
library(simstudy)
library(ordinal)


#### Set definitions

def1 <- defData(varname = "male", formula = 0.45, dist = "binary", id = "idG")
def1 <- defData(def1, varname="x", formula="-0.5;0.5", dist = "uniform")

defZ <- defDataAdd(varname = "z1", formula = "2*x", dist = "nonrandom")
defZ <- defDataAdd(defZ, varname = "z2", formula = "2*male", dist = "nonrandom")


#### Generate data

set.seed(20)

dx <- genData(10000, def1)
dx <- addColumns(defZ, dx)

probs<-c(0.40, 0.25, 0.15, 0.20)
dx <- genOrdCat(dx, adjVar = "z1", probs, catVar = "grp1")


m1 <- clm(grp1 ~ male + x, data = dx, link = "logit")
summary(m1)

log(cumsum(probs)/(1-cumsum(probs)))

cgrp <- cut(dx$x, breaks = seq(-.55, .55, .1), include.lowest = TRUE, labels = FALSE)
dx[ , cgrp := cgrp]

dtable <- dx[, table(grp1), keyby = cgrp]
dtable[, p := V1/sum(V1), keyby = cgrp]

library(ggplot2)

dx[, jcat:= jitter(as.numeric(grp1), factor = .5)]

ggplot(data = dx, aes(x = x, y = jcat, group = male)) +
  geom_point(color = "forestgreen") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggplot(data = dx, aes(x = x, y = jcat, group = male)) +
  geom_point(color = "forestgreen") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

