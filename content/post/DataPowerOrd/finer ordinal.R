library(ordinal)
library(simstudy)

genDT <- function(nobs, baseprobs, defA) {
  
  dT <- genData(nobs)
  dT <- trtAssign(dT, grpName = "exposed")
  dT <- addColumns(defA, dT)
  
  dT <- genOrdCat(dT, adjVar = "z", baseprobs, catVar = "r")
  dT[]
}

genDTnon <- function(nobs) {
  
  # define the data
  
  ps <- round(gtools::rdirichlet(2, rep(2, 8)), 2)
  
  p0 <- paste(ps[1, 1:7], collapse = ";")
  p1 <- paste(ps[2, 1:7], collapse = ";")
  
  defc <- defCondition(condition = "exposed == 0", formula = p0, dist = "categorical")
  defc <- defCondition(defc, condition = "exposed == 1", formula = p1, dist = "categorical")
  
  # generate the data
  
  dx <- genData(nobs)
  dx <- trtAssign(dx, grpName = "exposed")
  dx <- addCondition(defc, dx, "r")
  
  dx[, r := factor(r)]

  dx[]
}

analyzeDT <- function(nobs, baseprobs, defA) {

  dT <- genDT(nobs, baseprobs, defA)
  clmFit <- clm(r ~ exposed, data = dT)
  coef(summary(clmFit))
  
}


defA <- defDataAdd(varname = "z", formula = "-1.0 * exposed", dist = "nonrandom")

baseprobs <- list(
  c(0.0375*8, 0.0875*8),               
  c(rep(0.0375*8, 1), rep(0.0875*4, 2)), 
  c(rep(0.0375*4, 2), rep(0.0875*4, 2)),  
  c(rep(0.0375*4, 2), rep(0.0875*4, 1), rep(0.0875*2, 2)),
  c(rep(0.0375*2, 2), rep(0.0375*4, 1), rep(0.0875*4, 1), rep(0.0875*2, 2)),
  c(rep(0.0375*2, 2), rep(0.0375*4, 1), rep(0.0875*2, 4)),
  c(rep(0.0375*2, 4), rep(0.0875*2, 4)),
  c(rep(0.0375*2, 4), rep(0.0875*2, 3), rep(0.0875, 2)),
  c(rep(0.0375, 2), rep(0.0375*2, 3), rep(0.0875*2, 3), rep(0.0875, 2)),
  c(rep(0.0375, 2), rep(0.0375*2, 3), rep(0.0875*2, 2), rep(0.0875, 4)),
  c(rep(0.0375, 4), rep(0.0375*2, 2), rep(0.0875*2, 2), rep(0.0875, 4)),
  c(rep(0.0375, 4), rep(0.0375*2, 2), rep(0.0875*2, 1), rep(0.0875, 6)),
  c(rep(0.0375, 6), rep(0.0375*2, 1), rep(0.0875*2, 1), rep(0.0875, 6)),
  c(rep(0.0375, 6), rep(0.0375*2, 1), rep(0.0875, 8)),
  rep(c(.0375, 0.0875), each = 8)
)

#####

getRecCoords <- function(x) {
  
  lenx <- length(x)
  n1 <- floor(lenx/2)
  n2 <- ceiling(lenx/2)
  
  data.table(xstart = c(0, cumsum(x)[-lenx]), 
             xend = cumsum(x), 
             ymin = lenx - 2 + 0.3, 
             ymax = lenx - 1 - 0.3,
             grp = factor(rep(c(1, 2), times = c(n1, n2)))
  )
}
  
dx <- rbindlist(lapply(baseprobs, function(x) getRecCoords(x)))

ggplot(data = dx, aes(xmin= xstart, xmax = xend, ymin = ymin, ymax = ymax, 
                      fill = grp)) +
  geom_rect(color = "grey45", size = 1) +
  scale_fill_manual(values = c("#deb0ca", "#cadeb0")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "white"))

#####

iterate <- function(niter, nobs, baseprobs, defA) {
  res <- lapply(1:niter, function(x) analyzeDT(nobs, baseprobs, defA))
  mean( sapply(res, function(x) x["exposed",  "Pr(>|z|)"]) < 0.05) 
}

set.seed(1295)
power <- lapply(baseprobs, function(x) iterate(niter = 10000, nobs = 100, x, defA))
save(power, file = "Working/powerord.Rdata")

dp <- data.table(cats = unlist(lapply(baseprobs, function(x) length(x))),
           power = unlist(power)
)

ggplot(data = dp, aes(x = cats, y = power)) +
  geom_smooth(se = FALSE, size = 1) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = dp$cats, name = "number of categories") +
  scale_y_continuous(limits = c(.65, .85), breaks = c(.7, .75, .8), 
                     name = "estimated power") +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank())

