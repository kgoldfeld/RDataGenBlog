library(pscl)
library(MASS)
library(parallel)

defs <- function(hint, heffect, pint, peffect, nRes) {
  
  form.xHurdle <- gsub("hint", hint, "hint + heffect * rx")
  form.xHurdle <- gsub("heffect", heffect, form.xHurdle)
  
  form.xNon0Pois = gsub("pint", pint, "pint + peffect * rx + log(patDays)")
  form.xNon0Pois = gsub("peffect", peffect, form.xNon0Pois)
  
  def <- defDataAdd(varname = "nRes", formula = nRes, dist = "poisson")
  def <- defDataAdd(def, varname = "pDays", formula = 80, dist = "poisson")
  def <- defDataAdd(def, varname = "nDays", formula = "pmin(90, pDays)", 
                    dist = "nonrandom")
  def <- defDataAdd(def, varname = "patDays", formula = "nRes * nDays",
                    dist = "nonrandom")
  def <- defDataAdd(def, varname = "xHurdle", formula = form.xHurdle, 
                    dist = "binary")
  def <- defDataAdd(def, varname = "xNon0Pois", formula = form.xNon0Pois, 
                    dist = "noZeroPoisson", link = "log")
  def <- defDataAdd(def, varname = "y", formula = "xHurdle * xNon0Pois", 
                    dist = "nonrandom")
  
  return(def[])
}

gData <- function(n, def) {
  
  dx <- genData(n)
  dx <- trtAssign(dx, grpName = "rx")
  dx <- addColumns(def, dx)
  
  dx[]
  
}

estModel <- function(dx) {
  
  hfit <- hurdle(y ~ rx | rx, offset = log(patDays), data = dx, )
  hfit0 <- hurdle(y ~ 1 | 1, offset = log(patDays), data = dx)
  lrt <- -2*(logLik(hfit0) - logLik(hfit))

  data.table(p.zero = coef(summary(hfit))$zero["rx", "Pr(>|z|)"],
             p.count = coef(summary(hfit))$count["rx", "Pr(>|z|)"],
             X2 = 1 - pchisq(lrt, 2))
  
}

iter <- function(n, def) {
  
  dx <- gData(n, def)
  hfit <- estModel(dx)
  return(data.table(n = n, hfit))
  
}

iterp <- function(x) {
  
  def <- defs(x["hint"], x["heff"], x["pint"], x["peff"], x["nRes"])
  n <- x["nSites"]
  
  res <- rbindlist(mclapply(1:2500, function(x) iter(n, def), mc.cores = 4))
  power <- res[, mean(X2 < 0.05)]
  
  data.table(hint = x["hint"], heff = x["heff"], 
             pint = x["pint"], peff = x["peff"], 
             nRes = x["nRes"], nSites = n,
             power)
}

hint <- c(0.7)
heffect <- c(0.00, -0.10, -0.15)
pint <- c(log(10/4000), log(15/4000), log(20/4000), log(25/4000))
peffect <- c(log(0.90), log(0.85), log(0.80), log(0.75), log(0.70))
nRes <- c(50)
nSites <- c(100, 150, 200)

params <- expand.grid(hint = hint, heff = heffect, pint = pint, peff = peffect, 
                      nRes = nRes, nSites = nSites)
lparam <- asplit(params, 1)

set.seed(2872)
res.power.l <- lapply(lparam, function(x) iterp(x))

res.power <- rbindlist(res.power.l)

res.power[, infect := exp(pint)*1000]
res.power[, infect := factor(formatC(infect, 2, width = 4, format = "f"))]
res.power[, infect := factor(infect, levels = rev(levels(infect )))]

res.power[, rxinfect :=  .7 + heff]
res.power[, `P(# infections in rx > 0)` := factor(formatC(rxinfect, 2, width = 4, format = "f"))]
res.power[, `total # NHs` := factor(nSites)] 

save(res.power, file = "Data/respower.rda")


library(paletteer)
p1 <- ggplot(data = res.power, aes(x = exp(peff), y = power, group = pint)) +
  geom_hline(yintercept = 0.8, color = "grey50", lty = 3) +
  geom_line(aes(color = infect)) +
  facet_grid(`total # NHs` ~ `P(# infections in rx > 0)`, 
             labeller = label_both) +
  theme(panel.grid = element_blank()) +
  scale_x_reverse(name = "\nrate ratio",) +
  scale_color_paletteer_d(name = "control rate:\n# infections/1K resident-days", 
                          palette = "wesanderson::BottleRocket2") +
  scale_y_continuous(limits = c(.2, 1), breaks = c(.4, .6, .8, 1)) +
  theme(plot.caption = element_text(size = 10, hjust = 0, lineheight = .5),
        plot.title = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 8)) +
  labs(title = "Power analysis")

ggsave(file = "~/Local R Projects/RDataGenBlog/static/img/post-hurdle/power.png", 
       height = 5, width = 7, scale = 1.2)

