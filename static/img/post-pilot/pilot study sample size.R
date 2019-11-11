estVar <- function(n, sigma2) {
  x <- rnorm(n, 0, sigma2^0.5)
  var(x) 
}

vars <- sapply(1:100000, function(x) estVar(30, 16))
mean(vars)
hist(vars)
median(vars)

defd <- defDataAdd(varname = "y", formula = "rx * 50", variance = 100^2)

dd <- genData(n = 10*2)
dd <- trtAssign(dd, grpName = "rx")
dd <- addColumns(defd, dd)

estsd <- ANOVAreplication::pooled.sd(dd[, .(y, rx)])
estsd


library(pwr) 
nstar <- round(pwr.t.test(n = NULL, d =  50/estsd, sig.level = 0.05, power = 0.80, type = "two.sample")$n,0)

estpower <- function() {
  defx <- defDataAdd(varname = "y", formula = "rx * 50", variance = estsd^2)
  dx <- genData(n = nstar*2)
  dx <- trtAssign(dx, grpName = "rx")
  dx <- addColumns(defx, dx)
  t.test(y~rx, data=dx)$p.value
}

ps <- sapply(1:2500, function(x) estpower())
mean(ps <= 0.05)

pwr.t.test(n=nstar, d = 0.50, sig.level = 0.05, type = "two.sample")

estpower2 <- function() {
  defx <- defDataAdd(varname = "y", formula = "rx * 50", variance = 100^2)
  dx <- genData(n = nstar*2)
  dx <- trtAssign(dx, grpName = "rx")
  dx <- addColumns(defx, dx)
  t.test(y~rx, data=dx)$p.value
}

ps <- sapply(1:2500, function(x) estpower2())
mean(ps <= 0.05)

####

getPower <- function(ssize, esize, gamma = 0, use.est = FALSE) {
  
  estring <- paste0("rx * ", esize)
  defd <- defDataAdd(varname = "y", formula = estring, variance = 100^2)
  
  N <- ssize * 2
  
  dd <- genData(n = N)
  dd <- trtAssign(dd, grpName = "rx")
  dd <- addColumns(defd, dd)
  
  lmfit <- lm(y~rx, data = dd)

  estsd <- ANOVAreplication::pooled.sd(dd[, .(y, rx)])
  ucl <- sqrt( ( (N-2) * estsd^2 ) / qchisq(gamma, df = N - 2, lower.tail = FALSE) )

  pwr.sd <- estsd * (gamma == 0) + ucl * (gamma > 0)
  pwr.eff <- esize * (use.est == FALSE) + coef(lmfit)["rx"] * (use.est == TRUE)
  
  if (abs(pwr.eff/pwr.sd < 0.0002)) pwr.eff <- sign(pwr.eff) * .0002 * pwr.sd
  
  nstar <- round(pwr.t.test(n = NULL, d =  pwr.eff/pwr.sd, sig.level = 0.05, 
                            power = 0.80, type = "two.sample")$n,0)  
  
  
  power <- pwr.t.test(n=nstar, d = esize/100, sig.level = 0.05, type = "two.sample")$power
  
  return(data.table(ssize, esize, gamma, use.est,
    estsd = estsd, ucl = ucl, nstar, power = power,
    est = coef(lmfit)["rx"], 
    lcl.est = confint(lmfit)["rx",1] , 
    ucl.est = confint(lmfit)["rx",2])
  )
  
}

mean(sapply(1:2500, function(x) getPower(ssize = 10, esize = 40)) >= 0.80)
mean(sapply(1:2500, function(x) getPower(ssize = 50, esize = 40)) >= 0.80)
mean(sapply(1:2500, function(x) getPower(ssize = 100, esize = 40)) >= 0.80)
mean(sapply(1:2500, function(x) getPower(ssize = 100, esize = 75)) >= 0.80)


chk <- rbindlist(lapply(1:10000, function(x) getPower(ssize = 20, esize = 40)))
chk <- rbindlist(lapply(1:10000, function(x) getPower(ssize = 20, esize = 40, gamma = 0.95)))
chk <- rbindlist(lapply(1:10000, function(x) getPower(ssize = 20, esize = 40, gamma = 0.50)))
chk <- rbindlist(lapply(1:10000, function(x) getPower(ssize = 40, esize = 40, gamma = 0.50)))
chk <- rbindlist(lapply(1:10000, function(x) getPower(ssize = 40, esize = 75, gamma = 0.50)))
chk <- rbindlist(lapply(1:10000, function(x) getPower(ssize = 10, esize = 75, gamma = 0.50)))
chk <- rbindlist(lapply(1:10000, function(x) getPower(ssize = 10, esize = 75, gamma = 0.60)))
chk <- rbindlist(lapply(1:1000, function(x) getPower(ssize = 30, esize = 40, gamma = 0.60)))
chk <- rbindlist(lapply(1:20000, function(x) getPower(ssize = 100, esize = 75, gamma = 0.90)))


chk <- rbindlist(lapply(1:10000, 
                        function(x) getPower(ssize = 30, esize = 40, gamma = 0.60)))
chk <- rbindlist(lapply(1:10000, 
                        function(x) getPower(ssize = 30, esize = 40, gamma = 0.60, TRUE)))

eff <- 30
ss <- 60

chk[, mean(100 < ucl)]
chk[, mean(power >= 0.80)]

ggplot(data = chk, aes(x=power)) +
  geom_histogram(color = "grey92", fill = "grey50", 
                 breaks = seq(0, 1, .02)) +
  geom_rect(aes(xmin=qs[1], xmax=qs[2], ymin = 0, ymax = Inf), 
            fill = "pink", color = "grey92", alpha = .01) +
  geom_rect(aes(xmin=qs[2], xmax=qs[3], ymin = 0, ymax = Inf), 
            fill = "pink", color = "grey92", alpha = .01) +
  scale_x_continuous(breaks = seq(0, 1, .2)) + 
  scale_y_continuous(limits = c(0, 4100)) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey92"))

### 

set.seed(38273)

  dres <- data.table()
  
  for (i in c(30, 60)) {
    for (j in c(30, 75)) {
      for (k in c(0, .5, .7)) {
        for (l in c(FALSE, TRUE)) {
          dd <- rbindlist(lapply(1:5000, 
             function(x) getPower(ssize = i, esize = j, gamma = k, use.est = l)))
          dres <- rbind(dres, dd)
          cat(i, j, k, l, "\n")
        }
      }
    }
  }
  
  above80 <- dres[, .(x80 = mean(power >= 0.80)), keyby = .(ssize, esize, gamma, use.est)]
  above80[, l80 := scales::percent(x80, accuracy = 1)]
  
  g_labeller <- function(value) {
      paste0("gamma = ", value)
  }
  
  e_labeller <- function(value) {
    paste0("effect size = ", value)
  }
  
  p30 <- ggplot(data = dres[ssize == 30], 
         aes(x=factor(use.est, labels=c("'truth'", "pilot")), y=power)) +
    geom_hline(yintercept = 0.8, color = "white") +
    geom_boxplot(outlier.shape = NA, fill = "#9ba1cf", width = .4) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey92"),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 9, face = "bold")) +
    facet_grid(esize ~ gamma, labeller = labeller(gamma = g_labeller, esize = e_labeller)) +
    scale_x_discrete(name = "\n source of effect size used for power calculation") +
    scale_y_continuous(limits = c(0,1), breaks = c(0, .8),
                       name = "distribution of power estimates \n") +
    ggtitle("Distribution of power estimates (n = 30 per treatment arm)") +
    geom_text(data = above80[ssize == 30], 
              aes(label = l80), x=rep(c(0.63, 1.59), 6), y = 0.95,
              size = 2.5)
  
  p30
  
  p60 <- ggplot(data = dres[ssize == 60], 
                aes(x=factor(use.est, labels=c("'truth'", "pilot")), y=power)) +
    geom_hline(yintercept = 0.8, color = "white") +
    geom_boxplot(outlier.shape = NA, fill = "#9ba1cf", width = .4) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey92"),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 9, face = "bold")) +
    facet_grid(esize ~ gamma, labeller = labeller(gamma = g_labeller, esize = e_labeller)) +
    scale_x_discrete(name = "\n source of effect size used for power calculation") +
    scale_y_continuous(limits = c(0,1), breaks = c(0, .8),
                       name = "distribution of power estimates \n") +
    ggtitle("Distribution of power estimates (n = 60 per treatment arm)") +
    geom_text(data = above80[ssize == 60], 
              aes(label = l80), x=rep(c(0.63, 1.59), 6), y = 0.95,
              size = 2.5)
  
  p60

ggsave(filename =  "static/img/post-pilot/pilot60.png", p60, width = 6, height = 6)
ggsave(filename =  "static/img/post-pilot/pilot30.png", p30, width = 6, height = 6)

####

genRx <- function() {
  defx <- defDataAdd(varname = "y", formula = "rx * 40", variance = 100^2)
  dx <- genData(n = 30 * 2)
  dx <- trtAssign(dx, grpName = "rx")
  dx <- addColumns(defx, dx)
  lmfit <- lm(y~rx, data = dx)
  data.table(est = coef(lmfit)["rx"], 
             lcl = confint(lmfit)["rx",1] , 
             ucl = confint(lmfit)["rx",2])
}

estRes <- rbindlist(lapply(1:1000, function(x) genRx()))
estRes[, cover := (40 > lcl & 40 < ucl)]
estRes[, mean(cover)]

t.test(y~rx, data=dx)$p.value


