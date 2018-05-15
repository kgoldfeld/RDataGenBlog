dt <- genData(75)
dt <- trtAssign(dt, 2, grpName = "t")

def <- defDataAdd(varname = "y", formula = "0 + 0*t", variance = 1, dist = "normal")

results <- data.table()

for (i in 1:200) {
  
  dr <- addColumns(def, dt)
  
  lmfit <- lm(y~t, data = dr)
  
  halfWidth = qt(0.975, df = lmfit$df)
  se.t <- sqrt(diag(vcov(lmfit)))["t"]
  pointest.t <- coef(lmfit)["t"]
  
  ci.t <- pointest.t + c(-1, 1) * halfWidth * se.t
  
  ests <- data.table(i, est = pointest.t, l95 = ci.t[1], u95 = ci.t[2])
  results <- rbind(results, ests)
  
}

results[, inCI := (l95 <= 0 & u95 >= 0)]
results[, mean(inCI)]

ggplot(data = results, aes(x=i, y=est)) +
  geom_point(size = .5) +
  geom_segment(aes(x = i, xend = i, y = l95, yend=u95, color = factor(inCI))) +
  scale_color_manual(values = c("red", "grey60")) +
  geom_hline(yintercept = 0, color = "grey50") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() 


### Binary

dt <- genData(200)
dt <- trtAssign(dt, 2, grpName = "t")

def <- defDataAdd(varname = "y", formula = "0.5 + 0*t", 
                  dist = "binary", link = "logit")

results <- data.table()
resexp <- data.table()
resexp2 <- data.table()

for (i in 1:50000) {
  
  dr <- addColumns(def, dt)
  
  glmfit <- glm(y~t, family = binomial, data = dr)
  
  se.t <- sqrt(diag(vcov(glmfit)))["t"]
  pointest.t <- coef(glmfit)["t"]
  
  ci.t <- pointest.t + c(-1, 1) * 1.96 * se.t
  p.t <- coef(summary(glmfit))["t", "Pr(>|z|)"]
  
  ests <- data.table(i, est = pointest.t, l95 = ci.t[1], u95 = ci.t[2], p.t)
  results <- rbind(results, ests)
  
  # Wrong
  
  pointest.or <- exp(pointest.t)
  se.or <- exp(se.t)
  ci.or <- pointest.or + c(-1, 1) * 1.96 * se.or
  
  ests.or <- data.table(i, est = pointest.or, l95 = ci.or[1], 
                        u95 = ci.or[2])
  
  resexp <- rbind(resexp, ests.or)
  
  # Better?
  
  pointest.or <- exp(pointest.t)
  ci.or <- exp(ci.t)  
  
  ests.or <- data.table(i, est = pointest.or, l95 = ci.or[1], 
                        u95 = ci.or[2])
  
  resexp2 <- rbind(resexp2, ests.or)
  
}

results[, inCI := (l95 <= 0 & u95 >= 0)]
results[, mean(inCI)]
results[, mean(p.t <= .05)]

ggplot(data = results, aes(x=i, y=est)) +
  geom_point(size = .5) +
  geom_segment(aes(x = i, xend = i, y = l95, yend=u95, color = factor(inCI))) +
  scale_color_manual(values = c("red", "grey60")) +
  geom_hline(yintercept = 0, color = "grey50") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() 

resexp[, inCI := (l95 <= 1 & u95 >= 1)]
resexp[, mean(inCI)]

ggplot(data = resexp, aes(x=i, y=est)) +
  geom_point(size = .5) +
  geom_segment(aes(x = i, xend = i, y = l95, yend=u95, color = factor(inCI))) +
  scale_color_manual(values = c("red", "grey60")) +
  geom_hline(yintercept = 1, color = "grey50") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() 

resexp2[, inCI := (l95 <= 1 & u95 >= 1)]
resexp2[, mean(inCI)]

ggplot(data = resexp2, aes(x=i, y=est)) +
  geom_point(size = .5) +
  geom_segment(aes(x = i, xend = i, y = l95, yend=u95, color = factor(inCI))) +
  scale_color_manual(values = c("red", "grey60")) +
  geom_hline(yintercept = 1, color = "grey50") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() 

### Bootstrap

bsresults <- data.table()

for (j in 1: 250) {
  
  dr <- addColumns(def, dt)
  
  glmfit <- glm(y~t, family = binomial, data = dr)
  
  se.t <- sqrt(diag(vcov(glmfit)))["t"]
  pointest.t <- coef(glmfit)["t"]
  
  ci.t <- pointest.t + c(-1, 1) * 1.96 * se.t
  ci.or <- exp(ci.t)
  
  resultb <- data.table()
  
  for (i in 1:999) {
    
    ids <- dr[, .(ids = sample(id, .N, replace = TRUE)), keyby = t]
    bs <- dr[ids[,ids]]
    
    glmbs <- glm(y~t, family = binomial, data = bs)
    
    se.t <- sqrt(diag(vcov(glmbs)))["t"]
    pointest.t <- coef(glmbs)["t"]
    pointest.or <- exp(pointest.t)
    
    estsb <- data.table(i, est.or = pointest.or, est.lor = pointest.t)
    resultb <- rbind(resultb, estsb)
    
  }
  
  bsci.or = round(quantile(x = resultb$est.or, probs = c(0.025, 0.975)), 3)
  ci.or = round(ci.or, 3)
  
  bsci.lor = round(quantile(x = resultb$est.lor, probs = c(0.025, 0.975)), 3)
  ci.lor = round(ci.t, 3)
  
  dd <- data.table(l95bsor = bsci.or[1], u95bsor = bsci.or[2],
             l95or = ci.or[1], u95or = ci.or[2],
             l95bslor = bsci.lor[1], u95bslor = bsci.lor[2],
             l95lor = ci.lor[1], u95lor = ci.lor[2])
  
  bsresults <- rbind(bsresults, dd)
  
}


bsresults[, mean(l95bsor < 1 & u95bsor > 1)]
bsresults[, mean(l95or < 1 & u95or > 1)]

bsresults[, mean(l95bslor < 0 & u95bslor > 0)]
bsresults[, mean(l95lor < 0 & u95lor > 0)]

