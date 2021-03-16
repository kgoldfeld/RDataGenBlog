library(mice)
library(simstudy)
library(parallel)

theme_update(panel.grid = element_blank(),
             plot.title = element_text(size = 10, face = "bold"))


d1 <- defData(varname = "u", formula=0, variance = 1, dist = "normal")

d2 <- defDataAdd(varname = "y1", formula = "5 + a*2 + u*2", variance = 2)
d2 <- defDataAdd(d2, "y2", formula = "6 + a*3 + u*2", variance = 2)

dm <- defMiss(varname = "y2", formula = "-3.5 + 0.4*y1", logit.link = TRUE)


#--- 

s_generate <- function(n) {
  dd <- genData(n, d1)
  dd <- trtAssign(dd, grpName = "a")
  dd <- addColumns(d2, dd)
  
  dmiss <- genMiss(dd, dm, id = "id")
  dobs <- genObs(dd, dmiss, id = "id")
  
  return(list(dd, dobs))
}

s_replicate <- function(n) {
  dsets <- s_generate(n)
  est.complete <- coef(lm(y2 ~ a, data = dsets[[1]]))["a"]
  est.obs <- coef(lm(y2 ~ y1 + a, data = dsets[[2]]))["a"]
  
  imp <- mice(dsets[[2]][,-"id"], m=20, maxit=5, print=FALSE)
  
  fit <- with(imp, lm(y2 ~ a))
  pooled.ests <- summary(pool(fit))
  est.impute <- pooled.ests$estimate[2]
  
  diff.complete <- dsets[[1]][, .(avg = mean(y2)), keyby = a][ , avg - shift(avg)][2]    
  diff.obs<- dsets[[2]][!is.na(y2), .(avg = mean(y2)), keyby = a][ , avg - shift(avg)][2] 
  
  return(data.table(est.complete, diff.complete, est.obs, diff.obs, est.impute))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(28172)
results <- rbindlist(mclapply(1:2500, function(x) s_replicate(300), mc.cores = 4))

mu_c <- results[, mean(diff.complete)]
mu_o <- results[, mean(diff.obs)]

p1 <- ggplot(data = results, aes(x = diff.complete, y=diff.obs)) +
  geom_vline(xintercept = mu_c,color = "#b38f00", lty = 3) +
  geom_hline(yintercept = mu_o,color = "#b38f00", lty = 3) +
  geom_point(size = .5) +
  scale_x_continuous(limits = c(.7, 4.2), name = "Complete data") +
  scale_y_continuous(limits = c(.7, 4.2), name = "Observed data") +
  ggtitle("Difference in means")

mu_c <- results[, mean(est.complete)]
mu_o <- results[, mean(est.obs)]

p2 <- ggplot(data = results, aes(x = est.complete, y=est.obs)) +
  geom_vline(xintercept = mu_c,color = "#b38f00", lty = 3) +
  geom_hline(yintercept = mu_o,color = "#b38f00", lty = 3) +
  geom_point(size = .5) +
  scale_x_continuous(limits = c(.7, 4.2), name = "Complete data") +
  scale_y_continuous(limits = c(.7, 4.2), name = "Observed data") +
  ggtitle("Linear regression with adjustment")

mu_c <- results[, mean(est.complete)]
mu_i <- results[, mean(est.impute)]

p3 <- ggplot(data = results, aes(x = est.complete, y=est.impute)) +
  geom_vline(xintercept = mu_c,color = "#b38f00", lty = 3) +
  geom_hline(yintercept = mu_i,color = "#b38f00", lty = 3) +
  geom_point(size = .5) +
  scale_x_continuous(limits = c(1.8, 4.2), name = "Complete data") +
  scale_y_continuous(limits = c(1.8, 4.2), name = "Imputed data") +
  ggtitle("Linear regression with adjustment")

ggsave(
  filename = "content/post/2021-03-30-some-cases-where-imputing-missing-data-matters/img/full.png",
  gridExtra::grid.arrange(p1, p2, nrow = 1), width = 4, height = 2, scale = 1.7)

ggsave(
  filename = "content/post/2021-03-30-some-cases-where-imputing-missing-data-matters/img/impute.png",
  gridExtra::grid.arrange(p3, nrow = 1), width = 2, height = 2, scale = 2)

