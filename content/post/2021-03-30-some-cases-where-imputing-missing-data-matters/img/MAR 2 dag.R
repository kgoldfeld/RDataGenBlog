library(parallel)

theme_update(panel.grid = element_blank(),
             plot.title = element_text(size = 10, face = "bold"))


#---

d1 <- defData(varname = "x", formula=0.5, dist = "binary")
d2 <- defDataAdd(varname = "y", formula = "5 + 1*a + 5*x + 3*a*x", variance = 2)
dm <- defMiss(varname = "y", formula = "-3.5 + 2.3*x", logit.link = TRUE)

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
  
  diff.complete <- dsets[[1]][, .(avg = mean(y)), keyby = a][ , avg - shift(avg)][2]
  diff.obs <- dsets[[2]][!is.na(y), .(avg = mean(y)), keyby = a][ , avg - shift(avg)][2]
  
  coefs <- coef(lm(y ~ a*x, data = dsets[[2]]))
  est.obs_a <- coefs["a"]
  est.obs_ax <- sum(coefs[c("a", "a:x")])
  
  coefs <- coef(lm(y ~ a*x, data = dsets[[1]]))
  est.comp_a <- coefs["a"]
  est.comp_ax <- sum(coefs[c("a", "a:x")])
  
  results <- data.table(diff.complete, diff.obs, 
                        est.comp_a, est.obs_a, 
                        est.comp_ax, est.obs_ax)
  return(results)
}

RNGkind("L'Ecuyer-CMRG")
set.seed(18271)

results_sim <- rbindlist(mclapply(1:2500, function(x) s_replicate(300), mc.cores = 4))

results_sim[, .(mean(diff.complete), mean(diff.obs), mean(est.obs_a), mean(est.obs_ax)) ]
results_sim[, .(sd(diff.complete), sd(diff.obs), sd(est.obs)) ]

mu_c <- results_sim[, mean(diff.complete)]
mu_o <- results_sim[, mean(diff.obs)]

p1 <- ggplot(data = results_sim, aes(x = diff.complete, y=diff.obs)) +
  geom_vline(xintercept = mu_c,color = "#b38f00", lty = 3) +
  geom_hline(yintercept = mu_o,color = "#b38f00", lty = 3) +
  geom_point(size = .5) +
  scale_x_continuous(limits = c(1.2, 3.7), breaks = seq(1.5, 3.5, by = .5), name = "Complete data") +
  scale_y_continuous(limits = c(1.2, 3.7), breaks = seq(1.5, 3.5, by = .5), name = "Observed data") +
  ggtitle("Difference in means: unadjusted")

mu_c <- results_sim[, mean(est.comp_a)]
mu_o <- results_sim[, mean(est.obs_a)]

p2 <- ggplot(data = results_sim, aes(x = est.comp_a, y=est.obs_a)) +
  geom_vline(xintercept = mu_c,color = "#b38f00", lty = 3) +
  geom_hline(yintercept = mu_o,color = "#b38f00", lty = 3) +
  geom_point(size = .5) +
  scale_x_continuous(limits = c(-.2, 2.2), breaks = seq(0, 2, by = .5), name = "Complete data") +
  scale_y_continuous(limits = c(-.2, 2.2), breaks = seq(0, 2, by = .5), name = "Observed data") +
  ggtitle("Difference in means: x = 0")

mu_c <- results_sim[, mean(est.comp_ax)]
mu_o <- results_sim[, mean(est.obs_ax)]

p3 <- ggplot(data = results_sim, aes(x = est.comp_ax, y=est.obs_ax)) +
  geom_vline(xintercept = mu_c,color = "#b38f00", lty = 3) +
  geom_hline(yintercept = mu_o,color = "#b38f00", lty = 3) +  
  geom_point(size = .5) +
  scale_x_continuous(limits = c(2.8, 5.2), breaks =  seq(3, 5, by = .5), name = "Complete data") +
  scale_y_continuous(limits = c(2.8, 5.2), breaks =  seq(3, 5, by = .5), name = "Observed data") +
  ggtitle("Difference in means: x = 1")

ggsave(
  filename = "content/post/2021-03-30-some-cases-where-imputing-missing-data-matters/img/MAR_2_diff.png",
  p1, width = 2, height = 2, scale = 1.7)

ggsave(
  filename = "content/post/2021-03-30-some-cases-where-imputing-missing-data-matters/img/MAR_2_adj.png",
  gridExtra::grid.arrange(p2, p3, nrow = 1), width = 4, height = 2, scale = 1.5)
