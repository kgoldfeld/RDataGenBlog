library(parallel)

theme_update(panel.grid = element_blank(),
             plot.title = element_text(size = 10, face = "bold"))


#---

d1 <- defData(varname = "x", formula=0.5, dist = "binary")
d2 <- defDataAdd(varname = "y", formula = "5 + 2.5*a + 5*x", variance = 2)
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
  
  avg.complete <- dsets[[1]][, .(avg = mean(y)), keyby = a]$avg
  avg.obs <- dsets[[2]][, .(avg = mean(y, na.rm = TRUE)), keyby = a]$avg
  
  results <- data.table(diff.complete, diff.obs, t(avg.complete), t(avg.obs))
  setnames(results, c("diff.complete", "diff.obs", "avg0_c", "avg1_c","avg0_o", "avg1_o"))
  return(results)
}

RNGkind("L'Ecuyer-CMRG")
set.seed(18271)

results_sim <- rbindlist(mclapply(1:2500, function(x) s_replicate(300), mc.cores = 4))

mu_c <- results_sim[, mean(diff.complete)]
mu_o <- results_sim[, mean(diff.obs)]

p1 <- ggplot(data = results_sim, aes(x = diff.complete, y=diff.obs)) +
  geom_vline(xintercept = mu_c,color = "#b38f00", lty = 3) +
  geom_hline(yintercept = mu_o,color = "#b38f00", lty = 3) +
  geom_point(size = .5) +
  scale_x_continuous(limits = c(1.2, 3.7), breaks = seq(1.5, 3.5, by = .5), name = "Complete data") +
  scale_y_continuous(limits = c(1.2, 3.7), breaks = seq(1.5, 3.5, by = .5), name = "Observed data") +
  ggtitle("Difference in means: unadjusted")

mu_c <- results_sim[, mean(avg0_c)]
mu_o <- results_sim[, mean(avg0_o)]

p2 <- ggplot(data = results_sim, aes(x = avg0_c, y=avg0_o)) +
  geom_vline(xintercept = mu_c,color = "#b38f00", lty = 3) +
  geom_hline(yintercept = mu_o,color = "#b38f00", lty = 3) +
  geom_point(size = .5) +
  scale_x_continuous(limits = c(6.2, 8.4), breaks = seq(6.5, 8, by = .5), name = "Complete data") +
  scale_y_continuous(limits = c(6.2, 8.4), breaks = seq(6.5, 8, by = .5), name = "Observed data") +
  ggtitle("Average Y (control arm)")

mu_c <- results_sim[, mean(avg1_c)]
mu_o <- results_sim[, mean(avg1_o)]

p3 <- ggplot(data = results_sim, aes(x = avg1_c, y=avg1_o)) +
  geom_vline(xintercept = mu_c,color = "#b38f00", lty = 3) +
  geom_hline(yintercept = mu_o,color = "#b38f00", lty = 3) +  
  geom_point(size = .5) +
  scale_x_continuous(limits = c(8.8, 10.9), breaks = seq(9, 10.5, .5), name = "Complete data") +
  scale_y_continuous(limits = c(8.8, 10.9), breaks = seq(9, 10.5, .5), name = "Observed data") +
  ggtitle("Average Y (treatment arm)")

ggsave(
  filename = "content/post/2021-03-30-some-cases-where-imputing-missing-data-matters/img/MAR_1_diff.png",
  p1, width = 2, height = 2, scale = 2)

ggsave(
  filename = "content/post/2021-03-30-some-cases-where-imputing-missing-data-matters/img/MAR_1_y.png",
  gridExtra::grid.arrange(p2, p3, nrow = 1), width = 4, height = 2, scale = 1.5)
