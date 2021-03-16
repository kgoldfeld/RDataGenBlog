library(parallel)

theme_update(panel.grid = element_blank(),
             plot.title = element_text(size = 10, face = "bold"))


#---

d1 <- defData(varname = "x", formula=0, variance = 1, dist = "normal")

d2 <- defDataAdd(varname = "y", formula = "6 + 3*a + 1*x", variance = 3)

dm <- defMiss(varname = "y", formula = "-3.2 + 2*x", logit.link = TRUE)

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

results[, .(mean(diff.complete), mean(diff.obs), mean(avg0_c), mean(avg0_o), mean(avg1_c), mean(avg1_o)) ]
results[, .(sd(diff.complete), sd(diff.obs)) ]

mu_c <- results_sim[, mean(diff.complete)]
mu_o <- results_sim[, mean(diff.obs)]

p1 <- ggplot(data = results_sim, aes(x = diff.complete, y=diff.obs)) +
  geom_vline(xintercept = mu_c,color = "#b38f00", lty = 3) +
  geom_hline(yintercept = mu_o,color = "#b38f00", lty = 3) +
  geom_point(size = .5) +
  scale_x_continuous(limits = c(2, 4), breaks = c(2.5, 3, 3.5), name = "Complete data") +
  scale_y_continuous(limits = c(2, 4), breaks = c(2.5, 3, 3.5), name = "Observed data") +
  ggtitle("Effect size")

mu_c <- results_sim[, mean(avg0_c)]
mu_o <- results_sim[, mean(avg0_o)]

p2 <- ggplot(data = results_sim, aes(x = avg0_c, y=avg0_o)) +
  geom_vline(xintercept = mu_c,color = "#b38f00", lty = 3) +
  geom_hline(yintercept = mu_o,color = "#b38f00", lty = 3) +
  geom_point(size = .5) +
  scale_x_continuous(limits = c(5.2, 6.7), breaks = c(5.5, 6, 6.5), name = "Complete data") +
  scale_y_continuous(limits = c(5.2, 6.7), breaks = c(5.5, 6, 6.5), name = "Observed data") +
  ggtitle("Average Y (control arm)")

mu_c <- results_sim[, mean(avg1_c)]
mu_o <- results_sim[, mean(avg1_o)]

p3 <- ggplot(data = results_sim, aes(x = avg1_c, y=avg1_o)) +
  geom_vline(xintercept = mu_c,color = "#b38f00", lty = 3) +
  geom_hline(yintercept = mu_o,color = "#b38f00", lty = 3) +  
  geom_point(size = .5) +
  scale_x_continuous(limits = c(8.2, 9.7), breaks = c(8.5, 9, 9.5), name = "Complete data") +
  scale_y_continuous(limits = c(8.2, 9.7), breaks = c(8.5, 9, 9.5), name = "Observed data") +
  ggtitle("Average Y (treatment arm)")

ggsave(
  filename = "content/post/2021-03-30-some-cases-where-imputing-missing-data-matters/img/comp_obs.png",
  gridExtra::grid.arrange(p2, p3, p1, nrow = 1), width = 6.5, height = 2, scale = 1.5)
