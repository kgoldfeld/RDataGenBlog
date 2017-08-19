library(geepack)
genFromU <- function(U, dist, param) {

  if (dist == "binary") colVal <- qbinom(U, 1, param[1])
  else if (dist == "poisson") colVal <- qpois(U, param[1])

  return(colVal)

}

# need to make a little more flexible - to include covariates in mean model, and link

genCopulaCor <- function(n, nvars, params, dist, rho, corstr) {

  mu <- rep(0, nvars)
  dt <- genCorData(n, mu, sigma = 1, rho = rho, corstr = corstr )
  dtM <- melt(dt, id.vars = "id", variable.factor = TRUE, value.name = "Y", variable.name = "seq")
  dtM[, period := as.integer(seq) - 1]
  setkey(dtM, "id")
  dtM[, seqid := .I]
  dtM[, param1 := params[seq], keyby = seqid]
  dtM[, U := pnorm(Y)]
  dtM[, X := genFromU(U, dist, param1), keyby = seqid]

  return(dtM)

}

#### binary outcomes

dist <- "binary"
# corStr <- "cs"
corStr <- "ar1"
results <- data.table()
n = 5

for (i in seq(0.05, 0.95, by = 0.01)) {
  for (p in seq(0.1, 0.5, 0.1)) {

      params <- rep(p, n)

      dtBin <- genCopulaCor(n=1500, nvars = n, params = params, dist = dist, rho = i, corstr = corStr)
      geefit <- geepack::geeglm(X ~ period, data = dtBin, id = id, corstr = "exchangeable")

      results <- rbind(results, data.table(rho = i, p, brho = geefit$geese$alpha))
  }
}

results1 <- copy(results)
results2 <- copy(results)

p2 <- ggplot(data = results2, aes(x = rho, y = brho, group = p)) +
  geom_smooth(se = FALSE, aes(color=p)) +
  geom_abline(slope = 1, intercept = 0, lty = 3) +
  theme_ksg("grey95") +
  scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  xlab(expression(rho[Normal])) +
  ylab(expression(rho[Binary])) +
  theme(axis.title = element_text(size = 14)) +
  ggtitle("Binary (ar-1)")

png("static/img/post-copula/binaryAR1.png", height = 360, width = 400)
 p1
dev.off()

dist <- "poisson"
corStr = "cs"
corStr = "ar1"
results <- data.table()
n = 5

for (i in seq(0.05, 0.95, by = 0.01)) {
  for (p in seq(2, 20, by = 2)) {

    params <- rep(p, n)

    dtPois <- genCopulaCor(n=1500, nvars = n, params = params, dist = dist, rho = i, corstr = corStr)
    geefit <- geepack::geeglm(X ~ period, data = dtPois, id = id, corstr = corStr)

    results <- rbind(results, data.table(rho = i, p, brho = geefit$geese$alpha))
  }
}

results3 <- copy(results)
results4 <- copy(results)

p4 <- ggplot(data = results4[n !=1], aes(x = rho, y = brho, group = p)) +
  geom_smooth(se = FALSE, aes(color=p)) +
  geom_abline(slope = 1, intercept = 0, lty = 3) +
  theme_ksg("grey95") +
  scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  xlab(expression(rho[Normal])) +
  ylab(expression(rho[Poisson])) +
  theme(axis.title = element_text(size = 14)) +
  ggtitle("Poisson (ar-1)")

png("static/img/post-copula/dists.png", height = 360, width = 560)
gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()






