

coverage <- function(x, a, b) {
  range <- a + c(-1,1) * 1.96 * b
  (range[1] < x & x < range[2])
}

performance <- function(true, est, se) {
  
  bias = mean(est) - true
  rmse = sqrt(mean((true - est)^2))
  true.se = sd(est)
  est.se = mean(se)
  coverage = mean(mapply("coverage", true, est, se)) * 100

  return(list(bias = bias, rmse = rmse, true.se = true.se, 
              est.se = est.se, coverage = coverage))
}

load("sw_smooth_12.8.rda")
res <- rbindlist(res[sapply(res, function(x) length(x) == 6)])

res <- res[est.gam < 20]

rbind(
res[, performance(5, est.lmek, se.lmek)],
res[, performance(5, est.lmes, se.lmes)],
res[, performance(5, est.gam, se.gam)]
)



