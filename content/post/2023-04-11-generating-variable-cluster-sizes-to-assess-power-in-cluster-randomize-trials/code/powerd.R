library(simstudy)
library(data.table)
library(lmerTest)

s_define <- function() {
  
  d0 <- defData(varname = "a", formula = 0, variance = "..v")
  d0 <- defData(d0, varname = "rx", formula = "1;1", dist = "trtAssign")
  d0 <- defData(d0, varname = "n", formula = "..N", variance = "..d", dist = "clusterSize")
  
  d1 <- defDataAdd(varname = "y", formula = "20 + ..B1 * rx + a", 
          variance = "..var.e", dist = "normal")
  
  return(list(d0 = d0, d1 = d1)) 
}

s_generate <- function(list_of_defs, argsvec) {
  
  list2env(list_of_defs, envir = environment())
  list2env(as.list(argsvec), envir = environment())
  
  v <- iccRE(icc, "normal", varWithin = var.e)

  dc <- genData(nclusters, d0, "site")
  dd <- genCluster(dc, "site", "n", "id")
  dd <- addColumns(d1, dd)

  return(list(dc = dc, dd = dd))
}

s_model <- function(generated_data, argsvec) {
  
  list2env(as.list(argsvec), envir = environment())
  
  fit <- lmer(y ~ rx + (1 | site), data = generated_data$dd)
  
  model_results <- list(
    sim.params = data.table(icc = icc, d = d, N = N),
    size.stats = generated_data$dc[, .(avg = mean(n), sd = sd(n), cv = sd(n)/mean(n))],
    vars = data.table(
      re = VarCorr(fit)$site[1], 
      e = summary(fit)$sigma^2, 
      est.icc = VarCorr(fit)$site[1]/(VarCorr(fit)$site[1] + summary(fit)$sigma^2)),
    ests = data.table(t(coef(summary(fit))["rx",]))
  )
  
  return(model_results)
}

s_single_rep <- function(list_of_defs, argsvec) {
  
  generated_data <- s_generate(list_of_defs, argsvec)
  model_results <- s_model(generated_data, argsvec)
  
  return(model_results)
}

s_replicate <- function(argsvec, nsim) {
  
  list_of_defs <- s_define()
  
  model_results <- parallel::mclapply(
      X = 1 : nsim, 
      FUN = function(x) s_single_rep(list_of_defs, argsvec), 
      mc.cores = 4
  )
  
  
  #--- add summary statistics code ---#
  
  # return(summary_stats) # summary_stats is a data.table
}

#---- specify varying power-related parameters ---#

scenario_list <- function(...) {
  argmat <- expand.grid(...)
  return(asplit(argmat, MARGIN = 1))
}

B1 <- 1.6 # effect size
icc <- c(0.01, 0.02, 0.03, 0.04)
var.e <- 16
d <- seq(0, .5, by = .05)
N <- 500
nclusters <- 20

scenarios <- scenario_list(
  B1 = B1, icc = icc, var.e = var.e, d = d, N = N
)

#--- slurm ---#

job <- Slurm_lapply(
  X = scenarios, 
  FUN = s_replicate, 
  nsim = 10,
  njobs = length(scenarios), 
  mc.cores = 4,
  tmp_path = "/gpfs/data/troxellab/ksg/scratch",
  overwrite = TRUE,
  job_name = "i_d",
  sbatch_opt = list(time = "03:00:00", partition = "cpu_short"),
  export = c("s_define", "s_generate", "s_model", "s_single_rep"),
  plan = "wait")

est.list <- Slurm_collect(job)

getitems <- function(a) {
  
  getdatatable <- function(y) {
    if (is.list(y)) return(data.table(y$sim.params, y$ests))
  }
  
  rbindlist(lapply(a, function(x) getdatatable(x)))
}

reslist <- lapply(est.list, function(a) getitems(a))
res <- rbindlist(reslist)

save(res, file = "/gpfs/data/troxellab/ksg/data/powerbyd.rda")


