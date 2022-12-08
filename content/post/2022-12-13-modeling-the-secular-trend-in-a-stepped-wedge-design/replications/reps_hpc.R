library(simstudy)
library(ggplot2)
library(data.table)
library(mgcv)
library(lme4)
library(splines)
library(slurmR)

s_define <- function() {
  
  def <- defData(varname = "a", formula = 0, variance = 9)
  def <- defData(def, varname = "mu_b", formula = 0, dist = "nonrandom")
  def <- defData(def, varname = "s2_b", formula = 16, dist = "nonrandom")
  
  defOut <- defDataAdd(varname = "y", formula = "a + b + 5 * A", variance = 40)
  
  return(list(def = def, defOut = defOut)) 
}

s_generate <- function(list_of_defs) {
  
  list2env(list_of_defs, envir = environment())
  
  #--- add data generation code ---#
  
  ds <- genData(24, def, id = "site")
  ds <- addPeriods(ds, 25, "site", perName = "k")
  ds <- addCorGen(dtOld = ds, idvar = "site", 
                  rho = 0.7, corstr = "ar1",
                  dist = "normal", param1 = "mu_b", param2 = "s2_b", cnames = "b")
  ds <- trtStepWedge(ds, "site", nWaves = 24, lenWaves = 1, startPer = 1, 
                     grpName = "A", perName = "k")
  ds$site <- as.factor(ds$site)
  
  dd <- genCluster(ds, "timeID", numIndsVar = 30, level1ID = "id")
  dd <- addColumns(defOut, dd)
  dd[, normk := (k - min(k))/(max(k) - min(k))]
  
  return(dd) #  generated_data is a data.table
}

s_model <- function(dd) {
  
  fitlme_k <- lmer(y ~ A + factor(k) - 1 + (1|site)  , data = dd) # 4
  res_fitlme_k <- summary(fitlme_k)$coefficients["A", c("Estimate", "Std. Error")]
  
  knots <- c(.2, .4, .6, .8)
  
  fitlme_s <- lmer(y ~ A + ( splines::ns(normk, knots = knots) - 1 | site ), data = dd) 
  res_fitlme_s <- summary(fitlme_s)$coefficients["A", c("Estimate", "Std. Error")]
  
  fitgam <- gam(y ~ A + s(k, site, bs = "fs", k = 20), data = dd)  
  res_fitgam <- c(summary(fitgam)$p.coeff["A"], summary(fitgam)$se["A"])

  
  model_results <- data.table(t(res_fitlme_k), t(res_fitlme_s), t(res_fitgam))
  setnames(model_results, c("est.lmek", "se.lmek", "est.lmes", "se.lmes",
                            "est.gam", "se.gam"))

  return(model_results) # model_results is a data.table
}

s_single_rep <- function(args, list_of_defs) {
  
  generated_data <- s_generate(list_of_defs)
  model_results <- s_model(generated_data)
  
  return(model_results)
}

job <- Slurm_lapply(
  X = 1:5400,
  FUN = s_single_rep,
  list_of_defs = s_define(),
  njobs = 90,
  tmp_path = "/gpfs/data/troxellab/ksg/scratch",
  overwrite = TRUE,
  job_name = "i_sw",
  sbatch_opt = list(time = "03:00:00", partition = "cpu_short"),
  export = c("s_generate", "s_model", "s_single_rep"),
  plan = "wait"
)
  
res <- Slurm_collect(job)

save(res, file = "/gpfs/data/troxellab/ksg/r/sw_smooth.rda")

