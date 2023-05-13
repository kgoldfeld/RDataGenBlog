library(simstudy)
library(data.table)
library(survival)
library(coxme)
library(lmerTest)
library(parallel)
library(slurmR)

extract_coxme_table <- function (mod){
  beta <- mod$coefficients 
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z <- round(beta/se, 2)
  p <- signif(1 - pchisq((beta/se)^2, 1), 2)
  table=data.table(beta = beta, se = se, z = z, p = p)
  return(table)
}

s_def <- function() {
  
  defc <- defData(varname = "b", formula = 0, variance = "..s2")
  defc <- defData(defc, varname = "rx", formula = "1;1", dist = "trtAssign")

  defa <- defDataAdd(varname = "start_day", formula = "1;182", dist = "uniformInt")
  defa <- defDataAdd(defa, varname = "censor", 
    formula = "365 - start_day ", dist = "nonrandom")

  defs <- defSurv(varname = "ttc", formula = "-4.815 + 0.4 * rx + b", shape = 1.326)

  defa2 <- defDataAdd(varname = "event6", 
    formula = "1*(ttc <= 182)", dist = "nonrandom")
  
  return(list(defc = defc, defa = defa, defs = defs, defa2 = defa2))

}

s_generate <- function(argsvec, list_of_defs) {
  
  list2env(list_of_defs, envir = environment())
  list2env(as.list(argsvec), envir = environment())
  
  dc <- genData(nsites, defc, id = "site")
  dd <- genCluster(dc, "site", ninds, "id")
  dd <- addColumns(defa, dd)
  dd <- genSurv(dd, defs, digits = 0)
  dx <- addCompRisk(dd, events = c("ttc", "censor"), 
    timeName = "time", censorName = "censor", keepEvents = TRUE)
  dx <- addColumns(defa2, dx)
  
  dx[]
  
}

s_replicate <- function(argsvec, list_of_defs) {
  
  dx <- s_generate(argsvec, list_of_defs)

  coxfitm <-coxme(Surv(time, event) ~ rx + (1 | site), data = dx)
  lmefit <- glmer(event6 ~ rx + (1 | site), data = dx, family = binomial)
  m.model <- coxph(Surv(time, event) ~ rx, data = dx, cluster=site)
  
  list2env(as.list(argsvec), envir = environment())
  
  return(data.table(
    nsites = nsites,
    ninds = ninds,
    s2 = s2,
    est_s = fixef(coxfitm), 
    re.var_s = VarCorr(coxfitm)$site,
    p_s = extract_coxme_table(coxfitm)$p,
    p_l = coef(summary(lmefit))["rx", "Pr(>|z|)"],
    p_m = summary(m.model)$robscore["pvalue"]
  ))
  
}

s_scenarios <- function(argsvec, nreps) {
  
  list_of_defs <- s_def()
  
  rbindlist(
    parallel::mclapply(
      X = 1 : nreps, 
      FUN = function(x) s_replicate(argsvec, list_of_defs), 
      mc.cores = 4)
  )
  
}

scenario_list <- function(...) {
  argmat <- expand.grid(...)
  return(asplit(argmat, MARGIN = 1))
}

nsites <- c(18, 24, 30, 36)
ninds <- c(15, 18, 21)
s2 <- c(0.01, 0.02, 0.03, 0.04, 0.05)

scenarios <- scenario_list(nsites = nsites, ninds = ninds, s2 = s2)

# lapply(scenarios, function(a) s_scenarios(a, nrep = 5))

job <- Slurm_lapply(
  X = scenarios, 
  FUN = s_scenarios, 
  nreps = 5000,
  njobs = min(90, length(scenarios)),
  mc.cores = 4,
  tmp_path = "/gpfs/data/troxellab/ksg/scratch",
  overwrite = TRUE,
  job_name = "i_ss",
  sbatch_opt = list(time = "03:00:00", partition = "cpu_short"),
  export = c("extract_coxme_table", "s_def", "s_generate", "s_replicate"),
  plan = "wait")

job
res <- Slurm_collect(job)

res <- rbindlist(res)

resum <- res[, .(p_s = mean(p_s < 0.05), 
                 p_l = mean(p_l < 0.05), 
                 p_m = mean(p_m < 0.05),
                 est_s = mean(est_s),
                 est_s2 = mean(re.var_s)), 
             keyby=.(nsites, ninds, s2)]

save(resum, file = "/gpfs/data/troxellab/ksg/data/ss.rda")
