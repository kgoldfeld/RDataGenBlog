library(rstan)
library(slurmR)
library(simstudy)
library(data.table)
library(lme4)

#--- Stan model

print("compiling stan")

#--- Set up directories and files

lab_dir <- "/gpfs/data/troxellab/ksg"
stan_file <- file.path(lab_dir, "r", "binary.stan")

time_stamp <- gsub("-", "", Sys.Date())

temp_dir <- file.path(lab_dir, "scratch", time_stamp)
data_dir <- file.path(lab_dir, "data", time_stamp)

dir.create(file.path(temp_dir), showWarnings = FALSE)
dir.create(file.path(data_dir), showWarnings = FALSE)

# rt <- stanc("content/post/DataBayesMetaInterim/binary.stan")

rt <- stanc(file = stan_file)
sm <- stan_model(stanc_ret = rt, verbose=FALSE)



### simulations

print("running simulations")

iter <- function(iter, maxn, sm) {
  
  #--- Definitions
  
  def_s <- defDataAdd(varname = "b0", formula = 0, variance = 0.025)
  def_s <- defDataAdd(
    def_s, varname = "delta_k", 
  # formula = "(c_type==1) * 0.4 + (c_type==2) * 0.5 + (c_type==3) * 0.6", 
    formula = "(c_type==1) * 0.0 + (c_type==2) * 0.0 + (c_type==3) * 0.0", 
    dist = "nonrandom"
  )
  
  def_i <- defDataAdd(
    varname = "y", formula = "-1 + b0 + rx * delta_k", 
    dist = "binary", link = "logit")
  
  genDT <- function(n_i, def_s, def_i, n_study, n_clust) {
    
    dc <- genData(n_study, id = "c_type")
    
    ds <- genCluster(dc, "c_type", numIndsVar = n_clust, level1ID = "site")
    ds <- addColumns(def_s, ds)
    ds[, order := .SD[, .(order = .I)]$order, keyby = c_type] # for plotting
    
    di <- genCluster(ds, "site", n_i, "id")
    di <- trtAssign(di, 2, strata = "site", grp = "rx")
    di <- addColumns(def_i, di)
    
    di[]
    
  }
  
  prepData <- function(dx) {
    
    N <- nrow(dx) ;                       
    C <- dx[, length(unique(c_type))]     
    K <- dx[, length(unique(site))]       
    y <- as.numeric(dx$y)                 
    kk <- dx$site                         
    ctrl <- dx$rx                         
    cc <- dx[, .N, keyby = .(site, c_type)]$c_type  
    
    list(N=N, C=C, K=K, y=y, kk=kk, ctrl=ctrl, cc=cc, Delta_prior_sd = 2.5)
    
  }
  
  doMCMC <- function(iter, n, dxfull) {
    
    dx <- dxfull[, .SD[1:n], keyby = .(site)]
    
    sampdat <- prepData(dx)
    
    fit <-  sampling(sm, data=sampdat, 
                     iter = 3000, warmup = 500, show_messages = FALSE, 
                     chains = 4, cores = 4, refresh = 0)
    
    p.eff <- mean(extract(fit, pars = "OR")[[1]] > 1)
    p.clinic <- mean(extract(fit, pars = "OR")[[1]] > 1.2)
    
    data.table(iter, n, p.eff, p.clinic)
    
  }
  
  doGlmer <- function(iter, n, dxfull) {

      dx <- dxfull[, .SD[1:n], keyby = .(site)]
      glmerfit <- glmer(y ~ rx + (1 | site), data = dx, family=binomial)
      pval <- coef(summary(glmerfit))["rx", "Pr(>|z|)"]
      
      data.table(iter, n, pval)
    
  }
  
  cat("\f", iter, format(Sys.time(), "%X"))
  
  dxfull <- genDT(maxn, def_s, def_i, 3, 7)
  
  for (n in seq(40, maxn, 5)) {
    probs <- doMCMC(iter, n, dxfull)
    # success <- (probs$p.eff > .95 & probs$p.clinic > 0.50)
    success <- (probs$p.eff > .90)
    if ( success | n == maxn ) {
      probs[, success := as.integer(success)]
      break
    }
  }
  
  seqn <- seq(40, maxn, 20)
  bounds <- c(0.00305, 0.01625, 0.03070) # from ldbounds
  
  for (n in seq_along(seqn)) {
    pvals <- doGlmer(iter, seqn[n], dxfull)
    if ( pvals$pval < bounds[n] | seqn[n] == maxn ) {
      pvals$success <- as.integer(pvals$pval < bounds[n])
      break
    }
  }
  
  list(probs, pvals)
  
}

# lapply(1:4, function(x) iter(x, 30, sm))

job <- Slurm_lapply(
  X = 1:1050, 
  FUN = iter, 
  maxn = 80, 
  sm = sm,
  njobs = 75,  # 50, 100
  mc.cores = 4,  # 10, 5
  tmp_path = temp_dir,
  job_name = "iter",
  sbatch_opt = list(
    time = "03:00:00", 
    partition="cpu_short"),
  plan = "wait"
)

job

### Save data

print("saving data")

res <- Slurm_collect(job)

resBad <-  Filter((function(l) length(l) == 1), res)
if (length(resBad) > 0) print(resBad[[1]])

resOK <- Filter((function(l) length(l) > 1), res)
print(paste("Number of records with valid estimation:", length(resOK)))

res_stan <- rbindlist(lapply(resOK, function(x) x[[1]]))
res_glmer <- rbindlist(lapply(resOK, function(x) x[[2]]))

save(res_stan, res_glmer, file = file.path(data_dir, "interim.rda"))

