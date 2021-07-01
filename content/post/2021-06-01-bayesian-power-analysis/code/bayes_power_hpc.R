library(simstudy)
library(data.table)
library(ggplot2)
library(bayesplot)
library(posterior)
library(cmdstanr)
library(slurmR)
library(collapse)

s_define <- function() {
    
  defB <- defDataAdd(varname = "Y", formula = "-2 + ..lor * rx", 
                     dist = "binary", link="logit")
  
  return(list(defB = defB)) # list_of_defs is a list of simstudy data definitions
}

s_generate <- function(list_of_defs, argsvec) {
    
  list2env(list_of_defs, envir = environment())
  list2env(as.list(argsvec), envir = environment())
  
  #--- add data generation code ---#
  
  lor <- rnorm(1, mu.lor, sigma.lor)
  
  dT <- genData(nobs)
  dT <- trtAssign(dT, grpName = "rx")
  dT <- addColumns(defB, dT)
  
  return(dT[])
  
}

s_model <- function(generated_data, mod, argsvec) {
    
  list2env(as.list(argsvec), envir = environment())
    
  dt_to_list <- function(dx) {
    
    N <- nrow(dx)                 ## number of observations 
    y <- dx$Y                      ## individual outcome 
    x <- dx$rx                     ## treatment arm for individual 
    s <- t_sigma
    mu <- 0 # can be mu.lor
      
    list(N=N, y=y, x=x, s=s, mu = mu)
  }

  fit <- mod$sample(
    data = dt_to_list(generated_data),
    refresh = 0,
    chains = 4L,
    parallel_chains = 4L,
    iter_warmup = 1000,
    iter_sampling = 4000,
    step_size = 0.1,
    show_messages = FALSE
  )
    
  res <- data.table(fit$summary(variables = "beta"))[, .(median, sd, q95, len = q95-q5)]
    
  draws_array <- as_draws_array(fit$draws())
  betas <- data.table(beta = as.matrix(draws_array[,,"beta"]))
  res$p0 <- mean(betas$beta.V1 < 0)

  return(res) # model_results is a data.table
  }

s_single_rep <- function(list_of_defs, argsvec, mod) {
    
  set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
  
  list_of_defs <- s_define()
  generated_data <- s_generate(list_of_defs, argsvec)
  model_results <- s_model(generated_data, mod, argsvec)
  
  return(model_results)
}
  
s_replicate <- function(argsvec, nsim, mod) {
  
  list_of_defs <- s_define()
  
  model_results <- 
    lapply(
      X = 1 : nsim, 
      FUN = function(x) s_single_rep(list_of_defs, argsvec, mod)
    )
  
  #--- add summary statistics code ---#
  
  model_sums <- unlist2d(lapply(model_results, function(x) x), idcols = "replicate", DT = TRUE)
  summary_stats <- model_sums[ , 
    .(p_95 = mean(p0 >= 0.95), 
      p_len = mean(len <= 2),
      p_sd = mean(sd <= 0.5))
  ]
  
  model_ests <- data.table(t(argsvec), summary_stats)

  return(model_ests)
  
}

###

scenario_list <- function(...) {
  argmat <- expand.grid(...)
  return(asplit(argmat, MARGIN = 1))
}

mu.lor <- c(0, -0.5, -1.0, -1.5)
sigma.lor <- c(0.25)
nobs <- c(100, 150, 200, 250, 300, 350, 400)
t_sigma <- c(1, 5, 10)

scenarios <- scenario_list(mu.lor = mu.lor, sigma.lor = sigma.lor, 
                           nobs = nobs, t_sigma = t_sigma)

set_cmdstan_path(path = ".../cmdstan/2.25.0")
mod <- cmdstan_model("present.stan")

job <- Slurm_lapply(
  X = scenarios, 
  FUN = s_replicate, 
  mod = mod,
  nsim = 1200,
  njobs = min(length(scenarios), 90L), 
  mc.cores = 4L,
  job_name = "i_bp",
  tmp_path = "/gpfs/data/troxellab/ksg/scratch",
  plan = "wait",
  sbatch_opt = list(time = "03:00:00", partition = "cpu_short"),
  export = c("s_single_rep", "s_define", "s_generate", "s_model"),
  overwrite = TRUE
)

summary_stats <- Slurm_collect(job)
final_tab <- rbindlist(summary_stats)

save(final_tab, file = ".../bp.rda")

