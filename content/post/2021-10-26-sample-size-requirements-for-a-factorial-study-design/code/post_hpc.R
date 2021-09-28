library(cmdstanr)
library(simstudy)
library(data.table)
library(posterior)
library(slurmR)
library(glue)

s_define <- function() {
  
  #--- data definition code ---#
  
  f <- "..t_0 + ..t_a*a + ..t_b*b + ..t_c*c + 
      ..t_ab*a*b + ..t_ac*a*c + ..t_bc*b*c + ..t_abc*a*b*c"
  
  defY <- defDataAdd(varname = "y", formula = f, dist = "binary", link="logit")
  
  return(list(defY = defY)) 
  
}

s_generate <- function(list_of_defs, argsvec) {
  
  list2env(list_of_defs, envir = environment())
  list2env(as.list(argsvec), envir = environment())
  
  t_0 <- mu_int
  t_a <- rnorm(1, mu_a, .10)
  t_b <- rnorm(1, mu_b, .10)
  t_c <- rnorm(1, mu_c, .10)
  t_ab <- rnorm(1, mu_ab, .10)
  t_ac <- rnorm(1, mu_ac, .10)
  t_bc <- rnorm(1, mu_bc, .10)
  t_abc <- mu_abc
  
  dd <- genData(8 * n)
  dd <- addMultiFac(dd, nFactors = 3, colNames = c("a", "b", "c"))
  dd <- addColumns(defY, dd)
  
  return(dd)
  
}

s_model <- function(generated_data, mod) {
  
  dt_to_list <- function(dx) {
    
    N <- nrow(dx)                               ## number of observations 
    x_abc <- model.matrix(~a*b*c, data = dx)
    y <- dx[, y]
    
    list(N = N, x_abc = x_abc, y = y)
  }
  
  fit <- mod$sample(
    data = dt_to_list(generated_data),
    refresh = 0,
    chains = 4L,
    parallel_chains = 4L,
    iter_warmup = 500,
    iter_sampling = 2500,
    adapt_delta = 0.98,
    max_treedepth = 20,
    show_messages = FALSE
  )
  
  posterior <- data.frame(as_draws_rvars(fit$draws(variables = "lOR")))
  
  pcts <- c(.025, 0.25, .50, 0.75, .975)
  sumstats <- data.table(t(quantile(posterior$lOR, pcts)))
  setnames(sumstats, glue("p{pcts}"))
  sumstats$sd <- sd(posterior$lOR)
  sumstats$var <- glue("lOR[{1:7}]") 
  
  return(sumstats) # model_results is a data.table
  
}

s_replicate <- function(argsvec, mod) {
  
  set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
  
  list_of_defs <- s_define()
  generated_data <- s_generate(list_of_defs, argsvec)
  model_results <- s_model(generated_data, mod)
  
  #--- summary statistics ---#
  
  summary_stats <- data.table(t(argsvec), model_results)
  
  return(summary_stats) # summary_stats is a data.table
}

#--- Set arguments ---#

scenario_list <- function(...) {
  argmat <- expand.grid(...)
  return(asplit(argmat, MARGIN = 1))
}

n <- c(400, 450, 500, 550, 600, 650)

mu_int <- -1.4
mu_m <- 0.5
mu_x <- -0.3
mu_abc <- 0.3

scenarios <- scenario_list(n = n,
                           mu_int = mu_int, mu_a = mu_m, mu_b = mu_m, mu_c = mu_m, 
                           mu_ab = mu_x, mu_ac = mu_x, mu_bc = mu_x, mu_abc = mu_abc)

scenarios <- rep(scenarios, each = 250)

#--- run locally ---#

# scenarios <- rep(scenarios, each = 2)
# 
# if (file.exists("code/model_ind")) {
#   unlink("code/model_ind")
# }
# 
# mod <- cmdstan_model("code/model_ind.stan")
# summary_stats <- lapply(scenarios, function(a) s_replicate(a, mod))

#--- run on HPC ---#

set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
smodel <- cmdstan_model("/gpfs/data/troxellab/ksg/r/model_ind.stan")

job <- Slurm_lapply(
  X = scenarios, 
  FUN = s_replicate, 
  mod = smodel,
  njobs = min(90L, length(scenarios)), 
  mc.cores = 4L,
  job_name = "i_ss",
  tmp_path = "/gpfs/data/troxellab/ksg/scratch",
  plan = "wait",
  sbatch_opt = list(time = "12:00:00", partition = "cpu_short", `mem-per-cpu` = "4G"),
  export = c("s_define", "s_generate", "s_model"),
  overwrite = TRUE
)

job
res <- Slurm_collect(job)

save(res, file = "/gpfs/data/troxellab/ksg/r/post_ss.rda")

