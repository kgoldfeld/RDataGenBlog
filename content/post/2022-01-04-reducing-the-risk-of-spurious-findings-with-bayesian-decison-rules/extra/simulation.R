library(simstudy)
library(data.table)
library(rslurm)
library(posterior)
library(cmdstanr)

listdat <- function(dx, grpvar) {
  
  dx[, grp := factor(get(grpvar))]
  
  N <- dx[, .N]
  L <- dx[, nlevels(grp)]
  y <- dx[, y]
  rx <- dx[, rx]
  grp <-dx[, as.numeric(grp)]
  
  list(N = N, L = L, y = y, rx = rx, grp = grp)
}

fitmods <-function(dx, grpvar) {
  
  dat <- listdat(dx, grpvar)
  
  fitbayes <- mod$sample(
    data = dat,
    refresh = 0,
    chains = 4L,
    parallel_chains = 4L,
    iter_warmup = 1500,
    iter_sampling = 3000
  )
  
  fitglm <- glm(y ~ grp + rx:grp - 1, data = dx, family = "binomial")
  
  OR <- as_draws_rvars(fitbayes$draws(variables = "OR"))$OR
  bayes_res_2 <- any((mean(OR < 1) > 0.95) & (mean(OR < 0.80 ) > 0.5)) 
  bayes_res_ci <- any((quantile(OR, .025) > 1) | (quantile(OR, .975) < 1))
  
  
  pvals <- coef(summary(fitglm))[, "Pr(>|z|)"]
  lpval <- length(pvals)
  freq_res <- any(pvals[(lpval/2 + 1) : lpval] < 0.05)
  
  # data.table(var = grpvar, bayes_res, freq_res)
  
  return(data.table(var = grpvar, bayes_res_2, bayes_res_ci, freq_res))
}

s_replicate <- function(iter) {
  
  set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
  
  dd <- genData(1000, def)
  res <- rbindlist(lapply(paste0("g", 1:20), function(a) fitmods(dd, a)))
  # res <- rbindlist(lapply(paste0("g", 5), function(a) fitmods(dd, a)))
  
  first_true <- sapply(res[, c(2:4)], function(x) match(TRUE, x))
  
  data.table(iter, res[, .(bayes_2 = first_true[1], 
                           bayes_ci = first_true[2],
                           freq = first_true[3])])
  
}

### Data definitions

genRepeatDef <- function(nvars, prefix, formula, variance, dist, link = "identity") {
  varnames <- paste0(prefix, 1:nvars)
  data.table(varname = varnames, formula=formula, variance=variance, dist=dist, link=link)
}

def <- genRepeatDef(20, "g", "1/3;1/3;1/3", 0, "categorical")
def <- defData(def, "rx", formula = "1;1", dist = "trtAssign")
def <- defData(def, "y", formula = "-2", dist = "binary", link = "logit")
# def <- defData(def, "y", formula = "-2 - 0.9* rx * I(g5 %in% c(1,2))", dist = "binary", link = "logit")

### Compile stan code

set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
mod <- cmdstan_model("simulation.stan")

### Set rslurm arguments

sopts <- list(time = '12:00:00', partition = "cpu_short")
sobjs <- c("listdat", "fitmods", "mod", "def")

### Replicate over iterations

sjob <- slurm_apply(
  f = s_replicate, # the function
  params = data.frame(iter = 1:1050),
  jobname = 'sim_i',
  nodes = 75, 
  cpus_per_node = 4,
  processes_per_node = 4,
  slurm_options = sopts,
  global_objects = sobjs,
  submit = TRUE
)

### Collect the results and save them

res <- get_slurm_out(sjob, outtype = 'table', wait = TRUE)
save(res, file = "/gpfs/data/troxellab/ksg/data/simulation.rda")


