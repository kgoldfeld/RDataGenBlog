library(simstudy)
library(data.table)
library(parallel)
library(cmdstanr)
library(posterior)
library(rslurm)

freq_fit <- function(dx) {
  
  lmfit <- lm(y ~ rx, data = dx)
  coef(summary(lmfit))["rx", "Pr(>|t|)"]
  
}

bayes_fit <- function(dx, p_sigma, m_effect, decision_rule, x) {
  
  print(x)
  
  data_list <- list(N=nrow(dx), y=dx$y, rx=dx$rx, p_mu=0, p_sigma=p_sigma)

  fit <- mod$sample(
    data = data_list,
    refresh = 0,
    chains = 4L,
    parallel_chains = 4L,
    iter_warmup = 500,
    iter_sampling = 2500,
    step_size = 0.1,
    show_messages = FALSE
  )
  
  df <- data.frame(as_draws_rvars(fit$draws(variables = "beta")))
  
  if (decision_rule == 1) {
    return((mean(df$beta > 0) > 0.95))
  } else { # decision_rule == 2
    return( ((mean(df$beta > 0) > 0.95) & (mean(df$beta > m_effect ) > 0.5)) )  
  }
}

seq_mods <- function(dx, seq, iter, p_sigma, decision_rule, m_effect) {
  
  freq_ps <- sapply(seq(start, end, by = by), function(x) freq_fit(dx[1:x]))
  freq_effect <- any(freq_ps < 0.05)
  
  bayes_ci <- sapply(seq(start, end, by = by), 
    function(x) bayes_fit(dx[1:x], p_sigma, m_effect, decision_rule, x))
  bayes_effect <- any(bayes_ci)
  
  return(data.table(seq, iter, p_sigma, m_effect, decision_rule, freq_effect, bayes_effect))
  
}

s_replicate <- function(iter, p_sigma, decision_rule, m_effect, seq) {
  
  set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
  
  print(c(seq, iter, p_sigma, decision_rule, m_effect))
  
  def <- defData(varname = "rx", formula = "1;1", dist = "trtAssign")
  def <- defData(def, varname = "y", formula = 0, variance = 1, dist = "normal")
  
  dd <- genData(end, def)
  
  seq_mods(dd, seq, iter, p_sigma, decision_rule, m_effect)
  
}

set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.25.0")
mod <- cmdstan_model("multiple.stan")

###

scenario_dt <- function(...) {
  argdt <- data.table(expand.grid(...))
  argdt[, seq := .I]
  argdt[]
}

# start <- 80L
# end <- 160L
# by <- 20L

start <- 100L
end <- 1000L
by <- 100L

# decision_rule = 1
# iter <- c(1:3000)
# p_sigma <- c(1, 5, 10)
# m_effect <- c(0)

decision_rule = 2
iter <- c(1:1500)
p_sigma <- c(1, 5, 10)
m_effect <- c(0.2, 0.3, 0.4)

scenarios <- scenario_dt(
  iter = iter, 
  p_sigma = p_sigma, 
  decision_rule = decision_rule,
  m_effect = m_effect
)

###

sopts <- list(time = '12:00:00', partition = "cpu_short", `mem-per-cpu` = "5G")
sobjs <- c("seq_mods", "freq_fit", "bayes_fit", "mod", "start", "end", "by")

sjob <- slurm_apply(
  f = s_replicate, # the function
  params = scenarios, # a data frame
  jobname = 'mult_i',
  nodes = 50, 
  slurm_options = sopts,
  global_objects = sobjs,
  submit = TRUE
)

res <- data.table(get_slurm_out(sjob, outtype = 'table', wait = TRUE))

path <- "/gpfs/data/troxellab/ksg/data/"
fn <- paste0("mult_", end, "_d", decision_rule, ".rda")

save(res, file = paste0(path, fn))

cleanup_files(sjob, wait = TRUE)



