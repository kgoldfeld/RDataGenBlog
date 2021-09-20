library(simstudy)
library(data.table)
library(cmdstanr)
library(caret)
library(ggpubr)
library(posterior)
library(bayesplot)
library(ggdist)
library(glue)

t_0 <- 0
t_a <- c(-8, -1, 3, 6)
t_b <- c(-3, -1, 0, 4)

x <- c(4, 3) 
nox <- - sum(x) / (16 - length(x))

t_ab <- matrix(c(nox, nox, nox, nox,
                 nox,   4, nox, nox,
                 nox,   3, nox, nox,
                 nox, nox, nox, nox), nrow = 4, byrow = TRUE)

d1 <- defDataAdd(varname = "y", formula = "mu", variance = 16, dist = "normal")

set.seed(110)

dd <- genMultiFac(nFactors = 2, levels = 4, each = 30, colNames = c("a", "b"))
dd[, mu := t_0 + t_a[a] + t_b[b] + t_ab[a, b], keyby = id]
dd <- addColumns(d1, dd)

dsum <- dd[, .(yhat = mean(y)), keyby = .(a, b)]

ggplot(data = dsum, aes(x=b, y = yhat)) +
  geom_vline(aes(xintercept = b), color = "white", size = .25) +
  geom_line(color = "#06295e", size =1.25) +
  facet_grid(.~a, labeller = labeller(a = label_both)) +
  theme(panel.grid = element_blank()) 

dt_to_list <- function(dx) {
  
  dx[, a_f := factor(a)]
  dx[, b_f := factor(b)]
  
  dv <- dummyVars(~ b_f:a_f , data = dx, n = c(4, 4))
  dp <- predict(dv, dx )
  
  N <- nrow(dx)                               ## number of observations 
  I <- 2
  X2 <- 1
  
  main <- as.matrix(dx[,.(a,b)])

  ab <- as.vector(dp %*% c(1:16))  
  x <- as.matrix(ab, nrow = N, ncol = X2)
  
  y <- dx[, y]
  
  list(N=N, I=I, X2=X2, main=main, x=x, y=y)
}


dir <- "/Users/keith/git R projects/rdatagen/content/post/2021-09-28-analyzing-a-factorial-trial-with-a-bayesian-model"
if (getwd() != dir) setwd(dir)

if (file.exists("code/model_4_factors_alt")) unlink("code/model_4_factors_alt")
mod <- cmdstan_model("code/model_4_factors_alt.stan", force_recompile = TRUE)

fit <- mod$sample(
  data = dt_to_list(dd),
  refresh = 500,
  chains = 4L,
  parallel_chains = 4L,
  iter_warmup = 500,
  iter_sampling = 2500,
  adapt_delta = 0.99,
  step_size = .05,
  max_treedepth = 20,
  show_messages = FALSE
)

r <- as_draws_rvars(fit$draws(variables = c("t_0","t", "t_x")))

dnew <- data.frame(genMultiFac(nFactors = 2, levels = 4, each = 1, colNames = c("b", "a")))
dnew$yhat <- with(r, rep(t_0, 16) + rep(t[1, ], each = 4) + rep(t[2, ], times = 4) + t(t_x))

ggplot(data = dnew, aes(x=b, dist = yhat)) +
  geom_vline(aes(xintercept = b), color = "white", size = .25) +
  stat_dist_lineribbon() +
  geom_point(data = dsum, aes(y = yhat), color = "white", size = 2) +
  facet_grid(.~a, labeller = labeller(a = label_both)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())  + 
  scale_fill_brewer()

posterior <- as_draws_array(fit$draws())

mcmc_trace(posterior, pars = glue("t[{1},{1:4}]"))
mcmc_trace(posterior, pars = glue("t_x[1,{1:16}]"))

mcmc_intervals(posterior, pars = c(glue("t[{1},{1:4}]"), glue("t[{2},{1:4}]")))
mcmc_intervals(posterior, pars = glue("t_x[1,{1:16}]"))
