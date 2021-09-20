y <- rnorm(10)
x <- rep(1, 3)

mod <- cmdstan_model("code/test.stan", force_recompile = TRUE)

fit <- mod$sample(
  data = list(nrow = 0 , ncol = 0, y=y, x=x),
  refresh = 500,
  chains = 1L,
  parallel_chains = 1L,
  iter_warmup = 500,
  iter_sampling = 2500,
  adapt_delta = 0.99,
  step_size = .05,
  max_treedepth = 20,
  show_messages = FALSE
)


y <- rnorm(10)
x <- matrix(1, nrow = 3, ncol = 1)

mod <- cmdstan_model("code/test_matrix.stan", force_recompile = TRUE)

fit <- mod$sample(
  data = list(nrow = nrow(x) , ncol = ncol(x), y=y, x=x),
  refresh = 500,
  chains = 1L,
  parallel_chains = 1L,
  iter_warmup = 500,
  iter_sampling = 2500,
  adapt_delta = 0.99,
  step_size = .05,
  max_treedepth = 20,
  show_messages = FALSE
)
