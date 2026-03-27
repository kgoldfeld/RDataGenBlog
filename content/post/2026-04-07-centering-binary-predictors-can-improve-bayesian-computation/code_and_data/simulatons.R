library(simstudy)
library(data.table)
library(ggplot2)
library(rjags)
library(coda)
library(posterior)
library(parallel)
library(broom)

setwd("~/git R projects/rdatagen/content/post/2026-04-07-centering-binary-predictors-can-improve-bayesian-computation/code_and_data")

theme_set(theme_minimal(base_size = 13))
set.seed(824)

s_gen <- function(n = 2000,
                    alpha = -0.8,
                    beta_a = 0.5,
                    beta_b = 0.9,
                    beta_ab = -0.3) {
  
  def <- 
    defData(varname = "A", formula = 0.5, dist = "binary") |>
    defData(varname = "B", formula = 0.5, dist = "binary") |>
    defData(varname = "AB", formula = "A*B", dist = "nonrandom") |>
    defData(varname = "A_c", formula = "A - 0.5", dist = "nonrandom") |>
    defData(varname = "B_c", formula = "B - 0.5", dist = "nonrandom") |>
    defData(varname = "AB_c", formula = "A_c * B_c", dist = "nonrandom") |>
    defData(
      varname = "Y", 
      formula = "..alpha + ..beta_a * A + ..beta_b * B + ..beta_ab * AB",
      dist = "binary", link = "logit"
    )
    
  genData(n, def)
  
}

dd <- s_gen()

fit_01 <- glm(Y ~ A * B, data = dd, family = binomial)
fit_c  <- glm(Y ~ A_c * B_c, data = dd, family = binomial)

tidy(fit_01)
tidy(fit_c)

# Would like to show that lors are different - and point out that Standard Errors are smaller

model_01 <- "
model {
  for (i in 1:N) {
    Y[i] ~ dbern(p[i])
    logit(p[i]) <- alpha + beta_a * A[i] + beta_b * B[i] + beta_ab * AB[i]
  }
  
  alpha   ~ dnorm(0, 0.25)
  beta_a  ~ dnorm(0, 0.25)
  beta_b  ~ dnorm(0, 0.25)
  beta_ab ~ dnorm(0, 25)
}
"

model_c <- "
model {
  for (i in 1:N) {
    Y[i] ~ dbern(p[i])
    logit(p[i]) <- alpha + gamma_a * A_c[i] + gamma_b * B_c[i] + gamma_ab * AB_c[i]
  }
  
  alpha   ~ dnorm(0, 0.25)
  gamma_a  ~ dnorm(0, 0.25)
  gamma_b  ~ dnorm(0, 0.25)
  gamma_ab ~ dnorm(0, 25)
}
"

fit_jags <- function(dat, model_string, centered = FALSE,
                     n_chains = 3, burn = 2000, n_iter = 5000) {
  
  # jdat <- as.list(dat[, .(Y, A, B, AB, A_c, B_c, AB_c)])
  if (centered) {
    jdat <- as.list(dat[, .(Y, A_c, B_c, AB_c)])
    vars <- c("alpha", "gamma_a", "gamma_b", "gamma_ab")
  } else {
    jdat <- as.list(dat[, .(Y, A, B, AB)])
    vars <- c("alpha", "beta_a", "beta_b", "beta_ab")
  }
  jdat$N <- nrow(dat)
  
  mod <- jags.model(
    textConnection(model_string),
    data = jdat,
    n.chains = n_chains,
    quiet = TRUE
  )
  
  update(mod, burn, progress.bar = "none")
  
  samp <- coda.samples(
    mod,
    variable.names = vars,
    n.iter = n_iter,
    progress.bar = "none"
  )
  
  samp
}

samp_01 <- fit_jags(dd, model_01, centered = FALSE)
samp_c  <- fit_jags(dd, model_c, centered = TRUE)

diag_tbl <- function(samp, model_name) {
  post <- as_draws_df(samp)
  summ <- summarise_draws(post)
  out <- as.data.table(summ)
  out[, model := model_name]
  out[]
}

diag_01 <- diag_tbl(samp_01, "0/1-coded")
diag_c  <- diag_tbl(samp_c, "centered")

ddiag <- rbind(diag_01, diag_c)
num_cols <- names(ddiag)[sapply(ddiag, is.numeric)]
ddiag[, (num_cols) := lapply(.SD, round, 3), .SDcols = num_cols]

get_lor_summary <- function(samp, model_name) {
  dt <- as.data.table(as_draws_df(samp))
  
  if (model_name == "0/1-coded") {
    dt[, lOR_A := beta_a]
    dt[, lOR_B := beta_b]
    dt[, lOR_AB := beta_a + beta_b + beta_ab]
  } else {
    dt[, lOR_A := gamma_a - 0.5 * gamma_ab]
    dt[, lOR_B := gamma_b - 0.5 * gamma_ab]
    dt[, lOR_AB := gamma_a + gamma_b]
  }
  
  dt[, .(
    mean_A = mean(lOR_A),
    mean_B = mean(lOR_B),
    mean_AB = mean(lOR_AB),
    sd_A = sd(lOR_A),
    sd_B = sd(lOR_B),
    sd_AB = sd(lOR_AB)
  )]
}

lor_01 <- get_lor_summary(samp_01, "0/1-coded")
lor_c  <- get_lor_summary(samp_c,  "centered")

bayes_lors <- rbindlist(list(
  cbind(model = "0/1-coded", lor_01),
  cbind(model = "centered",  lor_c)
))

num_cols <- names(bayes_lors)[sapply(bayes_lors, is.numeric)]
bayes_lors[, (num_cols) := lapply(.SD, round, 3), .SDcols = num_cols]


save(ddiag, bayes_lors, file = "ddiag.Rdata")


post_01 <- as_draws_array(samp_01)
post_c  <- as_draws_array(samp_c)

plot_trace_dt <- function(draws, model_name, vars) {
  keep = c(".chain", ".iteration", vars)
  dt <- as.data.table(as_draws_df(draws))
  dt <- melt(dt[, ..keep],
             id.vars = c(".chain", ".iteration"),
             variable.name = "parameter",
             value.name = "value")
  dt[, model := model_name]
  
  dt[, parameter_plot := fcase(
    parameter %in% c("alpha"), "intercept",
    parameter %in% c("beta_a",  "gamma_a"),  "main effect of A",
    parameter %in% c("beta_b",  "gamma_b"),  "main effect of B",
    parameter %in% c("beta_ab", "gamma_ab"), "interaction",
    default = as.character(parameter)
  )]
  
  dt[, parameter_plot := factor(parameter_plot, 
    levels = c("intercept", "main effect of A", "main effect of B", "interaction"))]
  
  dt[]
}

tr_01 <- plot_trace_dt(
  post_01, "0/1-coded", 
  vars = c("alpha", "beta_a", "beta_b", "beta_ab")
)

tr_c  <- plot_trace_dt(
  post_c, "centered",
  vars = c("alpha", "gamma_a", "gamma_b", "gamma_ab")
)

tr_dt <- rbindlist(list(tr_01, tr_c))

trace_plot <- ggplot(tr_dt,
       aes(.iteration, value, group = interaction(.chain, model), color = factor(.chain))) +
  geom_line(alpha = 0.7, linewidth = 0.3) +
  facet_grid(parameter_plot ~ model, scales = "free_y") +
  labs(x = "Iteration", y = NULL, color = "Chain") +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 8),
        panel.spacing = unit(.15, "inches"))

ggsave(file = "trace_plot.png", trace_plot, width = 5, height = 3, scale = 1.9)

###

one_run <- function(
  n = 2000,
  truth = c(alpha = -0.8, beta_a = 0.5, beta_b = 0.9, beta_ab = -0.3),
  burn = 1000,
  n_iter = 3000) {
  
  dd <- s_gen(
    n = n,
    alpha = truth["alpha"],
    beta_a = truth["beta_a"],
    beta_b = truth["beta_b"],
    beta_ab = truth["beta_ab"]
  )
  
  samp_01 <- fit_jags(dd, model_01, centered = FALSE, burn = burn, n_iter = n_iter)
  samp_c  <- fit_jags(dd, model_c, centered = TRUE, burn = burn, n_iter = n_iter)
  
  get_metrics <- function(samp, model_name) {
    post <- as_draws_df(samp)
    summ <- as.data.table(summarise_draws(post))
    summ[, model := model_name]
    summ[]
  }
  
  out <- rbindlist(list(
    get_metrics(samp_01, "0/1-coded"),
    get_metrics(samp_c,  "centered")
  ))
  
  out[]
}

nsim <- 500

sim_res <- rbindlist(mclapply(seq_len(nsim), function(i) {
  out <- one_run()
  out[, sim := i]
  out[]
}, mc.cores = 5))

sim_diag <- sim_res[, .(
  mean_rhat = mean(rhat, na.rm = TRUE),
  median_rhat = median(rhat, na.rm = TRUE),
  mean_ess_bulk = mean(ess_bulk, na.rm = TRUE),
  median_ess_bulk = median(ess_bulk, na.rm = TRUE),
  mean_ess_tail = mean(ess_tail, na.rm = TRUE)
), by = .(model, variable)]

setkey(sim_diag, "variable", "model")
sim_diag[, .(variable, model, mean_rhat, mean_ess_bulk)]

sim_res[, variable := factor(
  variable,
  levels = c("alpha", "beta_a", "beta_b", "beta_ab", "gamma_a", "gamma_b", "gamma_ab")
)]

###


sim_res[, role := fcase(
  variable %in% c("alpha", "alpha_s"), "intercept",
  variable %in% c("beta_a", "gamma_a"), "A",
  variable %in% c("beta_b", "gamma_b"), "B",
  variable %in% c("beta_ab", "gamma_ab"), "AB",
  default = as.character(variable)
)]

sim_res[, role := factor(role, levels = c("intercept", "A", "B", "AB"))]

# save(sim_res, file = "sim.Rdata")
# load("sim.Rdata")

p_rhat <- ggplot(sim_res, aes(model, rhat)) +
  geom_jitter(width = 0.06, alpha = 0.15, size = 0.7, color = "grey50") +
  geom_boxplot(
    outliers = FALSE, color = "grey20", fill = "#E69F00",
    width = .2, coef = 0, median.linewidth = .25, alpha = .6
  ) +
  facet_wrap(
    ~ role,
    labeller = as_labeller(
      c(
        intercept = "Intercept",
        A = "A~effect",
        B = "B~effect",
        AB = "A %*% B~interaction"
      ),
      label_parsed
    )
  ) +
  labs(x = NULL, y = "R-hat") +
  scale_y_continuous(breaks = c(1, 1.0050, 1.0100, 1.0150)) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 8, face = "bold")
  )

p_rhat

###
  
ratio_dt <- dcast(
  sim_res,
  sim + role ~ model,
  value.var = "ess_bulk"
)

# ratio_dt[, variable := factor(
#   variable,
#   levels = c("alpha", "beta_a", "beta_b", "beta_ab")
# )]

ratio_dt[, ratio := centered / `0/1-coded`]

p_ess <- ggplot(ratio_dt, aes(x = role, y = ratio)) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  # geom_jitter(width = 0.08, alpha = 0.25, size = .8, color = "grey50") +
  geom_jitter(width = 0.08, alpha = 0.15, size = 0.7, color = "grey50") +
  geom_boxplot(
    outliers = FALSE, color = "grey20", fill = "#E69F00",
    width = .3, coef = 0, median.linewidth = .25, alpha = .6) + 
  scale_x_discrete(
    labels = c(
      alpha   = expression(alpha),
      beta_a  = expression(beta[a]),
      beta_b  = expression(beta[b]),
      beta_ab = expression(beta[ab])
    )
  ) +
  scale_y_continuous(breaks = c(1,3,5,7, 9, 11)) +
  coord_flip() +
  labs(x = NULL, y = "ESS ratio: centered / 0/1-coded") +
  theme(axis.text.y = element_text(size = 12),
        panel.grid = element_blank())

ggsave(get_versioned_filename(("plot_rhat")), plot = p_rhat, width = 4, height = 3, scale = 1.2)
ggsave(get_versioned_filename(("plot_ess")), plot = p_ess, width = 5, height = 3, scale = 1.2)

# save(fit_01, fit_c, post_01, post_c, tr_01, tr_c, tr_dt, file = "models.Rdata")

###

library(posterior)
library(GGally)

# convert mcmc.list → draws_df
post <- as_draws_df(samp_01)   # samp = your mcmc.list

labs_beta <- c(
  alpha   = "alpha",
  beta_a  = "beta[a]",
  beta_b  = "beta[b]",
  beta_ab = "beta[ab]"
)

p_cor_01 <- GGally::ggpairs(
  post[, c("alpha", "beta_a", "beta_b", "beta_ab")],
  lower = list(continuous = wrap(ggally_points, size = 0.1, color = "#FF9999")),
  diag  = list(continuous = "densityDiag"),
  upper = list(
    continuous = GGally::wrap("cor", size = 4, stars = FALSE)
  ),
  columnLabels = labs_beta,
  labeller = labeller(.cols = label_parsed, .rows = label_parsed)
  ) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 14, face = "bold")
  )

post <- as_draws_df(samp_c)   # samp = your mcmc.list

labs_gamma <- c(
  alpha    = "alpha",
  gamma_a  = "gamma[a]",
  gamma_b  = "gamma[b]",
  gamma_ab = "gamma[ab]"
)

p_cor_c <- GGally::ggpairs(
  post[, c("alpha", "gamma_a", "gamma_b", "gamma_ab")],
  lower = list(continuous = wrap(ggally_points, size = 0.1, color = "#FF9999")),
  diag  = list(continuous = "densityDiag"),
  upper = list(
    continuous = GGally::wrap("cor", size = 4, stars = FALSE)
  ),
  columnLabels = labs_gamma,
  labeller = labeller(.cols = label_parsed, .rows = label_parsed)
  ) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 14, face = "bold")
  )

ggsave(get_versioned_filename(("plot_cor_01")), plot = p_cor_01, width = 6, height = 4, scale = 1.2)
ggsave(get_versioned_filename(("plot_cor_c")), plot = p_cor_c, width = 6, height = 4, scale = 1.2)



