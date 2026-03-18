library(simstudy)
library(data.table)
library(ggplot2)

setwd("~/git R projects/rdatagen/content/post/2026-03-16-getting-to-the-bottom-of-tmle-targeting-in-action/code_and_output")

gen_dgp <- function(n, tau = 5) {
  
  def <-
    defData(varname = "x1", formula = 0.5, dist = "binary") |>
    defData(varname = "x2", formula = 0, variance = 1, dist = "normal") |>
    defData(
      varname = "a",
      formula = "-0.2 + 0.8 * x1 + 0.6 * x2",
      dist = "binary",
      link = "logit"
    ) |>
    defData(
      varname = "y",
      formula = "..tau * a + 1.0 * x1 + 1.0 * x2 + 1.5 * x1 * x2",
      variance = 1,
      dist = "normal"
    )
  
  genData(n, def)[]
}

fit_nuisance <- function(dt, scenario) {
  
  if (scenario %in% c("both_correct", "g_wrong")) {
    Q_fit <- lm(y ~ a + x1 + x2 + x1:x2, data = dt)
  } else {
    Q_fit <- lm(y ~ a + x1, data = dt)
  }
  
  if (scenario %in% c("both_correct", "Q_wrong")) {
    g_fit <- glm(a ~ x1 + x2, family = binomial(), data = dt)
  } else {
    g_fit <- glm(a ~ x1, family = binomial(), data = dt)
  }
  
  list(Q_fit = Q_fit, g_fit = g_fit)
}

predict_Q <- function(Q_fit, dt, a_val) {
  nd <- copy(dt)
  nd[, a := a_val]
  as.numeric(predict(Q_fit, newdata = nd))
}

predict_g <- function(g_fit, dt) {
  p <- as.numeric(predict(g_fit, newdata = dt, type = "response"))
  pmin(pmax(p, 0.01), 0.99)
}

phi_ate <- function(dt, Q1, Q0, g, psi) {
  A <- dt$a
  Y <- dt$y
  (Q1 - Q0 - psi) + A / g * (Y - Q1) - (1 - A) / (1 - g) * (Y - Q0)
}

tmle_update_gaussian <- function(dx, Q_fit, g_fit) {
  
  dt <- copy(dx)
  
  dt[, QAW := predict_Q(Q_fit, dt, dt$a)]
  dt[, Q1  := predict_Q(Q_fit, dt, 1)]
  dt[, Q0  := predict_Q(Q_fit, dt, 0)]
  dt[, g   := predict_g(g_fit, dt)]
  dt[, H   := a / g - (1 - a) / (1 - g)]
  
  # Estimate fluctuation parameter on this dataset unless supplied

  fluc_fit <- lm(y ~ -1 + offset(QAW) + H, data = dt)
  eps <- coef(fluc_fit)[["H"]]
  
  # Targeted updates
    
  dt[, QAW_star := QAW + eps * H]
  dt[, Q1_star  := Q1 + eps / g]
  dt[, Q0_star  := Q0 - eps / (1 - g)]
  
  list(
    eps = eps,
    Q1 = dt$Q1,
    Q0 = dt$Q0,
    g = dt$g,
    Q1_star = dt$Q1_star,
    Q0_star = dt$Q0_star
  )
}

s_estimate <- function(dx, scenario, tau, n, nfolds = 2) {
  
  dt <- copy(dx)
  dt[, fold := sample(rep(1:nfolds, length.out = .N))]

  # Storage for fold-specific population diagnostics
  
  p0_phi_plugin_vec <- numeric(nfolds)
  p0_phi_tmle_vec <- numeric(nfolds)
  eps_vec <- numeric(nfolds)
  
  for (k in 1:nfolds) {
    
    train <- dt[fold != k]
    test  <- copy(dt[fold == k])
    
    fits <- fit_nuisance(train, scenario)
    
    # Sample: fit on train, target/evaluate on held-out test fold
    
    tmle_k <- tmle_update_gaussian(test, fits$Q_fit, fits$g_fit)
    eps_vec[k] <- tmle_k$eps
    
    # Write held-out predictions back into main sample object
    
    dt[fold == k, `:=`(
      Q1 = tmle_k$Q1,
      Q0 = tmle_k$Q0,
      g = tmle_k$g,
      Q1_star = tmle_k$Q1_star,
      Q0_star = tmle_k$Q0_star
    )]
  }
  
  # Cross-fitted sample estimators
  
  psi_plugin <- mean(dt$Q1 - dt$Q0)
  psi_tmle   <- mean(dt$Q1_star - dt$Q0_star)
  
  pn_phi_plugin <- mean(
    phi_ate(
      dt,
      Q1 = dt$Q1,
      Q0 = dt$Q0,
      g  = dt$g,
      psi = psi_plugin
    )
  )
  
  pn_phi_tmle <- mean(
    phi_ate(
      dt,
      Q1 = dt$Q1_star,
      Q0 = dt$Q0_star,
      g  = dt$g,
      psi = psi_tmle
    )
  )
  
  data.table(
    tau = tau,
    n = n,
    scenario = scenario,
    psi_true = tau,
    psi_plugin = psi_plugin,
    psi_tmle = psi_tmle,
    abs_pn_phi_plugin = abs(pn_phi_plugin),
    abs_pn_phi_tmle = abs(pn_phi_tmle),
    abs_err_plugin = abs(psi_plugin - tau),
    abs_err_tmle = abs(psi_tmle - tau),
    eps = mean(eps_vec)
  )
}

s_simulate <- function(n, tau, scenarios) {
  
  dd <- gen_dgp(n, tau)
  rbindlist(lapply(scenarios, function(s) s_estimate(dd, s, tau, n)))
  
}

set.seed(2026)

tau <- 5
ns = rep(c(100, 250, 750, 1000, 2000), each = 1000)
scenarios <- c("both_correct", "Q_wrong", "g_wrong", "both_wrong")

res <- rbindlist(
  lapply(ns, function(x) s_simulate(x, tau, scenarios))
)

###

plot_dt <- melt(
  res,
  id.vars = c("n", "scenario"),
  measure.vars = c("abs_pn_phi_plugin", "abs_pn_phi_tmle"),
  variable.name = "stage",
  value.name = "abs_pn_phi"
)

plot_dt[, stage := factor(
  stage,
  levels = c("abs_pn_phi_plugin", "abs_pn_phi_tmle"),
  labels = c("Initial fit", "After targeting")
)]

plot_dt[, scenario := factor(
  scenario,
  levels = scenarios,
  labels = c("Both correct", "Q wrong", "g wrong", "Both wrong")
)]

EIF_plot <- ggplot(plot_dt, aes(x = factor(n), y = abs_pn_phi)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_jitter(width = .1, height = 0, size = .2, color = "grey68") +
  geom_boxplot(outliers = FALSE, color = "grey20", fill = "#E69F00",
               width = .4, coef = 0, median.linewidth = .5) +
  facet_grid(stage ~ scenario) +
  labs(
    x = "Sample size",
    y = expression("|" * P[n] * phi * "|"),
    title = "Magnitude of the empirical mean of the estimated EIF"
  ) +
  ylim(-0.5, 4) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold"))

###

err_dt <- melt(
  res,
  id.vars = c("n", "scenario"),
  measure.vars = c("abs_err_plugin", "abs_err_tmle"),
  variable.name = "estimator",
  value.name = "abs_error"
)

err_dt[, estimator := factor(
  estimator,
  levels = c("abs_err_plugin", "abs_err_tmle"),
  labels = c("Plug-in", "TMLE")
)]

err_dt[, scenario := factor(
  scenario,
  levels = scenarios,
  labels = c("Both correct", "Q wrong", "g wrong", "Both wrong")
)]

err_plot <- ggplot(err_dt, aes(x = factor(n), y = abs_error)) +
  geom_jitter(width = .1, size = .2, color = "grey68") +
  geom_boxplot(outliers = FALSE, color = "grey20", fill = "#E69F00",
               width = .4, coef = 0, median.linewidth = .5) +
  facet_grid(estimator ~ scenario) +
  labs(
    x = "Sample size",
    y = "Absolute error",
    title = "Absolute estimation error by estimator and nuisance scenario"
  ) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold"))

summ <- res[, .(
  mean_abs_pn_phi_plugin = mean(abs_pn_phi_plugin),
  mean_abs_pn_phi_tmle = mean(abs_pn_phi_tmle),
  bias_plugin = mean(psi_plugin - psi_true),
  bias_tmle   = mean(psi_tmle - psi_true),
  rmse_plugin = sqrt(mean((psi_plugin - psi_true)^2)),
  rmse_tmle   = sqrt(mean((psi_tmle - psi_true)^2))
), by = .(scenario, n)]

ggsave(filename = "EIF_plot.png", plot = EIF_plot, ,scale = 1.3, width = 6, height = 3.5)
ggsave(filename = "err_plot.png", plot = err_plot, ,scale = 1.3, width = 6, height = 3.5)

setwd("/Users/KSG/git R projects/rdatagen")
