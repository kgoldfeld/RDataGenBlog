library(simstudy)
library(data.table)
library(ggplot2)

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

tmle_update_gaussian <- function(dx, Q_fit, g_fit, eps = NULL) {
  
  dt <- copy(dx)
  
  dt[, QAW := predict_Q(Q_fit, dt, dt$a)]
  dt[, Q1  := predict_Q(Q_fit, dt, 1)]
  dt[, Q0  := predict_Q(Q_fit, dt, 0)]
  dt[, g   := predict_g(g_fit, dt)]
  dt[, H   := a / g - (1 - a) / (1 - g)]
  
  # Estimate fluctuation parameter on this dataset unless supplied
  if (is.null(eps)) {
    fluc_fit <- lm(y ~ -1 + offset(QAW) + H, data = dt)
    eps <- coef(fluc_fit)[["H"]]
  }
  
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

analyze_once <- function(dx, scenario, tau, n) {
  
  dt <- copy(dx)
  dt[, fold := sample(rep(1:2, length.out = .N))]
  
  # Storage for cross-fitted sample predictions
  dt[, `:=`(
    Q1 = NA_real_,
    Q0 = NA_real_,
    g = NA_real_,
    Q1_star = NA_real_,
    Q0_star = NA_real_
  )]
  
  # Storage for fold-specific population diagnostics
  p0_phi_init_vec <- numeric(2)
  p0_phi_tmle_vec <- numeric(2)
  eps_vec <- numeric(2)
  
  for (k in 1:2) {
    
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
    
    # Population diagnostics using the same fold-specific nuisance fit
    # and the same sample-estimated epsilon
    tmle_pop_k <- tmle_update_gaussian(
      pop_dt,
      fits$Q_fit,
      fits$g_fit,
      eps = tmle_k$eps
    )
    
    psi0_pop_k <- mean(tmle_pop_k$Q1 - tmle_pop_k$Q0)
    psitmle_pop_k <- mean(tmle_pop_k$Q1_star - tmle_pop_k$Q0_star)
    
    p0_phi_init_vec[k] <- mean(
      phi_ate(
        pop_dt,
        Q1 = tmle_pop_k$Q1,
        Q0 = tmle_pop_k$Q0,
        g  = tmle_pop_k$g,
        psi = psi0_pop_k
      )
    )
    
    p0_phi_tmle_vec[k] <- mean(
      phi_ate(
        pop_dt,
        Q1 = tmle_pop_k$Q1_star,
        Q0 = tmle_pop_k$Q0_star,
        g  = tmle_pop_k$g,
        psi = psitmle_pop_k
      )
    )
  }
  
  # Cross-fitted sample estimators
  psi_plugin <- mean(dt$Q1 - dt$Q0)
  psi_tmle   <- mean(dt$Q1_star - dt$Q0_star)
  
  pn_phi_init <- mean(
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
    n = n,
    scenario = scenario,
    psi_true = tau,
    psi_plugin = psi_plugin,
    psi_tmle = psi_tmle,
    pn_phi_init = pn_phi_init,
    pn_phi_tmle = pn_phi_tmle,
    p0_phi_init = mean(p0_phi_init_vec),
    p0_phi_tmle = mean(p0_phi_tmle_vec),
    abs_err_plugin = abs(psi_plugin - tau),
    abs_err_tmle = abs(psi_tmle - tau),
    eps = mean(eps_vec)
  )
}

run_sim <- function(n, tau, scenarios) {
  
  dd <- gen_dgp(n, tau)
  rbindlist(lapply(
    scenarios, 
    function(s) analyze_once(dd, s, tau, n)
  ))
  
}

set.seed(2026)

pop_dt <- gen_dgp(5e5)
tau <- 5
ns = rep(c(100, 250, 750, 1000, 2000), each = 500)
scenarios <- c("both_correct", "Q_wrong", "g_wrong", "both_wrong")

res <- rbindlist(
  lapply(ns, function(x) run_sim(x, tau, scenarios))
)

###

plot_dt <- melt(
  res,
  id.vars = c("n", "scenario"),
  measure.vars = c("pn_phi_init", "pn_phi_tmle"),
  variable.name = "stage",
  value.name = "pn_phi"
)

plot_dt[, stage := factor(
  stage,
  levels = c("pn_phi_init", "pn_phi_tmle"),
  labels = c("Initial fit", "After targeting")
)]

plot_dt[, scenario := factor(
  scenario,
  levels = scenarios,
  labels = c("Both correct", "Q wrong", "g wrong", "Both wrong")
)]

ggplot(plot_dt, aes(x = factor(n), y = pn_phi)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_boxplot(outlier.alpha = 0.15, outliers = FALSE) +
  facet_grid(stage ~ scenario) +
  labs(
    x = "Sample size",
    y = expression(P[n] * phi[hat(P)]),
    title = "Empirical mean of the estimated EIF before and after targeting"
  ) +
  # ylim(-0.5, 4) +
  theme(panel.grid.minor = element_blank())

###

plot_0_dt <- melt(
  res,
  id.vars = c("n", "scenario"),
  measure.vars = c("pn_phi_init", "p0_phi_init"),
  variable.name = "stage",
  value.name = "pn_phi"
)

plot_0_dt[, stage := factor(
  stage,
  levels = c("pn_phi_init", "p0_phi_init"),
  labels = c("P_n", "P_0")
)]

plot_0_dt[, scenario := factor(
  scenario,
  levels = scenarios,
  labels = c("Both correct", "Q wrong", "g wrong", "Both wrong")
)]

ggplot(plot_0_dt, aes(x = factor(n), y = pn_phi)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_boxplot(outlier.alpha = 0.15, outliers = FALSE) +
  facet_grid(stage ~ scenario) +
  labs(
    x = "Sample size",
    y = expression(P[n] * phi[hat(P)]),
    title = "Empirical mean of the estimated EIF before and after targeting"
  ) +
  # ylim(-0.5, 4) +
  theme(panel.grid.minor = element_blank())

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

ggplot(err_dt, aes(x = factor(n), y = abs_error)) +
  geom_boxplot(outlier.alpha = 0.15) +
  facet_grid(estimator ~ scenario) +
  labs(
    x = "Sample size",
    y = "Absolute error",
    title = "Absolute estimation error by estimator and nuisance scenario"
  ) +
  theme(panel.grid.minor = element_blank())


summ <- res[, .(
  bias_plugin = mean(psi_plugin - psi_true),
  bias_tmle   = mean(psi_tmle - psi_true),
  rmse_plugin = sqrt(mean((psi_plugin - psi_true)^2)),
  rmse_tmle   = sqrt(mean((psi_tmle - psi_true)^2)),
  mean_abs_pn_phi_init = mean(abs(pn_phi_init)),
  mean_abs_pn_phi_tmle = mean(abs(pn_phi_tmle)),
  mean_abs_p0_phi_init = mean(abs(p0_phi_init)),  # just added
  mean_abs_p0_phi_tmle = mean(abs(p0_phi_tmle))  # just added
  
), by = .(scenario, n)]

summ[]
