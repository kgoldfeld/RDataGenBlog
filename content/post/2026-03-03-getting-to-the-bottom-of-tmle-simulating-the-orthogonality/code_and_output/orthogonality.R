library(simstudy)
library(data.table)
library(ggplot2)
library(latex2exp)

gen_dgp <- function(n) {
  
  def <- 
    defData(varname = "x1", formula = .5, dist = "binary") |>
    defData(varname = "x2", formula = 0, variance = 1) |>
    defData(varname = "a", formula = "-0.2 + 0.8 * x1 + 0.6 * x2", dist = "binary", link="logit") |>
    defData(varname = "y", formula = "..tau * a + 1.0 * x1 + 1.0 * x2 + 1.5 * x1 * x2", variance = 1)
  
  genData(n, def)[]
  
}

fit_nuisance <- function(dt, scenario) {
  
  # Outcome regression Q(a,x)
  
  if (scenario %in% c("both_correct", "g_wrong")) {
    Q_fit <- lm(y ~ a + x1 + x2 + x1:x2, data = dt)  # correct
  } else {
    Q_fit <- lm(y ~ a + x1, data = dt)               # wrong on purpose
  }
  
  # Propensity model g(x)
  
  if (scenario %in% c("both_correct", "Q_wrong")) {
    g_fit <- glm(a ~ x1 + x2, data = dt, family = binomial())  # correct
  } else {
    g_fit <- glm(a ~ x1, data = dt, family = binomial())       # wrong on purpose
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
  pmin(pmax(p, 0.01), 0.99)  # simple stabilization
}

Q_true <- function(dt, a_val, tau) {
  tau * a_val + 1.0 * dt$x1 + 1.0 * dt$x2 + 1.5 * dt$x1 * dt$x2
}

g_true <- function(dt) {
  plogis(-0.2 + 0.8 * dt$x1 + 0.6 * dt$x2)
}

# EIF evaluator given Q1,Q0,g,psi
phi_ate <- function(dt, Q1, Q0, g, psi) {
  A <- dt$a
  Y <- dt$y
  (Q1 - Q0 - psi) + A/g * (Y - Q1) - (1 - A)/(1 - g) * (Y - Q0)
}

# Build phi_hat using your fitted nuisances (AIPW-style)
phi_hat_from_fits <- function(dt, Q_fit, g_fit, psi_hat) {
  Q1 <- predict_Q(Q_fit, dt, 1)
  Q0 <- predict_Q(Q_fit, dt, 0)
  g  <- predict_g(g_fit, dt)
  phi_ate(dt, Q1, Q0, g, psi_hat)
}

# Build phi0 from true nuisances
phi0_true <- function(dt, tau) {
  Q1 <- Q_true(dt, 1, tau)
  Q0 <- Q_true(dt, 0, tau)
  g  <- g_true(dt)
  psi0 <- tau
  phi_ate(dt, Q1, Q0, g, psi0)
}

# Cross-fitted AIPW estimate (just as a convenient psi_hat)
psi_hat_aipw_from_fits <- function(dt, Q_fit, g_fit) {
  Q1 <- predict_Q(Q_fit, dt, 1)
  Q0 <- predict_Q(Q_fit, dt, 0)
  g  <- predict_g(g_fit, dt)
  A <- dt$a; Y <- dt$y
  mean((Q1 - Q0) + A/g * (Y - Q1) - (1 - A)/(1 - g) * (Y - Q0))
}

est_2T <- function(scenario, dd, tau, dd_pop) {
  
  n <- nrow(dd)
  idx <- sample.int(n)
  I1 <- idx[1:floor(n/2)]
  I2 <- idx[(floor(n/2)+1):n]
  
  # fit nuisances on each training fold
  fits1 <- fit_nuisance(dd[I1], scenario)  # trained on fold 1
  fits2 <- fit_nuisance(dd[I2], scenario)  # trained on fold 2
  
  # cross-fitted psi_hat (evaluate each model on opposite fold, then average)
  psi1 <- psi_hat_aipw_from_fits(dd[I2], fits1$Q_fit, fits1$g_fit)  # train 1, eval 2
  psi2 <- psi_hat_aipw_from_fits(dd[I1], fits2$Q_fit, fits2$g_fit)  # train 2, eval 1
  psi_hat_cf <- 0.5 * (psi1 + psi2)
  
  # cross-fitted phi_hat on dd:
  # - for obs in fold 2, use fits1 (trained on fold 1)
  # - for obs in fold 1, use fits2 (trained on fold 2)
  phi_hat_dd <- numeric(n)
  phi_hat_dd[I2] <- phi_hat_from_fits(dd[I2], fits1$Q_fit, fits1$g_fit, psi_hat_cf)
  phi_hat_dd[I1] <- phi_hat_from_fits(dd[I1], fits2$Q_fit, fits2$g_fit, psi_hat_cf)
  
  dphi_dd <- phi_hat_dd - phi0_true(dd, tau)
  
  # approximate P0 expectation:
  # evaluate delta-phi under each fold-specific nuisance fit on independent pop,
  # then average them (since cross-fitting produces two fitted nuisance models)
  
  dphi_pop_1 <- phi_hat_from_fits(pop_dd, fits1$Q_fit, fits1$g_fit, psi_hat_cf) - phi0_true(pop_dd, tau)
  dphi_pop_2 <- phi_hat_from_fits(pop_dd, fits2$Q_fit, fits2$g_fit, psi_hat_cf) - phi0_true(pop_dd, tau)
  dphi_pop   <- 0.5 * (dphi_pop_1 + dphi_pop_2)
  
  T2 <- mean(dphi_dd) - mean(dphi_pop)
  
  data.table(scenario, n, T2)[]
}

run_sim <- function(n, tau, dd_pop, scenarios) {
  
  dd <- gen_dgp(n)
  rbindlist(lapply(
      scenarios, 
      function(s) est_2T(s, dd, tau, dd_pop)
    )
  )
  
}

# Example

set.seed(1)
tau <- 5
pop_dd <- gen_dgp(5e5)

n = rep(c(100, 250, 750, 1000, 2000), each = 500)
# scenarios = c("both_correct", "Q_wrong", "g_wrong", "both_wrong")
scenarios = c("both_correct", "Q_wrong")

res <- rbindlist(
  lapply(n, function(x) run_sim(x, tau, dd_pop, scenarios))
)

ressum <- res[, .(avg = round(mean(T2), 3), sd = round(sd(T2), 3), se=round(sd(T2)/sqrt(.N), 3)), keyby = .(scenario, n)]

###

p1 <- ggplot(data = res, aes(y = T2, x = factor(n))) +
  geom_hline(yintercept = 0, color = "white") +
  geom_jitter(size = .2, width = .2, color = "blue") +
  facet_wrap(~scenario,
             labeller = labeller(
               scenario = c(
                 "both_correct"    = "Both models correct",
                 "Q_wrong" = "Outcome (Q) model wrong"
               ))) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold")) +
  ylim(-5, 5) +
  xlab("sample size (n)") +
  ylab(TeX("$(P_n - P_0)(\\phi_{\\hat{P}} - \\phi_{P_0})$"))

ggsave(
  filename = "/Users/KSG/git R projects/rdatagen/content/post/2026-02-26-getting-to-the-bottom-of-tmle-simulating-the-orthogonality/code_and_output/ortho.png",
  height = 3, width = 7,
  scale = 1.2
)

