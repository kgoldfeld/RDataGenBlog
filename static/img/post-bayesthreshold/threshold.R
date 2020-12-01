library(simstudy)
library(bayesplot)
library(rstan)
library(glue)
library(gridExtra)
library(Rmpfr)

# Generate simulated data

d1 <- defData(varname = "antibody", formula = 0, variance = 1, dist = "normal")
d1 <- defData(d1, varname = "latent_status", formula = "-3 + 6 * (antibody > -0.7)",
              dist = "binary", link = "logit")
d1 <- defData(d1, varname = "y", formula = "0 + ..effect_size * latent_status", 
              variance = 1, dist = "normal")

d1 <- defData(varname = "antibody", formula = 0, variance = 1, dist = "normal")  
d1 <- defData(d1, varname = "y", formula = "antibody", variance = 1, dist = "normal")

dd <- genData(500, d1)

set.seed(387263)

effect_size = 3
dd3 <- genData(500, d1)

effect_size = 0
dd0 <- genData(500, d1)

pobs <- ggplot(data = dd3, aes(x = antibody, y = y)) +
  geom_point() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  ggtitle("Observed data") +
  scale_x_continuous(limits = c(-3.5, 3.5))

platent <- ggplot(data = dd3, aes(x = antibody, y = y)) +
  geom_point(aes(color = factor(latent_status))) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold")) +
  scale_color_manual(values = c("#bb5566", "#ddaa33")) +
  ggtitle("Latent threshold") +
  scale_x_continuous(limits = c(-3.5, 3.5))

p3 <- arrangeGrob(pobs, platent, nrow = 1)
grid.arrange(p3)

ggsave(p3, file = "static/img/post-bayesthreshold/p3.png", 
       height = 3, width = 6, scale = 1.75)

pobs <- ggplot(data = dd0, aes(x = antibody, y = y)) +
  geom_point() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  ggtitle("Observed data") +
  scale_x_continuous(limits = c(-3.5, 3.5))

platent <- ggplot(data = dd0, aes(x = antibody, y = y)) +
  geom_point(aes(color = factor(latent_status))) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold")) +
  scale_color_manual(values = c("#bb5566", "#ddaa33")) +
  ggtitle("Latent threshold") +
  scale_x_continuous(limits = c(-3.5, 3.5))

p0 <- arrangeGrob(pobs, platent, nrow = 1)
grid.arrange(p0)

ggsave(p0, file = "static/img/post-bayesthreshold/p0.png", 
       height = 3, width = 6, scale = 1.75)

# Fit with stan

rt <- stanc("static/img/post-bayesthreshold/threshold.stan");
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

N <- nrow(dd3)
y <- dd3[, y]
x <- dd3[, antibody] 

c <- seq(round(min(x), 1), round(max(x), 1), by = .1)
M <- length(c)

studydata3 <- list(N=N, x=x, y=y, M=M, c=c)
fit3 <-  sampling(sm, data = studydata3, iter = 3000, warmup = 500, 
                  cores = 4L, chains = 4, control = list(adapt_delta = 0.8))

posterior <- as.array(fit3) 
lp <- log_posterior(fit3)
np <- nuts_params(fit3)

color_scheme_set("mix-brightblue-gray")

p <- mcmc_trace(posterior, pars = c("alpha","beta", "sigma"), 
                facet_args = list(nrow = 3), np = np) + 
  xlab("Post-warmup iteration")
p

ggsave(p, file = "static/img/post-bayesthreshold/trace3.png", 
       height = 3, width = 4, scale = 1.75)

p <- mcmc_intervals(posterior, pars = c("alpha","beta", "sigma"), prob_outer = 0.95)
p + scale_x_continuous(limits = c(0, 3))

ggsave(p, file = "static/img/post-bayesthreshold/estimates3.png", 
       height = 2, width = 4, scale = 1.75)

a <- mpfr(exp(-100), precBits=64)

qs <- NULL
for(m in 1:M) {
  lp.i <- glue("lp[{m}]")
  le <- rstan::extract(fit3, pars = lp.i)[[1]]
  q <- a^(-le/100)
  qs[m] <- sum(q)
}

qss <- mpfr2array(qs, dim = M)
ps <- log(qss/sum(qss))
dps <- data.table(c, y=as.numeric(ps))

p <- ggplot(data = dps, aes(x = c, y = y)) +
  geom_vline(xintercept = -0.7, color = "red", lty = 3) +
  geom_line(color = "grey60") +
  geom_point(size = 1) +
  theme(panel.grid = element_blank()) +
  ylab("log(probability)") +
  xlab("threshold from low to not low") +
  scale_y_continuous(limits = c(-800, 0))

p
ggsave(p, file = "static/img/post-bayesthreshold/threshold3.png", 
       height = 2, width = 4, scale = 2.00)

###

rt <- stanc("static/img/post-bayesthreshold/no_threshold.stan");
sm_no <- stan_model(stanc_ret = rt, verbose=FALSE)

N <- nrow(dd3)
t <- dd3[, as.integer(antibody >= -0.7)]
y <- dd3[, y]

studydata3no <- list(N=N, t=t, y=y)
fit3no <-  sampling(sm_no, data = studydata3no, iter = 3000, warmup = 500, 
                  cores = 4L, chains = 4, control = list(adapt_delta = 0.8))

posterior <- as.array(fit3no) 
lp <- log_posterior(fit3no)
np <- nuts_params(fit3no)

color_scheme_set("mix-brightblue-gray")

p <- mcmc_trace(posterior, pars = c("alpha","beta", "sigma"), 
                facet_args = list(nrow = 3), np = np) + 
  xlab("Post-warmup iteration")
p

ggsave(p, file = "static/img/post-bayesthreshold/trace3.png", 
       height = 3, width = 4, scale = 1.75)

p <- mcmc_intervals(posterior, pars = c("alpha","beta", "sigma"), prob_outer = 0.95)
p + scale_x_continuous(limits = c(0, 3))

###

N <- nrow(dd0)
y <- dd0[, y]
x <- dd0[, antibody] 

c <- seq(round(min(x), 1), round(max(x), 1), by = .1)
M <- length(c)

studydata0 <- list(N=N, x=x, y=y, M=M, c=c)
fit0 <-  sampling(sm, data = studydata0, iter = 3000, warmup = 500, 
                  cores = 4L, chains = 4, control = list(adapt_delta = 0.8))


posterior <- as.array(fit0) 
lp <- log_posterior(fit0)
np <- nuts_params(fit0)

color_scheme_set("mix-brightblue-gray")

p <- mcmc_trace(posterior, pars = c("alpha","beta", "sigma"), 
                facet_args = list(nrow = 3), np = np) + 
  xlab("Post-warmup iteration")
p

ggsave(p, file = "static/img/post-bayesthreshold/trace0.png", 
       height = 3, width = 4, scale = 1.75)

p <- mcmc_intervals(posterior, pars = c("alpha","beta", "sigma"))
p

ggsave(p, file = "static/img/post-bayesthreshold/estimates0.png", 
       height = 2, width = 4, scale = 1.75)

a <- mpfr(exp(-100), precBits=64)

qs <- NULL
for(m in 1:M) {
  lp.i <- glue("lp[{m}]")
  le <- rstan::extract(fit0, pars = lp.i)[[1]]
  q <- a^(-le/100)
  qs[m] <- sum(q)
}

qss <- mpfr2array(qs, dim = M)
ps <- log(qss/sum(qss))
dps <- data.table(c, y=as.numeric(ps))

p <- ggplot(data = dps, aes(x = c, y = y)) +
  geom_vline(xintercept = -0.7, color = "red", lty = 3) +
  geom_line(color = "grey60") +
  geom_point(size = 1) +
  theme(panel.grid = element_blank()) +
  ylab("log(probability)") +
  xlab("threshold from low to not low") +
  scale_y_continuous(limits = c(-800, 0))

p
ggsave(p, file = "static/img/post-bayesthreshold/threshold0.png", 
       height = 2, width = 4, scale = 2.00)
###




