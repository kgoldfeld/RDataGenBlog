library(rstan)
library(simstudy)
library(bayesplot)

rt_c <- stanc("content/post/DataBayesMeta2/binary_center.stan")
sm_c <- stan_model(stanc_ret = rt_c, verbose=FALSE)

rt_nc <- stanc("content/post/DataBayesMeta2/binary_noncenter.stan")
sm_nc <- stan_model(stanc_ret = rt_nc, verbose=FALSE)

def_s <- defDataAdd(varname = "b0", formula = 0, variance = 0.025)
def_s <- defDataAdd(
  def_s, varname = "delta_k", 
  formula = "(c_type==1) * 0.4 + (c_type==2) * 0.5 + (c_type==3) * 0.6", 
  dist = "nonrandom"
)

def_i <- defDataAdd(
  varname = "y", formula = "-1 + b0 + rx * delta_k", 
  dist = "binary", link = "logit")

set.seed(19827)

dc <- genData(3, id = "c_type")

ds <- genCluster(dc, "c_type", numIndsVar = 7, level1ID = "site")
ds <- addColumns(def_s, ds)
ds[, order := .SD[, .(order = .I)]$order, keyby = c_type] # for plotting

di <- genCluster(ds, "site", 50, "id")
di <- trtAssign(di, 2, strata = "site", grp = "rx")
di <- addColumns(def_i, di)

N <- nrow(di) ;                       
C <- di[, length(unique(c_type))]     
K <- di[, length(unique(site))]       
y <- as.numeric(di$y)                 
kk <- di$site                         
ctrl <- di$rx                         
cc <- di[, .N, keyby = .(site, c_type)]$c_type  

sampdat <- list(N=N, C=C, K=K, y=y, kk=kk, ctrl=ctrl, cc=cc)

fit_c <-  sampling(
  sm_c, data = sampdat, iter = 3000, warmup = 500, 
  show_messages = FALSE, cores = 4, refresh = 0,
  control = list(adapt_delta = 0.8), seed = 27162
)

fit_nc <-  sampling(
  sm_nc, data = sampdat, iter = 3000, warmup = 500, 
  show_messages = FALSE, cores = 4, refresh = 0,
  control = list(adapt_delta = 0.8), seed = 27162
)


#---- Plots

posterior_c <- as.array(fit_c) 
lp_c <- log_posterior(fit_c)
np_c <- nuts_params(fit_c)

color_scheme_set("mix-brightblue-gray")

eta0_c <- mcmc_trace(posterior_c, pars = "eta_0", np = np_c) + 
  xlab("Post-warmup iteration")

eta0_c_window <- mcmc_trace(posterior_c, pars = "eta_0", np = np_c, window = c(700, 900)) + 
  xlab("Post-warmup iteration")

Delta_c <- mcmc_trace(posterior_c, pars = "Delta", np = np_c) + 
  xlab("Post-warmup iteration")

Delta_c_window <- mcmc_trace(posterior_c, pars = "Delta", np = np_c, window = c(700, 900)) + 
  xlab("Post-warmup iteration")

trace_c <- gridExtra::grid.arrange(Delta_c, Delta_c_window, eta0_c, eta0_c_window,
                        nrow = 2)

ggsave(trace_c, filename = "static/img/post-bayesdiag/trace_c.png", 
       width = 9, height = 6, scale = 1.4)

posterior_nc <- as.array(fit_nc) 
lp_nc <- log_posterior(fit_nc)
np_nc <- nuts_params(fit_nc)

eta0_nc <- mcmc_trace(posterior_nc, pars = "eta_0", np = np_nc) + 
  xlab("Post-warmup iteration")

eta0_nc_window <- mcmc_trace(posterior_nc, pars = "eta_0", np = np_nc, window = c(700, 900)) + 
  xlab("Post-warmup iteration")

Delta_nc <- mcmc_trace(posterior_nc, pars = "Delta", np = np_nc) + 
  xlab("Post-warmup iteration")

Delta_nc_window <- mcmc_trace(posterior_nc, pars = "Delta", np = np_nc, window = c(700, 900)) + 
  xlab("Post-warmup iteration")

trace_nc <- gridExtra::grid.arrange(Delta_nc, Delta_nc_window, eta0_nc, eta0_nc_window,
                                   nrow = 2)

ggsave(trace_nc, filename = "static/img/post-bayesdiag/trace_nc.png", 
       width = 9, height = 6, scale = 1.4)


color_scheme_set("darkgray")

parcoord_c <-mcmc_parcoord(posterior_c, np = np_c, 
                           pars = c("eta_0", "sigma_b", "alpha", 
                                    "delta_c[1]", "delta_c[2]",
                                    "delta_c[3]", "Delta"))

parc_c <- parcoord_c +
  scale_x_discrete(expand = c(.01,.01)) +
  theme(panel.background = element_rect(fill = "grey90")) +
  ylim(-3, 3) +
  ggtitle("Original model specification")

ggsave(parc_c, filename = "static/img/post-bayesdiag/parc_c.png", 
       width = 7, height = 4, scale = 1.3)

parcoord_nc <-mcmc_parcoord(posterior_nc, np = np_nc, 
                           pars = c("eta_0", "sigma_b", "alpha", 
                                    "delta_c[1]", "delta_c[2]",
                                    "delta_c[3]", "Delta"))

parc_nc <- parcoord_nc +
  scale_x_discrete(expand = c(.01,.01)) +
  theme(panel.background = element_rect(fill = "grey90")) +
  ylim(-3, 3) +
  ggtitle("Original model specification")

ggsave(parc_nc, filename = "static/img/post-bayesdiag/parc_nc.png", 
       width = 7, height = 4, scale = 1.3)

