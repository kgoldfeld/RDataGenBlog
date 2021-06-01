library(simstudy)
library(data.table)
library(distributional)
library(ggplot2)
library(bayesplot)
library(posterior)
library(ggdist)
library(cmdstanr)
library(ggpubr)
library(collapse)
library(paletteer)

s_define <- function() {
  
  defB <- defDataAdd(varname = "Y", formula = "-2 + ..lor * rx", 
                     dist = "binary", link="logit")
  
  return(list(defB = defB)) # list_of_defs is a list of simstudy data definitions
}

s_generate <- function(list_of_defs, argsvec) {
  
  list2env(list_of_defs, envir = environment())
  list2env(as.list(argsvec), envir = environment())
  
  #--- add data generation code ---#
  
  lor <- rnorm(1, mu.lor, sigma.lor)
  
  dT <- genData(nobs)
  dT <- trtAssign(dT, grpName = "rx")
  dT <- addColumns(defB, dT)
  
  return(dT[])
  
}

s_model <- function(generated_data, argsvec) {
  
  list2env(as.list(argsvec), envir = environment())
  
  dt_to_list <- function(dx) {
    
    N <- nrow(dx)                 ## number of observations 
    y <- dx$Y                      ## individual outcome 
    x <- dx$rx                     ## treatment arm for individual 
    s <- t_sigma
    mu <- 0
    
    list(N=N, y=y, x=x, s=s, mu = mu)
  }

  fit <- mod$sample(
    data = dt_to_list(generated_data),
    refresh = 0,
    chains = 4L,
    parallel_chains = 4L,
    iter_warmup = 500,
    iter_sampling = 2500,
    step_size = 0.1,
    show_messages = FALSE
  )
  
  res <- data.table(fit$summary(variables = "beta"))[, .(median, sd, q95, len = q95-q5)]
  
  draws_array <- as_draws_array(fit$draws())
  betas <- data.table(beta = as.matrix(draws_array[,,"beta"]))
  res$p0 <- mean(betas$beta.V1 < 0)

  d <- density(draws_array[,,"beta"])
  
  return(list(res, d)) # model_results is a data.table
}


s_single_rep <- function(list_of_defs, argsvec) {
  
  generated_data <- s_generate(list_of_defs, argsvec)
  model_results <- s_model(generated_data, argsvec)
  
  return(model_results)
}

s_replicate <- function(argsvec, nsim) {
  
  list_of_defs <- s_define()
  
  model_results <- 
    lapply(
      X = 1 : nsim, 
      FUN = function(x) s_single_rep(list_of_defs, argsvec)
  )
  
  #--- add summary statistics code ---#
  
  model_sums <- unlist2d(lapply(model_results, function(x) x[[1]]), idcols = "replicate", DT = TRUE)
  model_densities <- lapply(model_results, function(x) x[[2]])
  
  rm(model_results)
  
  summary_stats <- model_sums[, .(p95 = mean(p0 >= 0.95), estsd = mean(sd))]
  model_ests <- list(data.table(t(argsvec), summary_stats), model_densities, model_sums)
  # model_ests <- data.table(t(argsvec), summary_stats)
  
  return(model_ests) # summary_stats is a data.table
}

scenario_list <- function(...) {
  argmat <- expand.grid(...)
  return(asplit(argmat, MARGIN = 1))
}

mu.lor <- c(-0.5, -1, -1.5)
sigma.lor <- c(0.10)
nobs <- c(100, 200, 400)
t_sigma <- c(5)

scenarios <- scenario_list(mu.lor = mu.lor, sigma.lor = sigma.lor, 
                           nobs = nobs, t_sigma = t_sigma)


mod <- cmdstan_model("code/bayes_logistic.stan")


#--- run locally ---#

summary_stats <- lapply(scenarios, function(a) s_replicate(a, nsim = 20))

####

###
i <- 2

params <- summary_stats[[i]][[1]]
densities <- summary_stats[[i]][[2]]
dist_sum <- summary_stats[[i]][[3]]

title <- paste0("N: ", params$nobs, ", effect: ", params$mu.lor)
plot_points <- lapply(densities, function(x) as.data.table(x[c("x", "y")]))

dp <- unlist2d(plot_points, idcols = "replicate", DT=TRUE)
dp <- merge(dp, dist_sum, by = "replicate")

p <- ggplot(data = dp, aes(x = x, y = y, group = replicate )) +
  geom_line(aes(color = factor(p0 > 0.95))) +
  xlim(-10, 4) +
  # ylim(0, 1.5) +
  # ggtitle(title) +
  scale_color_manual(values = c("#e5b758", "#2561dd"),
                     name = "P(log(OR) < 0) > 95%") +    
  theme(panel.grid = element_blank(),
        plot.title = element_text(size=10)) +
  xlab("log(OR)") +
  ylab("density")

ggsave(file = "img/p95_single.png", height=3, width=6, scale = 1.5)

###


p <- list()
q <- list()
r <- list()

count <- 0

for (i in c(2, 5, 8)) {
#for (i in 1:length(summary_stats)) {
  
  count <- count + 1  
  
  params <- summary_stats[[i]][[1]]
  densities <- summary_stats[[i]][[2]]
  dist_sum <- summary_stats[[i]][[3]]
  
  title <- paste0("N: ", params$nobs)
  plot_points <- lapply(densities, function(x) as.data.table(x[c("x", "y")]))
  
  dp <- unlist2d(plot_points, idcols = "replicate", DT=TRUE)
  dp <- merge(dp, dist_sum, by = "replicate")
  
  p[[count]] <- ggplot(data = dp, aes(x = x, y = y, group = replicate )) +
    geom_line(aes(color = factor(p0 > 0.95))) +
    xlim(-4, 2) +
    ylim(0, 1.5) +
    ggtitle(title) +
    scale_color_manual(values = c("#e5b758", "#2561dd"),
                       name = "P(log(OR) < 0) > 95%") +
    theme(panel.grid = element_blank(),
      plot.title = element_text(size=10))
  
  q[[count]] <- ggplot(data = dp, aes(x = x, y = y, group = replicate )) +
    geom_line(aes(color = factor(len < 2))) +
    xlim(-4, 2) +
    ylim(0, 1.5) +
    ggtitle(title) +
    scale_color_paletteer_d(limits = c("FALSE","TRUE"),
      palette = "wesanderson::Moonrise2") +
    theme(panel.grid = element_blank(),
      plot.title = element_text(size = 10))
  
  r[[count]] <- ggplot(data = dp, aes(x = x, y = y, group = replicate )) +
    geom_line(aes(color = factor(sd < 0.5))) +
    xlim(-4, 2) +
    ylim(0, 1.5) +
    ggtitle(title) +
    scale_color_paletteer_d(limits = c("FALSE","TRUE"),
                            palette = "wesanderson::Moonrise2") +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 10))
  
}

ggsave(file = "img/p95.png",
       ggarrange(plotlist = p, common.legend = TRUE, legend = "bottom", nrow=1),
       width = 7, height = 3, scale=1.5)


ggarrange(plotlist = q, common.legend = TRUE, legend = "right")
ggarrange(plotlist = r, common.legend = TRUE, legend = "right")


ggsave("test.png")

psum <- rbindlist(lapply(summary_stats, function(x) x[[1]]))
psum



###



