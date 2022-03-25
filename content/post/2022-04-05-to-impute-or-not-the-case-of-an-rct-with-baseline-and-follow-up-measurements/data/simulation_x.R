library(simstudy)
library(data.table)
library(ggplot2)
library(mice)

s_define <- function() {
  
  def <- defData(varname = "y0", formula = 0, variance = 1)
  def <- defData(def, "a", formula = "1;1", dist = "trtAssign")
  def <- defData(def, "y1", 
                 formula = "y0 + ..delta * a - ..lambda * y0 * a", variance = 0.5)
  def <- defData(def, "y1_obs", formula = "y1", dist = "nonrandom")
  
  defM <- defMiss(
    varname = "y1_obs", formula = "-1.5  - ..alpha * y0 - ..beta * a - ..zeta * a *y0", 
    logit.link = TRUE
  )
  
  return(list(def = def, defM = defM))
}

s_generate <- function(list_of_defs, argsvec) {
  
  list2env(list_of_defs, envir = environment())
  list2env(as.list(argsvec), envir = environment())
  
  dd <- genData(200, def)
  dmiss <- genMiss(dd, defM, idvars = "id")
  dobs <- genObs(dd, dmiss, idvars = "id")
  
  return(dobs) #  generated_data is a data.table
}

s_model <- function(generated_data) {
  
  imp_dd <- generated_data[, -c("id", "y1")]
  imp <- mice(imp_dd, m=20, maxit=5, print=FALSE)
  
  ###
  
  a_all <- coef(lm(y1 ~ y0 + a, data = generated_data))["a"]
  a_missing <- coef(lm(y1_obs ~ y0 + a, data = generated_data))["a"]
  
  fit_imp <- with(imp, lm(y1_obs ~ y0 + a))
  a_imp <- summary(pool(fit_imp))[3, "estimate"]
  
  return(data.table(a_all, a_missing, a_imp)) # model_results is a data.table
}

s_single_rep <- function(list_of_defs, argsvec) {
  
  generated_data <- s_generate(list_of_defs, argsvec)
  model_results <- s_model(generated_data)
  
  return(model_results)
}

s_replicate <- function(argsvec, nsim) {
  
  list_of_defs <- s_define()
  
  model_results <- rbindlist(
    parallel::mclapply(
      X = 1 : nsim, 
      FUN = function(x) s_single_rep(list_of_defs, argsvec), 
      mc.cores = 4)
  )
  
  #--- add summary statistics code ---#
  
  summary_stats <- model_results[, .(
    mean_all = mean(a_all, na.rm = TRUE), 
    bias_all = mean(a_all - delta, na.rm = TRUE), 
    var_all = var(a_all, na.rm = TRUE), 
    
    mean_missing = mean(a_missing, na.rm = TRUE), 
    bias_missing = mean(a_missing - delta, na.rm = TRUE), 
    var_missing = var(a_missing, na.rm = TRUE),
    
    mean_imp = mean(a_imp, na.rm = TRUE), 
    bias_imp = mean(a_imp - delta, na.rm = TRUE), 
    var_imp = var(a_imp, na.rm = TRUE)
  )]
  
  summary_stats <- data.table(t(argsvec), summary_stats)
  
  return(summary_stats) # summary_stats is a data.table
}

#---- specify varying power-related parameters ---#

scenario_list <- function(...) {
  argmat <- expand.grid(...)
  return(asplit(argmat, MARGIN = 1))
}

delta <- 1
lambda <- c(0, 0.2, .4, .6, .8, 1)
alpha <- c(0, 0.5, 1)
beta <- c(0)
zeta <- c(0, 0.5, 1)

scenarios <- scenario_list(delta = delta, lambda = lambda, alpha = alpha, beta = beta, zeta = zeta)

summary_stats <- rbindlist(lapply(scenarios, function(a) s_replicate(a, nsim = 500)))

save(summary_stats, file = "data/summary_stats_x.rdata")


# save(summary_stats, file = "data/summary_stats_adj.rdata")

summary_stats[, mse_all := bias_all^2 + var_all]
summary_stats[, mse_missing := bias_missing^2 + var_missing]
summary_stats[, mse_imp := bias_imp^2 + var_imp]

dp <- melt(summary_stats, 
           measure.vars = c("bias_all", "bias_imp", "bias_missing"))
dp$variable <- factor(dp$variable, 
                      labels = c("No missing data", "Imputed data", "Complete cases only"))

ggplot(data = dp, aes(x = lambda, y = value  )) +
  geom_line(aes(group = variable, color=factor(variable))) +
  facet_grid(zeta ~ alpha) +
  scale_y_continuous(limits = c(-.25, .1), 
                     name = "bias", breaks = c(-.2, -.1, 0),
                     sec.axis = sec_axis(~ . , 
                                         name = expression(paste("missingness based on treatment and ", Y[0], " ", (zeta))), 
                                         breaks = NULL, labels = NULL)) +
  scale_x_continuous(
    name = expression(paste("interaction of treatment and ", Y[0], " ", (lambda))), 
    breaks = c(0, .2, .4, .6, .8, 1 ),
    sec.axis = sec_axis(~ . , 
                        name = expression(paste("missingness based on ", Y[0], " ", (alpha))), 
                        breaks = NULL, labels = NULL)) +
  scale_color_brewer(palette = "Dark2", ) +
  ggtitle(expression(paste("Bias of treatment effect estimate ", (hat(delta)) ))) +
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        panel.spacing = unit(1, "lines")) 

