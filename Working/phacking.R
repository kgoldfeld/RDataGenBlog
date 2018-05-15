library(data.table)
library(ggplot2)


phack <- function(nvars, ssize, plot = FALSE) {
  
  dt <- genCorData(n = ssize, mu = rep(0, nvars), sigma = 1, corstr = "ind")
  dt <- trtAssign(dt, nTrt = 2, grpName = "a")
  
  result <- data.table()
  
  for ( i in 1:nvars) {
    formula <- paste0("V",i, " ~ a")
    lmfit <- lm(formula, data = dt)
    pvalue <-  coef(summary(lmfit))["a", "Pr(>|t|)"]
    est <- coef(summary(lmfit))["a", "Estimate"]
    se <- coef(summary(lmfit))["a", "Std. Error"]
    
    resulti <- data.table(i, est, l95 = est - 1.96 * se, u95 = est + 1.96 * se, 
                          pvalue, sig = (pvalue<=0.025))
    
    result <- rbind(result, resulti)
  
  }
  
  falsesigs <- result[, sum(sig)]
  
  if (plot) {
    hackplot <- ggplot(data = result, aes(y = est, x = i)) +
      geom_hline(yintercept = 0, color = "white", size = 1) +
      geom_point(aes(color=sig), size = 1, show.legend = FALSE) +
      geom_errorbar(aes(ymin = l95, ymax = u95, color = sig),
                    size = .1, width = 0, show.legend = FALSE) +
      scale_color_manual(values=c("black", "red")) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      scale_y_continuous(name = "Estimate", limits = c(-0.5, 0.5), 
                         breaks = seq(-0.5, 0.5, 0.25)) +
      scale_x_continuous(name = "Variable")
  } else {
    hackplot <- NULL
  }
  
 
  return(list(falseSig = falsesigs, hackplot = hackplot))
  
}

ps <- phack(20, 500, plot = TRUE)

result <- rep(NA, 1000)

for (j in 1:2000) {
  ps <- phack(2, 500)
  result <- c(result, ps$falseSig)
}

prop.table(table(result))
