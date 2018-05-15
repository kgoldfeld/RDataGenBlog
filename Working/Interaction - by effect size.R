library(simstudy)
library(scales)
library(ggplot2)

def <- defData(varname = "disease", formula = .5, dist = "binary")

genFormula <- function(coefs, vars) {
  
  lcoef <- length(coefs)
  lvars <- length(vars)
  
  if ( !(lcoef == lvars | lcoef == lvars + 1) ) {
    stop("Coefficients or variables not properly specified")
  }
  
  if (lcoef != lvars) { # Intercept
    form <- paste0(coefs[1])
    coefs <- coefs[-1]
    lcoef <- lcoef - 1
  }
  
  for (i in 1:lcoef) {
    form <- paste(form, "+" , coefs[i], "*", vars[i])
  }
  
  return(form)
}

pvals <- data.table()

for (eff0 in seq(-0.5, 0.5, by = 0.10)) {
  for (eff1 in seq(-0.5, 0.5, by = 0.25)) {
    
    form0 <- genFormula(c(0.0, eff0), c("trt"))
    form1 <- genFormula(c(0.5, eff1), c("trt"))
    
    def2 <- defCondition(condition = "disease == 0", 
                         formula = form0, variance = 1,
                         dist = "normal")
    
    def2 <- defCondition(def2, condition = "disease == 1", 
                         formula = form1, variance = 1,
                         dist = "normal")
    
    for (i in 1: 1500) {
      
      dx <- genData(400, def)
      dx <- trtAssign(dx, nTrt = 2, balanced = TRUE, 
                      strata = "disease", grpName = "trt")
      
      
      dx <- addCondition(def2, dx, "y")
      
      lmMain <- lm(y ~ disease + trt, data = dx)
      lmInt <- lm(y ~ disease + trt + trt*disease, data = dx)
      lmDis <- lm(y ~ disease, data = dx)
      
      cM <- coef(summary(lmMain))["trt", 4]
      cI <- coef(summary(lmInt))["disease:trt", 4]
      fI <- anova(lmDis, lmInt)$`Pr(>F)`[2]
      
      pvals <- rbind(pvals, data.table(eff0, eff1, cM, cI, fI))
      
    }
  }
}

pvals[, mEffect := (cM <= .05)]

method1 <- pvals[, .(power = mean(mEffect), approach = "1"), 
              keyby = .(eff0, eff1)]

pvals[, mEffect := (cM <= .025)]
pvals[, iEffect := (cI <= .025)]

method2 <- pvals[, .(power = mean(iEffect) + mean((!iEffect) & mEffect), 
                  approach = "2"), 
              keyby = .(eff0, eff1)]


pvals[, fEffect := (fI <= 0.05)]
pvals[, iEffect := (cI <= .05)]
pvals[, mEffect := (cM <= .05)]

method3 <- pvals[, .(power = mean(fEffect & iEffect) + 
                         mean(fEffect & !(iEffect) & mEffect), 
                     approach = "3"), 
                 keyby = .(eff0, eff1)]

pow <- rbind(method1, method2, method3)

a_label <- function(textnum) {
  mynum <- as.numeric(textnum)
  paste0("Group 2 effect: ", sprintf("%1.2f", mynum))
}

p <- ggplot(data = pow[eff1 %in% c(-.5, 0, .5)], 
            aes(x=eff0, y = power, group = approach)) +
  geom_smooth(aes(color = approach), se = FALSE) +
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0.2, 0.8), name = "Power", 
                     label = percent) +
  xlab("Group 1 effect") +
  scale_color_manual(values=c("#7e42bc","#bcbb42", "#4280bc"), 
                     name = "First model",
                     labels = c("Main effects first", "Interaction first", "Global test")) +
  facet_grid(eff1 ~ ., labeller = labeller(eff1 = a_label)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.title = element_blank(),
        strip.text = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 12)
        )

png(filename = "static/img/post-interaction/Models.png", width = 600, height = 600)  
 p
dev.off()
