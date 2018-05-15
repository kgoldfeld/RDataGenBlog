library(simstudy)

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

form0 <- genFormula(c(0.0, 0), c("trt"))
form1 <- genFormula(c(0.5, 0), c("trt"))

def2 <- defCondition(condition = "disease == 0", 
                     formula = form0, variance = 1,
                     dist = "normal")

def2 <- defCondition(def2, condition = "disease == 1", 
                     formula = form1, variance = 1,
                     dist = "normal")

pvals <- data.table()
for (i in 1: 10000) {
  
  dx <- genData(400, def)
  dx <- trtAssign(dx, nTrt = 2, balanced = TRUE, 
                  strata = "disease", grpName = "trt")
  
  
  dx <- addCondition(def2, dx, "y")
  
  lmMain <- lm(y ~ disease + trt, data = dx)
  lmInt <- lm(y ~ disease + trt + trt*disease, data = dx)
  lmDis <- lm(y ~ disease, data = dx)
  
  cM <- coef(summary(lmMain))["trt", 4]
  cI <- coef(summary(lmInt))["disease:trt", 4]
  cIm <- coef(summary(lmInt))["trt", 4]
  
  fI <- summary(lmInt)$fstatistic
  fI <- pf(fI[1], fI[2], fI[3], lower.tail = FALSE)
  fDm <- anova(lmDis, lmInt)$`Pr(>F)`[2]
  
  pvals <- rbind(pvals, data.table(cM, cI, cIm, fI, fDm))
  
  # if (fDm > .05 & fI <= .05) break
  # if (fDm < .05 & cI > 0.05 & cIm > 0.05) break
  
}

pvals[, mEffect := (cM <= .05)]
pvals[, iEffect := (cI <= .05)]

pvals[, mean(mEffect & !iEffect)] + pvals[, mean(mEffect & iEffect)]
pvals[, mean(mEffect)]

pvals[, mEffect := (cM <= .025)]
pvals[, iEffect := (cI <= .025)]

pvals[, mean(iEffect)] + pvals[, mean((!iEffect) & mEffect)]

pvals[, fMeffect := (fDm <= 0.05)]
pvals[, iEffect := (cI <= .05)]
pvals[, mEffect := (cM <= .05)]

pvals[, mean(fMeffect & iEffect)] + pvals[, mean(fMeffect & !(iEffect) & mEffect)]
