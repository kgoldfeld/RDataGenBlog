function (formula = formula(data), id = id, data = parent.frame(), 
          subset, na.action, R = NULL, b = NULL, tol = 0.001, maxiter = 25, 
          family = gaussian, corstr = "independence", Mv = 1, silent = TRUE, 
          contrasts = NULL, scale.fix = FALSE, scale.value = 1, v4.4compat = FALSE) 
{
  message("Beginning Cgee S-function, @(#) geeformula.q 4.13 98/01/27")
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$R <- m$b <- m$tol <- m$maxiter <- m$link <- m$varfun <- m$corstr <- m$Mv <- m$silent <- m$contrasts <- m$family <- m$scale.fix <- m$scale.value <- m$v4.4compat <- NULL
  if (is.null(m$id)) 
    m$id <- as.name("id")
  if (!is.null(m$na.action) && m$na.action != "na.omit") {
    warning("Only 'na.omit' is implemented for gee\ncontinuing with 'na.action=na.omit'")
    m$na.action <- as.name("na.omit")
  }
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  y <- as.matrix(model.extract(m, "response"))
  x <- model.matrix(Terms, m, contrasts)
  Q <- qr(x)
  if (Q$rank < ncol(x)) 
    stop("rank-deficient model matrix")
  N <- rep(1, length(y))
  if (dim(y)[2] == 2) {
    N <- as.vector(y %*% c(1, 1))
    y <- y[, 1]
  }
  else {
    if (dim(y)[2] > 2) 
      stop("Only binomial response matrices (2 columns)")
  }
  offset <- model.extract(m, offset)
  id <- model.extract(m, id)
  if (is.null(id)) {
    stop("Id variable not found")
  }
  nobs <- nrow(x)
  p <- ncol(x)
  xnames <- dimnames(x)[[2]]
  if (is.null(xnames)) {
    xnames <- paste("x", 1:p, sep = "")
    dimnames(x) <- list(NULL, xnames)
  }
  if (is.character(family)) 
    family <- get(family)
  if (is.function(family)) 
    family <- family()
  if (!is.null(b)) {
    beta <- matrix(as.double(b), ncol = 1)
    if (nrow(beta) != p) {
      stop("Dim beta != ncol(x)")
    }
    message("user's initial regression estimate")
    print(beta)
  }
  else {
    message("running glm to get initial regression estimate")
    mm <- match.call(expand.dots = FALSE)
    mm$R <- mm$b <- mm$tol <- mm$maxiter <- mm$link <- mm$varfun <- mm$corstr <- mm$Mv <- mm$silent <- mm$contrasts <- mm$scale.fix <- mm$scale.value <- mm$id <- NULL
    mm[[1]] <- as.name("glm")
    beta <- eval(mm, parent.frame())$coef
    print(beta)
    beta <- as.numeric(beta)
  }
  if (length(id) != length(y)) 
    stop("Id and y not same length")
  maxclsz <- as.integer(max(unlist(lapply(split(id, id), "length"))))
  maxiter <- as.integer(maxiter)
  silent <- as.integer(silent)
  if (length(offset) <= 1) 
    offset <- rep(0, length(y))
  if (length(offset) != length(y)) 
    stop("offset and y not same length")
  offset <- as.double(offset)
  if (!is.null(R)) {
    Rr <- nrow(R)
    if (Rr != ncol(R)) 
      stop("R is not square!")
    if (Rr < maxclsz) 
      stop("R not big enough to accommodate some clusters.")
  }
  else {
    R <- matrix(as.double(rep(0, maxclsz * maxclsz)), nrow = maxclsz)
  }
  links <- c("identity", "log", "logit", "inverse", "probit", 
             "cloglog")
  fams <- c("gaussian", "poisson", "binomial", "Gamma", "quasi")
  varfuns <- c("constant", "mu", "mu(1-mu)", "mu^2")
  corstrs <- c("independence", "fixed", "stat_M_dep", "non_stat_M_dep", 
               "exchangeable", "AR-M", "unstructured")
  linkv <- as.integer(match(c(family$link), links, -1))
  famv <- match(family$family, fams, -1)
  if (famv < 1) 
    stop("unknown family")
  if (famv <= 4) 
    varfunv <- famv
  else varfunv <- match(family$varfun, varfuns, -1)
  varfunv <- as.integer(varfunv)
  corstrv <- as.integer(match(corstr, corstrs, -1))
  if (linkv < 1) 
    stop("unknown link.")
  if (varfunv < 1) 
    stop("unknown varfun.")
  if (corstrv < 1) 
    stop("unknown corstr.")
  naivvar <- matrix(rep(0, p * p), nrow = p)
  robvar <- matrix(rep(0, p * p), nrow = p)
  phi <- as.double(scale.value)
  scale.fix <- as.integer(scale.fix)
  errstate <- as.integer(1)
  tol <- as.double(tol)
  Mv <- as.integer(Mv)
  maxcl <- as.integer(0)
  if (!(is.double(x))) 
    x <- as.double(x)
  if (!(is.double(y))) 
    y <- as.double(y)
  if (!(is.double(id))) 
    id <- as.double(id)
  if (!(is.double(N))) 
    N <- as.double(N)
  modvec <- as.integer(c(linkv, varfunv, corstrv))
  if (v4.4compat) 
    compatflag <- 1
  else compatflag <- 0
  z <- .C(Cgee, x, y, id, N, offset, nobs, p, modvec, Mv, estb = beta, 
          nv = naivvar, rv = robvar, sc = phi, wcor = R, tol, mc = maxcl, 
          iter = maxiter, silent, err = errstate, scale.fix, as.integer(compatflag))
  if (z$err != 0) 
    warning("Cgee had an error (code= ", z$err, ").  Results suspect.")
  if (min(eigen(z$wcor)$values) < 0) {
    warning("Working correlation estimate not positive definite")
    z$err <- z$err + 1000
  }
  fit <- list()
  attr(fit, "class") <- c("gee", "glm")
  fit$title <- "GEE:  GENERALIZED LINEAR MODELS FOR DEPENDENT DATA"
  fit$version <- "gee S-function, version 4.13 modified 98/01/27 (1998)"
  links <- c("Identity", "Logarithm", "Logit", "Reciprocal", 
             "Probit", "Cloglog")
  varfuns <- c("Gaussian", "Poisson", "Binomial", "Gamma")
  corstrs <- c("Independent", "Fixed", "Stationary M-dependent", 
               "Non-Stationary M-dependent", "Exchangeable", "AR-M", 
               "Unstructured")
  fit$model <- list()
  fit$model$link <- links[linkv]
  fit$model$varfun <- varfuns[varfunv]
  fit$model$corstr <- corstrs[corstrv]
  if (!is.na(match(c(corstrv), c(3, 4, 6)))) 
    fit$model$M <- Mv
  fit$call <- call
  fit$terms <- Terms
  fit$formula <- as.vector(attr(Terms, "formula"))
  fit$contrasts <- attr(x, "contrasts")
  fit$nobs <- nobs
  fit$iterations <- z$iter
  fit$coefficients <- as.vector(z$estb)
  fit$nas <- is.na(fit$coefficients)
  names(fit$coefficients) <- xnames
  eta <- as.vector(x %*% fit$coefficients)
  fit$linear.predictors <- eta
  mu <- as.vector(family$linkinv(eta))
  fit$fitted.values <- mu
  fit$residuals <- y - mu
  fit$family <- family
  fit$y <- as.vector(y)
  fit$id <- as.vector(id)
  fit$max.id <- z$mc
  z$wcor <- matrix(z$wcor, ncol = maxclsz)
  fit$working.correlation <- z$wcor
  fit$scale <- z$sc
  fit$robust.variance <- z$rv
  fit$naive.variance <- z$nv
  fit$xnames <- xnames
  fit$error <- z$err
  dimnames(fit$robust.variance) <- list(xnames, xnames)
  dimnames(fit$naive.variance) <- list(xnames, xnames)
  fit
}