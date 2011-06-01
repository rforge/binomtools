# Data and model
Serum <- read.table('Data/Anti-pneumococcus serum.dat', header=T)
Serum.glm <- glm(cbind(y, n-y) ~ dose, family=binomial, data=Serum)
Serumlog.glm <- glm(cbind(y, n-y) ~ log(dose), family=binomial, data=Serum)
Serum.gau.glm <- glm(y ~ dose, family=gaussian, data=Serum)


## Bedre navn til funktion???
profBin <- function(object, which.par, alpha = 0.05, max.steps = 50,
               nsteps = 8, step.warn = 5, trace = F, ...)
{
  ## Match and test input arguments
  if(!any(class(object) == 'glm'))
    stop("Object must be of class glm")
  if(family(object)$family != 'binomial')
    stop("GLM object must be fitted with the binomial family")
### Må der fittes med quasi-binomial?
  stopifnot(is.numeric(alpha) && length(alpha) == 1 &&
            alpha > 0 && alpha < 1)
  stopifnot(round(max.steps) > round(nsteps))
  stopifnot(round(nsteps) > round(step.warn))
  stopifnot(round(nsteps) > 0 && round(step.warn) >= 0)
  max.steps <- round(max.steps)
  nsteps <- round(nsteps)
  step.warn <- round(step.warn)
  trace <- as.logical(trace)[1]
  ## Misc extraction from object
  mf <- model.frame(object)
  X <- model.matrix(object)
  y <- model.response(mf)
  n <- length(y)
  if(!is.null(dim(y))) n <- n/2
  if(!is.null(model.weights(mf))) w <- model.weights(mf)
  else w <- rep(1, n)
  if(!is.null(model.offset(mf))) o <- model.offset(mf)
  else o <- rep(0, n)
### Tjekke om ovenstående går godt uanset dim(y), w og o
  fam <- family(object)
  p <- length(pnames <- names(b0 <- coef(object)))
  if(missing(which.par)) which.par <- 1:p
  if(!is.vector(which.par, mode='integer') && !is.vector(which.par,
               mode='character')) 
    stop("'which.par' must be a character or integer vector")
  if(is.character(which.par)) which.par <- match(which.par, pnames, 0)
  if(any(is.na(which.par)))
    stop("Invalid parameter argument(s) in 'which.par'")
  stopifnot(length(which.par) > 0)
  std.err <- coef(summary(object))[,'Std. Error']
  orig.dev <- deviance(object)
  disp <- summary(object)$dispersion
### skal kun bruges hvis quasi er tilladt 
  ## Profile limit
  lroot.max <- qnorm(1 - alpha/2)
  ## Step length
  delta <- lroot.max / nsteps
  ## Result list
  prof.list <- vector('list', length=length(pnames[which.par]))
  names(prof.list) <- pnames[which.par]
  ## for each which.par move down and then up (from MLE), fit model
  ## with 'wp' fixed at ascending/descending values (using offset) and
  ## store signed likelihood root values and parameter values
  for(wp in which.par) {
    X.wp <- X[, -wp, drop = FALSE] ## Design matrix without wp
    par.wp <- b.wp <- b0[wp] ## Reset wp coefficient to initial value
    lroot.wp <- 0 ## lroot at MLE
    wp.name <- pnames[wp] ## Name of wp parameter
    for(direction in c(-1, 1)) { ## Move down and then up
      if(trace) {
        message("\nParameter '", wp.name, "' ",
               c('down','up')[(direction + 1)/2 + 1])
        utils::flush.console()
      }
      step <- 0
      lroot <- 0
      start <- b0[-wp]
      ## Increment wp (in offset) and refit model until threshold is reached
      while(step < max.steps && abs(lroot) <= lroot.max) {
        step <- step + 1
        bi <- b.wp + direction * step * delta * std.err[wp]
        oi <- o + X[, wp] * bi
        attributes(oi) <- NULL
        fm <- glm.fit(x = X.wp, y = y, weights = w, start = start,
                      offset = oi, family = fam, control = object$control)
### Tage højde for om intercept er med??
        ## Likelihood ratio statistic
        lrs <- (fm$deviance - orig.dev)/disp 
### Likelihood ratio statistic # er det rigtigt og gælder det
### generelt? - skal dispersion være med?
        ## Likelihood root statistic
        lroot <- direction * sqrt(lrs) 
        lroot.wp <- c(lroot.wp, lroot)
        par.wp <- c(par.wp, bi)    
        ## Update start values
        start <- coef(fm) 
      } ## end 'while step <= max.steps and lroot <= lroot.max'
      ## test that lroot.max is reached (e.g. max.steps not reached)
      if(abs(lroot) < lroot.max)
        warning("profile may be unreliable for '", wp.name,
                "' because lroot.max was not reached for ",
                wb, c(" down", " up")[(direction + 1)/2 + 1])
      ## test that enough steps are taken
      if(step <= step.warn)
        warning("profile may be unreliable for '", wp.name,
                "' because only ", step, "\n  steps were taken ",
                c("down", "up")[(direction + 1)/2 + 1])
    } ## end 'for direction down or up'
    par.order <- order(par.wp)
    prof.list[[wp.name]] <-
      structure(data.frame(par.wp[par.order]), names = "par.values")
    prof.list[[wp.name]]$lroot <- lroot.wp[par.order]

### Fikse warning    
###    if(!all(diff(lroot.wp[par.order, wp.name]) > 0))
###      warning("likelihood is not monotonically decreasing from maximum,\n",
###              "  so profile may be unreliable for ", wb.name)
    
  } ## end 'for wp in which.par'
  val <- structure(prof.list, summary.originalfit = summary(object)) 
  class(val) <- c('profBin.glm')
  return(val)
}
profBin(Serum.glm, trace=T)


## To-do next: confint og plot


  
# To consider...
# Hvad med dispersion parameter i originalt fit?
# profile deviance?
# restriktioner på link-typen??


# Options plot
# - likelihood eller likelihood root statistic
# - absolut eller relativ scale
# - indtegne kvadratisk approximation


