
profBin <- function(object, which.par, alpha = 0.001, max.steps = 50,
               nsteps = 8, step.warn = 5, trace = F, ...)
{
  ## Match and test input arguments
  call <- match.call()
  name.object <- call$object # Name of original fitted object
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
  call.object <- object$call # Original call
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
  prof.list <- vector('list', length=length(which.par))
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
  val <- structure(prof.list, summary.originalfit = summary(object),
  name.originalfit = name.object, call.originalfit = call.object) 
  class(val) <- c('profBin.glm')
  return(val)
}
pr1 <- profBin(Serum.glm)
pr1



#######################
## Plot of profile likelihood
#######################
# plot.profBin.glm - senere navn?
aa <- function(x, which.par, statistic = c("likelihood", "lroot"),
                             log = TRUE, approx = TRUE,
                             conf.int = TRUE, level = 0.95,
                             fig = TRUE)
                             # flere senere: n = ???, ylim = ...?, ...) 
{  
  ## Match and test arguments
  xnames <- names(x)
  if(missing(which.par)) which.par <- seq_along(xnames)
#  if(...) warning('should be character or numeric vector')
  if(is.character(which.par)) which.par <- match(which.par, xnames, 0)
#  if(!all...match...) warning('nogle matcher ikke')
#  if(which.par == 0) stop('ingen parametre matcher')
  p <- length(which.par)
  statistic <- match.arg(statistic)
  ## Misch extraction from original glm object
  orig.fit <- attr(x, 'summary.originalfit')
  name.originalfit <- attr(x, 'name.originalfit')
  call.originalfit <- attr(x, 'call.originalfit') 
  b0 <- coef(orig.fit)[, 'Estimate']
  std.err <- coef(orig.fit)[, 'Std. Error'] # brug vcov i stedet for
  ## Result list
  res.list <- c(name.originalfit, call.originalfit, vector('list', length=p))
  names(res.list) <- c('name.originalfit', 'call.originalfit', xnames[which.par])
  ## If values are plotted a device is opened
  if(fig) {
    rows <- ceiling(sqrt(p))
    op <- par(mfrow=c(rows, rows), oma = c(0,0,2,0))
    on.exit(par(op))
  } 
  ## For each wp calculate requisite values and plot them  
  for(wp in which.par) {
    par.wp <- x[[wp]]$par.values
    lik.wp <- x[[wp]]$lroot
    spline.wp <- spline(par.wp, lik.wp) # husk tilføj n=n
    wp.name <- xnames[wp]
    if(approx) approx.wp <- (par.wp - b0[wp.name])/std.err[wp.name]
### byttes om på par.wp og b0?
    if(conf.int) cutoff.wp <- qnorm((1-level)/2, lower.tail = F) 
    if(statistic == 'lroot') {
      ylab <- "Profile likelihood root"
      if(conf.int) cutoff.wp <- c(-1, 1) * cutoff.wp
    } ## end 'if statistic is lroot'
    else if(statistic == 'likelihood') {
      ylab <- "Profile log-likelihood"
      spline.wp$y <- -spline.wp$y^2/2
      if(approx) approx.wp <- -approx.wp^2/2
      if(conf.int) cutoff.wp <- -cutoff.wp^2/2
      if(!log)
        ylab <- "Profile likelihood"
        spline.wp$y <- exp(spline.wp$y)
        if(approx) approx.wp <- exp(approx.wp)
        if(conf.int) cutoff.wp <- exp(cutoff.wp)
    } ## end 'if statistic is likelihood'
    else
      stop("Type of statistic not recognized. Supported arguments are 'likelihood' and 'lroot'")  
    res.list[[wp.name]]$spline.vals <- spline.wp
    if(approx) res.list[[wp.name]]$approx.vals <- list(x = par.wp,
                  y = approx.wp) 
    if(conf.int) res.list[[wp.name]]$cutoff <- cutoff.wp

    if(fig) {
      plot(spline.wp$x, spline.wp$y, type='l', xlab=wp.name,
      ylab=ylab, col='black')
      if(statistic=='lroot') points(b0[wp.name], 0, pch=18)
      if(conf.int) abline(h = cutoff.wp)
      if(approx) {
        lines(par.wp, approx.wp, lty=2, col='steelblue')
        legend('topleft', legend=c(statistic, 'approx'), lty = 1:2,
        col=c('black','steelblue'), bty = 'n')
      }
      title(paste("\nProfile likelihood of parameter
      estimates from object '", name.originalfit, "'", sep=""), outer=T)
    } ## end 'if fig = T draw plot for each wp'
  } ## end 'for wp in which.par calculate values needed for plotting'
#    attr(spline.list, '')
    if(fig) invisible(res.list)
    else return(res.list)
}

cc <- aa(pr1, statistic='lik', conf.int=T, approx=T)
cc




