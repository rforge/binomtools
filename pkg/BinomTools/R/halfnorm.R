### Halfnorm points and optional envelope returned as a list
halfnormList <- function(object, resType, env, nsim)
{

  res <- resBin(object, type=resType)

  ## Absolute values of residuals in ascending order plotted against Phi^(-1){(i+n-1/8)/(2n+1/2)}
  no.res <- length(object$resid)
  expected <- qnorm((1:no.res + no.res - 1/8)/(2*no.res + 1/2))
  sorted <- sort(abs(res))

  EnvList <- halfnormEnv(object, resType, nsim)

  if(env) List <- c(list(residuals=sorted, expected=expected), EnvList)
  else List <- list(residuals=sorted, expected=expected, resType=resType) 

  List
}


### Simulated envelope
halfnormEnv <- function(object, resType, nsim) 
{
  
  mf <- model.frame(object)
  mr <- model.response(mf)
  groupSize <- apply(mr, 1, sum)
  if(!is.null(model.weights(mf))) mw <- model.weights(mf)
  else mw <- rep(1, nrow(mf))
  call <- update(object, formula=YY ~ ., data = simData, weights = mw, evaluate = F)

  ## Simulation of data
  simData <- object$model
  simData$mw <- mw
  resSim <- data.frame(rep(0, length(groupSize)))
  for(i in 1:nsim) {
    suc <- rbinom(dim(object$data)[1], groupSize, object$fitted)
    simData$YY <- cbind(suc, groupSize-suc)
    simObj <- eval(call)
    resSim[,i] <- sort(abs(resBin(simObj, type=resType)))
  }

  ## Obtaining min, mean and max values for the envelope
  minValues <- apply(resSim, 1, min)
  meanValues <- apply(resSim, 1, mean)
  maxValues <- apply(resSim, 1, max)

  envList <- c(list(minValues=minValues, meanValues=meanValues, maxValues=maxValues, resType=resType))
  envList
}


### Plot of points and optional envelope
halfplot <- function(list, env, type=NULL)
{

  xy <- xy.coords(list$expected, list$residuals); x <- xy$x; y <- xy$y; 
  minEnv <- list$minValues
  meanEnv <- list$meanValues
  maxEnv <- list$maxValues
  resType <- list$resType

  resType <- switch(resType, deviance='deviance', aLik='approximated likelihood', eLik='exact likelihood', pearson='Pearson')

  if(env) ylim <- c(0, max(maxEnv))
  else ylim <- c(0, max(y))

  op <- par(mar=c(5, 5.5, 4, 2) + 0.1)
  plot(x, y, pch=20, cex.axis=0.8, las=1, tcl=-0.3, mgp=c(3, 0.6, 0),
    ylim=ylim,
    xlab='Expected value of half-normal order statistic', 
    ylab=paste('Absolute value of standardized \n', resType, 'residuals'), 
    main='Half-normal plot of residuals', type=type)
  if(env) {
    lines(x, minEnv, lty=3)
    lines(x, meanEnv)
    lines(x, maxEnv, lty=3)
  }
  par(op)
}


### Function to identify potential outliers
idHalfplot <- function(list, n, env, pch=NULL, tolerance=0.25, ...)
{
  xy <- xy.coords(list$expected, list$residuals); x <- xy$x; y <- xy$y
  if(is.null(pch)) pch <- names(list$residuals)
  ind <- rep(FALSE, length(x)); obsList <- integer(0)

  while(sum(ind) < n) {
    obs <- identify(x[!ind], y[!ind], n=1, plot=FALSE, tolerance=tolerance, ...)
    obs <- which(!ind)[obs]
    obsList <- c(obsList, obs)
    ind[obs] <- TRUE
    halfplot(list, env, type='n')
    points(x[which(!ind)], y[which(!ind)], pch=20)
    text(x[obsList], y[obsList], labels=pch[obsList], cex=0.75)
  }
}


### Final halfnorm function
halfnorm <- function(object, resType=c('deviance', 'aLik', 'eLik', 'pearson'), env=T, nsim = 20, plot=T, identify=F, n=2)  {

  resType <- match.arg(resType)
  List <- halfnormList(object, resType, env, nsim)

  if(plot) {
    halfplot(List, env)
    if(identify)
      idHalfplot(List, n, env)}
  else 
    List

}







