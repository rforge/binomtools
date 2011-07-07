########################
### Exact likelihood residuals - leave one out
########################
exLikResid <- function(object) {

  if(is.null(object$model)) 
    stop('glm object has to be fitted with \'model=TRUE\'')

  start <- coef(object) # start values to be used in model fitting
  fit <- object$fitted
  devResidAbs <- abs(deviance(object))
  mf <- model.frame(object)
  N <- nrow(mf)

  if(!is.null(model.weights(mf))) mw <- model.weights(mf)
  else mw <- rep(1, nrow(object$model))

  resp <- model.response(mf)
  if(NCOL(resp) == 1L) {
    n <- rep(1, nrow(mf))
    resp <- cbind(resp, n-resp)
    y <- resp
  }
  else y <- resp[,1] / rowSums(resp)

  object$model[,1] <- resp*mw		### new response including the weights
  newCall <- update(object, subset=-i, start=start, evaluate=F)

  residLik <- rep(0, 5)
  for(i in 1:N) {
    objLeave1out <- eval(newCall)
    devResidLeave1out <- deviance(objLeave1out)
    residLik[i] <- sign(y[i]-fit[i]) * sqrt(devResidAbs-devResidLeave1out) 
  }

  residLik
}



########################
### Function to calculate standardised pearson, deviance and likelihood residuals
########################
# tjek om residualer er rigtige

resBin <- function(object, type=c('deviance', 'aLik','eLik','pearson')) { 

  type <- match.arg(type)
  h <- hatvalues(object)

  res <- switch(type, 
    deviance = rstandard(object), aLik = rstudent(object), eLik = exLikResid(object), 
    pearson = resid(object, type='pearson')/sqrt(1-h))

  return(res)
}

