\name{profile}
\alias{profile.glm}
\alias{plot.profile.glm}
\title{Profile likelihoods for parameters in binomial regression models}
\description{
  Generate and plot the profile likelihoods for each parameter in a
  binomial regression model
}
\usage{
\method{profile}{glm}(fitted, which.par, alpha = 0.005, max.steps = 50,
         nsteps = 8, step.warn = 5, trace = F, ...)  

\method{plot}{profile.glm}(x, which.par, likelihood = TRUE,
         log = FALSE, relative = TRUE, approx = TRUE, conf.int = TRUE,
         level = 0.95, n = 100, fig = TRUE, ylim = NULL, ...) 
}
\arguments{
  \item{fitted}{An object of class \code{glm} with a binomial family}
  \item{x}{An object of class \code{profile.glm}}
  \item{which.par}{A numeric or character vector with the parameters to
  be profiled. If missing all parameters are profiled}
  \item{alpha}{The likelihood is profiled in approximately the
  100*(1-alpha)\% confidence region} 
  \item{likelihood}{Logical for whether the profile likelihood or
  likelihood root should be plotted} 
  \item{log}{Logical for whether the profile likelihood should be
  plotted on log-scale. Ignored if \code{likelihood = FALSE}}
  \item{relative}{Logical for whether the profile likelihood or
  log-likelihood should be plotted on a relative or absolute
  scale. Ignored if \code{likelihood = FALSE}} 
  \item{approx}{Logical for whether a quadratic approximation should be
  included in the plot}
  \item{conf.int}{Logical for whether a confidence interval should be
  included in the plot}
  \item{level}{A scalar or numerical vector indicating the confindence
  level(s) to be included in the plot. Ignored if \code{conf.int =
  FALSE}}
  \item{n}{How many points to employ in the spline interpolation
  of the profile likelihood}
  \item{fig}{Logical for whether the profile likelihood should be
  plotted. If \code{fig = FALSE} the list of points from the spline
  interpolation is returned} 
  \item{ylim}{The limits of the y-axis in the plot} 
  \item{trace}{Logical for whether progress should be printed to the
  screen during the profiling process}
  \item{nsteps}{Number of profiling steps to take in each direction for
  each parameter. The number is approximate since the step size is
  determined according to a quadratic approximation to the profile
  log-likelihood, hence, the deviation of the value of \code{nsteps} to
  the actual number of steps performed is influenced by the degree of 
  irregularity of the profile likelihood}
  \item{max.steps}{The maximum number of profiling steps in each
  direction for each parameter. A warning is issued if the number of
  \code{max.steps} is reached}
  \item{step.warn}{A warning is issued if the the actual number of steps
  in either direction does not exceed the number of \code{step.warn}}
  \item{...}{Additional arguments passed to other methods}
}
\details{
  \code{lroot} returned by \code{profile} is the signed square-root
  of the usual profile likelihood
  \deqn{sgn(\theta-\hat{\theta})\sqrt{2(l(\hat{\theta})-l(\theta))}}
  where \eqn{\theta} is the
  parameter being profiled and \eqn{\hat{\theta}} is the maximum likelihood
  estimate of \eqn{\theta}. The appertaining \code{par.vals} is a vector
  of \eqn{\theta} values in an appropriate range around \eqn{\hat{\theta}}.
  
  The logical argument \code{likelihood} in \code{plot} controls if the
  profile likelihood or the likelihood root should be plotted. 
}
\value{
  For \code{profile}: a list of class profile.glm with a range of
  parameter values and lroot statistics for each parameter in
  \code{which.par} 

  For \code{plot}: if \code{fig = FALSE} a list with
  plotting points and confidence interval(s) for each parameter in
  \code{which.par} is returned. If \code{fig = TRUE} the list is
  returned invisibly.  
}
\references{
  Pawitan, Y. (2001). \emph{In All Likelihood: Statistical Modelling and
  Inference Using Likelihood}. Oxford University Press.
}
\note{
  The implementation of these functions are largely inspired by
  \code{\link[MASS]{profile.glm}} from the \code{MASS} package and
  \code{\link[ordinal]{profile.clm}} from the \code{ordinal}
  package. This work is a direct extension of \code{profile} from
  \code{MASS} with an extended set of warnings. The main
  difference, though, is in the plotting functionality, which enables
  plot of the usual profile likelihood and log-likelihood and the
  optional inclusion of confidence interval(s). 
}
\author{Merete K Hansen}
\seealso{
  \code{\link[stats]{glm}},
  \code{\link[MASS]{profile.glm}},
  \code{\link[MASS]{plot.profile}}
}
\examples{
data(serum)
serum.glm <- glm(cbind(y, n-y) ~ dose, family=binomial, data=serum)
pr <- profile(serum.glm)
plot(pr)
}
\keyword{regression}
\keyword{models}
