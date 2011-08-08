\name{halfnorm}
\alias{halfnorm}
\title{Half normal plot with simulated envelopes
}
\description{
  \code{halfnorm} produces a half normal plot of the residuals with
  simulated envelopes useful for model evaluation and detection of
  outliers 
}
\usage{
halfnorm(object, resType = c("approx.deletion", "exact.deletion",
         "standard.deviance", "standard.pearson", "deviance",
         "pearson", "working", "response", "partial"), env = T,
         nsim = 20, plot = T, identify = F, n = 2)   
}
\arguments{
  \item{object}{An object of class \code{glm} with a binomial family}
  \item{resType}{The type of residual used in the plot}
  \item{env}{Logical for whether envelopes are simulated}
  \item{nsim}{Number of simulations used for the envelopes}
  \item{plot}{Logical for whether the points should be plotted. If
    \code{plot = F} a list is returned} 
  \item{identify}{Logical for whether it should be possible to identify
    points interactively. Ignored if \code{plot = F}}
  \item{n}{How many points should be identified. Ignored if
    \code{identify = F}}
}
\details{
  Absolute values of the residuals are used in a half normal plot that
  otherwise corresponds to a regular normal probability plot.

  Residuals from a binomial glm are not necessarily uncorrelated and
  normally distributed and may accordingly deviate from a straight line
  even if the fitted model is true. If the fitted model is true the
  optional simulated envelopes are likely to contain the absolute
  residuals. 

  The different types of residuals are described in
  \code{\link[binomTools]{Residuals}} 
}
\value{
  If \code{plot = T} a plot is produced. Otherwise a list of the
  residuals and their expected values are returned
}
\references{
  Atkinson, A. C. (1981) Two graphical displays for outlying and
  influential observations in regression. \emph{Biometrika}, \bold{68},
  13-20.

  Collett, D. (2003) \emph{Modelling binary data}. Second
  edition. Chapman & Hall/CRC. 
}
\author{Merete K Hansen}
\seealso{
  \code{\link[binomTools]{Residuals}}, \code{\link[graphics]{identify}}
}
\examples{
## Halfnormal plot with simulated envelopes
data(beetles)
beetles.glm <- glm(cbind(y, n-y) ~ conc, family=binomial, data=beetles)
halfnorm(beetles.glm, resType='pearson')

\dontrun{
## Halfnormal plot with simulated envelopes
## Two points are interactively identified when they are selected with the mouse
halfnorm(beetles.glm, resType='deviance', identify = T, n = 2)
}
}
\keyword{ dplot }
\keyword{ iplot }
\keyword{ models }
