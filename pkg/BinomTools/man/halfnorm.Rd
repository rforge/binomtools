\name{halfnorm}
\alias{halfnorm}
\title{Half normal plot with simulated envelopes
}
\description{
  \code{halfnorm} produces a half normal plot of the residuals with
  simulated envelopes useful for model evaluation and detection of outliers
}
\usage{
halfnorm(object, resType = c("deviance", "aLik", "eLik", "pearson"), env
= T, nsim = 20, plot = T, identify = F, n = 2) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{glm} with a binomial family}
  \item{resType}{The type of residuals used in the plot}
  \item{env}{Logical for whether the envelopes are simulated}
  \item{nsim}{Number of simulations used for the envelope}
  \item{plot}{Logical for whether the points should be plotted. If
    \code{plot = F}} a list is returned} 
  \item{identify}{Logical for whether it should be possible to identify
    points interactively}
  \item{n}{How many points should be identified. Only meaningful if
    \code{identify = T}}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  
  %% Husk forklare de forskellige typer residualer og henvis evt. til
  %% exLikResid og resBin
}
\value{
  If \code{plot = T} a plot is produced. Otherwise a list of the
  residuals and their expected values are returned
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Merete Kjær Hansen
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, resType = c("deviance", "aLik", "eLik", "pearson"), 
    env = T, nsim = 20, plot = T, identify = F, n = 2) 
{
    resType <- match.arg(resType)
    List <- halfnormList(object, resType, env, nsim)
    if (plot) {
        halfplot(List, env)
        if (identify) 
            idHalfplot(List, n, env)
    }
    else List
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
