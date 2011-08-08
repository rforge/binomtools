\name{Rsq.glm}
\alias{Rsq}
\alias{Rsq.glm}
\alias{print.Rsq}
\alias{plot.Rsq}
\title{
  R-squared measures for binomial GLMs
}
\description{
  This function computes the R-squared measures for binomial GLMs
  proposed by Tjur (2010) "Coefficients of determination in
  logistic regression models - a new proposal: The coefficient of
  discrimination".
}
\usage{
\method{Rsq}{glm}(object, ...)

\method{print}{Rsq}(x, digits = getOption("digits"), ...)

\method{plot}{Rsq}(x, which=c("hist", "ecdf", "ROC"), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    a binomial \code{\link{glm}} object
  }
  \item{x}{
    an \code{Rsq} object
  }
  \item{which}{
    the desired plot: histograms, empirical cumulative distribution
    functions or ROC (receiver operating characteristic) curve
  }
  \item{digits}{
    the desired number of printed digits
  }
  \item{\dots}{
    currently not used
  }
}
\details{
  The plot method has the following options
  \describe{
  \item{\code{"hist"}}{
    Two histograms with ten bins of the fitted probabilities are plottet
    on top of each other; the upper one for \code{y = 0} and the lower
    one for \code{y = 1}.
  }
  \item{\code{"ecdf"}}{
    Two ecdf curves; one for \code{y = 0} and one for \code{y = 1}
  }
  \item{\code{"ROC"}}{
    The (empirical) ROC curve
  }
  }
}
\value{
  \code{Rsq.glm} returns an object of class \code{Rsq}. The \code{plot}
  and \code{print} methods returns the \code{Rsq} objects invisibly.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  Tjur, T. (2009) Coefficients of determination in
  logistic regression models - a new proposal: The coefficient of
  discrimination. \emph{The American Statistician}, \bold{63}(4), 366-372.
}
\author{
  Rune Haubo B Christensen
}
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
  A \code{\link[binomTools]{HLtest}} (Hosmer and Lemeshow test)
  method exists for \code{Rsq} objects.
}
\examples{
## Lifted from example(predict.glm):
ldose <- rep(0:5, 2)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
sex <- factor(rep(c("M", "F"), c(6, 6)))
SF <- cbind(numdead, numalive=20-numdead)
budworm.lg <- glm(SF ~ sex*ldose, family=binomial)
## summary(budworm.lg)

(Rsq.budworm <- Rsq(budworm.lg))

plot(Rsq.budworm, "hist") ## or simply 'plot(Rsq.budworm)'
plot(Rsq.budworm, "ecdf")
plot(Rsq.budworm, "ROC")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{tests}
