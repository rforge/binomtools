\name{HLtest.Rsq}
\alias{HLtest}
\alias{X2GOFtest}
\alias{HLtest.Rsq}
\alias{print.HLtest.Rsq}
\alias{X2GOFtest.Rsq}
\alias{print.X2GOFtest.Rsq}
\title{
  Goodness-of-fit tests for GLMs for binary data
}
\description{
  Goodness-of-fit tests for GLMs for binary data including the
  Hosmer-Lemeshow decile test and X-squared test with normal
  approximation.
}
\usage{

\method{HLtest}{Rsq}(object, method = c("deciles", "fixed"),
         decile.type = 8, ...)

\method{print}{HLtest.Rsq}(x, digits = getOption("digits"), ...)

\method{X2GOFtest}{Rsq}(x, ...)

\method{print}{X2GOFtest.Rsq}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    An \code{Rsq} object
  }
  \item{x}{
    An \code{HLtest.Rsq} or an \code{X2GOFtest.Rsq} object
  }
  \item{method}{
    The type of Hosmer-Lemeshow test to be performed. The
    \code{"deciles"} method should be more accurate (Hosmer et al, 1997)
  }
  \item{decile.type}{
    The quantile computation method; see \code{\link{quantile}} for
    details
  }
  \item{digits}{
    the desired number of printed digits
  }
  \item{\dots}{
    currently not used
}
}
\details{
  These tests are known to have very low power. They are only
  appropriate when the fitted frequencies are very low and when the
  covariate pattern dictates strictly binary observations.
}
\value{
  For \code{HLtest.Rsq} an object of class \code{HLtest.Rsq} with
  components
  \item{expected}{
    the expected frequencies in the 2 x 10 entries
  }
  \item{observed}{
    the observed frequencies in the 2 x 10 entries
  }
  \item{resid}{
    Pearson residuals
  }
  \item{X2}{
    the Pearson X-squared statistic
  }
  \item{p.value}{
    the p-value for the goodness-of-fit test
  }
  \item{method}{
    the method used for the test
  }

  For \code{X2GOFtest} an object of class \code{X2GOFtest} with
  components
  \item{p.value}{
    the p-value for the goodness-of-fit test
  }
  \item{z.score}{
    the standardized z-score for the goodness-of-fit test
  }
  \item{A}{
    A
  }
  \item{RSS}{
    the residual sums of squares term
  }
  \item{X2}{
    the pearson chi-squared statistic
  }
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  Hosmer et al (????) ...

  Some Germans that I don't remember right now.
}
\author{
  Rune Haubo B Christensen
}
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
## Lifted from example(predict.glm):
ldose <- rep(0:5, 2)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
sex <- factor(rep(c("M", "F"), c(6, 6)))
SF <- cbind(numdead, numalive=20-numdead)
budworm.lg <- glm(SF ~ sex*ldose, family=binomial)
## summary(budworm.lg)

(Rsq.budworm <- Rsq(budworm.lg))

HLtest(Rsq.budworm)
HLtest(Rsq.budworm, method="fixed")
X2GOFtest(Rsq.budworm)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tests}

