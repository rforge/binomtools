\name{group}
\alias{group}
\alias{group.glm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Group observations in a binomial glm
}
\description{
  This function groups the observations in a binomial glm based on the
  covariate structure. This can make it possible to assess
  goodness-of-fit in some models fitted to binary observations.
}
\usage{
\method{group}{glm}(object, eval = TRUE, ind = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    a binomial glm object
  }
  \item{eval}{
    should the new glm-model be evaluated?
  }
  \item{ind}{
    an indicator for which rows to keep(?). If this is no specified the
    grouping structure is based on the covariate structure in the
    model.
  }
  \item{\dots}{
    currently not used
  }
}
\details{
  The residual deviance and residual Pearson deviance are not meaningful
  measures of goodness-of-fit if the expected frequencies under the
  model are small (say less than five).

  if \code{eval = TRUE} it is tested whether the estimated coefficients
  are identical up to three significant digits and a warning is issued
  if this is not the case. This should be the case in well-behaved
  situations but may not happen in cases of complete separation.
}
\value{
  A list with components
  \item{newCall}{
    the new call
  }
  \item{newData}{
    a data frame with the aggregated data set
  }
  \item{oldData}{
    a data frame with the original data set
  }
  \item{oldN}{
    the number of rows (cases / observations) in the original data set
  }
  \item{newN}{
    the number of rows (cases / observations) in the aggregated data set
  }
  \item{oldObject}{
    the original fitted model
  }
  \item{newObject}{
    if \code{eval = TRUE} the new fitted model object, otherwise empty
  }
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  Collett (2001)
  Venables and Ripley (2002)
%% ~put references to the literature/web site here ~
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
## budworm.lg <- glm(SF ~ sex*ldose, family=binomial)
## summary(budworm.lg)
dat <- data.frame(SF=SF, sex, ldose)
dat[10, 1:2] <- rep(5, 2)
dat[13, ] <- dat[10, ]
rm(SF, sex, ldose)
SF <- as.matrix(dat[,1:2])
dat <- dat[,-(1:2)]
dat <- as.data.frame(cbind(SF, dat))

summary(m0 <- glm(SF ~ sex*ldose, binomial, dat))

## Various types of grouping:
(ind <- c(1:12, 10))
g <- group(m0, ind=ind, eval=TRUE)
g <- group(m0, eval=FALSE)
g <- group(m0, eval=TRUE)

## The correct GOF-test from the residual deviance is given by:
g$newObject

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{tests}
