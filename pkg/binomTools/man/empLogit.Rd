\name{empLogit}
\alias{empLogit}
\title{
  Calculates the empirical logit transform
}
\description{
  The empirical logit transform allows for a tolerence such that
  infinity is not returned when the argument is zero or one.
}
\usage{
empLogit(x, eps = 1e-3)
}
\arguments{
  \item{x}{
    numerical vector for which the empirical logit transform is desired
  }
  \item{eps}{
    numerical scalar; a tolerence to prevent infinite values
  }
}
\value{
  the empirical logit transform of \code{x}
}
\author{
  Rune Haubo B Christensen
}
\examples{

## The function is currently defined as
## function (x, eps = 1e-3) log((eps + x)/(1 - x + eps))

## Lifted from example(predict.glm):
ldose <- rep(0:5, 2)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
sex <- factor(rep(c("M", "F"), c(6, 6)))
SF <- cbind(numdead, numalive=20-numdead)
## budworm.lg <- glm(SF ~ sex*ldose, family=binomial)
## summary(budworm.lg)

empLogit(numdead/20)

## Possible usage:
## Explorative interaction plot:
interaction.plot(ldose, sex, empLogit(numdead/20))

}
\keyword{utilities}

