\name{empLogit}
\alias{empLogit}
\title{
  Calculates the empircal logit transform
}
\description{
  The empirical logit transform allows for a tolerence such that
  infinity is not returned when the argument is zero or one.
}
\usage{
empLogit(x, eps = 1e-3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    numerical vector for which the empirical logit transform is desired
  }
  \item{eps}{
    numerical scalar; a tolerence to prevent infinite values
  }
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{
  the empirical logit transform of \code{x}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
\author{
  Rune Haubo B Christensen
%%  ~~who you are~~
}
%% \note{
%% %%  ~~further notes~~
%% }
%%
%% %% ~Make other sections like Warning with \section{Warning }{....} ~
%%
%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

