\name{exact.deletion}
\alias{exact.deletion}
\title{Exact deletion residuals}
\description{
  Function to derive exact values of deletion (leave-one-out)
  residuals for binomial regression models}
\usage{
exact.deletion(object)
}
\arguments{
  \item{object}{An object of class \code{glm} with a binomial family}
}
\details{
  The \eqn{i}th deletion residual is calculated subtracting the
  deviances when fitting a linear logistic model to the full set of
  \eqn{n} observations and fitting the same model to a set of \eqn{n-1}
  observations excluding the \eqn{i}th observation, for \eqn{i =
    1,...,n}. This gives rise to \eqn{n+1} fitting processes and may be
  computationally heavy for large data sets.
  
  Approximations to the deletion residuals, as described in Williams
  (1987), are provided by \code{\link[stats]{rstudent}}.

  Inconsistency regarding the terminology implies that the deletion
  residuals are called different names in the litterature, including
  likelihood residuals, studentized residuals, externally studentized
  residuals, deleted studentized residuals and jack-knife
  residuals. Conversely, some of these terms refer to different types of
  residuals   
}
\value{A vector with exact deletion residuals}
\references{
  Collett, D. (2003) \emph{Modelling binary data}. Second
  edition. Chapman & Hall/CRC. 

  Fox, J. (2002) \emph{An R and S-Plus Companion to Applied
  Regression}. Sage Publ. 
  
  Williams, D. A. (1987) Generalized linear model diagnostics using the
  deviance and single case deletions. \emph{Applied Statistics}
  \bold{36}, 181-191. 
}
\author{Merete K Hansen}
\seealso{
  \code{\link[binomTools]{Residuals}}, \code{\link[stats]{rstudent}}
}
\examples{
data(beetles)
beetles.glm <- glm(cbind(y, n-y) ~ log(conc), family=binomial, data=beetles)
exact.deletion(beetles.glm)
}
\keyword{ models }
