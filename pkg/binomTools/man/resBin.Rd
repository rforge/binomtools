\name{resBin}
\alias{resBin}
\title{Derivation of residuals from a logistic regression model}
\description{
  Function to derive residuals from a logistic regression model} 
\usage{
resBin(object, type = c("deviance", "aLik", "eLik", "pearson"))
}
\arguments{
  \item{object}{An object of class \code{glm} with a binomial family}
  \item{type}{The type of residuals that should be returned}
}
\details{
  The \code{deviance} residuals are the standardized deviance residuals
  from \code{\link[stats]{rstandard}}. The approximated likelihood
  residuals (\code{aLik}) are from \code{\link[stats]{rstudent}} and the
  exact likelihood residuals (\code{eLik}) from
  \code{\link[binomTools]{exLikResid}}. The \code{pearson} residuals are
  standardized pearson residuals and are calculated dividing the pearson
  residuals (obtained from \code{\link[stats]{residuals}}) by
  \eqn{sqrt(1-h_i)} (\eqn{h_i} can be extracted with
  \code{\link[stats]{hatvalues}}).
}
\value{A vector of residuals}
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
\examples{
data(serum)
serum.glm <- glm(cbind(y, n-y) ~ log(dose), family = binomial, data = serum)
resBin(serum.glm, type='deviance')
}
\keyword{ models }
