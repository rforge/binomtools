\name{Residuals}
\alias{Residuals}
\alias{Resid}
\title{Residuals from a binomial regression model}
\description{
  Function to extract residuals from a binomial regression model} 
\usage{
Residuals(object, type = c("approx.deletion", "exact.deletion",
          "standard.deviance", "standard.pearson", "deviance",
          "pearson", "working", "response", "partial"))
}
\arguments{
  \item{object}{An object of class \code{glm} with a binomial family}
  \item{type}{The type of residuals to be returned. Default is
    \code{approx.deletion} residuals} 
}
\details{
A considerable terminology inconsistency regarding residuals is found in
the litterature, especially concerning the adjectives 
\emph{standardized} and \emph{studentized}. Here, we use the term
\emph{standardized} about residuals divided by \eqn{\sqrt(1-h_i)} and
avoid the term \emph{studentized} in favour of \emph{deletion} to
avoid confusion. See Hardin and Hilbe (2007) p. 52 for a short
discussion of this topic. 

The objective of \code{Residuals} is to enhance transparency of
residuals of binomial regression models in \R and to uniformise the
terminology. With the exception of \code{exact.deletion} all residuals
are extracted with a call to \code{\link[stats]{rstudent}},
\code{\link[stats]{rstandard}} and \code{\link[stats]{residuals}} from
the \code{stats} package (see the description of the individual
residuals below).

\itemize{
  \item{\code{response}: response residuals
    \deqn{y_i - \hat{y}_i}
    The response residuals are also called raw residuals

    The residuals are extracted with a call to \code{\link[stats]{residuals}}.
  }
  \item{\code{pearson}: Pearson residuals
    \deqn{X_i =
      \frac{y_i - n_i \hat{p}_i}{\sqrt{n_i\hat{p}_i(1-\hat{p}_i)}}}{X_i = y_i-n_i
      \hat{p}_i/\sqrt{n_i \hat{p}_i (1-\hat{p}_i)}}

    The residuals are extracted with a call to \code{\link[stats]{residuals}}.
  }
  \item{\code{standard.pearson}: standardized Pearson residuals
    \deqn{r_{P,i} = \frac{X_i}{\sqrt{1-h_i}} =
      \frac{y_i+n_i\hat{p}_i}{\sqrt{n_i\hat{p}_i(1-\hat{p}_i)(1-h_i)}}}{r_{P,i}
    = X_i/\sqrt(1-h_i)}
    where \eqn{X_i} are the Pearson residuals and \eqn{h_i} are the
    hatvalues obtainable with \code{\link[stats]{hatvalues}}. 
    
    The standardized Pearson residuals have many names including
    studentized Pearson residuals, standardized residuals, studentized
    residuals, internally studentized residuals.
    
    The residuals are extracted with a call to \code{\link[stats]{rstandard}}.
  }
  \item{\code{deviance}: deviance residual
    
    The deviance residuals are the signed square roots of the individual
    observations to the overall deviance
    \deqn{d_i = sgn(y_i-\hat{y}_i)
      \sqrt{2 y_i \log\left( \frac{y_i}{\hat{y}_i}\right) + 2(n_i-y_i)
      \log\left( \frac{n_i-y_i}{n_i-\hat{y}_i}\right)}}{d_i = sgn(y_i-\hat{y}_i)) \sqrt{2
      y_i log(y_i/ \hat{y}_i) + 2(n_i-y_i)
      log((n_i-y_i)/(n_i-\hat{y}_i))}}  

    The residuals are extracted with a call to \code{\link[stats]{residuals}}.
  }
  \item{\code{standard.deviance}: standardized deviance residuals
    \deqn{r_{D,i} = \frac{d_i}{\sqrt{1-h_i}}}{r_{D,i} = d_i/\sqrt(1-h_i)}
    where \eqn{d_i} are the deviance residuals and \eqn{h_i} are the
    hatvalues that can be obtained with \code{\link[stats]{hatvalues}}.
    
    The standardized deviance residuals are also called studentized
    deviance residuals. 
    
    The residuals are extracted with a call to \code{\link[stats]{rstandard}}.
  }
  \item{\code{approx.deletion}: approximate deletion residuals
    \deqn{sgn(y_i-\hat{y}_i)\sqrt{h_i r^2_{P,i}+(1-h_i)r^2_{D,i}}}{sgn(y_i-\hat{y}_i)
    \sqrt{h_i r_{P,i}^2 + (1-h_i) r_{D,i}^2}} 
    where \eqn{r_{P,i}} are the standardized Pearson residuals,
    \eqn{r_{D,i}} are the standardized deviance residuals and \eqn{h_i}
    are the hatvalues that is obtained with \code{\link[stats]{hatvalues}}  
    The approximate deletion residuals are approximations to the exact
    deletion residuals (see below) as suggested by Williams (1987). 

    The approximate deletion residuals are called many different names in the
    litterature including likelihood residuals, studentized residuals,
    externally studentized residuals, deleted studentized residuals and
    jack-knife residuals.  

    The residuals are extracted with a call to \code{\link[stats]{rstudent}}. 
  }
  \item{\code{exact.deletion}: exact deletion residuals
    
    The \eqn{i}th deletion residual is calculated subtracting the
    deviances when fitting a linear logistic model to the full set of
    \eqn{n} observations and fitting the same model to a set of \eqn{n-1}
    observations excluding the \eqn{i}th observation, for \eqn{i =
      1,...,n}. This gives rise to \eqn{n+1} fitting processes and may be
    computationally heavy for large data sets.
  }
  \item{\code{working}: working residuals
    
    The difference between the working response and the linear predictor
    at convergence
    \deqn{r_{W,i} = (y_i -
    \hat{y}_i)\frac{\partial\hat{\eta}_i}{\partial\hat{\mu}_i}}{r_{W,i} = (y_i-\hat{y}_i)d\hat{eta}_i/d\hat{mu}_i}

    The residuals are extracted with a call to \code{\link[stats]{residuals}}.
  }
  \item{\code{partial}: partial residuals
    \deqn{r_{W,i} + x_{ij} \hat{\beta}_j}{r_{W,i} + x_{ij} \hat{\beta}_j}
    where \eqn{j = 1,...,p} and \eqn{p} is the number of
    predictors. \eqn{x_{ij}} is the \eqn{i}th observation of the \eqn{j}th
    predictor and \eqn{\hat{\beta}_j} is the \eqn{j}th fitted coefficient.

    The residuals are useful for making partial residuals plots. They
    are extracted with a call to \code{\link[stats]{residuals}}
  }
}
}
\value{A vector of residuals}
\references{ 
  Collett, D. (2003) \emph{Modelling binary data}. Second
  edition. Chapman & Hall/CRC. 

  Fox, J. (2002) \emph{An R and S-Plus Companion to Applied
  Regression}. Sage Publ. 

  Hardin, J.W., Hilbe, J.M. (2007). \emph{Generalized Linear Models and
  Extensions}. Second edition. Stata Press.  

  Williams, D. A. (1987) Generalized linear model diagnostics using the
  deviance and single case deletions. \emph{Applied Statistics}
  \bold{36}, 181-191. 
}
\author{Merete K Hansen}
\examples{
data(serum)
serum.glm <- glm(cbind(y, n-y) ~ log(dose), family = binomial, data = serum)
Residuals(serum.glm, type='standard.deviance')
}
\keyword{ models }
