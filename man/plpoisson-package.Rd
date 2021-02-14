\name{plpoisson-package}
\Rdversion{1.1}
\alias{plpoisson-package}
\alias{plpoisson}
\docType{package}

\title{\packageTitle{plpoisson}}

\description{
\if{html}{\figure{logo.svg}{options: align='right'}}
Prediction limits for Poisson distribution are useful when quantifying the uncertainty associated with predicting the occurrences of real life phenomena. The \pkg{plpoisson} package provides a set of functions to compute prediction limits of the inferred Poisson distribution under both, frequentist and Bayesian frameworks.

For frequentist prediction a common approach is to estimate the parameter based on the observed data firstly, then, to predict based on the estimated parameter. Different from the common approach of frequentist prediction, this approach does not require the estimation of the parameter. In a Bayesian setting, the uniform, Jeffreys and gamma distributions are used as priors when deriving the predictive posterior distribution. 
}

\details{
\tabular{ll}{
  Package: \tab plpoisson\cr
  Type: \tab Package\cr
  Version: \tab 0.2.0\cr
  Date: \tab 2021-02-14\cr
  License: \tab GPL-3\cr
}

For a complete list of exported functions, use \code{library(help = "plpoisson")}.
}

\author{
Valbona Bejleri, Luca Sartore and Balgobin Nandram

Maintainer: Luca Sartore \email{drwolf85@gmail.com}
}

\references{
Bejleri, V., & Nandram, B. (2018). Bayesian and frequentist prediction limits for the Poisson distribution. \emph{Communications in Statistics-Theory and Methods}, \bold{47}(17), 4254-4271.

Bejleri, V. (2005). \emph{Bayesian Prediction Intervals for the Poisson Model, Noninformative Priors}, Ph.D. Dissertation, American University, Washington, DC.
}

\keyword{ package }
\keyword{ distribution }
\keyword{ poisson }
\keyword{ prediction }
\keyword{ interval }

\examples{
## Loading the package
library(plpoisson)

## Setting quantities of interest
xobs <- rpois(1, 50)    # Number of the observed occurrencies  
n <- 1                  # Total number of the time windows of
                        #   of size 's' observed in the past
s <- rgamma(1, 4, .567) # Fixed size of observed time windows
t <- rgamma(1, 3, .33)  # Future time window
a <- 5                  # Shape hyperparameter of a gamma prior
b <- 1.558              # Rate hyperparameter of a gamma prior

## Frequentist prediction limits
poiss(xobs, n, s, t)

## Bayesian prediction limits (with uniform prior)
poisUNIF(xobs, n, s, t)

## Bayesian prediction limits (with Jeffreys prior)
poisJEFF(xobs, n, s, t)

## Bayesian prediction limits (with gamma prior)
poisBayes(xobs, n, s, t, a, b)
}
