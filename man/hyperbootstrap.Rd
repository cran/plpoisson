\name{hyperbootstrap}
\Rdversion{1.1}
\alias{hyperbootstrap}

\title{Bootstrap Methods to Estimate Hyperparameters for a Gamma Prior}

\description{
The function provides three bootstrap implementations to estimate the hyperparameters of a gamma prior distribution.  The method of moments, maximum likelihood and chi-square approximation are implemented for studying the uncertainties associated with the choice of the hyperparameters \eqn{a}{a} (shape) and \eqn{b}{b} (rate).
}

\usage{hyperbootstrap(xvec, B = 1000L, 
               method = c("moments", "likelihood", "chisq"))}

\arguments{
  \item{xvec}{a numeric vector of data with the observed occurrencies (assumed to be Poisson distributed).}
  \item{B}{a numeric value representing the total number of bootstrap iterations.}
  \item{method}{a character string (or strings) with the name/s of the method/s chosen to obtain hyperparameter estiamtes.}
}

\details{The function performs a choosen number of iterations using either the method of momemnts (\code{method = "moments"}), the maximum likelihood (\code{method = "likelihood"}), or the chi-square approximation (\code{method = "chisq"}).
}

\value{
A list containing the following components:
\item{a}{A matrix of values for the shape hyperparameter of the gamma distribution. The results of each method are organized by column.}
\item{b}{A matrix of values for the rate hyperparameter of the gamma distribution. The results of each method are organized by column.}
}

\references{

Bejleri, V., Sartore, L. & Nandram, B. (2021). \emph{Asymptotic equivalence between frequentist and Bayesian prediction limits for the Poisson distribution}. Journal of the Korean Statistical Society \doi{10.1007/s42952-021-00157-x}

Bejleri, V. (2005). \emph{Bayesian Prediction Intervals for the hyperbootstrapon Model, Noninformative Priors}, Ph.D. Dissertation, American University, Washington, DC.

}

\author{Valbona Bejleri, Luca Sartore and Balgobin Nandram}

\seealso{\code{\link{poisBayes}}, \code{\link{poisJEFF}}, \code{\link{poisUNIF}}}

\examples{
# Loading the package
library(plpoisson)
set.seed(2021L)

# Number of observed time windows
n <- 555L

# Simulating a dataset
data <- cbind.data.frame(
    occ_obs = rpois(n, rgamma(n, 5.5, .5)),
    win_siz = rgamma(n, 1.44, .777)
) 

## Compute bootstrap estimates using all methods
hyperbootstrap(data$occ_obs, 10L) # only 10 iterations
}

\keyword{ distribution }% use one of  RShowDoc("KEYWORDS")
\keyword{ poisson }% __ONLY ONE__ keyword per line
\keyword{ prediction }
\keyword{ intervals }
\keyword{ frequentist }
