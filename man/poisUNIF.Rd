\name{poisUNIF}
\Rdversion{1.1}
\alias{poisUNIF}

\title{Bayesian Prediction Limits for Poisson Distribution (Uniform Prior)}

\description{
The function provides the Bayesian prediction limits of a Poisson random variable derived based on a uniform prior.  The resulting prediction bounds quantify the uncertainty associated to the predicted future number of occurences in a time windows of size \eqn{t}{t}.
}

\usage{poisUNIF(xobs, n, s, t, alpha = 0.05)}

\arguments{
  \item{xobs}{a numeric value denoting the number of the observed occurrencies.}
  \item{n}{a numeric value representing the total number of the time windows \code{s} in the past (observed time windows).}
  \item{s}{a numeric value corresponding to the fixed size (or average size) of the observed time windows.}
  \item{t}{a numeric value indicating the size of the future time window.}
  \item{alpha}{a numeric value associated to the credible probability.  By default \code{alpha = 0.05}, thus an prediction interval at 95\% will be returned.}
%  \item{epsilon}{optional numeric value representing a level of tollerance to adjust the estimates when \code{xobs = 0}.  By default \code{epsilon = NULL}, thus a value is automatically assigned depending on the machine's precision (see \code{\link{.Machine}}).}
}

\details{
The resulting limits are equivalent to those provided when running the function \code{poisBayes()} with arguments \code{a = 1} and \code{b = Inf}. 
}

\value{
A list containing the following components:
%\item{maxr}{The maximum value used as an upper limit to bound the support of the Poisson random variable.}
\item{lower}{An integer value representing the lower bound of the prediction limit.}
\item{upper}{An integer value representing the upper bound of the prediction limit.}
}

\references{
Bejleri, V., & Nandram, B. (2018). Bayesian and frequentist prediction limits for the Poisson distribution. \emph{Communications in Statistics-Theory and Methods}, \bold{47}(17), 4254-4271.

Bejleri, V. (2005). \emph{Bayesian Prediction Intervals for the Poisson Model, Noninformative Priors}, Ph.D. Dissertation, American University, Washington, DC.
}

\author{Valbona Bejleri, Luca Sartore and Balgobin Nandram}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{poiss}}, \code{\link{poisJEFF}}, \code{\link{poisBayes}}}

\examples{
# Loading the package
library(plpoisson)
set.seed(2020L)

# Number of observed time windows
n <- 555L

# Simulating a dataset
data <- cbind.data.frame(
    occ_obs = rpois(n, rgamma(n, 5.5, .5)),
    win_siz = rgamma(n, 1.44, .777)
) 

## Bayesian prediction limits 
##  (with uniform prior)
poisUNIF(sum(data$occ_obs), # Past occurrencies 
    nrow(data), # Total past time windows
    mean(data$win_siz), # Window size
    444) # Size of future window
}

\keyword{ distribution }% use one of  RShowDoc("KEYWORDS")
\keyword{ poisson }% __ONLY ONE__ keyword per line
\keyword{ prediction }
\keyword{ intervals }
\keyword{ bayesian }
\keyword{ uniform prior }
