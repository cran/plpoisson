\name{poisBayes}
\Rdversion{1.1}
\alias{poisBayes}

\title{Bayesian Prediction Limits for Poisson Distribution (Gamma Prior)}

\description{
The function provides the Bayesian prediction limits of a Poisson random variable derived based on a gamma prior.  The resulting prediction bounds quantify the uncertainty associated with the predicted future number of occurences in a time window of size \eqn{t}{t}.
}

\usage{poisBayes(xobs, n, s, t, a, b, alpha = 0.05)}

\arguments{
  \item{xobs}{a numeric value denoting the number of the observed occurrencies.}
  \item{n}{a numeric value representing the total number of the time windows \code{s} in the past (observed time windows).}
  \item{s}{a numeric value corresponding to the fixed size (or average size) of the observed time windows.}
  \item{t}{a numeric value indicating the size of the future time window.}
  \item{a}{a poisitive real number denoting the shape hyperparameter of a gamma prior distribution.}
  \item{b}{a poisitive real number representing the rate hyperparameter of a gamma prior distribution.}
  \item{alpha}{a numeric value associated to the credible probability.  By default \code{alpha = 0.05}, thus an prediction interval at 95\% will be returned.}
%  \item{epsilon}{optional numeric value representing a level of tollerance to adjust the estimates when \code{xobs = 0}.  By default \code{epsilon = NULL}, thus a value is automatically assigned depending on the machine's precision (see \code{\link{.Machine}}).}
}

\details{
When the argument \code{b = Inf}, one can obtain prediction limits with uniform prior by setting the argument \code{a = 1}.  Similarly, one can get the limits with a Jeffreys prior by setting the argument \code{a = 0}.

Prediction bounds are obtained through the binary search algorithm.
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

Davis, C. H. (1969). The binary search algorithm. \emph{American Documentation (pre-1986)}, \bold{20}(2), 167.
}

\author{Valbona Bejleri, Luca Sartore and Balgobin Nandram}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{poiss}}, \code{\link{poisJEFF}}, \code{\link{poisUNIF}}}

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
##  (with gamma prior)
poisBayes(sum(data$occ_obs), # Past occurrencies 
    nrow(data), # Total past time windows
    mean(data$win_siz), # Window size
    333, # Size of future window
    2, 2.22) # Hyper-parameters for gamma prior
}

\keyword{ distribution }% use one of  RShowDoc("KEYWORDS")
\keyword{ poisson }% __ONLY ONE__ keyword per line
\keyword{ prediction }
\keyword{ intervals }
\keyword{ bayesian }
\keyword{ gamma prior }
