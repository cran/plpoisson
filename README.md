# plpoisson: Prediction Limits for the Poisson Distribution


### Authors <img src="man/figures/logo.svg" align="right" alt="plpoisson logo" />
[Valbona Bejleri](mailto://valbona.bejleri@gmail.com)[<img alt="ORCID iD" src="https://cran.r-project.org/web/orcid.svg" width="16px" height="16px" style="margin-left:4px; margin-right:4px; vertical-align:middle">](https://orcid.org/0000-0001-9828-968X), [Luca Sartore](mailto://drwolf85@gmail.com)[<img alt="ORCID iD" src="https://cran.r-project.org/web/orcid.svg" width="16px" height="16px" style="width:16px; height:16px; margin-left:4px; margin-right:4px; vertical-align:middle">](https://orcid.org/0000-0002-0446-1328) and [Balgobin Nandram](mailto://balnan@wpi.edu)[<img alt="ORCID iD" src="https://cran.r-project.org/web/orcid.svg" width="16px" height="16px" style="width:16px; height:16px; margin-left:4px; margin-right:4px; vertical-align:middle">](https://orcid.org/0000-0002-3204-0301)

Maintainer: [Luca Sartore](mailto://drwolf85@gmail.com)

[![CRAN version](https://www.r-pkg.org/badges/version/plpoisson)](https://cran.r-project.org/package=plpoisson)
[![CRAN release](https://www.r-pkg.org/badges/ago/plpoisson)](https://cran.r-project.org/package=plpoisson)
[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-yellow.svg)](http://perso.crans.org/besson/LICENSE.html)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/plpoisson)](https://cran.r-project.org/package=plpoisson)
[![Total Downloads from CRAN RStudio mirror](https://cranlogs.r-pkg.org/badges/grand-total/plpoisson?color=orange)](https://CRAN.R-project.org/package=plpoisson)

## Features of the package
Prediction limits for Poisson distribution are useful when predicting the occurrences of some real life phenomena; in fact, these limits quantify the uncertainty associated with the predicted values.  The **plpoisson** package provides a set of functions to compute prediction limits of the inferred Poisson distribution under both frequentist and Bayesian frameworks.

For a complete list of exported functions, use `library(help = "plpoisson")` once the **plpoisson** package is installed (see the `inst/INSTALL.md` file for a detailed description of the setup process).

### Example
```R
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
```

## References

Bejleri, V. (2005). *Bayesian Prediction Intervals for the Poisson Model, Noninformative Priors*, Ph.D. Dissertation, American University, Washington, DC.

Bejleri, V., & Nandram, B. (2018). Bayesian and frequentist prediction limits for the Poisson distribution. *Communications in Statistics-Theory and Methods*, *47*(17), 4254-4271.

