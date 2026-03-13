\name{coxlogLik}
\alias{coxlogLik}
\title{Compute the Log-Likelihood for a Cox Proportional Hazards Model}
\description{
Computes either the partial log-likelihood or the full log-likelihood for 
a Cox proportional hazards (PH) model.

The function evaluates the log-likelihood given a design matrix, survival outcome, 
regression coefficients, and optional offset and baseline hazard function. 
The data may optionally be pre-sorted by time to improve performance when the 
function is called repeatedly with the same outcome data.
}

\usage{
coxlogLik(X, y, beta, offset = NULL, H0 = NULL, partial = TRUE, sorted = FALSE)
}

\arguments{
  \item{X}{
  A numeric design matrix of covariates with dimension \eqn{n \times p}.
  }

  \item{y}{
  Two-column numeric matrix containing survival data.
  The first column gives the observed time and the second column the
  event indicator (1 = event, 0 = censored).
  }

  \item{beta}{
  A numeric vector of regression coefficients of length \eqn{p}.
  }

  \item{offset}{
  Optional numeric vector of length \eqn{n} representing an offset added to the linear predictor. Default is \code{NULL}.
  }

  \item{H0}{
  Optional cumulative baseline hazard function. If provided, it should be a function that takes time values as input and returns the baseline hazard at those times. 
  If \code{NULL}, the baseline hazard is estimated using the Breslow-type estimator \eqn{h_t = d_t / S_0(t)} when computing the full likelihood.
  }

  \item{partial}{
  Logical value indicating whether to compute the Cox partial log-likelihood.  
  If \code{TRUE} (default), the partial likelihood is returned.  
  If \code{FALSE}, the full log-likelihood is computed.
  }

  \item{sorted}{
  Logical value indicating whether the observations are already sorted by decreasing time.  
  If \code{FALSE} (default), the function sorts the data internally.  
  If \code{TRUE}, sorting is skipped, which can improve performance when repeatedly evaluating the likelihood with the same outcome data.
  }
}

\details{
Let \eqn{\eta_i = x_i^\top \beta + \text{offset}_i} denote the linear predictor and
\eqn{\exp(\eta_i)} the relative risk.

For the partial likelihood, the function computes

\deqn{
\ell(\beta) = \sum_{i=1}^n d_i \log\left(\frac{\exp(\eta_i)}{S_0(t_i)}\right),
}

where \eqn{d_i} is the event indicator and
\eqn{S_0(t_i) = \sum_{j \in R(t_i)} \exp(\eta_j)} is the risk set sum.

For the full likelihood, the function evaluates

\deqn{
\ell(\beta) = \sum_{i=1}^n d_i (\log h(t_i) + \eta_i) - H(t_i)\exp(\eta_i),
}

where \eqn{h(t)} is the baseline hazard and \eqn{H(t)} is the cumulative baseline hazard.

If \code{H0} is not supplied, a Breslow-type estimate based on the observed events is used.
}

\value{
A numeric value representing the log-likelihood (partial or full) evaluated at \code{beta}.
}

\author{
Bingshu E. Chen
}

\examples{
set.seed(29)

n <- 100
p <- 3

X <- matrix(rnorm(n * p), n, p)
beta <- c(0.5, -0.3, 0.2)

time <- rexp(n)
event <- rbinom(n, 1, 0.7)

y <- cbind(time, event)

# Partial log-likelihood
coxlogLik(X, y, beta)

# Full log-likelihood with estimated baseline hazard
coxlogLik(X, y, beta, partial = FALSE)

# Example with offset
offset <- rnorm(n)
coxlogLik(X, y, beta, offset = offset)
}
