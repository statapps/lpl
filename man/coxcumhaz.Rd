\name{coxcumhaz}
\alias{coxcumhaz}
\title{Baseline Hazard, Cumulative Hazard and Survival Functions for a Cox Model}
\description{
Estimate the baseline hazard, cumulative hazard and survival functions for a
Cox proportional hazards model using the Breslow estimator. 
This function can replace \code{basehaz}, which is hard to use due to the 
confusion of centered = TRUE / FALSE concerns.
}

\usage{
coxcumhaz(y, linear.predictors = NULL, sorted = FALSE)
}

\arguments{
  \item{y}{
  A two-column numeric matrix containing the survival outcome.
  The first column gives the observed time and the second column the
  event indicator (1 = event, 0 = censored).
  }

  \item{linear.predictors}{
  Optional numeric vector of linear predictors \eqn{\eta_i = x_i^\top\beta}.
  If \code{NULL}, all relative risks are set to 1, corresponding to the
  baseline hazard estimate with no covariate effects.
  }

  \item{sorted}{
  Logical; if \code{TRUE}, the observations are assumed to be already
  sorted by increasing time. Use this to speed up the function when it is called 
  repeatedly for the same y. Otherwise the function sorts the data
  internally.
  }
}

\details{
Let \eqn{t_i} denote the observed time and \eqn{\delta_i} the event
indicator. Define the relative risk

\deqn{
r_i = \exp(\eta_i),
}

where \eqn{\eta_i} is the supplied linear predictor. If no linear predictors
are provided, \eqn{r_i = 1}.

Let \eqn{R(t)} denote the risk set just prior to time \eqn{t} and define

\deqn{
S_0(t) = \sum_{j \in R(t)} r_j.
}

The Breslow estimator of the baseline hazard increment at time \eqn{t_i} is

\deqn{
\hat{h}(t_i) = \frac{\delta_i}{S_0(t_i)}.
}

The cumulative baseline hazard estimate is

\deqn{
\hat{H}(t) = \sum_{t_i \le t} \hat{h}(t_i),
}

and the corresponding survival function is

\deqn{
\hat{S}(t) = \exp[-\hat{H}(t)].
}

The returned objects are interpolation functions created with
\code{approxfun}, allowing evaluation of \eqn{\hat{H}(t)} and
\eqn{\hat{S}(t)} at arbitrary time points.
}

\value{
A list containing:
\item{ht: }{
A function giving the estimated baseline hazard \eqn{\hat{H}(t)}.
}
\item{Ht: }{
A function giving the estimated baseline cumulative hazard \eqn{\hat{H}(t)}.
}
\item{St: }{
A function giving the estimated baseline survival function \eqn{\hat{S}(t)}.
}
}

\author{
Bingshu E. Chen
}

\examples{
set.seed(29)

n <- 20
time <- rexp(n)
event <- rbinom(n, 1, 0.7)
y <- cbind(time, event)

# Estimate baseline cumulative hazard and survival
fit <- coxcumhaz(y)

fit$ht(1)   # hazard at time 1
fit$Ht(1)   # cumulative hazard at time 1
fit$St(1)   # survival probability at time 1

# Example with linear predictors
lp  = rnorm(n)
chz = coxcumhaz(y, linear.predictors = lp)
Ht  = chz$Ht
}
