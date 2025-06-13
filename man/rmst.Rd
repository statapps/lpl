\name{rmst}
\alias{rmst}
\alias{rmstFit}
\alias{rmst.coxph}
\alias{rmst.default}
\alias{rmst.Surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The restricted mean survival time (RMST)
}
\description{
Calculate the restricted mean survival time (RMST) for Surv object, Cox proportional model and other survival objects.
}
\usage{
rmst(object, ...)
rmstFit(tau, h0 = NULL, H0 = function(x){x})
\method{rmst}{default}(object, ...)
\method{rmst}{coxph}(object, newdata = NULL, linear.predictors = NULL, tau=NULL, ...)
\method{rmst}{Surv}(object, tau = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{for rmst.Surv and rmst.default, it is a survival object created by Surv(time, event). 
        For others, it is a model object returned by coxph, lple.}
  \item{h0}{a hazard function to be used for restricted mean survival time calculation. If h0(t) is provided, then H0(t) will be ignored.}
  \item{H0}{a cumulative hazard function to be used for restricted mean survival time calculation. The default is H0(t) = t for t>0}
  \item{linear.predictors}{the linear predictor from the Cox PH model.}
  \item{newdata}{optional new data at which the RMST is calculated. If absent, RMST is for the dataframe used in the original model fit.}
  \item{tau}{the time point at which the restricted mean survival time is calculated.}
  \item{\dots}{additional arguments to be passed to the functions such as rmst.coxph, rmst.lple, rmst.Surv etc.}
}
\value{A value of the Brier score or integration of the Brier score is returned.}
\details{
  The restricted mean survival time (RMST) is the mean of the truncated survival time at some finite value tau.
  The RMST is defined as,

    RMST(tau) = E(min(T, tau)) = int_0 ^tau S(t)dt,

  where S(t) = P(T>t) is the survival function of the random variable T.

  rmstFit(tau, h0, H0) calculates the restricted mean survival time based on a hazard (or cumulative hazard) function. Only one function of either h0(t) or H0(t) is required. If h0(t) is provided, then H0(t) will be ignored. 
}
\author{Bingshu E. Chen}
\seealso{
  \code{\link[survival]{coxph}},
  \code{\link[survival]{Surv}}
}
%\references{}
\examples{
  set.seed(29)
  n     = 25
  time  = rexp(n, 1)
  event = rbinom(n, 1, 0.75)
  x     = rnorm(n)
  y     = Surv(time, event)

  ### calculate the restricted mean survival time at tau = 0.5
  rms   = rmst(y, tau = 0.5)

  ### calculate the integrated brier score
  #fit   = coxph(y~x)
  #RMST  = rmst(fit, tau = 2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Restricted Mean Survival Time }
\keyword{ Survival Analysis }
