\name{ibs}
\alias{brierScore}
\alias{ibs}
\alias{ibs.coxph}
\alias{ibs.default}
\alias{ibs.lple}
\alias{ibs.Surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Brier Score and Integrated Brier Score (IBS)
}
\description{
Calculate the Brier score and the integration of the Brier score (IBS) using the Inverse Probability of Censoring Weighting (IPCW) method.
}
\usage{
brierScore(object, St, tau)
\method{ibs}{default}(object, ...)
\method{ibs}{coxph}(object, newdata = NULL, newy = NULL, ...)
\method{ibs}{lple}(object, newdata = NULL, newy = NULL, ...)
\method{ibs}{Surv}(object, survProb, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{for ibs.Surv and ibs.default, it is a survival object created by Surv(time, event). 
        For others, it is a model object returned by coxph, lple.}
  \item{newdata}{optional new data at which the IBS is calculated. If absent, IBS is for the dataframe used in the original model fit.}
  \item{newy}{optional new survival object data. Default is NULL.}
  \item{St}{the predicted survival function at time tau to calcuate the Brier score.}
  \item{survProb}{the predicted survival function matrix. Row denotes each subject and column denotes each time points. 
        survProb[i,j] denotes the predicted survival probability of the ith subject at the time t[j].}
  \item{tau}{the time point at which the Brier score is calculated.}
  \item{\dots}{additional arguments to be passed to the functions such as ibs.coxph, ibs.lple, ibs.Surv etc.}
}
\value{A value of the Brier score or integration of the Brier score is returned.}
\details{
  The Brier score is the mean square difference between the true survival status and the predicted survival function.
  The Brier score is defined as,

    bs(tau) = 1/n*I(T_i>tau, delta_i = 1) S(t)^2/G(T_i) + (1-S(tau))^2/G(tau),

  where G = IPCW(Surv(time, event)), and \code{\link{IPCW}} is called to fit a KM model for the censoring time.

  The IBS is an integrated Brier Score over time. That is an integrated weighted squared distance 
  between the estimated survival function and the empirical survival function int_0 ^ 2 (I(T > t) - S(t))^2dt. 
  The inverse probability censoring weighting(\code{\link{IPCW}}) is used to adjust for censoring.
}
\author{Bingshu E. Chen}
\seealso{
 The \code{\link{IPCW}} method is used calculate the Brier score and the integrated Brier score. 
 A Cox proportional hazards (PH) model (\code{\link[survival]{coxph}}) shall be fitted to calculate Brier and IBS for the Cox PH model.

 The Brier score for the Cox model can also be calculated by \code{\link[survival]{brier}}.
}
\references{
  1. Brier, G. W. (1950). Verification of forecasts expressed in terms of probability. Monthly Weather Review, 78.

  2. Graf, Erika, Schmoor, Claudia, Sauerbrei, & Willi, et al. (1999). Assessment and comparison of prognostic classification schemes for survival data. Statistics in Medicine, 18, 2529-2545.
}
\examples{
  set.seed(29)
  n     = 25
  time  = rexp(n, 1)
  event = rbinom(n, 1, 0.75)

  ### calculate the Brier score at time tau
  tau   = 0.5
  St    = pexp(rep(tau, n), 1, lower.tail = FALSE)
  bs    = brierScore(Surv(time, event), St, tau)

  ### calculate the integrated Brier score
  #fit   = coxph(Surv(time, event)~1)
  #IBS  = ibs(fit)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IPCW }
\keyword{ Brier Score }
