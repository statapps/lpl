\name{survfit.lple}
\alias{survfit.lple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute a Survival Curve from a Local Linear Partial Likelihood Estimate. 
}
\description{
Computes the predicted survival function for a model fitted by (lple).
}
\usage{
## S3 method for class 'lple'
\method{survfit}{lple}(formula, se.fit=TRUE, conf.int=.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a fitted model from (lple) fit}
  \item{se.fit}{a logical value indicating whether standard errors shall be computed. Default is TRUE}
  \item{conf.int}{The level for a two-sided confidence interval on the survival curve. Default is 0.95}
  \item{...}{other arguments to the specific method}
}
\details{
  survfit.lple is called to compuate baseline survival function from the lple model \code{\link{lple}}.

  The default method, survfit has its own help page. Use methods("survfit") to get all the methods for the survfit generic.
}

\value{
%%  ~Describe the value returned
survfit.lple returns a list of predicted baseline survival function, cumulative hazard function and residuals. 

  \item{surv}{Predicted baseline survival function when beta(w) = 0.}
  \item{cumhaz}{Baseline cumulative hazard function, -log(surv).}
  \item{hazard}{Baseline hazard function.} 
  \item{varhaz}{Variance of the baseline hazard.}
  \item{residuals}{Martingale residuals of the (lple) model.} 
  \item{std.err}{Standard error for the cumulative hazard function, if se.fit = TRUE.}

See \code{\link[survival]{survfit}} for more detail about other output values such as upper, lower, conf.type. 
Confidence interval is based on log-transformation of survival function. 
}
\author{Bingshu E. Chen}
\seealso{
 The default method for survfit \code{\link[survival]{survfit}},
#\code{\link{survfit.lple}}
}
\examples{
#
# See example in lple
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survfit }

