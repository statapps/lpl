\name{predict.lple}
\alias{predict.lple}
\alias{residuals.lple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
predict a lple object
}
\description{
Compute fitted values and prediction error for a model fitted by lple
}
\usage{
## S3 method for class 'lple'
\method{predict}{lple}(object, newdata, newy=NULL, ...)
\method{residuals}{lple}(object, type=c("martingale", "deviance"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a model object from the lple fit}
  \item{newdata}{optional new data at which to do predictions. If absent, predictions are for the dataframe used in the original fit}
  \item{newy}{optional new response data. Default is NULL}
  \item{type}{type of residuals, the default is a martingale residual}
  \item{...}{additional arguments affecting the predictions produced}
}
\details{
  predict.lple is called to predict object from the lple model \code{\link{lple}}.

  The default method, predict has its own help page. Use methods("predict") to get all the methods for the predict generic.
}

\value{
%%  ~Describe the value returned
predict.lple returns a list of predicted values, prediction error and residuals. 

  \item{lp}{linear predictor of beta(w)*Z, where beta(w) is the fitted regression coefficient and Z is covariance matrix.}
  \item{risk}{risk score, exp(lp). When new y is provided, both lp and risk will be ordered by survival time of the new y.}
  \item{residuals}{martingale residuals of the prediction, if available.} 
  \item{pe.mres}{prediction error based on martingale residual, if both new data and new y is provided.}
  \item{cumhaz}{cumulative hzard function.}
  \item{time}{time for cumulative hazard function. Time from new y will be used is provided}
}
\author{Bingshu E. Chen}
\seealso{
 The default method for predict \code{\link{predict}},
  
 For the Cox model prediction: \code{\link[survival]{predict.coxph}}.
#\code{\link{survfit.lple}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ predict }

