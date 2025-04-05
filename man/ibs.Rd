\name{ibs}
\alias{brier}
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
Create the Inverse Probability of Censoring Weighting (IPCW) using the Kaplan-Meier (KM) method.
print are used to provide a short summary of lple outputs.
}
\usage{
brier(object, St, tau)
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
  \code{\link[survival]{survfit}} is called to fit a KM model for the censoring time.
}
\author{Bingshu E. Chen}
\seealso{
 The ibs function is used in \code{\link{brier}} to calculate the brier score and \code{\link{ibs}} to calculate the 
 integrated brier score.
}
\examples{
# See example in brier ibs
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IPCW }
\keyword{ Brier Score }
