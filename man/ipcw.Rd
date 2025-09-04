\name{IPCW}
\alias{IPCW}
\alias{ipcw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inverse probability of censoring weighting (IPCW)
}
\description{
Create the Inverse Probability of Censoring Weighting (IPCW) using the Kaplan-Meier (KM) method.
print are used to provide a short summary of lple outputs.
}
\usage{
IPCW(object)
ipcw(time, event)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a survival object created by Surv(time, event).}
  \item{time}{the survival time.}
  \item{event}{the status indicator, normally 0=alive, 1=dead.}
}
\value{A vector for the survival function of the censoring time is returned.}
\details{
  \code{\link[survival]{survfit}} is called to fit a KM model for the censoring time.
}
\author{Bingshu E. Chen}
\seealso{
 The IPCW function is used in \code{\link{brierScore}} to calculate the brier score and \code{\link{ibs}} to calculate the 
 integrated brier score.
}
\examples{
# See example in brier ibs
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IPCW }
\keyword{ Brier Score }
