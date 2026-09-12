\name{oddsRatio}
\alias{oddsRatio}
\title{Odds Ratios and Confidence Intervals for a Logistic Regression}
\description{
  Computes odds ratios, their confidence intervals, and the
  associated Wald test p-values for the coefficients of a fitted
  logistic regression model.
}
\usage{
oddsRatio(x, conf.int = 0.95, \dots)
}
\arguments{
  \item{x}{a fitted model object of class \code{"glm"}, typically the
    result of a call to \code{\link{glm}} with \code{family =
    binomial}.}
  \item{conf.int}{numeric; the confidence level for the confidence
    interval on each odds ratio. Defaults to \code{0.95} (a 95\%
    confidence interval).}
  \item{\dots}{additional arguments, currently unused.}
}
\details{
  \code{oddsRatio} exponentiates the estimated coefficients of
  \code{x} to obtain odds ratios, and exponentiates the corresponding
  Wald-type confidence limits from \code{\link{confint.default}} (that
  is, limits based on asymptotic normality, \code{estimate +/- z *
  SE}, rather than profile-likelihood limits from
  \code{\link{confint}}) to obtain confidence intervals on the odds
  ratio scale. The two-sided Wald test p-value for each coefficient,
  taken from \code{summary(x)$coefficients[, "Pr(>|z|)"]}, is included
  alongside the odds ratio and its confidence interval.

  The intercept, if present, is included in the output along with all
  other coefficients; its exponentiated value does not have an odds
  ratio interpretation in the usual sense and can be ignored or
  removed by the user if not needed.
}
\value{
  A numeric matrix with one row per coefficient in \code{x} (in the
  same order as \code{coef(x)}) and four columns:
  \item{OR}{the odds ratio, \code{exp(coef(x))}.}
  \item{CI_lower}{the lower confidence limit for the odds ratio, at
    the level given by \code{conf.int}.}
  \item{CI_upper}{the upper confidence limit for the odds ratio, at
    the level given by \code{conf.int}.}
  \item{p_value}{the Wald test p-value for the corresponding
    coefficient, from \code{summary(x)}.}
}
\note{
  Confidence intervals are computed with \code{\link{confint.default}}
  (Wald-based), not with the (generally preferred, but slower)
  profile-likelihood intervals from \code{\link{confint}}. For models
  with small samples or coefficients near a boundary, the two methods
  can disagree noticeably.
}
\seealso{
  \code{\link{glm}}, \code{\link{coef}}, \code{\link{confint.default}},
  \code{\link{confint}}, \code{\link{summary.glm}}
}
\examples{
## a simple logistic regression
fit = glm(am ~ mpg + wt, data = mtcars, family = binomial)
oddsRatio(fit)

## with a 90\% confidence interval instead of the 95\% default
oddsRatio(fit, conf.int = 0.90)
}
\keyword{models}
\keyword{regression}

