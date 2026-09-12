\name{csv}
\alias{csv}
\alias{csv.default}
\alias{csv.matrix}
\alias{csv.coxph}
\alias{csv.lm}
\alias{csv.survfit}
\title{Print an Object in Comma-Separated (CSV) Style}
\description{
  A generic function that prints an object to the console in a
  comma-separated style, so the output can be easily copied into a
  spreadsheet or a plain \code{.csv} file. Methods are provided for
  matrices, and for fitted model objects of class \code{"coxph"},
  \code{"lm"}, and \code{"survfit"}. There is no method for
  \code{"glm"} objects; such objects are handled by
  \code{csv.default}, which reports the class as unsupported.
}
\usage{
csv(x, \dots)

\method{csv}{default}(x, digits = 4, \dots)

\method{csv}{matrix}(x, digits = 4, \dots)

\method{csv}{coxph}(x, \dots)

\method{csv}{lm}(x, \dots)

\method{csv}{survfit}(x, time = NULL, digits = 3, title = FALSE, \dots)
}
\arguments{
  \item{x}{an object to be printed in CSV style. For
    \code{csv.default}, \code{csv.matrix}, \code{csv.coxph} and
    \code{csv.lm}, this is the object whose values (or, for model
    objects, whose coefficient/summary table) will be printed.}
  \item{digits}{integer; the number of decimal places to which
    numeric values are rounded before printing. Defaults to \code{4}
    for the \code{default} and \code{matrix} methods, and to \code{3}
    for the \code{survfit} method.}
  \item{time}{a numeric vector of time points at which the survival
    estimate, and its confidence limits, should be reported. If
    \code{NULL} (the default), the event/censoring times stored in
    \code{fit$time} are used.}
  \item{title}{logical; if \code{TRUE}, a header line
    (\code{"Time, Survival percentage, 95\% CI lower, 95\% CI upper"})
    is printed before the table of results. Defaults to \code{FALSE}.}
  \item{\dots}{additional arguments. For \code{csv.coxph} and
    \code{csv.lm}, these are passed on to \code{csv.matrix} (e.g.\
    \code{digits}). Currently unused by \code{csv.default} and
    \code{csv.matrix} beyond generic method dispatch.}
}
\details{
  \code{csv} is an S3 generic. Its intended use is quick, ad-hoc
  export of fitted-model summaries or numeric matrices to a
  comma-separated text format that prints directly to the console
  (via \code{\link{cat}}), rather than writing to a file.

  \code{csv.default} does not print \code{x} at all; it simply lists
  the classes supported by \code{csv} (\code{coxph}, \code{glm},
  \code{matrix}, \code{survfit}) and states that the class of the
  supplied object is not currently supported. Note that, despite this
  message, no explicit \code{csv.glm} method is defined, so a
  \code{"glm"} object will also be routed to \code{csv.default} unless
  a suitable method is added.

  \code{csv.matrix} rounds \code{x} to \code{digits} decimal places
  and prints it row by row, with values separated by commas. If
  \code{x} has column names, a header row of column names is printed
  first (with a leading comma to align with the row-name column), and
  row names default to \code{1:n} if none are present. If \code{x}
  has no column names, no header or row names are printed and only
  the data values are printed for each row.

  \code{csv.coxph} takes the coefficient table and confidence
  intervals from \code{summary(x)} for a fitted Cox proportional
  hazards model (as produced by \code{\link[survival]{coxph}}),
  column-binds them together, and passes the result to
  \code{csv.matrix} for printing.

  \code{csv.lm} takes the coefficient table from \code{summary(x)}
  for a fitted linear model (as produced by \code{\link{lm}}) and
  passes it to \code{csv.matrix} for printing.

  \code{csv.survfit} prints, for each requested time in \code{time},
  the last available Kaplan--Meier estimate at or before that time:
  the time, the survival estimate, and the lower and upper confidence
  limits, each rounded to \code{digits} decimal places and separated
  by commas. Unlike the other methods, \code{csv.survfit} both prints
  its output via \code{\link{cat}} and \strong{returns} a matrix
  (invisibly not enforced) with one column per requested time,
  containing the four rounded values described above.
}
\value{
  \code{csv.default}, \code{csv.matrix}, \code{csv.coxph}, and
  \code{csv.lm} are called for their side effect of printing to the
  console and do not return a useful value (they return \code{NULL},
  invisibly, as a side effect of \code{\link{cat}}).

  \code{csv.survfit} returns a numeric matrix, \code{sout}, with four
  rows (time, survival estimate, lower confidence limit, upper
  confidence limit) and one column for each element of \code{time},
  in addition to printing the same information to the console.
}
\note{
  These functions print directly to the console with \code{\link{cat}}
  rather than writing to disk. To save the output as an actual
  \code{.csv} file, redirect the printed output (e.g.\ with
  \code{\link{sink}}) or use \code{\link{write.csv}} instead.
}
\seealso{
  \code{\link{cat}}, \code{\link{write.csv}},
  \code{\link[survival]{coxph}}, \code{\link[survival]{survfit}},
  \code{\link{lm}}
}
\examples{
## matrix method
m = matrix(1:6, nrow = 2, dimnames = list(c("r1", "r2"), c("a", "b", "c")))
csv(m)

## default method: unsupported class
csv(1:5)

\dontrun{
## coxph method (requires the survival package)
library(survival)
fit = coxph(Surv(time, status) ~ age + sex, data = lung)
csv(fit)

## lm method
fit2 = lm(mpg ~ wt + hp, data = mtcars)
csv(fit2)

## survfit method
sfit = survfit(Surv(time, status) ~ 1, data = lung)
csv(sfit, time = c(100, 200, 300), title = TRUE)
}
}
\keyword{print}
\keyword{methods}
