\name{asymSCB}
\alias{asymSCB}
\title{Asymptotic Critical Value for Simultaneous Confidence Bands}
\description{
Computes an asymptotic critical value for simultaneous confidence
bands (SCBs) based on the bandwidth, kernel function, and desired
confidence level.
}
\usage{
asymSCB(h, kernel, conf.int = 0.95)
}
\arguments{
\item{h}{
A positive numeric bandwidth parameter. The bandwidth should be
sufficiently small for the asymptotic approximation to be
appropriate.
}

\item{kernel}{
Character string specifying the kernel function used to calculate
the asymptotic critical value. Supported choices are
\code{"epanechnikov"}, \code{"gaussian"}, \code{"rectangular"},
\code{"triangular"}, \code{"biweight"}, \code{"cosine"}, and
\code{"optcosine"}.
}

\item{conf.int}{
Numeric value specifying the confidence level of the simultaneous
confidence band. The default is \code{0.95}.
}
}
\value{
A numeric value pn giving the asymptotic critical value for the specified
bandwidth, kernel, and confidence level.
}
\details{
The function calculates an asymptotic critical value for simultaneous
confidence bands.

First, an intermediate quantity \eqn{d_n} is calculated using the
internal function \code{.dvn}. It is given by

\deqn{
d_n = \sqrt{-2\log(h)}
(-2\log(h))^{-1/2}
\log\left(\frac{v_1}{4\pi v_0}\right)
}{
d_n = sqrt(-2 log(h)) + (-2 log(h))^(-1/2)
log(v_1 / (4 pi v_0))
}
where \eqn{v_0} and \eqn{v_1} depend on the selected kernel function.

The final asymptotic critical value is calculated as

\deqn{
p_n = d_n +
\left[\log(2)-\log{-\log(\mathit{conf.int})}\right]
(-2\log(h))^{-1/2}
}{
p_n = d_n + [log(2) - log(-log(conf.int))]
(-2 log(h))^(-1/2)
}

The internal function \code{.dvn} calculates \eqn{v_0} and \eqn{v_1}
numerically. The kernel is evaluated on a grid from \eqn{-5} to
\eqn{5} with spacing \eqn{0.001}. The quantity \eqn{v_0} is obtained
from the numerical integral of the squared kernel, while \eqn{v_1}
is obtained from the numerical integral of the squared derivative of
the kernel.

The resulting value can be used as an asymptotic multiplier in the
construction of simultaneous confidence bands.
}
\examples{
#95 percent asymptotic critical value using the Epanechnikov kernel
asymSCB(h = 0.1, kernel = "epanechnikov")

#90 percent confidence level
asymSCB(h = 0.1, kernel = "epanechnikov", conf.int = 0.90)

#Compare different kernels
asymSCB(h = 0.1, kernel = "gaussian")
asymSCB(h = 0.1, kernel = "triangular")
asymSCB(h = 0.1, kernel = "biweight")

#Effect of the bandwidth
asymSCB(h = 0.05, kernel = "epanechnikov")
asymSCB(h = 0.10, kernel = "epanechnikov")
asymSCB(h = 0.20, kernel = "epanechnikov")
}
\seealso{
\code{\link{K_func}}
}
\keyword{confidence}
\keyword{nonparametric}
\keyword{smooth}
