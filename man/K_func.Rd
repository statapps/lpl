\name{K_func}
\alias{K_func}
\title{Kernel Function}
\description{
Computes a kernel function evaluated at the difference between a target
value \code{w} and evaluation point \code{u}, with bandwidth \code{h}.
Several commonly used kernel functions are supported, including
Epanechnikov, Gaussian, rectangular, triangular, biweight, cosine, and
optimal cosine kernels.
}
\usage{
K_func(w, u, h,
kernel = c("epanechnikov", "gaussian", "rectangular",
"triangular", "biweight", "cosine", "optcosine"))
}
\arguments{
\item{w}{
Numeric value or vector at which the kernel is evaluated.
}

\item{u}{
Numeric value or vector representing the evaluation point or center
of the kernel.
}

\item{h}{
A positive numeric bandwidth parameter controlling the scale of the
kernel.
}

\item{kernel}{
Character string specifying the kernel function to use. One of
\code{"epanechnikov"}, \code{"gaussian"}, \code{"rectangular"},
\code{"triangular"}, \code{"biweight"}, \code{"cosine"}, or
\code{"optcosine"}. The default is \code{"epanechnikov"}.
}
}
\value{
A numeric value or vector containing the kernel density weights
evaluated at \code{w - u}. Values outside the support of a compactly
supported kernel are assigned a small positive value
(\code{1e-40}) rather than zero.
}
\details{
The function computes
\code{x = w - u} and evaluates the selected kernel using bandwidth
\code{h}.

For the Gaussian kernel, the function uses a normal density with
standard deviation \code{h} and evaluates it only when
\code{|w-u| < 5h}. Outside this range, the value is set to
\code{1e-40}.

The other kernels have compact support on \code{|w-u| < h}. Their
kernel functions are:

\itemize{
\item \code{epanechnikov}: \eqn{K(x) = \frac{3}{4h}
\left[1-(x/h)^2\right]}.
\item \code{gaussian}: \eqn{K(x) =
\frac{1}{h}\phi(x/h)}, where \eqn{\phi} is the standard normal
density.
\item \code{rectangular}: \eqn{K(x) = \frac{1}{2h}}.
\item \code{triangular}: \eqn{K(x) =
\frac{1-|x|/h}{h}}.
\item \code{biweight}: \eqn{K(x) =
\frac{15}{16h}\left[1-(x/h)^2\right]^2}.
\item \code{cosine}: \eqn{K(x) =
\frac{1+\cos(\pi x/h)}{2h}}.
\item \code{optcosine}: \eqn{K(x) =
\frac{\pi}{4h}\cos(\pi x/(2h))}.
}

The implementation uses \code{match.arg} to validate the requested
kernel type.
}
\examples{

#Evaluate the Epanechnikov kernel

K_func(w = 0.5, u = 0, h = 1)

#Evaluate the Gaussian kernel

K_func(w = 0.5, u = 0, h = 1, kernel = "gaussian")

#Compare several kernel functions

x <- seq(-1, 1, length.out = 100)

K.epan <- K_func(x, 0, 1, kernel = "epanechnikov")
K.tri <- K_func(x, 0, 1, kernel = "triangular")
K.bi <- K_func(x, 0, 1, kernel = "biweight")

#Vectorized evaluation

w <- c(-0.5, 0, 0.5, 1)
K_func(w, u = 0, h = 1, kernel = "epanechnikov")
}
\seealso{
\code{\link[stats]{dnorm}},
\code{\link[stats]{density}}
}
\keyword{nonparametric}
\keyword{smooth}
\keyword{kernel}
