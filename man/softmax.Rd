\name{softmax}
\alias{softmax}
\title{Softmax Activation Function}
\description{
Computes the softmax transformation for a numeric vector or matrix.
The softmax function converts real-valued inputs into probabilities
that sum to 1. A numerically stable implementation is used by subtracting
the maximum value before exponentiation.

For matrix input, the softmax transformation is applied row-wise, so each
row is converted into a probability distribution whose elements sum to 1.
}
\usage{
softmax(x)
}
\arguments{
  \item{x}{
    A numeric vector or numeric matrix.
  }
}
\value{
If \code{x} is a vector, returns a numeric vector of the same length
containing softmax probabilities that sum to 1.

If \code{x} is a matrix, returns a numeric matrix of the same dimensions
where each row contains softmax probabilities summing to 1.
}
\details{
The softmax function is defined as:

\deqn{
\mathrm{softmax}(x_i) = \frac{e^{x_i}}{\sum_j e^{x_j}}
}

To improve numerical stability, the implementation subtracts the maximum
value from the input before exponentiation:

\deqn{
\frac{e^{x_i - \max(x)}}{\sum_j e^{x_j - \max(x)}}
}
}
\examples{
# Vector input
x <- c(1.2, 0.3, 2.4)
softmax(x)

# Matrix input
m <- matrix(c(1, 2, 3,
              0, 1, 0),
            nrow = 2,
            byrow = TRUE)

softmax(m)

# Row sums equal 1
rowSums(softmax(m))
}
\seealso{
\code{\link{exp}}
}
\author{
Bingshu E. Chen
}
\keyword{math}
