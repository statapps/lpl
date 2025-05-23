\name{coxScoreHess}
\alias{coxScoreHess}

\title{Calculate the Score vector / Hessian matrix for the Cox model}
\description{Calculate the Score vector or the Hessian matrix for the Cox proportional hazards model with inputs of covariates, survival outcomes and the relative risks}
\usage{
    coxScoreHess(X, y, exb, hess = FALSE)
}
\arguments{
    \item{X}{the covariate matrix from model.matrix, without the interecpt term.}
    \item{y}{y is a survival object, y = Surv(time, event).}
    \item{exb}{exb is the relative risks with exb = exp(X*beta).}
    \item{hess}{output the Hessian matrix, with hess = FALSE as the default, which outputs the score vector only,}
}
\value{An p by 1 vector of the score of the function calculated at the 
       point relative exp(X*beta). If hess = TRUE, then a list with the following three components is returned:

    \item{score}{a 1 x p score vector.}
    \item{Sigma}{a p x p matrix for the empirical varaince of the score.}
    \item{H}{a p x p hessian matrix.}
}

\details{The survival time shall be sorted from the largest to the smallest, an error will occur if y is not sorted.

     score = sum(event*(X - S1/S))

     Sigma = sum(S1*t(S1))

     H     = sum(event*(S2/S0 - S1*t(S1)/S0))

  the robust varaince can be calculated by inv(H)*Sigma*inv(H).

}

\seealso{
    \code{\link{numHessian}}
    \code{\link{numScore}}
    \code{\link{multiRoot}}
}

%this is a simple too, no example is needed here.
%\examples{}
