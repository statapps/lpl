\name{lpl-package}
\alias{lpl-package}
\alias{lpl-doc}
\docType{package}
\title{
Local Partial Likelihood Boostrap test
}
\description{
This package fits a multivariable local partial likelihood model for covariate-biomarker 
interaction with survival data.
}
\details{
"lpl" is a R package for multivariate covariate-biomarker interaction uisng local partial likelihood method.

Please use the following steps to install 'lpl' package:

1. First, you need to install the 'devtools' package. You can skip this step if you have 'devtools' installed in your R. Invoke R and then type

  install.packages("devtools")

2. Load the devtools package.

  library(devtools)

3. Install "lpl" package with R commond

  install_github("statapps/lpl")


"lpl" uses local partial likelihood to etimate covariate-biomarker interactions 
and bootstrap method to test the significance of the interactions.
}
\author{Siwei Zhang and Bingshu E. Chen

Maintainer: Bingshu E. Chen <bingshu.chen@queensu.ca>
}
\references{
1. Fan, J., Lin, H,, Zhou, Y. (2006). Local partial-likelihood estimation for lifetime data. The Annals of Statistics. 34, 290-325.

2. Liu, Y., Jiang, W. and Chen, B. E. (2015). Testing for treatment-biomarker interaction based on local partial-likelihood. Statistics in Medicine. 34, 3516-3530.

3. Zhang, S., Jiang, W. and Chen, B. E. (2016). Estimate and test of multivariate covariates and biomarker interactions for survival data based on local partial likelihood. Manuscript in preparation.
}
\keyword{lpl}
\keyword{biomarker}
\keyword{local linear model}
\seealso{
\code{coxph, }
\code{survival}
}
\examples{
# fit = lpl(y~trt+age+biomarker)
} 
