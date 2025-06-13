#library(lpl)
#set.seed(29)
rmst = function(object, ...) UseMethod("rmst")

rmst.default = function(object, ...) {
  message("Rrestricted mean survival time method for class '", class(object), "' is not available.") 
  message("Calculate the RMST for the KM estimate instead.")
  message("You can create function rmst.", class(object), "(object, ...) to implement method for this class.")
  return(rmst(object$y, ...))
}

### Restricted mean survival time based on a hazard (or cumulative hazard) function
rmstFit = function(tau, h0 = NULL, H0 = function(x){x}) {
  rms = integrate(psurv, 0, tau, h0=h0, H0 = H0, low.tail = FALSE)$value
  return(rms)
}

rmst.coxph = function(object, newdata=NULL, linear.predictors = NULL, tau=NULL, ...) {
  H    = basehaz(object, centered = FALSE)
  time = H$time
  chaz = H$hazard
  beta = object$coefficients
  p    = length(beta)

  rmsfunlp = function(lp, tm, chz, tau) {
    Ht = approxfun(tm, chz*exp(lp), rule = 2)
    return(rmstFit(tau, H0 = Ht))
  }

  if(is.null(newdata)) {
    if(is.null(linear.predictors)) return(rmst(object$y, tua = tau))
    else return(sapply(linear.predictors, rmsfunlp, tm = time, chz = chaz, tau=tau))
  }
  if(!is.null(linear.predictors)) stop("Only one of the newdata or the linear.predictors can be non-null.")
  if(is(newdata, "matrix") == FALSE) 
    newdata =as.matrix(newdata, ncol = p)
  else if (ncol(newdata) != p) stop("newdata must be a matrix with ", p, "columns.")

  if(is.null(tau)) tau = max(time)
  rmsfun = function(x, beta, tm, chz, tau) {
    #lp = sum(x*beta)
    return(rmsfunlp(sum(x*beta), tm, chz, tau))
  }
  rms = apply(newdata, 1, rmsfun, beta = beta, tm = time, chz = chaz, tau = tau) 
  return(rms)
}

rmst.Surv = function(object, tau = NULL, ...) {
  if(is.null(tau)) tau = max(object[, 1])
  sf = survfit(object ~ 1)
  St = approxfun(sf$time, exp(-sf$cumhaz), rule = 2)
  return(integrate(St, 0, tau)$value)
}

### psurv from dnn package, can be deleted in the future.
#.psurv = function(q, h0 = NULL, H0 = function(x){x}, low.tail=TRUE, log.p=FALSE) {
#  H = function(t) { integrate(h0, 0, t, subdivisions = 500L)$value }
#  if(!is.null(h0)) Ht = vapply(q, H, 1)
#  else Ht = vapply(q, H0, 1)
#
#  if(low.tail) {
#    if (log.p) return(Ht)
#    else return(1-exp(-Ht))
#  } else {
#    if (log.p) return(-Ht)
#    else return(exp(-Ht))
#  }
#}

### example for rmst
#n     = 29
#x     = rnorm(n)
#time  = rexp(n)
#event = rbinom(n, 1, 0.7)
#y     = Surv(time, event)
#fit   = coxph(y~x)
### create a test dataset
#xt    = as.matrix(runif(3), ncol = 1, nrow = 3) 

#print(rmst(y, tau = 5))
#print(rmst(fit, newdata = xt, tau = 2))
#print(rmst(fit, linear.predictors = c(1, 2, 3), tau = 2))
#print(rmst(fit, tau = 2))
#obj = list(y = y, x = x)
