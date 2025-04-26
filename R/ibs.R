ibs = function(object, ...) UseMethod("ibs")

ibs.default = function(object, ...) {
  message("IBS Method for class '", class(object), "' is not available. Calculate the IBS for the KM estimate instead.")
  y  = object$y
  sf = survfit(y~1)
  tm = sf$time
  St = exp(-sf$cumhaz)
  ds = diff(c(0, tm))
  n  = length(tm)
  bs = rep(0, n)
  for(i in 1:n) bs[i] = brierScore(y, rep(St[i], n), tm[i])
  return(sum(bs*ds))
}

ibs.Surv = function(object, survProb, ...) {
  ### sort by time
  idx    = order(object[, 1])
  St     = survProb[idx, ]
  y      =   object[idx, ]
  time   = y[, 1]
  status = y[, 2]

  n = length(time)
  Indicator = matrix(0, nrow = n, ncol = n)   ### at risk indicator, subject sorted by time
  Indicator[lower.tri(Indicator)] = 1         ### subject 2 will not be at risk at time 1

  Cens = matrix(status, nrow = n, ncol = n)   ### Censoring matrix, for i with event,
  Cens[lower.tri(Cens)] = 1                   ### I(T_i > t) will be observed
                                              ### For all i, Cens will be 1 for t < T_i
  G = IPCW(y)
  G_M = matrix(G, nrow = n, ncol = n)         ### G(t) for censoring, from t to tau
  G_M[lower.tri(G_M)] = 0
  G_M = G_M + (Indicator %*% diag(G))
  BSM = ((Indicator - St)^2) * Cens/G_M       ### (I(T>t) - S(t))^2 * Delta/G
  BSt = apply(BSM, 2, mean, na.rm = TRUE)
  dt  = diff(c(0, time))
  return(sum(BSt * dt))
}

ibs.coxph = function(object, newdata = NULL, newy = NULL, ...) {
  y = object$y
  if (is.null(newdata)) {
    newdata = object$x
    newy = y
  }
  else if (is.null(newy))
    stop("To calculate Brier score for newdata, newy cannot be NULL.")
  if (length(newdata[, 1]) != length(newy[, 1]))
    stop("New data and new y must have same number of observations.")

  ### sort the newdata and newy for the IBS calculation    
  idx = order(newy[, 1])
  Y   = newy[idx, ]
  X   = newdata[idx, ]

  hr  = exp(X%*%(object$coefficients))    ### hazard ratio = exp(linear predictor)
  HZ  = basehaz(object, centered = FALSE) ### to find the baseline cumhaz
  tau = max(HZ$time, Y[, 1])
  chz = approxfun(c(0, HZ$time, tau), c(0, HZ$hazard, max(HZ$hazard)))
  St  = exp(-hr %*% chz(Y[, 1]))
  return(ibs(Y, St))
}

ibs.lple = function(object, newdata=NULL, newy = NULL, ...) {
  # When one does not use new data to calculate integrated Brier score, 
  # original data and y will be used 
  y = object$y
  if (is.null(newdata)) {
    newdata = object$X
    newy = y
  } else if(is.null(newy)) 
    stop("To calculate Brier score for newdata, newy cannot be NULL.")

  if(length(newdata[, 1]) != length(newy[, 1])) 
    stop("New data and new y must have same number of observations.")
  
  sf = predict(object, newdata, newy)
  ## Matrix of survival function S(t) rows: subjects, columns: time
  ## where sf$risk = exp(beta(w)*z + g(w))
  S = exp(-sf$risk %*% t(sf$cumhaz))
  time = newy[, 1]
  status = newy[, 2]
  n = length(time)

  Indicator = matrix(0, nrow=n, ncol=n)
  Indicator[lower.tri(Indicator)]=1
  Cens = matrix(status, nrow=n, ncol=n)
  Cens[lower.tri(Cens)]=1

  # fit km curve for censoring
  G_fit = survfit(Surv(time, 1-status)~1)
  G = .appxf(G_fit$surv, x=G_fit$time, xout = time)

  G_M = matrix(G, nrow=n, ncol=n)
  G_M[lower.tri(G_M)]=0
  G_M = G_M+(Indicator %*% diag(G))

  BSM = ((Indicator-S)^2)*Cens/G_M
  BSt = apply(BSM, 2, mean, na.rm = TRUE)
    
  dt = diff(c(0, time))
  return(sum(BSt*dt))
}

### helper functions: IPCW, brierScore
ipcw = function(time, event) {IPCW(Surv(time, event))}
IPCW = function(object) {
  time = object[, 1]
  cnsr = 1 - object[, 2]
  Gf   = survfit(Surv(time, cnsr)~1)
  Gfun = approxfun(Gf$time, Gf$surv)
  G    = Gfun(time)    ### fix G(t) = 0 problem
  G    = ifelse(G>0, G, 1/length(time))
  return(G)
}

brierScore = function(object, St, tau) {
# formula I(T<=tau)*(0 - S(tau))^2*delta/G(T) + I(T>tau)*(1-S(tau))^2/G(tau)
# object: Surv(time, event), St: Survival function at tau.
  idx = order(object[, 1])
  St  = St[idx]
  y   = object[idx, ]

  time   = y[, 1]
  status = y[, 2]

  Itau = ifelse(time>tau, 1, 0)
  G    = IPCW(y)
  Gf   = approxfun(time, G)
  Gtau = Gf(tau)

  bs   = (1-Itau)*(St)^2*status/G + Itau*(1-St)^2/Gtau
  return(mean(bs))
}
