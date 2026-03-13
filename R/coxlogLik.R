### return a partial likelihood or full likelihood of a Cox PH model
### this is provided by Bingshu (November 27, 2025)
coxlogLik = function(X, y, beta, offset = NULL, H0 = NULL, h0 = NULL, 
                     partial = TRUE, sorted = FALSE) {
  ### sort the time and data
  ### if the function to be called multiple times for the same y,
  ### sort the data by y outside of this function and use sorted=TRUE
  ### to speed up the algorithm
  
  ## linear predictor
  lp = as.vector(X%*%beta)
  if(!is.null(offset)) lp = lp + offset
  
  ## sort data by decreasing time if needed
  if(!sorted) {
    time = y[, 1]
    idx  = order(time, decreasing = TRUE)
    time = time[idx]
    lp   = lp[idx]
    y    =  y[idx, ]
  }
  
  exb = exp(lp)
  S0 = cumsum(exb)
  event = y[, 2]
  
  ### Partial log-likelihood
  if(partial) {
    logLik = sum(event * (lp - log(S0)))
    return(logLik)
  }
  
  
  ## Full log-likelihood
  if(is.null(H0)) {
    ## Breslow baseline hazard estimate
    ht = event/S0
    Ht = rev(cumsum(rev(ht)))
  } else {
    ## use-supplied cumulative baseline hazard
    Ht = H0(time)
    if(is.null(h0)) stop("h0 cannot be null if H0 is provided. Try { coxcumhaz } to find h0\n")
    else ht = h0(time)

  }
  
  ## avoid log(0) when event = 0
  ht[event==0] <- 1
  
  logLik = sum(event*(log(ht)+lp) - Ht*exb)
  return(logLik)
}

### cumulative hazard function  for the cox model, to replace basehaz, which is hard to use due to centered = TRUE/FALSE
#n = 7
#y = Surv(runif(n), rbinom(n, 1, 0.7))
#x = runif(n)
#fit = coxph(y~x)
#print(basehaz(fit))
#print(coxcumhaz(y, predict(fit)))
coxcumhaz = function(y, linear.predictors = NULL, sorted = FALSE) {
  n = length(y[, 1])
  
  if(!is.null(linear.predictors)) exb = exp(linear.predictors)
  else exb = rep(1, n)
  
  if(!sorted) {
    idx = order(y[, 1])
    exb = exb[idx]
    y = y[idx, ]
    time = y[, 1]
  }
  event = y[, 2]
  events = sum(event)
  n.risk = n:1
  
  S0 = .rcumsum(exb)
  haz = event/S0
  cumhaz = cumsum((haz))
  surv = exp(-cumhaz)
  ht = approxfun(time, haz)
  Ht = approxfun(time, cumhaz)
  St = approxfun(time, surv)
  return(list(ht = ht, Ht = Ht, St = St))
}
