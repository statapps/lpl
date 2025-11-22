### score function for the log likelihood using numerical derivative
numScore = function(func, theta, h = 0.0001, ...) {
  p = length(theta)
  score = rep(NA, p)
  for(i in 1:p) {
    theta1 = theta; theta2 = theta
    theta1[i] = theta[i] - h
    theta2[i] = theta[i] + h
    score[i] = (func(theta2, ...) - func(theta1, ...))/(2*h)
  }
  return(score)
}

### numerical Jacobian function 
numJacobian = function (func, theta, m, h = 0.0001, ...) {
  p = length(theta)
  jacobian = matrix(NA, m, p)
  for (i in 1:p) {
    theta1 = theta
    theta2 = theta
    theta1[i] = theta[i] - h
    theta2[i] = theta[i] + h
    jacobian[ ,i] = (func(theta2, ...) - func(theta1, ...))/(2*h)
  }
  return(jacobian)
}

### Hessian matrix for the log pseudo likelihood using numerical derivative
numHessian = function(func, theta, h = 0.0001, method = c("fast", "easy"), ...) {
  p = length(theta)
  H = matrix(NA, p, p)
  method = match.arg(method)

  ### easy to understand, evaluate func(...) [4*p*p] times
  if(method == "easy") {
    message("Easy num hessian method\n")
    sc = function(theta, ...) {
      numScore(func = func, theta = theta, h = h, ...)
    }
    for(i in 1:p) {
      theta1 = theta; theta2 = theta
      theta1[i] = theta[i] - h
      theta2[i] = theta[i] + h
      H[i, ] = (sc(theta2, ...) - sc(theta1, ...))/(2*h)
    }
  }
  ### faster, evaluate func(theta, ...) [1 + p + p*p] times
  if (method == "fast") {
    #cat("Fast num hessian method\n")
    f0 = func(theta, ...)
    for(i in 1:p) {
      theta1 = theta; theta2 = theta
      theta1[i] = theta[i] - h
      theta2[i] = theta[i] + h
      H[i, i] = (func(theta2, ...) - 2*f0 + func(theta1, ...))/h^2
    }
    for(i in 1:(p-1)){
      for(j in (i+1):p) {
        theta1 = theta; theta2 = theta
        theta1[i] = theta[i] - h; theta1[j] = theta[j] - h
        theta2[i] = theta[i] + h; theta2[j] = theta[j] + h
        H[i, j] = (func(theta2, ...) - 2*f0 + func(theta1, ...))/(2*h^2) - (H[i, i] + H[j, j])/2
        H[j, i] = H[i, j]
      }
    }
  }
  return(H)
}

### find roots for the multiple non-linear equations.
multiRoot = function(func, theta,..., verbose = FALSE, maxIter = 50, 
        thetaUp = NULL, thetaLow = NULL,
        tol = .Machine$double.eps^0.25) {
  alpha = 0.0001
  rho = 0.5
  U = func(theta, ...)
  m = length(U)
  p = length(theta)
  if (m == p) mp = TRUE  ### for m = p
  else mp = FALSE
  convergence = 0
  mU1 = sum(U^2)
  i = 1

  while(i < maxIter) {
    J = numJacobian(func, theta, m=m, ...)
    if(mp) dtheta = solve(J, U)
    else {
      tJ = t(J)       ## m*1-(p*m) x (m*p) x (p*m) x (m*1)
      dtheta = solve(tJ%*%J, tJ%*%U)
    }
    theta0 = theta
    lambda = 1
    delta = 1
    Ud = sum(U*dtheta)
    ########Linear search
    while (delta > 0) {
      theta = theta0 - lambda*dtheta
      if(!is.null(thetaUp)) theta = ifelse(theta > thetaUp, thetaUp, theta)
      if(!is.null(thetaLow)) theta = ifelse(theta < thetaLow, thetaLow, theta)
      U = func(theta, ...)
      mU = sum(U^2)
      delta = mU-mU1-alpha*lambda*Ud
      lambda = rho*lambda
      #if(verbose) cat('delta = ', delta, '\n')
    }
    dU = abs(mU1 - mU)
    mU1 = mU
    i = i + 1
    if(verbose) cat("||U|| = ", mU, "dU = ", dU, '\n')
    if ((mU < tol) | (dU < tol)) {
      convergence = 1
      if(mU > tol) convergence = 0
      break
    }
  }
  return(list(root = theta, f.root = U, iter = i, convergence = convergence))
}

#reverse rcumsum can be avoided if scored largest time to smallest time
#rcumsum=function(x) rev(cumsum(rev(x))) # sum from last to first

coxScoreHess = function(X, y, exb, hess = FALSE) {
  ### exb = exp(X%*%beta)
  ### delta shall be sorted from largest to smallest to avoid using rcumsum.
  y1    = y[, 1]
  delta = y[, 2]
  if((y1[1] < y1[2]) | (y1[2] < y1[length(y1)])) stop("Sort survival time from the largest to the smallest")

  S0 = cumsum(exb)
  S1 = apply(exb*X, 2, cumsum)
  SX = delta * (X - S1/S0)
  score = colSums(SX)
  if(!hess) return(score)

  ### Sigma = Var(Score)
  Sigma = t(SX)%*%SX

  n = length(delta)
  p = ncol(X)
  SS1 = array(apply(S1, 1, function(x){return(x%*%t(x))}),c(p, p, n)) # ((p*p)*n)
  SS1 = aperm(array(SS1, c(p, p, n)), c(3, 1, 2))                     # (n*(p*p))

  Xt = apply(X, 1, function(x){return(x%*%t(x))})                    # X*t(X)
  X2 = array(Xt, c(p, p, n))
  ## multiply each X2(p, p, i) with exb[i],
  ## by change exb to a p*p*n array with each of ith p x p matrix = exb[i]
  X2eb = X2 * array(rep(exb, each = p*p), c(p, p, n))

  ## Sm is a upper triangular matrix of 1
  Sm = matrix(1, n, n)
  #Sm[lower.tri(Sm)] = 0
  Sm[upper.tri(Sm)] = 0

  ## calculate S2, a n*p*p array
  S2 = apply(X2eb, c(1, 2), function(x, y){return(y%*%x)}, Sm)
  H = colSums(delta*(S2/c(S0)-SS1/c(S0)^2), dims = 1)
  return(list(score = score, Sigma = Sigma, H = H))
}

### coxlogLik to calculate the logarithm of the partial likelihood for the Cox PH model
coxlogLik = function(X, y, beta, offset = NULL, sorted = FALSE) {
  ### sort the time and data
  ### if the function to be called multiple times for the same y,
  ### sort the data by y outside of this function and use sorted=TRUE
  ### to speed up the algorithm
  if(!sorted) {
    time = y[, 1]
    idx = order(time, decreasing = TRUE)
    X = X[idx, ]
    y = y[idx, ]
    if(!is.null(offset)) offset = offset[idx]
    #time = time[idx] ### once sorted, time does not contribute to coxlogLik
  }

  #####
  event = y[, 2]
  exb = exp(X%*%beta)
  if(!is.null(offset)) exb = exb*exp(offset)
  S0 = cumsum(exb)
  logLik = sum(event * log(exb/S0))
  return(logLik)
}

### survfit for the cox model, to replace basehaz, which is hard to use due to centered = TRUE/FALSE
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
  fit = cbind(time = time, n.risk = n.risk, hazard = haz, cumhaz = cumhaz, surv = surv)
  #class(fit) = 'survfit'
  return(fit)
}

