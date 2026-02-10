### csv class: print in a csv style
csv = function(x, ...) UseMethod("csv")

csv.default = function(x, digits = 4, ...) {
  cat("csv supports the following class: ")
  cat("coxph, glm, matrix, survfit\n\n")
  cat("the class", class(x), "is not supported now.\n")
}

csv.matrix = function(x, digits = 4, ...) {
  x = round(x, digits)
  n = nrow(x)
  p = ncol(x)
  rname = rownames(x)
  cname = colnames(x)
  if(!is.null(cname)) {
    for(j in 1:p) cat(", ", cname[j])
    cat("\n")
    if(is.null(rname)) rname = 1:n
  }
  for(i in 1:n) {
    if(!is.null(rname)) cat(rname[i], ', ', sep = "")
    for(j in 1:(p-1)) cat(x[i, j], ', ', sep = "")
    cat(x[i, p], "\n")
  }
}

csv.coxph = function(x, ...) {
  sf = summary(x)
  beta = sf$coefficients
  ci   = sf$conf.int
  out = cbind(beta, ci)
  csv(out, ...)
}

csv.lm = function(x, ...) {
  out = summary(x)$coefficients
  csv(out, ...)
}

csv.survfit = function(fit, time = NULL, digits = 3, title = FALSE) {
  if(title) cat("Time, Survival percentage, 95% CI lower, 95% CI upper\n")
  if(is.null(time)) time = fit$time
  
  sfit = summary(fit)
  sout = NULL
  out = cbind(sfit$time, sfit$surv, sfit$lower, sfit$upper)
  K = length(time)
  for(i in 1:K) {
    out1 = out[out[, 1]<= time[i], ,drop = FALSE]
    #print(out1)
    p = nrow(out1)
    #print(p)
    sv = out1[p, ]
    sv = round(sv, digits)
    #sv[1] = tx
    cat(round(sv[1], digits), ",", sv[2], ",", sv[3], ",", sv[4], "\n")
    sout = cbind(sout, sv)
  }
  return(sout)
}

### find odds ratio (OR) and 95% CI for a logistic regression
### oddsRatio(x, conf.int = 0.95, ...), where x = glm()
oddsRatio = function(x, conf.int = 0.95, ...) {
  # p-values (Wald tests)
  pvals <- summary(x)$coefficients[, "Pr(>|z|)"]
  
  results <- cbind(
    OR = exp(coef(x)),
    CI_lower = exp(confint.default(x, level = conf.int)[,1]),
    CI_upper = exp(confint.default(x, level = conf.int)[,2]), 
    p_value = pvals
  )
  results
}


