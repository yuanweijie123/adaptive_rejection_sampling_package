#supporting functions (might want to make one .R script per function in package)
#library(rootSolve)

# calculate derivative of a function
# instead of "grad"
cal_grad = function(x, f, lower, upper, ...) {
  eps <- (.Machine$double.eps)^(1/4)
  d <- numeric(0)
  if (x>=lower && x <= lower + eps) d <- (f(x + eps, ...) - f(x, ...))/eps
  else if (x>=upper-eps && x<=upper) d <-  (f(x, ...)- f (x - eps, ...))/eps
  else if (lower+eps <= x && x <= upper-eps) d <- (f(x + eps, ...)-f(x - eps, ...))/(2*eps)
  else stop("x is out of bounds.")
  
  # if limit doesn't exist then we need to stop
  if(is.na(d) | is.infinite(d) | is.nan(d)) stop("The derivative does not exist.")
  return(d)
}

# generate intersect zj
generate_intersect <- function(hk, dhk, xk, lb, ub){
  # length of hk
  k = length(hk)
  # various vectors
  xkl <- c(xk,0); xko <- c(0,xk)
  hkl <- c(hk,0); hko <- c(0,hk)
  dhkl <- c(dhk,0); dhko <- c(0,dhk)

  # compute zj
  zj <- (hkl-hko-xkl*dhkl+xko*dhko)/(dhko-dhkl)
  index = which(round(dhko-dhkl,8) == 0)
  zj[index] = NaN
  # set starting point and ending point
  zj[1] <- lb
  zj[k+1] <- ub
  return(zj)
}

# initialization
initialization_step <- function(g, lb, ub){
  
  # considering lb and ub is infinity
  # pre-set interval delta
  maxPoint <- try(optim(par=0.1, f = g, method = "L-BFGS-B", 
                    lower = lb, upper = ub, control=list(fnscale=-1))$par,silent=TRUE)
  if(class(maxPoint)=="try-error") stop("Error in finding maximum, please input log-concave density function.")
  
  if (lb==-Inf & ub==Inf){
    leftPoint = maxPoint-1
    rightPoint = maxPoint +1
    midPoint = maxPoint
  }
  else if (lb == -Inf & ub!=Inf){
    if (abs(maxPoint - ub)<1e-3){
      leftPoint = maxPoint-1
      rightPoint = maxPoint
      midPoint = maxPoint-1/2
    }
    else{
      delta = abs(maxPoint-ub)
      leftPoint = maxPoint-delta/2
      rightPoint = maxPoint+delta/2
      midPoint = maxPoint
    }
  }
  else if (lb != -Inf & ub == Inf){
    if (abs(maxPoint - lb)<1e-3){
      leftPoint = maxPoint
      rightPoint = maxPoint+1
      midPoint = maxPoint+1/2
    }
    else{
      delta = abs(maxPoint-lb)
      leftPoint = maxPoint-delta/2
      rightPoint = maxPoint+delta/2
      midPoint = maxPoint
    }
  }
  else{
    if (abs(maxPoint - ub) < 1e-3) {
      delta <- abs(maxPoint-lb)/2
      rightPoint <- maxPoint
      midPoint <- maxPoint - .5*delta
      left_point <- max - delta
    } else if (abs(maxPoint - lb) < 1e-3) {
      delta <- abs(maxPoint-ub)/2
      rightPoint <- maxPoint + delta
      midPoint <- maxPoint + .5*delta
      leftPoint <- maxPoint
    } else {
      delta <- min(abs(maxPoint-lb),abs(maxPoint-ub))/2
      rightPoint <- maxPoint + .5*delta
      leftPoint <- maxPoint - .5*delta
      midPoint <- maxPoint
    }
  }
  
  init <- c(leftPoint, midPoint, rightPoint)
  return(init)
}


# create upper hull in a vectorized fashion
upper_hull <- function(x,hk,dhk,xk,zk){
  i <- min(which(x < zk)-1)
  upperBound <- hk[i] + (x - xk[i]) * dhk[i]
  return(upperBound)
}

# take exponential values of upper hull
exp_upper_hull <- function(x,hk,dhk,xk,zk){
  expUpBound <- exp(upper_hull(x,hk,dhk,xk,zk))
  return(expUpBound)
}


# create lower hull in a vectorized fashion
lower_hull <- function(x,hk,dhk,xk){
  if(x < min(xk) | x > max(xk)){
    return(-Inf)
  }else{
    i <- max(which(x >= xk))
    lowerBound <- ((xk[i + 1] - x) * xk[i] + (x - xk[i]) * hk[i + 1])/(xk[i + 1] - xk[i])
    return(lowerBound)
  }
}


# sample from the envelope
draw_sample <- function(u, cumEnv, hk, xk, dhk, zk, portion){

  j <- max(which(u > cumEnv))

  if(dhk[j] == 0){
    x <- runif(1, zk[j],zk[j+1])
    return(x)
  } else {

    # Sample from uniform random
    w = runif(1)

    # Rescale sample value w to area of the selected segment
    wRescale = w/portion[j]*exp(hk[j] - xk[j]*dhk[j])*(exp(dhk[j]*zk[j+1]) -exp(dhk[j]*zk[j]))

    # Use inverse CDF of selected segment to generate a sample
    x = (1/dhk[j])*log(wRescale*portion[j]/(exp(hk[j] - xk[j]*dhk[j])) + exp(zk[j]*dhk[j]))
  }
  return(x)
}

# adaptive rejection test
rejection_test <- function(x, h, hk, dhk, xk, zk){

  # Generate random sample from uniform distribution
  w = runif(1)

  # squeeze and reject tests indicator for adding point in boolean form
  accept = FALSE
  add = FALSE

  # get rejection point for squeeze and accept test
  lowerTest = exp(lower_hull(x,hk,dhk,xk) - upper_hull(x,hk,dhk,xk,zk))

  # check if we need to keep sample points
  # check if we need to add new points to xk
  # kernel test
  if(w <= lowerTest){
    accept = TRUE
  } else {
    add = TRUE
    functionTest = exp(h(x) - upper_hull(x,hk,dhk,xk,zk))
    if(w <= functionTest){
      accept = TRUE
    }else{
      accept = FALSE
    }
  }

  # Return boolean indicator whether to accept candidate sample point
  # and update these points to xk
  return(list(acceptIndicator = accept , UpdateIndicator = add))
}

#function that checks if input is integer (as in, a whole number, not whatever is.integer() is doing...)
is.wholenumber <-	function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}

# check whether f is positive in range from var_lower to var_upper
# f is continuous
check_f_positive = function(f, lower, upper) {
  # choose a test point in interval
  # may be more efficient
  testPoint = min(upper, lower+1)
  if (f(testPoint) < 0) {
    return(FALSE)
  }

  # root
  fLower <- f(lower)
  fUpper <- f(upper)
  # if the sign of boundary values differ
  if(sign(fLower) != sign(fUpper)){
    return(FALSE)
  }

  roots <- try(rootSolve::uniroot.all(Vectorize(f), lower = lower, upper = upper))
  # if run error
  if(class(roots)=="try-error"){
    stop("Error in uniroot.all.")
  }

  # if no root in interval
  if(length(roots)==0){
    if(fLower>0) return(TRUE)
    if(fLower <= 0) return(FALSE)
  }
  else{
    rootP = roots + 1e-5
    rootM = roots - 1e-5
    if((sum(f(rootP)<0) == 0) & (sum(f(rootM)<0) == 0)){
      return(TRUE)
    }
    return(FALSE)
  }
}


# check h(x) is concave
check_concave = function(x, h) {
  # x is x_k in main function
  # h is log of g
  if (length(x) != length(h)) {
    stop("x and h should be of equal length.")
  }
  concave = TRUE
  for (i in 1:(length(x)-2)) {
    inter = h[i] + (x[i+1] - x[i]) * (h[i+2] - h[i])/(x[i+2] - x[i])
    if (inter > h[i+1]) {
      concave = FALSE
      break
    }
  }
  return(concave)
}

