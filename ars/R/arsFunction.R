#main adaptive rejection sampling function
#inputs:
# -g			    	original function
# -n						number of samples desired
# -lb						lower bound on x axis
# -ub						upper bound on x axis
# -batchSize		number of seeds for inverse CDF
# -randomState	seed for set.Seed()
#' Adaptive Rejection Sampling
#'
#' This function is a statistical algorithm for generating samples from 
#' a univariate, log-concave density.
#'

#' @param g original function
#' @param n number of samples
#' @param lb lower bound on x axis
#' @param ub upper bound on x axis
#' @param batchSize number of seeds for inverse CDF, default = 100
#' @param randomState seed for set.seed()
#' @return Samples from target density
#' @examples g <- function(x) dnorm(x,0,1) 
#' ars(g, 5000, 2, 6, 100)

ars <- function(g, n, lb=-Inf, ub=Inf, batchSize=100, randomState=1){

	#check inputs
  assertthat::assert_that(!missing(g), !missing(n), msg = "Missing input arguments")
  assertthat::assert_that(length(n) == 1, length(lb) == 1, length(ub) == 1, length(batchSize) == 1, msg = "Inputs for n, lb, ub, and batchSize must be single numeric values")
  assertthat::assert_that(is.function(g), msg = "g must be a function input")
  assertthat::assert_that(is.numeric(n), n > 0, is.wholenumber(n), msg = "n must be a positive integer value")
  assertthat::assert_that(is.numeric(lb), is.numeric(ub), msg = "Upper and lower bound must be numeric values")
  assertthat::assert_that(lb < ub, msg = "Lower bound must be smaller than upper bound")
  assertthat::assert_that(is.numeric(batchSize), batchSize > 0, is.wholenumber(batchSize), msg = "batchSize must be a positive integer value")

	#set random seed
	set.seed(randomState)
	
	# take log
	h <- function(x){
		return(log(g(x)))
	}

  #find starting xk
  xk <- initialization_step(g, lb, ub)

  #initialize output variable
  newSample <- NULL

  #iterate until sample size is satisfied
  while(length(newSample) < n){

    # calculate hk and derivative of hk
    hk <- h(xk)
    dhk <- try(sapply(xk,cal_grad,h,lb,ub), silent=TRUE)
    if(class(dhk)=="try-error") stop("Error in calculating derivative. Try reasonable lower bound or upper bound.")

    # check for log-concavity of function
    assertthat::assert_that(check_concave(xk, hk), msg = "The provided function appears to be non-log-concave in parts of the domain. ars() can only draw samples from log-concave functions.")
    
    #intersection points of upper envelope
    zk <- generate_intersect(hk,dhk,xk,lb,ub)
    
    # 
    indexNaN = c(which(is.na(zk)))
    xk[indexNaN] = NaN
    zk <- zk[complete.cases(zk)]
    xk <- xk[complete.cases(xk)]
    
    #cumulative envelope
    #Calculate areas under exponential upper bound function for normalization purposes
    portion <-  unlist(sapply(2:length(zk),
    													function(i){integrate(Vectorize(exp_upper_hull,vectorize.args =c("x")),
    																								zk[i-1],zk[i],hk,dhk,xk,zk)})[1,])
    
    cum <- sum(portion)
    
    # Normalize and cumulation
    envelop <- portion/cum
    cumEnv <- cumsum(envelop)
    cumEnv <- c(0,cumEnv)
    
    # Sampling: Generate seeds for Inverse CDF method
    # Generate random seeds
    seeds <- runif(batchSize)
    xSample <- sapply(seeds, draw_sample,
    									cumEnv = cumEnv, hk = hk,
    									xk = xk, dhk = dhk,
    									zk = zk,portion = portion)
  
    # drop nan
    # not affect random
    xSample <- na.omit(xSample)
    
    #Rejection testing
    testResult <- try(sapply(xSample, rejection_test,
                         h = h, hk = hk, xk = xk,
                         dhk = dhk, zk = zk), silent=TRUE)
    ### detect some error
    if(class(testResult)=="try-error") {
      stop("Error in calculating integrate. Try reasonable and log-concave density.")}
    

    keepSample <- testResult[1,]
    add2xk <- testResult[2,]

    # extract index of kept samples
    # update accpeted points to sample points
    keepX <- xSample[keepSample > 0]
    newSample <- c(keepX, newSample)

    
    # check for log-concavity of function
    assertthat::assert_that(check_concave(xk, hk), msg = "The provided function appears to be non-log-concave in parts of the domain. ars() can only draw samples from log-concave functions.")
    
    # update xk
    # sort xk for the next iteration
    xkUpdated = xSample[add2xk > 0]
    xk = sort(c(xk, xkUpdated))
  }
  return(newSample[1:n])
}