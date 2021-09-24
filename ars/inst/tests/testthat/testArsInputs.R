#tests for different variations of number of inputs and input types

#necessary variables
#functional inputs
g <- function(x) dnorm(x,0,1) #functional sample input for g
n <- 5000 #functional sample input for n
lb <- -10 #functional sample input for lb
ub <- 10 #functional sample input for ub
lbInf <- -Inf #functional sample input for lb
ubInf <- Inf #functional sample input for ub
batchSize <- 100 #functional sample input for batchSize
randomState <- 0 #functional sample input for randomState

#broken inputs
gBroken <- runif(100) #non-function input for g
nBroken1 <- seq_len(3) #vector input for n
nBroken2 <- 0.3 #non-integer input for n
nBroken3 <- -3 #negative input for n
nBroken4 <- "4" #character input for n
lbBroken1 <- "3" #character input for lb
lbBroken2 <- nBroken1 #vector input for ub
ubBroken1 <- "3" #character input for ub
ubBroken2 <- nBroken1 #vector input for ub
batchSizeBroken <- 

#faulty function calls -- number of inputs
testthat::test_that("ars() can't be called with wrong number of inputs", {
	testthat::expect_error(ars(), "Missing input arguments")
	testthat::expect_error(ars(g), "Missing input arguments")
	testthat::expect_error(ars(g, n, 0, 5, 100, 3, 6), "unused argument")
})

#acceptable function calls -- number of inputs
#test_that("ars() runs properly with correct number of inputs", {
#	expec
#})
#ars(g, 3, 0, 10)

#faulty function calls -- format of inputs
testthat::test_that("ars() cannot be called with inputs of the wrong format", {
	testthat::expect_error(ars(gBroken, n), "must be a function")
	testthat::expect_error(ars(g, nBroken1), "single numeric values")
	testthat::expect_error(ars(g, nBroken2), "positive integer value")
	testthat::expect_error(ars(g, nBroken3), "positive integer value")
	testthat::expect_error(ars(g, nBroken4), "positive integer value")
	testthat::expect_error(ars(g, n, "3"), "must be numeric values")
})


#faulty function calls -- nonsensical inputs
#test_that()