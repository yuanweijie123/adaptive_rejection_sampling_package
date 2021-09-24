#unit tests for supporting functions

#variables needed for tests
g <- function(x) dnorm(x,0,1) #normal distribution as function input
lb <- 0 #functional lower bound
lb2 <- 1 #another lower bound
lb3 <- -1 
ub <- 10 #functional upper bound
deriv0 <- 0 #functional derivative of 0
deriv1 <- 1 #functional derivative of 1

#test cal_grad()
testthat::test_that("cal_grad() returns an appropriate derivative", {
	testthat::expect_equal(cal_grad(0, g, lb3, ub), 0)
	testthat::expect_equal(cal_grad(0.5*pi, sin, lb, ub), 0)
})

#test generate_intersect()
testthat::test_that("generate_intersect returns appropriate intersects between points with known derivatives and has lower and upper bound as first and last values", {
	testthat::expect_equal(generate_intersect(1, deriv0, 2, lb, ub), c(lb, ub))
	testthat::expect_equal(generate_intersect(c(lb, lb+1), c(deriv0, deriv0), c(ub-1, ub), lb, ub), c(lb, NaN, ub))
	testthat::expect_equal(generate_intersect(c(lb, lb+1), c(deriv0, deriv1), c(ub-1, ub), lb, ub), c(lb, ub-1, ub))
})

#test initialization_step()
testthat::test_that("initialization_step() produces an initial step consisting of the lower bound, half of the difference between upper and lower bound, and evenly spaced points in between", {
	testthat::expect_equal(initialization_step(g, lb, ub), c(lb, (ub-lb)/4+lb, (ub-lb)/2+lb))
	testthat::expect_equal(initialization_step(g, lb2, ub), c(lb2, (ub-lb2)/4+lb2, (ub-lb2)/2+lb2))
})

#test upper_hull()


#test exp_upper_hull()


#test lower_hull()


#test draw_sample()


#test rejection_test()


#test is.wholenumber()
testthat::test_that("is.wholenumber() returns FALSE for non-integers", {
	testthat::expect_equal(is.wholenumber(2.3), FALSE)
	testthat::expect_equal(is.wholenumber(4), TRUE)
	testthat::expect_equal(is.wholenumber(3.0003, tol = 0.0001), FALSE)
	testthat::expect_equal(is.wholenumber(3.0003, tol = 0.001), TRUE)
})


#test check_f_positive()
testthat::test_that("check_f_positive() returns TRUE for positive function", {
	testthat::expect_equal(check_f_positive(exp,0,1), TRUE)
	testthat::expect_equal(check_f_positive(dnorm,-1e8,1e8), TRUE)
})

#test check_concave()
testthat::test_that("check_concave returns TRUE for non-concave", {
	testthat::expect_equal(check_concave(c(1,2,3),c(2,3,4)), TRUE)
	testthat::expect_equal(check_concave(c(-1,0,1),c(1,3,1)), TRUE)
	testthat::expect_equal(check_concave(c(-1,0,1),c(1,-3,1)), FALSE)
})

