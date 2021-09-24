#tests for correctness of output given known input distributions
testthat::test_that("test with normal distribution", {
	f <- dnorm
	samples <- ars(f, 1000,-Inf, Inf)
	p <- ks.test(samples, "pnorm")
	testthat::expect_gte(p$p.value, 0.1)
})

testthat::test_that("test with beta(2,2) distribution", {
	f <- function(x) dbeta(x, 2, 2)
	samples <- ars(f, 1000, 0, 1)
	p <- ks.test(samples, 'pbeta', 2, 2)
	testthat::expect_gte(p$p.value, 0.1)
})

testthat::test_that("test with gamma(2,2) distribution", {
	f <- function(x) dgamma(x, 2, 2)
	samples <- ars(f, 1000, 0, Inf)
	p <- ks.test(samples, 'pgamma', 2, 2)
	testthat::expect_gte(p$p.value, 0.1)
})

testthat::test_that("test with dchisq(3) distribution", {
	f <- function(x) dchisq(x, 3)
	samples <- ars(f, 1000, 0, Inf)
	p <- ks.test(samples, 'pchisq', 3)
	testthat::expect_gte(p$p.value, 0.1)
})


#tests to check that non-log-concavity is detected
testthat::test_that("test with t-distribution with 1 degree of freedom", {
	f <- function(x) dt(x, 1)
	testthat::expect_error(ars(f, 1000, 0, Inf), "non-log-concave")
})

testthat::test_that("test with default Cauchy distribution", {
	f <- function(x) dcauchy(x)
	testthat::expect_error(ars(f, 1000, 0, Inf), "non-log-concave")
})

testthat::test_that("test with default log-normal distribution", {
	f <- function(x) dlnorm(x)
	testthat::expect_error(ars(f, 1000, 0, 10), "non-log-concave")
})

testthat::test_that("test with chi-squared distribution with 1 degree of freedom", {
	f <- function(x) dchisq(x, 1)
	testthat::expect_error(ars(f, 1000, 1, 10), "non-log-concave")
})