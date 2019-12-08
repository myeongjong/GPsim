####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several tests for functions regarding covariance functions.
###
####################################################################################


test_that("Isotropic exponential covariance function", {
  expect_equal(cov_expo_iso(locs = matrix(c(0, log(2)/sqrt(2), log(2)/sqrt(2), 0), 2, 2, byrow = T), covparms = c(1, 1)), matrix(c(1, 0.5, 0.5, 1), 2, 2, byrow = T))

  expect_equal(cov_expo_iso(locs = matrix(c(0, log(2)/sqrt(2), log(2)/sqrt(2), 0), 2, 2, byrow = T), covparms = c(1, 0.5)), matrix(c(1, 0.25, 0.25, 1), 2, 2, byrow = T))

  expect_equal(cov_expo_iso(locs = matrix(c(0, log(2)/sqrt(2), log(2)/sqrt(2), 0), 2, 2, byrow = T), covparms = c(2, 0.5)), matrix(c(2, 0.5, 0.5, 2), 2, 2, byrow = T))
})


test_that("Isotropic exponential covariance function - example", {
  expect_equal(round(cov_expo_iso(locs = expand.grid(c(0.25, 0.75), c(0.25, 0.75)), covparms = c(1, 0.1)), 5),
               matrix(c(1.00000, 0.00674, 0.00674, 0.00085,
                        0.00674, 1.00000, 0.00085, 0.00674,
                        0.00674, 0.00085, 1.00000, 0.00674,
                        0.00085, 0.00674, 0.00674, 1.00000), 4, 4, byrow = T))
})


test_that("Anisotropic exponential covariance function", {
  expect_equal(cov_expo_aniso(locs = matrix(c(0, log(2)/sqrt(2), log(2)/sqrt(2)/2, 0), 2, 2, byrow = T), covparms = c(1, 1, 2)), matrix(c(1, 0.5, 0.5, 1), 2, 2, byrow = T))

  expect_equal(cov_expo_aniso(locs = matrix(c(0, log(2)/sqrt(2), log(2)/sqrt(2)/5, 0), 2, 2, byrow = T), covparms = c(1, 1, 5)), matrix(c(1, 0.5, 0.5, 1), 2, 2, byrow = T))

  expect_equal(cov_expo_aniso(locs = matrix(c(0, log(2)/sqrt(2), log(2)/sqrt(2)/5, 0), 2, 2, byrow = T), covparms = c(1, 0.5, 5)), matrix(c(1, 0.25, 0.25, 1), 2, 2, byrow = T))

  expect_equal(cov_expo_aniso(locs = matrix(c(0, log(2)/sqrt(2), log(2)/sqrt(2)/5, 0), 2, 2, byrow = T), covparms = c(2, 0.5, 5)), matrix(c(2, 0.5, 0.5, 2), 2, 2, byrow = T))
})


test_that("Anisotropic exponential covariance function - example", {
  expect_equal(round(cov_expo_aniso(locs = expand.grid(c(0.25, 0.75), c(0.25, 0.75)), covparms = c(1, 0.1, 10))^(0.1), 5),
               matrix(c(1.00000, 0.00674, 0.60653, 0.00657,
                        0.00674, 1.00000, 0.00657, 0.60653,
                        0.60653, 0.00657, 1.00000, 0.00674,
                        0.00657, 0.60653, 0.00674, 1.00000), 4, 4, byrow = T))
})






