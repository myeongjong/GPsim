####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several tests for simulation functions.
###
###   Description:
###
###
####################################################################################

test_that("locs, n, and p", {
  expect_equal(.checkargs_locsnp(locs = matrix(0, 3, 2), n = NULL, p = NULL),
               .checkargs_locsnp(locs = matrix(0, 3, 2), n = 3, p = 2))
  expect_equal(.checkargs_locsnp(locs = matrix(0, 3, 2), n = NULL, p = NULL),
               .checkargs_locsnp(locs = matrix(0, 3, 2), n = 3, p = NULL))
  expect_equal(.checkargs_locsnp(locs = matrix(0, 3, 2), n = NULL, p = NULL),
               .checkargs_locsnp(locs = matrix(0, 3, 2), n = NULL, p = 2))
  expect_equal(as.matrix(.checkargs_locsnp(locs = 'grid', n = 3^2, p = 2)$locs),
               as.matrix(.checkargs_locsnp(locs = expand.grid(c(0, 0.5, 1), c(0, 0.5, 1)), n = NULL, p = NULL)$locs))
  expect_identical(tryCatch(.checkargs_locsnp(locs = matrix(0, 3, 2), n = NULL, p = 1), error = function(e) TRUE), TRUE)
  expect_identical(tryCatch(.checkargs_locsnp(locs = matrix(0, 3, 2), n = 4, p = NULL), error = function(e) TRUE), TRUE)
  expect_identical(tryCatch(.checkargs_locsnp(locs = matrix(0, 3, 2), n = 3, p = 1), error = function(e) TRUE), TRUE)
  expect_identical(tryCatch(.checkargs_locsnp(locs = matrix(0, 3, 2), n = 4, p = 2), error = function(e) TRUE), TRUE)
})

test_that("meanmodel", {
  expect_equal(.checkargs_meanmodel(meanmodel = function(locs, meanparms) rep(0, nrow(locs)), meanparms = 0, locs = matrix(0, 3, 2), n = 3, p = 2),
               .checkargs_meanmodel(meanmodel = matrix(0, 3, 1), meanparms = 0, locs = matrix(0, 3, 2), n = 3, p = 2))
  expect_equal(.checkargs_meanmodel(meanmodel = function(locs, meanparms) rep(0, nrow(locs)), meanparms = 0, locs = matrix(0, 3, 2), n = 3, p = 2),
               .checkargs_meanmodel(meanmodel = as.data.frame(matrix(0, 3, 1)), meanparms = 0, locs = matrix(0, 3, 2), n = 3, p = 2))
  expect_equal(.checkargs_meanmodel(meanmodel = function(locs, meanparms) rep(0, nrow(locs)), meanparms = 0, locs = matrix(0, 3, 2), n = 3, p = 2),
               .checkargs_meanmodel(meanmodel = matrix(0, 1, 3), meanparms = 0, locs = matrix(0, 3, 2), n = 3, p = 2))
  expect_equal(.checkargs_meanmodel(meanmodel = function(locs, meanparms) rep(0, nrow(locs)), meanparms = 0, locs = matrix(0, 3, 2), n = 3, p = 2),
               .checkargs_meanmodel(meanmodel = as.data.frame(matrix(0, 1, 3)), meanparms = 0, locs = matrix(0, 3, 2), n = 3, p = 2))
  expect_equal(.checkargs_meanmodel(meanmodel = function(locs, meanparms) rep(0, nrow(locs)), meanparms = 0, locs = matrix(0, 3, 2), n = 3, p = 2),
               .checkargs_meanmodel(meanmodel = 0, meanparms = 0, locs = matrix(0, 3, 2), n = 3, p = 2))
  expect_identical(tryCatch(.checkargs_meanmodel(meanmodel = matrix(0, 2, 1), meanparms = 0, locs = matrix(0, 3, 2), n = 3, p = 2), error = function(e) TRUE), TRUE)
  expect_identical(tryCatch(.checkargs_meanmodel(meanmodel = matrix(0, 1, 2), meanparms = 0, locs = matrix(0, 3, 2), n = 3, p = 2), error = function(e) TRUE), TRUE)
  expect_identical(tryCatch(.checkargs_meanmodel(meanmodel = matrix(0, 2, 3), meanparms = 0, locs = matrix(0, 3, 2), n = 3, p = 2), error = function(e) TRUE), TRUE)
  expect_identical(tryCatch(.checkargs_meanmodel(meanmodel = 'a', meanparms = 0, locs = matrix(0, 3, 2), n = 3, p = 2), error = function(e) TRUE), TRUE)
})

test_that("covmodel", {
  expect_equal(.checkargs_covmodel(covmodel = function(locs, covparms) diag(1, nrow(locs)), covparms = 1, locs = matrix(runif(6), 3, 2), n = 3, p = 2),
               .checkargs_covmodel(covmodel = diag(1, 3), covparms = 1, locs = matrix(runif(6), 3, 2), n = 3, p = 2))
  expect_equal(.checkargs_covmodel(covmodel = function(locs, covparms) diag(1, nrow(locs)), covparms = 1, locs = matrix(runif(6), 3, 2), n = 3, p = 2),
               .checkargs_covmodel(covmodel = as.data.frame(diag(1, 3)), covparms = 1, locs = matrix(runif(6), 3, 2), n = 3, p = 2))
  expect_equal(.checkargs_covmodel(covmodel = function(locs, covparms) diag(1, nrow(locs)), covparms = 1, locs = matrix(runif(6), 3, 2), n = 3, p = 2),
               .checkargs_covmodel(covmodel = as.data.frame(diag(1, 3)), covparms = 1, locs = matrix(runif(6), 3, 2), n = 3, p = 2))
  expect_equal(.checkargs_covmodel(covmodel = function(locs, covparms) diag(1, nrow(locs)), covparms = 1, locs = matrix(runif(6), 3, 2), n = 3, p = 2),
               .checkargs_covmodel(covmodel = c(1, 1, 1), covparms = 1, locs = matrix(runif(6), 3, 2), n = 3, p = 2))
  expect_equal(.checkargs_covmodel(covmodel = function(locs, covparms) diag(1, nrow(locs)), covparms = 1, locs = matrix(runif(6), 3, 2), n = 3, p = 2),
               .checkargs_covmodel(covmodel = 1, covparms = 1, locs = matrix(runif(6), 3, 2), n = 3, p = 2))
  expect_equal(.checkargs_covmodel(covmodel = function(locs, covparms) diag(1, nrow(locs)), covparms = 1, locs = matrix(runif(6), 3, 2), n = 3, p = 2),
               .checkargs_covmodel(covmodel = "identity", covparms = 1, locs = matrix(runif(6), 3, 2), n = 3, p = 2))
  expect_identical(tryCatch(.checkargs_covmodel(covmodel = matrix(0, 3, 2), covparms = 1, locs = matrix(runif(6), 3, 2), n = 3, p = 2), error = function(e) TRUE), TRUE)
  expect_identical(tryCatch(.checkargs_covmodel(covmodel = matrix(0, 2, 3), covparms = 1, locs = matrix(runif(6), 3, 2), n = 3, p = 2), error = function(e) TRUE), TRUE)
  expect_identical(tryCatch(.checkargs_covmodel(covmodel = "a", covparms = 1, locs = matrix(runif(6), 3, 2), n = 3, p = 2), error = function(e) TRUE), TRUE)
})

test_that("availability of arguments", {
  expect_equal(simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = function(loc, meanparms) 0, meanparms = NULL, covmodel = function(loc1, loc2, covparms) ifelse(identical(loc1, loc2), 1, 0), covparms = NULL, seed = 1)$y,
               simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = function(loc, meanparms) 0, meanparms = NULL, covmodel = function(loc1, loc2, covparms) ifelse(identical(loc1, loc2), 1, 0), covparms = NULL, seed = 1)$y)
  expect_equal(simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = function(loc, meanparms) 0, meanparms = NULL, covmodel = function(loc1, loc2, covparms) ifelse(identical(loc1, loc2), 1, 0), covparms = NULL, seed = 1)$y,
               simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = matrix(0, 3, 1), meanparms = NULL, covmodel = diag(3), covparms = NULL, seed = 1)$y)
  expect_equal(simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = function(loc, meanparms) 0, meanparms = NULL, covmodel = function(loc1, loc2, covparms) ifelse(identical(loc1, loc2), 1, 0), covparms = NULL, seed = 1)$y,
               simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = matrix(0, 1, 3), meanparms = NULL, covmodel = diag(3), covparms = NULL, seed = 1)$y)
  expect_equal(simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = function(loc, meanparms) 0, meanparms = NULL, covmodel = function(loc1, loc2, covparms) ifelse(identical(loc1, loc2), 1, 0), covparms = NULL, seed = 1)$y,
               simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = as.data.frame(matrix(0, 3, 1)), meanparms = NULL, covmodel = as.data.frame(diag(3)), covparms = NULL, seed = 1)$y)
  expect_equal(simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = function(loc, meanparms) 0, meanparms = NULL, covmodel = function(loc1, loc2, covparms) ifelse(identical(loc1, loc2), 1, 0), covparms = NULL, seed = 1)$y,
               simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = as.data.frame(matrix(0, 1, 3)), meanparms = NULL, covmodel = as.data.frame(diag(3)), covparms = NULL, seed = 1)$y)
  expect_equal(simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = function(loc, meanparms) 0, meanparms = NULL, covmodel = function(loc1, loc2, covparms) ifelse(identical(loc1, loc2), 1, 0), covparms = NULL, seed = 1)$y,
               simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = matrix(0, 3, 1), meanparms = NULL, covmodel = rep(1, 3), covparms = NULL, seed = 1)$y)
  expect_equal(simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = function(loc, meanparms) 0, meanparms = NULL, covmodel = function(loc1, loc2, covparms) ifelse(identical(loc1, loc2), 1, 0), covparms = NULL, seed = 1)$y,
               simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = matrix(0, 1, 3), meanparms = NULL, covmodel = 'identity', covparms = NULL, seed = 1)$y)
})

