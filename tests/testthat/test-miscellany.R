####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several tests for small functions (only kldiv so far).
###
####################################################################################


test_that("Kullback Leiber (KL) divergence", {
  expect_equal(kldiv(covmat0 = diag(2), covmat1 = diag(2)), 0)

  expect_equal(kldiv(covmat0 = diag(2), covmat1 = diag(2), mu1 = c(1, 1)), 1)

  expect_equal(kldiv(covmat0 = diag(2), covmat1 = diag(2, 2)), 0.1931, tol = 1e-4)
})


test_that("Kullback Leiber (KL) divergence - example", {
  expect_equal(kldiv(covmat0 = diag(3), covmat1 = matrix(c(1, 0.9, 0.7, 0.9, 1, 0.4, 0.7, 0.4, 1), 3, 3)), 14.4382, tol = 1e-4)
})
