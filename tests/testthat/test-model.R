####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several tests for functions regarding baseline models.
###
####################################################################################


test_that("coordinate-based ordering", {

  expect_equal(GPsim:::.order_coordinate(matrix(c(0.25, 0.75, 0.5, 0.51, 0.75, 0.25), 3, 2, byrow = T), coordinate = c(1)), c(1, 2, 3))

  expect_equal(GPsim:::.order_coordinate(matrix(c(0.25, 0.75, 0.5, 0.51, 0.75, 0.25), 3, 2, byrow = T), coordinate = c(2)), c(3, 2, 1))

  expect_equal(GPsim:::.order_coordinate(matrix(c(0.25, 0.75, 0.5, 0.51, 0.75, 0.25), 3, 2, byrow = T), coordinate = c(1, 2)), c(1, 3, 2))

  expect_equal(GPsim:::.order_coordinate(matrix(c(0.25, 0.75, 0.5, 0.51, 0.75, 0.25), 3, 2, byrow = T), coordinate = NULL), c(1, 2, 3))
})


test_that("Maxmin ordering with Euclidean distance", { set.seed(2019)

  expect_equal(GPsim:::.order_maxmin_euclidean(matrix(round(runif(10), 3), 5, 2)), c(4, 5, 2, 3, 1))

  expect_equal(GPsim:::.order_maxmin_euclidean(matrix(round(runif(10), 3), 5, 2)), c(5, 3, 2, 4, 1))

  expect_equal(GPsim:::.order_maxmin_euclidean(matrix(round(runif(10), 3), 5, 2)), c(2, 1, 3, 5, 4))
})


test_that("coordinate-based ordering", { set.seed(2019)

  locs <- matrix(round(runif(10), 3), 5, 2)
  dinv <- GPsim:::.correlation(locs = locs, covmodel = cov_expo_aniso, covparms = c(1, 0.1, 10), def.dist = NULL)
  expect_equal(GPsim:::.order_maxmin_correlation(locs = locs, dinv = dinv, covmodel = cov_expo_aniso, covparms = c(1, 0.1, 10), initial.pt = NULL), c(2, 5, 3, 4, 1))

  locs <- matrix(round(runif(10), 3), 5, 2)
  dinv <- GPsim:::.correlation(locs = locs, covmodel = cov_expo_aniso, covparms = c(1, 0.1, 10), def.dist = NULL)
  expect_equal(GPsim:::.order_maxmin_correlation(locs = locs, dinv = dinv, covmodel = cov_expo_aniso, covparms = c(1, 0.1, 10), initial.pt = NULL), c(3, 1, 2, 5, 4))


  locs <- matrix(round(runif(10), 3), 5, 2)
  dinv <- GPsim:::.correlation(locs = locs, covmodel = cov_expo_aniso, covparms = c(1, 0.1, 10), def.dist = NULL)
  expect_equal(GPsim:::.order_maxmin_correlation(locs = locs, dinv = dinv, covmodel = cov_expo_aniso, covparms = c(1, 0.1, 10), initial.pt = NULL), c(4, 3, 2, 1, 5))
})


test_that("Nearest neighbor conditioning sets with an arbitrary distance matrix", {

  expect_equal(as.numeric(GPsim:::.conditioning_nn(m = 2, dist.matrix = as.matrix(dist(matrix(c(0.25, 0.75, 0.5, 0.51, 0.75, 0.25), 3, 2, byrow = T))))), c(1, 2, 3, NA, 1, 2, NA, NA, 1))
})
