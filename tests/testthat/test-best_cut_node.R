## Tests for split select

data(iris)
X <- as.matrix(iris[, 1:4])
y <- iris[[5]]
bestcut <- best.cut.node(X, y, split = "entropy")

test_that("Return three vectors with length 1, 1, and 'X' dimension", {
  expect_equal(unlist(lapply(bestcut, length)), c(BestCutVar = 1, BestCutVal = 1, BestIndex = ncol(X)))
})

test_that("Error if data X has Missing value, NA or NaN", {
  n <- 100
  y <- rbinom(n, 1, 0.5)
  X <- cbind(c(NA, NaN, rnorm(n - 2)), rnorm(n), c(rnorm(n - 3), NA, NaN, NA))
  expect_error(best.cut.node(X, y, split = "entropy"))
})

## weights
test_that("split select weights work", {
  expect_silent(best.cut.node(X, y, split = "entropy", weights = runif(length(y))))
})

test_that("split select weights work", {
  n <- 100
  y <- rnorm(n)
  X <- cbind(runif(n), rnorm(n), rbinom(n, 1, 0.5))
  expect_silent(best.cut.node(X, y, split = "mse", weights = runif(n)))
})

## criteria
test_that("split select 'g-classification' work", {
  expect_silent(best.cut.node(X, y, split = "gini"))
})

test_that("split select 'regression' work", {
  n <- 100
  y <- rnorm(n)
  X <- cbind(runif(n), rnorm(n), rbinom(n, 1, 0.5))
  expect_silent(best.cut.node(X, y, split = "mse"))
})

test_that("no warning if argument 'X' dimension equal to 1", {
  n <- 100
  X <- rnorm(100)
  y <- rbinom(100, 1, .5)
  expect_silent(best.cut.node(X, y, split = "entropy"))
})

test_that("Split points are at (A+B)/2 for numeric features, i-classification splitting", {
  y <- rbinom(100, 1, .5)
  X <- matrix(rbinom(5 * 100, 1, .5), 100, 5)
  best_cut <- best.cut.node(X, y, split = "entropy")
  if (best_cut$BestCutVar == -1) {
    expect_false(FALSE)
  } else {
    expect_equal(best_cut$BestCutVal, 0.5)
  }
})

test_that("Split points are at (A+B)/2 for numeric features, g-classification splitting", {
  y <- rbinom(100, 1, .5)
  X <- matrix(rbinom(5 * 100, 1, .5), 100, 5)
  best_cut <- best.cut.node(X, y, split = "gini")
  if (best_cut$BestCutVar == -1) {
    expect_false(FALSE)
  } else {
    expect_equal(best_cut$BestCutVal, 0.5)
  }
})

test_that("Split points are at (A+B)/2 for numeric features, regression", {
  y <- rnorm(100)
  X <- matrix(rbinom(5 * 100, 1, .5), 100, 5)
  best_cut <- best.cut.node(X, y, split = "mse")
  if (best_cut$BestCutVar == -1) {
    expect_false(FALSE)
  } else {
    expect_equal(best_cut$BestCutVal, 0.5)
  }
})
