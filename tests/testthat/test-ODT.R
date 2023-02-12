
test_that("classification seeds tree is of class ODRF with 10 elements", {
  data(body_fat, package = "ODRF")
  tree <- ODT(Density ~ ., body_fat, type = "mse")
  expect_s3_class(tree, "ODT")
  expect_length(tree, 10)
})

test_that("no warning if data.frame has two classes", {
  dat <- iris
  class(dat) <- c("data.frame", "data.table")
  expect_silent(ODT(Species ~ ., data = iris, type = "entropy"))
})

test_that("Error if subset is 0", {
  expect_error(ODT(Species ~ ., data = iris, type = "gini", subset = c(0)))
})

test_that("Error if argument 'X' dimension equal to 1", {
  n <- 100
  X <- matrix(runif(n))
  y <- rnorm(n)
  expect_error(ODT(X, y, type = "gini"))
})

test_that("Error if the number of factor levels of response variable equal to 1", {
  n <- 100
  X <- matrix(runif(10 * n), n, 10)
  y <- rep(1, n)
  expect_error(ODT(X, y, type = "gini"))
})

test_that("no warning if data.frame has Missing value, NA or NaN", {
  n <- 100
  y <- rnorm(n)
  X <- cbind(c(NA, NaN, rnorm(n - 2)), rnorm(n), c(rnorm(n - 3), NA, NaN, NA))
  dat <- data.frame(y = y, X)
  expect_silent(ODT(y ~ ., data = dat, type = "mse"))
})

test_that("regression splitting not working on classification data", {
  expect_error(ODT(Species ~ ., iris, type = "mse"))
})

test_that("as.factor() in the arguments formula and y automatically works for classification", {
  n <- 100
  X <- matrix(runif(10 * n), n, 10)
  y <- as.factor(rbinom(n, 1, 0.5))
  expect_warning(ODT(X, y))
})

test_that("holdout mode holding out data with 0 weight", {
  weights <- rbinom(nrow(iris), 1, 0.5)
  expect_silent(ODT(Species ~ ., data = iris, type = "entropy", weights = weights))
})


test_that("Split points are at (A+B)/2 for numeric features, i-classification splitting", {
  dat <- data.frame(y = rbinom(100, 1, .5), matrix(rbinom(5 * 100, 1, .5), 100, 5))
  tree <- ODT(y ~ ., data = dat, type = "entropy", NodeRotateFun = "RotMatRF")
  split_points <- tree$structure$nodeCutValue
  nsp <- which(split_points != 0)
  expect_equal(split_points[nsp], rep(0.5, length(nsp)))
})

test_that("Split points are at (A+B)/2 for numeric features, g-classification splitting", {
  dat <- data.frame(y = rbinom(100, 1, .5), matrix(rbinom(5 * 100, 1, .5), 100, 5))
  tree <- ODT(y ~ ., data = dat, type = "gini", NodeRotateFun = "RotMatRF")
  split_points <- tree$structure$nodeCutValue
  nsp <- which(split_points != 0)
  expect_equal(split_points[nsp], rep(0.5, length(nsp)))
})

test_that("Split points are at (A+B)/2 for numeric features, regression", {
  dat <- data.frame(y = rnorm(100, 1, .5), matrix(rbinom(5 * 100, 1, .5), 100, 5))
  tree <- ODT(y ~ ., data = dat, type = "mse", NodeRotateFun = "RotMatRF")
  split_points <- tree$structure$nodeCutValue
  nsp <- which(split_points != 0)
  expect_equal(split_points[nsp], rep(0.5, length(nsp)))
})

## Tests for using seeds
## Initialize the random forests
ind <- 1:150 %in% sample(150, 100)

set.seed(1)
tree1 <- ODT(Species ~ ., data = iris, type = "gini", NodeRotateFun = "RotMatRand")
pred1 <- predict(tree1, Xnew = iris[!ind, -5])

set.seed(1)
tree2 <- ODT(Species ~ ., data = iris, type = "gini", NodeRotateFun = "RotMatRand")
pred2 <- predict(tree2, Xnew = iris[!ind, -5])

set.seed(2)
tree3 <- ODT(Species ~ ., data = iris, type = "gini", NodeRotateFun = "RotMatRand")
pred3 <- predict(tree3, Xnew = iris[!ind, -5])

## Tests
test_that("same result with same seed", {
  expect_equal(pred1, pred2)
})

test_that("different result with different seed", {
  expect_false(identical(pred1, pred3))
})
