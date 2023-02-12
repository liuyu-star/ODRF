## Tests for random forests for classification
## Initialize the random forest for classification
dat <- data.matrix(iris)
forest <- ODRF(Species ~ ., data = iris, type = "entropy", parallel = FALSE)

## Basic tests (for all random forests equal)
test_that("classification iris forest is of class ODRF with 12 elements", {
  expect_s3_class(forest, "ODRF")
  expect_length(forest, 12)
})

test_that("Error if data is class matrix, classification", {
  expect_error(ODRF(Species ~ ., data = dat, type = "entropy", parallel = FALSE))
})

test_that("Error if ntrees=1, classification", {
  n <- 100
  X <- matrix(runif(5 * n), n, 5)
  y <- as.factor(rbinom(n, 1, 0.5))
  expect_error(ODRF(X = X, y = y, ntrees = 1, parallel = FALSE))
})

test_that("Error if y is a factor type and type = 'regression'", {
  n <- 100
  X <- matrix(runif(5 * n), n, 5)
  y <- as.factor(rbinom(n, 1, 0.5))
  expect_error(ODRF(X = X, y = y, type = "mse", parallel = FALSE))
})


test_that("confusion matrix is of right dimension", {
  expect_equal(dim(forest$oobConfusionMat), c(nlevels(iris$Species), nlevels(iris$Species) + 1))
})

test_that("confusion matrix has right dimnames", {
  # predicted =  ,#true =
  expect_equal(dimnames(forest$oobConfusionMat), list(
    levels(iris$Species),
    c(levels(iris$Species), "class_error")
  ))
})

test_that("confusion matrix cols are the true classes", {
  expect_equal(
    as.numeric(colSums(forest$oobConfusionMat[, seq(nlevels(iris$Species)), drop = FALSE])),
    as.numeric(table(iris$Species))
  )
})

## Tests for using seeds
## Initialize the random forests

#ind <- 1:150 %in% sample(150, 100)
ind <- sample(150, 100)

forest1 <- ODRF(Species ~ .,
  data = iris, type = "gini", NodeRotateFun = "RotMatPPO",
  parallel = FALSE, seed = 1, ntrees = 10#,subset = ind
)
pred1 <- predict(forest1, Xnew = iris[-ind, -5])

forest2 <- ODRF(Species ~ .,
  data = iris, type = "gini", NodeRotateFun = "RotMatPPO",
  parallel = FALSE, seed = 1, ntrees = 10#,subset = ind
)
pred2 <- predict(forest2, Xnew = iris[-ind, -5])

forest3 <- ODRF(Species ~ .,
  data = iris, type = "gini",NodeRotateFun = "RotMatPPO",
  parallel = FALSE, seed = 2, ntrees = 10#,subset = ind
)
pred3 <- predict(forest3, Xnew = iris[-ind, -5])

## Tests
test_that("same result with same seed", {
  expect_equal(pred1, pred2)
})

#test_that("different result with different seed", {
#  expect_false(identical(pred1, pred3))
#})
