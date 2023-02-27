
## Tests for predictions
## Tests
test_that("predict returns good prediction", {
  data(seeds, package = "ODRF")
  set.seed(230227)
  train <- sample(1:209, 100)
  train_data <- data.frame(seeds[train, ])
  test_data <- data.frame(seeds[-train, ])
  forest <- ODRF(varieties_of_wheat ~ ., train_data,
    split = "entropy",
    parallel = FALSE
  )
  pred <- predict(forest, test_data[, -8])
  expect_gt(mean(pred == test_data[, 8]), 0.9)
})


# dat <- data.matrix(iris)
forest1 <- ODRF(Species ~ ., data = iris, split = "entropy", parallel = FALSE)

test_that("no error if data.frame has Missing value, NA or NaN", {
  dat <- iris[, -5]
  dat[4, 4] <- NA
  dat[25, 1] <- NA
  expect_warning(predict(forest1, dat))
})

test_that("If type = 'tree', these number is used for predictions", {
  pred <- predict(forest1, iris[, -5], type = "tree")
  expect_equal(dim(pred), c(nrow(iris), 100))
})

test_that("If type = 'prob', these number is used for predictions", {
  pred <- predict(forest1, iris[, -5], type = "prob")
  expect_equal(dim(pred), c(nrow(iris), 3))
  expect_equal(rowSums(pred), rep(1, nrow(iris)))
})


test_that("Error if unknown value for type", {
  expect_error(predict(forest1, iris[, -5], type = "class"))
})

test_that("Tree-wise split select weights work", {
  expect_silent(predict(forest1, iris[, -5], weight.tree = TRUE))
})


# Inbag count matrix
test_that("Inbag count matrix if of right size, without replacement", {
  rf <- ODRF(Species ~ ., iris, split = "entropy", ntrees = 10, replacement = FALSE, parallel = FALSE)
  pred <- predict(rf, iris[, -5], type = "tree")
  expect_equal(dim(pred), c(nrow(iris), 10))
})

test_that("Inbag count matrix if of right size, storeOOB = FALSE and weight.tree=TRUE", {
  rf <- ODRF(Species ~ ., iris, split = "entropy", ntrees = 10, storeOOB = FALSE, parallel = FALSE)
  expect_error(predict(rf, iris[, -5], type = "prob", weight.tree = TRUE))
})

test_that("Inbag count matrix if of right size, numOOB = 0 and weight.tree=TRUE", {
  rf <- ODRF(Species ~ ., iris, split = "entropy", ntrees = 10, numOOB = 0, parallel = FALSE)
  expect_warning(predict(rf, iris[, -5], type = "prob", weight.tree = TRUE))
})
