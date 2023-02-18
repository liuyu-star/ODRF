data(seeds, package = "ODRF")
train <- sample(1:209, 100)
train_data <- data.frame(seeds[train, ])
test_data <- data.frame(seeds[-train, ])

test_that("predict returns good prediction", {
  tree <- ODT(varieties_of_wheat ~ ., train_data, split = "entropy")
  pred <- predict(tree, test_data[, -8])
  expect_gt(mean(pred == test_data[, 8]), 0.8)
})

test_that("Error if the dimensions of Xnew are not matched with the dimensions of the training data", {
  tree <- ODT(varieties_of_wheat ~ ., train_data, split = "gini")
  expect_error(predict(tree, test_data))
})

test_that("No error if category predictor exists", {
  Xcol1 <- sample(c("A", "B", "C"), 100, replace = TRUE)
  Xcol2 <- sample(c("1", "2", "3", "4", "5"), 100, replace = TRUE)
  Xcon <- matrix(rnorm(100 * 3), 100, 3)
  X <- data.frame(Xcol1, Xcol2, Xcon)
  Xcat <- c(1, 2)
  catLabel <- NULL
  y <- as.factor(sample(c(0, 1), 100, replace = TRUE))
  #options(warn = -1)
  #tree <- ODT(X,y, split = "entropy", Xcat = Xcat,catLabel=catLabel)
  #options(warn = 1)
  expect_warning(ODT(X,y, split = "entropy", Xcat = Xcat,catLabel=catLabel))
  #expect_silent(predict(tree, X))
  #expect_no_error(predict(tree, X))
  #expect_no_error
})


tree <- ODT(Species ~ ., iris, split = "entropy")

test_that("no error if data.frame has Missing value, NA or NaN", {
  dat <- iris[, -5]
  dat[4, 4] <- NA
  dat[25, 1] <- NA
  expect_warning(predict(tree, dat))
})


## Special tests for tree for classification
test_that("predict works for single observations, classification", {
  pred <- predict(tree, head(iris[, -5], 1))
  expect_equal(pred, as.character(iris[1, "Species"]))
})

test_that("Terminal nodes returned by predict are node ids, classification", {
  pred <- predict(tree, iris[, -5], leafnode = TRUE)
  expect_type(pred, "integer")
})

test_that("Terminal nodes returned by predict are node ids, regression", {
  data(body_fat, package = "ODRF")
  tree <- ODT(Density ~ ., body_fat, split = "mse")
  pred <- predict(tree, body_fat[, -1], leafnode = TRUE)
  expect_type(pred, "integer")
})
