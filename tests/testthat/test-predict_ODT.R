
test_that("predict returns good prediction", {
  data(seeds, package = "ODRF")
  train <- sample(1:209, 100)
  train_data <- data.frame(seeds[train, ])
  test_data <- data.frame(seeds[-train, ])

  tree <- ODT(varieties_of_wheat ~ ., train_data, type = "i-classification")
  pred <- predict(tree, test_data[, -8])
  expect_gt(mean(pred == test_data[, 8]), 0.8)
})

tree <- ODT(Species ~ ., iris, type = "i-classification")

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
  tree <- ODT(Density ~ ., body_fat, type = "regression")
  pred <- predict(tree, body_fat[, -1], leafnode = TRUE)
  expect_type(pred, "integer")
})
