## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----ODT----------------------------------------------------------------------
library(ODRF)
data(seeds, package = "ODRF")
set.seed(19)
train <- sample(1:209, 150)
train_data <- data.frame(seeds[train, ])
test_data <- data.frame(seeds[-train, ])
index <- seq(floor(1*nrow(train_data) / 2))

forest <- ODRF(varieties_of_wheat ~ ., train_data,
  split = "gini", parallel = FALSE
)
pred <- predict(forest, test_data[, -8])
e.forest <- mean(pred != test_data[, 8])
forest1 <- ODRF(varieties_of_wheat ~ ., train_data[index, ],
  split = "gini", parallel = FALSE
)
pred <- predict(forest1, test_data[, -8])
e.forest.1 <- mean(pred != test_data[, 8])
forest2 <- ODRF(varieties_of_wheat ~ ., train_data[-index, ],
  split = "gini", parallel = FALSE
)
pred <- predict(forest2, test_data[, -8])
e.forest.2 <- mean(pred != test_data[, 8])

forest.online <- online(
  forest1, train_data[-index, -8],
  train_data[-index, 8]
)
pred <- predict(forest.online, test_data[, -8])
e.forest.online <- mean(pred != test_data[, 8])
forest.prune <- prune(forest1, train_data[-index, -8],
  train_data[-index, 8],
  useOOB = FALSE
)
pred <- predict(forest.prune, test_data[, -8])
e.forest.prune <- mean(pred != test_data[, 8])
print(c(
  forest = e.forest, forest1 = e.forest.1, forest2 = e.forest.2,
  forest.online = e.forest.online, forest.prune = e.forest.prune
))

## ----ODRF---------------------------------------------------------------------
data(body_fat, package = "ODRF")
set.seed(42)
train <- sample(1:252, 150)
train_data <- data.frame(body_fat[train, ])
test_data <- data.frame(body_fat[-train, ])
index <- seq(floor(1*nrow(train_data) / 2))

tree <- ODT(Density ~ ., train_data, split = "mse")
pred <- predict(tree, test_data[, -1])
e.tree <- mean((pred - test_data[, 1])^2)
tree1 <- ODT(Density ~ ., train_data[index, ], split = "mse")
pred <- predict(tree1, test_data[, -1])
e.tree.1 <- mean((pred - test_data[, 1])^2)
tree2 <- ODT(Density ~ ., train_data[-index, ], split = "mse")
pred <- predict(tree2, test_data[, -1])
e.tree.2 <- mean((pred - test_data[, 1])^2)

tree.online <- online(tree1, train_data[-index, -1], train_data[-index, 1])
pred <- predict(tree.online, test_data[, -1])
e.tree.online <- mean((pred - test_data[, 1])^2)
tree.prune <- prune(tree1, train_data[-index, -1], train_data[-index, 1])
pred <- predict(tree.prune, test_data[, -1])
e.tree.prune <- mean((pred - test_data[, 1])^2)
print(c(
  tree = e.tree, tree1 = e.tree.1, tree2 = e.tree.2,
  tree.online = e.tree.online, tree.prune = e.tree.prune
))

## ----print--------------------------------------------------------------------
data(iris, package = "datasets")
tree <- ODT(Species ~ ., data = iris)
tree
forest <- ODRF(Species ~ ., data = iris, parallel = FALSE)
forest

## ----plot, fig.height=5.5,fig.width=7.2---------------------------------------
plot(tree)

