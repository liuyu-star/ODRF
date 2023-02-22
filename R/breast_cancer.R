#' Breast Cancer Dataset
#'
#' Breast cancer is the most common cancer amongst women in the world. It accounts for \eqn{25\%} of all cancer cases, and affected over 2.1 Million people in 2015 alone.
#' It starts when cells in the breast begin to grow out of control. These cells usually form tumors that can be seen via X-ray or felt as lumps in the breast area.
#' The key challenges against it's detection is how to classify tumors into malignant (cancerous) or benign(non cancerous).
#'
#' The actual linear program used to obtain the separating plane in the 3-dimensional space is that described in:
#' \itemize{
#' \item ID number
#' \item Diagnosis (M = malignant, B = benign)
#' \item Ten real-valued features are computed for each cell nucleus:
#' \itemize{
#' \item radius (mean of distances from center to points on the perimeter)
#' \item texture (standard deviation of gray-scale values)
#' \item perimeter
#' \item area
#' \item smoothness (local variation in radius lengths)
#' \item compactness (perimeter^2 / area - 1.0)
#' \item concavity (severity of concave portions of the contour)
#' \item concave points (number of concave portions of the contour)
#' \item symmetry
#' \item fractal dimension ("coastline approximation" - 1)
#' }}
#'
#' @docType data
#' @keywords datasets internal
#' @format A data frame with 569 rows and 30 covariate variables and 1 response variable
#' @source \url{https://www.kaggle.com/datasets/yasserh/breast-cancer-dataset?select=breast-cancer.csv} and \url{https://archive.ics.uci.edu/ml/datasets/breast+cancer+wisconsin+(diagnostic)}
#' @references Wolberg WH, Street WN, Mangasarian OL. Machine learning techniques to diagnose breast cancer from image-processed nuclear features of fine needle aspirates. Cancer Lett. 1994 Mar 15;77(2-3):163-71.
#' @name breast_cancer
#'
#' @seealso  \code{\link{body_fat}} \code{\link{seeds}}
#' @examples
#' data(breast_cancer)
#' set.seed(221212)
#' train <- sample(1:569, 80)
#' train_data <- data.frame(breast_cancer[train, -1])
#' test_data <- data.frame(breast_cancer[-train, -1])
#'
#' forest <- ODRF(diagnosis ~ ., train_data, split = "gini", parallel = FALSE, ntrees = 50)
#' pred <- predict(forest, test_data[, -1])
#' # classification error
#' (mean(pred != test_data[, 1]))
#'
#' tree <- ODT(diagnosis ~ ., train_data, split = "gini")
#' pred <- predict(tree, test_data[, -1])
#' # classification error
#' (mean(pred != test_data[, 1]))
NULL
