% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RotMat.R
\name{RotMatMake}
\alias{RotMatMake}
\title{Create rotation matrix used to determine linear combination of features.}
\usage{
RotMatMake(
  X = NULL,
  y = NULL,
  RotMatFun = "RotMatPPO",
  PPFun = "PPO",
  FunDir = getwd(),
  paramList = NULL,
  ...
)
}
\arguments{
\item{X}{An n by d numeric matrix (preferable) or data frame.}

\item{y}{A response vector of length n.}

\item{RotMatFun}{A self-defined projection matrix function name, which can also be \code{\link{RotMatRand}} and \code{\link{RotMatPPO}}. Note that \code{(,...)} is necessary.}

\item{PPFun}{A self-defined projection matrix function, which can also be \code{\link{PPO}}. Note that \code{(,...)} is necessary.}

\item{FunDir}{The path to the \code{function} of the user-defined \code{NodeRotateFun}. (default current Workspace)}

\item{paramList}{List of parameters used by the functions \code{RotMatFun} and \code{PPFun}. If left unchanged, default values will be used, for details see \code{\link[ODRF]{defaults}}.}

\item{...}{Used to handle superfluous arguments passed in using paramList.}
}
\value{
A random matrix to use in running \code{\link{ODT}}.
\itemize{
\item{Variable: Variables to be projected.}
\item{Number: Number of projections.}
\item{Coefficient: Coefficients of the projection matrix.}
}
}
\description{
Create any projection matrix with a self-defined projection matrix function and projection optimization model function
}
\details{
There are two ways for the user to define a projection direction function. The first way is to define a function directly, and just let the argument \code{RotMatFun} be the name of the defined function
and let the argument \code{paramList} be the arguments used in the defined function; the second way is to use the function \code{RotMatMake} following the way in the example below.
Note that the name of the defined function cannot be the name of an existing function in the ODRF package
}
\examples{
set.seed(220828)
X <- matrix(rnorm(1000), 100, 10)
y <- (rnorm(100) > 0) + 0
(RotMat <- RotMatMake(X, y, "RotMatRand", "PPO"))
library(nnet)
(RotMat <- RotMatMake(X, y, "RotMatPPO", "PPO", paramList = list(model = "Log")))

## Define projection matrix function and projection optimization model function.##
## Note that (,...) is necessary.
makeRotMat <- function(dimX, dimProj, numProj, ...) {
  RotMat <- matrix(1, dimProj * numProj, 3)
  for (np in seq(numProj)) {
    RotMat[(dimProj * (np - 1) + 1):(dimProj * np), 1] <-
      sample(1:dimX, dimProj, replace = FALSE)
    RotMat[(dimProj * (np - 1) + 1):(dimProj * np), 2] <- np
  }
  return(RotMat)
}

makePP <- function(dimProj, prob, ...) {
  pp <- sample(c(1L, -1L), dimProj, replace = TRUE, prob = c(prob, 1 - prob))
  return(pp)
}

RotMat <- RotMatMake(
  RotMatFun = "makeRotMat", PPFun = "makePP",
  paramList = list(dimX = 8, dimProj = 5, numProj = 4, prob = 0.5)
)
head(RotMat)
#>      Variable Number Coefficient
#> [1,]        6      1           1
#> [2,]        8      1           1
#> [3,]        1      1          -1
#> [4,]        4      1          -1
#> [5,]        5      1          -1
#> [6,]        6      2           1

# train ODT with defined projection matrix function
tree <- ODT(X, y,
  type = "i-classification", NodeRotateFun = "makeRotMat",
  paramList = list(dimX = ncol(X), dimProj = 5, numProj = 4)
)
# train ODT with defined projection matrix function and projection optimization model function
tree <- ODT(X, y,
  type = "i-classification", NodeRotateFun = "RotMatMake", paramList =
    list(
      RotMatFun = "makeRotMat", PPFun = "makePP",
      dimX = ncol(X), dimProj = 5, numProj = 4, prob = 0.5
    )
)

}
\seealso{
\code{\link{RotMatPPO}} \code{\link{RotMatRand}} \code{\link{RotMatRF}}
}
\keyword{rotation}