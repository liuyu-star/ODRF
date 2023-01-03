#' @export
ODT <- function(X, ...) {
  UseMethod("ODT")
  #formula X
}

#' @export
"ODRF" <- function(X, ...) {
  UseMethod("ODRF")
}

#' online structure learning for class \code{ODT} and \code{ODRF}.
#'
#' \code{\link{ODT}} and \code{\link{ODRF}} are constantly updated by multiple batches of data to optimize the model. \code{online} is a S3 method for class \code{ODT} and \code{ODRF}.
#'
#' @param obj an object of class \code{ODT} or \code{ODRF}.
#' @param ... For other parameters related to class \code{obj}, see \code{ODT} or \code{ODRF}.
#'
#' @return object of class \code{ODT} or \code{ODRF}.
#'
#' @seealso \code{\link{ODT}} \code{\link{ODRF}} \code{\link{online.ODT}} \code{\link{online.ODRF}}
#'
#' @export
online <- function(obj, ...) {
  UseMethod("online", obj)
}

#' prune \code{ODT} or \code{ODRF}
#'
#' Prune \code{ODT} or \code{ODRF} from bottom to top with validation data based on prediction error, and \code{prune} is a S3 method for class \code{ODT} and \code{ODRF}.
#'
#' @param obj An object of class \code{ODT} or \code{ODRF}.
#' @param ... For other parameters related to class \code{obj}, see \code{\link{ODT}} or \code{\link{ODRF}}.
#'
#' @return An object of class \code{ODT} and \code{prune.ODT}.
#'
#' @seealso \code{\link{ODT}} \code{\link{ODRF}} \code{\link{prune.ODT}} \code{\link{prune.ODRF}}
#'
#' @export
prune <- function(obj, ...) {
  UseMethod("prune", obj)
}
