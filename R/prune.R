#' prune \code{ODT} or \code{ODRF}
#'
#' Prune \code{ODT} or \code{ODRF} from bottom to top with validation data based on prediction error, and \code{prune} is a S3 method for class \code{ODT} and \code{ODRF}.
#'
#' @param obj An object of class \code{ODT} or \code{ODRF}.
#' @param ... For other parameters related to class \code{obj}, see \code{\link{ODT}} or \code{\link{ODRF}}.
#'
#' @return An object of class \code{ODT} and \code{prune.ODT}.
#'
#' @keywords internal
#'
#' @seealso \code{\link{ODT}} \code{\link{ODRF}} \code{\link{prune.ODT}} \code{\link{prune.ODRF}}
#'
#' @export
prune <- function(obj, ...) {
  UseMethod("prune", obj)
}
