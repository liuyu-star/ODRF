#' Pruning of class \code{ODT} or \code{ODRF}.
#' 
#' Prune \code{ODT} or \code{ODRF} from bottom to top with validation data based on prediction error. \code{prune} is a S3 method for class class \code{ODT} and \code{ODRF}.
#'
#' @param obj an object of class \code{ODT} or \code{ODRF}.
#' @param ... For other parameters related to class \code{obj}, see \code{\link{ODT}} or \code{\link{ODRF}}.
#' 
#' @return an object of class \code{ODT} and \code{prune.ODT}.
#'
#' @seealso \code{ODT} \code{\link{ODRF}} \code{prune.ODT} \code{\link{prune.ODRF}}
#' 
#' @export
prune <- function(obj, ...)
  UseMethod("prune",obj)
