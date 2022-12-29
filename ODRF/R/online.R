#' online structure learning for class \code{ODT} and \code{ODRF}.
#' 
#' The \code{\link{ODT}} and \code{\link{ODRF}} are constantly updated by multiple batches of data to optimize the model. \code{online} is a S3 method for class \code{ODT} and \code{ODRF}.
#'
#' @param obj an object of class \code{ODT} or \code{ODRF}.
#' @param ... For other parameters related to class \code{obj}, see \code{ODT} or \code{ODRF}.
#' 
#' @return object of class \code{ODT} or \code{ODRF}. 
#'
#' @seealso \code{\link{ODT}} \code{\link{ODRF}} \code{\link{online.ODT}} \code{\link{online.ODRF}}
#' 
#' @export
online <- function(obj, ...)
  UseMethod("online",obj) 