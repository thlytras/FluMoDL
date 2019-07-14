#' Get or set BLUP coefficients for a FluMoDL object
#'
#' This retrieves or sets the BLUP coefficients for a particular
#' \link[=fitFluMoDL]{FluMoDL} object.
#'
#' @param m An object of class \code{\link[=fitFluMoDL]{FluMoDL}}
#'
#' @return For \code{blup.FluMoDL}, the returned object of class
#'   \code{\link{summary.FluMoDL}} holding the BLUP coefficients associated
#'   with the FluMoDL object.
#'
#' @importFrom mvmeta blup
#'
#' @export
blup.FluMoDL <- function(m) {
  return(m$blup)
}


#' @export
`blup<-` <- function(m, value) {
  UseMethod("blup<-")
}


#' @rdname blup.FluMoDL
#' @export
`blup<-.FluMoDL` <- function(m, value) {
  if (!inherits(value, "summary.FluMoDL")) stop("argument should be of class 'summary.FluMoDL'")
  if (!hasRSV(m) && hasRSV(value)) {
    value$coef$proxyRSV <- NULL
    value$vcov$proxyRSV <- NULL
  }
  m$blup <- addPredictions(value, m)
  return(m)
}

