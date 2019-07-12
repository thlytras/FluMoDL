#' Summary method for FluMoDL objects
#'
#' This function creates a summarized version of a 'FluMoDL' object. It contains
#' the sets of coefficients and variance-covariance matrices for the incidence
#' proxy terms (for influenza, and for RSV if provided), and the predictions for these terms.
#'
#' @param m An object of class 'FluMoDL'
#'
#' @return An object of class 'summary.FluMoDL'. This is a list containing the following elements:
#'   \describe{
#'     \item{$type}{A string describing the meaning of the coefficients. Defaults to
#'     "summary", meaning a first-stage model summary. Alternatively, "BLUP" means
#'     Best Unbiased Linear Predictor coefficients, and "pooled" refers to coefficients
#'     pooled in the course of a multivariate meta-analysis. See ...}
#'
#'     \item{$coef}{A list of numeric vectors, with names 'proxyH1', 'proxyH3' and 'proxyB'
#'     (and 'proxyRSV' if provided in the function arguments), containing the model
#'     coefficients for these terms.}
#'
#'     \item{$vcov}{A list of variance-covariance matrices, with names 'proxyH1', 'proxyH3'
#'     and 'proxyB' (and 'proxyRSV' if provided in the function arguments), for the respective
#'     model coefficients.}
#'
#'     \item{$pred}{A list with names 'proxyH1', 'proxyH3' and 'proxyB' (and 'proxyRSV'
#'     if provided in the function arguments), containing
#'     predictions (in the form of \code{\link[dlnm]{crosspred}} objects) for each exposure.
#'     These can be plotted in both the exposure-response and lag-response dimensions, see
#'     \code{\link[dlnm]{crosspred}}, \code{\link[dlnm]{plot.crosspred}} and the example below.}
#'   }
#'
#' @details These summaries can be used to run a multivariate meta-analysis and calculate
#' pooled effect estimates and BLUP (Best Unbiased Linear Predictor) estimates
#' for influenza (and RSV if provided).
#'
#' @examples
#' data(greece) # Use example surveillance data from Greece
#' m <- with(greece, fitFluMoDL(deaths = daily$deaths,
#'     temp = daily$temp, dates = daily$date,
#'     proxyH1 = weekly$ILI * weekly$ppH1,
#'     proxyH3 = weekly$ILI * weekly$ppH3,
#'     proxyB = weekly$ILI * weekly$ppB,
#'     yearweek = weekly$yearweek))
#' summ <- summary(m)
#' summ
#'
#' # Plot the association between A(H1N1)pdm09 activity and mortality:
#' plot(summ$pred$proxyH1, "overall")
#'
#' @export
summary.FluMoDL <- function(m) {
  if (!("FluMoDL" %in% class(m))) stop("Argument `m` should be of class 'FluMoDL'.")
  if (is.null(m$pred)) stop("No 'pred' element found; object is corrupted.")
  res <- list(
    type = "summary",
    coef = lapply(m$pred[names(m$pred)[grep("proxy", names(m$pred))]], coef),
    vcov = lapply(m$pred[names(m$pred)[grep("proxy", names(m$pred))]], vcov),
    pred = m$pred[grep("proxy", names(m$pred))])
  class(res) <- "summary.FluMoDL"
  return(res)
}



#' @export
print.summary.FluMoDL <- function(s) {
  cat("\n** FluMoDL model summary **\n\n")
  if (!is.null(s$pred)) {
    cat("Object includes predictions objects (of class 'crosspred').\n")
    mid <- ceiling(length(s$pred$proxyH1$allfit)/2)
    cat(sprintf("Relative Risk for an indicative influenza incidence proxy of %s: (95%% CI)\n",
                names(s$pred$proxyH1$allfit)[mid]))
    cat(sprintf("Influenza A(H1N1)pdm09 = %.3f (%.3f - %.3f)\n",
                exp(s$pred$proxyH1$allfit[mid]),
                s$pred$proxyH1$allRRlow[mid], s$pred$proxyH1$allRRhigh[mid]))
    cat(sprintf("Influenza A(H3N2)      = %.3f (%.3f - %.3f)\n",
                exp(s$pred$proxyH3$allfit[mid]),
                s$pred$proxyH3$allRRlow[mid], s$pred$proxyH3$allRRhigh[mid]))
    cat(sprintf("Influenza B            = %.3f (%.3f - %.3f)\n",
                exp(s$pred$proxyB$allfit[mid]),
                s$pred$proxyB$allRRlow[mid], s$pred$proxyB$allRRhigh[mid]))
    if (!is.null(s$pred$proxyRSV)) {
      cat("Object contains a term for RSV (Respiratory Syncytial Virus)\n")
      cat(sprintf("RR for an indicative RSV incidence proxy of %s = %.3f (%.3f - %.3f)\n",
                  names(s$pred$proxyRSV$allfit)[mid],
                  exp(s$pred$proxyRSV$allfit[mid]),
                  s$pred$proxyRSV$allRRlow[mid],
                  s$pred$proxyRSV$allRRhigh[mid]))
    }

  } else {
    cat("Object does NOT include 'crosspred' objects.\n")
    if (!is.null(s$pred$proxyRSV)) {
      cat("Object contains a term for RSV (Respiratory Syncytial Virus)\n")
    }
    cat("Summary of coefficients per incidence proxy:")
    print(do.call("cbind", s$coef))
  }
}
