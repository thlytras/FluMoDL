#' Multivariate meta-analysis for FluMoDL objects
#'
#' This function runs multivariate meta-analysis (using package \code{\link{mvmeta}})
#' on the first-stage coefficients of influenza (and possibly RSV) incidence proxies
#' for multiple 'FluMoDL' object summaries.
#'
#' @param summaries A \emph{list} of objects of class \code{\link{summary.FluMoDL}}
#'   (at least two), representing the first-stage analyses. If the list is named, the
#'   names are kept in the output object and can be retrieved with \code{names()}, see below.
#'
#' @param par For which model terms (sets of coefficients) to run the meta-analysis?
#' Defaults to \code{c("H1","H3","B","RSV")}, which
#' indicates all three influenza proxies and RSV (for those summaries that have included
#' an RSV term). It is unlikely that you'll want to alter this default.
#'
#' @return Returns an object of class 'metaFluMoDL'. This is a list of objects of class
#'   \code{\link[mvmeta]{mvmeta}}, representing the results of the multivariate
#'   random-effects meta-analysis for the sets of coefficients corresponding to each
#'   term in argument \code{par}; they can be accessed directly using the \code{$}
#'   operator as \code{$proxyH1}, \code{$proxyH3} and \code{$proxyB} (and also
#'   \code{$proxyRSV} if there were RSV terms in at least two elements of
#'   \code{summaries} and \code{par} included "RSV" -- in which case,
#'   \code{\link[=hasRSV]{hasRSV()}} returns \code{TRUE} for objects of class 'metaFluMoDL').
#'
#'   \emph{However}, some methods have been redefined for class 'metaFluMoDL', and do not
#'   work the same as in simple lists. In particular: \code{\link[=length]{length()}}
#'   returns the number of summaries (number of "studies") meta-analyzed and
#'   \code{\link[=names]{names()}} returns the names of these summaries (if the list
#'   in \code{summaries} argument was named).
#'
#'   In addition, the \code{[[} and \code{[} operators have been redefined for class
#'   'metaFluMoDL', and now return the Best Linear Unbiased Predictor (BLUP)
#'   estimates for the selected summaries ("studies"), as objects of class
#'   \code{\link{summary.FluMoDL}}; selection can be made the usual way,
#'   with a logical or numeric index vector, or with the summary names
#'   (as provided by \code{names}). \code{[} returns a \emph{list} of
#'   \code{\link{summary.FluMoDL}} objects, whereas \code{[[} returns a single object.
#'   The returned objects contain the string "blup" in their \code{$type} element,
#'   to distinguish them from \code{\link[=summary.FluMoDL]{first-stage model summaries}}
#'   or \code{\link{pooled}}
#'   result summaries. In their \code{$description} element, they contain the
#'   name of the respective summary ("study") if a named list had been provided
#'   in the \code{summaries} argument of \code{metaFluMoDL()}. And finally, they contain
#'   no \code{$pred} element, as they are not associated with a particular dataset and
#'   cross-basis matrices (which is a prerequisite to create
#'   \code{\link[dlnm]{crosspred}} objects).
#'
#'   The pooled coefficients (for all three or four incidence proxies) can be obtained
#'   with function \code{\link[=pooled]{pooled()}}, which also returns an object of class
#'   \code{\link{summary.FluMoDL}} that you can further use.
#'
#' @seealso \code{\link{summary.FluMoDL}}, \code{\link{pooled}}
#'
#' @references \itemize{
#'  \item Gasparrini A, Armstrong B, Kenward MG. Multivariate meta-analysis for non-linear
#'   and other multi-parameter associations.
#'   \href{https://onlinelibrary.wiley.com/doi/full/10.1002/sim.5471}{Stat Med} 2012;31(29):3821–39.
#' }
#'
#' @export
metaFluMoDL <- function(summaries, par=c("H1","H3","B","RSV")) {
  nm <- names(summaries)
  par <- paste0("proxy", par)
  # Get the number of summaries with data for each parameter
    Nsum <- rowSums(sapply(summaries, function(s) par %in% names(s$coef)))
  # If a parameter has <2 summaries, drop it (without a warning; maybe add a warning in the future)
    par <- par[Nsum>=2]
  # Loop through the parameters and create the meta-analysis
  M <- lapply(par, function(p) {
    TE <- lapply(summaries, function(s) s$coef[[p]])
    seTE <- lapply(summaries, function(s) s$vcov[[p]])

    nmm <- nm[!sapply(seTE, is.null)]
    TE <- TE[!sapply(TE, is.null)]
    TE <- do.call("rbind", TE)
    seTE <- seTE[!sapply(seTE, is.null)]

    rownames(TE) <- nmm
    names(seTE) <- nmm

    mvmeta(TE, S=seTE)
  })
  names(M) <- par
  class(M) <- "metaFluMoDL"
  M
}



#' @export
names.metaFluMoDL <- function(m) {
  return(rownames(unclass(m)[[1]]$model))
}

#' @export
`names<-.metaFluMoDL` <- function(m, value) {
  return(m)
}


#' @export
length.metaFluMoDL <- function(m) {
  return(nrow(unclass(m)[[1]]$model))
}

#' @export
`length<-.metaFluMoDL` <- function(m, value) {
  return(m)
}



#' @export
`[[.metaFluMoDL` <- function(m, i) {
  if (length(i)>1) stop("attempt to select more than one element")
  nn <- names(m)
  if (is.character(i) && !(i %in% nn)) stop(sprintf("label '%s' not found in metaFluMoDL object", i))
  if (is.numeric(i)) {
    if (i>length(nn)) stop("index out of bounds")
    if (i<0) stop("index cannot be negative")
    i <- nn[i]
  }
  u <- unclass(m)
  par <- names(u)
  res0 <- lapply(par, function(p) {
    blup(u[[p]], vcov=TRUE)[[i]]
  })
  names(res0) <- par
  res0 <- res0[!sapply(res0, is.null)]
  res <- list(
    type = "blup",
    description = i,
    coef = lapply(res0, function(x) x$blup),
    vcov = lapply(res0, function(x) x$vcov),
    pred = NULL
  )
  class(res) <- "summary.FluMoDL"
  res
}


#' @export
`[.metaFluMoDL` <- function(m, i) {
  nn <- names(m)
  if (is.logical(i)) {
    i <- which(i[1:length(nn)])
  }
  if (is.character(i)) {
    if (sum(!i %in% nn)>0) stop("name(s) not found in object")
    i <- match(i, nn)
  }
  res <- lapply(i, function(ii) m[[ii]])
  names(res) <- nn[i]
  res
}


#' @export
`[[<-.metaFluMoDL` <- function(m, i) {
  return(m)
}


#' @export
`[<-.metaFluMoDL` <- function(m, i) {
  return(m)
}

