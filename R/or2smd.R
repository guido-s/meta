#' Conversion from log odds ratio to standardised mean difference
#' 
#' @description
#' Conversion from log odds ratio to standardised mean difference
#' using method by Hasselblad & Hedges (1995) or Cox (1970).
#' 
#' @param lnOR Log odds ratio(s) or meta-analysis object.
#' @param selnOR Standard error(s) of log odds ratio(s).
#' @param method A character string indicating which method is used to
#'   convert log odds ratios to standardised mean differences. Either
#'   \code{"HH"} or \code{"CS"}, can be abbreviated.
#' @param \dots Additional arguments (not considered if argument
#'   \code{lnOR} is a meta-analysis object).
#' 
#' @details
#' This function implements the following methods for the conversion
#' from log odds ratios to standardised mean difference:
#' \itemize{
#' \item Hasselblad & Hedges (1995) assuming logistic distributions
#'   (\code{method == "HH"})
#' \item Cox (1970) and Cox & Snell (1989) assuming normal
#'   distributions (\code{method == "CS"})
#' }
#' Internally, \code{\link{metagen}} is used to conduct a
#' meta-analysis with the standardised mean difference as summary
#' measure.
#' 
#' Argument \code{lnOR} can be either a vector of log odds ratios or a
#' meta-analysis object created with \code{\link{metabin}} or
#' \code{\link{metagen}} and the odds ratio as summary measure.
#'
#' Argument \code{selnOR} is mandatory if argument \code{lnOR} is a
#' vector and ignored otherwise. Additional arguments in \code{\dots}
#' are only passed on to \code{\link{metagen}} if argument \code{lnOR}
#' is a vector.
#' 
#' @return
#' An object of class \code{"meta"} and \code{"metagen"}; see
#' \code{\link{metagen}}.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#'
#' @seealso \code{\link{smd2or}}, \code{\link{metabin}},
#'   \code{\link{metagen}}, \code{\link{metacont}}
#' 
#' @references
#' Borenstein M, Hedges LV, Higgins JPT, Rothstein HR (2009):
#' \emph{Introduction to Meta-Analysis}.
#' Chichester: Wiley
#'
#' Cox DR (1970):
#' \emph{Analysis of Binary Data}.
#' London: Chapman and Hall / CRC
#'
#' Cox DR, Snell EJ (1989):
#' \emph{Analysis of Binary Data} (2nd edition).
#' London: Chapman and Hall / CRC
#'
#' Hasselblad V, Hedges LV (1995):
#' Meta-analysis of screening and diagnostic tests.
#' \emph{Psychological Bulletin},
#' \bold{117}, 167--78
#' 
#' @examples
#' # Example from Borenstein et al. (2010), Chapter 7
#' #
#' mb <- or2smd(0.9069, sqrt(0.0676))
#' # TE = standardised mean difference (SMD); seTE = standard error of SMD
#' data.frame(SMD = round(mb$TE, 4), varSMD = round(mb$seTE^2, 4))
#'
#' # Use dataset from Fleiss (1993)
#' #
#' data(Fleiss93)
#' m1 <- metabin(event.e, n.e, event.c, n.c,
#'               data = Fleiss93,
#'               studlab = paste(study, year),
#'               sm = "OR", comb.random = FALSE)
#' or2smd(m1)
#' 
#' @export or2smd


or2smd <- function(lnOR, selnOR, method = "HH", ...) {
  
  
  is.meta <- inherits(lnOR, "meta")
  ##
  if (is.meta) {
    if (lnOR$sm != "OR")
      stop("Effect measure must be equal to 'OR'.", call. = FALSE)
    else {
      mdat <- lnOR
      lnOR <- mdat$TE
      selnOR <- mdat$seTE
    }
  }
  
  
  method <- setchar(method, c("HH", "CS"))
  ##  
  if (method == "HH") {
    smd <- lnOR * sqrt(3) / pi
    se.smd <- sqrt(selnOR^2 * 3 / pi^2)
  }
  else if (method == "CS") {
    smd <- lnOR / 1.65
    se.smd <- sqrt(selnOR^2 / 1.65)
  }
  
  
  if (is.meta) {
    if (is.null(mdat$byvar))
      res <- metagen(smd, se.smd, sm = "SMD",
                     data = mdat,
                     studlab = mdat$studlab,
                     subset = mdat$subset, exclude = mdat$exclude,
                     level = mdat$level, level.comb = mdat$level.comb,
                     comb.fixed = mdat$comb.fixed,
                     comb.random = mdat$comb.random,
                     hakn = mdat$hakn, method.tau = mdat$method.tau,
                     tau.common = mdat$tau.common,
                     prediction = mdat$prediction,
                     level.predict = mdat$level.predict,
                     null.effect = 0,
                     method.bias = mdat$method.bias,
                     title = mdat$title, complab = mdat$complab,
                     outclab = mdat$outclab,
                     label.c = mdat$label.c, label.e = mdat$label.e,
                     label.left = mdat$label.left,
                     label.right = mdat$label.right,
                     control = mdat$control)
    else
      res <- metagen(smd, se.smd, sm = "SMD",
                     data = mdat,
                     studlab = mdat$studlab,
                     subset = mdat$subset, exclude = mdat$exclude,
                     level = mdat$level, level.comb = mdat$level.comb,
                     comb.fixed = mdat$comb.fixed,
                     comb.random = mdat$comb.random,
                     hakn = mdat$hakn, method.tau = mdat$method.tau,
                     tau.common = mdat$tau.common,
                     prediction = mdat$prediction,
                     level.predict = mdat$level.predict,
                     null.effect = 0,
                     method.bias = mdat$method.bias,
                     title = mdat$title, complab = mdat$complab,
                     outclab = mdat$outclab,
                     label.c = mdat$label.c,
                     label.e = mdat$label.e,
                     label.left = mdat$label.left,
                     label.right = mdat$label.right,
                     byvar = mdat$byvar, bylab = mdat$bylab,
                     print.byvar = mdat$print.byvar,
                     byseparator = mdat$byseparator,
                     control = mdat$control)
  }
  else {
    dat <- data.frame(lnOR, selnOR, OR = exp(lnOR), smd, se.smd)
    res <- metagen(smd, se.smd, data = dat, sm = "SMD", ...)
  }
  
  
  res
}
