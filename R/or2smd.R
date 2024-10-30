#' Conversion from log odds ratio to standardised mean difference
#' 
#' @description
#' Conversion from log odds ratio to standardised mean difference
#' using method by Hasselblad & Hedges (1995) or Cox (1970).
#' 
#' @param lnOR Log odds ratio(s) or meta-analysis object.
#' @param selnOR Standard error(s) of log odds ratio(s) (ignored if
#'   argument \code{lnOR} is a meta-analysis object).
#' @param studlab An optional vector with study labels (ignored if
#'   argument \code{lnOR} is a meta-analysis object).
#' @param data An optional data frame containing the study information
#'   (ignored if argument \code{lnOR} is a meta-analysis object).
#' @param subset An optional vector specifying a subset of studies to
#'   be used (ignored if argument \code{lnOR} is a meta-analysis
#'   object).
#' @param exclude An optional vector specifying studies to exclude
#'   from meta-analysis, however, to include in printouts and forest
#'   plots (ignored if argument \code{lnOR} is a meta-analysis
#'   object).
#' @param method A character string indicating which method is used to
#'   convert log odds ratios to standardised mean differences. Either
#'   \code{"HH"} or \code{"CS"}, can be abbreviated.
#' @param \dots Additional arguments passed on to
#'   \code{\link{metagen}} (ignored if argument \code{lnOR} is a
#'   meta-analysis object).
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
#' An object of class \code{c("metagen", "meta")} with corresponding
#' generic functions (see \code{\link{meta-object}}).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
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
#' # Example from Borenstein et al. (2009), Chapter 7
#' #
#' mb <- or2smd(0.9069, sqrt(0.0676))
#' # TE = standardised mean difference (SMD); seTE = standard error of SMD
#' data.frame(SMD = round(mb$TE, 4), varSMD = round(mb$seTE^2, 4))
#'
#' # Use dataset from Fleiss (1993)
#' #
#' data(Fleiss1993bin)
#' m1 <- metabin(d.asp, n.asp, d.plac, n.plac,
#'   data = Fleiss1993bin, studlab = paste(study, year),
#'   sm = "OR", random = FALSE)
#' or2smd(m1)
#' 
#' @export or2smd


or2smd <- function(lnOR, selnOR, studlab,
                   data = NULL, subset = NULL, exclude = NULL,
                   method = "HH", ...) {
  
  
  is.meta <- inherits(lnOR, "meta")
  ##
  if (is.meta) {
    lnOR <- updateversion(lnOR)
    if (lnOR$sm != "OR")
      stop("Effect measure must be equal to 'OR'.", call. = FALSE)
    else {
      mdat <- lnOR
      lnOR <- mdat$TE
      selnOR <- mdat$seTE
    }
  }
  else {
    ##
    ## Read data
    ##
    nulldata <- is.null(data)
    sfsp <- sys.frame(sys.parent())
    mc <- match.call()
    ##
    if (nulldata)
      data <- sfsp
    ##
    ## Catch 'lnOR' and 'selnOR' from data:
    ##
    lnOR <- catch("lnOR", mc, data, sfsp)
    ##
    selnOR <- catch("selnOR", mc, data, sfsp)
    ##
    k.All <- length(lnOR)
    chknull(lnOR)
    ##
    ## Catch 'studlab', 'subset', and 'exclude' from data:
    ##
    studlab <- catch("studlab", mc, data, sfsp)
    studlab <- setstudlab(studlab, k.All)
    print(studlab)
    ##
    missing.subset <- missing(subset)
    subset <- catch("subset", mc, data, sfsp)
    ##
    missing.exclude <- missing(exclude)
    exclude <- catch("exclude", mc, data, sfsp)
    ##
    ## Check length of essential variables
    ##
    arg <- "lnOR"
    ##
    chklength(selnOR, k.All, arg)
    chklength(studlab, k.All, arg)
    ##
    ## Subset and exclude studies
    ##
    if (!missing.subset)
      if ((is.logical(subset) & (sum(subset) > k.All)) ||
          (length(subset) > k.All))
        stop("Length of argument 'subset' is larger than number of studies.")
    ##
    if (!missing.exclude)
      if ((is.logical(exclude) & (sum(exclude) > k.All)) ||
          (length(exclude) > k.All))
        stop("Length of argument 'exclude' is larger than number of studies.")
  }
  
  
  method <- setchar(method, c("HH", "CS"))
  ##  
  if (method == "HH") {
    smd <- lnOR * sqrt(3) / pi
    se.smd <- selnOR * sqrt(3) / pi
  }
  else if (method == "CS") {
    smd <- lnOR / 1.65
    se.smd <- selnOR / 1.65
  }
  
  
  if (is.meta) {
    if (is.null(mdat$subgroup))
      res <- metagen(smd, se.smd, studlab = mdat$studlab,
                     data = mdat,
                     subset = mdat$subset, exclude = mdat$exclude,
                     cluster = mdat$cluster,
                     ##
                     sm = "SMD",
                     ##
                     level = mdat$level, level.ma = mdat$level.ma,
                     common = mdat$common,
                     random = mdat$random,
                     overall = mdat$overall,
                     overall.hetstat = mdat$overall.hetstat,
                     ##
                     method.random.ci = mdat$method.random.ci,
                     adhoc.hakn.ci = mdat$adhoc.hakn.ci,
                     ##
                     prediction = mdat$prediction,
                     method.predict = mdat$method.predict,
                     adhoc.hakn.pi = mdat$adhoc.hakn.pi,
                     level.predict = mdat$level.predict,
                     ##
                     method.tau = mdat$method.tau,
                     method.tau.ci = mdat$method.tau.ci,
                     level.hetstat = mdat$level.hetstat,
                     tau.common = mdat$tau.common,
                     detail.tau = mdat$detail.tau,
                     ##
                     null.effect = 0,
                     ##
                     method.bias = mdat$method.bias,
                     ##
                     text.common = mdat$text.common,
                     text.random = mdat$text.random,
                     text.predict = mdat$text.predict,
                     text.w.common = mdat$text.w.common,
                     text.w.random = mdat$text.w.random,
                     ##
                     title = mdat$title, complab = mdat$complab,
                     outclab = mdat$outclab,
                     #
                     label.e = mdat$label.e, label.c = mdat$label.c,
                     label.left = mdat$label.left,
                     label.right = mdat$label.right,
                     col.label.left = mdat$col.label.left,
                     col.label.right = mdat$col.label.right,
                     #
                     control = mdat$control)
    else
      res <- metagen(smd, se.smd, studlab = mdat$studlab,
                     data = mdat,
                     subset = mdat$subset, exclude = mdat$exclude,
                     cluster = mdat$cluster,
                     ##
                     sm = "SMD",
                     ##
                     level = mdat$level, level.ma = mdat$level.ma,
                     common = mdat$common,
                     random = mdat$random,
                     overall = mdat$overall,
                     overall.hetstat = mdat$overall.hetstat,
                     ##
                     method.random.ci = mdat$method.random.ci,
                     adhoc.hakn.ci = mdat$adhoc.hakn.ci,
                     ##
                     prediction = mdat$prediction,
                     method.predict = mdat$method.predict,
                     adhoc.hakn.pi = mdat$adhoc.hakn.pi,
                     level.predict = mdat$level.predict,
                     ##
                     method.tau = mdat$method.tau,
                     method.tau.ci = mdat$method.tau.ci,
                     level.hetstat = mdat$level.hetstat,
                     tau.common = mdat$tau.common,
                     detail.tau = mdat$detail.tau,
                     ##
                     null.effect = 0,
                     ##
                     method.bias = mdat$method.bias,
                     ##
                     text.common = mdat$text.common,
                     text.random = mdat$text.random,
                     text.predict = mdat$text.predict,
                     text.w.common = mdat$text.w.common,
                     text.w.random = mdat$text.w.random,
                     ##
                     title = mdat$title, complab = mdat$complab,
                     outclab = mdat$outclab,
                     #
                     label.e = mdat$label.e, label.c = mdat$label.c,
                     label.left = mdat$label.left,
                     label.right = mdat$label.right,
                     col.label.left = mdat$col.label.left,
                     col.label.right = mdat$col.label.right,
                     #
                     subgroup = mdat$subgroup,
                     subgroup.name = mdat$subgroup.name,
                     print.subgroup.name = mdat$print.subgroup.name,
                     sep.subgroup = mdat$sep.subgroup,
                     test.subgroup = mdat$test.subgroup,
                     prediction.subgroup = mdat$prediction.subgroup,
                     ##
                     control = mdat$control)
  }
  else {
    dat <- data.frame(lnOR, selnOR, OR = exp(lnOR), smd, se.smd, studlab)
    dat$subset <- subset
    dat$exclude <- exclude
    ##
    res <- metagen(smd, se.smd, studlab = studlab,
                   data = dat, subset = subset, exclude = exclude,
                   sm = "SMD", ...)
  }
  ##
  res$method.or2smd <- method
  
  
  res
}
