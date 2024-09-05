#' Conversion from standardised mean difference to log odds ratio
#' 
#' @description
#' Conversion from standardised mean difference to log odds ratio
#' using method by Hasselblad & Hedges (1995) or Cox (1970).
#' 
#' @param smd Standardised mean difference(s) (SMD) or meta-analysis
#'   object.
#' @param se.smd Standard error(s) of SMD (ignored if argument
#'   \code{smd} is a meta-analysis object).
#' @param studlab An optional vector with study labels (ignored if
#'   argument \code{smd} is a meta-analysis object).
#' @param data An optional data frame containing the study information
#'   (ignored if argument \code{smd} is a meta-analysis object).
#' @param subset An optional vector specifying a subset of studies to
#'   be used (ignored if argument \code{smd} is a meta-analysis
#'   object).
#' @param exclude An optional vector specifying studies to exclude
#'   from meta-analysis, however, to include in printouts and forest
#'   plots (ignored if argument \code{smd} is a meta-analysis object).
#' @param method A character string indicating which method is used to
#'   convert SMDs to log odds ratios. Either \code{"HH"} or
#'   \code{"CS"}, can be abbreviated.
#' @param backtransf A logical indicating whether odds ratios (if
#'   TRUE) or log odds ratios (if FALSE) should be shown in printouts
#'   and plots.
#' @param \dots Additional arguments passed on to
#'   \code{\link{metagen}} (ignored if argument \code{smd} is a
#'   meta-analysis object).
#' 
#' @details
#' This function implements the following methods for the conversion
#' from standardised mean difference to log odds ratio:
#' \itemize{
#' \item Hasselblad & Hedges (1995) assuming logistic distributions
#'   (\code{method == "HH"})
#' \item Cox (1970) and Cox & Snell (1989) assuming normal
#'   distributions (\code{method == "CS"})
#' }
#' Internally, \code{\link{metagen}} is used to conduct a
#' meta-analysis with the odds ratio as summary measure.
#' 
#' Argument \code{smd} can be either a vector of standardised mean
#' differences or a meta-analysis object created with
#' \code{\link{metacont}} or \code{\link{metagen}} and the
#' standardised mean difference as summary measure.
#'
#' Argument \code{se.smd} is mandatory if argument \code{smd} is a
#' vector and ignored otherwise. Additional arguments in \code{\dots}
#' are only passed on to \code{\link{metagen}} if argument \code{smd}
#' is a vector.
#' 
#' @return
#' An object of class \code{c("metagen", "meta")} with corresponding
#' generic functions (see \code{\link{meta-object}}).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @seealso \code{\link{or2smd}}, \code{\link{metacont}},
#'   \code{\link{metagen}}, \code{\link{metabin}}
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
#' mb <- smd2or(0.5, sqrt(0.0205), backtransf = FALSE)
#' # TE = log odds ratio; seTE = standard error of log odds ratio
#' data.frame(lnOR = round(mb$TE, 4), varlnOR = round(mb$seTE^2, 4))
#'
#' # Use dataset from Fleiss (1993)
#' #
#' data(Fleiss1993cont)
#' m1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, sm = "SMD",
#'                studlab = paste(study, year))
#' smd2or(m1)
#' 
#' @export smd2or


smd2or <- function(smd, se.smd, studlab,
                   data = NULL, subset = NULL, exclude = NULL,
                   method = "HH", backtransf = gs("backtransf"), ...) {
  
  
  is.meta <- inherits(smd, "meta")
  ##
  if (is.meta) {
    smd <- updateversion(smd)
    if (smd$sm != "SMD")
      stop("Effect measure must be equal to 'SMD'.", call. = FALSE)
    else {
      mdat <- smd
      smd <- mdat$TE
      se.smd <- mdat$seTE
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
    ## Catch 'smd' and 'se.smd' from data:
    ##
    smd <- catch("smd", mc, data, sfsp)
    ##
    se.smd <- catch("se.smd", mc, data, sfsp)
    ##
    k.All <- length(smd)
    chknull(smd)
    ##
    ## Catch 'studlab', 'subset', and 'exclude' from data:
    ##
    studlab <- catch("studlab", mc, data, sfsp)
    studlab <- setstudlab(studlab, k.All)
    ##
    missing.subset <- missing(subset)
    subset <- catch("subset", mc, data, sfsp)
    ##
    missing.exclude <- missing(exclude)
    exclude <- catch("exclude", mc, data, sfsp)
    ##
    ## Check length of essential variables
    ##
    arg <- "smd"
    ##
    chklength(se.smd, k.All, arg)
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
    lnOR <- smd * pi / sqrt(3)
    selnOR <- se.smd * pi / sqrt(3)
  }
  else if (method == "CS") {
    lnOR <- smd * 1.65
    selnOR <- se.smd * 1.65
  }
  ##
  chklogical(backtransf)
  
  
  if (is.meta) {
    if (is.null(mdat$subgroup))
      res <- metagen(lnOR, selnOR, studlab = mdat$studlab,
                     data = mdat,
                     subset = mdat$subset, exclude = mdat$exclude,
                     cluster = mdat$cluster,
                     ##
                     sm = "OR",
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
                     backtransf = backtransf,
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
      res <- metagen(lnOR, selnOR, studlab = mdat$studlab,
                     data = mdat,
                     subset = mdat$subset, exclude = mdat$exclude,
                     cluster = mdat$cluster,
                     ##
                     sm = "OR",
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
                     backtransf = backtransf,
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
    dat <- data.frame(smd, se.smd, lnOR, selnOR, OR = exp(lnOR), studlab)
    dat$subset <- subset
    dat$exclude <- exclude
    ##
    res <- metagen(lnOR, selnOR, studlab = studlab,
                   data = dat, subset = subset, exclude = exclude,
                   sm = "OR", backtransf = backtransf, ...)
  }
  ##
  res$method.smd2or <- method
  
  res
}
