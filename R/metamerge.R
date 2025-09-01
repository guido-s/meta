#' Merge results of two meta-analyses on the same data set
#' 
#' @description
#' This function can be used to merge results of two meta-analyses
#' into a single meta-analysis object if they are based on the same
#' data set. This is, for example, useful to produce a forest plot of
#' a random-effects meta-analysis with different estimates of the
#' between-study variance \eqn{\tau^2}.
#' 
#' @param meta1 First meta-analysis object (see Details).
#' @param meta2 Second meta-analysis object (see Details).
#' @param common1 A logical indicating whether results of common
#'   effect model should be considered for first meta-analysis.
#' @param random1 A logical indicating whether results of random
#'   effects model should be considered for first meta-analysis.
#' @param prediction1 A logical indicating whether prediction interval
#'   should be considered for first meta-analysis.
#' @param common2 A logical indicating whether results of common
#'   effect model should be considered for second meta-analysis.
#' @param random2 A logical indicating whether results of random
#'   effects model should be considered for second meta-analysis.
#' @param prediction2 A logical indicating whether prediction interval
#'   should be considered for second meta-analysis.
#' @param label1 Default setting for arguments 'label1.common',
#'   'label1.random', 'label1.predict' and 'label1.subgroup'.
#' @param label2 Default setting for arguments 'label2.common',
#'   'label2.random', 'label2.predict' and 'label2.subgroup'.
#' @param label1.common A character string to label the common effect estimate
#'   from the first meta-analysis.
#' @param label2.common A character string to label the common effect estimate
#'   from the second meta-analysis.
#' @param label1.random A character string to label the random effects estimate
#'   from the first meta-analysis (default label for arguments 'hetlabel1' and
#'   'taulabel1').
#' @param label2.random A character string to label the random effects estimate
#'   from the second meta-analysis (default label for arguments 'hetlabel2' and
#'   'taulabel2').
#' @param label1.predict A character string to label the prediction interval
#'   from the first meta-analysis.
#' @param label2.predict A character string to label the prediction interval
#'   from the second meta-analysis.
#' @param label1.subgroup A character string to label the subgroup results
#'   from the first meta-analysis.
#' @param label2.subgroup A character string to label the subgroup results
#'   from the second meta-analysis.
#' @param text.pooled1 A character string used in printouts and forest
#'   plot to label the results from the first meta-analysis.
#' @param text.pooled2 A character string used in printouts and forest
#'   plot to label the results from the second meta-analysis.
#' @param text.w.pooled1 A character string used to label weights of
#'   the first meta-analysis; can be of same length as the number of
#'   pooled estimates requested in argument \code{pooled1}.
#' @param text.w.pooled2 A character string used to label weights of
#'   the second meta-analysis; can be of same length as the number of
#'   pooled estimates requested in argument \code{pooled1}.
#' @param hetlabel1 A character string used to label heterogeneity
#'   statistics of the first meta-analysis.
#' @param hetlabel2 A character string used to label heterogeneity
#'   statistics of the second meta-analysis.
#' @param taulabel1 A character string used to label estimate of
#'   between-study variance of the first meta-analysis.
#' @param taulabel2 A character string used to label estimate of
#'   between-study variance of the second meta-analysis.
#' @param text.common1 A character string used in printouts and forest
#'   plot to label results for common effect models from the first
#'   meta-analysis.
#' @param text.common2 A character string used in printouts and forest
#'   plot to label results for common effect models from the second
#'   meta-analysis.
#' @param text.random1 A character string used in printouts and forest
#'   plot to label results for random effects models from the first
#'   meta-analysis.
#' @param text.random2 A character string used in printouts and forest
#'   plot to label results for random effects models from the second
#'   meta-analysis.
#' @param text.predict1 A character string used in printouts and
#'   forest plot to label prediction interval from the first
#'   meta-analysis.
#' @param text.predict2 A character string used in printouts and
#'   forest plot to label prediction interval from the second
#'   meta-analysis.
#' @param text.w.common1 A character string used to label common
#'   effect weights of the first meta-analysis; can be of same length
#'   as the number of common effect estimates.
#' @param text.w.common2 A character string used to label common
#'   effect weights of the second meta-analysis; can be of same length
#'   as the number of common effect estimates.
#' @param text.w.random1 A character string used to label random
#'   effects weights of the first meta-analysis; can be of same length
#'   as the number of random effects estimates.
#' @param text.w.random2 A character string used to label random
#'   effects weights of the second meta-analysis; can be of same
#'   length as the number of random effects estimates.
#' @param keep A logical indicating whether to keep additional
#'   information from second meta-analysis.
#' @param keep.Q A logical indicating whether heterogeneity statistic
#'   Q of second meta-analysis should be kept or ignored.
#' @param keep.I2 A logical indicating whether heterogeneity statistic
#'   I2 of second meta-analysis should be kept or ignored.
#' @param keep.w A logical indicating whether weights of the second
#'   meta-analysis should be kept or ignored.
#' @param common A logical indicating whether results of common effect
#'   meta-analyses should be reported.
#' @param random A logical indicating whether results of random
#'   effects meta-analyses should be reported.
#' @param overall A logical indicating whether overall summaries
#'   should be reported.
#' @param overall.hetstat A logical value indicating whether to print
#'   heterogeneity measures for overall treatment comparisons.
#' @param prediction A logical indicating whether prediction intervals
#'   should be reported.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=TRUE} (default), results for \code{sm="OR"} are
#'   printed as odds ratios rather than log odds ratios, for example.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param pooled1 Deprecated argument (replaced by 'common1',
#'   'random1', 'prediction1'). A character string indicating whether
#'   results of common effect or random effects model should be
#'   considered for first meta-analysis. Either \code{"both"},
#'   \code{"common"} or \code{"random"}, can be abbreviated.
#' @param pooled2 Deprecated argument (replaced by 'common2',
#'   'random2', 'prediction2'). A character string indicating whether
#'   results of common effect or random effects model should be
#'   considered for second meta-analysis. Either \code{"both"},
#'   \code{"common"} or \code{"random"}, can be abbreviated.
#' 
#' @details
#' In R package \bold{meta}, objects of class \code{"meta"} contain
#' results of both common effect and random effects
#' meta-analyses. This function enables the user to merge the results
#' of two meta-analysis object if they are based on the same data set.
#'
#' Applications of this function include printing and plotting results
#' of the common effect or random effects meta-analysis and the
#' \itemize{
#' \item trim-and-fill method (\code{\link{trimfill}}),
#' \item limit meta-analyis (\code{\link[metasens]{limitmeta}} from R
#'   package \bold{metasens}),
#' \item Copas selection model (\code{\link[metasens]{copas}} from R
#'   package \bold{metasens}),
#' \item robust variance meta-analysis model
#'   (\code{\link[robumeta]{robu}} from R package \bold{robumeta}).
#' }
#'
#' The first argument (\code{meta1}) must be an object created by a
#' meta-analysis function (see \code{\link{meta-object}}). If an
#' object created with \code{\link[metasens]{limitmeta}} or
#' \code{\link[metasens]{copas}} is provided as the first argument,
#' this object will be returned, i.e., argument \code{meta2} will be
#' ignored.
#'
#' The second meta-analysis could be an object created by a
#' meta-analysis function or with \code{\link{trimfill}},
#' \code{\link[metasens]{limitmeta}}, \code{\link[metasens]{copas}},
#' or \code{\link[robumeta]{robu}}.
#'
#' The created meta-analysis object only contains the study results,
#' i.e., estimated effects and confidence intervals, from the first
#' meta-analysis which are shown in printouts and forest plots. This
#' only makes a difference for meta-analysis methods where individual
#' study results differ, e.g., Mantel-Haenszel and Peto method for
#' binary outcomes (see \code{\link{metabin}}).
#'
#' R function \code{\link{metaadd}} can be used to add pooled results
#' from any (external) meta-analysis.
#'
#' R function \code{\link{metabind}} can be used to print and plot the
#' results of several meta-analyses without the restriction that the
#' same data set has to be used. Accordingly, individual study results
#' are ignored.
#' 
#' @return
#' An object of class \code{"meta"} and \code{"metamerge"} with
#' corresponding generic functions (see \code{\link{meta-object}}).
#' 
#' The following list elements have a different meaning:
#' \item{TE, seTE, studlab}{Treatment estimate, standard error, and
#'   study labels (first meta-analysis).}
#' \item{lower, upper}{Lower and upper confidence interval limits for
#'   individual studies (first meta-analysis).}
#' \item{statistic, pval}{Statistic and p-value for test of treatment
#'   effect for individual studies (first meta-analysis.}
#' \item{w.common}{Vector or matrix with common effect weights.}
#' \item{w.random}{Vector or matrix with random effects weights.}
#' \item{k}{Vector with number of estimates (same length as number of
#'   common effect and random effects estimates).}
#' \item{k.study}{Vector with number of studies (same length as
#'   number of common effect and random effects estimates).}
#' \item{k.all}{Vector with total number of studies (same length as
#'   number of common effect and random effects estimates).}
#' \item{k.TE}{Vector with number of studies with estimable effects
#'   (same length as number of common effect and random effects
#'   estimates).}
#' \item{k.MH}{Vector with number of studies combined with
#'   Mantel-Haenszel method (same length as number of common effect
#'   and random effects estimates).}
#' \item{TE.common}{Vector with common effect estimates.}
#' \item{seTE.common}{Vector with standard errors of common effect
#'   estimates.}
#' \item{lower.common}{Vector with lower confidence limits (common
#'   effect model).}
#' \item{upper.common}{Vector with upper confidence limits (common
#'   effect model).}
#' \item{statistic.common}{Vector with test statistics for test of
#'   overall effect (common effect model).}
#' \item{pval.common}{Vector with p-value of test for overall effect
#'   (common effect model).}
#' \item{TE.random}{Vector with random effects estimates.}
#' \item{seTE.random}{Vector with standard errors of random effects
#'   estimates.}
#' \item{lower.random}{Vector with lower confidence limits (random
#'   effects model).}
#' \item{upper.random}{Vector with upper confidence limits (random
#'   effects model).}
#' \item{statistic.random}{Vector with test statistics for test of
#'   overall effect (random effects model).}
#' \item{pval.random}{Vector with p-value of test for overall effect
#'   (random effects model).}
#' 
#' Furthermore, meta-analysis results of common effect or random
#' effects model are taken from first meta-analysis if only random
#' effects or common effects models are selected from both
#' meta-analyses (arguments \code{pooled1} and \code{pooled2}).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metagen}}, \code{\link{metabind}},
#'   \code{\link{metaadd}}
#' 
#' @examples
#' # Print results with more significant digits and do not show confidence
#' # intervals for tau^2 and tau
#' oldset <- settings.meta(digits = 6, digits.stat = 4, digits.pval = 6,
#'   digits.Q = 6, digits.I2 = 4, digits.H = 4,
#'   print.tau2.ci = FALSE, print.tau.ci = FALSE)
#' oldopts <- options(width = 120)
#' 
#' data(Fleiss1993bin)
#' 
#' # Mantel-Haenszel method
#' m1 <- metabin(d.asp, n.asp, d.plac, n.plac, data = Fleiss1993bin,
#'   studlab = paste(study, year), sm = "OR")
#' # Peto method
#' m2 <- update(m1, method = "Peto")
#' # Inverse variance method (only common effect model)
#' m3 <- update(m2, method = "Inverse", random = FALSE)
#' 
#' # Merge results from MH and Peto method
#' # - show individual results for MH method
#' #   (as this is the first meta-analysis)
#' # - keep all additional information from Peto meta-analysis (i.e.,
#' #   weights, Q statistic and I2 statistic)
#' m12 <- metamerge(m1, m2,
#'   label1 = "REML", label2 = "REML-Peto",
#'   label1.common = "MH", label2.common = "Peto", 
#'   text.common1 = "Mantel-Haenszel method",
#'   text.common2 = "Peto method",
#'   text.w.random1 = "REML", text.w.random2 = "REML-Peto",
#'   hetlabel1 = "MH/IV", hetlabel2 = "Peto",
#'   keep = TRUE)
#' 
#' # Add common effect results from inverse variance method
#' # - keep weights from IV meta-analysis
#' # - Q and I2 statistic are identical for sm = "MH" and sm = "Inverse"
#' #   as inverse variance method is used for sm = "MH" under random
#' #   effects model
#' m123 <- metamerge(m12, m3,
#'   label2 = "IV",
#'   text.common2 = "Inverse variance method",
#'   keep.w = TRUE)
#' summary(m123)
#' \dontrun{
#' forest(m123, digits = 6)
#' 
#' # Merge results (show individual results for Peto method)
#' m21 <- metamerge(m2, m1,
#'   label1 = "REML-Peto", label2 = "REML",
#'   label1.common = "Peto", label2.common = "MH", 
#'   hetlabel1 = "Peto", hetlabel2 = "MH/IV",
#'   text.common1 = "Peto method",
#'   text.common2 = "Mantel-Haenszel method",
#'   keep = TRUE)
#' 
#' # Add results from inverse variance method
#' # - keep weights from IV meta-analysis
#' # - Q and I2 statistic are identical for sm = "MH" and sm = "Inverse"
#' #   as inverse variance method is used for sm = "MH" under random
#' #   effects model
#' m213 <- metamerge(m21, m3,
#'   label2 = "IV",
#'   text.common2 = "Inverse variance method",
#'   keep.w = TRUE)
#' summary(m213)
#' 
#' # Random effects method using ML estimator for between-study variance tau2
#' m4 <- update(m1, common = FALSE, method.tau = "ML")
#' 
#' # Use DerSimonian-Laird estimator for tau2
#' m5 <- update(m4, method.tau = "DL")
#' 
#' # Use Paule-Mandel estimator for tau2
#' m6 <- update(m4, method.tau = "PM")
#' 
#' # Merge random effects results for ML and DL estimators
#' # - keep weights for DL estimator (which are different from ML)
#' m45 <- metamerge(m4, m5, label1 = "ML", label2 = "DL",
#'   text.w.random1 = "RE-ML", text.w.random2 = "RE-DL", keep.w = TRUE)
#' summary(m45)
#' 
#' # Add results for PM estimator
#' # - keep weights
#' m456 <- metamerge(m45, m6, label2 = "PM",
#'   text.w.random2 = "RE-PM", keep.w = TRUE)
#' summary(m456)
#' 
#' m123456 <- metamerge(m123, m456)
#' m123456
#' 
#' # Use Hartung-Knapp confidence intervals
#' # - do not keep information on Q, I2 and weights
#' m7 <- update(m4, method.random.ci = "HK",
#'   text.random = "Hartung-Knapp method")
#' m8 <- update(m5, method.random.ci = "HK",
#'   text.random = "Hartung-Knapp method")
#' m9 <- update(m6, method.random.ci = "HK",
#'   text.random = "Hartung-Knapp method")
#' 
#' # Merge results for Hartung-Knapp method (with REML and DL estimator)
#' # - RE weights for REML estimator are shown
#' m78 <- metamerge(m7, m8, label1 = "ML", label2 = "DL")
#' summary(m78)
#' 
#' m789 <- metamerge(m78, m9, label2 = "PM")
#' summary(m789)
#' 
#' # Merge everything
#' m1to9 <- metamerge(metamerge(m123, m456, keep.w = TRUE), m789)
#' summary(m1to9)
#' 
#' m10 <- update(m1, method = "GLMM")
#' 
#' m.all <- metamerge(m1to9, m10, keep.Q = TRUE,
#'   label2 = "GLMM", taulabel2 = "ML-GLMM")
#' summary(m.all)
#' 
#' forest(m.all, layout = "JAMA")
#' forest(m.all, details = TRUE)
#' }
#' 
#' settings.meta(oldset)
#' options(oldopts)
#' 
#' @export metamerge


metamerge <- function(meta1, meta2,
                      ##
                      common1 = meta1$common,
                      random1 = meta1$random,
                      prediction1 = meta1$prediction,
                      common2 = meta2$common,
                      random2 = meta2$random,
                      prediction2 = meta2$prediction,
                      ##
                      label1 = NULL, label2 = NULL,
                      label1.common = label1, label2.common = label2,
                      label1.random = label1, label2.random = label2,
                      label1.predict = label1, label2.predict = label2,
                      label1.subgroup = label1, label2.subgroup = label2,
                      ##
                      hetlabel1 = label1.random,
                      hetlabel2 = label2.random,
                      taulabel1 = label1.random,
                      taulabel2 = label2.random,
                      ##
                      text.pooled1 = NULL, text.pooled2 = NULL,
                      text.w.pooled1 = NULL, text.w.pooled2 = NULL,
                      ##
                      text.common1 = text.pooled1,
                      text.common2 = text.pooled2,
                      text.random1 = text.pooled1,
                      text.random2 = text.pooled2,
                      text.predict1 = text.pooled1,
                      text.predict2 = text.pooled2,
                      ##
                      text.w.common1 = text.w.pooled1,
                      text.w.common2 = text.w.pooled2,
                      text.w.random1 = text.w.pooled1,
                      text.w.random2 = text.w.pooled2,
                      ##
                      keep = FALSE,
                      keep.Q = keep, keep.I2 = keep.Q,
                      keep.w = keep,
                      ##
                      common = common1 | common2,
                      random = random1 | random2,
                      overall = common | random,
                      overall.hetstat = common | random,
                      prediction = prediction1 | prediction2,
                      ##
                      backtransf,
                      ##
                      warn.deprecated = gs("warn.deprecated"),
                      pooled1, pooled2) {
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  
  if (missing(meta1))
    stop("Argument 'meta1' must be provided.",
         call. = FALSE)
  ##
  if (missing(meta2)) {
    if (!inherits(meta1, "metamerge") & inherits(meta1, "copas"))
      return(metamerge(update(meta1$x), meta1,
                       label2 = if (is.null(label2)) "copas" else label2))
    ##
    else if (!inherits(meta1, "metamerge") & inherits(meta1, "limitmeta"))
      return(metamerge(update(meta1$x), meta1,
                       label2 = if (is.null(label2)) "limit" else label2))
    ##
    else if (!inherits(meta1, "metamerge") & inherits(meta1, "robu"))
      stop("Argument 'meta1' cannot be of class 'robu'.",
           call. = FALSE)
    ##
    else
      return(meta1)
  }
  else if (!inherits(meta1, "metamerge")) {
    if (inherits(meta1, "copas"))
      stop("Argument 'meta1' cannot be of class 'copas' ",
           "(use argument 'meta2').",
           call. = FALSE)
    else if (inherits(meta1, "limitmeta"))
      stop("Argument 'meta1' cannot be of class 'limitmeta' ",
           "(use argument 'meta2').",
           call. = FALSE)
    else if (inherits(meta1, "robu"))
      stop("Argument 'meta1' cannot be of class 'robu' ",
           "(use argument 'meta2').",
           call. = FALSE)
  }
  ##
  chkclass(meta1, "meta")
  chksuitable(meta1, "R function metamerge()", classes = "metaprop",
              check.mlm = FALSE)
  meta1 <- updateversion(meta1)
  #
  chkclass(meta2, c("meta", "limitmeta", "copas", "robu"))
  chksuitable(meta2, "R function metamerge()", classes = "metaprop",
              check.mlm = FALSE)
  #
  if (inherits(meta1, "netpairwise") | inherits(meta2, "netpairwise"))
    stop("R objects of class 'netpairwise' cannot be merged.",
         call. = FALSE)
  ##
  if (inherits(meta2, "meta"))
    meta2 <- updateversion(meta2)
  ##
  is.copas <- !inherits(meta2, "metamerge") & inherits(meta2, "copas")
  is.limit <- !inherits(meta2, "metamerge") & inherits(meta2, "limitmeta")
  is.robu  <- !inherits(meta2, "metamerge") & inherits(meta2, "robu")
  ##
  if (is.copas | is.limit | is.robu) {
    common2 <- FALSE
    random2 <- TRUE
    prediction2 <- FALSE
    meta2$three.level <- FALSE
    meta2$rho <- NA
    ##
    if (is.copas)
      meta2$k <- length(meta2$TE)
    ##
    if (is.limit)
      meta2$tau2 <- meta2$tau^2
    ##
    keep <- FALSE
    keep.Q <- FALSE
    keep.I2 <- FALSE
    keep.w <- FALSE
  }
  ##
  chklogical(warn.deprecated)
  ##
  missing.pooled1 <- missing(pooled1)
  missing.pooled2 <- missing(pooled2)
  ##
  missing.common1 <- missing(common1)
  missing.random1 <- missing(random1)
  ##
  missing.common2 <- missing(common2)
  missing.random2 <- missing(random2)
  ##
  deprecated2(common1, missing(common1),
              pooled1, missing.pooled1, warn.deprecated)
  deprecated2(random1, missing(random1),
              pooled1, missing.pooled1, warn.deprecated)
  ##
  if (missing.common1 & !missing(pooled1)) {
    pooled1 <- setchar(pooled1, c("both", "common", "random", "fixed"))
    pooled1[pooled1 == "fixed"] <- "common"
    common1 <- pooled1 == "common"
  }
  ##
  if (missing.random1 & !missing(pooled1)) {
    pooled1 <- setchar(pooled1, c("both", "common", "random", "fixed"))
    random1 <- pooled1 %in% c("common", "both")
  }
  ##
  chklogical(common1)
  chklogical(random1)
  ##
  deprecated2(common2, missing(common2),
              pooled2, missing.pooled2, warn.deprecated)
  deprecated2(random2, missing(random2),
              pooled2, missing.pooled2, warn.deprecated)
  ##
  if (missing.common2 & !missing(pooled2)) {
    pooled2 <- setchar(pooled2, c("both", "common", "random", "fixed"))
    pooled2[pooled2 == "fixed"] <- "common"
    common2 <- pooled2 == "common"
  }
  ##
  if (missing.random2 & !missing(pooled2)) {
    pooled2 <- setchar(pooled2, c("both", "common", "random", "fixed"))
    random2 <- pooled2 %in% c("common", "both")
  }
  ##
  chklogical(common2)
  chklogical(random2)
  ##
  chklogical(keep)
  chklogical(keep.Q)
  chklogical(keep.I2)
  chklogical(keep.w)
  ##
  chkchar(label1, length = 1, NULL.ok = TRUE)
  chkchar(label2, length = 1, NULL.ok = TRUE)
  ##
  chkchar(label1.common, length = 1, NULL.ok = TRUE)
  chkchar(label2.common, length = 1, NULL.ok = TRUE)
  ##
  chkchar(label1.random, length = 1, NULL.ok = TRUE)
  chkchar(label2.random, length = 1, NULL.ok = TRUE)
  ##
  chkchar(label1.predict, length = 1, NULL.ok = TRUE)
  chkchar(label2.predict, length = 1, NULL.ok = TRUE)
  ##
  chkchar(hetlabel1, length = 1, NULL.ok = TRUE)
  chkchar(taulabel1, length = 1, NULL.ok = TRUE)
  ##
  chkchar(hetlabel2, length = 1, NULL.ok = TRUE)
  chkchar(taulabel2, length = 1, NULL.ok = TRUE)
  ##
  chkchar(text.pooled1, length = 1, NULL.ok = TRUE)
  chkchar(text.common1, length = 1, NULL.ok = TRUE)
  chkchar(text.random1, length = 1, NULL.ok = TRUE)
  chkchar(text.predict1, length = 1, NULL.ok = TRUE)
  ##
  chkchar(text.common2, length = 1, NULL.ok = TRUE)
  chkchar(text.pooled2, length = 1, NULL.ok = TRUE)
  chkchar(text.random2, length = 1, NULL.ok = TRUE)
  chkchar(text.predict2, length = 1, NULL.ok = TRUE)
  ##
  chkchar(text.w.pooled1, length = 1, NULL.ok = TRUE)
  chkchar(text.w.common1, length = 1, NULL.ok = TRUE)
  chkchar(text.w.random1, length = 1, NULL.ok = TRUE)
  ##
  chkchar(text.w.pooled2, length = 1, NULL.ok = TRUE)
  chkchar(text.w.common2, length = 1, NULL.ok = TRUE)
  chkchar(text.w.random2, length = 1, NULL.ok = TRUE)
  ##
  chklogical(common)
  chklogical(random)
  chklogical(overall)
  chklogical(overall.hetstat)
  chklogical(prediction)
  ##
  if (!missing(backtransf))
    chklogical(backtransf)
  else {
    if (!is.null(meta1$backtransf) & !is.null(meta2$backtransf))
      backtransf <- meta1$backtransf | meta2$backtransf
    else if (!is.null(meta1$backtransf))
      backtransf <- meta1$backtransf
    else if (!is.null(meta2$backtransf))
      backtransf <- meta2$backtransf
    else
      backtransf <- FALSE
  }
  ##
  ## Check summary measures
  ##
  if (!is.robu)
    samesm(meta1, meta2)
  ##
  ## Check original data (if available)
  ##
  if (!is.robu)
    samedata(meta1, meta2)
  ##
  ## Check subgroup levels
  ##
  samesubgroups(meta1, meta2)
  ##
  ## Check levels of confidence intervals
  ##
  lvls <- sort(unique(c(meta1$level, meta2$level)))
  lvls.ma <- sort(unique(c(meta1$level.ma, meta2$level.ma)))
  ##
  if (length(lvls) != 1)
    stop("Different level for confidence intervals of individual studies: ",
         paste0(lvls, collapse = ", "),
         call. = FALSE)
  if (length(lvls.ma) != 1)
    stop("Different level for meta-analysis confidence intervals: ",
         paste0(lvls.ma, collapse = ", "),
         call. = FALSE)
  
  
  ##
  ##
  ## (2) Some assignments for trim-and-fill method (meta1 / meta2),
  ##     Copas selection model, limit meta-analysis and robust
  ##     variance meta-analysis (meta2)
  ##
  
  meta1 <- updateobj(meta1,
                     label1.common, label1.random, label1.predict,
                     hetlabel1, taulabel1, label1.subgroup,
                     text.common1, text.random1, text.predict1,
                     text.w.common1, text.w.random1)
  ##
  hetlabel1 <- meta1$hetlabel
  taulabel1 <- meta1$taulabel
  ##
  meta2 <- updateobj(meta2,
                     label2.common, label2.random, label2.predict,
                     hetlabel2, taulabel2, label2.subgroup,
                     text.common2, text.random2, text.predict2,
                     text.w.common2, text.w.random2)
  ##
  hetlabel2 <- meta2$hetlabel
  taulabel2 <- meta2$taulabel
  
  
  ##
  ##
  ## (3) Some more assignments
  ##
  ##
  
  meta1$detail.tau <- replaceNULL(meta1$detail.tau, "")
  meta2$detail.tau <- replaceNULL(meta2$detail.tau, "")
  ##
  meta1$method.tau <- replaceNULL(meta1$method.tau, "")
  meta2$method.tau <- replaceNULL(meta2$method.tau, "")
  ##
  meta1$method.tau.ci <- replaceNULL(meta1$method.tau.ci, "")
  meta2$method.tau.ci <- replaceNULL(meta2$method.tau.ci, "")
  ##
  #meta1$method.random.ci <- replaceNULL(meta1$method.random.ci, "")
  #meta2$method.random.ci <- replaceNULL(meta2$method.random.ci, "")
  ##
  meta1$df.random <- replaceNULL(meta1$df.random, NA)
  meta2$df.random <- replaceNULL(meta2$df.random, NA)
  
  
  ##
  ##
  ## (4) Remove results from first meta-analysis (if necessary)
  ##
  ##
  
  if (!common1 & common2)
    meta1 <- dropcommon(meta1)
  ##
  if (!random1 & random2)
    meta1 <- droprandom(meta1)
  ##
  if (!prediction1 & prediction2)
    meta1 <- droppredict(meta1)
  
  
  ##
  ##
  ## (5) Remove results from second meta-analysis (if necessary)
  ##
  ##
  
  if (!common2)
    meta2 <- dropcommon(meta2)
  ##
  if (!random2)
    meta2 <- droprandom(meta2)
  ##
  if (!prediction2)
    meta2 <- droppredict(meta2)
  
  
  ##
  ##
  ## (6) Merge results
  ##
  ##
  
  res <- meta1
  ##
  ncom1 <- length(meta1$TE.common)
  nran1 <- length(meta1$TE.random)
  ncom2 <- length(meta2$TE.common)
  nran2 <- length(meta2$TE.random)
  ##
  ## Individual study weights
  ##
  if (!common1 & common2) {
    res$w.common <- meta2$w.common
    res$text.w.common <- meta2$text.w.common
  }
  else if (!is.null(meta1$w.common) & !is.null(meta2$w.common) & keep.w) {
    res$text.w.common <- c(meta1$text.w.common, meta2$text.w.common)
    ##
    res$w.common <- cbind(meta1$w.common, meta2$w.common)
    colnames(res$w.common) <- res$text.w.common
    rownames(res$w.common) <- meta1$studlab
  }
  ##
  if (is.null(res$w.random) & !is.null(meta2$w.random)) {
    res$w.random <- meta2$w.random
    res$text.w.random <- meta2$text.w.random
  }
  else if (!is.null(res$w.random) & !is.null(meta2$w.random) & keep.w) {
    if (length(res$w.random) != length(meta2$w.random)) {
      warning("Argument 'keep.w' set to FALSE as number of weights differs ",
              "between meta-analyses.",
              call. = FALSE)
      keep.w <- FALSE
    }
    else {
      res$text.w.random <- c(meta1$text.w.random, meta2$text.w.random)
      ##
      res$w.random <- cbind(meta1$w.random, meta2$w.random)
      colnames(res$w.random) <- c(meta1$text.w.random, meta2$text.w.random)
      rownames(res$w.random) <- meta1$studlab
    }
  }
  ##
  res$hetlabel <- c(meta1$hetlabel, meta2$hetlabel)
  ##
  res$method <-
    expandmerge(meta1$method, meta2$method,
                nc1 = ncom1, nc2 = ncom2)
  res$method.random <-
    expandmerge(meta1$method.random, meta2$method.random,
                nr1 = nran1, nr2 = nran2)
  ##
  ## Number of studies
  ##
  if (!inherits(meta1, "metamerge") & (is.limit | is.copas | is.robu)) {
    res$k <- meta1$k
    res$k.all <- meta1$k.all
    res$k.MH <- meta1$k.MH
    res$k.study <- meta1$k.study
    res$k.TE <- meta1$k.TE
  }
  else {
    res$k <-
      expandmerge(meta1$k, meta2$k,
                  nc1 = ncom1, nc2 = ncom2, nr1 = nran1, nr2 = nran2)
    res$k.all <-
      expandmerge(meta1$k.all, meta2$k.all,
                  nc1 = ncom1, nc2 = ncom2, nr1 = nran1, nr2 = nran2)
    res$k.MH <-
      expandmerge(meta1$k.MH, meta2$k.MH,
                  nc1 = ncom1, nc2 = ncom2, nr1 = nran1, nr2 = nran2)
    res$k.study <-
      expandmerge(meta1$k.study, meta2$k.study,
                  nc1 = ncom1, nc2 = ncom2, nr1 = nran1, nr2 = nran2)
    res$k.TE <-
      expandmerge(meta1$k.TE, meta2$k.TE,
                  nc1 = ncom1, nc2 = ncom2, nr1 = nran1, nr2 = nran2)
    res$k0 <-
      expandmerge(meta1$k0, meta2$k0,
                  nc1 = ncom1, nc2 = ncom2, nr1 = nran1, nr2 = nran2)
  }
  ##
  ## Common effect model
  ##
  res$TE.common <- c(meta1$TE.common, meta2$TE.common)
  res$seTE.common <- c(meta1$seTE.common, meta2$seTE.common)
  res$statistic.common <- c(meta1$statistic.common, meta2$statistic.common)
  res$pval.common <- c(meta1$pval.common, meta2$pval.common)
  res$lower.common <- c(meta1$lower.common, meta2$lower.common)
  res$upper.common <- c(meta1$upper.common, meta2$upper.common)
  res$zval.common <- c(meta1$zval.common, meta2$zval.common)
  res$text.common <- c(meta1$text.common, meta2$text.common)
  ##
  ## Random effects model
  ##
  res$TE.random <- c(meta1$TE.random, meta2$TE.random)
  res$seTE.random <- c(meta1$seTE.random, meta2$seTE.random)
  res$statistic.random <- c(meta1$statistic.random, meta2$statistic.random)
  res$pval.random <- c(meta1$pval.random, meta2$pval.random)
  res$method.random.ci <-
    replaceNULL(c(meta1$method.random.ci, meta2$method.random.ci), "")
  res$df.random <- c(meta1$df.random, meta2$df.random)
  res$lower.random <- c(meta1$lower.random, meta2$lower.random)
  res$upper.random <- c(meta1$upper.random, meta2$upper.random)
  res$zval.random <- c(meta1$zval.random, meta2$zval.random)
  res$seTE.classic <- c(meta1$seTE.classic, meta2$seTE.classic)
  res$adhoc.hakn.ci <- c(meta1$adhoc.hakn.ci, meta2$adhoc.hakn.ci)
  res$df.hakn <- c(meta1$df.hakn, meta2$df.hakn)
  res$seTE.hakn.ci <- c(meta1$seTE.hakn.ci, meta2$seTE.hakn.ci)
  res$seTE.hakn.adhoc.ci <-
    c(meta1$seTE.hakn.adhoc.ci, meta2$seTE.hakn.adhoc.ci)
  res$df.kero <- c(meta1$df.kero, meta2$df.kero)
  res$seTE.kero <- c(meta1$seTE.kero, meta2$seTE.kero)
  res$text.random <- c(meta1$text.random, meta2$text.random)
  ##
  ## Prediction interval
  ##
  res$method.predict <- c(meta1$method.predict, meta2$method.predict)
  res$seTE.predict <- c(meta1$seTE.predict, meta2$seTE.predict)
  res$df.predict <- c(meta1$df.predict, meta2$df.predict)
  res$lower.predict <- c(meta1$lower.predict, meta2$lower.predict)
  res$upper.predict <- c(meta1$upper.predict, meta2$upper.predict)
  res$adhoc.hakn.pi <- c(meta1$adhoc.hakn.pi, meta2$adhoc.hakn.pi)
  res$seTE.hakn.pi <- c(meta1$seTE.hakn.pi, meta2$seTE.hakn.pi)
  res$seTE.hakn.adhoc.pi <-
    c(meta1$seTE.hakn.adhoc.pi, meta2$seTE.hakn.adhoc.pi)
  res$text.predict <- c(meta1$text.predict, meta2$text.predict)
  ##
  res$prediction.subgroup <-
    any(c(meta1$prediction.subgroup, meta2$prediction.subgroup))
  ##
  ## Heterogeneity statistics
  ##
  if (keep.Q) {
    res$Q <- c(meta1$Q, meta2$Q)
    res$df.Q <- c(meta1$df.Q, meta2$df.Q)
    res$pval.Q <- c(meta1$pval.Q, meta2$pval.Q)
  }
  ##
  if (keep.I2) {
    res$I2 <- c(meta1$I2, meta2$I2)
    res$lower.I2 <- c(meta1$lower.I2, meta2$upper.I2)
    res$upper.I2 <- c(meta1$upper.I2, meta2$upper.I2)
    ##
    res$H <- c(meta1$H, meta2$H)
    res$lower.H <- c(meta1$lower.H, meta2$upper.H)
    res$upper.H <- c(meta1$upper.H, meta2$upper.H)
    ##
    res$Rb <- c(meta1$Rb, meta2$Rb)
    res$lower.Rb <- c(meta1$lower.Rb, meta2$upper.Rb)
    res$upper.Rb <- c(meta1$upper.Rb, meta2$upper.Rb)
  }
  ##
  ## Trim-and-fill method
  ##
  res$left <-
    expandmerge(meta1$left, meta2$left,
                nc1 = ncom1, nc2 = ncom2, nr1 = nran1, nr2 = nran2)
  res$ma.common <-
    expandmerge(meta1$ma.common, meta2$ma.common,
                nc1 = ncom1, nc2 = ncom2, nr1 = nran1, nr2 = nran2)
  res$type <-
    expandmerge(meta1$type, meta2$type,
                nc1 = ncom1, nc2 = ncom2, nr1 = nran1, nr2 = nran2)
  res$n.iter.max <-
    expandmerge(meta1$n.iter.max, meta2$n.iter.max,
                nc1 = ncom1, nc2 = ncom2, nr1 = nran1, nr2 = nran2)
  res$n.iter <-
    expandmerge(meta1$n.iter, meta2$n.iter,
                nc1 = ncom1, nc2 = ncom2, nr1 = nran1, nr2 = nran2)
  ##
  ## Subgroup analyses
  ##
  replaceNULL(label1.subgroup, "")
  replaceNULL(label2.subgroup, "")
  ##
  res$TE.common.w <- c(meta1$TE.common.w, meta2$TE.common.w)
  res$seTE.common.w <- c(meta1$seTE.common.w, meta2$seTE.common.w)
  res$statistic.common.w <-
    c(meta1$statistic.common.w, meta2$statistic.common.w)
  res$pval.common.w <- c(meta1$pval.common.w, meta2$pval.common.w)
  res$lower.common.w <- c(meta1$lower.common.w, meta2$lower.common.w)
  res$upper.common.w <- c(meta1$upper.common.w, meta2$upper.common.w)
  ##
  res$w.common.w <- c(meta1$w.common.w, meta2$w.common.w)
  ##
  res$tau.w <- c(meta1$tau.w, meta2$tau.w)
  ##
  res$Q.w.common <- c(meta1$Q.w.common, meta2$Q.w.common)
  res$pval.Q.w.common <- c(meta1$pval.Q.w.common, meta2$pval.Q.w.common)
  ##
  res$Q.b.common <- c(meta1$Q.b.common, meta2$Q.b.common)
  res$df.Q.b.common <- c(meta1$df.Q.b.common, meta2$df.Q.b.common)
  res$pval.Q.b.common <- c(meta1$pval.Q.b.common, meta2$pval.Q.b.common)
  ##
  res$TE.random.w <- c(meta1$TE.random.w, meta2$TE.random.w)
  res$seTE.random.w <- c(meta1$seTE.random.w, meta2$seTE.random.w)
  res$statistic.random.w <-
    c(meta1$statistic.random.w, meta2$statistic.random.w)
  res$pval.random.w <- c(meta1$pval.random.w, meta2$pval.random.w)
  res$df.random.w <- c(meta1$df.random.w, meta2$df.random.w)
  res$lower.random.w <- c(meta1$lower.random.w, meta2$lower.random.w)
  res$upper.random.w <- c(meta1$upper.random.w, meta2$upper.random.w)
  res$df.hakn.w <- c(meta1$df.hakn.w, meta2$df.hakn.w)
  res$df.kero.w <- c(meta1$df.kero.w, meta2$df.kero.w)
  ##
  res$w.random.w <- c(meta1$w.random.w, meta2$w.random.w)
  ##
  res$Q.w.random <- c(meta1$Q.w.random, meta2$Q.w.random)
  res$pval.Q.w.random <- c(meta1$pval.Q.w.random, meta2$pval.Q.w.random)
  ##
  res$Q.b.random <- c(meta1$Q.b.random, meta2$Q.b.random)
  res$df.Q.b.random <- c(meta1$df.Q.b.random, meta2$df.Q.b.random)
  res$pval.Q.b.random <- c(meta1$pval.Q.b.random, meta2$pval.Q.b.random)
  ##
  res$seTE.predict.w <- c(meta1$seTE.predict.w, meta2$seTE.predict.w)
  res$df.predict.w <- c(meta1$df.predict.w, meta2$df.predict.w)
  res$lower.predict.w <- c(meta1$lower.predict.w, meta2$lower.predict.w)
  res$upper.predict.w <- c(meta1$upper.predict.w, meta2$upper.predict.w)
  
  
  ##
  ##
  ## (7) More settings
  ##
  ##
  if (is.null(taulabel1))
    taulabel1 <- names(meta1$tau)
  if (is.null(taulabel2))
    taulabel2 <- names(meta2$tau)
  ##
  names(meta1$tau) <- names(meta1$tau2) <- taulabel1
  names(meta2$tau) <- names(meta2$tau2) <- taulabel2
  ##
  if (!random1 & random2) {
    res$method.tau <- meta2$method.tau
    res$method.tau.ci <- meta2$method.tau.ci
    res$detail.tau <- meta2$detail.tau
    ##
    res$tau <- meta2$tau
    names(res$tau) <- taulabel2
    res$lower.tau <- meta2$lower.tau
    res$upper.tau <- meta2$upper.tau
    ##
    res$tau2 <- meta2$tau2
    names(res$tau2) <- taulabel2
    res$lower.tau2 <- meta2$lower.tau2
    res$upper.tau2 <- meta2$upper.tau2
    ##
    res$se.tau <- meta2$se.tau
    ##
    res$tau.preset <- meta2$tau.preset
    ##
    res$Q.Cochrane <- meta2$Q.Cochrane
  }
  else if (random1 & random2) {
    if (!inherits(meta1, "metamerge") && meta1$three.level) {
        meta1$tau2 <- sum(meta1$tau2)
        meta1$tau <- sqrt(meta1$tau2)
        meta1$lower.tau <- meta1$upper.tau <-
          meta1$lower.tau2 <- meta1$upper.tau2 <- NA
        meta1$sign.lower.tau <- meta1$sign.upper.tau <- ""
    }
    ##
    if (!inherits(meta2, "metamerge") && meta2$three.level) {
        meta2$tau2 <- sum(meta2$tau2)
        meta2$tau <- sqrt(meta2$tau2)
        meta2$lower.tau <- meta2$upper.tau <-
          meta2$lower.tau2 <- meta2$upper.tau2 <- NA
        meta2$sign.lower.tau <- meta2$sign.upper.tau <- ""
    }
    ##
    res$method.tau <-
      expandmerge(meta1$method.tau, meta2$method.tau,
                  nr1 = nran1, nr2 = nran2)
    res$method.tau.ci <-
      expandmerge(meta1$method.tau.ci, meta2$method.tau.ci,
                  nr1 = nran1, nr2 = nran2)
    ##
    res$tau <- expandmerge(meta1$tau, meta2$tau, nr1 = nran1, nr2 = nran2)
    res$lower.tau <-
      expandmerge(meta1$lower.tau, meta2$lower.tau, nr1 = nran1, nr2 = nran2)
    res$upper.tau <-
      expandmerge(meta1$upper.tau, meta2$upper.tau, nr1 = nran1, nr2 = nran2)
    ##
    res$tau2 <- expandmerge(meta1$tau2, meta2$tau2, nr1 = nran1, nr2 = nran2)
    res$lower.tau2 <- expandmerge(meta1$lower.tau2, meta2$lower.tau2,
                                  nr1 = nran1, nr2 = nran2)
    res$upper.tau2 <- expandmerge(meta1$upper.tau2, meta2$upper.tau2,
                                  nr1 = nran1, nr2 = nran2)
    ##
    res$se.tau <-
      expandmerge(meta1$se.tau, meta2$se.tau, nr1 = nran1, nr2 = nran2)
    ##
    res$detail.tau <-
      expandmerge(meta1$detail.tau, meta2$detail.tau,
                  nr1 = nran1, nr2 = nran2)
    ##
    res$tau.preset <-
      expandmerge(meta1$tau.preset, meta2$tau.preset, nr1 = nran1, nr2 = nran2)
    ##
    res$Q.Cochrane <-
      expandmerge(meta1$Q.Cochrane, meta2$Q.Cochrane, nr1 = nran1, nr2 = nran2)
  }
  ##
  res$sm <- if (meta1$sm != meta2$sm) "" else meta1$sm
  ##
  res$common <- any(common)
  res$random <- any(random)
  res$overall <- any(overall)
  res$overall.hetstat <- overall.hetstat
  res$prediction <- prediction
  ##
  res$backtransf <- backtransf
  ##
  ## Three-level model
  ##
  if (!inherits(meta1, "metamerge") & is.limit)
    res$three.level <- unique(res$three.level)
  else
    res$three.level <-
      expandmerge(meta1$three.level, meta2$three.level,
                  nr1 = nran1, nr2 = nran2)
  ##
  res$rho <-
    expandmerge(meta1$rho, meta2$rho, nr1 = nran1, nr2 = nran2)
  ##
  if (!random1 & random2)
    res$cluster <- meta2$cluster
  ##
  ## Additional arguments from metabin(), metainc(), metaprop()
  ##
  res$incr <- expandmerge(meta1$incr, meta2$incr,
                          ncom1, nran1, ncom2, nran2)
  res$method.incr <-
    expandmerge(meta1$method.incr, meta2$method.incr,
                ncom1, nran1, ncom2, nran2)
  res$sparse <-
    expandmerge(meta1$sparse, meta2$sparse,
                ncom1, nran1, ncom2, nran2)
  ##
  ## Additional arguments from metabin()
  ##
  res$allstudies <-
    expandmerge(meta1$allstudies, meta2$allstudies,
                ncom1, nran1, ncom2, nran2)
  res$doublezeros <-
    expandmerge(meta1$doublezeros, meta2$doublezeros,
                ncom1, nran1, ncom2, nran2)
  res$MH.exact <-
    expandmerge(meta1$MH.exact, meta2$MH.exact,
                ncom1, nran1, ncom2, nran2)
  res$RR.Cochrane <-
    expandmerge(meta1$RR.Cochrane, meta2$RR.Cochrane,
                ncom1, nran1, ncom2, nran2)
  res$model.glmm <-
    expandmerge(meta1$model.glmm, meta2$model.glmm,
                ncom1, nran1, ncom2, nran2)
  res$phi <-
    expandmerge(meta1$phi, meta2$phi,
                ncom1, nran1, ncom2, nran2)
  ##
  ## Additional arguments from metacont()
  ##
  res$pooledvar <- expandmerge(meta1$pooledvar, meta2$pooledvar,
                               ncom1, nran1, ncom2, nran2)
  res$method.smd <- expandmerge(meta1$method.smd, meta2$method.smd,
                                ncom1, nran1, ncom2, nran2)
  res$sd.glass <- expandmerge(meta1$sd.glass, meta2$sd.glass,
                              ncom1, nran1, ncom2, nran2)
  res$exact.smd <- expandmerge(meta1$exact.smd, meta2$exact.smd,
                               ncom1, nran1, ncom2, nran2)
  res$method.mean <- expandmerge(meta1$method.mean, meta2$method.mean,
                                 ncom1, nran1, ncom2, nran2)
  res$method.sd <- expandmerge(meta1$method.sd, meta2$method.sd,
                               ncom1, nran1, ncom2, nran2)
  
  
  ##
  ##
  ## (8) Backward compatibility
  ##
  ##
  res <- backward(res)
  
  
  ##
  ##
  ## (9) Set class
  ##
  ##
  ##  
  class(res) <- c(class(res), "metamerge")
  ##
  if (inherits(meta2, "trimfill"))
    class(res) <- c(class(res), "trimfill")
  ##
  res$call <- match.call()
  ##
  class(res) <- unique(class(res))
  
  
  res
}
