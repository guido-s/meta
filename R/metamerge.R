#' Merge pooled results of two meta-analyses
#' 
#' @description
#' This function can be used to merge pooled results of two
#' meta-analyses into a single meta-analysis object. This is, for
#' example, useful to produce a forest plot of a random-effects
#' meta-analysis with different estimates of the between-study
#' variance \eqn{\tau^2}.
#' 
#' @param meta1 First meta-analysis object (see Details).
#' @param meta2 Second meta-analysis object (see Details).
#' @param pooled1 A character string indicating whether results of
#'   common effect or random effects model should be considered for
#'   first meta-analysis. Either \code{"both"}, \code{"common"} or
#'   \code{"random"}, can be abbreviated.
#' @param pooled2 A character string indicating whether results of
#'   common effect or random effects model should be considered for
#'   second meta-analysis. Either \code{"both"}, \code{"common"} or
#'   \code{"random"}, can be abbreviated.
#' @param text.pooled1 A character string used in printouts and forest
#'   plot to label the results from the first meta-analysis.
#' @param text.pooled2 A character string used in printouts and forest
#'   plot to label the results from the second meta-analysis.
#' @param text.w.pooled1 A character string used to label weights of
#'   the first meta-analysis.
#' @param text.w.pooled2 A character string used to label weights of
#'   the second meta-analysis.
#' @param label1 A character string used to label estimate of
#'   between-study variance and heterogeneity statistics of the first
#'   meta-analysis.
#' @param label2 A character string used to label estimate of
#'   between-study variance and heterogeneity statistics of the second
#'   meta-analysis.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=TRUE} (default), results for \code{sm="OR"} are
#'   printed as odds ratios rather than log odds ratios, for example.
#' 
#' @details
#' In R package \bold{meta}, objects of class \code{"meta"} contain
#' results of both a common effect and random effects
#' meta-analysis. This function enables the user to keep the results
#' of one of these models and to add results from a second
#' meta-analysis or a sensitivity analysis.
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
#' The first argument must be an object created by a meta-analysis
#' function, e.g., \code{\link{metagen}} or \code{\link{metabin}}. It
#' is also possible to provide an object created with
#' \code{\link[metasens]{limitmeta}} or
#' \code{\link[metasens]{copas}}. In this case, arguments \code{meta2}
#' and \code{pooled2} will be ignored.
#'
#' The second meta-analysis could be an object created by a
#' meta-analysis function or with \code{\link{trimfill}},
#' \code{\link[metasens]{limitmeta}}, \code{\link[metasens]{copas}},
#' or \code{\link[robumeta]{robu}}.
#'
#' The created meta-analysis object only contains the study results
#' from the first meta-analysis which are shown in printouts and
#' forest plots. This only makes a difference for meta-analysis
#' methods where individual study results differ, e.g.,
#' Mantel-Haenszel and Peto method for binary outcomes (see
#' \code{\link{metabin}}).
#'
#' R function \code{\link{metabind}} can be used to print and plot the
#' results of more than two meta-analyses, however, without showing
#' individual study results.
#' 
#' @return
#' An object of class \code{"meta"} and \code{"metamerge"} with
#' corresponding generic functions (see \code{\link{meta-object}}).
#' 
#' The following list elements have a different meaning:
#' \item{TE, seTE, studlab}{Treatment estimate, standard error, and
#'   study labels (first meta-analyis).}
#' \item{lower, upper}{Lower and upper confidence interval limits for
#'   individual studies (first meta-analysis).}
#' \item{statistic, pval}{Statistic and p-value for test of treatment
#'   effect for individual studies (first meta-analysis.}
#' \item{w.common}{Weights of first common effect meta-analysis.}
#' \item{w.random}{Weights of first random effects meta-analysis.}
#' \item{k}{Number of studies combined in first meta-analysis.}
#'
#' Furthermore, meta-analysis results of common effect or random
#' effects model are taken from first meta-analysis if only random
#' effects or common effects models are selected from both
#' meta-analyses (arguments \code{pooled1} and \code{pooled2}).
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{metagen}}, \code{\link{metabind}}
#' 
#' @examples
#' data(Fleiss1993cont)
#' #
#' m1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, sm = "MD",
#'   text.random = "Random effects model (REML)", text.w.random = "DL")
#' #
#' # Use DerSimonian-Laird estimator of tau2
#' #
#' m2 <- update(m1, method.tau = "DL", common = FALSE,
#'   text.random = "Random effects model (DL)", text.w.random = "DL")
#' #
#' # Merge results of the two meta-analyses
#' #
#' m12 <- metamerge(m1, m2)
#' m12
#' forest(m12, rightcols = c("effect", "ci", "w.common"))
#'
#' # Show in addition the results for the Paule-Mandel estimate of
#' # between-study variance
#' #
#' m3 <- update(m1, method.tau = "PM",
#'   text.random = "Random effects moded (PM)", text.w.random = "PM")
#' #
#' m123 <- metamerge(m12, m3, pooled2 = "random")
#' m123
#'
#' data(Fleiss1993bin)
#' #
#' # Mantel-Haenszel method
#' #
#' m4 <- metabin(d.asp, n.asp, d.plac, n.plac, data = Fleiss1993bin,
#'   studlab = paste(study, year), sm = "OR", random = FALSE,
#'   text.common = "MH method", text.w.common = "MH")
#' #
#' # Peto method
#' #
#' m5 <- update(m4, method = "Peto", text.common = "Peto method",
#'   text.w.common = "Peto")
#' #
#' # Merge results (show individual results for MH method)
#' #
#' m45 <- metamerge(m4, m5)
#' summary(m45)
#' forest(m45, digits = 4)
#' #
#' # Merge results (show individual results for Peto method)
#' #
#' m54 <- metamerge(m5, m4)
#' summary(m54)
#' forest(m54)
#' 
#' @export metamerge


metamerge <- function(meta1, meta2, pooled1, pooled2,
                      text.pooled1, text.pooled2,
                      text.w.pooled1, text.w.pooled2,
                      label1, label2,
                      backtransf) {
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  
  if (missing(meta1))
    stop("Argument 'meta1' must be provided.",
         call. = FALSE)
  if (missing(meta2)) {
    warning("Argument 'meta2' not provided.",
            call. = FALSE)
    return(meta1)
  }
  ##
  if (inherits(meta1, "copas") |
      inherits(meta1, "limitmeta") |
      inherits(meta1, "robu"))
    stop("Argument 'meta1' cannot be of class '",
         class(meta1), "' (use argument 'meta2').",
         call. = FALSE)
  ##
  chkclass(meta1, "meta")
  meta1 <- updateversion(meta1)
  ##
  chkclass(meta2, c("meta", "limitmeta", "copas", "robu"))
  ##
  if (inherits(meta2, "meta"))
    meta2 <- updateversion(meta2)
  ##
  is.copas <- inherits(meta2, "copas")
  is.limit <- inherits(meta2, "limitmeta")
  is.robu <- inherits(meta2, "robu")
  ##
  if (!missing(pooled1)) {
    pooled1 <- setchar(pooled1, c("both", "common", "random", "fixed"))
    pooled1[pooled1 == "fixed"] <- "common"
  }
  else
    pooled1 <-
      if (meta1$common & meta1$random)
        "both"
      else if (meta1$common & !meta1$random)
        "common"
      else
        "random"
  ##
  if (!missing(pooled2)) {
    pooled2 <- setchar(pooled2, c("both", "common", "random", "fixed"))
    pooled2[pooled2 == "fixed"] <- "common"
  }
  else {
    if (is.copas | is.limit | is.robu)
      pooled2 <- "random"
    else
      pooled2 <-
        if (meta2$common & meta2$random)
          "both"
        else if (meta2$common & !meta2$random)
          "common"
        else
          "random"
  }
  ##
  missing.text.pooled1 <- missing(text.pooled1)
  if (!missing.text.pooled1)
    chkchar(text.pooled1, length = 1)
  ##
  missing.text.pooled2 <- missing(text.pooled2)
  if (!missing.text.pooled2)
    chkchar(text.pooled2, length = 1)
  ##
  missing.text.w.pooled1 <- missing(text.w.pooled1)
  if (!missing.text.w.pooled1)
    chkchar(text.w.pooled1, length = 1)
  ##
  missing.text.w.pooled2 <- missing(text.w.pooled2)
  if (!missing.text.w.pooled2)
    chkchar(text.w.pooled2, length = 1)
  ##
  missing.label1 <- missing(label1)
  if (!missing.label1)
    chkchar(label1, length = 1)
  ##
  missing.label2 <- missing(label2)
  if (!missing.label2)
    chkchar(label2, length = 1)
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
  samesm(meta1, meta2)
  ##
  ## Check original data (if available)
  ##
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
  
  meta1 <- updateobj(meta1, text.pooled1, missing.text.pooled1,
                     text.w.pooled1, missing.text.w.pooled1,
                     label1, missing.label1)
  ##
  meta2 <- updateobj(meta2, text.pooled2, missing.text.pooled2,
                     text.w.pooled2, missing.text.w.pooled2,
                     label2, missing.label2)
  ##
  if (is.null(meta1$label))
    meta1$label <- "meta2"
  ##
  if (is.null(meta2$label))
    meta2$label <- "meta2"
  
  
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
  meta1$method.random.ci <- replaceNULL(meta1$method.random.ci, "")
  meta2$method.random.ci <- replaceNULL(meta2$method.random.ci, "")
  ##
  meta1$df.random <- replaceNULL(meta1$df.random, NA)
  meta2$df.random <- replaceNULL(meta2$df.random, NA)
  
  
  ##
  ##
  ## (4) Remove results from first meta-analysis (if necessary)
  ##
  ##
  
  res <- meta1
  ##
  if (length(res$TE.random) == 1 &&
      length(res$TE.random) != length(res$seTE.random)) {
    res$TE.random <- rep(res$TE.random, length(res$seTE.random))
    names(res$TE.random) <- names(res$seTE.random)
  }
  ##
  if (pooled1 == "random" & pooled2 %in% c("common", "both"))
    res <- dropcommon(res)
  ##
  if (pooled1 == "common" & pooled2 %in% c("random", "both")) {
    res <- droprandom(res)
    res <- droppredict(res)
  }
  
  
  ##
  ##
  ## (5) Remove results from second meta-analysis (if necessary)
  ##
  ##
  
  if (length(meta2$TE.random) == 1 &&
      length(meta2$TE.random) != length(meta2$seTE.random)) {
    meta2$TE.random <- rep(meta2$TE.random, length(meta2$seTE.random))
    names(meta2$TE.random) <- names(meta2$seTE.random)
  }
  ##
  if (pooled2 == "random")
    meta2 <- dropcommon(meta2, !is.null(meta1$subgroup))
  ##
  if (pooled2 == "common")
    meta2 <- droprandom(meta2, !is.null(meta1$subgroup))
  
  
  ##
  ##
  ## (6) Merge results
  ##
  ##
  
  if (is.null(res$w.common) & !is.null(meta2$w.common))
    res$w.common <- meta2$w.common
  ##
  if (is.null(res$w.random) & !is.null(meta2$w.random))
    res$w.random <- meta2$w.random
  ##
  res$label <- c(res$label, meta2$label)
  ##
  res$method.tau.ci <- c(res$method.tau.ci, meta2$method.tau.ci)
  ##
  res$TE.common <- c(res$TE.common, meta2$TE.common)
  res$seTE.common <- c(res$seTE.common, meta2$seTE.common)
  res$statistic.common <- c(res$statistic.common, meta2$statistic.common)
  res$pval.common <- c(res$pval.common, meta2$pval.common)
  res$lower.common <- c(res$lower.common, meta2$lower.common)
  res$upper.common <- c(res$upper.common, meta2$upper.common)
  res$zval.common <- c(res$zval.common, meta2$zval.common)
  res$text.common <- c(res$text.common, meta2$text.common)
  ##
  res$Q <- c(res$Q, meta2$Q)
  res$df.Q <- c(res$df.Q, meta2$df.Q)
  res$pval.Q <- c(res$pval.Q, meta2$pval.Q)
  ##
  res$I2 <- c(res$I2, meta2$I2)
  res$lower.I2 <- c(res$lower.I2, meta2$lower.I2)
  res$upper.I2 <- c(res$upper.I2, meta2$upper.I2)
  ##
  res$H <- c(res$H, meta2$H)
  res$lower.H <- c(res$lower.H, meta2$lower.H)
  res$upper.H <- c(res$upper.H, meta2$upper.H)
  ##
  res$Rb <- c(res$Rb, meta2$Rb)
  res$lower.Rb <- c(res$lower.Rb, meta2$lower.Rb)
  res$upper.Rb <- c(res$upper.Rb, meta2$upper.Rb)
  ##
  if (!is.null(meta1$subgroup) | !is.null(meta2$subgroup)) {
    if (is.null(res$w.common.w) & !is.null(meta2$w.common.w))
      res$w.common.w <- meta2$w.common.w
    ##
    res$TE.common.w <- c(res$TE.common.w, meta2$TE.common.w)
    res$seTE.common.w <- c(res$seTE.common.w, meta2$seTE.common.w)
    res$statistic.common.w <-
      c(res$statistic.common.w, meta2$statistic.common.w)
    res$pval.common.w <- c(res$pval.common.w, meta2$pval.common.w)
    res$lower.common.w <- c(res$lower.common.w, meta2$lower.common.w)
    res$upper.common.w <- c(res$upper.common.w, meta2$upper.common.w)
    ##
    if (is.null(meta1$subgroup)) {
      res$Q.w.common <- meta2$Q.w.common
      res$pval.Q.w.common <- meta2$pval.Q.w.common
      ##
      res$Q.b.common <- meta2$Q.b.common
      res$df.Q.b.common <- meta2$df.Q.b.common
      res$pval.Q.b.common <- meta2$pval.Q.b.common
    }
  }
  ##
  res$TE.random <- c(res$TE.random, meta2$TE.random)
  res$seTE.random <- c(res$seTE.random, meta2$seTE.random)
  res$statistic.random <- c(res$statistic.random, meta2$statistic.random)
  res$pval.random <- c(res$pval.random, meta2$pval.random)
  res$method.random.ci <- c(res$method.random.ci, meta2$method.random.ci)
  res$df.random <- c(res$df.random, meta2$df.random)
  res$lower.random <- c(res$lower.random, meta2$lower.random)
  res$upper.random <- c(res$upper.random, meta2$upper.random)
  res$zval.random <- c(res$zval.random, meta2$zval.random)
  res$seTE.classic <- c(res$seTE.classic, meta2$seTE.classic)
  res$adhoc.hakn.ci <- c(res$adhoc.hakn.ci, meta2$adhoc.hakn.ci)
  res$df.hakn <- c(res$df.hakn, meta2$df.hakn)
  res$seTE.hakn.ci <- c(res$seTE.hakn.ci, meta2$seTE.hakn.ci)
  res$seTE.hakn.adhoc.ci <-
    c(res$seTE.hakn.adhoc.ci, meta2$seTE.hakn.adhoc.ci)
  res$df.kero <- c(res$df.kero, meta2$df.kero)
  res$seTE.kero <- c(res$seTE.kero, meta2$seTE.kero)
  res$text.random <- c(res$text.random, meta2$text.random)
  ##
  res$method.predict <- c(res$method.predict, meta2$method.predict)
  res$adhoc.hakn.pi <- c(res$adhoc.hakn.pi, meta2$adhoc.hakn.pi)
  res$seTE.predict <- c(res$seTE.predict, meta2$seTE.predict)
  res$df.predict <- c(res$df.predict, meta2$df.predict)
  res$lower.predict <- c(res$lower.predict, meta2$lower.predict)
  res$upper.predict <- c(res$upper.predict, meta2$upper.predict)
  res$seTE.hakn.pi <- c(res$seTE.hakn.pi, meta2$seTE.hakn.pi)
  res$seTE.hakn.adhoc.pi <-
    c(res$seTE.hakn.adhoc.pi, meta2$seTE.hakn.adhoc.pi)
  res$text.predict <- c(res$text.predict, meta2$text.predict)
  ##      
  if (!is.null(meta1$subgroup) | !is.null(meta2$subgroup)) {
    if (is.null(res$w.random.w) & !is.null(meta2$w.random.w))
      res$w.random.w <- meta2$w.random.w
    ##
    res$TE.random.w <- c(res$TE.random.w, meta2$TE.random.w)
    res$seTE.random.w <- c(res$seTE.random.w, meta2$seTE.random.w)
    res$statistic.random.w <-
      c(res$statistic.random.w, meta2$statistic.random.w)
    res$pval.random.w <- c(res$pval.random.w, meta2$pval.random.w)
    res$df.random.w <- c(res$df.random.w, meta2$df.random.w)
    res$lower.random.w <- c(res$lower.random.w, meta2$lower.random.w)
    res$upper.random.w <- c(res$upper.random.w, meta2$upper.random.w)
    res$w.random.w <- c(res$w.random.w, meta2$w.random.w)
    res$df.hakn.w <- c(res$df.hakn.w, meta2$df.hakn.w)
    res$df.kero.w <- c(res$df.kero.w, meta2$df.kero.w)
    ##
    res$seTE.predict.w <- c(res$seTE.predict.w, meta2$seTE.predict.w)
    res$df.predict.w <- c(res$df.predict.w, meta2$df.predict.w)
    res$lower.predict.w <- c(res$lower.predict.w, meta2$lower.predict.w)
    res$upper.predict.w <- c(res$upper.predict.w, meta2$upper.predict.w)
    ##
    if (is.null(meta1$subgroup)) {
      res$Q.w.random <- meta2$Q.w.random
      res$pval.Q.w.random <- meta2$pval.Q.w.random
    }
    ##
    res$Q.b.random <- c(res$Q.b.random, meta2$Q.b.random)
    res$df.Q.b.random <- c(res$df.Q.b.random, meta2$df.Q.b.random)
    res$pval.Q.b.random <- c(res$pval.Q.b.random, meta2$pval.Q.b.random)
  }
  ##
  ## Remove duplicates
  ##
  equalQ <- duplicated(res$Q) & duplicated(res$df.Q) &
    duplicated(res$pval.Q)
  if (any(equalQ)) {
    res$Q <- res$Q[!equalQ]
    res$df.Q <- res$df.Q[!equalQ]
    res$pval.Q <- res$pval.Q[!equalQ]
  }
  ##
  equalI2 <- duplicated(res$I2) & duplicated(res$lower.I2) &
    duplicated(res$upper.I2)
  if (any(equalI2)) {
    res$I2 <- res$I2[!equalI2]
    res$lower.I2 <- res$lower.I2[!equalI2]
    res$upper.I2 <- res$upper.I2[!equalI2]
  }
  ##
  equalH <- duplicated(res$H) & duplicated(res$lower.H) &
    duplicated(res$upper.H)
  if (any(equalH)) {
    res$H <- res$H[!equalH]
    res$lower.H <- res$lower.H[!equalH]
    res$upper.H <- res$upper.H[!equalH]
  }
  ##
  equalRb <- duplicated(res$Rb) & duplicated(res$lower.Rb) &
    duplicated(res$upper.Rb)
  if (any(equalRb)) {
    res$Rb <- res$Rb[!equalRb]
    res$lower.Rb <- res$lower.Rb[!equalRb]
    res$upper.Rb <- res$upper.Rb[!equalH]
  }
  
  
  ##
  ##
  ## (7) More settings
  ##
  ##
  
  res$method <- c(res$method, meta2$method)
  ##
  if (!is.null(meta1$Q.Cochrane) & !is.null(meta2$Q.Cochrane))
    res$Q.Cochrane <-
      if (meta1$Q.Cochrane == meta2$Q.Cochrane) meta1$Q.Cochrane else FALSE
  ##
  if (pooled1 == "common" & pooled2 %in% c("both", "random")) {
    res$method.tau <- meta2$method.tau
    res$method.tau.ci <- meta2$method.tau.ci    
    ##
    res$tau <- meta2$tau
    res$lower.tau <- meta2$lower.tau
    res$upper.tau <- meta2$upper.tau
    res$tau2 <- meta2$tau2
    res$lower.tau2 <- meta2$lower.tau2
    res$upper.tau2 <- meta2$upper.tau2
    res$se.tau <- meta2$se.tau
  }
  else {
    if (pooled1 == "common" & pooled2 == "common")
      res$overall.hetstat <- FALSE
    ##
    res$tau <- c(res$tau, meta2$tau)
    res$lower.tau <- c(res$lower.tau, meta2$lower.tau)
    res$upper.tau <- c(res$upper.tau, meta2$upper.tau)
    res$tau2 <- c(res$tau2, meta2$tau2)
    res$lower.tau2 <- c(res$lower.tau2, meta2$lower.tau2)
    res$upper.tau2 <- c(res$upper.tau2, meta2$upper.tau2)
    res$se.tau <- c(res$se.tau, meta2$se.tau)
    ##
    res$detail.tau <- c(res$detail.tau, meta2$detail.tau)
    res$method.tau <- c(res$method.tau, meta2$method.tau)
  }
  ##
  res$common <-
    pooled1 %in% c("both", "common") |
    pooled2 %in% c("both", "common")
  ##
  res$random <-
    pooled1 %in% c("both", "random") |
    pooled2 %in% c("both", "random")
  ##
  res$backtransf <- backtransf
  ##
  res$pooled1 <- pooled1
  res$pooled2 <- pooled2

  
  ##
  ##
  ## (8) Replace NULL with NA
  ##
  ##
  
  res$w.common <-
    replaceNULL(res$w.common, rep_len(0, length(res$studlab)))
  res$TE.common <- replaceNULL(res$TE.common)
  res$seTE.common <- replaceNULL(res$seTE.common)
  res$statistic.common <- replaceNULL(res$statistic.common)
  res$pval.common <- replaceNULL(res$pval.common)
  res$lower.common <- replaceNULL(res$lower.common)
  res$upper.common <- replaceNULL(res$upper.common)
  res$zval.common <- replaceNULL(res$zval.common)
  res$text.common <- replaceNULL(res$text.common, "")
  ##
  res$w.random <-
    replaceNULL(res$w.random, rep_len(0, length(res$studlab)))
  res$TE.random <- replaceNULL(res$TE.random)
  res$seTE.random <- replaceNULL(res$seTE.random)
  res$statistic.random <- replaceNULL(res$statistic.random)
  res$pval.random <- replaceNULL(res$pval.random)
  res$method.random.ci <- replaceNULL(res$method.random.ci, "")
  res$df.random <- replaceNULL(res$df.random)
  res$lower.random <- replaceNULL(res$lower.random)
  res$upper.random <- replaceNULL(res$upper.random)
  res$zval.random <- replaceNULL(res$zval.random)
  res$seTE.classic <- replaceNULL(res$seTE.classic)
  res$adhoc.hakn.ci <- replaceNULL(res$adhoc.hakn.ci, "")
  res$df.hakn <- replaceNULL(res$df.hakn)
  res$seTE.hakn.ci <- replaceNULL(res$seTE.hakn.ci)
  res$seTE.hakn.adhoc.ci <- replaceNULL(res$seTE.hakn.adhoc.ci)
  res$df.kero <- replaceNULL(res$df.kero)
  res$seTE.kero <- replaceNULL(res$seTE.kero)
  res$text.random <- replaceNULL(res$text.random, "")
  ##
  res$method.predict <- replaceNULL(res$method.predict, "")
  res$adhoc.hakn.pi <- replaceNULL(res$adhoc.hakn.pi, "")
  res$seTE.predict <- replaceNULL(res$seTE.predict)
  res$df.predict <- replaceNULL(res$df.predict)
  res$lower.predict <- replaceNULL(res$lower.predict)
  res$upper.predict <- replaceNULL(res$upper.predict)
  res$seTE.hakn.pi <- replaceNULL(res$seTE.hakn.pi)
  res$seTE.hakn.adhoc.pi <- replaceNULL(res$seTE.hakn.adhoc.pi)
  res$text.predict <- replaceNULL(res$text.predict, "")

  
  ##
  ##
  ## (9) Additional setting for trim-and-fill method
  ##
  ##
  
  if (inherits(meta2, "trimfill")) {
    class(res) <- c(class(res), "trimfill")
    ##
    if (is.null(res$k0)) {
      if (res$k == meta2$k - meta2$k0 &
          res$k.study == meta2$k.study - meta2$k0 &
          res$k.all == meta2$k.all - meta2$k0 &
          res$k.TE == meta2$k.TE - meta2$k0) {
        res$k <- meta2$k
        res$k.study <- meta2$k.study
        res$k.all <- meta2$k.all
        res$k.TE <- meta2$k.TE
        res$k0 <- meta2$k0
      }
      else {
        res$k <- c(res$k, meta2$k)
        res$k.study <- c(res$k.study, meta2$k.study)
        res$k.all <- c(res$k.all, meta2$k.all)
        res$k.TE <- c(res$k.TE, meta2$k.TE)
        res$k0 <- c(NA, meta2$k0)
      }
    }
    else if (res$k != meta2$k |
             res$k.study != meta2$k.study |
             res$k.all != meta2$k.all |
             res$k.TE != meta2$k.TE |
             res$k0 != meta2$k0) {
      res$k <- c(res$k, meta2$k)
      res$k.study <- c(res$k.study, meta2$k.study)
      res$k.all <- c(res$k.all, meta2$k.all)
      res$k.TE <- c(res$k.TE, meta2$k.TE)
      res$k0 <- c(res$k0, meta2$k0)
    }
  }

  
  ##
  ##
  ## (10) Backward compatibility
  ##
  ##
  
  res$comb.fixed <- res$fixed <- res$common
  res$comb.random <- res$random
  
  
  class(res) <- c(class(res), "metamerge")
  ##
  res
}
