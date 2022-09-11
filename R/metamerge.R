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
#'   plot to label the estimates from the first meta-analysis.
#' @param text.pooled2 A character string used in printouts and forest
#'   plot to label the estimates from the second meta-analysis.
#' @param text.w.pooled1 A character string used to label weights of
#'   the first meta-analysis.
#' @param text.w.pooled2 A character string used to label weights of
#'   the second meta-analysis.
#' @param detail.tau1 A character string used to label estimate of
#'   between-study variance of the first meta-analysis.
#' @param detail.tau2 A character string used to label estimate of
#'   between-study variance of the second meta-analysis.
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
#' \item Hartung-Knapp method (see argument \code{method.random.ci} in
#'   \code{\link{meta-object}}),
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
#' \item{w.common}{Weight of individual studies (first meta-analysis).}
#' \item{w.random}{Weight of individual studies (second
#'   meta-analysis).}
#' \item{TE.common, seTE.common}{Estimated overall treatment effect and
#'   standard error (first meta-analysis).}
#' \item{lower.common, upper.common}{Lower and upper confidence interval
#'   limits (first meta-analysis).}
#' \item{statistic.common, pval.common}{Statistic and p-value for test of
#'   overall treatment effect (first meta-analysis).}
#' \item{TE.random, seTE.random}{Estimated overall treatment effect and
#'   standard error (second meta-analysis).}
#' \item{lower.random, upper.random}{Lower and upper confidence interval
#'   limits (second meta-analysis).}
#' \item{statistic.random, pval.random}{Statistic and p-value for test of
#'   overall treatment effect (second meta-analysis).}
#' \item{lower.predict, upper.predict}{Lower and upper limits of
#'   prediction interval (related to first meta-analysis).}
#' \item{k}{Number of studies combined in first meta-analysis.}
#' \item{Q}{Heterogeneity statistic (first meta-analysis).}
#' \item{df.Q}{Degrees of freedom for heterogeneity statistic (first
#'   meta-analysis).}
#' \item{pval.Q}{P-value of heterogeneity test (first meta-analysis).}
#' \item{tau2}{Between-study variance(s) \eqn{\tau^2} (first and
#'   second meta-analysis).}
#' \item{lower.tau2, upper.tau2}{Lower and upper limit of confidence
#'   interval(s) for \eqn{\tau^2} (first and second meta-analysis).}
#' \item{tau}{Square-root of between-study variance(s) \eqn{\tau}
#'   (first and second meta-analysis).}
#' \item{lower.tau, upper.tau}{Lower and upper limit of confidence
#'   interval(s) for \eqn{\tau} (first and second meta-analysis).}
#' \item{text.common}{Label for the first meta-analysis.}
#' \item{text.random}{Label for the second meta-analysis.}
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
#' # Use Hartung-Knapp method
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
#' # Show results for DerSimonian-Laird and Paule-Mandel estimate of
#' # between-study variance
#' #
#' m3 <- update(m1, method.tau = "PM",
#'   text.random = "Random effects moded (PM)", text.w.random = "PM")
#' #
#' m34 <- metamerge(m2, m3)
#' m34
#'
#' data(Fleiss1993bin)
#' #
#' # Mantel-Haenszel method
#' #
#' m5 <- metabin(d.asp, n.asp, d.plac, n.plac, data = Fleiss1993bin,
#'   studlab = paste(study, year), sm = "OR", random = FALSE,
#'   text.common = "MH method", text.w.common = "MH")
#' #
#' # Peto method
#' #
#' m6 <- update(m5, method = "Peto", text.common = "Peto method",
#'   text.w.common = "Peto")
#' #
#' # Merge results (show individual results for MH method)
#' #
#' m56 <- metamerge(m5, m6)
#' m56
#' forest(m56, digits = 4)
#' #
#' # Merge results (show individual results for Peto method)
#' #
#' m65 <- metamerge(m6, m5)
#' m65
#' 
#' @export metamerge


metamerge <- function(meta1, meta2, pooled1, pooled2,
                      text.pooled1, text.pooled2,
                      text.w.pooled1, text.w.pooled2,
                      detail.tau1, detail.tau2,
                      backtransf) {
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chkclass(meta1, c("meta", "limitmeta", "copas"))
  if (inherits(meta1, "meta"))
    meta1 <- updateversion(meta1)
  else if (inherits(meta1, c("limitmeta", "copas"))) {
    if (!missing(meta2))
      warning("Argument 'meta2' ignored as argument 'meta1' is of class '",
              class(meta1), "'.",
              call. = FALSE)
    if (!missing(pooled2))
      warning("Argument 'pooled2' ignored as argument 'meta1' is of class '",
              class(meta1), "'.",
              call. = FALSE)
    meta2 <- meta1
    meta1 <- updateversion(meta1$x)
    pooled2 <- "random"
  }
  ##
  chkclass(meta2, c("meta", "limitmeta", "copas", "robu"))
  if (inherits(meta2, "meta"))
    meta2 <- updateversion(meta2)
  ##
  if (inherits(meta1, "metamerge"))
    stop("Argument 'meta1' already of class \"metameta\".",
         call. = FALSE)
  if (inherits(meta2, "metamerge"))
    stop("Argument 'meta2' already of class \"metameta\".",
         call. = FALSE)
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
  if (!missing(text.pooled1))
    chkchar(text.pooled1, length = 1)
  if (!missing(text.pooled2))
    chkchar(text.pooled2, length = 1)
  if (!missing(text.w.pooled1))
    chkchar(text.w.pooled1, length = 1)
  if (!missing(text.w.pooled2))
    chkchar(text.w.pooled2, length = 1)
  if (!missing(detail.tau1))
    chkchar(detail.tau1, length = 1)
  if (!missing(detail.tau2))
    chkchar(detail.tau2, length = 1)
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
  if (inherits(meta1, "metabin")) {
    if ((meta1$sm != meta2$sm) &
        any(c(meta1$sm, meta2$sm) %in% c("RD", "ASD")))
      stop("Summary measures do not fit.",
           call. = FALSE)
  }
  ##
  ## Check original data (if available)
  ##
  if (!is.null(meta1$data) & !is.null(meta2$data)) {
    ##
    if (nrow(meta1$data) != nrow(meta2$data))
      stop("Meta-analyses based on different data sets.",
           call. = FALSE)
    ##
    if (inherits(meta1, "metabin")) {
      if (any(meta1$data$.event.e != meta2$data$.event.e) |
          any(meta1$data$.n.e != meta2$data$.n.e) |
          any(meta1$data$.event.c != meta2$data$.event.c) |
          any(meta1$data$.n.c != meta2$data$.n.c))
        stop("Meta-analyses have different data.",
             call. = FALSE)
    }
  }
  ##
  ## Check subgroup levels
  ##
  if (!is.null(meta1$subgroup.levels) &
      !is.null(meta2$subgroup.levels)) {
    ##
    if (length(meta1$subgroup.levels) !=
        length(meta2$subgroup.levels))
      stop("Meta-analyses have different number of subgroups.",
           call. = FALSE)
    ##
    if (any(meta1$subgroup.levels != meta2$subgroup.levels))
      stop("Meta-analyses based on different subgroup-analyses.",
           call. = FALSE)
  }
  
  
  ##
  ##
  ## (2) Some assignments for Copas, Limit or Robu
  ##
  ##
  if (is.copas | is.limit | is.robu) {
    meta2$method <- "Inverse"
    ##
    meta2$method.random.ci <- "classic"
    meta2$adhoc.hakn.ci <- ""
    meta2$df.random <- Inf
  }
  if (is.copas) {
    meta2$method.tau <- "ML"
    meta2$detail.tau <- "copas"
    ##
    if (!missing(text.pooled2))
      meta2$text.random <- text.pooled2
    else
      meta2$text.random <- "Copas selection model"
    ##
    if (!missing(text.w.pooled2))
      meta2$text.w.random <- text.w.pooled2
    else
      meta2$text.w.random <- "Copas"
  }
  else if (is.limit) {
    meta2$method.tau <- meta1$method.tau
    meta2$detail.tau <- "limit"
    ##
    if (!missing(text.pooled2))
      meta2$text.random <- text.pooled2
    else
      meta2$text.random <- "Limit meta-analysis"
    ##
    if (!missing(text.w.pooled2))
      meta2$text.w.random <- text.w.pooled2
    else
      meta2$text.w.random <- "limit"
  }
  else if (is.robu) {
    meta2$method.tau <- "DL"
    meta2$detail.tau <- "RVE"
  }
  ##
  if (is.copas | is.limit) {
    ##
    meta2$TE.random <- meta2$TE.adjust
    meta2$seTE.random <- meta2$seTE.adjust
    meta2$lower.random <- meta2$lower.adjust
    meta2$upper.random <- meta2$upper.adjust
    meta2$statistic.random <- meta2$statistic.adjust
    meta2$pval.random <- meta2$pval.adjust
    ##
    meta2$w.random <- rep(0, length(meta2$w.random))
  }
  ##
  if (is.robu) {
    if (!missing(text.pooled2))
      meta2$text.random <- text.pooled2
    else
      meta2$text.random <- "RVE model"
    ##
    if (!missing(text.w.pooled2))
      meta2$text.w.random <- text.w.pooled2
    else
      meta2$text.w.random <- "RVE"
    ##
    meta2$TE.random <- meta2$reg_table$b.r[1]
    meta2$seTE.random <- meta2$reg_table$SE[1]
    meta2$lower.random <- meta2$reg_table$CI.L[1]
    meta2$upper.random <- meta2$reg_table$CI.U[1]
    meta2$statistic.random <- meta2$reg_table$t[1]
    meta2$pval.random <- meta2$reg_table$prob[1]
    ##
    meta2$w.random <- meta2$data.full$r.weights
  }
  
  
  ##
  ##
  ## (3) Some more assignments
  ##
  ##
  meta1$detail.tau <- replaceNULL(meta1$detail.tau, "")
  ##
  if (!missing(detail.tau1))
    meta1$detail.tau <- detail.tau1
  if (!missing(detail.tau2))
    meta2$detail.tau <- detail.tau2
  ##
  meta1$method.tau <- replaceNULL(meta1$method.tau, "")
  meta1$method.tau.ci <- replaceNULL(meta1$method.tau.ci, "")
  ##
  meta2$method.tau <- replaceNULL(meta2$method.tau, "")
  meta2$method.tau.ci <- replaceNULL(meta2$method.tau.ci, "")
  ##
  meta1$method.random.ci <- replaceNULL(meta1$method.random.ci, "")
  meta2$method.random.ci <- replaceNULL(meta2$method.random.ci, "")
  ##
  meta1$df.random <- replaceNULL(meta1$df.random, NA)
  meta2$df.random <- replaceNULL(meta2$df.random, NA)
  
  
  ##
  ##
  ## (4) Remove results from first meta-analysis
  ##
  ##
  res <- meta1
  ##
  if (pooled1 == "random") {
    res$w.common <- NULL
    res$TE.common <- NULL
    res$seTE.common <- NULL
    res$statistic.common <- NULL
    res$pval.common <- NULL
    res$lower.common <- NULL
    res$upper.common <- NULL
    res$zval.common <- NULL
    res$text.common <- NULL
    ##
    if (!is.null(meta1$subgroup)) {
      res$TE.common.w <- NULL
      res$seTE.common.w <- NULL
      res$statistic.common.w <- NULL
      res$pval.common.w <- NULL
      res$lower.common.w <- NULL
      res$upper.common.w <- NULL
      res$w.common.w <- NULL
      ##
      res$Q.w.common <- NULL
      res$pval.Q.w.common <- NULL
      ##
      res$Q.b.common <- NULL
      res$pval.Q.b.common <- NULL
    }
  }
  ##
  else if (pooled1 == "common") {
    res$w.random <- NULL
    res$TE.random <- NULL
    res$seTE.random <- NULL
    res$statistic.random <- NULL
    res$pval.random <- NULL
    res$method.random.ci <- NULL
    res$df.random <- NULL
    res$lower.random <- NULL
    res$upper.random <- NULL
    res$zval.random <- NULL
    res$seTE.classic <- NULL
    res$adhoc.hakn.ci <- NULL
    res$df.hakn <- NULL
    res$seTE.hakn.ci <- NULL
    res$seTE.hakn.adhoc.ci <- NULL
    res$df.kero <- NULL
    res$seTE.kero <- NULL
    res$text.random <- NULL
    ##
    res$method.predict <- NULL
    res$adhoc.hakn.pi <- NULL
    res$seTE.predict <- NULL
    res$df.predict <- NULL
    res$lower.predict <- NULL
    res$upper.predict <- NULL
    res$seTE.hakn.pi <- NULL
    res$seTE.hakn.adhoc.pi <- NULL
    #      
    if (!is.null(meta1$subgroup)) {
      res$TE.random.w <- NULL
      res$seTE.random.w <- NULL
      res$statistic.random.w <- NULL
      res$pval.random.w <- NULL
      res$df.random.w <- NULL
      res$lower.random.w <- NULL
      res$upper.random.w <- NULL
      res$w.random.w <- NULL
      res$df.hakn.w <- NULL
      res$df.kero.w <- NULL
      ##
      res$seTE.predict.w <- NULL
      res$df.predict.w <- NULL
      res$lower.predict.w <- NULL
      res$upper.predict.w <- NULL
      ##
      res$Q.w.random <- NULL
      res$pval.Q.w.random <- NULL
      ##
      res$Q.b.random <- NULL
      res$pval.Q.b.random <- NULL
    }
  }
  
  
  ##
  ##
  ## (5) Remove results from second meta-analysis
  ##
  ##
  if (pooled2 == "random") {
    meta2$w.common <- NULL
    meta2$TE.common <- NULL
    meta2$seTE.common <- NULL
    meta2$statistic.common <- NULL
    meta2$pval.common <- NULL
    meta2$lower.common <- NULL
    meta2$upper.common <- NULL
    meta2$zval.common <- NULL
    meta2$text.common <- NULL
    ##
    if (!is.null(meta1$subgroup)) {
      meta2$TE.common.w <- NULL
      meta2$seTE.common.w <- NULL
      meta2$statistic.common.w <- NULL
      meta2$pval.common.w <- NULL
      meta2$lower.common.w <- NULL
      meta2$upper.common.w <- NULL
      meta2$w.common.w <- NULL
      ##
      meta2$Q.w.common <- NULL
      meta2$pval.Q.w.common <- NULL
      ##
      meta2$Q.b.common <- NULL
      meta2$pval.Q.b.common <- NULL
    }
  }
  ##
  else if (pooled2 == "common") {
    meta2$w.random <- NULL
    meta2$TE.random <- NULL
    meta2$seTE.random <- NULL
    meta2$statistic.random <- NULL
    meta2$pval.random <- NULL
    meta2$method.random.ci <- NULL
    meta2$df.random <- NULL
    meta2$lower.random <- NULL
    meta2$upper.random <- NULL
    meta2$zval.random <- NULL
    meta2$seTE.classic <- NULL
    meta2$adhoc.hakn.ci <- NULL
    meta2$df.hakn <- NULL
    meta2$seTE.hakn.ci <- NULL
    meta2$seTE.hakn.adhoc.ci <- NULL
    meta2$df.kero <- NULL
    meta2$seTE.kero <- NULL
    meta2$text.random <- NULL
    ##
    meta2$method.predict <- NULL
    meta2$adhoc.hakn.pi <- NULL
    meta2$seTE.predict <- NULL
    meta2$df.predict <- NULL
    meta2$lower.predict <- NULL
    meta2$upper.predict <- NULL
    meta2$seTE.hakn.pi <- NULL
    meta2$seTE.hakn.adhoc.pi <- NULL
    #      
    if (!is.null(meta2$subgroup)) {
      meta2$TE.random.w <- NULL
      meta2$seTE.random.w <- NULL
      meta2$statistic.random.w <- NULL
      meta2$pval.random.w <- NULL
      meta2$df.random.w <- NULL
      meta2$lower.random.w <- NULL
      meta2$upper.random.w <- NULL
      meta2$w.random.w <- NULL
      meta2$df.hakn.w <- NULL
      meta2$df.kero.w <- NULL
      ##
      meta2$seTE.predict.w <- NULL
      meta2$df.predict.w <- NULL
      meta2$lower.predict.w <- NULL
      meta2$upper.predict.w <- NULL
      ##
      meta2$Q.w.random <- NULL
      meta2$pval.Q.w.random <- NULL
      ##
      meta2$Q.b.random <- NULL
      meta2$pval.Q.b.random <- NULL
    }
  }
  
  
  ##
  ##
  ## (6) Merge results
  ##
  ##
  if (is.null(res$w.common) & !is.null(meta2$w.common))
    res$w.common <- meta2$w.common
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
      res$pval.Q.b.common <- meta2$pval.Q.b.common
    }
  }
  ##
  if (is.null(res$w.random) & !is.null(meta2$w.random))
    res$w.random <- meta2$w.random
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
      ##
      res$Q.b.random <- meta2$Q.b.random
      res$pval.Q.b.random <- meta2$pval.Q.b.random
    }
  }
  
  
  ##
  ##
  ## (7) More settings
  ##
  ##
  if (!is.null(meta2$method))
    res$method <- if (meta1$method == meta2$method) meta1$method else ""
  else
    res$method <- meta1$method
  ##
  if (!is.null(meta1$Q.Cochrane) & !is.null(meta2$Q.Cochrane))
    res$Q.Cochrane <-
      if (meta1$Q.Cochrane == meta2$Q.Cochrane) meta1$Q.Cochrane else FALSE
  ##
  if (pooled1 == "common" & pooled2 == "common") {
    res$overall.hetstat <- FALSE
    ##
    res$method.tau <- ""
    res$method.tau.ci <- ""
    res$tau <- NA
    res$lower.tau <- NA
    res$upper.tau <- NA
    res$tau2 <- NA
    res$lower.tau2 <- NA
    res$upper.tau2 <- NA
    res$se.tau <- NA
  }
  ##
  else if (pooled1 == "common" & pooled2 %in% c("both", "random")) {
    res$method.tau <- meta2$method.tau
    res$method.tau.ci <- meta2$method.tau.ci    
    ##
    if (is.copas) {
      res$tau <- meta2$tau.adjust
      res$lower.tau <- NA
      res$upper.tau <- NA
      res$tau2 <- meta2$tau.adjust^2
      res$lower.tau2 <- NA
      res$upper.tau2 <- NA
      res$se.tau <- NA
    }
    else if (is.robu) {
      res$tau <- sqrt(meta2$mod_info$tau.sq)
      res$lower.tau <- NA
      res$upper.tau <- NA
      res$tau2 <- meta2$mod_info$tau.sq
      res$lower.tau2 <- NA
      res$upper.tau2 <- NA
      res$se.tau <- NA
    }
    else {
      res$tau <- meta2$tau
      res$lower.tau <- meta2$lower.tau
      res$upper.tau <- meta2$upper.tau
      res$tau2 <- meta2$tau2
      res$lower.tau2 <- meta2$lower.tau2
      res$upper.tau2 <- meta2$upper.tau2
      res$se.tau <- meta2$se.tau
    }
  }
  else if (pooled1 %in% c("both", "random") &
           pooled2 %in% c("both", "random")) {
    if (is.copas) {
      if (res$method.tau != "ML" & res$detail.tau == "") {
        res$detail.tau <- res$method.tau
        res$method.tau <- ""
      }
      ##
      res$tau <- c(res$tau, meta2$tau.adjust)
      res$lower.tau <- c(res$lower.tau, NA)
      res$upper.tau <- c(res$upper.tau, NA)
      res$tau2 <- c(res$tau2, meta2$tau.adjust^2)
      res$lower.tau2 <- c(res$lower.tau2, NA)
      res$upper.tau2 <- c(res$upper.tau2, NA)
      res$se.tau <- c(res$se.tau, NA)
      ##
      res$detail.tau <- c(res$detail.tau, meta2$detail.tau)
    }
    else if (is.robu) {
      res$detail.tau <- res$method.tau
      res$method.tau <- ""
      ##
      res$tau <- c(res$tau, sqrt(meta2$mod_info$tau.sq))
      res$lower.tau <- c(res$lower.tau, NA)
      res$upper.tau <- c(res$upper.tau, NA)
      res$tau2 <- c(res$tau2, meta2$mod_info$tau.sq)
      res$lower.tau2 <- c(res$lower.tau2, NA)
      res$upper.tau2 <- c(res$upper.tau2, NA)
      res$se.tau <- c(res$se.tau, NA)
      ##
      res$detail.tau <- c(res$detail.tau, meta2$detail.tau)
    }
    else {
      ##
      if (
        any(meta1$tau != meta2$tau) |
        (any(!is.na(meta1$lower.tau) & !is.na(meta2$lower.tau)) &&
         any(meta1$lower.tau != meta2$lower.tau))) {
        ##
        if (meta1$method.tau != meta2$method.tau) {
          if (res$detail.tau == "")
            res$detail.tau <- meta1$method.tau
          if (meta2$detail.tau == "")
            meta2$detail.tau <- meta2$method.tau
          res$method.tau <- ""
        }
        ##
        if (meta1$method.tau.ci != meta2$method.tau.ci) {
          if (res$detail.tau == "")
            res$detail.tau <- meta1$method.tau.ci
          if (meta2$detail.tau == "")
            meta2$detail.tau <- meta2$method.tau.ci
          res$method.tau.ci <- ""
        }
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
      }
    }
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
  ## Backward compatibility
  ##
  res$comb.fixed <- res$fixed <- res$common
  res$comb.random <- res$random
  ##
  class(res) <- c(class(res), "metamerge")
  ##
  res
}
