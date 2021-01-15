#' Merge pooled results of two meta-analyses
#' 
#' @description
#' This function can be used to merge pooled results of two
#' meta-analyses into a single meta-analysis object. This is, for
#' example, useful to produce a forest plot of a random-effects
#' meta-analysis with and without using the Hartung-Knapp method.
#' 
#' @param meta1 First meta-analysis object (of class \code{"meta"}).
#' @param meta2 Second meta-analysis object (see Details).
#' @param pooled1 A character string indicating whether results of
#'   fixed effect or random effects model should be considered for
#'   first meta-analysis. Either \code{"fixed"} or \code{"random"},
#'   can be abbreviated.
#' @param pooled2 A character string indicating whether results of
#'   fixed effect or random effects model should be considered for
#'   second meta-analysis. Either \code{"fixed"} or \code{"random"},
#'   can be abbreviated.
#' @param text.pooled1 A character string used in printouts and forest
#'   plot to label the estimate from the first meta-analysis.
#' @param text.pooled2 A character string used in printouts and forest
#'   plot to label the estimate from the second meta-analysis.
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
#' results of both a fixed effect and random effects
#' meta-analysis. This function enables the user to keep the results
#' of one of these models and to add results from a second
#' meta-analysis or a sensitivity analysis.
#'
#' Applications of this function include printing and plotting results
#' of the fixed effect or random effects meta-analysis and the
#' \itemize{
#' \item Hartung-Knapp method (see argument \code{hakn} in
#'   \code{\link{metagen}}),
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
#' function, e.g., \code{\link{metagen}} or \code{\link{metabin}}. The
#' second meta-analysis could also be an object created with
#' \code{\link{trimfill}}, \code{\link[metasens]{limitmeta}},
#' \code{\link[metasens]{copas}}, or \code{\link[robumeta]{robu}}.
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
#' corresponding \code{print}, \code{summary}, and \code{forest}
#' functions. The following list elements have a different meaning:
#' \item{TE, seTE, studlab}{Treatment estimate, standard error, and
#'   study labels (first meta-analyis).}
#' \item{lower, upper}{Lower and upper confidence interval limits for
#'   individual studies (first meta-analysis).}
#' \item{statistic, pval}{Statistic and p-value for test of treatment
#'   effect for individual studies (first meta-analysis.}
#' \item{w.fixed}{Weight of individual studies (first meta-analysis).}
#' \item{w.random}{Weight of individual studies (second
#'   meta-analysis).}
#' \item{TE.fixed, seTE.fixed}{Estimated overall treatment effect and
#'   standard error (first meta-analysis).}
#' \item{lower.fixed, upper.fixed}{Lower and upper confidence interval
#'   limits (first meta-analysis).}
#' \item{statistic.fixed, pval.fixed}{Statistic and p-value for test of
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
#' \item{text.fixed}{Label for the first meta-analysis.}
#' \item{text.random}{Label for the second meta-analysis.}
#'
#' See \code{\link{metagen}} for information on other list
#' elements.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{metagen}}, \code{\link{metabind}}
#' 
#' @examples
#' data(Fleiss1993cont)
#' #
#' m1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'                data = Fleiss1993cont, sm = "MD",
#'                comb.fixed = FALSE,
#'                text.random = "Classic random effects",
#'                text.w.random = "RE")
#' #
#' # Use Hartung-Knapp method
#' #
#' m2 <- update(m1, hakn = TRUE,
#'              text.random = "Hartung-Knapp method",
#'              text.w.random = "HK")
#' #
#' # Merge results of the two meta-analyses
#' #
#' m12 <- metamerge(m1, m2)
#' m12
#' forest(m12, rightcols = c("effect", "ci", "w.fixed"))
#'
#' # Show results for DerSimonian-Laird and REML estimate of
#' # between-study variance
#' #
#' m3 <- update(m1,
#'              text.random = "Random effects moded (DL)",
#'              text.w.random = "DL")
#' m4 <- update(m1, method.tau = "REML",
#'              text.random = "Random effects moded (REML)",
#'              text.w.random = "REML")
#' #
#' m34 <- metamerge(m3, m4)
#' m34
#'
#' data(Fleiss1993bin)
#' #
#' # Mantel-Haenszel method
#' #
#' m5 <- metabin(d.asp, n.asp, d.plac, n.plac, data = Fleiss1993bin,
#'               studlab = paste(study, year),
#'               sm = "OR", comb.random = FALSE,
#'               text.fixed = "MH method", text.w.fixed = "MH")
#' #
#' # Peto method
#' #
#' m6 <- update(m5, method = "Peto", text.fixed = "Peto method",
#'              text.w.fixed = "Peto")
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
  
  
  chkclass(meta1, "meta")
  chkclass(meta2, c("meta", "limitmeta", "copas", "robu"))
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
  if (!missing(pooled1))
    pooled1 <- setchar(pooled1, c("fixed", "random"))
  else
    pooled1 <- ifelse(meta1$comb.random, "random", "fixed")
  ##
  if (!missing(pooled2))
    pooled2 <- setchar(pooled2, c("fixed", "random"))
  else {
    if (is.copas | is.limit | is.robu)
      pooled2 <- "random"
    else
      pooled2 <- ifelse(meta2$comb.random, "random", "fixed")
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
  if (is.copas)
    meta2$detail.tau <- "copas"
  else if (is.limit)
    meta2$detail.tau <- "limit"
  else if (is.robu)
    meta2$detail.tau <- "RVE"
  ##
  if (is.null(meta1$detail.tau))
    meta1$detail.tau <- ""
  
  
  ##
  ## Check original data
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
  ## Check summary measures
  ##
  if (inherits(meta1, "metabin")) {
    if ((meta1$sm != meta2$sm) &
        any(c(meta1$sm, meta2$sm) %in% c("RD", "ASD")))
      stop("Summary measures do not fit.",
           call. = FALSE)
  }
  
  
  ##
  ## Some assignments
  ##
  if (!missing(detail.tau1))
    meta1$detail.tau <- detail.tau1
  if (!missing(detail.tau2))
    meta2$detail.tau <- detail.tau2
  
  
  ##
  ## Result of first meta-analysis is saved in list elements for fixed
  ## effect model
  ##
  res <- meta1
  ##
  if (pooled1 == "random") {
    if (!missing(text.pooled1))
      res$text.fixed <- text.pooled1
    else
      res$text.fixed <- meta1$text.random
    ##
    if (!missing(text.w.pooled1))
      res$text.w.fixed <- text.w.pooled1
    else
      res$text.w.fixed <- meta1$text.w.random
    ##
    res$detail.tau <- meta1$detail.tau
    ##
    res$TE.fixed <- meta1$TE.random
    res$seTE.fixed <- meta1$seTE.random
    res$lower.fixed <- meta1$lower.random
    res$upper.fixed <- meta1$upper.random
    res$statistic.fixed <- meta1$statistic.random
    res$pval.fixed <- meta1$pval.random
    res$w.fixed <- meta1$w.random
    ##
    if (!is.null(meta1$byvar)) {
      res$TE.fixed.w <- meta1$TE.random.w
      res$seTE.fixed.w <- meta1$seTE.random.w
      res$lower.fixed.w <- meta1$lower.random.w
      res$upper.fixed.w <- meta1$upper.random.w
      res$statistic.fixed.w <- meta1$statistic.random.w
      res$pval.fixed.w <- meta1$pval.random.w
      res$w.fixed.w <- meta1$w.random.w
      ##
      res$Q.w.fixed <- meta1$Q.w.random
      res$pval.Q.w.fixed <- meta1$pval.Q.w.random
      ##
      res$Q.b.fixed <- meta1$Q.b.random
      res$pval.Q.b.fixed <- meta1$pval.Q.b.random
    }
  }
  
  
  ##
  ## Merge results of second meta-analysis with first meta-analysis
  ##
  if (is.copas | is.limit) {
    if (!missing(text.pooled2))
      res$text.random <- text.pooled2
    else
      res$text.random <-
        if (is.limit) "Limit meta-analysis" else "Copas selection model"
    ##
    if (!missing(text.w.pooled2))
      res$text.w.random <- text.w.pooled2
    else
      res$text.w.random <-
        if (is.limit) "limit" else "Copas"
    ##
    res$TE.random <- meta2$TE.adjust
    res$seTE.random <- meta2$seTE.adjust
    res$lower.random <- meta2$lower.adjust
    res$upper.random <- meta2$upper.adjust
    res$statistic.random <- meta2$statistic.adjust
    res$pval.random <- meta2$pval.adjust
    ##
    res$w.random <- rep(0, length(res$w.random))
  }
  else if (is.robu) {
    if (!missing(text.pooled2))
      res$text.random <- text.pooled2
    else
      res$text.random <- "RVE model"
    ##
    if (!missing(text.w.pooled2))
      res$text.w.random <- text.w.pooled2
    else
      res$text.w.random <- "RVE"
    ##
    res$TE.random <- meta2$reg_table$b.r[1]
    res$seTE.random <- meta2$reg_table$SE[1]
    res$lower.random <- meta2$reg_table$CI.L[1]
    res$upper.random <- meta2$reg_table$CI.U[1]
    res$statistic.random <- meta2$reg_table$t[1]
    res$pval.random <- meta2$reg_table$prob[1]
    ##
    res$w.random <- meta2$data.full$r.weights
  }
  else if (pooled2 == "fixed") {
    if (!missing(text.pooled2))
      res$text.random <- text.pooled2
    else {
      if (inherits(meta2, "trimfill"))
        res$text.random <-
          paste(meta2$text.fixed, "(trim-and-fill)")
      else
        res$text.random <- meta2$text.fixed
    }
    ##
    if (!missing(text.w.pooled2))
      res$text.w.random <- text.w.pooled2
    else
      res$text.w.random <- meta2$text.w.fixed
    ##
    res$TE.random <- meta2$TE.fixed
    res$seTE.random <- meta2$seTE.fixed
    res$lower.random <- meta2$lower.fixed
    res$upper.random <- meta2$upper.fixed
    res$statistic.random <- meta2$statistic.fixed
    res$pval.random <- meta2$pval.fixed
    ##
    if (!inherits(meta1, "trimfill") & inherits(meta2, "trimfill"))
      res$w.random <- meta2$w.fixed[seq_along(res$w.random)]
    else if (inherits(meta1, "trimfill") & !inherits(meta2, "trimfill")) {
      res$w.random[res$w.random != 0] <- 0
      res$w.random[seq_along(meta2$w.fixed)] <-
        meta2$w.fixed
    }
    else
      res$w.random <- meta2$w.fixed
    ##
    if (!is.null(meta2$byvar)) {
      res$TE.random.w <- meta2$TE.fixed.w
      res$seTE.random.w <- meta2$seTE.fixed.w
      res$lower.random.w <- meta2$lower.fixed.w
      res$upper.random.w <- meta2$upper.fixed.w
      res$statistic.random.w <- meta2$statistic.fixed.w
      res$pval.random.w <- meta2$pval.fixed.w
      res$w.random.w <- meta2$w.fixed.w
      ##
      res$Q.w.random <- meta2$Q.w.fixed
      res$pval.Q.w.random <- meta2$pval.Q.w.fixed
      ##
      res$Q.b.random <- meta2$Q.b.fixed
      res$pval.Q.b.random <- meta2$pval.Q.b.fixed
    }
  }
  else {
    if (!missing(text.pooled2))
      res$text.random <- text.pooled2
    else {
      if (inherits(meta2, "trimfill"))
        res$text.random <-
          paste(meta2$text.random, "(trim-and-fill)")
      else
        res$text.random <- meta2$text.random
    }
    ##
    if (!missing(text.w.pooled2))
      res$text.w.random <- text.w.pooled2
    else
      res$text.w.random <- meta2$text.w.random
    ##
    res$TE.random <- meta2$TE.random
    res$seTE.random <- meta2$seTE.random
    res$lower.random <- meta2$lower.random
    res$upper.random <- meta2$upper.random
    res$statistic.random <- meta2$statistic.random
    res$pval.random <- meta2$pval.random
    ##
    if (!inherits(meta1, "trimfill") & inherits(meta2, "trimfill"))
      res$w.random <- meta2$w.random[seq_along(res$w.random)]
    else if (inherits(meta1, "trimfill") & !inherits(meta2, "trimfill")) {
      res$w.random[res$w.random != 0] <- 0
      res$w.random[seq_along(meta2$w.random)] <-
        meta2$w.random
    }
    else
      res$w.random <- meta2$w.random
    ##
    if (!is.null(meta2$byvar)) {
      res$TE.random.w <- meta2$TE.random.w
      res$seTE.random.w <- meta2$seTE.random.w
      res$lower.random.w <- meta2$lower.random.w
      res$upper.random.w <- meta2$upper.random.w
      res$statistic.random.w <- meta2$statistic.random.w
      res$pval.random.w <- meta2$pval.random.w
      res$w.random.w <- meta2$w.random.w
      ##
      res$Q.w.random <- meta2$Q.w.random
      res$pval.Q.w.random <- meta2$pval.Q.w.random
      ##
      res$Q.b.random <- meta2$Q.b.random
      res$pval.Q.b.random <- meta2$pval.Q.b.random
    }
  }
  
  
  ##
  ## Additional settings
  ##
  if (is.copas | is.limit | is.robu)
    meta2$method <- "Inverse"
  ##
  if (!is.null(meta2$method))
    res$method <- if (meta1$method == meta2$method) meta1$method else ""
  else
    res$method <- meta1$method
  ##
  if (is.copas)
    meta2$method.tau <- "ML"
  else if (is.limit)
    meta2$method.tau <- meta1$method.tau
  else if (is.robu)
    meta2$method.tau <- "DL"
  ##
  if (is.null(meta2$method.tau))
    meta2$method.tau <- ""
  ##
  if (is.null(meta2$method.tau.ci))
    meta2$method.tau.ci <- ""
  ##
  if (!is.null(meta2$hakn))
    res$hakn <-
      (pooled1 == "random" & meta1$hakn) |
      (pooled2 == "random" & meta2$hakn)
  else
    res$hakn <-
      (pooled1 == "random" & meta1$hakn)
  ##
  if (!is.null(meta1$Q.Cochrane) & !is.null(meta2$Q.Cochrane))
    res$Q.Cochrane <-
      if (meta1$Q.Cochrane == meta2$Q.Cochrane) meta1$Q.Cochrane else FALSE
  ##
  if (pooled1 == "fixed" & pooled2 == "fixed") {
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
  if (pooled1 == "fixed" & pooled2 == "random") {
    res$method.tau <- meta2$method.tau
    res$method.tau.ci <- meta2$method.tau.ci    
    ##
    if (is.limit) {
      res$tau <- meta2$tau
      res$lower.tau <- NA
      res$upper.tau <- NA
      res$tau2 <- meta2$tau^2
      res$lower.tau2 <- NA
      res$upper.tau2 <- NA
      res$se.tau <- NA
    }
    else if (is.copas) {
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
  ##
  if (pooled1 == "random" & pooled2 == "random") {
    if (is.limit) {
      res$tau <- c(res$tau, meta2$tau)
      res$lower.tau <- c(res$lower.tau, NA)
      res$upper.tau <- c(res$upper.tau, NA)
      res$tau2 <- c(res$tau2, meta2$tau^2)
      res$lower.tau2 <- c(res$lower.tau2, NA)
      res$upper.tau2 <- c(res$upper.tau2, NA)
      res$se.tau <- c(res$se.tau, NA)
      ##
      res$detail.tau <- c(res$detail.tau, meta2$detail.tau)
    }
    else if (is.copas) {
      if (res$method.tau != "ML" & res$detail.tau == "") {
        res$detail.tau <- res$method.tau
        res$method.tau <- ""
      }
      ##
      res$tau <- c(res$tau, meta2$tau.adjust)
      res$lower.tau <- c(res$lower.tau, NA)
      res$upper.tau <- c(res$upper.tau, NA)
      res$tau2 <- c(res$tau, meta2$tau.adjust^2)
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
    else if (
           any(meta1$tau != meta2$tau) |
           any(meta1$lower.tau != meta2$lower.tau)) {
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
  ##
  res$comb.fixed <- res$comb.random <- TRUE
  res$backtransf <- backtransf
  ##
  res$pooled1 <- pooled1
  res$pooled2 <- pooled2
  ##
  class(res) <- c(class(res), "metamerge")
  ##
  res
}
