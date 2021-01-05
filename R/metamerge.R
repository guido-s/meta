#' Merge two meta-analysis objects
#' 
#' @description
#' 
#' This function can be used to merge two meta-analysis objects.
#' 
#' @param meta1 First meta-analysis object.
#' @param meta2 Second meta-analysis object.
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
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=TRUE} (default), results for \code{sm="OR"} are
#'   printed as odds ratios rather than log odds ratios, for example.
#' 
#' @details
#' This function can be used to merge two meta-analysis objects and
#' is, for example, useful to print or plot results for the classic
#' random effects meta-analysis and the Hartung-Knapp method.
#' 
#' @return
#' An object of class \code{"meta"} with corresponding \code{print},
#' \code{summary}, and \code{forest} functions. See
#' \code{\link{metagen}} for more information on list elements.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{metagen}}, \code{\link{metabind}}
#' 
#' @examples
#' data(Fleiss1993cont)
#' 
#' m1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'                data = Fleiss1993cont, sm = "MD",
#'                comb.fixed = FALSE,
#'                text.random = "Classic random effects",
#'                text.w.random = "RE")
#'
#' # Use Hartung-Knapp method
#' #
#' m2 <- update(m1, hakn = TRUE,
#'              text.random = "Hartung-Knapp method",
#'              text.w.random = "HK")
#'
#' # Merge results of the two meta-analyses
#' #
#' m12 <- metamerge(m1, m2)
#' m12
#' forest(m12, rightcols = c("effect", "ci", "w.fixed"))
#'
#'
#' data(Fleiss1993bin)
#' 
#' # Mantel-Haenszel method
#' #
#' m3 <- metabin(d.asp, n.asp, d.plac, n.plac, data = Fleiss1993bin,
#'               studlab = paste(study, year),
#'               sm = "OR", comb.random = FALSE,
#'               text.fixed = "MH method", text.w.fixed = "MH")
#'
#' # Peto method
#' #
#' m4 <- update(m3, method = "Peto", text.fixed = "Peto method",
#'              text.w.fixed = "Peto")
#' 
#' # Merge results (show individual results for MH method)
#' #
#' m34 <- metamerge(m3, m4)
#' m34
#' forest(m34, digits = 4)
#' 
#' # Merge results (show individual results for Peto method)
#' #
#' m43 <- metamerge(m4, m3)
#' m43
#' 
#' @export metamerge


metamerge <- function(meta1, meta2, pooled1, pooled2,
                      text.pooled1, text.pooled2,
                      text.w.pooled1, text.w.pooled2,
                      backtransf) {


  replace.NULL <- function(x, val = NA) {
    if (is.null(x))
      res <- val
    else
      res <- x
    ##
    res
  }
  
  
  chkclass(meta1, "meta")
  chkclass(meta2, c("meta", "limitmeta", "copas", "robu"))
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
    else
      res$text.random <- meta2$text.fixed
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
    else
      res$text.random <- meta2$text.random
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
           any(meta1$lower.tau != meta2$.lower.tau)) {
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
    ##
    ##if (!is.null(meta2$method.tau))
    if (FALSE)
      if (meta1$method.tau != meta2$method.tau)
        warning("Between-study variance estimate from first meta-analysis ",
                "shown in printouts.", call. = FALSE)
  }
  ##
  res$comb.fixed <- res$comb.random <- TRUE
  res$backtransf <- backtransf
  ##
  res$pooled1 <- pooled1
  res$pooled2 <- pooled2
  
  
  res
}
