#' Add pooled results from external analysis to meta-analysis
#' 
#' @description
#' Add pooled results from external analysis to an existing
#' meta-analysis object. This is useful, for example, to add results
#' from a Bayesian meta-analysis which is not implemented in R package
#' \bold{meta}.
#' 
#' @param x Meta-analysis object.
#' @param type A character string or vector indicating whether added
#'   results are from common effect, random effects model or
#'   prediction interval. Either \code{"common"}, \code{"random"} or
#'   \code{"prediction"}, can be abbreviated.
#' @param TE Pooled estimate(s).
#' @param lower Lower limit(s) of confidence or prediction interval.
#' @param upper Upper limit(s) of confidence or prediction interval.
#' @param statistic Test statistic(s).
#' @param pval P-value(s).
#' @param text A character string or vector used in printouts and
#'   forest plot to label the added results.
#' @param data An optional data frame containing the new results or an
#'   object of class \code{meta}.
#' @param method.common A character string or vector to describe the
#'   common effect method(s).
#' @param method.random A character string or vector to describe the
#'   random effects method(s).
#' @param method.tau A character string or vector to describe the
#'   estimator(s) of the between-study variance.
#' @param method.random.ci A character string or vector to describe
#'   the method(s) to calculate confidence intervals under the random
#'   effects model.
#' @param method.predict A character string or vector to describe the
#'   method(s) used for prediction intervals.
#' @param transf A logical indicating whether inputs for arguments
#'   \code{TE}, \code{lower} and \code{upper} are already
#'   appropriately transformed to conduct the meta-analysis or on the
#'   original scale. If \code{transf = TRUE} (default), inputs are
#'   expected to be log odds ratios instead of odds ratios for
#'   \code{sm = "OR"} and Fisher's z transformed correlations instead
#'   of correlations for \code{sm = "ZCOR"}, for example.
#' 
#' @details
#' In R package \bold{meta}, objects of class \code{"meta"} contain
#' results of both common effect and random effects
#' meta-analyses. This function enables the user to add the pooled
#' results of an additional analysis to an existing meta-analysis
#' object. This is useful, for example, to add the result of a
#' Bayesian meta-analysis.
#'
#' If argument \code{data} is a meta-analysis object created with R
#' package \bold{meta}, arguments \code{TE}, \code{lower},
#' \code{upper}, \code{statistic} and \code{pval} are ignored as this
#' information is extracted from the meta-analysis.
#'
#' Otherwise, arguments \code{TE}, \code{lower} and \code{upper} have
#' to be provided if \code{type = "common"} or \code{type =
#' "random"}. For \code{type = "prediction"}, only arguments
#' \code{lower} and \code{upper} are mandatory.
#'
#' Note, R function \code{\link{metamerge}} can be used to add
#' meta-analysis results of another meta-analysis object (see
#' \code{\link{meta-object}}).
#' 
#' @return
#' An object of class \code{"meta"} with corresponding generic
#' functions (see \code{\link{meta-object}}).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metagen}}, \code{\link{metamerge}}
#' 
#' @examples
#' data(Fleiss1993bin)
#' 
#' # Common effect and random effects meta-analysis
#' m1 <- metabin(d.asp, n.asp, d.plac, n.plac, data = Fleiss1993bin,
#'   studlab = paste(study, year), sm = "OR")
#'
#' # Naive pooling
#' m2 <- metabin(sum(d.asp), sum(n.asp), sum(d.plac), sum(n.plac),
#'   data = Fleiss1993bin, sm = "OR", text.common = "Naive pooling")
#'
#' # Add results of second meta-analysis from common effect model
#' m12 <- metaadd(m1, data = m2, method.common = "Naive pooling")
#' m12
#'
#' forest(m12)
#' 
#' @export metaadd


metaadd <- function(x, type,
                    TE, lower, upper, statistic = NA, pval = NA,
                    text, data = NULL,
                    ##
                    method.common = "",
                    ##
                    method.random = "",
                    method.tau = "",
                    method.random.ci = "",
                    ##
                    method.predict = "",
                    ##
                    transf = gs("transf")) {
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  
  chkclass(x, "meta")
  res <- updateversion(x)
  ##
  missing.type <- missing(type)
  ##
  missing.TE <- missing(TE)
  missing.lower <- missing(lower)
  missing.upper <- missing(upper)
  missing.statistic <- missing(statistic)
  missing.pval <- missing(pval)
  missing.text <- missing(text)
  ##
  chklogical(transf)
  
  
  ##
  ##
  ## (2) Read data
  ##
  ##
  
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  if (nulldata)
    data <- sfsp
  ##
  ## Catch 'type', 'TE', 'lower', 'upper', 'statistic', and 'pval' from data:
  ##
  if (inherits(data, "meta")) {
    warnmeta <- function(miss, name)
      if (!miss)
        warning("Argument '", name, "' ignored as argument 'data' is a ",
                "meta-analysis object.",
                call. = FALSE)
    ##
    warnmeta(missing.type, "type")
    warnmeta(missing.TE, "TE")
    warnmeta(missing.lower, "lower")
    warnmeta(missing.upper, "upper")
    warnmeta(missing.statistic, "statistic")
    warnmeta(missing.pval, "pval")
    ##
    if (missing(method.common))
      method.common <- data$method
    ##
    if (missing(method.random))
      method.random <- data$method.random
    if (missing(method.tau))
      method.tau <- data$method.tau
    if (missing(method.random.ci))
      method.random.ci <- data$method.random.ci
    ##
    if (missing(method.predict))
      method.predict <- data$method.predict
    ##
    n.com <- length(data$lower.common)
    n.ran <- length(data$lower.random)
    n.prd <- length(data$lower.predict)
    ##
    type <- c(if (data$common)
                rep_len("common", n.com),
              if (data$random)
                rep_len("random", n.ran),
              if (data$prediction)
                rep_len("prediction", n.prd))
    ##
    if (length(type) == 0)
      type <- "common"
    ##
    res$tau <-
      c(if (res$random | !data$random) res$tau,
        if (data$random) data$tau)
    res$lower.tau <-
      c(if (res$random | !data$random) res$lower.tau,
        if (data$random) data$lower.tau)
    res$upper.tau <-
      c(if (res$random | !data$random) res$upper.tau,
        if (data$random) data$upper.tau)
    ##
    res$tau2 <-
      c(if (res$random | !data$random) res$tau2,
        if (data$random) data$tau2)
    res$lower.tau2 <-
      c(if (res$random | !data$random) res$lower.tau2,
        if (data$random) data$lower.tau2)
    res$upper.tau2 <-
      c(if (res$random | !data$random) res$upper.tau2,
        if (data$random) data$upper.tau2)
    ##
    res$method.tau <-
      c(if (res$random | !data$random) res$method.tau,
        if (data$random) data$method.tau)
    ##
    res$method.tau.ci <-
      c(if (res$random | !data$random) res$method.tau.ci,
        if (data$random) data$method.tau.ci)
    ##
    TE <- lower <- upper <- statistic <- pval <-
      vector("numeric", length(type))
    ##
    if (length(data$TE.common) == 1 & n.com > 1)
      data$TE.common <- rep_len(data$TE.common, n.com)
    ##
    if (length(data$TE.random) == 1 & n.ran > 1)
      data$TE.random <- rep_len(data$TE.random, n.ran)
    ##
    text <- vector("character", length(type))
    ##
    if (!missing.text) {
      text <- catch("text", mc, data, sfsp)
      ##
      if (length(text) == 1 & length(type) > 1)
        text <- rep_len(text, length(type))
      else
        chklength(text, length(type),
                  text =
                    paste0("Argument 'text' must be of length ",
                          length(type), "."))
    }
    ##
    j.c <- j.r <- j.p <- 0
    ##
    for (i in seq_along(type)) {
      if (type[i] == "common") {
        j.c <- j.c + 1
        ##
        TE[i] <- data$TE.common[j.c]
        lower[i] <- data$lower.common[j.c]
        upper[i] <- data$upper.common[j.c]
        statistic[i] <- data$statistic.common[j.c]
        pval[i] <- data$pval.common[j.c]
        ##
        if (missing.text)
          text[i] <- data$text.common[j.c]
      }
      else if (type[i] == "random") {
        j.r <- j.r + 1
        ##
        TE[i] <- data$TE.random[j.r]
        lower[i] <- data$lower.random[j.r]
        upper[i] <- data$upper.random[j.r]
        statistic[i] <- data$statistic.random[j.r]
        pval[i] <- data$pval.random[j.r]
        ##
        if (missing.text)
          text[i] <- data$text.random[j.r]
      }
      else if (type[i] == "prediction") {
        j.r <- j.r + 1
        ##
        TE[i] <- NA
        lower[i] <- data$lower.predict[j.p]
        upper[i] <- data$upper.predict[j.p]
        statistic[i] <- NA
        pval[i] <- NA
        ##
        if (missing.text)
          text[i] <- data$text.predict[j.p]
      }
    }
  }
  else {
    type <- catch("type", mc, data, sfsp)
    ##
    if (any(type %in% c("common", "random")))
      if (missing.TE)
        stop("Argument 'TE' must be provided for common effect or ",
             "random effects model.",
             call. = FALSE)
    ##
    if (missing.lower)
      stop("Argument 'lower' must be provided.",
           call. = FALSE)
    ##
    if (missing.upper)
      stop("Argument 'upper' must be provided.",
           call. = FALSE)
    ##
    TE <- catch("TE", mc, data, sfsp)
    lower <- catch("lower", mc, data, sfsp)
    upper <- catch("upper", mc, data, sfsp)
    statistic <- catch("statistic", mc, data, sfsp)
    pval <- catch("pval", mc, data, sfsp)
    text <- catch("text", mc, data, sfsp)
    ##
    type <- setchar(type, c("common", "random", "prediction"))
    ##
    if (length(type) == 1 & length(lower) > 1)
      type <- rep_len(type, length(lower))
    ##
    k.all <- length(type)
    ##
    if (!missing.TE) {
      chknumeric(TE)
      chklength(TE, k.all, name = "type")
    }
    ##
    chknumeric(lower)
    chklength(lower, k.all, name = "type")
    ##
    chknumeric(upper)
    chklength(upper, k.all, name = "type")
    ##
    if (!missing.statistic) {
      chknumeric(statistic)
      chklength(statistic, k.all, name = "type")
    }
    else
      statistic <- rep_len(NA, k.all)
    ##
    if (!missing.pval) {
      chknumeric(pval)
      chklength(pval, k.all, name = "type")
    }
    else
      pval <- rep_len(NA, k.all)
    ##
    if (!missing.text) {
      if (length(text) == 1 & k.all > 1)
        text <- rep_len(text, k.all)
      else
        chklength(text, k.all)
    }
    else
      text <- rep_len("Added result", k.all)
  }

  
  ##
  ##
  ## (3) Transform added results
  ##
  ##
  
  if (!transf) {
    TE <- transf(TE, x$sm, x$func.transf, x$args.transf)
    lower <- transf(lower, x$sm, x$func.transf, x$args.transf)
    upper <- transf(upper, x$sm, x$func.transf, x$args.transf)
  }
  
  
  ##
  ##
  ## (4) Add results
  ##
  ##
  
  j.c <- j.r <- j.p <- 0
  ##
  for (i in seq_along(type)) {
    if (type[i] == "common") {
      j.c <- j.c + 1
      ##
      if (!res$common)
        res$w.common[!is.na(res$w.common)] <- NA
      ##
      res$TE.common <-
        c(if (res$common) res$TE.common, TE[i])
      res$lower.common <-
        c(if (res$common) res$lower.common, lower[i])
      res$upper.common <-
        c(if (res$common) res$upper.common, upper[i])
      res$statistic.common <-
        c(if (res$common) res$statistic.common, statistic[i])
      res$pval.common <-
        c(if (res$common) res$pval.common, pval[i])
      res$text.common <-
        c(if (res$common) res$text.common, text[i])
      ##
      res$method <- c(if (res$common) res$method, method.common[j.c])
      ##
      res$common <- TRUE
      res$overall <- TRUE
    }
    ##
    if (type[i] == "random") {
      j.r <- j.r + 1
      ##
      if (!res$random)
        res$w.random[!is.na(res$w.random)] <- NA
      ##
      res$TE.random <-
        c(if (res$random) res$TE.random, TE[i])
      res$seTE.random <-
        c(if (res$random) res$seTE.random, NA)
      res$lower.random <-
        c(if (res$random) res$lower.random, lower[i])
      res$upper.random <-
        c(if (res$random) res$upper.random, upper[i])
      res$statistic.random <-
        c(if (res$random) res$statistic.random, statistic[i])
      res$pval.random <-
        c(if (res$random) res$pval.random, pval[i])
      res$text.random <-
        c(if (res$random) res$text.random, text[i])
      ##
      res$method.random <-
        c(if (res$random) res$method.random, method.random[j.r])
      res$method.tau <-
        c(if (res$random) res$method.tau, method.tau[j.r])
      res$method.random.ci <-
        c(if (res$random) res$method.random.ci, method.random.ci[j.r])
      res$df.random <-
        c(if (res$random) res$df.random, NA)
      ##
      res$random <- TRUE
      res$overall <- TRUE
    }
    ##
    if (type[i] == "prediction") {
      j.p <- j.p + 1
      ##
      res$lower.predict <-
        c(if (res$prediction) res$lower.predict, lower[i])
      res$upper.predict <-
        c(if (res$prediction) res$upper.predict, upper[i])
      res$text.predict <-
        c(if (res$prediction) res$text.predict, text[i])
      ##
      res$method.predict <-
        c(if (res$prediction) res$method.predict, method.predict[j.p])
      ##
      res$prediction <- TRUE
      res$overall <- TRUE
    }
  }
  ##
  class(res) <- c(class(res), "metaadd")
  
  res
}
