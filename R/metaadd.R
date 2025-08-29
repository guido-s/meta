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
#'   prediction interval. Either \code{"common"}, \code{"random"},
#'   \code{"prediction"}, or \code{"tau2"} can be abbreviated.
#' @param TE Pooled estimate(s) or between-study variance.
#' @param lower Lower limit(s) of confidence or prediction interval.
#' @param upper Upper limit(s) of confidence or prediction interval.
#' @param statistic Test statistic(s).
#' @param pval P-value(s).
#' @param df Degrees of freedom for confidence or prediction intervals.s
#' @param se Standard error(s).
#' @param method A character string or vector to describe the
#'   method used to get the pooled estimate(s), prediction interval(s) or
#'   between-study variance(s).
#' @param method.ci A character string or vector to describe the
#'   method used to get the confidence or prediction interval.
#' @param text A character string or vector used in printouts and
#'   forest plot to label the added results.
#' @param data An optional data frame containing the new results or an
#'   object of class \code{meta}.
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
#' m12 <- metaadd(m1, data = m2)
#' m12
#'
#' forest(m12)
#' 
#' @export metaadd


metaadd <- function(x, type = NULL,
                    #
                    TE = NA, lower = NA, upper = NA,
                    #
                    statistic = NA, pval = NA, df = NA, se = NA,
                    #
                    method = "", method.ci = "",
                    #
                    text = "Added result",
                    #
                    data = NULL,
                    #
                    transf = gs("transf")) {
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  
  chkclass(x, "meta")
  res <- updateversion(x)
  #
  if (res$k < 2)
    stop("Meta-analysis object provided in argument 'x' must combine ",
         "at least two studies.",
         call. = FALSE)
  #
  missing.type <- missing(type)
  ##
  missing.TE <- missing(TE)
  missing.lower <- missing(lower)
  missing.upper <- missing(upper)
  missing.statistic <- missing(statistic)
  missing.pval <- missing(pval)
  missing.df <- missing(df)
  missing.se <- missing(se)
  #
  missing.method <- missing(method)
  missing.method.ci <- missing(method.ci)
  #
  if (missing.type & missing.TE & missing.lower & missing.upper &
      missing.statistic & missing.pval & missing.df & missing.se)
    return(x)
  #
  missing.text <- missing(text)
  #
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
  #
  # Catch 'type', 'TE', 'lower', 'upper', 'statistic', 'pval', 'df', 'se',
  # 'method' and 'method.ci' from data:
  #
  if (inherits(data, "meta")) {
    if (inherits(data, "metaadd"))
      stop("Argument 'data' cannot be of class 'metaadd'.",
           call. = FALSE)
    #
    is.meta <- TRUE
    #
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
    warnmeta(missing.df, "df")
    warnmeta(missing.se, "se")
    #
    warnmeta(missing.method, "method")
    warnmeta(missing.method.ci, "method.ci")
    
    dat.c <- dat.r <- dat.t <- dat.p <- NULL
    #
    if (data$common | !(data$random | data$prediction)) {
      if (any(data$method == "GLMM") & is.null(res$model.glmm))
        res$model.glmm <- data$model.glmm
      #
      dat.c <- data.frame(type = "common",
                          #
                          TE = data$TE.common,
                          se = data$seTE.common,
                          lower = data$lower.common,
                          upper = data$upper.common,
                          statistic = data$statistic.common,
                          pval = data$pval.common,
                          df = NA,
                          #
                          method = data$method,
                          method.ci = data$method.common.ci,
                          #
                          text = data$text.common,
                          phi = NA)
    }
    #
    if (data$random) {
      if (any(data$method.random == "GLMM") & is.null(res$model.glmm))
        res$model.glmm <- data$model.glmm
      #
      dat.r <- data.frame(type = "random",
                          #
                          TE = data$TE.random,
                          se = data$seTE.random,
                          lower = data$lower.random,
                          upper = data$upper.random,
                          statistic = data$statistic.random,
                          pval = data$pval.random,
                          df = data$df.random,
                          #
                          method = data$method.random,
                          method.ci = data$method.random.ci,
                          #
                          text = data$text.random,
                          phi = NA)
      #
      dat.t <- data.frame(type = "tau2",
                          #
                          TE = data$tau2,
                          se = data$se.tau2,
                          lower = data$lower.tau2,
                          upper = data$upper.tau2,
                          statistic = NA,
                          pval = NA,
                          df = NA,
                          #
                          method = data$method.tau,
                          method.ci = data$method.tau.ci,
                          #
                          text = "",
                          #
                          phi = replaceNULL(data$phi))
    }
    #
    if (data$prediction) {
      dat.p <- data.frame(type = "prediction",
                          #
                          TE = NA,
                          se = data$seTE.predict,
                          lower = data$lower.predict,
                          upper = data$upper.predict,
                          statistic = NA,
                          pval = NA,
                          df = data$df.predict,
                          #
                          method = data$method.predict,
                          method.ci = data$method.predict.ci,
                          #
                          text = data$text.predict,
                          phi = NA)
      #
      n.prd <- nrow(dat.p)
    }
    #
    dat <- rbind(dat.c, dat.r, dat.t, dat.p)
    #
    if (nrow(dat) == 0) {
      warning("No pooled results in meta-analysis object.",
              call. = FALSE)  
      return(res)
    }
  }
  else {
    is.meta <- FALSE
    #
    if (!missing.type)
      type <- catch("type", mc, data, sfsp)
    type <- replaceNULL(type, "")
    type <- setchar(type, c("common", "random", "prediction", "tau2"))
    #
    if (any(type %in% c("common", "random", "tau2")))
      if (missing.TE)
        stop("Argument 'TE' must be provided for common effect, ",
             "random effects model, or between-study variance.",
             call. = FALSE)
    #
    if (missing.lower & any(type %in% c("common", "random", "prediction")))
      stop("Argument 'lower' must be provided.",
           call. = FALSE)
    #
    if (missing.upper & any(type %in% c("common", "random", "prediction")))
      stop("Argument 'upper' must be provided.",
           call. = FALSE)
    #
    if (!missing.TE) {
      TE <- catch("TE", mc, data, sfsp)
      chknumeric(TE)
    }
    if (!missing.lower) {
      lower <- catch("lower", mc, data, sfsp)
      chknumeric(lower)
    }
    if (!missing.upper) {
      upper <- catch("upper", mc, data, sfsp)
      chknumeric(upper)
    }
    #
    if (!missing.statistic) {
      statistic <- catch("statistic", mc, data, sfsp)
      chknumeric(statistic)
    }
    if (!missing.pval) {
      pval <- catch("pval", mc, data, sfsp)
      chknumeric(pval)
    }
    if (!missing.df) {
      df <- catch("df", mc, data, sfsp)
      chknumeric(df)
    }
    if (!missing.se) {
      se <- catch("se", mc, data, sfsp)
      chknumeric(se, min = 0)
    }
    #
    if (!missing.method) {
      method <- catch("method", mc, data, sfsp)
      chkchar(method)
    }
    if (!missing.method.ci) {
      method.ci <- catch("method.ci", mc, data, sfsp)
      chkchar(method.ci)
    }
    #
    if (!missing.text) {
      text <- catch("text", mc, data, sfsp)
      chkchar(text)
    }
    #
    dat <- data.frame(type,
                      TE, se, lower, upper, statistic, pval, df,
                      method, method.ci, text)
    #
    # Transform added results
    #
    for (i in seq_len(nrow(dat))) {
      if (!transf && dat$type[i] != "tau2") {
        dat$TE[i] <- transf(dat$TE[i], x$sm, x$func.transf, x$args.transf)
        dat$lower[i] <- transf(dat$lower[i], x$sm, x$func.transf, x$args.transf)
        dat$upper[i] <- transf(dat$upper[i], x$sm, x$func.transf, x$args.transf)
      }
    }
  }
  
  
  ##
  ##
  ## (3) Add results
  ##
  ##
  
  init.random <- res$random
  #
  j.t <- 0
  #
  for (i in seq_len(nrow(dat))) {
    if (dat$type[i] == "common") {
      if (!res$common)
        res$w.common[!is.na(res$w.common)] <- NA
      #
      res$TE.common <-
        c(if (res$common) res$TE.common, dat$TE[i])
      res$seTE.common <-
        c(if (res$common) res$seTE.common, dat$se[i])
      res$lower.common <-
        c(if (res$common) res$lower.common, dat$lower[i])
      res$upper.common <-
        c(if (res$common) res$upper.common, dat$upper[i])
      res$statistic.common <-
        c(if (res$common) res$statistic.common, dat$statistic[i])
      res$pval.common <-
        c(if (res$common) res$pval.common, dat$pval[i])
      #
      res$method <-
        c(if (res$common) res$method, dat$method[i])
      res$method.common.ci <-
        c(if (res$common) res$method.common.ci, dat$method.ci[i])
      #
      res$text.common <-
        c(if (res$common) res$text.common, dat$text[i])
      ##
      res$common <- TRUE
      res$overall <- TRUE
    }
    ##
    if (dat$type[i] == "random") {
      if (!res$random)
        res$w.random[!is.na(res$w.random)] <- NA
      #
      res$TE.random <-
        c(if (res$random) res$TE.random, dat$TE[i])
      res$seTE.random <-
        c(if (res$random) res$seTE.random, dat$se[i])
      res$lower.random <-
        c(if (res$random) res$lower.random, dat$lower[i])
      res$upper.random <-
        c(if (res$random) res$upper.random, dat$upper[i])
      res$statistic.random <-
        c(if (res$random) res$statistic.random, dat$statistic[i])
      res$pval.random <-
        c(if (res$random) res$pval.random, dat$pval[i])
      res$df.random <-
        c(if (res$random) res$df.random, dat$df[i])
      #
      res$method.random <-
        c(if (res$random) res$method.random, dat$method[i])
      res$method.random.ci <-
        c(if (res$random) res$method.random.ci, dat$method.ci[i])
      #
      res$text.random <-
        c(if (res$random) res$text.random, dat$text[i])
      #
      res$random <- TRUE
      res$overall <- TRUE
    }
    #
    if (dat$type[i] == "prediction") {
      res$seTE.predict <-
        c(if (res$prediction) res$seTE.predict, dat$se[i])
      res$lower.predict <-
        c(if (res$prediction) res$lower.predict, dat$lower[i])
      res$upper.predict <-
        c(if (res$prediction) res$upper.predict, dat$upper[i])
      res$df.prediction <-
        c(if (res$prediction) res$df.predict, dat$df[i])
      #
      res$method.predict <-
        c(if (res$prediction) res$method.predict, dat$method[i])
      res$method.predict.ci <-
        c(if (res$prediction) res$method.predict.ci, dat$method.ci[i])
      #
      res$text.predict <-
        c(if (res$prediction) res$text.predict, dat$text[i])
      #
      res$prediction <- TRUE
      res$overall <- TRUE
    }
    #
    if (dat$type[i] == "tau2") {
      j.t <- j.t + 1
      #
      sel.t.i <- res$random & ((is.meta & (init.random & j.t > 1)) | !is.meta)
      #
      res$tau2 <- c(if (sel.t.i) res$tau2, dat$TE[i])
      res$se.tau2 <- c(if (sel.t.i) res$se.tau2, dat$se[i])
      res$lower.tau2 <- c(if (sel.t.i) res$lower.tau2, dat$lower[i])
      res$upper.tau2 <- c(if (sel.t.i) res$upper.tau2, dat$upper[i])
      #
      res$method.tau <- c(if (sel.t.i) res$method.tau, dat$method[i])
      res$method.tau.ci <- c(if (sel.t.i) res$method.tau.ci, dat$method.ci[i])
      #
      if (!is.na(dat$phi[i]))
        res$phi <- c(if (sel.t.i) res$phi, dat$phi[i])
      #
      res$random <- TRUE
    }
  }
  #
  # Change results for between-study standard deviation
  #
  res$tau <- sqrt(res$tau2)
  res$lower.tau <- sqrt(res$lower.tau2)
  res$upper.tau <- sqrt(res$upper.tau2)
  #
  if (all(is.na(res$phi)))
    res$phi <- NULL
  #
  class(res) <- c(class(res), "metaadd")
  res
}
