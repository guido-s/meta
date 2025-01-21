#' Extract results from meta-analysis object
#' 
#' @description
#' Extract study and meta-analysis results from meta-analysis object
#' which can be stored in an Excel file.
#' 
#' @param x A meta-analysis object of class \code{meta}.
#' @param sortvar An optional vector used to sort the individual
#'   studies (must be of same length as \code{x$TE}).
#' @param study.results A logical indicating whether study results
#'   should be extracted.
#' @param common A logical indicating whether results of common effect
#'   meta-analysis should be extracted.
#' @param random A logical indicating whether results of random
#'   effects meta-analysis should be extracted.
#' @param prediction A logical indicating whether prediction interval
#'   should be extracted.
#' @param overall A logical indicating whether overall summaries
#'   should be extracted. This argument is useful in a meta-analysis
#'   with subgroups if overall results should not be extracted.
#' @param subgroup A logical indicating whether subgroup results
#'   should be extracted.
#' @param prediction.subgroup A single logical or logical vector
#'   indicating whether / which prediction intervals should be
#'   extracted for subgroups.
#' @param backtransf A logical indicating whether extracted results
#'   should be back transformed.
#' @param se A logical indicating whether standard errors should be
#'   extracted.
#' @param ci A logical indicating whether confidence / prediction
#'   interval should be extracted.
#' @param statistic A logical indicating whether to extract statistic
#'   of test for overall effect.
#' @param pval A logical indicating whether to extract p-value of test
#'   for overall effect.
#' @param n A logical indicating whether sample sizes should be
#'   extracted (if available).
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.se Minimal number of significant digits for standard
#'   errors, see \code{print.default}.
#' @param digits.stat Minimal number of significant digits for z- or
#'   t-statistic for test of overall effect, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall treatment effect, see \code{print.default}.
#' @param digits.tau Minimal number of significant digits for square
#'   root of between-study variance, see \code{print.default}.
#' @param text.tau2 Text printed to identify between-study variance
#'   \eqn{\tau^2}.
#' @param text.tau Text printed to identify \eqn{\tau}, the square
#'   root of the between-study variance \eqn{\tau^2}.
#' @param big.mark A character used as thousands separator.
#' @param details A logical specifying whether details on statistical
#'   methods should be printed.
#' @param writexl A logical indicating whether an Excel file should be
#'   created (R package \bold{writexl} must be available).
#' @param path A character string specifying the filename of the Excel
#'   file.
#' @param overwrite A logical indicating whether an existing Excel
#'   file should be overwritten.
#' @param \dots Additional arguments passed on to \code{prmatrix}.
#' 
#' @details
#' Extract study and meta-analysis results from meta-analysis object.
#' By default, a data frame with the results is
#' generated. Alternatively, an Excel file is created if argument
#' \code{writexl = TRUE}.
#'
#' The following information is extracted.
#'
#' \tabular{ll}{
#' \bold{Variable} \tab \bold{Content} \cr
#' studlab \tab Study label / descriptor of meta-analysis result \cr
#' estimate \tab (Back transformed) estimate for individual studies /
#'   meta-analysis \cr
#' se \tab Standard error \cr
#' lower \tab Lower (back transformed) confidence / prediction
#'   interval limit \cr
#' upper \tab Upper (back transformed) confidence / prediction
#'   interval limit \cr
#' k \tab Number of studies \cr
#' df \tab Degrees of freedom for confidence / prediction intervals \cr
#' statistic \tab Statistic for test of effect \cr
#' pval \tab P-value for test of effect \cr
#' n \tab Total sample size \cr
#' n.e \tab Sample size in first (experimental) group \cr
#' n.c \tab Sample size in second (control) group
#' }
#'
#' Some variables are only extracted if the corresponding logical
#' argument \code{se}, \code{ci}, \code{statistic} or \code{n} is
#' \code{TRUE}. Furthermore, (group) sample sizes are only extracted
#' if available in the meta-analysis object.
#'
#' The variables \code{estimate}, \code{lower} and \code{upper}
#' contain the back transformed effects and confidence interval
#' limits, e.g., odds ratio (argument \code{sm = "OR"} in
#' \code{\link{metabin}} or \code{\link{metagen}}) or correlation
#' (\code{sm = "ZCOR"} in \code{\link{metacor}} or
#' \code{\link{metagen}}), if argument \code{backtransf} is
#' \code{TRUE}. Otherwise, these variables contain the transformed
#' values, e.g., log odds ratio (\code{sm = "OR"}) or Fisher's Z
#' transformed correlations (\code{sm = "ZCOR"}). See
#' \code{\link{meta-sm}} for available summary measures and
#' \code{\link{meta-transf}} for the corresponding transformations and
#' back transformations.
#' 
#' @return
#' A data frame with additional class 'extract.meta' or an Excel file.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{meta-transf}}, \code{\link{meta-sm}},
#'   \code{\link{as.data.frame.meta}}
#' 
#' @examples
#' 
#' m1 <- metacor(c(0.85, 0.7, 0.95), c(20, 40, 10))
#' summary(m1)
#' estimates(m1)
#' estimates(m1, backtransf = FALSE)
#' estimates(update(m1, common = FALSE, random = FALSE))
#' estimates(update(m1, prediction = TRUE))
#' estimates(update(m1, prediction = TRUE,
#'   level.ma = 0.99, level.predict = 0.9))
#'
#' \dontrun{
#' # Create Excel file with extracted results
#' # (R package 'writexl' must be available)
#' if (requireNamespace("writexl", quietly = TRUE)) {
#'  fname1 <- tempfile(fileext = ".xlsx")
#'  estimates(m1, path = fname1)
#'  # An existing Excel file is not overwritten but a warning is printed
#'  estimates(m1, path = fname1)
#'  # Overwrite an existing Excel file
#'  estimates(m1, path = fname1, overwrite = TRUE)
#'  # Suppress message on file creation and overwrite existing file
#'  suppressMessages(estimates(m1, path = fname1, overwrite = TRUE))
#' }
#'
#' # Save the extracted results in a text file
#' fname2 <- tempfile(fileext = ".csv")
#' fname2
#' write.csv(estimates(m1), file = fname2, row.names = FALSE)
#' }
#'
#' @rdname estimates
#' @method estimates meta
#' @export


estimates.meta <- function(x,
                         sortvar,
                         ##
                         study.results = TRUE,
                         common = x$common,
                         random = x$random,
                         prediction = x$prediction,
                         overall = x$overall,
                         subgroup,
                         prediction.subgroup = x$prediction.subgroup,
                         ##
                         se = FALSE,
                         ci = TRUE,
                         statistic = FALSE,
                         pval = FALSE,
                         n = TRUE,
                         ##
                         backtransf = x$backtransf,
                         ##
                         digits = gs("digits"),
                         digits.se = gs("digits.se"),
                         digits.stat = gs("digits.stat"),
                         digits.pval = gs("digits.pval"),
                         ##
                         writexl = !missing(path),
                         path = "estimates.xlsx",
                         overwrite = FALSE,
                         ##
                         ...) {
  
  chkclass(x, "meta")
  x <- updateversion(x)
  chksuitable(x, method = "estimates", classes = "metabind", check.mlm = FALSE)
  ##
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  chklogical(study.results)
  chklogical(common)
  chklogical(random)
  chklogical(prediction)
  chklogical(overall)
  ##
  if (missing(subgroup)) {
    subgroup <- !is.null(x$subgroup)
    n.subgroups <- length(x$subgroup.levels)
  }
  else {
    chklogical(subgroup)
    if (is.null(x$subgroup))
      subgroup <- FALSE
    else
      n.subgroups <- length(x$subgroup.levels)
  }
  ##
  missing.prediction.subgroup <- missing(prediction.subgroup)
  ##
  if (!subgroup) {
    if (!missing.prediction.subgroup)
      warning("Argument 'prediction.subgroup' only considered for ",
              "meta-analysis with subgroups.",
              call. = FALSE)
    ##
    prediction.subgroup <- FALSE
  }
  else {
    prediction.subgroup <- catch("prediction.subgroup", mc, x, sfsp)
    prediction.subgroup <- replaceNULL(prediction.subgroup, FALSE)
    ##
    if (length(prediction.subgroup) == 1) {
      if (is.matrix(x$lower.predict.w))
        prediction.subgroup.logical <-
          prediction.subgroup &
          apply(x$lower.predict.w, 1, notallNA) &
          apply(x$upper.predict.w, 1, notallNA)
      else {
        prediction.subgroup.logical <-
          prediction.subgroup &
          notallNA(x$lower.predict.w) &
          notallNA(x$upper.predict.w)
        prediction.subgroup.logical <-
          rep(prediction.subgroup.logical, n.subgroups)
      }
    }
    else {
      chklength(prediction.subgroup, n.subgroups,
                text = paste("Length of argument 'prediction.subgroup' must be",
                             "equal to 1 or number of subgroups."))
    }
    ##
    chklogical(prediction.subgroup[1])
    ##
    if (length(prediction.subgroup) == 1 & n.subgroups > 1)
      prediction.subgroup <- rep(prediction.subgroup, n.subgroups)
  }
  ##
  chklogical(se)
  chklogical(ci)
  chklogical(statistic)
  chklogical(pval)
  chklogical(n)
  ##
  chklogical(backtransf)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 0, length = 1)
  ##
  chklogical(writexl)
  chkchar(path, length = 1)
  chklogical(overwrite)
  
  
  if (!study.results & !common & !random &
      !prediction & !any(prediction.subgroup)) {
    warning("Nothing to extract.")
    return(NULL)
  }
  ##
  avail.n <- n & !is.null(x$n)
  avail.n.e <- n & !is.null(x$n.e)
  avail.n.c <- n & !is.null(x$n.c)
  ##
  res.s <- NULL
  res.c <- res.c.w <- NULL
  res.r <- res.r.w <- NULL
  res.p <- res.p.w <- NULL
  ##
  n.s <- n.e.s <- n.c.s <- NULL
  n.c <- n.e.c <- n.c.c <- NULL
  n.c.w <- n.e.c.w <- n.c.c.w <- NULL
  n.r <- n.e.r <- n.c.r <- NULL
  n.r.w <- n.e.r.w <- n.c.r.w <- NULL
  n.p <- n.e.p <- n.c.p <- NULL
  n.p.w <- n.e.p.w <- n.c.p.w <- NULL
  ##
  if (study.results) {
    error <-
      try(sortvar <-
            catch("sortvar", mc, x, sfsp),
          silent = TRUE)
    if (inherits(error, "try-error")) {
      sortvar <- catch("sortvar", mc, x$data,  NULL)
      if (isCol(x$data, ".subset"))
        sortvar <- sortvar[x$data$.subset]
    }
    sort <- !is.null(sortvar)
    if (sort && (length(sortvar) != length(x$TE)))
      stop("Number of studies in object 'x' and ",
           "argument 'sortvar' have different length.")
    if (!sort)
      sortvar <- seq_along(x$TE)
    ##
    o <- order(sortvar)
    ##
    res.s <- data.frame(studlab = x$studlab,
                        subgroup = if (subgroup) x$subgroup else "",
                        estimate = x$TE,
                        se = x$seTE,
                        lower = x$lower,
                        upper = x$upper,
                        k = NA,
                        df = if (is.null(x$df)) NA else x$df,
                        statistic = x$statistic,
                        pval = x$pval)[o, ]
    ##
    if (avail.n)
      n.s <- x$n[o]
    else if (avail.n.e & avail.n.c)
      n.s <- x$n.e[o] + x$n.c[o]
    ##
    if (avail.n.e)
      n.e.s <- x$n.e[o]
    ##
    if (avail.n.c)
      n.c.s <- x$n.c[o]
  }
  ##
  if (common) {
    if (overall) {
      res.c <- data.frame(studlab = x$text.common,
                          subgroup = "",
                          estimate = x$TE.common,
                          se = x$seTE.common,
                          lower = x$lower.common,
                          upper = x$upper.common,
                          k = x$k,
                          df = NA,
                          statistic = x$statistic.common,
                          pval = x$pval.common)
      ##
      if (avail.n)
        n.c <- rep(sum(x$n), length(x$lower.common))
      else if (avail.n.e & avail.n.c)
        n.c <- rep(sum(x$n.e + x$n.c), length(x$lower.common))
      ##
      if (avail.n.e)
        n.e.c <- rep(sum(x$n.e), length(x$lower.common))
      ##
      if (avail.n.c)
        n.c.c <- rep(sum(x$n.c), length(x$lower.common))
    }
    ##
    if (subgroup) {
      rep.c <- rep(length(x$text.common), n.subgroups)
      ##
      res.c.w <- data.frame(studlab = rep(x$text.common, n.subgroups),
                            subgroup = rep(x$subgroup.levels, rep.c),
                            estimate = x$TE.common.w,
                            se = as.vector(t(x$seTE.common.w)),
                            lower = as.vector(t(x$lower.common.w)),
                            upper = as.vector(t(x$upper.common.w)),
                            k = rep(x$k.w, rep.c),
                            df = NA,
                            statistic = as.vector(t(x$statistic.common.w)),
                            pval = as.vector(t(x$pval.common.w)),
                            row.names =
                              seq_along(as.vector(t(x$lower.common.w))))
      ##
      if (avail.n)
        n.c.w <- rep(x$n.w, rep.c)
      else if (avail.n.e & avail.n.c)
        n.c.w <- rep(x$n.e.w + x$n.c.w, rep.c)
      ##
      if (avail.n.e)
        n.e.c.w <- rep(x$n.e.w, rep.c)
      ##
      if (avail.n.c)
        n.c.c.w <- rep(x$n.c.w, rep.c)
    }
  }
  ##
  if (random) {
    if (overall) {
      res.r <- data.frame(studlab = x$text.random,
                          subgroup = "",
                          estimate = x$TE.random,
                          se = x$seTE.random,
                          lower = x$lower.random,
                          upper = x$upper.random,
                          k = x$k,
                          df = x$df.random,
                          statistic = x$statistic.random,
                          pval = x$pval.random)
      ##
      if (avail.n)
        n.r <- rep(sum(x$n), length(x$lower.random))
      else if (avail.n.e & avail.n.c)
        n.r <- rep(sum(x$n.e + x$n.c), length(x$lower.random))
      ##
      if (avail.n.e)
        n.e.r <- rep(sum(x$n.e), length(x$lower.random))
      ##
      if (avail.n.c)
        n.c.r <- rep(sum(x$n.c), length(x$lower.random))
    }
    ##
    if (subgroup) {
      rep.r <- rep(length(x$text.random), n.subgroups)
      ##
      res.r.w <- data.frame(studlab = rep(x$text.random, n.subgroups),
                            subgroup = rep(x$subgroup.levels, rep.r),
                            estimate = x$TE.random.w,
                            se = as.vector(t(x$seTE.random.w)),
                            lower = as.vector(t(x$lower.random.w)),
                            upper = as.vector(t(x$upper.random.w)),
                            k = rep(x$k.w, rep.r),
                            df = as.vector(t(x$df.random.w)),
                            statistic = as.vector(t(x$statistic.random.w)),
                            pval = as.vector(t(x$pval.random.w)),
                            row.names =
                              seq_along(as.vector(t(x$lower.random.w))))
      ##
      if (avail.n)
        n.r.w <- rep(x$n.w, rep.r)
      else if (avail.n.e & avail.n.c)
        n.r.w <- rep(x$n.e.w + x$n.c.w, rep.r)
      ##
      if (avail.n.e)
        n.e.r.w <- rep(x$n.e.w, rep.r)
      ##
      if (avail.n.c)
        n.c.r.w <- rep(x$n.c.w, rep.r)
    }
  }
  ##
  if (prediction) {
    if (overall) {
      res.p <- data.frame(studlab = x$text.predict,
                          subgroup = "",
                          estimate = NA,
                          se = NA,
                          lower = x$lower.predict,
                          upper = x$upper.predict,
                          k = x$k,
                          df = x$df.predict,
                          statistic = NA,
                          pval = NA)
      ##
      if (avail.n)
        n.p <- rep(sum(x$n), length(x$lower.predict))
      else if (avail.n.e & avail.n.c)
        n.p <- rep(sum(x$n.e + x$n.c), length(x$lower.predict))
      ##
      if (avail.n.e)
        n.e.p <- rep(sum(x$n.e), length(x$lower.predict))
      ##
      if (avail.n.c)
        n.c.p <- rep(sum(x$n.c), length(x$lower.predict))
    }
  }
  ##
  if (any(prediction.subgroup)) {
    rep.p <- rep(length(x$text.predict), n.subgroups)
    sel.p <- rep(prediction.subgroup, rep.p)
    ##
    res.p.w <-
      data.frame(studlab = rep(x$text.predict, n.subgroups),
                 subgroup = rep(x$subgroup.levels, rep.p),
                 estimate = NA,
                 se = NA,
                 lower = as.vector(t(x$lower.predict.w)),
                 upper = as.vector(t(x$upper.predict.w)),
                 k = rep(x$k.w, rep.p),
                 df = as.vector(t(x$df.predict.w)),
                 statistic = NA,
                 pval = NA,
                 row.names =
                   seq_along(as.vector(t(x$lower.predict.w))))[sel.p, ]
    ##
    if (avail.n)
      n.p.w <- rep(x$n.w, rep.p)
    else if (avail.n.e & avail.n.c)
      n.p.w <- rep(x$n.e.w, rep.p) + rep(x$n.c.w, rep.p)
    ##
    if (avail.n.e)
      n.e.p.w <- rep(x$n.e.w, rep.p)
    ##
    if (avail.n.c)
      n.c.p.w <- rep(x$n.c.w, rep.p)
    ##
    n.p.w <- n.p.w[sel.p]
    n.e.p.w <- n.e.p.w[sel.p]
    n.c.p.w <- n.c.p.w[sel.p]
  }
  ##
  res <- rbind(res.s, res.c, res.c.w, res.r, res.r.w, res.p, res.p.w)
  ##
  if (n) {
    res$n <- c(n.s, n.c, n.c.w, n.r, n.r.w, n.p, n.p.w)
    res$n.e <- c(n.e.s, n.e.c, n.e.c.w, n.e.r, n.e.r.w, n.e.p, n.e.p.w)
    res$n.c <- c(n.c.s, n.c.c, n.c.c.w, n.c.r, n.c.r.w, n.c.p, n.c.p.w)
  }
  ##
  if (backtransf) {
    if (inherits(x, "metaprop")) {
      harmonic.mean <- 1 / mean(1 / x$n)
      ##
      nback <- c(if (study.results)
                   x$n[o],
                 if (common & overall)
                   rep(harmonic.mean, length(x$lower.common)),
                 if (common & subgroup)
                   rep(x$n.harmonic.mean.w, rep.c),
                 if (random & overall)
                   rep(harmonic.mean, length(x$lower.random)),
                 if (random & subgroup)
                   rep(x$n.harmonic.mean.w, rep.r),
                 if (prediction & overall)
                   rep(harmonic.mean, length(x$lower.predict)),
                 if (any(prediction.subgroup))
                   rep(x$n.harmonic.mean.w, rep.p)[sel.p]
                 )
    }
    else
      nback <- NULL
    ##
    if (inherits(x, "metarate")) {
      harmonic.mean <- 1 / mean(1 / x$time)
      ##
      timeback <- c(if (study.results)
                      x$time[o],
                    if (common & overall)
                      rep(harmonic.mean, length(x$lower.common)),
                    if (common & subgroup)
                      rep(x$t.harmonic.mean.w, rep.c),
                    if (random & overall)
                      rep(harmonic.mean, length(x$lower.random)),
                    if (random & subgroup)
                      rep(x$t.harmonic.mean.w, rep.r),
                    if (prediction & overall)
                      rep(harmonic.mean, length(x$lower.predict)),
                    if (any(prediction.subgroup))
                      rep(x$t.harmonic.mean.w, rep.p)[sel.p]
                    )
    }
    else
      timeback <- NULL
    ##
    res$estimate <- backtransf(res$estimate, x$sm, nback, timeback,
                               x$func.backtransf, x$args.backtransf)
    if (ci)
      res$lower <- backtransf(res$lower, x$sm, nback, timeback,
                              x$func.backtransf, x$args.backtransf)
    if (ci)
      res$upper <- backtransf(res$upper, x$sm, nback, timeback,
                              x$func.backtransf, x$args.backtransf)
    #
    if (inherits(x, "metaprop") && x$pscale != 1) {
      res$estimate <- x$pscale * res$estimate
      if (ci) {
        res$lower <- x$pscale * res$lower
        res$upper <- x$pscale * res$upper
      }
    }
    #
    if (inherits(x, "metarate") && x$irscale != 1) {
      res$estimate <- x$irscale * res$estimate
      if (ci) {
        res$lower <- x$irscale * res$lower
        res$upper <- x$irscale * res$upper
      }
    }
  }
  ##
  res$estimate <- round(res$estimate, digits = digits)
  ##
  if (se)
    res$se <- round(res$se, digits = digits.se)
  else
    res$se <- NULL
  ##
  if (ci)
    res$lower <- round(res$lower, digits = digits)
  else
    res$lower <- NULL
  ##
  if (ci)
    res$upper <- round(res$upper, digits = digits)
  else
    res$upper <- NULL
  ##
  if (statistic)
    res$statistic <- round(res$statistic, digits = digits)
  else
    res$statistic <- NULL
  ##
  if (pval)
    res$pval <- round(res$pval, digits = digits.pval)
  else
    res$pval <- NULL
  ##
  if (!any(is.finite(res$df)))
    res$df <- NULL
  else
    res$df <- round(ifelse(is.na(res$df), Inf, res$df), 2)
  ##
  if (!subgroup)
    res$subgroup <- NULL
  
  if (writexl) {
    if (!is_installed_package("writexl", stop = FALSE))
      stop(paste0("Package 'writexl' missing.",
                  "\n  ",
                  "Please use the following R command for installation:",
                  "\n  install.packages(\"writexl\")"),
           call. = FALSE)
    ##
    if (file.exists(path) & !overwrite)
      warning("File '", path, "' exists. ",
              "Use argument 'overwrite = TRUE' to overwrite file.",
              call. = FALSE)
    else {
      names(res)[names(res) == "estimate"] <-
        smlab(x$sm, backtransf, x$pscale, x$irscale)
      #
      writexl::write_xlsx(res, path = path, col_names = TRUE, ...)
      message(paste0("Extracted information saved in file '", path, "'."))
    }
    ##
    return(invisible(NULL))
  }
  
  
  attr(x, ".print.study.results.") <- study.results
  ##
  attr(res, "x") <- x
  ##
  attr(res, "study.results") <- study.results
  attr(res, "common") <- common
  attr(res, "random") <- random
  attr(res, "prediction") <- prediction
  attr(res, "overall") <- overall
  ##
  attr(res, "se") <- se
  attr(res, "ci") <- ci
  attr(res, "statistic") <- statistic
  attr(res, "pval") <- pval
  attr(res, "n") <- n
  ##
  attr(res, "backtransf") <- backtransf
  ##
  attr(res, "digits") <- digits
  attr(res, "digits.se") <- digits.se
  attr(res, "digits.stat") <- digits.stat
  attr(res, "digits.pval") <- digits.pval
  ##
  class(res) <- c("estimates.meta", "estimates", "data.frame")  
  
  res
}


#' @rdname estimates
#' @export estimates


estimates <- function(x, ...) 
  UseMethod("estimates")


#' @rdname estimates
#' @method print estimates.meta
#' @export

print.estimates.meta <- function(x,
                                 digits.tau = gs("digits.tau"),
                                 text.tau2 = gs("text.tau2"),
                                 text.tau = gs("text.tau"),
                                 big.mark = gs("big.mark"),
                                 details = TRUE, ...) {
  
  chkclass(x, "estimates.meta")
  
  meta <- attr(x, "x")
  ##
  digits <- attr(x, "digits")
  digits.se <- attr(x, "digits.se")
  digits.stat <- attr(x, "digits.stat")
  digits.pval <- attr(x, "digits.pval")
  ##
  chknumeric(digits.tau, min = 0, length = 1)
  
  common <- attr(x, "common")
  random <- attr(x, "random")
  prediction <- attr(x, "prediction")
  overall <- attr(x, "overall")
  ##
  backtransf <- attr(x, "backtransf")
  sm.lab <- smlab(meta$sm, backtransf, meta$pscale, meta$irscale)
  #
  ci.lab <- paste0(round(100 * meta$level, 1), "%-CI")
  #
  se <- attr(x, "se")
  ci <- attr(x, "ci")
  statistic <- attr(x, "statistic")
  pval <- attr(x, "pval")
  
  res <- x[, names(x) != "studlab"]
  ##
  res$estimate <- round(res$estimate, digits = digits)
  if (ci)
    res$lower <- round(res$lower, digits = digits)
  if (ci)
    res$upper <- round(res$upper, digits = digits)
  ##
  res$estimate <- formatN(res$estimate, digits, "", big.mark = big.mark)
  ##
  if (se)
    res$se <- formatN(res$se, digits, "", big.mark = big.mark)
  ##
  if (ci)
    res$lower <-
      ifelse(is.na(res$lower) & is.na(res$upper), "",
             formatCI(
               formatN(res$lower, digits, "", big.mark = big.mark),
               formatN(res$upper, digits, "", big.mark = big.mark)))
  ##
  if (statistic)
    res$statistic <-
      formatN(res$statistic, digits.stat, "", big.mark = big.mark)
  if (pval)
    res$pval <- formatPT(res$pval, digits = digits.pval,
                         lab.NA = "", big.mark = big.mark)
  ##
  if (all(is.na(res$k)))
    res$k <- NULL
  else
    res$k <- formatN(res$k, 0, "", big.mark = big.mark)
  ##
  if (!is.null(res$df))
    res$df <- ifelse(is.finite(res$df), res$df, "")
  ##
  res$upper <- NULL
  ##
  res <- as.matrix(res)
  ##
  rownames(res) <- x$studlab
  ##
  cn <- colnames(res)
  cn[cn == "estimate"] <- sm.lab
  cn[cn == "lower"] <- ci.lab
  colnames(res) <- cn
  
  crtitle(meta)
  ##
  prmatrix(res, quote = FALSE, right = TRUE, ...)
  ##
  if (details) {
    catmeth(meta,
            common, random, prediction, overall, random | prediction,
            #
            func.transf = NULL,
            backtransf = backtransf, func.backtransf = NULL,
            #
            big.mark = gs("big.mark"), digits = digits, digits.tau = digits.tau,
            text.tau = gs("text.tau"), text.tau2 = gs("text.tau2"))
    ##
    if ((common | random) && meta$level != meta$level.ma)
      cat(paste0("- ", round(100 * meta$level.ma, 1),
                 "%-CI calculated for meta-analysis results\n"))
    if (prediction && meta$level != meta$level.predict)
      if (meta$level != meta$level.ma)
        cat(paste0("- ", round(100 * meta$level.predict, 1),
                   "%-PI calculated\n"))
  }
  
  invisible(NULL)
}
