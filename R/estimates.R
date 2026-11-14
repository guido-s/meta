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
#' \code{\link{transf}} for the corresponding transformations and
#' back transformations.
#' 
#' @return
#' A data frame with additional class 'extract.meta' or an Excel file.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{transf}}, \code{\link{meta-sm}},
#'   \code{\link{as.data.frame.meta}}
#' 
#' @examples
#' ma <- metacor(c(0.85, 0.7, 0.95), c(20, 40, 10))
#' summary(ma)
#' estimates(ma)
#' estimates(ma, backtransf = FALSE)
#' estimates(update(ma, common = FALSE, random = FALSE))
#' estimates(update(ma, prediction = TRUE))
#' estimates(update(ma, prediction = TRUE,
#'   level.ma = 0.99, level.predict = 0.9))
#'
#' \dontrun{
#' # Create Excel file with extracted results
#' # (R package 'writexl' must be available)
#' if (requireNamespace("writexl", quietly = TRUE)) {
#'  fname1 <- tempfile(fileext = ".xlsx")
#'  estimates(ma, path = fname1)
#'  # An existing Excel file is not overwritten but a warning is printed
#'  estimates(ma, path = fname1)
#'  # Overwrite an existing Excel file
#'  estimates(ma, path = fname1, overwrite = TRUE)
#'  # Suppress message on file creation and overwrite existing file
#'  suppressMessages(estimates(ma, path = fname1, overwrite = TRUE))
#' }
#'
#' # Save the extracted results in a text file
#' fname2 <- tempfile(fileext = ".csv")
#' fname2
#' write.csv(estimates(ma), file = fname2, row.names = FALSE)
#' }
#'
#' @rdname estimates
#' @method estimates meta
#' @export

estimates.meta <- function(x,
                           sortvar,
                           #
                           study.results = TRUE,
                           common = x$common,
                           random = x$random,
                           prediction = x$prediction,
                           overall = x$overall,
                           subgroup,
                           prediction.subgroup = x$prediction.subgroup,
                           #
                           se = FALSE,
                           ci = TRUE,
                           statistic = FALSE,
                           pval = FALSE,
                           n = TRUE,
                           #
                           backtransf = x$backtransf,
                           #
                           digits = gs("digits"),
                           digits.se = gs("digits.se"),
                           digits.stat = gs("digits.stat"),
                           digits.pval = gs("digits.pval"),
                           #
                           writexl = !missing(path),
                           path = "estimates.xlsx",
                           overwrite = FALSE,
                           #
                           ...) {
  
  chkclass(x, "meta")
  chksuitable(x, method = "estimates", classes = "metabind", check.mlm = FALSE)
  x <- updateversion(x)
  #
  chklogical(study.results)
  chklogical(common)
  chklogical(random)
  chklogical(prediction)
  chklogical(overall)
  #
  if (missing(subgroup)) {
    subgroup <- !is.null(x$subgroup)
    n.subgroups <- length(x$subgroup.levels)
  }
  else {
    chklogical(subgroup)
    #
    if (is.null(x$subgroup))
      subgroup <- FALSE
    else
      n.subgroups <- length(x$subgroup.levels)
  }
  #
  if (subgroup) {
    for (model in c("common", "random", "predict")) {
      reference <- x[[paste0("lower.", model, ".w")]]
      #
      if (is.null(reference))
        next
      #
      vars <- if (model == "predict")
        c("lower", "upper", "df")
      else
        c("TE", "seTE", "lower", "upper", "statistic", "pval",
          if (model == "random") "df")
      #
      for (i in vars) {
        name <- paste0(i, ".", model, ".w")
        x[[name]] <- subgroup_result_matrix(x[[name]], reference)
      }
    }
  }
  #
  missing.prediction.subgroup <- missing(prediction.subgroup)
  #
  if (!subgroup) {
    if (!missing.prediction.subgroup)
      warning("Argument 'prediction.subgroup' only considered for ",
              "meta-analysis with subgroups.",
              call. = FALSE)
    #
    prediction.subgroup <- FALSE
  }
  else {
    prediction.subgroup <- catch(enquo(prediction.subgroup), x)
    prediction.subgroup <- replaceNULL(prediction.subgroup, FALSE)
    #
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
    #
    chklogical(prediction.subgroup[1])
    #
    if (length(prediction.subgroup) == 1 & n.subgroups > 1)
      prediction.subgroup <- rep(prediction.subgroup, n.subgroups)
  }
  #
  chklogical(se)
  chklogical(ci)
  chklogical(statistic)
  chklogical(pval)
  chklogical(n)
  #
  chklogical(backtransf)
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 0, length = 1)
  #
  chklogical(writexl)
  chkchar(path, length = 1)
  chklogical(overwrite)
  
  
  if (!study.results & !common & !random &
      !prediction & !any(prediction.subgroup)) {
    warning("Nothing to extract.")
    return(NULL)
  }
  #
  chars <- result_characteristics(x)
  #
  avail.n <- n & !is.null(x[["n"]])
  avail.n.e <- n & !is.null(x[["n.e"]])
  avail.n.c <- n & !is.null(x[["n.c"]])
  #
  res.s <- NULL
  res.c <- res.c.w <- NULL
  res.r <- res.r.w <- NULL
  res.p <- res.p.w <- NULL
  #
  study.characteristics <- NULL
  #
  if (study.results) {
    error <-
      try(sortvar <- catch(enquo(sortvar), x), silent = TRUE)
    if (inherits(error, "try-error")) {
      sortvar <- catch(enquo(sortvar), x$data)
      if (isCol(x$data, ".subset"))
        sortvar <- sortvar[x$data$.subset]
    }
    #
    sort <- !is.null(sortvar)
    if (sort && (length(sortvar) != length(x$TE)))
      stop("Number of studies in object 'x' and ",
           "argument 'sortvar' have different length.")
    #
    if (!sort)
      sortvar <- seq_along(x$TE)
    #
    o <- order(sortvar)
    #
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
    #
    NAs <- rep(NA_real_, nrow(res.s))
    #
    study.characteristics <-
      data.frame(
        k = NAs,
        n = if (avail.n)
          x[["n"]][o]
        else if (avail.n.e && avail.n.c)
          x[["n.e"]][o] + x[["n.c"]][o]
        else
          NAs,
        n.e = if (avail.n.e) x[["n.e"]][o] else NAs,
        n.c = if (avail.n.c) x[["n.c"]][o] else NAs,
        n.harmonic.mean =
          if (inherits(x, "metaprop")) x[["n"]][o] else NAs,
        t.harmonic.mean =
          if (inherits(x, "metarate")) x[["time"]][o] else NAs)
  }
  #
  if (common) {
    if (overall) {
      res.c <- data.frame(studlab = x$text.common,
                          subgroup = "",
                          estimate = x$TE.common,
                          se = x$seTE.common,
                          lower = x$lower.common,
                          upper = x$upper.common,
                          k = chars$common$overall$k,
                          df = NA,
                          statistic = x$statistic.common,
                          pval = x$pval.common)
    }
    #
    if (subgroup) {
      rep.c <- rep(length(x$text.common), n.subgroups)
      #
      TE.common.w <- x$TE.common.w
      #
      if (is.matrix(TE.common.w))
        TE.common.w <- as.vector(t(TE.common.w))
      else
        TE.common.w <- rep(TE.common.w, rep.c)
      #
      res.c.w <- data.frame(studlab = rep(x$text.common, n.subgroups),
                            subgroup = rep(x$subgroup.levels, rep.c),
                            estimate = TE.common.w,
                            se = as.vector(t(x$seTE.common.w)),
                            lower = as.vector(t(x$lower.common.w)),
                            upper = as.vector(t(x$upper.common.w)),
                            k = as.vector(t(chars$common$subgroup$k)),
                            df = NA,
                            statistic = as.vector(t(x$statistic.common.w)),
                            pval = as.vector(t(x$pval.common.w)),
                            row.names =
                              seq_along(as.vector(t(x$lower.common.w))))
    }
  }
  #
  if (random) {
    if (overall) {
      res.r <- data.frame(studlab = x$text.random,
                          subgroup = "",
                          estimate = x$TE.random,
                          se = x$seTE.random,
                          lower = x$lower.random,
                          upper = x$upper.random,
                          k = chars$random$overall$k,
                          df = x$df.random,
                          statistic = x$statistic.random,
                          pval = x$pval.random)
    }
    #
    if (subgroup) {
      rep.r <- rep(length(x$text.random), n.subgroups)
      #
      TE.random.w <- x$TE.random.w
      #
      if (is.matrix(TE.random.w))
        TE.random.w <- as.vector(t(TE.random.w))
      else
        TE.random.w <- rep(TE.random.w, rep.r)
      #
      res.r.w <- data.frame(
        studlab = rep(x$text.random, n.subgroups),
        subgroup = rep(x$subgroup.levels, rep.r),
        estimate = TE.random.w,
        se = as.vector(t(x$seTE.random.w)),
        lower = as.vector(t(x$lower.random.w)),
        upper = as.vector(t(x$upper.random.w)),
        k = as.vector(t(chars$random$subgroup$k)),
        df = as.vector(t(x$df.random.w)),
        statistic = as.vector(t(x$statistic.random.w)),
        pval = as.vector(t(x$pval.random.w)),
        row.names = seq_along(as.vector(t(x$lower.random.w))))
    }
  }
  #
  if (prediction) {
    if (overall) {
      res.p <- data.frame(studlab = x$text.predict,
                          subgroup = "",
                          estimate = NA,
                          se = NA,
                          lower = x$lower.predict,
                          upper = x$upper.predict,
                          k = chars$predict$overall$k,
                          df = x$df.predict,
                          statistic = NA,
                          pval = NA)
    }
  }
  #
  if (any(prediction.subgroup)) {
    rep.p <- rep(length(x$text.predict), n.subgroups)
    sel.p <- rep(prediction.subgroup, rep.p)
    #
    res.p.w <-
      data.frame(studlab = rep(x$text.predict, n.subgroups),
                 subgroup = rep(x$subgroup.levels, rep.p),
                 estimate = NA,
                 se = NA,
                 lower = as.vector(t(x$lower.predict.w)),
                 upper = as.vector(t(x$upper.predict.w)),
                 k = as.vector(t(chars$predict$subgroup$k)),
                 df = as.vector(t(x$df.predict.w)),
                 statistic = NA,
                 pval = NA,
                 row.names =
                   seq_along(as.vector(t(x$lower.predict.w))))[sel.p, ]
  }
  #
  res <- rbind(res.s, res.c, res.c.w, res.r, res.r.w, res.p, res.p.w)
  info <- rbind(
    study.characteristics,
    if (common && overall)
      chars$common$overall,
    if (common && subgroup)
      characteristic_rows(chars$common$subgroup),
    if (random && overall)
      chars$random$overall,
    if (random && subgroup)
      characteristic_rows(chars$random$subgroup),
    if (prediction && overall)
      chars$predict$overall,
    if (any(prediction.subgroup))
      characteristic_rows(chars$predict$subgroup)[sel.p, , drop = FALSE])
  #
  if (nrow(info) != nrow(res))
    stop("Result characteristics do not match the extracted rows.",
         call. = FALSE)
  #
  res$k <- info$k
  #
  if (n) {
    for (var in c("n", "n.e", "n.c")) {
      if (!all(is.na(info[[var]])))
        res[[var]] <- info[[var]]
    }
  }
  #
  if (backtransf) {
    if (inherits(x, "metaprop"))
      n_backtr <- info$n.harmonic.mean
    else
      n_backtr <- NULL
    #
    if (inherits(x, "metarate"))
      t_backtr <- info$t.harmonic.mean
    else
      t_backtr <- NULL
    #
    res$estimate <- backtransf(res$estimate, x$sm, n_backtr, t_backtr,
                               x$func.backtransf, x$args.backtransf)
    #
    if (ci) {
      res$lower <- backtransf(res$lower, x$sm, n_backtr, t_backtr,
                              x$func.backtransf, x$args.backtransf)
      res$upper <- backtransf(res$upper, x$sm, n_backtr, t_backtr,
                              x$func.backtransf, x$args.backtransf)
    }
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
  #
  res$estimate <- round(res$estimate, digits = digits)
  #
  if (se)
    res$se <- round(res$se, digits = digits.se)
  else
    res$se <- NULL
  #
  if (ci) {
    res$lower <- round(res$lower, digits = digits)
    res$upper <- round(res$upper, digits = digits)
    #
    if (backtransf) {
      wo <- !is.na(res$lower) & !is.na(res$upper) & res$lower > res$upper
      #
      tlower <- res$lower[wo]
      res$lower[wo] <- res$upper[wo]
      res$upper[wo] <- tlower
    }
  }
  else {
    res$lower <- NULL
    res$upper <- NULL
  }
  #
  if (statistic)
    res$statistic <- round(res$statistic, digits = digits)
  else
    res$statistic <- NULL
  #
  if (pval)
    res$pval <- round(res$pval, digits = digits.pval)
  else
    res$pval <- NULL
  #
  if (!any(is.finite(res$df)))
    res$df <- NULL
  else
    res$df <- round(ifelse(is.na(res$df), Inf, res$df), 2)
  #
  if (!subgroup)
    res$subgroup <- NULL
  
  if (writexl) {
    if (!is_installed_package("writexl", stop = FALSE))
      stop(paste0("Package 'writexl' missing.",
                  "\n  ",
                  "Please use the following R command for installation:",
                  "\n  install.packages(\"writexl\")"),
           call. = FALSE)
    #
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
    #
    return(invisible(NULL))
  }
  
  
  attr(x, ".print.study.results.") <- study.results
  #
  attr(res, "x") <- x
  #
  attr(res, "study.results") <- study.results
  attr(res, "common") <- common
  attr(res, "random") <- random
  attr(res, "prediction") <- prediction
  attr(res, "overall") <- overall
  #
  attr(res, "se") <- se
  attr(res, "ci") <- ci
  attr(res, "statistic") <- statistic
  attr(res, "pval") <- pval
  attr(res, "n") <- n
  #
  attr(res, "backtransf") <- backtransf
  #
  attr(res, "digits") <- digits
  attr(res, "digits.se") <- digits.se
  attr(res, "digits.stat") <- digits.stat
  attr(res, "digits.pval") <- digits.pval
  #
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
  #
  digits <- attr(x, "digits")
  digits.se <- attr(x, "digits.se")
  digits.stat <- attr(x, "digits.stat")
  digits.pval <- attr(x, "digits.pval")
  #
  chknumeric(digits.tau, min = 0, length = 1)
  
  common <- attr(x, "common")
  random <- attr(x, "random")
  prediction <- attr(x, "prediction")
  overall <- attr(x, "overall")
  #
  backtransf <- attr(x, "backtransf")
  sm.lab <- smlab(meta$sm, backtransf, meta$pscale, meta$irscale)
  sm.lab.se <- smlab(meta$sm, FALSE)
  #
  ci.lab <- percent_label(meta$level)
  #
  se <- attr(x, "se")
  ci <- attr(x, "ci")
  statistic <- attr(x, "statistic")
  pval <- attr(x, "pval")
  
  res <- x[, names(x) != "studlab"]
  #
  res$estimate <- round(res$estimate, digits = digits)
  #
  if (ci) {
    res$lower <- round(res$lower, digits = digits)
    res$upper <- round(res$upper, digits = digits)
  }
  #
  res$estimate <- formatN(res$estimate, digits, "", big.mark = big.mark)
  #
  if (se)
    res$se <- formatN(res$se, digits, "", big.mark = big.mark)
  #
  if (ci)
    res$lower <-
      ifelse(is.na(res$lower) & is.na(res$upper), "",
             formatCI(
               formatN(res$lower, digits, "", big.mark = big.mark),
               formatN(res$upper, digits, "", big.mark = big.mark)))
  #
  if (statistic)
    res$statistic <-
      formatN(res$statistic, digits.stat, "", big.mark = big.mark)
  if (pval)
    res$pval <- formatPT(res$pval, digits = digits.pval,
                         lab.NA = "", big.mark = big.mark)
  #
  if (all(is.na(res$k)))
    res$k <- NULL
  else
    res$k <- formatN(res$k, 0, "", big.mark = big.mark)
  #
  if (!is.null(res$df))
    res$df <- ifelse(is.finite(res$df), res$df, "")
  #
  res$upper <- NULL
  #
  res <- as.matrix(res)
  #
  rownames(res) <- x$studlab
  #
  cn <- colnames(res)
  cn[cn == "estimate"] <- sm.lab
  cn[cn == "se"] <- paste0("SE(", sm.lab.se, ")")
  cn[cn == "lower"] <- ci.lab
  colnames(res) <- cn
  
  crtitle(meta)
  #
  prmatrix(res, quote = FALSE, right = TRUE, ...)
  #
  if (details) {
    catmeth(meta,
            common, random, prediction, overall, random | prediction,
            #
            func.transf = NULL,
            backtransf = backtransf, func.backtransf = NULL,
            #
            big.mark = gs("big.mark"), digits = digits, digits.tau = digits.tau,
            text.tau = gs("text.tau"), text.tau2 = gs("text.tau2"))
    #
    if ((common | random) && meta$level != meta$level.ma)
      cat(paste0("- ", percent_label(meta$level.ma),
                 " calculated for meta-analysis results\n"))
    if (prediction && meta$level != meta$level.predict)
      if (meta$level != meta$level.ma)
        cat(paste0("- ", percent_label(meta$level.predict, "PI"),
                   " calculated\n"))
  }
  
  invisible(NULL)
}
