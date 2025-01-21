#' Calculate best linear unbiased predictor for \code{meta} object
#' 
#' @description
#' Calculate best linear unbiased predictors (BLUPs) for meta-analysis object
#' created with R package \bold{meta}.
#' 
#' @aliases blup blup.meta
#' 
#' @param x An object of class \code{meta}, \code{blup.meta}, or
#'   \code{estimates.blup.meta}.
#' @param level The level used to calculate prediction intervals for BLUPs.
#' @param backtransf A logical indicating whether BLUPs should be
#'   back transformed. If \code{backtransf = TRUE}, results for
#'   \code{sm = "OR"} will be odds ratios rather than log odds ratios,
#'    for example.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.se Minimal number of significant digits for standard errors.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance \eqn{\tau^2}, see \code{print.default}.
#' @param digits.tau Minimal number of significant digits for
#'   \eqn{\tau}, the square root of the between-study variance
#'   \eqn{\tau^2}.
#' @param big.mark A character used as thousands separator.
#' @param se A logical indicating whether standard errors should be
#'   printed / extracted.
#' @param print.tau2 A logical specifying whether between-study
#'   variance \eqn{\tau^2} should be printed.
#' @param print.tau A logical specifying whether \eqn{\tau}, the
#'   square root of the between-study variance \eqn{\tau^2}, should be
#'   printed.
#' @param details A logical specifying whether details on
#'   statistical methods should be printed.
#' @param writexl A logical indicating whether an Excel file should be
#'   created (R package \bold{writexl} must be available).
#' @param path A character string specifying the filename of the Excel
#'   file.
#' @param overwrite A logical indicating whether an existing Excel
#'   file should be overwritten.
#' @param \dots Additional arguments (passed on to \code{\link{prmatrix}}).
#' 
#' @return
#' Data frame with variables
#' \item{studlab}{Study label}
#' \item{blup}{estimated best linear unbiased predictor}
#' \item{se.blup}{standard error (only if argument \code{backtransf = FALSE})}
#' \item{lower}{lower prediction limits}
#' \item{upper}{upper prediction limits}
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}, \code{\link[metadat]{dat.bcg}}
#' 
#' @examples
#' m1 <- metabin(tpos, tpos + tneg, cpos, cpos + cneg,
#'   data = dat.bcg, studlab = paste(author, year), method = "Inverse")
#' summary(m1)
#' blup(m1)
#' 
#' \dontrun{
#' # Save estimates in Excel file (R package 'writexl' must be available)
#' if (requireNamespace("writexl", quietly = TRUE))
#'  estimates(blup(m1), path = "blup_m1.xlsx")
#' }
#' 
#' @method blup meta
#' @export

blup.meta <- function(x, level = x$level, backtransf = x$backtransf,
                      ...) {
  
  #
  # (1) Check for meta object and upgrade older meta object
  #
  chkclass(x, "meta")
  x <- updateversion(x)
  meta <- x
  
  #
  # (2) Check additional arguments
  #
  chklevel(level)
  #
  if (is_untransformed(x$sm))
    backtransf <- TRUE
  chklogical(backtransf)
  #
  by <- !is.null(x$subgroup)
  
  #
  # (3) Calculate BLUP in meta-analysis without subgroups
  #
  if (!by) {
    TE.random <- x$TE.random
    seTE.random <- x$seTE.random
    df.hakn.ci <- x$df.hakn.ci
    df.kero <- round(x$df.kero, 2)
    #
    if (length(x$method.random.ci) == 1) {
      names(TE.random) <- x$method.random.ci
      names(seTE.random) <- x$method.random.ci
      names(df.hakn.ci) <- x$method.random.ci
      names(df.kero) <- x$method.random.ci
    }
    #
    frac <- x$tau2 / (x$tau2 + x$seTE^2)
    res <- data.frame()
    #
    for (i in x$method.random.ci) {
      dat.i <- as.data.frame(x)
      dat.i$blup <- frac * x$TE + (1 - frac) * TE.random[i]
      dat.i$se.blup <- sqrt(frac * x$seTE^2 + (1 - frac)^2 * seTE.random[i]^2)
      dat.i$tau2 <- x$tau2
      dat.i$method <- i
      #
      dat.i$df.hakn <- x$df.hakn.ci[i]
      dat.i$df.kero <- round(x$df.kero, 2)
      #
      if (i %in% c("classic", "classic-KR")) {
        ci.blup.i <- ci(dat.i$blup, dat.i$se.blup, level = level)
      }
      else if (i == "HK") {
        ci.blup.i <-
          ci(dat.i$blup, dat.i$se.blup, level = level, df = dat.i$df.hakn)
      }
      else if (i == "KR") {
        ci.blup.i <-
          ci(dat.i$blup, dat.i$se.blup, level = level, df = dat.i$df.kero)
      }
      else {
        ci.blup.i <- list(lower = NA, upper = NA)
      }
      #
      res.i <- data.frame(studlab = x$studlab,
                          blup = dat.i$blup, se.blup = dat.i$se.blup,
                          lower = ci.blup.i$lower,
                          upper = ci.blup.i$upper)
      #
      if (length(x$method.random.ci) > 1)
        res.i$method <- i
      if (!any(x$method.random.ci == "HK"))
        res.i$df.hakn <- NULL
      if (!any(x$method.random.ci == "KR"))
        res.i$df.kero <- NULL
      #
      if (nrow(res) == 0)
        res <- res.i
      else
        res <- rbind(res, res.i)
    }
  }
  
  #
  # (4) Calculate BLUP in meta-analysis with subgroups
  #
  if (by) {
    TE.random <- x$TE.random.w
    seTE.random <- x$seTE.random.w
    #
    if (length(x$method.random.ci) == 1) {
      seTE.random <- matrix(seTE.random, ncol = 1)
      rownames(seTE.random) <- x$subgroup.levels
      colnames(seTE.random) <- x$method.random.ci
    }
    #
    res <- data.frame()
    #
    for (i in x$method.random.ci) {
      dat.i <- as.data.frame(x)
      dat.i$blup <- NA
      dat.i$se.blup <- NA
      dat.i$tau2 <- NA
      #
      for (g in x$subgroup.levels) {
        sel.g <- dat.i$subgroup == g
        tau2.g <- replaceNA(x$tau2.w[g], 0)
        #
        frac.g <- tau2.g / (tau2.g + dat.i$seTE[sel.g]^2)
        #
        dat.i$blup[sel.g] <-
          frac.g * dat.i$TE[sel.g] + (1 - frac.g) * TE.random[g]
        #
        dat.i$se.blup[sel.g] <-
          sqrt(frac.g * x$seTE[sel.g]^2 + (1 - frac.g)^2 * seTE.random[g, i]^2)
        #
        dat.i$df.hakn[sel.g] <- x$k.w[g] - 1
        dat.i$df.kero[sel.g] <- round(x$df.kero.w[g], 2)
        #
        dat.i$tau2[sel.g] <- tau2.g
      }
      #
      blup.i <- dat.i$blup
      se.blup.i <- dat.i$se.blup
      df.hakn <- dat.i$df.hakn
      df.kero <- dat.i$df.kero
      #
      if (i %in% c("classic", "classic-KR"))
        ci.blup.i <- ci(blup.i, se.blup.i, level = level)
      else if (i == "HK")
        ci.blup.i <- ci(blup.i, se.blup.i, level = level, df = df.hakn)
      else if (i == "KR")
        ci.blup.i <- ci(blup.i, se.blup.i, level = level, df = df.kero)
      else
        ci.blup.i <- list(lower = NA, upper = NA)
      #
      res.i <- data.frame(studlab = x$studlab,
                          subgroup = dat.i$subgroup,
                          blup = blup.i, se.blup = se.blup.i,
                          lower = ci.blup.i$lower,
                          upper = ci.blup.i$upper,
                          tau2 = dat.i$tau2,
                          df.hakn = dat.i$df.hakn,
                          df.kero = dat.i$df.kero)
      #
      if (length(x$method.random.ci) > 1)
        res.i$method <- i
      #
      if (nrow(res) == 0)
        res <- res.i
      else
        res <- rbind(res, res.i)
    }
  }
  
  #
  # (5) Return BLUPs
  #
  meta$data <- NULL
  meta$level <- level
  meta$backtransf <- backtransf
  #
  attr(res, "x") <- meta
  #
  class(res) <- c("blup.meta", "data.frame")
  #
  res
}





#' @rdname blup.meta
#' @method print blup.meta
#' @export

print.blup.meta <- function(x, backtransf = attr(x, "x")$backtransf,
                            #
                            digits = gs("digits"),
                            digits.se = gs("digits.se"),
                            digits.tau2 = gs("digits.tau2"),
                            digits.tau = gs("digits.tau"),
                            big.mark = gs("big.mark"),
                            #
                            se = FALSE,
                            print.tau2 = gs("print.tau2"),
                            print.tau = gs("print.tau"),
                            #
                            details = gs("details"),
                            ...) {
  
  #
  # (1) Check arguments
  #
  chkclass(x, "blup.meta")
  #
  meta <- attr(x, "x")
  #
  sm <- meta$sm
  pscale <- meta$pscale
  irscale <- meta$irscale
  #
  if (is_untransformed(sm))
    backtransf <- TRUE
  #
  chklogical(backtransf)
  chklogical(se)
  if (!isCol(x, "se.blup"))
    se <- FALSE
  #
  print.se <- se
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  chklogical(print.tau2)
  chklogical(print.tau)
  #
  chklogical(details)
  
  #
  # (2) Back-transform and round results
  #
  sm.lab <- smlab(sm, backtransf, pscale, irscale)
  #
  if (backtransf) {
    x$blup <- backtransf(x$blup, sm)
    x$lower <- backtransf(x$lower, sm)
    x$upper <- backtransf(x$upper, sm)
  }
  #
  blup <- round(x$blup, digits)
  lower <- round(x$lower, digits)
  upper <- round(x$upper, digits)
  #
  if (print.se) {
    se.blup <- round(x$se.blup, digits.se)
    #
    if (sm.lab != "")
      se.lab <- paste0("SE(", sm.lab, ")")
    else
      se.lab <- "SE"
  }
  #
  if (!isCol(x, "method"))
    x$method <- meta$method.random.ci
  #
  text.random <- meta$text.random
  if (length(text.random) == 1 & length(unique(x$method)) > 1)
    text.random <- paste0(text.random, "(", unique(x$method), ")")
  #
  by <- !is.null(meta$subgroup)
  
  #
  # (3) Print results
  #
  crtitle(meta)
  #
  j <- 0
  for (i in unique(x$method)) {
    j <- j + 1
    sel.i <- x$method == i
    #
    res.i <-
      cbind(blup = formatN(blup[sel.i], digits, "NA", big.mark = big.mark),
            if (print.se)
              se = formatPT(se.blup, digits = digits.se,
                            big.mark = big.mark,
                            lab = FALSE, lab.NA = "NA"),
            ci = formatCI(formatN(lower[sel.i], digits, "NA", 
                                  big.mark = big.mark),
                          formatN(upper[sel.i], digits, "NA",
                                  big.mark = big.mark))
      )
    #
    colnames(res.i) <- c(sm.lab, if (print.se) se.lab,
                         paste0(100 * meta$level, "%-PI"))
    #
    if (by)
      res.i <- cbind(res.i, subgroup = x$subgroup[sel.i])
    #
    if (print.tau2 & length(unique(x$tau2[sel.i])) > 1)
      res.i <- cbind(res.i, tau2 =
                       rmSpace(formatPT(x$tau2[sel.i], digits = digits.tau2,
                                        big.mark = big.mark,
                                        lab = FALSE, lab.NA = "NA")))
    #
    if (print.tau & length(unique(x$tau2[sel.i])) > 1)
      res.i <- cbind(res.i, tau =
                       rmSpace(formatPT(sqrt(x$tau2[sel.i]),
                                        digits = digits.tau,
                                        big.mark = big.mark,
                                        lab = FALSE, lab.NA = "NA")))
    #
    if (i == "HK" & (by | length(unique(x$df.hakn[sel.i])) > 1))
      res.i <- cbind(res.i, df.hakn = x$df.hakn[sel.i])
    #
    if (i == "KR" & (by | length(unique(x$df.kero[sel.i])) > 1))
      res.i <- cbind(res.i,
                     df.kero = formatN(round(x$df.kero[sel.i], 2),
                                       2, "NA", big.mark = big.mark))
    #
    rownames(res.i) <- x$studlab[sel.i]
    #
    if (length(unique(x$method)) > 1) {
      if (j > 1)
        cat("\n")
      cat(text.random[j], "\n")
    }
    prmatrix(res.i, quote = FALSE, right = TRUE, ...)
  }
  #
  if (details)
    details <-
    catmeth(meta,
            common = FALSE, random = TRUE, prediction = FALSE,
            overall = TRUE, overall.hetstat = TRUE,
            #
            func.transf = NULL,
            backtransf = backtransf, func.backtransf = NULL,
            #
            big.mark = big.mark, digits = digits, digits.tau = digits.tau,
            text.tau = gs("text.tau"), text.tau2 = gs("text.tau2"),
            #
            print.tau2 = print.tau2, print.tau2.ci = FALSE,
            print.tau = print.tau, print.tau.ci = FALSE,
            #
            print.df = !by)
  #
  invisible(details)
}





#' @rdname blup.meta
#' @method estimates blup.meta
#' @export

estimates.blup.meta <- function(x,
                                se = FALSE,
                                backtransf = attr(x, "x")$backtransf,
                                #
                                digits = gs("digits"),
                                digits.se = gs("digits.se"),
                                digits.tau2 = gs("digits.tau2"),
                                digits.tau = gs("digits.tau"),
                                #
                                writexl = !missing(path),
                                path = "estimates_blup.xlsx",
                                overwrite = FALSE,
                                #
                                ...) {
  
  #
  # (1) Check arguments
  #
  chkclass(x, "blup.meta")
  #
  chklogical(se)
  chklogical(backtransf)
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  #
  chklogical(writexl)
  chkchar(path, length = 1)
  chklogical(overwrite)
  #
  meta <- attr(x, "x")
  res <- x
  
  #
  # (2) Back-transform
  #
  if (backtransf) {
    res$blup <- backtransf(res$blup, meta$sm)
    res$lower <- backtransf(res$lower, meta$sm)
    res$upper <- backtransf(res$upper, meta$sm)
  }
  #
  sm.lab <- smlab(meta$sm, backtransf, meta$pscale, meta$irscale)
  #
  if (!any(meta$method.random.ci == "HK"))
    res$df.hakn <- NULL
  #
  if (!any(meta$method.random.ci == "KR"))
    res$df.kero <- NULL
  
  #
  # (3) Create Excel file
  #
  if (writexl) {
    if (!is_installed_package("writexl", stop = FALSE))
      stop(paste0("Package 'writexl' missing.",
                  "\n  ",
                  "Please use the following R command for installation:",
                  "\n  install.packages(\"writexl\")"),
           call. = FALSE)
    #
    res$blup <- round(res$blup, digits = digits)
    res$lower <- round(res$lower, digits = digits)
    res$upper <- round(res$upper, digits = digits)
    #
    if (isCol(res, "tau2"))
      res$tau2 <- round(res$tau2, digits = digits.tau2)
    #
    if (isCol(res, "tau"))
      res$tau <- round(res$tau, digits = digits.tau)
    #
    names(res)[names(res) == "blup"] <- sm.lab
    #
    if (se & isCol(res, "se.blup") &
        !(backtransf & is_relative_effect(meta$sm))) {
      res$se.blup <- round(res$se.blup, digits = digits.se)
      if (sm.lab != "")
        names(res)[names(res) == "se.blup"] <- paste0("SE(", sm.lab, ")")
      else
        names(res)[names(res) == "se.blup"]  <- "SE"
    }
    else
      res$se.blup <- NULL
    #
    if (file.exists(path) & !overwrite)
      warning("File '", path, "' exists. ",
              "Use argument 'overwrite = TRUE' to overwrite file.",
              call. = FALSE)
    else {
      writexl::write_xlsx(res, path = path, col_names = TRUE, ...)
      message(paste0("Extracted information saved in file '", path, "'."))
    }
    #
    return(invisible(NULL))
  }
  
  #
  # (4) Return estimates object
  #
  class(res) <- c("estimates.blup.meta", "data.frame")
  #
  attr(res, "digits") <- digits
  attr(res, "digits.se") <- digits.se
  attr(res, "digits.tau2") <- digits.tau2
  attr(res, "digits.tau") <- digits.tau
  #
  meta$backtransf <- backtransf
  attr(res, "x") <- meta
  #
  res
}





#' @rdname blup.meta
#' @method print estimates.blup.meta
#' @export

print.estimates.blup.meta <- function(x,
                                      big.mark = gs("big.mark"),
                                      details = gs("details"),
                                      ...) {
  
  chkclass(x, "estimates.blup.meta")
  
  meta <- attr(x, "x")
  #
  digits <- attr(x, "digits")
  digits.se <- attr(x, "digits.se")
  digits.tau2 <- attr(x, "digits.tau2")
  digits.tau <- attr(x, "digits.tau")
  #
  backtransf <- meta$backtransf
  #
  se <- isCol(x, "se.blup")
  
  sm.lab <- smlab(meta$sm, backtransf, meta$pscale, meta$irscale)
  #
  ci.lab <- paste0(round(100 * meta$level, 1), "%-PI")
  #
  x$blup <- round(x$blup, digits = digits)
  x$lower <- round(x$lower, digits = digits)
  x$upper <- round(x$upper, digits = digits)
  studlab <- x$studlab
  x$studlab <- NULL
  #
  if (se)
    x$se.blup <- round(x$se.blup, digits = digits.se)
  #
  x$blup <- formatN(x$blup, digits, "", big.mark = big.mark)
  #
  x$lower <-
    ifelse(is.na(x$lower) & is.na(x$upper), "",
           formatCI(
             formatN(x$lower, digits, "", big.mark = big.mark),
             formatN(x$upper, digits, "", big.mark = big.mark)))
  x$upper <- NULL
  #
  if (isCol(x, "tau2"))
    x$tau2 <- round(x$tau2, digits = digits.tau2)
  #
  if (isCol(x, "tau"))
    x$tau <- round(x$tau, digits = digits.tau)
  #
  x <- as.matrix(x)
  rownames(x) <- studlab
  colnames(x)[colnames(x) == "blup"] <- sm.lab
  colnames(x)[colnames(x) == "lower"] <- ci.lab
  if (se) {
    if (sm.lab != "")
      colnames(x)[colnames(x) == "se.blup"] <- paste0("SE(", sm.lab, ")")
    else
      colnames(x)[colnames(x) == "se.blup"] <- "SE"
  }
  
  crtitle(meta)
  #
  prmatrix(x, quote = FALSE, right = TRUE, ...)
  #
  if (details)
    details <-
      catmeth(meta,
              common = FALSE, random = TRUE, prediction = FALSE,
              overall = TRUE, overall.hetstat = TRUE,
              #
              func.transf = NULL,
              backtransf = backtransf, func.backtransf = NULL,
              #
              big.mark = big.mark, digits = digits, digits.tau = digits.tau,
              text.tau = gs("text.tau"), text.tau2 = gs("text.tau2"),
              #
              print.tau2 = TRUE, print.tau2.ci = FALSE,
              print.tau = FALSE, print.tau.ci = FALSE)
  #
  invisible(details)
}
