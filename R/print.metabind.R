#' Print results of several meta-analyses
#' 
#' @description
#' Print results of several meta-analyses
#' 
#' @aliases print.metabind
#' 
#' @param x An object of class \code{\link{metabind}}.
#' @param common A logical indicating whether common effect estimates
#'   should be plotted.
#' @param random A logical indicating whether random effects estimates
#'   should be plotted.
#' @param prediction A logical indicating whether prediction intervals
#'   should be plotted.
#' @param overall A logical indicating whether overall summaries
#'   should be plotted. This argument is useful in a meta-analysis
#'   with subgroups if summaries should only be plotted on group
#'   level.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If \code{backtransf=TRUE}
#'   (default), results for \code{sm="OR"} are printed as odds ratios rather
#'   than log odds ratios, for example.
#' @param print A named list specifying which information to print; see
#'   \code{\link{meta-formatting}}.
#' @param digits A named list specifying digits used for displayed results; a
#'   scalar value is also accepted. See \code{\link{meta-formatting}}.
#' @param text.tau2 Text printed to identify between-study variance
#'   \eqn{\tau^2}.
#' @param text.tau Text printed to identify \eqn{\tau}, the square
#'   root of the between-study variance \eqn{\tau^2}.
#' @param text.I2 Text printed to identify heterogeneity statistic
#'   I\eqn{^2}.
#' @param lab.NA A character string to label missing values.
#' @param big.mark A character used as thousands separator.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   overall effect should be printed according to JAMA reporting
#'   standards.
#' @param details.methods A logical specifying whether details on
#'   statistical methods should be printed.
#' @param \dots Additional arguments.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metabind}}, \code{\link{settings.meta}}
#' 
#' @examples
#' data(Fleiss1993bin)
#' ma1 <- metabin(d.asp, n.asp, d.plac, n.plac,
#'   data = Fleiss1993bin, studlab = study, sm = "RR", method = "I")
#' ma2 <- update(ma1, method = "MH", random = FALSE)
#' ma1
#' ma2
#' metabind(ma1, ma2)
#' 
#' @method print metabind
#' @export

print.metabind <- function(x,
                          #
                          common = !is.null(x$common),
                          random = !is.null(x$random),
                          prediction = !is.null(x$predict),
                          #
                          overall = common || random,
                          #
                          backtransf = x$backtransf,
                          #
                          print = NULL,
                          digits = NULL,
                          #
                          lab.NA = ".",
                          #
                          big.mark = gs("big.mark"),
                          scientific.pval = gs("scientific.pval"),
                          zero.pval = gs("zero.pval"),
                          JAMA.pval = gs("JAMA.pval"),
                          #
                          text.tau2 = gs("text.tau2"),
                          text.tau = gs("text.tau"),
                          text.I2 = gs("text.I2"),
                          #
                           details.methods = gs("details"),
                           ...) {
  # Get rid of warning 'no visible binding for global variable'
  ftr <- atr <- type <- name <- method <- method.interval <- TE <- pval <-
    seTE <- lower <- upper <- adhoc.hakn <- df <- Q <- upper.I2 <- pval.Q <-
    tau2 <- lower.tau2 <- upper.tau2 <- tau <- lower.tau <- upper.tau <-
    I2 <- lower.I2 <- H <- lower.H <- upper.H <-
    level.hetstat <- level.ma <- level.predict <- statistic <-
    id <- method.tau <- method.tau.ci <- method.I2 <- se.tau2 <- NULL
  
  #
  #
  # (1) Check and set arguments
  #
  #
  
  chkclass(x, "metabind")
  x <- updateversion(x)
  #
  # Drop internally used arguments from '...'
  #
  drop_args <- NULL
  #
  replace_args <- list()
  #
  args <- drop_from_dots(list(...), drop_args, replace_args)
  #
  with.subgroups <- !is.null(x$subgroup)
  #
  if (with.subgroups)
    return(invisible(NULL))
  #
  if (!with.subgroups) {
    overall <- FALSE
    overall.hetstat <- FALSE
  }
  #
  overall <- replaceNULL(overall, FALSE)
  overall.hetstat <- replaceNULL(overall.hetstat, FALSE)
  #
  missing.common <- missing(common)
  missing.random <- missing(random)
  missing.prediction <- missing(prediction)
  #
  chklogical(common)
  chklogical(random)
  chklogical(prediction)
  #
  if (is.null(x$common)) {
    if (!missing.common && common)
      warning("Argument 'common' set to FALSE as ",
              "no common effect results are available.",
              call. = FALSE)
    #
    common <- FALSE
  }
  #
  if (is.null(x$random)) {
    if (!missing.random && random)
      warning("Argument 'random' set to FALSE as ",
              "no random effects results are available.",
              call. = FALSE)
    #
    random <- FALSE
  }
  #
  if (is.null(x$predict)) {
    if (!missing.prediction && prediction)
      warning("Argument 'prediction' set to FALSE as ",
              "no predict intervals are available.",
              call. = FALSE)
    #
    prediction <- FALSE
  }
  #
  chklogical(overall)
  chklogical(overall.hetstat)
  #
  chklogical(backtransf)
  #
  if (!common)
    x$common <- NULL
  #
  if (!random)
    x$random <- NULL
  #
  if (!prediction)
    x$predict <- NULL
  #
  chkchar(big.mark, length = 1)
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  #
  chklogical(details.methods)
  #
  # Inputs related to printing
  #
  print.list <-
    merge_legacy(print, args, print_defaults(x, "metabind", func = "print"))
  #
  print <- print.list$res
  args <- print.list$args
  #
  print.class <- print$class
  print.method <- print$method
  print.method.interval <- print$method.interval
  print.k <- print$k
  print.sm <- print$sm
  print.se <- print$se
  print.pval <- print$pval
  print.stat <- print$stat
  print.level.ma <- print$level.ma
  print.level.predict <- print$level.predict
  print.level.hetstat <- print$level.hetstat
  print.null.effect <- print$null.effect
  print.I2 <- print$I2
  print.I2.ci <- print$I2.ci
  print.tau2 <- print$tau2
  print.tau2.ci <- print$tau2.ci
  print.tau <- print$tau
  print.tau.ci <- print$tau.ci
  print.Q <- print$Q
  print.pval.Q <- replaceNULL(print$pval.Q, print.Q) && print.Q
  print.H <- print$H
  print.H.ci <- print$H.ci
  print.pval.Q.b <- print$pval.Q.b
  print.subgroup.name <- print$subgroup.name
  #
  print.pval.Q.b <- replaceNULL(print.pval.Q.b, with.subgroups)
  print.tau2 <- replaceNULL(print.tau2, !with.subgroups && random)
  print.tau <- replaceNULL(print.tau, !with.subgroups && random)
  print.subgroup.name <- replaceNULL(print.subgroup.name, with.subgroups)
  #
  chklogical(print.class, NULL.ok = TRUE)
  chklogical(print.method, NULL.ok = TRUE)
  chklogical(print.method.interval, NULL.ok = TRUE)
  chklogical(print.k)
  chklogical(print.sm, NULL.ok = TRUE)
  chklogical(print.se)
  chklogical(print.pval)
  chklogical(print.stat)
  chklogical(print.level.ma, NULL.ok = TRUE)
  chklogical(print.level.predict, NULL.ok = TRUE)
  chklogical(print.level.hetstat, NULL.ok = TRUE)
  chklogical(print.null.effect)
  chklogical(print.I2)
  chklogical(print.I2.ci)
  chklogical(print.tau2)
  chklogical(print.tau2.ci)
  chklogical(print.tau)
  chklogical(print.tau.ci)
  chklogical(print.Q)
  chklogical(print.pval.Q)
  chklogical(print.H)
  chklogical(print.H.ci)
  chklogical(print.pval.Q.b)
  chklogical(print.subgroup.name)
  #
  print.tau2.ci <- print.tau2 && print.tau2.ci
  print.tau.ci <- print.tau && print.tau.ci
  print.I2.ci <- print.I2 && print.I2.ci
  print.H.ci <- print.H && print.H.ci
  #
  print.het <- print.Q || print.tau2 || print.tau || print.I2 || print.H
  #
  if (!common && !random && !prediction && !print.het) {
    warning("Nothing selected to print.", call. = FALSE)
    return(invisible(NULL))
  }
  #
  # Inputs related to rounding of results
  #
  digits <- as_named_list(digits)
  digits.list <-
    merge_legacy(digits, args, digits_defaults("metabind", func = "print"))
  #
  digitss <- digits.list$res
  args <- digits.list$args
  #
  digits <- digitss$digits
  digits.se <- digitss$se
  digits.stat <- digitss$stat
  digits.pval <- digitss$pval
  digits.pval.Q <- digitss$pval.Q
  digits.Q <- digitss$Q
  digits.tau2 <- digitss$tau2
  digits.tau <- digitss$tau
  digits.I2 <- digitss$I2
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  chknumeric(digits.pval.Q, min = 1, length = 1)
  chknumeric(digits.Q, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  
  #
  #
  # (2) Create column labels and more general settings
  #
  #
  
  pscale <- makeunique(x$overall$pscale, 1)
  irscale <- makeunique(x$overall$irscale, 1)
  irunit <- makeunique(x$overall$irunit, "units")
  #
  sm <- unique(x$overall$sm)
  same.sm <- length(sm) == 1
  #
  if (!same.sm && !backtransf) {
    warning("Argument 'backtransf' set to TRUE as ",
            "meta-analyses have different summary measures.",
            call. = FALSE)
    backtransf <- TRUE
  }
  if (same.sm) {
    sm <- sm[1]
    #
    fbt <- makeunique(x$overall$func.backtransf, "")
    is_relative <- is_relative_effect(sm) | (!is.null(fbt) && fbt == "exp")
    #
    sm.lab <- smlab(
      sm,
      backtransf = backtransf,
      pscale = pscale,
      irscale = irscale,
      func.backtransf = fbt, fallback = ".no_sm_provided")
    #
    lab.TE <- sm
    #
    if (is_relative || sm == "VE")
      lab.TE <- paste0(gs("log.prefix"), if (sm == "VE") "VR" else sm)
    else if (!is.null(ftr)) {
      lab.TE <-
        paste0(ftr, "(", sm,
               if (!is.null(atr) && length(names(atr)) >= 1)
                 paste0(", ",
                        paste0(names(atr), "=", paste(atr), collapse = ", ")),
               ")")
    }
    else if (sm == "")
      lab.TE <- "TE"
  }
  else {
    if (x$sm.family == "relative"){
      sm.lab <- "Effect"
    }
    else if (x$sm.family == "correlation") {
      sm.lab <- "Correlation"
    }
    else if (x$sm.family == "mean") {
      sm.lab <- "Mean"
    }
    else if (x$sm.family == "proportion") {
      if (pscale == 1)
        sm.lab <- "Proportion"
      else
        sm.lab <- "Events"
    }
    else if (x$sm.family == "rate") {
      if (irscale == 1)
        sm.lab <- "Rate"
      else
        sm.lab <- "Events"
    }
    else
      sm.lab <- "Effect"
    #
    lab.TE <- "TE"
  }
    
  #
  #
  # (3) Create data sets with results
  #
  #
  
  dat <- NULL
  #
  if (common || random || prediction) {
    dat <- as.data.frame(x, drop_all_NA = FALSE)
    #
    # Set certain information to missing
    #
    if (common && random) {
      sel.c <- dat$type == "Common effect model"
      #
      dat$tau2[sel.c] <- NA
      dat$lower.tau2[sel.c] <- NA
      dat$upper.tau2[sel.c] <- NA
      #
      dat$tau[sel.c] <- NA
      dat$lower.tau[sel.c] <- NA
      dat$upper.tau[sel.c] <- NA
      #
      dat$I2[sel.c] <- NA
      dat$lower.I2[sel.c] <- NA
      dat$upper.I2[sel.c] <- NA
    }
    #
    if (prediction && (common || random)) {
      sel.p <- dat$type == "Prediction interval"
      #
      dat$tau2[sel.p] <- NA
      dat$lower.tau2[sel.p] <- NA
      dat$upper.tau2[sel.p] <- NA
      #
      dat$tau[sel.p] <- NA
      dat$lower.tau[sel.p] <- NA
      dat$upper.tau[sel.p] <- NA
      #
      dat$I2[sel.p] <- NA
      dat$lower.I2[sel.p] <- NA
      dat$upper.I2[sel.p] <- NA
    }
    #
    # Act upon different summary measures
    #
    dat$TE.orig <- dat$TE
    dat$lower.orig <- dat$lower
    dat$upper.orig <- dat$upper
    #
    # Back transformation
    #
    if (backtransf) {
      for (i in seq_len(nrow(dat))) {
        if (is_prop(dat$sm[i]) || dat$sm[i] == "RD") {
          dat$TE[i] <- pscale *
            backtransf(dat$TE[i], dat$sm[i],
                       n = dat$n.harmonic.mean[i],
                       func = dat$func.backtransf[i],
                       args = dat$args.backtransf[i])
          #
          dat$lower[i] <- pscale *
            backtransf(dat$lower[i], dat$sm[i],
                       n = dat$n.harmonic.mean[i],
                       func = dat$func.backtransf[i],
                       args = dat$args.backtransf[i])
          #
          dat$upper[i] <- pscale *
            backtransf(dat$upper[i], dat$sm[i],
                       n = dat$n.harmonic.mean[i],
                       func = dat$func.backtransf[i],
                       args = dat$args.backtransf[i])
        }
        #
        else if (is_rate(dat$sm[i]) || dat$sm[i] == "IRD") {
          dat$TE[i] <- irscale *
            backtransf(dat$TE[i], dat$sm[i],
                       time = dat$t.harmonic.mean[i],
                       func = dat$func.backtransf[i],
                       args = dat$args.backtransf[i])
          #
          dat$lower[i] <- irscale *
            backtransf(dat$lower[i], dat$sm[i],
                       time = dat$t.harmonic.mean[i],
                       func = dat$func.backtransf[i],
                       args = dat$args.backtransf[i])
          #
          dat$upper[i] <- irscale *
            backtransf(dat$upper[i], dat$sm[i],
                       time = dat$t.harmonic.mean[i],
                       func = dat$func.backtransf[i],
                       args = dat$args.backtransf[i])
        }
        else {
          dat$TE[i] <- backtransf(dat$TE[i], dat$sm[i],
                                  n = dat$n.harmonic.mean[i],
                                  time = dat$t.harmonic.mean[i],
                                  func = dat$func.backtransf[i],
                                  args = dat$args.backtransf[i])
          #
          dat$lower[i] <- backtransf(dat$lower[i], dat$sm[i],
                                     n = dat$n.harmonic.mean[i],
                                     time = dat$t.harmonic.mean[i],
                                     func = dat$func.backtransf[i],
                                     args = dat$args.backtransf[i])
          #
          dat$upper[i] <- backtransf(dat$upper[i], dat$sm[i],
                                     n = dat$n.harmonic.mean[i],
                                     time = dat$t.harmonic.mean[i],
                                     func = dat$func.backtransf[i],
                                     args = dat$args.backtransf[i])
        }
      }
    }
    #
    # Results for common effect model
    #
    if (common) {
      dat.c <- dat %>%
        filter(type == "Common effect model") %>%
        select(name, class, sm, method, method.interval, level.ma, TE:pval) %>%
        #
        mutate(
          TE = formatN(TE, digits = digits, lab.NA = lab.NA,
                       big.mark = big.mark),
          seTE = formatN(seTE, digits = digits.se, lab.NA = lab.NA,
                         big.mark = big.mark, keep.allNA = TRUE),
          #
          lower = formatCI(lower, upper, digits = digits, lab.NA = lab.NA,
                           na.lab = TRUE),
          #
          statistic = formatN(statistic, digits = digits.stat, lab.NA = lab.NA,
                              keep.allNA = TRUE),
          #
          pval = formatPT(pval, digits = digits.pval, lab.NA = lab.NA,
                          scientific = scientific.pval,
                          zero = zero.pval, JAMA = JAMA.pval,
                          keep.allNA = TRUE),
          #
          level.ma = round(level.ma, 4)) %>%
        #
        select(-upper) %>%
        #
        rename(
          !!sm.lab := TE,
          interval = method.interval,
          level = level.ma) %>%
        #
        drop_allNAs()
      #
      same.level.ma.c <- length(unique(dat.c$level)) == 1
      print.level.ma.c <- replaceNULL(print.level.ma, !same.level.ma.c)
      #
      if (same.level.ma.c)
        names(dat.c)[names(dat.c) == "lower"] <-
        percent_label(dat.c$level[1])
      else
        names(dat.c)[names(dat.c) == "lower"] <- "CI"
      #
      if (!print.level.ma.c)
        dat.c$level <- NULL
      #
      if (isFALSE(print.class) ||
          (is.null(print.class) && length(unique(dat.c$class)) == 1))
        dat.c$class <- NULL
      #
      if (isFALSE(print.sm) ||
          (is.null(print.sm) && length(unique(dat.c$sm)) == 1))
        dat.c$sm <- NULL
      #
      if (isFALSE(print.method) ||
          (is.null(print.method) && length(unique(dat.c$method)) == 1))
        dat.c$method <- NULL
      #
      if (isFALSE(print.method.interval) ||
          (is.null(print.method.interval) &&
           length(unique(dat.c$interval)) == 1))
        dat.c$interval <- NULL
      #
      if (!is.null(dat.c$statistic) && print.stat)
        names(dat.c)[names(dat.c) == "statistic"] <- "z"
      else
        dat.c$statistic <- NULL
      #
      if (!is.null(dat.c$pval) && print.pval)
        names(dat.c)[names(dat.c) == "pval"] <- "P value"
      else
        dat.c$pval <- NULL
      #
      if (print.se)
        names(dat.c)[names(dat.c) == "seTE"] <- paste0("SE(", lab.TE, ")")
      else
        dat.c$seTE <- NULL
      #
      names(dat.c)[names(dat.c) == ".no_sm_provided"] <- ""
      #
      name.c <- dat.c$name
      dat.c$name <- NULL
    }
    #
    # Results for random effects model
    #
    if (random) {
      dat.r <- dat %>%
        filter(type == "Random effects model") %>%
        select(name, class, sm, method, method.interval,
               adhoc.hakn, level.ma, TE:pval, df) %>%
        #
        mutate(
          TE = formatN(TE, digits = digits, lab.NA = lab.NA,
                       big.mark = big.mark),
          seTE = formatN(seTE, digits = digits.se, lab.NA = lab.NA,
                         big.mark = big.mark, keep.allNA = TRUE),
          #
          lower = formatCI(lower, upper, digits = digits, lab.NA = lab.NA,
                           na.lab = TRUE),
          #
          statistic = formatN(statistic, digits = digits.stat, lab.NA = lab.NA,
                              keep.allNA = TRUE),
          #
          pval = formatPT(pval, digits = digits.pval, lab.NA = lab.NA,
                          scientific = scientific.pval,
                          zero = zero.pval, JAMA = JAMA.pval,
                          keep.allNA = TRUE),
          #
          df = round(df, 2),
          level.ma = round(level.ma, 4)) %>%
        #
        select(-upper) %>%
        #
        rename(
          !!sm.lab := TE,
          interval = method.interval,
          level = level.ma) %>%
        #
        drop_allNAs()
      #
      interval.r <- dat.r$interval
      same.level.ma.r <- length(unique(dat.r$level)) == 1
      print.level.ma.r <- replaceNULL(print.level.ma, !same.level.ma.r)
      #
      if (same.level.ma.r)
        names(dat.r)[names(dat.r) == "lower"] <-
        percent_label(dat.r$level[1])
      else
        names(dat.r)[names(dat.r) == "lower"] <- "CI"
      #
      if (!print.level.ma.r)
        dat.r$level <- NULL
      #
      if (isFALSE(print.class) ||
          (is.null(print.class) && length(unique(dat.r$class)) == 1))
        dat.r$class <- NULL
      #
      if (isFALSE(print.sm) ||
          (is.null(print.sm) && length(unique(dat.r$sm)) == 1))
        dat.r$sm <- NULL
      #
      if (all(is.infinite(dat.r$df)))
        dat.r$df <- NULL
      #
      if (isFALSE(print.method) ||
          (is.null(print.method) && length(unique(dat.r$method)) == 1))
        dat.r$method <- NULL
      #
      if (isFALSE(print.method.interval) ||
          (is.null(print.method.interval) &&
           length(unique(interval.r)) == 1))
        dat.r$interval <- NULL
      #
      if (is.null(dat.r$interval) ||
          !any(interval.r %in% "HK") ||
          all(is.na(dat.r$adhoc.hakn) | dat.r$adhoc.hakn == ""))
        dat.r$adhoc.hakn <- NULL
      #
      if (!is.null(dat.r$statistic) && print.stat) {
        is.t <- interval.r %in% c("HK", "KR", "CR0", "CR1", "CR2")
        #
        stat.lab <-
          if (all(is.t))
            "t"
        else if (all(!is.t))
          "z"
        else
          "z|t"
        #
        names(dat.r)[names(dat.r) == "statistic"] <- stat.lab
      }
      else
        dat.r$statistic <- NULL
      #
      if (!is.null(dat.r$pval) && print.pval)
        names(dat.r)[names(dat.r) == "pval"] <- "P value"
      else
        dat.r$pval <- NULL
      #
      if (print.se)
        names(dat.r)[names(dat.r) == "seTE"] <- paste0("SE(", lab.TE, ")")
      else
        dat.r$seTE <- NULL
      #
      names(dat.r)[names(dat.r) == ".no_sm_provided"] <- ""
      #
      name.r <- dat.r$name
      dat.r$name <- NULL
    }
    #
    # Prediction intervals
    #
    if (prediction) {
      dat.p <- dat %>%
        filter(type == "Prediction interval") %>%
        select(name, class, sm, method, method.interval,
               adhoc.hakn, level.predict, seTE:upper, df) %>%
        #
        mutate(
          seTE = formatN(seTE, digits = digits.se, lab.NA = lab.NA,
                         big.mark = big.mark, keep.allNA = TRUE),
          #
          lower = formatCI(lower, upper, digits = digits, lab.NA = lab.NA,
                           na.lab = TRUE),
          #
          df = round(df, 2),
          level.predict = round(level.predict, 4)) %>%
        #
        select(-upper) %>%
        #
        rename(
          interval = method.interval,
          level = level.predict) %>%
        #
        drop_allNAs()
      #
      same.level.predict <- length(unique(dat.p$level)) == 1
      print.level.predict <- replaceNULL(print.level.predict, !same.level.predict)
      #
      if (same.level.predict)
        names(dat.p)[names(dat.p) == "lower"] <-
        percent_label(dat.p$level[1], type = "PI")
      else
        names(dat.p)[names(dat.p) == "lower"] <- "PI"
      #
      if (!print.level.predict)
        dat.p$level <- NULL
      #
      if (isFALSE(print.class) ||
          (is.null(print.class) && length(unique(dat.p$class)) == 1))
        dat.p$class <- NULL
      #
      if (isFALSE(print.sm) ||
          (is.null(print.sm) && length(unique(dat.p$sm)) == 1))
        dat.p$sm <- NULL
      #
      if (all(is.infinite(dat.p$df)))
        dat.p$df <- NULL
      #
      if (isFALSE(print.method) ||
          (is.null(print.method) && (all(dat.p$method %in% "Inverse"))))
        dat.p$method <- NULL
      #
      if (isFALSE(print.method.interval) ||
          (is.null(print.method.interval) &&
           length(unique(dat.p$interval)) == 1))
        dat.p$interval <- NULL
      #
      if (is.null(dat.p$interval) ||
          !any(dat.p$interval %in% c("HK", "HK-PR")) ||
          all(is.na(dat.p$adhoc.hakn) | dat.p$adhoc.hakn == ""))
        dat.p$adhoc.hakn <- NULL
      #
      if (print.se)
        names(dat.p)[names(dat.p) == "seTE"] <- paste0("SE(", lab.TE, ")")
      else
        dat.p$seTE <- NULL
      #
      name.p <- dat.p$name
      dat.p$name <- NULL
    }
  }
  
  #
  #
  # (4) Create data set with heterogeneity statistics
  #
  #
  
  if (print.het) {
    if (common || random || prediction) {
      ids.het <-
        unique(c(
          if (common) x$common$id,
          if (random) x$random$id,
          if (prediction) x$predict$id))
    }
    else
      ids.het <- x$overall$id
    #
    dat.h <- x$overall %>%
      #
      filter(id %in% ids.het) %>%
      #
      select(name, class, sm, method.tau, method.tau.ci, method.I2,
             level.hetstat, Q:upper.I2) %>%
      #
      mutate(
        Q = formatN(Q, digits = digits.Q, lab.NA = lab.NA,
                    big.mark = big.mark, keep.allNA = TRUE),
        #
        pval.Q = formatPT(pval.Q, digits = digits.pval.Q, lab.NA = lab.NA,
                          scientific = scientific.pval,
                          zero = zero.pval, JAMA = JAMA.pval,
                          keep.allNA = TRUE),
        #
        tau2 = formatPT(tau2, digits = digits.tau2, lab.NA = lab.NA,
                        big.mark = big.mark, keep.allNA = TRUE),
        #
        lower.tau2 =
          if (allNA(lower.tau2) && allNA(upper.tau2))
            lower.tau2
        else
          formatCI(lower.tau2, upper.tau2,
                   digits = digits.tau2, lab.NA = lab.NA,
                   na.lab = TRUE, big.mark = big.mark),
        #
        tau = formatPT(tau, digits = digits.tau, lab.NA = lab.NA,
                       big.mark = big.mark, keep.allNA = TRUE),
        #
        lower.tau =
          if (allNA(lower.tau) && allNA(upper.tau))
            lower.tau
        else
          formatCI(lower.tau, upper.tau,
                   digits = digits.tau, lab.NA = lab.NA,
                   na.lab = TRUE, big.mark = big.mark),
        #
        I2 =
          if (allNA(I2))
            I2
        else
          percent_label(I2, "", digits.I2, formatN, lab.NA = lab.NA),
        #
        lower.I2 =
          if (allNA(lower.I2) && allNA(upper.I2))
            lower.I2
        else
          formatCI(lower.I2, upper.I2, digits = digits.I2,
                   scale = 100, suffix = "%",
                   lab.NA = lab.NA, na.lab = TRUE),
        #
        H = formatN(H, digits = digits.I2, lab.NA = lab.NA,
                    big.mark = big.mark, keep.allNA = TRUE),
        #
        lower.H =
          if (allNA(lower.H) && allNA(upper.H))
            lower.H
        else
          formatCI(lower.H, upper.H,
                   digits = digits.I2, lab.NA = lab.NA, na.lab = TRUE),
        #
        level.hetstat = round(level.hetstat, 4)) %>%
      #
      select(-se.tau2, -upper.tau2, -upper.tau, -upper.H, -upper.I2) %>%
      #
      rename(level = level.hetstat) %>%
      #
      drop_allNAs()
    #
    same.level.hetstat <- length(unique(dat.h$level)) == 1
    print.level.hetstat.h <-
      replaceNULL(print.level.hetstat, !same.level.hetstat)
    #
    if (same.level.hetstat)
      cilab.hetstat <- percent_label(dat.h$level[1])
    else
      cilab.hetstat <- "CI"
    #
    if (!print.level.hetstat.h)
      dat.h$level <- NULL
    #
    if (isFALSE(print.class) ||
        (is.null(print.class) && length(unique(dat.h$class)) == 1))
      dat.h$class <- NULL
    #
    if (isFALSE(print.sm) ||
        (is.null(print.sm) && length(unique(dat.h$sm)) == 1))
      dat.h$sm <- NULL
    #
    for (method.col in c("method.tau", "method.tau.ci", "method.I2")) {
      if (isFALSE(print.method) ||
          (is.null(print.method) &&
           length(unique(dat.h[[method.col]])) == 1))
        dat.h[[method.col]] <- NULL
    }
    #
    if (!(print.tau2 || print.tau))
      dat.h$method.tau <- NULL
    #
    if (!((print.tau2 && print.tau2.ci) ||
          (print.tau && print.tau.ci)))
      dat.h$method.tau.ci <- NULL
    #
    if (!print.I2)
      dat.h$method.I2 <- NULL
    #
    if (!print.Q)
      dat.h$Q <- NULL
    #
    if (print.Q && print.pval.Q) {
      names(dat.h)[names(dat.h) == "df.Q"] <- "df"
      names(dat.h)[names(dat.h) == "pval.Q"] <- "P value"
    }
    else {
      dat.h$df.Q <- NULL
      dat.h$pval.Q <- NULL
    }
    #
    if (print.tau2)
      names(dat.h)[names(dat.h) == "tau2"] <- text.tau2
    else
      dat.h$tau2 <- NULL
    #
    if (print.tau2 && print.tau2.ci)
      names(dat.h)[names(dat.h) == "lower.tau2"] <- cilab.hetstat
    else
      dat.h$lower.tau2 <- NULL
    #
    if (print.tau)
      names(dat.h)[names(dat.h) == "tau"] <- text.tau
    else
      dat.h$tau <- NULL
    #
    if (print.tau && print.tau.ci)
      names(dat.h)[names(dat.h) == "lower.tau"] <- cilab.hetstat
    else
      dat.h$lower.tau <- NULL
    #
    if (print.I2)
      names(dat.h)[names(dat.h) == "I2"] <- text.I2
    else
      dat.h$I2 <- NULL
    #
    if (print.I2 && print.I2.ci)
      names(dat.h)[names(dat.h) == "lower.I2"] <- cilab.hetstat
    else
      dat.h$lower.I2 <- NULL
    #
    if (!print.H)
      dat.h$H <- NULL
    #
    if (print.H && print.H.ci)
      names(dat.h)[names(dat.h) == "lower.H"] <- cilab.hetstat
    else
      dat.h$lower.H <- NULL
    #
    if (!(print.tau2 && print.tau2.ci) &&
        !(print.tau && print.tau.ci) &&
        !(print.I2 && print.I2.ci) &&
        !(print.H && print.H.ci))
      dat.h$level <- NULL
    #
    name.h <- dat.h$name
    dat.h$name <- NULL
  }
  
  #
  #
  # (5) Print results
  #
  #
  
  cat("Meta-analyses results\n")
  #
  if (common) {
    cat("\nCommon effect model\n\n")
    #
    prmatrix(dat.c, quote = FALSE, right = TRUE, rowlab = name.c)
  }
  #
  if (random) {
    cat("\nRandom effects model\n\n")
    #
    prmatrix(dat.r, quote = FALSE, right = TRUE, rowlab = name.r)
  }
  #
  if (prediction) {
    cat("\nPrediction interval\n\n")
    #
    prmatrix(dat.p, quote = FALSE, right = TRUE, rowlab = name.p)
  }
  #
  if (print.het) {
    cat("\nInformation on heterogeneity\n\n")
    #
    prmatrix(dat.h, quote = FALSE, right = TRUE, rowlab = name.h)
  }
  #
  if (details.methods) {
    catmeth_metabind(x,
                     #
                     common, random, prediction,
                     #
                     print.tau2 = print.tau2,
                     print.tau2.ci = print.tau2.ci,
                     print.tau = print.tau,
                     print.tau.ci = print.tau.ci,
                     print.I2 = print.I2,
                     print.pval = print.pval,
                     #
                     backtransf = backtransf,
                     #
                     digits = digits,
                     digits.tau2 = digits.tau2,
                     digits.tau = digits.tau,
                     #
                     pscale = pscale,
                     irscale = irscale, irunit = irunit,
                     #
                     big.mark = big.mark,
                     #
                     forest = FALSE)
  }
  #
  invisible(NULL)
}
