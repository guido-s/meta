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
  ftr <- atr <- model <- name <- method <- method.interval <- TE <- pval <-
    seTE <- lower <- upper <- adhoc.hakn <- df <- Q <- upper.I2 <- pval.Q <-
    tau2 <- lower.tau2 <- upper.tau2 <- tau <- lower.tau <- upper.tau <-
    I2 <- lower.I2 <- H <- lower.H <- upper.H <-
    level.hetstat <- level.ma <- level.predict <- NULL
  
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
  if (!common)
    x$common <- NULL
  #
  if (!random)
    x$random <- NULL
  #
  if (!prediction)
    x$predict <- NULL
  #
  # Print information
  #
  print.list <-
    merge_legacy(print, args, print_defaults(x, "metabind", func = "print"))
  #
  print <- print.list$res
  args <- print.list$args
  #
  print.k <- print$k
  print.sm <- replaceNULL(print$sm, !x$same.sm)
  print.se <- print$se
  print.pval <- print$pval
  print.stat <- print$stat
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
  chklogical(print.k)
  chklogical(print.sm)
  chklogical(print.se)
  chklogical(print.pval)
  chklogical(print.stat)
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
  # Digits
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
  chkchar(big.mark, length = 1)
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  #
  chklogical(details.methods)
  #
  print.het <- print.tau2 || print.tau || print.I2 || print.H
  #
  if (!common && !random && !prediction && !print.het) {
    warning("Nothing selected to print.", call. = FALSE)
    return(invisible(NULL))
  }
  
  #
  #
  # (2) Create data set
  #
  #
  
  dat <- as.data.frame(x, drop_all_NA = FALSE)
  #
  # Set certain information to missing
  #
  if (common && random) {
    sel.c <- dat$model == "Common effect model"
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
    sel.p <- dat$model == "Prediction interval"
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
  backtransf <- x$backtransf
  func.backtransf <- NULL
  #
  dat$TE.orig <- dat$TE
  dat$seTE.orig <- dat$seTE
  dat$lower.orig <- dat$lower
  dat$upper.orig <- dat$upper
  #
  pscale <- makeunique(x$overall$pscale, 1)
  irscale <- makeunique(x$overall$irscale, 1)
  irunit <- makeunique(x$overall$irunit, "units")
  #
  same.level.ma <- length(unique(x$overall$level.ma)) == 1
  same.level.predict <- length(unique(x$overall$level.predict)) == 1
  same.level.hetstat <- length(unique(x$overall$level.hetstat)) == 1
  #
  if (x$same.sm) {
    sm <- label.sm <- x$sm
    #
    fbt <- makeunique(x$overall$func.backtransf, "")
    is_relative <- is_relative_effect(sm) | (!is.null(fbt) && fbt == "exp")
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
    #
    for (i in seq_len(nrow(dat))) {
      dat$TE[i] <- backtransf(dat$TE[i], sm,
                              n = dat$n.harmonic.mean[i],
                              time = dat$t.harmonic.mean[i],
                              func = dat$func.backtransf[i],
                              args = dat$args.backtransf[i])
      #
      dat$lower[i] <- backtransf(dat$lower[i], sm,
                                 n = dat$n.harmonic.mean[i],
                                 time = dat$t.harmonic.mean[i],
                                 func = dat$func.backtransf[i],
                                 args = dat$args.backtransf[i])
      #
      dat$upper[i] <- backtransf(dat$upper[i], sm,
                                 n = dat$n.harmonic.mean[i],
                                 time = dat$t.harmonic.mean[i],
                                 func = dat$func.backtransf[i],
                                 args = dat$args.backtransf[i])
    }
    #
    if (is_prop(sm) || sm == "RD") {
      dat$TE <- pscale * dat$TE
      dat$lower <- pscale * dat$lower
      dat$upper <- pscale * dat$upper
    }
    #
    else if (is_rate(sm) || sm == "IRD") {
      dat$TE <- irscale * dat$TE
      dat$lower <- irscale * dat$lower
      dat$upper <- irscale * dat$upper
    }
  }
  else {
    if (x$sm.family == "relative"){
      label.sm <- "Effect"
      func.backtransf <- "exp"
    }
    else if (x$sm.family == "correlation") {
      label.sm <- "Correlation"
      backtransf <- FALSE
      #
      for (i in seq_len(nrow(dat))) {
        dat$TE[i] <- backtransf(dat$TE[i], dat$sm[i])
        dat$lower[i] <- backtransf(dat$lower[i], dat$sm[i])
        dat$upper[i] <- backtransf(dat$upper[i], dat$sm[i])
      }
    }
    else if (x$sm.family == "mean") {
      label.sm <- "Mean"
      backtransf <- FALSE
      #
      for (i in seq_len(nrow(dat))) {
        dat$TE[i] <- backtransf(dat$TE[i], dat$sm[i])
        dat$lower[i] <- backtransf(dat$lower[i], dat$sm[i])
        dat$upper[i] <- backtransf(dat$upper[i], dat$sm[i])
      }
    }
    else if (x$sm.family == "proportion") {
      if (pscale == 1)
        label.sm <- "Proportion"
      else
        label.sm <- "Events"
      #
      backtransf <- FALSE
      #
      for (i in seq_len(nrow(dat))) {
        dat$TE[i] <- pscale *
          backtransf(dat$TE[i], dat$sm[i], dat$n.harmonic.mean[i])
        dat$lower[i] <- pscale *
          backtransf(dat$lower[i], dat$sm[i], dat$n.harmonic.mean[i])
        dat$upper[i] <- pscale *
          backtransf(dat$upper[i], dat$sm[i], dat$n.harmonic.mean[i])
      }
    }
    else if (x$sm.family == "rate") {
      if (irscale == 1)
        label.sm <- "Rates"
      else
        label.sm <- "Events"
      #
      backtransf <- FALSE
      #
      for (i in seq_len(nrow(dat))) {
        dat$TE[i] <- irscale *
          backtransf(dat$TE[i], dat$sm[i], dat$n.harmonic.mean[i])
        dat$lower[i] <- irscale *
          backtransf(dat$lower[i], dat$sm[i], dat$n.harmonic.mean[i])
        dat$upper[i] <- irscale *
          backtransf(dat$upper[i], dat$sm[i], dat$n.harmonic.mean[i])
      }
    }
    else
      label.sm <- "Effect"
  }
  #
  # Results for common effect model
  #
  if (common) {
    dat.c <- dat %>% filter(model == "Common effect model") %>%
      select(name, text, sm, method, method.interval, level.ma, TE:pval) %>%
      mutate(
        TE = formatN(TE, digits = digits, lab.NA = lab.NA, big.mark = big.mark),
        seTE = formatN(seTE, digits = digits.se, lab.NA = lab.NA,
                       big.mark = big.mark),
        lower = 
          formatCI(lower, upper, digits = digits, lab.NA = lab.NA,
                   na.lab = TRUE),
        #
        level.ma = round(level.ma, 4)) %>%
      drop_allNAs()
    #
    if (is.null(dat.c$statistic) || allNA(dat.c$statistic))
      dat.c$statistic <- NULL
    else
      dat.c$statistic <-
        formatN(dat.c$statistic, digits = digits.stat, lab.NA = lab.NA)
    #
    if (is.null(dat.c$pval) || allNA(dat.c$pval))
      dat.c$pval <- NULL
    else
      dat.c$pval <- formatPT(dat.c$pval, digits = digits.pval, lab.NA = lab.NA,
                             scientific = scientific.pval,
                             zero = zero.pval, JAMA = JAMA.pval)
    #
    if (all(dat.c$text == "Common effect model"))
      dat.c$text <- NULL
    #
    names(dat.c)[names(dat.c) == "text"] <- ""
    #
    names(dat.c)[names(dat.c) == "TE"] <- label.sm
    #
    if (print.se)
      names(dat.c)[names(dat.c) == "seTE"] <- paste0("SE(", lab.TE, ")")
    else
      dat.c$seTE <- NULL
    #
    if (same.level.ma) {
      names(dat.c)[names(dat.c) == "lower"] <- percent_label(dat.c$level.ma[1])
      dat.c$level.ma <- NULL
    }
    else {
      names(dat.c)[names(dat.c) == "lower"] <- "CI"
      names(dat.c)[names(dat.c) == "level.ma"] <- "level"
    }
    #
    if (print.stat)
      names(dat.c)[names(dat.c) == "statistic"] <- "z"
    else
      dat.c$statistic <- NULL
    #
    if (print.pval)
      names(dat.c)[names(dat.c) == "pval"] <- "p-value"
    else
      dat.c$pval <- NULL
    #
    # Drop additional variables
    #
    dat.c$upper <- NULL
    #
    if (!print.sm)
      dat.c$sm <- NULL
    #
    names(dat.c)[names(dat.c) == "method.interval"] <- "interval"
    #
    if (all(dat.c$interval == "classic"))
      dat.c$interval <- NULL
    #
    if (length(unique(dat.c$text)) == nrow(dat.c))
      dat.c$method <- NULL
    #
    name.c <- dat.c$name
    dat.c$name <- NULL
  }
  #
  # Results for random effects model
  #
  if (random) {
    dat.r <- dat %>% filter(model == "Random effects model") %>%
      select(name, text, sm, method, method.interval, adhoc.hakn, level.ma,
             TE:pval, df) %>%
      mutate(
        TE = formatN(TE, digits = digits, lab.NA = lab.NA, big.mark = big.mark),
        seTE = formatN(seTE, digits = digits.se, lab.NA = lab.NA,
                       big.mark = big.mark),
        lower = 
          formatCI(lower, upper, digits = digits, lab.NA = lab.NA,
                   na.lab = TRUE),
        #
        df = round(df, 2),
        #
        level.ma = round(level.ma, 4)) %>%
      drop_allNAs()
    #
    if (is.null(dat.r$statistic) || allNA(dat.r$statistic))
      dat.r$statistic <- NULL
    else
      dat.r$statistic <-
        formatN(dat.r$statistic, digits = digits.stat, lab.NA = lab.NA)
    #
    if (is.null(dat.r$pval) || allNA(dat.r$pval))
      dat.r$pval <- NULL
    else
      dat.r$pval <- formatPT(dat.r$pval, digits = digits.pval, lab.NA = lab.NA,
                             scientific = scientific.pval,
                             zero = zero.pval, JAMA = JAMA.pval)
    #
    if (all(dat.r$text == "Random effects model"))
      dat.r$text <- NULL
    #
    names(dat.r)[names(dat.r) == "text"] <- ""
    #
    names(dat.r)[names(dat.r) == "TE"] <- label.sm
    #
    if (print.se)
      names(dat.r)[names(dat.r) == "seTE"] <- paste0("SE(", lab.TE, ")")
    else
      dat.r$seTE <- NULL
    #
    if (same.level.ma) {
      names(dat.r)[names(dat.r) == "lower"] <- percent_label(dat.r$level.ma[1])
      dat.r$level.ma <- NULL
    }
    else {
      names(dat.r)[names(dat.r) == "lower"] <- "CI"
      names(dat.r)[names(dat.r) == "level.ma"] <- "level"
    }
    #
    if (all(is.infinite(dat.r$df)))
      dat.r$df <- NULL
    #
    if (print.stat) {
      is.t <- dat.r$method.interval %in% c("HK", "KR", "CR0", "CR1", "CR2")
      #
      if (all(is.t))
        stat.lab <- "t"
      else if (all(!is.t))
        stat.lab <- "z"
      else
        stat.lab <- "z|t"
      #
      names(dat.r)[names(dat.r) == "statistic"] <- stat.lab
    }
    else
      dat.r$statistic <- NULL
    #
    if (print.pval)
      names(dat.r)[names(dat.r) == "pval"] <- "p-value"
    else
      dat.r$pval <- NULL
    #
    # Drop additional variables
    #
    dat.r$upper <- NULL
    #
    if (!print.sm)
      dat.r$sm <- NULL
    #
    names(dat.r)[names(dat.r) == "method.interval"] <- "interval"
    #
    if (all(dat.r$interval == "classic")) {
      dat.r$interval <- NULL
      dat.r$adhoc.hakn <- NULL
    }
    else if (!any(dat.r$interval == "HK"))
      dat.r$adhoc.hakn <- NULL
    #
    if (all(dat.r$adhoc.hakn == "")) {
      dat.r$adhoc.hakn <- NULL
    }
    #
    if (length(unique(dat.r$text)) == nrow(dat.r))
      dat.r$method <- NULL
    #
    name.r <- dat.r$name
    dat.r$name <- NULL
  }
  #
  # Prediction intervals
  #
  if (prediction) {
    dat.p <- dat %>% filter(model == "Prediction interval") %>%
      select(name, text, sm, method, method.interval, adhoc.hakn, level.predict,
             seTE:upper, df) %>%
      mutate(
        lower = 
          formatCI(lower, upper, digits = digits, lab.NA = lab.NA,
                   na.lab = TRUE),
        #
        df = round(df, 2),
        #
        level.predict = round(level.predict, 4)) %>%
      drop_allNAs()
    #
    if (all(dat.p$text == "Prediction interval"))
      dat.p$text <- NULL
    #
    names(dat.p)[names(dat.p) == "text"] <- ""
    #
    if (print.se)
      names(dat.p)[names(dat.p) == "seTE"] <- paste0("SE(", lab.TE, ")")
    else
      dat.p$seTE <- NULL
    #
    if (all(is.infinite(dat.p$df)))
      dat.p$df <- NULL
    #
    if (same.level.predict) {
      names(dat.p)[names(dat.p) == "lower"] <-
        percent_label(dat.p$level.predict[1])
      dat.p$level.predict <- NULL
    }
    else {
      names(dat.p)[names(dat.p) == "lower"] <- "PI"
      names(dat.p)[names(dat.p) == "level.predict"] <- "level"
    }
    #
    # Drop additional variables
    #
    dat.p$upper <- NULL
    #
    if (!print.sm)
      dat.p$sm <- NULL
    #
    names(dat.p)[names(dat.p) == "method.interval"] <- "interval"
    #
    if (all(dat.p$method == "Inverse")) {
      dat.p$method <- NULL
    }
    #
    if (!any(grepl("HK", dat.p$interval)))
      dat.p$adhoc.hakn <- NULL
    #
    if (all(dat.p$adhoc.hakn == "")) {
      dat.p$adhoc.hakn <- NULL
    }
    #
    if (length(unique(dat.p$text)) == nrow(dat.p))
      dat.p$method <- NULL
    #
    name.p <- dat.p$name
    dat.p$name <- NULL
  }
  #
  # Heterogeneity statistics
  #
  if (print.het) {
    dat.h <- dat %>% filter(model == "Random effects model") %>% 
      select(name, text, sm, method, level.hetstat, Q:upper.I2) %>%
      unique() %>%
      mutate(
        Q = formatN(Q, digits = digits.Q, lab.NA = lab.NA,
                    big.mark = big.mark),
        #
        pval.Q = formatPT(pval.Q, digits = digits.pval.Q, lab.NA = lab.NA,
                        scientific = scientific.pval,
                        zero = zero.pval, JAMA = JAMA.pval),
        #
        tau2 = formatPT(tau2, digits = digits.tau2, lab.NA = lab.NA,
                        big.mark = big.mark),
        #
        lower.tau2 = 
          formatCI(lower.tau2, upper.tau2,
                   digits = digits.tau2, lab.NA = lab.NA, na.lab = TRUE,
                   big.mark = big.mark),
        #
        tau = formatPT(tau, digits = digits.tau, lab.NA = lab.NA,
                       big.mark = big.mark),
        #
        lower.tau = 
          formatCI(lower.tau, upper.tau,
                   digits = digits.tau, lab.NA = lab.NA, na.lab = TRUE,
                   big.mark = big.mark),
        #
        I2 = percent_label(I2, "", digits.I2, formatN, lab.NA = lab.NA),
        #
        lower.I2 = 
          formatCI(lower.I2, upper.I2, digits = digits.I2,
                   scale = 100, suffix = "%",
                   lab.NA = lab.NA, na.lab = TRUE),
        #
        H = formatN(H, digits = digits.I2, lab.NA = lab.NA,
                    big.mark = big.mark),
        #
        lower.H = 
          formatCI(lower.H, upper.H,
                   digits = digits.I2, lab.NA = lab.NA, na.lab = TRUE)
      )
    #
    if (all(dat.h$text == "Random effects model"))
      dat.h$text <- NULL
    #
    dat.h$se.tau2 <- NULL
    dat.h$upper.tau2 <- NULL
    dat.h$upper.tau <- NULL
    dat.h$upper.H <- NULL
    dat.h$upper.I2 <- NULL
    #
    if (!print.sm)
      dat.h$sm <- NULL
    #
    if (length(unique(dat.h$method.interval)) == 1)
      dat.h$method.interval <- NULL
    #
    if (length(unique(dat.h$text)) == nrow(dat.h))
      dat.h$method <- NULL
    #
    names(dat.h)[names(dat.h) == "text"] <- ""
    #
    names(dat.h)[names(dat.h) == "TE"] <- label.sm
    #
    if (print.se)
      names(dat.h)[names(dat.h) == "seTE"] <- paste0("SE(", lab.TE, ")")
    else
      dat.h$seTE <- NULL
    #
    if (!print.Q) {
      dat.h$Q <- NULL
      dat.h$df.Q <- NULL
      dat.h$pval.Q <- NULL
    }
    #
    if (print.pval.Q) {
      names(dat.h)[names(dat.h) == "df.Q"] <- "df"
      names(dat.h)[names(dat.h) == "pval.Q"] <- "p-value"
    }
    else {
      dat.h$df.Q <- NULL
      dat.h$pval.Q <- NULL
    }
    #
    if (print.tau2)
      names(dat.h)[names(dat.h) == "tau2"] <- text.tau2
    else {
      dat.h$tau2 <- NULL
      dat.h$lower.tau2 <- NULL
    }
    #
    if (same.level.hetstat) {
      cilab.hetstat <- percent_label(dat.h$level.hetstat[1])
      dat.h$level.hetstat <- NULL
    }
    else
      cilab.hetstat <- "CI"
    #
    if (print.tau2.ci)
      names(dat.h)[names(dat.h) == "lower.tau2"] <- cilab.hetstat
    else
      dat.h$lower.tau2 <- NULL
    #
    if (print.tau)
      names(dat.h)[names(dat.h) == "tau"] <- text.tau
    else
      dat.h$tau <- NULL
    #
    if (print.tau.ci)
      names(dat.h)[names(dat.h) == "lower.tau"] <- cilab.hetstat
    else
      dat.h$lower.tau <- NULL
    #
    if (print.I2)
      names(dat.h)[names(dat.h) == "I2"] <- text.I2
    else
      dat.h$I2 <- NULL
    #
    if (print.I2.ci)
      names(dat.h)[names(dat.h) == "lower.I2"] <- cilab.hetstat
    else
      dat.h$lower.I2 <- NULL
    #
    if (!print.H)
      dat.h$H <- NULL
    #
    if (print.H.ci)
      names(dat.h)[names(dat.h) == "lower.H"] <- cilab.hetstat
    else
      dat.h$lower.H <- NULL
    #
    if (!print.tau2.ci & !print.tau.ci & !print.I2.ci & !print.H.ci)
      dat.h$level.hetstat <- NULL
    else
      names(dat.h)[names(dat.h) == "level.hetstat"] <- "level"
    #
    name.h <- dat.h$name
    dat.h$name <- NULL
  }
  
  #
  #
  # (4) Print results
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
                     backtransf = x$backtransf,
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
