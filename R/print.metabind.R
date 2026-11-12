#' Print results of several meta-analyses
#'
#' @description
#' Prints common effect estimates, random effects estimates, prediction
#' intervals, and heterogeneity information in separate tables.
#'
#' @aliases print.metabind
#'
#' @param x An object of class \code{\link{metabind}}.
#' @param common,random,prediction Logical values selecting common effect
#'   estimates, random effects estimates, and prediction intervals to print.
#' @param overall A logical indicating whether overall summaries should be
#'   printed. No additional overall summaries are printed for component
#'   results without source subgroup analyses.
#'   For source subgroup analyses, selected overall effects can be shared
#'   independently of heterogeneity. Shared heterogeneity is printed only
#'   when its values and settings also agree across components.
#' @param backtransf A logical indicating whether results should be
#'   back-transformed. Different summary measures require back-transformation.
#' @param print A named list controlling displayed information; see
#'   \code{\link{meta-formatting}}. Elements \code{func}, \code{sm},
#'   \code{method}, and \code{method.interval} accept \code{TRUE},
#'   \code{FALSE}, or \code{NULL} for automatic column selection within
#'   each result table. Element \code{func} identifies the originating
#'   meta-analysis function. Elements \code{level.ma}, \code{level.predict},
#'   and \code{level.hetstat} control confidence-level columns, which are
#'   shown automatically when the corresponding levels differ. Elements
#'   \code{backtransf} and \code{args.backtransf} control columns showing
#'   custom back-transformation functions and their named arguments, with
#'   headings \code{backtransf} and \code{args}. By default, these are
#'   shown when supplied; arguments are shown only with the function column.
#' @param digits A named list controlling rounding; see
#'   \code{\link{meta-formatting}}. A scalar value is also accepted.
#' @param text.tau2 Text printed to identify between-study variance
#'   \eqn{\tau^2}.
#' @param text.tau Text printed to identify \eqn{\tau}, the square
#'   root of the between-study variance \eqn{\tau^2}.
#' @param text.I2 Text printed to identify heterogeneity statistic
#'   I\eqn{^2}.
#' @param lab.NA Label for missing values.
#' @param big.mark Character used as thousands separator.
#' @param scientific.pval,zero.pval,JAMA.pval Logical values controlling
#'   the formatting of P values.
#' @param details.methods A logical specifying whether details on
#'   statistical methods should be printed.
#' @param \dots Additional formatting arguments, including individual
#'   \code{print.*} and \code{digits.*} arguments.
#'
#' @return Invisibly returns \code{NULL}.
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
  ftr <- atr <- type <- name <- func <- method <- method.interval <- TE <-
    pval <- seTE <- lower <- upper <- adhoc.hakn <- df <- Q <- upper.I2 <-
    pval.Q <- tau2 <- lower.tau2 <- upper.tau2 <-
    tau <- lower.tau <- upper.tau <- I2 <- lower.I2 <-
    H <- lower.H <- upper.H <- level.hetstat <- level.ma <- level.predict <-
    statistic <- id <- method.tau <- method.tau.ci <- method.I2 <- se.tau2 <-
    func.backtransf <- args.backtransf <- subgroup.levels <- NULL
  
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
      warning(
        "Argument 'random' set to FALSE as ",
        "no random effects results are available.",
        call. = FALSE
      )
    #
    random <- FALSE
  }
  #
  if (is.null(x$predict)) {
    if (!missing.prediction && prediction)
      warning(
        "Argument 'prediction' set to FALSE as ",
        "no predict intervals are available.",
        call. = FALSE
      )
    #
    prediction <- FALSE
  }
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
  with.subgroups <- !is.null(x$subgroup)
  #
  if (common || random || prediction) {
    dat <- as.data.frame(x, subgroup = FALSE, drop_all_NA = FALSE)
    #
    if (with.subgroups) {
      prepared <- prepare_metabind_subgroups(x)
      #
      if (missing(overall)) {
        overall <- prepared$shared.overall
        overall.hetstat <- prepared$shared.overall &&
          prepared$shared.heterogeneity
      }
      else
        overall.hetstat <- replaceNULL(overall, FALSE) &&
          prepared$shared.overall &&
          prepared$shared.heterogeneity
    }
  }
  else if (with.subgroups) {
    warning("Nothing selected to print for subgroup meta-analyses.",
            call. = FALSE)
    return(invisible(NULL))
  }
  #
  if (!with.subgroups) {
    overall <- FALSE
    overall.hetstat <- FALSE
  }
  #
  overall <- replaceNULL(overall, FALSE)
  overall.hetstat <- replaceNULL(overall.hetstat, FALSE)
  #
  chklogical(overall)
  chklogical(overall.hetstat)
  #
  chklogical(backtransf)
  #
  chkchar(big.mark, length = 1)
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  #
  chklogical(details.methods)
  #
  # Columns shown when printing
  #
  print.list <-
    merge_legacy(print, args, print_defaults(x, "metabind", func = "print"))
  #
  print <- print.list$res
  args <- print.list$args
  #
  print$pval.Q <- replaceNULL(print$pval.Q, print$Q) && print$Q
  print$pval.Q.b <- replaceNULL(print$pval.Q.b, with.subgroups)
  print$tau2 <- replaceNULL(print$tau2, !with.subgroups && random)
  print$tau <- replaceNULL(print$tau, !with.subgroups && random)
  print$subgroup.name <- replaceNULL(print$subgroup.name, with.subgroups)
  #
  # Check every option, retaining the partial matching used by '$'.
  #
  for (i in c(
    "func",
    "backtransf",
    "args.backtransf",
    "method",
    "method.interval",
    "k",
    "sm",
    "se",
    "pval",
    "stat",
    "level.ma",
    "level.predict",
    "level.hetstat",
    "null.effect",
    "I2",
    "I2.ci",
    "tau2",
    "tau2.ci",
    "tau",
    "tau.ci",
    "Q",
    "pval.Q",
    "H",
    "H.ci",
    "pval.Q.b",
    "subgroup.name"
  ))
  chklogical(
    print[[i, exact = FALSE]],
    name = paste0("print.", i),
    NULL.ok =
      i %in% c(
        "func",
        "backtransf",
        "args.backtransf",
        "method",
        "method.interval",
        "sm",
        "level.ma",
        "level.predict",
        "level.hetstat"
      )
  )
  #
  print$tau2.ci <- print$tau2 && print$tau2.ci
  print$tau.ci <- print$tau && print$tau.ci
  print$I2.ci <- print$I2 && print$I2.ci
  print$H.ci <- print$H && print$H.ci
  #
  print$het <- print$Q ||
    print$tau2 || print$tau || print$I2 || print$H
  #
  if (!common && !random && !prediction && !print$het) {
    warning("Nothing selected to print.", call. = FALSE)
    return(invisible(NULL))
  }
  #
  # Digits used to round results
  #
  digits <- as_named_list(digits)
  digits.list <-
    merge_legacy(digits, args, digits_defaults("metabind", func = "print"))
  #
  digits <- digits.list$res
  args <- digits.list$args
  #
  for (i in c("digits",
              "se",
              "stat",
              "pval",
              "pval.Q",
              "Q",
              "tau2",
              "tau",
              "I2"))
    chknumeric(
      digits[[i, exact = FALSE]],
      min = if (i %in% c("pval", "pval.Q"))
        1
      else
        0,
      length = 1,
      name = if (i == "digits")
        "digits"
      else
        paste0("digits.", i)
    )
  
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
    warning(
      "Argument 'backtransf' set to TRUE as ",
      "meta-analyses have different summary measures.",
      call. = FALSE
    )
    backtransf <- TRUE
  }
  #
  # Resolve forward-transformation settings for selected components
  #
  ids <- unique(c(x$common$id, x$random$id, x$predict$id))
  meta <- x$overall
  #
  if (length(ids) > 0)
    meta <- meta[meta$id %in% ids, , drop = FALSE]
  #
  transf.settings <- lapply(seq_len(nrow(meta)), function(i) {
    func <- if (is.null(meta[["func.transf"]]))
      NULL
    else
      meta[["func.transf"]][[i]]
    #
    args <- if (is.null(meta[["args.transf"]]))
      NULL
    else
      meta[["args.transf"]][[i]]
    #
    if (is.atomic(func) && length(func) == 1L &&
        (is.na(func) || identical(func, "")))
      func <- NULL
    #
    if (is.null(func) || length(args) == 0)
      args <- NULL
    #
    list(func = func, args = args)
  })
  #
  same.transf <- all(vapply(transf.settings, function(setting)
    identical(setting, transf.settings[[1]]), logical(1)))
  #
  ftr <- atr <- NULL
  #
  if (same.transf) {
    ftr <- transf.settings[[1]]$func
    atr <- transf.settings[[1]]$args
  }
  #
  if (same.sm) {
    sm <- sm[1]
    #
    fbt <- makeunique(x$overall$func.backtransf, "")
    is_relative <- is_relative_effect(sm) |
      (!is.null(fbt) && fbt == "exp")
    #
    sm.lab <- smlab(
      sm,
      backtransf = backtransf,
      pscale = pscale,
      irscale = irscale,
      func.backtransf = fbt,
      fallback = ".no_sm_provided"
    )
    #
    if (!same.transf)
      lab.TE <- "TE"
    else if (!is.null(ftr)) {
      func.lab <- if (is.character(ftr))
        ftr
      else
        deparse1(ftr)
      #
      lab.TE <- paste0(func.lab,
                       "(",
                       if (sm == "")
                         "TE"
                       else
                         sm,
                       if (length(atr) > 0)
                         paste0(", ", format_args(atr)),
                       ")")
    }
    else if (is_relative || sm == "VE")
      lab.TE <- paste0(gs("log.prefix"), if (sm == "VE")
        "VR"
        else
          sm)
    else if (sm == "")
      lab.TE <- "TE"
    else
      lab.TE <- sm
    #
    if (!backtransf)
      sm.lab <- lab.TE
  }
  else {
    if (x$sm.family == "relative") {
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
  
  if (common || random || prediction) {
    if (with.subgroups) {
      overall.hetstat <-
        prepared$shared.overall && prepared$shared.heterogeneity &&
        replaceNULL(overall.hetstat, TRUE)
      #
      pooled <- prepared$dat
      subdat <- prepared$subdat
      #
      # Parent heterogeneity must not appear as subgroup heterogeneity
      #
      for (i in c(
        "Q",
        "df.Q",
        "pval.Q",
        "tau2",
        "lower.tau2",
        "upper.tau2",
        "tau",
        "lower.tau",
        "upper.tau",
        "I2",
        "lower.I2",
        "upper.I2",
        "H",
        "lower.H",
        "upper.H"
      ))
      subdat[[i]] <- NA_real_
      #
      dat <- subdat
    }
    else {
      dat$subgroup.levels <- NA
      dat$subgroup.name <- NA
      #
      dat$Q.b <- NA
      dat$df.Q.b <- NA
      dat$pval.Q.b <- NA
    }

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
      dat <- backtransf_metabind(dat, pscale = pscale, irscale = irscale)
      #
      if (with.subgroups)
        pooled <- backtransf_metabind(pooled,
                                      pscale = pscale,
                                      irscale = irscale,
                                      suffix = ".pooled")
    }
    #
    # Results for common effect model
    #
    if (common) {
      dat.c <-
        format_metabind_results(
          dat,
          "common",
          print,
          digits,
          sm.lab,
          lab.TE,
          with.subgroups = with.subgroups,
          lab.NA = lab.NA,
          big.mark = big.mark,
          scientific.pval = scientific.pval,
          zero.pval = zero.pval,
          JAMA.pval = JAMA.pval
        )
    }
    #
    # Results for random effects model
    #
    if (random) {
      dat.r <-
        format_metabind_results(
          dat,
          "random",
          print,
          digits,
          sm.lab,
          lab.TE,
          with.subgroups = with.subgroups,
          lab.NA = lab.NA,
          big.mark = big.mark,
          scientific.pval = scientific.pval,
          zero.pval = zero.pval,
          JAMA.pval = JAMA.pval
        )
    }
    #
    # Prediction intervals
    #
    if (prediction) {
      dat.p <-
        format_metabind_results(
          dat,
          "prediction",
          print,
          digits,
          sm.lab,
          lab.TE,
          with.subgroups = with.subgroups,
          lab.NA = lab.NA,
          big.mark = big.mark,
          scientific.pval = scientific.pval,
          zero.pval = zero.pval,
          JAMA.pval = JAMA.pval
        )
    }
  }

  #
  #
  # (4) Create data set with heterogeneity statistics
  #
  #

  if (print$het) {
    if (common || random || prediction) {
      ids.het <-
        unique(c(if (common)
          x$common$id, if (random)
            x$random$id, if (prediction)
              x$predict$id))
    }
    else
      ids.het <- x$overall$id
    #
    dat.h <- format_metabind_heterogeneity(
      x$overall,
      print,
      digits,
      ids = ids.het,
      lab.NA = lab.NA,
      big.mark = big.mark,
      scientific.pval = scientific.pval,
      zero.pval = zero.pval,
      JAMA.pval = JAMA.pval,
      text.tau2 = text.tau2,
      text.tau = text.tau,
      text.I2 = text.I2
    )
  }

  #
  #
  # (5) Print results
  #
  #

  if (!with.subgroups) {
    cat("Meta-analyses results\n")
    #
    if (common) {
      cat("\nCommon effect model\n\n")
      #
      dat.c$id <- NULL
      #
      name.c <- dat.c$name
      dat.c$name <- NULL
      #
      prmatrix(dat.c,
               quote = FALSE,
               right = TRUE,
               rowlab = name.c)
    }
    #
    if (random) {
      cat("\nRandom effects model\n\n")
      #
      dat.r$id <- NULL
      #
      name.r <- dat.r$name
      dat.r$name <- NULL
      #
      prmatrix(dat.r,
               quote = FALSE,
               right = TRUE,
               rowlab = name.r)
    }
    #
    if (prediction) {
      cat("\nPrediction interval\n\n")
      #
      dat.p$id <- NULL
      #
      name.p <- dat.p$name
      dat.p$name <- NULL
      #
      prmatrix(dat.p,
               quote = FALSE,
               right = TRUE,
               rowlab = name.p)
    }
    #
    if (print$het) {
      cat("\nInformation on heterogeneity\n\n")
      #
      name.h <- dat.h$name
      dat.h$name <- NULL
      #
      prmatrix(dat.h,
               quote = FALSE,
               right = TRUE,
               rowlab = name.h)
    }
  }
  else {
    if (overall) {
      cat("Overall results\n\n")
      pooled  %<>%
        rename_with( ~ sub("\\.pooled$", "", .x), ends_with(".pooled"))
      #
      levels.pooled <-
        unique(round(c(
          pooled$level.ma[pooled$type %in%
                            c(if (common)
                              "Common effect model", if (random)
                                "Random effects model")], if (prediction)
                                  pooled$level.predict[pooled$type ==
                                                         "Prediction interval"]
        ), 6))
      #
      same.level.pooled <- length(levels.pooled) == 1 &&
        !anyNA(levels.pooled)
      #
      if (prediction && (common || random))
        type.interval <- "CI/PI"
      else if (prediction)
        type.interval <-  "PI"
      else
        type.interval <- "CI"
      #
      lab.interval.pooled <- if (same.level.pooled)
        percent_label(levels.pooled[1], type = type.interval)
      else
        type.interval
      #
      pooled.c <- pooled.r <- pooled.p <- NULL
      #
      if (common) {
        pooled.c <-
          format_metabind_results(
            pooled,
            "common",
            print,
            digits,
            sm.lab,
            lab.TE,
            with.subgroups = FALSE,
            lab.NA = lab.NA,
            big.mark = big.mark,
            scientific.pval = scientific.pval,
            zero.pval = zero.pval,
            JAMA.pval = JAMA.pval,
            keep = "text",
            lab.interval = lab.interval.pooled,
            show.level = replaceNULL(print[["level.ma"]], !same.level.pooled)
          )
      }
      #
      if (random) {
        pooled.r <-
          format_metabind_results(
            pooled,
            "random",
            print,
            digits,
            sm.lab,
            lab.TE,
            with.subgroups = FALSE,
            lab.NA = lab.NA,
            big.mark = big.mark,
            scientific.pval = scientific.pval,
            zero.pval = zero.pval,
            JAMA.pval = JAMA.pval,
            keep = "text",
            lab.interval = lab.interval.pooled,
            show.level = replaceNULL(print[["level.ma"]], !same.level.pooled)
          )
      }
      #
      if (prediction) {
        pooled.p <-
          format_metabind_results(
            pooled,
            "pred",
            print,
            digits,
            sm.lab,
            lab.TE,
            with.subgroups = FALSE,
            lab.NA = lab.NA,
            big.mark = big.mark,
            scientific.pval = scientific.pval,
            zero.pval = zero.pval,
            JAMA.pval = JAMA.pval,
            keep = "text",
            lab.interval = lab.interval.pooled,
            show.level =
              replaceNULL(print[["level.predict"]], !same.level.pooled)
          )
      }
      #
      pooled.all <- bind_rows(pooled.c, pooled.r, pooled.p)
      #
      names.pooled <- pooled.all$text
      pooled.all$text <- NULL
      #
      pooled.all[is.na(pooled.all)] <- lab.NA
      #
      prmatrix(pooled.all,
               quote = FALSE,
               right = TRUE,
               rowlab = names.pooled)
      #
      if (print$het && overall.hetstat) {
        het.pooled <- format_metabind_heterogeneity(
          x$overall,
          print,
          digits,
          pooled = TRUE,
          ids = unique(c(
            if (common)
              x$common$id, if (random)
                x$random$id, if (prediction)
                  x$predict$id
          )),
          lab.NA = lab.NA,
          big.mark = big.mark,
          scientific.pval = scientific.pval,
          zero.pval = zero.pval,
          JAMA.pval = JAMA.pval,
          text.tau2 = text.tau2,
          text.tau = text.tau,
          text.I2 = text.I2
        )
        #
        cat("\n")
        prmatrix(
          het.pooled,
          quote = FALSE,
          right = TRUE,
          rowlab = rep_len("", nrow(het.pooled))
        )
      }
      #
      cat("\n")
    }
    #
    cat("Subgroup meta-analysis results\n\n")
    #
    if (common) {
      cat("Common effect model\n")
      #
      j <- 0
      #
      for (i in unique(dat.c$id)) {
        j <- j + 1
        #
        dat.i <- dat.c %>% filter(id == i) %>%
          select(-id, -name) %>%
          rename(subgroup = subgroup.levels)
        #
        if (print[["subgroup.name"]])
          cat(paste0("\n", j, ". ", unique(dat.i$subgroup.name), "\n\n"))
        else
          cat(paste0("\n", j, ".\n\n"))
        #
        dat.i$subgroup.name <- NULL
        #
        if (isCol(dat.i, "level") &&
            length(unique(dat.i$level)) == 1 &&
            !anyNA(dat.i$level)) {
          names(dat.i)[names(dat.i) == "CI"] <-
            percent_label(dat.i$level[1])
        }
        #
        if (isCol(dat.i, "level"))
          dat.i$level <- formatN(
            dat.i$level,
            digits = 4,
            lab.NA = lab.NA,
            big.mark = big.mark
          )
        #
        prmatrix(
          dat.i,
          quote = FALSE,
          right = TRUE,
          rowlab = rep_len("", nrow(dat.i))
        )
      }
    }
    #
    if (random) {
      cat("\nRandom effects model\n")
      #
      j <- 0
      #
      for (i in unique(dat.r$id)) {
        j <- j + 1
        #
        dat.i <- dat.r %>% filter(id == i) %>%
          select(-id, -name) %>%
          rename(subgroup = subgroup.levels)
        #
        if (print[["subgroup.name"]])
          cat(paste0("\n", j, ". ", unique(dat.i$subgroup.name), "\n\n"))
        else
          cat(paste0("\n", j, ".\n\n"))
        #
        dat.i$subgroup.name <- NULL
        #
        if (isCol(dat.i, "level") &&
            length(unique(dat.i$level)) == 1 &&
            !anyNA(dat.i$level)) {
          names(dat.i)[names(dat.i) == "CI"] <-
            percent_label(dat.i$level[1])
        }
        #
        if (isCol(dat.i, "level"))
          dat.i$level <- formatN(
            dat.i$level,
            digits = 4,
            lab.NA = lab.NA,
            big.mark = big.mark
          )
        #
        prmatrix(
          dat.i,
          quote = FALSE,
          right = TRUE,
          rowlab = rep_len("", nrow(dat.i))
        )
      }
    }
    #
    if (prediction) {
      cat("\nPrediction interval\n")
      #
      j <- 0
      #
      for (i in unique(dat.p$id)) {
        j <- j + 1
        #
        dat.i <- dat.p %>% filter(id == i) %>%
          select(-id, -name) %>%
          rename(subgroup = subgroup.levels)
        #
        if (print[["subgroup.name"]])
          cat(paste0("\n", j, ". ", unique(dat.i$subgroup.name), "\n\n"))
        else
          cat(paste0("\n", j, ".\n\n"))
        #
        dat.i$subgroup.name <- NULL
        #
        if (isCol(dat.i, "level"))
          dat.i$level <- formatN(
            dat.i$level,
            digits = 4,
            lab.NA = lab.NA,
            big.mark = big.mark
          )
        #
        prmatrix(
          dat.i,
          quote = FALSE,
          right = TRUE,
          rowlab = rep_len("", nrow(dat.i))
        )
      }
    }
  }
  #
  if (details.methods) {
    catmeth_metabind(
      x,
      #
      common,
      random,
      prediction,
      #
      backtransf,
      print,
      digits,
      #
      pscale,
      irscale,
      irunit,
      #
      big.mark = big.mark,
      #
      forest = FALSE
    )
  }
  #
  invisible(NULL)
}
