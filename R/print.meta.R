print.meta <- function(x,
                       sortvar,
                       comb.fixed = x$comb.fixed,
                       comb.random = x$comb.random,
                       prediction = x$prediction,
                       details = FALSE, ma = TRUE,
                       backtransf = x$backtransf,
                       pscale = x$pscale,
                       irscale = x$irscale,
                       irunit = x$irunit,
                       digits = gs("digits"),
                       digits.se = gs("digits.se"),
                       digits.tau2 = gs("digits.tau2"),
                       digits.I2 = gs("digits.I2"),
                       digits.prop = gs("digits.prop"),
                       digits.weight = gs("digits.weight"),
                       big.mark = gs("big.mark"),
                       warn.backtransf = FALSE,
                       ...
                       ) {
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta object
  ##
  ##
  chkclass(x, "meta")
  x <- updateversion(x)
  ##
  k.all <- length(x$TE)


  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  mf <- match.call()
  error <- try(sortvar <- eval(mf[[match("sortvar", names(mf))]],
                               as.data.frame(x, stringsAsFactors = FALSE),
                               enclos = sys.frame(sys.parent())),
               silent = TRUE)
  if (class(error) == "try-error") {
    xd <- x$data
    sortvar <- eval(mf[[match("sortvar", names(mf))]],
                    xd, enclos = NULL)
    if (!is.null(x$data$.subset))
      sortvar <- sortvar[x$data$.subset]
  }
  sort <- !is.null(sortvar)
  if (sort && (length(sortvar) != k.all))
    stop("Number of studies in object 'x' and argument 'sortvar' have different length.")
  if (!sort)
    sortvar <- 1:k.all
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(prediction)
  chklogical(details)
  chklogical(ma)
  ##
  if (is.untransformed(x$sm))
    backtransf <- TRUE
  chklogical(backtransf)
  ##
  chklogical(warn.backtransf)
  ##
  if (!is.prop(x$sm) & x$sm != "RD")
    pscale <- 1
  if (!is.null(pscale))
    chknumeric(pscale, single = TRUE)
  else
    pscale <- 1
  if (!backtransf & pscale != 1 & !is.untransformed(x$sm)) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  if (!is.rate(x$sm) & x$sm != "IRD")
    irscale <- 1
  if (!is.null(irscale))
    chknumeric(irscale, single = TRUE)
  else
    irscale <- 1
  if (!backtransf & irscale != 1 & !is.untransformed(x$sm)) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  if (!is.null(irunit) && !is.na(irunit))
    chkchar(irunit)
  ##
  chknumeric(digits, min = 0, single = TRUE)
  chknumeric(digits.se, min = 0, single = TRUE)
  chknumeric(digits.tau2, min = 0, single = TRUE)
  chknumeric(digits.I2, min = 0, single = TRUE)
  chknumeric(digits.prop, min = 0, single = TRUE)
  ##
  ## Additional arguments / checks for metacont objects
  ##
  cl <- paste("update.meta() or ", class(x)[1], "()", sep = "")
  addargs <- names(list(...))
  ##
  fun <- "print.meta"
  ##
  warnarg("level", addargs, fun, cl)
  warnarg("level.comb", addargs, fun, cl)
  warnarg("level.predict", addargs, fun, cl)
  warnarg("logscale", addargs, fun, otherarg = "backtransf")
  ##
  level <- x$level
  level.comb <- x$level.comb
  level.predict <- x$level.predict
  
  
  ##
  ##
  ## (3) Some additional settings
  ##
  ##
  metainf.metacum <- inherits(x, "metainf") | inherits(x, "metacum")
  mb <- inherits(x, "metabind")
  ##  
  prediction <- prediction & x$k >= 3
  if (is.na(prediction))
    prediction <- FALSE
  ##  
  ci.lab <- paste(round(100 * level, 1), "%-CI", sep = "")
  ##  
  sm <- x$sm
  ##  
  sm.lab <- sm
  ##
  if (backtransf) {
    if (sm == "ZCOR")
      sm.lab <- "COR"
    else if (is.mean(sm))
      sm.lab <- "mean"
    else if (is.prop(sm)) {
      if (pscale == 1)
        sm.lab <- "proportion"
      else
        sm.lab <- "events"
    }
    else if (is.rate(sm)) {
      if (irscale == 1)
        sm.lab <- "rate"
      else
        sm.lab <- "events"
    }
  }
  else
    if (is.relative.effect(sm))
      sm.lab <- paste("log", sm, sep = "")
  
  
  ##
  ##
  ## (4) Print title and details
  ##
  ##
  crtitle(x)
  ##  
  if (details) {
    if (inherits(x, "metabin")) {      
      res <- data.frame(event.e = formatN(x$event.e, digits = 0,
                                          "NA", big.mark = big.mark),
                        n.e = formatN(x$n.e, digits = 0,
                                      "NA", big.mark = big.mark),
                        event.c = formatN(x$event.c, digits = 0,
                                          "NA", big.mark = big.mark),
                        n.c = formatN(x$n.c, digits = 0,
                                      "NA", big.mark = big.mark))
      ##
      if (pscale == 1) {
        res$p.e <- formatN(round(x$event.e / x$n.e, digits.prop),
                           big.mark = big.mark)
        res$p.c <- formatN(round(x$event.c / x$n.c, digits.prop),
                           big.mark = big.mark)
      }
      else {
        res$events.e <- formatN(round(pscale * x$event.e / x$n.e, digits),
                                digits,
                                "NA", big.mark = big.mark)
        res$events.c <- formatN(round(pscale * x$event.c / x$n.c, digits),
                                digits,
                                "NA", big.mark = big.mark)
      }
    }
    else if (inherits(x, "metacont")) {
      res <- data.frame(n.e = formatN(x$n.e, digits = 0,
                                      "NA", big.mark = big.mark),
                        mean.e = formatN(round(x$mean.e, digits), digits,
                                         "NA", big.mark = big.mark),
                        sd.e = formatN(round(x$sd.e, digits.se), digits.se,
                                       "NA", big.mark = big.mark),
                        n.c = formatN(x$n.c, digits = 0,
                                      "NA", big.mark = big.mark),
                        mean.c = formatN(round(x$mean.c, digits), digits,
                                         "NA", big.mark = big.mark),
                        sd.c = formatN(round(x$sd.c, digits.se), digits.se,
                                       "NA", big.mark = big.mark))
    }
    else if (inherits(x, "metacor")) {
      res <- data.frame(cor = x$cor,
                        n = formatN(x$n, digits = 0,
                                    "NA", big.mark = big.mark))
    }
    else if (inherits(x, "metagen")) {
      res <- data.frame(TE = formatN(round(x$TE, digits), digits,
                                     "NA", big.mark = big.mark),
                        seTE = formatN(round(x$seTE, digits.se), digits.se,
                                       "NA", big.mark = big.mark))
    }
    else if (inherits(x, "metainc")) {
      res <- data.frame(event.e = formatN(x$event.e, digits = 0,
                                          "NA", big.mark = big.mark),
                        time.e = formatN(round(x$time.e, digits), digits,
                                         "NA", big.mark = big.mark),
                        event.c = formatN(x$event.c, digits = 0,
                                          "NA", big.mark = big.mark),
                        time.c = formatN(round(x$time.c, digits), digits,
                                         "NA", big.mark = big.mark))
      ##
      if (irscale == 1) {
        res$rate.e <- formatN(round(x$event.e / x$time.e, digits.prop),
                              big.mark = big.mark)
        res$rate.c <- formatN(round(x$event.c / x$time.c, digits.prop),
                              big.mark = big.mark)
      }
      else {
        res$events.e <- formatN(round(irscale * x$event.e / x$n.e, digits),
                                digits,
                                "NA", big.mark = big.mark)
        res$events.c <- formatN(round(irscale * x$event.c / x$n.c, digits),
                                digits,
                                "NA", big.mark = big.mark)
      }
    }
    else if (inherits(x, "metaprop")) {
      res <- data.frame(event = formatN(x$event, digits = 0,
                                        "NA", big.mark = big.mark),
                        n = formatN(x$n, digits = 0,
                                    "NA", big.mark = big.mark))
      if (pscale == 1)
        res$p <- formatN(round(x$event / x$n, digits.prop), digits.prop,
                         "NA", big.mark = big.mark)
      else
        res$events <- formatN(round(pscale * x$event / x$n, digits), digits,
                              "NA", big.mark = big.mark)
    }
    else if (inherits(x, "metarate")) {
      res <- data.frame(event = formatN(x$event, digits = 0,
                                        "NA", big.mark = big.mark),
                        time = formatN(x$time, digits = digits,
                                       "NA", big.mark = big.mark))
      if (irscale == 1)
        res$rate <- formatN(round(x$event / x$time, digits.prop), digits.prop,
                            "NA", big.mark = big.mark)
      else
        res$events <- formatN(round(irscale * x$event / x$time, digits),
                              digits, "NA", big.mark = big.mark)
    }
    else {
      res <- data.frame(TE = formatN(round(x$TE, digits), digits,
                                     "NA", big.mark = big.mark),
                        seTE = formatN(round(x$seTE, digits), digits,
                                       "NA", big.mark = big.mark))
    }
    dimnames(res)[[1]] <- x$studlab
    prmatrix(res[order(sortvar), ], quote = FALSE, right = TRUE)
    cat("\n\n")
  }
  
  
  ##
  ##
  ## (5) Print results for individual studies
  ##
  ##
  if (k.all == 1 & !inherits(x, "metaprop")) {
    if (inherits(x, "metabin") & x$method == "MH")
      print(summary(metabin(x$event.e, x$n.e,
                            x$event.c, x$n.c,
                            sm = sm,
                            method = "Inverse",
                            studlab = x$studlab,
                            incr = x$incr,
                            allincr = x$allincr,
                            doublezeros = x$doublezeros,
                            MH.exact = x$MH.exact,
                            warn = FALSE, level.comb = level.comb)),
            header = FALSE,
            digits = digits,
            backtransf = backtransf, pscale = pscale,
            irscale = irscale, irunit = irunit, big.mark = big.mark,
            warn.backtransf = warn.backtransf,
            ...)
    else
      print(summary(x),
            header = FALSE,
            digits = digits,
            backtransf = backtransf, pscale = pscale,
            irscale = irscale, irunit = irunit, big.mark = big.mark,
            warn.backtransf = warn.backtransf,
            ...)
  }
  else {
    TE <- x$TE
    seTE <- x$seTE
    lowTE <- x$lower
    uppTE <- x$upper
    ##
    if (inherits(x, "metaprop") & !backtransf) {
      ciTE <- ci(TE, seTE, level = level)
      lowTE <- ciTE$lower
      uppTE <- ciTE$upper
      ##
      x$method.ci <- "NAsm"
    }
    ##
    if (backtransf) {
      ## Freeman-Tukey Arcsin transformation
      if (metainf.metacum) {
        if (sm == "IRFT")
          harmonic.mean <- x$t.harmonic.mean
        else
          harmonic.mean <- x$n.harmonic.mean
      }
      else {
        if (sm == "IRFT")
          harmonic.mean <- x$time
        else
          harmonic.mean <- x$n
      }
      ##
      if (inherits(x, "metaprop"))
        TE <- x$event / x$n
      else {
        TE    <- backtransf(   TE, sm, "mean",  harmonic.mean, warn.backtransf)
        lowTE <- backtransf(lowTE, sm, "lower", harmonic.mean, warn.backtransf)
        uppTE <- backtransf(uppTE, sm, "upper", harmonic.mean, warn.backtransf)
      }
      ##
      if (is.prop(sm) | sm == "RD") {
        TE <- pscale * TE
        lowTE <- pscale * lowTE
        uppTE <- pscale * uppTE
      }
      ##
      if (is.rate(sm) | sm == "IRD") {
        TE <- irscale * TE
        lowTE <- irscale * lowTE
        uppTE <- irscale * uppTE
      }
    }
    ##
    TE <- round(TE, digits)
    lowTE <- round(lowTE, digits)
    uppTE <- round(uppTE, digits)
    ##
    if (!metainf.metacum) {
      if (comb.fixed)
        if (!all(is.na(x$w.fixed)) && sum(x$w.fixed) > 0)
          w.fixed.p <- round(100 * x$w.fixed / sum(x$w.fixed, na.rm = TRUE),
                             digits.weight)
        else w.fixed.p <- x$w.fixed
      ##
      if (comb.random)
        if (!is.null(x$w.random) & !all(is.na(x$w.random)) &&
            sum(x$w.random) > 0)
          w.random.p <- round(100 * x$w.random / sum(x$w.random, na.rm = TRUE),
                              digits.weight)
        else w.random.p <- x$w.random
    }
    ##    
    if (metainf.metacum) {
      is.random <- x$pooled == "random"
      ##
      I2 <- formatN(round(100 * x$I2, digits.I2), digits.I2, "")
      ##
      sel <- is.na(x$p.value)
      p.value <- formatPT(x$p.value)
      p.value <- ifelse(sel, "", p.value)
      ##
      tau2 <- x$tau^2
      tau2 <- formatN(round(tau2, digits.tau2), digits.tau2, "",
                      big.mark = big.mark)
      ##
      res <- cbind(formatN(round(TE, digits), digits, "",
                           big.mark = big.mark),
                   formatCI(formatN(round(lowTE, digits), digits, "NA",
                                    big.mark = big.mark),
                            formatN(round(uppTE, digits), digits, "NA",
                                    big.mark = big.mark)),
                   p.value,
                   paste(" ", tau2, sep = ""),
                   paste(" ", I2, ifelse(I2 == "", "", "%"), sep = ""))
      dimnames(res) <- list(paste(x$studlab, "  ", sep = ""),
                            c(sm.lab, ci.lab, "p-value",
                              gs("text.tau2"), gs("text.I2")))
      ##
      if (inherits(x, "metainf")) {
        if (!is.random)
          cat("\nInfluential analysis (Fixed effect model)\n")
        else
          cat("\nInfluential analysis (Random effects model)\n")
      }
      if (inherits(x, "metacum")) {
        if (!is.random)
          cat("\nCumulative meta-analysis (Fixed effect model)\n")
        else
          cat("\nCumulative meta-analysis (Random effects model)\n")
      }
      cat("\n")
      prmatrix(res, quote = FALSE, right = TRUE, na.print = "--")
      ## Print information on summary method:
      catmeth(class = class(x),
              method = x$method,
              method.tau = if (is.random) x$method.tau else "",
              sm = sm,
              k.all = k.all,
              hakn = is.random & x$hakn,
              tau.preset = x$tau.preset,
              method.smd = x$method.smd,
              sd.glass = x$sd.glass,
              exact.smd = x$exact.smd,
              model.glmm = x$model.glmm,
              big.mark = big.mark,
              digits = digits, digits.tau2 = digits.tau2,
              text.tau2 = gs("text.tau2"))
    }
    else {
      res <- cbind(formatN(round(TE, digits), digits, "NA", big.mark = big.mark),
                   formatCI(formatN(round(lowTE, digits), digits, "NA",
                                    big.mark = big.mark),
                            formatN(round(uppTE, digits), digits, "NA",
                                    big.mark = big.mark)),
                   if (comb.fixed & !mb)
                     formatN(w.fixed.p, digits.weight,
                             big.mark = big.mark),
                   if (comb.random & !mb)
                     formatN(w.random.p, digits.weight,
                             big.mark = big.mark),
                   if (!is.null(x$byvar)) x$byvar,
                   if (!is.null(x$exclude))
                     ifelse(is.na(x$exclude), "",
                     ifelse(x$exclude, "*", "")))
      ## Printout for a single proportion:
      if (k.all == 1) {
        ##
        if (!is.null(x$method.ci)) {
          if  (x$method.ci == "CP")
            method.ci.details <- "Clopper-Pearson confidence interval:\n\n"
          else if (x$method.ci == "WS")
            method.ci.details <- "Wilson Score confidence interval:\n\n"
          else if (x$method.ci == "WSCC")
            method.ci.details <- "Wilson Score confidence interval with continuity correction:\n\n"
          else if (x$method.ci == "AC")
            method.ci.details <- "Agresti-Coull confidence interval:\n\n"
          else if (x$method.ci == "SA")
            method.ci.details <- "Simple approximation confidence interval:\n\n"
          else if (x$method.ci == "SACC")
            method.ci.details <- "Simple approximation confidence interval with continuity correction:\n\n"
          if (x$method.ci != "NAsm") {
            cat(method.ci.details)
            dimnames(res) <- list("", c(sm.lab, ci.lab,
                                        if (comb.fixed & !mb) "%W(fixed)",
                                        if (comb.random & !mb) "%W(random)",
                                        if (!is.null(x$byvar)) x$bylab,
                                        if (!is.null(x$exclude)) "exclude"))
            prmatrix(res, quote = FALSE, right = TRUE)
            cat("\n\n")
          }
        }
        cat("Normal approximation confidence interval:\n")
      }
      else {
        dimnames(res) <-
          list(x$studlab, c(sm.lab, ci.lab,
                            if (comb.fixed & !mb) "%W(fixed)",
                            if (comb.random & !mb) "%W(random)",
                            if (!is.null(x$byvar)) x$bylab,
                            if (!is.null(x$exclude)) "exclude"))
        prmatrix(res[order(sortvar),], quote = FALSE, right = TRUE)
      }
    }
    
    
    ##
    ##
    ## (6) Print result for meta-analysis
    ##
    ##
    if (ma & !metainf.metacum) {
      if (!is.na(x$k))
        cat("\n")
      ##
      print(summary(x, warn = FALSE),
            header = FALSE,
            digits = digits,
            comb.fixed = comb.fixed, comb.random = comb.random,
            prediction = prediction,
            backtransf = backtransf, pscale = pscale,
            irscale = irscale, irunit = irunit,
            digits.tau2 = digits.tau2, digits.I2 = digits.I2,
            big.mark = big.mark,
            warn.backtransf = warn.backtransf,
            ...)
    }
  }
  
  
  invisible(NULL)
}
