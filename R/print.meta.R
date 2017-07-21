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
                       digits.zval = gs("digits.zval"),
                       digits.pval = max(gs("digits.pval"), 2),
                       digits.pval.Q = max(gs("digits.pval.Q"), 2),
                       digits.Q = gs("digits.Q"),
                       digits.tau2 = gs("digits.tau2"),
                       digits.H = gs("digits.H"),
                       digits.I2 = gs("digits.I2"),
                       digits.prop = gs("digits.prop"),
                       digits.weight = gs("digits.weight"),
                       scientific.pval = gs("scientific.pval"),
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
  chklogical(backtransf)
  if (!(x$sm %in% c("PLOGIT", "PLN", "PRAW", "PAS", "PFT")))
    pscale <- 1
  if (!is.null(pscale))
    chknumeric(pscale, single = TRUE)
  else
    pscale <- 1
  if (!backtransf & pscale != 1) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  if (!(x$sm %in% c("IR", "IRLN", "IRS", "IRFT")))
    irscale <- 1
  if (!is.null(irscale))
    chknumeric(irscale, single = TRUE)
  else
    irscale <- 1
  if (!backtransf & irscale != 1) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  chknumeric(digits, min = 0, single = TRUE)
  chknumeric(digits.se, min = 0, single = TRUE)
  chknumeric(digits.tau2, min = 0, single = TRUE)
  chknumeric(digits.zval, min = 0, single = TRUE)
  chknumeric(digits.pval, min = 1, single = TRUE)
  chknumeric(digits.pval.Q, min = 1, single = TRUE)
  chknumeric(digits.Q, min = 0, single = TRUE)
  chknumeric(digits.H, min = 0, single = TRUE)
  chknumeric(digits.I2, min = 0, single = TRUE)
  chknumeric(digits.prop, min = 0, single = TRUE)
  chklogical(scientific.pval)
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
  ##  
  prediction <- prediction & x$k >= 3
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
    else if (sm %in% c("PLOGIT", "PLN", "PRAW", "PAS", "PFT")) {
      if (pscale == 1)
        sm.lab <- "proportion"
      else
        sm.lab <- "events"
    }
    else if (sm %in% c("IR", "IRLN", "IRS", "IRFT")) {
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
      res <- data.frame(event.e = x$event.e, n.e = x$n.e,
                        event.c = x$event.c, n.c = x$n.c,
                        p.e = format.NA(round(x$event.e / x$n.e, digits.prop)),
                        p.c = format.NA(round(x$event.c / x$n.c, digits.prop)))
    }
    else if (inherits(x, "metacont")) {
      res <- data.frame(n.e = x$n.e,
                        mean.e = format.NA(round(x$mean.e, digits), digits, "NA"),
                        sd.e = format.NA(round(x$sd.e, digits.se), digits.se, "NA"),
                        n.c = x$n.c,
                        mean.c = format.NA(round(x$mean.c, digits), digits, "NA"),
                        sd.c = format.NA(round(x$sd.c, digits.se), digits.se, "NA"))
    }
    else if (inherits(x, "metacor")) {
      res <- data.frame(cor = x$cor, n = x$n)
    }
    else if (inherits(x, "metagen")) {
      res <- data.frame(TE = format.NA(round(x$TE, digits), digits, "NA"),
                        seTE = format.NA(round(x$seTE, digits.se), digits.se, "NA"))
    }
    else if (inherits(x, "metainc")) {
      res <- data.frame(event.e = x$event.e,
                        time.e = format.NA(round(x$time.e, digits), digits, "NA"),
                        event.c = x$event.c,
                        time.c = format.NA(round(x$time.c, digits), digits, "NA"))
    }
    else if (inherits(x, "metaprop")) {
      res <- data.frame(event = x$event, n = x$n)
      if (pscale == 1)
        res$p <- format.NA(round(x$event / x$n, digits.prop), digits.prop, "NA")
      else
        res$events <- format.NA(round(pscale * x$event / x$n, digits.prop), digits.prop, "NA")
    }
    else if (inherits(x, "metarate")) {
      res <- data.frame(event = x$event, time = x$time)
      if (irscale == 1)
        res$rate <- format.NA(round(x$event / x$time, digits.prop), digits.prop, "NA")
      else
        res$events <- format.NA(round(irscale * x$event / x$time, digits.prop), digits.prop, "NA")
    }
    else {
      res <- data.frame(TE = format.NA(round(x$TE, digits), digits, "NA"),
                        seTE = format.NA(round(x$seTE, digits), digits, "NA"))
    }
    dimnames(res)[[1]] <- x$studlab
    prmatrix(res[order(sortvar),])
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
            digits = digits, header = FALSE,
            backtransf = backtransf, pscale = pscale,
            irscale = irscale, irunit = irunit,
            digits.zval = digits.zval, digits.pval = digits.pval,
            scientific.pval = scientific.pval)
    else
      print(summary(x),
            digits = digits, header = FALSE,
            backtransf = backtransf, pscale = pscale,
            irscale = irscale, irunit = irunit,
            digits.zval = digits.zval, digits.pval = digits.pval,
            scientific.pval = scientific.pval,
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
      if (sm %in% c("PLOGIT", "PLN", "PRAW", "PAS", "PFT")) {
        TE <- pscale * TE
        lowTE <- pscale * lowTE
        uppTE <- pscale * uppTE
      }
      ##
      if (sm %in% c("IR", "IRLN", "IRS", "IRFT")) {
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
          w.fixed.p <- round(100 * x$w.fixed / sum(x$w.fixed, na.rm = TRUE), digits.weight)
        else w.fixed.p <- x$w.fixed
      ##
      if (comb.random)
        if (!is.null(x$w.random) & !all(is.na(x$w.random)) && sum(x$w.random) > 0)
          w.random.p <- round(100 * x$w.random / sum(x$w.random, na.rm = TRUE), digits.weight)
        else w.random.p <- x$w.random
    }
    ##    
    if (metainf.metacum) {
      is.random <- x$pooled == "random"
      ##
      I2 <- format.NA(round(100 * x$I2, digits.I2), digits.I2, "")
      ##
      sel <- is.na(x$p.value)
      p.value <- format.p(x$p.value)
      p.value <- ifelse(sel, "", p.value)
      ##
      tau2 <- x$tau^2
      tau2 <- format.NA(round(tau2, digits.tau2), digits.tau2, "")
      ##
      res <- cbind(format.NA(round(TE, digits), digits, ""),
                   p.ci(format.NA(round(lowTE, digits), digits, "NA"),
                        format.NA(round(uppTE, digits), digits, "NA")),
                   p.value,
                   paste("  ", tau2, sep = ""),
                   paste("  ", I2, ifelse(I2 == "", "", "%"), sep = ""))
      dimnames(res) <- list(paste(x$studlab, "  ", sep = ""),
                            c(sm.lab, ci.lab, "p-value", "tau^2", "I^2"))
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
      catmeth(method = x$method,
              method.tau = if (is.random) x$method.tau else "",
              sm = sm,
              k.all = k.all,
              hakn = is.random & x$hakn,
              metaprop = inherits(x, "metaprop"),
              trimfill = inherits(x, "trimfill"),
              tau.preset = x$tau.preset,
              method.smd = x$method.smd,
              sd.glass = x$sd.glass,
              exact.smd = x$exact.smd,
              model.glmm = x$model.glmm)
    }
    else {
      res <- cbind(format.NA(round(TE, digits), digits, "NA"),
                   p.ci(format.NA(round(lowTE, digits), digits, "NA"),
                        format.NA(round(uppTE, digits), digits, "NA")),
                   if (comb.fixed) format.NA(w.fixed.p, digits.weight),
                   if (comb.random) format.NA(w.random.p, digits.weight),
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
                                        if (comb.fixed) "%W(fixed)",
                                        if (comb.random) "%W(random)",
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
                            if (comb.fixed) "%W(fixed)",
                            if (comb.random) "%W(random)",
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
      cat("\n")
      print(summary(x, warn = FALSE), digits = digits,
            comb.fixed = comb.fixed, comb.random = comb.random,
            prediction = prediction,
            header = FALSE,
            backtransf = backtransf, pscale = pscale,
            irscale = irscale, irunit = irunit,
            digits.tau2 = digits.tau2,
            digits.zval = digits.zval,
            digits.pval = digits.pval,
            digits.pval.Q = digits.pval.Q,
            digits.Q = digits.Q,
            digits.H = digits.H,
            digits.I2 = digits.I2,
            scientific.pval = scientific.pval)
    }
  }
  
  
  invisible(NULL)
}
