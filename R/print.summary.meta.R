print.summary.meta <- function(x,
                               comb.fixed = x$comb.fixed,
                               comb.random = x$comb.random,
                               prediction = x$prediction,
                               print.byvar = x$print.byvar,
                               byseparator = x$byseparator,
                               print.CMH = x$print.CMH,
                               header = TRUE,
                               backtransf = x$backtransf,
                               pscale = x$pscale,
                               irscale = x$irscale,
                               irunit = x$irunit,
                               bylab.nchar = 35,
                               digits = gs("digits"),
                               digits.zval = gs("digits.zval"),
                               digits.pval = max(gs("digits.pval"), 2),
                               digits.pval.Q = max(gs("digits.pval.Q"), 2),
                               digits.Q = gs("digits.Q"),
                               digits.tau2 = gs("digits.tau2"),
                               digits.H = gs("digits.H"),
                               digits.I2 = gs("digits.I2"),
                               scientific.pval = gs("scientific.pval"),
                               big.mark = gs("big.mark"),
                               print.I2 = gs("print.I2"),
                               print.H = gs("print.H"),
                               print.Rb = gs("print.Rb"),
                               text.tau2 = gs("text.tau2"),
                               text.I2 = gs("text.I2"),
                               text.Rb = gs("text.Rb"),
                               warn.backtransf = FALSE,
                               ...) {
  
  
  ##
  ##
  ## (1) Check for summary.meta object
  ##
  ##
  chkclass(x, "summary.meta")
  ##
  if (inherits(x, "metacum") | inherits(x, "metainf"))
    return(invisible(NULL))
  ##
  by <- !is.null(x$bylab)
  
  
  ##
  ##
  ## (2) Check and set other arguments
  ##
  ##
  chknumeric(digits, min = 0, single = TRUE)
  chknumeric(digits.tau2, min = 0, single = TRUE)
  chknumeric(digits.zval, min = 0, single = TRUE)
  chknumeric(digits.pval, min = 1, single = TRUE)
  chknumeric(digits.pval.Q, min = 1, single = TRUE)
  chknumeric(digits.Q, min = 0, single = TRUE)
  chknumeric(digits.H, min = 0, single = TRUE)
  chknumeric(digits.I2, min = 0, single = TRUE)
  chklogical(scientific.pval)
  ##
  if (is.untransformed(x$sm))
    backtransf <- TRUE
  chklogical(backtransf)
  ##
  chklogical(print.I2)
  chklogical(print.H)
  chklogical(print.Rb)
  chkchar(text.tau2)
  chkchar(text.I2)
  chkchar(text.Rb)
  chklogical(warn.backtransf)
  is.prop <- is.prop(x$sm)
  is.rate <- is.rate(x$sm)
  ##
  if (!is.prop & x$sm != "RD")
    pscale <- 1
  if (!is.null(pscale))
    chknumeric(pscale, single = TRUE)
  else
    pscale <- 1
  if (!backtransf & pscale != 1 & !is.untransformed(x$sm)) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  if (!is.rate & x$sm != "IRD")
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
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(prediction)
  if (by) {
    chklogical(print.byvar)
    chkchar(byseparator)
  }
  if (!is.null(print.CMH))
    chklogical(print.CMH)
  chklogical(header)
  chknumeric(bylab.nchar)
  ##
  ## Additional arguments / checks for metacont objects
  ##
  cl <- paste("update.meta() or ", class(x)[1], "()", sep = "")
  addargs <- names(list(...))
  ##
  fun <- "print.summary.meta"
  ##
  warnarg("logscale", addargs, fun, otherarg = "backtransf")
  
  
  ##
  ##
  ## (3) Some additional settings
  ##
  ##
  k.all <- length(x$study$TE)
  k <- x$k
  sm <- x$sm
  ##
  bip <- inherits(x, c("metabin", "metainc", "metaprop", "metarate"))
  null.given <- (inherits(x, c("metacor", "metagen", "metamean",
                               "metaprop", "metarate")) |
                 is.prop(sm) | is.rate(sm) | is.cor(sm) | is.mean(sm))
  ##
  null.effect <- x$null.effect
  ##
  if (null.given & !backtransf) {
    ##
    if (sm %in% c("PFT", "PAS"))
      null.effect <- asin(sqrt(null.effect))
    else if (is.log.effect(sm))
      null.effect <- log(null.effect)
    else if (sm == c("PLOGIT"))
      null.effect <- log(null.effect / (1 - null.effect))
    else if (sm %in% c("IRS", "IRFT"))
      null.effect <- sqrt(null.effect)
    else if (sm == "ZCOR")
      null.effect <- 0.5 * log((1 + null.effect) / (1 - null.effect))
  }
  ##
  prediction <- prediction & k >= 3
  if (is.na(prediction))
    prediction <- FALSE
  ##
  sm.lab <- sm
  ##
  if (backtransf) {
    if (sm == "ZCOR")
      sm.lab <- "COR"
    else if (is.mean(sm))
      sm.lab <- "mean"
    else if (is.prop) {
      if (pscale == 1)
        sm.lab <- "proportion"
      else
        sm.lab <- "events"
    }
    else if (is.rate) {
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
  if (length(x$tau.common) == 0)
    x$tau.common <- FALSE
  ##
  if (length(x$tau.common) == 0)
    x$tau.common <- FALSE
  ##
  if (by)
    bylevs <- ifelse(nchar(x$bylevs) > bylab.nchar,
                     paste(substring(x$bylevs, 1, bylab.nchar - 4), " ...", sep = ""),
                     x$bylevs)
  
  
  ##
  ##
  ## (4) Set and backtransform results of meta-analysis
  ##
  ##
  TE.fixed    <- x$fixed$TE
  lowTE.fixed <- x$fixed$lower
  uppTE.fixed <- x$fixed$upper
  ##
  TE.random    <- x$random$TE
  lowTE.random <- x$random$lower
  uppTE.random <- x$random$upper
  ##
  lowTE.predict <- x$predict$lower
  uppTE.predict <- x$predict$upper
  ##
  Q <- x$Q
  df.Q <- replaceNULL(x$df.Q, k - 1)
  pval.Q <- replaceNULL(x$pval.Q, pvalQ(Q, df.Q))
  ##
  if (!is.null(x$Q.CMH)) {
    Q.CMH <- x$Q.CMH
    df.Q.CMH <- replaceNULL(x$df.Q.CMH, 1)
    pval.Q.CMH <- replaceNULL(x$pval.Q.CMH, pvalQ(Q.CMH, df.Q.CMH))
  }
  ##
  if (x$method == "GLMM") {
    Q.LRT <- x$Q.LRT
    df.Q.LRT <- replaceNULL(x$df.Q.LRT, df.Q)
    pval.Q.LRT <- replaceNULL(x$pval.Q.LRT, pvalQ(Q.LRT, df.Q.LRT))
  }
  ##
  if (by) {
    TE.fixed.w     <- x$within.fixed$TE
    lowTE.fixed.w  <- x$within.fixed$lower
    uppTE.fixed.w  <- x$within.fixed$upper
    pval.fixed.w   <- x$within.fixed$p
    harmonic.mean.w <- x$within.fixed$harmonic.mean
    TE.random.w    <- x$within.random$TE
    lowTE.random.w <- x$within.random$lower
    uppTE.random.w <- x$within.random$upper
    pval.random.w   <- x$within.random$p
    ##
    Q.b.fixed <- x$Q.b.fixed
    Q.w.fixed <- x$Q.w.fixed
    Q.b.random <- x$Q.b.random
    Q.w.random <- x$Q.w.random
    ##
    Q.w <- x$Q.w
    ##
    k.w <- x$k.w
    ##
    df.Q.w <- replaceNULL(x$df.Q.w, sum((k.w - 1)[!is.na(x$Q.w)]))
    df.Q.b <- replaceNULL(x$df.Q.b, (k - 1) - sum((k.w - 1)[!is.na(x$Q.w)]))
    ##
    pval.Q.b.fixed  <- replaceNULL(x$pval.Q.b.fixed, pvalQ(Q.b.fixed, df.Q.b))
    pval.Q.w.fixed  <- replaceNULL(x$pval.Q.w.fixed, pvalQ(Q.w.fixed, df.Q.w))
    pval.Q.b.random <- replaceNULL(x$pval.Q.b.random, pvalQ(Q.b.random, df.Q.b))
    pval.Q.w.random <- replaceNULL(x$pval.Q.w.random, pvalQ(Q.w.random, df.Q.w))
  }
  ##
  if (backtransf) {
    if (sm %in% c("IR", "IRLN", "IRS", "IRFT"))
      harmonic.mean <- 1 / mean(1 / x$time)
    else
      harmonic.mean <- 1 / mean(1 / x$n)
    ##
    TE.fixed    <- backtransf(TE.fixed, sm, "mean",
                              harmonic.mean, warn = comb.fixed & warn.backtransf)
    lowTE.fixed <- backtransf(lowTE.fixed, sm, "lower",
                              harmonic.mean, warn = comb.fixed & warn.backtransf)
    uppTE.fixed <- backtransf(uppTE.fixed, sm, "upper",
                              harmonic.mean, warn = comb.fixed & warn.backtransf)
    ##
    TE.random <- backtransf(TE.random, sm, "mean",
                            harmonic.mean, warn = comb.random & warn.backtransf)
    lowTE.random <- backtransf(lowTE.random, sm, "lower",
                               harmonic.mean, warn = comb.random & warn.backtransf)
    uppTE.random <- backtransf(uppTE.random, sm, "upper",
                               harmonic.mean, warn = comb.random & warn.backtransf)
    ##
    lowTE.predict <- backtransf(lowTE.predict, sm, "lower",
                                harmonic.mean, warn = prediction & warn.backtransf)
    uppTE.predict <- backtransf(uppTE.predict, sm, "upper",
                                harmonic.mean, warn = prediction & warn.backtransf)
    ##
    if (by) {
      TE.fixed.w     <- backtransf(TE.fixed.w, sm, "mean",
                                   harmonic.mean.w, warn = comb.fixed & warn.backtransf)
      lowTE.fixed.w  <- backtransf(lowTE.fixed.w, sm, "lower",
                                   harmonic.mean.w, warn = comb.fixed & warn.backtransf)
      uppTE.fixed.w  <- backtransf(uppTE.fixed.w, sm, "upper",
                                   harmonic.mean.w, warn = comb.fixed & warn.backtransf)
      ##
      TE.random.w    <- backtransf(TE.random.w, sm, "mean",
                                   harmonic.mean.w, warn = comb.random & warn.backtransf)
      lowTE.random.w <- backtransf(lowTE.random.w, sm, "lower",
                                   harmonic.mean.w, warn = comb.random & warn.backtransf)
      uppTE.random.w <- backtransf(uppTE.random.w, sm, "upper",
                                   harmonic.mean.w, warn = comb.random & warn.backtransf)
    }
  }
  ##
  ## Apply argument 'pscale' to proportions / risk differences and
  ## 'irscale' to rates / incidence rate differences
  ##
  if (is.prop | sm == "RD" | is.rate | sm == "IRD") {
    if (is.prop | sm == "RD")
      scale <- pscale
    else if (is.rate | sm == "IRD")
      scale <- irscale
    ##
    TE.fixed    <- scale * TE.fixed
    lowTE.fixed <- scale * lowTE.fixed
    uppTE.fixed <- scale * uppTE.fixed
    ##
    TE.random    <- scale * TE.random
    lowTE.random <- scale * lowTE.random
    uppTE.random <- scale * uppTE.random
    ##
    lowTE.predict <- scale * lowTE.predict
    uppTE.predict <- scale * uppTE.predict
    ##
    if (by) {
      TE.fixed.w    <- scale * TE.fixed.w
      lowTE.fixed.w <- scale * lowTE.fixed.w
      uppTE.fixed.w <- scale * uppTE.fixed.w
      ##   
      TE.random.w    <- scale * TE.random.w
      lowTE.random.w <- scale * lowTE.random.w
      uppTE.random.w <- scale * uppTE.random.w
    }
  }
  ##
  ## Round and round ...
  ##
  TE.fixed    <- round(TE.fixed, digits)
  lowTE.fixed <- round(lowTE.fixed, digits)
  uppTE.fixed <- round(uppTE.fixed, digits)
  pTE.fixed <- x$fixed$p
  zTE.fixed <- round(x$fixed$z, digits.zval)
  ##
  TE.random    <- round(TE.random, digits)
  lowTE.random <- round(lowTE.random, digits)
  uppTE.random <- round(uppTE.random, digits)
  pTE.random <- x$random$p
  zTE.random <- round(x$random$z, digits.zval)
  ##
  lowTE.predict <- round(lowTE.predict, digits)
  uppTE.predict <- round(uppTE.predict, digits)
  ##
  if (by) {
    TE.fixed.w     <- round(TE.fixed.w, digits)
    lowTE.fixed.w  <- round(lowTE.fixed.w, digits)
    uppTE.fixed.w  <- round(uppTE.fixed.w, digits)
    ##
    TE.random.w    <- round(TE.random.w, digits)
    lowTE.random.w <- round(lowTE.random.w, digits)
    uppTE.random.w <- round(uppTE.random.w, digits)
    ##
    if (print.I2)
      I2.w <- round(100 * x$I2.w$TE, digits.I2)
    ##
    if (print.Rb)
      Rb.w <- round(100 * x$Rb.w$TE, digits.I2)
  }
  ##
  if (print.H) {
    H <- round(x$H$TE, digits.H)
    lowH <- round(x$H$lower, digits.H)
    uppH <- round(x$H$upper, digits.H)
  }
  ##
  if (print.I2) {
    I2 <- round(100 * x$I2$TE, digits.I2)
    lowI2 <- round(100 * x$I2$lower, digits.I2)
    uppI2 <- round(100 * x$I2$upper, digits.I2)
    print.ci.I2 <- ((Q > k & k >= 2) | (Q <= k & k > 2)) &
      !(is.na(lowI2) | is.na(uppI2))
    if (is.na(print.ci.I2))
      print.ci.I2 <- FALSE
  }
  else
    print.ci.I2 <- FALSE
  ##
  if (print.Rb) {
    Rb <- round(100 * x$Rb$TE, digits.I2)
    lowRb <- round(100 * x$Rb$lower, digits.I2)
    uppRb <- round(100 * x$Rb$upper, digits.I2)
  }
  
  
  ##
  ##
  ## (5) Print result for meta-analysis
  ##
  ##
  if (header)
    crtitle(x)
  ##
  if (is.na(k.all)) {
    ## Do nothing
  }
  else if (k.all == 1) {
    ##
    ## Print results for a single study
    ##
    res <- cbind(formatN(TE.fixed, digits, "NA",
                         big.mark = big.mark),
                 formatCI(formatN(lowTE.fixed, digits, "NA",
                                  big.mark = big.mark),
                          formatN(uppTE.fixed, digits, "NA",
                                  big.mark = big.mark)),
                 formatN(zTE.fixed, digits.zval, big.mark = big.mark),
                 formatPT(pTE.fixed, digits = digits.pval,
                          scientific = scientific.pval))
    dimnames(res) <- list("", c(sm.lab, x$ci.lab, "z", "p-value"))
    prmatrix(res, quote = FALSE, right = TRUE, ...)
    ## Print information on summary method:
    catmeth(class = class(x),
            method = x$method,
            sm = sm,
            k.all = k.all,
            sparse = ifelse(bip, x$sparse, FALSE),
            incr = if (bip) x$incr else FALSE,
            allincr = ifelse(bip, x$allincr, FALSE),
            addincr = ifelse(bip, x$addincr, FALSE),
            allstudies = x$allstudies,
            doublezeros = x$doublezeros,
            method.ci = x$method.ci,
            pooledvar = x$pooledvar,
            method.smd = x$method.smd,
            sd.glass = x$sd.glass,
            exact.smd = x$exact.smd,
            model.glmm = x$model.glmm,
            pscale = pscale,
            irscale = irscale,
            irunit = irunit,
            null.effect = if (null.given) null.effect else 0,
            big.mark = big.mark,
            digits = digits, digits.tau2 = digits.tau2)
  }
  else if (is.na(k)) {
    ## Do nothing
  }
  else {
    ##
    ## Print results for meta-analysis with more than one study
    ##
    if (comb.fixed | comb.random | prediction) {
      if (!inherits(x, "trimfill")) {
        if (x$method == "MH" &&
            (inherits(x, c("metabin", "metainc")) &
             comb.fixed & sm %in% c("RD", "IRD") &
             (!is.null(x$k.MH) == 1 && k != x$k.MH)))
          cat(paste("Number of studies combined:   k.MH = ", x$k.MH,
                    " (fixed effect), k = ", format(k, big.mark = big.mark),
                    " (random effects)\n\n", sep = ""))
        else
          cat(paste("Number of studies combined: k = ",
                    format(k, big.mark = big.mark), "\n\n", sep = ""))
      }
      else
        cat(paste("Number of studies combined: k = ",
                  format(k, big.mark = big.mark),
                  " (with ",
                  format(x$k0, big.mark = big.mark),
                  " added studies)\n\n", sep = ""))
      ##
      res <- cbind(formatN(c(if (comb.fixed) TE.fixed,
                             if (comb.random) TE.random,
                             if (prediction) NA),
                           digits, "NA",
                           big.mark = big.mark),
                   formatCI(formatN(c(if (comb.fixed) lowTE.fixed,
                                      if (comb.random) lowTE.random,
                                      if (prediction) lowTE.predict),
                                    digits, "NA", big.mark = big.mark),
                            formatN(c(if (comb.fixed) uppTE.fixed,
                                      if (comb.random) uppTE.random,
                                      if (prediction) uppTE.predict),
                                    digits, "NA", big.mark = big.mark)),
                   formatN(c(if (comb.fixed) zTE.fixed,
                             if (comb.random) zTE.random,
                             if (prediction) NA),
                           digits = digits.zval, big.mark = big.mark),
                   formatPT(c(if (comb.fixed) pTE.fixed,
                              if (comb.random) pTE.random,
                              if (prediction) NA),
                            digits = digits.pval,
                            scientific = scientific.pval))
      if (prediction)
        res[dim(res)[1], c(1,3:4)] <- ""
      if (!is.null(x$hakn) && x$hakn) {
        if (comb.fixed & comb.random)
          zlab <- "z|t"
        else if (comb.fixed & !comb.random)
          zlab <- "z"
        else if (!comb.fixed & comb.random)
          zlab <- "t"
      }
      else
        zlab <- "z"
      ##
      dimnames(res) <- list(c(if (comb.fixed) "Fixed effect model",
                              if (comb.random) "Random effects model",
                              if (prediction) "Prediction interval"),  
                            c(sm.lab, x$ci.lab, zlab, "p-value"))
      prmatrix(res, quote = FALSE, right = TRUE, ...)
      ##
      if (inherits(x, "metabin") && print.CMH) {
        Qdata <- cbind(formatN(round(Q.CMH, digits.Q), digits.Q, "NA",
                               big.mark = big.mark),
                       df.Q.CMH,
                       formatPT(pval.Q.CMH,
                                digits = digits.pval.Q,
                                scientific = scientific.pval))
        dimnames(Qdata) <- list("", c("Q", "d.f.", "p-value"))
        ##
        cat("\nCochran-Mantel-Haenszel (CMH) test for overall effect: \n")
        prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
      }
    }
    else
      cat(paste("Number of studies: k = ", k, "\n", sep = ""))
    ##
    ## Print information on heterogeneity
    ##
    if (k.all > 1)
      cat(paste("\nQuantifying heterogeneity:\n",
                ##
                formatPT(x$tau^2,
                         lab = TRUE, labval = text.tau2,
                         digits = digits.tau2,
                         lab.NA = "NA",
                         big.mark = big.mark),
                ##
                if (print.H)
                  paste("; H = ",
                        if (is.na(H)) "NA" else formatN(H, digits.H, "NA",
                                                        big.mark = big.mark),
                        ifelse(k > 2 & !(is.na(lowH) | is.na(uppH)),
                               paste(" ", formatCI(formatN(lowH, digits.H,
                                                           big.mark = big.mark),
                                                   formatN(uppH, digits.H,
                                                           big.mark = big.mark)),
                                     sep = ""),
                               ""),
                        sep = ""),
                ##
                if (print.I2)
                  paste("; ", text.I2, " = ",
                        if (is.na(I2)) "NA" else paste(formatN(I2, digits.I2),
                                                       "%", sep = ""),
                        if (print.ci.I2)
                          paste(" ",
                                formatCI(paste(formatN(lowI2, digits.I2),
                                               "%", sep = ""),
                                         paste(formatN(uppI2, digits.I2),
                                               "%", sep = "")),
                                sep = ""),
                        sep = ""),
                ##
                if (print.Rb)
                  paste("; ",
                        text.Rb, " = ",
                        if (is.na(Rb)) "NA" else paste(formatN(Rb, digits.I2,
                                                               big.mark = big.mark),
                                                       "%", sep = ""),
                        ifelse(k > 2 & !(is.na(lowRb) | is.na(uppRb)),
                               paste(" ",
                                     formatCI(paste(formatN(lowRb, digits.I2,
                                                            big.mark = big.mark),
                                                    "%", sep = ""),
                                              paste(formatN(uppRb, digits.I2,
                                                            big.mark = big.mark),
                                                    "%", sep = "")),
                                     sep = ""),
                               ""),
                        sep = ""),
                "\n", sep = "")
          )
    ##
    if (k.all > 1 & (comb.fixed|comb.random)) {
      if (k > 1) {
        if (x$method != "GLMM") {
          Qdata <- cbind(formatN(round(Q, digits.Q), digits.Q, "NA",
                                 big.mark = big.mark),
                         format(df.Q, big.mark = big.mark),
                         formatPT(pval.Q,
                                  digits = digits.pval.Q,
                                  scientific = scientific.pval))
          dimnames(Qdata) <- list("", c("Q", "d.f.", "p-value"))
        }
        else {
          Qdata <- cbind(formatN(round(c(Q, Q.LRT), digits.Q), digits.Q, "NA",
                                 big.mark = big.mark),
                         format(c(df.Q, df.Q.LRT), big.mark = big.mark),
                         formatPT(c(pval.Q, pval.Q.LRT),
                                  digits = digits.pval.Q,
                                  scientific = scientific.pval),
                         c("Wald-type", "Likelihood-Ratio"))
          dimnames(Qdata) <- list(rep("", 2),
                                  c("Q", "d.f.", "p-value", "Test"))
        }
        ##
        cat("\nTest of heterogeneity:\n")
        prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
      }
      ##
      if (by) {
        ##
        ## Print information for subgroup analysis
        ##
        if (comb.fixed) {
          ##
          ## Subgroup analysis based on fixed effect model
          ##
          Tdata <- cbind(format(k.w, big.mark = big.mark),
                         formatN(TE.fixed.w, digits, "NA",
                                 big.mark = big.mark),
                         formatCI(formatN(lowTE.fixed.w, digits, "NA",
                                          big.mark = big.mark),
                                  formatN(uppTE.fixed.w, digits, "NA",
                                          big.mark = big.mark)),
                         formatN(round(Q.w, digits.Q), digits.Q,
                                 big.mark = big.mark),
                         ifelse(k.w == 1, "--",
                                formatPT(x$tau.w^2,
                                         digits = digits.tau2,
                                         big.mark = big.mark,
                                         noblanks = TRUE)),
                         if (print.I2)
                           ifelse(is.na(I2.w),
                                  "--",
                                  paste(formatN(I2.w, digits.I2),
                                        "%", sep = "")),
                         if (print.Rb)
                           ifelse(is.na(Rb.w),
                                  "--",
                                  paste(formatN(Rb.w, digits.I2),
                                        "%", sep = ""))
                         )
          ##
          bylab <- bylabel(x$bylab, bylevs, print.byvar, byseparator,
                           big.mark = big.mark)
          ##
          dimnames(Tdata) <- list(bylab,
                                  c("  k", sm.lab, x$ci.lab,
                                    "Q", text.tau2,
                                    if (print.I2) text.I2,
                                    if (print.Rb) text.Rb)
                                  )
          if (inherits(x, "metabind"))
            cat("\nResults for meta-analyses (fixed effect model):\n")
          else
            cat("\nResults for subgroups (fixed effect model):\n")
          prmatrix(Tdata, quote = FALSE, right = TRUE, ...)
          ##
          if (!inherits(x, "metabind")) {
            cat("\nTest for subgroup differences (fixed effect model):\n")
            if (x$method == "MH") {
              Qdata <- cbind(formatN(round(Q.b.fixed, digits.Q), digits.Q, "NA",
                                     big.mark = big.mark),
                             format(df.Q.b, big.mark = big.mark),
                             formatPT(pval.Q.b.fixed,
                                      digits = digits.pval.Q,
                                      scientific = scientific.pval))
              dimnames(Qdata) <- list("Between groups  ",
                                      c("Q", "d.f.", "p-value"))
              prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
            }
            else {
              Qs  <- c(Q.b.fixed, Q.w.fixed)
              dfs <- c(df.Q.b, df.Q.w)
              pvals <- c(pval.Q.b.fixed, pval.Q.w.fixed)
              Qdata <- cbind(formatN(round(Qs, digits.Q), digits.Q, "NA",
                                     big.mark = big.mark),
                             format(dfs, big.mark = big.mark),
                             formatPT(pvals,
                                      digits = digits.pval.Q,
                                      scientific = scientific.pval))
              dimnames(Qdata) <- list(c("Between groups", "Within groups"),
                                      c("Q", "d.f.", "p-value"))
              prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
            }
          }
        }
        ##
        if (comb.random) {
          ##
          ## Subgroup analysis based on random effects model
          ##
          Tdata <- cbind(format(k.w, big.mark = big.mark),
                         formatN(TE.random.w, digits, "NA",
                                 big.mark = big.mark),
                         formatCI(formatN(lowTE.random.w, digits, "NA",
                                          big.mark = big.mark),
                                  formatN(uppTE.random.w, digits, "NA",
                                          big.mark = big.mark)),
                         formatN(round(Q.w, digits.Q), digits.Q,
                                 big.mark = big.mark),
                         ifelse(k.w == 1, "--",
                                formatPT(x$tau.w^2,
                                         digits = digits.tau2,
                                         big.mark = big.mark,
                                         noblanks = TRUE)),
                         if (print.I2)
                           ifelse(is.na(I2.w),
                                  "--",
                                  paste(formatN(I2.w, digits.I2),
                                        "%", sep = "")),
                         if (print.Rb)
                           ifelse(is.na(Rb.w),
                                  "--",
                                  paste(formatN(Rb.w, digits.I2,
                                                big.mark = big.mark),
                                        "%", sep = ""))
                         )
          ##
          bylab <- bylabel(x$bylab, bylevs, print.byvar, byseparator,
                           big.mark = big.mark)
          ##
          dimnames(Tdata) <- list(bylab,
                                  c("  k", sm.lab, x$ci.lab,
                                    "Q", text.tau2,
                                    if (print.I2) text.I2,
                                    if (print.Rb) text.Rb)
                                  )
          if (inherits(x, "metabind"))
            cat("\nResults for meta-analyses (random effects model):\n")
          else
            cat("\nResults for subgroups (random effects model):\n")
          prmatrix(Tdata, quote = FALSE, right = TRUE, ...)
          ##
          if (!inherits(x, "metabind")) {
            cat("\nTest for subgroup differences (random effects model):\n")
            if (is.na(Q.w.random)) {
              Qdata <- cbind(formatN(round(Q.b.random, digits.Q), digits.Q,
                                     "NA", big.mark = big.mark),
                             format(df.Q.b, big.mark = big.mark),
                             formatPT(pval.Q.b.random,
                                      digits = digits.pval.Q,
                                      scientific = scientific.pval))
              dimnames(Qdata) <- list("Between groups  ",
                                      c("Q", "d.f.", "p-value"))
            }
            else {
              Qs  <- c(Q.b.random, Q.w.random)
              dfs <- c(df.Q.b, df.Q.w)
              pvals <- c(pval.Q.b.random, pval.Q.w.random)
              Qdata <- cbind(formatN(round(Qs, digits.Q), digits.Q, "NA",
                                     big.mark = big.mark),
                             format(dfs, big.mark = big.mark),
                             formatPT(pvals,
                                      digits = digits.pval.Q,
                                      scientific = scientific.pval))
              dimnames(Qdata) <- list(c("Between groups", "Within groups"),
                                      c("Q", "d.f.", "p-value"))
            }
            prmatrix(Qdata, quote = FALSE, right = TRUE, ...)
          }
        }
      }
    }
  }
  
  
  if ((comb.fixed | comb.random | prediction) & (is.na(k.all) | k.all > 1)) {
    ##
    ## Print information on summary method:
    ##
    catmeth(class = class(x),
            method = x$method,
            method.tau = if (comb.random) x$method.tau else "",
            sm = sm,
            k.all = k.all,
            hakn = !is.null(x$hakn) && (x$hakn & comb.random),
            tau.common = by & x$tau.common,
            tau.preset = x$tau.preset,
            sparse = ifelse(bip, x$sparse, FALSE),
            incr = if (bip) x$incr else FALSE,
            allincr = ifelse(bip, x$allincr, FALSE),
            addincr = ifelse(bip, x$addincr, FALSE),
            allstudies = x$allstudies,
            doublezeros = x$doublezeros,
            MH.exact = ifelse(inherits(x, "metabin"), x$MH.exact, FALSE),
            method.ci = x$method.ci,
            pooledvar = x$pooledvar,
            method.smd = x$method.smd,
            sd.glass = x$sd.glass,
            exact.smd = x$exact.smd,
            model.glmm = x$model.glmm,
            pscale = pscale,
            irscale = irscale,
            irunit = irunit,
            null.effect = if (null.given) null.effect else 0,
            big.mark = big.mark,
            digits = digits, digits.tau2 = digits.tau2)
  }
  
  
  invisible(NULL)
}
