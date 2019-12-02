hetcalc <- function(TE, seTE, method.tau, method.tau.ci,
                    TE.tau, level.hetstats, byvar, control) {
  
  
  Ccalc <- function(x) {
    res <- (sum(x, na.rm = TRUE) -
            sum(x^2, na.rm = TRUE) /
            sum(x, na.rm = TRUE))
    ##
    res
  }
  
  
  by <- !missing(byvar)
  
  
  if (by) {
    if (is.numeric(byvar))
      byvar <- as.factor(byvar)
    ##
    mf1 <- rma.uni(yi = TE, sei = seTE, method = method.tau,
                   mods = ~ byvar, control = control)
    ##
    Q <- mf1$QE
    df.Q <- mf1$k - mf1$p
    pval.Q <- pvalQ(Q, df.Q)
    ##
    if (method.tau.ci == "BJ")
      ci1 <- confint.rma.uni(rma.uni(yi = TE, sei = seTE, weights = 1 / seTE^2,
                                     method = "GENQ",
                                     mods = ~ byvar, control = control))
    else if (method.tau.ci == "J")
      ci1 <- confint.rma.uni(rma.uni(yi = TE, sei = seTE, weights = 1 / seTE,
                                     method = "GENQ",
                                     mods = ~ byvar, control = control))
    else if (method.tau.ci == "QP")
      ci1 <- confint.rma.uni(mf1)
    ##
    Q <- mf1$QE
    df.Q <- mf1$k - mf1$p
    pval.Q <- pvalQ(Q, df.Q)
    ##
    tau2 <- mf1$tau2
    se.tau2 <- mf1$se.tau2
    tau <- sqrt(tau2)
    ##
    if (method.tau.ci %in% c("QP", "BJ", "J")) {
      lower.tau2 <- ci1$random["tau^2", "ci.lb"]
      upper.tau2 <- ci1$random["tau^2", "ci.ub"]
      ##
      lower.tau <- ci1$random["tau", "ci.lb"]
      upper.tau <- ci1$random["tau", "ci.ub"]
      ##
      sign.lower.tau <- ci1$lb.sign
      sign.upper.tau <- ci1$ub.sign
    }
    else {
      lower.tau2 <- upper.tau2 <- lower.tau <- upper.tau <- NA
      sign.lower.tau <- sign.upper.tau <- "" 
    }
    ##
    H.resid <- sqrt(mf1$H2)
    lower.H.resid <- upper.H.resid <- NA
    ##
    I2.resid <- mf1$I2 / 100
    lower.I2.resid <- upper.I2.resid <- NA
  }
  else {
    if (!(is.null(TE.tau)) & method.tau == "DL") {
      ##
      ## Mantel-Haenszel estimator to calculate Q and tau (like RevMan 5)
      ##
      w.fixed <- 1 / seTE^2
      w.fixed[is.na(w.fixed)] <- 0
      ##
      Q <- sum(w.fixed * (TE - TE.tau)^2, na.rm = TRUE)
      df.Q <- sum(!is.na(seTE)) - 1
      pval.Q <- pvalQ(Q, df.Q)
      ##
      if (df.Q == 0)
        tau2 <- NA
      else if (round(Q, digits = 18) <= df.Q)
        tau2 <- 0
      else
        tau2 <- (Q - df.Q) / Ccalc(w.fixed)
      ##
      se.tau2 <- lower.tau2 <- upper.tau2 <- NA
      tau <- sqrt(tau2)
      lower.tau <- upper.tau <- NA
      ##
      sign.lower.tau <- sign.upper.tau <- method.tau.ci <- ""
    }
    else {
      sel <- !(is.na(TE) | is.na(seTE))
      ##
      if (all(!sel) || sum(sel) < 2) {
        if (all(!sel))
          Q <- NA
        else
          Q <- 0
        ##
        df.Q <- 0
        pval.Q <- pvalQ(Q, df.Q)
        ##
        tau2 <- NA
        se.tau2 <- lower.tau2 <- upper.tau2 <- NA
        tau <- sqrt(tau2)
        lower.tau <- upper.tau <- NA
        ##
        sign.lower.tau <- sign.upper.tau <- method.tau.ci <- ""
      }
      else {
        mf2 <- rma.uni(yi = TE[sel], sei = seTE[sel],
                       method = method.tau, control = control)
        ##
        if (sum(sel) < 3)
          method.tau.ci <- ""
        ##
        if (method.tau.ci == "BJ")
          ci2 <- confint.rma.uni(rma.uni(yi = TE[sel], sei = seTE[sel],
                                         weights = 1 / seTE[sel]^2,
                                         method = "GENQ", control = control))
        else if (method.tau.ci == "J")
          ci2 <- confint.rma.uni(rma.uni(yi = TE[sel], sei = seTE[sel],
                                         weights = 1 / seTE[sel],
                                         method = "GENQ", control = control))
        else if (method.tau.ci == "QP")
          ci2 <- confint.rma.uni(mf2)
        ##
        Q <- mf2$QE
        df.Q <- mf2$k - mf2$p
        pval.Q <- pvalQ(Q, df.Q)
        ##
        tau2 <- mf2$tau2
        se.tau2 <- mf2$se.tau2
        tau <- sqrt(tau2)
        ##
        if (method.tau.ci %in% c("QP", "BJ", "J")) {
          lower.tau2 <- ci2$random["tau^2", "ci.lb"]
          upper.tau2 <- ci2$random["tau^2", "ci.ub"]
          ##
          lower.tau <- ci2$random["tau", "ci.lb"]
          upper.tau <- ci2$random["tau", "ci.ub"]
          ##
          sign.lower.tau <- ci2$lb.sign
          sign.upper.tau <- ci2$ub.sign
        }
        else {
          lower.tau2 <- upper.tau2 <- lower.tau <- upper.tau <- NA
          sign.lower.tau <- sign.upper.tau <- "" 
        }
      }
    }
  }
  
  
  H  <- calcH(Q, df.Q, level.hetstats)
  I2 <- isquared(Q, df.Q, level.hetstats)
  
  
  res <- list(tau2 = tau2,
              se.tau2 = se.tau2,
              lower.tau2 = lower.tau2,
              upper.tau2 = upper.tau2,
              ##
              tau = tau,
              lower.tau = lower.tau,
              upper.tau = upper.tau,
              ##
              method.tau.ci = method.tau.ci,
              sign.lower.tau = sign.lower.tau,
              sign.upper.tau = sign.upper.tau,
              ##
              Q = Q,
              df.Q = df.Q,
              pval.Q = pval.Q,
              ##
              H = H$TE,
              lower.H = H$lower,
              upper.H = H$upper,
              ##
              I2 = I2$TE,
              lower.I2 = I2$lower,
              upper.I2 = I2$upper,
              ##
              H.resid = if (by) H.resid else NA,
              lower.H.resid = if (by) lower.H.resid else NA,
              upper.H.resid = if (by) upper.H.resid else NA,
              ##
              I2.resid = if (by) I2.resid else NA,
              lower.I2.resid = if (by) lower.I2.resid else NA,
              upper.I2.resid = if (by) upper.I2.resid else NA
              )
  ##
  res
}
