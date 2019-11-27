hetcalc <- function(TE, seTE, method.tau, TE.tau,
                    level.hetstats, byvar, control) {
  
  
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
    ##
    tau2 <- mf1$tau2
    se.tau2 <- mf1$se.tau2
    ##
    H.resid <- sqrt(mf1$H2)
    lower.H.resid <- upper.H.resid <- NA
    ##
    I2.resid <- mf1$I2 / 100
    lower.I2.resid <- upper.I2.resid <- NA
  }
  else {
    if (!(is.null(TE.tau)) & (method.tau == "DL" | method.tau == "PM")) {
      w.fixed <- 1 / seTE^2
      w.fixed[is.na(w.fixed)] <- 0
      ##
      Q <- sum(w.fixed * (TE - TE.tau)^2, na.rm = TRUE)
      df.Q <- sum(!is.na(seTE)) - 1
      ##
      if (df.Q == 0)
        tau2 <- NA
      else if (round(Q, digits = 18) <= df.Q)
        tau2 <- 0
      else
        tau2 <- (Q - df.Q) / Ccalc(w.fixed)
      ##
      se.tau2 <- NULL
    }
    else {
      sel.NA <- is.na(TE) | is.na(seTE)
      ##
      if (sum(sel.NA) == length(TE)) {
        Q <- NA
        df.Q <- 0
        ##
        tau2 <- NA
        se.tau2 <- NULL
      }
      else {
        mf2 <- rma.uni(yi = TE[!sel.NA], sei = seTE[!sel.NA],
                       method = method.tau, control = control)
        Q <- mf2$QE
        df.Q <- mf2$k - mf2$p
        ##
        tau2 <- mf2$tau2
        se.tau2 <- mf2$se.tau2
      }
    }
  }
  
  
  H  <- calcH(Q, df.Q, level.hetstats)
  I2 <- isquared(Q, df.Q, level.hetstats)
  
  
  res <- list(tau = sqrt(tau2),
              se.tau2 = se.tau2,
              Q = Q,
              df.Q = df.Q,
              ##
              H = H$TE,
              lower.H = H$lower,
              upper.H = H$upper,
              ##
              I2 = I2$TE,
              lower.I2 = I2$lower,
              upper.I2 = I2$upper,
              ##
              H.resid = if (by) H.resid else NULL,
              lower.H.resid = if (by) lower.H.resid else NULL,
              upper.H.resid = if (by) upper.H.resid else NULL,
              ##
              I2.resid = if (by) I2.resid else NULL,
              lower.I2.resid = if (by) lower.I2.resid else NULL,
              upper.I2.resid = if (by) upper.I2.resid else NULL
              )
  ##
  res
}
