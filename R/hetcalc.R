hetcalc <- function(TE, seTE,
                    method.tau, method.tau.ci, TE.tau,
                    method.I2, level,
                    subgroup, control,
                    cluster = NULL, rho = 0) {
  
  Ccalc <- function(x) {
    res <- (sum(x, na.rm = TRUE) -
            sum(x^2, na.rm = TRUE) /
            sum(x, na.rm = TRUE))
    ##
    res
  }
  
  
  by <- !missing(subgroup)
  ##
  sel.noInf <- !is.infinite(TE) & !is.infinite(seTE)
  TE <- TE[sel.noInf]
  seTE <- seTE[sel.noInf]
  if (!is.null(cluster))
    cluster <- cluster[sel.noInf]
  if (by)
    subgroup <- subgroup[sel.noInf]
  ##
  sel.noNA <- !(is.na(TE) | is.na(seTE))
  TE <- TE[sel.noNA]
  seTE <- seTE[sel.noNA]
  if (!is.null(cluster))
    cluster <- cluster[sel.noNA]
  if (by)
    subgroup <- subgroup[sel.noNA]
  ##
  noHet <- all(!sel.noNA) || sum(sel.noNA) < 2
  allNA <- all(!sel.noNA)
  ##
  three.level <- FALSE
  ##
  ## Only conduct three-level meta-analysis if variable 'cluster'
  ## contains duplicate values after removing inestimable study
  ## results standard errors
  ##
  if (!is.null(cluster) &&
      length(unique(cluster)) != length(cluster))
    three.level <- TRUE
  
  
  ##
  ## Meta-analysis without subgroups
  ##
  if (!(is.null(TE.tau)) & method.tau == "DL") {
    ##
    ## Mantel-Haenszel estimator to calculate Q and tau (like RevMan 5)
    ##
    w.common <- 1 / seTE^2
    w.common[is.na(w.common)] <- 0
    ##
    Q <- sum(w.common * (TE - TE.tau)^2, na.rm = TRUE)
    df.Q <- sum(!is.na(seTE)) - 1
    pval.Q <- pvalQ(Q, df.Q)
    #
    if (df.Q == 0)
      tau2 <- NA
    else if (round(Q, digits = 18) <= df.Q)
      tau2 <- 0
    else
      tau2 <- (Q - df.Q) / Ccalc(w.common)
    ##
    se.tau2 <- lower.tau2 <- upper.tau2 <- NA
    tau <- sqrt(tau2)
    lower.tau <- upper.tau <- NA
    ##
    sign.lower.tau <- sign.upper.tau <- method.tau.ci <- ""
    #
    H  <- calcH(Q, df.Q, level)
    I2 <- isquared(Q, df.Q, level)
  }
  else {
    if (noHet) {
      if (allNA)
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
      #
      H  <- calcH(Q, df.Q, level)
      I2 <- isquared(Q, df.Q, level)
    }
    else {
      if (!three.level) {
        mf0 <- runNN(rma.uni,
                     list(yi = TE, sei = seTE, method = method.tau,
                          control = control))
        ##
        tau2 <- mf0$tau2
        se.tau2 <- mf0$se.tau2
        ## Calculate Cochran's Q
        w <- 1 / seTE^2
        Q <- sum(w * (TE - weighted.mean(TE, w))^2)
      }
      else {
        idx <- seq_along(TE)
        mf0 <-
          runNN(rma.mv,
                list(yi = TE,
                     V = vcalc(vi = seTE^2, cluster = cluster, obs = idx,
                               rho = rho),
                     method = method.tau,
                     random = as.call(~ 1 | cluster / idx),
                     control = control,
                     data = data.frame(cluster, idx)),
                warn = FALSE)
        ##
        tau2 <- mf0$sigma2
        se.tau2 <- NA
        ##
        Q <- mf0$QE
      }
      ##
      tau <- sqrt(tau2)
      ##
      df.Q <- mf0$k - mf0$p
      pval.Q <- pvalQ(Q, df.Q)
      #
      if (df.Q < 2)
        method.tau.ci <- ""
      else if (three.level & method.tau.ci != "")
        method.tau.ci <- "PL"
      ##
      ## Confidence interval for overall meta-analysis
      ##
      if (method.tau.ci == "BJ")
        ci0 <-
          confint.rma.uni(
            runNN(rma.uni,
                  list(yi = TE, sei = seTE, weights = 1 / seTE^2,
                       method = "GENQ", control = control)),
            level = 100 * level)
      else if (method.tau.ci == "J")
        ci0 <-
          confint.rma.uni(
            runNN(rma.uni,
                  list(yi = TE, sei = seTE, weights = 1 / seTE,
                       method = "GENQ", control = control)),
            level = 100 * level)
      else if (method.tau.ci == "QP")
        ci0 <- confint.rma.uni(mf0, level = 100 * level)
      else if (method.tau.ci == "PL")
        ci0 <- confint.rma.mv(mf0, level = 100 * level)
      #
      if (method.I2 == "Q") {
        H  <- calcH(Q, df.Q, level)
        I2 <- isquared(Q, df.Q, level)
      }
      else {
        H <- list(TE = sqrt(mf0$H2), lower = NA, upper = NA)
        I2 <- list(TE = mf0$I2 / 100, lower = NA, upper = NA)
      }
    }
  }
  
  
  ##
  ## Meta-analysis with subgroups
  ##
  useFE <- FALSE
  ##
  if (by) {
    if (is.numeric(subgroup))
      subgroup <- as.factor(subgroup)
    ##
    if (!three.level) {
      if (length(unique(subgroup)) == 1)
        mf1 <-
          runNN(rma.uni,
                list(yi = TE, sei = seTE, method = method.tau,
                     control = control))
      else {
        mf1 <-
          try(runNN(rma.uni,
                    list(yi = TE, sei = seTE, method = method.tau,
                         mods = as.call(~ subgroup), control = control,
                         data = data.frame(TE, seTE, subgroup))),
              silent = TRUE)
        ##
        if ("try-error" %in% class(mf1))
          if (grepl(paste0("Number of parameters to be estimated is ",
                           "larger than the number of observations"),
                    mf1)) {
            useFE <- TRUE
            mf1 <-
              runNN(rma.uni,
                    list(yi = TE, sei = seTE, method = "FE",
                         mods = as.call(~ subgroup), control = control,
                         data = data.frame(TE, seTE, subgroup)))
          }
          else
            stop(mf1)
      }
      ##
      tau2.resid <- mf1$tau2
      se.tau2.resid <- mf1$se.tau2
    }
    else {
      idx <- seq_along(TE)
      ##
      if (length(unique(subgroup)) == 1)
        mf1 <-
          runNN(rma.mv,
                list(yi = TE,
                     V = vcalc(vi = seTE^2, cluster = cluster, obs = idx,
                               rho = rho),
                     method = method.tau,
                     random = as.call(~ 1 | cluster / idx),
                     control = control,
                     data = data.frame(cluster, idx)),
                warn = FALSE)
      else {
        mf1 <-
          try(
            runNN(rma.mv,
                  list(yi = TE,
                       V = vcalc(vi = seTE^2, cluster = cluster, obs = idx,
                                 rho = rho),
                       method = method.tau,
                       random = as.call(~ 1 | cluster / idx),
                       mods = as.call(~ subgroup), control = control,
                       data = data.frame(TE, seTE, subgroup, cluster, idx)),
                  warn = FALSE),
            silent = TRUE)
        ##
        if ("try-error" %in% class(mf1))
          if (grepl(paste0("Number of parameters to be estimated is ",
                           "larger than the number of observations"),
                    mf1)) {
            useFE <- TRUE
            mf1 <-
              runNN(rma.mv,
                    list(yi = TE,
                         V = vcalc(vi = seTE^2, cluster = cluster, obs = idx,
                                   rho = rho),
                         method = "FE",
                         random = as.call(~ 1 | cluster / idx),
                         mods = as.call(~ subgroup), control = control,
                         data = data.frame(TE, seTE, subgroup, cluster, idx)),
                    warn = FALSE)
          }
          else
            stop(mf1)
      }
      ##
      tau2.resid <- mf1$sigma2
      se.tau2.resid <- NA
    }
    ##
    tau.resid <- sqrt(tau2.resid)
    ##
    Q.resid <- mf1$QE
    df.Q.resid <- mf1$k - mf1$p
    pval.Q.resid <- pvalQ(Q.resid, df.Q.resid)
    #
    if (method.I2 == "Q") {
      H.resid <- calcH(Q.resid, df.Q.resid, level)
      I2.resid <- isquared(Q.resid, df.Q.resid, level)
    }
    else {
      H.resid <- list(TE = sqrt(mf1$H2), lower = NA, upper = NA)
      I2.resid <- list(TE = mf1$I2 / 100, lower = NA, upper = NA)
    }
    ##
    if (df.Q < 2 || useFE)
      method.tau.ci <- ""
    else if (three.level & method.tau.ci != "")
      method.tau.ci <- "PL"
    ##
    ## Confidence interval for residual heterogeneity
    ##
    if (method.tau.ci == "BJ")
      ci1 <-
        confint.rma.uni(
          runNN(rma.uni,
                list(yi = TE, sei = seTE, weights = 1 / seTE^2,
                     method = "GENQ",
                     mods = as.call(~ subgroup), control = control,
                     data = data.frame(TE, seTE, subgroup))),
          level = 100 * level)
    else if (method.tau.ci == "J")
      ci1 <-
        confint.rma.uni(
          runNN(rma.uni,
                list(yi = TE, sei = seTE, weights = 1 / seTE,
                     method = "GENQ",
                     mods = as.call(~ subgroup), control = control,
                     data = data.frame(TE, seTE, subgroup))),
          level = 100 * level)
    else if (method.tau.ci == "QP")
      ci1 <- confint.rma.uni(mf1, level = 100 * level)
    else if (method.tau.ci == "PL")
      ci1 <- confint.rma.mv(mf1, level = 100 * level)
  }
  
  
  ##
  ## Confidence interval for tau2 and tau
  ##
  if (method.tau.ci %in% c("QP", "BJ", "J")) {
    lower.tau2 <- ci0$random["tau^2", "ci.lb"]
    upper.tau2 <- ci0$random["tau^2", "ci.ub"]
    ##
    lower.tau <- ci0$random["tau", "ci.lb"]
    upper.tau <- ci0$random["tau", "ci.ub"]
    ##
    sign.lower.tau <- ci0$lb.sign
    sign.upper.tau <- ci0$ub.sign
  }
  else if (method.tau.ci == "PL") {
    if (any(names(ci0) == "random")) {      
      lower.tau2 <- c(NA, ci0$random["sigma^2.2", "ci.lb"])
      upper.tau2 <- c(NA, ci0$random["sigma^2.2", "ci.ub"])
      ##
      lower.tau <- c(NA, ci0$random["sigma.2", "ci.lb"])
      upper.tau <- c(NA, ci0$random["sigma.2", "ci.ub"])
      ##
      sign.lower.tau <- c("", ci0$lb.sign)
      sign.upper.tau <- c("", ci0$ub.sign)
    }
    else {
      lower.tau2 <- c(ci0[[1]]$random["sigma^2.1", "ci.lb"],
                      ci0[[2]]$random["sigma^2.2", "ci.lb"])
      upper.tau2 <- c(ci0[[1]]$random["sigma^2.1", "ci.ub"],
                      ci0[[2]]$random["sigma^2.2", "ci.ub"])
      ##
      lower.tau <- c(ci0[[1]]$random["sigma.1", "ci.lb"],
                     ci0[[2]]$random["sigma.2", "ci.lb"])
      upper.tau <- c(ci0[[1]]$random["sigma.1", "ci.ub"],
                     ci0[[2]]$random["sigma.2", "ci.ub"])
      ##
      sign.lower.tau <- c(ci0[[1]]$lb.sign, ci0[[2]]$lb.sign)
      sign.upper.tau <- c(ci0[[1]]$ub.sign, ci0[[2]]$ub.sign)
    }
  }
  else {
    lower.tau2 <- upper.tau2 <- lower.tau <- upper.tau <- NA
    sign.lower.tau <- sign.upper.tau <- "" 
  }
  ##
  ## Confidence interval for tau2.resid and tau.resid
  ##
  if (by) {
    if (method.tau.ci %in% c("QP", "BJ", "J")) {
      lower.tau2.resid <- ci1$random["tau^2", "ci.lb"]
      upper.tau2.resid <- ci1$random["tau^2", "ci.ub"]
      ##
      lower.tau.resid <- ci1$random["tau", "ci.lb"]
      upper.tau.resid <- ci1$random["tau", "ci.ub"]
      ##
      sign.lower.tau.resid <- ci1$lb.sign
      sign.upper.tau.resid <- ci1$ub.sign
    }
    else if (method.tau.ci == "PL") {
      if (any(names(ci0) == "random")) {      
        lower.tau2.resid <- c(NA, ci0$random["sigma^2.2", "ci.lb"])
        upper.tau2.resid <- c(NA, ci0$random["sigma^2.2", "ci.ub"])
        ##
        lower.tau.resid <- c(NA, ci0$random["sigma.2", "ci.lb"])
        upper.tau.resid <- c(NA, ci0$random["sigma.2", "ci.ub"])
      }
      else {
        lower.tau2.resid <- c(ci1[[1]]$random["sigma^2.1", "ci.lb"],
                              ci1[[2]]$random["sigma^2.2", "ci.lb"])
        upper.tau2.resid <- c(ci1[[1]]$random["sigma^2.1", "ci.ub"],
                              ci1[[2]]$random["sigma^2.2", "ci.ub"])
        ##
        lower.tau.resid <- c(ci1[[1]]$random["sigma.1", "ci.lb"],
                             ci1[[2]]$random["sigma.2", "ci.lb"])
        upper.tau.resid <- c(ci1[[1]]$random["sigma.1", "ci.ub"],
                             ci1[[2]]$random["sigma.2", "ci.ub"])
      }
    }
    else {
      lower.tau2.resid <- upper.tau2.resid <-
        lower.tau.resid <- upper.tau.resid <- NA
    }
  }
  
  
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
              #
              level = level,
              #
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
              tau2.resid = if (by) tau2.resid else NA,
              se.tau2.resid = if (by) se.tau2.resid else NA,
              lower.tau2.resid = if (by) lower.tau2.resid else NA,
              upper.tau2.resid = if (by) upper.tau2.resid else NA,
              ##
              tau.resid = if (by) tau.resid else NA,
              lower.tau.resid = if (by) lower.tau.resid else NA,
              upper.tau.resid = if (by) upper.tau.resid else NA,
              ##
              Q.resid = if (by) Q.resid else NA,
              df.Q.resid = if (by) df.Q.resid else NA,
              pval.Q.resid = if (by) pval.Q.resid else NA,
              #
              H.resid = if (by) H.resid$TE else NA,
              lower.H.resid = if (by) H.resid$lower else NA,
              upper.H.resid = if (by) H.resid$upper else NA,
              ##
              I2.resid = if (by) I2.resid$TE else NA,
              lower.I2.resid = if (by) I2.resid$lower else NA,
              upper.I2.resid = if (by) I2.resid$upper else NA
              )
  ##
  res
}
