gm <- function(x, digits = 4, debug = FALSE) {
  
  if (inherits(x, c("metabin", "metainc", "metaprop", "metarate"))) {
    incr_not_0 <- sort(unique(x$incr[x$incr != 0]))
    #
    if (length(incr_not_0) == 1)
      incr <- incr_not_0
    else if (length(incr_not_0) == 0)
      incr <- 0
    else
      incr <- paste0("{", paste(incr_not_0, collapse = ", "), "}")
  }
  
  func <- if (debug) list else data.frame
  ## Get rid of warning 'Undefined global functions or variables'
  model <- NULL
  ##
  if (inherits(x, "metabind")) {
    res <-
      with(
        x$data,
        list(
          meth =
            func(
              model, method, three.level,
              k, k.all,
              k.MH = replaceNULL(x$k.MH, NA),
              k.study, k.TE,
              k0 = replaceNULL(x$k0, NA),
              method.tau, method.tau.ci,
              tau = round(tau, digits), tau.preset,
              method.I2,
              method.random.ci,
              df.random = replaceNULL(df.random),
              adhoc.hakn.ci,
              rho),
          pred =
            func(
              model, method.predict,
              df.predict = replaceNULL(df.predict), adhoc.hakn.pi)
        )
      )
    ##
    res$meth <- subset(res$meth, model %in% c("common", "random"))
    res$pred <- subset(res$pred, model == "predict")
    ##
    return(res)
  }
  else if (inherits(x, c("metacum", "metainf"))) {
    res <-
      with(
        x,
        list(
          meth =
            func(
              model = pooled,
              method =
                if (pooled == "common") method
                else method.random,
              three.level,
              k, k.all, k.MH = replaceNULL(x$k.MH),
              k.study, k.TE,
              k0 = replaceNULL(x$k0),
              method.tau,
              method.tau.ci =
                if (length(method.tau.ci[method.tau.ci != ""]) == 0) ""
                else method.tau.ci[method.tau.ci != ""],
              tau = NA, tau.preset = replaceNULL(tau.preset),
              method.I2,
              method.random.ci,
              df.random = replaceNULL(df.random), adhoc.hakn.ci,
              rho = replaceNULL(rho)),
          pred = data.frame()
        )
      )
    ##
    res$meth <- subset(res$meth, model %in% c("common", "random"))
    ##
    ## Add row with information on random effects model (to print
    ## information on tau2)
    ##
    if (all(res$meth$model == "common")) {
      res$meth <- rbind(res$meth, res$meth[nrow(res$meth), ])
      res$meth$model[nrow(res$meth)] <- "random"
    }
  }
  else {
    ##
    n.com <- length(x$lower.common)
    n.ran <- length(x$lower.random)
    n.prd <- length(x$lower.predict)
    ##
    x$method <- expandvar(x$method, n.com, 1)
    ##
    x$method.random <- expandvar(x$method.random, n.ran, 1)
    x$method.random.ci <- expandvar(x$method.random.ci, n.ran, 1)
    x$method.tau <- expandvar(x$method.tau, n.ran, 1)
    x$method.tau.ci <- expandvar(x$method.tau.ci, n.ran, 1)
    if (length(x$three.level) == 1 && x$three.level)
      x$tau <- sum(x$tau)
    x$three.level <- expandvar(x$three.level, n.ran, 1)
    x$tau <- expandvar(x$tau, n.ran, 1)
    x$adhoc.hakn.ci <- expandvar(x$adhoc.hakn.ci, n.ran, 1)
    x$rho <- expandvar(x$rho, n.ran, 1)
    ##
    x$method.predict <- expandvar(x$method.predict, n.prd, 1)
    ##
    res <-
      list(meth = 
             with(x,
                  func(
                    model = c(rep("common", n.com), rep("random", n.ran)),
                    method = c(method, method.random),
                    three.level = c(rep(FALSE, n.com), x$three.level),
                    ##
                    k = k, k.all = k.all, k.MH = NA, k.study = k.study,
                    k.TE = k.TE, k0 = NA,
                    ##
                    method.tau = c(rep("", n.com), x$method.tau),
                    method.tau.ci = c(rep("", n.com), x$method.tau.ci),
                    tau = round(c(rep(NA, n.com), x$tau), digits),
                    tau.preset = NA,
                    #
                    method.I2 = x$method.I2,
                    #
                    method.random.ci = c(rep("", n.com), method.random.ci),
                    df.random = c(rep(NA, n.com), unlist(df.random)),
                    adhoc.hakn.ci = c(rep("", n.com), adhoc.hakn.ci),
                    ##
                    rho = c(rep(0, n.com), x$rho))),
           pred =
             with(x,
                  data.frame(method.predict = method.predict))
           )
  }
  ##
  if (!is.null(x$tau.preset))
    res$meth$tau.preset <- round(c(rep(NA, n.com), x$tau.preset), digits)
  ##
  if (inherits(x, "metabin")) {
    res$meth$k.MH <- x$k.MH
    res$meth$incr <- incr
    res$meth$method.incr <- x$method.incr
    res$meth$sparse <- x$sparse
    res$meth$allstudies <- x$allstudies
    res$meth$doublezeros <- x$doublezeros
    res$meth$MH.exact <- x$MH.exact
    res$meth$RR.Cochrane <- x$RR.Cochrane
    if (length(x$Q.Cochrane) == 1)
      x$Q.Cochrane <- rep(x$Q.Cochrane, n.ran)
    res$meth$Q.Cochrane <- c(rep(FALSE, n.com), x$Q.Cochrane)
  }
  ##
  if (inherits(x, "metacont")) {
    res$meth$pooledvar <- x$pooledvar
    res$meth$method.smd <- x$method.smd
    res$meth$sd.glass <- x$sd.glass
    res$meth$exact.smd <- x$exact.smd
    res$meth$method.mean <- x$method.mean
    res$meth$method.sd <- x$method.sd
  }
  ##
  if (inherits(x, "metainc"))
    res$meth$k.MH <- x$k.MH
  ##
  if (inherits(x, c("metainc", "metaprop", "metarate"))) {
    res$meth$incr <- incr
    res$meth$method.incr <- x$method.incr
    res$meth$sparse <- x$sparse
  }
  if (!is.null(x$model.glmm))
    res$meth$model.glmm <- x$model.glmm
  #
  if (!is.null(x$phi))
    res$meth$phi <- x$phi
  #
  if (inherits(x, "trimfill")) {
    res$meth$k0 <- x$k0
    res$meth$left <- x$left
    res$meth$ma.common <- x$ma.common
    res$meth$type <- x$type
    res$meth$n.iter.max <- x$n.iter.max
    res$meth$n.iter <- x$n.iter
  }
  if (allNA(res$meth$k.MH))
    res$meth$k.MH <- NULL
  if (allNA(res$meth$k0))
    res$meth$k0 <- NULL
  ##
  if (!is.null(x$df.predict))
    res$pred$df.predict <- unlist(x$df.predict)
  if (!is.null(x$adhoc.hakn.pi))
    res$pred$adhoc.hakn.pi <- x$adhoc.hakn.pi
  ##
  res
}
