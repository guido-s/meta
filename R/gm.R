gm <- function(x, digits = 4) {
  ncom <- length(x$method)
  nran <- length(x$lower.random)
  nprd <- length(x$lower.predict)
  ##
  if (length(x$method.random) == 1 & nran > 1)
    x$method.random <- rep(x$method.random, nran)
  if (length(x$method.tau) == 1 & nran > 1)
    x$method.tau <- rep(x$method.tau, nran)
  if (length(x$method.tau.ci) == 1 & nran > 1)
    x$method.tau.ci <- rep(x$method.tau.ci, nran)
  if (length(x$tau) == 1 & nran > 1)
    x$tau <- rep(x$tau, nran)
  if (length(x$three.level) == 1 && x$three.level)
    x$tau <- sqrt(sum(x$tau^2))
  ##
  if (length(x$method.predict) == 1 && nprd > 1)
    x$method.predict <- rep(x$method.predict, nprd)
  ##
  res <-
    list(meth = 
           with(x,
                data.frame(
                  model = c(rep("common", length(method)),
                            rep("random", length(method.random))),
                  method = c(method, method.random),
                  three.level = three.level,
                  ##
                  k = k, k.all = k.all, k.MH = NA, k.study = k.study,
                  k.TE = k.TE, k0 = NA,
                  ##
                  method.tau = c(rep("", ncom), x$method.tau),
                  method.tau.ci = c(rep("", ncom), x$method.tau.ci),
                  tau = round(c(rep(NA, ncom), x$tau), digits),
                  tau.preset = NA,
                  ##
                  method.random.ci = c(rep("", ncom), method.random.ci),
                  df.random = c(rep(NA, ncom), unlist(df.random)),
                  adhoc.hakn.ci = c(rep("", ncom), adhoc.hakn.ci))),
         pred =
           with(x,
                data.frame(method.predict = method.predict))
         )
  ##
  if (!is.null(x$tau.preset))
    res$meth$tau.preset <- round(c(rep(NA, ncom), x$tau.preset), digits)
  ##
  if (inherits(x, "metabin")) {
    res$meth$k.MH <- x$k.MH
    res$meth$incr <- x$incr
    res$meth$method.incr <- x$method.incr
    res$meth$sparse <- x$sparse
    res$meth$allstudies <- x$allstudies
    res$meth$doublezeros <- x$doublezeros
    res$meth$MH.exact <- x$MH.exact
    res$meth$RR.Cochrane <- x$RR.Cochrane
    if (length(x$Q.Cochrane) == 1)
      x$Q.Cochrane <- rep(x$Q.Cochrane, nran)
    res$meth$Q.Cochrane <- c(rep(FALSE, ncom), x$Q.Cochrane)
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
    res$meth$incr <- x$incr
    res$meth$method.incr <- x$method.incr
    res$meth$sparse <- x$sparse
  }
  if (!is.null(x$model.glmm))
    res$meth$model.glmm <- x$model.glmm
  ##
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
