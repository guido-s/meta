trimfill.meta <- function(x, left = NULL, ma.fixed = TRUE,
                          type = "L", n.iter.max = 50,
                          level = x$level, level.comb = x$level.comb,
                          comb.fixed = FALSE, comb.random = TRUE,
                          hakn = x$hakn,
                          method.tau = x$method.tau,
                          prediction = x$prediction, level.predict = x$level.predict,
                          backtransf = x$backtransf, pscale = x$pscale,
                          irscale = x$irscale, irunit = x$irunit,
                          silent = TRUE, ...) {
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(x, "meta")
  if (inherits(x, "metacum"))
    stop("This function is not usable for an object of class \"metacum\"")
  if (inherits(x, "metainf"))
    stop("This function is not usable for an object of class \"metainf\"")
  x <- updateversion(x)
  
  
  ##
  ## Check arguments
  ##
  type <- setchar(type, c("L", "R"))
  ##
  chklevel(level)
  chklevel(level.comb)
  chklevel(level.predict)
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  ##
  chklogical(prediction)
  ##
  chklogical(backtransf)
  sm <- x$sm
  if (!is.prop(sm))
    pscale <- 1
  chknumeric(pscale, single = TRUE)
  if (!backtransf & pscale != 1) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  if (!is.rate(sm))
    irscale <- 1
  chknumeric(irscale, single = TRUE)
  if (!backtransf & irscale != 1) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  ##
  chklogical(silent)
  
  
  TE <- x$TE
  seTE <- x$seTE
  studlab <- x$studlab
  ##
  n.e <- x$n.e
  n.c <- x$n.c
  n <- x$n
  ##
  event.e <- x$event.e
  event.c <- x$event.c
  event <- x$event
  ##
  time.e <- x$time.e
  time.c <- x$time.c
  time <- x$time
  ##
  cor <- x$cor
  ##
  mean.e <- x$mean.e
  mean.c <- x$mean.c
  ##
  sd.e <- x$sd.e
  sd.c <- x$sd.c
  ##
  transf.null.effect <- null.effect <- x$null.effect
  ##
  if (sm %in% c("PFT", "PAS"))
    transf.null.effect <- asin(sqrt(null.effect))
  else if (is.log.effect(sm))
    transf.null.effect <- log(null.effect)
  else if (sm == c("PLOGIT"))
    transf.null.effect <- log(null.effect / (1 - null.effect))
  else if (sm %in% c("IRS", "IRFT"))
    transf.null.effect <- sqrt(null.effect)
  else if (sm == "ZCOR")
    transf.null.effect <- 0.5 * log((1 + null.effect) / (1 - null.effect))
  
  
  if(length(TE) != length(seTE))
    stop("length of argument TE and seTE must be equal")
  ##
  if(length(TE) != length(studlab))
    stop("length of argument TE and studlab must be equal")
  ##
  ## Exclude studies from meta-analysis
  ##
  if (!is.null(x$exclude)) {
    exclude <- x$exclude
    nomiss <- !is.na(TE) & !is.na(seTE)
    miss <- !nomiss & !exclude
    ##
    sel <- nomiss & !exclude
  }
  else {
    exclude <- exclude.na <- NULL
    nomiss <- !is.na(TE) & !is.na(seTE)
    miss <- !nomiss
    ##
    sel <- nomiss
  }
  ##
  if (any(miss))
    warning(paste(sum(miss),
                  "observation(s) dropped due to missing values"))
  
  
  TE <- TE[sel]
  seTE <- seTE[sel]
  studlab <- studlab[sel]
  ##
  if (!is.null(n.e))
    n.e <- n.e[sel]
  if (!is.null(n.c))
    n.c <- n.c[sel]
  if (!is.null(n))
    n <- n[sel]
  ##
  if (!is.null(event.e))
    event.e <- event.e[sel]
  if (!is.null(event.c))
    event.c <- event.c[sel]
  if (!is.null(event))
    event <- event[sel]
  ##
  if (!is.null(time.e))
    time.e <- time.e[sel]
  if (!is.null(time.c))
    time.c <- time.c[sel]
  if (!is.null(time))
    time <- time[sel]
  ##
  if (!is.null(cor))
    cor <- cor[sel]
  ##
  if (!is.null(mean.e))
    mean.e <- mean.e[sel]
  if (!is.null(mean.c))
    mean.c <- mean.c[sel]
  ##
  if (!is.null(sd.e))
    sd.e <- sd.e[sel]
  if (!is.null(sd.c))
    sd.c <- sd.c[sel]
  ##
  k <- length(TE)
  ##
  if (k <= 2) {
    warning("Minimal number of three studies for trim-and-fill method")
    return(invisible(NULL))
  }
  
  
  if (is.null(left))
    left <- as.logical(sign(metabias(TE, seTE, method = "linreg", k.min = 3)$estimate[1]) == 1)
  ##
  if (!left) TE <- -TE
  ##
  ord <- order(TE)
  TE <- TE[ord]
  seTE <- seTE[ord]
  studlab <- studlab[ord]
  ##
  if (!is.null(n.e))
    n.e <- n.e[ord]
  if (!is.null(n.c))
    n.c <- n.c[ord]
  if (!is.null(n))
    n <- n[ord]
  ##
  if (!is.null(event.e))
    event.e <- event.e[ord]
  if (!is.null(event.c))
    event.c <- event.c[ord]
  if (!is.null(event))
    event <- event[ord]
  ##
  if (!is.null(time.e))
    time.e <- time.e[ord]
  if (!is.null(time.c))
    time.c <- time.c[ord]
  if (!is.null(time))
    time <- time[ord]
  ##
  if (!is.null(cor))
    cor <- cor[ord]
  ##
  if (!is.null(mean.e))
    mean.e <- mean.e[ord]
  if (!is.null(mean.c))
    mean.c <- mean.c[ord]
  ##
  if (!is.null(sd.e))
    sd.e <- sd.e[ord]
  if (!is.null(sd.c))
    sd.c <- sd.c[ord]
  
  
  if (ma.fixed)
    TE.sum <- metagen(TE, seTE)$TE.fixed
  else
    TE.sum <- metagen(TE, seTE, method.tau = method.tau)$TE.random
  
  
  if (k == 1) {
    n.iter <- 0
    k0 <- -9
  }
  else {
    n.iter  <-  0
    k0.last <- -1
    k0      <-  0
    ##
    while (k0.last != k0 & k0 <= (k - 1) & n.iter < n.iter.max) {
      ##
      n.iter <- n.iter + 1
      ##
      k0.last <- k0
      ##
      sel <- 1:(k - k0)
      ##
      if (ma.fixed)
        TE.sum <- metagen(TE[sel], seTE[sel])$TE.fixed
      else
        TE.sum <- metagen(TE[sel], seTE[sel],
                          method.tau = method.tau)$TE.random
      ##
      trim1 <- estimate.missing(TE, TE.sum, type)
      ##
      if (!silent) {
        cat("n.iter = ", n.iter, "\n", sep = "")
        if (type == "L")
          cat("L0 = ", round(trim1$res0, 2), "\n\n", sep = "")
        if (type == "R")
          cat("R0 = ", round(trim1$res0 + 0.5, 2), "\n\n", sep = "")
      }
      ##
      k0 <- trim1$res0.plus
    }
  }
  
  
  if (k0 > (k - 1)) k0 <- k - 1
  ##
  if (k0 > 0) {
    TE.star   <- 2 * TE.sum - TE[(k - k0 + 1):k]
    seTE.star <- seTE[(k - k0 + 1):k]
    ##
    trimfill  <- c(rep(FALSE, length(TE)),
                   rep(TRUE, length(TE.star)))
    ##
    TE      <- c(TE[order(ord)], TE.star)
    seTE    <- c(seTE[order(ord)], seTE.star)
    studlab <- c(studlab[order(ord)],
                 paste("Filled:", studlab[(k - k0 + 1):k]))
    ##
    if (!is.null(n.e))
      n.e <- c(n.e[order(ord)], n.e[(k - k0 + 1):k])
    if (!is.null(n.c))
      n.c <- c(n.c[order(ord)], n.c[(k - k0 + 1):k])
    if (!is.null(n))
      n <- c(n[order(ord)], n[(k - k0 + 1):k])
    ##
    if (!is.null(event.e))
      event.e <- c(event.e[order(ord)], event.e[(k - k0 + 1):k])
    if (!is.null(event.c))
      event.c <- c(event.c[order(ord)], event.c[(k - k0 + 1):k])
    if (!is.null(event))
      event <- c(event[order(ord)], event[(k - k0 + 1):k])
    ##
    if (!is.null(time.e))
      time.e <- c(time.e[order(ord)], time.e[(k - k0 + 1):k])
    if (!is.null(time.c))
      time.c <- c(time.c[order(ord)], time.c[(k - k0 + 1):k])
    if (!is.null(time))
      time <- c(time[order(ord)], time[(k - k0 + 1):k])
    ##
    if (!is.null(cor))
      cor <- c(cor[order(ord)], cor[(k - k0 + 1):k])
    ##
    if (!is.null(mean.e))
      mean.e <- c(mean.e[order(ord)], mean.e[(k - k0 + 1):k])
    if (!is.null(mean.c))
      mean.c <- c(mean.c[order(ord)], mean.c[(k - k0 + 1):k])
    ##
    if (!is.null(sd.e))
      sd.e <- c(sd.e[order(ord)], sd.e[(k - k0 + 1):k])
    if (!is.null(sd.c))
      sd.c <- c(sd.c[order(ord)], sd.c[(k - k0 + 1):k])
  }
  else {
    TE.star   <- NA
    seTE.star <- NA
    trimfill  <- rep(FALSE, length(TE))
    TE        <- TE[order(ord)]
    seTE      <- seTE[order(ord)]
    studlab   <- studlab[order(ord)]
    ##
    if (!is.null(n.e))
      n.e <- n.e[order(ord)]
    if (!is.null(n.c))
      n.c <- n.c[order(ord)]
    if (!is.null(n))
      n <- n[order(ord)]
    ##
    if (!is.null(event.e))
      event.e <- event.e[order(ord)]
    if (!is.null(event.c))
      event.c <- event.c[order(ord)]
    if (!is.null(event))
      event <- event[order(ord)]
    ##
    if (!is.null(time.e))
      time.e <- time.e[order(ord)]
    if (!is.null(time.c))
      time.c <- time.c[order(ord)]
    if (!is.null(time))
      time <- time[order(ord)]
    ##
    if (!is.null(cor))
      cor <- cor[order(ord)]
    ##
    if (!is.null(mean.e))
      mean.e <- mean.e[order(ord)]
    if (!is.null(mean.c))
      mean.c <- mean.c[order(ord)]
    ##
    if (!is.null(sd.e))
      sd.e <- sd.e[order(ord)]
    if (!is.null(sd.c))
      sd.c <- sd.c[order(ord)]
  }
  
  
  if (!left)
    m <- metagen(-TE, seTE, studlab = studlab,
                 level = level, level.comb = level.comb,
                 hakn = hakn, method.tau = method.tau,
                 prediction = prediction, level.predict = level.predict,
                 null.effect = transf.null.effect)
  else
    m <- metagen(TE, seTE, studlab = studlab,
                 level = level, level.comb = level.comb,
                 hakn = hakn, method.tau = method.tau,
                 prediction = prediction, level.predict = level.predict,
                 null.effect = transf.null.effect)
  
  
  ##
  ## Calculate H, I-Squared, and Rb
  ##
  Hres  <- calcH(m$Q, m$df.Q, level.comb)
  I2res <- isquared(m$Q, m$df.Q, level.comb)
  Rbres <- with(m,
                Rb(seTE[!is.na(seTE)], seTE.random, tau^2, Q, df.Q, level.comb))
  
  
  ##
  ## Number of filled studies
  ##
  k0 <- sum(trimfill)
  
  
  if (!is.null(exclude) && any(exclude)) {
    exclude.na <- c(exclude, rep(NA, k0))
    exclude <- c(exclude, rep(FALSE, k0))
    TE.all <- seTE.all <- studlab.all <- rep(NA, length(exclude))
    ##
    TE.all[exclude] <- x$TE[exclude]
    TE.all[!exclude] <- TE
    ##
    seTE.all[exclude] <- x$seTE[exclude]
    seTE.all[!exclude] <- seTE
    ##
    studlab.all[exclude] <- x$studlab[exclude]
    studlab.all[!exclude] <- studlab
    ##
    if (!left)
      m.all <- metagen(-TE.all, seTE.all, studlab = studlab.all,
                       exclude = exclude, level = level,
                       null.effect = transf.null.effect)
    else
      m.all <- metagen(TE.all, seTE.all, studlab = studlab.all,
                       exclude = exclude, level = level,
                       null.effect = transf.null.effect)
  }
  else
    m.all <- m
  
  
  res <- list(studlab = m.all$studlab,
              TE = m.all$TE, seTE = m.all$seTE,
              lower = m.all$lower, upper = m.all$upper,
              zval = m.all$zval, pval = m.all$pval,
              w.fixed = m.all$w.fixed, w.random = m.all$w.random,
              exclude = exclude.na,
              ##
              TE.fixed = m$TE.fixed, seTE.fixed = m$seTE.fixed,
              lower.fixed = m$lower.fixed, upper.fixed = m$upper.fixed,
              zval.fixed = m$zval.fixed, pval.fixed = m$pval.fixed,
              ##
              TE.random = m$TE.random, seTE.random = m$seTE.random,
              lower.random = m$lower.random, upper.random = m$upper.random,
              zval.random = m$zval.random, pval.random = m$pval.random,
              ##
              seTE.predict = m$seTE.predict,
              lower.predict = m$lower.predict,
              upper.predict = m$upper.predict,
              level.predict = level.predict,
              ##
              k = m$k, Q = m$Q, df.Q = m$df.Q, tau = m$tau,
              ##
              H = Hres$TE,
              lower.H = Hres$lower,
              upper.H = Hres$upper,
              ##
              I2 = I2res$TE,
              lower.I2 = I2res$lower,
              upper.I2 = I2res$upper,
              ##
              Rb = Rbres$TE,
              lower.Rb = Rbres$lower,
              upper.Rb = Rbres$upper,
              ##
              sm = sm,
              method = m$method,
              ##
              call = match.call(),
              left = left,
              ma.fixed = ma.fixed,
              type = type,
              n.iter.max = n.iter.max,
              n.iter = n.iter,
              trimfill = trimfill,
              hakn = m$hakn,
              df.hakn = m$df.hakn,
              method.tau = m$method.tau,
              prediction = prediction,
              title = x$title,
              complab = x$complab,
              outclab = x$outclab,
              label.e = x$label.e,
              label.c = x$label.c,
              label.left = x$label.left,
              label.right = x$label.right,
              k0 = k0,
              level = level, level.comb = level.comb,
              comb.fixed = comb.fixed,
              comb.random = comb.random,
              ##
              n.e = n.e,
              n.c = n.c,
              n = n,
              ##
              event.e = event.e,
              event.c = event.c,
              event = event,
              ##
              time.e = time.e,
              time.c = time.c,
              time = time,
              ##
              cor = cor,
              ##
              mean.e = mean.e,
              mean.c = mean.c,
              ##
              sd.e = sd.e,
              sd.c = sd.c,
              ##
              null.effect = x$null.effect,
              ##
              class.x = class(x)[1]
              )
  
  res$backtransf <- backtransf
  res$pscale <- pscale
  res$irscale <- irscale
  res$irunit <- irunit
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metagen", "meta", "trimfill")
  ##
  res
}
