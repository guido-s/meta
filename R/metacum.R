#' Cumulative meta-analysis
#' 
#' @description
#' Performs a cumulative meta-analysis.
#' 
#' @param x An object of class \code{meta}.
#' @param pooled A character string indicating whether a common effect
#'   or random effects model is used for pooling. Either missing (see
#'   Details), \code{"common"}, or \code{"random"}, can be abbreviated.
#' @param sortvar An optional vector used to sort the individual
#'   studies (must be of same length as \code{x$TE}).
#' @param prediction A logical indicating whether to report prediction
#'   intervals.
#' @param no A numeric specifying which meta-analysis results to
#'   consider.
#' @param \dots Additional arguments (ignored).
#' 
#' @details
#' A cumulative meta-analysis is performed. Studies are included
#' sequentially as defined by \code{sortvar}.
#' 
#' Information from object \code{x} is utilised if argument
#' \code{pooled} is missing. A common effect model is assumed
#' (\code{pooled = "common"}) if argument \code{x$common} is
#' \code{TRUE}; a random effects model is assumed (\code{pooled =
#' "random"}) if argument \code{x$random} is \code{TRUE} and
#' \code{x$common} is \code{FALSE}.
#' 
#' @return
#' An object of class \code{"metacum"} and \code{"meta"} with
#' corresponding generic functions (see \code{\link{meta-object}}).
#' 
#' The following list elements have a different meaning:
#' \item{TE, seTE}{Estimated treatment effect and standard error of
#'   pooled estimate in cumulative meta-analyses.}
#' \item{lower, upper}{Lower and upper confidence interval limits.}
#' \item{statistic}{Statistic for test of overall effect.}
#' \item{pval}{P-value for test of overall effect.}
#' \item{studlab}{Study label describing addition of studies.}
#' \item{w}{Sum of weights from common effect or random effects model.}
#' \item{TE.common, seTE.common}{Value is \code{NA}.}
#' \item{TE.random, seTE.random}{Value is \code{NA}.}
#' \item{Q}{Value is \code{NA}.}
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}, \code{\link{metacont}},
#'   \code{\link{print.meta}}
#' 
#' @references
#' Cooper H & Hedges LV (1994):
#' \emph{The Handbook of Research Synthesis}.
#' Newbury Park, CA: Russell Sage Foundation
#'
#' @examples
#' data(Fleiss1993bin)
#' m1 <- metabin(d.asp, n.asp, d.plac, n.plac,
#'   data = Fleiss1993bin, studlab = study, sm = "RR", method = "I")
#' m1
#' metacum(m1)
#' metacum(m1, pooled = "random")
#' 
#' forest(metacum(m1))
#' forest(metacum(m1, pooled = "random"))
#' 
#' metacum(m1, sortvar = study)
#' metacum(m1, sortvar = 7:1)
#' 
#' m2 <- update(m1, title = "Fleiss1993bin meta-analysis", backtransf = FALSE)
#' metacum(m2)
#' 
#' data(Fleiss1993cont)
#' m3 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, sm = "SMD")
#' metacum(m3)
#' 
#' @rdname metacum
#' @method metacum meta
#' @export


metacum.meta <- function(x, pooled, sortvar, prediction, no = 1, ...) {
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(x, "meta")
  chksuitable(x, "Cumulative meta-analysis",
              c("trimfill", "metamerge", "netpairwise"))
  ##
  x <- updateversion(x)
  ##
  k.all <- length(x$TE)
  if (k.all < 2) {
    warning("Nothing calculated (minimum number of studies: 2).")
    return(invisible(NULL))
  }
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  if (!missing(pooled)) {
    pooled <- setchar(pooled, c("common", "random", "fixed"))
    pooled[pooled == "fixed"] <- "common"
  }
  else
    if (!x$common & x$random)
      pooled <- "random"
    else
      pooled <- "common"
  #
  if (missing(prediction))
    prediction <- pooled == "random" & x$prediction
  else
    chklogical(prediction)
  #
  mc <- match.call()
  error <-
    try(sortvar <-
          catch("sortvar", mc, x, sys.frame(sys.parent())),
        silent = TRUE)
  if (inherits(error, "try-error")) {
    sortvar <- catch("sortvar", mc, x$data,  NULL)
    if (isCol(x$data, ".subset"))
      sortvar <- sortvar[x$data$.subset]
  }
  sort <- !is.null(sortvar)
  if (sort && (length(sortvar) != k.all))
    stop("Number of studies in object 'x' and argument 'sortvar' ",
         "have different length.")
  if (!sort)
    sortvar <- seq_len(k.all)
  
  
  ##
  ##
  ## (3) Sort variables
  ##
  ##
  o <- order(sortvar)
  ##
  n.e <- x$n.e[o]
  n.c <- x$n.c[o]
  n   <- x$n[o]
  ##
  event.e <- x$event.e[o]
  event.c <- x$event.c[o]
  event   <- x$event[o]
  ##
  mean.e <- x$mean.e[o]
  mean.c <- x$mean.c[o]
  mean   <- x$mean[o]
  ##
  sd.e <- x$sd.e[o]
  sd.c <- x$sd.c[o]
  sd   <- x$sd[o]
  ##
  time.e <- x$time.e[o]
  time.c <- x$time.c[o]
  time   <- x$time[o]
  ##
  cor <- x$cor[o]
  ##
  TE <- x$TE[o]
  seTE <- x$seTE[o]
  #
  incr.e <- x$incr.e[o]
  incr.c <- x$incr.c[o]
  #
  if (length(x$incr) > 1)
    incr <- x$incr[o]
  else if (!is.null(x$incr))
    incr <- rep_len(x$incr, k.all)
  else
    incr <- x$incr
  ##
  ## Exclude studies from meta-analysis
  ##
  if (!is.null(x$exclude))
    exclude <- x$exclude[o]
  else
    exclude <- rep_len(FALSE, k.all)
  ##
  ncum <- cumsum(!exclude)
  ##
  studlab <- x$studlab[o]
  slab <- character(k.all)
  for (i in seq_len(k.all))
    slab[i] <- paste0("Adding ", studlab[i], " (k=", ncum[i], ")")
  ##
  chknumeric(no, min = 1, length = 1)
  ##
  ## Select a single common effect or random effects models
  ##
  if (pooled == "common") {
    if (no > length(x$seTE.common))
      stop("Argument 'no' must be smaller or equal to ",
           "number of common effect estimates.",
           call. = FALSE)
    ##
    no.c <- no
    no.r <- 1
  }
  else {
    if (no > length(x$seTE.random))
      stop("Argument 'no' must be smaller or equal to ",
           "number of random effects estimates.",
           call. = FALSE)
    ##
    no.c <- 1
    no.r <- no
  }
  ##
  x$TE.common <- x$TE.common[no.c]
  x$seTE.common <- x$seTE.common[no.c]
  x$statistic.common <- x$statistic.common[no.c]
  x$pval.common <- x$pval.common[no.c]
  x$lower.common <- x$lower.common[no.c]
  x$upper.common <- x$upper.common[no.c]
  x$zval.common <- x$zval.common[no.c]
  ##
  if (length(x$TE.random) == 1 &&
      length(x$TE.random) != length(x$seTE.random))
    x$TE.random <- rep_len(x$TE.random, length(x$seTE.random))
  x$TE.random <- x$TE.random[no.r]
  x$seTE.random <- x$seTE.random[no.r]
  x$df.random <- x$df.random[no.r]
  x$statistic.random <- x$statistic.random[no.r]
  x$pval.random <- x$pval.random[no.r]
  x$lower.random <- x$lower.random[no.r]
  x$upper.random <- x$upper.random[no.r]
  x$zval.random <- x$zval.random[no.r]
  x$seTE.hakn.adhoc.ci <- x$seTE.hakn.adhoc.ci[no.r]
  x$df.hakn.ci <- x$df.hakn.ci[no.r]
  ##
  x$text.random <- x$text.random[no.r]
  x$method.random.ci <- x$method.random.ci[no.r]
  x$adhoc.hakn.ci <- x$adhoc.hakn.ci[no.r]
  ##
  x$seTE.hakn.adhoc <- x$seTE.hakn.adhoc[no.r]
  x$df.hakn <- x$df.hakn[no.r]
  ##
  x$lower.predict <- x$lower.predict[1]
  x$upper.predict <- x$upper.predict[1]
  x$seTE.predict <- x$seTE.predict[1]
  x$df.predict <- x$df.predict[1]
  x$seTE.hakn.adhoc.pi <- x$seTE.hakn.adhoc.pi[1]
  x$df.hakn.pi <- x$df.hakn.pi[1]
  ##
  x$text.predict <- x$text.predict[1]
  x$method.predict <- x$method.predict[1]
  x$adhoc.hakn.pi <- x$adhoc.hakn.pi[1]
  
  
  ##
  ##
  ## (4) Do sensitivity analysis
  ##
  ##
  res.i <- matrix(NA, ncol = 28, nrow = k.all)
  add.i <- matrix(NA, ncol = 3, nrow = k.all)
  ##
  for (i in seq_len(k.all)) {
    sel <- 1:i
    ##
    if (length(incr) > 1)
      incr.i <- incr[sel]
    else
      incr.i <- incr
    ##
    if (inherits(x, "metabin"))
      m <- metabin(event.e[sel], n.e[sel], event.c[sel], n.c[sel],
                   ##
                   exclude = exclude[sel],
                   ##
                   method = x$method, sm = x$sm,
                   #
                   incr.e = incr.e[sel], incr.c = incr.c[sel],
                   allstudies = x$allstudies, MH.exact = x$MH.exact,
                   RR.Cochrane = x$RR.Cochrane, Q.Cochrane = x$Q.Cochrane,
                   model.glmm =
                     if (!is.null(x$model.glmm)) x$model.glmm else "UM.FS",
                   ##
                   method.tau = x$method.tau,
                   tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                   ##
                   level.ma = x$level.ma,
                   method.random.ci = x$method.random.ci,
                   adhoc.hakn.ci = x$adhoc.hakn.ci,
                   #
                   level.predict = x$level.predict,
                   method.predict = x$method.predict,
                   adhoc.hakn.pi = x$adhoc.hakn.pi,
                   #
                   keepdata = FALSE,
                   warn = FALSE,
                   ##
                   control = x$control)
    ##
    if (inherits(x, "metacont"))
      m <- metacont(n.e[sel], mean.e[sel], sd.e[sel],
                    n.c[sel], mean.c[sel], sd.c[sel],
                    ##
                    exclude = exclude[sel],
                    ##
                    sm = x$sm,
                    pooledvar = replaceNA(x$pooledvar, gs("pooledvar")),
                    ##
                    method.tau = x$method.tau,
                    tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                    ##
                    level.ma = x$level.ma,
                    method.random.ci = x$method.random.ci,
                    adhoc.hakn.ci = x$adhoc.hakn.ci,
                    #
                    level.predict = x$level.predict,
                    method.predict = x$method.predict,
                    adhoc.hakn.pi = x$adhoc.hakn.pi,
                    #
                    keepdata = FALSE,
                    warn = FALSE,
                    ##
                    control = x$control)
    ##
    if (inherits(x, "metacor"))
      m <- metacor(cor[sel], n[sel],
                   ##
                   exclude = exclude[sel],
                   ##
                   sm = x$sm, null.effect = x$null.effect,
                   ##
                   method.tau = x$method.tau,
                   tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                   ##
                   level.ma = x$level.ma,
                   method.random.ci = x$method.random.ci,
                   adhoc.hakn.ci = x$adhoc.hakn.ci,
                   #
                   level.predict = x$level.predict,
                   method.predict = x$method.predict,
                   adhoc.hakn.pi = x$adhoc.hakn.pi,
                   #
                   keepdata = FALSE,
                   ##
                   control = x$control)
    ##
    if (inherits(x, "metagen"))
      m <- metagen(TE[sel], seTE[sel],
                   ##
                   exclude = exclude[sel],
                   ##
                   sm = x$sm, null.effect = x$null.effect,
                   ##
                   method.tau = x$method.tau,
                   tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                   ##
                   level.ma = x$level.ma,
                   method.random.ci = x$method.random.ci,
                   adhoc.hakn.ci = x$adhoc.hakn.ci,
                   #
                   level.predict = x$level.predict,
                   method.predict = x$method.predict,
                   adhoc.hakn.pi = x$adhoc.hakn.pi,
                   #
                   keepdata = FALSE,
                   warn = FALSE,
                   ##
                   control = x$control)
    ##
    if (inherits(x,"metainc"))
      m <- metainc(event.e[sel], time.e[sel],
                   event.c[sel], time.c[sel],
                   ##
                   exclude = exclude[sel],
                   ##
                   method = x$method, sm = x$sm,
                   #
                   incr.e = incr.e[sel], incr.c = incr.c[sel],
                   model.glmm =
                     if (!is.null(x$model.glmm)) x$model.glmm else "UM.FS",
                   ##
                   method.tau = x$method.tau,
                   tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                   ##
                   level.ma = x$level.ma,
                   method.random.ci = x$method.random.ci,
                   adhoc.hakn.ci = x$adhoc.hakn.ci,
                   #
                   level.predict = x$level.predict,
                   method.predict = x$method.predict,
                   adhoc.hakn.pi = x$adhoc.hakn.pi,
                   #
                   keepdata = FALSE,
                   warn = FALSE,
                   ##
                   control = x$control)
    ##
    if (inherits(x, "metamean"))
      m <- metamean(n[sel], mean[sel], sd[sel],
                    ##
                    exclude = exclude[sel],
                    ##
                    sm = x$sm, null.effect = x$null.effect,
                    ##
                    method.tau = x$method.tau,
                    tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                    ##
                    level.ma = x$level.ma,
                    method.random.ci = x$method.random.ci,
                    adhoc.hakn.ci = x$adhoc.hakn.ci,
                    #
                    level.predict = x$level.predict,
                    method.predict = x$method.predict,
                    adhoc.hakn.pi = x$adhoc.hakn.pi,
                    #
                    keepdata = FALSE,
                    warn = FALSE,
                    ##
                    control = x$control)
    ##
    if (inherits(x, "metaprop"))
      m <- metaprop(event[sel], n[sel],
                    ##
                    exclude = exclude[sel],
                    ##
                    method = x$method, sm = x$sm, null.effect = x$null.effect,
                    ##
                    incr = incr.i, method.incr = x$method.incr,
                    method.ci = x$method.ci,
                    ##
                    method.tau = x$method.tau,
                    tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                    ##
                    level.ma = x$level.ma,
                    method.random.ci = x$method.random.ci,
                    adhoc.hakn.ci = x$adhoc.hakn.ci,
                    #
                    level.predict = x$level.predict,
                    method.predict = x$method.predict,
                    adhoc.hakn.pi = x$adhoc.hakn.pi,
                    #
                    keepdata = FALSE,
                    warn = FALSE,
                    ##
                    control = x$control)
    ##
    if (inherits(x, "metarate"))
      m <- metarate(event[sel], time[sel],
                    ##
                    exclude = exclude[sel],
                    ##
                    method = x$method, sm = x$sm, null.effect = x$null.effect,
                    ##
                    incr = incr.i, method.incr = x$method.incr,
                    ##
                    method.tau = x$method.tau,
                    tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                    ##
                    level.ma = x$level.ma,
                    method.random.ci = x$method.random.ci,
                    adhoc.hakn.ci = x$adhoc.hakn.ci,
                    #
                    level.predict = x$level.predict,
                    method.predict = x$method.predict,
                    adhoc.hakn.pi = x$adhoc.hakn.pi,
                    #
                    keepdata = FALSE,
                    warn = FALSE,
                    ##
                    control = x$control)
    ##
    sel.pft <- inherits(x, "metaprop") & x$sm == "PFT"
    sel.irft <- inherits(x, "metarate") & x$sm == "IRFT"
    ##
    add.i[i, ] <- c(m$method.tau.ci,  # 1
                    m$sign.lower.tau, # 2
                    m$sign.upper.tau  # 3
                    )
    ##
    if (pooled == "common") {
      res.i[i, ] <- c(m$TE.common,                                   #  1
                      m$seTE.common,                                 #  2
                      m$lower.common,                                #  3
                      m$upper.common,                                #  4
                      m$statistic.common,                            #  5
                      m$pval.common,                                 #  6
                      m$tau2,                                        #  7
                      m$lower.tau2,                                  #  8
                      m$upper.tau2,                                  #  9
                      m$se.tau2,                                     # 10
                      m$tau,                                         # 11
                      m$lower.tau,                                   # 12
                      m$upper.tau,                                   # 13
                      m$I2,                                          # 14
                      m$lower.I2,                                    # 15
                      m$upper.I2,                                    # 16
                      sum(m$w.common, na.rm = TRUE),                 # 17
                      if (sel.pft) 1 / mean(1 / n[sel]) else NA,     # 18
                      NA,                                            # 19
                      if (sel.irft) 1 / mean(1 / time[sel]) else NA, # 20
                      m$Rb,                                          # 21
                      NA,                                            # 22
                      NA,                                            # 23
                      m$k,                                           # 24
                      m$k.study,                                     # 25
                      m$k.all,                                       # 26
                      m$k.TE,                                        # 27
                      replaceNULL(m$k.MH)                            # 28
      )
    }
    ##
    else if (pooled == "random") {
      res.i[i, ] <- c(m$TE.random,                                   #  1
                      m$seTE.random,                                 #  2
                      m$lower.random,                                #  3
                      m$upper.random,                                #  4
                      m$statistic.random,                            #  5
                      m$pval.random,                                 #  6
                      m$tau2,                                        #  7
                      m$lower.tau2,                                  #  8
                      m$upper.tau2,                                  #  9
                      m$se.tau2,                                     # 10
                      m$tau,                                         # 11
                      m$lower.tau,                                   # 12
                      m$upper.tau,                                   # 13
                      m$I2,                                          # 14
                      m$lower.I2,                                    # 15
                      m$upper.I2,                                    # 16
                      sum(m$w.random, na.rm = TRUE),                 # 17
                      if (sel.pft) 1 / mean(1 / n[sel]) else NA,     # 18
                      if (x$method.random.ci %in% c("HK", "KR"))     #
                        m$df.random else NA,                         # 19
                      if (sel.irft) 1 / mean(1 / time[sel]) else NA, # 20
                      m$Rb,                                          # 21
                      m$lower.predict,                               # 22
                      m$upper.predict,                               # 23
                      m$k,                                           # 24
                      m$k.study,                                     # 25
                      m$k.all,                                       # 26
                      m$k.TE,                                        # 27
                      replaceNULL(m$k.MH)                            # 28
      )
    }
  }
  ##
  TE.i <- res.i[, 1]
  seTE.i <- res.i[, 2]
  lower.i <- res.i[, 3]
  upper.i <- res.i[, 4]
  statistic.i <- res.i[, 5]
  pval.i <- res.i[, 6]
  ##
  tau2.i <- res.i[, 7]
  lower.tau2.i <- res.i[, 8]
  upper.tau2.i <- res.i[, 9]
  se.tau2.i <- res.i[, 10]
  ##
  tau.i <- res.i[, 11]
  lower.tau.i <- res.i[, 12]
  upper.tau.i <- res.i[, 13]
  ##
  I2.i <- res.i[, 14]
  lower.I2.i <- res.i[, 15]
  upper.I2.i <- res.i[, 16]
  ##
  weight.i <- res.i[, 17]
  n.harmonic.mean.i <- res.i[, 18]
  if (pooled == "random" & x$method.random.ci %in% c("HK", "KR"))
    df.random.i <- res.i[, 19]
  t.harmonic.mean.i <- res.i[, 20]
  Rb.i <- res.i[, 21]
  #
  lower.predict.i <- res.i[, 22]
  upper.predict.i <- res.i[, 23]
  #
  k.i <- res.i[, 24]
  k.study.i <- res.i[, 25]
  k.all.i <- res.i[, 26]
  k.TE.i <- res.i[, 27]
  k.MH.i <- res.i[, 28]
  #
  method.tau.ci <- replaceNULL(unique(add.i[add.i[, 1] != "", 1]), "")
  sign.lower.tau.i <- replaceNULL(unique(add.i[add.i[, 2] != "", 2]), "")
  sign.upper.tau.i <- replaceNULL(unique(add.i[add.i[, 3] != "", 3]), "")
  #
  if (pooled == "common") {
    TE.pooled <- x$TE.common
    seTE.pooled <- x$seTE.common
    lower.pooled <- x$lower.common
    upper.pooled <- x$upper.common
    statistic.pooled <- x$statistic.common
    pval.pooled <- x$pval.common
    w.pooled <- sum(x$w.common, na.rm = TRUE)
    #
    lower.predict.pooled <- NA
    upper.predict.pooled <- NA
  }
  #
  else if (pooled == "random") {
    TE.pooled <- x$TE.random
    seTE.pooled <- x$seTE.random
    lower.pooled <- x$lower.random
    upper.pooled <- x$upper.random
    statistic.pooled <- x$statistic.random
    pval.pooled <- x$pval.random
    w.pooled <- sum(x$w.random, na.rm = TRUE)
    #
    lower.predict.pooled <- x$lower.predict
    upper.predict.pooled <- x$upper.predict
  }
  #
  df.random.pooled <- x$df.random
  #
  tau2.pooled <- x$tau2
  se.tau2.pooled <- x$se.tau2
  lower.tau2.pooled <- x$lower.tau2
  upper.tau2.pooled <- x$upper.tau2
  #
  tau.pooled <- x$tau
  lower.tau.pooled <- x$lower.tau
  upper.tau.pooled <- x$upper.tau
  #
  I2.pooled <- x$I2
  lower.I2.pooled <- x$lower.I2
  upper.I2.pooled <- x$upper.I2
  #
  Rb.pooled <- x$Rb
  #
  n.harmonic.mean.pooled <- 1 / mean(1 / n)
  t.harmonic.mean.pooled <- 1 / mean(1 / time)
  
  
  ##
  ##
  ## (5) Generate R object
  ##
  ##
  res <- list(studlab = slab,
              ##
              sm = x$sm,
              null.effect = x$null.effect,
              ##
              TE = TE.i,
              seTE = seTE.i,
              statistic = statistic.i,
              pval = pval.i,
              level = x$level.ma,
              lower = lower.i,
              upper = upper.i,
              lower.predict = lower.predict.i,
              upper.predict = upper.predict.i,
              w = weight.i,
              #
              tau2 = tau2.i,
              se.tau2 = se.tau2.i,
              lower.tau2 = lower.tau2.i,
              upper.tau2 = upper.tau2.i,
              #
              tau = tau.i,
              lower.tau = lower.tau.i,
              upper.tau = upper.tau.i,
              #
              sign.lower.tau.i = sign.lower.tau.i,
              sign.upper.tau.i = sign.upper.tau.i,
              #
              I2 = I2.i,
              lower.I2 = lower.I2.i,
              upper.I2 = upper.I2.i,
              #
              Rb = Rb.i,
              #
              n.harmonic.mean = n.harmonic.mean.i,
              t.harmonic.mean = t.harmonic.mean.i,
              #
              k = k.i,
              k.study = k.study.i,
              k.all = k.all.i,
              k.TE = k.TE.i,
              k.MH = k.MH.i,
              #
              three.level = x$three.level,
              cluster = x$cluster,
              rho = x$rho,
              #
              pooled = pooled,
              common = pooled == "common",
              random = pooled == "random",
              overall = TRUE,
              overall.hetstat = TRUE,
              prediction = prediction,
              method.predict = x$method.predict,
              adhoc.hakn.pi = x$adhoc.hakn.pi,
              df.predict = x$df.predict,
              backtransf = x$backtransf,
              func.backtransf = x$func.backtransf,
              #
              level.ma = x$level.ma,
              level.predict = x$level.predict,
              #
              method = x$method,
              method.random = x$method.random,
              #
              method.random.ci = x$method.random.ci,
              adhoc.hakn.ci = x$adhoc.hakn.ci,
              #
              method.tau = x$method.tau,
              method.tau.ci = method.tau.ci,
              #
              tau.preset = x$tau.preset,
              TE.tau = x$TE.tau,
              #
              method.I2 = x$method.I2,
              #
              k.pooled = x$k,
              k.study.pooled = x$k.study,
              k.all.pooled = x$k.all,
              k.TE.pooled = x$k.TE,
              k.MH.pooled = x$k.MH,
              #
              df.random =
                if (pooled == "random" &
                    x$method.random.ci %in% c("HK", "KR"))
                  df.random.i else NULL,
              #
              TE.pooled = TE.pooled,
              seTE.pooled = seTE.pooled,
              lower.pooled = lower.pooled,
              upper.pooled = upper.pooled,
              statistic.pooled = statistic.pooled,
              pval.pooled = pval.pooled,
              w.pooled = w.pooled,
              text.pooled = "Pooled estimate",
              #
              lower.predict.pooled = lower.predict.pooled,
              upper.predict.pooled = upper.predict.pooled,
              #
              df.random.pooled = df.random.pooled,
              #
              Q.pooled = x$Q,
              #
              tau2.pooled = tau2.pooled,
              se.tau2.pooled = se.tau2.pooled,
              lower.tau2.pooled = lower.tau2.pooled,
              upper.tau2.pooled = upper.tau2.pooled,
              #
              tau.pooled = tau.pooled,
              lower.tau.pooled = lower.tau.pooled,
              upper.tau.pooled = upper.tau.pooled,
              #
              I2.pooled = I2.pooled,
              lower.I2.pooled = lower.I2.pooled,
              upper.I2.pooled = upper.I2.pooled,
              #
              Rb.pooled = Rb.pooled,
              #
              n.harmonic.mean.pooled = n.harmonic.mean.pooled,
              t.harmonic.mean.pooled = t.harmonic.mean.pooled,
              #
              pscale = x$pscale,
              irscale = x$irscale, irunit = x$irunit,
              #
              label.e = x$label.e,
              label.c = x$label.c,
              #
              text.common = x$text.common, text.random = x$text.random,
              text.predict = x$text.predict,
              text.w.common = x$text.w.common, text.w.random = x$text.w.random,
              ##
              title = x$title, complab = x$complab,
              outclab = x$outclab,
              ##
              no = no,
              ##
              x = x,
              ##
              call = match.call())
  
  res$version <- packageDescription("meta")$Version
  ##
  res$x$common <- res$common
  res$x$random <- res$random
  #
  class(res) <- "metacum"
  #
  res
}





#' @rdname metacum
#' @export metacum


metacum <- function(x, ...) 
  UseMethod("metacum")





#' @rdname metacum
#' @method metacum default
#' @export


metacum.default <- function(x, ...)
  stop("Cumulative meta-analysis not available for an object of class '",
       class(x)[1], "'.",
       call. = FALSE)
