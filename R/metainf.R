#' Influence analysis in meta-analysis using leave-one-out method
#' 
#' @description
#' Performs an influence analysis. Pooled estimates are calculated
#' omitting one study at a time.
#' 
#' @param x An object of class \code{meta}.
#' @param pooled A character string indicating whether a common effect
#'   or random effects model is used for pooling. Either missing (see
#'   Details), \code{"common"} or \code{"random"}, can be abbreviated.
#' @param sortvar An optional vector used to sort the individual
#'   studies (must be of same length as \code{x$TE}).
#' @param no A numeric specifying which meta-analysis results to
#'   consider.
#' 
#' @details
#' Performs a influence analysis; pooled estimates are calculated
#' omitting one study at a time. Studies are sorted according to
#' \code{sortvar}.
#' 
#' Information from object \code{x} is utilised if argument
#' \code{pooled} is missing. A common effect model is assumed
#' (\code{pooled="common"}) if argument \code{x$common} is
#' \code{TRUE}; a random effects model is assumed
#' (\code{pooled="random"}) if argument \code{x$random} is
#' \code{TRUE} and \code{x$common} is \code{FALSE}.
#' 
#' @return
#' An object of class \code{"meta"} and \code{"metainf"} with
#' corresponding generic functions (see \code{\link{meta-object}}).
#' 
#' The following list elements have a different meaning:
#' \item{TE, seTE}{Estimated treatment effect and standard error of
#'   pooled estimate in influence analysis.}
#' \item{lower, upper}{Lower and upper confidence interval limits.}
#' \item{statistic}{Statistic for test of overall effect.}
#' \item{pval}{P-value for test of overall effect.}
#' \item{studlab}{Study label describing omission of studies.}
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
#' metainf(m1)
#' metainf(m1, pooled = "random")
#' 
#' forest(metainf(m1))
#' forest(metainf(m1), layout = "revman5")
#' forest(metainf(m1, pooled = "random"))
#' 
#' metainf(m1, sortvar = study)
#' metainf(m1, sortvar = 7:1)
#' 
#' m2 <- update(m1, title = "Fleiss1993bin meta-analysis", backtransf = FALSE)
#' metainf(m2)
#' 
#' data(Fleiss1993cont)
#' m3 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, sm = "SMD")
#' metainf(m3)
#' 
#' @export metainf


metainf <- function(x, pooled, sortvar, no = 1) {
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(x, "meta")
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
  ##
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
    sortvar <- 1:k.all
  
  
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
  ##
  if (length(x$incr) > 1)
    incr <- x$incr[o]
  else if (!is.null(x$incr))
    incr <- rep_len(x$incr, k.all)
  else
    incr <- x$incr
  ##
  studlab <- x$studlab[o]
  slab <- c(paste("Omitting", studlab), "Pooled estimate")
  studlab <- c(rev(rev(slab)[-1]), " ", rev(slab)[1])
  ##
  ## Exclude studies from meta-analysis
  ##
  if (!is.null(x$exclude))
    exclude <- x$exclude[o]
  else
    exclude <- rep_len(FALSE, k.all)
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
  res.i <- matrix(NA, ncol = 21, nrow = k.all)
  add.i <- matrix(NA, ncol = 3, nrow = k.all)
  ##
  for (i in 1:k.all) {
    sel <- -i
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
                   ##
                   incr = incr.i, method.incr = x$method.incr,
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
                   ##
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
                    sm = x$sm, pooledvar = x$pooledvar,
                    ##
                    method.tau = x$method.tau,
                    tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                    ##
                    level.ma = x$level.ma,
                    method.random.ci = x$method.random.ci,
                    adhoc.hakn.ci = x$adhoc.hakn.ci,
                    ##
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
                   ##
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
                   ##
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
                   incr = incr.i, method.incr = x$method.incr,
                   model.glmm =
                     if (!is.null(x$model.glmm)) x$model.glmm else "UM.FS",
                   ##
                   method.tau = x$method.tau,
                   tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                   ##
                   level.ma = x$level.ma,
                   method.random.ci = x$method.random.ci,
                   adhoc.hakn.ci = x$adhoc.hakn.ci,
                   ##
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
                    ##
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
                    ##
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
                    ##
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
                      m$Rb                                           # 21
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
                      m$Rb                                           # 21
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
  ##
  method.tau.ci <- unique(add.i[, 1])
  sign.lower.tau.i <- add.i[, 2]
  sign.upper.tau.i <- add.i[, 3]
  ##  
  if (pooled == "common") {
    TE.s <- x$TE.common
    seTE.s <- x$seTE.common
    lower.TE.s <- x$lower.common
    upper.TE.s <- x$upper.common
    statistic.s <- x$statistic.common
    pval.s <- x$pval.common
    w.s <- sum(x$w.common, na.rm = TRUE)
  }
  ##
  else if (pooled == "random") {
    TE.s <- x$TE.random
    seTE.s <- x$seTE.random
    lower.TE.s <- x$lower.random
    upper.TE.s <- x$upper.random
    statistic.s <- x$statistic.random
    pval.s <- x$pval.random
    w.s <- sum(x$w.random, na.rm = TRUE)
  }
  
  
  ##
  ##
  ## (5) Generate R object
  ##
  ##
  res <- list(studlab = studlab,
              ##
              sm = x$sm,
              null.effect = x$null.effect,
              ##
              TE = c(TE.i, NA, TE.s),
              seTE = c(seTE.i, NA, seTE.s),
              statistic = c(statistic.i, NA, statistic.s),
              pval = c(pval.i, NA, pval.s),
              level = x$level.ma,
              lower = c(lower.i, NA, lower.TE.s),
              upper = c(upper.i, NA, upper.TE.s),
              ##
              three.level = x$three.level,
              cluster = x$cluster,
              ##
              k = x$k, k.study = x$k.study, k.all = x$k.all, k.TE = x$k.TE,
              ##
              pooled = pooled,
              common = ifelse(pooled == "common", TRUE, FALSE),
              random = ifelse(pooled == "random", TRUE, FALSE),
              overall = TRUE,
              overall.hetstat = TRUE,
              prediction = FALSE,
              backtransf = x$backtransf,
              ##
              method = x$method,
              ##
              w = c(weight.i, NA, w.s),
              TE.common = NA, seTE.common = NA,
              TE.random = NA, seTE.random = NA,
              df.random =
                if (pooled == "random" &
                    x$method.random.ci %in% c("HK", "KR"))
                  c(df.random.i, NA, x$df.random) else NULL,
              level.ma = x$level.ma,
              method.random.ci = x$method.random.ci,
              adhoc.hakn.ci = x$adhoc.hakn.ci,
              ##
              Q = NA,
              ##
              method.tau = x$method.tau,
              method.tau.ci = method.tau.ci,
              tau2 = c(tau2.i, NA, x$tau2),
              se.tau2 = c(se.tau2.i, NA, x$se.tau2),
              lower.tau2 = c(lower.tau2.i, NA, x$lower.tau2),
              upper.tau2 = c(upper.tau2.i, NA, x$upper.tau2),
              tau = c(tau.i, NA, x$tau),
              lower.tau = c(lower.tau.i, NA, x$lower.tau),
              upper.tau = c(upper.tau.i, NA, x$upper.tau),
              tau.preset = x$tau.preset,
              TE.tau = x$TE.tau,
              sign.lower.tau.i = c(sign.lower.tau.i, NA, x$sign.lower.tau),
              sign.upper.tau.i = c(sign.upper.tau.i, NA, x$sign.upper.tau),
              ##
              I2 = c(I2.i, NA, x$I2),
              lower.I2 = c(lower.I2.i, NA, x$lower.I2),
              upper.I2 = c(upper.I2.i, NA, x$upper.I2),
              ##
              Rb = c(Rb.i, NA, x$Rb),
              ##
              n.harmonic.mean = c(n.harmonic.mean.i, NA, 1 / mean(1 / n)),
              t.harmonic.mean = c(t.harmonic.mean.i, NA, 1 / mean(1 / time)),
              ##
              pscale = x$pscale,
              irscale = x$irscale, irunit = x$irunit,
              ##
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
  ##
  ## Backward compatibility
  ##
  res$fixed <- res$x$fixed <- res$common
  res$TE.fixed <- res$TE.common
  res$seTE.fixed <- res$seTE.common
  ##
  res$text.fixed <- res$text.common
  res$text.w.fixed <- res$text.w.common
  
  
  class(res) <- c("metainf", "summary.meta", "meta")
  ##
  if (inherits(x, "trimfill"))
    class(res) <- c(class(res), "trimfill")
  
  res
}
