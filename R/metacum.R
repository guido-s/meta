#' Cumulative meta-analysis
#' 
#' @description
#' Performs a cumulative meta-analysis.
#' 
#' @param x An object of class \code{meta}.
#' @param pooled A character string indicating whether a fixed effect
#'   or random effects model is used for pooling. Either missing (see
#'   Details), \code{"fixed"}, or \code{"random"}, can be abbreviated.
#' @param sortvar An optional vector used to sort the individual
#'   studies (must be of same length as \code{x$TE}).
#' 
#' @details
#' A cumulative meta-analysis is performed. Studies are included
#' sequentially as defined by \code{sortvar}.
#' 
#' Information from object \code{x} is utilised if argument
#' \code{pooled} is missing. A fixed effect model is assumed
#' (\code{pooled = "fixed"}) if argument \code{x$comb.fixed} is
#' \code{TRUE}; a random effects model is assumed (\code{pooled =
#' "random"}) if argument \code{x$comb.random} is \code{TRUE} and
#' \code{x$comb.fixed} is \code{FALSE}.
#' 
#' @return
#' An object of class \code{c("metacum", "meta")} with corresponding
#' \code{print}, and \code{forest} functions. The object is a list
#' containing the following components:
#' \item{TE, seTE}{Estimated treatment effect and standard error of
#'   pooled estimate in cumulative meta-analyses.}
#' \item{lower, upper}{Lower and upper confidence interval limits.}
#' \item{zval}{z-value for test of overall effect.}
#' \item{pval}{P-value for test of overall effect.}
#' \item{studlab}{Study label describing addition of studies.}
#' \item{w}{Sum of weights from fixed effect or random effects model.}
#' \item{I2}{Heterogeneity statistic I\eqn{^2}.}
#' \item{Rb}{Heterogeneity statistic R\eqn{_b}.}
#' \item{tau}{Square-root of between-study variance.}
#' \item{df.hakn}{Degrees of freedom for test of treatment effect for
#'   Hartung-Knapp method (only if \code{hakn = TRUE}).}
#' \item{sm}{Summary measure.}
#' \item{method}{Method used for pooling.}
#' \item{k}{Number of studies combined in meta-analysis.}
#' \item{pooled}{As defined above.}
#' \item{comb.fixed}{A logical indicating whether analysis is based on
#'   fixed effect model.}
#' \item{comb.random}{A logical indicating whether analysis is based
#'   on random effects model.}
#' \item{TE.fixed, seTE.fixed}{Value is \code{NA}.}
#' \item{TE.random, seTE.random}{Value is \code{NA}.}
#' \item{Q}{Value is \code{NA}.}
#' \item{level.comb}{The level used to calculate confidence intervals
#'   for pooled estimates.}
#' \item{hakn}{A logical indicating whether the method by Hartung and
#'   Knapp is used to adjust test statistics and confidence
#'   intervals.}
#' \item{method.tau}{A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2}.}
#' \item{tau.preset}{Prespecified value for the square root of the
#'   between-study variance \eqn{\tau^2}.}
#' \item{TE.tau}{Overall treatment effect used to estimate the
#'   between-study variance \eqn{\tau^2}.}
#' \item{n.harmonic.mean}{Harmonic mean of number of observations (for
#'   back transformation of Freeman-Tukey Double arcsine
#'   transformation).}
#' \item{version}{Version of R package \bold{meta} used to create
#'   object.}
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
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
#' data(Fleiss93)
#' m1 <- metabin(event.e, n.e, event.c, n.c,
#'               data = Fleiss93, studlab = study,
#'               sm = "RR", method = "I")
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
#' m2 <- update(m1, title = "Fleiss93 meta-analysis",
#'              backtransf = FALSE)
#' metacum(m2)
#' 
#' data(Fleiss93cont)
#' m3 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c,
#'                data = Fleiss93cont, sm = "SMD")
#' metacum(m3)
#' 
#' @export metacum


metacum <- function(x, pooled, sortvar) {
  
  
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
  if (!missing(pooled))
    pooled <- setchar(pooled, c("fixed", "random"))
  else
    if (!x$comb.fixed & x$comb.random)
      pooled <- "random"
    else
      pooled <- "fixed"
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
    if (isCol(x$data, ".subset"))
      sortvar <- sortvar[x$data$.subset]
  }
  sort <- !is.null(sortvar)
  if (sort && (length(sortvar) != k.all))
    stop("Number of studies in object 'x' and argument 'sortvar' have different length.")
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
  for (i in 1:k.all)
    slab[i] <- paste("Adding ", studlab[i],
                     " (k=", ncum[i], ")", sep = "")
  slab <- c(slab, "Pooled estimate")
  studlab <- c(rev(rev(slab)[-1]), " ", rev(slab)[1])
  
  
  ##
  ##
  ## (4) Do sensitivity analysis
  ##
  ##
  res.i <- matrix(NA, ncol = 21, nrow = k.all)
  add.i <- matrix(NA, ncol = 3, nrow = k.all)
  ##
  for (i in 1:k.all) {
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
                   incr = incr.i, allincr = x$allincr, addincr = x$addincr,
                   allstudies = x$allstudies, MH.exact = x$MH.exact,
                   RR.Cochrane = x$RR.Cochrane, Q.Cochrane = x$Q.Cochrane,
                   model.glmm =
                     if (!is.null(x$model.glmm)) x$model.glmm else "UM.FS",
                   ##
                   level.comb = x$level.comb,
                   ##
                   hakn = x$hakn,
                   method.tau = x$method.tau,
                   tau.preset = x$tau.preset, TE.tau = x$TE.tau,
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
                    level.comb = x$level.comb,
                    ##
                    hakn = x$hakn,
                    method.tau = x$method.tau,
                    tau.preset = x$tau.preset, TE.tau = x$TE.tau,
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
                   level.comb = x$level.comb,
                   ##
                   hakn = x$hakn,
                   method.tau = x$method.tau,
                   tau.preset = x$tau.preset, TE.tau = x$TE.tau,
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
                   level.comb = x$level.comb,
                   ##
                   hakn = x$hakn,
                   method.tau = x$method.tau,
                   tau.preset = x$tau.preset, TE.tau = x$TE.tau,
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
                   incr = incr.i, allincr = x$allincr, addincr = x$addincr,
                   model.glmm =
                     if (!is.null(x$model.glmm)) x$model.glmm else "UM.FS",
                   ##
                   level.comb = x$level.comb,
                   ##
                   hakn = x$hakn, method.tau = x$method.tau,
                   tau.preset = x$tau.preset, TE.tau = x$TE.tau,
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
                    level.comb = x$level.comb,
                    ##
                    hakn = x$hakn,
                    method.tau = x$method.tau,
                    tau.preset = x$tau.preset, TE.tau = x$TE.tau,
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
                    incr = incr.i, allincr = x$allincr, addincr = x$addincr,
                    method.ci = x$method.ci,
                    ##
                    level.comb = x$level.comb,
                    ##
                    hakn = x$hakn,
                    method.tau = x$method.tau,
                    tau.preset = x$tau.preset, TE.tau = x$TE.tau,
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
                    incr = incr.i, allincr = x$allincr, addincr = x$addincr,
                    ##
                    level.comb = x$level.comb,
                    ##
                    hakn = x$hakn,
                    method.tau = x$method.tau,
                    tau.preset = x$tau.preset, TE.tau = x$TE.tau,
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
    if (pooled == "fixed") {
      res.i[i, ] <- c(m$TE.fixed,                                    #  1
                      m$seTE.fixed,                                  #  2
                      m$lower.fixed,                                 #  3
                      m$upper.fixed,                                 #  4
                      m$zval.fixed,                                  #  5
                      m$pval.fixed,                                  #  6
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
                      sum(m$w.fixed, na.rm = TRUE),                  # 17
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
                      m$zval.random,                                 #  5
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
                      if (x$hakn) m$df.hakn else NA,                 # 19
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
  zval.i <- res.i[, 5]
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
  if (pooled == "random" & x$hakn)
    df.hakn.i <- res.i[, 19]
  t.harmonic.mean.i <- res.i[, 20]
  Rb.i <- res.i[, 21]
  ##
  method.tau.ci <- unique(add.i[, 1])
  sign.lower.tau.i <- add.i[, 2]
  sign.upper.tau.i <- add.i[, 3]
  ##  
  if (pooled == "fixed") {
    TE.s <- x$TE.fixed
    seTE.s <- x$seTE.fixed
    lower.TE.s <- x$lower.fixed
    upper.TE.s <- x$upper.fixed
    zval.s <- x$zval.fixed
    pval.s <- x$pval.fixed
    w.s <- sum(x$w.fixed, na.rm = TRUE)
  }
  ##
  else if (pooled == "random") {
    TE.s <- x$TE.random
    seTE.s <- x$seTE.random
    lower.TE.s <- x$lower.random
    upper.TE.s <- x$upper.random
    zval.s <- x$zval.random
    pval.s <- x$pval.random
    w.s <- sum(x$w.random, na.rm = TRUE)
  }
  
  
  ##
  ##
  ## (5) Generate R object
  ##
  ##
  res <- list(TE = c(TE.i, NA, TE.s),
              seTE = c(seTE.i, NA, seTE.s),
              lower = c(lower.i, NA, lower.TE.s),
              upper = c(upper.i, NA, upper.TE.s),
              zval = c(zval.i, NA, zval.s),
              pval = c(pval.i, NA, pval.s),
              studlab = studlab,
              ##
              tau2 = c(tau2.i, NA, x$tau2),
              lower.tau2 = c(lower.tau2.i, NA, x$lower.tau2),
              upper.tau2 = c(upper.tau2.i, NA, x$upper.tau2),
              se.tau2 = c(se.tau2.i, NA, x$se.tau2),
              ##
              tau = c(tau.i, NA, x$tau),
              lower.tau = c(lower.tau.i, NA, x$lower.tau),
              upper.tau = c(upper.tau.i, NA, x$upper.tau),
              ##
              method.tau.ci = method.tau.ci,
              sign.lower.tau.i = c(sign.lower.tau.i, NA, x$sign.lower.tau),
              sign.upper.tau.i = c(sign.upper.tau.i, NA, x$sign.upper.tau),
              ##
              I2 = c(I2.i, NA, x$I2),
              lower.I2 = c(lower.I2.i, NA, x$lower.I2),
              upper.I2 = c(upper.I2.i, NA, x$upper.I2),
              ##
              Rb = c(Rb.i, NA, x$Rb),
              ##
              w = c(weight.i, NA, w.s),
              df.hakn = if (pooled == "random" & x$hakn) c(df.hakn.i, NA, x$df.hakn) else NULL,
              ##
              sm = x$sm, method = x$method, k = x$k,
              pooled = pooled,
              comb.fixed = ifelse(pooled == "fixed", TRUE, FALSE),
              comb.random = ifelse(pooled == "random", TRUE, FALSE),
              TE.fixed = NA, seTE.fixed = NA,
              TE.random = NA, seTE.random = NA,
              null.effect = x$null.effect,
              ##
              Q = NA,
              level.comb = x$level.comb,
              hakn = x$hakn,
              method.tau = x$method.tau,
              tau.preset = x$tau.preset,
              TE.tau = x$TE.tau,
              n.harmonic.mean = c(n.harmonic.mean.i, NA, 1 / mean(1 / n)),
              t.harmonic.mean = c(t.harmonic.mean.i, NA, 1 / mean(1 / time)),
              prediction = FALSE,
              ##
              backtransf = x$backtransf,
              pscale = x$pscale,
              irscale = x$irscale, irunit = x$irunit,
              title = x$title, complab = x$complab,
              outclab = x$outclab,
              ##
              call = match.call())
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metacum", "meta")
  if (inherits(x, "trimfill"))
    class(res) <- c(class(res), "trimfill")
  
  res
}
