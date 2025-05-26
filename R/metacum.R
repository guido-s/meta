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
#' @param overall A logical indicating whether overall results should be
#'   reported.
#' @param text.pooled A character string used in printouts and forest
#'   plots to label the pooled effect estimate.
#' @param no A numeric specifying which meta-analysis results to
#'   consider.
#' @param cid A numeric value or vector specifying clinically important
#'   differences (CID) / decision thresholds used to calculate expected
#'   proportions of clinically important benefit or harm
#'   (see \code{\link{cidprop}}).
#' @param cid.below.null A single numeric defining the decision threshold below
#'   the null effect to distinguish clinically important from not important
#'   effects (see \code{\link{cidprop}}).
#' @param cid.above.null A single numeric defining the decision threshold above
#'   the null effect to distinguish clinically important from not important
#'   effects (see \code{\link{cidprop}}).
#' @param small.values A character string specifying whether small
#'   treatment effects indicate a beneficial (\code{"desirable"}) or
#'   harmful (\code{"undesirable"}) effect, can be abbreviated
#'   (see \code{\link{cidprop}}).
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
#' An object of class \code{"metacum"} with dedicated print and forest
#' functions.
#' 
#' The following list elements provide results from meta-analyses, each
#' adding one study at a time (see \code{\link{meta-object}} for more
#' information on these list elements):
#' \tabular{l}{
#' \cr
#' studlab, TE, seTE, df.random, lower, upper, statistic, pval, \cr
#' lower.predict, upper.predict, df.predict, w (sum of weights), \cr
#' tau2, se.tau2, lower.tau2, upper.tau2, tau, lower.tau, upper.tau, \cr
#' I2, lower.I2, upper.I2, Rb, n.harmonic.mean, t.harmonic.mean, \cr
#' k, k.study, k.all, k.TE, k.MH.
#' }
#' 
#' The following list elements contain results of the original meta-analysis:
#' \tabular{l}{
#' \cr
#' TE.pooled, seTE.pooled, df.random.pooled, \cr
#' lower.pooled, upper.pooled, statistic.pooled, pval.pooled, \cr
#' lower.predict.pooled, upper.predict.pooled, \cr
#' df.predict.pooled, w.pooled, \cr
#' tau2.pooled, se.tau2.pooled, lower.tau2.pooled, upper.tau2.pooled, \cr
#' tau.pooled, lower.tau.pooled, upper.tau.pooled, \cr
#' I2.pooled, lower.I2.pooled, upper.I2.pooled, Rb.pooled, \cr
#' n.harmonic.mean.pooled, t.harmonic.mean.pooled, \cr
#' k.pooled, k.study.pooled, k.all.pooled, k.TE.pooled, k.MH.pooled.
#' }
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{forest.metacum}}, \code{\link{print.metacum}},
#'   \code{\link{cidprop}}
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

metacum.meta <- function(x, pooled, sortvar, prediction, overall = x$overall,
                         text.pooled, no = 1,
                         cid = NULL,
                         cid.below.null = NULL, cid.above.null = NULL,
                         small.values = "desirable",
                         ...) {
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  
  chkclass(x, "meta")
  chksuitable(x, "Cumulative meta-analysis",
              c("trimfill", "metamerge", "netpairwise"))
  #
  x <- updateversion(x)
  #
  if (!is.null(x$three.level) && x$three.level)
    stop("Cumulative meta-analysis not implemented for a ",
         "multi-level meta-analysis.",
         call. = FALSE)
  #
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
  chklogical(overall)
  #
  tdist_random <- pooled == "random" & x$method.random.ci %in% c("HK", "KR")
  tdist_predict <- !(x$method.predict %in% c("S", ""))
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
  #
  if (!is.null(cid))
    chknumeric(cid, length = 1)
  #
  if (!is.null(cid.below.null))
    chknumeric(cid.below.null, length = 1)
  #
  if (!is.null(cid.above.null))
    chknumeric(cid.above.null, length = 1)
  #
  avail.cid <- !is.null(cid) && !all(is.na(cid))
  avail.cid.below.null <-
    !is.null(cid.below.null) && !all(is.na(cid.below.null))
  avail.cid.above.null <-
    !is.null(cid.above.null) && !all(is.na(cid.above.null))
  #
  run_cidprop <- avail.cid | avail.cid.below.null | avail.cid.above.null
  
  
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
  if (!is.null(x$cluster))
    cluster <- x$cluster[o]
  else
    cluster <- NULL
  #
  if (!is.null(x$weights.common))
    weights.common <- x$weights.common[o]
  else
    weights.common <- NULL
  #
  if (!is.null(x$weights.random))
    weights.random <- x$weights.random[o]
  else
    weights.random <- NULL
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
  x$text.common <- x$text.common[no.c]
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
  
  res.i <- matrix(NA, ncol = 31, nrow = k.all)
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
                   #
                   cluster = cluster[sel], rho = x$rho,
                   #
                   weights.common = weights.common[sel],
                   weights.random = weights.random[sel],
                   #
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
                   method.common.ci = x$method.common.ci,
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
                    #
                    cluster = cluster[sel], rho = x$rho,
                    #
                    weights.common = weights.common[sel],
                    weights.random = weights.random[sel],
                    #
                    sm = x$sm,
                    pooledvar = replaceNA(x$pooledvar, gs("pooledvar")),
                    ##
                    method.tau = x$method.tau,
                    tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                    ##
                    level.ma = x$level.ma,
                    method.common.ci = x$method.common.ci,
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
                   #
                   cluster = cluster[sel], rho = x$rho,
                   #
                   weights.common = weights.common[sel],
                   weights.random = weights.random[sel],
                   #
                   sm = x$sm, null.effect = x$null.effect,
                   ##
                   method.tau = x$method.tau,
                   tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                   ##
                   level.ma = x$level.ma,
                   method.common.ci = x$method.common.ci,
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
                   #
                   cluster = cluster[sel], rho = x$rho,
                   #
                   weights.common = weights.common[sel],
                   weights.random = weights.random[sel],
                   #
                   sm = x$sm, null.effect = x$null.effect,
                   ##
                   method.tau = x$method.tau,
                   tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                   ##
                   level.ma = x$level.ma,
                   method.common.ci = x$method.common.ci,
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
                   #
                   cluster = cluster[sel], rho = x$rho,
                   #
                   weights.common = weights.common[sel],
                   weights.random = weights.random[sel],
                   #
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
                   method.common.ci = x$method.common.ci,
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
                    #
                    cluster = cluster[sel], rho = x$rho,
                    #
                    weights.common = weights.common[sel],
                    weights.random = weights.random[sel],
                    #
                    sm = x$sm, null.effect = x$null.effect,
                    ##
                    method.tau = x$method.tau,
                    tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                    ##
                    level.ma = x$level.ma,
                    method.common.ci = x$method.common.ci,
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
                    #
                    cluster = cluster[sel], rho = x$rho,
                    #
                    weights.common = weights.common[sel],
                    weights.random = weights.random[sel],
                    #
                    method = x$method, sm = x$sm, null.effect = x$null.effect,
                    ##
                    incr = incr.i, method.incr = x$method.incr,
                    method.ci = x$method.ci,
                    ##
                    method.tau = x$method.tau,
                    tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                    ##
                    level.ma = x$level.ma,
                    method.common.ci = x$method.common.ci,
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
                    #
                    cluster = cluster[sel], rho = x$rho,
                    #
                    weights.common = weights.common[sel],
                    weights.random = weights.random[sel],
                    #
                    method = x$method, sm = x$sm, null.effect = x$null.effect,
                    ##
                    incr = incr.i, method.incr = x$method.incr,
                    ##
                    method.tau = x$method.tau,
                    tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                    ##
                    level.ma = x$level.ma,
                    method.common.ci = x$method.common.ci,
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
                      #
                      m$tau2,                                        #  7
                      m$lower.tau2,                                  #  8
                      m$upper.tau2,                                  #  9
                      m$se.tau2,                                     # 10
                      #
                      m$tau,                                         # 11
                      m$lower.tau,                                   # 12
                      m$upper.tau,                                   # 13
                      #
                      m$I2,                                          # 14
                      m$lower.I2,                                    # 15
                      m$upper.I2,                                    # 16
                      #
                      m$Rb,                                          # 17
                      #
                      sum(m$w.common, na.rm = TRUE),                 # 18
                      if (sel.pft) 1 / mean(1 / n[sel]) else NA,     # 19
                      if (sel.irft) 1 / mean(1 / time[sel]) else NA, # 20
                      #
                      NA,                                            # 21
                      #
                      NA,                                            # 22
                      NA,                                            # 23
                      NA,                                            # 24
                      #
                      m$k,                                           # 25
                      m$k.study,                                     # 26
                      m$k.all,                                       # 27
                      m$k.TE,                                        # 28
                      replaceNULL(m$k.MH),                           # 29
                      NA,                                            # 30
                      NA                                             # 31
      )
    }
    ##
    else if (pooled == "random") {
      if (run_cidprop) {
        pp <-
          cidprop(m,
                   cid = cid,
                   cid.below.null = cid.below.null,
                   cid.above.null = cid.above.null,
                   small.values = small.values)
        #
        prop.cid.below.null <- pp$prop.cid.below.null
        prop.cid.above.null <- pp$prop.cid.above.null
      }
      else {
        prop.cid.below.null <- NA
        prop.cid.above.null <- NA
      }
      #
      res.i[i, ] <- c(m$TE.random,                                   #  1
                      m$seTE.random,                                 #  2
                      m$lower.random,                                #  3
                      m$upper.random,                                #  4
                      m$statistic.random,                            #  5
                      m$pval.random,                                 #  6
                      #
                      m$tau2,                                        #  7
                      m$lower.tau2,                                  #  8
                      m$upper.tau2,                                  #  9
                      m$se.tau2,                                     # 10
                      #
                      m$tau,                                         # 11
                      m$lower.tau,                                   # 12
                      m$upper.tau,                                   # 13
                      #
                      m$I2,                                          # 14
                      m$lower.I2,                                    # 15
                      m$upper.I2,                                    # 16
                      #
                      m$Rb,                                          # 17
                      #
                      sum(m$w.random, na.rm = TRUE),                 # 18
                      if (sel.pft) 1 / mean(1 / n[sel]) else NA,     # 19
                      if (sel.irft) 1 / mean(1 / time[sel]) else NA, # 20
                      #
                      if (tdist_random) m$df.random else NA,         # 21
                      #
                      m$lower.predict,                               # 22
                      m$upper.predict,                               # 23
                      if (tdist_predict) m$df.predict else NA,       # 24
                      #
                      m$k,                                           # 25
                      m$k.study,                                     # 26
                      m$k.all,                                       # 27
                      m$k.TE,                                        # 28
                      replaceNULL(m$k.MH),                           # 29
                      prop.cid.below.null,                           # 30
                      prop.cid.above.null                            # 31
      )
    }
  }
  #
  TE.i <- res.i[, 1]
  seTE.i <- res.i[, 2]
  lower.i <- res.i[, 3]
  upper.i <- res.i[, 4]
  statistic.i <- res.i[, 5]
  pval.i <- res.i[, 6]
  #
  tau2.i <- res.i[, 7]
  lower.tau2.i <- res.i[, 8]
  upper.tau2.i <- res.i[, 9]
  se.tau2.i <- res.i[, 10]
  #
  tau.i <- res.i[, 11]
  lower.tau.i <- res.i[, 12]
  upper.tau.i <- res.i[, 13]
  #
  I2.i <- res.i[, 14]
  lower.I2.i <- res.i[, 15]
  upper.I2.i <- res.i[, 16]
  #
  Rb.i <- res.i[, 17]
  #
  weight.i <- res.i[, 18]
  n.harmonic.mean.i <- res.i[, 19]
  t.harmonic.mean.i <- res.i[, 20]
  #
  df.random.i <- res.i[, 21]
  #
  lower.predict.i <- res.i[, 22]
  upper.predict.i <- res.i[, 23]
  df.predict.i <- res.i[, 24]
  #
  k.i <- res.i[, 25]
  k.study.i <- res.i[, 26]
  k.all.i <- res.i[, 27]
  k.TE.i <- res.i[, 28]
  k.MH.i <- res.i[, 29]
  #
  prop.cid.below.null <- res.i[, 30]
  prop.cid.above.null <- res.i[, 31]
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
  #
  if (pooled == "random" & run_cidprop) {
    pp.pooled <-
      cidprop(x,
               cid = cid,
               cid.below.null = cid.below.null,
               cid.above.null = cid.above.null,
               small.values = small.values)
    #
    prop.cid.below.null.pooled <- pp.pooled$prop.cid.below.null
    prop.cid.above.null.pooled <- pp.pooled$prop.cid.above.null
  }
  else {
    prop.cid.below.null.pooled <- NA
    prop.cid.above.null.pooled <- NA
  }
  
  
  ##
  ##
  ## (5) Generate R object
  ##
  ##
  
  if (missing(text.pooled))
    text.pooled <- if (pooled == "common") x$text.common else x$text.random
  #
  res <- list(studlab = slab,
              #
              TE = TE.i,
              seTE = seTE.i,
              df.random = df.random.i,
              lower = lower.i,
              upper = upper.i,
              statistic = statistic.i,
              pval = pval.i,
              #
              lower.predict = lower.predict.i,
              upper.predict = upper.predict.i,
              df.predict = df.predict.i,
              #
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
              sign.lower.tau = sign.lower.tau.i,
              sign.upper.tau = sign.upper.tau.i,
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
              prop.cid.below.null = prop.cid.below.null,
              prop.cid.above.null = prop.cid.above.null,
              cid.below.null = cid.below.null,
              cid.above.null = cid.above.null,
              small.values = small.values,
              #
              sm = x$sm,
              null.effect = x$null.effect,
              #
              pooled = pooled,
              common = pooled == "common",
              random = pooled == "random",
              overall = overall,
              overall.hetstat = FALSE,
              #
              prediction = prediction,
              method.predict = x$method.predict,
              adhoc.hakn.pi = x$adhoc.hakn.pi,
              #
              backtransf = x$backtransf,
              func.backtransf = x$func.backtransf,
              #
              level = x$level.ma,
              level.ma = x$level.ma,
              level.predict = x$level.predict,
              #
              method = x$method,
              method.random = x$method.random,
              #
              method.common.ci = x$method.common.ci,
              method.random.ci = x$method.random.ci,
              adhoc.hakn.ci = x$adhoc.hakn.ci,
              #
              method.tau = x$method.tau,
              method.tau.ci =
                if (length(method.tau.ci) > 0) method.tau.ci else "",
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
              TE.pooled = TE.pooled,
              seTE.pooled = seTE.pooled,
              lower.pooled = lower.pooled,
              upper.pooled = upper.pooled,
              df.random.pooled = df.random.pooled,
              statistic.pooled = statistic.pooled,
              pval.pooled = pval.pooled,
              w.pooled = w.pooled,
              text.pooled = text.pooled,
              #
              lower.predict.pooled = lower.predict.pooled,
              upper.predict.pooled = upper.predict.pooled,
              df.predict.pooled = x$df.predict,
              text.predict = x$text.predict,
              #
              prop.cid.below.null.pooled = prop.cid.below.null.pooled,
              prop.cid.above.null.pooled = prop.cid.above.null.pooled,
              #
              Q.pooled = x$Q,
              Q.Cochrane = x$Q.Cochrane,
              #
              tau2.pooled = tau2.pooled,
              se.tau2.pooled = se.tau2.pooled,
              lower.tau2.pooled = lower.tau2.pooled,
              upper.tau2.pooled = upper.tau2.pooled,
              #
              tau.pooled = tau.pooled,
              lower.tau.pooled = x$lower.tau,
              upper.tau.pooled = upper.tau.pooled,
              #
              sign.lower.tau.pooled = x$sign.lower.tau,
              sign.upper.tau.pooled = x$sign.upper.tau,
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
              label.left = x$label.left,
              label.right = x$label.right,
              #
              title = x$title, complab = x$complab,
              outclab = x$outclab,
              ##
              no = no,
              ##
              x = x,
              ##
              call = match.call())
  #
  if (!is.null(weights.common) & pooled == "common") {
    res$weights.common <- TRUE
    res$weights.random <- FALSE
  }
  #
  if (!is.null(weights.random) & pooled == "random") {
    res$weights.common <- FALSE
    res$weights.random <- TRUE
  }
  #
  if (run_cidprop) {
    res$cid.below.null <- pp$cid.below.null
    res$cid.above.null <- pp$cid.above.null
  }
  else {
    res$prop.cid.below.null <- NULL
    res$prop.cid.above.null <- NULL
    res$cid.below.null <- NULL
    res$cid.above.null <- NULL
    res$small.values <- NULL
    #
    res$prop.cid.below.null.pooled <- NULL
    res$prop.cid.above.null.pooled <- NULL
  }
  #
  res$version <- packageDescription("meta")$Version
  ##
  res$x$common <- res$common
  res$x$random <- res$random
  #
  res$classes <- class(x)[class(x) != "meta"]
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
