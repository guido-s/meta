#' Trim-and-fill method to adjust for bias in meta-analysis
#' 
#' @description
#' Trim-and-fill method for estimating and adjusting for the number
#' and outcomes of missing studies in a meta-analysis.
#' 
#' @aliases trimfill trimfill.meta trimfill.default
#' 
#' @param x An object of class \code{meta}, or estimated treatment
#'   effect in individual studies.
#' @param seTE Standard error of estimated treatment effect.
#' @param left A logical indicating whether studies are supposed to be
#'   missing on the left or right side of the funnel plot. If NULL,
#'   the linear regression test for funnel plot symmetry (i.e.,
#'   function \code{metabias(..., method="Egger")}) is used to
#'   determine whether studies are missing on the left or right side.
#' @param ma.common A logical indicating whether a common effect or
#'   random effects model is used to estimate the number of missing
#'   studies.
#' @param type A character indicating which method is used to estimate
#'   the number of missing studies. Either \code{"L"} or \code{"R"}.
#' @param n.iter.max Maximum number of iterations to estimate number
#'   of missing studies.
#' @param sm An optional character string indicating underlying
#'   summary measure, e.g., \code{"RD"}, \code{"RR"}, \code{"OR"},
#'   \code{"ASD"}, \code{"HR"}, \code{"MD"}, \code{"SMD"}, or
#'   \code{"ROM"}; ignored if \code{x} is of class \code{meta}.
#' @param studlab An optional vector with study labels; ignored if
#'   \code{x} is of class \code{meta}.
#' @param level The level used to calculate confidence intervals for
#'   individual studies. If existing, \code{x$level} is used as value
#'   for \code{level}; otherwise 0.95 is used.
#' @param level.ma The level used to calculate confidence interval for
#'   the pooled estimate. If existing, \code{x$level.ma} is used as
#'   value for \code{level.ma}; otherwise 0.95 is used.
#' @param common A logical indicating whether a common effect
#'   meta-analysis should be conducted.
#' @param random A logical indicating whether a random effects
#'   meta-analysis should be conducted.
#' @param method.common.ci A character string indicating which method
#'   is used to calculate confidence interval and test statistic for
#'   common effect estimate (see \code{\link{meta-package}}).
#' @param method.random.ci A character string indicating which method
#'   is used to calculate confidence interval and test statistic for
#'   random effects estimate (see \code{\link{meta-package}}).
#' @param adhoc.hakn.ci A character string indicating whether an
#'   \emph{ad hoc} variance correction should be applied in the case
#'   of an arbitrarily small Hartung-Knapp variance estimate (see
#'   \code{\link{meta-package}}).
#' @param method.predict A character string indicating which method is
#'   used to calculate a prediction interval (see
#'   \code{\link{meta-package}}).
#' @param adhoc.hakn.pi A character string indicating whether an
#'   \emph{ad hoc} variance correction should be applied for the
#'   prediction interval (see \code{\link{meta-package}}).
#' @param seed.predict A numeric value used as seed to calculate
#'   bootstrap prediction interval (see \code{\link{meta-package}}).
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2} and its
#'   square root \eqn{\tau} (see \code{\link{meta-package}}).
#' @param method.tau.ci A character string indicating which method is
#'   used to estimate the confidence interval of \eqn{\tau^2} and
#'   \eqn{\tau} (see \code{\link{meta-package}}).
#' @param level.hetstat The level used to calculate confidence intervals
#'   for heterogeneity statistics.
#' @param prediction A logical indicating whether a prediction
#'   interval should be printed.
#' @param level.predict The level used to calculate prediction
#'   interval for a new study.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=TRUE}, results for \code{sm="OR"} are printed as
#'   odds ratios rather than log odds ratios and results for
#'   \code{sm="ZCOR"} are printed as correlations rather than Fisher's
#'   z transformed correlations, for example.
#' @param pscale A numeric giving scaling factor for printing of
#'   single event probabilities or risk differences, i.e. if argument
#'   \code{sm} is equal to \code{"PLOGIT"}, \code{"PLN"},
#'   \code{"PRAW"}, \code{"PAS"}, \code{"PFT"}, or \code{"RD"}.
#' @param irscale A numeric defining a scaling factor for printing of
#'   single incidence rates or incidence rate differences, i.e. if
#'   argument \code{sm} is equal to \code{"IR"}, \code{"IRLN"},
#'   \code{"IRS"}, \code{"IRFT"}, or \code{"IRD"}.
#' @param irunit A character specifying the time unit used to
#'   calculate rates, e.g. person-years.
#' @param silent A logical indicating whether basic information on
#'   iterations shown.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
#'
#' @details
#' The trim-and-fill method (Duval, Tweedie 2000a, 2000b) can be used
#' for estimating and adjusting for the number and outcomes of missing
#' studies in a meta-analysis. The method relies on scrutiny of one
#' side of a funnel plot for asymmetry assumed due to publication
#' bias.
#' 
#' Three different methods have been proposed originally to estimate
#' the number of missing studies. Two of these methods (L- and
#' R-estimator) have been shown to perform better in simulations, and
#' are available in this R function (argument \code{type}).
#' 
#' A common effect or random effects model can be used to estimate the
#' number of missing studies (argument \code{ma.common}). Furthermore,
#' a common effect and/or random effects model can be used to
#' summaries study results (arguments \code{common} and
#' \code{random}). Simulation results (Peters et al. 2007) indicate
#' that the common-random model, i.e. using a common effect model to
#' estimate the number of missing studies and a random effects model
#' to summaries results, (i) performs better than the common-common
#' model, and (ii) performs no worse than and marginally better in
#' certain situations than the random-random model. Accordingly, the
#' common-random model is the default.
#' 
#' An empirical comparison of the trim-and-fill method and the Copas
#' selection model (Schwarzer et al. 2010) indicates that the
#' trim-and-fill method leads to excessively conservative inference in
#' practice. The Copas selection model is available in R package
#' \bold{metasens}.
#' 
#' The function \code{\link{metagen}} is called internally.
#' 
#' @return
#' An object of class \code{c("trimfill", "metagen", "meta")} with
#' corresponding generic functions (see \code{\link{meta-object}}).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metagen}}, \code{\link{metabias}},
#'   \code{\link{funnel}}
#' 
#' @references
#' Duval S & Tweedie R (2000a):
#' A nonparametric "Trim and Fill" method of accounting for
#' publication bias in meta-analysis.
#' \emph{Journal of the American Statistical Association},
#' \bold{95}, 89--98
#' 
#' Duval S & Tweedie R (2000b):
#' Trim and Fill: A simple funnel-plot-based method of testing and
#' adjusting for publication bias in meta-analysis.
#' \emph{Biometrics},
#' \bold{56}, 455--63
#' 
#' Peters JL, Sutton AJ, Jones DR, Abrams KR, Rushton L (2007):
#' Performance of the trim and fill method in the presence of
#' publication bias and between-study heterogeneity.
#' \emph{Statisics in Medicine},
#' \bold{10}, 4544--62
#' 
#' Schwarzer G, Carpenter J, Rücker G (2010):
#' Empirical evaluation suggests Copas selection model preferable to
#' trim-and-fill method for selection bias in meta-analysis
#' \emph{Journal of Clinical Epidemiology},
#' \bold{63}, 282--8
#' 
#' @examples
#' data(Fleiss1993bin)
#' m1 <- metabin(d.asp, n.asp, d.plac, n.plac, data = Fleiss1993bin, sm = "OR")
#' tf1 <- trimfill(m1)
#' tf1
#' funnel(tf1)
#' funnel(tf1, pch = ifelse(tf1$trimfill, 1, 16), level = 0.9, random = FALSE)
#' #
#' # Use log odds ratios on x-axis
#' #
#' funnel(tf1, backtransf = FALSE)
#' funnel(tf1, pch = ifelse(tf1$trimfill, 1, 16), level = 0.9, random = FALSE,
#'   backtransf = FALSE)
#' 
#' trimfill(m1$TE, m1$seTE, sm = m1$sm)
#' 
#' @rdname trimfill
#' @method trimfill meta
#' @export


trimfill.meta <- function(x, left = NULL, ma.common = TRUE,
                          type = "L", n.iter.max = 50,
                          ##
                          common = FALSE, random = TRUE,
                          prediction = x$prediction,
                          ##
                          backtransf = x$backtransf, pscale = x$pscale,
                          irscale = x$irscale, irunit = x$irunit,
                          silent = TRUE,
                          warn.deprecated = gs("warn.deprecated"),
                          ...) {
  
  ##
  ##
  ## (1) Check for (inadmissible) meta object
  ##
  ##
  chkclass(x, "meta")
  chksuitable(x, "Trim-and-fill method",
              c("metacum", "metainf", "metamerge", "netpairwise"))
  #
  if (!is.null(x$weights.common) | !is.null(x$weights.random))
    stop("Trim-and-fill method not implemented for user-specified weights.",
         call. = FALSE)
  #
  x <- updateversion(x)
  ##
  if (x$three.level)
    stop("Trim-and-fill method not implemented for three-level model.",
         call. = FALSE)
  
  
  ##
  ##
  ## (2) Check arguments
  ##
  ##
  args <- list(...)
  ##
  ma.common <- deprecated(ma.common, missing(ma.common), args, "ma.fixed",
                          warn.deprecated)
  chklogical(ma.common)
  ##
  type <- setchar(type, c("L", "R"))
  chknumeric(n.iter.max, min = 1, length = 1)
  ##
  common <- deprecated(common, missing(common), args, "fixed",
                       warn.deprecated)
  chklogical(common)
  chklogical(random)
  chklogical(prediction)
  chklogical(backtransf)
  ##
  sm <- x$sm
  if (!is_prop(sm))
    pscale <- 1
  chknumeric(pscale, length = 1)
  if (!backtransf & pscale != 1) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  if (!is_rate(sm))
    irscale <- 1
  chknumeric(irscale, length = 1)
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
  else if (is_log_effect(sm))
    transf.null.effect <- log(null.effect)
  else if (sm == c("PLOGIT"))
    transf.null.effect <- log(null.effect / (1 - null.effect))
  else if (sm %in% c("IRS", "IRFT"))
    transf.null.effect <- sqrt(null.effect)
  else if (sm == "ZCOR")
    transf.null.effect <- 0.5 * log((1 + null.effect) / (1 - null.effect))
  
  
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
    left <- as.logical(sign(metabias(TE, seTE, method = "Egger",
                                     k.min = 3)$estimate[1]) == 1)
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
  
  
  if (ma.common)
    TE.sum <-
      metagen(TE, seTE, method.tau = "DL", method.tau.ci = "")$TE.common
  else
    TE.sum <-
      metagen(TE, seTE, method.tau = x$method.tau, method.tau.ci = "")$TE.random
  
  
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
      if (ma.common)
        TE.sum <- metagen(TE[sel], seTE[sel],
                          method.tau = "DL",
                          method.tau.ci = "")$TE.common
      else
        TE.sum <- metagen(TE[sel], seTE[sel],
                          method.tau = x$method.tau,
                          method.tau.ci = "")$TE.random
      ##
      trim1 <- estimate.missing(TE, TE.sum, type)
      ##
      if (!silent) {
        cat(paste0("n.iter = ", n.iter, "\n"))
        if (type == "L")
          cat(paste0("L0 = ", round(trim1$res0, 2), "\n\n"))
        if (type == "R")
          cat(paste0("R0 = ", round(trim1$res0 + 0.5, 2), "\n\n"))
      }
      ##
      k0 <- trim1$res0.plus
    }
  }
  
  
  if (k0 > (k - 1))
    k0 <- k - 1
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
    m <- metagen(-TE, seTE, studlab = studlab, level = x$level,
                 ##
                 level.ma = x$level.ma,
                 method.common.ci = x$method.common.ci,
                 method.random.ci = x$method.random.ci,
                 adhoc.hakn.ci = x$adhoc.hakn.ci,
                 ##
                 method.tau = x$method.tau, method.tau.ci = x$method.tau.ci,
                 level.hetstat = x$level.hetstat,
                 #
                 method.I2 = replaceNULL(x$method.I2, gs("method.I2")),
                 #
                 method.predict = x$method.predict,
                 prediction = prediction, level.predict = x$level.predict,
                 adhoc.hakn.pi = x$adhoc.hakn.pi,
                 seed.predict = x$seed.predict,
                 ##
                 null.effect = transf.null.effect)
  else
    m <- metagen(TE, seTE, studlab = studlab, level = x$level,
                 ##
                 level.ma = x$level.ma,
                 method.common.ci = x$method.common.ci,
                 method.random.ci = x$method.random.ci,
                 adhoc.hakn.ci = x$adhoc.hakn.ci,
                 ##
                 method.tau = x$method.tau, method.tau.ci = x$method.tau.ci,
                 level.hetstat = x$level.hetstat,
                 #
                 method.I2 = replaceNULL(x$method.I2, gs("method.I2")),
                 #
                 method.predict = x$method.predict,
                 prediction = prediction, level.predict = x$level.predict,
                 adhoc.hakn.pi = x$adhoc.hakn.pi,
                 seed.predict = x$seed.predict,
                 ##
                 null.effect = transf.null.effect)
  
  
  ##
  ## Calculate H, I-Squared, and Rb
  ##
  Hres  <- list(TE = m$H, lower = m$lower.H, upper = m$upper.H)
  I2res <- list(TE = m$I2, lower = m$lower.I2, upper = m$upper.I2)
  Rbres <- list(TE = m$Rb, lower = m$lower.Rb, upper = m$upper.Rb)
  
  
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
      m <- metagen(-TE.all, seTE.all, studlab = studlab.all,
                   exclude = exclude, level = x$level,
                   ##
                   level.ma = x$level.ma,
                   method.common.ci = x$method.common.ci,
                   method.random.ci = x$method.random.ci,
                   adhoc.hakn.ci = x$adhoc.hakn.ci,
                   ##
                   method.tau = x$method.tau, method.tau.ci = x$method.tau.ci,
                   level.hetstat = x$level.hetstat,
                   #
                   method.I2 = replaceNULL(x$method.I2, gs("method.I2")),
                   #
                   method.predict = x$method.predict,
                   prediction = prediction,
                   level.predict = x$level.predict,
                   adhoc.hakn.pi = x$adhoc.hakn.pi,
                   seed.predict = x$seed.predict,
                   ##
                   null.effect = transf.null.effect)
    else
      m <- metagen(TE.all, seTE.all, studlab = studlab.all,
                   exclude = exclude, level = x$level,
                   ##
                   level.ma = x$level.ma,
                   method.common.ci = x$method.common.ci,
                   method.random.ci = x$method.random.ci,
                   adhoc.hakn.ci = x$adhoc.hakn.ci,
                   ##
                   method.tau = x$method.tau, method.tau.ci = x$method.tau.ci,
                   level.hetstat = x$level.hetstat,
                   #
                   method.I2 = replaceNULL(x$method.I2, gs("method.I2")),
                   #
                   method.predict = x$method.predict,
                   prediction = prediction,
                   level.predict = x$level.predict,
                   adhoc.hakn.pi = x$adhoc.hakn.pi,
                   seed.predict = x$seed.predict,
                   ##
                   null.effect = transf.null.effect)
  }
  
  
  res <- list(studlab = m$studlab,
              ##
              sm = sm,
              null.effect = x$null.effect,
              ##
              TE = m$TE, seTE = m$seTE,
              statistic = m$statistic, pval = m$pval,
              df = rep(NA, length(m$TE)),
              level = x$level,
              lower = m$lower, upper = m$upper,
              ##
              three.level = FALSE,
              cluster = NULL,
              ##
              k = m$k,
              k.study = m$k.study,
              k.all = m$k.all,
              k.TE = m$k.TE,
              k0 = k0,
              ##
              overall = common | random,
              overall.hetstat = TRUE,
              common = common,
              random = random,
              prediction = prediction,
              backtransf = backtransf,
              ##
              method = m$method,
              method.random = m$method.random,
              level = x$level,
              ##
              w.common = m$w.common,
              TE.common = m$TE.common, seTE.common = m$seTE.common,
              statistic.common = m$statistic.common,
              pval.common = m$pval.common,
              level.ma = x$level.ma,
              lower.common = m$lower.common,
              upper.common = m$upper.common,
              ##
              w.random = m$w.random,
              TE.random = m$TE.random, seTE.random = m$seTE.random,
              statistic.random = m$statistic.random,
              pval.random = m$pval.random,
              method.common.ci = x$method.common.ci,
              method.random.ci = x$method.random.ci,
              df.random = m$df.random,
              lower.random = m$lower.random, upper.random = m$upper.random,
              ##
              seTE.classic = m$seTE.classic,
              ##
              adhoc.hakn.ci = m$adhoc.hakn.ci,
              df.hakn = m$df.hakn,
              seTE.hakn.ci = m$seTE.hakn.ci,
              seTE.hakn.adhoc.ci = m$seTE.hakn.adhoc.ci,
              ##
              df.kero = m$df.kero,
              seTE.kero = m$seTE.kero,
              ##
              method.predict = m$method.predict,
              adhoc.hakn.pi = m$adhoc.hakn.pi,
              seTE.predict = m$seTE.predict,
              df.predict = m$df.predict,
              level.predict = x$level.predict,
              lower.predict = m$lower.predict,
              upper.predict = m$upper.predict,
              seTE.hakn.pi = m$seTE.hakn.pi,
              seTE.hakn.adhoc.pi = m$seTE.hakn.adhoc.pi,
              seed.predict = m$seed.predict,
              ##
              Q = m$Q, df.Q = m$df.Q, pval.Q = m$pval.Q,
              ##
              method.tau = m$method.tau, method.tau.ci = m$method.tau.ci,
              level.hetstat = m$level.hetstat,
              tau2 = m$tau2,
              se.tau2 = m$se.tau2,
              lower.tau2 = m$lower.tau2, upper.tau2 = m$upper.tau2,
              tau = m$tau,
              lower.tau = m$lower.tau, upper.tau = m$upper.tau,
              detail.tau = m$detail.tau,
              sign.lower.tau = m$sign.lower.tau,
              sign.upper.tau = m$sign.upper.tau,
              #
              method.I2 = m$method.I2,
              #
              H = Hres$TE, lower.H = Hres$lower, upper.H = Hres$upper,
              ##
              I2 = I2res$TE, lower.I2 = I2res$lower, upper.I2 = I2res$upper,
              ##
              Rb = Rbres$TE, lower.Rb = Rbres$lower, upper.Rb = Rbres$upper,
              ##
              method.bias = m$method.bias,
              ##
              text.common = x$text.common,
              text.random = x$text.random,
              text.predict = x$text.predict,
              text.w.common = x$text.w.common,
              text.w.random = x$text.w.random,
              ##
              title = x$title,
              complab = x$complab,
              outclab = x$outclab,
              label.e = x$label.e,
              label.c = x$label.c,
              label.left = x$label.left,
              label.right = x$label.right,
              ##
              keepdata = FALSE,
              exclude = exclude.na,
              ##
              left = left,
              ma.common = ma.common,
              type = type,
              n.iter.max = n.iter.max,
              n.iter = n.iter,
              ##
              trimfill = trimfill,
              ##
              class.x = class(x)[1]
              )
  ##
  res$pscale <- pscale
  res$irscale <- irscale
  res$irunit <- irunit
  ##
  res$n.e <- n.e
  res$n.c <- n.c
  res$n <- n
  ##
  res$event.e <- event.e
  res$event.c <- event.c
  res$event <- event
  ##
  res$time.e <- time.e
  res$time.c <- time.c
  res$time <- time
  ##
  res$cor <- cor
  ##
  res$mean.e <- mean.e
  res$mean.c <- mean.c
  ##
  res$sd.e <- sd.e
  res$sd.c <- sd.c
  ##
  call <- match.call()
  res$version <- packageDescription("meta")$Version
  ##
  ## Backward compatibility
  ##
  res <- backward(res)
  ##
  res$ma.fixed <- res$ma.common
  res$hakn <- m$hakn
  ##  
  class(res) <- c("metagen", "meta", "trimfill")
  
  
  res
}





#' @rdname trimfill
#' @method trimfill default
#' @export


trimfill.default <- function(x, seTE, left = NULL, ma.common = TRUE,
                             type = "L", n.iter.max = 50,
                             sm = "", studlab = NULL,
                             level = 0.95, level.ma = level,
                             common = FALSE, random = TRUE,
                             ##
                             method.common.ci = gs("method.common.ci"),
                             method.random.ci = gs("method.random.ci"),
                             adhoc.hakn.ci = gs("adhoc.hakn.ci"),
                             ##
                             method.tau = gs("method.tau"),
                             method.tau.ci =
                               if (method.tau == "DL") "J" else "QP",
                             level.hetstat = gs("level.hetstat"),
                             #
                             prediction = FALSE, level.predict = level,
                             method.predict = gs("method.predict"),
                             adhoc.hakn.pi = gs("adhoc.hakn.pi"),
                             seed.predict = NULL,
                             ##
                             backtransf = TRUE, pscale = 1,
                             irscale = 1, irunit = "person-years",
                             silent = TRUE, ...) {
  
  
  ##
  ##
  ## (1) Check essential arguments
  ##
  ##
  k.All <- length(x)
  ##
  chknumeric(x)
  chknumeric(seTE)
  chknull(sm)
  ##
  fun <- "trimfill"
  chklength(seTE, k.All, fun)
  ##
  if (!is.null(studlab))
    chklength(studlab, k.All, fun)
  else
    studlab <- seq(along = x)
  ##
  if (is.null(sm))
    sm <- ""
  
  
  ##
  ##
  ## (2) Do meta-analysis
  ##
  ##
  m <- metagen(x, seTE, studlab = studlab, sm = sm,
               level = level, level.ma = level.ma,
               common = common, random = random,
               ##
               method.common.ci = method.common.ci,
               method.random.ci = method.random.ci,
               adhoc.hakn.ci = adhoc.hakn.ci,
               ##
               method.tau = method.tau, method.tau.ci = method.tau.ci,
               level.hetstat = level.hetstat,
               #
               method.I2 = gs("method.I2"),
               #
               prediction = prediction,
               method.predict = method.predict,
               adhoc.hakn.pi = adhoc.hakn.pi,
               level.predict = level.predict,
               seed.predict = seed.predict,
               ##
               backtransf = backtransf, pscale = pscale,
               irscale = irscale, irunit = irunit,
               ...)
  
  
  ##
  ##
  ## (3) Run trim-and-fill method
  ##
  ##
  res <- trimfill(m, left = left, ma.common = ma.common,
                  type = type, n.iter.max = n.iter.max,
                  silent = silent, ...)
  
  
  res
}
