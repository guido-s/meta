#' Test for funnel plot asymmetry
#' 
#' @description
#' Test for funnel plot asymmetry, based on rank correlation or linear
#' regression method.
#' 
#' @aliases metabias metabias.meta metabias.default
#' 
#' @param x An object of class \code{meta} or estimated treatment
#'   effect in individual studies.
#' @param seTE Standard error of estimated treatment effect (mandatory
#'   if \code{x} not of class \code{meta}).
#' @param method.bias A character string indicating which test is to
#'   be used.  Either \code{"rank"}, \code{"linreg"}, \code{"mm"},
#'   \code{"count"}, \code{"score"}, or \code{"peters"}, can be
#'   abbreviated.
#' @param plotit A logical indicating whether a plot should be
#'   produced for method.bias \code{"rank"}, \code{"linreg"},
#'   \code{"mm"}, or \code{"score"}.
#' @param correct A logical indicating whether a continuity corrected
#'   statistic is used for rank correlation methods \code{"rank"} and
#'   \code{"count"}.
#' @param k.min Minimum number of studies to perform test for funnel
#'   plot asymmetry.
#' @param \dots Additional arguments (ignored at the moment).
#' 
#' @details
#' Functions to conduct rank correlation or linear regression tests
#' for funnel plot asymmetry.
#' 
#' Following recommendations by Sterne et al. (2011), by default, a
#' test for funnel plot asymmetry is only conducted if the number of
#' studies is ten or larger (argument \code{k.min = 10}). This
#' behaviour can be changed by setting a smaller value for argument
#' \code{k.min}. Note, the minimum number of studies is three.
#' 
#' If argument \code{method.bias} is \code{"rank"}, the test statistic
#' is based on the rank correlation between standardised treatment
#' estimates and variance estimates of estimated treatment effects;
#' Kendall's tau is used as correlation measure (Begg & Mazumdar,
#' 1994). The test statistic follows a standard normal
#' distribution. By default (if \code{correct} is FALSE), no
#' continuity correction is utilised (Kendall & Gibbons, 1990).
#' 
#' If argument \code{method.bias} is \code{"linreg"}, the test
#' statistic is based on a weighted linear regression of the treatment
#' effect on its standard error (Egger et al., 1997). The test
#' statistic follows a t distribution with \code{number of studies -
#' 2} degrees of freedom.
#' 
#' If argument \code{method.bias} is \code{"mm"}, the test statistic
#' is based on a weighted linear regression of the treatment effect on
#' its standard error using the method of moments estimator for the
#' additive between-study variance component (method 3a in Thompson,
#' Sharp, 1999). The test statistic follows a t distribution with
#' \code{number of studies - 2} degrees of freedom.
#' 
#' If argument \code{method.bias} is \code{"peters"}, the test
#' statistic is based on a weighted linear regression of the treatment
#' effect on the inverse of the total sample size using the variance
#' of the average event rate as weights (Peters et al., 2006). The
#' test statistic follows a t distribution with \code{number of
#' studies - 2} degrees of freedom. This test is available for
#' meta-analyses comparing two binary outcomes or combining single
#' proportions, i.e.  generated with functions \code{metabin} and
#' \code{metaprop}.
#' 
#' The following tests for funnel plot asymmetry are only available
#' for meta-analyses comparing two binary outcomes, i.e. meta-analyses
#' generated with the \code{metabin} function.
#' 
#' If argument \code{method.bias} is \code{"count"}, the test
#' statistic is based on the rank correlation between a standardised
#' cell frequency and the inverse of the variance of the cell
#' frequency; Kendall's tau is used as correlation measure (Schwarzer
#' et al., 2007). The test statistic follows a standard normal
#' distribution. By default (if \code{correct} is FALSE), no
#' continuity correction is utilised (Kendall & Gibbons, 1990).
#' 
#' If argument \code{method.bias} is \code{"score"}, the test
#' statistic is based on a weighted linear regression utilising
#' efficient score and score variance (Harbord et al., 2006,
#' 2009). The test statistic follows a t distribution with
#' \code{number of studies - 2} degrees of freedom.
#' 
#' In order to calculate an arcsine test for funnel plot asymmetry
#' (Rücker et al., 2008), one has to use the \code{metabin} function
#' with argument \code{sm = "ASD"} as input to the \code{metabias}
#' command. The three arcsine tests described in Rücker et al. (2008)
#' can be calculated by setting \code{method.bias} to \code{"rank"},
#' \code{"linreg"} and \code{"mm"}, respectively.
#' 
#' If argument \code{method.bias} is missing, the Harbord test
#' (\code{method.bias = "score"}) is used for the odds ratio as effect
#' measure and the Egger test (\code{method.bias = "linreg"}) for
#' other effect measures (Sterne et al., 2011).
#' 
#' No test for funnel plot asymmetry is conducted in meta-analyses
#' with subgroups.
#' @return
#' A list with class \code{htest} containing the following components
#' if a test for funnel plot asymmetry is conducted:
#' \item{estimate}{The estimated degree of funnel plot asymmetry, with
#'   name \code{"ks"} or \code{"bias"} corresponding to the method
#'   employed, i.e., rank correlation or regression method.}
#' \item{statistic}{The value of the test statistic.}
#' \item{parameters}{The degrees of freedom of the test statistic in
#'   the case that it follows a t distribution.}
#' \item{p.value}{The p-value for the test.}
#' \item{alternative}{A character string describing the alternative
#'   hypothesis.}
#' \item{method}{A character string indicating what type of test was
#'   used.}
#' \item{data.name}{A character string giving the names of the data.}
#'   \item{title}{Title of Cochrane review.}
#' \item{complab}{Comparison label.}
#' \item{outclab}{Outcome label.}
#' \item{version}{Version of R package \bold{meta} used to create
#'   object.}
#'
#' Or a list with the following elements if test is not conducted due
#' to the number of studies:
#' \item{k}{Number of studies in meta-analysis.}
#' \item{k.min}{Minimum number of studies to perform test for funnel
#'   plot asymmetry.}
#' \item{version}{Version of R package \bold{meta} used to create
#'   object.}
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{funnel}}, \code{\link{funnel.meta}},
#'   \code{\link{metabin}}, \code{\link{metacont}},
#'   \code{\link{metagen}}
#' 
#' @references
#' Begg CB & Mazumdar M (1994):
#' Operating characteristics of a rank correlation test for
#' publication bias.
#' \emph{Biometrics},
#' \bold{50}, 1088--101
#' 
#' Egger M, Smith GD, Schneider M & Minder C (1997):
#' Bias in meta-analysis detected by a simple, graphical test.
#' \emph{British Medical Journal},
#' \bold{315}, 629--34
#' 
#' Harbord RM, Egger M & Sterne J (2006):
#' A modified test for small-study effects in meta-analyses of
#' controlled trials with binary endpoints.
#' \emph{Statistics in Medicine},
#' \bold{25}, 3443--57
#'
#' Harbord RM, Harris RJ, Sterne JAC (2009):
#' Updated tests for small-study effects in meta–analyses.
#' \emph{The Stata Journal},
#' \bold{9}, 197--210
#' 
#' Kendall M & Gibbons JD (1990):
#' \emph{Rank Correlation Methods}.
#' London: Edward Arnold
#' 
#' Peters JL, Sutton AJ, Jones DR, Abrams KR & Rushton L (2006):
#' Comparison of two methods to detect publication bias in
#' meta-analysis.
#' \emph{Journal of the American Medical Association},
#' \bold{295}, 676--80
#' 
#' Rücker G, Schwarzer G, Carpenter JR (2008):
#' Arcsine test for publication bias in meta-analyses with binary
#' outcomes.
#' \emph{Statistics in Medicine},
#' \bold{27}, 746--63
#' 
#' Schwarzer G, Antes G & Schumacher M (2007):
#' A test for publication bias in meta-analysis with sparse binary
#' data.
#' \emph{Statistics in Medicine},
#' \bold{26}, 721--33
#' 
#' Sterne, JAC et al. (2011):
#' Recommendations for examining and interpreting funnel plot
#' asymmetry in meta-analyses of randomised controlled trials.
#' \emph{BMJ (Clinical research ed.)},
#' \bold{343}, 1
#' 
#' Thompson SG & Sharp, SJ (1999):
#' Explaining heterogeneity in meta-analysis: a comparison of methods,
#' \emph{Statistics in Medicine},
#' \bold{18}, 2693--708
#'
#' @keywords htest
#' 
#' @examples
#' data(Olkin95)
#' m1 <- metabin(event.e, n.e, event.c, n.c,
#'               data = Olkin95, subset = 1:10,
#'               sm = "RR", method = "I")
#' 
#' metabias(m1)
#' metabias(m1, plotit = TRUE)
#' 
#' metabias(m1, method.bias = "rank")
#' metabias(m1, method.bias = "rank", correct = TRUE)
#' 
#' metabias(m1, method.bias = "count")
#' metabias(m1, method.bias = "linreg")$p.value
#' 
#' # Arcsine test (based on linear regression)
#' #
#' m1.as <- update(m1, sm = "ASD")
#' metabias(m1.as)
#' # Same result (using function metabias.default)
#' metabias(m1.as$TE, m1.as$seTE)
#' 
#' # No test for funnel plot asymmetry calculated
#' #
#' m2 <- update(m1, subset = 1:5)
#' metabias(m2)
#' 
#' m3 <- update(m1, subset = 1:2)
#' metabias(m3)
#' 
#' # Test for funnel plot asymmetry calculated (use of argument k.min)
#' #
#' metabias(m2, k.min = 5)
#' 
#' @rdname metabias
#' @export metabias


metabias <- function(x, ...) 
  UseMethod("metabias")





#' @rdname metabias
#' @method metabias default
#' @export
#' @export metabias.default


metabias.default <- function(x, seTE,
                             method.bias = "linreg",
                             plotit = FALSE, correct = FALSE,
                             k.min = 10, ...) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  k.All <- length(x)
  ##
  chknumeric(x)
  chknumeric(seTE)
  ##
  fun <- "metabias"
  chklength(seTE, k.All, fun)
  ##
  method.bias <- setchar(method.bias, c("rank", "linreg", "mm"))
  
  
  ##
  ##
  ## (2) Do meta-analysis
  ##
  ##
  m <- metagen(x, seTE)
  
  
  ##
  ##
  ## (3) Conduct test for funnel plot asymmetry
  ##
  ##
  res <- metabias(m, method.bias = method.bias,
                  plotit = plotit, correct = correct,
                  k.min = k.min, ...)
  
  
  res
}





#' @rdname metabias
#' @method metabias meta
#' @export
#' @export metabias.meta


metabias.meta <- function(x, method.bias = x$method.bias,
                          plotit = FALSE, correct = FALSE,
                          k.min = 10, ...) {
  
  ##
  ##
  ## (1) Check for meta object
  ##
  ##
  
  chkclass(x, "meta")
  x.name <- deparse(substitute(x))
  ##  
  if (inherits(x, "metacum"))
    stop("Test for funnel plot asymmetry not meaningful for object of class \"metacum\".")
  ##
  if (inherits(x, "metainf"))
    stop("Test for funnel plot asymmetry not meaningful for object of class \"metainf\".")
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  ## --
  ##
  ## Set default for method.bias (see Sterne et al., BMJ 2011;343:d4002):
  ## - "score" for OR as effect measure
  ## - "linreg" otherwise
  ##
  if (length(method.bias) == 0)
    if (inherits(x, "metabin") & x$sm == "OR")
      method.bias <- "score"
    else
      method.bias <- "linreg"
  ##
  tests <- c("rank", "linreg", "mm", "count", "score", "peters")
  method.bias <- setchar(method.bias, tests)
  imeth <- charmatch(method.bias, tests)
  method <- c(paste("Rank correlation test of funnel plot asymmetry",
                    ifelse(correct == TRUE, " (with continuity correction)", ""),
                    sep = ""),
              "Linear regression test of funnel plot asymmetry",
              "Linear regression test of funnel plot asymmetry (methods of moment)",
              paste("Rank correlation test of funnel plot asymmetry (based on counts)",
                    ifelse(correct == TRUE, " (with continuity correction)", ""),
                    sep = ""),
              "Linear regression test of funnel plot asymmetry (efficient score)",
              "Linear regression test of funnel plot asymmetry (based on sample size)")[imeth]    ##
  chklogical(plotit)
  chklogical(correct)
  chknumeric(k.min, 1, single = TRUE)
  
  
  TE <- x$TE
  seTE <- x$seTE
  n.e <- x$n.e
  n.c <- x$n.c
  ##
  if (inherits(x, "metabin")) {
    event.e <- x$event.e
    event.c <- x$event.c
    ##
    if (x$method == "Peto") {
      warning(paste("Inverse variance method used for pooling",
                    "to perform test of funnel plot asymmetry."))
      m <- metabin(event.e, n.e, event.c, n.c, sm = x$sm, method = "Inverse",
                   incr = x$incr, allincr = x$allincr,
                   allstudies = x$allstudies, MH.exact = x$MH.exact,
                   RR.cochrane = x$RR.cochrane)
      TE <- m$TE
      seTE <- m$seTE
    }
  }
  ##
  if (inherits(x, "metaprop")) {
    event <- x$event
    n <- x$n
  }
  
  if(length(TE) != length(seTE))
    stop("length of argument TE and seTE must be equal")
  ##
  ## Exclude studies from meta-analysis
  ##
  if (!is.null(x$exclude)) {
    TE <- TE[!x$exclude]
    seTE <- seTE[!x$exclude]
    ##
    if (!is.null(n.e))
      n.e <- n.e[!x$exclude]
    if (!is.null(n.c))
      n.c <- n.c[!x$exclude]
    ##
    if (inherits(x, "metabin")) {
      event.e <- event.e[!x$exclude]
      event.c <- event.c[!x$exclude]
    }
    ##
    if (inherits(x, "metaprop")) {
      n <- n[!x$exclude]
      event <- event[!x$exclude]
    }
  }
  ##
  sel <- !is.na(TE) & !is.na(seTE)
  if (length(TE) != sum(sel))
    warning(paste(length(TE) - sum(sel),
                  "observation(s) dropped due to missing values"))
  ##
  TE <- TE[sel]
  seTE <- seTE[sel]
  ##
  if (inherits(x, "metabin")) {
    n.e <- n.e[sel]
    n.c <- n.c[sel]
    event.e <- event.e[sel]
    event.c <- event.c[sel]
  }
  ##
  k <- length(TE)
  
  
  if (k < k.min | k < 3)
    res <- list(k = k, k.min = k.min)
  else if (length(x$byvar) != 0)
    res <- list(subgroup = TRUE)
  else {
    if (method.bias == "rank") {
      if (length(unique(seTE)) == 1)
        stop("Test for small-study effects not feasible as all studies have same precision")
      ##
      ## Begg und Mazumdar (1994), Biometrics, 50, 1088-1101
      ##
      m <- metagen(TE, seTE)
      TE.fixed <- m$TE.fixed
      seTE.fixed <- m$seTE.fixed
      ##
      varTE.s <- seTE^2 - seTE.fixed^2
      TE.s <- (TE - TE.fixed) / sqrt(varTE.s)
      ##
      ktau <- kentau(TE.s, seTE^2, correct = correct)
      ##
      res <- list(estimate = c(ktau$ks, ktau$se.ks),
                  statistic = ktau$ks / ktau$se.ks,
                  p.value = ktau$p.value)
      names(res$statistic) <- "z"
      names(res$estimate) <- c("ks", "se.ks")
    }
    else if (method.bias == "linreg" | method.bias == "mm" | method.bias == "score" |
             method.bias == "peters") {
      
      if (method.bias == "linreg") {
        if (length(unique(seTE)) == 1)
          stop("Test for small-study effects not feasible as all studies have same precision")
        ##
        ## Egger, Smith, Schneider, Minder (1997), BMJ, 315, 629-34
        ##
        lreg <- linregcore(seTE, TE, 1 / seTE^2)
        se.bias <- lreg$se.slope
      }
      else if (method.bias == "mm") {
        if (length(unique(seTE)) == 1)
          stop("Test for small-study effects not feasible as all studies have same precision")
        ##
        ## Thompson und Sharp (1999), Stat Med, 18, 2693-2708
        ##
        fit1 <- linregcore(1 / seTE, TE / seTE)
        Q <- sum((TE / seTE - fit1$intercept - fit1$slope / seTE)^2)
        ##
        x <- seTE
        y <- TE
        w <- 1 / seTE^2
        ##
        tau2 <- (Q - (k - 2)) / (sum(w) - (sum(w^2) * sum(w * x^2) -
                                             2 * sum(w^2 * x) * sum(w * x) +
                                               sum(w) * sum(w^2 * x^2)) /
                                                 (sum(w) * sum(w * x^2) - (sum(w * x))^2))
        ##
        tau2 <- ifelse(tau2 < 0, 0, tau2)
        ##
        lreg <- linregcore(seTE, TE, 1 / (seTE^2 + tau2))
        se.bias <- lreg$se.slope / sqrt(lreg$MSE.w)
      }
      else if (method.bias == "score") {
        if (inherits(x, "metabin")) {
          if (x$sm == "RR") {
            ##
            ## Harbord et al. (2009), The Stata Journal
            ##
            Z <- (event.e * (n.e + n.c) - (event.e + event.c) * n.e) /
              (n.e - event.e + n.c - event.c)
            V <- n.e * n.c * (event.e + event.c) /  
              ((n.e + n.c) * (n.e - event.e + n.c - event.c))
          }
          else {
            if (x$sm != "OR")
              warning("Using odds ratio as effect measure in Harbord test.")
            ##
            ## Harbord et al. (2006), Statistics in Medicine
            ##
            Z <- event.e - (event.e + event.c) * n.e / (n.e + n.c)
            V <- n.e * n.c * (event.e + event.c) * 
              (n.e - event.e + n.c - event.c) /
              ((n.e + n.c)^2 * ((n.e + n.c) - 1))
          }
          ##
          TE.score <- Z / V
          seTE.score <- 1 / sqrt(V)
          ##
          lreg <- linregcore(seTE.score, TE.score, 1 / seTE.score^2)
          se.bias <- lreg$se.slope
        }
        else {
          stop("method.bias '", method.bias,
               "' only defined for meta-analysis ",
               "with binary outcome data (function 'metabin')")
        }
      }
      else if (method.bias == "peters") {
        ##
        ## Peters et al. (2006), JAMA
        ##
        if (inherits(x, "metabin")) {
          seTE.peters <- sqrt(1 / (event.e + event.c) +
                                1 / (n.e - event.e + n.c - event.c))
          ##
          lreg <- linregcore(1 / (n.e + n.c), TE, 1 / seTE.peters^2)
          se.bias <- lreg$se.slope
        }
        else if (inherits(x, "metaprop")) {
          seTE.peters <- sqrt(1 / event + 1 / (n - event))
          ##
          lreg <- linregcore(1 / n, TE, 1 / seTE.peters^2)
          se.bias <- lreg$se.slope
        }
        else {
          stop(paste("method.bias '", method.bias,
                     "' only defined for meta-analysis conducted with metabin() or metaprop()", sep = ""))
        }
      }
      ##
      ##
      bias <- lreg$slope
      df <- lreg$df
      slope <- lreg$intercept
      ##
      statistic <- bias / se.bias
      p.value <- 2 * pt(abs(statistic), df = df, lower.tail = FALSE)[[1]]
      ##
      res <- list(estimate = c(bias, se.bias, slope),
                  parameters = df, statistic = statistic,
                  p.value = p.value)
      
      names(res$statistic) <- "t"
      names(res$parameters) <- "df"
      names(res$estimate) <- c("bias", "se.bias", "slope")
    }
    else if (method.bias == "count") {
      ##
      ## Schwarzer, Antes, Schumacher (2006), Stat Med
      ##
      if (inherits(x, "metabin")) {
        TE.MH <- metabin(x$event.e, x$n.e,
                         x$event.c, x$n.c,
                         sm = "OR", method = "MH", warn = FALSE)$TE.fixed
        
        n.. <- n.e + n.c
        n11 <- event.e
        n1. <- n.e
        n.1 <- event.e + event.c
        
        n11.e <- rep(NA, k)
        n11.v <- rep(NA, k)
        ##
        for ( i in seq(1:k) ) {
          obj <- hypergeometric(n1.[i], n.1[i], n..[i], exp(TE.MH))
          n11.e[i] <- obj$mean()
          n11.v[i] <- obj$var()
        }
        
        ktau <- kentau((n11 - n11.e) / sqrt(n11.v), 1 / n11.v, correct = correct)
        ##
        res <- list(estimate = c(ktau$ks, ktau$se.ks),
                    statistic = ktau$ks / ktau$se.ks,
                    p.value = ktau$p.value)
        names(res$statistic) <- "z"
        names(res$estimate) <- c("ks", "se.ks")
      }
      else {
        stop(paste("method.bias '", method.bias, "' only defined for meta-analysis with binary outcome data (function 'metabin')", sep = ""))
      }
    }
    
    res$alternative <- "asymmetry in funnel plot"
    
    res$method <- method
    
    res$data.name <- x.name
    
    if (plotit) {
      ##
      if (method.bias == "linreg" | method.bias == "mm") {
        radial(TE, seTE, comb.fixed = FALSE)
        abline(lreg$slope, lreg$intercept)
      }
      else if (method.bias == "rank") {
        ##
        if (plotit) {
          plot(TE.s, seTE^2,
               xlab = "Standardised treatment effect",
               ylab = "Estimated variance of treatment estimate")
        }
      }
      else if (method.bias == "score") {
        ##
        radial(TE.score, seTE.score, comb.fixed = FALSE)
        abline(lreg$slope, lreg$intercept)
      }
    }

    if (inherits(x, "meta")) {
      if (!is.null(x$title))
        res$title <- x$title
      if (!is.null(x$complab))
        res$complab <- x$complab
      if (!is.null(x$outclab))
        res$outclab <- x$outclab
    }
  }
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metabias", "htest")
  
  res
}





#' @rdname metabias
#' @method print metabias
#' @export
#' @export print.metabias


print.metabias <- function(x, ...) {
  
  
  ##
  ##
  ## (1) Check for metabias object
  ##
  ##
  chkclass(x, "metabias")
  
  
  ## Print information for meta-analysis from Cochrane Review
  ##
  crtitle(x)
  
  class(x) <- "htest"
  ##
  if (length(x$p.value) != 0)
    print(x)
  else {
    ##
    ## Check whether number of studies is too small:
    ##
    if (length(x$k) != 0 & length(x$k.min != 0)) {
      if (x$k <= x$k.min) {
        if (x$k <= 2)
          warning("Number of studies (k=",  x$k,
                  ") too small to test for small study effects.")
        else
          warning("Number of studies (k=",  x$k,
                  ") too small to test for small study effects (k.min=",
                  x$k.min, "). Change argument 'k.min' if appropriate.\n")
      }
    }
    ##
    ## Check whether meta-analysis has subgroups:
    ##
    if (length(x$subgroup) != 0)
      warning("No test for small study effects conducted ",
              "for meta-analysis with subgroups.")
  }
  
  invisible(NULL)
}
