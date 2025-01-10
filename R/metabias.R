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
#'   be used (see Details), can be abbreviated.
#' @param plotit A logical indicating whether a plot should be
#'   produced (see Details).
#' @param correct A logical indicating whether a continuity corrected
#'   statistic is used for rank correlation tests.
#' @param k.min Minimum number of studies to perform test for funnel
#'   plot asymmetry.
#' @param digits Minimal number of significant digits for estimates,
#'   see \code{print.default}.
#' @param digits.stat Minimal number of significant digits for z- or
#'   t-value of test for test of funnel plot asymmetry, see
#'   \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of test for test of funnel plot asymmetry, see
#'   \code{print.default}.
#' @param digits.se Minimal number of significant digits for standard
#'   errors, see \code{print.default}.
#' @param digits.tau2 Minimal number of significant digits for
#'   residual heterogeneity variance, see \code{print.default}.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param big.mark A character used as thousands separator.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   overall effect should be printed according to JAMA reporting
#'   standards.
#' @param text.tau2 Text printed to identify residual heterogeneity
#'   variance \eqn{\tau^2}.
#' @param details.methods A logical specifying whether details on
#'   statistical methods should be printed.
#' @param \dots Additional arguments passed on to
#'   \code{\link[metafor]{rma.uni}}.
#' 
#' @details
#' Functions to conduct rank correlation or linear regression tests
#' for funnel plot asymmetry.
#'
#' \subsection{Classic generic tests}{
#' The following tests are generic tests for funnel plot asymmetry
#' which only require estimates of the treatment effect and
#' corresponding standard errors. Accordingly, these are the only
#' tests provided by R function \code{metabias.default}.
#' 
#' If argument \code{method.bias} is \code{"Begg"}, the test statistic
#' is based on the rank correlation between standardised treatment
#' estimates and variance estimates of estimated treatment effects;
#' Kendall's tau is used as correlation measure (Begg & Mazumdar,
#' 1994). The test statistic follows a standard normal
#' distribution. By default (if \code{correct} is FALSE), no
#' continuity correction is utilised (Kendall & Gibbons, 1990).
#' 
#' If argument \code{method.bias} is \code{"Egger"}, the test
#' statistic is based on a weighted linear regression of the treatment
#' effect on its standard error (Egger et al., 1997). The test
#' statistic follows a t distribution with \code{number of studies -
#' 2} degrees of freedom.
#' 
#' If argument \code{method.bias} is \code{"Thompson"}, the test
#' statistic is based on a weighted linear regression of the treatment
#' effect on its standard error using an additive between-study
#' variance component denoted as methods (3a) - (3d) in Thompson &
#' Sharp (1999). The test statistic follows a t distribution with
#' \code{number of studies - 2} degrees of freedom.  }
#' 
#' \subsection{Tests for meta-analysis with binary outcomes}{
#' The following tests for funnel plot asymmetry are only available
#' for meta-analyses comparing two binary outcomes, i.e. meta-analyses
#' generated with the \code{metabin} function. The only exception is
#' the test by Peters et al. (2006) which can also be used in a
#' meta-analysis of single proportions generated with \code{metaprop}.
#' 
#' If argument \code{method.bias} is \code{"Harbord"}, the test
#' statistic is based on a weighted linear regression utilising
#' efficient score and score variance (Harbord et al., 2006,
#' 2009). The test statistic follows a t distribution with
#' \code{number of studies - 2} degrees of freedom.
#' 
#' In order to calculate an arcsine test for funnel plot asymmetry
#' (Rücker et al., 2008), one has to use the \code{metabin} function
#' with argument \code{sm = "ASD"} as input to the \code{metabias}
#' command. The three arcsine tests described in Rücker et al. (2008)
#' can be calculated by setting \code{method.bias} to \code{"Begg"},
#' \code{"Egger"} and \code{"Thompson"}, respectively.
#' 
#' If argument \code{method.bias} is \code{"Macaskill"}, the test
#' statistic is based on a weighted linear regression of the treatment
#' effect on the total sample size with weights reciprocal to the
#' variance of the average event probability (Macaskill et al., 2001,
#' \emph{method FPV}). The test statistic follows a t distribution
#' with \code{number of studies - 2} degrees of freedom.
#' 
#' If argument \code{method.bias} is \code{"Peters"}, the test
#' statistic is based on a weighted linear regression of the treatment
#' effect on the inverse of the total sample size with weights
#' reciprocal to the variance of the average event probability (Peters
#' et al., 2006). The test statistic follows a t distribution with
#' \code{number of studies - 2} degrees of freedom. Note, this test is
#' a variant of Macaskill et al. (2001), \emph{method FPV}, using the
#' inverse sample size as covariate.
#' 
#' If argument \code{method.bias} is \code{"Schwarzer"}, the test
#' statistic is based on the rank correlation between a standardised
#' cell frequency and the inverse of the variance of the cell
#' frequency; Kendall's tau is used as correlation measure (Schwarzer
#' et al., 2007). The test statistic follows a standard normal
#' distribution. By default (if \code{correct} is FALSE), no
#' continuity correction is utilised (Kendall & Gibbons, 1990).
#' 
#' Finally, for meta-analysis of diagnostic test accuracy studies, if
#' argument \code{method.bias} is \code{"Deeks"}, the test statistic
#' is based on a weighted linear regression of the log diagnostic odds
#' ratio on the inverse of the squared effective sample size using the
#' effective sample size as weights (Deeks et al., 2005). The test
#' statistic follows a t distribution with \code{number of studies -
#' 2} degrees of freedom.
#' }
#'
#' \subsection{Test for the standardised mean difference}{
#' If argument \code{method.bias} is \code{"Pustejovsky"}, the test
#' statistic is based on a weighted linear regression of the treatment
#' effect on the square root of the sum of the inverse group sample
#' sizes using the treatment effect variance as weights (Pustejovsky &
#' Rodgers, 2019). The test statistic follows a t distribution with
#' \code{number of studies - 2} degrees of freedom.
#' }
#' 
#' \subsection{Recommendations and default settings}{
#' Following recommendations by Sterne et al. (2011), by default, a
#' test for funnel plot asymmetry is only conducted if the number of
#' studies is ten or larger (argument \code{k.min = 10}). This
#' behaviour can be changed by setting a smaller value for argument
#' \code{k.min}. Note, the minimum number of studies is three.
#' 
#' If argument \code{method.bias} is missing, the Harbord test
#' (\code{method.bias = "Harbord"}) is used in meta-analyses with a binary
#' outcome for the odds ratio and Deeks' test (\code{method.bias = "Deeks"})
#' for the diagnostic odds ratios. In all other settings, the Egger test
#' (\code{method.bias = "Egger"}) is used (Sterne et al., 2011).
#' 
#' No test for funnel plot asymmetry is conducted in meta-analyses
#' with subgroups.
#' }
#'
#' If argument \code{plotit = TRUE}, a scatter plot is shown if
#' argument \code{method.bias} is equal to \code{"Begg"},
#' \code{"Egger"}, \code{"Thompson"}, \code{"Harbord"}, or
#' \code{"Deeks"}.
#' 
#' @return
#' A list with class \code{metabias} containing the following
#' components if a test for funnel plot asymmetry is conducted:
#' \item{statistic}{Test statistic.}
#' \item{df}{The degrees of freedom of the test statistic in
#'   the case that it follows a t distribution.}
#' \item{pval}{The p-value for the test.}
#' \item{estimate}{Estimates used to calculate test statisic.}
#' \item{method}{A character string indicating what type of test was
#'   used.}
#' \item{title}{Title of Cochrane review.}
#' \item{complab}{Comparison label.}
#' \item{outclab}{Outcome label.}
#' \item{var.model}{A character string indicating whether none,
#'   multiplicative, or additive residual heterogeneity variance was
#'   assumed.}
#' \item{method.bias}{As defined above.}
#' \item{x}{Meta-analysis object.}
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
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
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
#' Deeks JJ, Macaskill P, Irwig L (2005):
#' The performance of tests of publication bias and other sample size
#' effects in systematic reviews of diagnostic test accuracy was
#' assessed.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{58}:882--93
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
#' Macaskill P, Walter SD, Irwig L (2001):
#' A comparison of methods to detect publication bias in
#' meta-analysis.
#' \emph{Statistics in Medicine},
#' \bold{20}, 641--54
#' 
#' Peters JL, Sutton AJ, Jones DR, Abrams KR & Rushton L (2006):
#' Comparison of two methods to detect publication bias in
#' meta-analysis.
#' \emph{Journal of the American Medical Association},
#' \bold{295}, 676--80
#'
#' Pustejovsky JE, Rodgers MA (2019):
#' Testing for funnel plot asymmetry of standardized mean differences.
#' \emph{Research Synthesis Methods},
#' \bold{10}, 57--71
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
#' @examples
#' data(Olkin1995)
#' m1 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
#'   data = Olkin1995, subset = 1:10, sm = "RR", method = "I")
#' 
#' metabias(m1)
#' metabias(m1, plotit = TRUE)
#' 
#' metabias(m1, method.bias = "Begg")
#' metabias(m1, method.bias = "Begg", correct = TRUE)
#' 
#' metabias(m1, method.bias = "Schwarzer")
#' metabias(m1, method.bias = "Egger")$pval
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
#' @method metabias meta
#' @export


metabias.meta <- function(x, method.bias = x$method.bias,
                          plotit = FALSE, correct = FALSE,
                          k.min = 10, ...) {
  
  ##
  ##
  ## (1) Check for meta object
  ##
  ##
  chkclass(x, "meta")
  chksuitable(x, "Test for funnel plot asymmetry", "metamerge")
  ##
  x <- updateversion(x)
  x.name <- deparse(substitute(x))
  
  
  ##
  ##
  ## (2) Check / set other arguments
  ##
  ##
  if (length(method.bias) == 0) {
    ## Set default for method.bias (see Sterne et al., BMJ 2011;343:d4002):
    ## - "Harbord" for OR as effect measure
    ## - "Egger" otherwise
    if (inherits(x, "metabin") & x$sm == "OR")
      method.bias <- "Harbord"
    else
      method.bias <- "Egger"
  }
  ##
  method.bias <- setmethodbias(method.bias)
  ##
  if (method.bias %in% c("Begg", "Schwarzer"))
    lab.method <-
      paste0("Rank correlation test of funnel plot asymmetry",
             if (correct == TRUE) " (with continuity correction)")
  else if (method.bias == "Deeks")
    lab.method <- "Funnel plot test for diagnostic odds ratios"
  else
    lab.method <- "Linear regression test of funnel plot asymmetry"
  ##
  chklogical(plotit)
  chklogical(correct)
  chknumeric(k.min, 1, length = 1)
  
  
  ##
  ##
  ## (3) Select studies for inclusion in test of funnel plot asymmetry
  ##
  ##
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
      warning("Inverse variance method used for pooling ",
              "to perform test of funnel plot asymmetry.",
              call = FALSE)
      m <- metabin(event.e, n.e, event.c, n.c, sm = x$sm, method = "Inverse",
                   incr = x$incr, allincr = x$allincr,
                   allstudies = x$allstudies, MH.exact = x$MH.exact,
                   RR.Cochrane = x$RR.Cochrane)
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
    stop("Length of argument TE and seTE must be equal.",
         call. = FALSE)
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
  if (inherits(x, c("metabin", "metacont"))) {
    n.e <- n.e[sel]
    n.c <- n.c[sel]
  }
  ##
  if (inherits(x, "metabin")) {
    event.e <- event.e[sel]
    event.c <- event.c[sel]
  }
  ##
  k <- length(TE)
  
  
  ##
  ##
  ## (4) Conduct test of funnel plot asymmetry
  ##
  ##
  if (k < k.min | k < 3)
    res <- list(k = k, k.min = k.min)
  else if (length(x$subgroup) != 0)
    res <- list(subgroup = TRUE)
  else {
    ##
    ## Check whether standard errors differ
    ##
    if (method.bias %in% c("Begg", "Egger", "Thompson",
                           "Macaskill", "Pustejovsky"))
      if (length(unique(seTE)) == 1)
        stop("Test for small-study effects not feasible as ",
             "all studies have the same precision.",
             call. = FALSE)
    ##
    ## (a) Rank correlation tests
    ##
    if (method.bias %in% c("Begg", "Schwarzer")) {
      if (method.bias == "Begg") {
        ##
        ## Begg und Mazumdar (1994), Biometrics, 50, 1088-1101
        ##
        m <- metagen(TE, seTE, method.tau = "DL", method.tau.ci = "")
        TE.common <- m$TE.common
        seTE.common <- m$seTE.common
        ##
        varTE.s <- seTE^2 - seTE.common^2
        TE.s <- (TE - TE.common) / sqrt(varTE.s)
        ##
        ktau <- kentau(TE.s, seTE^2, correct = correct)
      }
      else {
        ##
        ## Schwarzer, Antes, Schumacher (2006), Stat Med
        ##
        if (inherits(x, "metabin")) {
          TE.MH <- metabin(x$event.e, x$n.e, x$event.c, x$n.c,
                           sm = "OR", method = "MH",
                           method.tau = "DL", method.tau.ci = "",
                           warn = FALSE)$TE.common
          ##
          n.. <- n.e + n.c
          n11 <- event.e
          n1. <- n.e
          n.1 <- event.e + event.c
          ##
          n11.e <- rep(NA, k)
          n11.v <- rep(NA, k)
          ##
          for ( i in seq(1:k) ) {
            obj <- hypergeometric(n1.[i], n.1[i], n..[i], exp(TE.MH))
            n11.e[i] <- obj$mean()
            n11.v[i] <- obj$var()
          }
          ##
          ktau <- kentau((n11 - n11.e) / sqrt(n11.v),
                         1 / n11.v, correct = correct)
        }
      }
      ##
      res <- list(statistic = ktau$ks / ktau$se.ks,
                  pval = ktau$p.value,
                  estimate = c(ktau$ks, ktau$se.ks),
                  p.value = ktau$p.value)
      ##
      names(res$statistic) <- "z"
      names(res$estimate) <- c("ks", "se.ks")
    }
    ##
    ## (b) Linear regression tests
    ##
    else {
      if (method.bias == "Egger") {
        ##
        ## Egger, Smith, Schneider, Minder (1997), BMJ, 315, 629-34
        ##
        lreg <- linregcore(TE, seTE, ...)
      }
      else if (method.bias == "Thompson") {
        ##
        ## Thompson und Sharp (1999), Stat Med, 18, 2693-2708
        ##
        lreg <- linregcore(TE, seTE, model = "rma", method.tau = x$method.tau,
                           ...)
      }
      else if (method.bias == "Harbord") {
        if (inherits(x, "metabin")) {
          if (x$sm == "RR") {
            ##
            ## Harbord et al. (2009), The Stata Journal
            ##
            Z <- (event.e * (n.e + n.c) - (event.e + event.c) * n.e) /
              (n.e - event.e + n.c - event.c)
            V <- n.e * n.c * (event.e + event.c) /  
              ((n.e + n.c) * (n.e - event.e + n.c - event.c))
            V <- (n.e + n.c) / 4 *
              (event.e + event.c) / (n.e - event.e + n.c - event.c)
          }
          else {
            if (x$sm != "OR")
              warning("Using odds ratio as effect measure in Harbord test.",
                      call. = FALSE)
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
          lreg <- linregcore(TE.score, seTE.score, ...)
        }
        else
          stoponly("method.bias", method.bias, "metabin()")
      }
      else if (method.bias == "Macaskill") {
        ##
        ## Macaskill et al. (2001), Stat Med
        ##
        if (inherits(x, "metabin") || inherits(x, "metacont"))
          covar <- n.e + n.c
        else if (length(n.e) > 0 & length(n.c) > 0)
          covar <- n.e + n.c
        else if (inherits(x, "metaprop"))
          covar <- n
        else if (length(n.e) > 0)
          covar <- n.e
        else if (length(n.c) > 0)
          covar <- n.c
        else
          stop("No information on sample size available.")
        ##
        lreg <- linregcore(TE, seTE, covar, ...)
      }
      else if (method.bias == "Peters") {
        ##
        ## Peters et al. (2006), JAMA
        ##
        if (inherits(x, "metabin")) {
          seTE.Peters <- sqrt(1 / (event.e + event.c) +
                              1 / (n.e - event.e + n.c - event.c))
          covar <- 1 / (n.e + n.c)
        }
        else if (inherits(x, "metaprop")) {
          seTE.Peters <- sqrt(1 / event + 1 / (n - event))
          covar <- 1 / n
        }
        else
          stoponly("method.bias", method.bias, "metabin() or metaprop()")
        ##
        lreg <- linregcore(TE, seTE.Peters, covar, ...)
      }
      else if (method.bias == "Deeks") {
        ##
        ## Deeks et al. (2005), J Clin Epid
        ##
        if (inherits(x, "metabin")) {
          ESS <- 4 * n.e * n.c / (n.e + n.c)
          ##
          lreg <- linregcore(TE, sqrt(1 / ESS), 1 / sqrt(ESS), ...)
        }
        else
          stoponly("method.bias", method.bias, "metabin()")
      }
      else if (method.bias == "Pustejovsky") {
        ##
        ## Pustejovsky & Rodgers (2019)
        ##
        if (is.null(n.e))
          stop("Sample size in experimental group (n.e) missing.",
               call. = FALSE)
        if (is.null(n.c))
          stop("Sample size in control group (n.c) missing.",
               call. = FALSE)
        ##
        lreg <- linregcore(TE, seTE, sqrt(1 / n.e + 1 / n.c), ...)
      }
      ##
      res <- list(statistic = lreg$statistic, df = lreg$df, pval = lreg$pval,
                  estimate = c(lreg$slope, lreg$se.slope),
                  tau = lreg$tau, p.value = lreg$pval,
                  intercept = lreg$intercept,
                  se.intercept = lreg$se.intercept)
      names(res$estimate) <- c("bias", "se.bias")
    }
    
    
    ##
    ## Optional plot
    ##
    if (plotit) {
      ##
      if (method.bias == "Egger" | method.bias == "Thompson") {
        radial(TE, seTE, common = FALSE)
        abline(lreg$slope, lreg$intercept)
      }
      else if (method.bias == "Begg") {
        ##
        if (plotit) {
          plot(TE.s, seTE^2,
               xlab = "Standardised treatment effect",
               ylab = "Estimated variance of treatment estimate")
        }
      }
      else if (method.bias == "Harbord") {
        ##
        radial(TE.score, seTE.score, common = FALSE)
        abline(lreg$slope, lreg$intercept)
      }
    }
    
    
    res$method <- lab.method
    ##
    res$var.model <- ifelse(method.bias %in% c("Begg", "Schwarzer"),
                            "",
                     ifelse(method.bias == "Thompson",
                            "additive", "multiplicative"))
    ##
    res$method.bias <- method.bias
    
    
    if (inherits(x, "meta")) {
      if (!is.null(x$title))
        res$title <- x$title
      if (!is.null(x$complab))
        res$complab <- x$complab
      if (!is.null(x$outclab))
        res$outclab <- x$outclab
    }
  }
  
  
  res$x <- x
  res$version <- packageDescription("meta")$Version
  
  class(res) <- "metabias"
  
  res
}





#' @rdname metabias
#' @method print metabias
#' @export


print.metabias <- function(x,
                           digits = gs("digits"),
                           digits.stat = gs("digits.stat"),
                           digits.pval = max(gs("digits.pval"), 2),
                           digits.se = gs("digits.se"),
                           digits.tau2 = gs("digits.tau2"),
                           #
                           scientific.pval = gs("scientific.pval"),
                           big.mark = gs("big.mark"),
                           zero.pval = gs("zero.pval"),
                           JAMA.pval = gs("JAMA.pval"),
                           #
                           text.tau2 = gs("text.tau2"),
                           #
                           details.methods = gs("details"),
                           #
                           ...) {
  
  
  #
  #
  # (1) Check for metabias object
  #
  #
  
  chkclass(x, "metabias")
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  #
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  #
  chkchar(text.tau2, length = 1)
  chklogical(details.methods)
  
  # Update old metabias objects
  #
  if (inherits(x, "htest")) {
    x$pval <- x$p.value
    x$x <- list(method.tau = "DL")
    x$df <- x$parameters
    x$var.model <- "multiplicative"
    #
    if (grepl("Rank correlation test of funnel plot asymmetry",
              x$method)) {
      if (!grepl("counts", x$method))
        x$method.bias <- "Begg"
      else
        x$method.bias <- "Schwarzer"
      #
      x$method <- "Rank correlation test of funnel plot asymmetry"
      x$var.model <- ""
    }
    else if (grepl("Deeks", x$method)) {
      x$method.bias <- "Deeks"
      x$method <- "Funnel plot test for diagnostic odds ratios"
    }
    else {
      #
      if (grepl("score", x$method))
        x$method.bias <- "Harbord"
      else if (grepl("sample size", x$method))
        x$method.bias <- "Peters"
      else if (grepl("methods of moment", x$method)) {
        x$method.bias <- "Thompson"
      x$var.model <- "additive"
      }
      else
        x$method.bias <- "Egger"
      #
      x$method <- "Linear regression test of funnel plot asymmetry"
    }
  }
  
  # Print information for meta-analysis from Cochrane Review
  #
  crtitle(x)
  
  if (length(x$pval) != 0) {
    format.pvalue <- formatPT(x$pval, digits = digits.pval,
                              scientific = scientific.pval,
                              zero = zero.pval, JAMA = JAMA.pval)
    
    if (!grepl("<", format.pvalue))
      format.pvalue <- paste("=", format.pvalue)
    
    cat(paste0(x$method, "\n\n"))
    #
    cat(paste0("Test result: ",
               names(x$statistic), " = ",
               formatN(x$statistic, digits.stat, "NA",
                       big.mark = big.mark),
               ", ",
               if (!is.null(x$df))
                 paste0("df = ", x$df, ", "),
               "p-value ",
               format.pvalue,
               "\n",
               #
               "Bias estimate: ",
               formatN(x$estimate[1], digits, "NA", big.mark = big.mark),
               " (SE = ",
               formatN(x$estimate[2], digits.se, "NA", big.mark = big.mark),
               ")\n"
    )
    )
    #
    if (details.methods) {
      if (x$var.model != "") {
        cat("\nDetails:\n")
        cat(paste0("- ", x$var.model,
                   " residual heterogeneity variance (",
                   formatPT(x$tau^2,
                            lab = TRUE, labval = text.tau2,
                            digits = digits.tau2,
                            lab.NA = "NA",
                            big.mark = big.mark),
                   ")\n"))
        #
        if (x$var.model == "additive") {
          i.lab.method.tau <-
            charmatch(x$x$method.tau, c(gs("meth4tau"), ""), nomatch = NA)
          #
          lab.method.tau <-
            c("- DerSimonian-Laird estimator",
              "- Paule-Mandel estimator",
              "- restricted maximum-likelihood estimator",
              "- maximum-likelihood estimator",
              "- Hunter-Schmidt estimator",
              "- Sidik-Jonkman estimator",
              "- Hedges estimator",
              "- empirical Bayes estimator",
              ""
            )[i.lab.method.tau]
          #
          if (lab.method.tau != "")
            lab.method.tau <- paste(lab.method.tau, "for", text.tau2)
          #
          cat(paste0(lab.method.tau, "\n"))
        }
        #
        if (x$method.bias == "Egger")
          detail.predictor <- "standard error"
        else if (x$method.bias == "Thompson")
          detail.predictor <- "standard error"
        else if (x$method.bias == "Harbord")
          detail.predictor <- "standard error of score"
        else if (x$method.bias == "Macaskill")
          detail.predictor <- "total sample size"
        else if (x$method.bias == "Peters")
          detail.predictor <- "inverse of total sample size"
        else if (x$method.bias == "Deeks")
          detail.predictor <- "inverse of the squared effective sample size"
        else if (x$method.bias == "Pustejovsky")
          detail.predictor <-
          "square root of the sum of the inverse group sample sizes"
        #
        cat(paste0("- predictor: ", detail.predictor, "\n"))
        #
        if (x$method.bias == "Egger")
          detail.weights <- "inverse variance"
        else if (x$method.bias == "Thompson")
          detail.weights <- "inverse variance"
        else if (x$method.bias == "Harbord")
          detail.weights <- "inverse variance of score"
        else if (x$method.bias == "Macaskill")
          detail.weights <- "inverse variance of average event probability"
        else if (x$method.bias == "Peters")
          detail.weights <- "inverse variance of average event probability"
        else if (x$method.bias == "Deeks")
          detail.weights <- "effective sample size"
        else if (x$method.bias == "Pustejovsky")
          detail.weights <- "inverse variance"
        #
        cat(paste0("- weight:    ", detail.weights, "\n"))
      }
    }
    #
    if (x$method.bias == "Begg")
      detail.ref <- "Begg & Mazumdar (1993), Biometrics"
    else if (x$method.bias == "Egger")
      detail.ref <- "Egger et al. (1997), BMJ"
    else if (x$method.bias == "Thompson")
      detail.ref <- "Thompson & Sharp (1999), Stat Med"
    else if (x$method.bias == "Schwarzer")
      detail.ref <- "Schwarzer et al. (2006), Stat Med"
    else if (x$method.bias == "Harbord")
      detail.ref <- "Harbord et al. (2006), Stat Med"
    else if (x$method.bias == "Macaskill")
      detail.ref <- "Macaskill et al. (2001), Stat Med, method FPV"
    else if (x$method.bias == "Peters")
      detail.ref <- "Peters et al. (2006), JAMA"
    else if (x$method.bias == "Deeks")
      detail.ref <- "Deeks et al. (2005), J Clin Epid"
    else if (x$method.bias == "Pustejovsky")
      detail.ref <- "Pustejovsky & Rodgers (2019), RSM"
    else
      detail.ref <- ""
    #
    if (detail.ref != "")
      cat(
        paste0(
          if (details.methods && x$var.model != "")
            "- reference: "
          else
            "\nReference: ",
          detail.ref, "\n"))
  }
  else {
    #
    # Check whether number of studies is too small:
    #
    if (length(x$k) != 0 & length(x$k.min != 0)) {
      if (x$k <= x$k.min) {
        if (x$k <= 2)
          warning("Number of studies (k=",  x$k,
                  ") too small to test for small study effects.",
                  call. = FALSE)
        else
          warning("Number of studies (k=",  x$k,
                  ") too small to test for small study effects (k.min=",
                  x$k.min, "). Change argument 'k.min' if appropriate.",
                  call. = FALSE)
      }
    }
    #
    # Check whether meta-analysis has subgroups:
    #
    if (length(x$subgroup) != 0)
      warning("No test for small study effects conducted ",
              "for meta-analysis with subgroups.",
              call. = FALSE)
  }
  
  invisible(NULL)
}





#' @rdname metabias
#' @export metabias


metabias <- function(x, ...) 
  UseMethod("metabias")





#' @rdname metabias
#' @method metabias default
#' @export


metabias.default <- function(x, seTE,
                             method.bias = "Egger",
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
  method.bias <- setchar(method.bias, c("Begg", "Egger", "Thompson"))
  
  
  ##
  ##
  ## (2) Do meta-analysis
  ##
  ##
  m <- metagen(x, seTE, method.tau.ci = "")
  
  
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
