#' Meta-analysis of continuous outcome data
#' 
#' @description
#' Calculation of fixed and random effects estimates for meta-analyses
#' with continuous outcome data; inverse variance weighting is used
#' for pooling.
#' 
#' @param n.e Number of observations in experimental group.
#' @param mean.e Estimated mean in experimental group.
#' @param sd.e Standard deviation in experimental group.
#' @param n.c Number of observations in control group.
#' @param mean.c Estimated mean in control group.
#' @param sd.c Standard deviation in control group.
#' @param studlab An optional vector with study labels.
#' @param data An optional data frame containing the study
#'   information.
#' @param subset An optional vector specifying a subset of studies to
#'   be used.
#' @param exclude An optional vector specifying studies to exclude
#'   from meta-analysis, however, to include in printouts and forest
#'   plots.
#' @param level The level used to calculate confidence intervals for
#'   individual studies.
#' @param level.comb The level used to calculate confidence intervals
#'   for pooled estimates.
#' @param comb.fixed A logical indicating whether a fixed effect
#'   meta-analysis should be conducted.
#' @param comb.random A logical indicating whether a random effects
#'   meta-analysis should be conducted.
#' @param prediction A logical indicating whether a prediction
#'   interval should be printed.
#' @param level.predict The level used to calculate prediction
#'   interval for a new study.
#' @param hakn A logical indicating whether the method by Hartung and
#'   Knapp should be used to adjust test statistics and confidence
#'   intervals.
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2}. Either
#'   \code{"DL"}, \code{"PM"}, \code{"REML"}, \code{"ML"},
#'   \code{"HS"}, \code{"SJ"}, \code{"HE"}, or \code{"EB"}, can be
#'   abbreviated.
#' @param tau.preset Prespecified value for the square-root of the
#'   between-study variance \eqn{\tau^2}.
#' @param TE.tau Overall treatment effect used to estimate the
#'   between-study variance tau-squared.
#' @param tau.common A logical indicating whether tau-squared should
#'   be the same across subgroups.
#' @param method.bias A character string indicating which test is to
#'   be used.  Either \code{"rank"}, \code{"linreg"}, or \code{"mm"},
#'   can be abbreviated.  See function \code{\link{metabias}}
#' @param backtransf A logical indicating whether results for ratio of
#'   means (\code{sm="ROM"}) should be back transformed in printouts
#'   and plots. If TRUE (default), results will be presented as ratio
#'   of means; otherwise log ratio of means will be shown.
#' @param title Title of meta-analysis / systematic review.
#' @param complab Comparison label.
#' @param outclab Outcome label.
#' @param label.e Label for experimental group.
#' @param label.c Label for control group.
#' @param label.left Graph label on left side of forest plot.
#' @param label.right Graph label on right side of forest plot.
#' @param sm A character string indicating which summary measure
#'   (\code{"MD"}, \code{"SMD"}, or \code{"ROM"}) is to be used for
#'   pooling of studies.
#' @param pooledvar A logical indicating if a pooled variance should
#'   be used for the mean difference (\code{sm="MD"}).
#' @param method.smd A character string indicating which method is
#'   used to estimate the standardised mean difference
#'   (\code{sm="SMD"}). Either \code{"Hedges"} for Hedges' g
#'   (default), \code{"Cohen"} for Cohen's d, or \code{"Glass"} for
#'   Glass' delta, can be abbreviated.
#' @param sd.glass A character string indicating which standard
#'   deviation is used in the denominator for Glass' method to
#'   estimate the standardised mean difference. Either
#'   \code{"control"} using the standard deviation in the control
#'   group (\code{sd.c}) or \code{"experimental"} using the standard
#'   deviation in the experimental group (\code{sd.e}), can be
#'   abbreviated.
#' @param exact.smd A logical indicating whether exact formulae should
#'   be used in estimation of the standardised mean difference and its
#'   standard error (see Details).
#' @param byvar An optional vector containing grouping information
#'   (must be of same length as \code{n.e}).
#' @param bylab A character string with a label for the grouping
#'   variable.
#' @param print.byvar A logical indicating whether the name of the
#'   grouping variable should be printed in front of the group labels.
#' @param byseparator A character string defining the separator
#'   between label and levels of grouping variable.
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in meta object.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if studies are excluded from meta-analysis due to zero
#'   standard deviations).
#' @param control An optional list to control the iterative process to
#'   estimate the between-study variance tau^2. This argument is
#'   passed on to \code{\link[metafor]{rma.uni}}.
#' 
#' @details
#' Calculation of fixed and random effects estimates for meta-analyses
#' with continuous outcome data; inverse variance weighting is used
#' for pooling.
#' 
#' Three different types of summary measures are available for continuous
#' outcomes:
#' \itemize{
#' \item mean difference (argument \code{sm = "MD"})
#' \item standardised mean difference (\code{sm = "SMD"})
#' \item ratio of means (\code{sm = "ROM"})
#' }
#' 
#' Meta-analysis of ratio of means -- also called response ratios --
#' is described in Hedges et al. (1999) and Friedrich et al. (2008).
#' 
#' For the standardised mean difference three methods are implemented:
#' \itemize{
#' \item Hedges' g (default, \code{method.smd = "Hedges"}) - see
#'   Hedges (1981)
#' \item Cohen's d (\code{method.smd = "Cohen"}) - see Cohen (1988)
#' \item Glass' delta (\code{method.smd = "Glass"}) - see Glass (1976)
#' }
#'
#' Hedges (1981) calculated the exact bias in Cohen's d which is a
#' ratio of gamma distributions with the degrees of freedom,
#' i.e. total sample size minus two, as argument. By default (argument
#' \code{exact.smd = FALSE}), an accurate approximation of this bias
#' provided in Hedges (1981) is utilised for Hedges' g as well as its
#' standard error; these approximations are also used in RevMan
#' 5. Following Borenstein et al. (2009) these approximations are not
#' used in the estimation of Cohen's d. White and Thomas (2005) argued
#' that approximations are unnecessary with modern software and
#' accordingly promote to use the exact formulae; this is possible
#' using argument \code{exact.smd = TRUE}. For Hedges' g the exact
#' formulae are used to calculate the standardised mean difference as
#' well as the standard error; for Cohen's d the exact formula is only
#' used to calculate the standard error. In typical applications (with
#' sample sizes above 10), the differences between using the exact
#' formulae and the approximation will be minimal.
#' 
#' For Glass' delta, by default (argument \code{sd.glass =
#' "control"}), the standard deviation in the control group
#' (\code{sd.c}) is used in the denominator of the standard mean
#' difference. The standard deviation in the experimental group
#' (\code{sd.e}) can be used by specifying \code{sd.glass =
#' "experimental"}.
#' 
#' Calculations are conducted on the log scale for ratio of means
#' (\code{sm = "ROM"}). Accordingly, list elements \code{TE},
#' \code{TE.fixed}, and \code{TE.random} contain the logarithm of
#' ratio of means. In printouts and plots these values are back
#' transformed if argument \code{backtransf = TRUE}.
#' 
#' For several arguments defaults settings are utilised (assignments
#' using \code{\link{gs}} function). These defaults can be changed
#' using the \code{\link{settings.meta}} function.
#' 
#' Internally, both fixed effect and random effects models are
#' calculated regardless of values choosen for arguments
#' \code{comb.fixed} and \code{comb.random}. Accordingly, the estimate
#' for the random effects model can be extracted from component
#' \code{TE.random} of an object of class \code{"meta"} even if
#' argument \code{comb.random = FALSE}. However, all functions in R
#' package \bold{meta} will adequately consider the values for
#' \code{comb.fixed} and \code{comb.random}. E.g. function
#' \code{\link{print.meta}} will not print results for the random
#' effects model if \code{comb.random = FALSE}.
#' 
#' The function \code{metagen} is called internally to calculate
#' individual and overall treatment estimates and standard errors.
#' 
#' A prediction interval for the treatment effect of a new study is
#' calculated (Higgins et al., 2009) if arguments \code{prediction}
#' and \code{comb.random} are \code{TRUE}.
#' 
#' R function \code{\link{update.meta}} can be used to redo the
#' meta-analysis of an existing metacont object by only specifying
#' arguments which should be changed.
#' 
#' For the random effects, the method by Hartung and Knapp (2003) is
#' used to adjust test statistics and confidence intervals if argument
#' \code{hakn = TRUE}.
#' 
#' The DerSimonian-Laird estimate (1986) is used in the random effects
#' model if \code{method.tau = "DL"}. The iterative Paule-Mandel
#' method (1982) to estimate the between-study variance is used if
#' argument \code{method.tau = "PM"}.  Internally, R function
#' \code{paulemandel} is called which is based on R function
#' \code{mpaule.default} from R package \bold{metRology} from S.L.R.
#' Ellison <s.ellison at lgc.co.uk>.
#' 
#' If R package \bold{metafor} (Viechtbauer 2010) is installed, the
#' following methods to estimate the between-study variance
#' \eqn{\tau^2} (argument \code{method.tau}) are also available:
#' \itemize{
#' \item Restricted maximum-likelihood estimator (\code{method.tau =
#'   "REML"})
#' \item Maximum-likelihood estimator (\code{method.tau = "ML"})
#' \item Hunter-Schmidt estimator (\code{method.tau = "HS"})
#' \item Sidik-Jonkman estimator (\code{method.tau = "SJ"})
#' \item Hedges estimator (\code{method.tau = "HE"})
#' \item Empirical Bayes estimator (\code{method.tau = "EB"})
#' }
#' For these methods the R function \code{rma.uni} of R package
#' \bold{metafor} is called internally. See help page of R function
#' \code{rma.uni} for more details on these methods to estimate
#' between-study variance.
#' 
#' @return
#' An object of class \code{c("metacont", "meta")} with corresponding
#' \code{print}, \code{summary}, and \code{forest} functions. The
#' object is a list containing the following components:
#'
#' \item{n.e, mean.e, sd.e,}{As defined above.}
#' \item{n.c, mean.c, sd.c,}{As defined above.}
#' \item{studlab, exclude, sm, level, level.comb,}{As defined above.}
#' \item{comb.fixed, comb.random,}{As defined above.}
#' \item{pooledvar, method.smd, sd.glass,}{As defined above.}
#' \item{hakn, method.tau, tau.preset, TE.tau, method.bias,}{As
#'   defined above.}
#' \item{tau.common, title, complab, outclab,}{As defined above.}
#' \item{label.e, label.c, label.left, label.right,}{As defined
#'   above.}
#'
#' \item{byvar, bylab, print.byvar, byseparator}{As defined above.}
#' \item{TE, seTE}{Estimated treatment effect and standard error of
#'   individual studies.}
#' \item{lower, upper}{Lower and upper confidence interval limits for
#'   individual studies.}
#' \item{zval, pval}{z-value and p-value for test of treatment effect
#'   for individual studies.}
#' \item{w.fixed, w.random}{Weight of individual studies (in fixed and
#'   random effects model).}
#' \item{TE.fixed, seTE.fixed}{Estimated overall treatment effect and
#'   standard error (fixed effect model).}
#' \item{lower.fixed, upper.fixed}{Lower and upper confidence interval
#'   limits (fixed effect model).}
#' \item{zval.fixed, pval.fixed}{z-value and p-value for test of
#'   overall treatment effect (fixed effect model).}
#' \item{TE.random, seTE.random}{Estimated overall treatment effect
#'   and standard error (random effects model).}
#' \item{lower.random, upper.random}{Lower and upper confidence
#'   interval limits (random effects model).}
#' \item{zval.random, pval.random}{z-value or t-value and
#'   corresponding p-value for test of overall treatment effect
#'   (random effects model).}
#' \item{prediction, level.predict}{As defined above.}
#' \item{seTE.predict}{Standard error utilised for prediction
#'   interval.}
#' \item{lower.predict, upper.predict}{Lower and upper limits of
#'   prediction interval.}
#' \item{k}{Number of studies combined in meta-analysis.}
#' \item{Q}{Heterogeneity statistic Q.}
#' \item{df.Q}{Degrees of freedom for heterogeneity statistic.}
#' \item{pval.Q}{P-value of heterogeneity test.}
#' \item{tau}{Square-root of between-study variance.}
#' \item{se.tau}{Standard error of square-root of between-study
#'   variance.}
#' \item{C}{Scaling factor utilised internally to calculate common
#'   tau-squared across subgroups.}
#' \item{df.hakn}{Degrees of freedom for test of treatment effect for
#'   Hartung-Knapp method (only if \code{hakn = TRUE}).}
#' \item{method}{Pooling method: \code{"Inverse"}.}
#' \item{bylevs}{Levels of grouping variable - if \code{byvar} is not
#'   missing.}
#' \item{TE.fixed.w, seTE.fixed.w}{Estimated treatment effect and
#'   standard error in subgroups (fixed effect model) - if
#'   \code{byvar} is not missing.}
#' \item{lower.fixed.w, upper.fixed.w}{Lower and upper confidence
#'   interval limits in subgroups (fixed effect model) - if
#'   \code{byvar} is not missing.}
#' \item{zval.fixed.w, pval.fixed.w}{z-value and p-value for test of
#'   treatment effect in subgroups (fixed effect model) - if
#'   \code{byvar} is not missing.}
#' \item{TE.random.w, seTE.random.w}{Estimated treatment effect and
#'   standard error in subgroups (random effects model) - if
#'   \code{byvar} is not missing.}
#' \item{lower.random.w, upper.random.w}{Lower and upper confidence
#'   interval limits in subgroups (random effects model) - if
#'   \code{byvar} is not missing.}
#' \item{zval.random.w, pval.random.w}{z-value or t-value and
#'   corresponding p-value for test of treatment effect in subgroups
#'   (random effects model) - if \code{byvar} is not missing.}
#' \item{w.fixed.w, w.random.w}{Weight of subgroups (in fixed and
#'   random effects model) - if \code{byvar} is not missing.}
#' \item{df.hakn.w}{Degrees of freedom for test of treatment effect
#'   for Hartung-Knapp method in subgroups - if \code{byvar} is not
#'   missing and \code{hakn = TRUE}.}
#' \item{n.e.w}{Number of observations in experimental group in
#'   subgroups - if \code{byvar} is not missing.}
#' \item{n.c.w}{Number of observations in control group in subgroups -
#'   if \code{byvar} is not missing.}
#' \item{k.w}{Number of studies combined within subgroups - if
#'   \code{byvar} is not missing.}
#' \item{k.all.w}{Number of all studies in subgroups - if \code{byvar}
#'   is not missing.}
#' \item{Q.w.fixed}{Overall within subgroups heterogeneity statistic Q
#'   (based on fixed effect model) - if \code{byvar} is not missing.}
#' \item{Q.w.random}{Overall within subgroups heterogeneity statistic
#'   Q (based on random effects model) - if \code{byvar} is not
#'   missing (only calculated if argument \code{tau.common} is TRUE).}
#' \item{df.Q.w}{Degrees of freedom for test of overall within
#'   subgroups heterogeneity - if \code{byvar} is not missing.}
#' \item{pval.Q.w.fixed}{P-value of within subgroups heterogeneity
#'   statistic Q (based on fixed effect model) - if \code{byvar} is
#'   not missing.}
#' \item{pval.Q.w.random}{P-value of within subgroups heterogeneity
#'   statistic Q (based on random effects model) - if \code{byvar} is
#'   not missing.}
#' \item{Q.b.fixed}{Overall between subgroups heterogeneity statistic
#'   Q (based on fixed effect model) - if \code{byvar} is not
#'   missing.}
#' \item{Q.b.random}{Overall between subgroups heterogeneity statistic
#'   Q (based on random effects model) - if \code{byvar} is not
#'   missing.}
#' \item{df.Q.b}{Degrees of freedom for test of overall between
#'   subgroups heterogeneity - if \code{byvar} is not missing.}
#' \item{pval.Q.b.fixed}{P-value of between subgroups heterogeneity
#'   statistic Q (based on fixed effect model) - if \code{byvar} is
#'   not missing.}
#' \item{pval.Q.b.random}{P-value of between subgroups heterogeneity
#'   statistic Q (based on random effects model) - if \code{byvar} is
#'   not missing.}
#' \item{tau.w}{Square-root of between-study variance within subgroups
#'   - if \code{byvar} is not missing.}
#' \item{C.w}{Scaling factor utilised internally to calculate common
#'   tau-squared across subgroups - if \code{byvar} is not missing.}
#'   \item{H.w}{Heterogeneity statistic H within subgroups - if
#'   \code{byvar} is not missing.}
#' \item{lower.H.w, upper.H.w}{Lower and upper confidence limti for
#'   heterogeneity statistic H within subgroups - if \code{byvar} is
#'   not missing.}
#' \item{I2.w}{Heterogeneity statistic I2 within subgroups - if
#'   \code{byvar} is not missing.}
#' \item{lower.I2.w, upper.I2.w}{Lower and upper confidence limit for
#'   heterogeneity statistic I2 within subgroups - if \code{byvar} is
#'   not missing.}
#' \item{keepdata}{As defined above.}
#' \item{data}{Original data (set) used in function call (if
#'   \code{keepdata = TRUE}).}
#' \item{subset}{Information on subset of original data used in
#'   meta-analysis (if \code{keepdata = TRUE}).}
#' \item{call}{Function call.}
#' \item{version}{Version of R package \bold{meta} used to create
#'   object.}
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{update.meta}}, \code{\link{metabin}},
#'   \code{\link{metagen}}
#' 
#' @references
#' Borenstein M, Hedges LV, Higgins JPT, Rothstein HR (2009):
#' \emph{Introduction to Meta-Analysis}.
#' Chichester: Wiley
#' 
#' Cohen J (1988):
#' \emph{Statistical Power Analysis for the Behavioral Sciences
#'   (second ed.)}.
#' Lawrence Erlbaum Associates
#' 
#' Cooper H & Hedges LV (1994):
#' \emph{The Handbook of Research Synthesis}.
#' Newbury Park, CA: Russell Sage Foundation
#' 
#' DerSimonian R & Laird N (1986):
#' Meta-analysis in clinical trials.
#' \emph{Controlled Clinical Trials},
#' \bold{7}, 177--88
#' 
#' Friedrich JO, Adhikari NK, Beyene J (2008):
#' The ratio of means method as an alternative to mean differences for
#' analyzing continuous outcome variables in meta-analysis: A
#' simulation study.
#' \emph{BMC Medical Research Methodology},
#' \bold{8}, 32
#' 
#' Glass G (1976):
#' Primary, secondary, and meta-analysis of research.
#' \emph{Educational Researcher},
#' \bold{5}, 3--8
#' 
#' Hartung J & Knapp G (2001):
#' On tests of the overall treatment effect in meta-analysis with
#' normally distributed responses.
#' \emph{Statistics in Medicine},
#' \bold{20}, 1771--82
#' 
#' Hedges LV (1981):
#' Distribution theory for Glass's estimator of effect size and
#' related estimators.
#' \emph{Journal of Educational and Behavioral Statistics},
#' \bold{6}, 107--28
#' 
#' Hedges LV, Gurevitch J, Curtis PS (1999):
#' The meta-analysis of response ratios in experimental ecology.
#' \emph{Ecology},
#' \bold{80}, 1150--6
#' 
#' Higgins JPT, Thompson SG, Spiegelhalter DJ (2009):
#' A re-evaluation of random-effects meta-analysis.
#' \emph{Journal of the Royal Statistical Society: Series A},
#' \bold{172}, 137--59
#' 
#' Knapp G & Hartung J (2003):
#' Improved tests for a random effects meta-regression with a single
#' covariate.
#' \emph{Statistics in Medicine},
#' \bold{22}, 2693--710
#' 
#' Paule RC & Mandel J (1982):
#' Consensus values and weighting factors.
#' \emph{Journal of Research of the National Bureau of Standards},
#' \bold{87}, 377--85
#' 
#' \emph{Review Manager (RevMan)}
#' [Computer program]. Version 5.3.
#' Copenhagen: The Nordic Cochrane Centre, The Cochrane Collaboration, 2014
#' 
#' Viechtbauer W (2010):
#' Conducting Meta-Analyses in R with the Metafor Package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
#' 
#' White IR, Thomas J (2005):
#' Standardized mean differences in individually-randomized and
#' cluster-randomized trials, with applications to meta-analysis.
#' \emph{Clinical Trials},
#' \bold{2}, 141--51
#' 
#' @examples
#' data(Fleiss93cont)
#'
#' # Meta-analysis with Hedges' g as effect measure
#' #
#' m1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c,
#'                data = Fleiss93cont, sm = "SMD")
#' m1
#' forest(m1)
#' 
#' # Use Cohen's d instead of Hedges' g as effect measure
#' #
#' update(m1, method.smd = "Cohen")
#' 
#' # Use Glass' delta instead of Hedges' g as effect measure
#' #
#' update(m1, method.smd = "Glass")
#' 
#' # Use Glass' delta based on the standard deviation in the experimental group
#' #
#' update(m1, method.smd = "Glass", sd.glass = "experimental")
#' 
#' # Calculate Hedges' g based on exact formulae
#' #
#' update(m1, exact.smd = TRUE)
#' 
#' data(amlodipine)
#' m2 <- metacont(n.amlo, mean.amlo, sqrt(var.amlo),
#'                n.plac, mean.plac, sqrt(var.plac),
#'                data = amlodipine, studlab = study)
#' summary(m2)
#' 
#' # Use pooled variance
#' #
#' summary(update(m2, pooledvar = TRUE))
#' 
#' # Meta-analysis of response ratios (Hedges et al., 1999)
#' #
#' data(woodyplants)
#' m3 <- metacont(n.elev, mean.elev, sd.elev,
#' 		  n.amb, mean.amb, sd.amb,
#'                data = woodyplants, sm = "ROM")
#' summary(m3)
#' summary(m3, backtransf = FALSE)
#' 
#' @export metacont


metacont <- function(n.e, mean.e, sd.e, n.c, mean.c, sd.c, studlab,
                     ##
                     data = NULL, subset = NULL, exclude = NULL,
                     ##
                     sm = gs("smcont"),
                     ##
                     pooledvar = gs("pooledvar"),
                     method.smd = gs("method.smd"),
                     sd.glass = gs("sd.glass"),
                     exact.smd = gs("exact.smd"),
                     ##
                     level = gs("level"), level.comb = gs("level.comb"),
                     comb.fixed = gs("comb.fixed"),
                     comb.random = gs("comb.random"),
                     ##
                     hakn = gs("hakn"),
                     method.tau = gs("method.tau"),
                     tau.preset = NULL, TE.tau = NULL,
                     tau.common = gs("tau.common"),
                     ##
                     prediction = gs("prediction"),
                     level.predict = gs("level.predict"),
                     ##
                     method.bias = gs("method.bias"),
                     ##
                     backtransf = gs("backtransf"),
                     title = gs("title"), complab = gs("complab"),
                     outclab = "",
                     label.e = gs("label.e"), label.c = gs("label.c"),
                     label.left = gs("label.left"),
                     label.right = gs("label.right"),
                     ##
                     byvar, bylab, print.byvar = gs("print.byvar"),
                     byseparator = gs("byseparator"),
                     ##
                     keepdata = gs("keepdata"),
                     warn = gs("warn"),
                     ##
                     control = NULL
                     ) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chknull(sm)
  chklevel(level)
  chklevel(level.comb)
  chklogical(comb.fixed)
  chklogical(comb.random)
  ##
  chklogical(hakn)
  method.tau <- setchar(method.tau,
                        c("DL", "PM", "REML", "ML", "HS", "SJ", "HE", "EB"))
  chklogical(tau.common)
  ##
  chklogical(prediction)
  chklevel(level.predict)
  ##
  method.bias <- setchar(method.bias,
                         c("rank", "linreg", "mm", "count", "score", "peters"))
  ##
  chklogical(keepdata)
  ##
  ## Additional arguments / checks for metacont objects
  ##
  fun <- "metacont"
  sm <- setchar(sm, .settings$sm4cont)
  chklogical(pooledvar)
  method.smd <- setchar(method.smd, c("Hedges", "Cohen", "Glass"))
  sd.glass <- setchar(sd.glass, c("control", "experimental"))
  chklogical(warn)
  
  
  ##
  ##
  ## (2) Read data
  ##
  ##
  nulldata <- is.null(data)
  ##
  if (nulldata)
    data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Catch 'n.e', 'mean.e', 'sd.e', 'n.c', 'mean.c', and 'sd.c' from data:
  ##
  n.e <- eval(mf[[match("n.e", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  chknull(n.e)
  k.All <- length(n.e)
  ##
  mean.e <- eval(mf[[match("mean.e", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknull(mean.e)
  ##
  sd.e <- eval(mf[[match("sd.e", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  chknull(sd.e)
  ##
  n.c <- eval(mf[[match("n.c", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  chknull(n.c)
  ##
  mean.c <- eval(mf[[match("mean.c", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknull(mean.c)
  ##
  sd.c <- eval(mf[[match("sd.c", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  chknull(sd.c)
  ##
  ## Catch 'studlab', 'byvar', 'subset' and 'exclude' from data:
  ##
  studlab <- eval(mf[[match("studlab", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  studlab <- setstudlab(studlab, k.All)
  ##
  byvar <- eval(mf[[match("byvar", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  by <- !is.null(byvar)
  ##
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  missing.subset <- is.null(subset)
  ##
  exclude <- eval(mf[[match("exclude", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  missing.exclude <- is.null(exclude)
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  chklength(mean.e, k.All, fun)
  chklength(sd.e, k.All, fun)
  chklength(n.c, k.All, fun)
  chklength(mean.c, k.All, fun)
  chklength(sd.c, k.All, fun)
  chklength(studlab, k.All, fun)
  ##
  if (by)
    chklength(byvar, k.All, fun)
  ##
  ## Additional checks
  ##
  if (!by & tau.common) {
    warning("Value for argument 'tau.common' set to FALSE as argument 'byvar' is missing.")
    tau.common <- FALSE
  }
  if (by & !tau.common & !is.null(tau.preset)) {
    warning("Argument 'tau.common' set to TRUE as argument tau.preset is not NULL.")
    tau.common <- TRUE
  }
  
  
  ##
  ##
  ## (4) Subset, exclude studies, and subgroups
  ##
  ##
  if (!missing.subset)
    if ((is.logical(subset) & (sum(subset) > k.All)) ||
        (length(subset) > k.All))
      stop("Length of subset is larger than number of studies.")
  ##
  if (!missing.exclude) {
    if ((is.logical(exclude) & (sum(exclude) > k.All)) ||
        (length(exclude) > k.All))
      stop("Length of argument 'exclude' is larger than number of studies.")
    ##
    exclude2 <- rep(FALSE, k.All)
    exclude2[exclude] <- TRUE
    exclude <- exclude2
  }
  else
    exclude <- rep(FALSE, k.All)
  ##
  if (by) {
    chkmiss(byvar)
    byvar.name <- byvarname(mf[[match("byvar", names(mf))]])
    bylab <- if (!missing(bylab) && !is.null(bylab)) bylab else byvar.name
  }
  
  
  ##
  ##
  ## (5) Store complete dataset in list object data
  ##     (if argument keepdata is TRUE)
  ##
  ##
  if (keepdata) {
    if (nulldata)
      data <- data.frame(.n.e = n.e)
    else
      data$.n.e <- n.e
    ##
    data$.mean.e <- mean.e
    data$.sd.e <- sd.e
    data$.n.c <- n.c
    data$.mean.c <- mean.c
    data$.sd.c <- sd.c
    data$.studlab <- studlab
    ##
    if (by)
      data$.byvar <- byvar
    ##
    if (!missing.subset) {
      if (length(subset) == dim(data)[1])
        data$.subset <- subset
      else {
        data$.subset <- FALSE
        data$.subset[subset] <- TRUE
      }
    }
    ##
    if (!missing.exclude)
      data$.exclude <- exclude
  }
  
  
  ##
  ##
  ## (6) Use subset for analysis
  ##
  ##
  if (!missing.subset) {
    n.e <- n.e[subset]
    mean.e <- mean.e[subset]
    sd.e <- sd.e[subset]
    n.c <- n.c[subset]
    mean.c <- mean.c[subset]
    sd.c <- sd.c[subset]
    studlab <- studlab[subset]
    ##
    exclude <- exclude[subset]
    ##
    if (by)
      byvar <- byvar[subset]
  }
  ##
  ## Determine total number of studies
  ##
  k.all <- length(n.e)
  ##
  if (k.all == 0)
    stop("No studies to combine in meta-analysis.")
  ##
  ## No meta-analysis for a single study
  ##
  if (k.all == 1) {
    comb.fixed  <- FALSE
    comb.random <- FALSE
    prediction  <- FALSE
  }
  ##
  ## Check variable values
  ##
  chknumeric(n.e)
  chknumeric(mean.e)
  chknumeric(sd.e)
  chknumeric(n.c)
  chknumeric(mean.c)
  chknumeric(sd.c)
  ##
  ## Recode integer as numeric:
  ##
  n.e    <- int2num(n.e)
  mean.e <- int2num(mean.e)
  sd.e   <- int2num(sd.e)
  n.c    <- int2num(n.c)
  mean.c <- int2num(mean.c)
  sd.c   <- int2num(sd.c)
  
  
  ##
  ##
  ## (7) Calculate results for individual studies
  ##
  ##
  npn.n <- npn(n.e) | npn(n.c)
  ##
  N <- n.e + n.c
  if (sm == "MD" | sm == "ROM")
    var.pooled <- ((n.e - 1) * sd.e^2 + (n.c - 1) * sd.c^2) / (N - 2)
  ##
  if (any(npn.n) & warn)
    warning("Studies with non-positive values for n.e and / or n.c get no weight in meta-analysis.")
  ##
  if (sm == "MD") {
    TE <- ifelse(npn.n, NA, mean.e - mean.c)
    ##
    if (pooledvar)
      seTE <- ifelse(npn.n, NA,
                     sqrt(var.pooled * (1 / n.e + 1 / n.c)))
    else
      seTE <- ifelse(npn.n, NA,
                     sqrt(sd.e^2 / n.e + sd.c^2 / n.c))
    ##
    seTE[is.na(TE)] <- NA
  }
  else if (sm == "SMD") {
    J <- function(x) gamma(x / 2) / (sqrt(x / 2) * gamma((x - 1) / 2))
    K <- function(x) 1 - (x - 2) / (x * J(x)^2)
    ##
    if (method.smd %in% c("Hedges", "Cohen"))
      S.within <- sqrt(((n.e - 1) * sd.e^2 + (n.c - 1) * sd.c^2) / (N - 2))
    else
      S.within <- if (sd.glass == "control") sd.c else sd.e
    ##
    smd <- ifelse(npn.n, NA, (mean.e - mean.c) / S.within)
    ##
    if (method.smd == "Cohen") {
      ##
      ## Borenstein et al. (2009), p. 26-27;
      ## White and Thomas (2005), p. 143
      ##
      TE <- smd
      if (exact.smd) {
        J <- function(x) gamma(x / 2) / (sqrt(x / 2) * gamma((x - 1) / 2))
        K <- function(x) 1 - (x - 2) / (x * J(x)^2)
        seTE <- ifelse(npn.n, NA,
                       sqrt(N / (n.e * n.c) + (J(N - 2) * smd)^2 * K(N - 2)))
      }
      else
        seTE <- ifelse(npn.n, NA,
                       sqrt(N / (n.e * n.c) + TE^2 / (2 * N)))
    }
    else if (method.smd == "Hedges") {
      ##
      ## Hedges and Olkin (1985); White and Thomas (2005), p. 143;
      ## formulae used in RevMan 5 (exact.smd = FALSE)
      ##
      if (exact.smd) {
        J <- function(x) gamma(x / 2) / (sqrt(x / 2) * gamma((x - 1) / 2))
        K <- function(x) 1 - (x - 2) / (x * J(x)^2)
      }
      else {
        J <- function(x) 1 - 3 / (4 * x - 1)
        K <- function(x) 1 / (2 * (x - 1.94))
      }
      ##
      TE   <- J(N - 2) * smd
      seTE <- ifelse(npn.n, NA,
                     sqrt(N / (n.e * n.c) + TE^2 * K(N - 2)))
    }
    else if (method.smd == "Glass") {
      ##
      ## see Cooper & Hedges (1994), p. 238
      ##
      n.g  <- if (sd.glass == "control") n.c else n.e
      ##
      TE <- smd
      seTE <- ifelse(npn.n, NA,
                     sqrt(N / (n.e * n.c) + TE^2 / (2 * n.g - 1)))
    }
    ##
    seTE[is.na(TE)] <- NA
  }
  ##
  else if (sm == "ROM") {
    npn.mean <- npn(mean.e) | npn(mean.c)
    ##
    if (any(npn.mean) & warn)
      warning("Studies with negative or zero means get no weight in meta-analysis.")

    TE <- ifelse(npn.n | npn.mean, NA, log(mean.e / mean.c))
    ##
    if (pooledvar)
      seTE <- ifelse(npn.n, NA,
                     sqrt(var.pooled * (1 / (n.e * mean.e^2) + 1 / (n.c * mean.c^2))))
    else
      seTE <- ifelse(npn.n | npn.mean, NA,
                     sqrt(sd.e^2 / (n.e * mean.e^2) + sd.c^2 / (n.c * mean.c^2)))
    ##
    seTE[is.na(TE)] <- NA
  }
  ##
  ## Studies with non-positive variance get zero weight in meta-analysis
  ##
  sel <- sd.e <= 0 | sd.c <= 0
  ##
  if (any(sel, na.rm = TRUE) & warn)
    warning("Studies with non-positive values for sd.e or sd.c get no weight in meta-analysis.")
  ##
  seTE[sel] <- NA
  ##
  if (sm == "SMD")
    TE[sel] <- NA
  
  
  ##
  ##
  ## (8) Do meta-analysis
  ##
  ##
  m <- metagen(TE, seTE, studlab,
               exclude = if (missing.exclude) NULL else exclude,
               ##
               sm = sm,
               level = level,
               level.comb = level.comb,
               comb.fixed = comb.fixed,
               comb.random = comb.random,
               ##
               hakn = hakn,
               method.tau = method.tau,
               tau.preset = tau.preset,
               TE.tau = TE.tau,
               tau.common = FALSE,
               ##
               prediction = prediction,
               level.predict = level.predict,
               ##
               method.bias = method.bias,
               ##
               backtransf = backtransf,
               title = title, complab = complab, outclab = outclab,
               label.e = label.e, label.c = label.c,
               label.left = label.left, label.right = label.right,
               ##
               keepdata = FALSE,
               warn = warn,
               ##
               control = control)
  ##
  if (by & tau.common) {
    ## Estimate common tau-squared across subgroups
    hcc <- hetcalc(TE, seTE, method.tau, TE.tau,
                   level.comb, byvar, control)
  }
  
  
  ##
  ##
  ## (9) Generate R object
  ##
  ##
  res <- list(n.e = n.e, mean.e = mean.e, sd.e = sd.e,
              n.c = n.c, mean.c = mean.c, sd.c = sd.c,
              pooledvar = pooledvar,
              method.smd = method.smd, sd.glass = sd.glass,
              exact.smd = exact.smd)
  ##
  ## Add meta-analysis results
  ## (after removing unneeded list elements)
  ##
  m$n.e <- NULL
  m$n.c <- NULL
  ##
  res <- c(res, m)
  ##
  ## Add data
  ##
  res$call <- match.call()
  ##
  if (keepdata) {
    res$data <- data
    if (!missing.subset)
      res$subset <- subset
  }
  ##
  class(res) <- c(fun, "meta")
  ##
  ## Add results from subgroup analysis
  ##
  if (by) {
    res$byvar <- byvar
    res$bylab <- bylab
    res$print.byvar <- print.byvar
    res$byseparator <- byseparator
    res$tau.common <- tau.common
    ##
    if (!tau.common) {
      res <- c(res, subgroup(res))
      res$tau.resid <- NA
    }
    else if (!is.null(tau.preset)) {
      res <- c(res, subgroup(res, tau.preset))
      res$tau.resid <- NA
    }
    else {
      res <- c(res, subgroup(res, hcc$tau))
      res$Q.w.random <- hcc$Q
      res$df.Q.w.random <- hcc$df.Q
      res$tau.resid <- hcc$tau
    }
    ##
    if (!tau.common || method.tau == "DL") {
      ci.H.resid <- calcH(res$Q.w.fixed, res$df.Q.w, level.comb)
      ##
      res$H.resid <- ci.H.resid$TE
      res$lower.H.resid <- ci.H.resid$lower
      res$upper.H.resid <- ci.H.resid$upper
    }
    else {
      res$H.resid <- hcc$H.resid
      res$lower.H.resid <- hcc$lower.H.resid
      res$upper.H.resid <- hcc$upper.H.resid
    }
    ##
    if (!tau.common || method.tau == "DL") {
      ci.I2.resid <- isquared(res$Q.w.fixed, res$df.Q.w, level.comb)
      ##
      res$I2.resid <- ci.I2.resid$TE
      res$lower.I2.resid <- ci.I2.resid$lower
      res$upper.I2.resid <- ci.I2.resid$upper
    }
    else {
      res$I2.resid <- hcc$I2.resid
      res$lower.I2.resid <- hcc$lower.I2.resid
      res$upper.I2.resid <- hcc$upper.I2.resid
    }
    ##
    res$event.e.w <- NULL
    res$event.c.w <- NULL
    res$event.w <- NULL
    res$n.w <- NULL
    res$time.e.w <- NULL
    res$time.c.w <- NULL
  }
  ##
  class(res) <- c(fun, "meta")
  
  
  res
}
