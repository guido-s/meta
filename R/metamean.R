#' Meta-analysis of single means
#' 
#' @description
#' Calculation of an overall mean from studies reporting a single mean
#' using the inverse variance method for pooling; inverse variance
#' weighting is used for pooling.
#' 
#' @param n Number of observations.
#' @param mean Estimated mean.
#' @param sd Standard deviation.
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
#' @param null.effect A numeric value specifying the effect under the
#'   null hypothesis.
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
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots for \code{sm = "MLN"}. If
#'   TRUE (default), results will be presented as means; otherwise
#'   logarithm of means will be shown.
#' @param title Title of meta-analysis / systematic review.
#' @param complab Comparison label.
#' @param outclab Outcome label.
#' @param sm A character string indicating which summary measure
#'   (\code{"MRAW"} or \code{"MLN"}) is to be used for pooling of
#'   studies.
#' @param byvar An optional vector containing grouping information
#'   (must be of same length as \code{n}).
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
#' Fixed effect and random effects meta-analysis of single means to
#' calculate an overall mean; inverse variance weighting is used for
#' pooling. The following transformations of means are implemented to
#' calculate an overall mean:
#' \itemize{
#' \item Raw, i.e. untransformed, means (\code{sm = "MRAW"}, default)
#' \item Log transformed means (\code{sm = "MLN"})
#' }
#' 
#' Note, you should use R function \code{\link{metacont}} to compare
#' means of pairwise comparisons instead of using \code{metamean} for
#' each treatment arm separately which will break randomisation in
#' randomised controlled trials.
#' 
#' Calculations are conducted on the log scale if \code{sm =
#' "ROM"}. Accordingly, list elements \code{TE}, \code{TE.fixed}, and
#' \code{TE.random} contain the logarithm of means. In printouts and
#' plots these values are back transformed if argument
#' \code{backtransf = TRUE}.
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
#' meta-analysis of an existing metamean object by only specifying
#' arguments which should be changed.
#' 
#' For the random effects, the method by Hartung and Knapp (2001) /
#' Knapp and Hartung (2003) is used to adjust test statistics and
#' confidence intervals if argument \code{hakn = TRUE}.
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
#' An object of class \code{c("metamean", "meta")} with corresponding
#' \code{print}, \code{summary}, and \code{forest} functions. The
#' object is a list containing the following components:
#' \item{n, mean, sd,}{As defined above.}
#' \item{studlab, exclude, sm, level, level.comb,}{As defined above.}
#' \item{comb.fixed, comb.random,}{As defined above.}
#' \item{hakn, method.tau, tau.preset, TE.tau, method.bias,}{As
#'   defined above.}
#' \item{tau.common, title, complab, outclab,}{As defined above.}
#' \item{byvar, bylab, print.byvar, byseparator, warn}{As defined
#'   above.}
#' \item{TE, seTE}{Estimated effect (mean or log mean) and standard
#'   error of individual studies.}
#' \item{lower, upper}{Lower and upper confidence interval limits for
#'   individual studies.}
#' \item{zval, pval}{z-value and p-value for test of overall effect
#'   for individual studies.}
#' \item{w.fixed, w.random}{Weight of individual studies (in fixed and
#'   random effects model).}
#' \item{TE.fixed, seTE.fixed}{Estimated overall effect (mean or log
#'   mean) and standard error (fixed effect model).}
#' \item{lower.fixed, upper.fixed}{Lower and upper confidence interval
#'   limits (fixed effect model).}
#' \item{zval.fixed, pval.fixed}{z-value and p-value for test of
#'   overall effect (fixed effect model).}
#' \item{TE.random, seTE.random}{Estimated overall effect (mean or log
#'   mean) and standard error (random effects model).}
#' \item{lower.random, upper.random}{Lower and upper confidence
#'   interval limits (random effects model).}
#' \item{zval.random, pval.random}{z-value or t-value and
#'   corresponding p-value for test of overall effect (random effects
#'   model).}
#' \item{prediction, level.predict}{As defined above.}
#' \item{seTE.predict}{Standard error utilised for prediction
#'   interval.}
#' \item{lower.predict, upper.predict}{Lower and upper limits of
#'   prediction interval.}
#' \item{k}{Number of studies combined in meta-analysis.}
#' \item{Q}{Heterogeneity statistic.}
#' \item{tau}{Square-root of between-study variance.}
#' \item{se.tau}{Standard error of square-root of between-study
#'   variance.}
#' \item{C}{Scaling factor utilised internally to calculate common
#'   tau-squared across subgroups.}
#' \item{method}{Pooling method: \code{"Inverse"}.}
#' \item{df.hakn}{Degrees of freedom for test of treatment effect for
#'   Hartung-Knapp method (only if \code{hakn = TRUE}).}
#' \item{bylevs}{Levels of grouping variable - if \code{byvar} is not
#'   missing.}
#' \item{TE.fixed.w, seTE.fixed.w}{Estimated effect and standard error
#'   in subgroups (fixed effect model) - if \code{byvar} is not
#'   missing.}
#' \item{lower.fixed.w, upper.fixed.w}{Lower and upper confidence
#'   interval limits in subgroups (fixed effect model) - if
#'   \code{byvar} is not missing.}
#' \item{zval.fixed.w, pval.fixed.w}{z-value and p-value for test of
#'   treatment effect in subgroups (fixed effect model) - if
#'   \code{byvar} is not missing.}
#' \item{TE.random.w, seTE.random.w}{Estimated effect and standard
#'   error in subgroups (random effects model) - if \code{byvar} is
#'   not missing.}
#' \item{lower.random.w, upper.random.w}{Lower and upper confidence
#'   interval limits in subgroups (random effects model) - if
#'   \code{byvar} is not missing.}
#' \item{zval.random.w, pval.random.w}{z-value or t-value and
#'   corresponding p-value for test of effect in subgroups (random
#'   effects model) - if \code{byvar} is not missing.}
#' \item{w.fixed.w, w.random.w}{Weight of subgroups (in fixed and
#'   random effects model) - if \code{byvar} is not missing.}
#' \item{df.hakn.w}{Degrees of freedom for test of effect for
#'   Hartung-Knapp method in subgroups - if \code{byvar} is not
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
#' \item{H.w}{Heterogeneity statistic H within subgroups - if
#'   \code{byvar} is not missing.}
#' \item{lower.H.w, upper.H.w}{Lower and upper confidence limti for
#'   heterogeneity statistic H within subgroups - if \code{byvar} is
#'   not missing.}
#' \item{I2.w}{Heterogeneity statistic I2 within subgroups - if
#'   \code{byvar} is not missing.}
#' \item{lower.I2.w, upper.I2.w}{Lower and upper confidence limti for
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
#' @seealso \code{\link{update.meta}}, \code{\link{metamean}},
#'   \code{\link{metagen}}
#' 
#' @references
#' DerSimonian R & Laird N (1986):
#' Meta-analysis in clinical trials.
#' \emph{Controlled Clinical Trials},
#' \bold{7}, 177--88
#' 
#' Hartung J & Knapp G (2001):
#' On tests of the overall treatment effect in meta-analysis with
#' normally distributed responses.
#' \emph{Statistics in Medicine},
#' \bold{20}, 1771--82
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
#' Viechtbauer W (2010):
#' Conducting Meta-Analyses in R with the Metafor Package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
#' 
#' @examples
#' m1 <- metamean(rep(100, 3), 1:3, rep(1, 3))
#' m1
#' 
#' m2 <- update(m1, sm = "MLN")
#' m2
#' 
#' # With test for overall mean equal to 2
#' #
#' update(m1, null.effect = 2)
#' update(m2, null.effect = 2)
#' 
#' # Print results without back-transformation
#' #
#' update(m1, backtransf = FALSE)
#' update(m2, backtransf = FALSE)
#' update(m1, null.effect = 2, backtransf = FALSE)
#' update(m2, null.effect = 2, backtransf = FALSE)
#' 
#' @export metamean


metamean <- function(n, mean, sd, studlab,
                     ##
                     data = NULL, subset = NULL, exclude = NULL,
                     ##
                     sm = gs("smmean"),
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
                     null.effect = NA,
                     ##
                     method.bias = gs("method.bias"),
                     ##
                     backtransf = gs("backtransf"),
                     title = gs("title"), complab = gs("complab"),
                     outclab = "",
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
  chknumeric(null.effect, single = TRUE)
  ##
  method.bias <- setchar(method.bias,
                         c("rank", "linreg", "mm", "count", "score", "peters"))
  ##
  chklogical(keepdata)
  ##
  ## Additional arguments / checks for metamean objects
  ##
  fun <- "metamean"
  sm <- setchar(sm, .settings$sm4mean)
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
  ## Catch 'n', 'mean', and 'sd' from data:
  ##
  n <- eval(mf[[match("n", names(mf))]],
            data, enclos = sys.frame(sys.parent()))
  chknull(n)
  k.All <- length(n)
  ##
  mean <- eval(mf[[match("mean", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  chknull(mean)
  ##
  sd <- eval(mf[[match("sd", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  chknull(sd)
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
  chklength(mean, k.All, fun)
  chklength(sd, k.All, fun)
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
      data <- data.frame(.n = n)
    else
      data$.n <- n
    ##
    data$.mean <- mean
    data$.sd <- sd
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
    n <- n[subset]
    mean <- mean[subset]
    sd <- sd[subset]
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
  k.all <- length(n)
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
  chknumeric(n)
  chknumeric(mean)
  chknumeric(sd)
  ##
  ## Recode integer as numeric:
  ##
  n    <- int2num(n)
  mean <- int2num(mean)
  sd   <- int2num(sd)
  
  
  ##
  ##
  ## (7) Calculate results for individual studies
  ##
  ##
  npn.n <- npn(n)
  ##
  if (any(npn.n) & warn)
    warning("Studies with non-positive sample size get no weight in meta-analysis.")
  ##
  if (sm == "MRAW") {
    TE <- ifelse(npn.n, NA, mean)
    ##
    seTE <- ifelse(npn.n, NA, sqrt(sd^2 / n))
    ##
    seTE[is.na(TE)] <- NA
    ##
    transf.null.effect <- null.effect
  }
  ##
  else if (sm == "MLN") {
    npn.mean <- npn(mean)
    ##
    if (any(npn.mean) & warn)
      warning("Studies with negative or zero mean get no weight in meta-analysis.")

    TE <- ifelse(npn.n | npn.mean, NA, log(mean))
    ##
    seTE <- ifelse(npn.n | npn.mean, NA, sqrt(sd^2 / (n * mean^2)))
    ##
    seTE[is.na(TE)] <- NA
    ##
    transf.null.effect <- log(null.effect)
  }
  ##
  ## Studies with non-positive variance get zero weight in meta-analysis
  ##
  sel <- sd <= 0
  ##
  if (any(sel, na.rm = TRUE) & warn)
    warning("Studies with non-positive standard deviation get no weight in meta-analysis.")
  ##
  seTE[sel] <- NA
  
  
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
               null.effect = transf.null.effect,
               ##
               method.bias = method.bias,
               ##
               backtransf = backtransf,
               title = title, complab = complab, outclab = outclab,
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
  res <- list(n = n, mean = mean, sd = sd)
  ##
  ## Add meta-analysis results
  ## (after removing unneeded list elements)
  ##
  m$n.e <- NULL
  m$n.c <- NULL
  m$label.e <- ""
  m$label.c <- ""
  m$label.left <- ""
  m$label.right <- ""
  ##
  res <- c(res, m)
  res$null.effect <- null.effect
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
