#' Update a meta-analysis object
#' 
#' @description
#' Update an existing meta-analysis object.
#' 
#' @param object An object of class \code{meta}.
#' @param data Dataset.
#' @param subset Subset.
#' @param studlab Study label.
#' @param exclude An optional vector specifying studies to exclude
#'   from meta-analysis, however, to include in printouts and forest
#'   plots.
#' @param method A character string indicating which method is to be
#'   used for pooling of studies; see \code{\link{metabin}} and
#'   \code{\link{metainc}} function for admissible values.
#' @param sm A character string indicating which summary measure is
#'   used for pooling.
#' @param incr Either a numerical value or vector which can be added
#'   to each cell frequency for studies with a zero cell count or the
#'   character string \code{"TA"} which stands for treatment arm
#'   continuity correction.
#' @param allincr A logical indicating if \code{incr} is added to each
#'   cell frequency of all studies if at least one study has a zero
#'   cell count. If FALSE (default), \code{incr} is added only to each
#'   cell frequency of studies with a zero cell count.
#' @param addincr A logical indicating if \code{incr} is added to each
#'   cell frequency of all studies irrespective of zero cell counts.
#' @param allstudies A logical indicating if studies with zero or all
#'   events in both groups are to be included in the meta-analysis
#'   (applies only if \code{sm} is equal to \code{"RR"} or
#'   \code{"OR"}).
#' @param MH.exact A logical indicating if \code{incr} is not to be
#'   added to all cell frequencies for studies with a zero cell count
#'   to calculate the pooled estimate based on the Mantel-Haenszel
#'   method.
#' @param RR.Cochrane A logical indicating if 2*\code{incr} instead of
#'   1*\code{incr} is to be added to \code{n.e} and \code{n.c} in the
#'   calculation of the risk ratio (i.e., \code{sm="RR"}) for studies
#'   with a zero cell. This is used in RevMan 5, the program for
#'   preparing and maintaining Cochrane reviews.
#' @param Q.Cochrane A logical indicating if the Mantel-Haenszel
#'   estimate is used in the calculation of the heterogeneity
#'   statistic Q which is implemented in RevMan 5, the program for
#'   preparing and maintaining Cochrane reviews.
#' @param model.glmm A character string indicating which GLMM model
#'   should be used.
#' @param level The level used to calculate confidence intervals for
#'   individual studies.
#' @param level.comb The level used to calculate confidence intervals
#'   for pooled estimates.
#' @param comb.fixed A logical indicating whether a fixed effect
#'   meta-analysis should be conducted.
#' @param comb.random A logical indicating whether a random effects
#'   meta-analysis should be conducted.
#' @param overall A logical indicating whether overall summaries
#'   should be reported. This argument is useful in a meta-analysis
#'   with subgroups if overall results should not be reported.
#' @param overall.hetstat A logical value indicating whether to print
#'   heterogeneity measures for overall treatment comparisons. This
#'   argument is useful in a meta-analysis with subgroups if
#'   heterogeneity statistics should only be printed on subgroup
#'   level.
#' @param hakn A logical indicating whether the method by Hartung and
#'   Knapp should be used to adjust test statistics and confidence
#'   intervals.
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2} and its
#'   square root \eqn{\tau}. Either \code{"DL"}, \code{"PM"},
#'   \code{"REML"}, \code{"ML"}, \code{"HS"}, \code{"SJ"},
#'   \code{"HE"}, or \code{"EB"}, can be abbreviated. See function
#'   \code{\link{metagen}}.
#' @param method.tau.ci A character string indicating which method is
#'   used to estimate the confidence interval of \eqn{\tau^2} and
#'   \eqn{\tau}. Either \code{"QP"}, \code{"BJ"}, or \code{"J"}, can
#'   be abbreviated.
#' @param tau.preset Prespecified value for the square root of the
#'   between-study variance \eqn{\tau^2}.
#' @param TE.tau Overall treatment effect used to estimate the
#'   between-study variance \eqn{\tau^2}.
#' @param tau.common A logical indicating whether tau-squared should
#'   be the same across subgroups.
#' @param prediction A logical indicating whether a prediction
#'   interval should be printed.
#' @param level.predict The level used to calculate prediction
#'   interval for a new study.
#' @param null.effect A numeric value specifying the effect under the
#'   null hypothesis.
#' @param method.bias A character string indicating which test for
#'   funnel plot asymmetry is to be used. Either \code{"rank"},
#'   \code{"linreg"}, \code{"mm"}, \code{"count"}, \code{"score"}, or
#'   \code{"peters"}, can be abbreviated. See function
#'   \code{\link{metabias}}
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If \code{backtransf =
#'   TRUE}, results for \code{sm = "OR"} are printed as odds ratios
#'   rather than log odds ratios and results for \code{sm = "ZCOR"}
#'   are printed as correlations rather than Fisher's z transformed
#'   correlations, for example.
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
#' @param title Title of meta-analysis / systematic review.
#' @param complab Comparison label.
#' @param outclab Outcome label.
#' @param label.e Label for experimental group.
#' @param label.c Label for control group.
#' @param label.left Graph label on left side of forest plot.
#' @param label.right Graph label on right side of forest plot.
#' @param n.e Number of observations in experimental group. (only for
#'   metagen object)
#' @param n.c Number of observations in control group. (only for
#'   metagen object)
#' @param pooledvar A logical indicating if a pooled variance should
#'   be used for the mean difference (only for metacont object with
#'   \code{sm = "MD"}).
#' @param method.smd A character string indicating which method is
#'   used to estimate the standardised mean difference (only for
#'   metacont object with \code{sm = "SMD"}). Either \code{"Hedges"}
#'   for Hedges' g (default), \code{"Cohen"} for Cohen's d, or
#'   \code{"Glass"} for Glass' delta, can be abbreviated.
#' @param sd.glass A character string indicating which standard
#'   deviation is used in the denominator for Glass' method to
#'   estimate the standardised mean difference (only for metacont
#'   object with \code{sm = "SMD"}). Either \code{"control"} using the
#'   standard deviation in the control group (\code{sd.c}) or
#'   \code{"experimental"} using the standard deviation in the
#'   experimental group (\code{sd.e}), can be abbreviated.
#' @param exact.smd A logical indicating whether exact formulae should
#'   be used in estimation of the standardised mean difference and its
#'   standard error.
#' @param method.ci A character string indicating which method is used
#'   to calculate confidence intervals for individual studies. Either
#'   \code{"CP"}, \code{"WS"}, \code{"WSCC"}, \code{"AC"},
#'   \code{"SA"},, \code{"SACC"}, or \code{"NAsm"}, can be
#'   abbreviated. See function \code{\link{metaprop}}.
#' @param byvar An optional vector containing grouping information
#'   (must be of same length as \code{event.e}).
#' @param bylab A character string with a label for the grouping
#'   variable.
#' @param print.byvar A logical indicating whether the name of the
#'   grouping variable should be printed in front of the group labels.
#' @param byseparator A character string defining the separator
#'   between label and levels of grouping variable.
#' @param print.CMH A logical indicating whether result of the
#'   Cochran-Mantel-Haenszel test for overall effect should be
#'   printed.
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in meta object.
#' @param left A logical indicating whether studies are supposed to be
#'   missing on the left or right side of the funnel plot. If NULL,
#'   the linear regression test for funnel plot symmetry (i.e.,
#'   function \code{metabias(..., method = "linreg")}) is used to
#'   determine whether studies are missing on the left or right side.
#' @param ma.fixed A logical indicating whether a fixed effect or
#'   random effects model is used to estimate the number of missing
#'   studies.
#' @param type A character indicating which method is used to estimate
#'   the number of missing studies. Either \code{"L"} or \code{"R"}.
#' @param n.iter.max Maximum number of iterations to estimate number
#'   of missing studies.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if \code{incr} is added to studies with zero cell
#'   frequencies).
#' @param control An optional list to control the iterative process to
#'   estimate the between-study variance \eqn{\tau^2}. This argument
#'   is passed on to \code{\link[metafor]{rma.uni}} or
#'   \code{\link[metafor]{rma.glmm}}, respectively.
#' @param \dots Additional arguments (ignored at the moment).
#' 
#' @details
#' Wrapper function to update an existing meta-analysis object which
#' was created with R function \code{\link{metabin}},
#' \code{\link{metacont}}, \code{\link{metacor}},
#' \code{\link{metagen}}, \code{\link{metainc}},
#' \code{\link{metamean}}, \code{\link{metaprop}}, or
#' \code{\link{metarate}}. More details on function arguments are
#' available in help files of respective R functions
#' 
#' This function can also be used for objects of class 'trimfill',
#' 'metacum', and 'metainf'.
#' 
#' @return
#' An object of class \code{"meta"} and \code{"metabin"},
#' \code{"metacont"}, \code{"metacor"}, \code{"metainc"},
#' \code{"metagen"}, \code{"metamean"}, \code{"metaprop"}, or
#' \code{"metarate"}.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}, \code{\link{metacont}},
#'   \code{\link{metacor}}, \code{\link{metagen}},
#'   \code{\link{metainc}}, \code{\link{metamean}},
#'   \code{\link{metaprop}}, \code{\link{metarate}}
#' 
#' @examples
#' data(Fleiss93cont)
#' m1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c,
#'                data = Fleiss93cont, sm = "SMD", studlab = study)
#' m1
#' 
#' # Change summary measure (from 'SMD' to 'MD')
#' #
#' update(m1, sm = "MD")
#' 
#' # Restrict analysis to subset of studies
#' #
#' update(m1, subset = 1:2)
#' 
#' # Use different levels for confidence intervals
#' #
#' m2 <- update(m1, level = 0.66, level.comb = 0.99)
#' print(m2, digits = 2)
#' forest(m2)
#' 
#' @export
#' @export update.meta


update.meta <- function(object, 
                        data = object$data,
                        subset = object$subset,
                        studlab = object$data$.studlab,
                        exclude = object$data$.exclude,
                        method = object$method,
                        sm = object$sm,
                        incr,
                        allincr = object$allincr,
                        addincr = object$addincr,
                        allstudies = object$allstudies,
                        MH.exact = object$MH.exact,
                        RR.Cochrane = object$RR.Cochrane,
                        Q.Cochrane = object$Q.Cochrane,
                        model.glmm = object$model.glmm,
                        level = object$level,
                        level.comb = object$level.comb,
                        comb.fixed = object$comb.fixed,
                        comb.random = object$comb.random,
                        overall = object$overall,
                        overall.hetstat = object$overall.hetstat,
                        hakn = object$hakn,
                        method.tau = object$method.tau,
                        method.tau.ci = object$method.tau.ci,
                        tau.preset = object$tau.preset,
                        TE.tau = object$TE.tau,
                        tau.common = object$tau.common,
                        prediction = object$prediction,
                        level.predict = object$level.predict,
                        null.effect = object$null.effect,
                        method.bias = object$method.bias,
                        ##
                        backtransf = object$backtransf,
                        pscale = object$pscale,
                        irscale = object$irscale,
                        irunit = object$irunit,
                        title = object$title,
                        complab = object$complab,
                        outclab = object$outclab,
                        label.e = object$label.e,
                        label.c = object$label.c,
                        label.left = object$label.left,
                        label.right = object$label.right,
                        n.e = object$n.e,
                        n.c = object$n.c,
                        pooledvar = object$pooledvar,
                        method.smd = object$method.smd,
                        sd.glass = object$sd.glass,
                        exact.smd = object$exact.smd,
                        method.ci = object$method.ci,
                        byvar = object$byvar,
                        bylab = object$bylab,
                        print.byvar = object$print.byvar,
                        byseparator = object$byseparator,
                        print.CMH = object$print.CMH,
                        keepdata = TRUE,
                        ##
                        left = object$left,
                        ma.fixed = object$ma.fixed,
                        type = object$type,
                        n.iter.max = object$n.iter.max,
                        ##
                        warn = FALSE,
                        ##
                        control = object$control,
                        ...) {
  
  
  ##
  ##
  ## (1) Check for meta object
  ##
  ##
  chkclass(object, "meta")
  ##
  metabin  <- inherits(object, "metabin")
  metacont <- inherits(object, "metacont")
  metacor  <- inherits(object, "metacor")
  metagen  <- inherits(object, "metagen")
  metainc  <- inherits(object, "metainc")
  metamean <- inherits(object, "metamean")
  metaprop <- inherits(object, "metaprop")
  metarate <- inherits(object, "metarate")
  
  
  ##
  ##
  ## (2) Replace missing arguments with defaults
  ##
  ##
  replacemiss <- function(x, replace) {
    ##
    xnam <- deparse(substitute(x))
    ##
    if (is.null(x))
      if (missing(replace))
        res <- gs(xnam)
      else
        res <- replace
    else
      res <- x
    ##
    res
  }
  ##
  comb.fixed <- replacemiss(comb.fixed)
  comb.random <- replacemiss(comb.random)
  overall <- replaceNULL(overall, comb.fixed | comb.random)
  overall.hetstat <- replaceNULL(overall.hetstat, comb.fixed | comb.random)
  ##
  RR.Cochrane <- replacemiss(RR.Cochrane, object$RR.cochrane)
  Q.Cochrane <- replacemiss(Q.Cochrane, TRUE)
  if (Q.Cochrane & (!(sm %in% c("OR", "RR", "RD")) | method.tau != "DL"))
    Q.Cochrane <- FALSE
  ##
  model.glmm <- replacemiss(model.glmm)
  ##
  level <- replacemiss(level)
  level.comb <- replacemiss(level.comb)
  ##
  hakn <- replacemiss(hakn)
  method.tau <- replacemiss(method.tau)
  method.tau.ci <- replacemiss(method.tau.ci)
  tau.preset <- replacemiss(tau.preset, NULL)
  TE.tau <- replacemiss(TE.tau, NULL)
  null.effect <- replacemiss(null.effect, NA)
  method.bias <- replacemiss(method.bias)
  ##
  backtransf <- replacemiss(backtransf)
  label.left <- replacemiss(label.left)
  label.right <- replacemiss(label.right)
  ##
  tau.common <- replacemiss(tau.common)
  level.predict <- replacemiss(level.predict)
  prediction <- replacemiss(prediction)
  ##
  pscale  <- replacemiss(pscale, 1)
  irscale <- replacemiss(irscale, 1)
  irunit   <- replacemiss(irunit, 1)
  ##
  title <- replacemiss(title)
  complab <- replacemiss(complab)
  outclab <- replacemiss(outclab, "")
  label.e <- replacemiss(label.e)
  label.c <- replacemiss(label.c)
  ##
  print.byvar <- replacemiss(print.byvar)
  byseparator <- replacemiss(byseparator)
  ##
  warn <- replacemiss(warn)
  ##
  if (!backtransf & pscale != 1 & !is.untransformed(sm)) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  if (!backtransf & irscale != 1 & !is.untransformed(sm)) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  
  
  ##
  ##
  ## (3) Update trim-and-fill object
  ##
  ##
  if (inherits(object, "trimfill")) {
    ##
    rmfilled <- function(x) {
      ##
      if (!is.null(object[[x]]))
        res <- object[[x]][!object$trimfill]
      else
        res <- NULL
      ##
      res
    }
    ##
    tfnames <- c("TE", "seTE",
                 "studlab",
                 "n.e", "n.c",
                 "event.e", "event.c",
                 "mean.e", "mean.c", "sd.e", "sd.c",
                 "n", "event", "cor")
    ##
    for (i in tfnames)
      object[[i]] <- rmfilled(i)
    ##
    oldclass <- object$class.x
    ##
    res <- trimfill(object,
                    left = left, ma.fixed = ma.fixed,
                    type = type, n.iter.max = n.iter.max,
                    level = level, level.comb = level.comb,
                    comb.fixed = comb.fixed, comb.random = comb.random,
                    hakn = hakn,
                    method.tau = method.tau, method.tau.ci = method.tau.ci,
                    prediction = prediction, level.predict = level.predict,
                    silent = TRUE,
                    ...)
    ##
    res$call.object <- object$call
    res$call <- match.call()
    res$class.x <- oldclass
    ##
    return(res)
  }
  
  
  ##
  ##
  ## (4) Update metacum or metainf object
  ##
  ##
  if (inherits(object, "metacum") | inherits(object, "metainf")) {
    ##
    res <- object
    ##
    res$comb.fixed <- ifelse(res$pooled == "fixed", TRUE, FALSE)
    res$comb.random <- ifelse(res$pooled == "random", TRUE, FALSE)
    ##
    res$call.object <- object$call
    res$call <- match.call()
    res$version <- packageDescription("meta")$Version
    ##
    return(res)
  }
  
  
  ##
  ##
  ## (5) Prepare older meta object
  ##
  ##
  if (!(!is.null(object$version) &&
        as.numeric(unlist(strsplit(object$version, "-"))[1]) >= 3.2)) {
    ##
    ## Changes for meta objects with version < 3.2
    ##
    object$subset <- NULL
    ##
    object$data <- data.frame(.studlab = object$studlab,
                              .exclude = rep_len(FALSE, length(object$studlab)))
    ##
    if (!is.null(object$byvar))
      object$data$.byvar <- object$byvar
    ##
    if (metabin) {
      object$data$.event.e <- object$event.e
      object$data$.n.e <- object$n.e
      object$data$.event.c <- object$event.c
      object$data$.n.c <- object$n.c
    }
    ##
    if (metacont) {
      object$data$.n.e <- object$n.e
      object$data$.mean.e <- object$mean.e
      object$data$.sd.e <- object$sd.e
      object$data$.n.c <- object$n.c
      object$data$.mean.c <- object$mean.c
      object$data$.sd.c <- object$sd.c
    }
    ##
    if (metacor) {
      object$data$.cor <- object$cor
      object$data$.n <- object$n
    }
    ##
    if (metagen) {
      object$data$.TE <- object$TE
      object$data$.seTE <- object$seTE
    }
    ##
    if (metaprop) {
      object$data$.event <- object$event
      object$data$.n <- object$n
    }
  }
  ##
  if (!(!is.null(object$version) &&
        as.numeric(unlist(strsplit(object$version, "-"))[1]) >= 4.8)) {
    ##
    ## Changes for meta objects with version < 4.8
    ##
    if (metabin | metainc | metaprop | metarate)
      object$data$.incr <- object$incr
    ##
    if (metabin | metainc)
      if (object$method == "MH")
        object$k.MH <- sum(object$w.fixed > 0)
      else
        object$k.MH <- NA
  }
  ##
  if (is.null(object$data)) {
    warning("Necessary data not available. Please, recreate ",
            "meta-analysis object without option 'keepdata = FALSE'.")
    return(invisible(NULL))
  }
  ##
  missing.subset  <- missing(subset)
  missing.incr    <- missing(incr)
  missing.byvar   <- missing(byvar)
  missing.studlab <- missing(studlab)
  missing.exclude <- missing(exclude)
  ##
  mf <- match.call()
  ##
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  ##
  incr <- eval(mf[[match("incr", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  ##
  byvar <- eval(mf[[match("byvar", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  ##
  if (!missing.byvar) {
    byvar.name <- as.character(mf[[match("byvar", names(mf))]])
    if (length(byvar.name) > 1 & byvar.name[1] == "$")
      byvar.name <- byvar.name[length(byvar.name)]
    if (length(byvar.name) > 1)
      byvar.name <- "byvar"
    ##
    bylab <- if (!missing(bylab) && !is.null(bylab)) bylab else byvar.name
    ##
    data$.byvar <- byvar
  }
  ##
  studlab <- eval(mf[[match("studlab", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  ##
  exclude <- eval(mf[[match("exclude", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  ##
  if (missing.subset) {
    if (!is.null(object$subset))
      subset <- object$subset
    else if (isCol(object$data, ".subset"))
      subset <- object$data$.subset
  }
  ##
  if (missing.incr) {
    if (isCol(object$data, ".incr"))
      incr <- object$data$.incr
    else
      incr <- gs("incr")
  }
  ##
  if (missing.byvar & isCol(object$data, ".byvar"))
    byvar <- object$data$.byvar
  ##
  if (missing.studlab & isCol(object$data, ".studlab"))
    studlab <- object$data$.studlab
  ##
  if (missing.exclude & isCol(object$data, ".exclude"))
    exclude <- object$data$.exclude
  ##
  if (method == "GLMM")
    if (metabin & !missing(sm) & sm != "OR")
      warning("Summary measure 'sm = \"OR\" used as 'method = \"GLMM\".")
    else if (metainc & !missing(sm) & sm != "IRR")
      warning("Summary measure 'sm = \"IRR\" used as 'method = \"GLMM\".")
    else if (metaprop & !missing(sm) & sm != "PLOGIT")
      warning("Summary measure 'sm = \"PLOGIT\" used as 'method = \"GLMM\".")
    else if (metarate & !missing(sm) & sm != "IRLN")
      warning("Summary measure 'sm = \"IRLN\" used as 'method = \"GLMM\".")
  
  
  ##
  ##
  ## (6) Update meta object
  ##
  ##
  if (metabin) {
    sm <- setchar(sm, .settings$sm4bin)
    method <- setchar(method, c("Inverse", "MH", "Peto", "GLMM"))
    if (sm == "ASD") {
      if (!missing.incr)
        warning("Note, no continuity correction considered for arcsine difference (sm = \"ASD\").")
      incr <- 0
      object$data$.incr <- 0
    }
    ##
    if (method == "Peto") {
      if (!missing.incr)
        warning("Note, no continuity correction considered for method = \"Peto\".")
      incr <- 0
      object$data$.incr <- 0
    }
    ##
    m <- metabin(event.e = object$data$.event.e,
                 n.e = object$data$.n.e,
                 event.c = object$data$.event.c,
                 n.c = object$data$.n.c,
                 ##
                 studlab = studlab,
                 exclude = exclude,
                 ##
                 data = data, subset = subset,
                 ##
                 method = method,
                 sm = ifelse(method == "GLMM", "OR", sm),
                 incr = incr,
                 allincr = allincr, addincr = addincr, allstudies = allstudies,
                 MH.exact = MH.exact, RR.Cochrane = RR.Cochrane,
                 Q.Cochrane = Q.Cochrane, model.glmm = model.glmm,
                 ##
                 level = level, level.comb = level.comb,
                 comb.fixed = comb.fixed, comb.random = comb.random,
                 overall = overall, overall.hetstat = overall.hetstat,
                 ##
                 hakn = hakn,
                 method.tau = ifelse(method == "GLMM", "ML", method.tau),
                 method.tau.ci = method.tau.ci,
                 tau.preset = tau.preset, TE.tau = TE.tau,
                 tau.common = tau.common,
                 ##
                 prediction = prediction, level.predict = level.predict,
                 ##
                 method.bias = method.bias,
                 ##
                 backtransf = backtransf, pscale = pscale,
                 title = title, complab = complab, outclab = outclab,
                 label.e = label.e, label.c = label.c,
                 label.right = label.right, label.left = label.left,
                 ##
                 byvar = byvar, bylab = bylab, print.byvar = print.byvar,
                 byseparator = byseparator,
                 print.CMH = print.CMH,
                 ##
                 keepdata = keepdata,
                 warn = warn,
                 ##
                 control = control,
                 ...)
  }
  ##
  if (metacont)
    m <- metacont(n.e = object$data$.n.e,
                  mean.e = object$data$.mean.e,
                  sd.e = object$data$.sd.e,
                  n.c = object$data$.n.c,
                  mean.c = object$data$.mean.c,
                  sd.c = object$data$.sd.c,
                  ##
                  studlab = studlab,
                  exclude = exclude,
                  ##
                  data = data, subset = subset,
                  ##
                  sm = sm, pooledvar = pooledvar,
                  method.smd = method.smd, sd.glass = sd.glass, exact.smd = exact.smd,
                  ##
                  level = level, level.comb = level.comb,
                  comb.fixed = comb.fixed, comb.random = comb.random,
                  overall = overall, overall.hetstat = overall.hetstat,
                  ##
                  hakn = hakn,
                  method.tau = method.tau, method.tau.ci = method.tau.ci,
                  tau.preset = tau.preset, TE.tau = TE.tau,
                  tau.common = tau.common,
                  ##
                  prediction = prediction, level.predict = level.predict,
                  ##
                  method.bias = method.bias,
                  ##
                  title = title, complab = complab, outclab = outclab,
                  label.e = label.e, label.c = label.c,
                  label.right = label.right, label.left = label.left,
                  ##
                  byvar = byvar, bylab = bylab, print.byvar = print.byvar,
                  byseparator = byseparator,
                  ##
                  keepdata = keepdata,
                  warn = warn,
                  ##
                  control = control)
  ##
  if (metacor)
    m <- metacor(cor = object$data$.cor,
                 n = object$data$.n,
                 ##
                 studlab = studlab,
                 exclude = exclude,
                 ##
                 data = data, subset = subset,
                 ##
                 sm = sm,
                 ##
                 level = level, level.comb = level.comb,
                 comb.fixed = comb.fixed, comb.random = comb.random,
                 overall = overall, overall.hetstat = overall.hetstat,
                 ##
                 hakn = hakn,
                 method.tau = method.tau, method.tau.ci = method.tau.ci,
                 tau.preset = tau.preset, TE.tau = TE.tau,
                 tau.common = tau.common,
                 ##
                 prediction = prediction, level.predict = level.predict,
                 ##
                 null.effect = null.effect,
                 ##
                 method.bias = method.bias,
                 ##
                 backtransf = backtransf,
                 title = title, complab = complab, outclab = outclab,
                 byvar = byvar, bylab = bylab, print.byvar = print.byvar,
                 byseparator = byseparator,
                 ##
                 keepdata = keepdata,
                 ##
                 control = control)
  ##
  if (metagen) {
    data.m <- data
    add.e <- FALSE
    add.c <- FALSE
    ##
    if ("n.e" %in% names(data)) {
      add.e <- TRUE
      data.m <- data.m[, names(data.m) != "n.e"]
    }
    if ("n.c" %in% names(data)) {
      add.c <- TRUE
      data.m <- data.m[, names(data.m) != "n.c"]
    }
    ##
    m <- metagen(TE = object$data$.TE,
                 seTE = object$data$.seTE,
                 ##
                 studlab = studlab,
                 exclude = exclude,
                 ##
                 data = data.m, subset = subset,
                 ##
                 sm = sm,
                 ##
                 level = level, level.comb = level.comb,
                 comb.fixed = comb.fixed, comb.random = comb.random,
                 overall = overall, overall.hetstat = overall.hetstat,
                 ##
                 hakn = hakn,
                 method.tau = method.tau, method.tau.ci = method.tau.ci,
                 tau.preset = tau.preset, TE.tau = TE.tau,
                 tau.common = tau.common,
                 ##
                 prediction = prediction, level.predict = level.predict,
                 ##
                 method.bias = method.bias,
                 ##
                 n.e = n.e, n.c = n.c,
                 ##
                 backtransf = backtransf,
                 title = title, complab = complab, outclab = outclab,
                 label.e = label.e, label.c = label.c,
                 label.right = label.right, label.left = label.left,
                 ##
                 byvar = byvar, bylab = bylab, print.byvar = print.byvar,
                 byseparator = byseparator,
                 ##
                 keepdata = keepdata,
                 warn = warn,
                 ##
                 control = control)
    if (add.e)
      m$data$n.e <- data$n.e
    if (add.c)
      m$data$n.c <- data$n.c
    if (add.e | add.c)
      m$data <- m$data[, names(data)]
  }
  ##
  if (metainc) {
    data.m <- data
    add.e <- FALSE
    add.c <- FALSE
    ##
    if ("n.e" %in% names(data)) {
      add.e <- TRUE
      data.m <- data.m[, names(data.m) != "n.e"]
    }
    if ("n.c" %in% names(data)) {
      add.c <- TRUE
      data.m <- data.m[, names(data.m) != "n.c"]
    }
    ##
    m <- metainc(event.e = object$data$.event.e,
                 time.e = object$data$.time.e,
                 event.c = object$data$.event.c,
                 time.c = object$data$.time.c,
                 ##
                 studlab = studlab,
                 exclude = exclude,
                 ##
                 data = data, subset = subset,
                 ##
                 method = method,
                 sm = ifelse(method == "GLMM", "IRR", sm),
                 incr = incr,
                 allincr = allincr, addincr = addincr,
                 model.glmm = model.glmm,
                 ##
                 level = level, level.comb = level.comb,
                 comb.fixed = comb.fixed, comb.random = comb.random,
                 overall = overall, overall.hetstat = overall.hetstat,
                 ##
                 hakn = hakn,
                 method.tau = ifelse(method == "GLMM", "ML", method.tau),
                 method.tau.ci = method.tau.ci,
                 tau.preset = tau.preset, TE.tau = TE.tau,
                 tau.common = tau.common,
                 ##
                 prediction = prediction, level.predict = level.predict,
                 ##
                 method.bias = method.bias,
                 ##
                 n.e = n.e, n.c = n.c,
                 ##
                 backtransf = backtransf, irscale = irscale, irunit = irunit,
                 title = title, complab = complab, outclab = outclab,
                 label.e = label.e, label.c = label.c,
                 label.right = label.right, label.left = label.left,
                 ##
                 byvar = byvar, bylab = bylab, print.byvar = print.byvar,
                 byseparator = byseparator,
                 ##
                 keepdata = keepdata,
                 warn = warn,
                 ##
                 control = control,
                 ...)
    if (add.e)
      m$data$n.e <- data$n.e
    if (add.c)
      m$data$n.c <- data$n.c
    if (add.e | add.c)
      m$data <- m$data[, names(data)]
  }
  ##
  if (metamean)
    m <- metamean(n = object$data$.n,
                  mean = object$data$.mean,
                  sd = object$data$.sd,
                  ##
                  studlab = studlab,
                  exclude = exclude,
                  ##
                  data = data, subset = subset,
                  ##
                  sm = sm,
                  ##
                  level = level, level.comb = level.comb,
                  comb.fixed = comb.fixed, comb.random = comb.random,
                  overall = overall, overall.hetstat = overall.hetstat,
                  ##
                  hakn = hakn,
                  method.tau = method.tau, method.tau.ci = method.tau.ci,
                  tau.preset = tau.preset, TE.tau = TE.tau,
                  tau.common = tau.common,
                  ##
                  prediction = prediction, level.predict = level.predict,
                  ##
                  null.effect = null.effect,
                  ##
                  method.bias = method.bias,
                  ##
                  backtransf = backtransf,
                  title = title, complab = complab, outclab = outclab,
                  ##
                  byvar = byvar, bylab = bylab, print.byvar = print.byvar,
                  byseparator = byseparator,
                  ##
                  keepdata = keepdata,
                  warn = warn,
                  ##
                  control = control)
  ##
  if (metaprop)
    m <- metaprop(event = object$data$.event,
                  n = object$data$.n,
                  ##
                  studlab = studlab,
                  exclude = exclude,
                  ##
                  data = data, subset = subset, method = method,
                  ##
                  sm = ifelse(method == "GLMM", "PLOGIT", sm),
                  incr = incr,
                  allincr = allincr, addincr = addincr,
                  method.ci = ifelse(is.null(method.ci), "CP", method.ci),
                  ##
                  level = level, level.comb = level.comb,
                  comb.fixed = comb.fixed, comb.random = comb.random,
                  overall = overall, overall.hetstat = overall.hetstat,
                  ##
                  hakn = hakn,
                  method.tau = ifelse(method == "GLMM", "ML", method.tau),
                  method.tau.ci = method.tau.ci,
                  tau.preset = tau.preset, TE.tau = TE.tau,
                  tau.common = tau.common,
                  ##
                  prediction = prediction, level.predict = level.predict,
                  ##
                  null.effect = null.effect,
                  ##
                  method.bias = method.bias,
                  ##
                  backtransf = backtransf, pscale = pscale,
                  title = title, complab = complab, outclab = outclab,
                  byvar = byvar, bylab = bylab, print.byvar = print.byvar,
                  byseparator = byseparator,
                  ##
                  keepdata = keepdata,
                  warn = warn,
                  ##
                  control = control,
                  ...)
  ##
  if (metarate)
    m <- metarate(event = object$data$.event,
                  time = object$data$.time,
                  ##
                  studlab = studlab,
                  exclude = exclude,
                  ##
                  data = data, subset = subset, method = method,
                  ##
                  sm = ifelse(method == "GLMM", "IRLN", sm),
                  incr = incr,
                  allincr = allincr, addincr = addincr,
                  ##
                  level = level, level.comb = level.comb,
                  comb.fixed = comb.fixed, comb.random = comb.random,
                  overall = overall, overall.hetstat = overall.hetstat,
                  ##
                  hakn = hakn,
                  method.tau = ifelse(method == "GLMM", "ML", method.tau),
                  method.tau.ci = method.tau.ci,
                  tau.preset = tau.preset, TE.tau = TE.tau,
                  tau.common = tau.common,
                  ##
                  prediction = prediction, level.predict = level.predict,
                  ##
                  null.effect = null.effect,
                  ##
                  method.bias = method.bias,
                  ##
                  backtransf = backtransf, irscale = irscale, irunit = irunit,
                  title = title, complab = complab, outclab = outclab,
                  byvar = byvar, bylab = bylab, print.byvar = print.byvar,
                  byseparator = byseparator,
                  ##
                  keepdata = keepdata,
                  warn = warn,
                  ##
                  control = control,
                  ...)
  ##  
  m$call.object <- object$call
  m$call <- match.call()
  
  
  m
}
