#' Print and change default meta-analysis settings in R package \bold{meta}
#' 
#' @description
#' Print and change default settings to conduct and print or plot
#' meta-analyses in R package \bold{meta}. The following general
#' settings are available: \emph{Review Manager 5}, \emph{Journal of
#' the American Medical Association}.
#' 
#' @param ... Arguments to change default settings.
#' @param quietly A logical indicating whether information on settings
#'   should be printed.
#' 
#' @details
#' This function can be used to define defaults for several arguments
#' (i.e., assignments using \code{\link{gs}}) of the following R
#' functions: \code{\link{metabin}}, \code{\link{metacont}},
#' \code{\link{metacor}}, \code{\link{metacr}}, \code{\link{metagen}},
#' \code{\link{metainc}}, \code{\link{metaprop}},
#' \code{\link{metarate}}
#' 
#' Furthermore, some of these settings are considered to print
#' meta-analysis results and to produce forest plots.
#' 
#' The function can be used to either change individual settings (see
#' Examples) or use one of the following general settings:
#' \itemize{
#' \item \code{settings.meta("RevMan5")}
#' \item \code{settings.meta("BMJ")}
#' \item \code{settings.meta("JAMA")}
#' \item \code{settings.meta("IQWiG5")}
#' \item \code{settings.meta("IQWiG6")}
#' \item \code{settings.meta("geneexpr")}
#' \item \code{settings.meta("IVhet")}
#' \item \code{settings.meta("meta4")}
#' \item \code{settings.meta("meta7")}
#' }
#'
#' The first command can be used to reproduce meta-analyses from
#' Cochrane reviews conducted with \emph{Review Manager 5} (RevMan 5)
#' and specifies to use a RevMan 5 layout in forest plots.
#'
#' The second command can be used to generate forest plots in BMJ layout.
#'
#' The third command can be used to generate forest plots following
#' instructions for authors of the \emph{Journal of the American
#' Medical Association}. Study labels according to JAMA guidelines can be
#' generated using \code{\link{labels.meta}}.
#'
#' The next commands implement the recommendations of the Institute
#' for Quality and Efficiency in Health Care, Germany (IQWiG)
#' accordinging to General Methods 5 and 6, respectively
#' (\url{https://www.iqwig.de/en/about-us/methods/methods-paper/}).
#'
#' The setting \code{"geneexpr"} can be used to print p-values in
#' scientific notation and to suppress the calculation of confidence
#' intervals for the between-study variance.
#'
#' The setting \code{"IVhet"} can be used for the inverse variance
#' heterogeneity model (Henmi & Copas, 2010; Doi et al., 2015).
#'
#' The last settings use the default settings of R package
#' \bold{meta}, version 4 and 7.0-0, respectively, or below.
#' 
#' RevMan 5 settings, in detail:
#' \tabular{lll}{
#' \bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
#' \code{method.random.ci} \tab "classic" \tab only available method
#'   in RevMan 5 \cr
#' \code{method.tau} \tab "DL" \tab only available method in RevMan 5
#'   \cr
#' \code{method.I2} \tab "Q" \tab only available method in RevMan 5 \cr
#' \code{tau.common} \tab FALSE \tab common between-study variance in
#'   subgroups \cr
#' \code{MH.exact} \tab FALSE \tab exact Mantel-Haenszel method \cr
#' \code{RR.Cochrane} \tab TRUE \tab calculation of risk ratios \cr
#' \code{Q.Cochrane} \tab TRUE \tab calculation of heterogeneity
#'   statistic \cr
#' \code{exact.smd} \tab FALSE \tab exact formulae for Hedges' g and
#'   Cohen's d \cr
#' \code{layout} \tab "RevMan5" \tab layout for forest plots \cr
#' \code{prediction} \tab FALSE \tab no prediction interval \cr
#' \code{test.overall} \tab TRUE \tab print information on test of
#'   overall effect \cr
#' \code{test.subgroup} \tab TRUE \tab print information on test for
#'   subgroup differences \cr
#' \code{test.effect.subgroup} \tab TRUE \tab print information on
#'   test for effect in subgroups \cr
#' \code{forest.I2} \tab TRUE \tab show heterogeneity statistic I2 in
#'   forest plots \cr
#' \code{forest.tau2} \tab TRUE \tab show between-study heterogeneity \cr
#'   \tab \tab variance in forest plots \cr
#' \code{forest.tau} \tab FALSE \tab do not show between-study heterogeneity \cr
#'   \tab \tab standard deviation in forest plots \cr
#' \code{forest.Q} \tab TRUE \tab show heterogeneity statistic Q in
#'   forest plots
#'   \cr
#' \code{forest.pval.Q} \tab TRUE \tab show p-value of test for heterogeneity
#'   in forest plots
#'   \cr
#' \code{forest.Rb} \tab FALSE \tab do not show heterogeneity statistic Rb in
#'   forest plots
#'   \cr
#' \code{digits.tau2} \tab 3 \tab number of digits for tau-squared \cr
#' \code{digits.tau} \tab 4 \tab number of digits for square root of
#'   tau-squared \cr
#' \code{digits.I2} \tab 0 \tab number of digits for I-squared measure
#'   \cr
#' \code{CIbracket}, \tab "[" \tab \cr
#' \code{CIseparator} \tab ", " \tab print confidence intervals as
#'   "\code{[., .]}" \cr
#' \code{header.line}, \tab TRUE \tab print header line
#' }
#'
#' BMJ settings:
#' \tabular{lll}{
#' \bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
#' \code{layout} \tab "BMJ" \tab layout for forest plots \cr
#' \code{test.overall} \tab TRUE \tab print information on test of
#'   overall effect \cr
#' \code{test.subgroup} \tab FALSE \tab print information on test for
#'   subgroup differences \cr
#' \code{test.effect.subgroup} \tab FALSE \tab print information on
#'   test for effect in subgroups \cr
#' \code{forest.I2} \tab TRUE \tab show heterogeneity statistic I2 in
#'   forest plots
#'   \cr
#' \code{forest.tau2} \tab TRUE \tab show between-study heterogeneity \cr
#'   \tab \tab variance in forest plots \cr
#' \code{forest.tau} \tab FALSE \tab do not show between-study heterogeneity \cr
#'   \tab \tab standard deviation in forest plots \cr
#' \code{forest.Q} \tab TRUE \tab show heterogeneity statistic Q in
#'   forest plots
#'   \cr
#' \code{forest.pval.Q} \tab TRUE \tab show p-value of test for heterogeneity
#'   in forest plots
#'   \cr
#' \code{forest.Rb} \tab FALSE \tab do not show heterogeneity statistic Rb in
#'   forest plots
#'   \cr
#' \code{digits.I2} \tab 0 \tab number of digits for I-squared measure
#'   \cr
#' \code{digits.pval} \tab 2 \tab number of digits for p-values \cr
#' \code{CIbracket}, \tab "(" \tab \cr
#' \code{CIseparator} \tab " to " \tab print confidence intervals as
#'   "\code{(. to .)}" \cr
#' \code{hetlab}, \tab \tab "Test for heterogeneity: " \cr
#' \code{header.line}, \tab TRUE \tab print header line
#' }
#'
#' JAMA settings:
#' \tabular{lll}{
#' \bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
#' \code{layout} \tab "JAMA" \tab layout for forest plots \cr
#' \code{test.overall} \tab TRUE \tab print information on test of
#'   overall effect \cr
#' \code{test.subgroup} \tab FALSE \tab print information on test for
#'   subgroup differences \cr
#' \code{test.effect.subgroup} \tab FALSE \tab print information on
#'   test for effect in subgroups \cr
#' \code{forest.I2} \tab TRUE \tab show heterogeneity statistic I2 in
#'   forest plots
#'   \cr
#' \code{forest.tau2} \tab FALSE \tab do not show between-study heterogeneity \cr
#'   \tab \tab variance in forest plots \cr
#' \code{forest.tau} \tab FALSE \tab do not show between-study heterogeneity \cr
#'   \tab \tab standard deviation in forest plots \cr
#' \code{forest.Q} \tab TRUE \tab show heterogeneity statistic Q in
#'   forest plots
#'   \cr
#' \code{forest.pval.Q} \tab TRUE \tab show p-value of test for heterogeneity
#'   in forest plots
#'   \cr
#' \code{forest.Rb} \tab FALSE \tab do not show heterogeneity statistic Rb in
#'   forest plots
#'   \cr
#' \code{digits.I2} \tab 0 \tab number of digits for I-squared measure
#'   \cr
#' \code{digits.pval} \tab 3 \tab number of digits for p-values \cr
#' \code{CIbracket}, \tab "(" \tab \cr
#' \code{CIseparator} \tab "-" \tab print confidence intervals as
#'   "\code{(.-.)}" \cr
#' \code{zero.pval}, \tab FALSE \tab print p-values with leading zero
#' \cr
#' \code{JAMA.pval}, \tab TRUE \tab round p-values to three digits
#'   (for 0.001 < p \eqn{\le} 0.01) \cr
#'   \tab \tab or two digits (p > 0.01) \cr
#' \code{header.line}, \tab TRUE \tab print header line
#' }
#' 
#' IQWiG, General Methods 5 settings:
#' \tabular{lll}{
#' \bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
#' \code{method.random.ci} \tab "HK" \tab Hartung-Knapp method \cr
#' \code{prediction} \tab TRUE \tab Prediction interval \cr
#' }
#' 
#' IQWiG, General Methods 6 settings:
#' \tabular{lll}{
#' \bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
#' \code{method.random.ci} \tab "HK" \tab Hartung-Knapp method \cr
#' \code{adhoc.hakn.ci} \tab "IQWiG6" \tab \emph{ad hoc} variance correction \cr
#' \code{method.tau} \tab "PM" \tab Paule-Mandel estimator for
#'   between-study variance \cr
#' \code{prediction} \tab TRUE \tab Prediction interval \cr
#' }
#' 
#' Settings for gene expression data:
#' \tabular{lll}{
#' \bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
#' \code{scientific.pval} \tab TRUE \tab Scientific notation for p-values \cr
#' \code{method.tau.ci} \tab FALSE \tab
#'   no confidence interval for between-study \cr
#'  \tab \tab heterogeneity variance \cr
#' }
#'
#' IVhet settings:
#' \tabular{lll}{
#' \bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
#' \code{method.common.ci} \tab "IVhet" \tab inverse variance heterogeneity \cr
#' \code{text.common} \tab "IVhet model" \tab  \cr
#' \code{text.w.common} \tab "IVhet" \tab  \cr
#' }
#' 
#' Settings for \bold{meta}, version 4 or below:
#' \tabular{lll}{
#' \bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
#' \code{method.tau} \tab "DL" \tab DerSimonian-Laird estimator \cr
#' \code{method.I2} \tab "Q" \tab Use Q to calculate I-squared \cr
#' \code{method.predict} \tab "HTS" \tab prediction interval with \emph{k-2}
#'   degrees \cr
#'   \tab \tab of freedom \cr
#' \code{exact.smd} \tab FALSE \tab Use exact formula for standardised mean \cr
#'   \tab \tab difference (White and Thomas, 2005) \cr
#' \code{text.common} \tab "Fixed effect model" \tab \cr
#' \code{text.w.common} \tab "fixed" \tab \cr
#' \code{warn.deprecated} \tab FALSE \tab Do not print warnings for deprecated
#'   \cr
#'   \tab \tab arguments
#' }
#' 
#' Settings for \bold{meta}, version 7.0-0 or below:
#' \tabular{lll}{
#' \bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
#' \code{method.tau} \tab "REML" \tab REML estimator \cr
#' \code{method.I2} \tab "Q" \tab Use Q to calculate I-squared \cr
#' \code{method.predict} \tab "HTS" \tab prediction interval with \emph{k-2}
#'   degrees \cr
#'   \tab \tab of freedom \cr
#' \code{exact.smd} \tab TRUE \tab Use exact formula for standardised mean \cr
#'   \tab \tab difference (White and Thomas, 2005) \cr
#' \code{text.common} \tab "Common effect model" \tab \cr
#' \code{text.w.common} \tab "common" \tab \cr
#' \code{warn.deprecated} \tab FALSE \tab Do not print warnings for deprecated
#'   \cr
#'   \tab \tab arguments
#' }
#' 
#' A list of all arguments with current settings is printed using the
#' command \code{settings.meta()}.
#' 
#' In order to reset all settings of R package \bold{meta} the command
#' \code{settings.meta(reset = TRUE)} can be used.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{gs}}, \code{\link{forest.meta}},
#'   \code{\link{print.meta}}, \code{\link{labels.meta}}
#'
#' @references
#' Doi SAR, Barendregt JJ, Khan S, Thalib L, Williams GM (2015):
#' Advances in the meta-analysis of heterogeneous clinical trials I:
#' The inverse variance heterogeneity model.
#' \emph{Contemporary Clinical Trials},
#' \bold{45}, 130--8
#' 
#' Henmi M & Copas JB (2010):
#' Confidence intervals for random effects meta-analysis and
#' robustness to publication bias.
#' \emph{Statistics in Medicine},
#' \bold{29}, 2969--83
#' 
#' White IR, Thomas J (2005):
#' Standardized mean differences in individually-randomized and
#' cluster-randomized trials, with applications to meta-analysis.
#' \emph{Clinical Trials},
#' \bold{2}, 141--51
#' 
#' @examples
#' # Get listing of current settings
#' #
#' settings.meta()
#' 
#' # Meta-analyses using default settings
#' #
#' metabin(10, 20, 15, 20)
#' metaprop(4, 20)
#' metabin(10, 20, 15, 20, sm = "RD")
#' metaprop(4, 20, sm = "PLN")
#'
#' # Change summary measure for R functions metabin and metaprop
#' # and store old settings
#' #
#' oldset <- settings.meta(smbin = "RD", smprop = "PLN")
#' #
#' metabin(10, 20, 15, 20)
#' metaprop(4, 20)
#'
#' # Use old settings
#' #
#' settings.meta(oldset)
#' 
#' # Change level used to calculate confidence intervals
#' # (99%-CI for studies, 99.9%-CI for pooled effects)
#' #
#' metagen(1:3, 2:4 / 10, sm = "MD")
#' settings.meta(level = 0.99, level.ma = 0.999)
#' metagen(1:3, 2:4 / 10, sm = "MD")
#' 
#' # Always print a prediction interval
#' #
#' settings.meta(prediction = TRUE)
#' metagen(1:3, 2:4 / 10, sm = "MD")
#' metagen(4:6, 4:2 / 10, sm = "MD")
#' 
#' # Try to set unknown argument results in a warning
#' #
#' try(settings.meta(unknownarg = TRUE))
#' 
#' # Reset to default settings of R package meta
#' #
#' settings.meta("reset")
#' metabin(10, 20, 15, 20)
#' metaprop(4, 20)
#' metagen(1:3, 2:4 / 10, sm = "MD")
#' 
#' # Do not back transform results (e.g. print log odds ratios instead
#' # of odds ratios, print transformed correlations / proportions
#' # instead of correlations / proportions)
#' #
#' settings.meta(backtransf = FALSE)
#' metabin(10, 20, 15, 20)
#' metaprop(4, 20)
#' metacor(c(0.85, 0.7, 0.95), c(20, 40, 10))
#' 
#' # Forest plot using RevMan 5 style
#' #
#' settings.meta("RevMan5")
#' forest(metagen(1:3, 2:4 / 10, sm = "MD", common = FALSE),
#'   label.left = "Favours A", label.right = "Favours B",
#'   colgap.studlab = "2cm", colgap.forest.left = "0.2cm")
#' 
#' # Forest plot using JAMA style
#' #
#' settings.meta("JAMA")
#' forest(metagen(1:3, 2:4 / 10, sm = "MD", common = FALSE),
#'   label.left = "Favours A", label.right = "Favours B",
#'   colgap.studlab = "2cm", colgap.forest.left = "0.2cm")
#'
#' # Use slightly different layout for confidence intervals
#' # (especially useful if upper confidence limit can be negative)
#' #
#' settings.meta(CIseparator = " - ")
#' forest(metagen(-(1:3), 2:4 / 10, sm = "MD", common = FALSE),
#'   label.left = "Favours A", label.right = "Favours B",
#'   colgap.studlab = "2cm", colgap.forest.left = "0.2cm")
#' 
#' # Use old settings
#' #
#' settings.meta(oldset)
#' 
#' @export settings.meta


settings.meta <- function(..., quietly = TRUE) {
  
  ##
  ## Check argument
  ##
  missing.quietly <- missing(quietly)
  chklogical(quietly)
  
  
  ##
  ## Save object with current settings
  ##
  oldset <- .settings
  oldset$argslist <- oldset$argslist.internal <- NULL
  
  
  ##
  ## Set internal variables
  ##
  settings <- c("BMJ", "JAMA", "RevMan5", 
                "IQWiG5", "IQWiG6",
                "geneexpr", "IVhet",
                "meta4", "meta7")
  layouts <- c(settings[1:2], "meta")
  ##
  print.settings <- FALSE
  reset.settings <- FALSE
  specific.settings <- FALSE
  ##
  args  <- list(...)
  
  
  ##
  ## Print settings if no argument is provided
  ##
  if (length(args) == 0) {
    if (missing.quietly || !quietly)
      settings.meta("print", quietly = FALSE)
    return(invisible(oldset))
  }
  
  
  ##
  ## Check whether first argument is a list. In this case only use
  ## this list as input.
  ##
  warn.depr <- TRUE
  if (length(args) > 0 && is.list(args[[1]])) {
    if (!is.null(names(args))) {
      print(names(args))
      warning("Additional arguments ignored as first argument is a list.",
              call. = FALSE)
    }
    warn.depr <- FALSE
    args <- args[[1]]
  }
  
  
  ##
  ## Unnamed first (and only) argument must be character string or a
  ## logical
  ##
  if (length(args) == 1 & is.null(names(args))) {
    if (is.character(unlist(args)))
      action <- setchar(unlist(args), c(settings, "reset", "print"),
                        stop.at.error = FALSE)
    else
      action <- unlist(args)
    ##
    if (is.null(action))
      stop("First argument can be one of the following character strings:",
           "\n reset, print, ", paste(settings, collapse = ", "),
           call. = FALSE)
    else if (action == "reset")
      settings.meta(reset = TRUE, quietly = quietly)
    else if (action == "print" | (is.logical(action) && action))
      settings.meta(print = TRUE, quietly = FALSE)
    else
      settings.meta(setting = action, quietly = quietly)
    ##
    return(invisible(oldset))
  }
  ##
  else if (length(args) > 1 & names(args)[1] == "") {
    if (is.character(unlist(args[[1]])))
      action <- setchar(unlist(args[[1]]), c(settings, "reset", "print"),
                        stop.at.error = FALSE)
    else
      action <- unlist(args[[1]])
    ##
    if (is.null(action))
      stop("First argument can be one of the following character strings:",
           "\n reset, print, ", paste(settings, collapse = ", "),
           call. = FALSE)
    else if (action == "reset")
      settings.meta(reset = TRUE, quietly = quietly)
    else if (action == "print")
      settings.meta(print = TRUE, quietly = FALSE)
    else
      settings.meta(setting = action, quietly = quietly)
  }
  
  
  ##
  ## Check and warn about deprecated arguments
  ##
  names.all <- names(args)
  ##
  if (warn.depr) {
    chkdeprecated(names.all, "level.ma", "level.comb")
    chkdeprecated(names.all, "common", "comb.fixed")
    chkdeprecated(names.all, "common", "fixed")
    chkdeprecated(names.all, "random", "comb.random")
    ##
    chkdeprecated(names.all, "method.random.ci", "hakn")
    chkdeprecated(names.all, "adhoc.hakn.ci", "adhoc.hakn")
    ##
    chkdeprecated(names.all, "digits.stat", "digits.zval")
    chkdeprecated(names.all, "print.subgroup.name", "print.byvar")
    chkdeprecated(names.all, "sep.subgroup", "byseparator")
    ##
    chkdeprecated(names.all, "method.incr", "addincr")
    chkdeprecated(names.all, "method.incr", "allincr")
    #
    chkdeprecated(names.all, "cid.below.null", "lower.equi")
    chkdeprecated(names.all, "cid.above.null", "upper.equi")
    chkdeprecated(names.all, "lty.cid", "lty.equi")
    chkdeprecated(names.all, "col.cid", "col.equi")
  }
  
  
  #
  # Check for new argument names
  #
  newargs <- vector("character", 0)
  #
  for (i in seq_along(names.all)) {
    if (is.null(
      setchar(names.all[i],
              c(.settings$argslist.meta, "reset", "print", "setting", ""),
              stop.at.error = FALSE)))
      newargs <- c(newargs, names.all[i])
  }
  #
  for (i in seq_len(length(newargs))) {
    setOption("argslist", c(.settings$argslist, newargs[i]))
    setOption(newargs[i], args[[newargs[i]]])
  }
  #
  names <- names.all[!(names.all %in% .settings$argslist.internal)]
  
  
  ##
  ## Check argument names
  ##
  if (length(names) != length(unique(names)))
    stop("Arguments must be unique.")
  
  
  ##
  ## Determine whether to print, reset or use specific settings
  ##
  if (any(names == "print") && args[["print"]]) {
    print.settings <- TRUE
    quietly <- FALSE
  }
  if (any(names == "reset") && args[["reset"]])
    reset.settings <- TRUE
  if (any(names == "setting")) {
    setting <- setchar(args[["setting"]], settings)
    specific.settings <- TRUE
  }
  
  
  ##
  ## Reset settings
  ##
  if (reset.settings) {
    if (!quietly)
      cat("\n** Reset all meta-analysis settings (R package meta). **\n\n")
    ##
    ## General settings
    ##
    setOption("level", 0.95)
    setOption("level.ma", 0.95)
    ##
    setOption("common", TRUE)
    setOption("fixed", TRUE)
    setOption("comb.fixed", TRUE)
    setOption("method.common.ci", "classic")
    setOption("random", TRUE)
    setOption("comb.random", TRUE)
    setOption("method.random.ci", "classic")
    setOption("hakn", FALSE)
    setOption("adhoc.hakn", "")
    setOption("adhoc.hakn.ci", "")
    setOption("adhoc.hakn.pi", "")
    setOption("method.tau", "REML")
    setOption("method.tau.ci", NULL)
    setOption("level.hetstat", 0.95)
    setOption("tau.common", FALSE)
    setOption("method.I2", "Q")
    setOption("prediction", FALSE)
    setOption("level.predict", 0.95)
    setOption("method.predict", "V")
    setOption("test.subgroup", TRUE)
    setOption("prediction.subgroup", FALSE)
    setOption("method.bias", "Egger")
    setOption("tool.rob", NULL)
    setOption("overall.hetstat", NULL)
    setOption("text.common", "Common effect model")
    setOption("text.random", "Random effects model")
    setOption("text.predict", "Prediction interval")
    setOption("text.w.common", "common")
    setOption("text.w.random", "random")
    setOption("title", "")
    setOption("complab", "")
    setOption("CIbracket", "[")
    setOption("CIseparator", "; ")
    setOption("CIlower.blank", TRUE)
    setOption("CIupper.blank", TRUE)
    setOption("print.subgroup.name", TRUE)
    setOption("print.byvar", TRUE)
    setOption("sep.subgroup", " = ")
    setOption("byseparator", " = ")
    setOption("keepdata", TRUE)
    setOption("keeprma", FALSE)
    setOption("warn", TRUE)
    setOption("warn.deprecated", TRUE)
    setOption("transf", TRUE)
    setOption("backtransf", TRUE)
    setOption("digits", 4)
    setOption("digits.mean", 2)
    setOption("digits.sd", 4)
    setOption("digits.se", 4)
    setOption("digits.stat", 2)
    setOption("digits.zval", 2)
    setOption("digits.Q", 2)
    setOption("digits.tau2", 4)
    setOption("digits.tau", 4)
    setOption("digits.H", 2)
    setOption("digits.I2", 1)
    setOption("digits.prop", 4)
    setOption("digits.weight", 1)
    setOption("digits.pval", 4)
    setOption("digits.pval.Q", 4)
    setOption("scientific.pval", FALSE)
    setOption("big.mark", "")
    setOption("zero.pval", TRUE)
    setOption("JAMA.pval", FALSE)
    setOption("digits.df", 4)
    setOption("digits.cid", 4)
    ##
    setOption("details", TRUE)
    ##
    setOption("print.tau2", TRUE)
    setOption("print.tau2.ci", TRUE)
    setOption("print.tau", TRUE)
    setOption("print.tau.ci", TRUE)
    setOption("print.I2", TRUE)
    setOption("print.I2.ci", TRUE)
    setOption("print.H", TRUE)
    setOption("print.Rb", FALSE)
    ##
    setOption("text.tau2", "tau^2")
    setOption("text.tau", "tau")
    setOption("text.I2", "I^2")
    setOption("text.Rb", "Rb")
    ##
    setOption("print.Q", TRUE)
    ##
    ## R function metabin
    ##
    setOption("smbin", "RR")
    setOption("method", "MH")
    setOption("incr", 0.5)
    setOption("method.incr", "only0")
    setOption("allincr", FALSE)
    setOption("addincr", FALSE)
    setOption("allstudies", FALSE)
    setOption("MH.exact", FALSE)
    setOption("RR.Cochrane", FALSE)
    setOption("Q.Cochrane", TRUE)
    setOption("model.glmm", "UM.FS")
    setOption("print.CMH", FALSE)
    ##
    ## R function metacont
    ##
    setOption("smcont", "MD")
    setOption("pooledvar", FALSE)
    setOption("method.smd", "Hedges")
    setOption("sd.glass", "control")
    setOption("exact.smd", TRUE)
    setOption("method.ci.cont", "z")
    ##
    ## R function metaprop
    ##
    setOption("smprop", "PLOGIT")
    setOption("method.ci.prop", "CP")
    ##
    ## R function metarate
    ##
    setOption("smrate", "IRLN")
    setOption("method.ci.rate", "NAsm")
    ##
    ## Other meta-analysis functions
    ##
    setOption("smcor", "ZCOR")
    setOption("sminc", "IRR")
    setOption("smmean", "MRAW")
    ##
    ## R functions comparing two treatments
    ##
    setOption("label.e", "Experimental")
    setOption("label.c", "Control")
    setOption("label.left", "")
    setOption("label.right", "")
    ##
    ## R function forest.meta
    ##
    setOption("layout", "meta")
    setOption("forest.details", FALSE)
    setOption("test.overall", FALSE)
    setOption("test.effect.subgroup", FALSE)
    setOption("digits.forest", 2)
    setOption("digits.TE.forest", 4)
    ##
    setOption("lty.common", 2)
    setOption("lty.random", 3)
    setOption("col.common", "black")
    setOption("col.random", "black")
    ##
    setOption("sort.subgroup", FALSE)
    ##
    setOption("pooled.events", FALSE)
    setOption("pooled.times", FALSE)
    setOption("study.results", TRUE)
    ##
    setOption("cid", NA)
    setOption("cid.below.null", NA)
    setOption("cid.above.null", NA)
    setOption("lty.cid", 1)
    setOption("col.cid", "blue")
    setOption("fill.cid", "transparent")
    setOption("cid.pooled.only", FALSE)
    #
    setOption("fill", "transparent")
    setOption("fill.equi", "transparent")
    ##
    setOption("leftcols", NULL)
    setOption("rightcols", NULL)
    setOption("leftlabs", NULL)
    setOption("rightlabs", NULL)
    ##
    setOption("label.e.attach", NULL)
    setOption("label.c.attach", NULL)
    ##
    setOption("bottom.lr", TRUE)
    ##
    setOption("lab.NA", ".")
    setOption("lab.NA.effect", NULL)
    setOption("lab.NA.weight", ".")
    ##
    setOption("lwd", 1)
    setOption("lwd.square", 1)
    setOption("lwd.diamond", 1)
    ##
    setOption("arrow.type", "open")
    setOption("arrow.length", 0.05)
    ##
    setOption("type.study", "square")
    setOption("type.common", "diamond")
    ##
    setOption("col.study", "black")
    setOption("col.square", "gray")
    setOption("col.square.lines", "gray")
    setOption("col.circle", "royalblue")
    setOption("col.inside", "white")
    setOption("col.diamond", "gray")
    setOption("col.diamond.lines", "black")
    setOption("col.predict", "red")
    setOption("col.predict.lines", "black")
    setOption("col.subgroup", "black")
    setOption("col.label.right", "black")
    setOption("col.label.left", "black")
    ##
    setOption("col.lines", "black")
    setOption("col.label", "black")
    ##
    setOption("hetlab", "Heterogeneity: ")
    setOption("resid.hetstat", NULL)
    setOption("resid.hetlab", "Residual heterogeneity: ")
    ##
    setOption("forest.I2", NULL)
    setOption("forest.I2.ci", FALSE)
    setOption("forest.tau2", NULL)
    setOption("forest.tau2.ci", FALSE)
    setOption("forest.tau", FALSE)
    setOption("forest.tau.ci", FALSE)
    setOption("forest.Q", FALSE)
    setOption("forest.pval.Q", NULL)
    setOption("forest.Rb", FALSE)
    setOption("forest.Rb.ci", FALSE)
    ##
    setOption("text.subgroup.nohet", "not applicable")
    ##
    setOption("LRT", FALSE)
    ##
    setOption("forest.stat", TRUE)
    setOption("forest.Q.subgroup", TRUE)
    ##
    setOption("header.line", FALSE)
    ##
    setOption("fontsize", 12)
    setOption("fontfamily", NULL)
    setOption("fs.common", NULL)
    setOption("fs.random", NULL)
    setOption("fs.predict", NULL)
    setOption("fs.common.labels", NULL)
    setOption("fs.random.labels", NULL)
    setOption("fs.predict.labels", NULL)
    setOption("fs.hetstat", NULL)
    setOption("fs.test.overall", NULL)
    setOption("fs.test.subgroup", NULL)
    setOption("fs.test.effect.subgroup", NULL)
    setOption("fs.addline", NULL)
    ##
    setOption("ff.heading", "bold")
    setOption("ff.common", NULL)
    setOption("ff.random", NULL)
    setOption("ff.predict", NULL)
    setOption("ff.common.labels", NULL)
    setOption("ff.random.labels", NULL)
    setOption("ff.predict.labels", NULL)
    setOption("ff.study", "plain")
    setOption("ff.hetstat", NULL)
    setOption("ff.test.overall", NULL)
    setOption("ff.test.subgroup", NULL)
    setOption("ff.test.effect.subgroup", NULL)
    setOption("ff.addline", NULL)
    setOption("ff.axis", "plain")
    setOption("ff.smlab", "bold")
    setOption("ff.xlab", "plain")
    setOption("ff.lr", "plain")
    ##
    setOption("colgap", "2mm")
    setOption("colgap.forest", "2mm")
    ##
    setOption("width", NULL)
    ##
    setOption("calcwidth.predict", FALSE)
    setOption("calcwidth.hetstat", FALSE)
    setOption("calcwidth.tests", FALSE)
    setOption("calcwidth.subgroup", FALSE)
    setOption("calcwidth.addline", FALSE)
    ##
    setOption("just.studlab", "left")
    setOption("just.addcols", "center")
    ##
    setOption("spacing", 1)
    setOption("addrow", NULL)
    setOption("addrow.overall", NULL)
    setOption("addrow.subgroups", NULL)
    setOption("addrows.below.overall", NULL)
  }
  
  
  ##
  ## Specific settings
  ##
  if (specific.settings) {
    #
    if (setting == "BMJ") {
      specificSettings(
        args = c("layout", "test.overall",
                 "test.subgroup", "test.effect.subgroup",
                 #
                 "forest.I2", "forest.tau2", "forest.tau",
                 "forest.Q", "forest.pval.Q", "forest.Rb",
                 #
                 "digits.I2", "digits.tau2", "digits.tau", "digits.pval",
                 "digits.forest", "digits.mean", "digits.sd",
                 #
                 "CIbracket", "CIseparator",
                 #
                 "colgap.forest",
                 #
                 "col.common", "col.random", "col.subgroup",
                 "col.study", "col.square", "col.square.lines",
                 "col.diamond", "col.diamond.lines",
                 #
                 "col.lines",
                 #
                 "lwd.square", "lwd.diamond",
                 #
                 "arrow.type",
                 #
                 "ff.lr",
                 "zero.pval", "JAMA.pval",
                 #
                 "hetlab", "header.line"),
        new = list("BMJ",
                   replaceNULL(args[["test.overall"]], TRUE),
                   replaceNULL(args[["test.subgroup"]], TRUE),
                   replaceNULL(args[["test.effect.subgroup"]], TRUE),
                   #
                   replaceNULL(args[["forest.I2"]], TRUE),
                   replaceNULL(args[["forest.tau2"]], TRUE),
                   replaceNULL(args[["forest.tau"]], FALSE),
                   replaceNULL(args[["forest.Q"]], TRUE),
                   replaceNULL(args[["forest.pval.Q"]], TRUE),
                   replaceNULL(args[["forest.Rb"]], FALSE),
                   #
                   replaceNULL(args[["digits.I2"]], 0),
                   replaceNULL(args[["digits.tau2"]], 2),
                   replaceNULL(args[["digits.tau"]], 2),
                   replaceNULL(args[["digits.pval"]], 2),
                   replaceNULL(args[["digits.forest"]], 2),
                   replaceNULL(args[["digits.mean"]], 1),
                   replaceNULL(args[["digits.sd"]], 1),
                   #
                   replaceNULL(args[["CIbracket"]], "("),
                   replaceNULL(args[["CIseparator"]], " to "),
                   #
                   replaceNULL(args[["colgap.forest"]], "5mm"),
                   #
                   replaceNULL(args[["col.common"]], "#6b58a6"),
                   replaceNULL(args[["col.random"]], "#6b58a6"),
                   replaceNULL(args[["col.subgroup"]], "black"),
                   replaceNULL(args[["col.study"]], "#6b58a6"),
                   replaceNULL(args[["col.square"]], "#6b58a6"),
                   replaceNULL(args[["col.square.lines"]], "white"),
                   replaceNULL(args[["col.diamond"]], "#6b58a6"),
                   replaceNULL(args[["col.diamond.lines"]], "white"),
                   #
                   replaceNULL(args[["col.lines"]], "#a7a9ac"),
                   #
                   replaceNULL(args[["lwd.square"]], 0.5),
                   replaceNULL(args[["lwd.diamond"]], 0.5),
                   #
                   replaceNULL(args[["arrow.type"]], "closed"),
                   #
                   replaceNULL(args[["ff.lr"]], "bold"),
                   replaceNULL(args[["zero.pval"]], TRUE),
                   replaceNULL(args[["JAMA.pval"]], TRUE),
                   #
                   replaceNULL(args[["hetlab"]], "Test for heterogeneity: "),
                   replaceNULL(args[["header.line"]], TRUE)),
        setting = "BMJ settings",
        quietly = quietly)
    }
    #
    else if (setting == "JAMA") {
      specificSettings(
        args = c("layout", "test.overall",
                 "test.subgroup", "test.effect.subgroup",
                 #
                 "forest.I2", "forest.tau2", "forest.tau",
                 "forest.Q", "forest.pval.Q", "forest.Rb",
                 #
                 "digits.I2", "digits.pval",
                 "CIbracket", "CIseparator",
                 "zero.pval", "JAMA.pval",
                 "hetlab", "header.line"),
        new = list("JAMA",
                   replaceNULL(args[["test.overall"]], TRUE),
                   replaceNULL(args[["test.subgroup"]], FALSE),
                   replaceNULL(args[["test.effect.subgroup"]], FALSE),
                   #
                   replaceNULL(args[["forest.tau2"]], FALSE),
                   replaceNULL(args[["forest.tau"]], FALSE),
                   replaceNULL(args[["forest.I2"]], TRUE),
                   replaceNULL(args[["forest.Q"]], TRUE),
                   replaceNULL(args[["forest.pval.Q"]], TRUE),
                   replaceNULL(args[["forest.Rb"]], FALSE),
                   #
                   replaceNULL(args[["digits.I2"]], 0),
                   replaceNULL(args[["digits.pval"]], 3),
                   replaceNULL(args[["CIbracket"]], "("),
                   replaceNULL(args[["CIseparator"]], "-"),
                   replaceNULL(args[["zero.pval"]], FALSE),
                   replaceNULL(args[["JAMA.pval"]], TRUE),
                   replaceNULL(args[["hetlab"]], "Heterogeneity: "),
                   replaceNULL(args[["header.line"]], TRUE)),
        setting = "JAMA settings",
        quietly = quietly)
    }
    #
    else if (setting == "RevMan5") {
      specificSettings(
        args = c("method.random.ci", "method.tau", "method.I2",
                 "tau.common",
                 "MH.exact", "RR.Cochrane", "Q.Cochrane",
                 "exact.smd",
                 "layout", "prediction", "test.overall",
                 "test.subgroup", "test.effect.subgroup",
                 "col.subgroup",
                 #
                 "forest.I2", "forest.tau2", "forest.tau",
                 "forest.Q", "forest.pval.Q", "forest.Rb",
                 #
                 "digits", "digits.I2", "digits.tau2", "digits.tau",
                 #
                 "CIbracket", "CIseparator",
                 "zero.pval", "JAMA.pval",
                 "text.common", "text.w.common",
                 "hetlab", "header.line"),
        new = list(replaceNULL(args[["method.random.ci"]], "classic"),
                   replaceNULL(args[["method.tau"]], "DL"),
                   replaceNULL(args[["method.I2"]], "Q"),
                   replaceNULL(args[["tau.common"]], FALSE),
                   replaceNULL(args[["MH.exact"]], FALSE),
                   replaceNULL(args[["RR.Cochrane"]], TRUE),
                   replaceNULL(args[["Q.Cochrane"]], TRUE),
                   replaceNULL(args[["exact.smd"]], FALSE), "RevMan5",
                   replaceNULL(args[["prediction"]], FALSE),
                   replaceNULL(args[["test.overall"]], TRUE),
                   replaceNULL(args[["test.subgroup"]], TRUE),
                   replaceNULL(args[["test.effect.subgroup"]], TRUE),
                   replaceNULL(args[["col.subgroup"]], "black"),
                   #
                   replaceNULL(args[["forest.I2"]], TRUE),
                   replaceNULL(args[["forest.tau2"]], TRUE),
                   replaceNULL(args[["forest.tau"]], FALSE),
                   replaceNULL(args[["forest.Q"]], TRUE),
                   replaceNULL(args[["forest.pval.Q"]], TRUE),
                   replaceNULL(args[["forest.Rb"]], FALSE),
                   #
                   replaceNULL(args[["digits"]], 2),
                   replaceNULL(args[["digits.I2"]], 0),
                   replaceNULL(args[["digits.tau2"]], 3),
                   replaceNULL(args[["digits.tau"]], 4),
                   #
                   replaceNULL(args[["CIbracket"]], "["),
                   replaceNULL(args[["CIseparator"]], ", "),
                   replaceNULL(args[["zero.pval"]], TRUE),
                   replaceNULL(args[["JAMA.pval"]], FALSE),
                   replaceNULL(args[["text.common"]], "Fixed effect model"),
                   replaceNULL(args[["text.w.common"]], "fixed"),
                   replaceNULL(args[["hetlab"]], "Heterogeneity: "),
                   replaceNULL(args[["header.line"]], TRUE)
                   ),
        setting = "RevMan 5 settings",
        quietly = quietly)
    }
    ##
    else if (setting == "IQWiG5") {
      specificSettings(args = c("method.random.ci", "prediction"),
                       new = list(replaceNULL(args[["method.random.ci"]],
                                              "HK"),
                                  replaceNULL(args[["prediction"]], TRUE)),
                       setting = "IQWiG 5 settings",
                       quietly = quietly)
    }
    ##
    else if (setting == "IQWiG6") {
      specificSettings(args = c("method.random.ci", "adhoc.hakn.ci",
                                "method.tau", "prediction"),
                       new = list(replaceNULL(args[["method.random.ci"]],
                                              "HK"),
                                  replaceNULL(args[["adhoc.hakn.ci"]],
                                              "IQWiG6"),
                                  replaceNULL(args[["method.tau"]], "PM"),
                                  replaceNULL(args[["prediction"]], TRUE)),
                       setting = "IQWiG 6 settings",
                       quietly = quietly)
    }
    ##
    else if (setting == "meta4") {
      specificSettings(args = c("method.tau", "method.I2", "method.predict",
                                "exact.smd",
                                "text.common", "text.w.common",
                                "warn.deprecated"),
                       new = list(replaceNULL(args[["method.tau"]], "DL"),
                                  replaceNULL(args[["method.I2"]], "Q"),
                                  replaceNULL(args[["method.predict"]], "HTS"),
                                  replaceNULL(args[["exact.smd"]], FALSE),
                                  replaceNULL(args[["text.common"]],
                                              "Fixed effect model"),
                                  replaceNULL(args[["text.w.common"]],
                                              "fixed"),
                                  replaceNULL(args[["warn.deprecated"]],
                                              FALSE)),
                       setting =
                         "settings from meta, version 4 or below",
                       quietly = quietly)
    }
    ##
    else if (setting == "meta7") {
      specificSettings(args = c("method.tau", "method.I2", "method.predict",
                                "exact.smd",
                                "text.common", "text.w.common",
                                "warn.deprecated"),
                       new = list(replaceNULL(args[["method.tau"]], "REML"),
                                  replaceNULL(args[["method.I2"]], "Q"),
                                  replaceNULL(args[["method.predict"]], "HTS"),
                                  replaceNULL(args[["exact.smd"]], TRUE),
                                  replaceNULL(args[["text.common"]],
                                              "Common effect model"),
                                  replaceNULL(args[["text.w.common"]],
                                              "common"),
                                  replaceNULL(args[["warn.deprecated"]],
                                              FALSE)),
                       setting =
                         "settings from meta, version 7.0-0 or below",
                       quietly = quietly)
    }
    ##
    else if (setting == "geneexpr") {
      specificSettings(args = c("scientific.pval", "method.tau.ci"),
                       new = list(replaceNULL(args[["scientific.pval"]], TRUE),
                                  replaceNULL(args[["method.tau.ci"]], "")),
                       setting = "Settings for gene expression data",
                       quietly = quietly)
    }
    #
    else if (setting == "IVhet") {
      specificSettings(args = c("method.common.ci",
                                "text.common", "text.w.common"),
                       new =
                         list(replaceNULL(args[["method.common.ci"]], "IVhet"),
                              replaceNULL(args[["text.common"]], "IVhet model"),
                              replaceNULL(args[["text.w.common"]], "IVhet")),
                       setting = "Settings for IVhet model",
                       quietly = quietly)
    }
  }
  
  
  ##
  ## Print settings
  ##
  if (print.settings & !quietly) {
    cat(paste0("\n** Settings for meta-analysis method (R package meta, ",
               "version ", utils::packageDescription("meta")$Version,
               ") **\n\n"))
    ##
    cat(paste0("* General settings *\n"))
    catarg("level              ")
    catarg("level.ma           ")
    catarg("common             ")
    catarg("method.common.ci   ")
    catarg("random             ")
    catarg("method.random.ci   ")
    catarg("adhoc.hakn.ci      ")
    catarg("adhoc.hakn.pi      ")
    catarg("method.tau         ")
    catarg("method.tau.ci      ")
    catarg("level.hetstat      ")
    catarg("tau.common         ")
    catarg("method.I2          ")
    catarg("prediction         ")
    catarg("level.predict      ")
    catarg("method.predict     ")
    catarg("test.subgroup      ")
    catarg("prediction.subgroup")
    catarg("method.bias        ")
    catarg("tool.rob           ")
    catarg("overall.hetstat    ")
    catarg("cid                ")
    catarg("cid.below.null     ")
    catarg("cid.above.null     ")
    catarg("text.common        ")
    catarg("text.random        ")
    catarg("text.predict       ")
    catarg("text.w.common      ")
    catarg("text.w.random      ")
    catarg("title              ")
    catarg("complab            ")
    catarg("CIbracket          ")
    catarg("CIseparator        ")
    catarg("CIlower.blank      ")
    catarg("CIupper.blank      ")
    catarg("print.subgroup.name")
    catarg("sep.subgroup       ")
    catarg("keepdata           ")
    catarg("keeprma            ")
    catarg("warn               ")
    catarg("warn.deprecated    ")
    catarg("backtransf         ")
    catarg("digits             ")
    catarg("digits.mean        ")
    catarg("digits.sd          ")
    catarg("digits.se          ")
    catarg("digits.stat        ")
    catarg("digits.Q           ")
    catarg("digits.tau2        ")
    catarg("digits.tau         ")
    catarg("digits.H           ")
    catarg("digits.I2          ")
    catarg("digits.prop        ")
    catarg("digits.weight      ")
    catarg("digits.pval        ")
    catarg("digits.pval.Q      ")
    catarg("scientific.pval    ")
    catarg("big.mark           ")
    catarg("zero.pval          ")
    catarg("JAMA.pval          ")
    catarg("digits.df          ")
    catarg("digits.cid         ")
    catarg("details            ")
    catarg("print.tau2         ")
    catarg("print.tau2.ci      ")
    catarg("print.tau          ")
    catarg("print.tau.ci       ")
    catarg("print.I2           ")
    catarg("print.I2.ci        ")
    catarg("print.H            ")
    catarg("print.Rb           ")
    catarg("text.tau2          ")
    catarg("text.tau           ")
    catarg("text.I2            ")
    catarg("text.Rb            ")
    catarg("print.Q            ")
    ##
    cat(paste("\n* Default summary measure (argument 'sm' in",
               "corresponding function) *\n"))
    cat("- metabin():  ")
    catarg("smbin ", newline = FALSE)
    cat("- metacont(): ")
    catarg("smcont", newline = FALSE)
    cat("- metacor():  ")
    catarg("smcor ", newline = FALSE)
    cat("- metainc():  ")
    catarg("sminc ", newline = FALSE)
    cat("- metamean(): ")
    catarg("smmean", newline = FALSE)
    cat("- metaprop(): ")
    catarg("smprop", newline = FALSE)
    cat("- metarate(): ")
    catarg("smrate", newline = FALSE)
    ##
    cat(paste("\n* Additional settings for metabin(), metainc(),",
              "metaprop(), and metarate() *\n"))
    catarg("incr       ")
    catarg("method.incr")
    #catarg("allincr")
    #catarg("addincr")
    ##
    cat("\n* Additional settings for metabin() *\n")
    catarg("method     ")
    catarg("allstudies ")
    catarg("MH.exact   ")
    catarg("RR.Cochrane")
    catarg("Q.Cochrane ")
    catarg("model.glmm ")
    catarg("print.CMH  ")
    ##
    cat("\n* Additional settings for metacont() *\n")
    catarg("pooledvar ")
    catarg("method.smd")
    catarg("sd.glass  ")
    catarg("exact.smd ")
    catarg("method.ci.cont")
    ##
    cat("\n* Additional settings for metagen() *\n")
    catarg("transf ")
    ##
    cat("\n* Additional setting for metaprop() *\n")
    catarg("method.ci.prop")
    ##
    cat("\n* Additional setting for metarate() *\n")
    catarg("method.ci.rate")
    ##
    cat("\n* Settings for R functions comparing two treatments *\n")
    catarg("label.e    ")
    catarg("label.c    ")
    catarg("label.left ")
    catarg("label.right")
    ##
    cat("\n* Settings for forest.meta() *\n")
    catarg("layout                 ")
    catarg("forest.details         ")
    catarg("test.overall           ")
    catarg("test.effect.subgroup   ")
    catarg("digits.forest          ",
           end = "\n  (argument 'digits' in forest.meta())")
    catarg("digits.TE.forest          ",
           end = "\n  (argument 'digits.TE' in forest.meta())")
    ##
    catarg("lty.common             ")
    catarg("lty.random             ")
    catarg("col.common             ")
    catarg("col.random             ")
    ##
    catarg("sort.subgroup          ")
    ##
    catarg("pooled.events          ")
    catarg("pooled.times           ")
    catarg("study.results          ")
    ##
    catarg("lty.cid                ")
    catarg("col.cid                ")
    catarg("fill.cid               ")
    catarg("cid.pooled.only        ")
    catarg("fill.equi              ")
    ##
    catarg("fill                   ")
    ##
    catarg("leftcols               ")
    catarg("rightcols              ")
    catarg("leftlabs               ")
    catarg("rightlabs              ")
    ##
    catarg("label.e.attach         ")
    catarg("label.c.attach         ")
    ##
    catarg("bottom.lr              ")
    ##
    catarg("lab.NA                 ")
    catarg("lab.NA.effect          ")
    catarg("lab.NA.weight          ")
    ##
    catarg("lwd                    ")
    catarg("lwd.square             ")
    catarg("lwd.diamond            ")
    ##
    catarg("arrow.type             ")
    catarg("arrow.length           ")
    ##
    catarg("type.study             ")
    catarg("type.common            ")
    ##
    catarg("col.study              ")
    catarg("col.square             ")
    catarg("col.square.lines       ")
    catarg("col.circle             ")
    catarg("col.inside             ")
    catarg("col.diamond            ")
    catarg("col.diamond.lines      ")
    catarg("col.predict            ")
    catarg("col.predict.lines      ")
    catarg("col.subgroup           ")
    catarg("col.label.right        ")
    catarg("col.label.left         ")
    ##
    catarg("col.lines              ")
    catarg("col.label              ")
    ##
    catarg("hetlab                 ")
    catarg("resid.hetstat          ")
    catarg("resid.hetlab           ")
    ##
    catarg("forest.I2              ")
    catarg("forest.I2.ci           ")
    catarg("forest.tau2            ")
    catarg("forest.tau2.ci         ")
    catarg("forest.tau             ")
    catarg("forest.tau.ci          ")
    catarg("forest.Q               ")
    catarg("forest.pval.Q          ")
    catarg("forest.Rb              ")
    catarg("forest.Rb.ci           ")
    ##
    catarg("text.subgroup.nohet    ")
    ##
    catarg("LRT                    ")
    ##
    catarg("forest.stat            ")
    catarg("forest.Q.subgroup      ")
    ##
    catarg("header.line            ")
    ##
    catarg("fontsize               ")
    catarg("fontfamily             ")
    catarg("fs.common              ")
    catarg("fs.random              ")
    catarg("fs.predict             ")
    catarg("fs.common.labels       ")
    catarg("fs.random.labels       ")
    catarg("fs.predict.labels      ")
    catarg("fs.hetstat             ")
    catarg("fs.test.overall        ")
    catarg("fs.test.subgroup       ")
    catarg("fs.test.effect.subgroup")
    catarg("fs.addline             ")
    ##
    catarg("ff.heading             ")
    catarg("ff.common              ")
    catarg("ff.random              ")
    catarg("ff.predict             ")
    catarg("ff.common.labels       ")
    catarg("ff.random.labels       ")
    catarg("ff.predict.labels      ")
    catarg("ff.study               ")
    catarg("ff.hetstat             ")
    catarg("ff.test.overall        ")
    catarg("ff.test.subgroup       ")
    catarg("ff.test.effect.subgroup")
    catarg("ff.addline             ")
    catarg("ff.axis                ")
    catarg("ff.smlab               ")
    catarg("ff.xlab                ")
    catarg("ff.lr                  ")
    ##
    catarg("colgap                 ")
    catarg("colgap.forest          ")
    ##
    catarg("width                  ")
    ##
    catarg("calcwidth.predict      ")
    catarg("calcwidth.hetstat      ")
    catarg("calcwidth.tests        ")
    catarg("calcwidth.subgroup     ")
    catarg("calcwidth.addline      ")
    ##
    catarg("just.studlab           ")
    catarg("just.addcols           ")
    ##
    catarg("spacing                ")
    catarg("addrow                 ")
    catarg("addrow.overall         ")
    catarg("addrow.subgroups       ")
    catarg("addrows.below.overall  ")
  }
  
  
  ##
  ## Set individual settings
  ##
  if (any(!(names %in% c("reset", "print", "setting")))) {
    ##
    ## General settings
    ##
    setlevel("level", args)
    #
    setOptionDepr(args, "level.ma", "level.comb", setlevel)
    #
    setOptionDepr(args, "common", "comb.fixed", setlogical)
    setOptionDepr(args, "common", "fixed", setlogical)
    #
    setOptionDepr(args, "random", "comb.random", setlogical)
    #
    setcharacter("method.common.ci", args, gs("meth4common.ci"))
    #
    na <- is.na(setcharacter("method.random.ci", args, gs("meth4random.ci")))
    depr <- setlogical("hakn", args)
    if (na & !is.na(depr)) {
      if (gs("hakn"))
        setOption("method.random.ci", "HK")
      else
        setOption("method.random.ci", "classic")
    }
    #
    setOptionDepr(args, "adhoc.hakn.ci", "adhoc.hakn", setcharacter,
                  set = gs("adhoc4hakn.ci"))
    #
    setcharacter("adhoc.hakn.pi", args, gs("adhoc4hakn.pi"))
    setcharacter("method.tau", args, gs("meth4tau"))
    setcharacter("method.tau.ci", args, c("J", "BJ", "QP", "PL", ""))
    setlevel("level.hetstat", args)
    setlogical("tau.common", args)
    setcharacter("method.I2", args, gs("meth4i2"))
    setlogical("prediction", args)
    setlevel("level.predict", args)
    setcharacter("method.predict", args, gs("meth4pi"))
    setcharacter("method.bias", args, gs("meth4bias"))
    setcharacter("tool.rob", args, gs("tool4rob"), NULL.ok = TRUE)
    setlogical("overall.hetstat", args, NULL.ok = TRUE)
    setcharacter("text.common", args)
    setcharacter("text.random", args)
    setcharacter("text.predict", args)
    setcharacter("text.w.common", args)
    setcharacter("text.w.random", args)
    setcharacter("title", args)
    setcharacter("complab", args)
    setcharacter("CIbracket", args, c("[", "(", "{", ""))
    setcharacter("CIseparator", args)
    setlogical("CIlower.blank", args)
    setlogical("CIupper.blank", args)
    #
    setOptionDepr(args, "print.subgroup.name", "print.byvar", setlogical)
    #
    setOptionDepr(args, "sep.subgroup", "byseparator", setcharacter)
    #
    setlogical("keepdata", args)
    setlogical("keeprma", args)
    setlogical("warn", args)
    setlogical("warn.deprecated", args)
    setlogical("backtransf", args)
    setnumeric("digits", args)
    setnumeric("digits.mean", args)
    setnumeric("digits.sd", args)
    setnumeric("digits.se", args)
    #
    setOptionDepr(args, "digits.stat", "digits.zval", setnumeric)
    #
    setnumeric("digits.Q", args) 
    setnumeric("digits.tau2", args)
    setnumeric("digits.tau", args)
    setnumeric("digits.H", args) 
    setnumeric("digits.I2", args)
    setnumeric("digits.prop", args)
    setnumeric("digits.weight", args)
    setnumeric("digits.pval", args)
    setnumeric("digits.pval.Q", args)
    setlogical("scientific.pval", args)
    setcharacter("big.mark", args)
    setlogical("zero.pval", args)
    setlogical("JAMA.pval", args)
    setnumeric("digits.df", args)
    setnumeric("digits.cid", args)
    ##
    setlogical("details", args)
    setlogical("print.tau2", args)
    setlogical("print.tau2.ci", args)
    setlogical("print.tau", args)
    setlogical("print.tau.ci", args)
    setlogical("print.I2", args)
    setlogical("print.I2.ci", args)
    setlogical("print.H", args)
    setlogical("print.Rb", args)
    ##
    setcharacter("text.tau2", args)
    setcharacter("text.tau", args)
    setcharacter("text.I2", args)
    setcharacter("text.Rb", args)
    ##
    setlogical("print.Q", args)
    #
    setnumeric("cid", args)
    setOptionDepr(args, "cid.below.null", "lower.equi", setnumeric)
    setOptionDepr(args, "cid.above.null", "upper.equi", setnumeric)
    setOptionDepr(args, "lty.cid", "lty.equi", setnumeric)
    setOptionDepr(args, "col.cid", "col.equi", setcolor)
    setcolor("fill.cid", args)
    setcolor("fill.equi", args)
    setlogical("cid.pooled.only", args)
    ##
    ## R function metabin
    ##
    setcharacter("smbin", args, gs("sm4bin"))
    setcharacter("method", args, gs("meth4bin"))
    setcharacter("model.glmm", args, c("UM.FS", "UM.RS", "CM.EL", "CM.AL"))
    idincr <- argid(names.all, "incr")
    if (!is.na(idincr)) {
      incr <- args[[idincr]]
      if (!is.numeric(incr))
        incr <- setchar(incr, "TACC",
                        "should be numeric or the character string \"TACC\"")
      setOption("incr", incr)
    }
    ##
    na <- is.na(setcharacter("method.incr", args, gs("meth4incr")))
    depr1 <- setlogical("allincr", args)
    depr2 <- setlogical("addincr", args)
    if (na & !is.na(depr1) & gs("allincr"))
      setOption("method.incr", "if0all")
    if (na & !is.na(depr2) & gs("addincr"))
      setOption("method.incr", "all")
    ##
    setlogical("allstudies", args)
    setlogical("MH.exact", args)
    setlogical("RR.Cochrane", args)
    setlogical("Q.Cochrane", args)
    setlogical("print.CMH", args)
    ##
    ## R function metacont
    ##
    setcharacter("smcont", args, gs("sm4cont"))
    setlogical("pooledvar", args)
    setcharacter("method.smd", args, c("Hedges", "Cohen", "Glass"))
    setcharacter("sd.glass", args, c("control", "experimental"))
    setlogical("exact.smd", args)
    setcharacter("method.ci.cont", args, gs("ci4cont"))
    ##
    ## R function metagen
    ##
    setlogical("transf", args)
    ##
    ## R function metaprop
    ##
    setcharacter("smprop", args, gs("sm4prop"))
    setcharacter("method.ci.prop", args, gs("ci4prop"))
    ##
    ## R function metarate
    ##
    setcharacter("smrate", args, gs("sm4rate"))
    setcharacter("method.ci.rate", args, gs("ci4rate"))
    ##
    ## Other meta-analysis functions
    ##
    setcharacter("smcor" , args, gs("sm4cor"))
    setcharacter("sminc" , args, gs("sm4inc"))
    setcharacter("smmean", args, gs("sm4mean"))
    ##
    ## R functions comparing two treatments
    ##
    setcharacter("label.e", args)
    setcharacter("label.c", args)
    setcharacter("label.left", args)
    setcharacter("label.right", args)
    ##
    ## R function forest.meta
    ##
    setcharacter("layout", args, layouts)
    setlogical("forest.details", args)
    setlogical("test.overall", args)
    setlogical("test.subgroup", args)
    setlogical("prediction.subgroup", args)
    setlogical("test.effect.subgroup", args)
    setnumeric("digits.forest", args)
    setnumeric("digits.TE.forest", args)
    ##
    setnumeric("lty.common", args)
    setnumeric("lty.random", args)
    setcolor("col.common", args)
    setcolor("col.random", args)
    ##
    setlogical("sort.subgroup", args)
    ##
    setlogical("pooled.events", args)
    setlogical("pooled.times", args)
    setlogical("study.results", args)
    ##
    setcolor("fill", args)
    ##
    setcharacter("leftcols", args, length = 0,
                 NULL.ok = TRUE, logical.ok = TRUE)
    setcharacter("rightcols", args, length = 0,
                 NULL.ok = TRUE, logical.ok = TRUE)
    setcharacter("leftlabs", args, length = 0, NULL.ok = TRUE)
    setcharacter("rightlabs", args, length = 0, NULL.ok = TRUE)
    ##
    setcharacter("label.e.attach", args, NULL.ok = TRUE)
    setcharacter("label.c.attach", args, NULL.ok = TRUE)
    ##
    setlogical("bottom.lr", args)
    ##
    setcharacter("lab.NA", args)
    setcharacter("lab.NA.effect", args, NULL.ok = TRUE)
    setcharacter("lab.NA.weight", args)
    ##
    setnumeric("lwd", args)
    setnumeric("lwd.square", args)
    setnumeric("lwd.diamond", args)
    ##
    setcharacter("arrow.type", args)
    setnumeric("arrow.length", args)
    ##
    setcharacter("type.study", args, c("square", "diamond", "predict"))
    setcharacter("type.common", args, c("square", "diamond", "predict"))
    ##
    setcolor("col.study", args)
    setcolor("col.square", args)
    setcolor("col.square.lines", args)
    setcolor("col.circle", args)
    setcolor("col.inside", args)
    setcolor("col.diamond", args)
    setcolor("col.diamond.lines", args)
    setcolor("col.predict", args)
    setcolor("col.predict.lines", args)
    setcolor("col.subgroup", args)
    setcolor("col.label.right", args)
    setcolor("col.label.left", args)
    ##
    setcolor("col.lines", args)
    setcolor("col.label", args)
    ##
    setcharacter("hetlab", args)
    setlogical("resid.hetstat", args, TRUE)
    setcharacter("resid.hetlab", args)
    ##
    setlogical("forest.I2", args, TRUE)
    setlogical("forest.I2.ci", args)
    setlogical("forest.tau2", args, TRUE)
    setlogical("forest.tau2.ci", args)
    setlogical("forest.tau", args)
    setlogical("forest.tau.ci", args)
    setlogical("forest.Q", args)
    setlogical("forest.pval.Q", args, TRUE)
    setlogical("forest.Rb", args)
    setlogical("forest.Rb.ci", args)
    ##
    setcharacter("text.subgroup.nohet", args)
    ##
    setlogical("LRT", args)
    ##
    setlogical("forest.stat", args)
    setlogical("forest.Q.subgroup", args)
    ##
    setlogical("header.line", args, ignore.other = TRUE)
    setcharacter("header.line", args, c("both", "below", ""),
                 ignore.other = TRUE)
    ##
    setnumeric("fontsize", args)
    setcharacter("fontfamily", args, NULL.ok = TRUE)
    setnumeric("fs.common", args, TRUE)
    setnumeric("fs.random", args, TRUE)
    setnumeric("fs.predict", args, TRUE)
    setnumeric("fs.common.labels", args, TRUE)
    setnumeric("fs.random.labels", args, TRUE)
    setnumeric("fs.predict.labels", args, TRUE)
    setnumeric("fs.hetstat", args, TRUE)
    setnumeric("fs.test.overall", args, TRUE)
    setnumeric("fs.test.subgroup", args, TRUE)
    setnumeric("fs.test.effect.subgroup", args, TRUE)
    setnumeric("fs.addline", args, TRUE)
    ##
    setcharacter("ff.heading", args)
    setcharacter("ff.common", args, NULL.ok = TRUE)
    setcharacter("ff.random", args, NULL.ok = TRUE)
    setcharacter("ff.predict", args, NULL.ok = TRUE)
    setcharacter("ff.common.labels", args, NULL.ok = TRUE)
    setcharacter("ff.random.labels", args, NULL.ok = TRUE)
    setcharacter("ff.predict.labels", args, NULL.ok = TRUE)
    setcharacter("ff.study", "plain")
    setcharacter("ff.hetstat", args, NULL.ok = TRUE)
    setcharacter("ff.test.overall", args, NULL.ok = TRUE)
    setcharacter("ff.test.subgroup", args, NULL.ok = TRUE)
    setcharacter("ff.test.effect.subgroup", args, NULL.ok = TRUE)
    setcharacter("ff.addline", args, NULL.ok = TRUE)
    setcharacter("ff.axis", args)
    setcharacter("ff.smlab", args)
    setcharacter("ff.xlab", args)
    setcharacter("ff.lr", args)
    ##
    setcharacter("colgap", args)
    setcharacter("colgap.forest", args)
    ##
    setnumeric("width", args, TRUE)
    ##
    setlogical("calcwidth.predict", args)
    setlogical("calcwidth.hetstat", args)
    setlogical("calcwidth.tests", args)
    setlogical("calcwidth.subgroup", args)
    setlogical("calcwidth.addline", args)
    ##
    setcharacter("just.studlab", args, c("right", "center", "left"))
    setcharacter("just.addcols", args, c("right", "center", "left"))
    ##
    setnumeric("spacing", args)
    setlogical("addrow", args, TRUE)
    setlogical("addrow.overall", args, TRUE)
    setlogical("addrow.subgroups", args, TRUE)
    setnumeric("addrows.below.overall", args)
  }
  
  
  if (any(gs("method.random.ci") == "KR") & gs("method.tau") != "REML")
    warning("Default settings for arguments 'method.tau' and ",
            "'method.random.ci' do not match:\n",
            "- Kenward-Roger method (method.random.ci = \"KR\") ",
            "can only be used with\n  ",
            "REML estimator (method.tau = \"REML\").",
            call. = FALSE)
  #
  if (any(gs("method.predict") == "KR") & gs("method.tau") != "REML")
    warning("Default settings for arguments 'method.tau' and ",
            "'method.predict' do not match:\n",
            "- Kenward-Roger method (method.predict = \"KR\") ",
            "can only be used with\n  ",
            "REML estimator (method.tau = \"REML\").",
            call. = FALSE)
  #
  if (any(gs("method.predict") == "KR-PR") & gs("method.tau") != "REML")
    warning("Default settings for arguments 'method.tau' and ",
            "'method.predict' do not match:\n",
            "- Kenward-Roger method (method.predict = \"KR-PR\") ",
            "can only be used with\n  ",
            "REML estimator (method.tau = \"REML\").",
            call. = FALSE)
  
  
  invisible(oldset)
}
