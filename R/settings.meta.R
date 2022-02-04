#' Print and change default settings to conduct and print or plot
#' meta-analyses in R package \bold{meta}.
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
#' \item \code{settings.meta("revman5")}
#' \item \code{settings.meta("jama")}
#' \item \code{settings.meta("iqwig5")}
#' \item \code{settings.meta("iqwig6")}
#' \item \code{settings.meta("geneexpr")}
#' \item \code{settings.meta("meta4")}
#' }
#'
#' The first command can be used to reproduce meta-analyses from
#' Cochrane reviews conducted with \emph{Review Manager 5} (RevMan 5,
#' \url{https://training.cochrane.org/online-learning/core-software-cochrane-reviews/revman})
#' and specifies to use a RevMan 5 layout in forest plots.
#'
#' The second command can be used to generate forest plots following
#' instructions for authors of the \emph{Journal of the American
#' Medical Association}
#' (\url{https://jamanetwork.com/journals/jama/pages/instructions-for-authors/}).Study
#' labels according to JAMA guidelines can be generated using
#' \code{\link{labels.meta}}.
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
#' The last setting uses the default settings of R package
#' \bold{meta}, version 4 or below.
#' 
#' RevMan 5 settings, in detail:
#' \tabular{lll}{
#' \bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
#' \code{hakn} \tab FALSE \tab method not available in RevMan 5 \cr
#' \code{method.tau} \tab "DL" \tab only available method in RevMan 5
#'   \cr
#' \code{tau.common} \tab FALSE \tab common between-study variance in
#'   subgroups \cr
#' \code{MH.exact} \tab FALSE \tab exact Mantel-Haenszel method \cr
#' \code{RR.Cochrane} \tab TRUE \tab calculation of risk ratios \cr
#' \code{Q.Cochrane} \tab TRUE \tab calculation of heterogeneity statistic \cr
#' \code{layout} \tab "RevMan5" \tab layout for forest plots \cr
#' \code{test.overall} \tab TRUE \tab print information on test of
#'   overall effect \cr
#' \code{digits.I2} \tab 0 \tab number of digits for I-squared measure
#'   \cr
#' \code{digits.tau2} \tab 2 \tab number of digits for tau-squared \cr
#' \code{digits.tau} \tab 4 \tab number of digits for square root of
#'   tau-squared \cr
#' \code{CIbracket}, \tab "[" \tab \cr
#' \code{CIseparator} \tab ", " \tab print confidence intervals as
#'   "\code{[., .]}"
#' }
#'
#' JAMA settings:
#' \tabular{lll}{
#' \bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
#' \code{layout} \tab "JAMA" \tab layout for forest plots \cr
#' \code{test.overall} \tab TRUE \tab print information on test of
#'   overall effect \cr
#' \code{digits.I2} \tab 0 \tab number of digits for I-squared measure
#'   \cr
#' \code{CIbracket}, \tab "(" \tab \cr
#' \code{CIseparator} \tab "-" \tab print confidence intervals as
#'   "\code{(.-.)}" \cr
#' \code{zero.pval}, \tab TRUE \tab print p-values with leading zero
#' \cr
#' \code{JAMA.pval}, \tab TRUE \tab round p-values to three digits
#'   (for 0.001 < p \eqn{\le} 0.01) or two digits (p > 0.01)
#' }
#' 
#' IQWiG, General Methods 5 settings:
#' \tabular{lll}{
#' \bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
#' \code{hakn} \tab TRUE \tab Hartung-Knapp method \cr
#' \code{prediction} \tab TRUE \tab Prediction interval \cr
#' }
#' 
#' IQWiG, General Methods 6 settings:
#' \tabular{lll}{
#' \bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
#' \code{hakn} \tab TRUE \tab Hartung-Knapp method \cr
#' \code{adhoc.hakn} \tab "ci" \tab \emph{ad hoc} variance correction \cr
#' \code{method.tau} \tab "PM" \tab Paule-Mandel estimator for
#'   between-study variance \cr
#' \code{prediction} \tab TRUE \tab Prediction interval \cr
#' }
#' 
#' Settings for gene expression data:
#' \tabular{lll}{
#' \bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
#' \code{scientific.pval} \tab TRUE \tab Scientific notation for p-values \cr
#' \code{method.tau.ci} \tab FALSE \tab no confidence interval for \cr
#'  \tab between-study heterogeneity variance \cr
#' }
#' 
#' Settings for \bold{meta}, version 4 or below:
#' \tabular{lll}{
#' \bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
#' \code{method.tau} \tab "DL" \tab DerSimonian-Laird estimator \cr
#' }
#' 
#' A list of all arguments with current settings is printed using the
#' command \code{settings.meta("print")}.
#' 
#' In order to reset all settings of R package \bold{meta} the command
#' \code{settings.meta("reset")} or \code{settings.meta(reset = TRUE)}
#' can be used.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{gs}}, \code{\link{forest.meta}},
#'   \code{\link{print.meta}}, \code{\link{labels.meta}}
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
#' settings.meta("revman5")
#' forest(metagen(1:3, 2:4 / 10, sm = "MD", fixed = FALSE),
#'   label.left = "Favours A", label.right = "Favours B",
#'   colgap.studlab = "2cm", colgap.forest.left = "0.2cm")
#' 
#' # Forest plot using JAMA style
#' #
#' settings.meta("jama")
#' forest(metagen(1:3, 2:4 / 10, sm = "MD", fixed = FALSE),
#'   label.left = "Favours A", label.right = "Favours B",
#'   colgap.studlab = "2cm", colgap.forest.left = "0.2cm")
#'
#' # Use slightly different layout for confidence intervals
#' # (especially useful if upper confidence limit can be negative)
#' #
#' settings.meta(CIseparator = " - ")
#' forest(metagen(-(1:3), 2:4 / 10, sm="MD", fixed=FALSE),
#'   label.left="Favours A", label.right="Favours B",
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
  settings <- c("RevMan5", "JAMA", "IQWiG5", "IQWiG6", "geneexpr", "meta4")
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
    chkdeprecated(names.all, "fixed", "comb.fixed")
    chkdeprecated(names.all, "random", "comb.random")
    chkdeprecated(names.all, "digits.stat", "digits.zval")
    ##
    chkdeprecated(names.all, "print.subgroup.name", "print.byvar")
    chkdeprecated(names.all, "sep.subgroup", "byseparator")
  }
  ##  
  names <- names.all[!(names.all %in% .settings$argslist.internal)]
  
  
  ##
  ## Check argument names
  ##
  for (i in seq_along(names))
    names[i] <- setchar(names[i],
                        c(.settings$argslist, "reset", "print", "setting", ""),
                        "unmatched",
                        name = names[i])
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
    setOption("level", 0.95)
    setOption("level.ma", 0.95)
    setOption("fixed", TRUE)
    setOption("random", TRUE)
    setOption("hakn", FALSE)
    setOption("adhoc.hakn", "")
    setOption("method.tau", "REML")
    setOption("method.tau.ci", NULL)
    setOption("tau.common", FALSE)
    setOption("prediction", FALSE)
    setOption("level.predict", 0.95)
    setOption("test.subgroup", TRUE)
    setOption("prediction.subgroup", FALSE)
    setOption("method.bias", "Egger")
    setOption("text.fixed", "Common effect model")
    setOption("text.random", "Random effects model")
    setOption("text.predict", "Prediction interval")
    setOption("text.w.fixed", "common")
    setOption("text.w.random", "random")
    setOption("title", "")
    setOption("complab", "")
    setOption("CIbracket", "[")
    setOption("CIseparator", "; ")
    setOption("CIlower.blank", TRUE)
    setOption("CIupper.blank", TRUE)
    setOption("print.subgroup.name", TRUE)
    setOption("sep.subgroup", " = ")
    setOption("keepdata", TRUE)
    setOption("warn", TRUE)
    setOption("warn.deprecated", TRUE)
    setOption("backtransf", TRUE)
    setOption("digits", 4)
    setOption("digits.se", 4)
    setOption("digits.stat", 2)
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
    setOption("print.I2", TRUE)
    setOption("print.H", TRUE)
    setOption("print.Rb", FALSE)
    setOption("text.tau2", "tau^2")
    setOption("text.tau", "tau")
    setOption("text.I2", "I^2")
    setOption("text.Rb", "Rb")
    ##
    setOption("method", "MH")
    setOption("incr", 0.5)
    setOption("allincr", FALSE)
    setOption("addincr", FALSE)
    setOption("allstudies", FALSE)
    setOption("MH.exact", FALSE)
    setOption("RR.Cochrane", FALSE)
    setOption("Q.Cochrane", TRUE)
    setOption("model.glmm", "UM.FS")
    setOption("print.CMH", FALSE)
    ##
    setOption("smbin", "RR")
    setOption("smcont", "MD")
    setOption("smcor", "ZCOR")
    setOption("sminc", "IRR")
    setOption("smmean", "MRAW")
    setOption("smprop", "PLOGIT")
    setOption("smrate", "IRLN")
    ##
    setOption("pooledvar", FALSE)
    setOption("method.smd", "Hedges")
    setOption("sd.glass", "control")
    setOption("exact.smd", TRUE)
    setOption("method.ci.cont", "z")
    ##
    setOption("method.ci.prop", "CP")
    ##
    setOption("label.e", "Experimental")
    setOption("label.c", "Control")
    setOption("label.left", "")
    setOption("label.right", "")
    ##
    ## Forest plots
    ##
    setOption("layout", "meta")
    setOption("test.overall", FALSE)
    setOption("test.effect.subgroup", FALSE)
    setOption("digits.forest", 2)
  }
  
  
  ##
  ## Specific settings
  ##
  if (specific.settings) {
    ##
    ## Remember:
    ## settings <- c("RevMan5", "JAMA", "IQWiG5", "IQWiG6",
    ##               "geneexpr", "meta4")
    ##
    if (setting == "RevMan5") {
      specificSettings(args = c("hakn", "method.tau", "tau.common",
                                "MH.exact", "RR.Cochrane", "Q.Cochrane",
                                "layout", "test.overall",
                                "test.subgroup", "test.effect.subgroup",
                                "digits.I2", "digits.tau2", "digits.tau",
                                "CIbracket", "CIseparator"),
                       new = list(replaceNULL(args[["hakn"]], FALSE),
                                  replaceNULL(args[["method.tau"]], "DL"),
                                  replaceNULL(args[["tau.common"]], FALSE),
                                  replaceNULL(args[["MH.exact"]], FALSE),
                                  replaceNULL(args[["RR.Cochrane"]], TRUE),
                                  replaceNULL(args[["Q.Cochrane"]], TRUE),
                                  "RevMan5",
                                  replaceNULL(args[["test.overall"]], TRUE),
                                  replaceNULL(args[["test.subgroup"]], TRUE),
                                  replaceNULL(args[["test.effect.subgroup"]], TRUE),
                                  replaceNULL(args[["digits.I2"]], 0),
                                  replaceNULL(args[["digits.tau2"]], 3),
                                  replaceNULL(args[["digits.tau"]], 4),
                                  replaceNULL(args[["CIbracket"]], "["),
                                  replaceNULL(args[["CIseparator"]], ", ")
                                  ),
                       setting = "RevMan 5 settings",
                       quietly = quietly)
    }
    ##
    else if (setting == "JAMA") {
      specificSettings(args = c("layout", "test.overall",
                                "test.subgroup", "test.effect.subgroup",
                                "digits.I2", "digits.pval",
                                "CIbracket", "CIseparator",
                                "zero.pval", "JAMA.pval"),
                       new = list("JAMA", TRUE,
                                  FALSE, FALSE,
                                  0, 3,
                                  "(", "-",
                                  FALSE, TRUE),
                       setting = "JAMA settings",
                       quietly = quietly)
    }
    ##
    else if (setting == "IQWiG5") {
      specificSettings(args = c("hakn", "prediction"),
                       new = list(TRUE, TRUE),
                       setting = "IQWiG 5 settings",
                       quietly = quietly)
    }
    ##
    else if (setting == "IQWiG6") {
      specificSettings(args = c("hakn", "adhoc.hakn",
                                "method.tau", "prediction"),
                       new = list(TRUE, "iqwig6", "PM", TRUE),
                       setting = "IQWiG 6 settings",
                       quietly = quietly)
    }
    ##
    else if (setting == "meta4") {
      specificSettings(args = c("method.tau", "exact.smd",
                                "text.fixed", "text.w.fixed",
                                "warn.deprecated"),
                       new = list("DL", FALSE,
                                  "Fixed effect model", "fixed",
                                  FALSE),
                       setting =
                         "settings from meta, version 4 or below",
                       quietly = quietly)
    }
    ##
    else if (setting == "geneexpr") {
      specificSettings(args = c("scientific.pval", "method.tau.ci"),
                       new = list(TRUE, ""),
                       setting = "Settings for gene expression data",
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
    catarg("fixed              ")
    catarg("random             ")
    catarg("hakn               ")
    catarg("adhoc.hakn         ")
    catarg("method.tau         ")
    catarg("method.tau.ci      ")
    catarg("tau.common         ")
    catarg("prediction         ")
    catarg("level.predict      ")
    catarg("test.subgroup      ")
    catarg("prediction.subgroup")
    catarg("method.bias        ")
    catarg("text.fixed         ")
    catarg("text.random        ")
    catarg("text.predict       ")
    catarg("text.w.fixed       ")
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
    catarg("warn               ")
    catarg("warn.deprecated    ")
    catarg("backtransf         ")
    catarg("digits             ")
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
    catarg("print.I2           ")
    catarg("print.H            ")
    catarg("print.Rb           ")
    catarg("text.tau2          ")
    catarg("text.tau           ")
    catarg("text.I2            ")
    catarg("text.Rb            ")
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
    catarg("incr   ")
    catarg("allincr")
    catarg("addincr")
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
    cat("\n* Additional setting for metaprop() *\n")
    catarg("method.ci.prop")
    ##
    cat("\n* Settings for R functions comparing two treatments *\n")
    catarg("label.e    ")
    catarg("label.c    ")
    catarg("label.left ")
    catarg("label.right")
    ##
    cat("\n* Settings for forest.meta() *\n")
    catarg("layout              ")
    catarg("test.overall        ")
    catarg("test.effect.subgroup")
    catarg("digits.forest       ",
           end = "\n  (argument 'digits' in forest.meta())")
  }
  
  
  ##
  ## Set individual settings
  ##
  if (any(!(names %in% c("reset", "print", "setting")))) {
    idlevel <- argid(names.all, "level")
    idlevel.ma <- argid(names.all, "level.ma")
    idfixed <- argid(names.all, "fixed")
    idrandom <- argid(names.all, "random")
    idhakn <- argid(names.all, "hakn")
    idadhoc.hakn <- argid(names.all, "adhoc.hakn")
    idmethod.tau <- argid(names.all, "method.tau")
    idmethod.tau.ci <- argid(names.all, "method.tau.ci")
    idtau.common <- argid(names.all, "tau.common")
    idprediction <- argid(names.all, "prediction")
    idlevel.predict <- argid(names.all, "level.predict")
    idmethod.bias <- argid(names.all, "method.bias")
    idtext.fixed <- argid(names.all, "text.fixed")
    idtext.random <- argid(names.all, "text.random")
    idtext.predict <- argid(names.all, "text.predict")
    idtext.w.fixed <- argid(names.all, "text.w.fixed")
    idtext.w.random <- argid(names.all, "text.w.random")
    idtitle <- argid(names.all, "title")
    idcomplab <- argid(names.all, "complab")
    idCIbracket <- argid(names.all, "CIbracket")
    idCIseparator <- argid(names.all, "CIseparator")
    idCIlower.blank <- argid(names.all, "CIlower.blank")
    idCIupper.blank <- argid(names.all, "CIupper.blank")
    idprint.subgroup.name <- argid(names.all, "print.subgroup.name")
    idsep.subgroup <- argid(names.all, "sep.subgroup")
    idkeepdata <- argid(names.all, "keepdata")
    idwarn <- argid(names.all, "warn")
    idwarn.deprecated <- argid(names.all, "warn.deprecated")
    idbacktransf <- argid(names.all, "backtransf")
    iddigits <- argid(names.all, "digits")
    iddigits.se <- argid(names.all, "digits.se")
    iddigits.stat <- argid(names.all, "digits.stat")
    iddigits.Q <- argid(names.all, "digits.Q") 
    iddigits.tau2 <- argid(names.all, "digits.tau2")
    iddigits.tau <- argid(names.all, "digits.tau")
    iddigits.H <- argid(names.all, "digits.H") 
    iddigits.I2 <- argid(names.all, "digits.I2")
    iddigits.prop <- argid(names.all, "digits.prop")
    iddigits.weight <- argid(names.all,"digits.weight")
    iddigits.pval <- argid(names.all, "digits.pval")
    iddigits.pval.Q <- argid(names.all, "digits.pval.Q")
    idscientific.pval <- argid(names.all, "scientific.pval")
    idbig.mark <- argid(names.all, "big.mark")
    idzero.pval <- argid(names.all, "zero.pval")
    idJAMA.pval <- argid(names.all, "JAMA.pval")
    idprint.I2 <- argid(names.all, "print.I2")
    idprint.H <- argid(names.all, "print.H")
    idprint.Rb <- argid(names.all, "print.Rb")
    idtext.tau2 <- argid(names.all, "text.tau2")
    idtext.tau <- argid(names.all, "text.tau")
    idtext.I2 <- argid(names.all, "text.I2")
    idtext.Rb <- argid(names.all, "text.Rb")
    ##
    idsmbin <- argid(names.all, "smbin")
    idmethod <- argid(names.all, "method")
    idincr <- argid(names.all, "incr")
    idallincr <- argid(names.all, "allincr")
    idaddincr <- argid(names.all, "addincr")
    idallstudies <- argid(names.all, "allstudies")
    idMH.exact <- argid(names.all, "MH.exact")
    idRR.Cochrane <- argid(names.all, "RR.Cochrane")
    idQ.Cochrane <- argid(names.all, "Q.Cochrane")
    idmodel.glmm <- argid(names.all, "model.glmm")
    idprint.CMH <- argid(names.all, "print.CMH")
    ##
    idsmcont <- argid(names.all, "smcont")
    idsmcor <- argid(names.all, "smcor")
    idsminc <- argid(names.all, "sminc")
    idsmmean <- argid(names.all, "smmean")
    idsmprop <- argid(names.all, "smprop")
    idsmrate <- argid(names.all, "smrate")
    ##
    idpooledvar <- argid(names.all, "pooledvar")
    idmethod.smd <- argid(names.all, "method.smd")
    idsd.glass <- argid(names.all, "sd.glass")
    idexact.smd <- argid(names.all, "exact.smd")
    idmethod.ci.cont <- argid(names.all, "method.ci.cont")
    ##
    idmethod.ci.prop <- argid(names.all, "method.ci.prop")
    ##
    idlabel.e <- argid(names.all, "label.e")
    idlabel.c <- argid(names.all, "label.c")
    idlabel.left <- argid(names.all, "label.left")
    idlabel.right <- argid(names.all, "label.right")
    ##
    idlayout <- argid(names.all, "layout")
    idtest.overall <- argid(names.all, "test.overall")
    idtest.subgroup <- argid(names.all, "test.subgroup")
    idprediction.subgroup <- argid(names.all, "prediction.subgroup")
    idtest.effect.subgroup <- argid(names.all, "test.effect.subgroup")
    iddigits.forest <- argid(names.all, "digits.forest")
    ##
    ## General settings
    ##
    if (!is.na(idlevel)) {
      level <- args[[idlevel]]
      chklevel(level)
      setOption("level", level)
    }
    if (!is.na(idlevel.ma)) {
      level.ma <- args[[idlevel.ma]]
      chklevel(level.ma)
      setOption("level.ma", level.ma)
    }
    if (!is.na(idfixed)) {
      fixed <- args[[idfixed]]
      chklogical(fixed)
      setOption("fixed", fixed)
    }
    if (!is.na(idrandom)) {
      random <- args[[idrandom]]
      chklogical(random)
      setOption("random", random)
    }
    if (!is.na(idhakn)) {
      hakn <- args[[idhakn]]
      chklogical(hakn)
      setOption("hakn", hakn)
    }
    if (!is.na(idadhoc.hakn)) {
      adhoc.hakn <- args[[idadhoc.hakn]]
      adhoc.hakn <- setchar(adhoc.hakn, gs("adhoc4hakn"))
      setOption("adhoc.hakn", adhoc.hakn)
    }
    if (!is.na(idmethod.tau)) {
      method.tau <- args[[idmethod.tau]]
      method.tau <- setchar(method.tau, gs("meth4tau"))
      setOption("method.tau", method.tau)
    }
    if (!is.na(idmethod.tau.ci)) {
      method.tau.ci <- args[[idmethod.tau.ci]]
      method.tau.ci <- setchar(method.tau.ci, c("J", "BJ", "QP", "PL", ""))
      setOption("method.tau.ci", method.tau.ci)
    }
    if (!is.na(idtau.common)) {
      tau.common <- args[[idtau.common]]
      chklogical(tau.common)
      setOption("tau.common", tau.common)
    }
    if (!is.na(idprediction)) {
      prediction <- args[[idprediction]]
      chklogical(prediction)
      setOption("prediction", prediction)
    }
    if (!is.na(idlevel.predict)) {
      level.predict <- args[[idlevel.predict]]
      chklevel(level.predict)
      setOption("level.predict", level.predict)
    }
    if (!is.na(idmethod.bias)) {
      method.bias <- args[[idmethod.bias]]
      method.bias <- setchar(method.bias, gs("meth4bias"))
      setOption("method.bias", method.bias)
    }
    if (!is.na(idtext.fixed)) {
      text.fixed <- args[[idtext.fixed]]
      if (length(text.fixed) != 1)
        stop("Argument 'text.fixed' must be a character string.")
      ##
      setOption("text.fixed", text.fixed)
    }
    if (!is.na(idtext.random)) {
      text.random <- args[[idtext.random]]
      if (length(text.random) != 1)
        stop("Argument 'text.random' must be a character string.")
      ##
      setOption("text.random", text.random)
    }
    if (!is.na(idtext.predict)) {
      text.predict <- args[[idtext.predict]]
      if (length(text.predict) != 1)
        stop("Argument 'text.predict' must be a character string.")
      ##
      setOption("text.predict", text.predict)
    }
    if (!is.na(idtext.w.fixed)) {
      text.w.fixed <- args[[idtext.w.fixed]]
      if (length(text.w.fixed) != 1)
        stop("Argument 'text.w.fixed' must be a character string.")
      ##
      setOption("text.w.fixed", text.w.fixed)
    }
    if (!is.na(idtext.w.random)) {
      text.w.random <- args[[idtext.w.random]]
      if (length(text.w.random) != 1)
        stop("Argument 'text.w.random' must be a character string.")
      ##
      setOption("text.w.random", text.w.random)
    }
    if (!is.na(idtitle)) {
      title <- args[[idtitle]]
      if (length(title) != 1)
        stop("Argument 'title' must be a character string.")
      ##
      setOption("title", title)
    }
    if (!is.na(idcomplab)) {
      complab <- args[[idcomplab]]
      if (length(complab) != 1)
        stop("Argument 'complab' must be a character string.")
      ##
      setOption("complab", complab)
    }
    if (!is.na(idCIbracket)) {
      CIbracket <- args[[idCIbracket]]
      CIbracket <- setchar(CIbracket,
                           c("[", "(", "{", ""))
      setOption("CIbracket", CIbracket)
    }
    if (!is.na(idCIseparator)) {
      CIseparator <- args[[idCIseparator]]
      if (length(CIseparator) != 1)
        stop("Argument 'CIseparator' must be a character string.")
      ##
      setOption("CIseparator", CIseparator)
    }
    if (!is.na(idCIlower.blank)) {
      CIlower.blank <- args[[idCIlower.blank]]
      chklogical(CIlower.blank)
      setOption("CIlower.blank", CIlower.blank)
    }
    if (!is.na(idCIupper.blank)) {
      CIupper.blank <- args[[idCIupper.blank]]
      chklogical(CIupper.blank)
      setOption("CIupper.blank", CIupper.blank)
    }
    if (!is.na(idprint.subgroup.name)) {
      print.subgroup.name <- args[[idprint.subgroup.name]]
      chklogical(print.subgroup.name)
      setOption("print.subgroup.name", print.subgroup.name)
    }
    if (!is.na(idsep.subgroup)) {
      sep.subgroup <- args[[idsep.subgroup]]
      if (length(sep.subgroup) != 1)
        stop("Argument 'sep.subgroup' must be a character string.")
      ##
      setOption("sep.subgroup", sep.subgroup)
    }
    if (!is.na(idkeepdata)) {
      keepdata <- args[[idkeepdata]]
      chklogical(keepdata)
      setOption("keepdata", keepdata)
    }
    if (!is.na(idwarn)) {
      warn <- args[[idwarn]]
      chklogical(warn)
      setOption("warn", warn)
    }
    if (!is.na(idwarn.deprecated)) {
      warn.deprecated <- args[[idwarn.deprecated]]
      chklogical(warn.deprecated)
      setOption("warn.deprecated", warn.deprecated)
    }
    if (!is.na(idbacktransf)) {
      backtransf <- args[[idbacktransf]]
      chklogical(backtransf)
      setOption("backtransf", backtransf)
    }
    if (!is.na(iddigits)) {
      digits <- args[[iddigits]]
      chknumeric(digits, min = 0, length = 1)
      setOption("digits", digits)
    }
    if (!is.na(iddigits.se)) {
      digits.se <- args[[iddigits.se]]
      chknumeric(digits.se, min = 0, length = 1)
      setOption("digits.se", digits.se)
    }
    if (!is.na(iddigits.stat)) {
      digits.stat <- args[[iddigits.stat]]
      chknumeric(digits.stat, min = 0, length = 1)
      setOption("digits.stat", digits.stat)
    }
    if (!is.na(iddigits.Q)) {
      digits.Q <- args[[iddigits.Q]]
      chknumeric(digits.Q, min = 0, length = 1)
      setOption("digits.Q", digits.Q)
    }
    if (!is.na(iddigits.tau2)) {
      digits.tau2 <- args[[iddigits.tau2]]
      chknumeric(digits.tau2, min = 0, length = 1)
      setOption("digits.tau2", digits.tau2)
    }
    if (!is.na(iddigits.tau)) {
      digits.tau <- args[[iddigits.tau]]
      chknumeric(digits.tau, min = 0, length = 1)
      setOption("digits.tau", digits.tau)
    }
    if (!is.na(iddigits.H)) {
      digits.H <- args[[iddigits.H]]
      chknumeric(digits.H, min = 0, length = 1)
      setOption("digits.H", digits.H)
    }
    if (!is.na(iddigits.I2)) {
      digits.I2 <- args[[iddigits.I2]]
      chknumeric(digits.I2, min = 0, length = 1)
      setOption("digits.I2", digits.I2)
    }
    if (!is.na(iddigits.prop)) {
      digits.prop <- args[[iddigits.prop]]
      chknumeric(digits.prop, min = 0, length = 1)
      setOption("digits.prop", digits.prop)
    }
    if (!is.na(iddigits.weight)) {
      digits.weight <- args[[iddigits.weight]]
      chknumeric(digits.weight, min = 0, length = 1)
      setOption("digits.weight", digits.weight)
    }
    if (!is.na(iddigits.pval)) {
      digits.pval <- args[[iddigits.pval]]
      chknumeric(digits.pval, min = 0, length = 1)
      setOption("digits.pval", digits.pval)
    }
    if (!is.na(iddigits.pval.Q)) {
      digits.pval.Q <- args[[iddigits.pval.Q]]
      chknumeric(digits.pval.Q, min = 0, length = 1)
      setOption("digits.pval.Q", digits.pval.Q)
    }
    if (!is.na(idscientific.pval)) {
      scientific.pval <- args[[idscientific.pval]]
      chklogical(scientific.pval)
      setOption("scientific.pval", scientific.pval)
    }
    if (!is.na(idbig.mark)) {
      big.mark <- args[[idbig.mark]]
      if (length(big.mark) != 1)
        stop("Argument 'big.mark' must be a character string.")
      ##
      setOption("big.mark", big.mark)
    }
    if (!is.na(idzero.pval)) {
      zero.pval <- args[[idzero.pval]]
      chklogical(zero.pval)
      setOption("zero.pval", zero.pval)
    }
    if (!is.na(idJAMA.pval)) {
      JAMA.pval <- args[[idJAMA.pval]]
      chklogical(JAMA.pval)
      setOption("JAMA.pval", JAMA.pval)
    }
    if (!is.na(idprint.I2)) {
      print.I2 <- args[[idprint.I2]]
      chklogical(print.I2)
      setOption("print.I2", print.I2)
    }
    if (!is.na(idprint.H)) {
      print.H <- args[[idprint.H]]
      chklogical(print.H)
      setOption("print.H", print.H)
    }
    if (!is.na(idprint.Rb)) {
      print.Rb <- args[[idprint.Rb]]
      chklogical(print.Rb)
      setOption("print.Rb", print.Rb)
    }
    if (!is.na(idtext.tau2)) {
      text.tau2 <- args[[idtext.tau2]]
      if (length(text.tau2) != 1)
        stop("Argument 'text.tau2' must be a character string.")
      ##
      setOption("text.tau2", text.tau2)
    }
    if (!is.na(idtext.tau)) {
      text.tau <- args[[idtext.tau]]
      if (length(text.tau) != 1)
        stop("Argument 'text.tau' must be a character string.")
      ##
      setOption("text.tau", text.tau)
    }
    if (!is.na(idtext.I2)) {
      text.I2 <- args[[idtext.I2]]
      if (length(text.I2) != 1)
        stop("Argument 'text.I2' must be a character string.")
      ##
      setOption("text.I2", text.I2)
    }
    if (!is.na(idtext.Rb)) {
      text.Rb <- args[[idtext.Rb]]
      if (length(text.Rb) != 1)
        stop("Argument 'text.Rb' must be a character string.")
      ##
      setOption("text.Rb", text.Rb)
    }
    ##
    ## R function metabin
    ##
    if (!is.na(idsmbin)) {
      smbin <- args[[idsmbin]]
      smbin <- setchar(smbin, gs("sm4bin"))
      setOption("smbin", smbin)
    }
    if (!is.na(idmethod)) {
      method <- args[[idmethod]]
      method <- setchar(method, gs("meth4bin"))
      setOption("method", method)
    }
    if (!is.na(idmodel.glmm)) {
      model.glmm <- args[[idmodel.glmm]]
      model.glmm <- setchar(model.glmm, c("UM.FS", "UM.RS", "CM.EL", "CM.AL"))
      setOption("model.glmm", model.glmm)
    }
    ##
    if (!is.na(idincr)) {
      incr <- args[[idincr]]
      if (!is.numeric(incr))
        incr <- setchar(incr, "TACC",
                        "should be numeric or the character string \"TACC\"")
      setOption("incr", incr)
    }
    if (!is.na(idallincr)) {
      allincr <- args[[idallincr]]
      chklogical(allincr)
      setOption("allincr", allincr)
    }
    if (!is.na(idaddincr)) {
      addincr <- args[[idaddincr]]
      chklogical(addincr)
      setOption("addincr", addincr)
    }
    if (!is.na(idallstudies)) {
      allstudies <- args[[idallstudies]]
      chklogical(allstudies)
      setOption("allstudies", allstudies)
    }
    if (!is.na(idMH.exact)) {
      MH.exact <- args[[idMH.exact]]
      chklogical(MH.exact)
      setOption("MH.exact", MH.exact)
    }
    if (!is.na(idRR.Cochrane)) {
      RR.Cochrane <- args[[idRR.Cochrane]]
      chklogical(RR.Cochrane)
      setOption("RR.Cochrane", RR.Cochrane)
    }
    if (!is.na(idQ.Cochrane)) {
      Q.Cochrane <- args[[idQ.Cochrane]]
      chklogical(Q.Cochrane)
      setOption("Q.Cochrane", Q.Cochrane)
    }
    if (!is.na(idprint.CMH)) {
      print.CMH <- args[[idprint.CMH]]
      chklogical(print.CMH)
      setOption("print.CMH", print.CMH)
    }
    ##
    ## R function metacont
    ##
    if (!is.na(idsmcont)) {
      smcont <- args[[idsmcont]]
      smcont <- setchar(smcont, gs("sm4cont"))
      setOption("smcont", smcont)
    }
    ##
    ## R function metacor
    ##
    if (!is.na(idsmcor)) {
      smcor <- args[[idsmcor]]
      smcor <- setchar(smcor, gs("sm4cor"))
      setOption("smcor", smcor)
    }
    ##
    ## R function metacont
    ##
    if (!is.na(idpooledvar)) {
      pooledvar <- args[[idpooledvar]]
      chklogical(pooledvar)
      setOption("pooledvar", pooledvar)
    }
    if (!is.na(idmethod.smd)) {
      method.smd <- args[[idmethod.smd]]
      method.smd <- setchar(method.smd, c("Hedges", "Cohen", "Glass"))
      setOption("method.smd", method.smd)
    }
    if (!is.na(idsd.glass)) {
      sd.glass <- args[[idsd.glass]]
      sd.glass <- setchar(sd.glass, c("control", "experimental"))
      setOption("sd.glass", sd.glass)
    }
    if (!is.na(idexact.smd)) {
      exact.smd <- args[[idexact.smd]]
      chklogical(exact.smd)
      setOption("exact.smd", exact.smd)
    }
    ##
    if (!is.na(idmethod.ci.cont)) {
      method.ci.cont <- args[[idmethod.ci.cont]]
      method.ci.cont <- setchar(method.ci.cont, gs("ci4cont"))
      setOption("method.ci.cont", method.ci.cont)
    }
    ##
    ## R function metainc
    ##
    if (!is.na(idsminc)) {
      sminc <- args[[idsminc]]
      sminc <- setchar(sminc, gs("sm4inc"))
      setOption("sminc", sminc)
    }
    ##
    ## R function metamean
    ##
    if (!is.na(idsmmean)) {
      smmean <- args[[idsmmean]]
      smmean <- setchar(smmean, gs("sm4mean"))
      setOption("smmean", smmean)
    }
    ##
    ## R function metaprop
    ##
    if (!is.na(idsmprop)) {
      smprop <- args[[idsmprop]]
      smprop <- setchar(smprop, gs("sm4prop"))
      setOption("smprop", smprop)
    }
    ##
    if (!is.na(idmethod.ci.prop)) {
      method.ci.prop <- args[[idmethod.ci.prop]]
      method.ci.prop <- setchar(method.ci.prop, gs("ci4prop"))
      setOption("method.ci.prop", method.ci.prop)
    }
    ##
    ## R function metarate
    ##
    if (!is.na(idsmrate)) {
      smrate <- args[[idsmrate]]
      smrate <- setchar(smrate, gs("sm4rate"))
      setOption("smrate", smrate)
    }
    ##
    ## R functions comparing two treatments
    ##
    if (!is.na(idlabel.e)) {
      label.e <- args[[idlabel.e]]
      if (length(label.e) != 1)
        stop("Argument 'label.e' must be a character string.")
      ##
      setOption("label.e", label.e)
    }
    if (!is.na(idlabel.c)) {
      label.c <- args[[idlabel.c]]
      if (length(label.c) != 1)
        stop("Argument 'label.c' must be a character string.")
      ##
      setOption("label.c", label.c)
    }
    if (!is.na(idlabel.left)) {
      label.left <- args[[idlabel.left]]
      if (length(label.left) != 1)
        stop("Argument 'label.left' must be a character string.")
      ##
      setOption("label.left", label.left)
    }
    if (!is.na(idlabel.right)) {
      label.right <- args[[idlabel.right]]
      if (length(label.right) != 1)
        stop("Argument 'label.right' must be a character string.")
      ##
      setOption("label.right", label.right)
    }
    ##
    ## R function forest.meta
    ##
    if (!is.na(idlayout)) {
      layout <- args[[idlayout]]
      layout <- setchar(layout, layouts)
      setOption("layout", layout)
    }
    if (!is.na(idtest.overall)) {
      test.overall <- args[[idtest.overall]]
      chklogical(test.overall)
      setOption("test.overall", test.overall)
    }
    if (!is.na(idtest.subgroup)) {
      test.subgroup <- args[[idtest.subgroup]]
      chklogical(test.subgroup)
      setOption("test.subgroup", test.subgroup)
    }
    if (!is.na(idprediction.subgroup)) {
      prediction.subgroup <- args[[idprediction.subgroup]]
      chklogical(prediction.subgroup)
      setOption("prediction.subgroup", prediction.subgroup)
    }
    if (!is.na(idtest.effect.subgroup)) {
      test.effect.subgroup <- args[[idtest.effect.subgroup]]
      chklogical(test.effect.subgroup)
      setOption("test.effect.subgroup", test.effect.subgroup)
    }
    if (!is.na(iddigits.forest)) {
      digits.forest <- args[[iddigits.forest]]
      chknumeric(digits.forest, min = 0, length = 1)
      setOption("digits.forest", digits.forest)
    }
  }
  
  
  invisible(oldset)
}
