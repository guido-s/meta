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
#' meta-analysis results using \code{\link{print.meta}} and
#' \code{\link{print.summary.meta}}, and to produce forest plots using
#' \code{\link{forest.meta}}.
#' 
#' The function can be used to either change individual settings (see
#' Examples) or use one of the following general settings:
#' \itemize{
#' \item \code{settings.meta("revman5")}
#' \item \code{settings.meta("jama")}
#' \item \code{settings.meta("iqwig5")}
#' \item \code{settings.meta("iqwig6")}
#' \item \code{settings.meta("geneexpr")}
#' }
#'
#' The first command can be used to reproduce meta-analyses from
#' Cochrane reviews conducted with \emph{Review Manager 5} (RevMan 5,
#' \url{https://training.cochrane.org/online-learning/core-software-cochrane-reviews/revman})
#' and specifies to use a RevMan 5 layout in forest plots. The second
#' command can be used to generate forest plots following instructions
#' for authors of the \emph{Journal of the American Medical
#' Association}
#' (\url{https://jamanetwork.com/journals/jama/pages/instructions-for-authors/}). The
#' next commands implement the recommendations of the Institute for
#' Quality and Efficiency in Health Care, Germany (IQWiG) accordinging
#' to General Methods 5 and 6, respectively
#' (\url{https://www.iqwig.de/en/about-us/methods/methods-paper/}). The
#' last setting can be used to print p-values in scientific notation
#' and to suppress the calculation of confidence intervals for the
#' between-study variance.
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
#' A list of all arguments with current settings is printed using the
#' command \code{settings.meta("print")}.
#' 
#' In order to reset all settings of R package \bold{meta} the command
#' \code{settings.meta("reset")} can be used.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{gs}}, \code{\link{forest.meta}}
#' 
#' @examples
#' # Get listing of current settings
#' #
#' settings.meta("print")
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
#' settings.meta(level = 0.99, level.comb = 0.999)
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
#' forest(metagen(1:3, 2:4 / 10, sm = "MD", comb.fixed = FALSE),
#'        label.left = "Favours A", label.right = "Favours B",
#'        colgap.studlab = "2cm",
#'        colgap.forest.left = "0.2cm")
#' 
#' # Forest plot using JAMA style
#' #
#' settings.meta("jama")
#' forest(metagen(1:3, 2:4 / 10, sm = "MD", comb.fixed = FALSE),
#'        label.left = "Favours A", label.right = "Favours B",
#'        colgap.studlab = "2cm",
#'        colgap.forest.left = "0.2cm")
#'
#' # Use slightly different layout for confidence intervals
#' # (especially useful if upper confidence limit can be negative)
#' #
#' settings.meta(CIseparator = " - ")
#' forest(metagen(-(1:3), 2:4 / 10, sm="MD", comb.fixed=FALSE),
#'        label.left="Favours A", label.right="Favours B",
#'        colgap.studlab = "2cm",
#'        colgap.forest.left = "0.2cm")
#' 
#' # Use old settings
#' #
#' settings.meta(oldset)
#' 
#' @export settings.meta


settings.meta <- function(...) {
  
  ## Return current settings
  ##
  res <- .settings
  res$metafor <- NULL
  res$sm4bin <- res$sm4cont <- res$sm4cor <- res$sm4inc <-
    res$sm4mean <- res$sm4prop <- res$sm4rate <- NULL
  res$ci4cont <- res$ci4prop <- NULL
  res$meth4bin <- res$meth4inc <- res$meth4prop <- res$meth4rate <- NULL
  res$meth4tau <- res$meth4tau.ci <- NULL
  res$adhoc4hakn <- res$meth4bias <- NULL
  res$argslist <- NULL
  res$Wan2014.Table1 <- res$Wan2014.Table2 <- NULL
  res$digits.zval <- NULL
  
  
  catarg <- function(x, newline = TRUE, end = "") {
    xname <- x
    x <- gsub(" ", "", x)
    ##
    if (newline)
      cat("- ")
    ##
    if (is.null(.settings[[x]]))
      cat(paste0(xname, ' = NULL', end, '\n'))
    else if (is.character(.settings[[x]]))
      cat(paste0(xname, ' = "', .settings[[x]], '"', end, '\n'))
    else
      cat(paste0(xname, ' = ', .settings[[x]], end, "\n"))
    invisible(NULL)
  }
  
  
  specificSettings <- function(args, new, setting) {
    isnull.old <- as.vector(unlist(lapply(.settings[args], is.null)))
    ischar.old <- as.vector(unlist(lapply(.settings[args], is.character)))
    old <- as.vector(unlist(.settings[args]))
    ##
    ischar.new <- as.vector(unlist(lapply(new, is.character)))
    new <- as.vector(unlist(new))
    ##
    label.old <- ifelse(isnull.old, "NULL",
                        ifelse(ischar.old, paste0("\"", old, "\""), old))
    label.new <- ifelse(ischar.new, paste0("\"", new, "\""), new)
    ##
    sel <- new != old
    if (any(sel)) {
      tdata <- data.frame(argument = c("Argument",
                                       "--------",
                                       args[sel]),
                          space1 = rep("  ", along = c(1:2, sel)),
                          new = c("New value",
                                  "---------",
                                  label.new[sel]),
                          space2 = rep("  ", along = c(1:2, sel)),
                          previous = c("Previous value",
                                       "--------------",
                                       label.old[sel]))
      
      names(tdata) <- c("--------", "", "---------",
                        "", "--------------")
      ##
      cat(paste0("\n** Use ", setting, " (R package meta) **\n\n"))
      prmatrix(tdata, quote = FALSE, right = FALSE,
               rowlab = rep_len("", 2 + sum(sel)))
      ##
      for (i in seq(along = args)) {
        new.i <- new[i]
        if (!ischar.new[i]) {
          if (new.i %in% c("TRUE", "FALSE"))
            new.i <- as.logical(new.i)
          else
            new.i <- as.numeric(new.i)
        }
        setOption(args[i], new.i)
      }
    }
    else {
      if (substring(setting, 1, 1) == "s")
        setting <- paste0("S", substring(setting, 2))
      cat(paste0("\n** ", setting, " already in used (R package meta). **\n\n"))
    }
  }
  
  
  settings <- c("RevMan5", "JAMA", "IQWiG5", "IQWiG6", "geneexpr", "meta4")
  layouts <- c(settings[1:2], "meta")
  
  
  args  <- list(...)
  ## Check whether first argument is a list. In this case only use
  ## this list as input.
  if (length(args) > 0 && is.list(args[[1]]))
    args <- args[[1]]
  ##
  names <- names(args)
  
  
  for (i in seq_along(names))
    names[i] <- setchar(names[i],
                        c(.settings$argslist, "reset", "print", "setting", ""),
                        "unmatched",
                        name = names[i])
  ##
  if (length(names) != length(unique(names)))
    stop("Arguments must be unique.")
  
  
  if (length(args) > 1 && any(names == "reset")) {
    cat("To reset settings in R package meta use command 'settings.meta(\"reset\")'\n")
    return(invisible(res))
  }
  ##
  if (length(args) > 1 && any(names == "setting")) {
    cat("Argument 'setting' can only be used without other arguments (R package meta)\n")
    return(invisible(res))
  }
  
  
  print.settings <- FALSE
  reset.settings <- FALSE
  specific.settings <- FALSE
  ##
  if (length(args) == 1) {
    if (!is.null(names)) {
      if (names == "print") {
        chklogical(args[[1]], "print")
        print.settings <- args[[1]]
      }
      ##
      else if (names == "reset") {
        chklogical(args[[1]], "reset")
        if (args[[1]])
          reset.settings <- TRUE
        else {
          cat("To reset settings in R package meta use command 'settings.meta(\"reset\")'\n")
          return(invisible(res))
        }
      }
      ##
      else if (names == "setting") {
        setting <- setchar(args[[1]], settings, name = "setting")
        specific.settings <- TRUE
      }
    }
    else if (is.null(names) & is.character(args[[1]])) {
      ##
      idx <- charmatch(tolower(args[[1]]), "print", nomatch = NA)
      if (!(anyNA(idx) || any(idx == 0)))
        if (nchar(args[[1]]) >= 1)
          print.settings <- TRUE
      ##
      idx <- charmatch(tolower(args[[1]]), "reset", nomatch = NA)
      if (!(anyNA(idx) || any(idx == 0)))
        if (nchar(args[[1]]) >= 3)
          reset.settings <- TRUE
      ##
      if (!print.settings & !reset.settings) {
        setting <- setchar(args[[1]], settings, name = "setting")
        specific.settings <- TRUE
      }
    }
  }
  
  
  if (print.settings) {
    cat(paste0("\n** Settings for meta-analysis method (R package meta, ",
               "version ", utils::packageDescription("meta")$Version,
               ") **\n\n"))
    ##
    cat(paste0("* General settings *\n"))
    catarg("level          ")
    catarg("level.comb     ")
    catarg("comb.fixed     ")
    catarg("comb.random    ")
    catarg("hakn           ")
    catarg("adhoc.hakn     ")
    catarg("method.tau     ")
    catarg("method.tau.ci  ")
    catarg("tau.common     ")
    catarg("prediction     ")
    catarg("level.predict  ")
    catarg("method.bias    ")
    catarg("text.fixed     ")
    catarg("text.random    ")
    catarg("text.predict   ")
    catarg("text.w.fixed   ")
    catarg("text.w.random  ")
    catarg("title          ")
    catarg("complab        ")
    catarg("CIbracket      ")
    catarg("CIseparator    ")
    catarg("print.byvar    ")
    catarg("byseparator    ")
    catarg("keepdata       ")
    catarg("warn           ")
    catarg("backtransf     ")
    catarg("digits         ")
    catarg("digits.se      ")
    catarg("digits.stat    ")
    catarg("digits.Q       ")
    catarg("digits.tau2    ")
    catarg("digits.tau     ")
    catarg("digits.H       ")
    catarg("digits.I2      ")
    catarg("digits.prop    ")
    catarg("digits.weight  ")
    catarg("digits.pval    ")
    catarg("digits.pval.Q  ")
    catarg("scientific.pval")
    catarg("big.mark       ")
    catarg("zero.pval      ")
    catarg("JAMA.pval      ")
    catarg("print.I2       ")
    catarg("print.H        ")
    catarg("print.Rb       ")
    catarg("text.tau2      ")
    catarg("text.tau       ")
    catarg("text.I2        ")
    catarg("text.Rb        ")
    ##
    cat("\n* Default summary measure (argument 'sm' in corresponding function) *\n")
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
    cat("\n* Additional settings for metabin(), metainc(), metaprop(), and metarate() *\n")
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
    catarg("test.subgroup       ")
    catarg("test.effect.subgroup")
    catarg("digits.forest       ",
           end = "\n  (argument 'digits' in forest.meta())")
  }
  else if (reset.settings) {
    cat("\n** Reset all meta-analysis settings (R package meta). **\n\n")
    ##
    setOption("level", 0.95)
    setOption("level.comb", 0.95)
    setOption("comb.fixed", TRUE)
    setOption("comb.random", TRUE)
    setOption("hakn", FALSE)
    setOption("adhoc.hakn", "")
    setOption("method.tau", "DL")
    setOption("method.tau.ci", NULL)
    setOption("tau.common", FALSE)
    setOption("prediction", FALSE)
    setOption("level.predict", 0.95)
    setOption("method.bias", "linreg")
    setOption("text.fixed", "Fixed effect model")
    setOption("text.random", "Random effects model")
    setOption("text.predict", "Prediction interval")
    setOption("text.w.fixed", "fixed")
    setOption("text.w.random", "random")
    setOption("title", "")
    setOption("complab", "")
    setOption("CIbracket", "[")
    setOption("CIseparator", "; ")
    setOption("print.byvar", TRUE)
    setOption("byseparator", " = ")
    setOption("keepdata", TRUE)
    setOption("warn", TRUE)
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
    setOption("exact.smd", FALSE)
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
    setOption("test.subgroup", FALSE)
    setOption("test.effect.subgroup", FALSE)
    setOption("digits.forest", 2)
  }
  else if (specific.settings) {
    ##
    ## Remember:
    ## settings <- c("RevMan5", "JAMA", "IQWiG5", "IQWiG6", "geneexpr", "meta4")
    ##
    if (setting == "RevMan5") {
      specificSettings(args = c("hakn", "method.tau", "tau.common",
                                "MH.exact", "RR.Cochrane", "Q.Cochrane",
                                "layout", "test.overall",
                                "test.subgroup", "test.effect.subgroup",
                                "digits.I2", "digits.tau2", "digits.tau",
                                "CIbracket", "CIseparator"),
                       new = list(FALSE, "DL", FALSE,
                                  FALSE, TRUE, TRUE,
                                  "RevMan5", TRUE,
                                  TRUE, TRUE,
                                  0, 2, 4,
                                  "[", ", "),
                       setting = "RevMan 5 settings")
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
                       setting = "JAMA settings")
    }
    ##
    else if (setting == "IQWiG5") {
      specificSettings(args = c("hakn", "prediction"),
                       new = list(TRUE, TRUE),
                       setting = "IQWiG 5 settings")
    }
    ##
    else if (setting == "IQWiG6") {
      specificSettings(args = c("hakn", "adhoc.hakn",
                                "method.tau", "prediction"),
                       new = list(TRUE, "iqwig6", "PM", TRUE),
                       setting = "IQWiG 6 settings")
    }
    ##
    else if (setting == "meta4") {
      specificSettings(args = c("hakn", "method.tau",
                                "RR.Cochrane", "Q.Cochrane",
                                "CIbracket", "CIseparator"),
                       new = list(FALSE, "DL",
                                  FALSE, TRUE,
                                  "[", "; "),
                       setting = "Settings from R package meta (version 4.y-z and older)")
    }
    ##
    else if (setting == "geneexpr") {
      specificSettings(args = c("scientific.pval", "method.tau.ci"),
                       new = list(TRUE, ""),
                       setting = "Settings for gene expression data")
    }
  }
  else {
    argid <- function(x, value) {
      if (any(names == value))
        res <- seq(along = x)[names == value]
      else
        res <- NA
      res
    }
    ##
    idlevel <- argid(names, "level")
    idlevel.comb <- argid(names, "level.comb")
    idcomb.fixed <- argid(names, "comb.fixed")
    idcomb.random <- argid(names, "comb.random")
    idhakn <- argid(names, "hakn")
    idadhoc.hakn <- argid(names, "adhoc.hakn")
    idmethod.tau <- argid(names, "method.tau")
    idmethod.tau.ci <- argid(names, "method.tau.ci")
    idtau.common <- argid(names, "tau.common")
    idprediction <- argid(names, "prediction")
    idlevel.predict <- argid(names, "level.predict")
    idmethod.bias <- argid(names, "method.bias")
    idtext.fixed <- argid(names, "text.fixed")
    idtext.random <- argid(names, "text.random")
    idtext.predict <- argid(names, "text.predict")
    idtext.w.fixed <- argid(names, "text.w.fixed")
    idtext.w.random <- argid(names, "text.w.random")
    idtitle <- argid(names, "title")
    idcomplab <- argid(names, "complab")
    idCIbracket <- argid(names, "CIbracket")
    idCIseparator <- argid(names, "CIseparator")
    idprint.byvar <- argid(names, "print.byvar")
    idbyseparator <- argid(names, "byseparator")
    idkeepdata <- argid(names, "keepdata")
    idwarn <- argid(names, "warn")
    idbacktransf <- argid(names, "backtransf")
    iddigits <- argid(names, "digits")
    iddigits.se <- argid(names, "digits.se")
    iddigits.stat <- argid(names, "digits.stat")
    iddigits.Q <- argid(names, "digits.Q") 
    iddigits.tau2 <- argid(names, "digits.tau2")
    iddigits.tau <- argid(names, "digits.tau")
    iddigits.H <- argid(names, "digits.H") 
    iddigits.I2 <- argid(names, "digits.I2")
    iddigits.prop <- argid(names, "digits.prop")
    iddigits.weight <- argid(names,"digits.weight")
    iddigits.pval <- argid(names, "digits.pval")
    iddigits.pval.Q <- argid(names, "digits.pval.Q")
    idscientific.pval <- argid(names, "scientific.pval")
    idbig.mark <- argid(names, "big.mark")
    idzero.pval <- argid(names, "zero.pval")
    idJAMA.pval <- argid(names, "JAMA.pval")
    idprint.I2 <- argid(names, "print.I2")
    idprint.H <- argid(names, "print.H")
    idprint.Rb <- argid(names, "print.Rb")
    idtext.tau2 <- argid(names, "text.tau2")
    idtext.tau <- argid(names, "text.tau")
    idtext.I2 <- argid(names, "text.I2")
    idtext.Rb <- argid(names, "text.Rb")
    ##
    idsmbin <- argid(names, "smbin")
    idmethod <- argid(names, "method")
    idincr <- argid(names, "incr")
    idallincr <- argid(names, "allincr")
    idaddincr <- argid(names, "addincr")
    idallstudies <- argid(names, "allstudies")
    idMH.exact <- argid(names, "MH.exact")
    idRR.Cochrane <- argid(names, "RR.Cochrane")
    idQ.Cochrane <- argid(names, "Q.Cochrane")
    idmodel.glmm <- argid(names, "model.glmm")
    idprint.CMH <- argid(names, "print.CMH")
    ##
    idsmcont <- argid(names, "smcont")
    idsmcor <- argid(names, "smcor")
    idsminc <- argid(names, "sminc")
    idsmmean <- argid(names, "smmean")
    idsmprop <- argid(names, "smprop")
    idsmrate <- argid(names, "smrate")
    ##
    idpooledvar <- argid(names, "pooledvar")
    idmethod.smd <- argid(names, "method.smd")
    idsd.glass <- argid(names, "sd.glass")
    idexact.smd <- argid(names, "exact.smd")
    idmethod.ci.cont <- argid(names, "method.ci.cont")
    ##
    idmethod.ci.prop <- argid(names, "method.ci.prop")
    ##
    idlabel.e <- argid(names, "label.e")
    idlabel.c <- argid(names, "label.c")
    idlabel.left <- argid(names, "label.left")
    idlabel.right <- argid(names, "label.right")
    ##
    idlayout <- argid(names, "layout")
    idtest.overall <- argid(names, "test.overall")
    idtest.subgroup <- argid(names, "test.subgroup")
    idtest.effect.subgroup <- argid(names, "test.effect.subgroup")
    iddigits.forest <- argid(names, "digits.forest")
    ##
    ## General settings
    ##
    if (!is.na(idlevel)) {
      level <- args[[idlevel]]
      chklevel(level)
      setOption("level", level)
    }
    if (!is.na(idlevel.comb)) {
      level.comb <- args[[idlevel.comb]]
      chklevel(level.comb)
      setOption("level.comb", level.comb)
    }
    if (!is.na(idcomb.fixed)) {
      comb.fixed <- args[[idcomb.fixed]]
      chklogical(comb.fixed)
      setOption("comb.fixed", comb.fixed)
    }
    if (!is.na(idcomb.random)) {
      comb.random <- args[[idcomb.random]]
      chklogical(comb.random)
      setOption("comb.random", comb.random)
    }
    if (!is.na(idhakn)) {
      hakn <- args[[idhakn]]
      chklogical(hakn)
      setOption("hakn", hakn)
    }
    if (!is.na(idadhoc.hakn)) {
      adhoc.hakn <- args[[idadhoc.hakn]]
      adhoc.hakn <- setchar(adhoc.hakn, .settings$adhoc4hakn)
      setOption("adhoc.hakn", adhoc.hakn)
    }
    if (!is.na(idmethod.tau)) {
      method.tau <- args[[idmethod.tau]]
      method.tau <- setchar(method.tau, .settings$meth4tau)
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
      method.bias <- setchar(method.bias, .settings$meth4bias)
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
    if (!is.na(idprint.byvar)) {
      print.byvar <- args[[idprint.byvar]]
      chklogical(print.byvar)
      setOption("print.byvar", print.byvar)
    }
    if (!is.na(idbyseparator)) {
      byseparator <- args[[idbyseparator]]
      if (length(byseparator) != 1)
        stop("Argument 'byseparator' must be a character string.")
      ##
      setOption("byseparator", byseparator)
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
    if (!is.na(argid(names, "digits.zval"))) {
      warning("Argument 'digits.zval' is deprecated; ",
              "use instead argument 'digits.stat'.", call. = FALSE)
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
      smbin <- setchar(smbin, .settings$sm4bin)
      setOption("smbin", smbin)
    }
    if (!is.na(idmethod)) {
      method <- args[[idmethod]]
      method <- setchar(method, .settings$meth4bin)
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
      smcont <- setchar(smcont, .settings$sm4cont)
      setOption("smcont", smcont)
    }
    ##
    ## R function metacor
    ##
    if (!is.na(idsmcor)) {
      smcor <- args[[idsmcor]]
      smcor <- setchar(smcor, .settings$sm4cor)
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
      method.ci.cont <- setchar(method.ci.cont, .settings$ci4cont)
      setOption("method.ci.cont", method.ci.cont)
    }
    ##
    ## R function metainc
    ##
    if (!is.na(idsminc)) {
      sminc <- args[[idsminc]]
      sminc <- setchar(sminc, .settings$sm4inc)
      setOption("sminc", sminc)
    }
    ##
    ## R function metamean
    ##
    if (!is.na(idsmmean)) {
      smmean <- args[[idsmmean]]
      smmean <- setchar(smmean, .settings$sm4mean)
      setOption("smmean", smmean)
    }
    ##
    ## R function metaprop
    ##
    if (!is.na(idsmprop)) {
      smprop <- args[[idsmprop]]
      smprop <- setchar(smprop, .settings$sm4prop)
      setOption("smprop", smprop)
    }
    ##
    if (!is.na(idmethod.ci.prop)) {
      method.ci.prop <- args[[idmethod.ci.prop]]
      method.ci.prop <- setchar(method.ci.prop, .settings$ci4prop)
      setOption("method.ci.prop", method.ci.prop)
    }
    ##
    ## R function metarate
    ##
    if (!is.na(idsmrate)) {
      smrate <- args[[idsmrate]]
      smrate <- setchar(smrate, .settings$sm4rate)
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
  
  invisible(res)
}
