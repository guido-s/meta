settings.meta <- function(...) {
  
  ## Return current settings
  ##
  res <- .settings
  res$metafor <- NULL
  res$argslist <- NULL
  
  
  catarg <- function(x, newline = TRUE, end = "") {
    xname <- x
    x <- gsub(" ", "", x)
    ##
    if (newline)
      cat("- ")
    ##
    if (is.character(.settings[[x]]))
      cat(xname, ' = "', .settings[[x]], '"', end, '\n', sep = "")
    else
      cat(xname, ' = ', .settings[[x]], end, "\n", sep = "")
    invisible(NULL)
  }
  
  
  specificSetting <- function(args, new, setting) {
    old <- as.vector(unlist(.settings[args]))
    ischar <- as.vector(unlist(lapply(new, is.character)))
    new <- as.vector(unlist(new))
    ##
    label.old <- ifelse(ischar, paste("\"", old, "\"", sep = ""), old)
    label.new <- ifelse(ischar, paste("\"", new, "\"", sep = ""), new)
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
      cat(paste("\n** Use ", setting, " (R package meta) **\n\n", sep = ""))
      prmatrix(tdata, quote = FALSE, right = FALSE,
               rowlab = rep_len("", 2 + sum(sel)))
      ##
      for (i in seq(along = args)) {
        new.i <- new[i]
        if (!ischar[i]) {
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
        setting <- paste("S", substring(setting, 2), sep = "")
      cat(paste("\n** ", setting, " already in used (R package meta). **\n\n", sep = ""))
    }
  }
  
  
  settings <- c("RevMan5", "JAMA", "IQWiG4.2", "meta4")
  layouts <- c(settings[1:2], "meta")
  
  
  args  <- list(...)
  ## Check whether first argument is a list. In this case only use
  ## this list as input.
  if (length(args) > 0 && is.list(args[[1]]))
    args <- args[[1]]
  ##
  names <- names(args)
  
  
  if (length(names) != length(unique(names)))
    stop("Arguments must be unique.")
  
  
  unknown <- !(names %in% c(.settings$argslist, "reset", "print", "setting", ""))
  ##
  if (sum(unknown) == 1)
    warning(paste("Argument '", names[unknown], "' unknown.", sep = ""))
  else if (sum(unknown) > 1)
    warning(paste("Unknown arguments: ", 
                  paste(paste("'", names[unknown], "'", sep = ""),
                        collapse = " - "), sep = ""))
  
  
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
      if (!(any(is.na(idx)) || any(idx == 0)))
        if (nchar(args[[1]]) >= 1)
          print.settings <- TRUE
      ##
      idx <- charmatch(tolower(args[[1]]), "reset", nomatch = NA)
      if (!(any(is.na(idx)) || any(idx == 0)))
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
    cat(paste("\n** Settings for meta-analysis method (R package meta, version ",
              utils::packageDescription("meta")$Version, ") **\n\n", sep = ""))
    ##
    cat("* General settings *\n", sep = "")
    catarg("level        ")
    catarg("level.comb   ")
    catarg("comb.fixed   ")
    catarg("comb.random  ")
    catarg("hakn         ")
    catarg("method.tau   ")
    catarg("tau.common   ")
    catarg("prediction   ")
    catarg("level.predict")
    catarg("method.bias  ")
    catarg("title        ")
    catarg("complab      ")
    catarg("CIbracket    ")
    catarg("CIseparator  ")
    catarg("print.byvar  ")
    catarg("byseparator  ")
    catarg("keepdata     ")
    catarg("warn         ")
    catarg("backtransf   ")
    catarg("digits       ")
    catarg("digits.se    ")
    catarg("digits.zval  ")
    catarg("digits.Q     ")
    catarg("digits.tau2  ")
    catarg("digits.H     ")
    catarg("digits.I2    ")
    catarg("digits.prop  ")
    catarg("digits.weight")
    catarg("digits.pval  ")
    catarg("digits.pval.Q")
    catarg("scientific.pval")
    catarg("print.I2")
    catarg("print.H")
    catarg("print.Rb")
    catarg("text.tau2")
    catarg("text.I2")
    catarg("text.Rb")
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
    catarg("RR.cochrane")
    catarg("model.glmm ")
    catarg("print.CMH  ")
    ##
    cat("\n* Additional settings for metacont() *\n")
    catarg("pooledvar ")
    catarg("method.smd")
    catarg("sd.glass  ")
    catarg("exact.smd ")
    ##
    cat("\n* Additional setting for metaprop() *\n")
    catarg("method.ci")
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
    setOption("method.tau", "DL")
    setOption("tau.common", FALSE)
    setOption("prediction", FALSE)
    setOption("level.predict", 0.95)
    setOption("method.bias", "linreg")
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
    setOption("digits.zval", 2)
    setOption("digits.Q", 2)
    setOption("digits.tau2", 4)
    setOption("digits.H", 2)
    setOption("digits.I2", 1)
    setOption("digits.prop", 4)
    setOption("digits.weight", 1)
    setOption("digits.pval", 4)
    setOption("digits.pval.Q", 4)
    setOption("scientific.pval", FALSE)
    setOption("print.I2", TRUE)
    setOption("print.H", TRUE)
    setOption("print.Rb", FALSE)
    setOption("text.tau2", "tau^2")
    setOption("text.I2", "I^2")
    setOption("text.Rb", "Rb")
    ##
    setOption("method", "MH")
    setOption("incr", 0.5)
    setOption("allincr", FALSE)
    setOption("addincr", FALSE)
    setOption("allstudies", FALSE)
    setOption("MH.exact", FALSE)
    setOption("RR.cochrane", FALSE)
    setOption("model.glmm", "UM.FS")
    setOption("print.CMH", FALSE)
    ##
    setOption("smbin", "RR")
    setOption("smcont", "MD")
    setOption("smcor", "ZCOR")
    setOption("sminc", "IRR")
    setOption("smprop", "PLOGIT")
    setOption("smrate", "IRLN")
    ##
    setOption("pooledvar", FALSE)
    setOption("method.smd", "Hedges")
    setOption("sd.glass", "control")
    setOption("exact.smd", FALSE)
    ##
    setOption("method.ci", "CP")
    ##
    setOption("label.e", "Experimental")
    setOption("label.c", "Control")
    setOption("label.left", "")
    setOption("label.right", "")
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
    ## settings <- c("RevMan5", "JAMA", "IQWiG4.2", "meta4")
    ##
    if (setting == "RevMan5") {
      specificSetting(args = c("hakn", "method.tau", "tau.common",
                               "MH.exact", "RR.cochrane",
                               "layout", "test.overall",
                               "test.subgroup", "test.effect.subgroup",
                               "digits.I2", "digits.tau2",
                               "CIbracket", "CIseparator"),
                      new = list(FALSE, "DL", FALSE,
                                 FALSE, TRUE,
                                 "RevMan5", TRUE,
                                 TRUE, TRUE,
                                 0, 2,
                                 "[", ", "),
                      setting = "RevMan 5 settings")
    }
    ##
    else if (setting == "JAMA") {
      specificSetting(args = c("layout", "test.overall",
                               "test.subgroup", "test.effect.subgroup",
                               "digits.I2",
                               "CIbracket", "CIseparator"),
                      new = list("JAMA", TRUE,
                                 FALSE, FALSE,
                                 0,
                                 "(", "-"),
                      setting = "JAMA settings")
    }
    ##
    else if (setting == "IQWiG4.2") {
      specificSetting(args = c("hakn", "method.tau", "RR.cochrane"),
                      new = list(TRUE, "PM", FALSE),
                      setting = "IQWiG 4.2 settings")
    }
    ##
    else if (setting == "meta4") {
      specificSetting(args = c("hakn", "method.tau", "RR.cochrane",
                               "CIbracket", "CIseparator"),
                      new = list(FALSE, "DL", FALSE,
                                 "[", "; "),
                      setting = "settings from R package meta (version 4.y-z and older)")
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
    idmethod.tau <- argid(names, "method.tau")
    idtau.common <- argid(names, "tau.common")
    idprediction <- argid(names, "prediction")
    idlevel.predict <- argid(names, "level.predict")
    idmethod.bias <- argid(names, "method.bias")
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
    iddigits.zval <- argid(names, "digits.zval")
    iddigits.Q <- argid(names, "digits.Q") 
    iddigits.tau2 <- argid(names, "digits.tau2")
    iddigits.H <- argid(names, "digits.H") 
    iddigits.I2 <- argid(names, "digits.I2")
    iddigits.prop <- argid(names, "digits.prop")
    iddigits.weight <- argid(names,"digits.weight")
    iddigits.pval <- argid(names, "digits.pval")
    iddigits.pval.Q <- argid(names, "digits.pval.Q")
    idscientific.pval <- argid(names, "scientific.pval")
    idprint.I2 <- argid(names, "print.I2")
    idprint.H <- argid(names, "print.H")
    idprint.Rb <- argid(names, "print.Rb")
    idtext.tau2 <- argid(names, "text.tau2")
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
    idRR.cochrane <- argid(names, "RR.cochrane")
    idmodel.glmm <- argid(names, "model.glmm")
    idprint.CMH <- argid(names, "print.CMH")
    ##
    idsmcont <- argid(names, "smcont")
    idsmcor <- argid(names, "smcor")
    idsminc <- argid(names, "sminc")
    idsmprop <- argid(names, "smprop")
    idsmrate <- argid(names, "smrate")
    ##
    idpooledvar <- argid(names, "pooledvar")
    idmethod.smd <- argid(names, "method.smd")
    idsd.glass <- argid(names, "sd.glass")
    idexact.smd <- argid(names, "exact.smd")
    ##
    idmethod.ci <- argid(names, "method.ci")
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
    if (!is.na(idmethod.tau)) {
      method.tau <- args[[idmethod.tau]]
      method.tau <- setchar(method.tau,
                            c("DL", "PM", "REML", "ML", "HS", "SJ", "HE", "EB"))
      if (method.tau %in% c("REML", "ML", "HS", "SJ", "HE", "EB"))
        is.installed.package("metafor", chksettings = TRUE,
                             argument = "method.tau", value = method.tau)
      setOption("method.tau", method.tau)
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
      method.bias <- setchar(method.bias,
                            c("rank", "linreg", "mm", "count", "score", "peters"))
      setOption("method.bias", method.bias)
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
      chknumeric(digits, min = 0, single = TRUE)
      setOption("digits", digits)
    }
    if (!is.na(iddigits.se)) {
      digits.se <- args[[iddigits.se]]
      chknumeric(digits.se, min = 0, single = TRUE)
      setOption("digits.se", digits.se)
    }
    if (!is.na(iddigits.zval)) {
      digits.zval <- args[[iddigits.zval]]
      chknumeric(digits.zval, min = 0, single = TRUE)
      setOption("digits.zval", digits.zval)
    }
    if (!is.na(iddigits.Q)) {
      digits.Q <- args[[iddigits.Q]]
      chknumeric(digits.Q, min = 0, single = TRUE)
      setOption("digits.Q", digits.Q)
    }
    if (!is.na(iddigits.tau2)) {
      digits.tau2 <- args[[iddigits.tau2]]
      chknumeric(digits.tau2, min = 0, single = TRUE)
      setOption("digits.tau2", digits.tau2)
    }
    if (!is.na(iddigits.H)) {
      digits.H <- args[[iddigits.H]]
      chknumeric(digits.H, min = 0, single = TRUE)
      setOption("digits.H", digits.H)
    }
    if (!is.na(iddigits.I2)) {
      digits.I2 <- args[[iddigits.I2]]
      chknumeric(digits.I2, min = 0, single = TRUE)
      setOption("digits.I2", digits.I2)
    }
    if (!is.na(iddigits.prop)) {
      digits.prop <- args[[iddigits.prop]]
      chknumeric(digits.prop, min = 0, single = TRUE)
      setOption("digits.prop", digits.prop)
    }
    if (!is.na(iddigits.weight)) {
      digits.weight <- args[[iddigits.weight]]
      chknumeric(digits.weight, min = 0, single = TRUE)
      setOption("digits.weight", digits.weight)
    }
    if (!is.na(iddigits.pval)) {
      digits.pval <- args[[iddigits.pval]]
      chknumeric(digits.pval, min = 0, single = TRUE)
      setOption("digits.pval", digits.pval)
    }
    if (!is.na(iddigits.pval.Q)) {
      digits.pval.Q <- args[[iddigits.pval.Q]]
      chknumeric(digits.pval.Q, min = 0, single = TRUE)
      setOption("digits.pval.Q", digits.pval.Q)
    }
    if (!is.na(idscientific.pval)) {
      scientific.pval <- args[[idscientific.pval]]
      chklogical(scientific.pval)
      setOption("scientific.pval", scientific.pval)
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
      smbin <- setchar(smbin, c("OR", "RD", "RR", "ASD"))
      setOption("smbin", smbin)
    }
    if (!is.na(idmethod)) {
      method <- args[[idmethod]]
      method <- setchar(method, c("Inverse", "MH", "Peto", "GLMM"))
      if (method == "GLMM") {
        is.installed.package("metafor", chksettings = TRUE,
                             argument = "method", value = method)
        is.installed.package("lme4", chksettings = TRUE,
                             argument = "method", value = method)
        is.installed.package("numDeriv", chksettings = TRUE,
                             argument = "method", value = method)
      }
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
    if (!is.na(idRR.cochrane)) {
      RR.cochrane <- args[[idRR.cochrane]]
      chklogical(RR.cochrane)
      setOption("RR.cochrane", RR.cochrane)
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
      smcont <- setchar(smcont, c("MD", "SMD", "ROM"))
      setOption("smcont", smcont)
    }
    ##
    ## R function metacor
    ##
    if (!is.na(idsmcor)) {
      smcor <- args[[idsmcor]]
      smcor <- setchar(smcor, c("ZCOR", "COR"))
      setOption("smcor", smcor)
    }
    ##
    ## R function metainc
    ##
    if (!is.na(idsminc)) {
      sminc <- args[[idsminc]]
      sminc <- setchar(sminc, c("IRR", "IRD"))
      setOption("sminc", sminc)
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
    ## R function metaprop
    ##
    if (!is.na(idsmprop)) {
      smprop <- args[[idsmprop]]
      smprop <- setchar(smprop, c("PFT", "PAS", "PRAW", "PLN", "PLOGIT"))
      setOption("smprop", smprop)
    }
    ##
    if (!is.na(idmethod.ci)) {
      method.ci <- args[[idmethod.ci]]
      method.ci <- setchar(method.ci,
                          c("CP", "WS", "WSCC", "AC", "SA", "SACC", "NAsm"))
      setOption("method.ci", method.ci)
    }
    ##
    ## R function metarate
    ##
    if (!is.na(idsmrate)) {
      smrate <- args[[idsmrate]]
      smrate <- setchar(smrate, c("IR", "IRLN", "IRS", "IRFT"))
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
      chknumeric(digits.forest, min = 0, single = TRUE)
      setOption("digits.forest", digits.forest)
    }
  }
  
  invisible(res)
}
