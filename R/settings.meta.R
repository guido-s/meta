settings.meta <- function(...){
  
  ## Return current settings
  ##
  res <- .settings
  res$CIbracket <- NULL
  res$CIseparator <- NULL
  res$argslist <- NULL
  
  
  catarg <- function(x, lab, newline=TRUE){
    if (!missing(lab))
      xname <- lab
    else
      xname <- x
    ##
    if (newline)
      cat("- ")
    ##
    if (is.character(.settings[[x]]))
      cat(xname, '="', .settings[[x]], '"\n', sep="")
    else
      cat(xname, '=', .settings[[x]], "\n", sep="")
    invisible(NULL)
  }
  
  
  specificSetting <- function(args, new, ischar, title){
    old <- as.vector(unlist(.settings[args]))
    label.old <- ifelse(ischar, paste("\"", old, "\"", sep=""), old)
    label.new <- ifelse(ischar, paste("\"", new, "\"", sep=""), new)
    ##
    sel <- new != old
    if (any(sel)){
      tdata <- data.frame(argument=args[sel],
                          new=label.new[sel],
                          previous=label.old[sel])
      names(tdata) <- c("Argument", "New value", "Previous value")
      ##
      cat(paste("\n** ", title, " **\n\n", sep=""))
      prmatrix(tdata, quote=FALSE, right=FALSE,
               rowlab=rep("", length(args)))
      ##
      for (i in seq(along=args))
        setOption(args[i], new[i])
    }
    else
      cat("No setting changed.\n")
  }
  
  
  args  <- list(...)
  ## Check whether first argument is a list. In this case only use
  ## this list as input.
  if (length(args)>0 && is.list(args[[1]]))
    args <- args[[1]]
  ##
  names <- names(args)
  
  
  if (length(names) != length(unique(names)))
    stop("Arguments must be unique.")
  
  
  unknown <- !(names %in% c(.settings$argslist, "reset", "print", "setting", ""))
  ##
  if (sum(unknown)==1)
    warning(paste("Argument '", names[unknown], "' unknown.", sep=""))
  else if (sum(unknown)>1)
    warning(paste("Unknown arguments: ", 
                  paste(paste("'", names[unknown], "'", sep=""),
                        collapse=" - "), sep=""))
  
  
  if (length(args)>1 && any(names=="reset")){
    cat("To reset all settings use a single argument 'reset=TRUE' (R package meta)\n")
    return(invisible(res))
  }
  ##
  if (length(args)>1 && any(names=="setting")){
    cat("Argument 'setting' can only be used without other arguments (R package meta)\n")
    return(invisible(res))
  }
  
  
  print.settings <- FALSE
  reset.settings <- FALSE
  specific.settings <- FALSE
  ##
  if (length(args)==1){
    if (!is.null(names)){
      if (names=="print"){
        chklogical(args[[1]], "print")
        print.settings <- args[[1]]
      }
      ##
      else if (names=="reset"){
        chklogical(args[[1]], "reset")
        if (args[[1]])
          reset.settings <- TRUE
        else{
          cat("To reset all settings use argument 'reset=TRUE' (R package meta)\n")
          return(invisible(res))
        }
      }
      ##
      else if (names=="setting"){
        setting <- setchar(args[[1]], c("RevMan5", "IQWiG", "meta4"), name="setting")
      specific.settings <- TRUE
      }
    }
    else if (is.null(names) & is.character(args[[1]])){
      setting <- setchar(args[[1]], c("RevMan5", "IQWiG", "meta4"), name="setting")
      specific.settings <- TRUE
    }
  }
  
  
  if (print.settings){
    cat(paste("\n*** Settings for meta-analysis method (R package meta, version ",
              utils::packageDescription("meta")$Version, ") ***\n\n", sep=""))
    ##
    cat("General settings:\n", sep="")
    catarg("level")
    catarg("level.comb")
    catarg("comb.fixed")
    catarg("comb.random")
    catarg("hakn")
    catarg("method.tau")
    catarg("tau.common")
    catarg("prediction")
    catarg("level.predict")
    catarg("method.bias")
    catarg("title")
    catarg("complab")
    catarg("print.byvar")
    catarg("keepdata")
    catarg("warn")
    catarg("backtransf")
    ##
    cat("\nDefault summary measure (argument 'sm' in corresponding function):\n")
    cat("- metabin:  ")
    catarg("smbin", newline=FALSE)
    cat("- metacont: ")
    catarg("smcont", newline=FALSE)
    cat("- metacor:  ")
    catarg("smcor", newline=FALSE)
    cat("- metainc:  ")
    catarg("sminc", newline=FALSE)
    cat("- metaprop: ")
    catarg("smprop", newline=FALSE)
    ##
    cat("\nSettings for R functions metabin, metainc, and metaprop:\n")
    catarg("incr")
    catarg("allincr")
    catarg("addincr")
    ##
    cat("\nAdditional settings for R function metabin:\n")
    catarg("method")
    catarg("allstudies")
    catarg("MH.exact")
    catarg("RR.cochrane")
    catarg("print.CMH")
    ##
    cat("\nAdditional setting for R function metacont:\n")
    catarg("pooledvar")
    catarg("method.smd")
    catarg("sd.glass")
    catarg("exact.smd")
    ##
    cat("\nAdditional setting for R function metaprop:\n")
    catarg("method.ci")
    ##
    cat("\nSettings for R functions comparing two treatments:\n")
    catarg("label.e")
    catarg("label.c")
    catarg("label.left")
    catarg("label.right")
    ##
    cat("\nSettings for R function forest.meta:\n")
    catarg("test.overall")
    catarg("test.subgroup")
  }
  else if (reset.settings){
    cat("Reset all settings back to default (R package meta).\n")
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
    setOption("print.byvar", TRUE)
    setOption("keepdata", TRUE)
    setOption("warn", TRUE)
    setOption("backtransf", TRUE)
    ##
    setOption("method", "MH")
    setOption("incr", 0.5)
    setOption("allincr", FALSE)
    setOption("addincr", FALSE)
    setOption("allstudies", FALSE)
    setOption("MH.exact", FALSE)
    setOption("RR.cochrane", FALSE)
    setOption("print.CMH", FALSE)
    ##
    setOption("smbin", "RR")
    setOption("smcont", "MD")
    setOption("smcor", "ZCOR")
    setOption("sminc", "IRR")
    setOption("smprop", "PLOGIT")
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
    setOption("test.overall", FALSE)
    setOption("test.subgroup", FALSE)
  }
  else if (specific.settings){
    if (setting=="RevMan5")
      specificSetting(args=c("RR.cochrane"),
                      new=TRUE,
                      ischar=FALSE,
                      title="Use RevMan 5 settings")
    ##
    else if (setting=="IQWiG")
      specificSetting(args=c("hakn", "method.tau"),
                      new=c(TRUE, "PM"),
                      ischar=c(FALSE, TRUE),
                      title="Use IQWiG settings")
    ##
    else if (setting=="meta4")
      specificSetting(args=c("hakn", "method.tau", "RR.cochrane"),
                      new=c(FALSE, "DL", FALSE),
                      ischar=c(FALSE, TRUE, FALSE),
                      title="Use settings from R package meta (version 4.3-1 and older)")
  }
  else{
    argid <- function(x, value){
      if (any(names==value))
        res <- seq(along=x)[names==value]
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
    idprint.byvar <- argid(names, "print.byvar")
    idkeepdata <- argid(names, "keepdata")
    idwarn <- argid(names, "warn")
    idbacktransf <- argid(names, "backtransf")
    ##
    idsmbin <- argid(names, "smbin")
    idmethod <- argid(names, "method")
    idincr <- argid(names, "incr")
    idallincr <- argid(names, "allincr")
    idaddincr <- argid(names, "addincr")
    idallstudies <- argid(names, "allstudies")
    idMH.exact <- argid(names, "MH.exact")
    idRR.cochrane <- argid(names, "RR.cochrane")
    idprint.CMH <- argid(names, "print.CMH")
    ##
    idsmcont <- argid(names, "smcont")
    idsmcor <- argid(names, "smcor")
    idsminc <- argid(names, "sminc")
    idsmprop <- argid(names, "smprop")
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
    idtest.overall <- argid(names, "test.overall")
    idtest.subgroup <- argid(names, "test.subgroup")
    ##
    ## General settings
    ##
    if (!is.na(idlevel)){
      level <- args[[idlevel]]
      chklevel(level)
      setOption("level", level)
    }
    if (!is.na(idlevel.comb)){
      level.comb <- args[[idlevel.comb]]
      chklevel(level.comb)
      setOption("level.comb", level.comb)
    }
    if (!is.na(idcomb.fixed)){
      comb.fixed <- args[[idcomb.fixed]]
      chklogical(comb.fixed)
      setOption("comb.fixed", comb.fixed)
    }
    if (!is.na(idcomb.random)){
      comb.random <- args[[idcomb.random]]
      chklogical(comb.random)
      setOption("comb.random", comb.random)
    }
    if (!is.na(idhakn)){
      hakn <- args[[idhakn]]
      chklogical(hakn)
      setOption("hakn", hakn)
    }
    if (!is.na(idmethod.tau)){
      method.tau <- args[[idmethod.tau]]
      method.tau <- setchar(method.tau,
                           c("DL", "PM", "REML", "ML", "HS", "SJ", "HE", "EB"))
      setOption("method.tau", method.tau)
    }
    if (!is.na(idtau.common)){
      tau.common <- args[[idtau.common]]
      chklogical(tau.common)
      setOption("tau.common", tau.common)
    }
    if (!is.na(idprediction)){
      prediction <- args[[idprediction]]
      chklogical(prediction)
      setOption("prediction", prediction)
    }
    if (!is.na(idlevel.predict)){
      level.predict <- args[[idlevel.predict]]
      chklevel(level.predict)
      setOption("level.predict", level.predict)
    }
    if (!is.na(idmethod.bias)){
      method.bias <- args[[idmethod.bias]]
      method.bias <- setchar(method.bias,
                            c("rank", "linreg", "mm", "count", "score", "peters"))
      setOption("method.bias", method.bias)
    }
    if (!is.na(idtitle)){
      title <- args[[idtitle]]
      if (length(title)!=1)
        stop("Argument 'title' must be a character string.")
      ##
      setOption("title", title)
    }
    if (!is.na(idcomplab)){
      complab <- args[[idcomplab]]
      if (length(complab)!=1)
        stop("Argument 'complab' must be a character string.")
      ##
      setOption("complab", complab)
    }
    if (!is.na(idprint.byvar)){
      print.byvar <- args[[idprint.byvar]]
      chklogical(print.byvar)
      setOption("print.byvar", print.byvar)
    }
    if (!is.na(idkeepdata)){
      keepdata <- args[[idkeepdata]]
      chklogical(keepdata)
      setOption("keepdata", keepdata)
    }
    if (!is.na(idwarn)){
      warn <- args[[idwarn]]
      chklogical(warn)
      setOption("warn", warn)
    }
    if (!is.na(idbacktransf)){
      backtransf <- args[[idbacktransf]]
      chklogical(backtransf)
      setOption("backtransf", backtransf)
    }
    ##
    ## R function metabin
    ##
    if (!is.na(idsmbin)){
      smbin <- args[[idsmbin]]
      smbin <- setchar(smbin, c("OR", "RD", "RR", "ASD"))
      setOption("smbin", smbin)
    }
    if (!is.na(idmethod)){
      method <- args[[idmethod]]
      method <- setchar(method, c("Inverse", "MH", "Peto"))
      setOption("method", method)
    }
    if (!is.na(idincr)){
      incr <- args[[idincr]]
      if (!is.numeric(incr))
        incr <- setchar(incr, "TACC",
                       "should be numeric or the character string \"TACC\"")
      setOption("incr", incr)
    }
    if (!is.na(idallincr)){
      allincr <- args[[idallincr]]
      chklogical(allincr)
      setOption("allincr", allincr)
    }
    if (!is.na(idaddincr)){
      addincr <- args[[idaddincr]]
      chklogical(addincr)
      setOption("addincr", addincr)
    }
    if (!is.na(idallstudies)){
      allstudies <- args[[idallstudies]]
      chklogical(allstudies)
      setOption("allstudies", allstudies)
    }
    if (!is.na(idMH.exact)){
      MH.exact <- args[[idMH.exact]]
      chklogical(MH.exact)
      setOption("MH.exact", MH.exact)
    }
    if (!is.na(idRR.cochrane)){
      RR.cochrane <- args[[idRR.cochrane]]
      chklogical(RR.cochrane)
      setOption("RR.cochrane", RR.cochrane)
    }
    if (!is.na(idprint.CMH)){
      print.CMH <- args[[idprint.CMH]]
      chklogical(print.CMH)
      setOption("print.CMH", print.CMH)
    }
    ##
    ## R function metacont
    ##
    if (!is.na(idsmcont)){
      smcont <- args[[idsmcont]]
      smcont <- setchar(smcont, c("MD", "SMD"))
      setOption("smcont", smcont)
    }
    ##
    ## R function metacor
    ##
    if (!is.na(idsmcor)){
      smcor <- args[[idsmcor]]
      smcor <- setchar(smcor, c("ZCOR", "COR"))
      setOption("smcor", smcor)
    }
    ##
    ## R function metainc
    ##
    if (!is.na(idsminc)){
      sminc <- args[[idsminc]]
      sminc <- setchar(sminc, c("IRR", "IRD"))
      setOption("sminc", sminc)
    }
    ##
    ## R function metacont
    ##
    if (!is.na(idpooledvar)){
      pooledvar <- args[[idpooledvar]]
      chklogical(pooledvar)
      setOption("pooledvar", pooledvar)
    }
    if (!is.na(idmethod.smd)){
      method.smd <- args[[idmethod.smd]]
      method.smd <- setchar(method.smd, c("Hedges", "Cohen", "Glass"))
      setOption("method.smd", method.smd)
    }
    if (!is.na(idsd.glass)){
      sd.glass <- args[[idsd.glass]]
      sd.glass <- setchar(sd.glass, c("control", "experimental"))
      setOption("sd.glass", sd.glass)
    }
    if (!is.na(idexact.smd)){
      exact.smd <- args[[idexact.smd]]
      chklogical(exact.smd)
      setOption("exact.smd", exact.smd)
    }
    ##
    ## R function metaprop
    ##
    if (!is.na(idsmprop)){
      smprop <- args[[idsmprop]]
      smprop <- setchar(smprop, c("PFT", "PAS", "PRAW", "PLN", "PLOGIT"))
      setOption("smprop", smprop)
    }
    ##
    if (!is.na(idmethod.ci)){
      method.ci <- args[[idmethod.ci]]
      method.ci <- setchar(method.ci,
                          c("CP", "WS", "WSCC", "AC", "SA", "SACC", "NAsm"))
      setOption("method.ci", method.ci)
    }
    ##
    ## R functions comparing two treatments
    ##
    if (!is.na(idlabel.e)){
      label.e <- args[[idlabel.e]]
      if (length(label.e)!=1)
        stop("Argument 'label.e' must be a character string.")
      ##
      setOption("label.e", label.e)
    }
    if (!is.na(idlabel.c)){
      label.c <- args[[idlabel.c]]
      if (length(label.c)!=1)
        stop("Argument 'label.c' must be a character string.")
      ##
      setOption("label.c", label.c)
    }
    if (!is.na(idlabel.left)){
      label.left <- args[[idlabel.left]]
      if (length(label.left)!=1)
        stop("Argument 'label.left' must be a character string.")
      ##
      setOption("label.left", label.left)
    }
    if (!is.na(idlabel.right)){
      label.right <- args[[idlabel.right]]
      if (length(label.right)!=1)
        stop("Argument 'label.right' must be a character string.")
      ##
      setOption("label.right", label.right)
    }
    ##
    ## R function forest.meta
    ##
    if (!is.na(idtest.overall)){
      test.overall <- args[[idtest.overall]]
      chklogical(test.overall)
      setOption("test.overall", test.overall)
    }
    if (!is.na(idtest.subgroup)){
      test.subgroup <- args[[idtest.subgroup]]
      chklogical(test.subgroup)
      setOption("test.subgroup", test.subgroup)
    }
  }
  
  invisible(res)
}
