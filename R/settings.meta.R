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
  
  
  args  <- list(...)
  ## Check whether first argument is a list. In this case only use
  ## this list as input.
  if (length(args)>0 && is.list(args[[1]]))
    args <- args[[1]]
  ##
  names <- names(args)
  
  
  if (length(names) != length(unique(names)))
    stop("Arguments must be unique.")
  
  
  unknown <- !(names %in% c(.settings$argslist, "reset", "print"))
  ##
  if (sum(unknown)==1)
    warning(paste("Argument '", names[unknown], "' unknown.", sep=""))
  else if (sum(unknown)>1)
    warning(paste("Unknown arguments: ", 
                  paste(paste("'", names[unknown], "'", sep=""),
                        collapse=" - "), sep=""))
  

  if (length(args)==1 && names=="print"){
    cat(paste("\n*** Settings for meta-analysis method (R package meta, version ",
              utils::packageDescription("meta")$Version, ") ***\n\n", sep=""))
    ##cat("General settings (i.e. for R functions metabin, metacont,\n",
    ##    "                  metacor, metagen, metainc, metaprop):\n", sep="")
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
    ##
    cat("\nAdditional setting for R function metaprop:\n")
    catarg("method.ci")
    ##
    cat("\nSettings for R functions comparing two treatments:\n")
    catarg("label.e")
    catarg("label.c")
    catarg("label.left")
    catarg("label.right")
  }
  else if (length(args)==1 && names=="reset"){
    if (is.logical(args[[1]]) && args[[1]]==TRUE){
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
      ##
      setOption("method.ci", "CP")
      ##
      setOption("label.e", "Experimental")
      setOption("label.c", "Control")
      setOption("label.left", "")
      setOption("label.right", "")
    }
    else
      cat("To reset all settings use argument 'reset=TRUE' (R package meta)\n")
  }
  else if (length(args)>1 && any(names=="reset"))
    cat("To reset all settings use a single argument 'reset=TRUE' (R package meta)\n")
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
    ##
    idmethod.ci <- argid(names, "method.ci")
    ##
    idlabel.e <- argid(names, "label.e")
    idlabel.c <- argid(names, "label.c")
    idlabel.left <- argid(names, "label.left")
    idlabel.right <- argid(names, "label.right")
    ##
    ## General settings
    ##
    if (!is.na(idlevel)){
      level <- args[[idlevel]]
      if (!is.numeric(level) | length(level)!=1)
        stop("Argument 'level' must be a numeric of length 1.")
      if (level <= 0 | level >= 1)
        stop("Argument 'level': no valid level for confidence interval")
      ##
      setOption("level", level)
    }
    if (!is.na(idlevel.comb)){
      level.comb <- args[[idlevel.comb]]
      if (!is.numeric(level.comb) | length(level.comb)!=1)
        stop("Argument 'level.comb' must be a numeric of length 1.")
      if (level.comb <= 0 | level.comb >= 1)
        stop("Argument 'level.comb': no valid level for confidence interval")
      ##
      setOption("level.comb", level.comb)
    }
    if (!is.na(idcomb.fixed)){
      comb.fixed <- args[[idcomb.fixed]]
      if (length(comb.fixed)!= 1 || !is.logical(comb.fixed))
        stop("Argument 'comb.fixed' must be a logical.")
      ##
      setOption("comb.fixed", comb.fixed)
    }
    if (!is.na(idcomb.random)){
      comb.random <- args[[idcomb.random]]
      if (length(comb.random)!= 1 || !is.logical(comb.random))
        stop("Argument 'comb.random' must be a logical.")
      ##
      setOption("comb.random", comb.random)
    }
    if (!is.na(idhakn)){
      hakn <- args[[idhakn]]
      if (length(hakn)!= 1 || !is.logical(hakn))
        stop("Argument 'hakn' must be a logical.")
      ##
      setOption("hakn", hakn)
    }
    if (!is.na(idmethod.tau)){
      method.tau <- args[[idmethod.tau]]
      ##
      imethod.tau <- charmatch(tolower(method.tau),
                               c("dl", "pm", "reml", "ml", "hs", "sj", "he", "eb"), nomatch = NA)
      ##
      if (is.na(imethod.tau) || imethod.tau==0)
        stop('Argument \'method.tau\' should be "DL", "PM", "REML", "ML", "HS", "SJ", "HE", or "EB".')
      ##
      method.tau <- c("DL", "PM", "REML", "ML", "HS", "SJ", "HE", "EB")[imethod.tau]
      ##
      setOption("method.tau", method.tau)
    }
    if (!is.na(idtau.common)){
      tau.common <- args[[idtau.common]]
      if (length(tau.common)!= 1 || !is.logical(tau.common))
        stop("Argument 'tau.common' must be a logical.")
      ##
      setOption("tau.common", tau.common)
    }
    if (!is.na(idprediction)){
      prediction <- args[[idprediction]]
      if (length(prediction)!= 1 || !is.logical(prediction))
        stop("Argument 'prediction' must be a logical.")
      ##
      setOption("prediction", prediction)
    }
    if (!is.na(idlevel.predict)){
      level.predict <- args[[idlevel.predict]]
      if (!is.numeric(level.predict) | length(level.predict)!=1)
        stop("Argument 'level.predict' must be a numeric of length 1.")
      if (level.predict <= 0 | level.predict >= 1)
        stop("Argument 'level.predict': no valid level for confidence interval")
      ##
      setOption("level.predict", level.predict)
    }
    if (!is.na(idmethod.bias)){
      method.bias <- args[[idmethod.bias]]
      ##
      imethod.bias <- charmatch(method.bias,
                         c("rank", "linreg", "mm", "count", "score", "peters"),
                         nomatch = NA)
      if(is.na(imethod.bias) | imethod.bias==0)
        stop("method.bias should be \"rank\", \"linreg\", \"mm\", \"count\", \"score\", or \"peters\"")
      ##
      method.bias <- c("rank", "linreg", "mm", "count", "score", "peters")[imethod.bias]
      ##
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
      if (length(print.byvar)!= 1 || !is.logical(print.byvar))
        stop("Argument 'print.byvar' must be a logical.")
      ##
      setOption("print.byvar", print.byvar)
    }
    if (!is.na(idkeepdata)){
      keepdata <- args[[idkeepdata]]
      if (length(keepdata)!= 1 || !is.logical(keepdata))
        stop("Argument 'keepdata' must be a logical.")
      ##
      setOption("keepdata", keepdata)
    }
    if (!is.na(idwarn)){
      warn <- args[[idwarn]]
      if (length(warn)!= 1 || !is.logical(warn))
        stop("Argument 'warn' must be a logical.")
      ##
      setOption("warn", warn)
    }
    if (!is.na(idbacktransf)){
      backtransf <- args[[idbacktransf]]
      if (length(backtransf)!= 1 || !is.logical(backtransf))
        stop("Argument 'backtransf' must be a logical.")
      ##
      setOption("backtransf", backtransf)
    }
    ##
    ## R function metabin
    ##
    if (!is.na(idsmbin)){
      smbin <- args[[idsmbin]]
      ##
      if (match(smbin, c("OR", "RD", "RR", "AS"), nomatch=0) == 0)
        stop("possible summary measures are \"OR\", \"RD\", \"RR\", and \"AS\"")
      ##
      setOption("smbin", smbin)
    }
    if (!is.na(idmethod)){
      method <- args[[idmethod]]
      ##
      imeth <- charmatch(tolower(method),
                         c("inverse", "mh", "peto"), nomatch = NA)
      ##
      if(is.na(imeth))
        stop("Argument 'method' should be \"Inverse\", \"MH\", or \"Peto\".")
      ##
      method <- c("Inverse", "MH", "Peto")[imeth]
      ##
      setOption("method", method)
    }
    if (!is.na(idincr)){
      incr <- args[[idincr]]
      ##
      if (!is.numeric(incr)){
        ##
        iincr <- charmatch(tolower(incr),
                           c("tacc"), nomatch = NA)
        ##
        if(is.na(iincr))
          stop("incr should be numeric or the character string \"TACC\"")
        ##
        incr <- c("TACC")[iincr]
      }
      ##
      setOption("incr", incr)
    }
    if (!is.na(idallincr)){
      allincr <- args[[idallincr]]
      if (length(allincr)!= 1 || !is.logical(allincr))
        stop("Argument 'allincr' must be a logical.")
      ##
      setOption("allincr", allincr)
    }
    if (!is.na(idcomb.fixed)){
      addincr <- args[[idaddincr]]
      if (length(addincr)!= 1 || !is.logical(addincr))
        stop("Argument 'addincr' must be a logical.")
      ##
      setOption("addincr", addincr)
    }
    if (!is.na(idallstudies)){
      allstudies <- args[[idallstudies]]
      if (length(allstudies)!= 1 || !is.logical(allstudies))
        stop("Argument 'allstudies' must be a logical.")
      ##
      setOption("allstudies", allstudies)
    }
    if (!is.na(idMH.exact)){
      MH.exact <- args[[idMH.exact]]
      if (length(MH.exact)!= 1 || !is.logical(MH.exact))
        stop("Argument 'MH.exact' must be a logical.")
      ##
      setOption("MH.exact", MH.exact)
    }
    if (!is.na(idRR.cochrane)){
      RR.cochrane <- args[[idRR.cochrane]]
      if (length(RR.cochrane)!= 1 || !is.logical(RR.cochrane))
        stop("Argument 'RR.cochrane' must be a logical.")
      ##
      setOption("RR.cochrane", RR.cochrane)
    }
    if (!is.na(idprint.CMH)){
      print.CMH <- args[[idprint.CMH]]
      if (length(print.CMH)!= 1 || !is.logical(print.CMH))
        stop("Argument 'print.CMH' must be a logical.")
      ##
      setOption("print.CMH", print.CMH)
    }
    ##
    ## R function metacont
    ##
    if (!is.na(idsmcont)){
      smcont <- args[[idsmcont]]
      ##
      if (smcont == "WMD"|smcont=="wmd"){
        if (warn)
          warning("Effect measure '", smcont, "' renamed as 'MD'.")
        smcont <- "MD"
      }
      ##
      if (match(smcont, c("MD", "SMD"), nomatch=0) == 0)
        stop("Possible summary measures are \"MD\" and \"SMD\".")
      ##
      setOption("smcont", smcont)
    }
    ##
    ## R function metacor
    ##
    if (!is.na(idsmcor)){
      smcor <- args[[idsmcor]]
      ##
      imeth <- charmatch(tolower(smcor), c("zcor", "cor"), nomatch = NA)
      ##
      if(is.na(imeth) || imeth==0)
        stop("Argument 'smcr' should be \"ZCOR\" or \"COR\".")
      ##
      smcor <- c("ZCOR", "COR")[imeth]
      ##
      setOption("smcor", smcor)
    }
    ##
    ## R function metainc
    ##
    if (!is.na(idsminc)){
      sminc <- args[[idsminc]]
      ##
      imeth <- charmatch(tolower(sminc),
                         c("irr", "ird"), nomatch = NA)
      ##
      if(is.na(imeth) || imeth==0)
        stop("Argument 'sminc' should be \"IRR\", \"IRD\".")
      ##
      sminc <- c("IRR", "IRD")[imeth]
      ##
      setOption("sminc", sminc)
    }
    ##
    ## R function metacont
    ##
    if (!is.na(idpooledvar)){
      pooledvar <- args[[idpooledvar]]
      if (length(pooledvar)!= 1 || !is.logical(pooledvar))
        stop("Argument 'pooledvar' must be a logical.")
      ##
      setOption("pooledvar", pooledvar)
    }
    ##
    ## R function metaprop
    ##
    if (!is.na(idsmprop)){
      smprop <- args[[idsmprop]]
      ##
      imeth <- charmatch(tolower(smprop),
                         c("pft", "pas", "praw", "pln", "plogit"), nomatch = NA)
      ##
      if(is.na(imeth) || imeth==0)
        stop("Argument 'smprop' should be \"PLOGIT\", \"PLN\", \"PFT\", \"PAS\", or \"PRAW\".")
      ##
      smprop <- c("PFT", "PAS", "PRAW", "PLN", "PLOGIT")[imeth]
      ##
      setOption("smprop", smprop)
    }
    ##
    if (!is.na(idmethod.ci)){
      method.ci <- args[[idmethod.ci]]
      ##
      imci <- charmatch(tolower(method.ci),
                        c("cp", "ws", "wscc", "ac", "sa", "sacc", "nasm"), nomatch = NA)
      ##
      if(is.na(imci) || imci==0)
        stop("Argument 'method.ci' should be \"CP\", \"WS\", \"WSCC\", \"AC\", \"SA\", \"SACC\", or \"NAsm\".")
      ##
      method.ci <- c("CP", "WS", "WSCC", "AC", "SA", "SACC", "NAsm")[imci]
      ##
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
  }
  
  invisible(res)
}
