update.meta <- function(object, 
                        data=object$data,
                        subset=object$subset,
                        studlab=object$data$.studlab,
                        method=object$method,
                        sm=object$sm,
                        incr=object$incr,
                        allincr=object$allincr,
                        addincr=object$addincr,
                        allstudies=object$allstudies,
                        MH.exact=object$MH.exact,
                        RR.cochrane=object$RR.cochrane,
                        level=object$level,
                        level.comb=object$level.comb,
                        comb.fixed=object$comb.fixed,
                        comb.random=object$comb.random,
                        hakn=object$hakn,
                        method.tau=object$method.tau,
                        tau.preset=object$tau.preset,
                        TE.tau=object$TE.tau,
                        tau.common=object$tau.common,
                        prediction=object$prediction,
                        level.predict=object$level.predict,
                        method.bias=object$method.bias,
                        ##
                        backtransf=object$backtransf,
                        title=object$title,
                        complab=object$complab,
                        outclab=object$outclab,
                        label.e=object$label.e,
                        label.c=object$label.c,
                        label.left=object$label.left,
                        label.right=object$label.right,
                        n.e=object$n.e,
                        n.c=object$n.c,
                        pooledvar=object$pooledvar,
                        method.smd=object$method.smd,
                        sd.glass=object$sd.glass,
                        exact.smd=object$exact.smd,
                        method.ci=object$method.ci,
                        byvar=object$byvar,
                        bylab=object$bylab,
                        print.byvar=object$print.byvar,
                        print.CMH=object$print.CMH,
                        keepdata=TRUE,
                        ##
                        left=object$left,
                        ma.fixed=object$ma.fixed,
                        type=object$type,
                        n.iter.max=object$n.iter.max,
                        ##
                        warn=object$warn,
                        ...){
  
  
  ##
  ##
  ## (1) Check for meta object
  ##
  ##
  chkclass(object, "meta")
  
  
  ##
  ##
  ## (2) Replace missing arguments with defaults
  ##
  ##
  replacemiss <- function(x, replacement){
    ##
    xnam <- deparse(substitute(x))
    ##
    if (is.null(x))
      if (missing(replacement))
        res <- .settings[[xnam]]
      else
        res <- replacement
    else
      res <- x
    ##
    res
  }
  ##
  comb.fixed <- replacemiss(comb.fixed)
  comb.random <- replacemiss(comb.random)
  ##
  level <- replacemiss(level)
  level.comb <- replacemiss(level.comb)
  ##
  hakn <- replacemiss(hakn)
  method.tau <- replacemiss(method.tau)
  tau.preset <- replacemiss(tau.preset)
  TE.tau <- replacemiss(TE.tau)
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
  title <- replacemiss(title)
  complab <- replacemiss(complab)
  outclab <- replacemiss(outclab, "")
  label.e <- replacemiss(label.e)
  label.c <- replacemiss(label.c)
  ##
  print.byvar <- replacemiss(print.byvar)
  ##
  warn <- replacemiss(warn)
  
  
  ##
  ##
  ## (3) Update trim-and-fill object
  ##
  ##
  if (inherits(object, "trimfill")){
    ##
    rmfilled <- function(x){
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
                    left=left, ma.fixed=ma.fixed,
                    type=type, n.iter.max=n.iter.max,
                    level=level, level.comb=level.comb,
                    comb.fixed=comb.fixed, comb.random=comb.random,
                    hakn=hakn,
                    method.tau=method.tau,
                    prediction=prediction, level.predict=level.predict,
                    silent=TRUE,
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
  if (inherits(object, "metacum") | inherits(object, "metainf")){
    ##
    res <- object
    ##
    res$comb.fixed <- ifelse(res$pooled=="fixed", TRUE, FALSE)
    res$comb.random <- ifelse(res$pooled=="random", TRUE, FALSE)
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
        as.numeric(unlist(strsplit(object$version, "-"))[1]) >= 3.2)){
    ## Some additional changes for meta objects with version < 3.2
    object$subset <- NULL
    ##
    object$data <- data.frame(.studlab=object$studlab)
    ##
    if (!is.null(object$byvar))
      object$data$.byvar <- object$byvar
    ##
    if (inherits(object,"metabin")){
      object$data$.event.e <- object$event.e
      object$data$.n.e <- object$n.e
      object$data$.event.c <- object$event.c
      object$data$.n.c <- object$n.c
    }
    if (inherits(object,"metacont")){
      object$data$.n.e <- object$n.e
      object$data$.mean.e <- object$mean.e
      object$data$.sd.e <- object$sd.e
      object$data$.n.c <- object$n.c
      object$data$.mean.c <- object$mean.c
      object$data$.sd.c <- object$sd.c
    }
    if (inherits(object,"metagen")){
      object$data$.TE <- object$TE
      object$data$.seTE <- object$seTE
    }
    if (inherits(object,"metaprop")){
      object$data$.event <- object$event
      object$data$.n <- object$n
    }
    if (inherits(object,"metacor")){
      object$data$.cor <- object$cor
      object$data$.n <- object$n
    }
  }
  ##  
  if (is.null(object$data)){
    warning("Necessary data not available. Please, recreate meta-analysis object without option 'keepdata=FALSE'.")
    return(invisible(NULL))
  }
  ##  
  missing.subset  <- missing(subset)
  missing.byvar   <- missing(byvar)
  missing.studlab <- missing(studlab)
  ##  
  mf <- match.call()
  ##
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  ##
  byvar <- eval(mf[[match("byvar", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  ##
  if (!missing.byvar){
    byvar.name <- as.character(mf[[match("byvar", names(mf))]])
    if (length(byvar.name)>1 & byvar.name[1]=="$")
      byvar.name <- byvar.name[length(byvar.name)]
    if (length(byvar.name)>1)
      byvar.name <- "byvar"
    ##
    bylab <- if (!missing(bylab) && !is.null(bylab)) bylab else byvar.name
  }
  ##
  studlab <- eval(mf[[match("studlab", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  ##
  if (missing.subset){
    if (!is.null(object$subset))
      subset <- object$subset
    else if (!is.null(object$data$.subset))
      subset <- object$data$.subset
  }
  ##
  if (missing.byvar & !is.null(object$data$.byvar))
    byvar <- object$data$.byvar
  ##
  if (missing.studlab & !is.null(object$data$.studlab))
    studlab <- object$data$.studlab
  
  
  ##
  ##
  ## (6) Update meta object
  ##
  ##
  if (inherits(object,"metabin"))
    m <- metabin(event.e=object$data$.event.e,
                 n.e=object$data$.n.e,
                 event.c=object$data$.event.c,
                 n.c=object$data$.n.c,
                 studlab=studlab,
                 ##
                 data=data, subset=subset,
                 ##
                 method=method,
                 sm=sm,
                 incr=incr, allincr=allincr, addincr=addincr, allstudies=allstudies,
                 MH.exact=MH.exact, RR.cochrane=RR.cochrane,
                 ##
                 level=level, level.comb=level.comb,
                 comb.fixed=comb.fixed, comb.random=comb.random,
                 ##
                 hakn=hakn, method.tau=method.tau,
                 tau.preset=tau.preset, TE.tau=TE.tau, tau.common=tau.common,
                 ##
                 prediction=prediction, level.predict=level.predict,
                 ##
                 method.bias=method.bias,
                 ##
                 backtransf=backtransf,
                 title=title, complab=complab, outclab=outclab,
                 label.e=label.e, label.c=label.c,
                 label.right=label.right, label.left=label.left,
                 ##
                 byvar=byvar, bylab=bylab, print.byvar=print.byvar,
                 print.CMH=print.CMH,
                 ##
                 keepdata=keepdata,
                 warn=warn)
  ##
  if (inherits(object,"metacont"))
    m <- metacont(n.e=object$data$.n.e,
                  mean.e=object$data$.mean.e,
                  sd.e=object$data$.sd.e,
                  n.c=object$data$.n.c,
                  mean.c=object$data$.mean.c,
                  sd.c=object$data$.sd.c,
                  studlab=studlab,
                  ##
                  data=data, subset=subset,
                  ##
                  sm=sm, pooledvar=pooledvar,
                  method.smd=method.smd, sd.glass=sd.glass, exact.smd=exact.smd,
                  ##
                  level=level, level.comb=level.comb,
                  comb.fixed=comb.fixed, comb.random=comb.random,
                  ##
                  hakn=hakn, method.tau=method.tau,
                  tau.preset=tau.preset, TE.tau=TE.tau, tau.common=tau.common,
                  ##
                  prediction=prediction, level.predict=level.predict,
                  ##
                  method.bias=method.bias,
                  ##
                  title=title, complab=complab, outclab=outclab,
                  label.e=label.e, label.c=label.c,
                  label.right=label.right, label.left=label.left,
                  ##
                  byvar=byvar, bylab=bylab, print.byvar=print.byvar,
                  ##
                  keepdata=keepdata,
                  warn=warn)
  ##
  if (inherits(object,"metacor"))
    m <- metacor(cor=object$data$.cor,
                 n=object$data$.n,
                 studlab=studlab,
                 ##
                 data=data, subset=subset,
                 ##
                 sm=sm,
                 ##
                 level=level, level.comb=level.comb,
                 comb.fixed=comb.fixed, comb.random=comb.random,
                 ##
                 hakn=hakn, method.tau=method.tau,
                 tau.preset=tau.preset, TE.tau=TE.tau, tau.common=tau.common,
                 ##
                 prediction=prediction, level.predict=level.predict,
                 ##
                 method.bias=method.bias,
                 ##
                 backtransf=backtransf,
                 title=title, complab=complab, outclab=outclab,
                 byvar=byvar, bylab=bylab, print.byvar=print.byvar,
                 ##
                 keepdata=keepdata)
  ##
  if (inherits(object,"metagen"))
    m <- metagen(TE=object$data$.TE,
                 seTE=object$data$.seTE,
                 studlab=studlab,
                 ##
                 data=data, subset=subset,
                 ##
                 sm=sm,
                 ##
                 level=level, level.comb=level.comb,
                 comb.fixed=comb.fixed, comb.random=comb.random,
                 ##
                 hakn=hakn, method.tau=method.tau,
                 tau.preset=tau.preset, TE.tau=TE.tau, tau.common=tau.common,
                 ##
                 prediction=prediction, level.predict=level.predict,
                 ##
                 method.bias=method.bias,
                 ##
                 n.e=n.e, n.c=n.c,
                 ##
                 backtransf=backtransf,
                 title=title, complab=complab, outclab=outclab,
                 label.e=label.e, label.c=label.c,
                 label.right=label.right, label.left=label.left,
                 ##
                 byvar=byvar, bylab=bylab, print.byvar=print.byvar,
                 ##
                 keepdata=keepdata,
                 warn=warn)
  ##
  if (inherits(object,"metainc"))
    m <- metainc(event.e=object$data$.event.e,
                 time.e=object$data$.time.e,
                 event.c=object$data$.event.c,
                 time.c=object$data$.time.c,
                 studlab=studlab,
                 ##
                 data=data, subset=subset,
                 ##
                 method=method,
                 sm=sm,
                 incr=incr, allincr=allincr, addincr=addincr,
                 ##
                 level=level, level.comb=level.comb,
                 comb.fixed=comb.fixed, comb.random=comb.random,
                 ##
                 hakn=hakn, method.tau=method.tau,
                 tau.preset=tau.preset, TE.tau=TE.tau, tau.common=tau.common,
                 ##
                 prediction=prediction, level.predict=level.predict,
                 ##
                 method.bias=method.bias,
                 ##
                 n.e=n.e, n.c=n.c,
                 ##
                 backtransf=backtransf,
                 title=title, complab=complab, outclab=outclab,
                 label.e=label.e, label.c=label.c,
                 label.right=label.right, label.left=label.left,
                 ##
                 byvar=byvar, bylab=bylab, print.byvar=print.byvar,
                 ##
                 keepdata=keepdata,
                 warn=warn)
  ##
  if (inherits(object,"metaprop"))
    m <- metaprop(event=object$data$.event,
                  n=object$data$.n,
                  studlab=studlab,
                  ##
                  data=data, subset=subset,
                  ##
                  sm=sm,
                  incr=incr, allincr=allincr, addincr=addincr,
                  method.ci=ifelse(is.null(method.ci), "CP", method.ci),
                  ##
                  level=level, level.comb=level.comb,
                  comb.fixed=comb.fixed, comb.random=comb.random,
                  ##
                  hakn=hakn, method.tau=method.tau,
                  tau.preset=tau.preset, TE.tau=TE.tau, tau.common=tau.common,
                  ##
                  prediction=prediction, level.predict=level.predict,
                  ##
                  method.bias=method.bias,
                  ##
                  backtransf=backtransf,
                  title=title, complab=complab, outclab=outclab,
                  byvar=byvar, bylab=bylab, print.byvar=print.byvar,
                  ##
                  keepdata=keepdata,
                  warn=warn)
  ##  
  m$call.object <- object$call
  m$call <- match.call()
  
  
  m
}
