update.meta <- function(object, 
                        data=object$data,
                        subset=object$subset,
                        studlab=object$data$studlab,
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
                        title=object$title,
                        complab=object$complab,
                        outclab=object$outclab,
                        label.e=object$label.e,
                        label.c=object$label.c,
                        label.left=object$label.left,
                        label.right=object$label.right,
                        n.e=object$n.e,
                        n.c=object$n.c,
                        byvar=object$byvar,
                        bylab=object$bylab,
                        print.byvar=object$print.byvar,
                        print.CMH=object$print.CMH,
                        keepdata=TRUE,
                        warn=object$warn,
                        ...){
  
  
  if (!inherits(object,"meta"))
    stop("Argument 'object' must be an object of class \"meta\"")
  
  if (is.null(object$version) ||
      as.numeric(strsplit(object$version, "-")[[1]][1]) < 3.0){
    warning("Function not applicable to objects from older versions of meta (< 3.0).")
    return(invisible(NULL))
  }
  
  if (is.null(object$data)){
    warning("Necessary data not available. Please, recreate meta-analysis object without option 'keepdata=FALSE'.")
    return(invisible(NULL))
  }
  
  ## Extract byvar from dataset in meta object
  ##
  if (!missing(byvar)){
    byname <- deparse(substitute(byvar))
    byname <- gsub("\"", "", byname)
    if (length(grep("\\$", byname))==0)
      if (match(byname, names(data), nomatch=-1)!=-1)
        byvar <- data[,byname]
    if (is.null(byvar))
      stop("Arguments 'byvar' contains no sensible information (value: NULL)")
    if (missing(bylab))
      bylab <- byname
  }
  
  
  if (inherits(object,"metabin"))
    m <- metabin(event.e=object$data$event.e,
                 n.e=object$data$n.e,
                 event.c=object$data$event.c,
                 n.c=object$data$n.c,
                 studlab=studlab,
                 data=data, subset=subset,
                 method=method,
                 sm=sm,
                 incr=incr, allincr=allincr, addincr=addincr, allstudies=allstudies,
                 MH.exact=MH.exact, RR.cochrane=RR.cochrane,
                 level=level, level.comb=level.comb,
                 comb.fixed=comb.fixed, comb.random=comb.random,
                 hakn=hakn, method.tau=method.tau,
                 tau.preset=tau.preset, TE.tau=TE.tau, tau.common=tau.common,
                 prediction=prediction, level.predict=level.predict,
                 method.bias=method.bias,
                 title=title, complab=complab, outclab=outclab,
                 label.e=label.e, label.c=label.c,
                 label.right=label.right, label.left=label.left,
                 byvar=byvar, bylab=bylab, print.byvar=print.byvar,
                 print.CMH=print.CMH,
                 keepdata=keepdata,
                 warn=warn)
  ##
  if (inherits(object,"metacont"))
    m <- metacont(n.e=object$data$n.e,
                  mean.e=object$data$mean.e,
                  sd.e=object$data$sd.e,
                  n.c=object$data$n.c,
                  mean.c=object$data$mean.c,
                  sd.c=object$data$sd.c,
                  studlab=studlab,
                  data=data, subset=subset,
                  sm=sm,
                  level=level, level.comb=level.comb,
                  comb.fixed=comb.fixed, comb.random=comb.random,
                  hakn=hakn, method.tau=method.tau,
                  tau.preset=tau.preset, TE.tau=TE.tau, tau.common=tau.common,
                  prediction=prediction, level.predict=level.predict,
                  method.bias=method.bias,
                  title=title, complab=complab, outclab=outclab,
                  label.e=label.e, label.c=label.c,
                  label.right=label.right, label.left=label.left,
                  byvar=byvar, bylab=bylab, print.byvar=print.byvar,
                  keepdata=keepdata,
                  warn=warn)
  ##
  if (inherits(object,"metagen"))
    m <- metagen(TE=object$data$TE,
                 seTE=object$data$seTE,
                 studlab=studlab,
                 data=data, subset=subset,
                 sm=sm,
                 level=level, level.comb=level.comb,
                 comb.fixed=comb.fixed, comb.random=comb.random,
                 hakn=hakn, method.tau=method.tau,
                 tau.preset=tau.preset, TE.tau=TE.tau, tau.common=tau.common,
                 prediction=prediction, level.predict=level.predict,
                 method.bias=method.bias,
                 n.e=n.e, n.c=n.c,
                 title=title, complab=complab, outclab=outclab,
                 label.e=label.e, label.c=label.c,
                 label.right=label.right, label.left=label.left,
                 byvar=byvar, bylab=bylab, print.byvar=print.byvar,
                 keepdata=keepdata,
                 warn=warn)
  ##
  if (inherits(object,"metaprop"))
    m <- metaprop(event=object$data$event,
                  n=object$data$n,
                  studlab=studlab,
                  data=data, subset=subset,
                  sm=sm,
                  incr=incr, allincr=allincr, addincr=addincr,
                  level=level, level.comb=level.comb,
                  comb.fixed=comb.fixed, comb.random=comb.random,
                  hakn=hakn, method.tau=method.tau,
                  tau.preset=tau.preset, TE.tau=TE.tau, tau.common=tau.common,
                  prediction=prediction, level.predict=level.predict,
                  method.bias=method.bias,
                  title=title, complab=complab, outclab=outclab,
                  byvar=byvar, bylab=bylab, print.byvar=print.byvar,
                  keepdata=keepdata,
                  warn=warn)
  ##
  if (inherits(object,"metacor"))
    m <- metacor(cor=object$data$cor,
                 n=object$data$n,
                 studlab=studlab,
                 data=data, subset=subset,
                 sm=sm,
                 level=level, level.comb=level.comb,
                 comb.fixed=comb.fixed, comb.random=comb.random,
                 hakn=hakn, method.tau=method.tau,
                 tau.preset=tau.preset, TE.tau=TE.tau, tau.common=tau.common,
                 prediction=prediction, level.predict=level.predict,
                 method.bias=method.bias,
                 title=title, complab=complab, outclab=outclab,
                 byvar=byvar, bylab=bylab, print.byvar=print.byvar,
                 keepdata=keepdata)
  
  
  m
}

