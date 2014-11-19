forest.meta <- function(x,
                        sortvar,
                        studlab=TRUE,
                        ##
                        comb.fixed=x$comb.fixed,
                        comb.random=x$comb.random,
                        overall=TRUE,
                        text.fixed=if (x$level!=x$level.comb) paste("Fixed effect model (",
                                         round(x$level.comb*100), "%-CI)", sep="") else "Fixed effect model",
                        text.random=if (x$level!=x$level.comb) paste("Random effects model (",
                                          round(x$level.comb*100), "%-CI)", sep="") else "Random effects model",
                        lty.fixed=2, lty.random=3,
                        ##
                        prediction=x$prediction,
                        text.predict=if (!(length(x$level.predict)==0) && x$level!=x$level.predict) paste("Prediction interval (",
                                           round(x$level.predict*100), "%)", sep="") else "Prediction interval",
                        ##
                        bylab=x$bylab,
                        print.byvar=x$print.byvar,
                        text.fixed.w=text.fixed,
                        text.random.w=text.random,
                        bysort=FALSE,
                        ##
                        pooled.totals=comb.fixed|comb.random, pooled.events=FALSE,
                        ##
                        xlab="", xlab.pos=ref,
                        smlab=NULL, smlab.pos=ref, xlim="symmetric",
                        ##
                        allstudies=TRUE,
                        weight,
                        pscale=1,
                        ##
                        ref=ifelse(backtransf & is.relative.effect(x$sm), 1, 0),
                        ##
                        leftcols=NULL, rightcols=NULL,
                        leftlabs=NULL, rightlabs=NULL,
                        ##
                        lab.e=x$label.e,
                        lab.c=x$label.c,
                        ##
                        lab.e.attach.to.col=NULL,
                        lab.c.attach.to.col=NULL,
                        ##
                        label.right=x$label.right,
                        label.left=x$label.left,
                        ##
                        lab.NA=".",
                        ##
                        lwd=1,
                        ##
                        at=NULL,
                        label=TRUE,
                        ##
                        col.i="black",
                        col.i.inside.square="white",
                        col.square="gray",
                        col.square.lines=col.square,
                        ##
                        col.diamond="gray",
                        col.diamond.fixed=col.diamond,
                        col.diamond.random=col.diamond,
                        col.diamond.lines="black",
                        col.diamond.fixed.lines=col.diamond.lines,
                        col.diamond.random.lines=col.diamond.lines,
                        ##
                        col.predict="red",
                        col.predict.lines="black",
                        ##
                        col.by="darkgray",
                        ##
                        print.I2=comb.fixed|comb.random,
                        print.tau2=comb.fixed|comb.random,
                        print.Q=FALSE,
                        print.pval.Q=comb.fixed|comb.random,
                        hetstat=print.I2|print.tau2|print.Q|print.pval.Q,
                        overall.hetstat=overall&hetstat,
                        hetlab="Heterogeneity: ",
                        text.I2="I-squared",
                        text.tau2="tau-squared",
                        ##
                        fontsize=12,
                        fs.heading=fontsize,
                        fs.fixed=fontsize,
                        fs.random=fs.fixed,
                        fs.predict=fs.fixed,
                        fs.study=fontsize,
                        fs.fixed.labels=fs.fixed,
                        fs.random.labels=fs.random,
                        fs.predict.labels=fs.predict,
                        fs.study.labels=fs.study,
                        fs.hetstat=fontsize-2,
                        fs.axis=fontsize,
                        fs.smlab=fontsize,
                        fs.xlab=fontsize,
                        fs.lr=fontsize,
                        ##
                        ff.heading="bold",
                        ff.fixed="bold",
                        ff.random=ff.fixed,
                        ff.predict=ff.fixed,
                        ff.study="plain",
                        ff.fixed.labels=ff.fixed,
                        ff.random.labels=ff.random,
                        ff.predict.labels=ff.predict,
                        ff.study.labels=ff.study,
                        ff.hetstat="bold.italic",
                        ff.axis="plain",
                        ff.smlab="bold",
                        ff.xlab="plain",
                        ff.lr="plain",
                        ##
                        squaresize=0.8,
                        ##
                        plotwidth=grid::unit(6, "cm"),
                        colgap=grid::unit(2, "mm"),
                        colgap.left=colgap,
                        colgap.right=colgap,
                        colgap.forest=colgap,
                        colgap.forest.left=colgap.forest,
                        colgap.forest.right=colgap.forest,
                        ##
                        just="right",
                        just.studlab="left",
                        just.addcols="center",
                        ##
                        addspace=TRUE,
                        ##
                        new=TRUE,
                        ##
                        backtransf=x$backtransf,
                        digits=2, ...){
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(x, "meta")
  x <- updateversion(x)
  ##
  K.all <- length(x$TE)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  mf <- match.call()
  error <- try(sortvar <- eval(mf[[match("sortvar", names(mf))]],
                               as.data.frame(x, stringsAsFactors=FALSE),
                               enclos = sys.frame(sys.parent())),
               silent=TRUE)
  if (class(error)=="try-error"){
    xd <- x$data
    sortvar <- eval(mf[[match("sortvar", names(mf))]],
                    xd, enclos = NULL)
    if (!is.null(x$data$.subset))
      sortvar <- sortvar[x$data$.subset]
  }
  sort <- !is.null(sortvar)
  if (sort && (length(sortvar) != K.all))
    stop("Number of studies in object 'x' and argument 'sortvar' have different length.")
  if (!sort)
    sortvar <- 1:K.all
  ##
  slab <- TRUE
  if (length(studlab) == 1 & is.logical(studlab))
    if (studlab == FALSE){
      studlab <- rep("", K.all)
      slab <- FALSE
    }
    else studlab <- x$studlab
  ##
  if (length(studlab) != K.all)
    stop("Number of studies in object 'x' and argument 'studlab' have different length.")
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(overall)
  chknumeric(lty.fixed)
  chknumeric(lty.random)
  chklogical(prediction)
  if (!is.null(print.byvar))
    chklogical(print.byvar)
  chklogical(bysort)
  chklogical(pooled.totals)
  chklogical(pooled.events)
  ## chknumeric(xlab.pos) ??
  ## chknumeric(smlab.pos) ??
  chklogical(allstudies)
  chknumeric(pscale, single=TRUE)
  chknumeric(ref)
  if (!is.null(at))
    chknumeric(at)
  chkchar(col.diamond)
  chkchar(col.diamond.fixed)
  chkchar(col.diamond.random)
  chkchar(col.diamond.lines)
  chkchar(col.diamond.fixed.lines)
  chkchar(col.diamond.random.lines)
  chkchar(col.predict)
  chkchar(col.predict.lines)
  chklogical(print.I2)
  chklogical(print.tau2)
  chklogical(print.Q)
  chklogical(print.pval.Q)
  chklogical(hetstat)
  chklogical(overall.hetstat)
  chknumeric(fontsize, single=TRUE)
  chknumeric(fs.heading, single=TRUE)
  chknumeric(fs.fixed, single=TRUE)
  chknumeric(fs.random, single=TRUE)
  chknumeric(fs.predict, single=TRUE)
  chknumeric(fs.study, single=TRUE)
  chknumeric(fs.fixed.labels, single=TRUE)
  chknumeric(fs.random.labels, single=TRUE)
  chknumeric(fs.predict.labels, single=TRUE)
  chknumeric(fs.study.labels, single=TRUE)
  chknumeric(fs.hetstat, single=TRUE)
  chknumeric(fs.axis, single=TRUE)
  chknumeric(fs.smlab, single=TRUE)
  chknumeric(fs.xlab, single=TRUE)
  chknumeric(fs.lr, single=TRUE)
  chknumeric(squaresize, single=TRUE)
  just.cols <- setchar(just, c("right", "center", "left"))
  just.studlab <- setchar(just.studlab, c("right", "center", "left"))
  just.addcols <- setchar(just.addcols, c("right", "center", "left"))
  ##
  if (missing(weight))
    weight <- ifelse(comb.random & !comb.fixed, "random", "fixed")
  weight <- setchar(weight, c("same", "fixed", "random"))
  ##
  cl <- class(x)[1]
  addargs <- names(list(...))
  ##
  fun <- "forest.meta"
  ##
  warnarg("byvar", addargs, fun, cl)
  warnarg("level", addargs, fun, cl)
  warnarg("level.comb", addargs, fun, cl)
  warnarg("level.predict", addargs, fun, cl)
  
  
  ##
  ##
  ## (3) Check length of variables
  ##
  ##
  if (length(col.i)==1)
    col.i <- rep(col.i, K.all)
  else
    chklength(col.i, K.all, fun)
  ##
  if (length(col.i.inside.square)==1)
    col.i.inside.square <- rep(col.i.inside.square, K.all)
  else
    chklength(col.i.inside.square, K.all, fun)
  ##
  if (length(col.square)==1)
    col.square <- rep(col.square, K.all)
  else
    chklength(col.square, K.all, fun)
  ##
  if (length(col.square.lines)==1)
    col.square.lines <- rep(col.square.lines, K.all)
  else
    chklength(col.square.lines, K.all, fun)
  
  
  ##
  ##
  ## (4) Some assignments and additional checks
  ##
  ##
  metaprop <- inherits(x, "metaprop")
  metacor <- inherits(x, "metacor")
  metainf.metacum <- inherits(x, "metainf") | inherits(x, "metacum")
  ##  
  if (metainf.metacum){
    overall.hetstat <- FALSE
    hetstat <- FALSE
    prediction <- FALSE
  }
  ##
  prediction <- prediction & comb.random & x$k >= 3
  ##  
  x.name <- deparse(substitute(x))
  ##
  byvar <- x$byvar
  level <- x$level
  level.comb <- x$level.comb
  level.predict <- x$level.predict
  ##
  chklevel(level)
  chklevel(level.comb)
  if (prediction)
    chklevel(level.predict)
  ##  
  if (is.logical(leftcols)){
    warning("Logical value not possible for argument 'leftcols', set to 'NULL'.")
    leftcols <- NULL
  }
  ##
  notmiss.xlim <- !missing(xlim)
  ##
  if (just.studlab=="left")
    xpos.studlab <- 0
  else if (just.studlab=="center")
      xpos.studlab <- 0.5
  else if (just.studlab=="right")
    xpos.studlab <- 1
  ##
  sm <- x$sm
  ##
  log.xaxis <- FALSE
  ##
  if (sm %in% c("PLOGIT", "PLN", "PRAW", "PAS", "PFT"))
    ref <- NA
  ##
  if (backtransf & is.relative.effect(sm)){
    ref <- log(ref)
    log.xaxis <- TRUE
  }
  ##  
  if (!backtransf & pscale!=1){
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  ##  
  if (is.null(xlab))
    xlab <- xlab(sm, backtransf)
  ##
  if (is.null(smlab))
    smlab <- xlab(sm, backtransf)
  ##
  if (is.null(label.right))
    label.right <- ""
  if (is.null(label.left))
    label.left <- ""
  ##
  by <- !is.null(byvar)
  
  
  ##
  ##
  ## (5) Determine columns on left and right side of forest plot
  ##
  ##
  ## Determine whether to print columns on right and/or left side
  ## of forest plot
  ##
  rsel <- !(is.logical(rightcols) && length(rightcols)==1 && !rightcols)
  ##
  if (!rsel)
    rightcols <- NULL
  ##
  ## Check for duplicate columns
  ##
  if (length(c(rightcols, leftcols))>0 &&
      any(duplicated(c(rightcols, leftcols))))
    stop("Duplicate entries in 'leftcols' and 'rightcols'.")
  ##
  ## Predefined columns and labels
  ##
  sm.lab <- sm
  ##
  if (backtransf){
    if (sm=="ZCOR")
      sm.lab <- "COR"
    else if (sm %in% c("PLOGIT", "PLN", "PRAW", "PAS", "PFT")){
      if (pscale==100)
        sm.lab <- "Prop (in %)"
      else if (pscale==1)
        sm.lab <- "Proportion"
      else
        sm.lab <- sm
    }
    else if (sm=="proportion")
      sm.lab <- "Proportion"
  }
  else 
    if (is.relative.effect(sm))
      sm.lab <- paste("log", sm, sep="")
  ##
  colnames <- c("studlab", "TE", "seTE",
                "n.e", "n.c",
                "event.e", "event.c",
                "mean.e", "mean.c",
                "sd.e", "sd.c",
                "cor",
                "time.e", "time.c",
                "effect", "ci",
                "w.fixed", "w.random")
  ##
  labnames <- c("Study", "TE", "seTE",
                "Total", "Total", "Events", "Events",
                "Mean", "Mean", "SD", "SD",
                "Cor",
                "Time", "Time",
                sm.lab, paste(100*level, "%-CI", sep=""),
                "W(fixed)", "W(random)")
  ##
  ## Identify and process columns in addition to columns
  ## defined above in variables 'colnames' and 'labnames'
  ##
  colnames.new <- c(rightcols, leftcols)[!c(rightcols, leftcols) %in% colnames]
  ##
  newcols <- length(colnames.new)>0
  ##
  if (newcols){
    if (is.null(x$data))
      dataset <- as.data.frame(x)
    else{
      if (!is.null(x$subset))
        dataset <- x$data[x$subset,]
      else
        dataset <- x$data
    }
    ##
    ## Check whether additional variables are
    ## part of meta-object
    ##
    for (i in colnames.new)
      if (length(x[[i]]) == 0 & length(dataset[[i]]) == 0)
        stop("Variable '", i, "' not available in '", x.name, "'.")
    ##
    rightcols.new <- rightcols[! rightcols %in% colnames]
    leftcols.new  <- leftcols[! leftcols %in% colnames]
    ##
    ## Determine label for new columns
    ## 1. Use column name as label if no label is given
    ##    argument right|left|labs
    ## 2. Otherwise use corresponding entry from
    ##    argument right|left|labs
    ##
    if (length(rightcols.new)>0){
      pos.rightcols.new <- match(rightcols.new, rightcols)
      ##
      if (missing(rightlabs))
        rightlabs.new <- rightcols.new
      else if (length(rightcols.new) == length(rightlabs))
        rightlabs.new <- rightlabs
      else if (max(pos.rightcols.new) <= length(rightlabs))
        rightlabs.new <- rightlabs[pos.rightcols.new]
      else if (max(pos.rightcols.new) > length(rightlabs))
        stop("Too few labels defined for argument 'rightcols'.")
    }
    if (length(leftcols.new)>0){
      pos.leftcols.new <- match(leftcols.new, leftcols)
      ##
      if (missing(leftlabs))
        leftlabs.new <- leftcols.new
      else if (length(leftcols.new) == length(leftlabs))
        leftlabs.new <- leftlabs
      else if (max(pos.leftcols.new) <= length(leftlabs))
        leftlabs.new <- leftlabs[pos.leftcols.new]
      else if (max(pos.leftcols.new) > length(leftlabs))
        stop("Too few labels defined for argument 'leftcols'.")
    }
  }
  ##
  ## Default set of columns if argument leftcols and/or
  ## rightcols not specified
  ##
  if (is.null(leftcols)){
    if (inherits(x, "metabin"))
      leftcols <- c("studlab",
                    "event.e", "n.e",
                    "event.c", "n.c")
    ##
    if (inherits(x, "metacont"))
      leftcols <- c("studlab",
                    "n.e", "mean.e", "sd.e",
                    "n.c", "mean.c", "sd.c")
    ##
    if (inherits(x, "metagen"))
      leftcols <- c("studlab",
                    "TE", "seTE")
    ##
    if (metainf.metacum)
      leftcols <- "studlab"
    ##
    if (metaprop)
      leftcols <- c("studlab",
                    "event.e", "n.e")
    ##
    if (metacor)
      leftcols <- c("studlab",
                    "n.e")
    ##
    if (inherits(x, "metainc"))
      leftcols <- c("studlab",
                    "event.e", "time.e",
                    "event.c", "time.c")
  }
  ##
  if (is.null(rightcols) & rsel){
    rightcols <- c("effect", "ci")
    ##
    if (!metainf.metacum){
      if (comb.fixed & overall)
        rightcols <- c(rightcols, "w.fixed")
      if (comb.random & overall)
        rightcols <- c(rightcols, "w.random")
    }
  }
  
  
  ##
  ##
  ## (6) Select data for forest plot
  ##
  ##
  if (metaprop){
    x$event.e <- x$event
    x$n.e <- x$n
    ##
    if (!is.null(rightcols)){
      if (any(rightcols=="n"))
        rightcols[rightcols=="n"] <- "n.e"
      if (any(rightcols=="event"))
        rightcols[rightcols=="event"] <- "event.e"
    }
    ##
    if (!is.null(leftcols)){
      if (any(leftcols=="n"))
        leftcols[leftcols=="n"] <- "n.e"
      if (any(leftcols=="event"))
        leftcols[leftcols=="event"] <- "event.e"
    }
  }
  ##
  if (metacor){
    x$n.e <- x$n
  }
  ##  
  if (metainf.metacum){
    ##
    x$TE.fixed    <- rev(x$TE)[1]
    x$seTE.fixed  <- rev(x$seTE)[1]
    x$lower.fixed <- rev(x$lower)[1]
    x$upper.fixed <- rev(x$upper)[1]
    ##
    x$TE.random   <- rev(x$TE)[1]
    x$seTE.random <- rev(x$seTE)[1]
    x$lower.random <- rev(x$lower)[1]
    x$upper.random <- rev(x$upper)[1]
    ##
    x$n.harmonic.mean.ma <- rev(x$n.harmonic.mean)[1]
    ##
    x$w.all <- rev(x$w)[1]
    ##
    x$TE <- rev(rev(x$TE)[-(1:2)])
    x$seTE <- rev(rev(x$seTE)[-(1:2)])
    x$studlab <- rev(rev(x$studlab)[-(1:2)])
    ##
    x$lower <- rev(rev(x$lower)[-(1:2)])
    x$upper <- rev(rev(x$upper)[-(1:2)])
    ##
    x$w.fixed <- rev(rev(x$w)[-(1:2)])
    x$w.random <- rev(rev(x$w)[-(1:2)])
    ##
    x$n.harmonic.mean <- rev(rev(x$n.harmonic.mean)[-(1:2)])
    ##
    if (overall & x$pooled=="fixed"){
      comb.fixed <- TRUE
      comb.random <- FALSE
      if (weight != "same")
        weight <- "fixed"
    }
    else if (overall & x$pooled=="random"){
      comb.fixed <- FALSE
      comb.random <- TRUE
      if (weight != "same")
        weight <- "random"
    }
  }
  ## Total number of studies to plot (*not* number of studies combined)
  k.all <- length(x$TE)
  ##
  if (allstudies)
    n.stud <- k.all # all studies
  else
    n.stud <- x$k   # number of studies combined in meta-analysis
  ##
  if (!by)
    byvar <- rep(1, k.all)
  ##
  if (by & any(is.na(byvar)))
    stop("Missing values in 'byvar'")
  if (allstudies)
    sel <- 1:length(x$TE)
  else
    sel <- !is.na(x$TE)
  ##
  if (n.stud != sum(sel>0))
    warning("n.stud != sum(sel)")
  ##
  x$n.e <- x$n.e[sel]
  x$n.c <- x$n.c[sel]
  ##
  x$event.e <- x$event.e[sel]
  x$event.c <- x$event.c[sel]
  ##
  x$mean.e <- x$mean.e[sel]
  x$mean.c <- x$mean.c[sel]
  ##
  x$sd.e <- x$sd.e[sel]
  x$sd.c <- x$sd.c[sel]
  ##
  x$cor <- x$cor[sel]
  ##
  x$time.e <- x$time.e[sel]
  x$time.c <- x$time.c[sel]
  ##
  x$TE   <- x$TE[sel]
  x$seTE <- x$seTE[sel]
  ##
  x$lower <- x$lower[sel]
  x$upper <- x$upper[sel]
  ##
  x$w.fixed  <- x$w.fixed[sel]
  x$w.random <- x$w.random[sel]
  studlab  <- studlab[sel]
  ##
  x$n.harmonic.mean <- x$n.harmonic.mean[sel]
  ##
  byvar   <- byvar[sel]
  sortvar <- sortvar[sel]
  ##
  col.i <- col.i[sel]
  col.square <- col.square[sel]
  col.square.lines <- col.square.lines[sel]
  ##
  col.i.inside.square <- col.i.inside.square[sel]
  ##  
  if (sort | by){
    if (bysort)
      bylevs <- sort(x$bylevs)
    else
      bylevs <- x$bylevs
    ##
    byvar.factor <- factor(byvar, levels=bylevs)
    o <- order(byvar.factor, sortvar)
    ##
    x$n.e <- x$n.e[o]
    x$n.c <- x$n.c[o]
    ##
    x$event.e <- x$event.e[o]
    x$event.c <- x$event.c[o]
    ##
    x$mean.e <- x$mean.e[o]
    x$mean.c <- x$mean.c[o]
    ##
    x$sd.e <- x$sd.e[o]
    x$sd.c <- x$sd.c[o]
    ##
    x$cor <- x$cor[o]
    ##
    x$time.e <- x$time.e[o]
    x$time.c <- x$time.c[o]
    ##
    x$TE   <- x$TE[o]
    x$seTE <- x$seTE[o]
    ##
    x$lower <- x$lower[o]
    x$upper <- x$upper[o]
    ##
    x$w.fixed  <- x$w.fixed[o]
    x$w.random <- x$w.random[o]
    studlab  <- studlab[o]
    ##
    x$n.harmonic.mean <- x$n.harmonic.mean[o]
    ##
    byvar   <- byvar[o]
    sortvar <- sortvar[o]
    ##
    col.i <- col.i[o]
    col.square <- col.square[o]
    col.square.lines <- col.square.lines[o]
    ##
    col.i.inside.square <- col.i.inside.square[o]
    ##
    if (newcols)
      dataset <- dataset[o,]
  }
  ##
  if (by)
    n.by <- length(bylevs)
  else
    n.by <- 0
  ##
  if (metainf.metacum){
    TE    <- x$TE
    seTE  <- x$seTE
    lowTE <- x$lower
    uppTE <- x$upper
    ##
    TE.fixed    <- x$TE.fixed
    lowTE.fixed <- x$lower.fixed
    uppTE.fixed <- x$upper.fixed
    ##
    TE.random    <- x$TE.random
    lowTE.random <- x$lower.random
    uppTE.random <- x$upper.random
    ##
    lowTE.predict <- NA
    uppTE.predict <- NA
    ##
    Q    <- NA
    df   <- NA
    I2   <- NA
    tau2 <- NA
  }
  else{
    TE <- x$TE
    seTE <- x$seTE
    lowTE <- x$lower
    uppTE <- x$upper
    ##
    if (metaprop & !backtransf){
      ciTE <- ci(TE, seTE, level=level)
      lowTE <- ciTE$lower
      uppTE <- ciTE$upper
    }
    ##
    TE.fixed <- x$TE.fixed
    lowTE.fixed <- x$lower.fixed
    uppTE.fixed <- x$upper.fixed
    ##
    TE.random <- x$TE.random
    lowTE.random <- x$lower.random
    uppTE.random <- x$upper.random
    ##
    lowTE.predict <- x$lower.predict
    uppTE.predict <- x$upper.predict
    ##
    Q    <- x$Q
    df   <- x$df.Q
    tau2 <- x$tau^2
    ##
    I2 <- x$I2
  }
  ##
  if (overall.hetstat){
    dummy <- FALSE
    ##
    hetstat.overall <- hetlab
    ##
    if (print.I2){
      hetstat.overall <- paste(hetstat.overall,
                               text.I2, "=",
                               round(100*I2, 1), "%",
                               sep="")
      dummy <- TRUE
    }
    ##
    if (print.tau2){
      hetstat.overall <- paste(hetstat.overall,
                               if (dummy) ", ",
                               if (tau2==0) paste(text.tau2, "=0", sep="")
                               else format.tau(tau2, noblanks=TRUE,
                                               lab=TRUE, labval=text.tau2),
                               sep="")
      dummy <- TRUE
    }
    ##
    if (print.Q){
      hetstat.overall <- paste(hetstat.overall,
                               if (dummy) ", ",
                               "Q=", round(Q, 1),
                               ", df=", df,
                               sep="")
      dummy <- TRUE
    }
    ##
    if (print.pval.Q){
      hetstat.overall <- paste(hetstat.overall,
                               if (dummy) ", ",
                               format.p(1-pchisq(Q, df),
                                        lab=TRUE, noblanks=TRUE),
                               sep="")
    }
  }
  else
    hetstat.overall <- ""
  
  
  ##
  ##
  ## (7) Prepare data for subgroup analysis
  ##
  ##
  if (by){
    o <- order(factor(x$bylevs, levels=bylevs))
    k.w <- x$k.w[o]
    TE.fixed.w <- x$TE.fixed.w[o]
    lower.fixed.w <- x$lower.fixed.w[o]
    upper.fixed.w <- x$upper.fixed.w[o]
    TE.random.w <- x$TE.random.w[o]
    lower.random.w <- x$lower.random.w[o]
    upper.random.w <- x$upper.random.w[o]
    Q.w        <- x$Q.w[o]
    I2.w       <- x$I2.w[o]
    tau.w      <- x$tau.w[o]
    w.fixed.w  <- x$w.fixed.w[o]
    w.random.w <- x$w.random.w[o]
    e.e.w <- if (metaprop) x$event.w[o] else x$event.e.w[o]
    n.e.w <- if (metacor | metaprop) x$n.w[o] else x$n.e.w[o]
    e.c.w <- x$event.c.w[o]
    n.c.w <- x$n.c.w[o]
    n.harmonic.mean.w <- x$n.harmonic.mean.w[o]
    k.all.w <- x$k.all.w[o]
    ##
    if (allstudies)
      sel <- 1:length(k.w)
    else
      sel <- k.w>0
    k.w <- k.w[sel]
    TE.fixed.w <- TE.fixed.w[sel]
    lower.fixed.w <- lower.fixed.w[sel]
    upper.fixed.w <- upper.fixed.w[sel]
    TE.random.w <- TE.random.w[sel]
    lower.random.w <- lower.random.w[sel]
    upper.random.w <- upper.random.w[sel]
    Q.w        <- Q.w[sel]
    I2.w       <- I2.w[sel]
    tau.w      <- tau.w[sel]
    w.fixed.w  <- w.fixed.w[sel]
    w.random.w <- w.random.w[sel]
    e.e.w <- e.e.w[sel]
    n.e.w <- n.e.w[sel]
    e.c.w <- e.c.w[sel]
    n.c.w <- n.c.w[sel]
    n.harmonic.mean.w <- n.harmonic.mean.w[sel]
    k.all.w <- k.all.w[sel]
    bylevs <- bylevs[sel]
    ##
    if (comb.fixed){
      if (sum(w.fixed.w)>0)
        w.fixed.w.p <- 100*round(w.fixed.w/sum(w.fixed.w, na.rm=TRUE), 3)
      else
        w.fixed.w.p <- w.fixed.w
    }
    else{
      TE.fixed.w <- lower.fixed.w <- upper.fixed.w <- rep(NA, n.by)
      ##
      w.fixed.w.p <- rep(NA, n.by)
      if (missing(text.fixed.w))
        text.fixed.w <- rep("Overall", n.by)
    }
    ##
    if (comb.random){
      if (sum(w.random.w)>0)
        w.random.w.p <- 100*round(w.random.w/sum(w.random.w, na.rm=TRUE), 3)
      else
        w.random.w.p <- w.random.w
    }
    else{
      TE.random.w <- lower.random.w <- upper.random.w <- rep(NA, n.by)
      ##
      w.random.w.p <- rep(NA, n.by)
      text.random.w <- rep("", n.by)
    }
    ##
    if (hetstat){
      dummy <- FALSE
      ##
      hetstat.w <- hetlab
      ##
      if (print.I2){
        hetstat.w <- paste(hetstat.w,
                           text.I2, "=",
                           round(100*I2.w, 1), "%",
                           sep="")
        dummy <- TRUE
      }
      ##
      if (print.tau2){
        hetstat.w <- paste(hetstat.w,
                           if (dummy) ", ",
                           ifelse(tau.w==0, paste(text.tau2, "=0", sep=""),
                                  format.tau(tau.w^2, noblanks=TRUE,
                                             lab=TRUE, labval=text.tau2)),
                           sep="")
        dummy <- TRUE
      }
      ##
      if (print.Q){
        hetstat.w <- paste(hetstat.w,
                           if (dummy) ", ",
                           "Q=", round(Q.w, 1),
                           ", df=", k.w-1,
                           sep="")
        dummy <- TRUE
      }
      ##
      if (print.pval.Q){
        hetstat.w <- paste(hetstat.w,
                           if (dummy) ", ",
                           format.p(1-pchisq(Q.w, k.w-1),
                                    lab=TRUE, noblanks=TRUE),
                           sep="")
      }
      ##
      sel <- k.w==1
      hetstat.w[sel] <- paste(hetlab, "not applicable for a single study", sep="")
    }
    else
      hetstat.w <- rep("", n.by)
    ##
    TE.w <- c(TE.fixed.w, TE.random.w, rep(NA, n.by))
    lowTE.w <- c(lower.fixed.w, lower.random.w, rep(NA, n.by))
    uppTE.w <- c(upper.fixed.w, upper.random.w, rep(NA, n.by))
    n.harmonic.mean.w <- c(n.harmonic.mean.w, n.harmonic.mean.w, rep(NA, n.by))
    weight.w.p <- c(w.fixed.w.p, w.random.w.p, rep(NA, n.by))
  }
  
  
  ##
  ##
  ## (8) Backtransform data
  ##
  ##
  TE.orig <- TE
  ##
  if (backtransf){
    ##
    ## Freeman-Tukey Arcsin transformation (sm="PFT")
    ##
    if (metainf.metacum){
      npft <- x$n.harmonic.mean
      npft.ma <- x$n.harmonic.mean.ma
    }
    else{
      npft <- x$n
      npft.ma <- 1/mean(1/x$n)
    }
    ##
    ## Individual study results
    ##
    if (metaprop){
      TE <- x$event/x$n
    }
    ## Relative effect measures will be back transformed later
    else if (!is.relative.effect(sm)){
      TE <- backtransf(TE, sm, "mean", npft)
      lowTE <- backtransf(lowTE, sm, "lower", npft)
      uppTE <- backtransf(uppTE, sm, "upper", npft)
    }
    ##
    ## Results of meta-analysis
    ##
    if (!is.relative.effect(sm)){
      TE.fixed    <- backtransf(TE.fixed, sm, "mean",
                                npft.ma, warn=comb.fixed)
      lowTE.fixed <- backtransf(lowTE.fixed, sm, "lower",
                                npft.ma, warn=comb.fixed)
      uppTE.fixed <- backtransf(uppTE.fixed, sm, "upper",
                                npft.ma, warn=comb.fixed)
      ##
      TE.random <- backtransf(TE.random, sm, "mean",
                              npft.ma, warn=comb.random)
      lowTE.random <- backtransf(lowTE.random, sm, "lower",
                                 npft.ma, warn=comb.random)
      uppTE.random <- backtransf(uppTE.random, sm, "upper",
                                 npft.ma, warn=comb.random)
      ##
      if (!metainf.metacum){
        lowTE.predict <- backtransf(lowTE.predict, sm, "lower",
                                    npft.ma, warn=prediction)
        uppTE.predict <- backtransf(uppTE.predict, sm, "upper",
                                    npft.ma, warn=prediction)
      }
      ##
      if (by){
        npft.w <- n.harmonic.mean.w
        ##
        TE.w     <- backtransf(TE.w, sm, "mean", npft.w)
        lowTE.w  <- backtransf(lowTE.w, sm, "lower", npft.w)
        uppTE.w  <- backtransf(uppTE.w, sm, "upper", npft.w)
      }
    }
    ##
    ## Apply argument 'pscale' to proportions
    ##
    if (sm %in% c("PLOGIT", "PLN", "PRAW", "PAS", "PFT")){
      TE <- pscale*TE
      lowTE <- pscale*lowTE
      uppTE <- pscale*uppTE
      ##
      TE.fixed    <- pscale*TE.fixed
      lowTE.fixed <- pscale*lowTE.fixed
      uppTE.fixed <- pscale*uppTE.fixed
      ##
      TE.random    <- pscale*TE.random
      lowTE.random <- pscale*lowTE.random
      uppTE.random <- pscale*uppTE.random
      ##
      lowTE.predict <- pscale*lowTE.predict
      uppTE.predict <- pscale*uppTE.predict
      ##
      if (by){
        TE.w    <- pscale*TE.w
        lowTE.w <- pscale*lowTE.w
        uppTE.w <- pscale*uppTE.w
      }
    }
  }
  ##    
  if (!comb.fixed){
    TE.fixed    <- NA
    lowTE.fixed <- NA
    uppTE.fixed <- NA
  }
  ##
  if (!comb.random){
    TE.random    <- NA
    lowTE.random <- NA
    uppTE.random <- NA
  }
  ##
  if (!prediction){
    lowTE.predict <- NA
    uppTE.predict <- NA
  }
  ##  
  if (sum(x$w.fixed)>0)
    w.fixed.p <- 100*round(x$w.fixed/sum(x$w.fixed, na.rm=TRUE), 3)
  else
    w.fixed.p <- x$w.fixed
  ##
  if (sum(x$w.random)>0)
    w.random.p <- 100*round(x$w.random/sum(x$w.random, na.rm=TRUE), 3)
  else
    w.random.p <- x$w.random
  ##
  if (metainf.metacum){
    if (sum(x$w.fixed)>0)
      w.fixed.p <- 100*round(x$w.fixed/x$w.all, 3)
    else
      w.fixed.p <- x$w.fixed
    ##
    if (sum(x$w.random)>0)
      w.random.p <- 100*round(x$w.random/x$w.all, 3)
    else
      w.random.p <- x$w.random
  }
  
  
  ##
  ##
  ## (9) Determine column labels
  ##
  ##
  labs <- list()
  ##
  if (missing(leftlabs) || length(leftcols) != length(leftlabs)){
    for (i in seq(along=leftcols)){
      j <- match(leftcols[i], colnames)
      if (!is.na(j))
        labs[[paste("lab.", leftcols[i], sep="")]] <- labnames[j]
    }
  }
  else if (length(leftcols) == length(leftlabs)){
    for (i in seq(along=leftcols)){
      j <- match(leftcols[i], colnames)
      if (!is.na(leftlabs[i]))
        labs[[paste("lab.", leftcols[i], sep="")]] <- leftlabs[i]
      else
        if (!is.na(j))
          labs[[paste("lab.", leftcols[i], sep="")]] <- labnames[j]
    }
  }
  ##
  if (missing(rightlabs) || length(rightcols) != length(rightlabs)){
    for (i in seq(along=rightcols)){
      j <- match(rightcols[i], colnames)
      if (!is.na(j))
        labs[[paste("lab.", rightcols[i], sep="")]] <- labnames[j]
    }
  }
  else if (length(rightcols) == length(rightlabs)){
    for (i in seq(along=rightcols)){
      j <- match(rightcols[i], colnames)
      if (!is.na(rightlabs[i]))
        labs[[paste("lab.", rightcols[i], sep="")]] <- rightlabs[i]
      else
        if (!is.na(j))
          labs[[paste("lab.", rightcols[i], sep="")]] <- labnames[j]
    }
  }
  ##
  if (!slab)
    labs[["lab.studlab"]] <- ""
  
  
  ##
  ##
  ## (10) Define columns in forest plot as well as x- and y-limits
  ##
  ##
  if (by){
    ##
    if (print.byvar){
      if (length(bylab)==0 || bylab=="")
        bylab <- format(bylevs)
      else
        bylab <- paste(bylab,
                       " = ",
                       format(bylevs), sep="")
    }
    else
      bylab <- format(bylevs)
    ##
    if (length(text.fixed.w)==1&n.by>1)
      text.fixed.w <- rep(text.fixed.w, n.by)
    if (length(text.random.w)==1&n.by>1)
      text.random.w <- rep(text.random.w, n.by)
    ##
    modlabs <- c(text.fixed, text.random, hetstat.overall, text.predict,
                 bylab, text.fixed.w, text.random.w, hetstat.w,
                 studlab)
    ##
    TEs    <- c(TE.fixed, TE.random, NA, TE.w, TE)
    lowTEs <- c(lowTE.fixed, lowTE.random, lowTE.predict, lowTE.w, lowTE)
    uppTEs <- c(uppTE.fixed, uppTE.random, uppTE.predict, uppTE.w, uppTE)
    ##
    TEs.study <- c("", "", "", rep("", 3*n.by),
                   ifelse(is.na(TE.orig), lab.NA,
                          format(round(TE.orig, digits),
                                 scientific=FALSE))
                   )
    seTEs.study <- c("", "", "", rep("", 3*n.by),
                     ifelse(is.na(seTE), lab.NA,
                            format(round(seTE, 4),
                                   scientific=FALSE))
                     )
    ##
    w.fixeds  <- c(NA, NA, NA, rep(NA, length(weight.w.p)), w.fixed.p)
    w.randoms <- c(NA, NA, NA, rep(NA, length(weight.w.p)), w.random.p)
    ##
    w.fixeds.text  <- c(100, "--", "", format(c(weight.w.p, w.fixed.p), scientific=FALSE))
    w.randoms.text <- c("--", 100, "", format(c(weight.w.p, w.random.p), scientific=FALSE))
    ##
    sel.fixed  <- w.fixeds.text=="--"
    sel.random <- w.randoms.text=="--"
    ##
    sel.fixed[3+n.by+1:n.by] <- TRUE
    sel.random[3+1:n.by] <- TRUE
    ##
    col.diamond <- c(col.diamond.fixed, col.diamond.random, col.predict,
                     rep(col.diamond.fixed, n.by),
                     rep(col.diamond.random, n.by))
    col.diamond.lines <- c(col.diamond.fixed.lines, col.diamond.random.lines,
                           col.predict.lines,
                           rep(col.diamond.fixed.lines, n.by),
                           rep(col.diamond.random.lines, n.by))
  }
  else{
    modlabs <- c(text.fixed, text.random, hetstat.overall, text.predict, studlab)
    ##
    TEs    <- c(TE.fixed, TE.random, NA, TE)
    lowTEs <- c(lowTE.fixed, lowTE.random, lowTE.predict, lowTE)
    uppTEs <- c(uppTE.fixed, uppTE.random, uppTE.predict, uppTE)
    ##
    TEs.study <- c("", "", "",
                   ifelse(is.na(TE.orig), lab.NA,
                          format(round(TE.orig, digits),
                                 scientific=FALSE))
                   )
    seTEs.study <- c("", "", "",
                     ifelse(is.na(seTE), lab.NA,
                            format(round(seTE, 4),
                                   scientific=FALSE))
                     )
    ##
    w.fixeds  <- c(NA, NA, NA, w.fixed.p)
    w.randoms <- c(NA, NA, NA, w.random.p)
    ##
    w.fixeds.text  <- c(100, "--", "", format(w.fixed.p, scientific=FALSE))
    w.randoms.text <- c("--", 100, "", format(w.random.p, scientific=FALSE))
    ##
    sel.fixed <- w.fixeds.text=="--"
    sel.random <- w.randoms.text=="--"
    ##
    col.diamond <- c(col.diamond.fixed, col.diamond.random, col.predict)
    col.diamond.lines <- c(col.diamond.fixed.lines, col.diamond.random.lines,
                           col.predict.lines)
  }
  ##
  ## Treatment effect and confidence interval
  ##
  if (backtransf & is.relative.effect(sm)){
    effect.format <- format(round(exp(TEs), digits), scientific=FALSE)
    ci.format <- p.ci(format(round(exp(lowTEs), digits), scientific=FALSE),
                      format(round(exp(uppTEs), digits), scientific=FALSE))
  }
  else{
    effect.format <- format(round(TEs, digits), scientific=FALSE)
    ci.format <- p.ci(format(round(lowTEs, digits), scientific=FALSE),
                      format(round(uppTEs, digits), scientific=FALSE))
  }
  ##
  effect.format <- gsub("NA", "  ", effect.format)
  ##
  ## Weights of fixed and random effects model
  ##
  w.fixed.format  <- paste(w.fixeds.text, "%", sep="")
  w.random.format <- paste(w.randoms.text, "%", sep="")
  ##
  w.fixed.format[w.fixed.format=="%"] <- ""
  w.random.format[w.random.format=="%"] <- ""
  ##
  w.fixed.format[sel.fixed] <- "--"
  if (by)
    w.fixed.format[3+2*n.by+1:n.by] <- ""
  w.random.format[sel.random] <- "--"
  if (by)
    w.random.format[3+2*n.by+1:n.by] <- ""
  ##
  ## Treatment estimate and its standard error
  ##
  TE.format <- TEs.study
  seTE.format <- seTEs.study
  ##
  ## Number of Events and patients and person time
  ##
  sum.n.e <- sum(x$n.e, na.rm=TRUE)
  sum.n.c <- sum(x$n.c, na.rm=TRUE)
  sum.e.e <- sum(x$event.e, na.rm=TRUE)
  sum.e.c <- sum(x$event.c, na.rm=TRUE)
  ##
  if (by){
    if (pooled.totals){
      Ne <- c(sum.n.e, sum.n.e, NA, n.e.w, n.e.w, rep(NA, n.by), x$n.e)
      Nc <- c(sum.n.c, sum.n.c, NA, n.c.w, n.c.w, rep(NA, n.by), x$n.c)
    }
    else{
      Ne <- c(NA, NA, NA, rep(NA, 3*n.by), x$n.e)
      Nc <- c(NA, NA, NA, rep(NA, 3*n.by), x$n.c)
    }
    if (pooled.events){
      Ee <- c(sum.e.e, sum.e.e, NA, e.e.w, e.e.w, rep(NA, n.by), x$event.e)
      Ec <- c(sum.e.c, sum.e.c, NA, e.c.w, e.c.w, rep(NA, n.by), x$event.c)
    }
    else{
      Ee <- c(NA, NA, NA, rep(NA, 3*n.by), x$event.e)
      Ec <- c(NA, NA, NA, rep(NA, 3*n.by), x$event.c)
    }
    Te <- c(NA, NA, NA, rep(NA, 3*n.by), x$time.e)
    Tc <- c(NA, NA, NA, rep(NA, 3*n.by), x$time.c)
  }
  else{
    if (pooled.totals){
      Ne <- c(sum.n.e, sum.n.e, NA, x$n.e)
      Nc <- c(sum.n.c, sum.n.c, NA, x$n.c)
    }
    else{
      Ne <- c(NA, NA, NA, x$n.e)
      Nc <- c(NA, NA, NA, x$n.c)
    }
    if (pooled.events){
      Ee <- c(sum.e.e, sum.e.e, NA, x$event.e)
      Ec <- c(sum.e.c, sum.e.c, NA, x$event.c)
    }
    else{
      Ee <- c(NA, NA, NA, x$event.e)
      Ec <- c(NA, NA, NA, x$event.c)
    }
    Te <- c(NA, NA, NA, x$time.e)
    Tc <- c(NA, NA, NA, x$time.c)
  }
  ##
  Ne.format <- ifelse(is.na(Ne), "", format(Ne, scientific=FALSE))
  Nc.format <- ifelse(is.na(Nc), "", format(Nc, scientific=FALSE))
  Ee.format <- ifelse(is.na(Ee), "", format(Ee, scientific=FALSE))
  Ec.format <- ifelse(is.na(Ec), "", format(Ec, scientific=FALSE))
  Te.format <- ifelse(is.na(Te), "", format(Te, scientific=FALSE))
  Tc.format <- ifelse(is.na(Tc), "", format(Tc, scientific=FALSE))
  ##
  if (comb.fixed & comb.random){
    Ne.format[2] <- ""
    Nc.format[2] <- ""
    Ee.format[2] <- ""
    Ec.format[2] <- ""
    Te.format[2] <- ""
    Tc.format[2] <- ""
    if (by){
      Ne.format[3+n.by+1:n.by] <- ""
      Nc.format[3+n.by+1:n.by] <- ""
      Ee.format[3+n.by+1:n.by] <- ""
      Ec.format[3+n.by+1:n.by] <- ""
      Te.format[3+n.by+1:n.by] <- ""
      Tc.format[3+n.by+1:n.by] <- ""
    }
  }
  ##
  ## Mean and standard deviation
  ##
  if (by){
    Me <- c("", "", "", rep("", 3*length(n.e.w)), format(x$mean.e, scientific=FALSE))
    Mc <- c("", "", "", rep("", 3*length(n.c.w)), format(x$mean.c, scientific=FALSE))
    Se <- c("", "", "", rep("", 3*length(n.e.w)), format(x$sd.e, scientific=FALSE))
    Sc <- c("", "", "", rep("", 3*length(n.c.w)), format(x$sd.c, scientific=FALSE))
  }
  else{
    Me <- c("", "", "", format(x$mean.e, scientific=FALSE))
    Mc <- c("", "", "", format(x$mean.c, scientific=FALSE))
    Se <- c("", "", "", format(x$sd.e, scientific=FALSE))
    Sc <- c("", "", "", format(x$sd.c, scientific=FALSE))
  }
  ##
  Me.format <- Me
  Mc.format <- Mc
  Se.format <- Se
  Sc.format <- Sc
  ##
  Me.format[Me.format=="NA"] <- lab.NA
  Mc.format[Mc.format=="NA"] <- lab.NA
  Se.format[Se.format=="NA"] <- lab.NA
  Sc.format[Sc.format=="NA"] <- lab.NA
  ##
  ## Correlation
  ##
  if (by)
    cor <- c("", "", rep("", length(n.e.w)), format(x$cor, scientific=FALSE))
  else
    cor <- c("", "", format(x$cor, scientific=FALSE))
  ##
  cor.format <- cor
  ##
  ## y-axis:
  ##
  if (any(rightcols %in% c("n.e", "n.c")) |
      any(leftcols  %in% c("n.e", "n.c")) |
      (inherits(x, "metainc") &
       (any(rightcols %in% c("time.e", "time.c")) |
        any(leftcols  %in% c("time.e", "time.c")))
       )
      ){
    yHead <- 2
    yHeadadd <- 1
  }
  else{
    yHead <- 1
    yHeadadd <- NA
  }
  ##
  if (!by){
    N <- n.stud
    yTE <- 1:N
    yTE <- yTE
  }
  else{
    ##
    j <- 1
    k <- 0
    yBylab <- rep(NA, n.by)
    yTE <- rep(NA, n.stud)
    yTE.w.fixed <- yBylab
    yTE.w.random <- yBylab
    yTE.w.hetstat <- yBylab
    ##
    for (i in 1:n.by){
      ##
      k.i <- k.all.w[i]
      k <- k+k.i
      ##
      yBylab[i] <- j
      j <- j+1
      ##
      yTE[(k-k.i+1):k] <- j:(j+k.i-1)
      j <- j+k.i
      ##
      ## Fixed effect model
      ##
      if (comb.fixed){
        yTE.w.fixed[i] <- j
        j <- j+1
      }
      else
        yTE.w.fixed[i] <- NA
      ##
      ## Random effect model
      ##
      if (comb.random){
        yTE.w.random[i] <- j
        j <- j+1
      }
      else
        yTE.w.random[i] <- NA
      ##
      ## Only pooled totals
      ##
      if (pooled.totals&!(comb.fixed|comb.random)){
        yTE.w.fixed[i] <- j
        j <- j+1
      }
      ##
      ## Heterogeneity statistics
      ##
      if (hetstat){
        yTE.w.hetstat[i] <- j
        j <- j+1
      }
      else
        yTE.w.hetstat[i] <- NA
      ##
      j <- j+1
    }
    ##
    yTE.w <- c(yTE.w.fixed, yTE.w.random, yTE.w.hetstat)
  }
  ##
  ## x-axis:
  ##
  if (notmiss.xlim && is.numeric(xlim[1]))
    if (is.relative.effect(sm))
      xlim <- log(xlim)
  ##
  if (is.null(xlim)){
    if (metaprop){
      xlim <- c(min(c(lowTE, lowTE.predict), na.rm=TRUE),
                max(c(uppTE, uppTE.predict), na.rm=TRUE))
      ##
      if (!is.na(ref) && ref < xlim[1])
        xlim[1] <- ref
      if (!is.na(ref) && ref > xlim[2])
        xlim[2] <- ref
    }
    else{
      sel.low <- is.finite(lowTE)
      sel.upp <- is.finite(uppTE)
      xlim <- c(min(c(lowTE[sel.low], lowTE.predict), na.rm=TRUE),
                max(c(uppTE[sel.upp], uppTE.predict), na.rm=TRUE))
      ##
      if (!is.na(ref) && ref < xlim[1])
        xlim[1] <- ref
      if (!is.na(ref) && ref > xlim[2])
        xlim[2] <- ref
    }
  }
  ##
  symmetric <- FALSE
  ##
  if (!is.null(xlim) && is.character(xlim[1])){
    ##
    xlim <- setchar(xlim, "symmetric",
                    "should be a numeric vector (min, max) or the character string \"symmetric\"")
    symmetric <- TRUE
    ##
    if (metaprop){
      xlim <- c(min(c(lowTE, lowTE.predict), na.rm=TRUE),
                max(c(uppTE, uppTE.predict), na.rm=TRUE))
    }
    else{
      sel.low <- is.finite(lowTE)
      sel.upp <- is.finite(uppTE)
      minTE <- min(c(lowTE[sel.low], lowTE.predict), na.rm=TRUE)
      maxTE <- max(c(uppTE[sel.upp], uppTE.predict), na.rm=TRUE)
      ##
      if (minTE < 0 & maxTE < 0)
        xlim <- c(minTE, -minTE)
      else if (minTE > 0 & maxTE > 0)
        xlim <- c(-maxTE, maxTE)
      else
        xlim <- c(-max(abs(c(minTE, maxTE))), max(abs(c(minTE, maxTE))))
    }
  }
  ##
  if (!is.na(ref) &&
      round(xlim[2]-ref, 6) == round(ref-xlim[1], 6))
    symmetric <- TRUE
  ##  
  if (by)
    max.yTE <- max(c(yTE, yTE.w), na.rm=TRUE)
  else
    max.yTE <- max(yTE, na.rm=TRUE)
  ##  
  if (!is.na(ref) & missing(xlab.pos))
    if (ref <= min(xlim) | ref >= max(xlim))
      xlab.pos <- mean(xlim)
  ##
  if (!is.na(ref) & missing(smlab.pos))
    if (ref <= min(xlim) | ref >= max(xlim))
      smlab.pos <- mean(xlim)
  ##
  if (is.null(xlab.pos) || is.na(xlab.pos))
    xlab.pos <- mean(xlim)
  ##
  if (is.null(smlab.pos) || is.na(smlab.pos))
    smlab.pos <- mean(xlim)
  ##
  yTE.fixed   <- NA
  yTE.random  <- NA
  yTE.predict <- NA
  ##
  if (comb.fixed & comb.random & overall){
    yTE.fixed  <- max.yTE+2
    yTE.random <- max.yTE+3
  }
  ##
  if (comb.fixed & !comb.random & overall){
    yTE.fixed <- max.yTE+2
  }
  ##
  if (!comb.fixed & comb.random & overall){
    yTE.random <- max.yTE+2
  }
  ##
  if (!comb.fixed & !comb.random & pooled.totals & overall){
    yTE.fixed  <- max.yTE+2
    if (missing(text.fixed))
        text.fixed <- "Overall"
  }
  ##
  if (!is.na(yTE.random) & prediction)
    yTE.predict <- yTE.random+1
  ##
  if (overall.hetstat)
    if (is.na(yTE.fixed) & is.na(yTE.random))
      yTE.hetstat <- max.yTE+2
    else
      yTE.hetstat <- max(max.yTE, yTE.fixed, yTE.random, yTE.predict, na.rm=TRUE)+1
  else if (!overall.hetstat & addspace)
      yTE.hetstat <- max(max.yTE, yTE.fixed, yTE.random, yTE.predict, na.rm=TRUE)+1
  else
    yTE.hetstat <- NA
  ##
  if (!comb.fixed & !pooled.totals) text.fixed <- ""
  if (!comb.random) text.random <- ""
  if (!prediction) text.predict <- ""
  ##  
  yTE         <- yHead + yTE + addspace
  yTE.fixed   <- yHead + yTE.fixed + addspace
  yTE.random  <- yHead + yTE.random + addspace
  yTE.predict <- yHead + yTE.predict + addspace
  yTE.hetstat <- yHead + yTE.hetstat + addspace
  ##
  if (by){
    yBylab <- yHead + yBylab + addspace
    yTE.w  <- yHead + yTE.w + addspace
  }
  ##  
  if (by){
    yLab <- c(yHead,
              yTE.fixed, yTE.random, yTE.hetstat, yTE.predict,
              yBylab, yTE.w,
              yTE)
    ##
    yS <- c(yHead, yTE.fixed, yTE.random, yTE.predict, yTE.w, yTE)
  }
  else{
    yLab <- c(yHead, yTE.fixed, yTE.random, yTE.hetstat, yTE.predict, yTE)
    yS   <- c(yHead, yTE.fixed, yTE.random, yTE.predict, yTE)
  }
  
  
  ##
  ##
  ## (11) Format columns in forest plot
  ##
  ##
  formatcol <- function(x, y, rows, just="right"){
    if (just=="left")
      xpos <- 0
    if (just=="center")
      xpos <- 0.5
    if (just=="right")
      xpos <- 1
    ##
    res <- list(labels=
                lapply(c(x, as.list(y)),
                       textGrob, x=xpos, just=just,
                       gp=gpar(
                         fontsize=fs.study,
                         fontface=ff.study)
                       ),
                rows=rows)
    ##
    ## Study label:
    ##
    res$labels[[1]] <- textGrob(x,
                                x=xpos, just=just,
                                gp=gpar(
                                  fontsize=fs.heading,
                                  fontface=ff.heading)
                                )
    ##
    ## Fixed effect estimate:
    ##
    res$labels[[2]] <- textGrob(y[1],
                                x=xpos, just=just,
                                gp=gpar(
                                  fontsize=fs.fixed,
                                  fontface=ff.fixed)
                                )
    ##
    ## Random effects estimate:
    ##
    res$labels[[3]] <- textGrob(y[2],
                                x=xpos, just=just,
                                gp=gpar(
                                  fontsize=fs.random,
                                  fontface=ff.random)
                                )
    ##
    ## Prediction interval:
    ##
    res$labels[[4]] <- textGrob(y[3],
                                x=xpos, just=just,
                                gp=gpar(
                                  fontsize=fs.predict,
                                  fontface=ff.predict)
                                )
    ##
    if (by)
      for (i in 1:n.by){
        ##
        ## Fixed effect estimates:
        ##
        res$labels[[4+i]] <- textGrob(y[3+i],
                                      x=xpos, just=just,
                                      gp=
                                      gpar(
                                           fontsize=fs.fixed,
                                           fontface=ff.fixed,
                                           col=col.by)
                                      )
        ##
        ## Random effects estimates:
        ##
        res$labels[[4+n.by+i]] <- textGrob(y[3+n.by+i],
                                           x=xpos, just=just,
                                           gp=
                                           gpar(
                                                fontsize=fs.random,
                                                fontface=ff.random,
                                                col=col.by)
                                           )
      }
    ##
    res
  }
  ##
  col.studlab <- list(labels=
                      lapply(as.list(c(labs[["lab.studlab"]], modlabs)),
                             textGrob, x=xpos.studlab, just=just.studlab,
                             gp=gpar(
                               fontsize=fs.study.labels,
                               fontface=ff.study.labels)),
                      rows=yLab
                      )
  ##
  ## Study label:
  ##
  col.studlab$labels[[1]] <- textGrob(labs[["lab.studlab"]],
                                      x=xpos.studlab, just=just.studlab,
                                      gp=gpar(
                                        fontsize=fs.heading,
                                        fontface=ff.heading)
                                      )
  ##
  ## Fixed effect estimate:
  ##
  col.studlab$labels[[2]] <- textGrob(text.fixed,
                                      x=xpos.studlab, just=just.studlab,
                                      gp=gpar(
                                        fontsize=fs.fixed.labels,
                                        fontface=ff.fixed.labels)
                                      )
  ##
  ## Random effects estimate:
  ##
  col.studlab$labels[[3]] <- textGrob(text.random,
                                      x=xpos.studlab, just=just.studlab,
                                      gp=gpar(
                                        fontsize=fs.random.labels,
                                        fontface=ff.random.labels)
                                      )
  ##
  ## Heterogeneity statistics:
  ##
  col.studlab$labels[[4]] <- textGrob(hetstat.overall,
                                      x=xpos.studlab, just=just.studlab,
                                      gp=gpar(
                                        fontsize=fs.hetstat,
                                        fontface=ff.hetstat)
                                      )
  ##
  ## Prediction interval:
  ##
  col.studlab$labels[[5]] <- textGrob(text.predict,
                                      x=xpos.studlab, just=just.studlab,
                                      gp=gpar(
                                        fontsize=fs.predict.labels,
                                        fontface=ff.predict.labels)
                                      )
  ##
  if (by){
    for (i in 1:n.by){
      ##
      ## Subgroup labels:
      ##
      col.studlab$labels[[5+i]] <- textGrob(bylab[i],
                                            x=xpos.studlab, just=just.studlab,
                                            gp=
                                            gpar(
                                                 fontsize=fs.heading,
                                                 fontface=ff.heading,
                                                 col=col.by)
                                            )
      ##
      ## Fixed effect estimates:
      ##
      col.studlab$labels[[5+n.by+i]] <- textGrob(text.fixed.w[i],
                                                 x=xpos.studlab, just=just.studlab,
                                                 gp=
                                                 gpar(
                                                      fontsize=fs.fixed.labels,
                                                      fontface=ff.fixed.labels,
                                                      col=col.by)
                                                 )
      ##
      ## Random effects estimates:
      ##
      col.studlab$labels[[5+2*n.by+i]] <- textGrob(text.random.w[i],
                                                   x=xpos.studlab, just=just.studlab,
                                                   gp=
                                                   gpar(
                                                        fontsize=fs.random.labels,
                                                        fontface=ff.random.labels,
                                                        col=col.by)
                                                   )
      ##
      ## Heterogeneity statistics:
      ##
      col.studlab$labels[[5+3*n.by+i]] <- textGrob(hetstat.w[i],
                                                   x=xpos.studlab, just=just.studlab,
                                                   gp=
                                                   gpar(
                                                        fontsize=fs.hetstat,
                                                        fontface=ff.hetstat,
                                                        col=col.by)
                                                   )
    }
  }
  ##
  col.effect <- formatcol(labs[["lab.effect"]], effect.format, yS, just.cols)
  ##
  col.ci <- formatcol(labs[["lab.ci"]], ci.format, yS, just.cols)
  ##
  col.w.fixed  <- formatcol(labs[["lab.w.fixed"]], w.fixed.format, yS, just.cols)
  col.w.random <- formatcol(labs[["lab.w.random"]], w.random.format, yS, just.cols)
  ##
  col.TE <- formatcol(labs[["lab.TE"]], TE.format, yS, just.cols)
  col.seTE <- formatcol(labs[["lab.seTE"]], seTE.format, yS, just.cols)
  ##  
  col.n.e <- formatcol(labs[["lab.n.e"]], Ne.format, yS, just.cols)
  col.n.c <- formatcol(labs[["lab.n.c"]], Nc.format, yS, just.cols)
  ##
  col.event.e <- formatcol(labs[["lab.event.e"]], Ee.format, yS, just.cols)
  col.event.c <- formatcol(labs[["lab.event.c"]], Ec.format, yS, just.cols)
  ##
  col.mean.e <- formatcol(labs[["lab.mean.e"]], Me.format, yS, just.cols)
  col.mean.c <- formatcol(labs[["lab.mean.c"]], Mc.format, yS, just.cols)
  ##
  col.sd.e <- formatcol(labs[["lab.sd.e"]], Se.format, yS, just.cols)
  col.sd.c <- formatcol(labs[["lab.sd.c"]], Sc.format, yS, just.cols)
  ##
  col.cor  <- formatcol(labs[["lab.cor"]], cor.format, yS, just.cols)
  ##
  col.time.e <- formatcol(labs[["lab.time.e"]], Te.format, yS, just.cols)
  col.time.c <- formatcol(labs[["lab.time.c"]], Tc.format, yS, just.cols)
  ##  
  col.forest <- list(eff=TEs,
                     low=lowTEs,
                     upp=uppTEs,
                     rows=yS[-1],
                     ##
                     ## "p" means prediction, "s" means summary, "n" means normal
                     ##
                     type=c("s", "s", "p", rep("s", length(TEs)-length(TE)-3), rep("n", length(TE))),
                     col=c(rep("", length(TEs)-length(TE)), col.i),
                     col.square=c(rep("", length(TEs)-length(TE)), col.square),
                     col.square.lines=c(rep("", length(TEs)-length(TE)), col.square.lines),
                     col.i.inside.square=c(rep("", length(TEs)-length(TE)), col.i.inside.square),
                     col.diamond=c(col.diamond, rep("", length(TE))),
                     col.diamond.lines=c(col.diamond.lines, rep("", length(TE)))
                     )
  ##
  ## Sizes of squares
  ##
  if (weight=="same"){
    information <- rep(0.9, length(TEs))
  }
  else{
    ##
    if (weight=="fixed")
      information <- sqrt(w.fixeds)
    else if (weight=="random")
      information <- sqrt(w.randoms)
    ## Square height equal to 1 for most precise study result
    information <- information/max(information, na.rm=TRUE)
    ## Same/maximum polygon height for all meta-analytical results
    ## (both overall and subgroup results)
    information[is.na(information)] <- 1
  }
  ##
  col.forest$sizes <- information
  col.forest$sizes <- col.forest$sizes * squaresize
  ##
  ## Width of column 3
  col.forestwidth <- plotwidth
  ##
  ## Range on the x-axis for column 3
  col.forest$range <- xlim
  ##  
  cols <- list(col.studlab=col.studlab,
               col.effect=col.effect,
               col.ci=col.ci,
               col.w.fixed=col.w.fixed,
               col.w.random=col.w.random,
               col.TE=col.TE,
               col.seTE=col.seTE)
  ##
  cols[["col.n.e"]] <- col.n.e
  cols[["col.n.c"]] <- col.n.c
  cols[["col.event.e"]] <- col.event.e
  cols[["col.event.c"]] <- col.event.c
  ##
  cols[["col.mean.e"]] <- col.mean.e
  cols[["col.mean.c"]] <- col.mean.c
  cols[["col.sd.e"]] <- col.sd.e
  cols[["col.sd.c"]] <- col.sd.c
  ##
  cols[["col.cor"]] <- col.cor
  ##
  cols[["col.time.e"]] <- col.time.e
  cols[["col.time.c"]] <- col.time.c
  ##
  if (newcols){
    if (by){
      for (i in seq(along=rightcols.new)){
        tname <- paste("col.", rightcols.new[i], sep="")
        if (length(dataset[[rightcols.new[i]]])!=0)
          tmp.r <- dataset[[rightcols.new[i]]]
        else if (length(x[[rightcols.new[i]]])!=0)
          tmp.r <- x[[rightcols.new[i]]]
        if (is.factor(tmp.r))
          tmp.r <- as.character(tmp.r)
        tmp.r <- ifelse(is.na(tmp.r), "", tmp.r)
        cols[[tname]] <- formatcol(rightlabs.new[i],
                                   c("", "", "",
                                     rep("", length(TE.w)),
                                     tmp.r),
                                   yS,
                                   just=just.addcols)
      }
      for (i in seq(along=leftcols.new)){
        tname <- paste("col.", leftcols.new[i], sep="")
        if (length(dataset[[leftcols.new[i]]])!=0)
          tmp.l <- dataset[[leftcols.new[i]]]        
        else if (length(x[[leftcols.new[i]]])!=0)
          tmp.l <- x[[leftcols.new[i]]]
        if (is.factor(tmp.l))
          tmp.l <- as.character(tmp.l)
        tmp.l <- ifelse(is.na(tmp.l), "", tmp.l)
        cols[[tname]] <- formatcol(leftlabs.new[i],
                                   c("", "", "",
                                     rep("", length(TE.w)),
                                     tmp.l),
                                   yS,
                                   just=just.addcols)
      }
    }
    else{
      for (i in seq(along=rightcols.new)){
        tname <- paste("col.", rightcols.new[i], sep="")
        if (length(dataset[[rightcols.new[i]]])!=0)
          tmp.r <- dataset[[rightcols.new[i]]]
        else if (length(x[[rightcols.new[i]]])!=0)
          tmp.r <- x[[rightcols.new[i]]]
        if (is.factor(tmp.r))
          tmp.r <- as.character(tmp.r)
        tmp.r <- ifelse(is.na(tmp.r), "", tmp.r)
        cols[[tname]] <- formatcol(rightlabs.new[i],
                                   c("", "", "", tmp.r),
                                   yS,
                                   just=just.addcols)
      }
      for (i in seq(along=leftcols.new)){
        tname <- paste("col.", leftcols.new[i], sep="")
        if (length(dataset[[leftcols.new[i]]])!=0)
          tmp.l <- dataset[[leftcols.new[i]]]        
        else if (length(x[[leftcols.new[i]]])!=0)
          tmp.l <- x[[leftcols.new[i]]]
        if (is.factor(tmp.l))
          tmp.l <- as.character(tmp.l)
        tmp.l <- ifelse(is.na(tmp.l), "", tmp.l)
        cols[[tname]] <- formatcol(leftlabs.new[i],
                                   c("", "", "", tmp.l),
                                   yS,
                                   just=just.addcols)
      }
    }
  }
  ##  
  col.lab.e <- list(labels=list(textGrob(lab.e,
                      x=1, just=just.cols,
                      gp=gpar(
                        fontsize=fs.heading,
                        fontface=ff.heading)
                      )),
                    rows=1)
  ##
  col.lab.c <- list(labels=list(textGrob(lab.c,
                      x=1, just=just.cols,
                      gp=gpar(
                        fontsize=fs.heading,
                        fontface=ff.heading)
                      )),
                    rows=1)
  ##  
  leftcols  <- paste("col.", leftcols, sep="")
  rightcols <- paste("col.", rightcols, sep="")
  
  
  ##
  ##
  ## (12) Definition of auxiliary plot functions
  ##
  ##
  drawLabelCol <- function(col, j) {
    ##
    ## Function to draw a cell in a text column
    ##
    for (i in 1:length(col$rows)) {
      if (!is.na(col$rows[i])){
        pushViewport(viewport(layout.pos.row=col$rows[i], layout.pos.col=j))
        ##
        ## Labels are grobs containing their location so just
        ## have to grid.draw() them
        ##
        grid.draw(col$labels[[i]])
        popViewport()
      }
    }
  }
  ##
  drawNormalCI <- function(low, eff, upp, size, min, max,
                           col, col.square, col.square.lines,
                           col.i.inside.square){
    ##
    ## Function to draw a non-summary rect-plus-CI
    ##
    ## NOTE the use of "native" units to position relative to
    ## the x-axis scale, and "snpc" units to size relative to
    ## the height of the row
    ## ("snpc" stands for "square normalised parent coordinates"
    ##  which means that the value is calculated as a proportion
    ##  of the width and height of the current viewport and the
    ##  physically smaller of these is used)
    ##
    if (!is.na(eff)){
      ##
      ## Draw lines in colour "col.i.inside.square" if totally inside rect
      ##
      if (!is.na(size)){
        TElineCol <- if (size > 0 &&
                       (convertX(unit(eff, "native") + unit(0.5*size, "lines"),
                                 "native", valueOnly=TRUE) > upp) &&
                       (convertX(unit(eff, "native") - unit(0.5*size, "lines"),
                                 "native", valueOnly=TRUE) < low))
          col.i.inside.square
        else
          col
      }
    }
    ##
    if (!is.na(eff) && (eff >= min & eff <= max)){
      if (!is.na(size) && size > 0){
        grid.rect(x=unit(eff, "native"),
                  width=unit(size, "snpc"),
                  height=unit(size, "snpc"),
                  gp=gpar(fill=col.square, col=col.square.lines))
        grid.lines(x=unit(c(eff, eff), "native"),
                   y = unit(c(0.4, 0.6), "npc"),
                   gp=gpar(col=TElineCol, lwd=lwd))
      }
      else
        grid.lines(x=unit(c(eff, eff), "native"),
                   y = unit(c(0.4, 0.6), "npc"),
                   gp=gpar(col=TElineCol, lwd=lwd))
    }
    if (!is.na(eff)){
      ##
      ## Draw lines in colour "col.i.inside.square" if totally inside rect
      ##
      if (!is.na(size)){
        lineCol <- if (size > 0 &&
                       (convertX(unit(eff, "native") + unit(0.5*size, "lines"),
                                 "native", valueOnly=TRUE) > upp) &&
                       (convertX(unit(eff, "native") - unit(0.5*size, "lines"),
                                 "native", valueOnly=TRUE) < low))
          col.i.inside.square
        else
          col
        ##
        ## Draw arrow if exceed col range
        ## convertX() used to convert between coordinate systems
        ##
        if (!is.na(low) && !is.na(upp) &&
            (low >= min & upp <= max))
          grid.lines(x=unit(c(low, upp), "native"), y=0.5,
                     gp=gpar(col=lineCol, lwd=lwd))
        ##
        if (!is.na(low) && !is.na(upp) &&
            (low < min & upp > max))
          grid.lines(x=unit(c(min, max), "native"), y=0.5,
                     gp=gpar(col=lineCol, lwd=lwd))
        ##
        if (!is.na(low) && !is.na(upp) &&
            (low < min & (upp <= max & upp > min)))
          grid.lines(x=unit(c(min, upp), "native"), y=0.5,
                     gp=gpar(col=lineCol, lwd=lwd))
        ##
        if (!is.na(low) && !is.na(upp) &&
            ((low >= min & low < max) & upp > max))
          grid.lines(x=unit(c(low, max), "native"), y=0.5,
                     gp=gpar(col=lineCol, lwd=lwd))
        ##
        if (!is.na(low) && low < min)
          grid.lines(x=unit(c(min-0.00001, min), "native"), y=0.5,
                     gp=gpar(col=lineCol, lwd=lwd),
                     arrow=arrow(ends="first", length=unit(0.05, "inches")))
        if (!is.na(upp) && upp > max)
          grid.lines(x=unit(c(max, max+0.00001), "native"), y=0.5,
                     gp=gpar(col=lineCol, lwd=lwd),
                     arrow=arrow(ends="last", length=unit(0.05, "inches")))
      }
    }
  }
  ##
  drawSummaryCI <- function(low, eff, upp, size, min, max,
                            col.diamond, col.diamond.lines) {
    ##
    ## Function to draw a summary "diamond"
    ##
    ## Not sure how to calc the heights of the diamonds so
    ## I'm just using half the height of the equivalent rect
    ##
    if (!is.na(eff) &&
        ((min <= eff & eff <= max) |
         (min <= low & low <= max) |
         (min <= upp & upp <= max))
        )
      grid.polygon(x=unit(c(low, eff, upp, eff), "native"),
                   y=unit(0.5 + c(0, 0.3*size, 0, -0.3*size), "npc"),
                   gp=gpar(fill=col.diamond, col=col.diamond.lines))
  }
  ##
  drawPredictionCI <- function(low, upp, size, min, max,
                               col.predict, col.predict.lines) {
    ##
    ## Function to draw a prediction interval
    ##
    if (!(is.na(low) | is.na(upp))){
      ## Plot prediction interval only within plotting range
      if ((min <= low & low <= max) |
          (min <= upp & upp <= max))
        grid.polygon(x=unit(c(low, low, upp, upp), "native"),
                     y=unit(0.5 + c(-0.1*size, 0.1*size, 0.1*size, -0.1*size), "npc"),
                     gp=gpar(fill=col.predict, col=col.predict.lines))
    }
  }
  ##
  drawDataCol <- function(col, j) {
    ##
    ## Function to draw a "data" column
    ##
    pushViewport(viewport(layout.pos.col=j, xscale=col$range))
    ##
    ymax.line <- max(yS, na.rm=TRUE)-1
    ##
    ## Reference line:
    ##
    if (!is.na(ref) && (col$range[1] <= ref & ref <= col$range[2]))
      grid.lines(x=unit(ref, "native"),
                 y=unit(c(0, ymax.line), "lines"),
                 gp=gpar(lwd=lwd))
    if (comb.fixed & overall)
      if (col$range[1] <= TE.fixed & TE.fixed <= col$range[2])
        if (!is.null(lty.fixed))
          grid.lines(x=unit(TE.fixed, "native"),
                     y=unit(c(0, ymax.line), "lines"),
                     gp=gpar(lty=lty.fixed, lwd=lwd)) # lty="dashed"
    if (comb.random & overall)
      if (col$range[1] <= TE.random & TE.random <= col$range[2])
        if (!is.null(lty.random))
          grid.lines(x=unit(TE.random, "native"),
                     y=unit(c(0, ymax.line), "lines"),
                     gp=gpar(lty=lty.random, lwd=lwd))
    if (log.xaxis){
      if (is.null(at)){
        x1000 <- c(0.001, 0.1, 1,  10, 1000)
        x100  <- c(0.01 , 0.1, 1,  10, 100)
        x10   <- c(0.1  , 0.5, 1,   2, 10)
        x5    <- c(0.2  , 0.5, 1,   2, 5)
        x2    <- c(0.5  , 1, 2)
        x1.5  <- c(0.75 , 1, 1.5)
        x1.25 <- c(0.8  , 1, 1.25)
        x1    <- c(0.9  , 1, 1.1)
        ##
        tval.min <- min(exp(col$range[1]), 1)
        tval.max <- max(exp(col$range[2]), 1)
        ##
        if (all(x1000 >= tval.min) &
            all(x1000 <= tval.max))
          label <- x1000
        else if (all(x100 >= tval.min) &
            all(x100 <= tval.max))
          label <- x100
        else if (all(x10 >= tval.min) &
            all(x10 <= tval.max))
          label <- x10
        else if (all(x5 >= tval.min) &
            all(x5 <= tval.max))
          label <- x5
        else if (all(x2 >= tval.min) &
            all(x2 <= tval.max))
          label <- x2
        else if (all(x1.5 >= tval.min) &
            all(x1.5 <= tval.max))
          label <- x1.5
        else if (all(x1.25 >= tval.min) &
            all(x1.25 <= tval.max))
          label <- x1.25
        else if (all(x1 >= tval.min) &
            all(x1 <= tval.max))
          label <- x1
        else
          label <- 1
        ##
        if (notmiss.xlim && is.numeric(xlim[1])){
          if (exp(min(xlim)) < min(label))
            label <- c(exp(min(xlim)), label)
          if (exp(max(xlim)) > max(label))
            label <- c(label, exp(max(xlim)))
        }
        at <- log(label)
        label <- round(label, 2)
      }
      else{
        if (length(label)==1 && is.logical(label) && label)
          label <- round(at, 2)
        at <- log(at)
      }
      grid.xaxis(at=at, label=label,
                 gp=gpar(fontsize=fs.axis, fontface=ff.axis, lwd=lwd))
    }
    else{
      if (is.null(at))
        grid.xaxis(gp=gpar(fontsize=fs.axis, fontface=ff.axis, lwd=lwd))
      else
        if ((length(label)==1 && is.logical(label) && label) |
            (length(label)>=1 & !is.logical(label)))
          grid.xaxis(at=at, label=label,
                     gp=gpar(fontsize=fs.axis, fontface=ff.axis, lwd=lwd))
        else
          grid.xaxis(at=at,
                     gp=gpar(fontsize=fs.axis, fontface=ff.axis, lwd=lwd))
    }
    ##
    ## sm-Label on top:
    ##
    grid.text(smlab,
              x=unit(smlab.pos, "native"),
              y=unit(1, "npc"),
              just=c("center", "top"),
              gp=gpar(fontsize=fs.smlab, fontface=ff.smlab))
    ##
    ## Label on x-axis:
    ##
    grid.text(xlab,
              x=unit(xlab.pos, "native"),
              ##y=-0.2,
              y=unit(-2.5-1*(symmetric && (label.left!="" | label.right!="")), "lines"),
              just=c("center", "top"),
              gp=gpar(fontsize=fs.xlab, fontface=ff.xlab))
    ##
    ## Left and right label on x-axis:
    ##
    if (symmetric && (label.left!="" | label.right!="")){
      grid.text(label.left,
                x=unit(xlab.pos-(xlim[2]-xlim[1])/30, "native"),
                y=unit(-2.5, "lines"),
                just=c("right", "center"),
                gp=gpar(fontsize=fs.lr, fontface=ff.lr))
      grid.text(label.right,
                x=unit(xlab.pos+(xlim[2]-xlim[1])/30, "native"),
                y=unit(-2.5, "lines"),
                just=c("left", "center"),
                gp=gpar(fontsize=fs.lr, fontface=ff.lr))
    }
    ##
    popViewport()
    for (i in 1:length(col$rows)) {
      if (!is.na(col$rows[i])){
        pushViewport(viewport(layout.pos.row=col$rows[i], layout.pos.col=j,
                              xscale=col$range))
        if (col$type[i] == "n")
          drawNormalCI(low=col$low[i], eff=col$eff[i], upp=col$upp[i],
                       size=col$sizes[i],
                       min=col$range[1], max=col$range[2],
                       col=col$col[i], col.square=col$col.square[i],
                       col.square.lines=col$col.square.lines[i],
                       col.i.inside.square=col$col.i.inside.square[i])
        else if (col$type[i] == "s")
          drawSummaryCI(low=col$low[i], eff=col$eff[i], upp=col$upp[i],
                        size=col$sizes[i],
                        min=col$range[1], max=col$range[2],
                        col.diamond=col.diamond[i],
                        col.diamond.lines=col.diamond.lines[i])
        else if (col$type[i] == "p"){
          drawPredictionCI(low=col$low[i], upp=col$upp[i],
                           size=col$sizes[i],
                           min=col$range[1], max=col$range[2],
                           col.predict=col.diamond[i],
                           col.predict.lines=col.diamond.lines[i])
        }
        popViewport()
      }
    }
  }
  
  
  ##
  ##
  ## (13) Calculate width of columns in forest plot
  ##
  ##
  wcalc <- function(x)
    max(unit(rep(1, length(x)), "grobwidth", x))
  ##
  for (i in seq(along=leftcols)){
    if (i==1)
      if (leftcols[[i]]=="col.studlab")
        x1 <- unit.c(wcalc(cols[[leftcols[i]]]$labels[-(c(4, 5+3*n.by+1:n.by))]))
      else
        x1 <- unit.c(wcalc(cols[[leftcols[i]]]$labels))
    else
      if (leftcols[[i]]=="col.studlab")
        x1 <- unit.c(x1,
                     colgap.left,
                     wcalc(cols[[leftcols[i]]]$labels[-(c(4, 5+3*n.by+1:n.by))]))
      else
        x1 <- unit.c(x1,
                     colgap.left,
                     wcalc(cols[[leftcols[i]]]$labels))
  }
  ##
  x1 <- unit.c(x1, colgap.forest.left, col.forestwidth)
  ##
  if (rsel){
    for (i in seq(along=rightcols)){
      x1 <- unit.c(x1,
                   if (i==1) colgap.forest.right else colgap.right,
                   wcalc(cols[[rightcols[i]]]$labels))
    }
  }
  
  
  ##
  ##
  ## (14) Generate forest plot
  ##
  ##  
  if (new)
    grid.newpage()
  ##
  if (by)
    nrow <- max(c(yTE, yTE.fixed, yTE.random, yTE.hetstat, yTE.predict, yTE.w), na.rm=TRUE)
  else
    nrow <- max(c(yTE, yTE.fixed, yTE.random, yTE.hetstat, yTE.predict), na.rm=TRUE)
  ##
  pushViewport(viewport(layout=grid.layout(
                          nrow,
                          length(x1),
                          widths=x1,
                          heights=unit(rep(1, max(yTE)), "lines"))))
  ##
  ##
  j <- 1
  ##
  for (i in seq(along=leftcols)){
    drawLabelCol(cols[[leftcols[i]]], j)
    ##
    if (!is.na(yHeadadd)){
      if (!is.null(lab.e.attach.to.col)){
        if (leftcols[i]==paste("col.", lab.e.attach.to.col, sep=""))
          drawLabelCol(col.lab.e, j)
      }
      else if (inherits(x, "metabin")){
        if (leftcols[i]=="col.n.e" & just.cols=="right")
          drawLabelCol(col.lab.e, j)
        else if (leftcols[i]=="col.event.e" & just.cols %in% c("left", "center"))
          drawLabelCol(col.lab.e, j)
      }
      else if (inherits(x, "metacont")){
        if (leftcols[i]=="col.sd.e" & just.cols=="right")
          drawLabelCol(col.lab.e, j)
        else if (leftcols[i]=="col.mean.e" & just.cols %in% c("left", "center"))
          drawLabelCol(col.lab.e, j)
      }
      else if (inherits(x, "metainc")){
        if (leftcols[i]=="col.time.e" & just.cols=="right")
          drawLabelCol(col.lab.e, j)
        else if (leftcols[i]=="col.event.e" & just.cols %in% c("left", "center"))
          drawLabelCol(col.lab.e, j)
      }
      ##
      if (!is.null(lab.c.attach.to.col)){
        if (leftcols[i]==paste("col.", lab.c.attach.to.col, sep=""))
          drawLabelCol(col.lab.c, j)
      }
      else if (inherits(x, "metabin")){
        if (leftcols[i]=="col.n.c" & just.cols=="right")
          drawLabelCol(col.lab.c, j)
        else if (leftcols[i]=="col.event.c" & just.cols %in% c("left", "center"))
          drawLabelCol(col.lab.c, j)
      }
      else if (inherits(x, "metacont")){
        if (leftcols[i]=="col.sd.c" & just.cols=="right")
          drawLabelCol(col.lab.c, j)
        else if (leftcols[i]=="col.mean.c" & just.cols %in% c("left", "center"))
          drawLabelCol(col.lab.c, j)
      }
      else if (inherits(x, "metainc")){
        if (leftcols[i]=="col.time.c" & just.cols=="right")
          drawLabelCol(col.lab.c, j)
        else if (leftcols[i]=="col.event.c" & just.cols %in% c("left", "center"))
          drawLabelCol(col.lab.c, j)
      }
    }
    ##
    j <- j+2
  }
  ##
  drawDataCol(col.forest, j)
  j <- j+2
  ##
  if (rsel){
    for (i in seq(along=rightcols)){
      drawLabelCol(cols[[rightcols[i]]], j)
      ##
      if (!is.na(yHeadadd)){
        if (!is.null(lab.e.attach.to.col)){
          if (rightcols[i]==paste("col.", lab.e.attach.to.col, sep=""))
            drawLabelCol(col.lab.e, j)
        }
        else if (inherits(x, "metabin")){
          if (rightcols[i]=="col.n.e" & just.cols=="right")
            drawLabelCol(col.lab.e, j)
          else if (rightcols[i]=="col.event.e" & just.cols %in% c("left", "center"))
            drawLabelCol(col.lab.e, j)
        }
        else if (inherits(x, "metacont")){
          if (rightcols[i]=="col.sd.e" & just.cols=="right")
            drawLabelCol(col.lab.e, j)
          else if (rightcols[i]=="col.mean.e" & just.cols %in% c("left", "center"))
            drawLabelCol(col.lab.e, j)
        }
        else if (inherits(x, "metainc")){
          if (rightcols[i]=="col.time.e" & just.cols=="right")
            drawLabelCol(col.lab.e, j)
          else if (rightcols[i]=="col.event.e" & just.cols %in% c("left", "center"))
            drawLabelCol(col.lab.e, j)
        }
        ##
        if (!is.null(lab.c.attach.to.col)){
          if (rightcols[i]==paste("col.", lab.c.attach.to.col, sep=""))
            drawLabelCol(col.lab.c, j)
        }
        else if (inherits(x, "metabin")){
          if (rightcols[i]=="col.n.c" & just.cols=="right")
            drawLabelCol(col.lab.c, j)
          else if (rightcols[i]=="col.event.c" & just.cols %in% c("left", "center"))
            drawLabelCol(col.lab.c, j)
        }
        else if (inherits(x, "metacont")){
          if (rightcols[i]=="col.sd.c" & just.cols=="right")
            drawLabelCol(col.lab.c, j)
          else if (rightcols[i]=="col.mean.c" & just.cols %in% c("left", "center"))
            drawLabelCol(col.lab.c, j)
        }
        else if (inherits(x, "metainc")){
          if (rightcols[i]=="col.time.c" & just.cols=="right")
            drawLabelCol(col.lab.c, j)
          else if (rightcols[i]=="col.event.c" & just.cols %in% c("left", "center"))
            drawLabelCol(col.lab.c, j)
        }
      }
      ##
      j <- j+2
    }
  }
  ##
  popViewport()
  
  
  invisible(NULL)
}
