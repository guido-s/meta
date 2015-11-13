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
                        layout="meta",
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
                        lab.NA=".", lab.NA.effect="",
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
                        print.I2.ci=FALSE,
                        print.tau2=comb.fixed|comb.random,
                        print.Q=FALSE,
                        print.pval.Q=comb.fixed|comb.random,
                        hetstat=print.I2|print.tau2|print.Q|print.pval.Q,
                        overall.hetstat=overall&hetstat,
                        hetlab="Heterogeneity: ",
                        text.I2="I-squared",
                        text.tau2="tau-squared",
                        ##
                        test.overall=.settings$test.overall,
                        test.overall.fixed=comb.fixed&overall&test.overall,
                        test.overall.random=comb.random&overall&test.overall,
                        label.test.overall.fixed=paste("Test for overall effect",
                          if (comb.fixed & comb.random) " (fixed effect)", ": ", sep=""),
                        label.test.overall.random=paste("Test for overall effect",
                          if (comb.fixed & comb.random) " (random effects)", ": ", sep=""),
                        ##
                        test.subgroup=.settings$test.subgroup,
                        test.subgroup.fixed=if (missing(test.subgroup)) FALSE else test.subgroup,
                        test.subgroup.random=if (missing(test.subgroup)) !is.null(x$byvar)&comb.random&test.subgroup else test.subgroup,
                        print.Q.subgroup=print.Q,
                        label.test.subgroup.fixed="Test for subgroup differences (fixed effect): ",
                        label.test.subgroup.random=paste("Test for subgroup differences",
                          if (test.subgroup.fixed | comb.fixed) " (random effects)", ": ", sep=""),
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
                        fs.test.overall=fs.hetstat,
                        fs.test.subgroup=fs.hetstat,
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
                        ff.test.overall=ff.hetstat,
                        ff.test.subgroup=ff.hetstat,
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
                        digits=2,
                        digits.se=4,
                        digits.tau2=4,
                        digits.pval=4,
                        digits.pval.Q=digits.pval,
                        digits.Q=1,
                        digits.I2=1,
                        ...){
  
  
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
  layout <- setchar(layout, c("meta", "revman5"))
  chkchar(lab.NA)
  chkchar(lab.NA.effect)
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
  chklogical(print.I2.ci)
  chklogical(print.tau2)
  chklogical(print.Q)
  chklogical(print.pval.Q)
  chklogical(hetstat)
  chklogical(overall.hetstat)
  chklogical(test.overall.fixed)
  chklogical(test.overall.random)
  chklogical(test.subgroup.fixed)
  chklogical(test.subgroup.random)
  chklogical(print.Q.subgroup)
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
  chknumeric(digits, single=TRUE)
  chknumeric(digits.se, single=TRUE)
  chknumeric(digits.tau2, single=TRUE)
  chknumeric(digits.pval, single=TRUE)
  chknumeric(digits.pval.Q, single=TRUE)
  chknumeric(digits.Q, single=TRUE)
  chknumeric(digits.I2, single=TRUE)
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
    test.overall.fixed   <- FALSE
    test.overall.random  <- FALSE
    test.subgroup.fixed  <- FALSE
    test.subgroup.random <- FALSE
  }
  ##
  if (metaprop){
    test.overall.fixed  <- FALSE
    test.overall.random <- FALSE
  }
  ##
  prediction <- prediction & x$k >= 3
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
  if (layout=="revman5")
    rsel <- FALSE
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
    dataset2 <- as.data.frame(x)
    ##
    if (is.null(x$data))
      dataset1 <- dataset2
    else
      dataset1 <- x$data
    ##
    if (!is.null(x$subset))
      dataset1 <- dataset1[x$subset,]
    ##
    ## Check whether additional variables are
    ## part of meta-object
    ##
    for (i in colnames.new)
      if (length(dataset1[[i]]) == 0 & length(dataset2[[i]]) == 0)
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
      ##
      if ( (metacor|metaprop) & any(rightcols.new=="n"))
        rightlabs.new[rightlabs.new=="n"] <- "Total"
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
      ##
      if ( (metacor|metaprop) & any(leftcols.new=="n"))
        leftlabs.new[leftlabs.new=="n"] <- "Total"
    }
  }
  ##
  ## Default set of columns if argument leftcols and/or
  ## rightcols not specified
  ##
  if (is.null(leftcols)){
    if (layout=="meta"){
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
      ##
      if (metainf.metacum)
        leftcols <- "studlab"
    }
    else if (layout=="revman5"){
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
      ##
      if (metainf.metacum)
        leftcols <- "studlab"
      ##
      if (!metainf.metacum){
        if (comb.fixed & overall)
          leftcols <- c(leftcols, "w.fixed")
        if (comb.random & overall)
          leftcols <- c(leftcols, "w.random")
      }
      ##
      leftcols <- c(leftcols, "effect", "ci")
    }
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
  ##
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
    if (newcols){
      dataset1 <- dataset1[o,]
      dataset2 <- dataset2[o,]
    }
  }
  ##
  if (by){
    n.by <- length(bylevs)
    sel.by.fixed  <- 3 + 0*n.by + 1:n.by
    sel.by.random <- 3 + 1*n.by + 1:n.by
    sel.by.het    <- 3 + 2*n.by + 1:n.by
    sel.by <- c(sel.by.fixed, sel.by.random, sel.by.het)
 }
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
    lowI2 <- NA
    uppI2 <- NA
    ##
    Q.b.fixed  <- NA
    Q.b.random <- NA
    df.Q.b     <- NA
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
    lowI2 <- x$lower.I2
    uppI2 <- x$upper.I2
    ##
    Q.b.fixed  <- x$Q.b.fixed
    Q.b.random <- x$Q.b.random
    df.Q.b     <- x$df.Q.b
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
                               round(100*I2, digits.I2), "%",
                               sep="")
      if (print.I2.ci & x$k>2)
        hetstat.overall <- paste(hetstat.overall,
                                 " ",
                                 p.ci(paste(round(100*lowI2, digits.I2), "%", sep=""),
                                      paste(round(100*uppI2, digits.I2), "%", sep="")),
                                 sep="")
      dummy <- TRUE
    }
    ##
    if (print.tau2){
      hetstat.overall <- paste(hetstat.overall,
                               if (dummy) ", ",
                               if (tau2==0) paste(text.tau2, "=0", sep="")
                               else format.tau(tau2, noblanks=TRUE,
                                               lab=TRUE, labval=text.tau2,
                                               digits=digits.tau2),
                               sep="")
      dummy <- TRUE
    }
    ##
    if (print.Q){
      hetstat.overall <- paste(hetstat.overall,
                               if (dummy) ", ",
                               "Q=", round(Q, digits.Q),
                               ", df=", df,
                               sep="")
      dummy <- TRUE
    }
    ##
    if (print.pval.Q){
      hetstat.overall <- paste(hetstat.overall,
                               if (dummy) ", ",
                               format.p(1-pchisq(Q, df),
                                        lab=TRUE, noblanks=TRUE,
                                        digits=digits.pval.Q),
                               sep="")
    }
  }
  else
    hetstat.overall <- ""
  ##
  ## Text of test for subgroup differences
  ##
  if (!by){
    test.subgroup.fixed  <- FALSE
    test.subgroup.random <- FALSE
    Q.b.fixed  <- NA
    Q.b.random <- NA
    df.Q.b     <- NA
  }
  ##
  ## Label of test for overall effect
  ##
  pvals.overall <- format.p(c(x$pval.fixed, x$pval.random),
                            lab=TRUE, noblanks=TRUE,
                            digits=digits.pval)
  ##
  if (test.overall.fixed)
    text.overall.fixed  <- paste(label.test.overall.fixed, pvals.overall[1])
  else
    text.overall.fixed <- ""
  ##
  if (test.overall.random)
    text.overall.random  <- paste(label.test.overall.random, pvals.overall[2])
  else
    text.overall.random <- ""
  ##
  ##
  ## Label of test for subgroup differences
  ##
  if (!test.subgroup.fixed)
    label.test.subgroup.fixed <- ""
  if (!test.subgroup.random)
    label.test.subgroup.random <- ""
  ##
  dummy.het <- FALSE
  ##
  Q.bs <- c(Q.b.fixed, Q.b.random)
  Q.bs.format <- gsub(" ", "", format(round(Q.bs, digits.Q)))
  pval.Qbs <- format.p(1-pchisq(Q.bs, df.Q.b), lab=TRUE, noblanks=TRUE,
                       digits=digits.pval.Q)
  ##
  if (print.Q.subgroup){
    label.test.subgroup.fixed  <- paste(label.test.subgroup.fixed,
                                        "Q=", Q.bs.format[1],
                                        ", df=", df.Q.b, sep="")
    label.test.subgroup.random <- paste(label.test.subgroup.random,
                                        "Q=", Q.bs.format[2],
                                        ", df=", df.Q.b, sep="")
    dummy.het <- TRUE
  }
  ##
  text.subgroup.fixed <- paste(label.test.subgroup.fixed,
                               if (dummy.het) ", ",
                               pval.Qbs[1],
                               sep="")
  ##
  text.subgroup.random <- paste(label.test.subgroup.random,
                                if (dummy.het) ", ",
                                pval.Qbs[2],
                                sep="")
  
  
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
    lowI2.w    <- x$lower.I2.w[o]
    uppI2.w    <- x$upper.I2.w[o]
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
    lowI2.w    <- lowI2.w[sel]
    uppI2.w    <- uppI2.w[sel]
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
                           round(100*I2.w, digits.I2), "%",
                           sep="")
        if (print.I2.ci)
          hetstat.w <- paste(hetstat.w,
                             ifelse(k.w>2,
                                    paste(" ",
                                          p.ci(paste(round(100*lowI2.w, digits.I2), "%", sep=""),
                                               paste(round(100*uppI2.w, digits.I2), "%", sep="")),
                                          sep=""),
                                    ""),
                             sep="")
        dummy <- TRUE
      }
      ##
      if (print.tau2){
        hetstat.w <- paste(hetstat.w,
                           if (dummy) ", ",
                           ifelse(tau.w==0, paste(text.tau2, "=0", sep=""),
                                  format.tau(tau.w^2, noblanks=TRUE,
                                             lab=TRUE, labval=text.tau2,
                                             digits=digits.tau2)),
                           sep="")
        dummy <- TRUE
      }
      ##
      if (print.Q){
        hetstat.w <- paste(hetstat.w,
                           if (dummy) ", ",
                           "Q=", round(Q.w, digits.Q),
                           ", df=", k.w-1,
                           sep="")
        dummy <- TRUE
      }
      ##
      if (print.pval.Q){
        hetstat.w <- paste(hetstat.w,
                           if (dummy) ", ",
                           format.p(1-pchisq(Q.w, k.w-1),
                                    lab=TRUE, noblanks=TRUE,
                                    digits=digits.pval.Q),
                           sep="")
      }
      ##
      sel <- k.w==0 | k.w==1
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
      TE <- x$event.e/x$n.e
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
    modlabs <- c(text.fixed, text.random, text.predict,
                 hetstat.overall,
                 text.overall.fixed, text.overall.random,
                 text.subgroup.fixed, text.subgroup.random,
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
                            format(round(seTE, digits.se),
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
    sel.fixed[sel.by.random] <- TRUE
    sel.random[sel.by.fixed] <- TRUE
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
    modlabs <- c(text.fixed, text.random, text.predict,
                 hetstat.overall,
                 text.overall.fixed, text.overall.random,
                 text.subgroup.fixed, text.subgroup.random, studlab)
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
                            format(round(seTE, digits.se),
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
    effect.format <- ifelse(is.na(TEs), lab.NA.effect,
                            format(round(exp(TEs), digits), scientific=FALSE))
    ci.format <- ifelse(is.na(lowTEs) | is.na(uppTEs), lab.NA.effect,
                        p.ci(format(round(exp(lowTEs), digits), scientific=FALSE),
                             format(round(exp(uppTEs), digits), scientific=FALSE)))
  }
  else{
    effect.format <- ifelse(is.na(TEs), lab.NA.effect,
                            format(round(TEs, digits), scientific=FALSE))
    ci.format <- ifelse(is.na(lowTEs) | is.na(uppTEs), lab.NA.effect,
                        p.ci(format(round(lowTEs, digits), scientific=FALSE),
                             format(round(uppTEs, digits), scientific=FALSE)))
  }
  ##
  effect.format[3] <- ""
  if (!prediction)
    ci.format[3] <- ""
  if (by){
    effect.format[sel.by.het] <- ""
    ci.format[sel.by.het] <- ""
  }
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
    w.fixed.format[sel.by.het] <- ""
  w.random.format[sel.random] <- "--"
  if (by)
    w.random.format[sel.by.het] <- ""
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
  Ne.format <- ifelse(is.na(Ne), lab.NA, format(Ne, scientific=FALSE))
  Nc.format <- ifelse(is.na(Nc), lab.NA, format(Nc, scientific=FALSE))
  Ee.format <- ifelse(is.na(Ee), lab.NA, format(Ee, scientific=FALSE))
  Ec.format <- ifelse(is.na(Ec), lab.NA, format(Ec, scientific=FALSE))
  Te.format <- ifelse(is.na(Te), lab.NA, format(Te, scientific=FALSE))
  Tc.format <- ifelse(is.na(Tc), lab.NA, format(Tc, scientific=FALSE))
  ##
  ## Do not print treatment estimates for summary results
  ## (this information is provided in column 'effect')
  ##
  Te.format[1:3] <- Tc.format[1:3] <- ""
  ##
  ## Print nothing in line with prediction interval
  ##
  Ne.format[3] <- Nc.format[3] <- Ee.format[3] <- Ec.format[3] <- ""
  ##
  if (by){
    ##
    ## Do not print treatment estimates for summary results
    ## (this information is provided in column 'effect')
    ##
    Te.format[sel.by] <- Tc.format[sel.by] <- ""
    ##
    ## Print nothing in lines with heterogeneity results for subgroups
    ##
    Ne.format[sel.by.het] <- Nc.format[sel.by.het] <- ""
    Ee.format[sel.by.het] <- Ec.format[sel.by.het] <- ""
  }
  ##
  if (comb.fixed & comb.random){
    ##
    ## Print nothing in lines with results for random effects model
    ##
    Ne.format[2] <- Nc.format[2] <- Ee.format[2] <- Ec.format[2] <- ""
    ##
    if (by){
      Ne.format[sel.by.random] <- Nc.format[sel.by.random] <- ""
      Ee.format[sel.by.random] <- Ec.format[sel.by.random] <- ""
    }
  }
  ##
  ## Only print total number of events if pooled.events is TRUE
  ##
  if (!pooled.events){
    Ee.format[1:2] <- Ec.format[1:2] <- ""
    ##
    if (by){
      Ee.format[sel.by.fixed]  <- Ec.format[sel.by.fixed] <- ""
      Ee.format[sel.by.random] <- Ec.format[sel.by.random] <- ""
    }
  }
  ##
  ## Mean and standard deviation
  ##
  if (by){
    Me <- c(NA, NA, NA, rep(NA, 3*n.by), x$mean.e)
    Mc <- c(NA, NA, NA, rep(NA, 3*n.by), x$mean.c)
    Se <- c(NA, NA, NA, rep(NA, 3*n.by), x$sd.e)
    Sc <- c(NA, NA, NA, rep(NA, 3*n.by), x$sd.c)
  }
  else{
    Me <- c(NA, NA, NA, x$mean.e)
    Mc <- c(NA, NA, NA, x$mean.c)
    Se <- c(NA, NA, NA, x$sd.e)
    Sc <- c(NA, NA, NA, x$sd.c)
  }
  ##
  Me.format <- ifelse(is.na(Me), lab.NA, format(Me, scientific=FALSE))
  Mc.format <- ifelse(is.na(Mc), lab.NA, format(Mc, scientific=FALSE))
  Se.format <- ifelse(is.na(Se), lab.NA, format(Se, scientific=FALSE))
  Sc.format <- ifelse(is.na(Sc), lab.NA, format(Sc, scientific=FALSE))
  ##
  ## Print nothing for lines with summary results
  ##
  Me.format[1:3] <- Mc.format[1:3] <- Se.format[1:3] <- Sc.format[1:3] <- ""
  ##
  if (by){
    Me.format[sel.by] <- Mc.format[sel.by] <- ""
    Se.format[sel.by] <- Sc.format[sel.by] <- ""
  }
  ##
  ## Correlation
  ##
  if (by)
    cor <- c(NA, NA, NA, rep(NA, 3*n.by), x$cor)
  else
    cor <- c(NA, NA, NA, x$cor)
  ##
  cor.format <- ifelse(is.na(cor), lab.NA, format(cor, scientific=FALSE))
  ##
  ## Print nothing for lines with summary results
  ##
  cor.format[1:3] <- ""
  ##
  if (by)
    cor.format[sel.by] <- ""
  ##
  ##
  ## y-axis:
  ##
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
  ##
  ## x-axis:
  ##
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
  yNext <- max.yTE + 2
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
  yTE.fixed  <- NA
  yTE.random <- NA
  yPredict   <- NA
  yHetstat <- NA
  yOverall.fixed  <- NA
  yOverall.random <- NA
  ySubgroup.fixed  <- NA
  ySubgroup.random <- NA
  ##
  if (comb.fixed & comb.random & overall){
    yTE.fixed  <- yNext
    yTE.random <- yNext + 1
    yNext      <- yNext + 2
  }
  ##
  else if (comb.fixed & !comb.random & overall){
    yTE.fixed <- yNext
    yNext     <- yNext + 1
  }
  ##
  else if (!comb.fixed & comb.random & overall){
    yTE.random <- yNext
    yNext      <- yNext + 1
  }
  ##
  else if (!comb.fixed & !comb.random & pooled.totals & overall){
    yTE.fixed  <- yNext
    yNext      <- yNext + 1
    if (missing(text.fixed))
      text.fixed <- "Overall"
  }
  ##
  if (prediction){
    yPredict <- yNext
    yNext    <- yNext + 1
  }
  ##
  if (overall.hetstat){
    yHetstat <- yNext
    yNext    <- yNext + 1
  }
  ##
  if (test.overall.fixed){
    yOverall.fixed <- yNext
    yNext          <- yNext + 1
  }
  ##
  if (test.overall.random){
    yOverall.random <- yNext
    yNext           <- yNext + 1
  }
  ##
  if (test.subgroup.fixed){
    ySubgroup.fixed <- yNext
    yNext           <- yNext + 1
  }
  ##
  if (test.subgroup.random){
    ySubgroup.random <- yNext
    yNext            <- yNext + 1
  }
  ##
  if (!comb.fixed & !pooled.totals) text.fixed <- ""
  if (!comb.random) text.random <- ""
  if (!prediction) text.predict <- ""
  ##  
  yTE        <- yHead + yTE + addspace
  yTE.fixed  <- yHead + yTE.fixed + addspace
  yTE.random <- yHead + yTE.random + addspace
  yPredict   <- yHead + yPredict + addspace
  yHetstat <- yHead + yHetstat + addspace
  yOverall.fixed  <- yHead + yOverall.fixed + addspace
  yOverall.random <- yHead + yOverall.random + addspace
  ySubgroup.fixed  <- yHead + ySubgroup.fixed + addspace
  ySubgroup.random <- yHead + ySubgroup.random + addspace
  ##
  if (by){
    yBylab <- yHead + yBylab + addspace
    yTE.w  <- yHead + yTE.w + addspace
  }
  ##  
  if (by){
    yLab <- c(yHead,
              yTE.fixed, yTE.random, yPredict,
              yHetstat,
              yOverall.fixed, yOverall.random,
              ySubgroup.fixed, ySubgroup.random,
              yBylab, yTE.w,
              yTE)
    ##
    yS <- c(yHead, yTE.fixed, yTE.random, yPredict, yTE.w, yTE)
  }
  else{
    yLab <- c(yHead, yTE.fixed, yTE.random, yPredict,
              yHetstat,
              yOverall.fixed, yOverall.random,
              ySubgroup.fixed, ySubgroup.random,
              yTE)
    yS   <- c(yHead, yTE.fixed, yTE.random, yPredict, yTE)
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
  ## Prediction interval:
  ##
  col.studlab$labels[[4]] <- textGrob(text.predict,
                                      x=xpos.studlab, just=just.studlab,
                                      gp=gpar(
                                        fontsize=fs.predict.labels,
                                        fontface=ff.predict.labels)
                                      )
  ##
  ## Heterogeneity statistics:
  ##
  col.studlab$labels[[5]] <- textGrob(hetstat.overall,
                                      x=xpos.studlab, just=just.studlab,
                                      gp=gpar(
                                        fontsize=fs.hetstat,
                                        fontface=ff.hetstat)
                                      )
  ##
  ## Test for overall effect (fixed effect model):
  ##
  col.studlab$labels[[6]] <- textGrob(text.overall.fixed,
                                      x=xpos.studlab, just=just.studlab,
                                      gp=gpar(
                                        fontsize=fs.test.overall,
                                        fontface=ff.test.overall)
                                      )
  ##
  ## Test for overall effect (random effects model):
  ##
  col.studlab$labels[[7]] <- textGrob(text.overall.random,
                                      x=xpos.studlab, just=just.studlab,
                                      gp=gpar(
                                        fontsize=fs.test.overall,
                                        fontface=ff.test.overall)
                                      )
  ##
  ## Test for subgroup differences (fixed effect model):
  ##
  col.studlab$labels[[8]] <- textGrob(text.subgroup.fixed,
                                      x=xpos.studlab, just=just.studlab,
                                      gp=gpar(
                                        fontsize=fs.test.subgroup,
                                        fontface=ff.test.subgroup)
                                      )
  ##
  ## Test for subgroup differences (random effects model):
  ##
  col.studlab$labels[[9]] <- textGrob(text.subgroup.random,
                                      x=xpos.studlab, just=just.studlab,
                                      gp=gpar(
                                        fontsize=fs.test.subgroup,
                                        fontface=ff.test.subgroup)
                                      )
  ##
  n.summaries <- 9
  ##
  if (by){
    for (i in 1:n.by){
      ##
      ## Subgroup labels:
      ##
      col.studlab$labels[[n.summaries+i]] <- textGrob(bylab[i],
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
      col.studlab$labels[[n.summaries+n.by+i]] <- textGrob(text.fixed.w[i],
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
      col.studlab$labels[[n.summaries+2*n.by+i]] <- textGrob(text.random.w[i],
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
      col.studlab$labels[[n.summaries+3*n.by+i]] <- textGrob(hetstat.w[i],
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
        if (length(dataset1[[rightcols.new[i]]])!=0)
          tmp.r <- dataset1[[rightcols.new[i]]]
        else if (length(dataset2[[rightcols.new[i]]])!=0)
          tmp.r <- dataset2[[rightcols.new[i]]]
        if (is.factor(tmp.r))
          tmp.r <- as.character(tmp.r)
        tmp.r <- ifelse(is.na(tmp.r), lab.NA, tmp.r)
        cols[[tname]] <- formatcol(rightlabs.new[i],
                                   c("", "", "",
                                     rep("", length(TE.w)),
                                     tmp.r),
                                   yS,
                                   just=just.addcols)
      }
      for (i in seq(along=leftcols.new)){
        tname <- paste("col.", leftcols.new[i], sep="")
        if (length(dataset1[[leftcols.new[i]]])!=0)
          tmp.l <- dataset1[[leftcols.new[i]]]        
        else if (length(dataset2[[leftcols.new[i]]])!=0)
          tmp.l <- dataset2[[leftcols.new[i]]]
        if (is.factor(tmp.l))
          tmp.l <- as.character(tmp.l)
        tmp.l <- ifelse(is.na(tmp.l), lab.NA, tmp.l)
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
        if (length(dataset1[[rightcols.new[i]]])!=0)
          tmp.r <- dataset1[[rightcols.new[i]]]
        else if (length(dataset2[[rightcols.new[i]]])!=0)
          tmp.r <- dataset2[[rightcols.new[i]]]
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
        if (length(dataset1[[leftcols.new[i]]])!=0)
          tmp.l <- dataset1[[leftcols.new[i]]]        
        else if (length(dataset2[[leftcols.new[i]]])!=0)
          tmp.l <- dataset2[[leftcols.new[i]]]
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
        if (!is.null(lty.random) & !is.na(TE.random))
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
        x1 <- unit.c(wcalc(cols[[leftcols[i]]]$labels[-(c(4:n.summaries, n.summaries+3*n.by+1:n.by))]))
      else
        x1 <- unit.c(wcalc(cols[[leftcols[i]]]$labels))
    else
      if (leftcols[[i]]=="col.studlab")
        x1 <- unit.c(x1,
                     colgap.left,
                     wcalc(cols[[leftcols[i]]]$labels[-(c(4:n.summaries, n.summaries+3*n.by+1:n.by))]))
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
  if (by){
    addline <- addspace*(!any(c(test.overall.fixed, test.overall.random,
                                hetstat,
                                test.subgroup.fixed, test.subgroup.random)))
    ##
    nrow <- max(addline+c(yTE, yTE.fixed, yTE.random, yPredict,
                          yHetstat,
                          yOverall.fixed, yOverall.random,
                          ySubgroup.fixed, ySubgroup.random, yTE.w), na.rm=TRUE)
  }
  else{
    addline <- addspace*(!any(c(test.overall.fixed, test.overall.random, hetstat)))
    ##
    nrow <- max(addline+c(yTE, yTE.fixed, yTE.random, yPredict,
                          yHetstat,
                          yOverall.fixed, yOverall.random), na.rm=TRUE)
  }
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
