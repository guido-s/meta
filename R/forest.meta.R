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
                        plotwidth=unit(6, "cm"),
                        colgap=unit(2, "mm"),
                        colgap.left=colgap,
                        colgap.right=colgap,
                        colgap.forest=colgap,
                        colgap.forest.left=colgap.forest,
                        colgap.forest.right=colgap.forest,
                        ##
                        just="center",
                        just.studlab="left",
                        ##
                        addspace=TRUE,
                        ##
                        new=TRUE,
                        ##
                        backtransf=x$backtransf,
                        digits=2, ...){
  
  
  if (!inherits(x, "meta"))
    stop("Argument 'x' must be an object of class \"meta\"")
  
  
  ##
  ## Definition of auxiliary functions
  ##
  wcalc <- function(x)
    max(unit(rep(1, length(x)), "grobwidth", x))
  ##
  ##
  testchar <- function(x)
    if (length(x)!=1 || !is.character(x))
      stop(paste("Parameter '", deparse(substitute(x)),
                 "' must be a character string", sep=""))
  ##
  ##
  repchar <- function(x, len){
    if (length(x)==1)
      res <- rep(x, len)
    else if (length(x)!=len)
      stop(paste("Parameter '", deparse(substitute(x)),
                 "' has different length than number of studies in meta-analysis",
                 sep=""))
    else
      res <- x
    res
  }
  ##
  ##
  formatcol <- function(x, y, rows, just="right"){
    ##
    ##
    ##
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
      ##if ((!is.na(size) && size > 0) |
      ##    (!is.na(size) && (size == 0 & print.infinite))){
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
        ##grid.circle(x=unit(eff, "native"),
        ##            r=unit(0.15, "snpc"),
        ##            gp=gpar(fill=col, col=col))
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
      ##if ((!is.na(size) && size > 0) |
      ##    (!is.na(size) && (size == 0 & print.infinite))){
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
        ##if (convertX(unit(upp, "native"), "npc", valueOnly=TRUE) > 1)
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
                   ##y=unit(0.5 + c(0, 0.25*size, 0, -0.25*size), "npc"),
                   y=unit(0.5 + c(0, 0.3*size, 0, -0.3*size), "npc"),
                   gp=gpar(fill=col.diamond, col=col.diamond.lines))
  }
  ##
  ##
  drawPredictionCI <- function(low, upp, size, min, max,
                               col.predict, col.predict.lines) {
    ##
    ## Function to draw a prediction interval
    ##
    ##if (!(is.na(low) | is.na(upp)) &&
    ##    ((min <= low & low <= max) |
    ##     (min <= upp & upp <= max)))
    if (!(is.na(low) | is.na(upp))){
      ## Plot prediction interval only within plotting range
      ##if (low < min) low <- min
      ##if (upp > max) upp <- max
      if ((min <= low & low <= max) |
          (min <= upp & upp <= max))
        grid.polygon(x=unit(c(low, low, upp, upp), "native"),
                     y=unit(0.5 + c(-0.1*size, 0.1*size, 0.1*size, -0.1*size), "npc"),
                     gp=gpar(fill=col.predict, col=col.predict.lines))
    }
  }
  ##
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
    ##
    ## 
    ##
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
    ##
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
    ##    grid.text("Favours xyz   ",
    ##              x=unit(ref, "native"),
    ##              ##y=-0.2,
    ##              y=unit(-3.35, "lines"),
    ##              just=c("right", "top"),
    ##              gp=gpar(fontsize=fs.xlab, fontface=ff.xlab, col="darkgrey"))
    ##    ##
    ##    grid.text("   Favours zyx",
    ##              x=unit(ref, "native"),
    ##              ##y=-0.2,
    ##              y=unit(-3.35, "lines"),
    ##              just=c("left", "top"),
    ##              gp=gpar(fontsize=fs.xlab, fontface=ff.xlab, col="darkgrey"))
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
  
  
  if (new)
    grid.newpage()
  
  
  ## Upgrade meta objects created with older versions of meta
  ##
  if (!(!is.null(x$version) &&
        as.numeric(unlist(strsplit(x$version, "-"))[1]) >= 3.8))
    x <- update(x, warn=FALSE)
  
  
  metainf.metacum <- inherits(x, "metainf") | inherits(x, "metacum")
  
  
  if (metainf.metacum){
    hetstat <- FALSE
    prediction <- FALSE
  }
  
  x.name <- deparse(substitute(x))
  
  
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
  byvar <- x$byvar
  level <- x$level
  level.comb <- x$level.comb
  level.predict <- x$level.predict
  
  
  if (length(comb.fixed)==0)
    comb.fixed <- FALSE
  ##
  if (length(comb.random)==0)
    comb.random <- FALSE
  ##
  if (length(prediction)==0)
    prediction <- FALSE
  ##
  if (length(print.byvar)==0)
    print.byvar <- TRUE
  ##
  if (length(lab.e)==0)
    lab.e <- "Experimental"
  ##
  if (length(lab.c)==0)
    lab.c <- "Control"
  
  
  prediction <- prediction & comb.random & x$k>=3
  
  
  if (length(level)==0){
    warning("level set to 0.95")
    level <- 0.95
  }
  ##
  if (length(level.comb)==0){
    if (comb.fixed | comb.random)
      warning("level.comb set to 0.95")
    level.comb <- 0.95
  }
  ##
  if (length(level.predict)==0){
    if (prediction)
      warning("level.predict set to 0.95")
    level.predict <- 0.95
  }
  
  
  ##
  ## Check for levels of confidence interval
  ##
  if (!is.numeric(level) | length(level)!=1)
    stop("list object 'level' must be a numeric of length 1")
  if (level <= 0 | level >= 1)
    stop("list object 'level': no valid level for confidence interval")
  ##
  if (!is.numeric(level.comb) | length(level.comb)!=1)
    stop("list object 'level.comb' must be a numeric of length 1")
  if (level.comb <= 0 | level.comb >= 1)
    stop("list object 'level.comb': no valid level for confidence interval")
  ##
  if (!is.numeric(level.predict) | length(level.predict)!=1)
    stop("list object 'level.predict' must be a numeric of length 1")
  if (level.predict <= 0 | level.predict >= 1)
    stop("list object 'level.predict': no valid level for confidence interval")
  
  
  if (is.logical(leftcols)){
    warning("Logical value not possible for parameter 'leftcols', set to 'NULL'.")
    leftcols <- NULL
  }
  
  
  if (missing(weight))
    weight <- ifelse(comb.random & !comb.fixed, "random", "fixed")
  
  
  notmiss.xlim <- !missing(xlim)
  
  
  ijs <- charmatch(tolower(just.studlab),
                   c("right", "center", "left"), nomatch = NA)
  ##
  if(is.na(ijs))
    stop("Argument 'just.studlab' should be \"right\", \"center\", or \"left\"")
  ##
  just.studlab <- c("right", "center", "left")[ijs]
  ##
  if (just.studlab=="left")
    xpos.studlab <- 0
  else if (just.studlab=="center")
      xpos.studlab <- 0.5
  else if (just.studlab=="right")
    xpos.studlab <- 1
  
  
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
  
  
  if (!backtransf & pscale!=1){
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  
  
  if (is.null(xlab))
    xlab <- xlab(sm, backtransf)
  
  if (is.null(smlab))
    smlab <- xlab(sm, backtransf)
  
  if (is.null(label.right))
    label.right <- ""
  if (is.null(label.left))
    label.left <- ""
  
  
  testchar(col.diamond)
  testchar(col.diamond.fixed)
  testchar(col.diamond.random)
  testchar(col.diamond.lines)
  testchar(col.diamond.fixed.lines)
  testchar(col.diamond.random.lines)
  testchar(col.predict)
  testchar(col.predict.lines)
  
  
  iweight <- charmatch(tolower(weight),
                       c("same", "fixed", "random"), nomatch = NA)
  ##
  if(is.na(iweight))
    stop("weight should be \"same\", \"fixed\", or \"random\"")
  ##
  weight <- c("same", "fixed", "random")[iweight]
  
  
  ijust <- charmatch(tolower(just),
                       c("right", "center", "left"), nomatch = NA)
  ##
  if(is.na(ijust))
    stop("just should be \"right\", \"center\", or \"left\"")
  ##
  just <- c("right", "center", "left")[ijust]
  
  
  if (inherits(x, "metaprop")){
    x$event.e <- x$event
    x$n.e <- x$n
  }
  ##
  if (inherits(x, "metacor")){
    x$n.e <- x$n
  }
  
  
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
  
  
  ##
  ## Determine whether to print columns on right and/or left side
  ## of forest plot
  ##
  rsel <- !(is.logical(rightcols) && length(rightcols)==1 && !rightcols)
  ##
  if (!rsel)
    rightcols <- NULL
  ##
  ## lsel <- !(is.logical(leftcols) && length(leftcols)==1 && !leftcols)
  ## ##
  ## if (!lsel)
  ##   leftcols <- NULL
  
  
  ##
  ## Check for duplicate columns
  ##
  if (length(c(rightcols, leftcols))>0 &&
      any(duplicated(c(rightcols, leftcols))))
    stop("Duplicate entries in 'leftcols' and 'rightcols'")
  
  
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
    for (i in colnames.new){
      if (length(x[[i]]) == 0 & length(dataset[[i]]) == 0)
        stop(paste("variable '", i,
                   "' not available in '",
                   x.name, "'", sep=""))
    }
    ##
    rightcols.new <- rightcols[! rightcols %in% colnames]
    leftcols.new  <- leftcols[! leftcols %in% colnames]
    ##
    ## Determine label for new columns
    ## 1. Use column name as label if no label is given
    ##    parameter right|left|labs
    ## 2. Otherwise use corresponding entry from
    ##    parameter right|left|labs
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
        stop("Too few labels defined for parameter 'rightcols'")
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
        stop("Too few labels defined for parameter 'leftcols'")
    }
  }
  
  
  ##
  ## Default set of columns if parameter leftcols and/or
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
    if (inherits(x, "metaprop"))
      leftcols <- c("studlab",
                    "event.e", "n.e")
    ##
    if (inherits(x, "metacor"))
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
    ##    else{
    ##      if (x$pooled=="fixed" & overall)
    ##        rightcols <- c(rightcols, "w.fixed")
    ##      if (x$pooled=="random" & overall)
    ##        rightcols <- c(rightcols, "w.random")
    ##    }
  }
  
  
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
  ## total number of studies to plot (*not* number of studies combined)
  ##
  k.all <- length(x$TE)


  if (allstudies)
    n.stud <- k.all # all studies
  else
    n.stud <- x$k   # number of studies combined in meta-analysis
  
  
  by <- length(byvar)>0
  sort <- !missing(sortvar)
  
  
  mf <- match.call()
  ##
  error <- try(sortvar <- eval(mf[[match("sortvar", names(mf))]],
                               as.data.frame(x, stringsAsFactors=FALSE),
                               enclos = sys.frame(sys.parent())),
               silent=TRUE)
  ##
  if (class(error)=="try-error"){
    xd <- x$data
    sortvar <- eval(mf[[match("sortvar", names(mf))]],
                    xd, enclos = NULL)
    ##
    if (!is.null(x$data$.subset))
      sortvar <- sortvar[x$data$.subset]
  }
  
  
  if (by & metainf.metacum)
    stop("Use of 'byvar' not possible for 'metainf' or 'metacum' object.") 
  ##  
  if (!by) byvar <- rep(1, k.all)
  if (!sort) sortvar <- rep(1, k.all)
  ##
  if (sort & length(sortvar) != k.all)
    stop("'x' and 'sortvar' have different length")
  ##
  if (by & length(byvar) != k.all)
    stop("'x' and 'byvar' have different length")
  ##
  if (by & any(is.na(byvar)))
    stop("Missing values in 'byvar'")
  
  
  if (length(studlab) == 1 & is.logical(studlab))
    if (studlab == FALSE){
      studlab <- rep("", k.all)
      labs[["lab.studlab"]] <- ""
    }
    else studlab <- x$studlab
  ##
  if (length(studlab) != k.all)
    stop("'x' and 'studlab' have different length")
  
  
  col.i               <- repchar(col.i, length(x$TE))
  col.i.inside.square <- repchar(col.i.inside.square, length(x$TE))
  col.square          <- repchar(col.square, length(x$TE))
  col.square.lines    <- repchar(col.square.lines, length(x$TE))
  
  
  if (allstudies)
    sel <- 1:length(x$TE)
  else
    sel <- !is.na(x$TE)
  ##
  if (n.stud != sum(sel>0))
    warning("n.stud != sum(sel)")
  ##
  ##
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

  
  if (sort | by){
    if (bysort){
      ##
      o <- order(byvar, sortvar)
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
    }
    else{
      ##
      byvar.factor <- factor(byvar, levels=unique(byvar))
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
    }
  }
  
  
  by.levs <- unique(byvar)
  ##
  if (by)
    n.by <- length(by.levs)
  else
    n.by <- 0
  
  
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
    if (inherits(x, "metaprop") & !backtransf){
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
  
  
  if (overall.hetstat){
    dummy <- FALSE
    ##
    hetstat.overall <- hetlab
    ##
    if (print.I2){
      hetstat.overall <- paste(hetstat.overall,
                               "I-squared=",
                               round(100*I2, 1), "%",
                               sep="")
      dummy <- TRUE
    }
    ##
    if (print.tau2){
      hetstat.overall <- paste(hetstat.overall,
                               if (dummy) ", ",
                               if (tau2==0) "tau-squared=0"
                               else format.tau(tau2, noblanks=TRUE,
                                               lab=TRUE, labval="tau-squared"),
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

  
  
  if (by){
    
    if (x$tau.common)
      x$tau.preset <- x$tau
    
    res.w <- matrix(NA, ncol=16, nrow=n.by)
    j <- 0
    ##
    for ( i in by.levs){
      j <- j+1
      sel <- byvar == i
      ##
      if (inherits(x, "metabin"))
        m.w <- metabin(x$event.e[sel], x$n.e[sel],
                       x$event.c[sel], x$n.c[sel],
                       method=x$method, sm=sm,
                       incr=x$incr, allincr=x$allincr,
                       addincr=x$addincr, allstudies=x$allstudies,
                       MH.exact=x$MH.exact, RR.cochrane=x$RR.cochrane,
                       hakn=x$hakn, method.tau=x$method.tau,
                       tau.preset=x$tau.preset, TE.tau=x$TE.tau,
                       warn=x$warn)
      ##
      if (inherits(x, "metacont"))
        m.w <- metacont(x$n.e[sel], x$mean.e[sel], x$sd.e[sel],
                        x$n.c[sel], x$mean.c[sel], x$sd.c[sel],
                        sm=sm,
                        hakn=x$hakn, method.tau=x$method.tau,
                        tau.preset=x$tau.preset, TE.tau=x$TE.tau)
      ##
      if (inherits(x, "metagen"))
        if (!is.null(x$tau.preset))
          m.w <- metagen(x$TE[sel], x$seTE[sel], sm=sm,
                         hakn=x$hakn, method.tau=x$method.tau,
                         tau.preset=x$tau.preset, TE.tau=x$TE.tau)
        else
          m.w <- metagen(x$TE[sel], x$seTE[sel], sm=sm,
                         hakn=x$hakn, method.tau=x$method.tau,
                         TE.tau=x$TE.tau)
      ##
      if (inherits(x, "metaprop")){
        m.w <- metaprop(x$event.e[sel], x$n.e[sel], sm=sm,
                        incr=x$incr, allincr=x$allincr,
                        addincr=x$addincr,
                        hakn=x$hakn, method.tau=x$method.tau,
                        tau.preset=x$tau.preset, TE.tau=x$TE.tau,
                        warn=x$warn)
        m.w$event.e <- m.w$event
        m.w$n.e <- m.w$n
      }
      ##
      if (inherits(x, "metacor")){
        m.w <- metacor(x$cor[sel], x$n[sel], sm=sm,
                       hakn=x$hakn, method.tau=x$method.tau,
                       tau.preset=x$tau.preset, TE.tau=x$TE.tau)
        m.w$n.e <- m.w$n
      }
      ##
      if (inherits(x, "metainc"))
        m.w <- metainc(x$event.e[sel], x$time.e[sel],
                       x$event.c[sel], x$time.c[sel],
                       method=x$method, sm=sm,
                       incr=x$incr, allincr=x$allincr,
                       addincr=x$addincr,
                       hakn=x$hakn, method.tau=x$method.tau,
                       tau.preset=x$tau.preset, TE.tau=x$TE.tau,
                       warn=x$warn)
      ##
      if (is.null(summary(m.w)$fixed$harmonic.mean))
        harmmean <- NA
      else
        harmmean <- summary(m.w)$fixed$harmonic.mean
      ##
      res.w[j,] <- c(m.w$TE.fixed, m.w$seTE.fixed,
                     m.w$TE.random, m.w$seTE.random,
                     m.w$Q, m.w$k, length(m.w$TE),
                     summary(m.w)$I2$TE, m.w$tau^2,
                     sum(x$w.fixed[sel]), sum(x$w.random[sel]),
                     sum(m.w$event.e), sum(m.w$n.e),
                     sum(m.w$event.c), sum(m.w$n.c),
                     harmmean)
    }
    ##
    k.w     <- res.w[,6]
    k.all.w <- res.w[,7]
    ##
    if (allstudies) sel <- 1:length(k.w)
    else sel <- k.w>0
    ##
    TE.fixed.w    <- res.w[sel,1]
    seTE.fixed.w  <- res.w[sel,2]
    TE.random.w   <- res.w[sel,3]
    seTE.random.w <- res.w[sel,4]
    ##
    Q.w           <- res.w[sel,5]
    I2.w          <- res.w[sel,8]
    tau2.w        <- res.w[sel,9]
    w.fixed.w     <- res.w[sel,10]
    w.random.w    <- res.w[sel,11]
    ##
    e.e.w         <- res.w[sel,12]
    n.e.w         <- res.w[sel,13]
    e.c.w         <- res.w[sel,14]
    n.c.w         <- res.w[sel,15]
    ##
    harmonic.mean.w <- res.w[sel,16]
    ##
    by.levs <- by.levs[sel]
    k.all.w <- k.all.w[sel]
    ##
    if (comb.fixed){
      meta.fixed.w <- ci(TE.fixed.w, seTE.fixed.w, level)
      if (sum(w.fixed.w)>0)
        w.fixed.w.p <- 100*round(w.fixed.w/sum(w.fixed.w, na.rm=TRUE), 3)
      else
        w.fixed.w.p <- w.fixed.w
    }
    else{
      meta.fixed.w <- list(TE=rep(NA, n.by),
                           lower=rep(NA, n.by),
                           upper=rep(NA, n.by))
      w.fixed.w.p <- rep(NA, n.by)
      if (missing(text.fixed.w))
        text.fixed.w <- rep("Overall", n.by)
    }
    ##
    if (comb.random){
      if (!is.null(x$hakn) && x$hakn)
        meta.random.w <- ci(TE.random.w, seTE.random.w, level, df=k.w-1)
      else
        meta.random.w <- ci(TE.random.w, seTE.random.w, level)
      ##
      if (sum(w.random.w)>0)
        w.random.w.p <- 100*round(w.random.w/sum(w.random.w, na.rm=TRUE), 3)
      else
        w.random.w.p <- w.random.w
    }
    else{
      meta.random.w <- list(TE=rep(NA, n.by),
                            lower=rep(NA, n.by),
                            upper=rep(NA, n.by))
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
                           "I-squared=",
                           round(100*I2.w, 1), "%",
                           sep="")
        dummy <- TRUE
      }
      ##
      if (print.tau2){
        hetstat.w <- paste(hetstat.w,
                           if (dummy) ", ",
                           ifelse(tau2.w==0, "tau-squared=0",
                                  format.tau(tau2.w, noblanks=TRUE,
                                             lab=TRUE, labval="tau-squared")),
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
    }
    else
      hetstat.w <- rep("", n.by)
    ##
    TE.w <- c(meta.fixed.w$TE, meta.random.w$TE, rep(NA, n.by))
    lowTE.w <- c(meta.fixed.w$lower, meta.random.w$lower, rep(NA, n.by))
    uppTE.w <- c(meta.fixed.w$upper, meta.random.w$upper, rep(NA, n.by))
    harmonic.mean.w <- c(harmonic.mean.w, harmonic.mean.w, rep(NA, n.by))
    ##
    weight.w.p <- c(w.fixed.w.p, w.random.w.p, rep(NA, n.by))
  }
  
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
    if (inherits(x, "metaprop")){
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
        npft.w <- harmonic.mean.w
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
  
  
  if (sum(x$w.fixed)>0)
    w.fixed.p <- 100*round(x$w.fixed/sum(x$w.fixed, na.rm=TRUE), 3)
  else
    w.fixed.p <- x$w.fixed
  ##
  if (sum(x$w.random)>0)
    w.random.p <- 100*round(x$w.random/sum(x$w.random, na.rm=TRUE), 3)
  else
    w.random.p <- x$w.random

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
  
  
  if (print.byvar){
    if (length(bylab)==0 || bylab=="")
      bylab <- format(by.levs)
    else
      bylab <- paste(bylab,
                     " = ",
                     format(by.levs), sep="")
  }
  else
    bylab <- format(by.levs)
  
  
  ##
  ## Format the columns on right and left
  ## hand side of the funnel plot
  ##
  if (by){
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
    w.fixeds  <- c(100, "--", "", format(c(weight.w.p, w.fixed.p), scientific=FALSE))
    w.randoms <- c("--", 100, "", format(c(weight.w.p, w.random.p), scientific=FALSE))
    ##
    sel.fixed  <- w.fixeds=="--"
    sel.random <- w.randoms=="--"
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
    w.fixeds  <- c(100, "--", "", format(w.fixed.p, scientific=FALSE))
    w.randoms <- c("--", 100, "", format(w.random.p, scientific=FALSE))
    ##
    sel.fixed <- w.fixeds=="--"
    sel.random <- w.randoms=="--"
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
  w.fixed.format  <- paste(w.fixeds, "%", sep="")
  w.random.format <- paste(w.randoms, "%", sep="")
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
    if (inherits(x, "metaprop")){
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
    ixlim <- charmatch(tolower(xlim), c("symmetric"), nomatch = NA)
    ##
    if (is.na(ixlim))
      stop("xlim should be a numeric vector (min, max) or the character string \"symmetric\"")
    ##
    symmetric <- TRUE
    xlim <- c("symmetric")[ixlim]
    ##
    if (inherits(x, "metaprop")){
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
  
  
  if (by)
    max.yTE <- max(c(yTE, yTE.w), na.rm=TRUE)
  else
    max.yTE <- max(yTE, na.rm=TRUE)
  
  
  if (!is.na(ref) & missing(xlab.pos))
    if (ref <= min(xlim) | ref >= max(xlim))
      xlab.pos <- mean(xlim)
  
  if (!is.na(ref) & missing(smlab.pos))
    if (ref <= min(xlim) | ref >= max(xlim))
      smlab.pos <- mean(xlim)
  
  if (is.null(xlab.pos) || is.na(xlab.pos))
    xlab.pos <- mean(xlim)
  
  if (is.null(smlab.pos) || is.na(smlab.pos))
    smlab.pos <- mean(xlim)
  
  
  ##
  ##
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
  
  
  ##if (!comb.fixed) text.fixed <- ""
  if (!comb.fixed & !pooled.totals) text.fixed <- ""
  if (!comb.random) text.random <- ""
  if (!prediction) text.predict <- ""
  
  
  ##
  ##
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
  
  
  ##  print(yLab)
  ##  print(yS)
  
  
  ##  cat("yHead:\n")
  ##  print(yHead)
  ##  ##
  ##  cat("yTE:\n")
  ##  print(yTE)
  ##  ##
  ##  cat("yTE.fixed:\n")
  ##  print(yTE.fixed)
  ##  ##
  ##  cat("yTE.random:\n")
  ##  print(yTE.random)
  ##  ##
  ##  if (by){
  ##    cat("yBylab:\n")
  ##    print(yBylab)
  ##    ##
  ##    cat("yTE.w:\n")
  ##    print(yTE.w)
  ##  }
  ##  ##
  ##  cat("TEs:\n")
  ##  print(round(TEs, 2))
  
  col.studlab <- list(labels=
                      lapply(as.list(c(labs[["lab.studlab"]], modlabs)),
                      ##lapply(as.list(c(lab.studlab, modlabs)),
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
  ##col.studlab$labels[[1]] <- textGrob(lab.studlab,
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
  
  
  col.effect <- formatcol(labs[["lab.effect"]], effect.format, yS)
  ##
  col.ci <- formatcol(labs[["lab.ci"]],
                      ci.format, yS)
  ##
  col.w.fixed  <- formatcol(labs[["lab.w.fixed"]], w.fixed.format, yS)
  col.w.random <- formatcol(labs[["lab.w.random"]], w.random.format, yS)
  ##
  col.TE <- formatcol(labs[["lab.TE"]], TE.format, yS)
  col.seTE <- formatcol(labs[["lab.seTE"]], seTE.format, yS)
  
  
  col.n.e <- formatcol(labs[["lab.n.e"]], Ne.format, yS)
  col.n.c <- formatcol(labs[["lab.n.c"]], Nc.format, yS)
  ##
  col.event.e <- formatcol(labs[["lab.event.e"]], Ee.format, yS)
  col.event.c <- formatcol(labs[["lab.event.c"]], Ec.format, yS)
  ##
  col.mean.e <- formatcol(labs[["lab.mean.e"]], Me.format, yS)
  col.mean.c <- formatcol(labs[["lab.mean.c"]], Mc.format, yS)
  ##
  col.sd.e <- formatcol(labs[["lab.sd.e"]], Se.format, yS)
  col.sd.c <- formatcol(labs[["lab.sd.c"]], Sc.format, yS)
  ##
  col.cor  <- formatcol(labs[["lab.cor"]], cor.format, yS)
  ##
  col.time.e <- formatcol(labs[["lab.time.e"]], Te.format, yS)
  col.time.c <- formatcol(labs[["lab.time.c"]], Tc.format, yS)
  
  
  ##  print(length(TEs))
  ##  print(length(lowTEs))
  ##  print(length(uppTEs))
  ##  print(length(yS[-1]))
  ##  ##
  ##  print(round(TEs, 2))
  ##  print(round(lowTEs, 2))
  ##  print(round(uppTEs, 2))
  ##  print(yS[-1])
  
  
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
    if (weight=="fixed"){
      tmp.w.fixeds <- w.fixeds
      tmp.w.fixeds[c(1:2, grep("NA", w.fixeds))] <- ""
      tmp.w.fixeds <- as.numeric(tmp.w.fixeds)
      ##
      information <- sqrt(tmp.w.fixeds)
    }
    if (weight=="random"){
      tmp.w.randoms <- w.randoms
      tmp.w.randoms[c(1:2, grep("NA", w.randoms))] <- ""
      tmp.w.randoms <- as.numeric(tmp.w.randoms)
      ##
      information <- sqrt(tmp.w.randoms)
    }
    information <- information/max(information, na.rm=TRUE)
    information[is.na(information)] <- 1
    ##information <- sqrt(1 / ((col.forest$upp - col.forest$eff)/1.96))
  }
  ##
  col.forest$sizes <- information
  col.forest$sizes <- col.forest$sizes * squaresize
  ##
  ## Width of column 3
  col.forestwidth <- plotwidth
  ##
  ## Range on the x-axis for column 3
  col.forest$range <- xlim # c(0, 2)
  
  
  if (by)
    nrow <- max(c(yTE, yTE.fixed, yTE.random, yTE.hetstat, yTE.predict, yTE.w), na.rm=TRUE)
  else
    nrow <- max(c(yTE, yTE.fixed, yTE.random, yTE.hetstat, yTE.predict), na.rm=TRUE)
  
  
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
                                     tmp.r[o]),
                                   yS,
                                   just=just)
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
                                     tmp.l[o]),
                                   yS,
                                   just=just)
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
                                   c("", "", "",
                                     if (sort) tmp.r[o] else tmp.r
                                     ),
                                   yS,
                                   just=just)
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
                                     if (sort) tmp.l[o] else tmp.l
                                     ),
                                   yS,
                                   just=just)
      }
    }
  }
  
  
  col.lab.e <- list(labels=list(textGrob(lab.e,
                      x=1, just="right",
                      gp=gpar(
                        fontsize=fs.heading,
                        fontface=ff.heading)
                      )),
                    rows=1)
  ##
  col.lab.c <- list(labels=list(textGrob(lab.c,
                      x=1, just="right",
                      gp=gpar(
                        fontsize=fs.heading,
                        fontface=ff.heading)
                      )),
                    rows=1)
  
  
  leftcols  <- paste("col.", leftcols, sep="")
  rightcols <- paste("col.", rightcols, sep="")
  
  
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
  
  
  pushViewport(viewport(layout=grid.layout(
                          nrow,
                          length(x1),
                          widths=x1,
                          heights=unit(rep(1, max(yTE)), "lines"))))
  ##
  ##
  j <- 1
  ##
  ##if (lsel){
  for (i in seq(along=leftcols)){
    drawLabelCol(cols[[leftcols[i]]], j)
    ##
    if (!is.na(yHeadadd)){
      if (!is.null(lab.e.attach.to.col)){
        if (leftcols[i]==paste("col.", lab.e.attach.to.col, sep=""))
          drawLabelCol(col.lab.e, j)
      }
      else if (inherits(x, "metabin")){
        if (leftcols[i]=="col.n.e")
          drawLabelCol(col.lab.e, j)
      }
      else if (inherits(x, "metacont")){
        if (leftcols[i]=="col.sd.e")
          drawLabelCol(col.lab.e, j)
      }
      else if (inherits(x, "metainc")){
        if (leftcols[i]=="col.time.e")
          drawLabelCol(col.lab.e, j)
      }
      ##
      if (!is.null(lab.c.attach.to.col)){
        if (leftcols[i]==paste("col.", lab.c.attach.to.col, sep=""))
          drawLabelCol(col.lab.c, j)
      }
      else if (inherits(x, "metabin")){
        if (leftcols[i]=="col.n.c")
          drawLabelCol(col.lab.c, j)
      }
      else if (inherits(x, "metacont")){
        if (leftcols[i]=="col.sd.c")
          drawLabelCol(col.lab.c, j)
      }
      else if (inherits(x, "metainc")){
        if (leftcols[i]=="col.time.c")
          drawLabelCol(col.lab.c, j)
      }
    }
    ##    if (!is.na(yHeadadd)){
    ##      ##if (inherits(x, "metabin")){
    ##      if (TRUE){
    ##        if (leftcols[i]=="col.n.e")
    ##          drawLabelCol(col.lab.e, j)
    ##        if (leftcols[i]=="col.n.c")
    ##          drawLabelCol(col.lab.c, j)
    ##      }
    ##      if (inherits(x, "metacont")){
    ##        if (leftcols[i]=="col.sd.e")
    ##          drawLabelCol(col.lab.e, j)
    ##        if (leftcols[i]=="col.sd.c")
    ##          drawLabelCol(col.lab.c, j)
    ##      }
    ##    }
    ##
    j <- j+2
  }
  ##}
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
          if (rightcols[i]=="col.n.e")
            drawLabelCol(col.lab.e, j)
        }
        else if (inherits(x, "metacont")){
          if (rightcols[i]=="col.sd.e")
            drawLabelCol(col.lab.e, j)
        }
        else if (inherits(x, "metainc")){
          if (rightcols[i]=="col.time.e")
            drawLabelCol(col.lab.e, j)
        }
        ##
        if (!is.null(lab.c.attach.to.col)){
          if (rightcols[i]==paste("col.", lab.c.attach.to.col, sep=""))
            drawLabelCol(col.lab.c, j)
        }
        else if (inherits(x, "metabin")){
          if (rightcols[i]=="col.n.c")
            drawLabelCol(col.lab.c, j)
        }
        else if (inherits(x, "metacont")){
          if (rightcols[i]=="col.sd.c")
            drawLabelCol(col.lab.c, j)
        }
        else if (inherits(x, "metainc")){
          if (rightcols[i]=="col.time.c")
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
