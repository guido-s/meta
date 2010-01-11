forest <- function(x,
                   byvar=x$byvar,
                   bylab=x$bylab,
                   print.byvar=x$print.byvar,
                   sortvar,
                   studlab=TRUE,
                   ##
                   level=x$level,
                   level.comb=x$level.comb,
                   ##
                   comb.fixed=x$comb.fixed, comb.random=x$comb.random,
                   overall=TRUE,
                   text.fixed="Fixed effect model",
                   text.random="Random effects model",
                   lty.fixed=2, lty.random=3,
                   ##
                   xlab=NULL, xlab.pos=ref, xlim,
                   ##
                   allstudies=TRUE,
                   weight,
                   ##
                   ref=ifelse(x$sm %in% c("RR", "OR", "HR"), 1, 0),
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
                   lwd=1,
                   ##
                   at=NULL,
                   label=TRUE,
                   ##
                   fontsize=12,
                   boxsize=0.8,
                   ##
                   plotwidth=unit(6, "cm"),
                   colgap=unit(2, "mm"),
                   ##
                   col.i="black",
                   col.by="darkgray",
                   ##
                   digits=2){
  
  
  if (!inherits(x, "meta"))
    stop("Argument 'x' must be an object of class \"meta\"")
  
  
  x.name <- deparse(substitute(x))
  
  byvar.name <- deparse(substitute(byvar))
  ##
  if (!missing(byvar) & length(x$byvar)==0){
    if (!is.null(x[[byvar.name]]))
      byvar <- x[[byvar.name]]
  }
  
  
  if (length(comb.fixed)==0){
    comb.fixed <- FALSE
  }
  ##
  if (length(comb.random)==0){
    comb.random <- FALSE
  }
  ##
  if (length(print.byvar)==0){
    print.byvar <- TRUE
  }
  ##
  if (length(lab.e)==0){
    lab.e <- "Experimental"
  }
  ##
  if (length(lab.c)==0){
    lab.c <- "Control"
  }
  
  
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

  if (missing(weight))
    weight <- ifelse(comb.random & !comb.fixed, "random", "fixed")
  
  wcalc <- function(x)
    max(unit(rep(1, length(x)), "grobwidth", x))
  
  
  formatcol <- function(x, y, rows){
    ##
    ##
    ##
    res <- list(labels=
                lapply(c(x, as.list(y)),
                       textGrob, x=1, just="right",
                       gp=gpar(
                         fontsize=fontsize)
                       ),
                rows=rows)
    ##
    res$labels[[1]] <- textGrob(x,
                                x=1, just="right",
                                gp=gpar(
                                  fontsize=fontsize,
                                  fontface="bold")
                                )
    ##
    res$labels[[2]] <- textGrob(y[1],
                                x=1, just="right",
                                gp=gpar(
                                  fontsize=fontsize,
                                  fontface="bold")
                                )
    ##
    res$labels[[3]] <- textGrob(y[2],
                                x=1, just="right",
                                gp=gpar(
                                  fontsize=fontsize,
                                  fontface="bold")
                                )
    ##
    if (by)
      for (i in 1:n.by)
        res$labels[[3+i]] <- textGrob(y[2+i],
                                      x=1, just="right",
                                      gp=
                                      gpar(
                                           fontsize=fontsize,
                                           fontface="bold",
                                           col=col.by)
                                      )
    ##
    res
  }
  
  
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
  
  
  drawNormalCI <- function(low, eff, upp, size, min, max, col){
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
    if (!is.na(eff) && (eff >= min & eff <= max)){
      if (!is.na(size) && size > 0){
        grid.rect(x=unit(eff, "native"),
                  width=unit(size, "snpc"),
                  height=unit(size, "snpc"),
                  gp=gpar(fill=col, col=col))
      }
      else
        grid.lines(x=unit(c(eff, eff), "native"),
                   y = unit(c(0.4, 0.6), "npc"),
                   gp=gpar(col=col, lwd=lwd))
    }
    if (!is.na(eff)){
      ##
      ## Draw line white if totally inside rect
      ##
      ##if ((!is.na(size) && size > 0) |
      ##    (!is.na(size) && (size == 0 & print.infinite))){
      if (!is.na(size)){
        lineCol <- if (size > 0 &&
                       (convertX(unit(eff, "native") + unit(0.5*size, "lines"),
                                 "native", valueOnly=TRUE) > upp) &&
                       (convertX(unit(eff, "native") - unit(0.5*size, "lines"),
                                 "native", valueOnly=TRUE) < low))
          "white"
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
  
  drawSummaryCI <- function(low, eff, upp, size, min, max) {
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
                   gp=gpar(fill=col.by))
  }
  
  
  drawDataCol <- function(col, j) {
    ##
    ## Function to draw a "data" column
    ##
    pushViewport(viewport(layout.pos.col=j, xscale=col$range))
    ##
    ## Reference line:
    ##
    if (!is.na(ref) && (col$range[1] <= ref & ref <= col$range[2]))
      grid.lines(x=unit(ref, "native"), y=0:1,
                 gp=gpar(lwd=lwd))
    ##
    ## 
    ##
    if (comb.fixed & overall)
      if (col$range[1] <= TE.fixed & TE.fixed <= col$range[2] )
        grid.lines(x=unit(TE.fixed, "native"), y=0:1,
                   gp=gpar(lty=lty.fixed, lwd=lwd)) # lty="dashed"
    if (comb.random & overall)
      if (col$range[1] <= TE.random & TE.random <= col$range[2] )
        grid.lines(x=unit(TE.random, "native"), y=0:1,
                   gp=gpar(lty=lty.random, lwd=lwd))
    ##
    if (log){
      if (is.null(at)){
        x1000 <- c(0.001, 0.1, 1,  10, 1000)
        x100  <- c(0.01 , 0.1, 1,  10, 100)
        x10   <- c(0.1  , 0.5, 1,   2, 10)
        x5    <- c(0.2  , 0.5, 1,   2, 5)
        x2    <- c(0.5  , 1, 2)
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
        at <- log(label)
        label <- round(label, 2)
      }
      else{
        if (label)
          label <- round(at, 2)
        at <- log(at)
      }
      grid.xaxis(at=at, label=label,
                 gp=gpar(fontsize=fontsize, lwd=lwd))
    }
    else{
      grid.xaxis(gp=gpar(fontsize=fontsize, lwd=lwd))
    }
    ##
    ## Label on x-axis:
    ##
    grid.text(xlab,
              x=unit(xlab.pos, "native"),
              y=unit(-2.35, "lines"),
              just=c("center", "top"),
              gp=gpar(fontsize=fontsize))
    ##
    popViewport()
    for (i in 1:length(col$rows)) {
      if (!is.na(col$rows[i])){
        pushViewport(viewport(layout.pos.row=col$rows[i], layout.pos.col=j,
                              xscale=col$range))
        if (col$type[i] == "n")
          drawNormalCI(col$low[i], col$eff[i], col$upp[i], col$sizes[i],
                       col$range[1], col$range[2], col$col[i])
        else
          drawSummaryCI(col$low[i], col$eff[i], col$upp[i], col$sizes[i],
                        col$range[1], col$range[2])
        popViewport()
      }
    }
  }
  
  
  if (is.null(xlab)){
    if      (x$sm=="OR" ) xlab <- "Odds Ratio"
    else if (x$sm=="RD" ) xlab <- "Risk Difference"
    else if (x$sm=="RR" ) xlab <- "Relative Risk"
    else if (x$sm=="SMD") xlab <- "Standardised mean difference"
    else if (x$sm=="WMD"|x$sm=="MD") xlab <- "Mean difference"
    else if (x$sm=="HR" ) xlab <- "Hazard Ratio"
    else if (x$sm=="AS" ) xlab <- "Arcus Sinus Transformation"
    else if (x$sm=="proportion" ) xlab <- "Proportion"
    else xlab <- x$sm
  }
  
  
  iweight <- charmatch(tolower(weight),
                       c("same", "fixed", "random"), nomatch = NA)
  ##
  if(is.na(iweight))
    stop("weight should be \"same\", \"fixed\", or \"random\"")
  ##
  weight <- c("same", "fixed", "random")[iweight]
  
  
  if (inherits(x, "metaprop")){
    x$event.e <- x$event
    x$n.e <- x$n
  }
  
  
  if (inherits(x, "metainf")|inherits(x, "metacum")){
    ##
    x$TE.fixed    <- rev(x$TE)[1]
    x$seTE.fixed  <- rev(x$seTE)[1]
    x$TE.random   <- rev(x$TE)[1]
    x$seTE.random <- rev(x$seTE)[1]
    x$w.all <- rev(x$w)[1]
    ##
    x$TE <- rev(rev(x$TE)[-(1:2)])
    x$seTE <- rev(rev(x$seTE)[-(1:2)])
    x$studlab <- rev(rev(x$studlab)[-(1:2)])
    ##
    x$w.fixed <- rev(rev(x$w)[-(1:2)])
    x$w.random <- rev(rev(x$w)[-(1:2)])
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
  ## Check for duplicate columns
  ##
  if (length(c(rightcols, leftcols))>0 &&
      any(duplicated(c(rightcols, leftcols))))
    stop("Duplicate entries in 'leftcols' and 'rightcols'")
  
  
  ##
  ## Predefined columns and labels
  ##
  colnames <- c("studlab", "TE", "seTE",
                "n.e", "n.c",
                "event.e", "event.c",
                "mean.e", "mean.c",
                "sd.e", "sd.c",
                "effect", "ci",
                "w.fixed", "w.random")
  ##
  labnames <- c("Study", "TE", "seTE",
                "Total", "Total", "Events", "Events",
                "Mean", "Mean", "SD", "SD",
                x$sm, paste(100*level, "%-CI", sep=""),
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
    ##
    ## Check whether additional variables are
    ## part of meta-object
    ##
    for (i in colnames.new){
      if (length(x[[i]]) == 0)
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
    if (inherits(x, "metacont"))
      leftcols <- c("studlab",
                    "n.e", "mean.e", "sd.e",
                    "n.c", "mean.c", "sd.c")
    if (inherits(x, "metagen"))
      leftcols <- c("studlab",
                    "TE", "seTE")
    if (inherits(x, "metainf") | inherits(x, "metacum"))
      leftcols <- "studlab"
    if (inherits(x, "metaprop"))
      leftcols <- c("studlab",
                    "event.e", "n.e")
  }
  ##
  if (is.null(rightcols) & rsel){
    rightcols <- c("effect", "ci")
    ##
    if (!(inherits(x, "metainf") |
          inherits(x, "metacum"))){
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
      labs[[paste("lab.", leftcols[i], sep="")]] <- leftlabs[i]
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
      labs[[paste("lab.", rightcols[i], sep="")]] <- rightlabs[i]
    }
  }
  
  
  ##
  ## total number of trials to plot (*not* number of trials combined)
  ##
  k.all <- length(x$TE)


  if (allstudies)
    n.stud <- k.all # all trials
  else
    n.stud <- x$k   # number of trials combined in meta-analysis

  
  by <- length(byvar)>0
  sort <- !missing(sortvar)
  
  ##  if (by)
  ##    if (!is.null(x[[byvar.name]]))
  ##      byvar <- x[[byvar.name]]
  
  if (by & (inherits(x, "metainf")| inherits(x, "metacum")))
    stop("Use of 'byvar' not possible for 'metainf' or 'metacum' object.") 
  
  if (!by) byvar <- rep(1, k.all)
  if (!sort) sortvar <- rep(1, k.all)

  if (sort & length(sortvar) != k.all)
    stop("'x' and 'sortvar' have different length")
  if (by & length(byvar) != k.all)
    stop("'x' and 'byvar' have different length")
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
  
  
  if (length(col.i)==1)
    col.i <- rep(col.i, length(x$TE))
  else if (length(col.i)!=length(x$TE))
    stop("Parameter 'col.i' has different length than number of studies in meta-analysis")
  
  
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
  x$n   <- x$n[sel]
  x$n.e <- x$n.e[sel]
  x$n.c <- x$n.c[sel]
  ##
  x$event   <- x$event[sel]
  x$event.e <- x$event.e[sel]
  x$event.c <- x$event.c[sel]
  ##
  x$mean.e <- x$mean.e[sel]
  x$mean.c <- x$mean.c[sel]
  ##
  x$sd.e <- x$sd.e[sel]
  x$sd.c <- x$sd.c[sel]
  ##
  x$TE   <- x$TE[sel]
  x$seTE <- x$seTE[sel]
  ##
  x$w.fixed  <- x$w.fixed[sel]
  x$w.random <- x$w.random[sel]
  studlab  <- studlab[sel]
  ##
  byvar   <- byvar[sel]
  sortvar <- sortvar[sel]
  ##
  col.i <- col.i[sel]

  
  if (sort | by){
    ##
    o <- order(byvar, sortvar)
    ##
    x$n   <- x$n[o]
    x$n.e <- x$n.e[o]
    x$n.c <- x$n.c[o]
    ##
    x$event   <- x$event[o]
    x$event.e <- x$event.e[o]
    x$event.c <- x$event.c[o]
    ##
    x$mean.e <- x$mean.e[o]
    x$mean.c <- x$mean.c[o]
    ##
    x$sd.e <- x$sd.e[o]
    x$sd.c <- x$sd.c[o]
    ##
    x$TE   <- x$TE[o]
    x$seTE <- x$seTE[o]
    ##
    x$w.fixed  <- x$w.fixed[o]
    x$w.random <- x$w.random[o]
    studlab  <- studlab[o]
    ##
    byvar   <- byvar[o]
    sortvar <- sortvar[o]
    ##
    col.i <- col.i[o]
  }
  
  
  by.levs <- unique(byvar)
  
  
  sum.meta <- summary(x, level=level, level.comb=level.comb, warn=FALSE)
  ##
  TE    <- sum.meta$study$TE
  seTE  <- sum.meta$study$seTE
  lowTE <- sum.meta$study$lower
  uppTE <- sum.meta$study$upper
  ##
  TE.fixed    <- sum.meta$fixed$TE
  lowTE.fixed <- sum.meta$fixed$lower
  uppTE.fixed <- sum.meta$fixed$upper
  ##
  TE.random    <- sum.meta$random$TE
  lowTE.random <- sum.meta$random$lower
  uppTE.random <- sum.meta$random$upper
  
  
  if (by){

    res.w <- matrix(NA, ncol=14, nrow=length(by.levs))
    j <- 0
    ##
    for ( i in by.levs){
      j <- j+1
      sel <- byvar == i
      ##
      if (inherits(x, "metabin")){
        meta1 <- metabin(x$event.e[sel], x$n.e[sel],
                         x$event.c[sel], x$n.c[sel],
                         method=x$method, sm=x$sm,
                         incr=x$incr, allincr=x$allincr,
                         addincr=x$addincr, allstudies=x$allstudies,
                         MH.exact=x$MH.exact, RR.cochrane=x$RR.cochrane,
                         warn=x$warn)
      }
      ##
      if (inherits(x, "metacont")){
        meta1 <- metacont(x$n.e[sel], x$mean.e[sel], x$sd.e[sel],
                          x$n.c[sel], x$mean.c[sel], x$sd.c[sel],
                          sm=x$sm)
      }
      ##
      if (inherits(x, "metagen")){
        meta1 <- metagen(x$TE[sel], x$seTE[sel], sm=x$sm)
      }
      ##
      if (inherits(x, "metaprop")){
        meta1 <- metaprop(event=x$event[sel], n=x$n[sel],
                          freeman.tukey=x$freeman.tukey)
        meta1$event.e <- meta1$event
        meta1$n.e <- meta1$n
      }
      ##
      res.w[j,] <- c(meta1$TE.fixed, meta1$seTE.fixed,
                     meta1$TE.random, meta1$seTE.random,
                     meta1$Q, meta1$k, length(meta1$TE),
                     summary(meta1)$I2$TE,
                     sum(meta1$w.fixed), sum(meta1$w.random),
                     sum(meta1$event.e), sum(meta1$n.e),
                     sum(meta1$event.c), sum(meta1$n.c))
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
    Q.w           <- res.w[sel,5]
    I2.w          <- res.w[sel,8]
    w.fixed.w     <- res.w[sel,9]
    w.random.w    <- res.w[sel,10]
    e.e.w         <- res.w[sel,11]
    n.e.w         <- res.w[sel,12]
    e.c.w         <- res.w[sel,13]
    n.c.w         <- res.w[sel,14]
    ##
    by.levs <- by.levs[sel]
    k.all.w <- k.all.w[sel]
    ##
    if (comb.fixed | !(comb.fixed|comb.random)){
      if (comb.random)
        warning("Estimate from fixed effect model used in groups defined by 'byvar'")
      ##
      comb.w <- "fixed"
      meta.w <- ci(TE.fixed.w, seTE.fixed.w, level)
      weight.w <- w.fixed.w
      text.w <- text.fixed
    }
    if (!comb.fixed & comb.random){
      comb.w <- "random"
      meta.w <- ci(TE.random.w, seTE.random.w, level)
      weight.w <- w.random.w
      text.w <- text.random
    }
    ##
    TE.w <- meta.w$TE
    lowTE.w <- meta.w$lower
    uppTE.w <- meta.w$upper
    ##
    if (sum(weight.w)>0)
      weight.w.p <- 100*round(weight.w/sum(weight.w, na.rm=TRUE), 3)
    else
      weight.w.p <- weight.w
  }
  
  
  if (inherits(x, "metaprop")){
    ref <- NA
    ##
    denum <- 1 + x$freeman.tukey
    ##
    TE.fixed <- sin(TE.fixed/denum)^2
    lowTE.fixed <- sin(lowTE.fixed/denum)^2
    uppTE.fixed <- sin(uppTE.fixed/denum)^2
    ##
    TE.random <- sin(TE.random/denum)^2
    lowTE.random <- sin(lowTE.random/denum)^2
    uppTE.random <- sin(uppTE.random/denum)^2
    ##
    if (by){
      TE.w <- sin(TE.w/denum)^2
      lowTE.w <- sin(lowTE.w/denum)^2
      uppTE.w <- sin(uppTE.w/denum)^2
    }
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

  if (inherits(x, "metainf")|inherits(x, "metacum")){
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
  
  
  if (by)
    n.by <- length(by.levs)
  else
    n.by <- 0
  
  
  if (print.byvar){
    if (missing(bylab))
      bylab <- paste(byvar.name,
                     " = ",
                     format(by.levs), sep="")
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
    modlabs <- c(text.fixed, text.random, bylab,
              rep(text.w, n.by), studlab)
    ##
    TEs    <- c(TE.fixed, TE.random, TE.w, TE)
    lowTEs <- c(lowTE.fixed, lowTE.random, lowTE.w, lowTE)
    uppTEs <- c(uppTE.fixed, uppTE.random, uppTE.w, uppTE)
    ##
    TEs.study <- c(NA, NA, rep(NA, n.by), TE)
    seTEs.study <- c(NA, NA, rep(NA, n.by), seTE)
    ##
    w.fixeds  <- c(100, "--", weight.w.p, w.fixed.p)
    w.randoms <- c("--", 100, weight.w.p, w.random.p)
    ##
    sel.fixed  <- w.fixeds=="--"
    sel.random <- w.randoms=="--"
    ##
    if (comb.w=="fixed")
      sel.random[2+1:n.by] <- TRUE
    if (comb.w=="random")
      sel.fixed[2+1:n.by] <- TRUE
  }
  else{
    modlabs <- c(text.fixed, text.random, studlab)
    ##
    TEs    <- c(TE.fixed, TE.random, TE)
    lowTEs <- c(lowTE.fixed, lowTE.random, lowTE)
    uppTEs <- c(uppTE.fixed, uppTE.random, uppTE)
    ##
    TEs.study <- c(NA, NA, TE)
    seTEs.study <- c(NA, NA, seTE)
    ##
    w.fixeds  <- c(100, "--", w.fixed.p)
    w.randoms <- c("--", 100, w.random.p)
    ##
    sel.fixed <- w.fixeds=="--"
    sel.random <- w.randoms=="--"
  }
  ##
  ## Treatment effect and confidence interval
  ##
  if (x$sm %in% c("HR", "OR", "RR")){
    effect.format <- format(round(exp(TEs), digits))
    ci.format <- p.ci(format(round(exp(lowTEs), digits)),
                      format(round(exp(uppTEs), digits)))
    ref <- log(ref)
    log <- TRUE
  }
  else{
    effect.format <- format(round(TEs, digits))
    ci.format <- p.ci(format(round(lowTEs, digits)),
                      format(round(uppTEs, digits)))
    log <- FALSE
  }
  ##
  effect.format <- gsub("NA", "  ", effect.format)
  ##
  ## Weights of fixed and random effects model
  ##
  w.fixed.format  <- paste(formatC(w.fixeds, 4), "%", sep="")
  w.random.format <- paste(formatC(w.randoms, 4), "%", sep="")
  ##
  w.fixed.format[sel.fixed] <- "--"
  w.random.format[sel.random] <- "--"
  ##
  ## Treatment estimate and its standard error
  ##
  TE.format <- ifelse(is.na(TEs.study), "",
                      format(round(TEs.study, digits)))
  seTE.format <- ifelse(is.na(seTEs.study), "",
                        format(round(seTEs.study, 4)))
  ##
  ## Number of Events and patients
  ##
  sum.n.e <- sum(x$n.e, na.rm=TRUE)
  sum.n.c <- sum(x$n.c, na.rm=TRUE)
  sum.e.e <- sum(x$event.e, na.rm=TRUE)
  sum.e.c <- sum(x$event.c, na.rm=TRUE)
  ##
  if (by){
    Ne <- c(sum.n.e, sum.n.e, n.e.w, x$n.e)
    Nc <- c(sum.n.c, sum.n.c, n.c.w, x$n.c)
    Ee <- c(sum.e.e, sum.e.e, e.e.w, x$event.e)
    Ec <- c(sum.e.c, sum.e.c, e.c.w, x$event.c)
  }
  else{
    Ne <- c(sum.n.e, sum.n.e, x$n.e)
    Nc <- c(sum.n.c, sum.n.c, x$n.c)
    Ee <- c(sum.e.e, sum.e.e, x$event.e)
    Ec <- c(sum.e.c, sum.e.c, x$event.c)
  }
  ##
  Ne.format <- format(Ne)
  Nc.format <- format(Nc)
  Ee.format <- format(Ee)
  Ec.format <- format(Ec)
  ##
  if (comb.fixed & comb.random){
    Ne.format[2] <- ""
    Nc.format[2] <- ""
    Ee.format[2] <- ""
    Ec.format[2] <- ""
  }
  ##
  ## Mean and standard deviation
  ##
  if (by){
    Me <- c("", "", rep("", length(n.e.w)), x$mean.e)
    Mc <- c("", "", rep("", length(n.c.w)), x$mean.c)
    Se <- c("", "", rep("", length(n.e.w)), x$sd.e)
    Sc <- c("", "", rep("", length(n.c.w)), x$sd.c)
  }
  else{
    Me <- c("", "", x$mean.e)
    Mc <- c("", "", x$mean.c)
    Se <- c("", "", x$sd.e)
    Sc <- c("", "", x$sd.c)
  }
  ##
  Me.format <- formatC(Me)
  Mc.format <- formatC(Mc)
  Se.format <- formatC(Se)
  Sc.format <- formatC(Sc)
  
  
  if (!comb.fixed) text.fixed <- ""
  if (!comb.random) text.random <- ""
  
  
  ##
  ## y-axis:
  ##
  if (any(rightcols %in% c("n.e", "n.c")) |
      any(leftcols  %in% c("n.e", "n.c"))){
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
  }
  else{
    ##
    if (comb.fixed|comb.random){
      rows.up   <- 1
      rows.down <- 2
      ##
      N <- n.stud + (rows.up+rows.down)*n.by
      ##
      sel <- rep(rep(c(FALSE, TRUE, FALSE), n.by),
                 as.vector(rbind(rep(rows.up, n.by),
                                 k.all.w,
                                 rep(rows.down, n.by))))
      ##
      yTE <- (1:N)[sel]
      ##
      yBylab <- (1:N)[!sel][rep(c(TRUE, FALSE, FALSE), n.by)]
      ##
      yTE.w <- (1:N)[!sel][rep(c(FALSE, TRUE, FALSE), n.by)]
    }
    else{
      rows.up   <- 1
      rows.down <- 1
      ##
      N <- n.stud + (rows.up+rows.down)*n.by
      ##
      sel <- rep(rep(c(FALSE, TRUE, FALSE), n.by),
                 as.vector(rbind(rep(rows.up, n.by),
                                 k.all.w,
                                 rep(rows.down, n.by))))
      ##
      yTE <- (1:N)[sel]
      ##
      yBylab <- (1:N)[!sel][rep(c(TRUE, FALSE), n.by)]
      ##
      yTE.w <- rep(NA, n.by)
    }
  }
  
  
  ##
  ## x-axis:
  ##
  if (!missing(xlim))
    if (x$sm %in% c("HR", "OR", "RR"))
      xlim <- log(xlim)
  ##
  if (missing(xlim)){
    if (inherits(x, "metaprop"))
      xlim <- c(0,1)
    else{
      sel.low <- is.finite(lowTE)
      sel.upp <- is.finite(uppTE)
      xlim <- c(min(lowTE[sel.low], na.rm=TRUE),
              max(uppTE[sel.upp], na.rm=TRUE))
    }
  }
  
  if (by)
    max.yTE <- max(c(yTE, yTE.w), na.rm=TRUE)
  else
    max.yTE <- max(yTE, na.rm=TRUE)
  
  
  ##
  ##
  ##
  yTE.fixed  <- NA
  yTE.random <- NA
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
  ##
  ##
  yTE        <- yHead + yTE
  yTE.fixed  <- yHead + yTE.fixed
  yTE.random <- yHead + yTE.random
  ##
  if (by){
    yBylab <- yHead + yBylab
    yTE.w  <- yHead + yTE.w
  }
  
  
  if (by){
    yLab <- c(yHead,
              yTE.fixed, yTE.random,
              yBylab, yTE.w,
              yTE)
    ##
    yS <- c(yHead, yTE.fixed, yTE.random, yTE.w, yTE)
  }
  else{
    yLab <- c(yHead, yTE.fixed, yTE.random, yTE)
    yS   <- c(yHead, yTE.fixed, yTE.random, yTE)
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
                             textGrob, x=0, just="left",
                             gp=gpar(
                               fontsize=fontsize)),
                      rows=yLab
                      )
  ##
  ##
  col.studlab$labels[[1]] <- textGrob(labs[["lab.studlab"]],
  ##col.studlab$labels[[1]] <- textGrob(lab.studlab,
                                      x=0, just="left",
                                      gp=gpar(
                                        fontsize=fontsize,
                                        fontface="bold")
                                      )
  ##
  col.studlab$labels[[2]] <- textGrob(text.fixed,
                                      x=0, just="left",
                                      gp=gpar(
                                        fontsize=fontsize,
                                        fontface="bold")
                                      )
  ##
  col.studlab$labels[[3]] <- textGrob(text.random,
                                      x=0, just="left",
                                      gp=gpar(
                                        fontsize=fontsize,
                                        fontface="bold")
                                      )
  ##
  if (by){
    for (i in 1:n.by){
      col.studlab$labels[[3+i]] <- textGrob(bylab[i],
                                            x=0, just="left",
                                            gp=
                                            gpar(
                                                 fontsize=fontsize,
                                                 fontface="bold",
                                                 col=col.by)
                                            )
      ##
      col.studlab$labels[[3+n.by+i]] <- textGrob(text.w,
                                                 x=0, just="left",
                                                 gp=
                                                 gpar(
                                                      fontsize=fontsize,
                                                      fontface="bold",
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
                     ## "s" means summary, "n" means normal
                     ##
                     type=c(rep("s", length(TEs)-length(TE)),
                       rep("n", length(TE))),
                     col=c(rep("", length(TEs)-length(TE)),
                       col.i)
                     )
  ##
  ## Sizes of boxes
  ##
  if (weight=="same"){
    information <- rep(0.9, length(TEs))
  }
  else{
    if (weight=="fixed")
      information <- sqrt(as.numeric(w.fixeds[-(1:2)]))
    if (weight=="random")
      information <- sqrt(as.numeric(w.randoms[-(1:2)]))
    ##
    information <- information/max(information, na.rm=TRUE)
    information <- c(max(information), max(information),
                     information)
    ##information <- sqrt(1 / ((col.forest$upp - col.forest$eff)/1.96))
  }
  ##
  col.forest$sizes <- information
  col.forest$sizes <- col.forest$sizes * boxsize
  ##
  ## Width of column 3
  col.forestwidth <- plotwidth
  ##
  ## Range on the x-axis for column 3
  col.forest$range <- xlim # c(0, 2)
  
  
  if (by)
    nrow <- max(c(yTE, yTE.fixed, yTE.random, yTE.w), na.rm=TRUE)
  else
    nrow <- max(c(yTE, yTE.fixed, yTE.random), na.rm=TRUE)
  
  
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
  ##                                 )
  if (newcols){
    if (by){
      for (i in seq(along=rightcols.new)){
        tname <- paste("col.", rightcols.new[i], sep="")
        ##cols[[tname]] <- formatcol(ifelse(missing(rightlabs),
        ##                                  rightcols.new[[i]],
        ##                                  rightlabs[i]),
        cols[[tname]] <- formatcol(rightlabs.new[i],
                                   c("", "",
                                     rep("", length(TE.w)),
                                     x[[rightcols.new[i]]][o]),
                                   yS)
      }
      for (i in seq(along=leftcols.new)){
        tname <- paste("col.", leftcols.new[i], sep="")
        ##cols[[tname]] <- formatcol(ifelse(missing(leftlabs),
        ##                                  leftcols.new[[i]],
        ##                                  leftlabs[i]),
        cols[[tname]] <- formatcol(leftlabs.new[i],
                                   c("", "",
                                     rep("", length(TE.w)),
                                     x[[leftcols.new[i]]][o]),
                                   yS)
      }
    }
    else{
      for (i in seq(along=rightcols.new)){
        tname <- paste("col.", rightcols.new[i], sep="")
        cols[[tname]] <- formatcol(rightlabs.new[i],
                                   c("", "",
                                     if (sort) x[[rightcols.new[i]]][o]
                                     else x[[rightcols.new[i]]]
                                     ),
                                   yS)
      }
      for (i in seq(along=leftcols.new)){
        tname <- paste("col.", leftcols.new[i], sep="")
        cols[[tname]] <- formatcol(leftlabs.new[i],
                                   c("", "",
                                     if (sort) x[[leftcols.new[i]]][o]
                                     else x[[leftcols.new[i]]]
                                     ),
                                   yS)
      }
    }
  }
  
  
  col.lab.e <- list(labels=list(textGrob(lab.e,
                      x=1, just="right",
                      gp=gpar(
                        fontsize=fontsize,
                        fontface="bold")
                      )),
                    rows=1)
  ##
  col.lab.c <- list(labels=list(textGrob(lab.c,
                      x=1, just="right",
                      gp=gpar(
                        fontsize=fontsize,
                        fontface="bold")
                      )),
                    rows=1)
  
  
  leftcols  <- paste("col.", leftcols, sep="")
  rightcols <- paste("col.", rightcols, sep="")
  
  
  for (i in seq(along=leftcols)){
    if (i==1)
      x1 <- unit.c(wcalc(cols[[leftcols[i]]]$labels))
    else
      x1 <- unit.c(x1,
                   colgap,
                     wcalc(cols[[leftcols[i]]]$labels))
  }
  ##
  x1 <- unit.c(x1, colgap, col.forestwidth)
  ##
  if (rsel)
    for (i in seq(along=rightcols)){
      x1 <- unit.c(x1,
                   colgap,
                   wcalc(cols[[rightcols[i]]]$labels))
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
      }
      ##
      j <- j+2
    }
  }
  ##
  popViewport()
  
  
  invisible(NULL)
}
