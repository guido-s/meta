metacr <- function(x, comp.no=1, outcome.no=1,
                   sm, method, comb.fixed, comb.random,
                   swap.events, logscale,
                   method.tau="DL",
                   hakn=FALSE,
                   title, complab, outclab,
                   warn=FALSE){
  ##
  if (!inherits(x, "rm5"))
    stop("Argument 'x' must be an object of class \"rm5\"")
  ##
  sel <- x$comp.no==comp.no & x$outcome.no==outcome.no
  ##
  if (sum(sel)==0){
    warning("No data available for comp.no=", comp.no,
            " and outcome.no=", outcome.no)
    return(NULL)
  }
  ##
  x$sel <- sel
  ##
  if (missing(title))
    title   <- attributes(x)$title
  ##
  if (missing(complab))
    complab <- unique(x$complab[sel])
  ##
  if (missing(outclab))
    outclab <- unique(x$outclab[sel])
  ##
  label.e <- unique(x$label.e[sel])
  label.c <- unique(x$label.c[sel])
  ##
  label.left  <- unique(x$label.left[sel])
  label.right <- unique(x$label.right[sel])
  ##
  type <- unique(x$type[sel])
  ##
  if (missing(sm))
    sm <- unique(x$sm[sel])
  ##
  if (missing(method))
    method <- unique(x$method[sel])
  ##
  if (missing(comb.fixed))
    comb.fixed  <- unique(x$comb.fixed[sel])
  ##
  if (missing(comb.random))
    comb.random <- unique(x$comb.random[sel])
  
  
  if (!all(is.na(x$logscale[sel]))){
    if (!unique(x$logscale[sel]))
      x$TE[sel] <- log(x$TE[sel])
  }
  else{
    if (!missing(logscale)){
      if (!logscale)
        x$TE[sel] <- log(x$TE[sel])
    }
    else
      if ((type=="I" & method!="Peto") & sm %in% c("HR", "RR", "OR"))
        warning("Assuming that values for 'TE' are on log scale. Please use argument 'logscale=FALSE' if values are on natural scale.")
  }
  
  ##
  if (type=="D"){
    ##
    if (missing(swap.events))
      swap.events <- unique(x$swap.events[sel])
    if (!is.na(swap.events) && swap.events){
      event.e  <- x$n.e[sel] - x$event.e[sel]
      event.c  <- x$n.c[sel] - x$event.c[sel]
    }
    else{
      event.e <- x$event.e[sel]
      event.c <- x$event.c[sel]
    }
    m1 <- metabin(event.e, x$n.e[sel], event.c, x$n.c[sel],
                  sm=sm, method=method, studlab=x$studlab[sel],
                  comb.fixed=comb.fixed, comb.random=comb.random,
                  method.tau=method.tau, hakn=hakn,
                  title=title,
                  complab=complab, outclab=outclab,
                  label.e=label.e, label.c=label.c,
                  label.left=label.left, label.right=label.right,
                  RR.cochrane=TRUE, warn=warn)
  }
  ##
  if (type=="C")
    m1 <- metacont(x$n.e[sel], x$mean.e[sel], x$sd.e[sel],
                   x$n.c[sel], x$mean.c[sel], x$sd.c[sel],
                   sm=sm, studlab=x$studlab[sel],
                   comb.fixed=comb.fixed, comb.random=comb.random,
                   method.tau=method.tau, hakn=hakn,
                   title=title,
                   complab=complab, outclab=outclab,
                   label.e=label.e, label.c=label.c,
                   label.left=label.left, label.right=label.right)
  ##
  if (type=="P")
    m1 <- metagen(x$O.E[sel]/x$V[sel], sqrt(1/x$V[sel]),
                  sm=sm, studlab=x$studlab[sel],
                  comb.fixed=comb.fixed, comb.random=comb.random,
                  method.tau=method.tau, hakn=hakn,
                  title=title,
                  complab=complab, outclab=outclab,
                  label.e=label.e, label.c=label.c,
                  label.left=label.left, label.right=label.right)
  ##
  if (type=="I" & method!="Peto")
    m1 <- metagen(x$TE[sel], x$seTE[sel],
                  sm=sm, studlab=x$studlab[sel],
                  comb.fixed=comb.fixed, comb.random=comb.random,
                  method.tau=method.tau, hakn=hakn,
                  n.e=x$n.e[sel],
                  n.c=x$n.c[sel],
                  title=title,
                  complab=complab, outclab=outclab,
                  label.e=label.e, label.c=label.c,
                  label.left=label.left, label.right=label.right)
  ##
  if (type=="I" & method=="Peto")
    m1 <- metagen(x$O.E[sel]/x$V[sel], sqrt(1/x$V[sel]),
                  sm=sm, studlab=x$studlab[sel],
                  comb.fixed=comb.fixed, comb.random=comb.random,
                  method.tau=method.tau, hakn=hakn,
                  title=title,
                  complab=complab, outclab=outclab,
                  label.e=label.e, label.c=label.c,
                  label.left=label.left, label.right=label.right)
  ##
  if (length(unique(x$group.no[sel]))>1){
    m1$byvar <- x$grplab[sel]
    m1$bylab <- "grp"
    m1$print.byvar <- FALSE
  }
  ##
  if (sm=="OTHER"){
    warning('Meta-analysis not possible for sm="OTHER".')
    res <- NULL
  }
  else
    res <- m1
  res
}
