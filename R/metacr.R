metacr <- function(x, comp.no=1, outcome.no=1,
                   smother="", logscale=TRUE){
  ##
  if (!inherits(x, "rm5"))
    stop("Argument 'x' must be an object of class \"rm5\"")
  ##
  sel <- x$comp.no==comp.no & x$outcome.no==outcome.no
  ##
  if (sum(sel)==0)
    stop("No data available for comp.no=", comp.no,
         " and outcome.no=", outcome.no)
  ##
  x$sel <- sel
  ##
  title   <- attributes(x)$title
  complab <- unique(x$complab[sel])
  outclab <- unique(x$outclab[sel])
  ##
  label.e <- unique(x$label.e[sel])
  label.c <- unique(x$label.c[sel])
  ##
  type <- unique(x$type[sel])
  meth <- unique(x$meth[sel])
  totals <- unique(x$totals[sel])
  ##
  sm <- unique(x$sm[sel])
  RR.cochrane <- unique(x$RR.cochrane[sel])
  ##
  comb.fixed  <- unique(x$comb.fixed[sel])
  comb.random <- unique(x$comb.random[sel])
  ##
  if (sm=="OTHER" & smother!=""){
    sm <- smother
    if (logscale)
      x$TE[sel] <- log(x$TE[sel])
  }
  ##
  if (type=="D")
    m1 <- metabin(x$event.e[sel], x$n.e[sel],
                  x$event.c[sel], x$n.c[sel],
                  sm=sm, meth=meth, studlab=x$studlab[sel],
                  comb.fixed=comb.fixed, comb.random=comb.random,
                  title=title,
                  complab=complab, outclab=outclab,
                  label.e=label.e, label.c=label.c,
                  RR.cochrane=RR.cochrane)
  ##
  if (type=="C")
    m1 <- metacont(x$n.e[sel], x$mean.e[sel], x$sd.e[sel],
                   x$n.c[sel], x$mean.c[sel], x$sd.c[sel],
                   sm=sm, studlab=x$studlab[sel],
                   comb.fixed=comb.fixed, comb.random=comb.random,
                   title=title,
                   complab=complab, outclab=outclab,
                   label.e=label.e, label.c=label.c)
  ##
  if (type=="P")
    m1 <- metagen(x$O.E[sel]/x$V[sel], sqrt(1/x$V[sel]),
                  sm=sm, studlab=x$studlab[sel],
                  comb.fixed=comb.fixed, comb.random=comb.random,
                  title=title,
                  complab=complab, outclab=outclab,
                  label.e=label.e, label.c=label.c)
  ##
  if (type=="I" & meth!="Peto")
    m1 <- metagen(x$TE[sel], x$seTE[sel],
                  sm=sm, studlab=x$studlab[sel],
                  comb.fixed=comb.fixed, comb.random=comb.random,
                  title=title,
                  complab=complab, outclab=outclab,
                  label.e=label.e, label.c=label.c)
  ##
  if (type=="I" & meth=="Peto")
    m1 <- metagen(x$O.E[sel]/x$V[sel], sqrt(1/x$V[sel]),
                  sm=sm, studlab=x$studlab[sel],
                  comb.fixed=comb.fixed, comb.random=comb.random,
                  title=title,
                  complab=complab, outclab=outclab,
                  label.e=label.e, label.c=label.c)
  ##
  if (length(unique(x$group.no[sel]))>1){
    m1$byvar <- x$grplab[sel]
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
