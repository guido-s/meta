read.rm5 <- function(file, sep=",", quote = "\"",
                     title, numbers.in.labels=TRUE){
  ##
  selvar <- function(x, sel, value=NA){
    res <-
      if (!is.null(x))
        x[sel]
      else value[sel]
    res
  }
  ##
  numchar <- function(x){
    res <- as.numeric(as.character(x))
    res
  }
  ##
  nachar <- function(x){
    res <- x
    res[is.na(res)] <- ""
    res
  }
  ##
  if (missing(title)){
    title <- strsplit(file, "\\.csv$")[[1]]
    tmp <- strsplit(title, "\\/")
    title <- tmp[[1]][length(tmp[[1]])]
  }
  ##
  tdata <- read.table(file, header=TRUE,
                      sep=sep, quote=quote,
                      comment.char="")
  ##
  nam <- names(tdata)
  ##
  ##print(nam)
  ##
  if (!all(c("Comparison.Number", "Outcome.Number", "Subgroup.Number") %in% nam))
    stop("Mandatory fields 'Comparison Number', 'Outcome Number',",
         " and 'Subgroup Number' not included in export file ",
         deparse(substitute(file)),
         " (see help page of function read.rm5).")
  ##
  nam[nam=="Comparison.Number"] <- "comp.no"
  nam[nam=="Outcome.Number"] <- "outcome.no"
  nam[nam=="Subgroup.Number"] <- "group.no"
  nam[nam=="Name"] <- "author"
  ##
  nam[nam=="Data.Type"] <- "type"
  nam[nam=="Statistical.Method"] <- "method"
  nam[nam=="Effect.Measure"] <- "sm"
  nam[nam=="Analysis.Model"] <- "model"
  ##
  nam[nam=="Events.1"] <- "event.e"
  nam[nam=="Mean.1"] <- "mean.e"
  nam[nam=="SD.1"] <- "sd.e"
  nam[nam=="Total.1"] <- "n.e"
  ##
  nam[nam=="Events.2"] <- "event.c"
  nam[nam=="Mean.2"] <- "mean.c"
  nam[nam=="SD.2"] <- "sd.c"
  nam[nam=="Total.2"] <- "n.c"
  ##
  nam[nam=="O.E"] <- "O.E"
  nam[nam=="Var"] <- "V"
  nam[nam=="Effect.Estimate"] <- "TE"
  nam[nam=="SE"] <- "seTE"
  ##
  nam[nam=="CI.Start"] <- "lower"
  nam[nam=="CI.End"]   <- "upper"
  ##
  nam[nam=="Weight"]   <- "weight"
  ##
  nam[nam=="P.Q."]  <- "pval.Q"
  nam[nam=="I..Q."] <- "I2"
  nam[nam=="Tau."]  <- "tau2"
  ##
  nam[nam=="P.Z."]     <- "pval.TE"
  nam[nam=="P.Qint."]  <- "pval.Qint"
  nam[nam=="I..Qint."] <- "I2.Qint"
  ##
  nam[nam=="Group.Label.1"]     <- "label.e"
  nam[nam=="Group.Label.2"]     <- "label.c"
  nam[nam=="Left.Graph.Label"]  <- "label.left"
  nam[nam=="Right.Graph.Label"] <- "label.right"
  ##
  nam[nam=="Year.of.study"]      <- "year"
  nam[nam=="User.defined.order"] <- "order"
  ##
  ## Fancy coding to set variable names of I^2 and tau^2
  ## as "I2" and "tau2" (not necessary for Linux)
  ## I don't know how to do this differently with Windows ...
  ##
  if (!all(c("I2", "tau2") %in% nam)){
    pos1 <- seq(along=nam)[nam=="pval.Q"]
    pos2 <- seq(along=nam)[nam=="Z"]
    if ((pos2-pos1)==3)
      nam[pos1+1:2] <- c("I2", "tau2")
  }
  if (!"I2.Qint" %in% nam){
    pos3 <- seq(along=nam)[nam=="pval.Qint"]
    pos4 <- seq(along=nam)[nam=="df"]
    if ((pos4-pos3)==2)
      nam[pos1+1] <- c("I2.Qint")
  }
  ##
  names(tdata) <- nam
  ##
  tdata$author <- as.character(tdata$author)
  tdata$type <- as.character(tdata$type)
  tdata$method <- as.character(tdata$method)
  tdata$method[tdata$method=="IV"] <- "Inverse"
  tdata$sm <- as.character(tdata$sm)
  tdata$model <- as.character(tdata$model)
  tdata$comb.fixed <- tdata$model=="Fixed"
  tdata$comb.random <- tdata$model=="Random"
  ##
  sel.oe <- tdata$method=="EXP_O_E_VAR"
  tdata$method[sel.oe] <- "Peto"
  tdata$sm[sel.oe] <- "OR"
  ##
  sel.peto <- tdata$method=="PETO" & tdata$sm=="PETO_OR"
  tdata$method[sel.peto] <- "Peto"
  tdata$sm[sel.peto] <- "OR"
  ##
  tdata$RR.cochrane <- rep(NA, length(tdata$sm))
  ##tdata$RR.cochrane[tdata$method=="MH" & tdata$sm=="RR"] <- TRUE
  tdata$RR.cochrane[tdata$sm=="RR"] <- TRUE
  ##
  diffcomp <- c(1, diff(tdata$comp.no))
  diffoutc <- c(1, diff(tdata$outcome.no))
  diffgrp  <- c(1, diff(tdata$group.no))
  ##
  sel.comp  <- diffcomp != 0
  sel.outc  <- diffoutc != 0 & !sel.comp
  sel.grp   <- diffgrp  != 0 & !sel.comp & !sel.outc
  sel.study <- !sel.comp & !sel.outc & !sel.grp
  ##
  ##print(data.frame(comp.no=tdata$comp.no,
  ##                 outcome.no=tdata$outcome.no,
  ##                 group.no=tdata$group.no,
  ##                 sel.comp=sel.comp,
  ##                 sel.outc=sel.outc,
  ##                 sel.grp=sel.grp,
  ##                 sel.study=sel.study))
  ##
  ##print(names(tdata))
  ##
  res <- list(title=title,
              comparison=
              data.frame(comp.no= selvar(tdata$comp.no, sel.comp),
                         complab= selvar(tdata$author, sel.comp)
                         ),
              outcome=
              data.frame(comp.no    = selvar(tdata$comp.no, sel.outc),
                         outcome.no = selvar(tdata$outcome.no, sel.outc),
                         ##
                         type        = selvar(tdata$type, sel.outc),
                         method      = selvar(tdata$method, sel.outc),
                         sm          = selvar(tdata$sm, sel.outc),
                         model       = selvar(tdata$model, sel.outc),
                         comb.fixed  = selvar(tdata$comb.fixed, sel.outc),
                         comb.random = selvar(tdata$comb.random, sel.outc),
                         outclab     = selvar(tdata$author, sel.outc),
                         ##
                         k = selvar(tdata$df, sel.outc)+1,
                         ##
                         event.e.pooled = selvar(tdata$event.e, sel.outc),
                         n.e.pooled     = selvar(tdata$n.e, sel.outc),
                         event.c.pooled = selvar(tdata$event.c, sel.outc),
                         n.c.pooled     = selvar(tdata$n.c, sel.outc),
                         ##
                         TE.pooled    = selvar(tdata$TE, sel.outc),
                         lower.pooled = selvar(tdata$lower, sel.outc),
                         upper.pooled = selvar(tdata$upper, sel.outc),
                         ##
                         weight.pooled = selvar(tdata$weight, sel.outc),
                         ##
                         Z.pooled       = selvar(tdata$Z, sel.outc),
                         pval.TE.pooled = selvar(tdata$pval.TE, sel.outc),
                         ##
                         Q      = selvar(tdata$Q, sel.outc),
                         pval.Q = selvar(tdata$pval.Q, sel.outc),
                         I2     = selvar(tdata$I2, sel.outc),
                         tau2   = selvar(tdata$tau2, sel.outc),
                         ##
                         Q.w      = selvar(tdata$Qint, sel.outc),
                         pval.Q.w = selvar(tdata$pval.Qint, sel.outc),
                         I2.w     = selvar(tdata$I2.Qint, sel.outc),
                         ##
                         label.e     = selvar(tdata$label.e, sel.outc),
                         label.c     = selvar(tdata$label.c, sel.outc),
                         label.left  = selvar(tdata$label.left, sel.outc),
                         label.right = selvar(tdata$label.right, sel.outc),
                         RR.cochrane = selvar(tdata$RR.cochrane, sel.outc)
                         ),
              group=
              if (sum(sel.grp) > 0 ){
                data.frame(comp.no    = selvar(tdata$comp.no, sel.grp),
                           outcome.no = selvar(tdata$outcome.no, sel.grp),
                           group.no   = selvar(tdata$group.no, sel.grp),
                           grplab     = selvar(tdata$author, sel.grp)
                           )}
              else {
                NA
              },
              study=
              data.frame(comp.no    = selvar(tdata$comp.no, sel.study),
                         outcome.no = selvar(tdata$outcome.no, sel.study),
                         group.no   = selvar(tdata$group.no, sel.study),
                         studlab    = selvar(tdata$author, sel.study),
                         year       = selvar(tdata$year, sel.study),
                         ##
                         event.e = selvar(tdata$event.e, sel.study),
                         n.e     = selvar(tdata$n.e, sel.study),
                         event.c = selvar(tdata$event.c, sel.study),
                         n.c     = selvar(tdata$n.c, sel.study),
                         ##
                         mean.e = selvar(tdata$mean.e, sel.study),
                         sd.e   = selvar(tdata$sd.e, sel.study),
                         mean.c = selvar(tdata$mean.c, sel.study),
                         sd.c   = selvar(tdata$sd.c, sel.study),
                         ##
                         O.E = selvar(tdata$O.E, sel.study),
                         V   = selvar(tdata$V, sel.study),
                         ##
                         TE   = selvar(tdata$TE, sel.study),
                         seTE = selvar(tdata$seTE, sel.study),
                         ##
                         lower.TE = selvar(tdata$lower, sel.study),
                         upper.TE = selvar(tdata$upper, sel.study),
                         weight   = selvar(tdata$weight, sel.study),
                         ##
                         order = selvar(tdata$order, sel.study)
                         )
              )
  ##
  if (is.data.frame(res$group))
    res2 <- merge(res$study, res$group,
                  by=c("comp.no", "outcome.no", "group.no"),
                  all.x=TRUE)
  else
    res2 <- res$study
  res2 <- merge(res2, res$outcome,
                by=c("comp.no", "outcome.no"))
  res2 <- merge(res2, res$comparison,
                by=c("comp.no"))
  ##
  sel.nocont <- (res2$mean.e==0 & res2$sd.e==0 &
                 res2$mean.c==0 & res2$sd.c==0)
  res2$mean.e[sel.nocont] <- NA
  res2$sd.e[sel.nocont]   <- NA
  res2$mean.c[sel.nocont] <- NA
  res2$sd.c[sel.nocont]   <- NA
  ##
  if (is.factor(res2$O.E))
    res2$O.E <- numchar(res2$O.E)
  if (is.factor(res2$V))
    res2$V <- numchar(res2$V)
  sel.noOE <- res2$O.E==0 & res2$V==0
  res2$O.E[sel.noOE] <- NA
  res2$V[sel.noOE]   <- NA
  ##
  res2$complab <- as.character(res2$complab)
  res2$outclab <- nachar(as.character(res2$outclab))
  res2$studlab <- nachar(as.character(res2$studlab))
  if (is.data.frame(res$group))
    res2$grplab <- nachar(as.character(res2$grplab))
  ##
  res2$sm     <- as.character(res2$sm)
  res2$method <- as.character(res2$method)
  res2$type   <- substring(res2$type, 1, 1)
  res2$model  <- as.character(res2$model)
  ##
  res2$label.e     <- nachar(as.character(res2$label.e))
  res2$label.c     <- nachar(as.character(res2$label.c))
  res2$label.left  <- nachar(as.character(res2$label.left))
  res2$label.right <- nachar(as.character(res2$label.right))
  ##
  if (numbers.in.labels & length(res2$comp.no)>0){
    res2$complab <- paste(res2$comp.no, " ", res2$complab, sep="")
    res2$outclab <- paste(res2$comp.no, ".", res2$outcome.no, " ",
                          res2$outclab, sep="")
    if (is.data.frame(res$group))
      res2$grplab <- paste(res2$comp.no, ".", res2$outcome.no, ".",
                           res2$group.no, " ", res2$grplab, sep="")
  }
  ##
  res2 <- res2[order(res2$comp.no,
                     res2$outcome.no,
                     res2$group.no),]
  ##
  attr(res2, "title") <- res$title
  ##
  if (length(res2$comp.no)>0){
    for (i in unique(res2$comp.no)){
      for (j in unique(res2$outcome.no[res2$comp.no==i])){
        sel2 <- res2$comp.no==i&res2$outcome.no==j
        if (unique(res2$sm[sel2])=="OTHER")
          warning("Summary measure unclear for outcome '",
                  ##paste(unique(res2$comp.no[sel2]), ".",
                  ##      unique(res2$outcome.no[sel2]), sep=""),
                  unique(res2$outclab[sel2]),
                  "'. Please use parameter 'sm' in function metacr to choose adequate summary measure.")
      }
    }
  }
  ##
  class(res2) <- c("data.frame", "rm5")
  res2
}


##  res2$event.e[res2$type == "I" & res2$event.e==0 & res2$n.e==1] <- NA
##  res2$n.e[res2$type == "I" & res2$n.e==1] <- NA
##  res2$event.c[res2$type == "I" & res2$event.c==0 & res2$n.c==1] <- NA
##  res2$n.c[res2$type == "I" & res2$n.c==1] <- NA
##  res2$mean.c[res2$type == "I" & res2$mean.c==0 & res2$sd.c==0] <- NA
##  res2$sd.c[res2$type == "I" & res2$sd.c==0] <- NA
##  res2$O.E[res2$type == "I" & res2$O.E==0 & res2$V==0] <- NA
##  res2$V[res2$type == "I" & res2$V==0] <- NA
##  ##
##  res2$mean.e[res2$type == "D" & res2$mean.e==0 & res2$sd.e==0] <- NA
##  res2$sd.e[res2$type == "D" & res2$sd.e==0] <- NA
##  res2$mean.c[res2$type == "D" & res2$mean.c==0 & res2$sd.c==0] <- NA
##  res2$sd.c[res2$type == "D" & res2$sd.c==0] <- NA
##  res2$O.E[res2$type == "D" & res2$O.E==0 & res2$V==0] <- NA
##  res2$V[res2$type == "D" & res2$V==0] <- NA
##  ##
##  res2$event.e[res2$type == "C" & res2$event.e==0] <- NA
##  res2$event.c[res2$type == "C" & res2$event.c==0] <- NA
##  res2$O.E[res2$type == "C" & res2$O.E==0 & res2$V==0] <- NA
##  res2$V[res2$type == "C" & res2$V==0] <- NA
##  ##
##  res2$mean.e[res2$type == "P" & res2$mean.e==0 & res2$sd.e==0] <- NA
##  res2$sd.e[res2$type == "P" & res2$sd.e==0] <- NA
##  res2$mean.c[res2$type == "P" & res2$mean.c==0 & res2$sd.c==0] <- NA
##  res2$sd.c[res2$type == "P" & res2$sd.c==0] <- NA
