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
  ## Check number of fields in input data file
  ##
  if (length(unique(count.fields(file, sep=sep, quote=quote, comment.char=""))) != 1)
    stop("Input file has different number of elements. Probably your export data file contains field 'Risk of bias tables'. Please create a new export data file in RevMan 5 without this field.")
  ##
  tdata <- read.table(file, header=TRUE,
                      sep=sep, quote=quote,
                      comment.char="")
  ##
  nam <- names(tdata)
  
  ##
  ## Overcome problems to import files with UTF-8, byte order mark encoding
  ## (see http://en.wikipedia.org/wiki/Byte_order_mark)
  ##
  nam[grep("Comparison.Number$", nam)] <- "Comparison.Number"
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
  ## Variable 'Totals' not included in data frame 
  ## Variable 'Study.Confidence.Interval' not included in data frame 
  ## Variable 'Total.Confidence.Interval' not included in data frame 
  ## Variable 'Test.for.subgroup.differences' not included in data frame 
  ##
  nam[nam=="Swap.event.and.non.event"] <- "swap.events"
  nam[nam=="Entered.data.are.on.log.scale..Generic.Inverse.Variance.only."] <- "logscale"
  nam[nam=="Enter.number.of.participants..Generic.Inverse.Variance..for.display.only."] <- "enter.n"
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
  ##
  nam[nam=="Effect.Estimate"] <- "TE"
  nam[nam=="SE"] <- "seTE"
  ##
  nam[nam=="CI.Start"] <- "lower"
  nam[nam=="CI.End"]   <- "upper"
  ##
  nam[nam=="Weight"]   <- "weight"
  ##
  ## Variable 'Q' is Q
  nam[nam=="P.Q."]  <- "pval.Q"
  nam[nam=="I..Q."] <- "I2"
  nam[nam=="Tau."]  <- "tau2"
  ##
  ## Variable 'Z' is Z
  nam[nam=="P.Z."]     <- "pval.TE"
  ## Variable 'Qint' is Qint
  nam[nam=="P.Qint."]  <- "pval.Qint"
  nam[nam=="I..Qint."] <- "I2.Qint"
  ## Variable 'df' is df
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
    if ((pos2-pos1) == 3)
      nam[pos1+1:2] <- c("I2", "tau2")
  }
  if (!"I2.Qint" %in% nam){
    pos3 <- seq(along=nam)[nam=="pval.Qint"]
    pos4 <- seq(along=nam)[nam=="df"]
    if ((pos4-pos3) == 2)
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
  tdata$comb.fixed  <- tdata$model == "Fixed"
  tdata$comb.random <- tdata$model == "Random"
  ##
  tdata$sm[tdata$sm=="Odds Ratio"] <- "OR"
  tdata$sm[tdata$sm=="Odds Ratio (Non-event)"] <- "OR"
  tdata$sm[tdata$sm=="Risk Ratio"] <- "RR"
  tdata$sm[tdata$sm=="Risk Difference"] <- "RD"
  tdata$sm[tdata$sm=="Mean Difference"] <- "MD"
  tdata$sm[tdata$sm=="Standardized Mean Difference"] <- "SMD"
  tdata$sm[tdata$sm=="Std. Mean Difference"] <- "SMD"
  tdata$sm[tdata$sm=="Hazard Ratio"] <- "HR"
  ##
  sel.oe <- tdata$method == "EXP_O_E_VAR"
  tdata$method[sel.oe] <- "Peto"
  tdata$sm[sel.oe] <- "OR"
  ##
  sel.peto <- tdata$method == "PETO" & tdata$sm == "PETO_OR"
  tdata$method[sel.peto] <- "Peto"
  tdata$sm[sel.peto] <- "OR"
  ##
  tdata$event.e <- as.numeric(gsub(",", "", tdata$event.e))
  tdata$n.e <- as.numeric(gsub(",", "", tdata$n.e))
  tdata$event.c <- as.numeric(gsub(",", "", tdata$event.c))
  tdata$n.c <- as.numeric(gsub(",", "", tdata$n.c))
  tdata$mean.e <- as.numeric(gsub(",", "", tdata$mean.e))
  tdata$sd.e <- as.numeric(gsub(",", "", tdata$sd.e))
  tdata$mean.c <- as.numeric(gsub(",", "", tdata$mean.c))
  tdata$sd.c <- as.numeric(gsub(",", "", tdata$sd.c))
  if (!is.null(tdata$O.E))
    tdata$O.E <- as.numeric(gsub(",", "", tdata$O.E))
  if (!is.null(tdata$V))
    tdata$V <- as.numeric(gsub(",", "", tdata$V))
  tdata$TE <- as.numeric(gsub(",", "", tdata$TE))
  ## No warning concerning infinite values for
  ## seTE, lower and upper CI bound
  oldopts <- options(warn=-1)
  tdata$seTE <- as.numeric(gsub(",", "", tdata$seTE))
  tdata$lower <- as.numeric(gsub(",", "", tdata$lower))
  tdata$upper <- as.numeric(gsub(",", "", tdata$upper))
  options(oldopts)
  ##
  tdata$weight <- as.numeric(gsub(",", "", tdata$weight))
  tdata$Q <- as.numeric(gsub(",", "", tdata$Q))
  tdata$tau2 <- as.numeric(gsub(",", "", tdata$tau2))
  tdata$Qint <- as.numeric(gsub(",", "", tdata$Qint))
  ##
  diffcomp <- c(1, diff(tdata$comp.no))
  diffoutc <- c(1, diff(tdata$outcome.no))
  diffgrp  <- c(1, diff(tdata$group.no))
  ##
  sel.comp  <- diffcomp != 0
  sel.outc  <- diffoutc != 0 & !sel.comp
  sel.grp   <- tdata$type != "" & diffgrp  != 0 & !sel.comp & !sel.outc
  sel.study <- !sel.comp & !sel.outc & !sel.grp
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
                         swap.events = selvar(tdata$swap.events, sel.outc),
                         enter.n     = selvar(tdata$enter.n, sel.outc),
                         logscale    = selvar(tdata$logscale, sel.outc),
                         ##
                         label.e     = selvar(tdata$label.e, sel.outc),
                         label.c     = selvar(tdata$label.c, sel.outc),
                         label.left  = selvar(tdata$label.left, sel.outc),
                         label.right = selvar(tdata$label.right, sel.outc)
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
  sel.nocont <- (res2$mean.e == 0 & res2$sd.e == 0 &
                 res2$mean.c == 0 & res2$sd.c == 0)
  res2$mean.e[sel.nocont] <- NA
  res2$sd.e[sel.nocont]   <- NA
  res2$mean.c[sel.nocont] <- NA
  res2$sd.c[sel.nocont]   <- NA
  ##
  if (is.factor(res2$O.E))
    res2$O.E <- numchar(res2$O.E)
  if (is.factor(res2$V))
    res2$V <- numchar(res2$V)
  sel.noOE <- res2$O.E == 0 & res2$V == 0
  res2$O.E[sel.noOE] <- NA
  res2$V[sel.noOE]   <- NA
  ##
  res2$complab <- as.character(res2$complab)
  res2$outclab <- nachar(as.character(res2$outclab))
  res2$studlab <- nachar(as.character(res2$studlab))
  if (is.data.frame(res$group)){
    res2$grplab <- nachar(as.character(res2$grplab))
    ## Replace some XML code:
    res2$grplab <- gsub("\\&lt;", "<", res2$grplab)
    res2$grplab <- gsub("\\&gt;", ">", res2$grplab)
  }
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
        sel2 <- res2$comp.no == i & res2$outcome.no == j
        if (unique(res2$sm[sel2]) == "OTHER")
          warning("Summary measure unclear for outcome '",
                  ##paste(unique(res2$comp.no[sel2]), ".",
                  ##      unique(res2$outcome.no[sel2]), sep=""),
                  unique(res2$outclab[sel2]),
                  "'. Please use argument 'sm' in function metacr to choose adequate summary measure.")
      }
    }
  }
  
  
  ##res2$TE[res2$sm=="HR"] <- log(res2$TE[res2$sm=="HR"])
  ##res2$sm[res2$sm=="log HR"] <- "HR"

  iswap.events <- charmatch(tolower(as.character(res2$swap.events)), c("yes", "no"), nomatch=NA)
  res2$swap.events <- ifelse(iswap.events==1, TRUE, ifelse(iswap.events==2, FALSE, NA))
  ##
  ienter.n <- charmatch(tolower(as.character(res2$enter.n)), c("yes", "no"), nomatch=NA)
  res2$enter.n <- ifelse(ienter.n==1, TRUE, ifelse(ienter.n==2, FALSE, NA))
  ##
  ilogscale <- charmatch(tolower(as.character(res2$logscale)), c("no", "yes"), nomatch=NA)
  res2$logscale <- ifelse(ilogscale==1, TRUE, ifelse(ilogscale==2, FALSE, NA))
  
  attr(res2, "version") <- packageDescription("meta")$Version
  
  class(res2) <- c("rm5", "data.frame")
  res2
}
