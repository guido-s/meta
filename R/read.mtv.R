read.mtv <- function(file){
  ##
  line <- scan(file,
               what="character",
               blank.lines.skip=FALSE,
               comment.char="",
               sep="\n")
  ##
  ##line <- line[substring(line, 1, 10) != "START DATA"]
  ##line <- line[substring(line, 1,  8) != "END DATA"]
  ##line <- line[substring(line, 1,  8) != "STUDIES:"]
  ##
  sel <- substring(line, 1, 9) == "METAVIEW:"
  title <- rmSpace(substring(line[sel], 10), end=TRUE)
  ##line <- line[!sel]
  ##
  ##
  comp  <- line[substring(line, 1, 1) == 3]
  outc  <- line[substring(line, 1, 1) == 2]
  sgrp  <- line[substring(line, 1, 1) == 1]
  study <- line[substring(line, 1, 1) == 0]
  ##
  ##
  res <- list(title=title,
              ##
              comparison=
              data.frame(comp.no = I(substring(comp, 3, 4)),
                         complab = I(rmSpace(substring(comp, 9), end=TRUE))
                         ),
              ##
              outcome=
              data.frame(#id         = substring(outc, 5, 9),
                         comp.no    = I(substring(outc, 5, 6)),
                         outcome.no = I(substring(outc, 8, 9)),
                         type       = I(substring(outc, 3, 3)),
                         totals     = I(substring(outc, 11, 11)),
                         outclab    = I(rmSpace(substring(outc,  26, 132), end=TRUE)),
                         graph.exp  = I(rmSpace(substring(outc, 134, 153), end=TRUE)),
                         graph.cont = I(rmSpace(substring(outc, 155, 174), end=TRUE)),
                         label.exp  = I(rmSpace(substring(outc, 176, 195), end=TRUE)),
                         label.cont = I(rmSpace(substring(outc, 197, 216), end=TRUE))),
                         #details    = I(rmSpace(substring(outc,   5,  24), end=TRUE))),
              ##
              group=
              data.frame(#id         = I(substring(sgrp, 3, 10)),
                         comp.no    = I(substring(sgrp, 3,  4)),
                         outcome.no = I(substring(sgrp, 6,  7)),
                         group.no   = I(substring(sgrp, 9, 10)),
                         grplab     = I(rmSpace(substring(sgrp, 15), end=TRUE))),
              ##
              study=
              data.frame(#id         = I(substring(study, 3, 10)),
                         comp.no    = I(substring(study, 3,  4)),
                         outcome.no = I(substring(study, 6,  7)),
                         group.no   = I(substring(study, 9, 10)),
                         ##
                         studlab    = I(rmSpace(substring(study, 12, 31), end=TRUE)),
                         year       = as.numeric(substring(study, 33, 36)),
                         ##
                         event.e    = as.numeric(substring(study, 39, 44)),
                         n.e        = as.numeric(substring(study, 46, 51)),
                         event.c    = as.numeric(substring(study, 73, 78)),
                         n.c        = as.numeric(substring(study, 80, 85)),
                         ##
                         mean.e     = as.numeric(substring(study, 53,  61)),
                         sd.e       = as.numeric(substring(study, 63,  71)),
                         mean.c     = as.numeric(substring(study, 87,  95)),
                         sd.c       = as.numeric(substring(study, 97, 105)),
                         ##
                         O.E        = as.numeric(substring(study, 107, 115)),
                         V          = as.numeric(substring(study, 117, 125)),
                         ##
                         order      = as.numeric(substring(study, 127, 130)),
                         conceal    = I(substring(study, 132, 132)))
              )
  ##
  if (dim(res$group)[[1]] > 0)
    res2 <- merge(res$study, res$group, by=c("comp.no", "outcome.no", "group.no"), all.x=TRUE)
  else
    res2 <- res$study
  res2 <- merge(res2, res$outcome, by=c("comp.no", "outcome.no"))
  res2 <- merge(res2, res$comparison, by=c("comp.no"))
  ##
  res2$event.e[res2$type == "I" & res2$event.e==0 & res2$n.e==1] <- NA
  res2$n.e[res2$type == "I" & res2$n.e==1] <- NA
  res2$event.c[res2$type == "I" & res2$event.c==0 & res2$n.c==1] <- NA
  res2$n.c[res2$type == "I" & res2$n.c==1] <- NA
  res2$mean.c[res2$type == "I" & res2$mean.c==0 & res2$sd.c==0] <- NA
  res2$sd.c[res2$type == "I" & res2$sd.c==0] <- NA
  res2$O.E[res2$type == "I" & res2$O.E==0 & res2$V==0] <- NA
  res2$V[res2$type == "I" & res2$V==0] <- NA
  ##
  res2$mean.e[res2$type == "D" & res2$mean.e==0 & res2$sd.e==0] <- NA
  res2$sd.e[res2$type == "D" & res2$sd.e==0] <- NA
  res2$mean.c[res2$type == "D" & res2$mean.c==0 & res2$sd.c==0] <- NA
  res2$sd.c[res2$type == "D" & res2$sd.c==0] <- NA
  res2$O.E[res2$type == "D" & res2$O.E==0 & res2$V==0] <- NA
  res2$V[res2$type == "D" & res2$V==0] <- NA
  ##
  res2$event.e[res2$type == "C" & res2$event.e==0] <- NA
  res2$event.c[res2$type == "C" & res2$event.c==0] <- NA
  res2$O.E[res2$type == "C" & res2$O.E==0 & res2$V==0] <- NA
  res2$V[res2$type == "C" & res2$V==0] <- NA
  ##
  res2$mean.e[res2$type == "P" & res2$mean.e==0 & res2$sd.e==0] <- NA
  res2$sd.e[res2$type == "P" & res2$sd.e==0] <- NA
  res2$mean.c[res2$type == "P" & res2$mean.c==0 & res2$sd.c==0] <- NA
  res2$sd.c[res2$type == "P" & res2$sd.c==0] <- NA
  ##
  attr(res2, "title") <- res$title
  
  attr(res2, "version") <- packageDescription("meta")$Version
  
  res2
}
