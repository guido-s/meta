#' Import RevMan 4 data files (.mtv)
#' 
#' @description
#' Reads a file created with RevMan 4 and creates a data frame from
#' it.
#' 
#' @param file The name of a file to read data values from.
#' 
#' @details
#' Reads a file created with RevMan 4 (Menu: "File" - "Export" -
#' "Analysis data file...") and creates a data frame from it.
#' 
#' @return
#' A data frame containing the following components:
#' \item{comp.no}{Comparison number.}
#' \item{outcome.no}{Outcome number.}
#' \item{group.no}{Group number.}
#' \item{studlab}{Study label.}
#' \item{year}{Year of publication.}
#' \item{event.e}{Number of events in experimental group.}
#' \item{n.e}{Number of observations in experimental group.}
#' \item{event.c}{Number of events in control group.}
#' \item{n.c}{Number of observations in control group.}
#' \item{mean.e}{Estimated mean in experimental group.}
#' \item{sd.e}{Standard deviation in experimental group.}
#' \item{mean.c}{Estimated mean in control group.}
#' \item{sd.c}{Standard deviation in control group.}
#' \item{O.E}{Observed minus expected (IPD analysis).}
#' \item{V}{Variance of \code{O.E} (IPD analysis).}
#' \item{order}{Ordering of studies.}
#' \item{conceal}{Concealment of treatment allocation.}
#' \item{grplab}{Group label.}
#' \item{type}{Type of outcome. D = dichotomous, C = continuous, P =
#'   IPD.}
#' \item{outclab}{Outcome label.}
#' \item{graph.exp}{Graph label for experimental group.}
#' \item{graph.cont}{Graph label for control group.}
#' \item{label.exp}{Label for experimental group.}
#' \item{label.cont}{Label for control group.}
#' \item{complab}{Comparison label.}
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}, \code{\link{metacont}},
#'   \code{\link{metagen}}
#' 
#' @references
#' \emph{Review Manager (RevMan)} [Computer program]. Version 4.2.
#' Copenhagen: The Nordic Cochrane Centre, The Cochrane Collaboration, 2003
#' 
#' @keywords datagen
#' 
#' @examples
#' # Locate MTV-data file "FLEISS1993.MTV" in sub-directory of R package
#' # meta
#' #
#' filename <- system.file("extdata/FLEISS1993.MTV", package = "meta")
#' fleiss1933.cc <- read.mtv(filename)
#' 
#' # Same result as R Command example(Fleiss1993bin):
#' #
#' metabin(event.e, n.e, event.c, n.c,
#'   data = fleiss1933.cc, subset = type == "D",
#'   studlab = paste(studlab, year))
#' 
#' # Same result: example(Fleiss1993cont)
#' #
#' metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c,
#'   data = fleiss1933.cc, subset = type == "C",
#'   studlab = paste(studlab, year))
#' 
#' @export read.mtv


read.mtv <- function(file) {
  ##
  line <- scan(file,
               what = "character",
               blank.lines.skip = FALSE,
               comment.char = "",
               sep = "\n")
  ##
  ##line <- line[substring(line, 1, 10) != "START DATA"]
  ##line <- line[substring(line, 1,  8) != "END DATA"]
  ##line <- line[substring(line, 1,  8) != "STUDIES:"]
  ##
  sel <- substring(line, 1, 9) == "METAVIEW:"
  title <- rmSpace(substring(line[sel], 10), end = TRUE)
  ##line <- line[!sel]
  ##
  ##
  comp  <- line[substring(line, 1, 1) == 3]
  outc  <- line[substring(line, 1, 1) == 2]
  sgrp  <- line[substring(line, 1, 1) == 1]
  study <- line[substring(line, 1, 1) == 0]
  ##
  ##
  res <- list(title = title,
              ##
              comparison =
                data.frame(comp.no = I(substring(comp, 3, 4)),
                           complab = I(rmSpace(substring(comp, 9), end = TRUE))
                           ),
              ##
              outcome =
                data.frame(#id         = substring(outc, 5, 9),
                  comp.no    = I(substring(outc, 5, 6)),
                  outcome.no = I(substring(outc, 8, 9)),
                  type       = I(substring(outc, 3, 3)),
                  totals     = I(substring(outc, 11, 11)),
                  outclab    = I(rmSpace(substring(outc,  26, 132), end = TRUE)),
                  graph.exp  = I(rmSpace(substring(outc, 134, 153), end = TRUE)),
                  graph.cont = I(rmSpace(substring(outc, 155, 174), end = TRUE)),
                  label.exp  = I(rmSpace(substring(outc, 176, 195), end = TRUE)),
                  label.cont = I(rmSpace(substring(outc, 197, 216), end = TRUE))),
              ##  details    = I(rmSpace(substring(outc,   5,  24), end = TRUE))),
              ##
              group =
                data.frame(#id         = I(substring(sgrp, 3, 10)),
                  comp.no    = I(substring(sgrp, 3,  4)),
                  outcome.no = I(substring(sgrp, 6,  7)),
                  group.no   = I(substring(sgrp, 9, 10)),
                  grplab     = I(rmSpace(substring(sgrp, 15), end = TRUE))),
              ##
              study =
                data.frame(#id         = I(substring(study, 3, 10)),
                  comp.no    = I(substring(study, 3,  4)),
                  outcome.no = I(substring(study, 6,  7)),
                  group.no   = I(substring(study, 9, 10)),
                  ##
                  studlab    = I(rmSpace(substring(study, 12, 31), end = TRUE)),
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
    res2 <- merge(res$study, res$group, by = c("comp.no", "outcome.no", "group.no"), all.x = TRUE)
  else
    res2 <- res$study
  res2 <- merge(res2, res$outcome, by = c("comp.no", "outcome.no"))
  res2 <- merge(res2, res$comparison, by = c("comp.no"))
  ##
  res2$event.e[res2$type == "I" & res2$event.e == 0 & res2$n.e == 1] <- NA
  res2$n.e[res2$type == "I" & res2$n.e == 1] <- NA
  res2$event.c[res2$type == "I" & res2$event.c == 0 & res2$n.c == 1] <- NA
  res2$n.c[res2$type == "I" & res2$n.c == 1] <- NA
  res2$mean.c[res2$type == "I" & res2$mean.c == 0 & res2$sd.c == 0] <- NA
  res2$sd.c[res2$type == "I" & res2$sd.c == 0] <- NA
  res2$O.E[res2$type == "I" & res2$O.E == 0 & res2$V == 0] <- NA
  res2$V[res2$type == "I" & res2$V == 0] <- NA
  ##
  res2$mean.e[res2$type == "D" & res2$mean.e == 0 & res2$sd.e == 0] <- NA
  res2$sd.e[res2$type == "D" & res2$sd.e == 0] <- NA
  res2$mean.c[res2$type == "D" & res2$mean.c == 0 & res2$sd.c == 0] <- NA
  res2$sd.c[res2$type == "D" & res2$sd.c == 0] <- NA
  res2$O.E[res2$type == "D" & res2$O.E == 0 & res2$V == 0] <- NA
  res2$V[res2$type == "D" & res2$V == 0] <- NA
  ##
  res2$event.e[res2$type == "C" & res2$event.e == 0] <- NA
  res2$event.c[res2$type == "C" & res2$event.c == 0] <- NA
  res2$O.E[res2$type == "C" & res2$O.E == 0 & res2$V == 0] <- NA
  res2$V[res2$type == "C" & res2$V == 0] <- NA
  ##
  res2$mean.e[res2$type == "P" & res2$mean.e == 0 & res2$sd.e == 0] <- NA
  res2$sd.e[res2$type == "P" & res2$sd.e == 0] <- NA
  res2$mean.c[res2$type == "P" & res2$mean.c == 0 & res2$sd.c == 0] <- NA
  res2$sd.c[res2$type == "P" & res2$sd.c == 0] <- NA
  ##
  attr(res2, "title") <- res$title
  
  attr(res2, "version") <- packageDescription("meta")$Version
  
  res2
}
