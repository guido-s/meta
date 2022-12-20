#' Transform data from pairwise comparisons to long arm-based format
#' 
#' @description
#' This function transforms data from pairwise comparisons to a long
#' arm-based format, i.e., two rows for a pairwise comparison.
#' 
#' @param treat1 Either label for first treatment or a meta-analysis
#'   or pairwise object (see Details).
#' @param treat2 Label for second treatment.
#' @param event1 Number of events (first treatment).
#' @param n1 Number of observations (first treatment).
#' @param event2 Number of events (second treatment).
#' @param n2 Number of observations (second treatment)
#' @param mean1 Estimated mean (first treatment).
#' @param sd1 Standard deviation (first treatment).
#' @param mean2 Estimated mean (second treatment).
#' @param sd2 Standard deviation (second treatment).
#' @param time1 Person time at risk (first treatment)
#' @param time2 Person time at risk (second treatment)
#' @param data An optional data frame containing the study
#'   information.
#' @param studlab A vector with study labels (optional).
#' @param append A logical indicating if data frame provided in
#'   argument 'data' should be returned.
#' @param keep.duplicated A logical indicating if duplicated rows
#'   should be returned (see Details).
#' @param keep.internal A logical indicating if variables generated
#'   internally should be returned (typically only relevant for data
#'   checking).
#' 
#' @details
#' This function transforms data given as one pairwise comparison per
#' row to a long arm-based format with one row per treatment arm. The
#' long arm-based format is, for example, the required input format
#' for WinBUGS.
#' 
#' The function can be used to transform data with a binary,
#' continuous or count outcome. The corresponding meta-analysis
#' functions are \code{\link{metabin}}, \code{\link{metacont}} and
#' \code{\link{metainc}}. Accordingly, a meta-analysis object created
#' with one of these functions can be provided as argument
#' \code{treat1}. It is also possible to use the longarm function with
#' an R objected created with \code{\link[netmeta]{pairwise}} from R
#' package \bold{netmeta}.
#'
#' Otherwise, arguments \code{treat1} and \code{treat2} are mandatory
#' to identify the individual treatments and, depending on the
#' outcome, the following additional arguments are mandatory:
#' 
#' \itemize{
#' \item event1, n1, event2, n2 (binary outcome);
#' \item n1, mean1, sd1, n2, mean2, sd2 (continuous outcome);
#' \item time1, n1, time2, n2 (count outcome).
#' }
#' 
#' Argument \code{studlab} must be provided if several pairwise
#' comparisons come from a single study with more than two treatments.
#'
#' The following variables will be returned:
#'
#' \tabular{rl}{
#' \bold{\emph{studlab}}\tab study label \cr
#' \bold{\emph{treat}}\tab treatment label \cr
#' \bold{\emph{n}}\tab group sample size (count outcome only if provided) \cr
#' \bold{\emph{events}}\tab number of events (binary or count outcome) \cr
#' \bold{\emph{nonevents}}\tab number of non-events (binary outcome) \cr
#' \bold{\emph{mean}}\tab estimated mean (continuous outcome) \cr
#' \bold{\emph{sd}}\tab standard deviation (continuous outcome) \cr
#' \bold{\emph{time}}\tab person time at risk (count outcome) \cr
#' }
#'
#' In addition, the data set provided in argument \code{data} will be
#' returned if argument \code{append = TRUE} (default).
#' 
#' Argument \code{keep.duplicated} can be used to keep duplicated rows
#' from the data set. Duplicated rows can occur, for example, in a
#' three-arm study comparing treatments A and B with placebo. In this
#' situation, the placebo arm will be returned twice in the data set
#' in long arm-based format if \code{keep.duplicated = TRUE}. By
#' default, duplicated rows with not be kept in the data set.
#'
#' @note
#' R function \code{\link[metafor]{to.long}} from R package
#' \bold{metafor} is called internally.
#' 
#' @return
#' A data frame in long arm-based format.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}, \code{\link{metacont}},
#'   \code{\link{metainc}}, \code{\link[netmeta]{pairwise}}
#' 
#' @keywords datagen
#' 
#' @examples
#' # Artificial example with three studies
#' m <- metabin(1:3, 100:102, 4:6, 200:202, studlab = LETTERS[1:3])
#' # Transform data to long arm-based format
#' longarm(m)
#' # Keep internal variables
#' longarm(m, keep.internal = TRUE)
#' 
#' @export longarm


longarm <- function(treat1, treat2,
                    event1, n1, event2, n2,
                    mean1, sd1, mean2, sd2,
                    time1, time2,
                    data = NULL, studlab,
                    append = TRUE,
                    keep.duplicated = FALSE,
                    keep.internal = FALSE) {  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chklogical(append)
  chklogical(keep.duplicated)
  chklogical(keep.internal)
  
  
  ##
  ##
  ## (2) Read data
  ##
  ##
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  if (nulldata)
    data <- sfsp
  ##
  ## Catch treat1 and check for pairwise object
  ##
  missing.treat1 <- missing(treat1)
  missing.treat2 <- missing(treat2)
  ##
  if (missing.treat1)
    stop("Argument 'treat1' mandatory.")
  ##
  treat1 <- catch("treat1", mc, data, sfsp)
  ##
  ignore <- function(miss, name, func)
    if (!miss)
      warning("Argument '", name,
              "' ignored for object created with ",
              func, "().",
              call. = FALSE)
  ##
  if (inherits(treat1, c("metabin", "metacont", "metainc"))) {
    cls <- class(treat1)[1]
    ##
    ignore(missing.treat2, "treat2", cls)
    ignore(missing(event1), "event1", cls)
    ignore(missing(event2), "event2", cls)
    ignore(missing(n1), "n1", cls)
    ignore(missing(n2), "n2", cls)
    ignore(missing(mean1), "mean1", cls)
    ignore(missing(mean2), "mean2", cls)
    ignore(missing(sd1), "sd1", cls)
    ignore(missing(sd2), "sd2", cls)
    ignore(missing(time1), "time1", cls)
    ignore(missing(time2), "time2", cls)
    ##
    event1 <- event2 <- n1 <- n2 <-
      mean1 <- mean2 <- sd1 <- sd2 <-
        time1 <- time2 <- NULL
    ##
    if (!is.null(treat1$data)) {
      data <- treat1$data
      if (isCol(data, ".subset"))
        data <- data[data$.subset, , keep = TRUE]
      if (isCol(data, ".exclude"))
        data <- data[!data$.exclude, , keep = TRUE]
    }
    else
      data <- as.data.frame(treat1)
    ##
    data$.treat1 <- rep_len(treat1$label.e, nrow(data))
    data$.treat2 <- rep_len(treat1$label.c, nrow(data))
    ##
    if (cls == "metabin") {     
      event1 <- data$.event.e
      event2 <- data$.event.c
      n1 <- data$.n.e
      n2 <- data$.n.c
      ##
      type <- "binary"
    }
    else if (cls == "metacont") {
      n1 <- data$.n.e
      n2 <- data$.n.c
      mean1 <- data$.mean.e
      mean2 <- data$.mean.c
      sd1 <- data$.sd.e
      sd2 <- data$.sd.c
      ##
      type <- "continuous"
    }
    else if (cls == "metainc") {
      event1 <- data$.event.e
      event2 <- data$.event.c
      time1 <- data$.time.e
      time2 <- data$.time.c
      n1 <- data$.n.e
      n2 <- data$.n.c
      ##
      type <- "count"
    }
    ##
    studlab <- data$.studlab
    treat2 <- data$.treat2
    treat1 <- data$.treat1
  }
  else if (is.data.frame(treat1) & !is.null(attr(treat1, "pairwise"))) {
    if (!nulldata)
      warning("Argument 'data' ignored for object created with pairwise().",
              call. = FALSE)
    ##
    ignore(missing.treat2, "treat2", "pairwise")
    ignore(missing(event1), "event1", "pairwise")
    ignore(missing(event2), "event2", "pairwise")
    ignore(missing(n1), "n1", "pairwise")
    ignore(missing(n2), "n2", "pairwise")
    ignore(missing(mean1), "mean1", "pairwise")
    ignore(missing(mean2), "mean2", "pairwise")
    ignore(missing(sd1), "sd1", "pairwise")
    ignore(missing(sd2), "sd2", "pairwise")
    ignore(missing(time1), "time1", "pairwise")
    ignore(missing(time2), "time2", "pairwise")
    ##
    event1 <- event2 <- n1 <- n2 <-
      mean1 <- mean2 <- sd1 <- sd2 <-
        time1 <- time2 <- NULL
    ##
    if (attr(treat1, "type") == "binary") {
      event1 <- treat1$event1
      event2 <- treat1$event2
      n1 <- treat1$n1
      n2 <- treat1$n2
      ##
      type <- "binary"
    }
    else if (attr(treat1, "type") == "continuous") {
      n1 <- treat1$n1
      n2 <- treat1$n2
      mean1 <- treat1$mean1
      mean2 <- treat1$mean2
      sd1 <- treat1$sd1
      sd2 <- treat1$sd2
      ##
      type <- "continuous"
    }
    else if (attr(treat1, "type") == "count") {
      event1 <- treat1$event1
      event2 <- treat1$event2
      time1 <- treat1$time1
      time2 <- treat1$time2
      n1 <- treat1$n1
      n2 <- treat1$n2
      ##
      type <- "count"
    }
    else
      stop("Function cannot be used with generic outcome.")
    ##
    studlab <- treat1$studlab
    treat2 <- treat1$treat2
    data <- treat1
    treat1 <- treat1$treat1
    ##
    data$.treat1 <- treat1
    data$.treat2 <- treat2
    data$.event1 <- event1
    data$.event2 <- event2
    data$.n1 <- n1
    data$.n2 <- n2
    data$.mean1 <- mean1
    data$.mean2 <- mean2
    data$.sd1 <- sd1
    data$.sd2 <- sd2
    data$.time1 <- time1
    data$.time2 <- time2
    data$.studlab <- studlab
  }
  else {
    ##
    ## Catch studlab, treat2, event1, event2, n1, n2, mean1, mean2,
    ## sd1, sd2, time1, time2 from data:
    ##
    studlab <- catch("studlab", mc, data, sfsp)
    treat2 <- catch("treat2", mc, data, sfsp)
    event1 <- catch("event1", mc, data, sfsp)
    event2 <- catch("event2", mc, data, sfsp)
    n1 <- catch("n1", mc, data, sfsp)
    n2 <- catch("n2", mc, data, sfsp)
    mean1 <- catch("mean1", mc, data, sfsp)
    mean2 <- catch("mean2", mc, data, sfsp)
    sd1 <- catch("sd1", mc, data, sfsp)
    sd2 <- catch("sd2", mc, data, sfsp)
    time1 <- catch("time1", mc, data, sfsp)
    time2 <- catch("time2", mc, data, sfsp)
    if (missing.treat2)
      stop("Argument 'treat2' mandatory.")
    ##
    if (!is.null(event1))
      chknumeric(event1)
    if (!is.null(event2))
      chknumeric(event2)
    ##
    if (!is.null(n1))
      chknumeric(n1)
    if (!is.null(n2))
      chknumeric(n2)
    ##
    if (!is.null(mean1))
      chknumeric(mean1)
    if (!is.null(mean2))
      chknumeric(mean2)
    ##
    if (!is.null(sd1))
      chknumeric(sd1)
    if (!is.null(sd2))
      chknumeric(sd2)
    ##
    if (!is.null(time1))
      chknumeric(time1)
    if (!is.null(time2))
      chknumeric(time2)
    ##
    if (!is.null(event1) & !is.null(time1) &
        !is.null(event2) & !is.null(time2) &
        is.null(mean1) & is.null(sd1) &
        is.null(mean2) & is.null(sd2))
      type <- "count"
    else if (!is.null(event1) & !is.null(n1) &
             !is.null(event2) & !is.null(n2) &
             is.null(mean1) & is.null(sd1) &
             is.null(mean2) & is.null(sd2))
      type <- "binary"
    else if (!is.null(n1) & !is.null(n2) &
             !is.null(mean1) & !is.null(mean2) &
             !is.null(sd1) & !is.null(sd2))
      type <- "continuous"
    else
      stop("Type of outcome unclear. Please provide the necessary ",
           "information:\n  ",
           "- event1, n1, event2, n2 (binary outcome)\n  ",
           "- n1, mean1, sd1, n2, mean2, sd2 (continuous outcome)\n  ",
           "- event1, time1, event2, time2 (incidence rates).")
    ##
    if (nulldata) {
      data <-
        data.frame(.studlab = studlab,
                   .treat1 = treat1, .treat2 = treat2,
                   .n1 = n1, .n2 = n2)
      data$.mean1 <- mean1
      data$.mean2 <- mean2
      data$.sd1 <- sd1
      data$.sd2 <- sd2
      data$.time1 <- time1
      data$.time2 <- time2
    }
    else {
      data$.treat1 <- treat1
      data$.treat2 <- treat2
      data$.event1 <- event1
      data$.event2 <- event2
      data$.n1 <- n1
      data$.n2 <- n2
      data$.mean1 <- mean1
      data$.mean2 <- mean2
      data$.sd1 <- sd1
      data$.sd2 <- sd2
      data$.time1 <- time1
      data$.time2 <- time2
      data$.studlab <- studlab
    }
  }
  
  
  ##
  ##
  ## (3) Check length of variables
  ##
  ##
  k.all <- length(treat1)
  chklength(treat2, k.all, name = "treat1")
  ##
  if (!is.null(event1))
    chklength(event1, k.all, name = "treat1")
  if (!is.null(event2))
    chklength(event2, k.all, name = "treat1")
  ##
  if (!is.null(n1))
      chklength(n1, k.all, name = "treat1")
  if (!is.null(n2))
    chklength(n2, k.all, name = "treat1")
  ##
  if (!is.null(mean1))
    chklength(mean1, k.all, name = "treat1")
  if (!is.null(mean2))
    chklength(mean2, k.all, name = "treat1")
  ##
  if (!is.null(sd1))
    chklength(sd1, k.all, name = "treat1")
  if (!is.null(sd2))
    chklength(sd2, k.all, name = "treat1")
  ##
  if (!is.null(time1))
    chklength(time1, k.all, name = "treat1")
  if (!is.null(time2))
    chklength(time2, k.all, name = "treat1")
  ##
  if (is.null(studlab)) {
    studlab <- seq_len(k.all)
    data$.studlab <- studlab
  }
  
  
  ##
  ##
  ## (4) Create data set in long arm-based format
  ##
  ##
  if (type == "binary") {
    dat.l <-
      to.long("RD",
              ai = event1, bi = n1 - event1,
              ci = event2, di = n2 - event2,
              slab = studlab,
              data = data,
              var.names = c(".id.m4", ".grp.m4",
                            "events", "nonevents"))
    dat.l$n <- dat.l$events + dat.l$nonevents
    nam <- c("studlab", "treat", "n", "events", "nonevents")
  }
  else if (type == "continuous") {
    dat.l <-
      to.long("MD",
              n1i = n1, m1i = mean1, sd1i = sd1,
              n2i = n2, m2i = mean2, sd2i = sd2,
              slab = studlab,
              data = data,
              var.names = c(".id.m4", ".grp.m4",
                            "mean", "sd", "n"))
    nam <- c("studlab", "treat", "n", "mean", "sd")
  }
  else if (type == "count") {
    dat.l <-
      to.long("IRD",
              x1i = event1, t1i = time1,
              x2i = event2, t2i = time2,
              slab = studlab,
              data = data,
              var.names = c(".id.m4", ".grp.m4", "events", "time"))
    ##
    n.given <- FALSE
    if (!is.null(dat.l$.n1) & !is.null(dat.l$.n2)) {
      n.given <- TRUE
      dat.l$n <- ifelse(dat.l$.grp.m4 == 1, dat.l$.n1, dat.l$.n2)
    }
    else if (!is.null(dat.l$.n.e) & !is.null(dat.l$.n.c)) {
      n.given <- TRUE
      dat.l$n <- ifelse(dat.l$.grp.m4 == 1, dat.l$.n.e, dat.l$.n.c)
    } 
    ##
    nam <- c("studlab", "treat", if (n.given) "n", "events", "time")
  }
  ##
  dat.l$treat <- ifelse(dat.l$.grp.m4 == 1, dat.l$.treat1, dat.l$.treat2)
  if (!isCol(dat.l, "studlab") & isCol(dat.l, ".studlab"))
    dat.l$studlab <- dat.l$.studlab
  ##
  if (!keep.duplicated)
    dat.l <- dat.l[!duplicated(dat.l[, nam]), ]
  ##
  if (append) {
    allnames <- names(dat.l)
    other <- allnames[!(allnames %in% nam)]
    if (!keep.internal) {
      nam.internal <-
        c(".id.m4", ".grp.m4",
          ".treat1", ".treat2",
          ".event1", ".event2", ".n1", ".n2",
          ".mean1", ".mean2", ".sd1", ".sd2",
          ".time1", ".time2",
          ".studlab",
          ".event.e", ".event.c", ".n.e", ".n.c",
          ".mean.e", ".mean.c", ".sd.e", ".sd.c",
          ".time.e", ".time.c",
          ".incr", ".subset", ".exclude", ".subgroup")
      other <- other[!(other %in% nam.internal)]
    }
    res <- dat.l[, c(nam, other)]
  }
  else
    res <- dat.l[, nam]
  
  
  attr(res, "type") <- type
  attr(res, "longarm") <- TRUE
  ##
  res
}
