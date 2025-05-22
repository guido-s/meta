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
#' @param n2 Number of observations (second treatment).
#' @param mean1 Estimated mean (first treatment).
#' @param sd1 Standard deviation (first treatment).
#' @param mean2 Estimated mean (second treatment).
#' @param sd2 Standard deviation (second treatment).
#' @param time1 Person time at risk (first treatment).
#' @param time2 Person time at risk (second treatment).
#' @param agent1 Agent (first treatment).
#' @param agent2 Agent (second treatment).
#' @param dose1 Dose (first treatment).
#' @param dose2 Dose (second treatment).
#' @param sep.ag A character used as separator between agent and dose to
#'   create treatment labels.
#' @param data An optional data frame containing the study
#'   information.
#' @param studlab A vector with study labels (optional).
#' @param id1 Last character(s) of variable names for additional
#'   variables with group specific information for first treatment.
#' @param id2 Last character(s) of variable names for additional
#'   variables with group specific information for second treatment.
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
#' an R objected created with \code{\link{pairwise}}.
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
#'   \code{\link{metainc}}, \code{\link{pairwise}}
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
                    agent1, agent2,
                    dose1, dose2,
                    sep.ag = "*",
                    data = NULL, studlab,
                    id1 = NULL,
                    id2 = NULL,
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
  #
  dosres <-
    !(missing(agent1) | missing(agent2) |missing(dose1) | missing(dose2))
  #
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
  else if (inherits(treat1, "meta")) {
    cls <- class(treat1)
    cls <- cls[cls != "meta"]
    #
    stop("R function longarm() cannot be used with ",
         "meta-analysis object of class", if (length(cls) > 1) "es", " ",
         paste(paste0("'", cls, "'"), collapse = " - "), ".",
         call. = FALSE)
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
    ignore(missing(agent1), "agent1", cls)
    ignore(missing(agent2), "agent2", cls)
    ignore(missing(dose1), "dose1", cls)
    ignore(missing(dose2), "dose2", cls)
    ##
    event1 <- event2 <- n1 <- n2 <-
      mean1 <- mean2 <- sd1 <- sd2 <-
      time1 <- time2 <-
      agent1 <- agent2 <- dose1 <- dose2 <- NULL
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
    #
    dosres <- !is.null(treat1$agent1) & !is.null(treat1$agent2) &
      !is.null(treat1$dose1) & !is.null(treat1$dose2)
    #
    if (dosres) {
      agent1 <- treat1$agent1
      agent2 <- treat1$agent2
      dose1 <- treat1$dose1
      dose2 <- treat1$dose2
    }
    #
    data <- treat1
    treat1 <- treat1$treat1
    ##
    data$.treat1 <- treat1
    data$.treat2 <- treat2
    #
    if (dosres) {
      data$.agent1 <- agent1
      data$.agent2 <- agent2
      data$.dose1 <- dose1
      data$.dose2 <- dose2
    }
    #
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
    #
    agent1 <- catch("agent1", mc, data, sfsp)
    agent2 <- catch("agent2", mc, data, sfsp)
    dose1 <- catch("dose1", mc, data, sfsp)
    dose2 <- catch("dose2", mc, data, sfsp)
    #
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
    #
    if (type == "count") {
      k.all <-
        max(length(event1), length(time1), length(event2), length(time2))
    }
    else if (type == "binary") {
      k.all <-
        max(length(event1), length(n1), length(event2), length(n2))
    }
    else if (type == "continuous") {
      k.all <-
        max(length(n1), length(mean1), length(sd1),
            length(n2), length(mean2), length(sd2))
    }
    #
    if (missing.treat1 & missing.treat2) {
      if (dosres) {
        treat1 <- paste(agent1, dose1, sep = sep.ag)
        treat2 <- paste(agent2, dose2, sep = sep.ag)
        #
        if (length(treat1) == 1)
          treat1 <- rep(treat1, k.all)
        #
        if (length(treat2) == 1)
          treat2 <- rep(treat2, k.all)
      }
      else {
        treat1 <- rep("B", k.all)
        treat2 <- rep("A", k.all)
      }
      #
      missing.treat2 <- FALSE
    }
    else if (!missing.treat1 && !missing.treat2 &&
             (length(treat1) == 1 | length(treat2) == 1)) {
      if (length(treat1) == 1)
        treat1 <- rep(treat1, k.all)
      if (length(treat2) == 1)
        treat2 <- rep(treat2, k.all)
    }
    #
    if (dosres && (length(agent1) == 1 | length(dose1) == 1)) {
      if (length(agent1) == 1)
        agent1 <- rep(agent1, k.all)
      if (length(dose1) == 1)
        dose1 <- rep(dose1, k.all)
    }
    #
    if (dosres && (length(agent2) == 1 | length(dose2) == 1)) {
      if (length(agent2) == 1)
        agent2 <- rep(agent2, k.all)
      if (length(dose2) == 1)
        dose2<- rep(dose2, k.all)
    }
    #
    if (missing.treat2)
      stop("Argument 'treat2' mandatory.")
    #
    studlab <- replaceNULL(studlab, seq_along(treat1))
    ##
    ## Keep data set
    ##
    if (nulldata) {
      data <-
        data.frame(.studlab = studlab,
                   .treat1 = treat1, .treat2 = treat2,
                   .agent1 = if (dosres) agent1 else NA,
                   .dose1 = if (dosres) dose1 else NA,
                   .agent2 = if (dosres) agent2 else NA,
                   .dose2 = if (dosres) dose2 else NA,
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
      #
      if (dosres) {
        data$.agent1 <- agent1
        data$.agent2 <- agent2
        data$.dose1 <- dose1
        data$.dose2 <- dose2
      }
      #
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
    nam <- c("studlab", "treat",
             if (dosres) "agent", if (dosres) "dose",
             "n", "events", "nonevents")
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
    nam <- c("studlab", "treat",
             if (dosres) "agent", if (dosres) "dose",
             "n", "mean", "sd")
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
    nam <- c("studlab", "treat",
             if (dosres) "agent", if (dosres) "dose",
             if (n.given) "n", "events", "time")
  }
  ##
  dat.l$treat <- if_else(dat.l$.grp.m4 == 1, dat.l$.treat1, dat.l$.treat2)
  #
  if (dosres) {
    dat.l$agent <- if_else(dat.l$.grp.m4 == 1, dat.l$.agent1, dat.l$.agent2)
    dat.l$dose <- if_else(dat.l$.grp.m4 == 1, dat.l$.dose1, dat.l$.dose2)
  }
  #
  if (!isCol(dat.l, "studlab") & isCol(dat.l, ".studlab"))
    dat.l$studlab <- dat.l$.studlab
  ##
  if (!keep.duplicated)
    dat.l <- dat.l[!duplicated(dat.l[, nam]), ]
  ##
  ## Catch additional variables with group specific information
  ##
  if (!is.null(id1) & !is.null(id2) & !is.null(data)) {
    chklength(id1, 1)
    chklength(id2, 1)
    ##
    ext1 <- paste0(id1, "$")
    ext2 <- paste0(id2, "$")
    ##
    vars1 <- gsub(ext1, "", names(data)[grepl(ext1, names(data))])
    vars2 <- gsub(ext2, "", names(data)[grepl(ext2, names(data))])
    ##
    j <- 0
    both <- character(0)
    for (i in seq_along(vars1)) {
      if (vars1[i] %in% vars2) {
        j <- j + 1
        both[j] <- vars1[i]
      }
    }
    ##
    if (length(both) > 0) {
      bothlist <- list()
      for (i in seq_along(both)) {
        bothlist[[i]] <-
          data.frame(var1 = data[[paste0(both[i], id1)]],
                     var2 = data[[paste0(both[i], id2)]])
        names(bothlist)[[i]] <- both[i]
      }
      ##
      for (var.i in both) {
        if (!(var.i %in% names(dat.l)))
          dat.l[[var.i]] <- addvars2long(bothlist[[var.i]])
      }
    }
  }
  ##
  if (append) {
    allnames <- names(dat.l)
    other <- allnames[!(allnames %in% nam)]
    #
    nam.internal <-
      c(".studlab", ".treat1", ".treat2",
        #
        if (dosres) c(".agent1", ".agent2", ".dose1", ".dose2"),
        #
        ".event1", ".time1", ".n1", ".incr1", ".mean1", ".sd1",
        ".event2", ".time2", ".n2", ".incr2", ".mean2", ".sd2",
        #
        ".event.e", ".time.e", ".n.e", ".incr.e", ".mean.e", ".sd.e",
        ".event.c", ".time.c", ".n.c", ".incr.c", ".mean.c", ".sd.c",
        #
        ".incr",
        #
        ".subset", ".exclude", ".subgroup",
        #
        ".id.m4", ".grp.m4")
    #
    nam.internal <- nam.internal[nam.internal %in% other]
    other <- other[!(other %in% nam.internal)]
    #
    res <- dat.l[, c(nam, other, if (keep.internal) nam.internal)]
  }
  else
    res <- dat.l[, nam]
  
  
  attr(res, "type") <- type
  #
  class(res) <- unique(c("longarm", class(res)))
  #
  res
}
