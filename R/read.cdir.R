#' Import data of Cochrane intervention review
#' 
#' @description
#' Reads Cochrane data package (version 1) of a Cochrane intervention
#' review and creates a data frame from it.
#' 
#' @aliases read.cdir
#' 
#' @param file The name of a file to read data values from.
#' @param title Title of Cochrane review.
#' @param exdir The directory to extract files to (the equivalent of
#'   ‘unzip -d’). It will be created if necessary.
#' @param numbers.in.labels A logical indicating whether comparison
#'   number and outcome number should be printed at the beginning of
#'   the comparison (argument \code{complab}) and outcome label
#'   (argument \code{outclab}); this is the default in RevMan Web.
#' @param rob A logical indicating whether risk of bias (RoB)
#'   assessment should be considered in meta-analyses.
#' @param tool Risk of bias (RoB) tool.
#' @param categories Possible RoB categories.
#' @param col Colours for RoB categories.
#' @param symbols Corresponding symbols for RoB categories.
#' @param keep.orig A logical indicating whether to return the
#'   original data files.
#' @param \dots Additional arguments (passed on to
#'   \code{\link[utils]{unzip}})
#' @param x An object of class "cdir".
#' 
#' @details
#' RevMan Web is the current software used for preparing and maintaining
#' Cochrane reviews. RevMan Web includes the ability to write systematic
#' reviews of interventions or diagnostic test accuracy reviews.
#' 
#' This function provides the ability to read the Cochrane data
#' package from a Cochrane intervention review created with RevMan
#' Web. The ZIP-file is extracted with \code{\link[utils]{unzip}}.
#' 
#' Argument \code{title} can be used to overwrite the title of the
#' Cochrane review.
#'
#' Information on the risk of bias (RoB) assessment can be provided
#' with arguments \code{tool}, \code{categories}, \code{col} and
#' \code{symbols}. This is only useful if (i) all outcomes are based
#' on the same RoB categories and (ii) an overall RoB assessment has
#' not been done. If no overall RoB assessment was conducted, R
#' function \code{\link{metacr}} can be used to provide the RoB
#' information for a single outcome. R function \code{\link{rob}} is
#' the most flexible way to add RoB information to a meta-analysis
#' object.
#' 
#' \subsection{Creation of Cochrane data package}{
#' 
#' Two possible ways exist to create the ZIP-file.
#'
#' In RevMan Web, press the "Export" button at the bottom of the
#' \emph{Default view} website. After a couple of seconds, the data package will
#' be shown at the bottom of the \emph{Default view} website under "Downloads".
#'
#' In the Cochrane Library, press on "Download statistical data" in the
#' Contents menu to download an rm5-file. This file can be converted to a data
#' package in RevMan Web using \emph{Help} - \emph{Convert a RevMan 5 file}.
#' }
#' 
#' @return
#' A list consisting of a data frame 'data' with the study data and
#' (if available) a data frame 'rob' with information on the risk of
#' bias assessment. If \code{keep.orig = TRUE}, an additional list
#' 'orig' is returned containing elements 'settings', 'datarows',
#' 'subgroup' and 'rob' (if available).
#'
#' The data frame 'data' contains the following variables:
#' 
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
#' \item{TE, seTE}{Estimated treatment effect and standard error of
#'   individual studies.}
#' \item{lower, upper}{Lower and upper limit of 95\% confidence
#'   interval for treatment effect in individual studies.}
#' \item{weight}{Weight of individual studies (according to
#'   meta-analytical method used in respective meta-analysis - see
#'   details).}
#' \item{order}{Ordering of studies.}
#' \item{grplab}{Group label.}
#' \item{type}{Type of outcome. D = dichotomous, C = continuous, P =
#'   IPD.}
#' \item{method}{A character string indicating which method has been
#'   used for pooling of studies. One of \code{"Inverse"},
#'   \code{"MH"}, or \code{"Peto"}.}
#' \item{sm}{A character string indicating which summary measure has
#'   been used for pooling of studies.}
#' \item{model}{A character string indicating which meta-analytical
#'   model has been used (either \code{"Fixed"} or \code{"Random"}).}
#' \item{common}{A logical indicating whether common effect
#'   meta-analysis has been used in respective meta-analysis (see
#'   details).}
#' \item{random}{A logical indicating whether random effects
#'   meta-analysis has been used in respective meta-analysis (see
#'   details).}
#' \item{title}{Title of Cochrane review.}
#' \item{complab}{Comparison label.}
#' \item{outclab}{Outcome label.}
#' \item{label.e}{Label for experimental group.}
#' \item{label.c}{Label for control group.}
#' \item{label.left}{Graph label on left side of forest plot.}
#' \item{label.right}{Graph label on right side of forest plot.}
#'
#' The data frame 'rob' contains the following variables:
#'
#' \item{studlab}{Study label.}
#' \item{D1, D2, \dots}{Risk of bias domain 1, 2, \dots}
#' \item{D1.details, D2.details, \dots}{Details on risk of bias domain
#'   1, 2, \dots}
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metacr}}, \code{\link{rob}},
#'   \code{\link{read.rm5}}
#' 
#' @references
#' 
#' \url{https://documentation.cochrane.org/revman-kb/data-package-user-guide-243761660.html}
#' 
#' \url{https://documentation.cochrane.org/revman-kb/data-package-specification-249561249.html}
#' 
#' @keywords datagen
#' 
#' @examples
#' # Locate file "Fleiss1993.zip" with Cochrane data package in
#' # sub-directory of R package meta
#' #
#' filename <- system.file("extdata/Fleiss1993.zip", package = "meta")
#' Fleiss1993_CR <- read.cdir(filename)
#' Fleiss1993_CR
#' 
#' # Same result as R Command example(Fleiss1993bin):
#' #
#' metacr(Fleiss1993_CR)
#' 
#' # Same result as R Command example(Fleiss1993cont):
#' #
#' metacr(Fleiss1993_CR, 1, 2)
#' 
#' @rdname read.cdir
#' @export read.cdir


read.cdir <- function(file, title = "Cochrane Review of Interventions",
                      exdir = tempdir(),
                      numbers.in.labels = TRUE,
                      ##
                      rob =
                        !missing(tool) | !missing(categories) |
                        !missing(col) | !missing(symbols),
                      tool = NULL,
                      categories = NULL,
                      col = NULL,
                      symbols = NULL,
                      ##
                      keep.orig = FALSE,
                      ...) {
  
  ##
  ##
  ## (1) Check argument and extract filename and Cochrane data package
  ##     identifier (CD)
  ##
  ##
  
  if (missing(file))
    stop("File name must be provided", call. = FALSE)
  ##
  chkchar(file, length = 1)
  ##
  if (!file.exists(file))
    stop("File '", file, "' does not exists.",
         call. = FALSE)
  ##
  path <- file
  file <- basename(file)
  if (substring(file, 1, 2) == "CD")
    CD <- unlist(strsplit(file, "-"))[1]
  else
    CD <- NULL
  ##
  data.package <- grepl("dataPackage.zip$", file)
  ##
  chkchar(title, length = 1)
  chklogical(numbers.in.labels)
  chkchar(exdir, length = 1)
  chklogical(keep.orig)
  chklogical(rob)
  ##
  ## Get rid of warnings 'no visible binding for global variable'
  ##
  common <- comp.no <- complab <- eligibility <-
    event.c <- event.e <- footnotes <- group.no <- label.c <- label.e <-
      label.left <- label.right <- level <- level.ma <- logscale <-
        lower.TE <- mean.c <- mean.e <- method <- model <- n.c <- n.e <-
          O.E <- outclab <- outclab2 <- outcome.no <- overall <- random <-
            sd.c <- sd.e <- seTE <- show.subgroup <- sm <- studlab <-
              Study <- swap.events <- TE <- test.subgroup <- type <-
                upper.TE <- V <- weight <- year <- NULL
  
  
  ##
  ##
  ## (2) Unzip Cochrane data package with subdirectories
  ##
  ##
  
  dir.create(exdir, showWarnings = FALSE, recursive = TRUE)
  unzip(path, exdir = exdir, ...)
  ##
  if (is.null(CD)) {
    lf <- list.files(exdir)
    ##
    if (any(grepl("-files$", lf)))
      CD <- unlist(strsplit(lf[grep("-files$", lf)], "-"))[1]
    else {
      ld <- list.dirs(exdir)
      ld <- ld[ld != exdir]
      exdir <- ld[1]
      lf <- list.files(exdir)
      CD <- unlist(strsplit(lf[grep("-files$", lf)], "-"))[1]
    }
  }
  ##
  if (data.package)
    wd <- exdir
  else
    wd <- paste0(exdir, "/", CD, "-files")  
  ##
  if (!dir.exists(wd))
    stop("Data directory '", wd, "' does not exist.")
  ##
  if (!data.package) {
    datazip <- paste0(wd, "/", CD, "-data.zip")
    ##
    if (!file.exists(datazip))
      stop("File '", datazip, "' does not exists.",
           call. = FALSE)
    else
      unzip(datazip, exdir = wd)
  }
  
  
  ##
  ##
  ## (3) Read files from directory CD-analysis-data
  ##
  ##
  
  dir.analysis <- paste0(wd, "/", CD, "-analysis-data")
  ##
  if (!dir.exists(dir.analysis))
    stop("Directory '", dir.analysis, "' does not exist.",
         call. = FALSE)
  ##
  ## Check for file specific to Cochrane review of DTA studies
  ##
  if (file.exists(paste0(dir.analysis, "/", CD, "-parameters.csv")))
    stop("Cochrane data package is from diagnostic test accuracy review.",
         call. = FALSE)
  ##
  file.settings <-
    paste0(dir.analysis, "/", CD, "-overall-estimates-and-settings.csv")
  file.datarows <-
    paste0(dir.analysis, "/", CD, "-data-rows.csv")
  file.subgroup <-
    paste0(dir.analysis, "/", CD, "-subgroup-estimates.csv")
  ##
  if (file.exists(file.settings))
    settings <- read_csv(file.settings, col_types = cols())
  else
    stop("File '", file.settings, "' does not exists.",
         call. = FALSE)
  ##
  if (file.exists(file.datarows))
    datarows <- read_csv(file.datarows, col_types = cols())
  else
    stop("File '", file.datarows, "' does not exists.",
         call. = FALSE)
  ##
  if (file.exists(file.subgroup))
    subgroup <- read_csv(file.subgroup, col_types = cols())
  else
    stop("File '", file.subgroup, "' does not exists.",
         call. = FALSE)
  ##
  if (keep.orig) {
    settings.orig <- settings
    datarows.orig <- datarows
    subgroup.orig <- subgroup
  }
  ##
  ## Check whether meta-analyses have been conducted
  ##
  if (nrow(datarows) == 0) {
    warning("No meta-analyses have been conducted.",
            call. = FALSE)
    return(NULL)
  }
  
  
  ##
  ##
  ## (4) Create study data set (with information on analyses)
  ##
  ##
  
  settings %<>%
    rename(
      comp.no = "Analysis group",
      outcome.no = "Analysis number",
      type = "Data type",
      method = "Statistical method",
      sm = "Effect measure",
      level = "Study CI",
      level.ma = "Estimate CI",
      unit = "Unit of effect measure",
      model = "Analysis model",
      logscale = "Log-scale data",
      ##
      label.e = "Experimental group label",
      label.c = "Control group label",
      ##
      overall = "Overall estimates",
      test.subgroup = "Test for subgroup differences",
      show.subgroup = "Subgroup estimates",
      ##
      swap.events = "Swap event and non-event",
      ##
      complab = "Analysis group name",
      outclab = "Analysis name",
      ##
      data.source = "Data source",
      eligibility = "Data source eligibility"
    ) %>%
  mutate(
    comp.no = as.numeric(comp.no),
    outcome.no = as.numeric(outcome.no),
    type = if_else(type == "Contrast level", "I",
                   if_else(type == "O-E and variance", "P",
                           substring(type, 1, 1))),
    common = model == "Fixed effect",
    random = model == "Random effect",
    method = sm2meth(method),
    sm = em2sm(sm, type),
    eligibility = if_else(eligibility == "#N/A", NA, eligibility),
    model = rmSpace(substring(model, 1, 6), end = TRUE),
    level = as.numeric(gsub("%", "", level)) / 100,
    level.ma = as.numeric(gsub("%", "", level.ma)) / 100,
    label.e = if_else(is.na(label.e), "Experimental", label.e),
    label.c = if_else(is.na(label.c), "Control", label.c),
    label.left = "",
    label.right = ""
  ) %>%
  select(
    comp.no, outcome.no, complab, outclab,
    type, method, sm, model, common, random, swap.events, logscale, unit,
    level, level.ma,
    label.e, label.c, label.left, label.right,
    overall, test.subgroup, show.subgroup, swap.events) %>%
  as.data.frame()
  ##
  datarows %<>%
    rename(
      comp.no = "Analysis group",
      outcome.no = "Analysis number",
      ##
      outclab2 = "Analysis name",
      ##
      subgroup = "Subgroup",
      studlab = "Study",
      year = "Study year",
      ##
      TE = "GIV Mean",
      seTE = "GIV SE",
      lower.TE = "CI start",
      upper.TE = "CI end",
      ##
      event.e = "Experimental cases",
      n.e = "Experimental N",
      mean.e = "Experimental mean",
      sd.e = "Experimental SD",
      ##
      event.c = "Control cases",
      n.c = "Control N",
      mean.c = "Control mean",
      sd.c = "Control SD",
      ##
      O.E = "O-E",
      V = "Variance",
      weight = "Weight",
      ##
      footnotes = "Footnotes"
    ) %>%
  mutate(
    comp.no = as.numeric(comp.no),
    outcome.no = as.numeric(outcome.no),
    lower.TE = if_else(is.na(TE), NA, lower.TE),
    upper.TE = if_else(is.na(TE), NA, upper.TE)
    ) %>%
    select(
      comp.no, outcome.no, outclab2, subgroup,
      studlab, year,
      event.e, n.e, event.c, n.c, mean.e, sd.e, mean.c, sd.c,
      O.E, V, TE, seTE, lower.TE, upper.TE,
      weight, footnotes) %>%
  as.data.frame()
  ##
  subgroup %<>%
    rename(
      comp.no = "Analysis group",
      outcome.no = "Analysis number",
      group.no = "Subgroup number",
      subgroup = "Subgroup"
    ) %>%
    mutate(
      comp.no = as.numeric(comp.no),
      outcome.no = as.numeric(outcome.no)
    ) %>%
    select(
      comp.no, outcome.no, group.no, subgroup
    ) %>%
    as.data.frame()
  ##
  ## Full Cochrane data set
  ##
  data <- merge(settings, datarows,
                by = c("comp.no", "outcome.no"), all = TRUE)
  ##
  data <- merge(data, subgroup,
                by = c("comp.no", "outcome.no", "subgroup"), all.x = TRUE)
  ##
  suppressWarnings(
    data %<>%
    mutate(
      lower.TE = if_else(logscale & !is.na(lower.TE), log(lower.TE), lower.TE),
      upper.TE = if_else(logscale & !is.na(upper.TE), log(upper.TE), upper.TE),
      order = seq_len(nrow(data))
    )
  )
  ##
  data %<>%
    select(
      comp.no, outcome.no, group.no, studlab, year,
      event.e, n.e, event.c, n.c, mean.e, sd.e, mean.c, sd.c,
      O.E, V, TE, seTE, lower.TE, upper.TE,
      weight, order,
      type, method, sm, model, common, random, swap.events, logscale,
      level, level.ma,
      label.e, label.c, label.left, label.right, unit,
      complab, outclab, footnotes
    )
  ##
  if (numbers.in.labels) {
    data$complab <- paste(data$comp.no, data$complab)
    data$outclab <- paste0(data$comp.no, ".", data$outcome.no, " ", data$outclab)
  }
  ##
  attr(data, "title") <- title
  
  
  ##
  ##
  ## (5) Extract risk of bias (if available)
  ##
  ##
  
  dir.study <- paste0(wd, "/", CD, "-study-data")
  ##
  robdata <- NULL
  ##
  if (dir.exists(dir.study)) {
    ##
    file.robdata <- paste0(dir.study, "/", CD, "-risk-of-bias.csv")
    ##
    rob <- rob && file.exists(file.robdata)
    ##
    if (rob) {
      robdata <- as.data.frame(read_csv(file.robdata, col_types = cols()))
      rob <- rob && nrow(robdata) >= 1
    }
    ##
    if (rob) {
      if (keep.orig)
        rob.orig <- robdata
      ##
      if (nrow(robdata) > 0) {
        robdata %<>%
          rename(studlab = Study)
        ##
        vars <- names(robdata)
        j <- 0
        k <- 0
        domains <- vector("character", 0)
        ##
        for (i in seq_along(vars)) {
          if (substring(vars[i], 1, 18) == "Domain (judgement)") {
            j <- j + 1
            robdata[paste0("D", j)] <- robdata[vars[i]]
            domains <- c(domains, substring(vars[i], 21))
          }
          ##
          if (substring(vars[i], 1, 16) == "Domain (support)") {
            k <- k + 1
            robdata[paste0("D", k, ".details")] <- robdata[vars[i]]
          }
        }
        ##
        selvars <- "studlab"
        ##
        if (j > 0)
          selvars <- c(selvars, paste0("D", seq_len(j)))
        if (k > 0)
          selvars <- c(selvars, paste0("D", seq_len(j), ".details"))
        ##
        robdata %<>% select(all_of(selvars))
        ##
        domains <- gsub(": All outcomes", "", domains)
        ##
        attr(robdata, "rob") <- rob
        attr(robdata, "domains") <- domains
        attr(robdata, "tool") <- tool
        attr(robdata, "categories") <- categories
        attr(robdata, "col") <- col
        attr(robdata, "symbols") <- symbols
      }
    }
  }
  else
    rob <- FALSE
  
  
  ##
  ##
  ## (6) Extract risk of bias (if available)
  ##
  ##
  
  res <- list(data = data, rob = robdata)
  ##
  if (keep.orig)
    res$orig <-
      list(settings = settings.orig, datarows = datarows.orig,
           subgroup = subgroup.orig, if (rob) rob = rob.orig)
  ##
  class(res) <- "cdir"
  res
}





#' @rdname read.cdir
#' @method print cdir
#' @export


print.cdir <- function(x, ...) {
  
  chkclass(x, "cdir")
  ##
  tl <- options()$width - 12
  newline <- FALSE
  ##
  title <- attr(x$data, "title")
  ##
  if (!is.null(title)) {
    if (title != "") {
      newline <- TRUE
      if (nchar(title) <= tl)
        cat(paste0("Review:     ", title, "\n"))
      else
        cat(paste0("Review:     ", substring(title, 1, tl - 4), " ...\n"))
    }
  }
  ##
  if (!is.null(x$data$complab)) {
    cat("Available comparisons:\n")
    ##
    cat(paste0(unique(x$data$complab), "\n", collapse = ""))
  }
  ##
  invisible(NULL)
}
