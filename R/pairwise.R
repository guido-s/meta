#' Transform meta-analysis data from two arm-based formats into contrast-based
#' format
#' 
#' @description
#' This function transforms data that are given in wide or long arm-based
#' format (e.g. input format for WinBUGS) to a contrast-based format that is
#' needed as input to R functions \code{\link{metabin}}, \code{\link{metacont}},
#' \code{\link{metainc}}, \code{\link{metagen}}, or
#' \code{\link[netmeta]{netmeta}} from R package \bold{netmeta}. The function
#' can transform data with binary, continuous, or generic outcomes as well as
#' incidence rates from arm-based to contrast-based format.
#' 
#' @param treat A list or vector with treatment information for
#'   individual treatment arms (see Details).
#' @param event A list or vector with information on number of events
#'   for individual treatment arms (see Details).
#' @param n A list or vector with information on number of
#'   observations for individual treatment arms (see Details).
#' @param mean A list or vector with estimated means for individual
#'   treatment arms (see Details).
#' @param sd A list or vector with information on the standard
#'   deviation for individual treatment arms (see Details).
#' @param TE A list or vector with estimated treatment effects for
#'   individual treatment arms (see Details).
#' @param seTE A list or vector with standard errors of estimated
#'   treatment effect for individual treatment arms (see Details).
#' @param time A list or vector with information on person time at
#'   risk for individual treatment arms (see Details).
#' @param agent A list or vector with agent information for
#'   individual treatment arms (see Details).
#' @param dose A list or vector with dose information for
#'   individual treatment arms (see Details).
#' @param data An optional data frame containing the study
#'   information.
#' @param studlab A vector with study labels (optional).
#' @param method A character string indicating which method is to be
#'   used to calculate treatment estimates (see Details).
#' @param sm A character string indicating which summary measure is to be
#'   used to calculate treatment estimates (see Details).
#' @param incr A numerical value which is added to cell frequencies
#'   for studies with a zero cell count, see Details.
#' @param method.incr A character string indicating which continuity
#'   correction method should be used (\code{"only0"},
#'   \code{"if0all"}, or \code{"all"}), see \code{\link{metabin}}.
#' @param allstudies A logical indicating if studies with zero or all
#'   events in two treatment arms are to be included in the
#'   meta-analysis (applies only if \code{sm} is equal to \code{"RR"}
#'   or \code{"OR"}).
#' @param reference.group Reference treatment (first treatment is used
#'   if argument is missing).
#' @param sep.ag A character used as separator between agent and dose to
#'   create treatment labels.
#' @param keep.all.comparisons A logical indicating whether all
#'   pairwise comparisons or only comparisons with the study-specific
#'   reference group should be kept ('basic parameters').
#' @param varnames Character vector of length 2 with the variable names for the
#'   treatment estimate and its standard error; by default, "TE" and "seTE".
#' @param append Either a logical indicating whether variables from the dataset
#'   provided in argument \code{data} are appended to the dataset with
#'   pairwise comparisons or a character vector with variable names to append to
#'   the dataset.
#' @param allincr Deprecated argument (replaced by 'method.incr');
#'   see \code{\link{metabin}}.
#' @param addincr Deprecated argument (replaced by 'method.incr');
#'   see \code{\link{metabin}}.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if studies are excluded due to only providing a single
#'   treatment arm).
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments passed-through to the functions
#'   to calculate effects.
#' 
#' @details
#' The pairwise function transforms data given in (wide or long)
#' arm-based format into the contrast-based format which consists of
#' \emph{pairwise} comparisons and which is needed as input to R functions
#' \code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metainc}},
#' \code{\link{metagen}}, or \code{\link[netmeta]{netmeta}} from R package
#' \bold{netmeta}.
#' 
#' The pairwise function can transform data with binary outcomes
#' continuous outcomes (\code{\link{metacont}} function), incidence
#' rates (\code{\link{metainc}} function), and generic outcomes
#' (\code{\link{metagen}} function). Depending on the outcome, the
#' following arguments are mandatory:
#' \itemize{
#' \item treat, event, n (see \code{\link{metabin}});
#' \item treat, n, mean, sd (see \code{\link{metacont}});
#' \item treat, event, time (see \code{\link{metainc}});
#' \item treat, TE, seTE (see \code{\link{metagen}}).
#' }
#' 
#' (using the \code{\link{metabin}} function from R package meta),
#'
#' Admissible values for arguments \code{method} and \code{sm} are outcome
#' specific; see help pages of R functions \code{\link{metabin}},
#' \code{\link{metacont}}, \code{\link{metainc}}, and \code{\link{metagen}}.
#' 
#' Argument \code{treat} is mandatory to identify the individual
#' treatments. The other arguments contain outcome specific
#' data. These arguments must be either lists (wide arm-based format,
#' i.e., one row per study) or vectors (long arm-based format,
#' i.e., multiple rows per study) of the same length.
#' 
#' For the wide arm-based format, each list consists of as many
#' vectors of the same length as the multi-arm study with the largest
#' number of treatments. If a single multi-arm study has five arms,
#' five vectors have to be provided for each lists. Two-arm studies
#' have entries with \code{NA} for the third and subsequent
#' vectors. Each list entry is a vector with information for each
#' individual study; i.e., the length of this vector corresponds to
#' the total number of studies incorporated in the network
#' meta-analysis. Typically, list elements are part of a data frame
#' (argument \code{data}, optional); see Examples. An optional vector
#' with study labels can be provided which can be part of the data
#' frame.
#' 
#' In the long arm-based format, argument \code{studlab} is mandatory
#' to identify rows contributing to individual studies.
#' 
#' Additional arguments for meta-analysis functions can be provided
#' using argument '\dots'; see help pages of R functions
#' \code{\link{metabin}}, \code{\link{metacont}},
#' \code{\link{metainc}}, and \code{\link{metagen}}.
#'
#' For standardised mean differences (argument \code{sm = "SMD"}),
#' equations (4) and (5) in Crippa & Orsini (2016) are used to
#' calculated SMDs and standard errors. These equations guarantee
#' consistent SMDs and standard errors for multi-arm studies. Note,
#' the summary measure is actually Cohen's d as Hedges' g is not
#' consistent in multi-arm studies.
#' 
#' For binary outcomes, 0.5 is added to all cell frequencies (odds
#' ratio) or only the number of events (risk ratio) for studies with a
#' zero cell count. For odds ratio and risk ratio, treatment estimates
#' and standard errors are only calculated for studies with zero or
#' all events in both groups if \code{allstudies} is \code{TRUE}. This
#' continuity correction is used both to calculate individual study
#' results with confidence limits and to conduct meta-analysis based
#' on the inverse variance method. For the risk difference, 0.5 is
#' only added to all cell frequencies to calculate the standard error.
#' 
#' For incidence rates, 0.5 is added to all cell frequencies for the
#' incidence rate ratio as summary measure. For the incidence risk
#' difference, 0.5 is only added to all cell frequencies to calculate
#' the standard error.
#' 
#' The value of pairwise is a data frame with as many rows as there
#' are pairwise comparisons. For each study with \emph{p} treatments,
#' \emph{p*(p-1) / 2} contrasts are generated. Each row contains the
#' treatment effect (\code{TE}), its standard error (\code{seTE}), the
#' treatments compared ((\code{treat1}), (\code{treat2})) and the
#' study label ((\code{studlab})). Further columns are added according
#' to type of data.
#' 
#' All variables from the original dataset are also part of the output
#' dataset if argument \code{append = TRUE}. If data are provided in the long
#' arm-based format, the value of a variable can differ between treatment arms;
#' for example, the mean age or percentage of women in the treatment arm. In
#' this situation, two variables instead of one variable will be included
#' in the output dataset. The values "1" and "2" are added to the
#' names for these variables, e.g. "mean.age1" and "mean.age2" for the
#' mean age.
#' 
#' In general, any variable names in the original dataset that are
#' identical to the main variable names (i.e., "TE", "seTE", ...) will
#' be renamed to variable names with ending ".orig".
#'
#' A reduced dataset with basic comparisons (Rücker & Schwarzer,
#' 2014) can be generated using argument \code{keep.all.comparisons =
#' FALSE}. Furthermore, the reference group for the basic comparisons
#' can be specified with argument \code{reference.group}.
#'
#' \subsection{Use in network meta-analysis}{
#' R function \code{\link[netmeta]{netmeta}} expects data in a
#' \bold{contrast-based format}, where each row corresponds to a
#' comparison of two treatments and contains a measure of the
#' treatment effect comparing two treatments with standard error,
#' labels for the two treatments and an optional study label.  In
#' contrast-based format, a three-arm study contributes three rows
#' with treatment comparison and corresponding standard error for
#' pairwise comparison \emph{A} vs \emph{B}, \emph{A} vs \emph{C}, and
#' \emph{B} vs \emph{C} whereas a four-arm study contributes six rows
#' / pairwise comparisons: \emph{A} vs \emph{B}, \emph{A} vs \emph{C},
#' \dots{}, \emph{C} vs \emph{D}.
#' 
#' Other programs for network meta-analysis in WinBUGS and Stata
#' require data in an \emph{arm-based} format, i.e. treatment estimate
#' for each treatment arm instead of a difference of two treatments. A
#' common \bold{(wide) arm-based format} consists of one data row per
#' study, containing treatment and other necessary information for all
#' study arms. For example, a four-arm study contributes one row with
#' four treatment estimates and corresponding standard errors for
#' treatments \emph{A}, \emph{B}, \emph{C}, and \emph{D}.  Another
#' possible arm-based format is a long format where each row
#' corresponds to a single study arm. Accordingly, in the \bold{long
#' arm-based format} a study contributes as many rows as treatments
#' considered in the study.
#' }
#' 
#' @note
#' This function must not be confused with \code{\link[netmeta]{netpairwise}}
#' which can be used to conduct pairwise meta-analyses for all
#' comparisons with direct evidence in a network meta-analysis.
#' 
#' @return
#' A data frame with the following columns:
#' \item{TE}{Treatment estimate comparing treatment 'treat1' and
#'   'treat2'.}
#' \item{seTE}{Standard error of treatment estimate.}
#' \item{studlab}{Study labels.}
#' \item{treat1}{First treatment in comparison.}
#' \item{treat2}{Second treatment in comparison.}
#' \item{event1}{Number of events for first treatment arm (for metabin
#'   and metainc).}
#' \item{event2}{Number of events for second treatment arm (for
#'   metabin and metainc).}
#' \item{n1}{Number of observations for first treatment arm (for
#'   metabin and metacont).}
#' \item{n2}{Number of observations for second treatment arm (for
#'   metabin and metacont).}
#' \item{mean1}{Estimated mean for first treatment arm (for
#'   metacont).}
#' \item{mean2}{Estimated mean for second treatment arm (for
#'   metacont).}
#' \item{sd1}{Standard deviation for first treatment arm (for
#'   metacont).}
#' \item{sd2}{Standard deviation for second treatment arm (for
#'   metacont).}
#' \item{TE1}{Estimated treatment effect for first treatment arm (for
#'   metagen).}
#' \item{TE2}{Estimated treatment effect for second treatment arm (for
#'   metagen).}
#' \item{seTE1}{Standard error of estimated treatment effect for first
#'   treatment arm (for metagen).}
#' \item{seTE2}{Standard error of estimated treatment effect for
#'   second treatment arm (for metagen).}
#' \item{time1}{Person time at risk for first treatment arm (for
#'   metainc).}
#' \item{time2}{Person time at risk for second treatment arm (for
#'   metainc).}
#' \item{agent1}{First agent in comparison.}
#' \item{agent2}{Second agent in comparison.}
#' \item{dose1}{Dose of first agent in comparison.}
#' \item{dose2}{Dose of second agent in comparison.}
#' 
#' All variables from the original dataset are also part of the output
#' dataset; see Details.
#' 
#' @author Gerta Rücker\email{gerta.ruecker@@uniklinik-freiburg.de}, Guido
#'   Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{longarm}}, \code{\link{metabin}},
#'   \code{\link{metacont}}, \code{\link{metagen}},
#'   \code{\link{metainc}}, \code{\link[netmeta]{netmeta}},
#'   \code{\link[netmeta]{netgraph.netmeta}},
#'   \code{\link[metadat]{dat.senn2013}},
#'   \code{\link[metadat]{dat.franchini2012}},
#'   \code{\link[metadat]{dat.franchini2012}}
#'
#' @references
#' Crippa A, Orsini N (2016):
#' Dose-response meta-analysis of differences in means.
#' \emph{BMC Medical Research Methodology},
#' \bold{16}:91.
#' 
#' @keywords datagen
#' 
#' @examples
#' pw0 <- pairwise(studlab = study, treat = treatment,
#'   n = ni, mean = mi, sd = sdi, data = dat.senn2013,
#'   append = c("study", "comment"))
#' head(pw0)
#' # Meta-analysis of studies comparing metformin to placebo
#' metagen(pw0, subset = treat1 == "metformin" & treat2 == "placebo")
#' 
#' \dontrun{
#' # Use pairwise() to run network meta-analyses
#' # (R package 'netmeta' must be available)
#' if (requireNamespace("netmeta", quietly = TRUE)) {
#'  # Example using continuous outcomes (internal call of function
#'  # metacont)
#'  #
#'  Franchini2012 <- dat.franchini2012
#'  # Transform data from arm-based format to contrast-based format
#'  pw1 <- pairwise(list(Treatment1, Treatment2, Treatment3),
#'    n = list(n1, n2, n3),
#'    mean = list(y1, y2, y3), sd = list(sd1, sd2, sd3),
#'    data = Franchini2012, studlab = Study)
#'  pw1
#' 
#'  # Conduct network meta-analysis
#'  library("netmeta")
#'  #
#'  net1 <- netmeta(pw1)
#'  net1
#' 
#'  # Draw network graphs
#'  #
#'  netgraph(net1, points = TRUE, cex.points = 3, cex = 1.5,
#'    thickness = "se.common")
#'  netgraph(net1, points = TRUE, cex.points = 3, cex = 1.5,
#'    plastic = TRUE, thickness = "se.common",
#'    iterate = TRUE)
#'  netgraph(net1, points = TRUE, cex.points = 3, cex = 1.5,
#'    plastic = TRUE, thickness = "se.common",
#'    iterate = TRUE, start = "eigen")
#' 
#'  # Example using generic outcomes (internal call of function
#'  # metagen)
#'  #
#'  # Calculate standard error for means y1, y2, y3
#'  Franchini2012$se1 <- with(Franchini2012, sqrt(sd1^2 / n1))
#'  Franchini2012$se2 <- with(Franchini2012, sqrt(sd2^2 / n2))
#'  Franchini2012$se3 <- with(Franchini2012, sqrt(sd3^2 / n3))
#'  # Transform data from arm-based format to contrast-based format
#'  # using means and standard errors (note, argument 'sm' has to be
#'  # used to specify that argument 'TE' is a mean difference)
#'  pw2 <- pairwise(list(Treatment1, Treatment2, Treatment3),
#'    TE = list(y1, y2, y3), seTE = list(se1, se2, se3),
#'    n = list(n1, n2, n3),
#'    data = Franchini2012, studlab = Study,
#'    sm = "MD")
#'  pw2
#' 
#'  # Compare pairwise objects pw1 (based on continuous outcomes) and pw2
#'  # (based on generic outcomes)
#'  #
#'  all.equal(
#'    pw1[, c("TE", "seTE", "studlab", "treat1", "treat2")],
#'    pw2[, c("TE", "seTE", "studlab", "treat1", "treat2")])
#' 
#'  # Same result as network meta-analysis based on continuous outcomes
#'  # (object net1)
#'  net2 <- netmeta(pw2)
#'  net2
#' 
#'  # Example with binary data
#'  #
#'  data(smokingcessation)
#'  # Transform data from arm-based format to contrast-based format
#'  # (internal call of metabin function). Argument 'sm' has to be used
#'  # for odds ratio as risk ratio (sm = "RR") is default of metabin
#'  # function.
#'  #
#'  pw3 <- pairwise(list(treat1, treat2, treat3),
#'    list(event1, event2, event3), list(n1, n2, n3),
#'    data = smokingcessation,
#'    sm = "OR")
#'  pw3
#' 
#'  # Conduct network meta-analysis
#'  #
#'  net3 <- netmeta(pw3)
#'  net3
#' 
#'  # Example with incidence rates
#'  #
#'  data(dietaryfat)
#' 
#'  # Transform data from arm-based format to contrast-based format
#'  #
#'  pw4 <- pairwise(list(treat1, treat2, treat3),
#'    list(d1, d2, d3), time = list(years1, years2, years3),
#'    studlab = ID,
#'    data = dietaryfat)
#'  pw4
#' 
#'  # Conduct network meta-analysis using incidence rate ratios (sm =
#'  # "IRR"). Note, the argument 'sm' is not necessary as this is the
#'  # default in R function metainc called internally.
#'  #
#'  net4 <- netmeta(pw4, sm = "IRR")
#'  summary(net4)
#' 
#'  # Example with long data format
#'  #
#'  # Transform data from long arm-based format to contrast-based
#'  # format Argument 'sm' has to be used for odds ratio as summary
#'  # measure; by default the risk ratio is used in the metabin
#'  # function called internally.
#'  #
#'  pw5 <- pairwise(treatment, event = r, n = N,
#'    studlab = author, data = dat.woods2010, sm = "OR")
#'  pw5
#' 
#'  # Conduct network meta-analysis
#'  net5 <- netmeta(pw5)
#'  net5
#' }
#' }
#' 
#' @export pairwise

pairwise <- function(treat,
                     event, n, mean, sd, TE, seTE, time,
                     agent, dose,
                     data = NULL, studlab,
                     #
                     method = "Inverse",
                     sm = NULL,
                     incr = gs("incr"),
                     method.incr = gs("method.incr"),
                     allstudies = gs("allstudies"),
                     #
                     reference.group,
                     keep.all.comparisons,
                     #
                     sep.ag = "*",
                     #
                     varnames = c("TE", "seTE"),
                     #
                     append = !is.null(data),
                     #
                     addincr = gs("addincr"),
                     allincr = gs("allincr"),
                     #
                     warn = FALSE, warn.deprecated = gs("warn.deprecated"),
                     ...) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  
  missing.event <- missing(event)
  missing.n <- missing(n)
  missing.mean <- missing(mean)
  missing.sd <- missing(sd)
  missing.TE <- missing(TE)
  missing.seTE <- missing(seTE)
  missing.time <- missing(time)
  missing.agent <- missing(agent)
  missing.dose <- missing(dose)
  missing.studlab <- missing(studlab)
  #
  missing.incr <- missing(incr)
  missing.method.incr <- missing(method.incr)
  missing.allstudies <- missing(allstudies)
  missing.allincr <- missing(allincr)
  missing.addincr <- missing(addincr)
  #
  missing.append <- missing(append)
  missing.sep.ag <- missing(sep.ag)
  missing.reference.group <- missing(reference.group)
  missing.keep.all.comparisons <- missing(keep.all.comparisons)
  missing.varnames <- missing(varnames)
  #
  chknumeric(incr, min = 0, length = 1)
  chklogical(allstudies)
  #
  chklogical(warn.deprecated)
  allincr <-
    deprecated2(method.incr, missing.method.incr, allincr, missing.allincr,
                warn.deprecated)
  addincr <-
    deprecated2(method.incr, missing.method.incr, addincr, missing.addincr,
                warn.deprecated)
  if (missing.method.incr) {
    if (is.logical(addincr) && addincr)
      method.incr <- "all"
    else if (is.logical(allincr) && allincr)
      method.incr <- "if0all"
  }
  #
  method.incr <-
    setchar(method.incr, gs("meth4incr")[gs("meth4incr") != "user"])
  #
  if (!is.character(append)) {
    chklogical(append, text = "or vector with variable names")
    append.logical <- append
  }
  else {
    if (is.null(data))
      append.logical <- NULL
    else {
      append <- setchar(append, names(data),
                        pre = "a logical or a character vector with ")
      #
      append.logical <- !is.null(append)
      #
      append <- c(append, paste0(append, 1), paste0(append, 2))
    }
  }
  #
  chklogical(warn)
  #
  chkchar(sep.ag)
  #
  chkchar(varnames, length = 2)
  #
  args <- list(...)
  nam.args <- names(args)
  
  
  ##
  ## Auxiliary functions
  ##
  sumzero <- function(x)
    sum(x[!is.na(x)] == 0)
  ##
  anytrue <- function(x)
    any(x == TRUE, na.rm = TRUE)
  
  
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
  #
  ##
  ## Catch studlab, treat, agent, dose, event, n, mean, sd, time from data:
  ##
  treat <- catch("treat", mc, data, sfsp)
  ##
  if (inherits(treat, "pairwise")) {
    is.pairwise <- TRUE
    #
    res.attr <- attributes(treat)
    res <- treat
    #
    if (any(replaceNULL(res.attr$varnames, c("TE", "seTE")) !=
            c("TE", "seTE"))) {
      names(res)[names(res) == res.attr$varnames[1]] <- "TE"
      names(res)[names(res) == res.attr$varnames[2]] <- "seTE"
      #
      if (missing.varnames)
        varnames <- res.attr$varnames
    }
    #
    if (missing.append) {
      append <- res.attr$append
      append.logical <- res.attr$append.logical
    }
    #
    if (missing.reference.group)
      if (!is.null(attributes(treat)$reference.group)) {
        reference.group <- attributes(treat)$reference.group
        missing.reference.group <- FALSE
      }
    #
    if (missing.keep.all.comparisons) {
      if (!is.null(attributes(treat)$keep.all.comparisons))
        keep.all.comparisons <- attributes(treat)$keep.all.comparisons
      else
        keep.all.comparisons <- TRUE
    }
    #
    txt.ignore <- "as first argument is a pairwise object"
    #
    warn_ignore_input(event, !missing.event, txt.ignore)
    warn_ignore_input(n, !missing.n, txt.ignore)
    warn_ignore_input(mean, !missing.mean, txt.ignore)
    warn_ignore_input(sd, !missing.sd, txt.ignore)
    warn_ignore_input(TE, !missing.TE, txt.ignore)
    warn_ignore_input(seTE, !missing.seTE, txt.ignore)
    warn_ignore_input(time, !missing.time, txt.ignore)
    warn_ignore_input(agent, !missing.agent, txt.ignore)
    warn_ignore_input(dose, !missing.dose, txt.ignore)
    warn_ignore_input(data, !nulldata, txt.ignore)
    warn_ignore_input(studlab, !missing.studlab, txt.ignore)
    warn_ignore_input(incr, !missing.incr, txt.ignore)
    warn_ignore_input(method.incr, !missing.method.incr, txt.ignore)
    warn_ignore_input(allincr, !missing.allincr, txt.ignore)
    warn_ignore_input(addincr, !missing.addincr, txt.ignore)
    warn_ignore_input(allstudies, !missing.allstudies, txt.ignore)
    #
    type <- attributes(res)$type
    #
    avail.treat <- TRUE
    avail.agent <- avail.dose <- FALSE
  }
  else {
    is.pairwise <- FALSE
    #
    if (missing.keep.all.comparisons)
      keep.all.comparisons <- TRUE
    chklogical(keep.all.comparisons)
    #
    agent <- catch("agent", mc, data, sfsp)
    dose <- catch("dose", mc, data, sfsp)
    #
    studlab <- catch("studlab", mc, data, sfsp)
    #
    event <- catch("event", mc, data, sfsp)
    n <- catch("n", mc, data, sfsp)
    mean <- catch("mean", mc, data, sfsp)
    sd <- catch("sd", mc, data, sfsp)
    TE <- catch("TE", mc, data, sfsp)
    seTE <- catch("seTE", mc, data, sfsp)
    time <- catch("time", mc, data, sfsp)
    #
    avail.studlab <- !is.null(studlab)
    #
    avail.event <- !is.null(event)
    avail.n <- !is.null(n)
    avail.mean <- !is.null(mean)
    avail.sd <- !is.null(sd)
    avail.TE <- !is.null(TE)
    avail.seTE <- !is.null(seTE)
    avail.time <- !is.null(time)
    #
    avail.treat <- !is.null(treat)
    avail.agent <- !is.null(agent)
    avail.dose <- !is.null(dose)
    
    #if (!(avail.treat || (avail.agent & avail.dose)))
    #  stop("Mandatory argument(s) are 'treat' or 'agent' and 'dose'.")
    #
    if (avail.treat & avail.agent & avail.dose)
      stop("Either provide argument 'treat' or arguments 'agent' and 'dose'.")
    #
    #if (!avail.treat & (avail.agent + avail.dose != 2))
    #  stop("Mandatory arguments 'agent' and 'dose' for dose-response data.")
    #
    if (avail.treat)
      if (is.list(treat))
        chklist(treat)
    #
    if (avail.agent)
      if (is.list(agent))
        chklist(agent)
    #
    if (avail.dose)
      if (is.list(dose))
        chklist(dose)
    ##
    if (avail.event)
      if (is.list(event))
        chklist(event)
      else
        chknumeric(event)
    ##
    if (avail.n)
      if (is.list(n))
        chklist(n)
      else
        chknumeric(n)
    ##
    if (avail.mean)
      if (is.list(mean))
        chklist(mean)
      else
        chknumeric(mean)
    ##
    if (avail.sd)
      if (is.list(sd))
        chklist(sd)
      else
        chknumeric(sd)
    ##
    if (avail.TE)
      if (is.list(TE))
        chklist(TE)
      else
        chknumeric(TE)
    ##
    if (avail.seTE)
      if (is.list(seTE))
        chklist(seTE)
      else
        chknumeric(seTE)
    ##
    if (avail.time)
      if (is.list(time))
        chklist(time)
      else
        chknumeric(time)
    
    
    if (avail.TE & avail.seTE)
      type <- "generic"
    else if (avail.event & avail.time & !avail.mean & !avail.sd)
      type <- "count"
    else if (avail.event & avail.n & !avail.mean & !avail.sd)
      type <- "binary"
    else if (avail.n & avail.mean & avail.sd)
      type <- "continuous"
    else if (!avail.event & !avail.mean & !avail.sd &
             !avail.TE & !avail.seTE & !avail.time)
      type <- "onlytreat"
    else
      stop("Type of outcome unclear. Please provide the necessary ",
           "information:\n  - event, n (binary outcome)\n  - n, ",
           "mean, sd (continuous outcome)\n  - TE, seTE (generic outcome)\n",
           "  - event, time (incidence rates).")
    
    
    #
    # Determine whether data is in long, wide or comparison-based format
    #
    if (type == "generic") {
      if (is.list(TE) & is.list(seTE)) {
        if (length(TE) == 1 & length(seTE) == 1 &
            (length(treat) == 1 || (length(agent) == 1 & length(dose) == 1))) {
          data.format <- "long"
          #
          TE <- unlist(TE)
          seTE <- unlist(seTE)
          #
          if (avail.treat)
            treat <- unlist(treat)
          else {
            agent <- unlist(agent)
            dose <- unlist(dose)
            #
            sep.ag <- setsep(agent, sep.ag, type = "agent")
            #
            treat <- paste(agent, dose, sep = sep.ag)
          }
        }
        else if (length(TE) == 2 & length(seTE) == 2 &
                 (length(treat) == 2 ||
                  (length(agent) == 2 & length(dose) == 2) ||
                  (!avail.treat & !avail.agent & !avail.dose))) {
          data.format <- "comparison"
          #
          if (!avail.treat & !avail.agent & !avail.dose) {
            treat <- list(expand("trt1", TE[[1]]), expand("trt2", TE[[2]]))
            avail.treat <- TRUE
          }
          #
          dat.st <- data.frame(studlab)
          #
          if (avail.treat) {
            dat.st$treat1 <- treat[[1]]
            dat.st$treat2 <- treat[[2]]
          }
          else {
            dat.st$agent1 <- agent[[1]]
            dat.st$dose1 <- dose[[1]]
            #
            dat.st$agent2 <- agent[[2]]
            dat.st$dose2 <- dose[[2]]
            #
            sep.ag <-
              setsep(c(dat.st$agent1, dat.st$agent2), sep.ag, type = "agent")
            #
            dat.st$treat1 <- paste(dat.st$agent1, dat.st$dose1, sep = sep.ag)
            dat.st$treat2 <- paste(dat.st$agent2, dat.st$dose2, sep = sep.ag)
          }
        }
        else {
          data.format <- "wide"
          #
          if (!avail.treat) {
            sep.ag <- setsep(unlist(agent), sep.ag, type = "agent")
            agent.without.NA <- lapply(agent, replaceNA, replace = "")
            dose.without.NA <- lapply(dose, replaceNA, replace = "")
            #
            treat.without.NA <-
              mapply(paste,
                     agent.without.NA, dose.without.NA,
                     MoreArgs = list(sep = sep.ag),
                     SIMPLIFY = FALSE)
            #
            treat <- lapply(treat.without.NA, replaceVal,
                            old = sep.ag, new = NA)
          }
        }
      }
      else {
        data.format <- "long"
        #
        if (!avail.treat) {
          sep.ag <- setsep(agent, sep.ag, type = "agent")
          treat <- paste(agent, dose, sep = sep.ag)
        }
      }
    }
    #
    else if (type == "binary") {
      if (is.list(event) & is.list(n)) {
        if (length(event) == 1 & length(n) == 1 &
            (length(treat) == 1 || (length(agent) == 1 & length(dose) == 1))) {
          data.format <- "long"
          #
          event <- unlist(event)
          n <- unlist(n)
          #
          if (avail.treat)
            treat <- unlist(treat)
          else {
            agent <- unlist(agent)
            dose <- unlist(dose)
            #
            sep.ag <- setsep(agent, sep.ag, type = "agent")
            #
            treat <- paste(agent, dose, sep = sep.ag)
          }
        }
        else if (length(event) == 2 & length(n) == 2 & 
                 (length(treat) == 2 ||
                  (length(agent) == 2 & length(dose) == 2) ||
                  (!avail.treat & !avail.agent & !avail.dose))) {
          data.format <- "comparison"
          #
          if (!avail.treat & !avail.agent & !avail.dose) {
            treat <- list(expand("trt1", n[[1]]), expand("trt2", n[[2]]))
            avail.treat <- TRUE
          }
          #
          dat.st <- data.frame(studlab)
          #
          if (avail.treat) {
            dat.st$treat1 <- treat[[1]]
            dat.st$treat2 <- treat[[2]]
          }
          else {
            dat.st$agent1 <- agent[[1]]
            dat.st$dose1 <- dose[[1]]
            #
            dat.st$agent2 <- agent[[2]]
            dat.st$dose2 <- dose[[2]]
            #
            sep.ag <-
              setsep(c(dat.st$agent1, dat.st$agent2), sep.ag, type = "agent")
            #
            dat.st$treat1 <- paste(dat.st$agent1, dat.st$dose1, sep = sep.ag)
            dat.st$treat2 <- paste(dat.st$agent2, dat.st$dose2, sep = sep.ag)
          }
        }
        else {
          data.format <- "wide"
          #
          if (!avail.treat) {
            sep.ag <- setsep(unlist(agent), sep.ag, type = "agent")
            agent.without.NA <- lapply(agent, replaceNA, replace = "")
            dose.without.NA <- lapply(dose, replaceNA, replace = "")
            #
            treat.without.NA <-
              mapply(paste,
                     agent.without.NA, dose.without.NA,
                     MoreArgs = list(sep = sep.ag),
                     SIMPLIFY = FALSE)
            #
            treat <- lapply(treat.without.NA, replaceVal,
                            old = sep.ag, new = NA)
          }
        }
      }
      else {
        data.format <- "long"
        #
        if (!avail.treat) {
          sep.ag <- setsep(agent, sep.ag, type = "agent")
          treat <- paste(agent, dose, sep = sep.ag)
        }
      }
    }
    #
    else if (type == "continuous") {
      if (is.list(n) & is.list(mean) & is.list(sd)) {
        if (length(n) == 1 & length(mean) == 1 & length(sd) == 1 &
            (length(treat) == 1 || (length(agent) == 1 & length(dose) == 1))) {
          data.format <- "long"
          #
          n <- unlist(n)
          mean <- unlist(mean)
          sd <- unlist(sd)
          #
          if (avail.treat)
            treat <- unlist(treat)
          else {
            agent <- unlist(agent)
            dose <- unlist(dose)
            #
            sep.ag <- setsep(agent, sep.ag, type = "agent")
            #
            treat <- paste(agent, dose, sep = sep.ag)
          }
        }
        else if (length(n) == 2 & length(mean) == 2 & length(sd) == 2 &
                 (length(treat) == 2 ||
                  (length(agent) == 2 & length(dose) == 2) ||
                  (!avail.treat & !avail.agent & !avail.dose))) {
          data.format <- "comparison"
          #
          if (!avail.treat & !avail.agent & !avail.dose) {
            treat <- list(expand("trt1", n[[1]]), expand("trt2", n[[2]]))
            avail.treat <- TRUE
          }
          #
          dat.st <- data.frame(studlab)
          #
          if (avail.treat) {
            dat.st$treat1 <- treat[[1]]
            dat.st$treat2 <- treat[[2]]
          }
          else {
            dat.st$agent1 <- agent[[1]]
            dat.st$dose1 <- dose[[1]]
            #
            dat.st$agent2 <- agent[[2]]
            dat.st$dose2 <- dose[[2]]
            #
            sep.ag <-
              setsep(c(dat.st$agent1, dat.st$agent2), sep.ag, type = "agent")
            #
            dat.st$treat1 <- paste(dat.st$agent1, dat.st$dose1, sep = sep.ag)
            dat.st$treat2 <- paste(dat.st$agent2, dat.st$dose2, sep = sep.ag)
          }
        }
        else {
          data.format <- "wide"
          #
          if (!avail.treat) {
            sep.ag <- setsep(unlist(agent), sep.ag, type = "agent")
            agent.without.NA <- lapply(agent, replaceNA, replace = "")
            dose.without.NA <- lapply(dose, replaceNA, replace = "")
            #
            treat.without.NA <-
              mapply(paste,
                     agent.without.NA, dose.without.NA,
                     MoreArgs = list(sep = sep.ag),
                     SIMPLIFY = FALSE)
            #
            treat <- lapply(treat.without.NA, replaceVal,
                            old = sep.ag, new = NA)
          }
        }
      }
      else {
        data.format <- "long"
        #
        if (!avail.treat) {
          sep.ag <- setsep(agent, sep.ag, type = "agent")
          treat <- paste(agent, dose, sep = sep.ag)
        }
      }
    }
    #
    else if (type == "count") {
      if (is.list(event) & is.list(time)) {
        if (length(event) == 1 & length(time) == 1 &
            (length(treat) == 1 || (length(agent) == 1 & length(dose) == 1))) {
          data.format <- "long"
          #
          event <- unlist(event)
          time <- unlist(time)
          #
          if (avail.n && is.list(n) && length(n) == 1)
            n <- unlist(n)
          #
          if (avail.treat)
            treat <- unlist(treat)
          else {
            agent <- unlist(agent)
            dose <- unlist(dose)
            #
            sep.ag <- setsep(agent, sep.ag, type = "agent")
            #
            treat <- paste(agent, dose, sep = sep.ag)
          }
        }
        else if (length(event) == 2 & length(time) == 2 &
                 (length(treat) == 2 ||
                  (length(agent) == 2 & length(dose) == 2) ||
                  (!avail.treat & !avail.agent & !avail.dose))) {
          data.format <- "comparison"
          #
          if (!avail.treat & !avail.agent & !avail.dose) {
            treat <- list(expand("trt1", time[[1]]), expand("trt2", time[[2]]))
            avail.treat <- TRUE
          }
          #
          dat.st <- data.frame(studlab)
          #
          if (avail.treat) {
            dat.st$treat1 <- treat[[1]]
            dat.st$treat2 <- treat[[2]]
          }
          else {
            dat.st$agent1 <- agent[[1]]
            dat.st$dose1 <- dose[[1]]
            #
            dat.st$agent2 <- agent[[2]]
            dat.st$dose2 <- dose[[2]]
            #
            sep.ag <-
              setsep(c(dat.st$agent1, dat.st$agent2), sep.ag, type = "agent")
            #
            dat.st$treat1 <- paste(dat.st$agent1, dat.st$dose1, sep = sep.ag)
            dat.st$treat2 <- paste(dat.st$agent2, dat.st$dose2, sep = sep.ag)
          }
        }
        else {
          data.format <- "wide"
          #
          if (!avail.treat) {
            sep.ag <- setsep(unlist(agent), sep.ag, type = "agent")
            agent.without.NA <- lapply(agent, replaceNA, replace = "")
            dose.without.NA <- lapply(dose, replaceNA, replace = "")
            #
            treat.without.NA <-
              mapply(paste,
                     agent.without.NA, dose.without.NA,
                     MoreArgs = list(sep = sep.ag),
                     SIMPLIFY = FALSE)
            #
            treat <- lapply(treat.without.NA, replaceVal,
                            old = sep.ag, new = NA)
          }
        }
      }
      else {
        data.format <- "long"
        #
        if (!avail.treat) {
          sep.ag <- setsep(agent, sep.ag, type = "agent")
          treat <- paste(agent, dose, sep = sep.ag)
        }
      }
    }
    #
    else if (type == "onlytreat") {
      if (is.list(treat)) {
        if ((length(treat) == 1 || (length(agent) == 1 & length(dose) == 1))) {
          data.format <- "long"
          #
          if (avail.treat)
            treat <- unlist(treat)
          else {
            agent <- unlist(agent)
            dose <- unlist(dose)
            #
            sep.ag <- setsep(agent, sep.ag, type = "agent")
            #
            treat <- paste(agent, dose, sep = sep.ag)
          }
        }
        else if (length(treat) == 2 ||
                 (length(agent) == 2 & length(dose) == 2)) {
          data.format <- "comparison"
          #
          dat.st <- data.frame(studlab)
          #
          if (avail.treat) {
            dat.st$treat1 <- treat[[1]]
            dat.st$treat2 <- treat[[2]]
          }
          else {
            dat.st$agent1 <- agent[[1]]
            dat.st$dose1 <- dose[[1]]
            #
            dat.st$agent2 <- agent[[2]]
            dat.st$dose2 <- dose[[2]]
            #
            sep.ag <-
              setsep(c(dat.st$agent1, dat.st$agent2), sep.ag, type = "agent")
            #
            dat.st$treat1 <- paste(dat.st$agent1, dat.st$dose1, sep = sep.ag)
            dat.st$treat2 <- paste(dat.st$agent2, dat.st$dose2, sep = sep.ag)
          }
        }
        else {
          data.format <- "wide"
          #
          if (!avail.treat) {
            sep.ag <- setsep(unlist(agent), sep.ag, type = "agent")
            agent.without.NA <- lapply(agent, replaceNA, replace = "")
            dose.without.NA <- lapply(dose, replaceNA, replace = "")
            #
            treat.without.NA <-
              mapply(paste,
                     agent.without.NA, dose.without.NA,
                     MoreArgs = list(sep = sep.ag),
                     SIMPLIFY = FALSE)
            #
            treat <- lapply(treat.without.NA, replaceVal,
                            old = sep.ag, new = NA)
          }
        }
      }
      else {
        data.format <- "long"
        #
        if (!avail.treat) {
          sep.ag <- setsep(agent, sep.ag, type = "agent")
          treat <- paste(agent, dose, sep = sep.ag)
        }
      }
    }
    
     
    #
    # Use longarm() to transform outcome variables from comparison-based
    # to long arm-based format
    #
    if (data.format == "comparison") {
      if (!avail.studlab)
        stop("Argument 'studlab' mandatory for comparison-based format.")
      #
      if (type == "binary") {
        ldat <- longarm(studlab = unlist(studlab),
                        treat1 = expand(treat[[1]], n[[1]]),
                        treat2 = expand(treat[[2]], n[[2]]),
                        event1 = event[[1]], event2 = event[[2]],
                        n1 = n[[1]], n2 = n[[2]])
        #
        studlab <- ldat$studlab
        treat <- ldat$treat
        event <- ldat$event
        n <- ldat$n
      }
      #
      else if (type == "continuous") {
        ldat <- longarm(studlab = unlist(studlab),
                        treat1 = expand(treat[[1]], n[[1]]),
                        treat2 = expand(treat[[2]], n[[2]]),
                        n1 = n[[1]], n2 = n[[2]],
                        mean1 = mean[[1]], mean2 = mean[[2]],
                        sd1 = sd[[1]], sd2 = sd[[2]])
        #
        studlab <- ldat$studlab
        treat <- ldat$treat
        n <- ldat$n
        mean <- ldat$mean
        sd <- ldat$sd
      }
      #
      else if (type == "count") {
        avail.n.comp <- avail.n && is.list(n) && length(n) == 2
        #
        ldat <- longarm(studlab = unlist(studlab),
                        treat1 = expand(treat[[1]], time[[1]]),
                        treat2 = expand(treat[[2]], time[[2]]),
                        event1 = event[[1]], event2 = event[[2]],
                        time1 = time[[1]], time2 = time[[2]],
                        n1 = if (avail.n.comp) n[[1]] else NULL,
                        n2 = if (avail.n.comp) n[[2]] else NULL)
        #
        studlab <- ldat$studlab
        treat <- ldat$treat
        event <- ldat$event
        time <- ldat$time
        #
        if (avail.n.comp)
          n <- ldat$n
      }
    }
    
    
    #
    # Transform outcome variables from long arm-based to list format
    #
    if (data.format %in% c("comparison", "long")) {
      if (!avail.studlab)
        stop("Argument 'studlab' mandatory for long arm-based format.")
      ##
      studlab <- as.character(studlab)
      ##
      treat <- as.character(treat)
      ##
      ttab <- table(studlab, treat)
      n.arms <- apply(ttab, 1, sum)
      max.arms <- max(n.arms)
      ##
      treat.list <- vector("list", max.arms)
      event.list <- vector("list", max.arms)
      n.list     <- vector("list", max.arms)
      mean.list  <- vector("list", max.arms)
      sd.list    <- vector("list", max.arms)
      TE.list    <- vector("list", max.arms)
      seTE.list  <- vector("list", max.arms)
      time.list  <- vector("list", max.arms)
      ##
      if (!nulldata)
        adddata <- vector("list", max.arms)
      ##
      if (type == "binary") {
        ##
        ## Generate lists
        ##
        tdat <- data.frame(studlab, treat, event, n,
                           .order = seq_along(treat), stringsAsFactors = FALSE)
        ##
        if (!nulldata & data.format == "long") {
          tdat <- cbind(tdat, data)
          dupl <- duplicated(names(tdat))
          if (any(dupl))
            names(tdat)[dupl] <- paste(names(tdat)[dupl], "orig", sep = ".")
        }
        ##
        studlab <- names(n.arms)
        dat.studlab <- data.frame(studlab, stringsAsFactors = FALSE)
        ##
        for (i in 1:max.arms) {
          sel.i <- !duplicated(tdat$studlab)
          tdat.i <- merge(dat.studlab, tdat[sel.i, ],
                          by = "studlab", all.x = TRUE)
          ##
          treat.list[[i]] <- tdat.i$treat
          event.list[[i]] <- tdat.i$event
          n.list[[i]]     <- tdat.i$n
          ##
          tdat.i$event <- NULL
          tdat.i$n     <- NULL
          ##
          if (!nulldata)
            adddata[[i]] <- tdat.i
          ##
          tdat <- tdat[!sel.i, ]
        }
        ##
        treat <- treat.list
        event <- event.list
        n     <- n.list
      }
      ##
      else if (type == "continuous") {
        ##
        ## Generate lists
        ##
        tdat <- data.frame(studlab, treat, n, mean, sd,
                           .order = seq_along(treat), stringsAsFactors = FALSE)
        ##
        if (!nulldata & data.format == "long") {
          tdat <- cbind(tdat, data)
          dupl <- duplicated(names(tdat))
          if (any(dupl))
            names(tdat)[dupl] <- paste(names(tdat)[dupl], "orig", sep = ".")
        }
        ##
        studlab <- names(n.arms)
        dat.studlab <- data.frame(studlab, stringsAsFactors = FALSE)
        ##
        for (i in 1:max.arms) {
          sel.i <- !duplicated(tdat$studlab)
          tdat.i <- merge(dat.studlab, tdat[sel.i, ],
                          by = "studlab", all.x = TRUE)
          ##
          treat.list[[i]] <- tdat.i$treat
          n.list[[i]]     <- tdat.i$n
          mean.list[[i]]  <- tdat.i$mean
          sd.list[[i]]    <- tdat.i$sd
          ##
          tdat.i$n    <- NULL
          tdat.i$mean <- NULL
          tdat.i$sd   <- NULL
          ##
          if (!nulldata)
            adddata[[i]] <- tdat.i
          ##
          tdat <- tdat[!sel.i, ]
        }
        ##
        treat <- treat.list
        n     <- n.list
        mean  <- mean.list
        sd    <- sd.list
      }
      ##
      else if (type == "count") {
        ##
        ## Generate lists
        ##
        tdat <- data.frame(studlab, treat, event, time,
                           .order = seq_along(treat), stringsAsFactors = FALSE)
        ##
        if (avail.n)
          tdat$n <- n
        ##
        if (!nulldata & data.format == "long") {
          tdat <- cbind(tdat, data)
          dupl <- duplicated(names(tdat))
          if (any(dupl))
            names(tdat)[dupl] <- paste(names(tdat)[dupl], "orig", sep = ".")
        }
        ##
        studlab <- names(n.arms)
        dat.studlab <- data.frame(studlab, stringsAsFactors = FALSE)
        ##
        for (i in 1:max.arms) {
          sel.i <- !duplicated(tdat$studlab)
          tdat.i <- merge(dat.studlab, tdat[sel.i, ],
                          by = "studlab", all.x = TRUE)
          ##
          treat.list[[i]] <- tdat.i$treat
          event.list[[i]] <- tdat.i$event
          time.list[[i]]  <- tdat.i$time
          #
          if (avail.n)
            n.list[[i]] <- tdat.i$n
          ##
          tdat.i$event <- NULL
          tdat.i$time  <- NULL
          ##
          if (!nulldata)
            adddata[[i]] <- tdat.i
          ##
          tdat <- tdat[!sel.i, ]
        }
        #
        treat <- treat.list
        event <- event.list
        time  <- time.list
        #
        if (avail.n)
          n  <- n.list
      }
      ##
      else if (type == "generic") {
        ##
        ## Generate lists
        ##
        tdat <- data.frame(studlab, treat, TE, seTE,
                           .order = seq_along(treat), stringsAsFactors = FALSE)
        ##
        if (avail.n)
          tdat$n <- n
        ##
        if (avail.event)
          tdat$event <- event
        ##
        if (!nulldata & data.format == "long") {
          tdat <- cbind(tdat, data)
          dupl <- duplicated(names(tdat))
          if (any(dupl))
            names(tdat)[dupl] <- paste(names(tdat)[dupl], "orig", sep = ".")
        }
        ##
        studlab <- names(n.arms)
        dat.studlab <- data.frame(studlab, stringsAsFactors = FALSE)
        ##
        for (i in 1:max.arms) {
          sel.i <- !duplicated(tdat$studlab)
          tdat.i <- merge(dat.studlab, tdat[sel.i, ],
                          by = "studlab", all.x = TRUE)
          ##
          treat.list[[i]] <- tdat.i$treat
          TE.list[[i]]    <- tdat.i$TE
          seTE.list[[i]]  <- tdat.i$seTE
          ##
          tdat.i$TE   <- NULL
          tdat.i$seTE <- NULL
          ##
          if (!nulldata)
            adddata[[i]] <- tdat.i
          ##
          tdat <- tdat[!sel.i, ]
        }
        ##
        treat <- treat.list
        TE    <- TE.list
        seTE  <- seTE.list
      }
      ##
      else if (type == "onlytreat") {
        ##
        ## Generate lists
        ##
        tdat <- data.frame(studlab, treat,
                           .order = seq_along(treat), stringsAsFactors = FALSE)
        ##
        if (avail.n)
          tdat$n <- n
        ##
        if (!nulldata) {
          if (data.format == "long")
            tdat <- cbind(tdat, data)
          dupl <- duplicated(names(tdat))
          if (any(dupl))
            names(tdat)[dupl] <- paste(names(tdat)[dupl], "orig", sep = ".")
        }
        ##
        studlab <- names(n.arms)
        dat.studlab <- data.frame(studlab, stringsAsFactors = FALSE)
        ##
        for (i in 1:max.arms) {
          sel.i <- !duplicated(tdat$studlab)
          tdat.i <- merge(dat.studlab, tdat[sel.i, ],
                          by = "studlab", all.x = TRUE)
          ##
          treat.list[[i]] <- tdat.i$treat
          ##
          if (!nulldata)
            adddata[[i]] <- tdat.i
          ##
          tdat <- tdat[!sel.i, ]
        }
        ##
        treat <- treat.list
      }
    }

    
    
    
    ##
    ## Check and set study labels
    ##
    if (!avail.studlab)
      studlab <- seq(along = treat[[1]])
    ##
    if (length(treat) != 2 && length(studlab) != length(unique(studlab)))
      stop("Study labels must all be distinct.")
    ##
    levs <- unique(studlab)


    narms <- length(unique(treat))
    nstud <- length(unique(studlab))


    ##
    ##
    ## Generate dataset with variables from original dataset
    ##
    ##
    if (!nulldata & data.format == "long") {
      names.adddata <- names(adddata[[1]])
      ##
      notunique <- matrix(NA,
                          ncol = length(names.adddata),
                          nrow = narms * (narms - 1) / 2)
      colnames(notunique) <- names.adddata
      ##
      oneNA <- matrix(NA,
                      ncol = length(names.adddata),
                      nrow = narms * (narms - 1) / 2)
      ##
      allNA <- matrix(NA,
                      ncol = length(names.adddata),
                      nrow = narms * (narms - 1) / 2)
      colnames(oneNA) <- names.adddata
      ##
      n.ij <- 0
      ##
      for (i in 1:(narms - 1)) {
        for (j in (i + 1):narms) {
          n.ij <- n.ij + 1
          notunique[n.ij, ] <-
            apply(adddata[[i]] != adddata[[j]], 2, anytrue)
          ##
          allNA[n.ij, ] <- apply(is.na(adddata[[j]]), 2, all)
          ##
          oneNA[n.ij, ] <-
            apply(is.na(adddata[[i]]) & !is.na(adddata[[j]]), 2, anytrue)
        }
      }
      ##
      notunique <- apply(notunique, 2, anytrue)
      oneNA <- apply(oneNA, 2, anytrue)
      allNA <- apply(allNA, 2, anytrue)
      notunique <- apply(rbind(notunique, oneNA, allNA), 2, anytrue)
      ## print(notunique)
      ##
      for (i in 1:(narms - 1)) {
        for (j in (i + 1):narms) {
          dat.i <- adddata[[i]]
          dat.j <- adddata[[j]]
          ##
          if (any(!notunique))
            dat.ij <- dat.i[, names.adddata[!notunique], drop = FALSE]
          else
            stop("Study label must be unique for single treatment arm.")
          ##
          for (nam in names.adddata[notunique]) {
            dat.ij[, paste0(nam, 1)] <- adddata[[i]][nam]
            dat.ij[, paste0(nam, 2)] <- adddata[[j]][nam]
          }
          ##
          if (i == 1 & j == 2)
            newdata <- dat.ij
          else
            newdata <- rbind(newdata, dat.ij)
        }
      }
      ##
      names.basic <- c("studlab", "treat1", "treat2")
      names.newdata <- names(newdata)
      ##
      newdata <- newdata[, c(names.basic,
                             names.newdata[!(names.newdata %in% names.basic)])]
      newdata <- newdata[!is.na(newdata$treat1) & !is.na(newdata$treat2), ]
    }
    
    
    if (type == "binary") {
      #
      method <- setchar(method, gs("meth4bin"))
      #
      if (length(event) != narms)
        stop("Different length of lists 'treat' and 'event'.")
      if (length(n) != narms)
        stop("Different length of lists 'treat' and 'n'.")
      ##
      ## Determine increment for individual studies
      ##
      n.zeros <- apply(matrix(unlist(event), ncol = length(event)), 1, sumzero)
      n.all   <- apply(matrix(unlist(n), ncol = length(event)) -
                       matrix(unlist(event), ncol = length(event)),
                       1, sumzero)
      ##
      incr.study <- rep(0, length(n.zeros))
      ##
      if (is.null(sm))
        sm <- if (method == "Peto") "OR" else gs("smbin")
      else
        sm <- setchar(sm, gs("sm4bin"))
      #
      addincr <- allincr <- FALSE
      #
      if (!(sm == "ASD" | method == "Peto")) {
        if (method.incr == "all")
          addincr <- TRUE
        else if (method.incr == "if0all")
          allincr <- TRUE
      }
      #
      sparse <- switch(sm,
                       OR = (n.zeros > 0) | (n.all > 0),
                       RD = (n.zeros > 0) | (n.all > 0),
                       RR = (n.zeros > 0) | (n.all > 0),
                       ASD = rep(FALSE, length(n.zeros)))
      ##
      if (!allincr & !addincr)
        incr.study[sparse] <- incr
      else if (addincr)
        incr.study[] <- incr
      else {
        if (any(n.zeros > 0))
          incr.study[] <- incr
        else
          incr.study[] <- 0
      }
      ##
      for (i in 1:(narms - 1)) {
        ##
        if (i == 1 & (length(treat[[i]]) != length(event[[i]])))
          stop("Different length of element ", i, " of ",
               "lists 'treat' and 'event'.")
        if (i == 1 & (length(event[[i]]) != length(n[[i]])))
          stop("Different length of element ", i, " of ",
               "lists 'event' and 'n'.")
        ##
        for (j in (i + 1):narms) {
          ##
          if (length(treat[[j]]) != length(event[[j]]))
            stop("Different length of element ", j, " of ",
                 "lists 'treat' and 'event'.")
          if (length(event[[j]]) != length(n[[j]]))
            stop("Different length of element ", j, " of ",
                 "lists 'event' and 'n'.")
          ##
          dat <- data.frame(studlab, treat1 = treat[[i]], treat2 = treat[[j]],
                            #
                            TE = NA, seTE = NA,
                            #
                            event1 = event[[i]], n1 = n[[i]],
                            event2 = event[[j]], n2 = n[[j]],
                            #
                            incr1 = incr.study, incr2 = incr.study,
                            #
                            .order = seq_along(studlab),
                            stringsAsFactors = FALSE, row.names = NULL)
          ##
          if (!nulldata & data.format == "wide") {
            dat <- cbind(dat, data, stringsAsFactors = FALSE)
            dupl <- duplicated(names(dat))
            if (any(dupl))
              names(dat)[dupl] <- paste(names(dat)[dupl], "orig", sep = ".")
          }
          ##
          dat <- dat[!(is.na(dat$event1) & is.na(dat$n1)), ]
          dat <- dat[!(is.na(dat$event2) & is.na(dat$n2)), ]
          ##
          if (nrow(dat) > 0) {
            m1 <- metabin(dat$event1, dat$n1, dat$event2, dat$n2,
                          #
                          method = method, sm = sm,
                          #
                          incr.e = dat$incr1, incr.c = dat$incr2,
                          method.incr = "user", allstudies = allstudies,
                          #
                          method.tau = "DL", method.tau.ci = "",
                          #
                          warn = warn,
                          warn.deprecated = FALSE, ...)
            ##
            dat$TE <- m1$TE
            dat$seTE <- m1$seTE
            ##
            dat$TE[is.infinite(dat$TE)] <- NA
            dat$seTE[is.infinite(dat$seTE)] <- NA
            #
            dat$incr1 <- m1$incr.e
            dat$incr2 <- m1$incr.c
            #
            dat.NAs <- dat[is.na(dat$TE) | is.na(dat$seTE) | dat$seTE <= 0, ]
            ##
            if (i == 1 & j == 2) {
              res <- dat
              res.NAs <- dat.NAs
            }
            else {
              res <- rbind(res, dat)
              res.NAs <- rbind(res.NAs, dat.NAs)
            }
          }
          else
            if (i == 1 & j == 2)
              stop("No studies available for comparison of ",
                   "first and second treatment.")
        }
      }
    }
    #
    else if (type == "continuous") {
      if (length(n) != narms)
        stop("Different length of lists 'treat' and 'n'.")
      if (length(mean) != narms)
        stop("Different length of lists 'treat' and 'mean'.")
      if (length(sd) != narms)
        stop("Different length of lists 'treat' and 'sd'.")
      #
      if (is.null(sm))
        sm <- gs("smcont")
      else
        sm <- setchar(sm, gs("sm4cont"))
      #
      for (i in seq_len(narms)) {
        ##
        if (length(treat[[i]]) != length(n[[i]]))
          stop("Different length of element ", i, " of ",
               "lists 'treat' and 'n'.",
               call. = FALSE)
        if (length(treat[[i]]) != length(mean[[i]]))
          stop("Different length of element ", i, " of ",
               "lists 'treat' and 'mean'.",
               call. = FALSE)
        if (length(treat[[i]]) != length(sd[[i]]))
          stop("Different length of element ", i, " of ",
               "lists 'treat' and 'sd'.",
               call. = FALSE)
        if (length(treat[[i]]) != nstud)
          stop("Different length of study labels and ",
               "element ", i, " of list 'treat'.",
               call. = FALSE)
      }
      ##
      ## For standardised mean difference, calculate pooled standard
      ## deviation for multi-arm studies
      ##
      if (sm == "SMD" & narms > 2) {
        pooled.sd <- function(sd, n) {
          sel <- !is.na(sd) & !is.na(n)
          ##
          if (any(sel))
            res <- sqrt(sum((n[sel] - 1) * sd[sel]^2) / sum(n[sel] - 1))
          else
            res <- NA
          ##
          res
        }
        ##
        N <- matrix(unlist(n), ncol = narms, nrow = nstud, byrow = FALSE)
        M <- matrix(unlist(mean), ncol = narms, nrow = nstud, byrow = FALSE)
        S <- matrix(unlist(sd), ncol = narms, nrow = nstud, byrow = FALSE)
        ##
        sel.n <- apply(!is.na(N) & N > 0, 1, sum) > 2
        sel.mean <- apply(!is.na(M), 1, sum) > 2
        sel.sd <- apply(!is.na(S) & S > 0, 1, sum) > 2
        sel <- sel.n & sel.mean & sel.sd
        ##
        if (any(sel)) {
          N <- N[sel, , drop = FALSE]
          S <- S[sel, , drop = FALSE]
          sd.p <- rep_len(NA, nrow(N))
          ##
          for (i in seq_len(nrow(N)))
            sd.p[i] <- pooled.sd(S[i, ], N[i, ])
        }
        ##
        for (i in seq_len(narms))
          sd[[i]][sel] <- ifelse(is.na(sd[[i]][sel]), NA, sd.p)
      }
      ##
      for (i in 1:(narms - 1)) {
        for (j in (i + 1):narms) {
          dat <- data.frame(studlab, treat1 = treat[[i]], treat2 = treat[[j]],
                            #
                            TE = NA, seTE = NA,
                            n1 = n[[i]], mean1 = mean[[i]], sd1 = sd[[i]],
                            n2 = n[[j]], mean2 = mean[[j]], sd2 = sd[[j]],
                            #
                            .order = seq_along(studlab),
                            stringsAsFactors = FALSE, row.names = NULL)
          ##
          if (data.format == "wide") {
            dat <- cbind(dat, data, stringsAsFactors = FALSE)
            dupl <- duplicated(names(dat))
            if (any(dupl))
              names(dat)[dupl] <- paste(names(dat)[dupl], "orig", sep = ".")
          }
          ##
          dat <- dat[!(is.na(dat$n1) & is.na(dat$mean1) & is.na(dat$sd1)), ]
          dat <- dat[!(is.na(dat$n2) & is.na(dat$mean2) & is.na(dat$sd2)), ]
          ##
          if (nrow(dat) > 0) {
            m1 <- metacont(dat$n1, dat$mean1, dat$sd1,
                           dat$n2, dat$mean2, dat$sd2,
                           #
                           sm = sm,
                           method.tau = "DL", method.tau.ci = "",
                           method.smd = "Cohen",
                           #
                           warn = warn,
                           warn.deprecated = FALSE, ...)
            ##
            dat$TE <- m1$TE
            dat$seTE <- m1$seTE
            ##
            dat$TE[is.infinite(dat$TE)] <- NA
            dat$seTE[is.infinite(dat$seTE)] <- NA
            ##
            dat.NAs <- dat[is.na(dat$TE) | is.na(dat$seTE) | dat$seTE <= 0, ]
            ##
            if (i == 1 & j == 2) {
              res <- dat
              res.NAs <- dat.NAs
            }
            else {
              res <- rbind(res, dat)
              res.NAs <- rbind(res.NAs, dat.NAs)
            }
          }
          else
            if (i == 1 & j == 2)
              stop("No studies available for comparison of ",
                   "first and second treatment.",
                   call. = FALSE)
        }
      }
    }
    #
    else if (type == "generic") {
      if (length(TE) != narms)
        stop("Different length of lists 'treat' and 'TE'.",
             call. = FALSE)
      if (length(seTE) != narms)
        stop("Different length of lists 'treat' and 'seTE'.",
             call. = FALSE)
      #
      if (is.null(sm))
        sm <- ""
      #
      for (i in 1:(narms - 1)) {
        ##
        if (i == 1 & (length(treat[[i]]) != length(TE[[i]])))
          stop("Different length of element ", i, " of ",
               "lists 'treat' and 'TE'.",
               call. = FALSE)
        if (i == 1 & (length(treat[[i]]) != length(seTE[[i]])))
          stop("Different length of element ", i, " of ",
               "lists 'treat' and 'seTE'.",
               call. = FALSE)
        ##
        for (j in (i + 1):narms) {
          ##
          if (length(treat[[j]]) != length(TE[[j]]))
            stop("Different length of element ", j, " of ",
                 "lists 'treat' and 'TE'.",
                 call. = FALSE)
          if (length(treat[[j]]) != length(seTE[[j]]))
            stop("Different length of element ", j, " of ",
                 "lists 'treat' and 'seTE'.",
                 call. = FALSE)
          ##
          dat <- data.frame(studlab, treat1 = treat[[i]], treat2 = treat[[j]],
                            #
                            TE = NA, seTE = NA,
                            TE1 = TE[[i]], seTE1 = seTE[[i]],
                            TE2 = TE[[j]], seTE2 = seTE[[j]],
                            #
                            .order = seq_along(studlab),
                            stringsAsFactors = FALSE, row.names = NULL)
          ##
          if (avail.event) {
            dat$event1 <- event[[i]]
            dat$event2 <- event[[j]]
          }
          ##
          if (avail.n) {
            dat$n1 <- n[[i]]
            dat$n2 <- n[[j]]
          }
          ##
          if (avail.mean) {
            dat$mean1 <- mean[[i]]
            dat$mean2 <- mean[[j]]
          }
          ##
          if (avail.sd) {
            dat$sd1 <- sd[[i]]
            dat$sd2 <- sd[[j]]
          }
          ##
          if (avail.time) {
            dat$time1 <- time[[i]]
            dat$time2 <- time[[j]]
          }
          ##
          if (data.format == "wide") {
            dat <- cbind(dat, data, stringsAsFactors = FALSE)
            dupl <- duplicated(names(dat))
            if (any(dupl))
              names(dat)[dupl] <- paste(names(dat)[dupl], "orig", sep = ".")
          }
          ##
          dat <- dat[!(is.na(dat$TE1) & is.na(dat$seTE1)), ]
          dat <- dat[!(is.na(dat$TE2) & is.na(dat$seTE2)), ]
          ##
          if (nrow(dat) > 0) {
            m1 <- metagen(dat$TE1 - dat$TE2,
                          sqrt(dat$seTE1^2 + dat$seTE2^2),
                          #
                          sm = sm,
                          method.tau = "DL", method.tau.ci = "",
                          #
                          warn = warn,
                          warn.deprecated = FALSE, ...)
            ##
            dat$TE <- m1$TE
            dat$seTE <- m1$seTE
            ##
            dat$TE[is.infinite(dat$TE)] <- NA
            dat$seTE[is.infinite(dat$seTE)] <- NA
            ##
            dat.NAs <- dat[is.na(dat$TE) | is.na(dat$seTE) | dat$seTE <= 0, ]
            ##
            if (i == 1 & j == 2) {
              res <- dat
              res.NAs <- dat.NAs
            }
            else {
              res <- rbind(res, dat)
              res.NAs <- rbind(res.NAs, dat.NAs)
            }
          }
          else
            if (i == 1 & j == 2)
              stop("No studies available for comparison of ",
                   "first and second treatment.",
                   call. = FALSE)
        }
      }
    }
    #
    else if (type == "count") {
      #
      method <- setchar(method, gs("meth4inc"))
      #
      if (length(event) != narms)
        stop("Different length of lists 'treat' and 'event'.",
             call. = FALSE)
      if (length(time) != narms)
        stop("Different length of lists 'treat' and 'time'.",
             call. = FALSE)
      #
      if (is.null(sm))
        sm <- gs("sminc")
      else
        sm <- setchar(sm, gs("sm4inc"))
      #
      addincr <- allincr <- FALSE
      #
      if (method.incr == "all")
        addincr <- TRUE
      else if (method.incr == "if0all")
        allincr <- TRUE
      ##
      ## Determine increment for individual studies
      ##
      n.zeros <- apply(matrix(unlist(event), ncol = length(event)), 1, sumzero)
      ##
      incr.study <- rep(0, length(n.zeros))
      ##
      sparse <- n.zeros > 0
      ##
      if (!allincr & !addincr)
        incr.study[sparse] <- incr
      else if (addincr)
        incr.study[] <- incr
      else {
        if (any(n.zeros > 0))
          incr.study[] <- incr
        else
          incr.study[] <- 0
      }
      ##
      for (i in 1:(narms - 1)) {
        #
        if (i == 1 & (length(treat[[i]]) != length(event[[i]])))
          stop("Different length of element ", i, " of ",
               "lists 'treat' and 'event'.",
               call. = FALSE)
        #
        if (i == 1 & (length(treat[[i]]) != length(time[[i]])))
          stop("Different length of element ", i, " of ",
               "lists 'treat' and 'time'.",
               call. = FALSE)
        #
        if (avail.n)
          if (i == 1 & (length(treat[[i]]) != length(n[[i]])))
            stop("Different length of element ", i, " of ",
                 "lists 'treat' and 'n'.",
                 call. = FALSE)
        ##
        for (j in (i + 1):narms) {
          ##
          if (length(treat[[j]]) != length(event[[j]]))
            stop("Different length of element ", j, " of ",
                 "lists 'treat' and 'event'.",
                 call. = FALSE)
          #
          if (length(treat[[j]]) != length(time[[j]]))
            stop("Different length of element ", j, " of ",
                 "lists 'treat' and 'time'.",
                 call. = FALSE)
          #
          if (avail.n)
            if (length(treat[[j]]) != length(n[[j]]))
              stop("Different length of element ", j, " of ",
                   "lists 'treat' and 'n'.",
                   call. = FALSE)
          ##
          dat <- data.frame(studlab, treat1 = treat[[i]], treat2 = treat[[j]],
                            #
                            TE = NA, seTE = NA,
                            #
                            event1 = event[[i]], time1 = time[[i]],
                            event2 = event[[j]], time2 = time[[j]],
                            #
                            incr1 = incr.study, incr2 = incr.study,
                            #
                            n1 = if (avail.n) n[[i]] else NA,
                            n2 = if (avail.n) n[[j]] else NA,
                            #
                            .order = seq_along(studlab),
                            stringsAsFactors = FALSE, row.names = NULL)
          ##
          if (data.format == "wide") {
            dat <- cbind(dat, data, stringsAsFactors = FALSE)
            dupl <- duplicated(names(dat))
            if (any(dupl))
              names(dat)[dupl] <- paste(names(dat)[dupl], "orig", sep = ".")
          }
          ##
          dat <- dat[!(is.na(dat$event1) & is.na(dat$time1)), ]
          dat <- dat[!(is.na(dat$event2) & is.na(dat$time2)), ]
          ##
          if (nrow(dat) > 0) {
            m1 <- metainc(dat$event1, dat$time1,
                          dat$event2, dat$time2,
                          #
                          n.e = if (avail.n) dat$n1 else NULL,
                          n.c = if (avail.n) dat$n2 else NULL,
                          #
                          method = method, sm = sm,
                          #
                          incr.e = dat$incr1, incr.c = dat$incr2,
                          method.incr = "user",
                          #
                          method.tau = "DL", method.tau.ci = "",
                          #
                          warn = warn, warn.deprecated = FALSE, ...)
            ##
            dat$TE <- m1$TE
            dat$seTE <- m1$seTE
            ##
            dat$TE[is.infinite(dat$TE)] <- NA
            dat$seTE[is.infinite(dat$seTE)] <- NA
            #
            dat$incr1 <- m1$incr.e
            dat$incr2 <- m1$incr.c
            #
            dat.NAs <- dat[is.na(dat$TE) | is.na(dat$seTE) | dat$seTE <= 0, ]
            ##
            if (i == 1 & j == 2) {
              res <- dat
              res.NAs <- dat.NAs
            }
            else {
              res <- rbind(res, dat)
              res.NAs <- rbind(res.NAs, dat.NAs)
            }
          }
          else
            if (i == 1 & j == 2)
              stop("No studies available for comparison of ",
                   "first and second treatment.",
                   call. = FALSE)
        }
      }
    }
    #
    else if (type == "onlytreat") {
      if (is.null(sm))
        sm <- ""
      #
      for (i in 1:(narms - 1)) {
        for (j in (i + 1):narms) {
          ##
          dat <- data.frame(studlab, treat1 = treat[[i]], treat2 = treat[[j]],
                            .order = seq_along(studlab),
                            stringsAsFactors = FALSE, row.names = NULL)
          ##
          if (avail.n) {
            dat$n1 <- n[[i]]
            dat$n2 <- n[[j]]
          }
          ##
          if (data.format == "wide") {
            dat <- cbind(dat, data, stringsAsFactors = FALSE)
            dupl <- duplicated(names(dat))
            if (any(dupl))
              names(dat)[dupl] <- paste(names(dat)[dupl], "orig", sep = ".")
          }
          ##
          dat <- dat[!is.na(dat$treat2), ]
          ##
          if (nrow(dat) > 0) {
            dat.NAs <- data.frame()
            ##
            if (i == 1 & j == 2) {
              res <- dat
              res.NAs <- dat.NAs
            }
            else {
              res <- rbind(res, dat)
              res.NAs <- rbind(res.NAs, dat.NAs)
            }
          }
          else
            if (i == 1 & j == 2)
              stop("No studies available for comparison of ",
                   "first and second treatment.",
                   call. = FALSE)
        }
      }
    }
    #
    method <- if (type == "onlytreat") "" else m1$method
    #
    if (!nulldata & data.format == "long")
      res <- merge(res, newdata,
                   by = c("studlab", "treat1", "treat2"),
                   suffixes = c("",".orig"),
                   all.x = TRUE)
    #
    if (!nulldata & data.format == "comparison") {
      res$.order <- NULL
      res2 <- data.frame()
      for (i in seq_len(nrow(dat.st))) {
        sel.i <-
          res$studlab == dat.st$studlab[i] &
          res$treat1 == dat.st$treat1[i] &
          res$treat2 == dat.st$treat2[i]
        #
        sel.wo.i <-
          res$studlab == dat.st$studlab[i] &
          res$treat1 == dat.st$treat2[i] &
          res$treat2 == dat.st$treat1[i]
        #
        if (any(sel.i))
          res2 <- rbind(res2, res[sel.i, ])
        else if (any(sel.wo.i)) {
          res2.i <- res[sel.wo.i, ]
          #
          res2.i$TE <- -res2.i$TE
          #
          t1 <- res2.i$treat1
          res2.i$treat1 <- res2.i$treat2
          res2.i$treat2 <- t1
          #
          if (type == "binary") {
            ev1 <- res2.i$event1
            res2.i$event1 <- res2.i$event2
            res2.i$event2 <- ev1
            #
            tn1 <- res2.i$n1
            res2.i$n1 <- res2.i$n2
            res2.i$n2 <- tn1
          }
          #
          else if (type == "continuous") {
            tn1 <- res2.i$n1
            res2.i$n1 <- res2.i$n2
            res2.i$n2 <- tn1
            #
            me1 <- res2.i$mean1
            res2.i$mean1 <- res2.i$mean2
            res2.i$mean2 <- me1
            #
            ts1 <- res2.i$sd1
            res2.i$sd1 <- res2.i$sd2
            res2.i$sd2 <- ts1
          }
          #
          else if (type == "count") {
            ev1 <- res2.i$event1
            res2.i$event1 <- res2.i$event2
            res2.i$event2 <- ev1
            #
            tt1 <- res2.i$time1
            res2.i$time1 <- res2.i$time2
            res2.i$time2 <- tt1
            #
            if (avail.n) {
              tn1 <- res2.i$n1
              res2.i$n1 <- res2.i$n2
              res2.i$n2 <- tn1
            }
          }
          res2 <- rbind(res2, res2.i)
        }
      }
      #
      res <- res2
      #
      if (!nulldata) {
        res <- cbind(res, data)
        dupl <- duplicated(names(res))
        if (any(dupl))
          names(res)[dupl] <- paste(names(res)[dupl], "orig", sep = ".")
      }
    }
    
    
    ##
    ## Additional checks
    ##
    ##
    ## a) Duplicate treatments ?
    ##
    sel.treat <- as.character(res$treat1) == as.character(res$treat2)
    ##
    if (any(sel.treat)) {
      sel.stud <- unique(sort(res$studlab[sel.treat]))
      ##
      stop(paste0("Identical treatments for the following stud",
                  if (length(sel.stud) == 1) "y: " else "ies:\n  ",
                  paste0(paste0("'", sel.stud, "'"),
                         collapse = " - "),
                  "\n  Please check dataset."),
           call. = FALSE)
    }
    ##
    ## b) Studies missing ?
    ##
    sel.study <- !(studlab %in% unique(as.character(res$studlab)))
    ##
    if (any(sel.study) & warn)
      warning(paste0("The following stud",
                     if (sum(sel.study) == 1) "y is " else "ies are ",
                     "excluded from the analysis\n  ",
                     "(due to a single study arm or missing values):",
                     if (sum(sel.study) == 1) " " else "\n  ",
                     paste0(paste0("'", studlab[sel.study], "'"),
                            collapse = " - ")),
              call. = FALSE)
    ##
    ## c) Missing treatment estimates or standard errors?
    ##
    if (type != "onlytreat" && nrow(res.NAs) > 0 & warn) {
      warning("Comparison",
              if (nrow(res.NAs) > 1) "s",
              " with missing TE / seTE or zero seTE",
              " will not be considered in network meta-analysis.",
              call. = FALSE)
      cat("Comparison",
          if (nrow(res.NAs) > 1) "s",
          " will not be considered in network meta-analysis:\n",
          sep = "")
      ##
      res.NAs$.order <- NULL
      res.NAs$.order1 <- NULL
      res.NAs$.order2 <- NULL
      ##
      prmatrix(res.NAs,
               quote = FALSE, right = TRUE, na.print = "NA",
               rowlab = rep("", nrow(res.NAs)))
    }
    
    
    ## Calculate standard error for SMDs
    ##
    if (type == "continuous" && m1$sm == "SMD") {
      res$.seTE <- res$seTE
      ##
      for (i in unique(res$studlab)) {
        sel.i <- res$studlab == i
        dat.i <- res[sel.i, ]
        ##
        ## Calculate total sample size in study i
        ##
        ndat.i <- rbind(data.frame(treat = dat.i$treat1, n = dat.i$n1),
                        data.frame(treat = dat.i$treat2, n = dat.i$n2))
        n.i <- sum(ndat.i[!duplicated(ndat.i), ]$n)
        ##
        ## Crippa & Orsini (2016), BMC Med Res Meth, equation (5)
        ##
        varTE.i <-
          1 / res$n1[sel.i] + 1 / res$n2[sel.i] + res$TE[sel.i]^2 / (2 * n.i)
        ##
        res$seTE[sel.i] <- sqrt(varTE.i)
      }
    }
    
     
    if (data.format != "comparison") {
      if (!is.null(res$.order1)) {
        res <- res[order(res$.order1), ]
        res$.order1 <- NULL
        res$.order2 <- NULL
        res$.order <- NULL
        res <- unique(res)
      }
      else if (!is.null(res$.order)) {
        res <- res[order(res$.order), ]
        res$.order <- NULL
        res$.order.orig <- NULL
        res <- unique(res)
      }
      else {
        res <- res[order(factor(res$studlab, levels = levs),
                         res$treat1, res$treat2), ]
      }
    }
    #
    rownames(res) <- 1:nrow(res)
  }
  #
  # Only keep core variables
  #
  corevars <- c("TE", "seTE", "studlab", "treat1", "treat2",
                "event1", "event2", "n1", "n2",
                "mean1", "mean2", "sd1", "sd2",
                "TE1", "TE2", "seTE1", "seTE2",
                "time1", "time2", "incr1", "incr2",
                ".seTE")
  #
  if (!append.logical)
    res <- res[, names(res) %in% corevars]
  else if (is.character(append))
    res <- res[, names(res) %in% c(corevars, append)]
  #
  # Drop columns 'n1' and 'n2' if argument 'n' is missing 
  #
  if (type == "count" & missing.n) {
    res$n1 <- NULL
    res$n2 <- NULL
  }
  
  
  ##
  ## Use first treatment with estimable effect as reference if
  ## argument is missing
  ##
  labels <- unique(sort(c(res$treat1, res$treat2)))
  ##
  if (missing.reference.group) {
    if (type == "onlytreat")
      reference.group <- ""
    else {
      go.on <- TRUE
      i <- 0
      while (go.on) {
        i <- i + 1
        sel.i <-
          !is.na(res$TE) & !is.na(res$seTE) &
          (res$treat1 == labels[i] | res$treat2 == labels[i])
        if (sum(sel.i) > 0) {
          go.on <- FALSE
          reference.group <- labels[i]
          }
        else if (i == length(labels)) {
          go.on <- FALSE
          reference.group <- ""
        }
      }
    }
  }
  ##
  if (is.factor(reference.group))
    reference.group <- as.character(reference.group)
  ##
  if (is.numeric(reference.group))
    chknumeric(reference.group, length = 1)
  else
    chkchar(reference.group, length = 1)
  ##
  reference.group <- setchar(reference.group, c(labels, ""))
  ##
  if (!keep.all.comparisons) {
    ##
    drop <- logical(0)
    ##
    for (i in unique(res$studlab)) {
      d.i <- res[res$studlab == i, , drop = FALSE]
      trts.i <- unique(sort(c(d.i$treat1, d.i$treat2)))
      ##
      ## Keep comparisons with reference group or first treatment if
      ## reference treatment is missing in study
      ##
      if (reference.group %in% trts.i)
        ref.i <- reference.group
      else
        ref.i <- rev(trts.i)[1]
      ##
      drop.i <- d.i$treat1 != ref.i & d.i$treat2 != ref.i
      ##
      drop <- c(drop, drop.i)
    }
    ##
    if (!keep.all.comparisons)
      res <- res[!drop, , drop = FALSE]
  }
  
  if (!avail.treat) {
    if (isCol(res, "agent1"))
      names(res)[names(res) == "agent1"] <- "agent1.orig"
    if (isCol(res, "dose1"))
      names(res)[names(res) == "dose1"] <- "dose1.orig"
    if (isCol(res, "agent2"))
      names(res)[names(res) == "agent2"] <- "agent2.orig"
    if (isCol(res, "dose2"))
      names(res)[names(res) == "dose2"] <- "dose2.orig"
    #
    res$agent1 <- sapply(compsplit(res$treat1, sep.ag), first)
    res$dose1 <- sapply(compsplit(res$treat1, sep.ag), second)
    res$agent2 <- sapply(compsplit(res$treat2, sep.ag), first)
    res$dose2 <- sapply(compsplit(res$treat2, sep.ag), second)
    #
    res$dose1 <- as.numeric(res$dose1)
    res$dose2 <- as.numeric(res$dose2)
    #
    nam.res <- names(res)
    nam1 <- c("TE", "seTE", "studlab")
    nam2 <- c("agent1", "dose1", "agent2", "dose2")
    #
    res <- res[, c(nam1, nam2, nam.res[!nam.res %in% c(nam1, nam2)])]
  }
  
  if (any(varnames != c("TE", "seTE"))) {
    names(res)[names(res) == "TE"] <- varnames[1]
    names(res)[names(res) == "seTE"] <- varnames[2]
  }
  
  attr(res, "pairwise") <- TRUE
  attr(res, "reference.group") <- reference.group
  attr(res, "keep.all.comparisons") <- keep.all.comparisons
  attr(res, "type") <- type
  #
  if (is.pairwise) {
    attr(res, "sm") <- res.attr$sm
    attr(res, "method") <- res.attr$method
    #
    attr(res, "incr") <- res.attr$incr
    attr(res, "method.incr") <- res.attr$method.incr
    #
    attr(res, "allstudies") <- res.attr$allstudies
  }
  else {
    attr(res, "sm") <- if (type != "onlytreat") replaceNULL(sm, "") else ""
    attr(res, "method") <- if (type != "onlytreat") method else ""
    #
    if (type %in% c("binary", "count")) {
      attr(res, "incr") <- incr
      attr(res, "method.incr") <- method.incr
      #
      if (type == "binary")
        attr(res, "allstudies") <- allstudies
    }
  }
  #
  attr(res, "varnames") <- varnames
  attr(res, "append") <- append
  attr(res, "append.logical") <- append.logical
  attr(res, "version") <- packageDescription("meta")$Version
  #
  attr(res, "args") <- args
  
  class(res) <- unique(c("pairwise", class(res)))
  
  res
}
