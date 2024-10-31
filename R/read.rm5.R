#' Import RevMan 5 analysis data
#' 
#' @description
#' Reads analysis data from Cochrane intervention review created with
#' RevMan 5 and creates a data frame from it.
#' 
#' @aliases read.rm5 Fleiss1993_CR
#' 
#' @param file The name of a file to read data values from.
#' @param sep The field separator character (only considered for
#'   CSV-files). Values on each line of the file are separated by this
#'   character. The comma is the default field separator character in
#'   RevMan 5.
#' @param quote The set of quoting characters (only considered for
#'   CSV-files). In RevMan 5 a "\"" is the default quoting character.
#' @param title Title of Cochrane review.
#' @param numbers.in.labels A logical indicating whether comparison
#'   number and outcome number should be printed at the beginning of
#'   the comparison (argument \code{complab}) and outcome label
#'   (argument \code{outclab}); this is the default in RevMan 5.
#' @param debug An integer between 0 and 3 indicating whether to print
#'   debug messages (only considered for RM5-files).
#' 
#' @details
#' Review Manager 5 (RevMan 5) was the software used for preparing and
#' maintaining Cochrane reviews. RevMan 5 includes the ability to write
#' systematic reviews of interventions, diagnostic test accuracy reviews,
#' methodology reviews and overviews of reviews.
#' 
#' This function provides the ability to read the analysis data from a
#' Cochrane intervention review created with RevMan 5; a data frame is
#' created from it. Cochrane intervention reviews are based on
#' comparisons of two interventions.
#' 
#' By default in RevMan 5, the name of the exported CSV data file is
#' the title of the Cochrane review. Furthermore, the title is part of
#' the RM5-file. Argument \code{title} can be used to overwrite the
#' title of the Cochrane review.
#'
#' \subsection{Import RM5-file}{
#' 
#' A RM5-file (which is in a specific XML format) can be used directly
#' to import the analysis dataset.
#'
#' If the import fails, use argument \code{debug = 3} for more details.
#' }
#'
#' \subsection{Import CSV-file}{
#' 
#' In the past, the following (rather complicated) procedure based on
#' a CSV-file generated within RevMan 5 was necessary - which is only
#' described here for backward compatibility.
#' 
#' In order to generate a data analysis file in RevMan 5 use the
#' following Menu points: \code{"File"} - \code{"Export"} -
#' \code{"Data and analyses"}. It is mandatory to include the
#' following fields in the exported data file by selecting them with
#' the mouse cursor in the Export Analysis Data Wizard: (i) Comparison
#' Number, (ii) Outcome Number, (iii) Subgroup Number. When these
#' fields are not selected a corresponding error message will be
#' printed in R.  It is recommended to include all fields in the
#' exported data file except for the last field
#' "Risk of bias tables". For example, in order to redo the
#' meta-analysis in R for the RevMan 5 data type
#' \code{"O-E and Variance"} the fields \code{"O-E"} and
#' \code{"Variance"} have to be selected in the Export Analysis Data
#' Wizard. If the last field "Risk of bias tables" is selected the
#' import in R fails with an error message
#' "line X did not have Y elements".
#' }
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
#' \item{outclab}{Outcome label.}
#' \item{k}{Total number of studies combined in respective
#'   meta-analysis).}
#' \item{event.e.pooled}{Number of events in experimental group in
#'   respective meta-analysis (see details).}
#' \item{n.e.pooled}{Number of observations in experimental group in
#'   respective meta-analysis (see details).}
#' \item{event.c.pooled}{Number of events in control group in
#'   respective meta-analysis (see details).}
#' \item{n.c.pooled}{Number of observations in control group in
#'   respective meta-analysis (see details).}
#' \item{TE.pooled}{Estimated treatment effect in respective
#'   meta-analysis (see details).}
#' \item{lower, upper}{Lower and upper limit of 95\% confidence
#'   interval for treatment effect in respective meta-analysis (see
#'   details).}
#' \item{weight.pooled}{Total weight in respective meta-analysis (see
#'   details).}
#' \item{Z.pooled}{Z-score for test of overall treatment effect in
#'   respective meta-analysis (see details).}
#' \item{pval.pooled}{P-value for test of overall treatment effect in
#'   respective meta-analysis (see details).}
#' \item{Q}{Heterogeneity statistic Q in respective meta-analysis (see
#'   details).}
#' \item{pval.Q}{P-value of heterogeneity statistic Q in respective
#'   meta-analysis (see details).}
#' \item{I2}{Heterogeneity statistic I\eqn{^2} in respective meta-analysis
#'   (see details).}
#' \item{tau2}{Between-study variance (moment estimator of
#'   DerSimonian-Laird) in respective meta-analysis.}
#' \item{Q.w}{Heterogeneity statistic Q within groups in respective
#'   meta-analysis (see details).}
#' \item{pval.Q.w}{P-value of heterogeneity statistic Q within groups
#'   in respective meta-analysis (see details).}
#' \item{I2.w}{Heterogeneity statistic I\eqn{^2} within groups in respective
#'   meta-analysis (see details).}
#' \item{label.e}{Label for experimental group.}
#' \item{label.c}{Label for control group.}
#' \item{label.left}{Graph label on left side of forest plot.}
#' \item{label.right}{Graph label on right side of forest plot.}
#' \item{complab}{Comparison label.}
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{summary.rm5}}, \code{\link{metabias.rm5}},
#'   \code{\link{metabin}}, \code{\link{metacont}},
#'   \code{\link{metagen}}, \code{\link{metacr}},
#'   \code{\link{print.rm5}}
#' 
#' @references
#' \emph{Review Manager (RevMan)} [Computer program]. Version 5.4.
#' The Cochrane Collaboration, 2020
#' 
#' @keywords datagen
#' 
#' @examples
#' # Locate export data file "Fleiss1993_CR.csv"
#' # in sub-directory of package "meta"
#' #
#' filename <- system.file("extdata", "Fleiss1993_CR.csv", package = "meta")
#' Fleiss1993_CR <- read.rm5(filename)
#' 
#' # Same result as R command example(Fleiss1993bin):
#' #
#' metacr(Fleiss1993_CR)
#' 
#' # Same result as R command example(Fleiss1993cont):
#' #
#' metacr(Fleiss1993_CR, 1, 2)
#'
#' \dontrun{
#' # Locate file "Fleiss1993.rm5" in sub-directory of R package meta
#' #
#' filename <- system.file("extdata/Fleiss1993.rm5", package = "meta")
#' Fleiss1993_CR <- read.cdir(filename)
#' Fleiss1993_CR
#' 
#' # Same result as R Command example(Fleiss1993bin):
#' #
#' metacr(Fleiss1993_CR)
#' }
#'
#' @importFrom xml2 as_xml_document xml_attr xml_find_all xml_text
#' 
#' @rdname read.rm5
#' @export read.rm5


read.rm5 <- function(file, sep = ",", quote = "\"",
                     title, numbers.in.labels = TRUE,
                     debug = 0) {
  ##
  if (missing(file))
    stop("File name must be provided", call. = FALSE)
  ##
  chkchar(file, length = 1)
  ##
  chkchar(sep, length = 1)
  chkchar(quote, length = 1)
  ##
  missing.title <- missing(title)
  if (!missing.title)
    chkchar(title, length = 1)
  ##
  chklogical(numbers.in.labels)
  ##
  if (!is.logical(debug))
    chknumeric(debug)
  if (!(debug %in% 0:3))
    stop("Argument 'debug' must be 0, 1, 2, or 3.", call. = FALSE)
  
  
  nc <- nchar(file)
  ##
  if (tolower(substring(file, nc - 2, nc)) == "csv") {
    if (missing.title)
      res <- read.rm5.csv(file, sep = sep, quote = quote,
                          numbers.in.labels = numbers.in.labels)
    else
      res <- read.rm5.csv(file, sep = sep, quote = quote,
                          title = title, numbers.in.labels = numbers.in.labels)
  }
  else if (tolower(substring(file, nc - 2, nc)) == "rm5") {
    if (missing.title)
      res <- read.rm5.rm5(file,
                          numbers.in.labels = numbers.in.labels,
                          debug = debug)
    else
      res <- read.rm5.rm5(file, title = title,
                          numbers.in.labels = numbers.in.labels,
                          debug = debug)
  }
  else
    stop("Argument 'file' must refer to a RM5- or CSV-file.", call. = FALSE)
  ##
  if (!is.null(res))
    attr(res, "filename") <- file
  ##
  res
}
