\name{read.rm5}
\alias{read.rm5}
\alias{Fleiss93_CR}
\title{Import RevMan 5 data files (.csv)}

\description{
  Reads data file from Cochrane Intervention review created with
  RevMan 5 and creates a data frame from it.
}

\usage{
read.rm5(file, sep=",", quote = "\"", title,
         numbers.in.labels=TRUE)
}

\arguments{
  \item{file}{The name of a file to read data values from.}
  \item{sep}{The field separator character. Values on each line of the
    file are separated by this character. The comma is the default field
    separator character in RevMan 5.}
  \item{quote}{The set of quoting characters. In RevMan 5 a "\""
    is the default quoting character.}
  \item{title}{Title of Cochrane review.}
  \item{numbers.in.labels}{A logical indicating whether comparision
    number and outcome number should be printed at the beginning of the
    comparison (argument \code{complab}) and outcome label (argument
    \code{outclab}); this is the default in RevMan 5.}
}

\details{
  Review Manager 5 (RevMan 5) is the current software used for
  preparing and maintaining Cochrane Reviews
  (\url{http://tech.cochrane.org/revman/}). RevMan 5 includes the ability
  to write Systematic reviews of interventions, Diagnostic test accuracy
  reviews, Methodology reviews and Overviews of reviews.
  
  This function provides the ability to read a data file from a
  Cochrane Intervention review created with RevMan 5; a data frame
  is created from it. Cochrane Intervention reviews are based on the
  comparison of two interventions.

  In order to generate a data analysis file in RevMan 5 use the
  following Menu points: \code{"File"} - \code{"Export"} - \code{"Data
  and analyses"}. It is mandatory to include the following fields in
  the exported data file by selecting them with the mouse cursor in
  the Export Analysis Data Wizard: (i) Comparison Number, (ii) Outcome
  Number, (iii) Subgroup Number. When these fields are not selected a
  corresponding error message will be printed in R. It is recommended
  to include all fields in the exported data file except for the last
  field "Risk of bias tables". For example, in order to redo the
  meta-analysis in R for the RevMan 5 data type \code{"O-E and
  Variance"} the fields \code{"O-E"} and \code{"Variance"} have to be
  selected in the Export Analysis Data Wizard. If the last field "Risk
  of bias tables" is selected the import in R fails with an error
  message "line X did not have Y elements".

  By default in RevMan 5, the name of the exported data file is the
  title of the Cochrane Review. Accordingly, information on the title is
  extracted from the name of the exported data file (argument:
  \code{file}) if argument \code{title} is missing (default).
  
  Each respective meta-analysis for arguments \code{event.e.pooled} --
  \code{df.pooled} is defined by values for \code{"comp.no"} and
  \code{"outcome.no"}, and \code{"grp.no"}.  }

\value{
  A data frame containing the following components:
  \item{comp.no}{Comparison number.}
  \item{outcome.no}{Outcome number.}
  \item{group.no}{Group number.}
  \item{studlab}{Study label.}
  \item{year}{Year of publication.}
  \item{event.e}{Number of events in experimental group.}
  \item{n.e}{Number of observations in experimental group.}
  \item{event.c}{Number of events in control group.}
  \item{n.c}{Number of observations in control group.}
  \item{mean.e}{Estimated mean in experimental group.}
  \item{sd.e}{Standard deviation in experimental group.}
  \item{mean.c}{Estimated mean in control group.}
  \item{sd.c}{Standard deviation in control group.}
  \item{O.E}{Observed minus expected (IPD analysis).}
  \item{V}{Variance of \code{O.E} (IPD analysis).}
  \item{TE, seTE}{Estimated treatment effect and standard error of individual studies.}
  \item{lower, upper}{Lower and upper limit of 95\% confidence
    interval for treatment effect in individual studies.}
  \item{weight}{Weight of individual studies (according to
    meta-analytical method used in respective meta-analysis - see below
    for details).}
  \item{order}{Ordering of studies.}
  \item{grplab}{Group label.}
  \item{type}{Type of outcome. D = dichotomous, C = continuous, P = IPD.}
  \item{method}{A character string indicating which method has been used
    for pooling of studies. One of \code{"Inverse"}, \code{"MH"}, or
    \code{"Peto"}.}
  \item{sm}{A character string indicating which summary measure
    has been used for pooling of studies.}
  \item{model}{A character string indicating which meta-analytical model
    has been used (either \code{"Fixed"} or \code{"Random"}).}
  \item{comb.fixed}{A logical indicating whether fixed effect
    meta-analysis has been used in respective meta-analysis (see below
    for details).}
  \item{comb.random}{A logical indicating whether random effects
    meta-analysis has been used in respective meta-analysis (see below
    for details).}
  \item{outclab}{Outcome label.}
  \item{k}{Total number of studies combined in respective
    meta-analysis).}
  \item{event.e.pooled}{Number of events in experimental group in
    respective meta-analysis (see below for details).}
  \item{n.e.pooled}{Number of observations in experimental group in
    respective meta-analysis (see below for details).}
  \item{event.c.pooled}{Number of events in control group in
    respective meta-analysis (see below for details).}
  \item{n.c.pooled}{Number of observations in control group in
    respective meta-analysis (see below for details).}
  \item{TE.pooled}{Estimated treatment effect in respective
    meta-analysis (see below for details).}
  \item{lower, upper}{Lower and upper limit of 95\% confidence
    interval for treatment effect in respective meta-analysis (see below
    for details).}
  \item{weight.pooled}{Total weight in respective meta-analysis (see below
    for details).}
  \item{Z.pooled}{Z-score for test of overall treatment effect in
    respective meta-analysis (see below for details).}
  \item{pval.pooled}{P-value for test of overall treatment effect in
    respective meta-analysis (see below for details).}
  \item{Q}{Heterogeneity statistic Q in respective meta-analysis
    (see below for details).}
  \item{pval.Q}{P-value of heterogeneity statistic Q in
    respective meta-analysis (see below for details).}
  \item{I2}{Heterogeneity statistic I2 in respective
    meta-analysis (see below for details).}
  \item{tau2}{Between-study variance (moment estimator of
    DerSimonian-Laird) in respective meta-analysis.}
  \item{Q.w}{Heterogeneity statistic Q within groups in
    respective meta-analysis (see below for details).}
  \item{pval.Q.w}{P-value of heterogeneity statistic Q within groups in
    respective meta-analysis (see below for details).}
  \item{I2.w}{Heterogeneity statistic I2 within groups in respective
    meta-analysis (see below for details).}
  \item{label.e}{Label for experimental group.}
  \item{label.c}{Label for control group.}
  \item{label.left}{Graph label on left side of forest plot.}
  \item{label.right}{Graph label on right side of forest plot.}
  %  \item{graph.exp}{Graph label for experimental group.}
  %  \item{graph.cont}{Graph label for control group.}
  %  \item{label.exp}{Label for experimental group.}
  %  \item{label.cont}{Label for control group.}
  \item{RR.cochrane}{A logical indicating if 2*\code{incr} instead of
    1*\code{incr} is to be added to \code{n.e} and \code{n.c} in the
    calculation of the risk ratio (i.e., \code{sm="RR"}) for studies
    with a zero cell. This is used in RevMan 5.}
  \item{complab}{Comparison label.}
}

\references{  
  \emph{Review Manager (RevMan)} [Computer program]. Version
  5.3. Copenhagen: The Nordic Cochrane Centre, The Cochrane
  Collaboration, 2014.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metagen}}, \code{\link{metacr}}}

\examples{
# Locate export data file "Fleiss93_CR.csv"
# in sub-directory of package "meta"
#
filename <- system.file("data/Fleiss93_CR.csv.gz", package = "meta")
#
Fleiss93_CR <- read.rm5(filename)

# Same result as R command example(Fleiss93):
#
metacr(Fleiss93_CR)

# Same result as R command example(Fleiss93cont):
#
metacr(Fleiss93_CR, 1, 2)
}
\keyword{datagen}
