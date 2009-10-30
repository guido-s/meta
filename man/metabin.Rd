\name{metabin}

\alias{metabin}

\title{Meta-analysis of binary outcome data}

\description{
  Calculation of fixed and random effects estimates (relative risk, odds
  ratio, risk difference or arcsine difference) for meta-analyses with
  binary outcome data. Mantel-Haenszel, inverse variance and Peto method
  are available for pooling.
}

\usage{
metabin(event.e, n.e, event.c, n.c, studlab,
        data = NULL, subset = NULL, method = "MH",
        sm = ifelse(!is.na(charmatch(method, c("Peto", "peto"), nomatch = NA)), "OR", "RR"),
        incr = 0.5, allincr = FALSE, addincr = FALSE, allstudies = FALSE,
        MH.exact = FALSE, RR.cochrane = FALSE,
        level = 0.95, level.comb = level,
        comb.fixed=TRUE, comb.random=TRUE,
        title="", complab="", outclab="",
        label.e="Experimental", label.c="Control",
        byvar, bylab, print.byvar=TRUE,
        warn = TRUE)
}

\arguments{
  \item{event.e}{Number of events in experimental group.}
  \item{n.e}{Number of observations in experimental group.}
  \item{event.c}{Number of events in control group.}
  \item{n.c}{Number of observations in control group.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information,
    i.e., event.e, n.e, event.c, and n.c.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{method}{A character string indicating which method is to be used
    for pooling of studies. One of \code{"Inverse"}, \code{"MH"}, or
    \code{"Peto"}, can be abbreviated.}
  \item{sm}{A character string indicating which summary measure
    (\code{"RR"}, \code{"OR"}, \code{"RD"}, or \code{"AS"}) is to be used
    for pooling of studies, see Details.}
  \item{incr}{Could be either a numerical value which is added to each cell
  frequency for studies with a zero cell count or the character string
  \code{"TA"} which stands for treatment arm continuity correction, see
  Details.}
  \item{allincr}{A logical indicating if \code{incr} is added to each
    cell frequency of all studies if at least one study has a zero cell
    count. If false, \code{incr} is added only to each cell frequency of
    studies with a zero cell count.}
  \item{addincr}{A logical indicating if \code{incr} is added to each cell
    frequency of all studies irrespective of zero cell counts.}
  \item{allstudies}{A logical indicating if studies with zero or all
    events in both groups are to be included in the meta-analysis
    (applies only if sm = \code{"RR"} or \code{"OR"}).}
  \item{MH.exact}{A logical indicating if \code{incr} is not to be added
    to all cell frequencies for studies with a zero cell count to
    calculate the pooled estimate based on the Mantel-Haenszel method.}
  \item{RR.cochrane}{A logical indicating if 2*\code{incr} instead of
    1*\code{incr} is to be added to \code{n.e} and \code{n.c} in the
    calculation of the relative risk (i.e., \code{sm="RR"}) for studies
    with a zero cell. This is used in RevMan 5, the
    Cochrane Collaboration's program for preparing and maintaining
    Cochrane reviews.}
  \item{level}{The level used to calculate confidence intervals for
    individual studies.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{comb.fixed}{A logical indicating whether a fixed effect
    meta-analysis should be conducted.}
  \item{comb.random}{A logical indicating whether a random effects
    meta-analysis should be conducted.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{complab}{Comparison label.}
  \item{outclab}{Outcome label.}
  \item{label.e}{Label for experimental group.}
  \item{label.c}{Label for control group.}
  \item{byvar}{An optional vector containing grouping information (must
    be of same length as \code{event.e}).}
  \item{bylab}{A character string with a label for the grouping variable.}
  \item{print.byvar}{A logical indicating whether the name of the grouping
    variable should be printed in front of the group labels.}
  \item{warn}{A logical indicating whether the addition of \code{incr}
    to studies with zero cell frequencies should result in a warning.}
}

\details{
  Treatment estimates and standard errors are calculated for each
  study. The following measures of treatment effect are available:
  relative risk (if \code{sm="RR"}), odds ratio (\code{sm="OR"}), risk
  difference (\code{sm="RD"}), and arcsine difference (\code{sm="AS"}).

  For studies with a zero cell count, by default, 0.5 is added to
  all cell frequencies of these studies; if \code{incr} is
  \code{"TA"} a treatment arm continuity correction is used instead
  (Sweeting et al., 2004; Diamond et al., 2007). Treatment estimates
  and standard errors are only calculated for studies with zero or all
  events in both groups if \code{allstudies} is \code{TRUE}.
  
  Internally, both fixed effect and random effects models are calculated
  regardless of values choosen for arguments \code{comb.fixed} and
  \code{comb.random}. Accordingly, the estimate for the random effects
  model can be extracted from component \code{TE.random} of an object
  of class \code{"meta"} even if \code{comb.random=FALSE}. However, all
  functions in R package \code{meta} will adequately consider the values
  for \code{comb.fixed} and \code{comb.random}. E.g. function
  \code{\link{print.meta}} will not print results for the random effects
  model if \code{comb.random=FALSE}.

  By default, both fixed effect and random effects models are considered
  (arguments \code{comb.fixed=TRUE} and \code{comb.random=TRUE}). If
  \code{method} is \code{"MH"} (default), the Mantel-Haenszel method is
  used to calculate the fixed effect estimate; if \code{method} is
  \code{"Inverse"}, inverse variance weighting is used for
  pooling; finally, if \code{method} is \code{"Peto"}, the Peto method
  is used for pooling. The DerSimonian-Laird estimate is used in the
  random effects model. For the Peto method, Peto's log odds ratio,
  i.e. \code{(O-E)/V} and its standard error \code{sqrt(1/V)} with
  \code{O-E} and \code{V} denoting "Observed minus Expected" and
  "V", are utilised in the random effects model. Accordingly, results of
  a random effects model using \code{sm="Peto"} can be (slightly)
  different to results from a random effects model using \code{sm="MH"} or
  \code{sm="Inverse"}.
  
  For the Mantel-Haenszel method, by default (if \code{MH.exact} is
  FALSE), 0.5 is added to all cell frequencies of a study with a zero cell
  count in the calculation of the pooled estimate. This approach is also
  used in other software, e.g. RevMan 5 and the Stata procedure metan.
  According to Fleiss (in Cooper & Hedges, 1994), there is no need to
  add 0.5 to a cell frequency of zero to calculate the Mantel-Haenszel
  estimate and he advocates the exact method
  (\code{MH.exact}=TRUE). Note, the estimate based on the exact method
  is not defined if the number of events is zero in all studies either
  in the experimental or control group.   
}

\value{
  An object of class \code{c("metabin", "meta")} with corresponding
  \code{print}, \code{summary}, \code{plot} function. The object is a
  list containing the following components:
  \item{event.e, n.e, event.c, n.c, studlab,}{}
  \item{sm, method, incr, allincr, addincr, }{}
  \item{allstudies, MH.exact, RR.cochrane, warn,}{}
  \item{level, level.comb, comb.fixed, comb.random,}{}
  \item{byvar, bylab, print.byvar}{As defined above.}
  \item{TE, seTE}{Estimated treatment effect and standard error of individual studies.}
  \item{w.fixed, w.random}{Weight of individual studies (in fixed and
    random effects model).}
  \item{TE.fixed, seTE.fixed}{Estimated overall treatment effect and
    standard error (fixed effect model).}
  \item{TE.random, seTE.random}{Estimated overall treatment effect and
    standard error (random effects model).}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{Q}{Heterogeneity statistic Q.}
  \item{tau}{Square-root of between-study variance (moment estimator of
    DerSimonian-Laird).}
  \item{Q.CMH}{Cochrane-Mantel-Haenszel heterogeneity statistic.}
  \item{incr.e, incr.c}{Increment added to cells in the experimental and
    control group, respectively}
  \item{sparse}{Logical flag indicating if any study included in
    meta-analysis has any zero cell frequencies.}
  \item{call}{Function call.}
}

\references{
  Cooper H & Hedges LV (1994),
  \emph{The Handbook of Research Synthesis}.
  Newbury Park, CA: Russell Sage Foundation.

  Diamond GA, Bax L, Kaul S (2007,
  Uncertain Effects of Rosiglitazone on the Risk for Myocardial
  Infarction and Cardiovascular Death.
  \emph{Annals of Internal Medicine}, \bold{147}, 578--581.

  DerSimonian R & Laird N (1986),
  Meta-analysis in clinical trials. \emph{Controlled Clinical Trials},
  \bold{7}, 177--188.

  Fleiss JL (1993),
  The statistical basis of meta-analysis.
  \emph{Statistical Methods in Medical Research}, \bold{2}, 121--145.

  Greenland S & Robins JM (1985),
  Estimation of a common effect parameter from sparse follow-up data.
  \emph{Biometrics}, \bold{41}, 55--68.
  
  \emph{Review Manager (RevMan)} [Computer program]. Version 5.0. The
  Nordic Cochrane Centre, The Cochrane Collaboration, 2008.
  
  Ruecker G, Schwarzer G, Carpenter JR (2008)
  Arcsine test for publication bias in meta-analyses with binary
  outcomes. \emph{Statistics in Medicine}, \bold{27},746--763.
  
  StataCorp. 2009.
  \emph{Stata Statistical Software: Release 11}. College Station, TX:
  Stata Corporation.

  Sweeting MJ, Sutton AJ, Lambert PC (2004),
  What to add to nothing? Use and avoidance of continuity corrections
  in meta-analysis of sparse data.
  \emph{Statistics in Medicine}, \bold{23}, 1351--1375.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{funnel}}, \code{\link{metabias}}, \code{\link{metacont}}, \code{\link{metagen}}, \code{\link{print.meta}}}

\examples{
metabin(10, 20, 15, 20, sm="OR")

##
## Different results:
##
metabin(0, 10, 0, 10, sm="OR")
metabin(0, 10, 0, 10, sm="OR", allstudies=TRUE)


data(Olkin95)

meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Olkin95, subset=c(41,47,51,59),
                 sm="RR", meth="I")
summary(meta1)
funnel(meta1)

meta2 <- metabin(event.e, n.e, event.c, n.c,
                 data=Olkin95, subset=Olkin95$year<1970,
                 sm="RR", meth="I")
summary(meta2)
}
\keyword{htest}
