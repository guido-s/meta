#' meta: Brief overview of methods and general hints
#' 
#' @description
#' R package \bold{meta} is a user-friendly general package providing
#' standard methods for meta-analysis and supporting Schwarzer et
#' al. (2015),
#' \url{https://link.springer.com/book/10.1007/978-3-319-21416-0}.
#' 
#' @details
#' R package \bold{meta} (Schwarzer, 2007; Balduzzi et al., 2019)
#' provides the following statistical methods for meta-analysis.
#' \enumerate{
#' \item Fixed effect and random effects model:
#' \itemize{
#'  \item Meta-analysis of continuous outcome data (\code{\link{metacont}})
#'  \item Meta-analysis of binary outcome data (\code{\link{metabin}})
#'  \item Meta-analysis of incidence rates (\code{\link{metainc}})
#'  \item Generic inverse variance meta-analysis (\code{\link{metagen}})
#'  \item Meta-analysis of single correlations (\code{\link{metacor}})
#'  \item Meta-analysis of single means (\code{\link{metamean}})
#'  \item Meta-analysis of single proportions (\code{\link{metaprop}})
#'  \item Meta-analysis of single incidence rates (\code{\link{metarate}})
#' }
#' \item Several plots for meta-analysis:
#' \itemize{
#'  \item Forest plot (\code{\link{forest.meta}}, \code{\link{forest.metabind}})
#'  \item Funnel plot (\code{\link{funnel.meta}})
#'  \item Galbraith plot / radial plot (\code{\link{radial.meta}})
#'  \item L'Abbe plot for meta-analysis with binary outcome data
#'   (\code{\link{labbe.metabin}}, \code{\link{labbe.default}})
#'  \item Baujat plot to explore heterogeneity in meta-analysis
#'   (\code{\link{baujat.meta}})
#'  \item Bubble plot to display the result of a meta-regression
#'   (\code{\link{bubble.metareg}})
#' }
#' \item Statistical tests for funnel plot asymmetry
#'  (\code{\link{metabias.meta}}, \code{\link{metabias.rm5}}) and
#'  trim-and-fill method (\code{\link{trimfill.meta}},
#'  \code{\link{trimfill.default}}) to evaluate bias in meta-analysis
#' \item Cumulative meta-analysis (\code{\link{metacum}}) and
#'   leave-one-out meta-analysis (\code{\link{metainf}})
#' \item Meta-regression (\code{\link{metareg}})
#' \item Import data from Review Manager 5 (\code{\link{read.rm5}});
#'   see also \code{\link{metacr}} to conduct meta-analysis for a
#'   single comparison and outcome from a Cochrane review
#' \item Prediction interval for the treatment effect of a new study
#'  (Higgins et al., 2009); see argument \code{prediction} in
#'  meta-analysis functions, e.g., \code{\link{metagen}}
#' \item Hartung-Knapp method for random effects meta-analysis
#'  (Hartung & Knapp, 2001a,b); see argument \code{hakn} in
#'  meta-analysis functions, e.g., \code{\link{metagen}}
#' \item Various estimators for the between-study variance
#'  \eqn{\tau^2} in a random effects model (Veroniki et al., 2016);
#'  see argument \code{method.tau} in meta-analysis functions, e.g.,
#'  \code{\link{metagen}}
#' \item Generalised linear mixed models (\code{\link{metabin}},
#'   \code{\link{metainc}}, \code{\link{metaprop}}, and
#'   \code{\link{metarate}})
#' }
#' 
#' The following more advanced statistical methods are provided by
#' add-on R packages:
#' \itemize{
#' \item Frequentist methods for network meta-analysis (R package
#'   \bold{netmeta})
#' \item Advanced methods to model and adjust for bias in
#'   meta-analysis (R package \bold{metasens})
#' }
#' 
#' Results of several meta-analyses can be combined with
#' \code{\link{metabind}}. This is, for example, useful to generate a
#' forest plot with results of subgroup analyses.
#' 
#' See \code{\link{settings.meta}} to learn how to print and specify
#' default meta-analysis methods used during your R session. For
#' example, the function can be used to specify general settings:
#' \itemize{
#' \item \code{settings.meta("revman5")}
#' \item \code{settings.meta("jama")}
#' \item \code{settings.meta("iqwig5")}
#' \item \code{settings.meta("iqwig6")}
#' \item \code{settings.meta("geneexpr")}
#' }
#' 
#' The first command can be used to reproduce meta-analyses from
#' Cochrane reviews conducted with \emph{Review Manager 5} (RevMan 5,
#' \url{https://training.cochrane.org/online-learning/core-software-cochrane-reviews/revman})
#' and specifies to use a RevMan 5 layout in forest plots.
#'
#' The second command can be used to generate forest plots following
#' instructions for authors of the \emph{Journal of the American
#' Medical Association}
#' (\url{https://jamanetwork.com/journals/jama/pages/instructions-for-authors/}). Study
#' labels according to JAMA guidelines can be generated using
#' \code{\link{labels.meta}}.
#'
#' The next two commands implement the recommendations of the
#' Institute for Quality and Efficiency in Health Care (IQWiG),
#' Germany accordinging to General Methods 5 and 6, respectively
#' (\url{https://www.iqwig.de/en/about-us/methods/methods-paper/}).
#'
#' The last setting can be used to print p-values in scientific
#' notation and to suppress the calculation of confidence intervals
#' for the between-study variance.
#' 
#' In addition, \code{\link{settings.meta}} can be used to change
#' individual settings. For example, the following R command specifies
#' the use of the Hartung-Knapp and Paule-Mandel methods, and the
#' printing of prediction intervals in the current R session for any
#' meta-analysis generated after execution of this command:
#' \itemize{
#' \item \code{settings.meta(hakn=TRUE, method.tau="PM", prediction=TRUE)}
#' }
#' 
#' Type \code{help(package = "meta")} for a listing of R functions and
#' datasets available in \bold{meta}.
#' 
#' Balduzzi et al. (2019) is the preferred citation in publications
#' for \bold{meta}. Type \code{citation("meta")} for a BibTeX entry of
#' this publication.
#' 
#' To report problems and bugs
#' \itemize{
#' \item type \code{bug.report(package = "meta")} if you do not use
#'   RStudio,
#' \item send an email to Guido Schwarzer
#'   \email{sc@imbi.uni-freiburg.de} if you use RStudio.
#' }
#' 
#' The development version of \bold{meta} is available on GitHub
#' \url{https://github.com/guido-s/meta/}.
#'
#' @note
#' R package \bold{meta} imports R functions from \bold{metafor}
#' (Viechtbauer, 2010) to
#' \itemize{
#' \item estimate the between-study variance \eqn{\tau^2},
#' \item conduct meta-regression,
#' \item estimate generalised linear mixed models.
#' }
#' 
#' @name meta-package
#' 
#' @aliases meta-package meta
#' 
#' @docType package
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @references
#' Balduzzi S, Rücker G, Schwarzer G (2019):
#' How to perform a meta-analysis with R: a practical tutorial.
#' \emph{Evidence-Based Mental Health},
#' \bold{22}, 153--160
#' 
#' Hartung J, Knapp G (2001a):
#' On tests of the overall treatment effect in meta-analysis with
#' normally distributed responses.
#' \emph{Statistics in Medicine},
#' \bold{20}, 1771--82
#' 
#' Hartung J, Knapp G (2001b):
#' A refined method for the meta-analysis of controlled clinical
#' trials with binary outcome.
#' \emph{Statistics in Medicine},
#' \bold{20}, 3875--89
#' 
#' Higgins JPT, Thompson SG, Spiegelhalter DJ (2009):
#' A re-evaluation of random-effects meta-analysis.
#' \emph{Journal of the Royal Statistical Society: Series A},
#' \bold{172}, 137--59
#' 
#' Schwarzer G (2007):
#' meta: An R package for meta-analysis.
#' \emph{R News},
#' \bold{7}, 40--5
#' 
#' Schwarzer G, Carpenter JR and Rücker G (2015):
#' \emph{Meta-Analysis with R (Use-R!)}.
#' Springer International Publishing, Switzerland
#'
#' Veroniki AA, Jackson D, Viechtbauer W, Bender R, Bowden J, Knapp G,
#' et al. (2016):
#' Methods to estimate the between-study variance and its uncertainty
#' in meta-analysis.
#' \emph{Research Synthesis Methods},
#' \bold{7}, 55--79 
#' 
#' Viechtbauer W (2010):
#' Conducting Meta-Analyses in R with the metafor Package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
#'
#' @keywords package
#'
#' @importFrom grid arrow gpar grid.draw grid.layout grid.lines
#'   grid.newpage grid.polygon grid.rect grid.text grid.xaxis textGrob
#'   popViewport pushViewport viewport unit unit.c convertX
#'
#' @importFrom grDevices gray gray.colors
#'
#' @importFrom graphics abline axis box mtext lines par plot points
#'   polygon text
#' 
#' @importFrom stats as.formula binom.test coef cor lm pchisq pnorm pt
#'   qlogis qnorm qt runif update var weighted.mean weights
#'
#' @importFrom utils count.fields read.table assignInNamespace
#'   getFromNamespace packageDescription packageVersion head tail
#'
#' @importFrom metafor forest funnel funnel.default baujat labbe
#'   radial trimfill rma.uni rma.glmm rma.mv predict.rma
#'   confint.rma.uni confint.rma.mv escalc regtest to.long
#'
#' @importFrom lme4 glmer
#'
#' @importFrom CompQuadForm farebrother
#'
#' @export forest funnel baujat labbe radial trimfill


NULL
