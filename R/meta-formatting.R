#' Formatting arguments for printouts and forest plots
#'
#' @description
#' Formatting arguments for printouts and forest plots.
#'
#' @details
#' Formatting arguments for printouts can be provided in the following named
#' lists, with availability shown below:
#' 
#' \code{print}, \code{digits}.
#' 
#' Furthermore, formatting arguments for forest plots can be
#' provided in the following named lists:
#' 
#' \code{test}, \code{label}, \code{rob}, \code{fs}, \code{ff},
#' \code{just}, \code{fill}, \code{col}, \code{print}, \code{type},
#' \code{pooled}, \code{text}, \code{lty}, \code{lwd}, \code{digits},
#' \code{calcwidth}, \code{colgap}, \code{cid}, \code{addrow}.
#' 
#' List elements available in these named lists are described in the following
#' subsections.
#' 
#' The corresponding individual arguments are also supported. Their names are
#' formed by combining the list name and element name, such as \code{fs.main}
#' or \code{col.diamond}. Supplying the same setting both as a list element
#' and as an individual argument is an error, even if the values agree.
#'
#' In \code{\link{forest.meta}}, the primary value can also be supplied directly
#' without wrapping it in a named list for \code{addrow}, \code{cid},
#' \code{colgap}, \code{digits}, \code{fill}, \code{just}, and \code{lwd}.
#' The \code{rob} argument also accepts a risk-of-bias assessment directly.
#'
#' In \code{\link{forest.metacum}} or \code{\link{forest.metainf}}, the primary
#' value can also be supplied directly without wrapping it in a named list for
#' \code{col} and \code{digits}.
#' A scalar \code{digits} in print methods sets \code{digits$digits}, except
#' in \code{\link{print.cidprop}}, where it sets \code{digits$cid}.
#'
#' Most individual arguments can also be set globally through
#' \code{\link{settings.meta}}. The corresponding defaults are then used by
#' the named lists.
#'
#' Function keys in the tables:
#' \describe{
#' \item{M}{\code{\link{forest.meta}}.}
#' \item{B}{\code{\link{forest.metabind}}.}
#' \item{C}{\code{\link{forest.metacum}} and \code{\link{forest.metainf}}.}
#' \item{P}{\code{\link{print.metabind}}.}
#' \item{N}{\code{\link{print.meta}} and \code{\link{print.summary.meta}}.}
#' \item{S}{\code{\link{print.summary.meta}} only.}
#' \item{R}{\code{\link{print.metacum}} and \code{\link{print.metainf}}.}
#' \item{D}{\code{\link{print.cidprop}}.}
#' \item{L}{\code{\link{print.blup.meta}}; also
#'   \code{\link{print.estimates.blup.meta}} for \code{digits}.}
#' }
#' A star, as in \bold{B*} or \bold{C*}, denotes forwarding through
#' \code{...} to \code{forest.meta}. Unmarked keys indicate direct processing;
#' C uses the same implementation for both methods.
#'
#' Availability does not imply that every setting affects every plot:
#' forwarded arguments remain subject to the wrapper's documented restrictions
#' and the results displayed. In particular, C supplies study symbols and
#' disables overall-effect tests; B does not provide tests between input
#' meta-analyses and fixes \code{pooled.totals = FALSE}.
#'
#' Each subsection first describes elements and availability, then gives
#' one table of general defaults. For elements supported by M, native
#' \code{forest.meta} defines the baseline; other elements have the general
#' defaults shown. Individual function documentation describes exceptions.
#' Named-list and individual arguments use the same defaults.
#'
#' \code{gs("name")} denotes the current
#' setting from \code{\link{settings.meta}}, \code{x} the object supplied,
#' and other names the corresponding function arguments. Defaults depending
#' on the object, layout, or displayed results are indicated as such;
#' \code{NULL} can request automatic selection rather than suppress output.
#' In default tables, \code{*} stands for the element name in a setting or
#' object field. CI denotes a confidence interval; CID denotes a clinically
#' important difference.
#' 
#' \subsection{Tests (test.*)}{
#' Logical values control whether the
#' corresponding tests are printed.
#'
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{overall} \tab Overall effects \tab M, B \cr
#' \code{overall.common} \tab Common effect \tab M, B \cr
#' \code{overall.random} \tab Random effects \tab M, B \cr
#' \code{subgroup} \tab Subgroup differences \tab M, B, C*, N \cr
#' \code{subgroup.common} \tab Common-effect differences \tab M, B, C*, N \cr
#' \code{subgroup.random} \tab Random-effects differences \tab M, B, C*, N \cr
#' \code{effect.subgroup} \tab Effects within subgroups \tab M, B, C* \cr
#' \code{effect.subgroup.common} \tab Common subgroup effects \tab M, B, C* \cr
#' \code{effect.subgroup.random} \tab Random subgroup effects \tab M, B, C*
#' }
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{General default} \cr
#' \code{overall} \tab Automatic (see below) \cr
#' \code{overall.common} \tab Automatic (see below) \cr
#' \code{overall.random} \tab Automatic (see below) \cr
#' \code{subgroup} \tab Automatic (see below) \cr
#' \code{subgroup.common} \tab Automatic (see below) \cr
#' \code{subgroup.random} \tab Automatic (see below) \cr
#' \code{effect.subgroup} \tab \code{gs("test.effect.subgroup")} \cr
#' \code{effect.subgroup.common} \tab \code{NULL} (automatic) \cr
#' \code{effect.subgroup.random} \tab \code{NULL} (automatic)
#' }
#' The overall-test default is \code{gs("test.overall")}; if NULL, it is TRUE
#' for BMJ, RevMan5, JAMA, and meta9 layouts, and FALSE otherwise. Model-specific
#' defaults also depend on \code{common}, \code{random}, and \code{overall}.
#' Subgroup tests default to \code{x$test.subgroup}, or TRUE if NULL;
#' model-specific defaults also depend on \code{common} and \code{random}.
#' C fixes overall-effect tests to FALSE; do not supply these elements.
#' In \code{\link{print.meta}}, \code{test} supports only \code{subgroup},
#' \code{subgroup.common}, and \code{subgroup.random}.
#' These settings can also be passed through \code{\link{print.summary.meta}}.
#' }
#'
#' \subsection{Labels (label.*)}{
#' Labels identify groups, attached study-data columns, graph regions, column
#' headings, and test results.
#'
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{e} \tab Experimental group \tab M, B, C* \cr
#' \code{c} \tab Control group \tab M, B, C* \cr
#' \code{e.attach} \tab Attached experimental column \tab M, B, C* \cr
#' \code{c.attach} \tab Attached control column \tab M, B, C* \cr
#' \code{left}, \code{right} \tab Either side of null effect \tab M, B, C* \cr
#' \code{n}, \code{events} \tab Sample sizes / events \tab M, B, C* \cr
#' \code{mean}, \code{sd} \tab Mean / standard deviation \tab M, B, C* \cr
#' \code{cor} \tab Correlation \tab M, B, C* \cr
#' \code{time}, \code{pval} \tab Time / p-value \tab M, B, C* \cr
#' \code{tau2}, \code{tau}, \code{I2} \tab Heterogeneity \tab M, B, C* \cr
#' \code{cluster} \tab Cluster variable \tab M, B, C* \cr
#' \code{cycles} \tab Cycle count \tab M, B, C* \cr
#' \code{test.overall.common} \tab Overall common-effect test \tab M, B, C* \cr
#' \code{test.overall.random} \tab Overall random-effects test \tab M, B, C* \cr
#' \code{test.subgroup.common} \tab Common subgroup test \tab M, B, C* \cr
#' \code{test.subgroup.random} \tab Random subgroup test \tab M, B, C* \cr
#' \code{test.effect.subgroup.common} \tab Common subgroup-effect test \tab M, B, C* \cr
#' \code{test.effect.subgroup.random} \tab Random subgroup-effect test \tab M, B, C*
#' }
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{General default} \cr
#' \code{e}, \code{c} \tab Corresponding \code{x$label.*} value \cr
#' \code{left}, \code{right} \tab Corresponding \code{x$label.*} value \cr
#' \code{e.attach}, \code{c.attach} \tab Corresponding \code{gs("label.*")} setting \cr
#' \code{n}, \code{events} \tab Corresponding \code{gs("label.*")} setting \cr
#' \code{mean}, \code{sd} \tab Corresponding \code{gs("label.*")} setting \cr
#' \code{cor}, \code{time} \tab Corresponding \code{gs("label.*")} setting \cr
#' \code{pval}, \code{tau2} \tab Corresponding \code{gs("label.*")} setting \cr
#' \code{tau}, \code{I2} \tab Corresponding \code{gs("label.*")} setting \cr
#' \code{cluster}, \code{cycles} \tab Corresponding \code{gs("label.*")} setting \cr
#' All \code{test.*} elements \tab \code{NULL} (test-specific labels)
#' }
#' }
#'
#' \subsection{Risk-of-bias settings (rob.*)}{
#' In \code{\link{forest.meta}}, settings for risk-of-bias assessments can be
#' supplied in the named list \code{rob}.
#'
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{rob} \tab Risk-of-bias assessment \tab M, B*, C* \cr
#' \code{text} \tab Table heading \tab M, B*, C* \cr
#' \code{legend} \tab Domain legend \tab M, B*, C* \cr
#' \code{only} \tab Risk of bias alone on right \tab M, B*, C*
#' }
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{General default} \cr
#' \code{rob} \tab \code{x$rob} \cr
#' \code{text} \tab \code{"Risk of Bias"} \cr
#' \code{legend} \tab \code{TRUE} \cr
#' \code{only} \tab \code{FALSE}
#' }
#' }
#'
#' \subsection{Font sizes (fs.*) and font faces (ff.*)}{
#' Numeric font sizes, in points. The corresponding elements can be supplied
#' in the named list \code{fs}. Font faces accepted by \code{\link[grid]{gpar}},
#' such as \code{"plain"}, \code{"bold"}, and \code{"italic"}. The corresponding
#' elements can be supplied in the named list \code{ff}.
#'
#' Both lists can contain the following elements:
#'
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{main}, \code{heading} \tab Title / column headings \tab M, B*, C* \cr
#' \code{common}, \code{random} \tab Pooled results \tab M, B*, C* \cr
#' \code{predict} \tab Prediction intervals \tab M, B*, C* \cr
#' \code{common.labels} \tab Common-effect labels \tab M, B*, C* \cr
#' \code{random.labels} \tab Random-effects labels \tab M, B*, C* \cr
#' \code{predict.labels} \tab Prediction-interval labels \tab M, B*, C* \cr
#' \code{study} \tab Study results \tab M, B*, C* \cr
#' \code{study.labels} \tab Study labels \tab M, B*, C* \cr
#' \code{hetstat} \tab Heterogeneity statistics \tab M, B*, C* \cr
#' \code{test.overall} \tab Overall-effect tests \tab M, B*, C* \cr
#' \code{test.subgroup} \tab Subgroup-difference tests \tab M, B*, C* \cr
#' \code{test.effect.subgroup} \tab Subgroup-effect tests \tab M, B*, C* \cr
#' \code{addline} \tab Additional lines \tab M, B*, C* \cr
#' \code{axis} \tab Axis text \tab M, B*, C* \cr
#' \code{smlab} \tab Summary measure \tab M, B*, C* \cr
#' \code{xlab} \tab x-axis label \tab M, B*, C* \cr
#' \code{lr} \tab Left / right labels \tab M, B*, C* \cr
#' \code{rob} \tab Risk-of-bias text \tab M, B*, C* \cr
#' \code{rob.symbols} \tab Risk-of-bias symbols \tab M, B*, C* \cr
#' \code{details} \tab Method details \tab M, B*, C*
#' }
#'
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Size default (fs)} \tab \bold{Face default (ff)} \cr
#' \code{heading} \tab \code{fontsize} \tab \code{"bold"} \cr
#' \code{study}, \code{study.labels} \tab \code{fontsize} \tab \code{"plain"} \cr
#' \code{smlab}, \code{xlab} \tab \code{fontsize} \tab Corresponding \code{gs("ff.*")} setting \cr
#' \code{lr} \tab \code{fontsize} \tab \code{gs("ff.lr")} \cr
#' \code{axis} \tab \code{fontsize - 1} \tab \code{gs("ff.axis")} \cr
#' \code{details} \tab \code{fontsize - 1} \tab \code{"plain"} \cr
#' \code{rob} \tab \code{gs("fs.rob")} \tab \code{"plain"} \cr
#' \code{rob.symbols} \tab \code{gs("fs.rob.symbols")} \tab \code{"bold"} \cr
#' Remaining elements \tab \code{gs("fs.*")} \tab \code{gs("ff.*")}
#' }
#' }
#'
#' \subsection{Justifications (just.*)}{
#' Possible values are \code{"left"}, \code{"right"}, and \code{"center"}.
#'
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{just} \tab Ordinary columns \tab M, B*, C* \cr
#' \code{main} \tab Main title \tab M, B*, C* \cr
#' \code{studlab} \tab Study labels \tab M, B*, C* \cr
#' \code{addcols} \tab Additional columns \tab M, B* \cr
#' \code{addcols.left} \tab Left additional columns \tab M, B*, C* \cr
#' \code{addcols.right} \tab Right additional columns \tab M, B*, C* \cr
#' \code{label.e}, \code{label.c} \tab Group headings \tab M, B*, C* \cr
#' \code{rob} \tab Risk-of-bias columns \tab M, B*, C*
#' }
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{General default} \cr
#' \code{just} \tab \code{"left"} for JAMA, otherwise \code{"right"} \cr
#' \code{label.e}, \code{label.c} \tab \code{"left"} for JAMA, otherwise \code{"right"} \cr
#' \code{main}, \code{studlab} \tab Corresponding \code{gs("just.*")} setting \cr
#' \code{addcols} \tab \code{gs("just.addcols")} \cr
#' \code{addcols.left}, \code{addcols.right} \tab \code{gs("just.addcols")} \cr
#' \code{rob} \tab \code{"left"}
#' }
#' C fixes \code{just.addcols = "right"}; do not supply \code{just$addcols}.
#' }
#'
#' \subsection{Fill settings (fill.*)}{
#' Colours used to fill plot regions and symbols.
#'
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{fill} \tab General fill \tab M, B*, C* \cr
#' \code{cid} \tab CID regions \tab M, B*, C* \cr
#' \code{cid.below.null} \tab CID below null effect \tab M, B*, C* \cr
#' \code{cid.above.null} \tab CID above null effect \tab M, B*, C* \cr
#' \code{equi} \tab Between equivalence limits \tab M, B*, C* \cr
#' \code{lower.equi} \tab Lower limit to null effect \tab M, B*, C* \cr
#' \code{upper.equi} \tab Null effect to upper limit \tab M, B*, C*
#' }
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{General default} \cr
#' \code{fill} \tab \code{gs("fill")} \cr
#' \code{cid}, \code{cid.below.null} \tab \code{gs("fill.cid")} \cr
#' \code{cid.above.null} \tab \code{rev(gs("fill.cid"))} \cr
#' \code{equi}, \code{lower.equi} \tab \code{gs("fill.equi")} \cr
#' \code{upper.equi} \tab \code{rev(gs("fill.equi"))}
#' }
#' }
#'
#' \subsection{Colours (col.*)}{
#' Colour values used for plot elements. Vectors can be supplied when colours
#' differ between subgroups or rows.
#'
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{main} \tab Main title \tab M, B* \cr
#' \code{common} \tab Common effect \tab M, B* \cr
#' \code{random} \tab Random effects \tab M, B* \cr
#' \code{cid} \tab Decision lines \tab M, B* \cr
#' \code{label} \tab Axis labels \tab M, B* \cr
#' \code{label.left}, \code{label.right} \tab Left / right labels \tab M, B* \cr
#' \code{study}, \code{study.text} \tab Study text \tab M, B* \cr
#' \code{square}, \code{square.lines} \tab Square fill / border \tab M, B \cr
#' \code{circle}, \code{circle.lines} \tab Circle fill / border \tab M, B \cr
#' \code{inside} \tab Inside study symbols \tab M, B* \cr
#' \code{inside.common} \tab Inside common symbols \tab M, B* \cr
#' \code{inside.random} \tab Inside random symbols \tab M, B* \cr
#' \code{diamond} \tab Diamond fill \tab M, B \cr
#' \code{diamond.common} \tab Common diamond fill \tab M, B \cr
#' \code{diamond.random} \tab Random diamond fill \tab M, B \cr
#' \code{diamond.common.subgroup} \tab Common subgroup fill \tab M, B* \cr
#' \code{diamond.random.subgroup} \tab Random subgroup fill \tab M, B* \cr
#' \code{diamond.lines} \tab Diamond border \tab M, B \cr
#' \code{diamond.lines.common} \tab Common diamond border \tab M, B \cr
#' \code{diamond.lines.random} \tab Random diamond border \tab M, B \cr
#' \code{diamond.lines.common.subgroup} \tab Common subgroup border \tab M, B* \cr
#' \code{diamond.lines.random.subgroup} \tab Random subgroup border \tab M, B* \cr
#' \code{predict} \tab Prediction-interval fill \tab M, B \cr
#' \code{predict.lines} \tab Prediction-interval border \tab M, B \cr
#' \code{predict.subgroup} \tab Subgroup prediction fill \tab M, B* \cr
#' \code{predict.lines.subgroup} \tab Subgroup prediction border \tab M, B* \cr
#' \code{subgroup} \tab Subgroup text \tab M, B* \cr
#' \code{lines}, \code{header.line} \tab Plot / header lines \tab M, B* \cr
#' \code{jama.line} \tab JAMA lines \tab M, B* \cr
#' \code{col} \tab Square results \tab C \cr
#' \code{bg} \tab Symbol fill \tab C \cr
#' \code{border} \tab Symbol border \tab C \cr
#' \code{bg.predict} \tab Prediction-interval fill \tab C \cr
#' \code{border.predict} \tab Prediction-interval border \tab C
#' }
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{General default} \cr
#' \code{main}, \code{common} \tab Corresponding \code{gs("col.*")} setting \cr
#' \code{random}, \code{cid} \tab Corresponding \code{gs("col.*")} setting \cr
#' \code{label}, \code{study} \tab Corresponding \code{gs("col.*")} setting \cr
#' \code{study.text}, \code{square} \tab Corresponding \code{gs("col.*")} setting \cr
#' \code{square.lines}, \code{circle} \tab Corresponding \code{gs("col.*")} setting \cr
#' \code{circle.lines} \tab \code{gs("col.circle")} \cr
#' \code{inside} \tab \code{gs("col.inside")} \cr
#' \code{inside.common}, \code{inside.random} \tab \code{gs("col.inside")} \cr
#' \code{diamond}, \code{diamond.common} \tab \code{gs("col.diamond")} \cr
#' \code{diamond.random} \tab \code{gs("col.diamond")} \cr
#' \code{diamond.common.subgroup} \tab \code{gs("col.diamond")} \cr
#' \code{diamond.random.subgroup} \tab \code{gs("col.diamond")} \cr
#' \code{diamond.lines} \tab \code{gs("col.diamond.lines")} \cr
#' \code{diamond.lines.common} \tab \code{gs("col.diamond.lines")} \cr
#' \code{diamond.lines.random} \tab \code{gs("col.diamond.lines")} \cr
#' \code{diamond.lines.common.subgroup} \tab \code{gs("col.diamond.lines")} \cr
#' \code{diamond.lines.random.subgroup} \tab \code{gs("col.diamond.lines")} \cr
#' \code{predict}, \code{predict.subgroup} \tab \code{gs("col.predict")} \cr
#' \code{predict.lines}, \code{predict.lines.subgroup} \tab \code{gs("col.predict.lines")} \cr
#' \code{subgroup}, \code{jama.line} \tab \code{gs("col.subgroup")} \cr
#' \code{lines}, \code{header.line} \tab \code{gs("col.lines")} \cr
#' \code{label.left} \tab \code{x$col.label.left} \cr
#' \code{label.right} \tab \code{x$col.label.right} \cr
#' \code{col}, \code{bg} \tab \code{NULL} \cr
#' \code{border} \tab \code{NULL} \cr
#' \code{bg.predict}, \code{border.predict} \tab \code{NULL}
#' }
#' B processes its own colour elements and forwards other native elements;
#' C uses only its alternative elements.
#' }
#'
#' \subsection{Printing results (print.*)}{
#' Logical settings controlling displayed results, columns, and labels.
#'
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{I2}, \code{I2.ci} \tab I-squared / CI \tab M, B, C, P, N, R \cr
#' \code{tau2} \tab Tau-squared \tab M, B, C, P, N, R, L \cr
#' \code{tau2.ci} \tab Tau-squared CI \tab M, B, C, P, N, R \cr
#' \code{tau} \tab Tau \tab M, B, C, P, N, R, L \cr
#' \code{tau.ci} \tab Tau CI \tab M, B, C, P, N, R \cr
#' \code{Q} \tab Heterogeneity statistic Q \tab M, B, C, P, N \cr
#' \code{pval.Q} \tab Heterogeneity p-value \tab M, B, C, P \cr
#' \code{Rb} \tab Rb statistic \tab M, N \cr
#' \code{Rb.ci} \tab Rb CI \tab M \cr
#' \code{stat} \tab Effect z- or t-statistic \tab M, B, C, P, R \cr
#' \code{pval} \tab Effect p-value \tab B, C, P \cr
#' \code{Q.subgroup} \tab Subgroup-difference Q \tab M \cr
#' \code{tau2.ci.subgroup} \tab Subgroup tau-squared CI \tab M \cr
#' \code{tau.ci.subgroup} \tab Subgroup tau CI \tab M \cr
#' \code{subgroup.labels} \tab Subgroup labels \tab M \cr
#' \code{subgroup.name} \tab Subgroup variable name \tab M, B, P, N \cr
#' \code{CMH} \tab Cochran-Mantel-Haenszel test \tab N \cr
#' \code{prob} \tab CID benefit / harm probability \tab C, R \cr
#' \code{studlab} \tab Meta-analysis labels \tab B \cr
#' \code{type} \tab Result type \tab B \cr
#' \code{func} \tab Analysis function \tab B, P \cr
#' \code{method} \tab Pooling method \tab B, P \cr
#' \code{method.interval} \tab Interval method \tab B, P \cr
#' \code{method.tau} \tab Tau-squared method \tab B \cr
#' \code{backtransf} \tab Back-transformation function \tab B, P \cr
#' \code{args.backtransf} \tab Back-transformation arguments \tab B, P \cr
#' \code{k} \tab Number of studies \tab B, P \cr
#' \code{sm} \tab Summary measure \tab B, P \cr
#' \code{se} \tab Standard error \tab B, P \cr
#' \code{level.ma} \tab Pooled-effect CI level \tab B, P \cr
#' \code{level.predict} \tab Prediction level \tab B, P \cr
#' \code{level.hetstat} \tab Heterogeneity CI level \tab B, P \cr
#' \code{null.effect} \tab Null effect \tab B, P \cr
#' \code{H} \tab H statistic \tab B, P, N \cr
#' \code{H.ci} \tab H CI \tab B, P \cr
#' \code{Q.b} \tab Subgroup-test statistic \tab B, P \cr
#' \code{pval.Q.b} \tab Subgroup-test p-value \tab B, P
#' }
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{General default} \cr
#' \code{I2}, \code{I2.ci} \tab Corresponding \code{gs("forest.*")} setting \cr
#' \code{tau2}, \code{tau2.ci} \tab Corresponding \code{gs("forest.*")} setting \cr
#' \code{tau}, \code{tau.ci} \tab Corresponding \code{gs("forest.*")} setting \cr
#' \code{Q}, \code{pval.Q} \tab Corresponding \code{gs("forest.*")} setting \cr
#' \code{Rb}, \code{Rb.ci} \tab Corresponding \code{gs("forest.*")} setting \cr
#' \code{stat}, \code{Q.subgroup} \tab Corresponding \code{gs("forest.*")} setting \cr
#' \code{tau2.ci.subgroup} \tab \code{gs("print.tau2.ci.subgroup")} \cr
#' \code{tau.ci.subgroup} \tab \code{gs("print.tau.ci.subgroup")} \cr
#' \code{subgroup.labels} \tab \code{TRUE} \cr
#' \code{subgroup.name} \tab \code{x$print.subgroup.name} \cr
#' \code{CMH} \tab \code{x$print.CMH} \cr
#' \code{studlab}, \code{k} \tab \code{TRUE} \cr
#' \code{pval} \tab \code{TRUE} \cr
#' \code{se}, \code{null.effect} \tab \code{FALSE} \cr
#' \code{H}, \code{H.ci} \tab \code{FALSE} \cr
#' \code{prob} \tab \code{FALSE} \cr
#' \code{type}, \code{func} \tab \code{NULL} (automatic) \cr
#' \code{method}, \code{method.interval} \tab \code{NULL} (automatic) \cr
#' \code{method.tau}, \code{backtransf} \tab \code{NULL} (automatic) \cr
#' \code{args.backtransf}, \code{sm} \tab \code{NULL} (automatic) \cr
#' \code{level.ma}, \code{level.predict} \tab \code{NULL} (automatic) \cr
#' \code{level.hetstat} \tab \code{NULL} (automatic) \cr
#' \code{Q.b}, \code{pval.Q.b} \tab \code{NULL} (automatic)
#' }
#' Confidence-interval switches require the corresponding statistic.
#' B and C use their own print lists; native-only elements are not a shared
#' wrapper interface.
#' \code{Q.b} and \code{pval.Q.b} refer to subgroup differences within each
#' meta-analysis.
#' In P, \code{pval.Q.b} controls the whole subgroup-test block, including
#' the statistic and degrees of freedom; \code{Q.b} has no separate effect.
#' }
#'
#' \subsection{Graphical-type (type.*)}{
#' Character strings specifying how treatment effects and confidence intervals
#' are plotted, such as \code{"square"}, \code{"diamond"}, or
#' \code{"circle"}.
#'
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{study} \tab Individual study results \tab M \cr
#' \code{common} \tab Overall common-effect results \tab M, B \cr
#' \code{random} \tab Overall random-effects results \tab M, B \cr
#' \code{type} \tab General result symbol \tab B \cr
#' \code{predict} \tab Prediction-interval symbol \tab B \cr
#' \code{subgroup} \tab Subgroup results \tab M \cr
#' \code{subgroup.common} \tab Common-effect subgroup results \tab M \cr
#' \code{subgroup.random} \tab Random-effects subgroup results \tab M
#' }
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{General default} \cr
#' \code{study} \tab \code{gs("type.study")} \cr
#' \code{common}, \code{random} \tab \code{gs("type.common")} for both \cr
#' \code{subgroup} \tab Automatic (see below) \cr
#' \code{subgroup.common} \tab Automatic (see below) \cr
#' \code{subgroup.random} \tab Automatic (see below) \cr
#' \code{type}, \code{predict} \tab \code{NULL} (automatic)
#' }
#' Subgroup symbols default to \code{"diamond"} if \code{study.results} is
#' TRUE, and \code{"square"} otherwise.
#' C uses an ordinary character argument \code{type}, not a named list.
#' }
#'
#' \subsection{Pooled-results (pooled.*)}{
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{totals} \tab Sample sizes / totals \tab M, C* \cr
#' \code{events} \tab Event counts \tab M, B*, C* \cr
#' \code{times} \tab Person-time totals \tab M, B*, C*
#' }
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{General default} \cr
#' \code{totals} \tab \code{common | random} \cr
#' \code{events} \tab \code{gs("pooled.events")} \cr
#' \code{times} \tab \code{gs("pooled.times")}
#' }
#' B fixes \code{pooled.totals = FALSE}; do not supply \code{pooled$totals}.
#' }
#'
#' \subsection{Text settings (text.*)}{
#' Character strings used for labels and explanatory text.
#'
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{common} \tab Common-effect label \tab M, B, C* \cr
#' \code{random} \tab Random-effects label \tab M, B, C* \cr
#' \code{w.common} \tab Common weight label \tab M, B, C* \cr
#' \code{w.random} \tab Random weight label \tab M, B, C* \cr
#' \code{predict} \tab Prediction-interval label \tab M, B, C* \cr
#' \code{subgroup} \tab Subgroup label \tab M, B, C* \cr
#' \code{common.w} \tab Common subgroup label \tab M, B, C* \cr
#' \code{random.w} \tab Random subgroup label \tab M, B, C* \cr
#' \code{predict.w} \tab Subgroup prediction label \tab M, B, C* \cr
#' \code{subgroup.nohet} \tab No-heterogeneity text \tab M, B, C* \cr
#' \code{addline1}, \code{addline2} \tab Additional lines \tab M, B, C*
#' }
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{General default} \cr
#' \code{common}, \code{random} \tab Corresponding \code{x$text.*} value \cr
#' \code{w.common}, \code{w.random} \tab Corresponding \code{x$text.*} value \cr
#' \code{predict} \tab \code{x$text.predict} \cr
#' \code{common.w}, \code{random.w} \tab \code{x$text.common}, \code{x$text.random} \cr
#' \code{predict.w} \tab \code{x$text.predict} \cr
#' \code{subgroup} \tab \code{NULL} \cr
#' \code{addline1}, \code{addline2} \tab \code{NULL} \cr
#' \code{subgroup.nohet} \tab \code{gs("text.subgroup.nohet")}
#' }
#' }
#'
#' \subsection{Line-type settings (lty.*)}{
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{common} \tab Common-effect lines \tab M, B*, C* \cr
#' \code{random} \tab Random-effects lines \tab M, B*, C* \cr
#' \code{cid} \tab CID lines \tab M, B*, C*
#' }
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{General default} \cr
#' \code{common}, \code{random} \tab Corresponding \code{gs("lty.*")} setting \cr
#' \code{cid} \tab \code{gs("lty.cid")}
#' }
#' }
#'
#' \subsection{Line-width settings (lwd.*)}{
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{lwd} \tab General line width \tab M, B*, C* \cr
#' \code{square} \tab Square borders \tab M, B*, C* \cr
#' \code{diamond} \tab Diamond borders \tab M, B*, C*
#' }
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{General default} \cr
#' \code{lwd} \tab \code{gs("lwd")} \cr
#' \code{square}, \code{diamond} \tab \code{gs("lwd.square")}, \code{gs("lwd.diamond")}
#' }
#' }
#'
#' \subsection{Digits (digits.*)}{
#' The values specify the minimum number of significant digits.
#'
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{digits} \tab Treatment effects \tab M, B, C, P, N, R, L \cr
#' \code{se} \tab Standard errors \tab M, B, P, S, L \cr
#' \code{stat} \tab z- or t-statistics \tab M, B, C, P, N, R \cr
#' \code{pval} \tab Effect p-values \tab M, B, C, P, N, R \cr
#' \code{pval.Q} \tab Heterogeneity p-values \tab M, B, P, N \cr
#' \code{Q} \tab Heterogeneity statistic Q \tab M, B, P, N \cr
#' \code{df} \tab Degrees of freedom \tab N \cr
#' \code{H} \tab H statistic \tab N \cr
#' \code{tau2} \tab Between-study variance \tab M, B, C, P, N, R, L \cr
#' \code{tau} \tab Between-study SD \tab M, B, C, P, N, R, L \cr
#' \code{I2} \tab I-squared statistic \tab M, B, C, P, N, R \cr
#' \code{weight} \tab Study weights \tab M, B*, S \cr
#' \code{prop} \tab Proportions \tab S \cr
#' \code{mean}, \code{sd} \tab Means / standard deviations \tab M, B* \cr
#' \code{cor}, \code{time} \tab Correlations / follow-up \tab M, B* \cr
#' \code{n}, \code{event} \tab Sample sizes / events \tab M, B* \cr
#' \code{TE} \tab Additional treatment estimates \tab M, B* \cr
#' \code{addcols} \tab Additional columns \tab M, B* \cr
#' \code{addcols.left} \tab Left additional columns \tab M, B* \cr
#' \code{addcols.right} \tab Right additional columns \tab M, B* \cr
#' \code{cid} \tab CID limits \tab C, R, D \cr
#' \code{percent} \tab Percentages \tab C, R, D
#' }
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{General default} \cr
#' \code{digits}, \code{cor} \tab \code{gs("digits.forest")} \cr
#' \code{time}, \code{addcols} \tab \code{gs("digits.forest")} \cr
#' \code{addcols.left}, \code{addcols.right} \tab \code{gs("digits.forest")} \cr
#' \code{se}, \code{stat} \tab Corresponding \code{gs("digits.*")} setting \cr
#' \code{pval}, \code{pval.Q} \tab Corresponding \code{gs("digits.*")} setting \cr
#' \code{Q}, \code{tau2} \tab Corresponding \code{gs("digits.*")} setting \cr
#' \code{df}, \code{H} \tab Corresponding \code{gs("digits.*")} setting \cr
#' \code{tau}, \code{I2} \tab Corresponding \code{gs("digits.*")} setting \cr
#' \code{weight}, \code{mean} \tab Corresponding \code{gs("digits.*")} setting \cr
#' \code{sd} \tab \code{gs("digits.sd")} \cr
#' \code{prop} \tab \code{gs("digits.prop")} \cr
#' \code{n}, \code{event} \tab \code{0} \cr
#' \code{TE} \tab \code{gs("digits.TE.forest")} \cr
#' \code{cid} \tab \code{gs("digits.cid")} \cr
#' \code{percent} \tab \code{1}
#' }
#' }
#'
#' \subsection{Width calculation (calcwidth.*)}{
#' Logical arguments controlling whether generated text is included when
#' calculating column widths.
#'
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{pooled} \tab Pooled results \tab M, B*, C* \cr
#' \code{common} \tab Common-effect text \tab M, B*, C* \cr
#' \code{random} \tab Random-effects text \tab M, B*, C* \cr
#' \code{predict} \tab Prediction-interval text \tab M, B*, C* \cr
#' \code{hetstat} \tab Heterogeneity statistics \tab M, B*, C* \cr
#' \code{tests} \tab Overall / subgroup tests \tab M, B*, C* \cr
#' \code{subgroup} \tab Subgroup labels \tab M, B*, C* \cr
#' \code{addline} \tab Additional lines \tab M, B*, C* \cr
#' \code{details} \tab Method details \tab M, B*, C*
#' }
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{General default} \cr
#' \code{pooled} \tab Automatic (see below) \cr
#' \code{common}, \code{random} \tab Automatic (see below) \cr
#' \code{predict}, \code{hetstat} \tab Corresponding \code{gs("calcwidth.*")} setting \cr
#' \code{tests}, \code{subgroup} \tab Corresponding \code{gs("calcwidth.*")} setting \cr
#' \code{addline}, \code{details} \tab Corresponding \code{gs("calcwidth.*")} setting
#' }
#' Width calculation for pooled, common, and random text is enabled when
#' either \code{common} or \code{random} is TRUE and either \code{overall}
#' is TRUE or \code{subgroup} is non-NULL.
#' }
#'
#' \subsection{Column-gaps (colgap.*)}{
#' Gaps can be specified as character strings or \code{grid::unit} objects.
#'
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{colgap} \tab Between columns \tab M, B*, C* \cr
#' \code{left}, \code{right} \tab Between left / right columns \tab M, B*, C* \cr
#' \code{studlab} \tab After study labels \tab M, B*, C* \cr
#' \code{forest} \tab Beside forest plot \tab M, B*, C* \cr
#' \code{forest.left} \tab Before forest plot \tab M, B*, C* \cr
#' \code{forest.right} \tab After forest plot \tab M, B*, C* \cr
#' \code{rob} \tab Between risk-of-bias columns \tab M, B*, C* \cr
#' \code{rob.overall} \tab Before overall risk of bias \tab M, B*, C*
#' }
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{General default} \cr
#' \code{colgap}, \code{left} \tab \code{gs("colgap")} \cr
#' \code{right}, \code{studlab} \tab \code{gs("colgap")} \cr
#' \code{forest} \tab \code{gs("colgap.forest")} \cr
#' \code{forest.left}, \code{forest.right} \tab \code{gs("colgap.forest")} \cr
#' \code{rob} \tab \code{"1mm"} \cr
#' \code{rob.overall} \tab \code{"2mm"}
#' }
#' }
#'
#' \subsection{Clinically important difference settings (cid.*)}{
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{cid} \tab CID / decision thresholds \tab M, B*, C* \cr
#' \code{below.null} \tab Limits below null effect \tab M, B*, C* \cr
#' \code{above.null} \tab Limits above null effect \tab M, B*, C* \cr
#' \code{pooled.only} \tab Pooled estimates only \tab M, B*, C*
#' }
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{General default} \cr
#' \code{cid} \tab \code{gs("cid")} \cr
#' \code{below.null}, \code{above.null} \tab Corresponding \code{gs("cid.*")} setting \cr
#' \code{pooled.only} \tab \code{gs("cid.pooled.only")}
#' }
#' }
#'
#' \subsection{Additional-row settings (addrow.*)}{
#' \tabular{lll}{
#' \bold{Element} \tab \bold{Description} \tab \bold{Functions} \cr
#' \code{addrow} \tab Empty row above studies \tab M, B*, C* \cr
#' \code{overall} \tab Empty row above overall results \tab M, B*, C* \cr
#' \code{subgroup} \tab Empty row between subgroups \tab M, B*, C*
#' }
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{General default} \cr
#' \code{addrow} \tab \code{gs("addrow")} \cr
#' \code{overall} \tab \code{gs("addrow.overall")} \cr
#' \code{subgroup} \tab \code{gs("addrow.subgroups")}
#' }
#' }
#'
#' @name meta-formatting
#' @seealso \code{\link{forest.meta}}, \code{\link{settings.meta}},
#'   \code{\link[grid]{gpar}}
#' 
#' @examples
#' data(Olkin1995)
#' ma1 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
#'   data = Olkin1995, subset = c(41, 47, 51, 59),
#'   sm = "RR", method = "I",
#'   studlab = paste(author, year))
#' 
#' # Use individual arguments label.left to col.label.right
#' #
#' forest(ma1, layout = "RevMan5", common = FALSE, prediction = TRUE,
#'   label.left = "Favours experimental", label.right = "Favours control",
#'   col.label.left = "green", col.label.right = "red",
#'   main = "Use individual arguments")
#' 
#' \dontrun{
#' # Same forest plot using named lists label and col
#' #
#' forest(ma1, layout = "RevMan5", common = FALSE, prediction = TRUE,
#'   label = list(left = "Favours experimental", right = "Favours control"),
#'   col = list(label.left = "green", label.right = "red"),
#'   main = "Use named lists 'label' and 'col'")
#'
#' # Same forest plot with a mixture
#' #
#' forest(ma1, layout = "RevMan5", common = FALSE, prediction = TRUE,
#'   label.left = "Favours experimental",
#'   col = list(label.left = "green"),
#'   label = list(right = "Favours control"),
#'   col.label.right = "red",
#'   main = "Use mixture of individual arguments and named lists")
#' }

NULL
