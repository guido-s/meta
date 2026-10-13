#' Formatting arguments for forest plots
#'
#' @description
#' Formatting arguments for forest plots created with \code{\link{forest.meta}}.
#'
#' @details
#' In \code{\link{forest.meta}}, the preferred interfaces for specifying related
#' formatting arguments are the following named lists:
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
#' or \code{col.diamond}.
#'
#' For \code{addrow}, \code{cid}, \code{colgap}, \code{digits}, \code{fill},
#' \code{just}, and \code{lwd}, the primary value can also be supplied directly
#' without wrapping it in a named list; see \code{\link{forest.meta}}.
#' The \code{rob} argument also accepts a risk-of-bias assessment directly.
#'
#' Most individual arguments can also be set globally through
#' \code{\link{settings.meta}}. The corresponding defaults are then used by
#' the named lists.
#' 
#' \subsection{Tests (test.*)}{
#' Logical values control whether the corresponding tests are printed.
#'
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{overall} \tab Tests for overall effects. \cr
#' \code{overall.common} \tab Overall common-effect test. \cr
#' \code{overall.random} \tab Overall random-effects test. \cr
#' \code{subgroup} \tab Tests for subgroup differences. \cr
#' \code{subgroup.common} \tab Common-effect subgroup-difference test. \cr
#' \code{subgroup.random} \tab Random-effects subgroup-difference test. \cr
#' \code{effect.subgroup} \tab Tests for effects in subgroups. \cr
#' \code{effect.subgroup.common} \tab Common-effect subgroup-effect tests. \cr
#' \code{effect.subgroup.random} \tab Random-effects subgroup-effect tests.
#' }
#' }
#'
#' \subsection{Labels (label.*)}{
#' Label arguments set group labels, attached study-data columns, graph labels,
#' column headings, and labels preceding test results.
#'
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{e} \tab Experimental group label. \cr
#' \code{c} \tab Control group label. \cr
#' \code{e.attach} \tab Column attached to the experimental group label. \cr
#' \code{c.attach} \tab Column attached to the control group label. \cr
#' \code{left} \tab Label to the left of the null effect. \cr
#' \code{right} \tab Label to the right of the null effect. \cr
#' \code{n} \tab Sample-size label. \cr
#' \code{events} \tab Event-count label. \cr
#' \code{mean} \tab Mean label. \cr
#' \code{sd} \tab Standard-deviation label. \cr
#' \code{cor} \tab Correlation label. \cr
#' \code{time} \tab Time label. \cr
#' \code{pval} \tab P-value label. \cr
#' \code{tau2} \tab Tau-squared label. \cr
#' \code{tau} \tab Tau label. \cr
#' \code{I2} \tab I-squared label. \cr
#' \code{cluster} \tab Cluster-variable label. \cr
#' \code{cycles} \tab Cycle-count label. \cr
#' \code{test.overall.common} \tab Label preceding overall common-effect test. \cr
#' \code{test.overall.random} \tab Label preceding overall random-effects test. \cr
#' \code{test.subgroup.common} \tab Label preceding common-effect subgroup test. \cr
#' \code{test.subgroup.random} \tab Label preceding random-effects subgroup test. \cr
#' \code{test.effect.subgroup.common} \tab Label preceding common-effect subgroup-effect test. \cr
#' \code{test.effect.subgroup.random} \tab Label preceding random-effects subgroup-effect test.
#' }
#' }
#'
#' \subsection{Risk-of-bias settings (rob.*)}{
#' Settings for risk-of-bias assessments can be supplied in the named list
#' \code{rob}.
#'
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{rob} \tab Risk-of-bias assessment. \cr
#' \code{text} \tab Column heading for the risk-of-bias table. \cr
#' \code{legend} \tab Print a legend with risk-of-bias domains. \cr
#' \code{only} \tab Print risk of bias as the only right-side information.
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
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{main}, \code{heading} \tab Main title and column headings. \cr
#' \code{common} \tab Common-effect results. \cr
#' \code{random} \tab Random-effects results. \cr
#' \code{predict} \tab Prediction-interval results. \cr
#' \code{common.labels} \tab Labels for common-effect pooled results. \cr
#' \code{random.labels} \tab Labels for random-effects pooled results. \cr
#' \code{predict.labels} \tab Labels for prediction intervals. \cr
#' \code{study}, \code{study.labels} \tab Study results and study labels. \cr
#' \code{hetstat} \tab Heterogeneity statistics. \cr
#' \code{test.overall} \tab Tests for overall effects. \cr
#' \code{test.subgroup} \tab Tests for subgroup differences. \cr
#' \code{test.effect.subgroup} \tab Tests for effects in subgroups. \cr
#' \code{addline}\tab Additional lines. \cr
#' \code{axis} \tab Axis text. \cr
#' \code{smlab} \tab Summary measure. \cr
#' \code{xlab} \tab x-axis label. \cr
#' \code{lr} \tab Left / right labels. \cr
#' \code{rob} \tab Risk-of-bias text. \cr
#' \code{rob.symbols} \tab Risk-of-bias symbols. \cr
#' \code{details} \tab Statistical-method details.
#' }
#' }
#'
#' \subsection{Justifications (just.*)}{
#' Possible values are \code{"left"}, \code{"right"}, and \code{"center"}.
#'
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{just} \tab Justification for ordinary columns. \cr
#' \code{main} \tab Justification for the main title. \cr
#' \code{studlab} \tab Justification for study labels. \cr
#' \code{addcols} \tab Justification for additional columns. \cr
#' \code{addcols.left} \tab Justification for left additional columns. \cr
#' \code{addcols.right} \tab Justification for right additional columns. \cr
#' \code{label.e} \tab Justification for the experimental-group heading. \cr
#' \code{label.c} \tab Justification for the control-group heading. \cr
#' \code{rob} \tab Justification for risk-of-bias columns.
#' }
#' }
#'
#' \subsection{Fill settings (fill.*)}{
#' Colours used to fill plot regions and symbols.
#'
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{fill} \tab General fill colour. \cr
#' \code{cid} \tab Fill colour for clinically important differences. \cr
#' \code{cid.below.null} \tab CID fill below the null effect. \cr
#' \code{cid.above.null} \tab CID fill above the null effect. \cr
#' \code{equi} \tab Fill between equivalence limits. \cr
#' \code{lower.equi} \tab Fill between the lower limit and the null effect. \cr
#' \code{upper.equi} \tab Fill between the null effect and the upper limit.
#' }
#' }
#'
#' \subsection{Colours (col.*)}{
#' Colour values used for plot elements. Vectors can be supplied when colours
#' differ between subgroups or rows.
#'
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{main} \tab Main title colour. \cr
#' \code{common} \tab Common-effect colour. \cr
#' \code{random} \tab Random-effects colour. \cr
#' \code{cid} \tab Decision-line colour. \cr
#' \code{label} \tab Axis-label colour. \cr
#' \code{label.left} \tab Left-label colour. \cr
#' \code{label.right} \tab Right-label colour. \cr
#' \code{study}, \code{study.text} \tab Study text colours. \cr
#' \code{square}, \code{square.lines} \tab Square fill and border colours. \cr
#' \code{circle}, \code{circle.lines} \tab Circle fill and border colours. \cr
#' \code{inside} \tab Colour inside study symbols. \cr
#' \code{inside.common} \tab Colour inside common-effect symbols. \cr
#' \code{inside.random} \tab Colour inside random-effects symbols. \cr
#' \code{diamond} \tab General diamond fill colour. \cr
#' \code{diamond.common} \tab Common-effect diamond fill colour. \cr
#' \code{diamond.common.subgroup} \tab Common-effect subgroup diamond fill. \cr
#' \code{diamond.random} \tab Random-effects diamond fill colour. \cr
#' \code{diamond.random.subgroup} \tab Random-effects subgroup diamond fill. \cr
#' \code{diamond.lines} \tab General diamond border colour. \cr
#' \code{diamond.lines.common} \tab Common-effect diamond border colour. \cr
#' \code{diamond.lines.common.subgroup} \tab Common-effect subgroup diamond border. \cr
#' \code{diamond.lines.random} \tab Random-effects diamond border colour. \cr
#' \code{diamond.lines.random.subgroup} \tab Random-effects subgroup diamond border. \cr
#' \code{predict} \tab Prediction-interval fill colour. \cr
#' \code{predict.subgroup} \tab Subgroup prediction-interval fill. \cr
#' \code{predict.lines} \tab Prediction-interval border colour. \cr
#' \code{predict.lines.subgroup} \tab Subgroup prediction-interval border. \cr
#' \code{subgroup} \tab Subgroup text colour. \cr
#' \code{lines}, \code{header.line}, \code{jama.line} \tab Plot-line colours.
#' }
#' }
#'
#' \subsection{Printing results (print.*)}{
#' Logical arguments controlling which results are printed.
#'
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{subgroup.labels} \tab Subgroup labels. \cr
#' \code{subgroup.name} \tab Subgroup variable name. \cr
#' \code{I2} \tab I-squared. \cr
#' \code{I2.ci} \tab Confidence interval for I-squared. \cr
#' \code{tau2} \tab Tau-squared. \cr
#' \code{tau2.ci} \tab Confidence interval for tau-squared. \cr
#' \code{tau} \tab Tau. \cr
#' \code{tau.ci} \tab Confidence interval for tau. \cr
#' \code{Q} \tab Heterogeneity statistic Q. \cr
#' \code{pval.Q} \tab P-value for heterogeneity statistic Q. \cr
#' \code{Rb} \tab Rb statistic. \cr
#' \code{Rb.ci} \tab Confidence interval for Rb. \cr
#' \code{stat} \tab The z- or t-statistic for treatment effect. \cr
#' \code{Q.subgroup} \tab Q for subgroup differences. \cr
#' \code{tau2.ci.subgroup} \tab Subgroup confidence interval for tau-squared. \cr
#' \code{tau.ci.subgroup} \tab Subgroup confidence interval for tau.
#' }
#' }
#'
#' \subsection{Graphical-type (type.*)}{
#' Character strings specifying how treatment effects and confidence intervals
#' are plotted, such as \code{"square"}, \code{"diamond"}, or
#' \code{"circle"}.
#'
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{study} \tab Individual study results. \cr
#' \code{common} \tab Overall common-effect results. \cr
#' \code{random} \tab Overall random-effects results. \cr
#' \code{subgroup} \tab Subgroup results. \cr
#' \code{subgroup.common} \tab Common-effect subgroup results. \cr
#' \code{subgroup.random} \tab Random-effects subgroup results.
#' }
#' }
#'
#' \subsection{Pooled-results (pooled.*)}{
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{totals} \tab Print pooled sample sizes or totals. \cr
#' \code{events} \tab Print pooled event counts. \cr
#' \code{times} \tab Print pooled person-time totals.
#' }
#' }
#'
#' \subsection{Text settings (text.*)}{
#' Character strings used for labels and explanatory text.
#'
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{common} \tab Overall common-effect label. \cr
#' \code{random} \tab Overall random-effects label. \cr
#' \code{w.common} \tab Common-effect weight label. \cr
#' \code{w.random} \tab Random-effects weight label. \cr
#' \code{predict} \tab Prediction-interval label. \cr
#' \code{subgroup} \tab Subgroup label. \cr
#' \code{common.w} \tab Common-effect subgroup label. \cr
#' \code{random.w} \tab Random-effects subgroup label. \cr
#' \code{predict.w} \tab Subgroup prediction-interval label. \cr
#' \code{subgroup.nohet} \tab Text for subgroups without heterogeneity. \cr
#' \code{addline1} \tab First additional line. \cr
#' \code{addline2} \tab Second additional line.
#' }
#' }
#'
#' \subsection{Line-type settings (lty.*)}{
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{common} \tab Line type for common-effect estimates. \cr
#' \code{random} \tab Line type for random-effects estimates. \cr
#' \code{cid} \tab Line type for CID lines.
#' }
#' }
#'
#' \subsection{Line-width settings (lwd.*)}{
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{lwd} \tab General line width. \cr
#' \code{square} \tab Border width of squares. \cr
#' \code{diamond} \tab Border width of diamonds.
#' }
#' }
#'
#' \subsection{Digits (digits.*)}{
#' The values specify the minimum number of significant digits.
#'
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{digits} \tab Treatment effects. \cr
#' \code{se} \tab Standard errors. \cr
#' \code{stat} \tab z- or t-statistics. \cr
#' \code{pval} \tab P-values for treatment effects. \cr
#' \code{pval.Q} \tab P-values for heterogeneity tests. \cr
#' \code{Q} \tab Heterogeneity statistic Q. \cr
#' \code{tau2} \tab Between-study variance. \cr
#' \code{tau} \tab Between-study standard deviation. \cr
#' \code{I2} \tab I-squared statistic. \cr
#' \code{weight} \tab Study weights. \cr
#' \code{mean} \tab Means. \cr
#' \code{sd} \tab Standard deviations. \cr
#' \code{cor} \tab Correlations. \cr
#' \code{time} \tab Follow-up times. \cr
#' \code{n} \tab Sample sizes. \cr
#' \code{event} \tab Event counts. \cr
#' \code{TE} \tab Treatment estimates in additional data. \cr
#' \code{addcols} \tab Additional columns. \cr
#' \code{addcols.left} \tab Left additional columns. \cr
#' \code{addcols.right} \tab Right additional columns.
#' }
#' }
#'
#' \subsection{Width calculation (calcwidth.*)}{
#' Logical arguments controlling whether generated text is included when
#' calculating column widths.
#'
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{pooled} \tab Pooled results. \cr
#' \code{common} \tab Common-effect text. \cr
#' \code{random} \tab Random-effects text. \cr
#' \code{predict} \tab Prediction-interval text. \cr
#' \code{hetstat} \tab Heterogeneity statistics. \cr
#' \code{tests} \tab Overall-effect and subgroup-difference tests. \cr
#' \code{subgroup} \tab Subgroup labels. \cr
#' \code{addline} \tab Additional lines. \cr
#' \code{details} \tab Statistical-method details.
#' }
#' }
#'
#' \subsection{Column-gaps (colgap.*)}{
#' Gaps can be specified as character strings or \code{grid::unit} objects.
#'
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{colgap} \tab Default gap between columns. \cr
#' \code{left} \tab Gap between columns on the left side. \cr
#' \code{right} \tab Gap between columns on the right side. \cr
#' \code{studlab} \tab Gap after the study-label column. \cr
#' \code{forest} \tab Gap adjacent to the forest plot. \cr
#' \code{forest.left} \tab Gap between left columns and the forest plot. \cr
#' \code{forest.right} \tab Gap between the forest plot and right columns. \cr
#' \code{rob} \tab Gap between risk-of-bias columns. \cr
#' \code{rob.overall} \tab Gap before the overall risk-of-bias column.
#' }
#' }
#'
#' \subsection{Clinically important difference settings (cid.*)}{
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{cid} \tab Clinically important differences or decision thresholds. \cr
#' \code{below.null} \tab CID limits below the null effect. \cr
#' \code{above.null} \tab CID limits above the null effect. \cr
#' \code{pooled.only} \tab Restrict CID regions to pooled estimates.
#' }
#' }
#'
#' \subsection{Additional-row settings (addrow.*)}{
#' \tabular{ll}{
#' \bold{List element} \tab \bold{Description} \cr
#' \code{addrow} \tab Print an empty row above study results. \cr
#' \code{overall} \tab Print an empty row above overall results. \cr
#' \code{subgroup} \tab Print an empty row between subgroup results.
#' }
#' }
#'
#' @name forest-formatting
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
