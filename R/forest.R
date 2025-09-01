#' Forest plot to display the result of a meta-analysis
#' 
#' @description
#' Draw a forest plot (using grid graphics system) in the active
#' graphics window or store the forest plot in a file.
#' 
#' @aliases forest forest.meta
#' 
#' @param x An object of class \code{meta}.
#' @param sortvar An optional vector used to sort the individual
#'   studies (must be of same length as \code{x$TE}).
#' @param studlab A logical indicating whether study labels should be
#'   printed in the graph. A vector with study labels can also be
#'   provided (must be of same length as the vector with estimates \code{x$TE}).
#' @param layout A character string specifying the layout of the
#'   forest plot (see Details).
#' @param common A logical indicating whether common effect estimate
#'   should be plotted.
#' @param random A logical indicating whether random effects estimate
#'   should be plotted.
#' @param overall A logical indicating whether overall summaries
#'   should be plotted. This argument is useful in a meta-analysis
#'   with subgroups if summaries should only be plotted on group
#'   level.
#' @param text.common A character string used in the plot to label the
#'   pooled common effect estimate.
#' @param text.random A character string used in the plot to label the
#'   pooled random effects estimate.
#' @param lty.common Line type of pooled common effect estimate.
#' @param lty.random Line type of pooled random effects estimate.
#' @param col.common Line colour of pooled common effect estimate.
#' @param col.random Line colour of pooled random effects estimate.
#' @param text.w.common A character string used to label weights of
#'   common effect model.
#' @param text.w.random A character string used to label weights of
#'   random effects model.
#' @param prediction A logical indicating whether a prediction
#'   interval should be printed.
#' @param text.predict A character string used in the plot to label
#'   the prediction interval.
#' @param subgroup A single logical or logical vector indicating
#'   whether / which subgroup results should be shown in forest
#'   plot. This argument is useful in a meta-analysis with subgroups
#'   if summaries should not be plotted for (some) subgroups.
#' @param subgroup.hetstat A single logical or logical vector
#'   indicating whether / which information on heterogeneity in
#'   subgroups should be shown in forest plot. This argument is useful
#'   in a meta-analysis with subgroups if heterogeneity statistics
#'   should not be printed for (some) subgroups.
#' @param print.subgroup.labels A logical indicating whether subgroup
#'   label should be printed.
#' @param subgroup.name A character string with a label for the
#'   grouping variable.
#' @param print.subgroup.name A logical indicating whether the name of
#'   the grouping variable should be printed in front of the group
#'   labels.
#' @param sep.subgroup A character string defining the separator
#'   between label and levels of grouping variable.
#' @param text.common.w A character string to label the pooled common
#'   effect estimate within subgroups, or a character vector of same
#'   length as number of subgroups with corresponging labels.
#' @param text.random.w A character string to label the pooled random
#'   effect estimate within subgroups, or a character vector of same
#'   length as number of subgroups with corresponging labels.
#' @param text.predict.w A character string to label the prediction
#'   interval within subgroups, or a character vector of same length
#'   as number of subgroups with corresponging labels.
#' @param sort.subgroup A logical indicating whether groups should be
#'   ordered alphabetically.
#' @param pooled.totals A logical indicating whether total number of
#'   observations should be given in the figure.
#' @param pooled.events A logical indicating whether total number of
#'   events should be given in the figure.
#' @param pooled.times A logical indicating whether total person time
#'   at risk should be given in the figure.
#' @param study.results A logical indicating whether results for
#'   individual studies should be shown in the figure (useful to only
#'   plot subgroup results).
#' @param rob Risk of bias (RoB) assessment.
#' @param rob.text Column heading for RoB table.
#' @param rob.xpos A numeric specifying the horizontal position of the
#'   risk of bias label in RoB table heading. The value is a so called
#'   normalised parent coordinate in the horizontal direction (see
#'   \code{\link[grid]{unit}}).
#' @param rob.legend A logical specifying whether a legend with RoB
#'   domains should be printed.
#' @param rob.only A logical indicating whether the risk of bias
#'   assessment is the only information printed on the right side of
#'   the forest plot.
#' @param xlab A label for the x-axis.
#' @param xlab.pos A numeric specifying the center of the label on the
#'   x-axis.
#' @param smlab A label for the summary measure (printed at top of
#'   figure).
#' @param smlab.pos A numeric specifying the center of the label for
#'   the summary measure.
#' @param xlim The x limits (min,max) of the plot, or the character
#'   string "symmetric" to produce symmetric forest plots.
#' @param allstudies A logical indicating whether studies with
#'   inestimable treatment effects should be included in the forest
#'   plot.
#' @param weight.study A character string indicating weighting used to
#'   determine size of squares or diamonds (argument
#'   \code{type.study}) to plot individual study results. One of
#'   missing, \code{"same"}, \code{"common"}, or \code{"random"}, can
#'   be abbreviated. Plot symbols have the same size for all studies
#'   or represent study weights from common effect or random effects
#'   model.
#' @param pscale A numeric giving scaling factor for printing of
#'   single event probabilities or risk differences, i.e. if argument
#'   \code{sm} is equal to \code{"PLOGIT"}, \code{"PLN"},
#'   \code{"PRAW"}, \code{"PAS"}, \code{"PFT"}, or \code{"RD"}.
#' @param irscale A numeric defining a scaling factor for printing of
#'   single incidence rates or incidence rate differences, i.e. if
#'   argument \code{sm} is equal to \code{"IR"}, \code{"IRLN"},
#'   \code{"IRS"}, \code{"IRFT"}, or \code{"IRD"}.
#' @param irunit A character specifying the time unit used to
#'   calculate rates, e.g., person-years.
#' @param file File name.
#' @param width Width of graphics file.
#' @param rows.gr Additional rows in forest plot to change height of
#'   graphics file (e.g., in order to add a title at the top of the
#'   forest plot).
#' @param func.gr Name of graphics function, e.g., \code{\link{pdf}}.
#' @param args.gr List with additional graphical parameters passed on
#'   to graphics function (argument 'height' cannot be provided as the
#'   height is calculated internally; use instead argument 'rows.gr').
#' @param dev.off A logical to specify whether current graphics device
#'   should be shut down, i.e., whether file should be stored.
#' @param ref A numerical giving the reference value to be plotted as
#'   a line in the forest plot. No reference line is plotted if
#'   argument \code{ref} is equal to \code{NA}.
#' @param cid A numeric value or vector specifying clinically important
#'   differences (CID) / decision thresholds used to calculate probabilities
#'   of clinically important benefit or harm, or not important effects
#'   (see Details).
#' @param cid.below.null A numeric value or vector specifying CID limits below
#'   the null effect (see Details).
#' @param cid.above.null A numeric value or vector specifying CID limits above
#'   the null effect (see Details).
#' @param lty.cid Line type for CID lines.
#' @param col.cid Line colour for CID lines.
#' @param fill.cid Colour(s) for regions below or above CID limits.
#' @param fill.cid.below.null Colour of CID regions below null effect /
#'   reference value. Can be equal to the number of lower limits or
#'   the number of limits plus 1 (in this case the region between
#'   minimum and smallest limit is also filled).
#' @param fill.cid.above.null Colour of CID regions above null effect /
#'   reference value. Can be equal to the number of upper limits or the
#'   number of limits plus 1 (in this case the region between largest
#'   limit and maximum is also filled).
#' @param cid.pooled.only A logical indicating whether CID regions should only
#'   be visible for pooled estimates or also individual studies.
#' @param fill Colour for background of confidence interval plot (also used
#'   as colour for region between CID limits if argument \code{fill.equi} was
#'   not provided).
#' @param fill.equi Colour(s) for region between limits of equivalence defined
#'   by arguments \code{cid}, \code{cid.lower} or \code{cid.upper}.
#' @param fill.lower.equi Colour of region between lower limit(s) and
#'   reference value. Can be equal to the number of lower limits or
#'   the number of limits plus 1 (in this case the the region between
#'   minimum and smallest limit is also filled).
#' @param fill.upper.equi Colour of region between reference value and
#'   upper limit(s). Can be equal to the number of upper limits or the
#'   number of limits plus 1 (in this case the region between largest
#'   limit and maximum is also filled).
#' @param leftcols A character vector specifying (additional) columns
#'   to be printed on the left side of the forest plot or a logical
#'   value (see Details).
#' @param rightcols A character vector specifying (additional) columns
#'   to be printed on the right side of the forest plot or a logical
#'   value (see Details).
#' @param leftlabs A character vector specifying labels for
#'   (additional) columns on left side of the forest plot (see
#'   Details).
#' @param rightlabs A character vector specifying labels for
#'   (additional) columns on right side of the forest plot (see
#'   Details).
#' @param label.e Label to be used for experimental group in table
#'   heading.
#' @param label.c Label to be used for control group in table heading.
#' @param label.e.attach A character string or vector specifying the column
#'   name(s) where label \code{label.e} should be attached to in table heading.
#' @param label.c.attach A character string or vector specifying the column
#'   name(s) where label \code{label.c} should be attached to in table heading.
#' @param label.left Graph label on left side of null effect.
#' @param label.right Graph label on right side of null effect.
#' @param bottom.lr A logical indicating whether labels on right and
#'   left side should be printed at bottom or top of forest plot.
#' @param lab.NA A character string to label missing values.
#' @param lab.NA.effect A character string to label missing values in
#'   individual treatment estimates and confidence intervals.
#' @param lab.NA.weight A character string to label missing weights.
#' @param lwd The line width, see \code{\link{par}}.
#' @param at The points at which tick-marks are to be drawn, see
#'   \code{grid.xaxis}.
#' @param label A logical value indicating whether to draw the labels
#'   on the tick marks, or an expression or character vector which
#'   specify the labels to use. See \code{\link{grid.xaxis}}.
#' @param col.label The colour of labels on the x-axis.
#' @param type.study A character string or vector specifying how to
#'   plot treatment effects and confidence intervals for individual
#'   studies (see Details).
#' @param type.common A character string specifying how to plot
#'   treatment effect and confidence interval for common effect
#'   meta-analysis (see Details).
#' @param type.random A character string specifying how to plot
#'   treatment effect and confidence interval for random effects
#'   meta-analysis (see Details).
#' @param type.subgroup A character string specifying how to plot
#'   treatment effect and confidence interval for subgroup results
#'   (see Details).
#' @param type.subgroup.common A character string specifying how to
#'   plot treatment effect and confidence interval for subgroup
#'   results (common effect model).
#' @param type.subgroup.random A character string specifying how to
#'   plot treatment effect and confidence interval for subgroup
#'   results (random effects model).
#' @param col.study The colour for individual study results and
#'   confidence limits.
#' @param col.inside The colour for individual study results and
#'   confidence limits if confidence limits are completely within
#'   squares.
#' @param col.inside.common The colour for result of common effect
#'   meta-analysis if confidence limit lies completely within square.
#' @param col.inside.random The colour for result of random effects
#'   meta-analysis if confidence limit lies completely within square.
#' @param col.square The colour for squares reflecting study's weight
#'   in the meta-analysis.
#' @param col.square.lines The colour for the outer lines of squares
#'   reflecting study's weight in the meta-analysis.
#' @param col.circle The colour for circles reflecting study weights
#'   in the meta-analysis.
#' @param col.circle.lines The colour for the outer lines of circles
#'   reflecting study's weight in the meta-analysis.
#' @param col.diamond The colour of diamonds representing the results
#'   for common effect and random effects models.
#' @param col.diamond.common The colour(s) of diamonds for common effect
#'   estimates.
#' @param col.diamond.random The colour(s) of diamonds for random effects
#'   estimates.
#' @param col.diamond.lines The colour of the outer lines of diamonds
#'   representing the results for common effect and random effects
#'   models.
#' @param col.diamond.lines.common The colour(s) of the outer lines of
#'   diamond for common effect estimate.
#' @param col.diamond.lines.random The colour(s) of the outer lines of
#'   diamond for random effects estimate.
#' @param col.predict Background colour(s) of prediction intervals.
#' @param col.predict.lines Colour(s) of outer lines of prediction
#'   intervals.
#' @param col.subgroup The colour to print information on subgroups.
#' @param col.label.left The colour of the label on the left side of the null
#'   effect.
#' @param col.label.right The colour of the label on the right side of the null
#'   effect.
#' @param hetstat Either a logical value indicating whether to print
#'   results for heterogeneity measures at all or a character string
#'   (see Details).
#' @param overall.hetstat A logical value indicating whether to print
#'   heterogeneity measures for overall treatment comparisons. This
#'   argument is useful in a meta-analysis with subgroups if
#'   heterogeneity statistics should only be printed on subgroup
#'   level.
#' @param hetlab Label printed in front of results for heterogeneity
#'   measures.
#' @param resid.hetstat A logical value indicating whether to print
#'   measures of residual heterogeneity in a meta-analysis with
#'   subgroups.
#' @param resid.hetlab Label printed in front of results for residual
#'   heterogeneity measures.
#' @param print.I2 A logical value indicating whether to print the
#'   value of the I-squared statistic.
#' @param print.I2.ci A logical value indicating whether to print the
#'   confidence interval of the I-squared statistic.
#' @param print.tau2 A logical value indicating whether to print the
#'   value of the between-study variance \eqn{\tau^2}.
#' @param print.tau2.ci A logical value indicating whether to print
#'   the confidence interval of \eqn{\tau^2}.
#' @param print.tau A logical value indicating whether to print
#'   \eqn{\tau}, the square root of the between-study variance
#'   \eqn{\tau^2}.
#' @param print.tau.ci A logical value indicating whether to print the
#'   confidence interval of \eqn{\tau}.
#' @param print.Q A logical value indicating whether to print the
#'   value of the heterogeneity statistic Q.
#' @param print.pval.Q A logical value indicating whether to print the
#'   p-value of the heterogeneity statistic Q.
#' @param print.Rb A logical value indicating whether to print the
#'   value of the I-squared statistic.
#' @param print.Rb.ci A logical value indicating whether to print the
#'   confidence interval of the I-squared statistic.
#' @param text.subgroup.nohet A logical value or character string
#'   which is printed to indicate subgroups with less than two studies
#'   contributing to meta-analysis (and thus without
#'   heterogeneity). If FALSE, heterogeneity statistics are printed
#'   (with NAs).
#' @param LRT A logical value indicating whether to report
#'   Likelihood-Ratio or Wald-type test of heterogeneity for
#'   generalised linear mixed models.
#' @param test.overall A logical value indicating whether to print
#'   results of test for overall effect.
#' @param test.overall.common A logical value indicating whether to
#'   print results of test for overall effect (common effect model).
#' @param test.overall.random A logical value indicating whether to
#'   print results of test for overall effect (random effects model).
#' @param label.test.overall.common Label printed in front of results
#'   of test for overall effect (common effect model).
#' @param label.test.overall.random Label printed in front of results
#'   of test for overall effect (random effects model).
#' @param print.stat A logical value indicating whether z- or t-value
#'   for test of treatment effect should be printed.
#' @param test.subgroup A logical value indicating whether to print
#'   results of test for subgroup differences.
#' @param test.subgroup.common A logical value indicating whether to
#'   print results of test for subgroup differences (common effect
#'   model).
#' @param test.subgroup.random A logical value indicating whether to
#'   print results of test for subgroup differences (random effects
#'   model).
#' @param common.subgroup A single logical or logical vector
#'   indicating whether / which common effect estimates should be
#'   printed for subgroups.
#' @param random.subgroup A single logical or logical vector
#'   indicating whether / which random effects estimates should be
#'   printed for subgroups.
#' @param prediction.subgroup A single logical or logical vector
#'   indicating whether / which prediction intervals should be printed
#'   for subgroups.
#' @param print.Q.subgroup A logical value indicating whether to print
#'   the value of the heterogeneity statistic Q (test for subgroup
#'   differences).
#' @param label.test.subgroup.common Label printed in front of results
#'   of test for subgroup differences (common effect model).
#' @param label.test.subgroup.random Label printed in front of results
#'   of test for subgroup differences (random effects model).
#' @param test.effect.subgroup A single logical or logical vector
#'   indicating whether / which tests for effect in subgroups should
#'   be printed.
#' @param test.effect.subgroup.common A single logical or logical
#'   vector indicating whether / which tests for effect in subgroups
#'   should be printed (common effect model).
#' @param test.effect.subgroup.random A single logical or logical
#'   vector indicating whether / which tests for effect in subgroups
#'   should be printed (random effects model).
#' @param label.test.effect.subgroup.common Label printed in front of
#'   results of test for effect in subgroups (common effect model).
#' @param label.test.effect.subgroup.random Label printed in front of
#'   results of test for effect in subgroups (random effects model).
#' @param text.addline1 Text for first additional line (below
#'   meta-analysis results).
#' @param text.addline2 Text for second additional line (below
#'   meta-analysis results).
#' @param details A logical specifying whether details on statistical
#'   methods should be printed.
#' @param col.lines The colour of lines.
#' @param header.line A logical value indicating whether to print a
#'   header line or a character string ("both", "below", "").
#' @param col.header.line Colour of the header line(s).
#' @param col.jama.line Colour of the additional JAMA lines.
#' @param data.pooled Data set with information for line(s) with pooled
#'   results (see Details).
#' @param fontsize The size of text (in points), see
#'   \code{\link{gpar}}.
#' @param fontfamily The font family, see \code{\link{gpar}}.
#' @param fs.heading The size of text for column headings, see
#'   \code{\link{gpar}}.
#' @param fs.common The size of text for results of common effect
#'   model, see \code{\link{gpar}}.
#' @param fs.random The size of text for results of random effects
#'   model, see \code{\link{gpar}}.
#' @param fs.predict The size of text for results of prediction
#'   interval, see \code{\link{gpar}}.
#' @param fs.common.labels The size of text for label of common effect
#'   model, see \code{\link{gpar}}.
#' @param fs.random.labels The size of text for label of random
#'   effects model, see \code{\link{gpar}}.
#' @param fs.predict.labels The size of text for label of prediction
#'   interval, see \code{\link{gpar}}.
#' @param fs.study The size of text for results of individual studies,
#'   see \code{\link{gpar}}.
#' @param fs.study.labels The size of text for labels of individual
#'   studies, see \code{\link{gpar}}.
#' @param fs.hetstat The size of text for heterogeneity measures, see
#'   \code{\link{gpar}}.
#' @param fs.test.overall The size of text of test for overall effect,
#'   see \code{\link{gpar}}.
#' @param fs.test.subgroup The size of text of test of subgroup
#'   differences, see \code{\link{gpar}}.
#' @param fs.test.effect.subgroup The size of text of test of effect
#'   in subgroups, see \code{\link{gpar}}.
#' @param fs.addline The size of text for additional lines, see
#'   \code{\link{gpar}}.
#' @param fs.axis The size of text on x-axis, see \code{\link{gpar}}.
#' @param fs.smlab The size of text of label for summary measure, see
#'   \code{\link{gpar}}.
#' @param fs.xlab The size of text of label on x-axis, see
#'   \code{\link{gpar}}.
#' @param fs.lr The size of text of label on left and right side of
#'   forest plot, see \code{\link{gpar}}.
#' @param fs.rob The size of text of risk of bias items in the legend,
#'   see \code{\link{gpar}}.
#' @param fs.rob.symbols The size of risk of bias symbols, see
#'   \code{\link{gpar}}.
#' @param fs.details The size of text for details on (meta-analysis)
#'   methods, see \code{\link{gpar}}.
#' @param ff.heading The fontface for column headings, see
#'   \code{\link{gpar}}.
#' @param ff.common The fontface of text for results of common effect
#'   model, see \code{\link{gpar}}.
#' @param ff.random The fontface of text for results of random effects
#'   model, see \code{\link{gpar}}.
#' @param ff.predict The fontface of text for results of prediction
#'   interval, see \code{\link{gpar}}.
#' @param ff.common.labels The fontface of text for label of common
#'   effect model, see \code{\link{gpar}}.
#' @param ff.random.labels The fontface of text for label of random
#'   effects model, see \code{\link{gpar}}.
#' @param ff.predict.labels The fontface of text for label of
#'   prediction interval, see \code{\link{gpar}}.
#' @param ff.study The fontface of text for results of individual
#'   studies, see \code{\link{gpar}}.
#' @param ff.study.labels The fontface of text for labels of
#'   individual studies, see \code{\link{gpar}}.
#' @param ff.hetstat The fontface of text for heterogeneity measures,
#'   see \code{\link{gpar}}.
#' @param ff.test.overall The fontface of text of test for overall
#'   effect, see \code{\link{gpar}}.
#' @param ff.test.subgroup The fontface of text for test of subgroup
#'   differences, see \code{\link{gpar}}.
#' @param ff.test.effect.subgroup The fontface of text for test of
#'   effect in subgroups, see \code{\link{gpar}}.
#' @param ff.addline The fontface of text for additional lines, see
#'   \code{\link{gpar}}.
#' @param ff.axis The fontface of text on x-axis, see
#'   \code{\link{gpar}}.
#' @param ff.smlab The fontface of text of label for summary measure,
#'   see \code{\link{gpar}}.
#' @param ff.xlab The fontface of text of label on x-axis, see
#'   \code{\link{gpar}}.
#' @param ff.lr The fontface of text of label on left and right side
#'   of forest plot, see \code{\link{gpar}}.
#' @param ff.rob The fontface of text of risk of bias items, see
#'   \code{\link{gpar}}.
#' @param ff.rob.symbols The fontface of risk of bias symbols, see
#'   \code{\link{gpar}}.
#' @param ff.details The fontface for details on (meta-analysis)
#'   methods, see \code{\link{gpar}}.
#' @param squaresize A numeric used to increase or decrease the size
#'   of squares in the forest plot.
#' @param lwd.square The line width of the border around squares.
#' @param lwd.diamond The line width of the border around diamonds.
#' @param arrow.type A character string indicating whether arrows
#'   printed for results outside the forest plot should be
#'   \code{"open"}, or \code{"closed"}, can be abbreviated.
#' @param arrow.length The length of arrows in inches.
#' @param plotwidth Either a character string, e.g., "8cm", "60mm", or
#'   "3inch", or a \code{\link[grid]{unit}} object specifying width of
#'   the forest plot.
#' @param colgap Either a character string or a
#'   \code{\link[grid]{unit}} object specifying gap between columns
#'   printed on left and right side of forest plot.
#' @param colgap.left Either a character string or a
#'   \code{\link[grid]{unit}} object specifying gap between columns
#'   printed on left side of forest plot.
#' @param colgap.right Either a character string or a
#'   \code{\link[grid]{unit}} object specifying gap between columns
#'   printed on right side of forest plot.
#' @param colgap.studlab Either a character string or a
#'   \code{\link[grid]{unit}} object specifying gap between column
#'   with study labels and subsequent column.
#' @param colgap.forest Either a character string or a
#'   \code{\link[grid]{unit}} object specifying gap between column
#'   adjacent to forest plot and the forest plot.
#' @param colgap.forest.left Either a character string or a
#'   \code{\link[grid]{unit}} object specifying gap between column on
#'   the left side of forest plot and the forest plot.
#' @param colgap.forest.right Either a character string or a
#'   \code{\link[grid]{unit}} object specifying gap between column on
#'   the right side of forest plot and the forest plot.
#' @param colgap.rob Either a character string or a
#'   \code{\link[grid]{unit}} object specifying gap between risk of
#'   bias columns.
#' @param colgap.rob.overall Either a character string or a
#'   \code{\link[grid]{unit}} object specifying gap before column with
#'   overall risk of bias assessment.
#' @param calcwidth.pooled A logical indicating whether text for
#'   common effect and random effects model should be considered to
#'   calculate width of the column with study labels.
#' @param calcwidth.common A logical indicating whether text given in
#'   arguments \code{text.common} and \code{text.common.w} should be
#'   considered to calculate width of the column with study labels.
#' @param calcwidth.random A logical indicating whether text given in
#'   arguments \code{text.random} and \code{text.random.w} should be
#'   considered to calculate width of the column with study labels.
#' @param calcwidth.predict A logical indicating whether text given in
#'   argument \code{text.predict} and \code{text.predict.w} should be
#'   considered to calculate width of the column with study labels.
#' @param calcwidth.hetstat A logical indicating whether text for
#'   heterogeneity statistics should be considered to calculate width
#'   of the column with study labels.
#' @param calcwidth.tests A logical indicating whether text for tests
#'   of overall effect or subgroup differences should be considered to
#'   calculate width of the column with study labels.
#' @param calcwidth.subgroup A logical indicating whether text with
#'   subgroup labels should be considered to calculate width of the
#'   column with study labels.
#' @param calcwidth.addline A logical indicating whether text for
#'   additional lines should be considered to calculate width of the
#'   column with study labels.
#' @param just Justification of text in all columns but columns with
#'   study labels and additional variables (possible values: "left",
#'   "right", "center").
#' @param just.studlab Justification of text for study labels
#'   (possible values: "left", "right", "center").
#' @param just.addcols Justification of text for additional columns
#'   (possible values: "left", "right", "center").
#' @param just.addcols.left Justification of text for additional
#'   columns on left side of forest plot (possible values: "left",
#'   "right", "center"). Can be of same length as number of additional
#'   columns on left side of forest plot.
#' @param just.addcols.right Justification of text for additional
#'   columns on right side of forest plot (possible values: "left",
#'   "right", "center"). Can be of same length as number of additional
#'   columns on right side of forest plot.
#' @param just.label.e Justification of text for experimental group in table
#'   heading (possible values: "left", "right", "center").
#' @param just.label.c Justification of text for control group in table
#'   heading (possible values: "left", "right", "center").
#' @param bmj.text A character string used in the plot with BMJ layout
#'   to label the group specific information.
#' @param bmj.xpos A numeric specifying the horizontal position of the
#'   BMJ label. The value is a so called normalised parent coordinate
#'   in the horizontal direction (see \code{\link[grid]{unit}}).
#' @param bmj.sep A character string used to separate sample sizes
#'   from number of events or means / standard deviations.
#' @param spacing A numeric determining line spacing in a forest plot.
#' @param addrow A logical value indicating whether an empty row is
#'   printed above study results.
#' @param addrow.overall A logical value indicating whether an empty
#'   row is printed above overall meta-analysis results.
#' @param addrow.subgroups A logical value indicating whether an empty
#'   row is printed between results for subgroups.
#' @param addrows.below.overall A numeric value indicating how many
#'   empty rows are printed between meta-analysis results and
#'   heterogeneity statistics and test results.
#' @param new A logical value indicating whether a new figure should
#'   be printed in an existing graphics window.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in forest plots. If \code{backtransf = TRUE},
#'   results for \code{sm = "OR"} are presented as odds ratios rather
#'   than log odds ratios and results for \code{sm = "ZCOR"} are
#'   presented as correlations rather than Fisher's z transformed
#'   correlations, for example.
#' @param digits Minimal number of significant digits for treatment
#'   effects, see \code{print.default}.
#' @param digits.se Minimal number of significant digits for standard
#'   errors.
#' @param digits.stat Minimal number of significant digits for z- or
#'   t-statistic for test of overall effect.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance.
#' @param digits.tau Minimal number of significant digits for square
#'   root of between-study variance.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall treatment effect.
#' @param digits.pval.Q Minimal number of significant digits for
#'   p-value of heterogeneity test.
#' @param digits.Q Minimal number of significant digits for
#'   heterogeneity statistic Q.
#' @param digits.I2 Minimal number of significant digits for I-squared
#'   statistic.
#' @param digits.weight Minimal number of significant digits for
#'   weights.
#' @param digits.mean Minimal number of significant digits for means;
#'   only applies to \code{\link{metacont}} objects.
#' @param digits.sd Minimal number of significant digits for standard
#'   deviations; only applies to \code{\link{metacont}} objects.
#' @param digits.cor Minimal number of significant digits for
#'   correlations; only applies to \code{\link{metacor}} objects.
#' @param digits.time Minimal number of significant digits for times;
#'   only applies to \code{\link{metainc}} and \code{\link{metarate}}
#'   objects.
#' @param digits.n Minimal number of significant digits for sample
#'   sizes.
#' @param digits.event Minimal number of significant digits for event
#'   numbers.
#' @param digits.TE Minimal number of significant digits for list
#'   element 'TE'.
#' @param digits.addcols A vector or scalar with minimal number of
#'   significant digits for additional columns.
#' @param digits.addcols.left A vector or scalar with minimal number
#'   of significant digits for additional columns on left side of
#'   forest plot.
#' @param digits.addcols.right A vector or scalar with minimal number
#'   of significant digits for additional columns on right side of
#'   forest plot.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param big.mark A character used as thousands separator.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   overall effect should be printed according to JAMA reporting
#'   standards.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional graphical arguments.
#' 
#' @details
#' A forest plot, also called confidence interval plot, is drawn in
#' the active graphics window. The forest functions in R package
#' \bold{meta} are based on the grid graphics system. Resize the
#' graphics windows if the forest plot is too large or too small for
#' the graphics window. Alternatively, save the forest plot in a file.
#' 
#' \subsection{Saving forest plots}{
#' 
#' A forest plot can be directly stored in a file using argument
#' \code{file} or specifying the R function for the graphics device
#' driver using argument \code{func.gr}, e.g., \code{\link{pdf}}. If
#' only the filename is provided, the extension is checked and matched
#' against the most common graphics device drivers.
#' 
#' \tabular{ll}{
#' \bold{Extension} \tab \bold{Graphics device} \cr
#' \code{.pdf} \tab \code{\link{pdf}} \cr
#' \code{.ps} \tab \code{\link{postscript}} \cr
#' \code{.svg} \tab \code{\link{svg}} \cr
#' \code{.bmp} \tab \code{\link{bmp}} \cr
#' \code{.jpg} / \code{.jpeg} \tab \code{\link{jpeg}} \cr
#' \code{.png} \tab \code{\link{png}} \cr
#' \code{.tif} / \code{.tiff} \tab \code{\link{tiff}}
#' }
#'
#' The height of the graphics device is automatically determined if
#' the forest plot is saved to a file. Argument \code{rows.gr} can be
#' used to increase or decrease the number of rows shown in the forest
#' plot (either to show missing information or to remove
#' whitespace). The width of the graphics device can be specified with
#' argument \code{width}, see, for example, \code{\link{pdf}} or
#' \code{\link{jpeg}}. Other arguments of graphics device functions
#' can be provided as a list in argument \code{args.gr}.
#'
#' Alternatively, the (resized) graphics window can be stored to a
#' file using either \code{\link{dev.copy2eps}} or
#' \code{\link{dev.copy2pdf}}. It is also possible to manually create
#' a file using, for example, \code{\link{pdf}}, \code{\link{png}}, or
#' \code{\link{svg}} and to specify the width and height of the
#' graphic (see Examples).
#' }
#' 
#' \subsection{Default layout for studies and pooled effects}{
#' 
#' By default, treatment estimates and confidence intervals are
#' plotted in the following way:
#' \itemize{
#' \item For an individual study, a square with treatment estimate in
#'   the center and confidence interval as line extending either side
#'   of the square (\code{type.study = "square"})
#' \item For meta-analysis results, a diamond with treatment estimate
#'   in the center and right and left side corresponding to lower and
#'   upper confidence limits (\code{type.common = "diamond"},
#'   \code{type.random = "diamond"}, and \code{type.subgroup = "diamond"})
#' }
#' 
#' In a forest plot, size of the squares typically reflects the precision of
#' individual treatment estimates based either on the common effect
#' (\code{weight.study = "common"}) or random effects meta-analysis
#' (\code{weight.study = "random"}). Information from meta-analysis object
#' \code{x} is utilised if argument \code{weight.study} is missing. Weights
#' from the common effect model are used if argument \code{x$common} is
#' \code{TRUE}; weights from the random effects model are used if argument
#' \code{x$random} is \code{TRUE} and \code{x$common} is \code{FALSE}.
#' The same square sizes are used if \code{weight.study = "same"}.
#' 
#' A prediction interval for treatment effect of a new study (Higgins
#' et al., 2009) is given in the forest plot if arguments
#' \code{prediction} and \code{random} are \code{TRUE}. For
#' graphical presentation of prediction intervals the approach by
#' Guddat et al. (2012) is used.
#' }
#' 
#' \subsection{Columns printed on left side of forest plot}{
#' 
#' Argument \code{leftcols} can be used to specify columns which are
#' printed on the left side of the forest plot. By default, i.e. if
#' argument \code{leftcols} is \code{NULL} and \code{layout = "meta"},
#' and depending on the class of the meta-analysis object (which is
#' defined by the R function used to generate the object) a different
#' set of \emph{\bold{columns}} is printed \emph{\bold{on the left
#' side of the forest plot}}:
#' \tabular{ll}{
#' \bold{Function} \tab \bold{Value of argument leftcols} \cr
#' \code{\link{metabin}} \tab \code{c("studlab", "event.e", "n.e",
#'   "event.c", "n.c")} \cr
#' \code{\link{metacont}} \tab \code{c("studlab", "n.e", "mean.e",
#'   "sd.e", "n.c", "mean.c", "sd.c")} \cr
#' \code{\link{metacor}} \tab \code{c("studlab", "n")} \cr
#' \code{\link{metagen}} \tab \code{c("studlab", "TE", "seTE")} \cr
#' \code{\link{metainc}} \tab \code{c("studlab", "event.e", "time.e",
#'   "event.c", "time.c")} \cr
#' \code{\link{metamean}} \tab \code{c("studlab", "n", "mean", "sd")}
#'   \cr
#' \code{\link{metaprop}} \tab \code{c("studlab", "event", "n")} \cr
#' \code{\link{metarate}} \tab \code{c("studlab", "event", "time", "n")}
#' }
#'
#' For three-level models, the cluster variable is printed next to the
#' study labels (value \code{"cluster"} in argument \code{leftcols}).
#'
#' By default, study labels and labels for pooled estimates and
#' heterogeneity statistics will be printed in the first column on the
#' left side of the forest plot. The character string \code{"studlab"}
#' is used to identify study labels as this is the name of the list
#' element of a meta-analysis object.
#'
#' If the character string \code{"studlab"} is not provided in
#' \code{leftcols} and \code{rightcols}, the first \emph{additional}
#' variable specified by the user is used as study labels (and labels
#' for pooled estimates are printed in this column). Additional
#' variables are any variables not mentioned in the section on
#' predefined column names below. For example, \code{leftcols =
#' "studlab"} and \code{leftcols = "study"} would result in the same
#' forest plot if the variable \code{"study"} was used in the command
#' to conduct the meta-analysis. If no additional variable is provided
#' by the user, no study labels will be printed.
#' }
#' 
#' \subsection{Overlapping information on left side of forest plot}{
#'
#' Depending on the number of columns printed on the left side of the
#' forest plot, information on heterogeneity measures or statistical
#' tests (see below) can be overlapping with the x-axis. Argument
#' \code{addrows.below.overall} can be used to specify the number of
#' empty rows that are printed between meta-analysis results and
#' information on heterogeneity measures and statistical tests. By
#' default, no additional rows are added to the forest plot. If
#' \code{addrows.below.overall = NULL}, the function tries to add a
#' sufficient number of empty rows to prevent overlapping
#' text. Another possibility is to manually increase the space between
#' the columns on the left side (argument \code{colgap.left}) or
#' between the columns on the left side and the forest plot (argument
#' \code{colgap.forest.left}).
#' }
#' 
#' \subsection{Columns printed on right side of forest plot}{
#' 
#' Argument \code{rightcols} can be used to
#' specify columns which are printed on the right side of the
#' forest plot. If argument \code{rightcols} is
#' \code{FALSE}, no columns will be printed on the right side. By
#' default, i.e. if argument \code{rightcols} is
#' \code{NULL} and \code{layout = "meta"}, the following
#' \emph{\bold{columns}} will be printed \emph{\bold{on the right side
#' of the forest plot}}:
#' \tabular{ll}{
#' \bold{Meta-analysis results} \tab \bold{Value of argument
#'   rightcols} \cr
#' No summary \tab \code{c("effect", "ci")} \cr
#' Only common effect model \tab \code{c("effect", "ci", "w.common")}
#'   \cr
#' Only random effects model \tab \code{c("effect", "ci", "w.random")}
#'   \cr
#' Both models \tab \code{c("effect", "ci", "w.common", "w.random")}
#' }
#'
#' By default, estimated treatment effect and corresponding confidence
#' interval will be printed. Depending on arguments \code{common} and
#' \code{random}, weights of the common effect and/or random effects
#' model will be given too.
#' }
#'
#' \subsection{Predefined columns and column labels}{
#' 
#' The arguments \code{leftlabs} and \code{rightlabs} can be used to
#' specify column headings which are printed on the left or right side of
#' the forest plot. For certain columns predefined labels exist which
#' are used by default, i.e., if arguments \code{leftlabs} and
#' \code{rightlabs} are \code{NULL}:
#' \tabular{rcccccc}{
#' Column: \tab \code{studlab} \tab \code{TE} \tab \code{seTE} \tab
#'   \code{cluster} \tab \code{n.e} \tab \code{n.c} \cr 
#' Label: \tab "Study" \tab "TE" \tab "seTE" \tab "Cluster" \tab
#'   "Total" \tab "Total" \cr
#' \cr
#' Column: \tab \code{n} \tab \code{event.e} \tab \code{event.c} \tab
#'   \code{event} \tab \code{mean.e} \tab \code{mean.c} \cr
#' Label: \tab "Total" \tab "Events" \tab "Events" \tab "Events" \tab
#'   "Mean" \tab "Mean" \cr
#' \cr
#' Column: \tab \code{sd.e} \tab \code{sd.c} \tab \code{time.e}
#'   \tab \code{time.c} \tab \code{effect} \tab \cr
#' Label: \tab "SD" \tab "SD" \tab "Time" \tab "Time" \tab
#'   \code{x$sm} \tab \cr
#' \cr
#' Column: \tab \code{ci} \tab \code{effect.ci} \tab
#'   \code{w.common} \tab \code{w.random} \tab \code{cycles} \tab \cr
#' Label: \tab \code{x$level}"\%-CI" \tab \emph{effect+ci} \tab
#'   "W(common)" \tab "W(random)" \tab "Cycles" \tab \cr
#' \cr
#' Column: \tab \code{pval} \tab \code{tau2} \tab
#'   \code{tau} \tab \tab \tab \cr
#' Label: \tab "P-value" \tab "Tau2" \tab "Tau" \tab \tab \tab
#' }
#'
#' For other columns, the column name will be used as a label if no
#' column label is defined. It is possible to only provide labels for
#' new columns (see Examples). Otherwise the length of \code{leftlabs}
#' and \code{rightlabs} must be the same as the number of printed
#' columns. The value \code{NA} can be used to specify columns which
#' should use default labels (see Examples).
#'
#' In pairwise meta-analysis comparing two groups (i.e.,
#' \code{\link{metabin}}, \code{\link{metacont}},
#' \code{\link{metainc}}, and \code{\link{metagen}} depending on the
#' outcome), arguments \code{label.e} and \code{label.c} are used to
#' label columns belonging to the two treatment groups. By default,
#' labels defined in the meta-analysis object are used. The columns
#' where treatment labels are attached can be changed using arguments
#' \code{label.e.attach} and \code{label.c.attach}.
#' }
#'
#' \subsection{Risk of bias assessment}{
#' 
#' A risk of bias (RoB) assessment can be shown in the forest plot by
#' either using a meta-analysis object with an RoB assessment as main
#' input or providing a suitable object created with
#' \code{\link{rob}}. Argument \code{rob = FALSE} can be used to
#' suppress the print of the risk of bias information.
#'
#' RoB assessments are shown as the only information on the right side
#' of the forest plot. Thus, arguments \code{rightcols} and
#' \code{rightlabs} should not be used. Predefined columns shown by
#' default on the right side of a forest plot will be moved to the
#' left side.
#' }
#' 
#' \subsection{Information on heterogeneity and statistical tests}{
#' 
#' Argument \code{hetstat} can be a character string to specify where
#' to print heterogeneity information:
#' \itemize{
#' \item row with results for common effect model (\code{hetstat =
#' "common"}),
#' \item row with results for random effects model (\code{hetstat =
#' "random"}).
#' }
#' 
#' Otherwise, information on heterogeneity measures is printed below
#' the meta-analysis results if argument \code{overall.hetstat = TRUE}
#' (default). The heterogeneity measures to print can be specified
#' (see list of arguments following \code{overall.hetstat}).
#' 
#' In addition, the following arguments can be used to print results
#' for various statistical tests:
#' \tabular{ll}{
#' \bold{Argument} \tab \bold{Statistical test} \cr
#' \code{test.overall.common} \tab Test for overall effect (common
#'   effect model) \cr
#' \code{test.overall.random} \tab Test for overall effect (random
#'   effects model) \cr
#' \code{test.effect.subgroup.common} \tab Test for effect in subgroup
#'   (CE model) \cr
#' \code{test.effect.subgroup.random} \tab Test for effect in subgroup
#'   (RE model) \cr
#' \code{test.subgroup.common} \tab Test for subgroup differences (CE
#'   model) \cr
#' \code{test.subgroup.random} \tab Test for subgroup differences (RE
#'   model)
#' }
#' 
#' By default, these arguments are \code{FALSE} with exception of
#' tests for subgroup differences which are \code{TRUE}. R function
#' \code{\link{settings.meta}} can be used to change this default for
#' the entire R session. For example, use the following command to
#' always print results of tests for an overall effect:
#' \code{settings.meta(test.overall = TRUE)}.
#' }
#' 
#' \subsection{Highlight regions corresponding to minimal clinically important
#' differences}{
#'
#' Regions corresponding to minimal clinically important differences can be
#' added to the forest plot using either argument \code{cid} or
#' \code{cid.below.null} and \code{cid.above.null}. Input for the later arguments will
#' be ignored if argument \code{cid} was specified. In this case, the values
#' of \code{cid.below.null} and \code{cid.above.null} will be equal to
#' \itemize{
#' \item \code{cid} and \code{1 / cid} for ratio measures,
#' \item \code{cid} and \code{-cid} for difference measures.
#' }
#' 
#' Thresholds based on argument \code{cid} will always be symmetric. Asymmetric
#' thresholds can be defined using arguments \code{cid.below.null} and
#' \code{cid.above.null}.
#' }
#'
#' \subsection{Flexible printing of subgroup results}{
#'
#' Argument \code{subgroup} determines whether summary results are
#' printed for subgroups. A logical vector of length equal to the
#' number of subgroups can be provided to determine which subgroup
#' summaries are printed. By default, only subgroup results based on
#' at least two studies are printed which is identical to use argument
#' \code{subgroup = k.w > 1}. The order of the logical vector
#' corresponds to the order of subgroups in list element 'subgroup.levels' of a
#' meta-analysis object. Argument \code{subgroup = k.w >= 1} can be
#' used to show results for all subgroups (including those with a
#' single study).
#'
#' The following arguments can be used in a similar way:
#'
#' \itemize{
#' \item \code{subgroup.hetstat} (heterogeneity statistic in
#'   subgroups),
#' \item \code{common.subgroup} (common effect estimates in
#'   subgroups),
#' \item \code{random.subgroup} (random effects estimates in
#'   subgroups),
#' \item \code{prediction.subgroup} (prediction interval in
#'   subgroups),
#' \item \code{test.effect.subgroup} (test for effect in subgroups),
#' \item \code{test.effect.subgroup.common} (test for effect in
#'   subgroups, common effect model),
#' \item \code{test.effect.subgroup.random} (test for effect in
#'   subgroups, random effects model).
#' }
#' }
#'
#' \subsection{Data set with information to print with overall results}{
#' 
#' Argument \code{data.pooled} can be used to provide information printed in
#' the line(s) with overall results, i.e., for the common effect or random
#' effects model. Input must be a data frame with variable names equal to those
#' provided in arguments \code{leftcols} or \code{rightcols}. Only variables
#' for additional variables are considered.
#' 
#' It is possible to provide a row in data set \code{data.pooled} for each
#' common effect or random effects estimate. The order in the data set
#' corresponds to the order of common effect and random effects estimates in
#' the forest plot, i.e., common effect followed by random effects estimates.
#' If the data set contains a single row, the value provided for a variable is
#' considered for all printed common effect and random effects estimates 
#' 
#' In meta-analyses with subgroups, a row must be provided in data set
#' \code{data.pooled} for each overall common effect or random effects estimate,
#' followed by common effect or random effects estimates within subgroups. The
#' order for subgroup results is the same as see in the forest plot.
#' }
#' 
#' \subsection{Additional general settings}{
#' 
#' Arguments \code{text.common}, \code{text.random}, and
#' \code{text.predict} can be used to change the label to identify
#' overall results (common effect and random effects model as well as
#' prediction interval). By default the following text is printed:
#' \itemize{
#' \item "Common effect model" (argument \code{text.common})
#' \item "Random effects model" (\code{text.random})
#' \item "Prediction interval" (\code{text.predict})
#' }
#'
#' If confidence interval levels are different for individual studies,
#' meta-analysis, and prediction interval (arguments \code{level},
#' \code{level.ma}, \code{level.predict} in meta-analysis functions,
#' e.g., \code{\link{metabin}}), additional information is printed,
#' e.g., " (99\%-CI)" for a 99\% confidence interval in the
#' meta-analysis.
#' 
#' Argument \code{pscale} can be used to rescale single proportions or
#' risk differences, e.g., \code{pscale = 1000} means that proportions
#' are expressed as events per 1000 observations. This is useful in
#' situations with (very) low event probabilities.
#' 
#' Argument \code{irscale} can be used to rescale single rates or rate
#' differences, e.g., \code{irscale = 1000} means that rates are
#' expressed as events per 1000 time units, e.g., person-years. This is
#' useful in situations with (very) low rates. Argument \code{irunit}
#' can be used to specify the time unit used in individual studies
#' (default: "person-years"). This information is printed in summaries
#' and forest plots if argument \code{irscale} is not equal to 1.
#' }
#' 
#' \subsection{Forest plots in RevMan5 layout}{
#' 
#' If argument \code{layout = "RevMan5"} (and arguments \code{leftcols} and
#' \code{rightcols} are \code{NULL}), the layout for forest plots used for
#' Cochrane reviews (which were generated with Review Manager 5) is reproduced:
#' \enumerate{
#' \item All columns are printed on the left side of the forest plot
#'   (see arguments \code{leftcols} and \code{rightcols})
#' \item Tests for overall effect and subgroup differences are printed
#'   (\code{test.overall}, \code{test.effect.subgroup},
#'   \code{test.subgroup})
#' \item Diamonds representing meta-analysis results are printed in
#'   black (\code{diamond.common}, \code{diamond.random})
#' \item Colour of squares depends on the meta-analysis object
#'   (\code{col.square}, \code{col.square.lines})
#' \item Information on effect measure and meta-analysis method is
#'   printed above the forest plot (\code{smlab})
#' \item Label "Study or Subgroup" is printed for meta-analysis with
#'   subgroups (\code{leftlabs})
#' }
#' }
#'
#' \subsection{Forest plots in JAMA layout}{
#' 
#' If argument \code{layout = "JAMA"} (and arguments \code{leftcols} and
#' \code{rightcols} are \code{NULL}), instructions for authors of the
#' \emph{Journal of the American Medical Association} are taken into account:
#' \enumerate{
#' \item Graph labels on right and left side are printed in bold font
#'   at top of forest plot (see arguments \code{bottom.lr} and
#'   \code{ff.lr})
#' \item Information on effect measure and level of confidence
#'   interval is printed at bottom of forest plot (\code{xlab})
#' \item Tests for overall effect are printed (\code{test.overall})
#' \item Diamonds representing meta-analysis results are printed in
#'   lightblue (\code{diamond.common}, \code{diamond.random})
#' \item Squares representing individual study results are printed in
#'   darkblue (\code{col.square}, \code{col.square.lines})
#' \item Between-study variance \eqn{\tau^2} is not printed
#' \item Empty rows are omitted (\code{addrow}, \code{addrow.overall},
#'   \code{addrow.subgroups})
#' \item Label "Source" is printed instead of "Study" (\code{leftlabs})
#' \item P-values are printed without leading zeros (\code{zero.pval})
#' \item P-values are rounded to three digits (for 0.001 < p \eqn{\le}
#'   0.01) or two digits (p > 0.01) (\code{JAMA.pval})
#' }
#' Study labels according to JAMA guidelines can be generated using
#' \code{\link{labels.meta}}.
#' }
#'
#' \subsection{Forest plots showing results of subgroups}{
#' 
#' The following changes are conducted if argument
#' \code{layout = "subgroup"} (and arguments \code{leftcols} and
#' \code{rightcols} are \code{NULL}) and a subgroup analysis was
#' conducted:
#' \enumerate{
#' \item Individual study results are omitted (see argument
#'   \code{study.results})
#' \item Total number of observations is not printed
#'   (\code{pooled.totals})
#' \item Label "Subgroup" is printed instead of "Study"
#'   (\code{leftlabs})
#' }
#' } 
#'
#' @note
#' R function \code{.forestArgs} generates a character vector with all
#' arguments of \code{forest.meta}.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}, \code{\link{metacont}},
#'   \code{\link{metagen}}, \code{\link{forest.metabind}},
#'   \code{\link{settings.meta}}, \code{\link{labels.meta}}
#' 
#' @references
#' Guddat C, Grouven U, Bender R, Skipka G (2012):
#' A note on the graphical presentation of prediction intervals in
#' random-effects meta-analyses.
#' \emph{Systematic Reviews},
#' \bold{1}, 34
#' 
#' Higgins JPT, Thompson SG, Spiegelhalter DJ (2009): 
#' A re-evaluation of random-effects meta-analysis.
#' \emph{Journal of the Royal Statistical Society: Series A},
#' \bold{172}, 137-59
#' 
#' @keywords hplot
#' 
#' @examples
#' data(Olkin1995)
#' m1 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
#'   data = Olkin1995, subset = c(41, 47, 51, 59),
#'   sm = "RR", method = "I",
#'   studlab = paste(author, year))
#' 
#' 
#' \dontrun{
#' # Do standard (symmetric) forest plot
#' #
#' forest(m1)
#' }
#' 
#' # Layout of forest plot similar to Review Manager 5
#' #
#' # Furthermore, add labels on both sides of forest plot and
#' # prediction interval
#' #
#' forest(m1, layout = "RevMan5", common = FALSE,
#'   label.left = "Favours experimental", col.label.left = "green",
#'   label.right = "Favours control", col.label.right = "red",
#'   prediction = TRUE)
#' 
#' 
#' \dontrun{
#' # Create PDF files with the forest plot
#' #
#' # - specify filename (R function pdf() is used due to extension .pdf)
#' # - height of the figure is automatically determined
#' # - width is set to 10 inches
#' forest(m1, file = "forest-m1-1.pdf", width = 10)
#' #
#' # - specify graphics device function
#' #   (filename "Rplots.pdf" used, see help page of R function pdf())
#' # - height of the figure is automatically determined
#' # - width is set to 10 inches
#' # - set title for PDF file
#' # - set background of forest plot
#' forest(m1, func.gr = pdf, width = 10,
#'   args.gr = list(title = "My Forest Plot", bg = "green"))
#' #
#' # - manually specify the height of the figure
#' pdf("forest-m1-2.pdf", width = 10, height = 3)
#' forest(m1)
#' dev.off()
#' 
#' # Define equivalence limits: 0.75 and 1 / 0.75
#' #
#' forest(m1, layout = "RevMan5", common = FALSE, cid = 0.75,
#'   fill = "lightgray",
#'   fill.cid = "white")
#' 
#' # Fill regions with beneficial and detrimental effects
#' #
#' forest(m1, layout = "RevMan5", common = FALSE, cid = 0.75,
#'   fill = "lightgray",
#'   fill.cid.below.null = "green",
#'   fill.cid.above.null = "red")
#' 
#' # Define thresholds for small, moderate and large effects
#' # and use hcl.colors() to define colours to fill regions
#' #
#' thresholds <- c(0.25, 0.5, 0.75)
#' n.cols <- length(thresholds)
#' forest(m1, layout = "RevMan5", common = FALSE,
#'   label.left = "Desirable effect", 
#'   label.right = "Undesirable effect", 
#'   lty.cid = 3, col.cid = "darkgray",
#'   cid.below.null = thresholds, cid.above.null = 1 / rev(thresholds),
#'   fill.cid.below.null =
#'     hcl.colors(n.cols, palette = "Blues 2", alpha = 0.6),
#'   fill.cid.above.null =
#'     hcl.colors(n.cols, palette = "Oranges", alpha = 0.6, rev = TRUE))
#' 
#' # Conduct subgroup meta-analysis
#' #
#' m2 <- update(m1,
#'   subgroup = ifelse(year < 1987, "Before 1987", "1987 and later"),
#'   print.subgroup.name = FALSE)
#' 
#' # Show summary results for subgroups with at least two studies
#' #
#' forest(m2, sortvar = -TE, random = FALSE)
#' 
#' # Show results for all subgroups
#' #
#' forest(m2, sortvar = -TE, random = FALSE, subgroup = k.w >= 1)
#' 
#' # Forest plot specifying argument xlim
#' #
#' forest(m1, xlim = c(0.01, 10))
#' 
#' # Print results of test for overall effect
#' #
#' forest(m1, test.overall.common = TRUE, test.overall.random = TRUE)
#' 
#' # Forest plot with 'classic' layout used in R package meta,
#' # version < 1.6-0
#' #
#' forest(m1, col.square = "black", hetstat = FALSE)
#' 
#' # Change set of columns printed on left side of forest plot
#' # (resulting in overlapping text)
#' #
#' forest(m1, random = FALSE, leftcols = "studlab")
#' 
#' # Use argument 'calcwidth.hetstat' to consider text for heterogeneity
#' # measures in width of column with study labels
#' #
#' forest(m1, random = FALSE, leftcols = "studlab",
#'   calcwidth.hetstat = TRUE)
#' 
#' # Use argument 'addrows.below.overall' to manually add two empty
#' # rows
#' #
#' forest(m1, random = FALSE, leftcols = "studlab", addrows = 2)
#' 
#' # Do not print columns on right side of forest plot
#' #
#' forest(m1, rightcols = FALSE)
#' 
#' # Change study label to "Author"
#' #
#' forest(m1, random = FALSE, leftlabs = c("Author", NA, NA, NA, NA))
#' 
#' # Just give effect estimate and 95% confidence interval on right
#' # side of forest plot (in one column)
#' #
#' forest(m1, rightcols = "effect.ci")
#' 
#' # Just give effect estimate and 95% confidence interval on right
#' # side of forest plot
#' #
#' forest(m1, rightcols = c("effect", "ci"))
#' 
#' # 1. Change order of columns on left side
#' # 2. Attach labels to columns 'event.e' and 'event.c' instead of
#' #    columns 'n.e' and 'n.c'
#' #
#' forest(m1,
#'   leftcols = c("studlab", "n.e", "event.e", "n.c", "event.c"),
#'   label.e.attach = "event.e", label.c.attach = "event.c")
#' 
#' # Specify column labels only for variables 'year' and 'author'
#' # (and define digits for additional variables)
#' #
#' forest(m1,
#'   leftcols = c("studlab", "event.e", "n.e", "event.c", "n.c", "author", "year"),
#'   leftlabs = c("Author", "Year of Publ"))
#' 
#' # Center text in all columns
#' #
#' forest(m1,
#'   leftcols = c("studlab", "event.e", "n.e", "event.c", "n.c",
#'                "author", "year"),
#'   leftlabs = c("Author", "Year of Publ"), hetstat = FALSE,
#'   just = "center", just.addcols = "center", just.studlab = "center")
#' 
#' # Same result
#' #
#' forest(m1,
#'   leftcols = c("studlab", "event.e", "n.e", "event.c", "n.c",
#'              "author", "year"),
#'   leftlabs = c("Author", "Year of Publ"), hetstat = FALSE,
#'   just = "c", just.addcols = "c", just.studlab = "c")
#' 
#' # Change some fontsizes and fontfaces
#' #
#' forest(m1,
#'   fs.study = 10, ff.study = "italic",
#'   fs.study.label = 11, ff.study.label = "bold",
#'   fs.axis = 5, ff.axis = "italic",
#'   ff.smlab = "bold.italic",
#'   ff.common = "plain", ff.hetstat = "plain")
#' 
#' # Change some colours
#' #
#' forest(m1,
#'   col.diamond = "green", col.diamond.lines = "red",
#'   col.study = c("green", "blue", "red", "orange"),
#'   col.square = "pink", col.square.lines = "black")
#' 
#' # Sort by weight in common effect model
#' #
#' forest(m1, sortvar = w.common, random = FALSE)
#' 
#' # Sort by decreasing weight in common effect model
#' #
#' forest(m1, sortvar = -w.common, random = FALSE)
#' 
#' # Sort by size of treatment effect
#' #
#' forest(m1, sortvar = TE, random = FALSE)
#' 
#' # Sort by size of treatment effect
#' #
#' forest(m1, sortvar = -TE, random = FALSE)
#' 
#' # Sort by decreasing year of publication
#' #
#' forest(m1, sortvar = -year, random = FALSE)
#' 
#' # Print results of test for subgroup differences (random effects
#' # model)
#' #
#' forest(m2, sortvar = -TE, common = FALSE)
#' 
#' # Print only subgroup results
#' #
#' forest(m2, layout = "subgroup")
#' 
#' # Print only subgroup results (and consider text for tests of
#' # subgroup differences in width of subgroup column)
#' #
#' forest(m2, layout = "subgroup", calcwidth.tests = TRUE)
#' 
#' # Print only subgroup results (and consider text for heterogeneity
#' # in width of subgroup column)
#' #
#' forest(m2, layout = "subgroup", calcwidth.hetstat = TRUE)
#' }
#'
#' @method forest meta
#' @export


forest.meta <- function(x,
                        sortvar,
                        studlab = TRUE,
                        #
                        layout = gs("layout"),
                        #
                        common = x$common,
                        random = x$random,
                        overall = x$overall,
                        text.common = x$text.common,
                        text.random = x$text.random,
                        lty.common = gs("lty.common"),
                        lty.random = gs("lty.random"),
                        col.common = gs("col.common"),
                        col.random = gs("col.random"),
                        text.w.common = x$text.w.common,
                        text.w.random = x$text.w.random,
                        #
                        prediction = x$prediction,
                        text.predict = x$text.predict,
                        #
                        subgroup = TRUE,
                        subgroup.hetstat =
                          subgroup & (is.character(hetstat) || hetstat),
                        print.subgroup.labels = TRUE,
                        #
                        subgroup.name = x$subgroup.name,
                        print.subgroup.name = x$print.subgroup.name,
                        sep.subgroup = x$sep.subgroup,
                        text.common.w = text.common,
                        text.random.w = text.random,
                        text.predict.w = text.predict,
                        sort.subgroup = gs("sort.subgroup"),
                        #
                        pooled.totals = common | random,
                        pooled.events = gs("pooled.events"),
                        pooled.times = gs("pooled.times"),
                        #
                        study.results = gs("study.results"),
                        #
                        rob = x$rob,
                        rob.text = "Risk of Bias",
                        rob.xpos = 0,
                        rob.legend = TRUE,
                        rob.only = FALSE,
                        #
                        xlab = "", xlab.pos,
                        smlab = NULL, smlab.pos, xlim,
                        #
                        allstudies = TRUE,
                        weight.study = NULL,
                        pscale = x$pscale,
                        irscale = x$irscale, irunit = x$irunit,
                        #
                        file = NULL, width = gs("width"), rows.gr = NULL,
                        func.gr = NULL, args.gr = NULL,
                        dev.off = NULL,
                        #
                        ref,
                        #
                        cid = gs("cid"),
                        cid.below.null = gs("cid.below.null"),
                        cid.above.null = gs("cid.above.null"),
                        lty.cid = gs("lty.cid"),
                        col.cid = gs("col.cid"),
                        fill.cid = gs("fill.cid"),
                        fill.cid.below.null = fill.cid,
                        fill.cid.above.null = rev(fill.cid),
                        cid.pooled.only = gs("cid.pooled.only"),
                        #
                        fill = gs("fill"),
                        #
                        fill.equi = gs("fill.equi"),
                        fill.lower.equi = fill.equi,
                        fill.upper.equi = rev(fill.equi),
                        #
                        leftcols = gs("leftcols"),
                        rightcols = gs("rightcols"),
                        leftlabs = gs("leftlabs"),
                        rightlabs = gs("rightlabs"),
                        #
                        label.e = x$label.e,
                        label.c = x$label.c,
                        #
                        label.e.attach = gs("label.e.attach"),
                        label.c.attach = gs("label.c.attach"),
                        #
                        label.left = x$label.left,
                        label.right = x$label.right,
                        bottom.lr = gs("bottom.lr"),
                        #
                        lab.NA = gs("lab.NA"),
                        lab.NA.effect = gs("lab.NA.effect"),
                        lab.NA.weight = gs("lab.NA.weight"),
                        #
                        lwd = gs("lwd"),
                        #
                        at = NULL,
                        label = TRUE,
                        col.label = gs("col.label"),
                        #
                        type.study = gs("type.study"),
                        type.common = gs("type.common"),
                        type.random = type.common,
                        type.subgroup =
                          ifelse(study.results, "diamond", "square"),
                        type.subgroup.common = type.subgroup,
                        type.subgroup.random = type.subgroup,
                        #
                        col.study = gs("col.study"),
                        col.square = gs("col.square"),
                        col.square.lines = gs("col.square.lines"),
                        col.circle = gs("col.circle"),
                        col.circle.lines = col.circle,
                        #
                        col.inside = gs("col.inside"),
                        col.inside.common = col.inside,
                        col.inside.random = col.inside,
                        #
                        col.diamond = gs("col.diamond"),
                        col.diamond.common = col.diamond,
                        col.diamond.random = col.diamond,
                        col.diamond.lines = gs("col.diamond.lines"),
                        col.diamond.lines.common = col.diamond.lines,
                        col.diamond.lines.random = col.diamond.lines,
                        #
                        col.predict = gs("col.predict"),
                        col.predict.lines = gs("col.predict.lines"),
                        #
                        col.subgroup = gs("col.subgroup"),
                        #
                        col.label.left = x$col.label.left,
                        col.label.right = x$col.label.right,
                        #
                        hetstat = common | random | overall.hetstat,
                        overall.hetstat =
                          x$overall.hetstat & !inherits(x, "metamerge"),
                        hetlab = gs("hetlab"),
                        resid.hetstat = gs("resid.hetstat"),
                        resid.hetlab = gs("resid.hetlab"),
                        #
                        print.I2 = gs("forest.I2"),
                        print.I2.ci = gs("forest.I2.ci"),
                        print.tau2 = gs("forest.tau2"),
                        print.tau2.ci = gs("forest.tau2.ci"),
                        print.tau = gs("forest.tau"),
                        print.tau.ci = gs("forest.tau.ci"),
                        print.Q = gs("forest.Q"),
                        print.pval.Q = gs("forest.pval.Q"),
                        print.Rb = gs("forest.Rb"),
                        print.Rb.ci = gs("forest.Rb.ci"),
                        text.subgroup.nohet = gs("text.subgroup.nohet"),
                        #
                        LRT = gs("LRT"),
                        #
                        test.overall = gs("test.overall"),
                        test.overall.common = common & overall & test.overall,
                        test.overall.random =
                          random & overall & test.overall,
                        label.test.overall.common,
                        label.test.overall.random,
                        #
                        print.stat = gs("forest.stat"),
                        #
                        test.subgroup = x$test.subgroup,
                        test.subgroup.common = test.subgroup & common,
                        test.subgroup.random = test.subgroup & random,
                        #
                        common.subgroup = common,
                        random.subgroup = random,
                        prediction.subgroup = x$prediction.subgroup,
                        #
                        print.Q.subgroup = gs("forest.Q.subgroup"),
                        label.test.subgroup.common,
                        label.test.subgroup.random,
                        #
                        test.effect.subgroup = gs("test.effect.subgroup"),
                        test.effect.subgroup.common,
                        test.effect.subgroup.random,
                        label.test.effect.subgroup.common,
                        label.test.effect.subgroup.random,
                        #
                        text.addline1,
                        text.addline2,
                        #
                        details = gs("forest.details"),
                        #
                        col.lines = gs("col.lines"),
                        header.line,
                        col.header.line = col.lines,
                        col.jama.line = col.subgroup,
                        #
                        data.pooled = NULL,
                        #
                        fontsize = gs("fontsize"),
                        fontfamily = gs("fontfamily"),
                        fs.heading = fontsize,
                        fs.common = gs("fs.common"),
                        fs.random = gs("fs.random"),
                        fs.predict = gs("fs.predict"),
                        fs.common.labels = gs("fs.common.labels"),
                        fs.random.labels = gs("fs.random.labels"),
                        fs.predict.labels = gs("fs.predict.labels"),
                        fs.study = fontsize,
                        fs.study.labels = fs.study,
                        fs.hetstat = gs("fs.hetstat"),
                        fs.test.overall = gs("fs.test.overall"),
                        fs.test.subgroup = gs("fs.test.subgroup"),
                        fs.test.effect.subgroup = gs("fs.test.effect.subgroup"),
                        fs.addline = gs("fs.addline"),
                        fs.axis = fontsize,
                        fs.smlab = fontsize,
                        fs.xlab = fontsize,
                        fs.lr = fontsize,
                        fs.rob = fontsize,
                        fs.rob.symbols = fontsize,
                        fs.details = fontsize,
                        #
                        ff.heading = "bold",
                        ff.common = gs("ff.common"),
                        ff.random = gs("ff.random"),
                        ff.predict = gs("ff.predict"),
                        ff.common.labels = gs("ff.common.labels"),
                        ff.random.labels = gs("ff.random.labels"),
                        ff.predict.labels = gs("ff.predict.labels"),
                        ff.study = "plain",
                        ff.study.labels = ff.study,
                        ff.hetstat = gs("ff.hetstat"),
                        ff.test.overall = gs("ff.test.overall"),
                        ff.test.subgroup = gs("ff.test.subgroup"),
                        ff.test.effect.subgroup = gs("ff.test.effect.subgroup"),
                        ff.addline = gs("ff.addline"),
                        ff.axis = gs("ff.axis"),
                        ff.smlab = gs("ff.smlab"),
                        ff.xlab = gs("ff.xlab"),
                        ff.lr = gs("ff.lr"),
                        ff.rob = "plain",
                        ff.rob.symbols = "bold",
                        ff.details = "plain",
                        #
                        squaresize =
                          if (layout == "BMJ") 0.9 / spacing else 0.8 / spacing,
                        lwd.square = gs("lwd.square"),
                        lwd.diamond = gs("lwd.diamond"),
                        #
                        arrow.type = gs("arrow.type"),
                        arrow.length = gs("arrow.length"),
                        #
                        plotwidth =
                          if (layout %in% c("BMJ", "JAMA")) "8cm" else "6cm",
                        colgap = gs("colgap"),
                        colgap.left = colgap,
                        colgap.right = colgap,
                        colgap.studlab = colgap.left,
                        colgap.forest = gs("colgap.forest"),
                        colgap.forest.left = colgap.forest,
                        colgap.forest.right = colgap.forest,
                        colgap.rob = "1mm",
                        colgap.rob.overall = "2mm",
                        #
                        calcwidth.pooled =
                          (common | random) &
                          (overall | !is.null(x$subgroup)),
                        calcwidth.common = calcwidth.pooled,
                        calcwidth.random = calcwidth.pooled,
                        calcwidth.predict = gs("calcwidth.predict"),
                        calcwidth.hetstat = gs("calcwidth.hetstat"),
                        calcwidth.tests  = gs("calcwidth.tests"),
                        calcwidth.subgroup = gs("calcwidth.subgroup"),
                        calcwidth.addline = gs("calcwidth.addline"),
                        #
                        just = if (layout == "JAMA") "left" else "right",
                        just.studlab = gs("just.studlab"),
                        just.addcols = gs("just.addcols"),
                        just.addcols.left = just.addcols,
                        just.addcols.right = just.addcols,
                        just.label.e = just,
                        just.label.c = just,
                        #
                        bmj.text = NULL,
                        bmj.xpos = 0,
                        bmj.sep = " / ",
                        #
                        spacing = gs("spacing"),
                        addrow = gs("addrow"),
                        addrow.overall = gs("addrow.overall"),
                        addrow.subgroups = gs("addrow.subgroups"),
                        addrows.below.overall = gs("addrows.below.overall"),
                        #
                        new = TRUE,
                        #
                        backtransf = x$backtransf,
                        digits = gs("digits.forest"),
                        digits.se = gs("digits.se"),
                        digits.stat = gs("digits.stat"),
                        digits.pval = gs("digits.pval"),
                        digits.pval.Q = gs("digits.pval.Q"),
                        digits.Q = gs("digits.Q"),
                        digits.tau2 = gs("digits.tau2"),
                        digits.tau = gs("digits.tau"),
                        digits.I2 = gs("digits.I2"),
                        digits.weight = gs("digits.weight"),
                        #
                        digits.mean = gs("digits.mean"),
                        digits.sd = gs("digits.sd"),
                        digits.cor = digits,
                        digits.time = digits,
                        #
                        digits.n = 0,
                        digits.event = 0,
                        #
                        digits.TE = gs("digits.TE.forest"),
                        #
                        digits.addcols = digits,
                        digits.addcols.right = digits.addcols,
                        digits.addcols.left = digits.addcols,
                        #
                        scientific.pval = gs("scientific.pval"),
                        big.mark = gs("big.mark"),
                        zero.pval =
                          if (layout == "JAMA") FALSE else gs("zero.pval"),
                        JAMA.pval =
                          if (layout == "JAMA") TRUE else gs("JAMA.pval"),
                        #
                        warn.deprecated = gs("warn.deprecated"),
                        ...) {
  
  
  #
  #
  # (1) Check for meta object and upgrade older meta objects
  #
  #
  chkclass(x, "meta")
  x.name <- deparse(substitute(x))
  x <- updateversion(x)
  #
  K.all <- length(x$TE)
  #
  sm <- x$sm
  #
  metabin <- inherits(x, "metabin")
  metacont <- inherits(x, "metacont")
  metacor <- inherits(x, "metacor")
  metagen <- inherits(x, "metagen")
  metainc <- inherits(x, "metainc")
  metamean <- inherits(x, "metamean")
  metaprop <- inherits(x, "metaprop")
  metarate <- inherits(x, "metarate")
  metabind <- inherits(x, "is.metabind")
  metacum <- inherits(x, "metacum")
  metainf <- inherits(x, "metainf")
  #
  metamerge <- inherits(x, "metamerge")
  #
  meta <- !metabind &&
    (metabin | metacont | metacor | metagen | metainc | metamean |
       metaprop | metarate)
  #
  ftr <- x$func.transf
  atr <- x$args.transf
  fbt <- x$func.backtransf
  abt <- x$args.backtransf
  #
  args <- list(...)
  # Check whether first argument is a list. In this case only use
  # this list as input.
  if (length(args) > 0 && is.list(args[[1]]))
    args <- args[[1]]
  nam.args <- names(args)
  #
  # Logical variables for missing arguments
  #
  missing.calcwidth.common <- missing(calcwidth.common)
  missing.col.circle <- missing(col.circle) | is.null(col.circle)
  missing.col.circle.lines <-
    missing(col.circle.lines) | is.null(col.circle.lines)
  missing.col.common <- missing(col.common)
  missing.col.diamond <- missing(col.diamond)
  missing.col.diamond.common <-
    missing(col.diamond.common) | is.null(col.diamond.common)
  missing.col.diamond.fixed <- is.na(argid(nam.args, "col.diamond.fixed"))
  missing.col.diamond.fixed.lines <-
    is.na(argid(nam.args, "col.diamond.fixed.lines"))
  missing.col.diamond.lines.common <- missing(col.diamond.lines.common)
  missing.col.diamond.lines.fixed <-
    is.na(argid(nam.args, "col.diamond.lines.fixed"))
  missing.col.diamond.lines.random <- missing(col.diamond.lines.random)
  missing.col.diamond.random <-
    missing(col.diamond.random) | is.null(col.diamond.random)
  missing.col.inside <- missing(col.inside)
  missing.col.inside.common <- missing(col.inside.common)
  missing.col.square <- missing(col.square) | is.null(col.square)
  missing.col.square.lines <-
    missing(col.square.lines) | is.null(col.square.lines)
  missing.col.subgroup <- missing(col.subgroup)
  missing.common <- missing(common)
  missing.common.subgroup <- missing(common.subgroup)
  missing.digits.TE <- missing(digits.TE)
  missing.digits.addcols <- missing(digits.addcols)
  missing.digits.addcols.left <- missing(digits.addcols.left)
  missing.digits.addcols.right <- missing(digits.addcols.right)
  missing.digits.cor <- missing(digits.cor)
  missing.digits.mean <- missing(digits.mean)
  missing.digits.sd <- missing(digits.sd)
  missing.digits.stat <- missing(digits.stat)
  missing.digits.time <- missing(digits.time)
  missing.ff.addline <- missing(ff.addline)
  missing.ff.common <- missing(ff.common)
  missing.ff.common.labels <- missing(ff.common.labels)
  missing.ff.fixed <- is.na(argid(nam.args, "ff.fixed"))
  missing.ff.fixed.labels <- is.na(argid(nam.args, "ff.fixed.labels"))
  missing.ff.hetstat <- missing(ff.hetstat)
  missing.ff.lr <- missing(ff.lr)
  missing.ff.predict <- missing(ff.predict)
  missing.ff.predict.labels <- missing(ff.predict.labels)
  missing.ff.random <- missing(ff.random)
  missing.ff.random.labels <- missing(ff.random.labels)
  missing.ff.test.effect.subgroup <- missing(ff.test.effect.subgroup)
  missing.ff.test.overall <- missing(ff.test.overall)
  missing.ff.test.subgroup <- missing(ff.test.subgroup)
  missing.fs.addline <- missing(fs.addline)
  missing.fs.common <- missing(fs.common)
  missing.fs.common.labels <- missing(fs.common.labels)
  missing.fs.fixed <- is.na(argid(nam.args, "fs.fixed"))
  missing.fs.fixed.labels <- is.na(argid(nam.args, "fs.fixed.labels"))
  missing.fs.hetstat <- missing(fs.hetstat)
  missing.fs.predict <- missing(fs.predict)
  missing.fs.predict.labels <- missing(fs.predict.labels)
  missing.fs.random <- missing(fs.random)
  missing.fs.random.labels <- missing(fs.random.labels)
  missing.fs.test.effect.subgroup <- missing(fs.test.effect.subgroup)
  missing.fs.test.overall <- missing(fs.test.overall)
  missing.fs.test.subgroup <- missing(fs.test.subgroup)
  missing.header.line <- missing(header.line)
  missing.hetstat <- missing(hetstat)
  missing.irscale <- missing(irscale)
  missing.just <- missing(just)
  missing.label <- missing(label)
  missing.label.c <- missing(label.c)
  missing.label.c.attach <- missing(label.c.attach)
  missing.label.e <- missing(label.e)
  missing.label.e.attach <- missing(label.e.attach)
  missing.label.test.effect.subgroup.common <-
    missing(label.test.effect.subgroup.common)
  missing.label.test.effect.subgroup.fixed <-
    is.na(argid(nam.args, "label.test.effect.subgroup.fixed"))
  missing.label.test.effect.subgroup.random <-
    missing(label.test.effect.subgroup.random)
  missing.label.test.overall.common <- missing(label.test.overall.common)
  missing.label.test.overall.fixed <-
    is.na(argid(nam.args, "label.test.overall.fixed"))
  missing.label.test.overall.random <- missing(label.test.overall.random)
  missing.label.test.subgroup.common <- missing(label.test.subgroup.common)
  missing.label.test.subgroup.fixed <-
    is.na(argid(nam.args, "label.test.subgroup.fixed"))
  missing.label.test.subgroup.random <- missing(label.test.subgroup.random)
  missing.leftcols <- missing(leftcols)
  missing.leftlabs <- missing(leftlabs)
  missing.lty.common <- missing(lty.common)
  missing.overall.hetstat <- missing(overall.hetstat)
  missing.pooled.totals <- missing(pooled.totals)
  missing.prediction.subgroup <- missing(prediction.subgroup)
  missing.print.I2 <- missing(print.I2)
  missing.print.Q <- missing(print.Q)
  missing.print.Rb <- missing(print.Rb)
  missing.print.pval.Q <- missing(print.pval.Q)
  missing.print.stat <- missing(print.stat)
  missing.print.tau <- missing(print.tau)
  missing.print.tau2 <- missing(print.tau2)
  missing.pscale <- missing(pscale)
  missing.random.subgroup <- missing(random.subgroup)
  missing.ref <- missing(ref)
  missing.resid.hetstat <- missing(resid.hetstat)
  missing.rightcols <- missing(rightcols)
  missing.rightlabs <- missing(rightlabs)
  missing.smlab.pos <- missing(smlab.pos)
  missing.sort.subgroup <- missing(sort.subgroup)
  missing.studlab <- missing(studlab)
  missing.study.results <- missing(study.results)
  missing.subgroup <- missing(subgroup)
  missing.subgroup.hetstat <- missing(subgroup.hetstat)
  missing.test.effect.subgroup <- missing(test.effect.subgroup)
  missing.test.effect.subgroup.common <- missing(test.effect.subgroup.common)
  missing.test.effect.subgroup.random <- missing(test.effect.subgroup.random)
  missing.test.overall.common <- missing(test.overall.common)
  missing.test.subgroup.common <- missing(test.subgroup.common)
  missing.text.addline1 <- missing(text.addline1)
  missing.text.addline2 <- missing(text.addline2)
  missing.text.common <- missing(text.common)
  missing.text.common.w <- is.na(argid(nam.args, "text.common.w"))
  missing.text.common.w <- missing(text.common.w)
  missing.text.fixed <- is.na(argid(nam.args, "text.fixed"))
  missing.text.predict <- missing(text.predict)
  missing.text.random <- missing(text.random)
  missing.text.w.common <- missing(text.w.common)
  missing.type.common <- missing(type.common)
  missing.type.subgroup <- missing(type.subgroup)
  missing.type.subgroup.common <- missing(type.subgroup.common)
  missing.type.subgroup.random <- missing(type.subgroup.random)
  missing.weight.study <- missing(weight.study)
  missing.xlab.pos <- missing(xlab.pos)
  #
  notavail.digits.addcols.left <-
    missing.digits.addcols & missing.digits.addcols.left
  notavail.digits.addcols.right <-
    missing.digits.addcols & missing.digits.addcols.right
  #
  avail.xlim <- !missing(xlim)
  avail.fill.equi <- !missing(fill.equi) & !is.null(fill.equi)
  
  
  #
  #
  # (2) Extract data on subgroup analysis
  #
  #
  by <- !is.null(x$subgroup)
  #
  if (by) {
    n.by <- length(x$subgroup.levels)
    #
    TE.common.w <- collapsemat(x$TE.common.w)
    statistic.common.w <- collapsemat(x$statistic.common.w)
    pval.common.w <- collapsemat(x$pval.common.w)
    lower.common.w <- collapsemat(x$lower.common.w)
    upper.common.w <- collapsemat(x$upper.common.w)
    #
    if (!is.matrix(x$TE.random.w) & is.matrix(x$seTE.random.w)) {
      TE.random.w <-
        collapsemat(matrix(x$TE.random.w,
                           nrow = nrow(x$seTE.random.w),
                           ncol = ncol(x$seTE.random.w),
                           dimnames = list(rownames(x$seTE.random.w),
                                           colnames(x$seTE.random.w))))
    }
    else
      TE.random.w <- collapsemat(x$TE.random.w)
    #
    lower.random.w <- collapsemat(x$lower.random.w)
    upper.random.w <- collapsemat(x$upper.random.w)
    statistic.random.w <- collapsemat(x$statistic.random.w)
    pval.random.w <- collapsemat(x$pval.random.w)
    #
    lower.predict.w <- collapsemat(x$lower.predict.w)
    upper.predict.w <- collapsemat(x$upper.predict.w)
  }
  else
    n.by <- 0
  
  
  #
  #
  # (3) Determine columns on left and right side of forest plot
  #
  #
  layout <- setchar(layout, c("meta", "BMJ", "RevMan5", "JAMA", "subgroup"))
  if (layout == "subgroup" & is.null(x$subgroup)) {
    warning("Argument 'layout' set to \"meta\" (default) as ",
            "no subgroup analysis was conducted.")
    layout <- "meta"
  }
  if (layout == "subgroup") {
    if (missing.type.subgroup)
      type.subgroup <- "square"
    if (missing.type.subgroup.common)
      type.subgroup.common <- "square"
    if (missing.type.subgroup.random)
      type.subgroup.random <- "square"
    #
    if (missing.pooled.totals)
      pooled.totals <- FALSE
  }
  bmj <- layout == "BMJ"
  revman5 <- layout == "RevMan5"
  jama <- layout == "JAMA"
  revman5.jama <- revman5 | jama
  bmj.revman5 <- bmj | revman5
  bmj.revman5.jama <- bmj | revman5 | jama
  #
  lsel <- TRUE
  if (is.logical(leftcols)) {
    if (!leftcols)
      lsel <- FALSE
    #
    leftcols <- NULL
  }
  #
  if (revman5.jama)
    rsel <- FALSE
  else
    rsel <- !(is.logical(rightcols) && length(rightcols) == 1 && !rightcols)
  #
  chkchar(rob.text, length = 1)
  chknumeric(rob.xpos, length = 1)
  chklogical(rob.legend)
  chklogical(rob.only)
  #
  text.rob <- ""
  RoB.available <- !is.null(rob) && !(is.logical(rob) && !rob)
  RoB.legend <- RoB.available & rob.legend
  #
  if (RoB.available) {
    if (!inherits(rob, "rob"))
      stop("Argument 'rob' must be of class \"rob\".", call. = FALSE)
    #
    rob.labels <-
      colnames(rob)[!(colnames(rob) %in% c("Study", "Weight"))]
    rob.labels[rob.labels == "Overall"] <- "O"
    #
    rob.categories <- attr(rob, "categories")
    rob.symbols <- attr(rob, "symbols")
    rob.col <- attr(rob, "col")
    #
    text.rob <- c("Risk of bias legend", catleg(rob))
    #
    rob <- rob[, !(colnames(rob) %in% c("Study", "Weight")), drop = FALSE]
    #
    colnames(rob)[colnames(rob) == "Overall"] <- "O"
    colnames(rob) <- gsub(" ", "_", paste0("RoB.", rob.labels))
    #
    rightcols.rob <- colnames(rob)
    #
    text.rob <- text.rob[text.rob != ""]
  }
  else
    rightcols.rob <- NULL
  #
  if (!rsel)
    rightcols <- NULL
  #
  # Check for duplicate columns
  #
  if (length(c(rightcols, rightcols.rob, leftcols)) > 0 &&
      any(duplicated(c(rightcols, rightcols.rob, leftcols))))
    stop("Duplicate entries in 'leftcols' and 'rightcols'.")
  #
  # Must be of same length as labnames!!!
  #
  colnames <- c("studlab",
                "TE", "seTE",
                "cluster", "cycles",
                #
                "n.e", "n.c", "event.e", "event.c",
                "event.n.e", "event.n.c", "event.n",
                "mean.e", "mean.c", "sd.e", "sd.c",
                "mean.sd.n.e", "mean.sd.n.c", "mean.sd.n",
                #
                "cor",
                "time.e", "time.c",
                #
                "effect", "ci",
                "effect.ci",
                #
                "w.fixed", "w.common", "w.random")
  #
  # If any of the following list elements is NULL, these 'special'
  # variable names are searched for in original data set (i.e., list
  # element x$data)
  #
  colnames.notNULL <- colnames
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "cluster")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "cycles")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "n.e")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "n.c")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "event.e")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "event.c")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "mean.e")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "mean.c")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "sd.e")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "sd.c")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "cor")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "time.e")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "time.c")
  #
  ev.n.bin <-
    metabin && (bmj ||
                  "event.n.e" %in% c(leftcols, rightcols) ||
                  "event.n.c" %in% c(leftcols, rightcols))
  #
  m.s.n.cont <-
    metacont && (bmj ||
                   "mean.sd.n.e" %in% c(leftcols, rightcols) ||
                   "mean.sd.n.c" %in% c(leftcols, rightcols) )
  #
  ev.n.prop <- metaprop && (bmj || "event.n" %in% c(leftcols, rightcols))
  #
  # Identify and process columns in addition to columns
  # defined above in variable 'colnames'
  #
  colnames.new <-
    c(leftcols, rightcols, rightcols.rob)[
      !c(leftcols, rightcols, rightcols.rob) %in% colnames.notNULL]
  #
  newcols <- length(colnames.new) > 0
  #
  withstudlab <-
    (lsel & missing.leftcols) |
    (lsel & !missing.leftcols & "studlab" %in% leftcols) |
    (rsel & !missing.rightcols & "studlab" %in% rightcols)
  #
  if (!withstudlab || newcols) {
    dataset2 <- as.data.frame(x)
    #
    if (!is.null(x$data))
      dataset1 <- x$data
    else if (!is.null(x$x$data))
      dataset1 <- x$x$data
    else
      dataset1 <- dataset2
    #
    if (!is.null(x$subset)) {
      dataset1 <- dataset1[x$subset, ]
    }
    #
    if (!withstudlab && newcols) {
      #
      # Check whether first additional variable is part of
      # meta-object
      #
      firstvar <- colnames.new[1]
      colnames.new <- colnames.new[-1]
      newcols <- length(colnames.new) > 0
      #
      if (length(dataset1[[firstvar]]) == 0 &
          length(dataset2[[firstvar]]) == 0)
        stop("Variable '", firstvar, "' not available in '",
             x.name, "'.",
             call. = FALSE)
      else {        
        if (length(dataset1[[firstvar]]) != 0)
          studlab.new <- dataset1[[firstvar]]
        else if (length(dataset2[[firstvar]]) != 0)
          studlab.new <- dataset2[[firstvar]]
        #
        if (length(x$studlab) != length(studlab.new))
          stop("Variable '", firstvar, "' has different length ",
               "than list element '", x.name, "' with study labels.",
               call. = FALSE)
        #
        x$studlab <- studlab.new
        switch.l <- leftcols == firstvar
        switch.r <- rightcols == firstvar
        if (any(switch.l))
          leftcols[switch.l] <- "studlab"
        if (any(switch.r))
          rightcols[switch.r] <- "studlab"
      }
    }
  }
  #
  if (newcols) {
    #
    # Determine labels for new columns
    # 1. Use column name as label if no label is given
    # 2. Otherwise use specified labels
    #
    rightcols.new <- c(rightcols, rightcols.rob)[
      !(c(rightcols, rightcols.rob) %in% colnames.notNULL)]
    leftcols.new <- leftcols[!leftcols %in% colnames.notNULL]
  }
  
  #
  #
  # (3) Check other arguments
  #
  #
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  error <-
    try(sortvar <- catch("sortvar", mc, x, sfsp),
        silent = TRUE)
  if (inherits(error, "try-error") || is.function(error)) {
    sortvar <- catch("sortvar", mc, x$data, sfsp)
    if (isCol(x$data, ".subset"))
      sortvar <- sortvar[x$data$.subset]
  }
  #
  sort <- !is.null(sortvar)
  if (sort && (length(sortvar) != K.all))
    stop("Number of studies in object 'x' and argument 'sortvar' ",
         "have different length.")
  if (!sort)
    sortvar <- 1:K.all
  #
  if (!by) {
    if (!missing.subgroup)
      warning("Argument 'subgroup' only considered for ",
              "meta-analysis with subgroups.",
              call. = FALSE)
    if (!missing.subgroup.hetstat)
      warning("Argument 'subgroup.hetstat' only considered for ",
              "meta-analysis with subgroups.",
              call. = FALSE)
    if (!missing.common.subgroup)
      warning("Argument 'common.subgroup' only considered for ",
              "meta-analysis with subgroups.",
              call. = FALSE)
    if (!missing.random.subgroup)
      warning("Argument 'random.subgroup' only considered for ",
              "meta-analysis with subgroups.",
              call. = FALSE)
    if (!missing.prediction.subgroup)
      warning("Argument 'prediction.subgroup' only considered for ",
              "meta-analysis with subgroups.",
              call. = FALSE)
    if (!missing.test.effect.subgroup)
      warning("Argument 'test.effect.subgroup' only considered for ",
              "meta-analysis with subgroups.",
              call. = FALSE)
    if (!missing.test.effect.subgroup.common)
      warning("Argument 'test.effect.subgroup.common' only considered for ",
              "meta-analysis with subgroups.",
              call. = FALSE)
    if (!missing.test.effect.subgroup.random)
      warning("Argument 'test.effect.subgroup.random' only considered for ",
              "meta-analysis with subgroups.",
              call. = FALSE)
  }
  #
  if (!missing.studlab) {
    error <-
      try(studlab <- catch("studlab", mc, x, sfsp),
          silent = TRUE)
    if (inherits(error, "try-error")) {
      studlab <- catch("studlab", mc, x$data, NULL)
      if (isCol(x$data, ".subset"))
        studlab <- studlab[x$data$.subset]
    }
  }
  #
  slab <- TRUE
  if (length(studlab) == 1 & is.logical(studlab)) {
    if (studlab == FALSE) {
      studlab <- rep("", K.all)
      slab <- FALSE
    }
    else studlab <- x$studlab
  }
  #
  if (missing.studlab && K.all == 1 && studlab == "")
    studlab <- "1"
  #
  if (length(studlab) != K.all)
    stop("Number of studies in object 'x' and argument 'studlab' have ",
         "different length.")
  #
  # Additional arguments in '...'
  #
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  chklogical(random)
  overall <- replaceNULL(overall, common | random)
  chklogical(overall)
  #
  lty.common <- deprecated(lty.common, missing.lty.common, args, "lty.fixed",
                           warn.deprecated)
  if (!is.null(lty.common))
    chknumeric(lty.common, length = 1)
  if (!is.null(lty.random))
    chknumeric(lty.random, length = 1)
  col.common <- deprecated(col.common, missing.col.common, args, "col.fixed",
                           warn.deprecated)
  chkcolor(col.common, length = 1)
  chkcolor(col.random, length = 1)
  chklogical(prediction)
  #
  chklogical(print.subgroup.labels)
  #
  if (!is.null(print.subgroup.name))
    chklogical(print.subgroup.name)
  #
  if (!is.null(sep.subgroup))
    chkchar(sep.subgroup)
  #
  sort.subgroup <- deprecated(sort.subgroup, missing.sort.subgroup,
                              args, "bysort", warn.deprecated)
  chklogical(sort.subgroup)
  #
  chklogical(pooled.totals)
  chklogical(pooled.events)
  chklogical(pooled.times)
  chklogical(study.results)
  # chknumeric(xlab.pos) ??
  # chknumeric(smlab.pos) ??
  chklogical(allstudies)
  #
  chklogical(backtransf)
  #
  if (!is.null(pscale))
    chknumeric(pscale, length = 1)
  else
    pscale <- 1
  #
  if (!is.null(irscale))
    chknumeric(irscale, length = 1)
  else
    irscale <- 1
  #
  if (!is.null(irunit) && !is.na(irunit))
    chkchar(irunit)
  #
  if (!is.null(file))
    chkchar(file, length = 1)
  #
  if (!is.null(width))
    chknumeric(width, min = 0, zero = TRUE, length = 1)
  #
  if (!is.null(rows.gr))
    chknumeric(rows.gr, length = 1)
  else
    rows.gr <- 0
  #
  if (!is.null(args.gr))
    chklist(args.gr)
  #
  if (!is.null(dev.off))
    chklogical(dev.off)
  else if (!is.null(file) | !is.null(func.gr))
    dev.off <- TRUE
  else
    dev.off <- FALSE
  #
  # Use logarithmic x-axis?
  # (for back-transformed relative effect measures)
  is_relative <- is_relative_effect(sm) | (!is.null(fbt) && fbt == "exp")
  #
  log.xaxis <- backtransf & is_relative
  #
  if (missing.ref) {
    if (is_prop(sm) | is_rate(sm) | is_mean(sm)) {
      if (!is.null(x$null.effect))
        ref <- x$null.effect
      else
        ref <- NA
    }
    else if (log.xaxis)
      ref <- 1
    else
      ref <- 0
  }
  else
    chknumeric(ref, length = 1)
  #
  avail.cid <- !missing(cid) & !is.null(cid) & !all(is.na(cid))
  chknumeric(cid)
  #
  cid.below.null <-
    deprecated(cid.below.null, missing(cid.below.null),
               args, "lower.equi", warn.deprecated)
  #
  cid.above.null <-
    deprecated(cid.above.null, missing(cid.above.null),
               args, "upper.equi", warn.deprecated)
  #
  lty.cid <-
    deprecated(lty.cid, missing(lty.cid), args, "lty.equi", warn.deprecated)
  chknumeric(lty.cid)
  #
  col.cid <-
    deprecated(col.cid, missing(col.cid), args, "col.equi", warn.deprecated)
  chkcolor(col.cid)
  #
  avail.cid.below.null <- !is.null(cid.below.null) & !all(is.na(cid.below.null))
  avail.cid.above.null <- !is.null(cid.above.null) & !all(is.na(cid.above.null))
  #
  if (avail.cid) {
    if (any(is.na(cid)))
      stop("Missing values not allows in argument 'cid'.",
           call. = FALSE)
    #
    if (avail.cid.below.null + avail.cid.above.null == 2)
      warning("Arguments 'cid.below.null' and 'cid.above.null' ignored as ",
              "argument 'cid' is provided.",
              call. = FALSE)
    else if (avail.cid.below.null)
      warning(
        "Argument 'cid.below.null' ignored as argument 'cid' is provided.",
        call. = FALSE)
    else if (avail.cid.above.null)
      warning(
        "Argument 'cid.above.null' ignored as argument 'cid' is provided.",
        call. = FALSE)
    #
    if (any(diff(cid) <= 0))
      stop("Values for argument 'cid' must be increasing.",
           call. = FALSE)
    #
    if (any(cid < ref) & any(cid > ref))
      stop("All values provided for argument 'cid' must be either ",
           "smaller or larger than reference value of ", ref, ".",
           call. = FALSE)
    #
    if (all(cid < ref)) {
      cid.below.null <- cid
      #
      if (log.xaxis)
        cid.above.null <- rev(1 / cid)
      else
        cid.above.null <- rev(-cid)
    }
    else {
      cid.above.null <- cid
      #
      if (log.xaxis)
        cid.below.null <- rev(1 / cid)
      else
        cid.below.null <- rev(-cid)
    }
    #
    avail.cid.below.null <- TRUE
    avail.cid.above.null <- TRUE
  }
  #
  chknumeric(cid.below.null)
  chknumeric(cid.above.null)
  #
  n.cid.below.null <- sum(!is.na(cid.below.null))
  n.cid.above.null <- sum(!is.na(cid.above.null))
  #
  if (n.cid.below.null > 0) {
    max.cid.below.null <- max(cid.below.null, na.rm = TRUE)
    #
    if (!is.na(ref) && any(cid.below.null[!is.na(cid.below.null)] > ref))
      stop("All values provided for argument 'cid.below.null' must be ",
           "smaller than reference value of ", ref, ".",
           call. = FALSE)
  }
  else
    max.cid.below.null <- NA
  #
  if (n.cid.above.null > 0) {
    min.cid.above.null <- min(cid.above.null, na.rm = TRUE)
    #
    if (!is.na(ref) && any(cid.above.null[!is.na(cid.above.null)] < ref))
      stop("All values provided for argument 'cid.above.null' must be ",
           "larger than reference value of ", ref, ".",
           call. = FALSE)
  }
  else
    min.cid.above.null <- NA
  #
  if (n.cid.below.null > 0 && n.cid.above.null > 0 &&
      max.cid.below.null > min.cid.above.null)
    stop("Value", if (length(cid.below.null) > 1) "s " else " ",
         "of 'cid.below.null' must be smaller than 'cid.above.null'.",
         call. = FALSE)
  #
  if (any(cid.below.null != sort(cid.below.null), na.rm = TRUE))
    stop("Values of 'cid.below.null' must be increasing.",
         call. = FALSE)
  #
  if (any(cid.above.null != sort(cid.above.null), na.rm = TRUE))
    stop("Values of 'cid.above.null' must be increasing.",
         call. = FALSE)
  #
  # Define colours for CID regions
  #
  missing.fill.cid <- missing(fill.cid)
  missing.fill.cid.below.null <- missing(fill.cid.below.null)
  missing.fill.cid.above.null <- missing(fill.cid.above.null)
  #
  if (n.cid.below.null > 0) {
    if (length(fill.cid.below.null) == 1 & n.cid.below.null > 1)
      fill.cid.below.null <- rep(fill.cid.below.null, n.cid.below.null)
    else if (all(length(fill.cid.below.null) != n.cid.below.null))
      stop("Number of fill colours must be equal to the number of values ",
           "for 'cid.below.null'.",
           call. = FALSE)
  }
  #
  fill.cid.below.null <-
    c(fill.cid.below.null, if (avail.fill.equi) fill.equi else fill)
  #
  if (n.cid.above.null > 0) {
    if (length(fill.cid.above.null) == 1 & n.cid.above.null > 1)
      fill.cid.above.null <- rep(fill.cid.above.null, n.cid.above.null)
    else if (all(length(fill.cid.above.null) != n.cid.above.null))
      stop("Number of fill colours must be equal to the number of values ",
           "for 'cid.above.null'.",
           call. = FALSE)
  }
  #
  fill.cid.above.null <-
    c(if (avail.fill.equi) fill.equi else fill, fill.cid.above.null)
  #
  # Define colours based on equivalence limits
  #
  missing.fill.lower.equi <- missing(fill.lower.equi)
  missing.fill.upper.equi <- missing(fill.upper.equi)
  #
  if (!missing.fill.lower.equi | !missing.fill.upper.equi) {
    if (!missing.fill.cid |
        !missing.fill.cid.above.null |
        !missing.fill.cid.above.null) {
      warning("Input for arguments 'fill.lower.equi' and 'fill.upper.equi' ",
              "ignored if argument 'fill.cid', 'fill.cid.below.null' or ",
              "'fill.cid.above.null' is provided.",
              call. = FALSE)
    }
    else {
      if (n.cid.below.null > 0) {
        if (length(fill.lower.equi) == 1 & n.cid.below.null > 1)
          fill.lower.equi <- rep(fill.lower.equi, n.cid.below.null)
        else if (all(length(fill.lower.equi) != n.cid.below.null + 0:1))
          stop("Number of fill colours must be equal to the number of values ",
               "for 'cid.below.null' or +1.",
               call. = FALSE)
        #
        fill.cid.below.null <- fill.lower.equi
      }
      #
      if (n.cid.above.null > 0) {
        if (length(fill.upper.equi) == 1 & n.cid.above.null > 1)
          fill.upper.equi <- rep(fill.upper.equi, n.cid.above.null)
        else if (all(length(fill.upper.equi) != n.cid.above.null + 0:1))
          stop("Number of fill colours must be equal to the number of values ",
               "for 'cid.above.null' or +1.",
               call. = FALSE)
        #
        fill.cid.above.null <- fill.upper.equi
      }
    }
  }
  #
  chklogical(cid.pooled.only)
  #
  if (bmj)
    type.study <- "squarediamond"
  else
    type.study <- setchar(type.study,
                          c("square", "diamond", "predict", "circle",
                            "squarediamond"))
  #
  type.common <-
    setchar(type.common,
            c("square", "diamond", "predict", "circle", "squarediamond"))
  #
  type.random <-
    setchar(type.random,
            c("square", "diamond", "predict", "circle", "squarediamond"))
  #
  type.subgroup <-
    setchar(type.subgroup,
            c("square", "diamond", "predict", "circle", "squarediamond"))
  type.subgroup.common <-
    setchar(type.subgroup.common,
            c("square", "diamond", "predict", "circle", "squarediamond"))
  #
  type.subgroup.random <-
    setchar(type.subgroup.random,
            c("square", "diamond", "predict", "circle", "squarediamond"))
  #
  chklogical(bottom.lr)
  #
  chkchar(lab.NA)
  if (is.null(lab.NA.effect))
      lab.NA.effect <- ""
  #
  chkchar(lab.NA.effect)
  chkchar(lab.NA.weight)
  if (!is.null(at))
    chknumeric(at)
  #
  chkcolor(col.diamond, length = 1)
  chkcolor(col.diamond.random)
  chkcolor(col.predict)
  chkcolor(col.predict.lines)
  #
  chkcolor(col.lines)
  #
  col.subgroup <- deprecated(col.subgroup, missing.col.subgroup, args,
                             "col.by", warn.deprecated)
  #
  overall.hetstat <- replaceNULL(overall.hetstat, TRUE)
  chklogical(overall.hetstat)
  #
  col.label.left <- replaceNULL(col.label.left, gs("col.label.left"))
  col.label.right <- replaceNULL(col.label.right, gs("col.label.right"))
  #
  if (is.null(print.I2)) {
    if (is.character(hetstat) || hetstat || overall.hetstat)
      print.I2 <- TRUE
    else
      print.I2 <- FALSE
  }
  else
    chklogical(print.I2)
  #
  chklogical(print.I2.ci)
  #
  if (is.null(print.tau2)) {
    if (is.character(hetstat) || hetstat || overall.hetstat)
      print.tau2 <- TRUE
    else
      print.tau2 <- FALSE
  }
  else
    chklogical(print.tau2)
  #
  chklogical(print.tau2.ci)
  #
  chklogical(print.tau)
  #
  chklogical(print.tau.ci)
  print.tau2.tau <- print.tau2 | print.tau
  if (print.tau2 & print.tau)
    print.tau2 <- FALSE
  if (print.tau2.ci & print.tau.ci)
    print.tau2.ci <- FALSE
  #
  chklogical(print.Q)
  #
  if (is.null(print.pval.Q)) {
    if (is.character(hetstat) || hetstat || overall.hetstat)
      print.pval.Q <- TRUE
    else
      print.pval.Q <- FALSE
  }
  else
    chklogical(print.pval.Q)
  #
  chklogical(print.Rb)
  chklogical(print.Rb.ci)
  #
  if (bmj & overall.hetstat) {
    if (!missing.print.I2 && !print.I2)
      warning("Heterogeneity statistic I2 printed for BMJ layout.")
    if (!missing.print.Q && !print.Q)
      warning("Heterogeneity statistic Q printed for BMJ layout.")
    if (!missing.print.pval.Q && !print.pval.Q)
      warning("P-value of test for heterogeneity printed for BMJ layout.")
    if (!missing.print.Rb && print.Rb)
      warning("Heterogeneity statistic Rb not printed for BMJ layout.")
  }
  #
  if (jama & overall.hetstat) {
    if (!missing.print.I2 && !print.I2)
      warning("Heterogeneity statistic I2 printed for JAMA layout.")
    if (!missing.print.Q && !print.Q)
      warning("Heterogeneity statistic Q printed for JAMA layout.")
    if (!missing.print.pval.Q && !print.pval.Q)
      warning("P-value of test for heterogeneity printed for JAMA layout.")
    if (!missing.print.Rb && print.Rb)
      warning("Heterogeneity statistic Rb not printed for JAMA layout.")
  }
  #
  if (revman5 && overall.hetstat) {
    if ((!missing.print.tau2 & !print.tau2 & !print.tau) |
        (!missing.print.tau  & !print.tau2 & !print.tau))
      warning(paste("Information on between-study variance printed for",
                    "RevMan5 layout"))
    if (!missing.print.I2 && !print.I2)
      warning("Heterogeneity statistic I2 printed for RevMan5 layout.")
    if (!missing.print.Q && !print.Q)
      warning("Heterogeneity statistic Q printed for RevMan5 layout.")
    if (!missing.print.pval.Q && !print.pval.Q)
      warning("P-value of test for heterogeneity printed for RevMan5 layout.")
    if (!missing.print.Rb && print.Rb)
      warning("Heterogeneity statistic Rb not printed for RevMan5 layout.")
  }    
  #
  if (!is.logical(text.subgroup.nohet))
    chkchar(text.subgroup.nohet)
  else if (text.subgroup.nohet)
    text.subgroup.nohet <- "not applicable"
  #
  hetstat.pooled <- ""
  if (is.character(hetstat)) {
    hetstat <- setchar(hetstat, c("common", "random", "study", "fixed"))
    hetstat[hetstat == "fixed"] <- "common"
    hetstat.pooled <- hetstat
    if (!metabind & hetstat.pooled == "study") {
      warning("Argument 'hetstat = \"study\"' ",
              "only considered for 'metabind' objects.")
      hetstat <- print.I2 | print.tau2.tau | print.Q | print.pval.Q | print.Rb
    }
  }
  else
    chklogical(hetstat)
  #
  if (hetstat.pooled == "common") {
    common <- TRUE
    overall <- TRUE
  }
  if (hetstat.pooled == "random") {
    if (length(x$TE.random) == 1) {
      random <- TRUE
      overall <- TRUE
    }
    else
      hetstat.pooled <- ""
  }
  #
  if (missing.overall.hetstat) {
    if (!missing.hetstat) {
      if (is.character(hetstat))
        overall.hetstat <- FALSE
      else
        overall.hetstat <- hetstat
    }
  }
  else
    chklogical(overall.hetstat)
  #
  chklogical(LRT)
  if (LRT & all(x$method != "GLMM")) {
    warning("Likelihood-Ratio test of heterogeneity only ",
            "available for generalised linear mixed models.")
    LRT <- FALSE
  }
  #
  chkchar(hetlab)
  if (!is.null(resid.hetstat))
    chklogical(resid.hetstat)
  else {
    if (overall && (is.character(hetstat) || hetstat) && !LRT &&
        !is.null(x$tau.common) && x$tau.common)
      resid.hetstat <- TRUE
    else
      resid.hetstat <- FALSE
  }
  chkchar(resid.hetlab)
  #
  test.overall.common <-
    deprecated(test.overall.common, missing.test.overall.common,
               args, "test.overall.fixed",
               warn.deprecated)
  chklogical(test.overall.common)
  chklogical(test.overall.random)
  #
  test.subgroup <- replaceNULL(test.subgroup, TRUE)
  chklogical(test.subgroup)
  test.subgroup.common <-
    deprecated(test.subgroup.common, missing.test.subgroup.common,
               args, "test.subgroup.fixed",
               warn.deprecated)
  chklogical(test.subgroup.common)
  chklogical(test.subgroup.random)
  #
  chklogical(print.Q.subgroup)
  #
  if (missing.header.line) {
    if (is.character(gs("header.line")))
      header.line <- gs("header.line")
    else
      header.line <- gs("header.line") | bmj.revman5.jama
  }
  #
  if (is.character(header.line)) {
    header.line.pos <- setchar(header.line, c("below", "both", ""))
    header.line <- header.line.pos != ""
  }
  else {
    chklogical(header.line)
    if (header.line)
      header.line.pos <- "below"
    else
      header.line.pos <- ""
  }
  #
  chknumeric(fontsize, length = 1)
  #
  if (!is.null(fontfamily)) {
    chkchar(fontfamily, length = 1)
    monospaced <- fontfamily %in% "Courier New"
  }
  else
    monospaced <- FALSE
  #
  chknumeric(fs.heading, length = 1)
  #
  if (!missing.fs.fixed) {
    fs.common <-
      deprecated(fs.common, missing.fs.common, args, "fs.fixed",
                 warn.deprecated)
    missing.fs.common <- FALSE
  }
  #
  if (!missing.fs.common)
    chknumeric(fs.common, length = 1)
  #
  if (!is.null(fs.random))
    chknumeric(fs.random, length = 1)
  if (!is.null(fs.predict))
    chknumeric(fs.predict, length = 1)
  #
  if (!missing.fs.fixed.labels) {
    fs.common.labels <-
      deprecated(fs.common.labels, missing.fs.common.labels,
                 args, "fs.fixed.labels",
                 warn.deprecated)
    missing.fs.common.labels <- FALSE
  }
  if (!missing.fs.common.labels)
    chknumeric(fs.common.labels, length = 1)
  #
  if (!is.null(fs.random.labels))
    chknumeric(fs.random.labels, length = 1)
  if (!is.null(fs.predict.labels))
    chknumeric(fs.predict.labels, length = 1)
  if (!is.null(fs.study))
    chknumeric(fs.study, length = 1)
  if (!is.null(fs.study.labels))
    chknumeric(fs.study.labels, length = 1)
  if (!is.null(fs.hetstat))
    chknumeric(fs.hetstat, length = 1)
  if (!is.null(fs.test.overall))
    chknumeric(fs.test.overall, length = 1)
  if (!is.null(fs.test.subgroup))
    chknumeric(fs.test.subgroup, length = 1)
  if (!is.null(fs.test.effect.subgroup))
    chknumeric(fs.test.effect.subgroup, length = 1)
  if (!is.null(fs.addline))
    chknumeric(fs.addline, length = 1)
  chknumeric(fs.axis, length = 1)
  chknumeric(fs.smlab, length = 1)
  chknumeric(fs.xlab, length = 1)
  chknumeric(fs.lr, length = 1)
  #
  if (!missing.ff.fixed) {
    ff.common <-
      deprecated(ff.common, missing.ff.common, args, "ff.fixed",
                 warn.deprecated)
    missing.ff.common <- FALSE
  }
  #
  if (!missing.ff.fixed.labels) {
    ff.common.labels <-
      deprecated(ff.common.labels, missing.ff.common.labels,
                 args, "ff.fixed.labels",
                 warn.deprecated)
    missing.ff.common.labels <- FALSE
  }
  #
  chknumeric(squaresize, length = 1)
  chknumeric(lwd, length = 1)
  chknumeric(lwd.square, length = 1)
  chknumeric(lwd.diamond, length = 1)
  #
  arrow.type <- setchar(arrow.type, c("open", "closed"))
  chknumeric(arrow.length, min = 0, zero = FALSE)
  #
  chklogical(calcwidth.pooled)
  calcwidth.common <-
    deprecated(calcwidth.common, missing.calcwidth.common,
               args, "calcwidth.fixed",
               warn.deprecated)
  chklogical(calcwidth.common)
  chklogical(calcwidth.random)
  chklogical(calcwidth.predict)
  chklogical(calcwidth.hetstat)
  chklogical(calcwidth.tests)
  chklogical(calcwidth.subgroup)
  chklogical(calcwidth.addline)
  #
  if (missing.just && bmj)
    just <- "center"
  just.cols <- setchar(just, c("right", "center", "left"))
  just.studlab <- setchar(just.studlab, c("right", "center", "left"))
  just.addcols <- setchar(just.addcols, c("right", "center", "left"))
  just.addcols.left <- setchar(just.addcols.left, c("right", "center", "left"))
  just.addcols.right <-
    setchar(just.addcols.right, c("right", "center", "left"))
  just.label.e <- setchar(just.label.e, c("right", "center", "left"))
  just.label.c <- setchar(just.label.c, c("right", "center", "left"))
  #
  if (!is.null(bmj.text))
    chkchar(bmj.text, length = 1)
  chknumeric(bmj.xpos)
  chkchar(bmj.sep, length = 1)
  #
  chknumeric(spacing, length = 1)
  #
  if (!missing.text.addline1) {
    chkchar(text.addline1)
    if (text.addline1 == "")
      missing.text.addline1 <- TRUE
  }
  else
    text.addline1 <- ""
  #
  if (!missing.text.addline2) {
    chkchar(text.addline2)
    if (text.addline2 == "")
      missing.text.addline2 <- TRUE
  }
  else
    text.addline2 <- ""
  #
  # Check and set additional empty rows in forest plot
  #
  if (!is.null(addrow))
    chklogical(addrow)
  else
    addrow <- !jama
  if (!is.null(addrow.overall)) {
    chklogical(addrow.overall)
    if (!(overall & (common | random | prediction)))
      addrow.overall <- FALSE
  }
  else
    addrow.overall <- !jama & overall & (common | random | prediction)
  if (!is.null(addrow.subgroups))
    chklogical(addrow.subgroups)
  else
    addrow.subgroups <- !jama
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  chknumeric(digits.pval.Q, min = 1, length = 1)
  chknumeric(digits.Q, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  if (!missing.digits.mean)
    chknumeric(digits.mean, min = 0, length = 1)
  if (!missing.digits.sd)
    chknumeric(digits.sd, min = 0, length = 1)
  if (!missing.digits.cor)
    chknumeric(digits.cor, min = 0, length = 1)
  if (!missing.digits.time)
    chknumeric(digits.time, min = 0, length = 1)
  chknumeric(digits.n, min = 0, length = 1)
  chknumeric(digits.event, min = 0, length = 1)
  if (!missing.digits.TE)
    chknumeric(digits.TE, min = 0, length = 1)
  if (!missing.digits.addcols)
    chknumeric(digits.addcols, min = 0)
  if (!missing.digits.addcols.right)
    chknumeric(digits.addcols.right, min = 0)
  if (!missing.digits.addcols.left)
    chknumeric(digits.addcols.left, min = 0)
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  #
  # Check for deprecated arguments in '...'
  #
  weight.study <-
    deprecated(weight.study, missing.weight.study, args, "weight",
               warn.deprecated)
  if (missing.weight.study || is.null(weight.study))
    weight.study <- ifelse(random & !common, "random", "common")
  weight.study <- setchar(weight.study, c("same", "common", "random", "fixed"))
  weight.study[weight.study == "fixed"] <- "common"
  #
  digits.stat <-
    deprecated(digits.stat, missing.digits.stat, args, "digits.zval",
               warn.deprecated)
  chknumeric(digits.stat, min = 0, length = 1)
  #
  print.stat <-
    deprecated(print.stat, missing.print.stat, args, "print.zval",
               warn.deprecated)
  chklogical(print.stat)
  #
  label <- deprecated(label, missing.label, args, "labels",
                      warn.deprecated)
  #
  label.e <- deprecated(label.e, missing.label.e, args, "lab.e",
                        warn.deprecated)
  label.c <- deprecated(label.c, missing.label.c, args, "lab.c",
                        warn.deprecated)
  #
  label.e.attach <- deprecated(label.e.attach, missing.label.e.attach,
                               args, "lab.e.attach.to.col",
                               warn.deprecated)
  label.c.attach <- deprecated(label.c.attach, missing.label.c.attach,
                               args, "lab.c.attach.to.col",
                               warn.deprecated)
  #
  col.inside <-
    deprecated(col.inside, missing.col.inside, args, "col.i.inside.square",
               warn.deprecated)
  chkcolor(col.inside)
  chkcolor(col.inside.common)
  chkcolor(col.inside.random)
  #
  if (!missing.col.diamond.lines.fixed) {
    col.diamond.lines.common <-
      deprecated(col.diamond.lines.common, missing.col.diamond.lines.common,
                 args, "col.diamond.lines.fixed",
                 warn.deprecated)
  }
  else if (!missing.col.diamond.fixed.lines) {
    col.diamond.lines.common <-
      deprecated(col.diamond.lines.common, missing.col.diamond.lines.common,
                 args, "col.diamond.fixed.lines",
                 warn.deprecated)
  }
  chkcolor(col.diamond.lines)
  chkcolor(col.diamond.lines.common)
  #
  col.diamond.lines.random <-
    deprecated(col.diamond.lines.random, missing.col.diamond.lines.random,
               args, "col.diamond.random.lines",
               warn.deprecated)
  chkcolor(col.diamond.lines.random)
  #
  # Additional assignments
  #
  if (bmj) {
    if (missing.ff.common)
      ff.common <- "plain"
    if (missing.ff.random)
      ff.random <- ff.common
    if (missing.ff.predict)
      ff.predict <- ff.common
    if (missing.ff.common.labels)
      ff.common.labels <- ff.common
    if (missing.ff.random.labels)
      ff.random.labels <- ff.random
    if (missing.ff.predict.labels)
      ff.predict.labels <- ff.predict
    #
    if (missing.fs.common)
      fs.common <- fontsize
    if (missing.fs.random)
      fs.random <- fs.common
    if (missing.fs.predict)
      fs.predict <- fs.common
    if (missing.fs.common.labels)
      fs.common.labels <- fs.common
    if (missing.fs.random.labels)
      fs.random.labels <- fs.random
    if (missing.fs.predict.labels)
      fs.predict.labels <- fs.predict
  }
  else if (jama) {
    if (missing.ff.common)
      ff.common <- "plain"
    if (missing.ff.random)
      ff.random <- ff.common
    if (missing.ff.predict)
      ff.predict <- ff.common
    if (missing.ff.common.labels)
      ff.common.labels <- ff.common
    if (missing.ff.random.labels)
      ff.random.labels <- ff.random
    if (missing.ff.predict.labels)
      ff.predict.labels <- ff.predict
    #
    if (missing.fs.common)
      fs.common <- fontsize
    if (missing.fs.random)
      fs.random <- fs.common
    if (missing.fs.predict)
      fs.predict <- fs.common
    if (missing.fs.common.labels)
      fs.common.labels <- fs.common
    if (missing.fs.random.labels)
      fs.random.labels <- fs.random
    if (missing.fs.predict.labels)
      fs.predict.labels <- fs.predict
  }
  else {
    if (missing.ff.common)
      ff.common <- "bold"
    if (missing.ff.random)
      ff.random <- ff.common
    if (missing.ff.predict)
      ff.predict <- ff.common
    if (missing.ff.common.labels)
      ff.common.labels <- ff.common
    if (missing.ff.random.labels)
      ff.random.labels <- ff.random
    if (missing.ff.predict.labels)
      ff.predict.labels <- ff.predict
    #
    if (missing.fs.common)
      fs.common <- fontsize
    if (missing.fs.random)
      fs.random <- fs.common
    if (missing.fs.predict)
      fs.predict <- fs.common
    if (missing.fs.common.labels)
      fs.common.labels <- fs.common
    if (missing.fs.random.labels)
      fs.random.labels <- fs.random
    if (missing.fs.predict.labels)
      fs.predict.labels <- fs.predict
  }
  hetseparator <- " = "
  #
  if (bmj) {
    if (missing.ff.hetstat)
      ff.hetstat <- "plain"
    if (missing.ff.test.overall)
      ff.test.overall <- ff.hetstat
    if (missing.ff.test.subgroup)
      ff.test.subgroup <- ff.hetstat
    if (missing.ff.test.effect.subgroup)
      ff.test.effect.subgroup <- ff.hetstat
    if (missing.ff.addline)
      ff.addline <- ff.hetstat
    #
    if (missing.fs.hetstat)
      fs.hetstat <- fontsize - 1
    if (missing.fs.test.overall)
      fs.test.overall <- fs.hetstat
    if (missing.fs.test.subgroup)
      fs.test.subgroup <- fs.hetstat
    if (missing.fs.test.effect.subgroup)
      fs.test.effect.subgroup <- fs.hetstat
    if (missing.fs.addline)
      fs.addline <- fs.hetstat
  }
  else if (revman5) {
    if (missing.ff.hetstat)
      ff.hetstat <- "plain"
    if (missing.ff.test.overall)
      ff.test.overall <- ff.hetstat
    if (missing.ff.test.subgroup)
      ff.test.subgroup <- ff.hetstat
    if (missing.ff.test.effect.subgroup)
      ff.test.effect.subgroup <- ff.hetstat
    if (missing.ff.addline)
      ff.addline <- ff.hetstat
    #
    if (missing.fs.hetstat)
      fs.hetstat <- fontsize - 1
    if (missing.fs.test.overall)
      fs.test.overall <- fs.hetstat
    if (missing.fs.test.subgroup)
      fs.test.subgroup <- fs.hetstat
    if (missing.fs.test.effect.subgroup)
      fs.test.effect.subgroup <- fs.hetstat
    if (missing.fs.addline)
      fs.addline <- fs.hetstat
  }
  else if (jama) {
    if (missing.ff.hetstat)
      ff.hetstat <- "plain"
    if (missing.ff.test.overall)
      ff.test.overall <- ff.hetstat
    if (missing.ff.test.subgroup)
      ff.test.subgroup <- ff.hetstat
    if (missing.ff.test.effect.subgroup)
      ff.test.effect.subgroup <- ff.hetstat
    if (missing.ff.addline)
      ff.addline <- ff.hetstat
    #
    if (missing.fs.hetstat)
      fs.hetstat <- fontsize - 1
    if (missing.fs.test.overall)
      fs.test.overall <- fs.hetstat
    if (missing.fs.test.subgroup)
      fs.test.subgroup <- fs.hetstat
    if (missing.fs.test.effect.subgroup)
      fs.test.effect.subgroup <- fs.hetstat
    if (missing.fs.addline)
      fs.addline <- fs.hetstat
  }
  else {
    if (missing.ff.hetstat)
      ff.hetstat <- "plain"
    if (missing.ff.test.overall)
      ff.test.overall <- ff.hetstat
    if (missing.ff.test.subgroup)
      ff.test.subgroup <- ff.hetstat
    if (missing.ff.test.effect.subgroup)
      ff.test.effect.subgroup <- ff.hetstat
    if (missing.ff.addline)
      ff.addline <- ff.hetstat
    #
    if (missing.fs.hetstat)
      fs.hetstat <- fontsize - 1
    if (missing.fs.test.overall)
      fs.test.overall <- fs.hetstat
    if (missing.fs.test.subgroup)
      fs.test.subgroup <- fs.hetstat
    if (missing.fs.test.effect.subgroup)
      fs.test.effect.subgroup <- fs.hetstat
    if (missing.fs.addline)
      fs.addline <- fs.hetstat
  }
  #
  chklogical(details)
  #
  # Set colours for JAMA and RevMan5 layouts
  #
  if (jama) {
    if (missing.col.square)
      col.square <- rep("darkblue", K.all)
    if (missing.col.square.lines)
      col.square.lines <- rep("darkblue", K.all)
    #
    if (missing.col.circle)
      col.circle <- rep("darkblue", K.all)
    if (missing.col.circle.lines)
      col.circle.lines <- rep("darkblue", K.all)
    #
    if (missing.col.diamond.common)
      col.diamond.common <- "lightblue"
    if (missing.col.diamond.random)
      col.diamond.random <- "lightblue"
  }
  else if (revman5) {
    if (missing.col.square) {
      if (metacont | metamean)
        col.square <- rep("green", K.all)
      else if (metabin)
        col.square <- rep("blue", K.all)
      else
        col.square <- rep("red", K.all)
    }
    if (missing.col.square.lines) {
      if (metacont | metamean)
        col.square.lines <- rep("green", K.all)
      else if (metabin)
        col.square.lines <- rep("darkblue", K.all)
      else
        col.square.lines <- rep("red", K.all)
    }
    #
    if (missing.col.circle) {
      if (metacont | metamean)
        col.circle <- rep("green", K.all)
      else if (metabin)
        col.circle <- rep("blue", K.all)
      else
        col.circle <- rep("red", K.all)
    }
    if (missing.col.circle.lines) {
      if (metacont | metamean)
        col.circle.lines <- rep("green", K.all)
      else if (metabin)
        col.circle.lines <- rep("darkblue", K.all)
      else
        col.circle.lines <- rep("red", K.all)
    }
    #
    if (missing.col.diamond.common)
      col.diamond.common <- "black"
    if (missing.col.diamond.random)
      col.diamond.random <- "black"
  }
  
  
  #
  #
  # (4) Check length of variables
  #
  #
  fun <- "forest.meta"
  #
  if (length(col.study) == 1)
    col.study <- rep(col.study, K.all)
  else
    chklength(col.study, K.all, fun)
  #
  if (length(col.inside) == 1)
    col.inside <- rep(col.inside, K.all)
  else
    chklength(col.inside, K.all, fun)
  #
  if (length(col.square) == 1)
    col.square <- rep(col.square, K.all)
  else
    chklength(col.square, K.all, fun)
  #
  if (length(col.square.lines) == 1)
    col.square.lines <- rep(col.square.lines, K.all)
  else
    chklength(col.square.lines, K.all, fun)
  #
  if (length(col.circle) == 1)
    col.circle <- rep(col.circle, K.all)
  else
    chklength(col.circle, K.all, fun)
  #
  if (length(col.circle.lines) == 1)
    col.circle.lines <- rep(col.circle.lines, K.all)
  else
    chklength(col.circle.lines, K.all, fun)
  
  
  #
  #
  # (5) Some assignments and additional checks
  #
  #
  n.com <- max(length(x$TE.common), 1)
  n.ran <- max(length(x$lower.random), 1)
  n.prd <- max(length(x$lower.predict), 1)
  #
  prediction <- prediction &
    any(!is.na(x$lower.predict) & !is.na(x$upper.predict))
  #
  if (by) {
    if (!missing.subgroup & !metabind)
      subgroup <- catch("subgroup", mc, x, sfsp)
    #
    if (length(subgroup) == 1)
      subgroup.logical <- rep(subgroup, n.by) &
        (x$k.w > 1 | layout == "subgroup")
    else {
      chklength(subgroup, n.by,
                text = paste("Length of argument 'subgroup' must be",
                             "equal to 1 or number of subgroups."))
      subgroup.logical <- subgroup
    }
    chklogical(subgroup[1])
    #
    if (!missing.subgroup.hetstat & !metabind)
      subgroup.hetstat <- catch("subgroup.hetstat", mc, x, sfsp)
    #
    if (length(subgroup.hetstat) == 1 & is.character(subgroup.hetstat))
      subgroup.hetstat.logical <- rep(TRUE, n.by)
    else if (length(subgroup.hetstat) == 1) {
      chklogical(subgroup.hetstat)
      subgroup.hetstat.logical <- subgroup.hetstat &
        (x$k.w > 1 | layout == "subgroup")
    }
    else {
      chklength(subgroup.hetstat, n.by,
                text = paste("Length of argument 'subgroup.hetstat' must be",
                             "equal to 1 or number of subgroups."))
      chklogical(subgroup.hetstat[1])
      subgroup.hetstat.logical <- subgroup.hetstat
    }
    #
    if (!missing.common.subgroup)
      common.subgroup <- catch("common.subgroup", mc, x, sfsp)
    #
    if (!missing.random.subgroup)
      random.subgroup <- catch("random.subgroup", mc, x, sfsp)
    #
    if (!missing.prediction.subgroup & !metabind)
      prediction.subgroup <- catch("prediction.subgroup", mc, x, sfsp)
    #
    common.subgroup <- replaceNULL(common.subgroup, FALSE)
    random.subgroup <- replaceNULL(random.subgroup, FALSE)
    prediction.subgroup <- replaceNULL(prediction.subgroup, FALSE)
    #
    # if (length(prediction.subgroup) == 1) {
    #   if (is.matrix(x$lower.predict.w))
    #     prediction.subgroup.logical <-
    #       prediction.subgroup &
    #       apply(x$lower.predict.w, 1, notallNA) &
    #       apply(x$upper.predict.w, 1, notallNA)
    #   else {
    #     prediction.subgroup.logical <-
    #       prediction.subgroup &
    #       notallNA(x$lower.predict.w) &
    #       notallNA(x$upper.predict.w)
    #     prediction.subgroup.logical <-
    #       rep(prediction.subgroup.logical, n.by)
    #   }
    # }
    # else {
    #   chklength(prediction.subgroup, n.by,
    #             text = paste("Length of argument 'prediction.subgroup' must be",
    #                          "equal to 1 or number of subgroups."))
    #   prediction.subgroup.logical <-
    #     prediction.subgroup
    # }
    #
    common.subgroup.logical <-
      show_subgroup_results(common.subgroup, n.by,
                            x$lower.common.w, x$upper.common.w)
    #
    random.subgroup.logical <-
      show_subgroup_results(random.subgroup, n.by,
                            x$lower.random.w, x$upper.random.w)
    #
    prediction.subgroup.logical <-
      show_subgroup_results(prediction.subgroup, n.by,
                            x$lower.predict.w, x$upper.predict.w)
    #
    chklogical(common.subgroup[1])
    chklogical(random.subgroup[1])
    chklogical(prediction.subgroup[1])
    #
    if (!missing.test.effect.subgroup) {
      test.effect.subgroup <-
        catch("test.effect.subgroup", mc, x, sfsp)
      test.effect.subgroup <- replaceNULL(test.effect.subgroup, FALSE)
      #
      if (length(test.effect.subgroup) == 1) {
        chklogical(test.effect.subgroup)
        #
        test.effect.subgroup.logical <-
          rep(test.effect.subgroup, n.by) &
          (x$k.w > 1 | layout == "subgroup")
      }
      else {
        chklength(test.effect.subgroup, n.by,
                  text = paste("Length of argument 'test.effect.subgroup'",
                               "must be equal to 1 or number of subgroups."))
        chklogical(test.effect.subgroup[1])
        test.effect.subgroup.logical <- test.effect.subgroup
      }
    }
    else
      test.effect.subgroup.logical <-
        rep(test.effect.subgroup, n.by) & subgroup.logical
    #
    if (missing.test.effect.subgroup.common)
      test.effect.subgroup.common.logical <-
        common & test.effect.subgroup.logical
    else {
      test.effect.subgroup.common <-
        catch("test.effect.subgroup.common", mc, x, sfsp)
      test.effect.subgroup.common <-
        replaceNULL(test.effect.subgroup.common, FALSE)
      #
      if (length(test.effect.subgroup.common) == 1) {
        chklogical(test.effect.subgroup.common)
        #
        test.effect.subgroup.common.logical <-
          rep(test.effect.subgroup.common, n.by) &
          (x$k.w > 1 | layout == "subgroup")
      }
      else {
        chklength(test.effect.subgroup.common, n.by,
                  text = paste("Length of argument",
                               "'test.effect.subgroup.common'",
                               "must be equal to 1 or number of subgroups."))
        chklogical(test.effect.subgroup.common[1])
        test.effect.subgroup.common.logical <- test.effect.subgroup.common
      }
    }
    #
    if (missing.test.effect.subgroup.random)
      test.effect.subgroup.random.logical <-
        random & test.effect.subgroup.logical
    else {
      test.effect.subgroup.random <-
        catch("test.effect.subgroup.random", mc, x, sfsp)
      test.effect.subgroup.random <-
        replaceNULL(test.effect.subgroup.random, FALSE)
      #
      if (length(test.effect.subgroup.random) == 1) {
        chklogical(test.effect.subgroup.random)
        #
        test.effect.subgroup.random.logical <-
          rep(test.effect.subgroup.random, n.by) &
          (x$k.w > 1 | layout == "subgroup")
      }
      else {
        chklength(test.effect.subgroup.random, n.by,
                  text = paste("Length of argument",
                               "'test.effect.subgroup.random'",
                               "must be equal to 1 or number of subgroups."))
        chklogical(test.effect.subgroup.random[1])
        test.effect.subgroup.random.logical <- test.effect.subgroup.random
      }
    }
  }
  #
  level <- x$level
  level.ma <- x$level.ma
  level.predict <- x$level.predict
  #
  if (is.null(label.right))
    label.right <- ""
  if (is.null(label.left))
    label.left <- ""
  #
  subgroup <- x$subgroup
  #
  three.level <- !is.null(x$three.level) && any(x$three.level)
  n_of_1 <- !is.null(x$cycles)
  #
  if (!by) {
    common.random <- common & random
    #
    test.subgroup.common <- FALSE
    test.subgroup.random <- FALSE
  }
  else
    common.random <-
      any(common | test.subgroup.common | test.effect.subgroup.common.logical) &
      any(random | test.subgroup.random | test.effect.subgroup.random.logical)
  #
  if (layout == "subgroup") {
    if (!missing.study.results & study.results)
      warning("Argument 'study.results' set to FALSE as ",
              "argument 'layout' is \"subgroup\".")
    study.results <- FALSE
  }
  #
  if (!missing.text.fixed) {
    text.common <-
      deprecated(text.common, missing.text.common, args, "text.fixed",
                 warn.deprecated)
    missing.text.common <- FALSE
  }
  #
  if (missing.text.common | is.null(text.common)) {
    if (is.null(text.common) || length(text.common) == 1) {
      if (study.results & (x$level != x$level.ma | bmj.revman5)) {
        if (revman5.jama)
          text.common <- paste0("Total (",
                                if (common.random)
                                  paste0(gs("text.w.common"), " effect, "),
                                round(x$level.ma * 100), "% CI)")
        else if (bmj)
          text.common <- paste0("Total (",
                                round(x$level.ma * 100), "% CI)",
                                if (common.random)
                                  paste0(", ", gs("text.w.common"))
                                )
        else if (!is.null(text.common))
          text.common <- paste0(text.common, " (",
                                round(x$level.ma * 100), "%-CI)")
        else
          text.common <- paste0(gs("text.common"), " (",
                                round(x$level.ma * 100), "%-CI)")
      }
      else {
        if (bmj.revman5.jama) {
          text.common <- "Total"
          if (common.random)
            text.common <- paste0(text.common, " (", gs("text.w.common"),
                                  " effect)")
        }
        else if (is.null(text.common))
          text.common <- gs("text.common")
      }
    }
  }
  #
  if (missing.text.random | is.null(text.random)) {
    if (is.null(text.random) || length(text.random) == 1) {
      if (study.results & (x$level != x$level.ma | bmj.revman5)) {
        if (bmj.revman5.jama) {
          if (revman5.jama) {
            text.random <- paste0("Total (",
                                  if (common.random)
                                    paste0(gs("text.w.random"), " effect, "))
            #
            if (length(x$lower.random) > 1) {
              meth.r <- gsub("classic", "", x$method.random.ci)
              text.random <-
                paste0(text.random, meth.r,
                       ifelse(meth.r == "HK" & x$adhoc.hakn.ci != "",
                              paste0("-",
                                     toupper(substring(x$adhoc.hakn.ci, 1, 2))),
                              ""),
                       ifelse(meth.r != "", ", ", ""))
            }
            text.random <- paste0(text.random,
                                  round(x$level.ma * 100), "% CI)")
          }
          else if (bmj)
            text.random <- paste0("Total (",
                                  round(x$level.ma * 100), "% CI)",
                                  if (common.random)
                                    paste0(", ", gs("text.w.random"))
                                  )
        }
        else if (!is.null(text.random))
          text.random <- paste0(text.random, " (",
                                round(x$level.ma * 100), "%-CI)")
        else
          text.random <- paste0(gs("text.random"), " (",
                                round(x$level.ma * 100), "%-CI)")
      }
      else {
        if (bmj.revman5.jama) {
          text.random <- "Total"
          if (common.random || length(x$lower.random) > 1) {
            text.random <- paste0(text.random, " (", gs("text.w.random"),
                                  " effect")
            if (length(x$lower.random) > 1) {
              meth.r <- gsub("classic", "", x$method.random.ci)
              text.random <-
                paste0(text.random,
                       ifelse(meth.r != "", ", ", ""),
                       meth.r,
                       ifelse(meth.r == "HK" & x$adhoc.hakn.ci != "",
                              paste0("-",
                                     toupper(substring(x$adhoc.hakn.ci, 1, 2))),
                              ""))
            }
            text.random <- paste0(text.random, ")")
          }
        }
        else if (is.null(text.random)) {
          text.random <- gs("text.random")
          if (length(x$lower.random) > 1) {
            meth.r <- gsub("classic", "", x$method.random.ci)
            text.random <-
              paste0(text.random,
                     ifelse(meth.r != "", "(", ""),
                     meth.r,
                     ifelse(meth.r == "HK" & x$adhoc.hakn.ci != "",
                            paste0("-",
                                   toupper(substring(x$adhoc.hakn.ci, 1, 2))),
                            ""),
                     ifelse(meth.r != "", ")", ""))
          }
        }
      }
    }
  }
  #
  if (length(text.random) == 1 & n.ran > 1)
    text.random <- rep(text.random, n.ran)
  #
  if (missing.text.predict | is.null(text.predict)) {
    if (is.null(text.predict))
      text.predict <- rep("Prediction interval", n.prd)
    if (!(length(x$level.predict) == 0) &&
        (study.results & (x$level != x$level.predict |
                          x$level.ma != x$level.predict)))
      text.predict <- paste0(text.predict, " (",
                             round(x$level.predict * 100), "%-PI)")
  }
  #
  if (is.null(x$null.effect) || is.na(x$null.effect)) {
    test.overall.common <- FALSE
    test.overall.random <- FALSE
  }
  #
  if (!overall) {
    if (test.overall)
      test.overall <- FALSE
    if (test.overall.common)
      test.overall.common <- FALSE
    if (test.overall.random)
      test.overall.random <- FALSE
  }
  #
  # Add space for heterogeneity statistics (if needed)
  #
  if (is.null(addrows.below.overall)) {
    addrows.below.overall <- 0
    #
    if (jama)
      addrows.below.overall <- 3
    else if (layout == "meta" & (metacor | metagen) &
             missing.leftcols &
             (overall.hetstat | test.overall.common |
              test.overall.random) &
             !(calcwidth.hetstat | calcwidth.tests) &
             !(details | RoB.legend))
      addrows.below.overall <- 2
    else if (layout %in% c("meta", "subgroup") &
             (test.subgroup.common & test.subgroup.random) &
             !calcwidth.tests &
             !(details | RoB.legend)) {
      addrows.below.overall <- 1
      #
      if (layout == "meta")
        addrows.below.overall <-
          addrows.below.overall +
          as.numeric((label.left != "" | label.right != "")) +
          as.numeric(xlab != "")
    }
  }
  #
  chknumeric(addrows.below.overall, min = 0, length = 1, integer = TRUE)
  #
  if (!avail.xlim) {
    mrm <- c("metaprop", "metarate", "metamean")
    #
    if (metaprop || metarate || metamean ||
        (metabind && sm %in% c(gs("sm4prop"), gs("sm4rate"), gs("sm4mean"))) ||
        (metacum && any(x$classes %in% mrm)) ||
        (metainf && any(x$classes %in% mrm))) {
      xlim <- NULL
      #
      avail.xlim <- FALSE
    }
    else
      xlim <- "symmetric"
  }
  #
  if (just.studlab == "left")
    xpos.s <- 0
  else if (just.studlab == "center")
    xpos.s <- 0.5
  else if (just.studlab == "right")
    xpos.s <- 1
  #
  if (just.cols == "left")
    xpos.c <- 0
  else if (just.cols == "center")
    xpos.c <- 0.5
  else if (just.cols == "right")
    xpos.c <- 1
  #
  if (just.label.e == "left")
    xpos.label.e <- 0
  else if (just.label.e == "center")
    xpos.label.e <- 0.5
  else if (just.label.e == "right")
    xpos.label.e <- 1
  #
  if (just.label.c == "left")
    xpos.label.c <- 0
  else if (just.label.c == "center")
    xpos.label.c <- 0.5
  else if (just.label.c == "right")
    xpos.label.c <- 1
  #
  if (log.xaxis) {
    ref <- log(ref)
    cid.below.null <- log(cid.below.null)
    cid.above.null <- log(cid.above.null)
  }
  #
  if (!backtransf & !missing.pscale & pscale != 1 & !is_untransformed(sm)) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  #
  if (!backtransf & pscale != 1)
    pscale <- 1
  #
  if (!backtransf & !missing.irscale & irscale != 1 & !is_untransformed(sm)) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  #
  if (!backtransf & irscale != 1)
    irscale <- 1
  #
  if (is.null(xlab))
    xlab <- xlab_meta(sm, backtransf, newline = revman5.jama, revman5 = revman5,
                 big.mark = big.mark)
  #
  scale <- 1
  if (pscale != 1 || irscale != 1) {
    if (pscale != 1 && irscale != 1)
      stop("Provide either arguments 'pscale' or 'irscale' but not ",
           "both arguments.",
           call. = FALSE)
    if (pscale != 1)
      scale <- pscale
    else
      scale <- irscale
  }
  #
  smlab.null <- is.null(smlab)
  if (smlab.null) {
    if (is_rate(sm))
      smlab <- xlab_meta(sm, backtransf, irscale = irscale, irunit = irunit,
                    newline = !bmj.revman5.jama, revman5 = revman5,
                    big.mark = big.mark)
    else if (is_prop(sm))
      smlab <- xlab_meta(sm, backtransf, pscale = pscale,
                    newline = !bmj.revman5.jama, revman5 = revman5,
                    big.mark = big.mark)
    else
      smlab <- xlab_meta(sm, backtransf, pscale = pscale,
                    irscale = irscale, irunit = irunit,
                    newline = !bmj.revman5.jama, revman5 = revman5,
                    big.mark = big.mark)
  }
  #
  print.label <- (label.left != "" | label.right != "") & !is.na(ref)
  if (print.label & !bottom.lr) {
    if (!smlab.null)
      warning("Argument 'smlab' ignored as argument 'bottom.lr' is FALSE.")
    smlab <- ""
  }
  #
  if (!by)
    addrow.subgroups <- FALSE
  #
  if (resid.hetstat &&
      (!by || (by && is.null(x$tau.common) || !x$tau.common))) {
    if (!missing.resid.hetstat)
      warning("Information on residual heterogeneity only added to ",
              "forest plot of meta-analysis with subgroups ",
              "assuming common estimator for between-study heterogeneity ",
              "(argument 'tau.common = TRUE' in meta-analysis functions)",
              call. = FALSE)
    resid.hetstat <- FALSE
  }
  #
  plotwidth <- setunit(plotwidth)
  colgap <- setunit(colgap)
  colgap.left <- setunit(colgap.left)
  colgap.right <- setunit(colgap.right)
  colgap.studlab <- setunit(colgap.studlab)
  colgap.forest <- setunit(colgap.forest)
  colgap.forest.left <- setunit(colgap.forest.left)
  colgap.forest.right <- setunit(colgap.forest.right)
  colgap.rob <- setunit(colgap.rob)
  colgap.rob.overall <- setunit(colgap.rob.overall)
  #
  if (!missing.label.test.overall.fixed) {
    label.test.overall.common <-
      deprecated(label.test.overall.common, missing.label.test.overall.common,
                 args, "label.test.overall.fixed",
                 warn.deprecated)
    missing.label.test.overall.common <- FALSE
  }
  if (missing.label.test.overall.common)
    label.test.overall.common <-
      paste0("Test for overall effect",
             if (common.random)
               paste0(" (", gs("text.w.common"), " effect)"),
             ": ")
  if (missing.label.test.overall.random)
    label.test.overall.random <-
      paste0("Test for overall effect",
             if (common.random)
               paste0(" (", gs("text.w.random"), " effects)"),
             ": ")
  #
  if (!missing.label.test.subgroup.fixed) {
    label.test.subgroup.common <-
      deprecated(label.test.subgroup.common,
                 missing.label.test.subgroup.common,
                 args, "label.test.subgroup.fixed",
                 warn.deprecated)
    missing.label.test.subgroup.common <- FALSE
  }
  if (missing.label.test.subgroup.common)
    label.test.subgroup.common <-
      paste0("Test for subgroup differences",
             if (common.random)
               paste0(" (", gs("text.w.common"), " effect)"),
             ": ")
  if (missing.label.test.subgroup.random)
    label.test.subgroup.random <-
      paste0("Test for subgroup differences",
             if (common.random)
               paste0(" (", gs("text.w.random"), " effects)"),
             ": ")
  #
  if (!missing.label.test.effect.subgroup.fixed) {
    label.test.effect.subgroup.common <-
      deprecated(label.test.effect.subgroup.common,
                 missing.label.test.effect.subgroup.common,
                 args, "label.test.effect.subgroup.fixed",
                 warn.deprecated)
    missing.label.test.effect.subgroup.common <- FALSE
  }
  if (missing.label.test.effect.subgroup.common)
    label.test.effect.subgroup.common <-
      paste0(if (bmj.revman5.jama)
               "Test for overall effect"
             else
               "Test for effect in subgroup",
             if (common.random)
               paste0(" (", gs("text.w.common"), " effect)"),
             ": ")
  if (missing.label.test.effect.subgroup.random)
    label.test.effect.subgroup.random <-
      paste0(if (bmj.revman5.jama)
               "Test for overall effect"
             else
               "Test for effect in subgroup",
             if (common.random)
               paste0(" (", gs("text.w.random"), " effects)"),
             ": ")
  #
  fs.head <- fs.heading
  ff.head <- ff.heading
  #
  just.c <- just.cols
  just.s <- just.studlab
  
  
  #
  #
  # (6) Labels for columns on left and right side of forest plot
  #
  #
  sm.lab <- sm
  #
  if (backtransf) {
    if (sm == "ZCOR")
      sm.lab <- "COR"
    else if (is_prop(sm)) {
      if (pscale == 1)
        sm.lab <- "Proportion"
      else
        sm.lab <- "Events"
    }
    else if (is_rate(sm)) {
      if (irscale == 1)
        sm.lab <- "Rate"
      else
        sm.lab <- "Events"
    }
    else if (sm == "proportion")
      sm.lab <- "Proportion"
    else if (sm == "MLN")
      sm.lab <- "Mean"
  }
  else if (is_relative_effect(sm) | sm == "VE" |
           (!is.null(fbt) && fbt == "exp"))
    sm.lab <- paste0("log", if (sm == "VE") "VR" else sm)
  #
  sel.studlab <-
    pmatch(layout, c("meta", "BMJ", "RevMan5", "JAMA", "subgroup"))
  lab.studlab <-
    c("Study", "Study or\nsubgroup", "Study",
      "Source", "Subgroup")[sel.studlab]
  if (revman5 & by)
    lab.studlab <- c("Study or\nSubgroup")
  #
  if (bmj.revman5.jama)
    cisep <- " "
  else
    cisep <- "-"
  #
  if (study.results)
    ci.lab <- paste0(100 * level, "%", cisep, "CI")
  else
    ci.lab <- paste0(100 * level.ma, "%", cisep, "CI")
  #
  if (!missing.col.diamond.fixed) {
    col.diamond.common <-
      deprecated(col.diamond.common, missing.col.diamond.common,
                 args, "col.diamond.fixed",
                 warn.deprecated)
    missing.col.diamond.common <- FALSE
  }
  chkcolor(col.diamond.common)
  #
  if (jama) {
    if (missing.ff.lr)
      ff.lr <- "bold"
    #
    if (xlab == "")
      xlab <- paste0(sm.lab, " (", ci.lab, ")")
    #
    smlab <- ""
    bottom.lr <- FALSE
  }
  else {
    if (revman5) {
      sel.method <- pmatch(x$method, c("Inverse", "MH", "Peto", "GLMM"))
      lab.method <- c("IV", "MH", "Peto", "GLMM")[sel.method]
      #
      if (common.random)
        lab.model <- "Fixed + Random, "
      else if (common)
        lab.model <- "Fixed, "
      else if (random)
        lab.model <- "Random, "
      else
        lab.model <- ""
      #
      if (smlab.null) {
        if (smlab != "")
          smlab <- paste0(smlab, "\n")
        smlab <- paste0(smlab, lab.method[1], ", ", lab.model, ci.lab)
      }
    }
    else if (bmj) {
      sel.method <- pmatch(x$method, c("Inverse", "MH", "Peto", "GLMM"))
      lab.method <- c("IV", "MH", "Peto", "GLMM")[sel.method]
      #
      if (common.random)
        lab.model <- "common + random "
      else if (common)
        lab.model <- "common "
      else if (random)
        lab.model <- "random "
      else
        lab.model <- ""
      #
      if (smlab.null) {
        if (smlab != "")
          smlab <- paste0(smlab, ", ")
        smlab <- paste0(smlab, lab.method[1], ",\n", lab.model,
                        "(", ci.lab, ")")
      }
    }
  }
  #
  if (jama | gs("CIbracket") == "(")
    ci.lab.bracket <- paste0("(", ci.lab, ")")
  else if (gs("CIbracket") == "[")
    ci.lab.bracket <- paste0("[", ci.lab, "]")
  else if (gs("CIbracket") == "{")
    ci.lab.bracket <- paste0("{", ci.lab, "}")
  else if (gs("CIbracket") == "")
    ci.lab.bracket <- ci.lab
  #
  if (!common.random) {
    text.w.common <- paste0("Weight", if (bmj) "\n(%)")
    text.w.random <- paste0("Weight", if (bmj) "\n(%)")
  }
  else {
    text.w.common <-
      deprecated(text.w.common, missing.text.w.common, args, "text.w.fixed",
                 warn.deprecated)
    #
    if (bmj) {
      if (is.null(text.w.common))
        text.w.common <- paste0("Weight (%),\n", gs("text.w.common"))
      else
        text.w.common <- paste0("Weight (%),\n", text.w.common)
      #
      if (is.null(text.w.random))
        text.w.random <- paste0("Weight (%),\n", gs("text.w.random"))
      else
        text.w.random <- paste0("Weight (%),\n", text.w.random)
    }
    else {
      if (is.null(text.w.common))
        text.w.common <- paste0("Weight\n(", gs("text.w.common"), ")")
      else
        text.w.common <- paste0("Weight\n(", text.w.common, ")")
      #
      if (is.null(text.w.random))
        text.w.random <- paste0("Weight\n(", gs("text.w.random"), ")")
      else
        text.w.random <- paste0("Weight\n(", text.w.random, ")")
    }
  }
  #
  lab.TE <- sm
  #
  if (is_relative)
    lab.TE <- paste0("log", if (sm == "VE") "VR" else sm)
  else if (!is.null(ftr)) {
    lab.TE <-
      paste0(ftr, "(", sm,
             if (!is.null(atr) && length(names(atr)) >= 1)
               paste0(", ",
                      paste0(names(atr), "=", paste(atr), collapse = ", ")),
             ")")
  }
  else if (sm == "")
    lab.TE <- "TE"
  #
  # Must be of same length as colnames!!!
  #
  labnames <- c(lab.studlab,
                lab.TE,
                if (revman5) "SE" else paste0("SE(", lab.TE, ")"),
                "Cluster", "Cycles",
                #
                "Total", "Total", "Events", "Events",
                label.e, label.c, label.e,
                "Mean", "Mean", "SD", "SD",
                label.e, label.c, label.e,
                #
                "Cor",
                "Time", "Time",
                #
                sm.lab,
                ci.lab,
                if (bmj.revman5 & smlab.null)
                  smlab
                else
                  paste(sm.lab, ci.lab.bracket),
                #
                text.w.common,
                text.w.common,
                text.w.random,
                #
                "P-value")
  #
  if (newcols) {
    #
    if (length(rightcols.new) > 0) {
      if (missing.rightlabs) {
        rightlabs.new <- rightcols.new
        if (RoB.available)
          rightlabs.new[rightlabs.new %in% colnames(rob)] <- rob.labels
        #
        if ((metacor | metaprop | metamean | metarate) &
            any(rightcols.new == "n"))
          rightlabs.new[rightlabs.new == "n"] <- "Total"
        #
        if (metamean & any(rightcols.new == "mean"))
          rightlabs.new[rightlabs.new == "mean"] <- "Mean"
        #
        if (metamean & any(rightcols.new == "sd"))
          rightlabs.new[rightlabs.new == "sd"] <- "SD"
        #
        if (metarate & any(rightcols.new == "time"))
          rightlabs.new[rightlabs.new == "time"] <- "Time"
        #
        if (any(rightcols.new == "pval"))
          rightlabs.new[rightlabs.new == "pval"] <- "P-value"
        #
        if (any(rightcols.new == "tau2"))
          rightlabs.new[rightlabs.new == "tau2"] <- "Tau2"
        #
        if (any(rightcols.new == "tau"))
          rightlabs.new[rightlabs.new == "tau"] <- "Tau"
        #
        if (any(rightcols.new == "I2"))
          rightlabs.new[rightlabs.new == "I2"] <- "I2"
        #
        if (three.level & any(rightcols.new == "cluster"))
          rightlabs.new[rightlabs.new == "cluster"] <- "Cluster"
        #
        if (n_of_1 & any(rightcols.new == "cycles"))
          rightlabs.new[rightlabs.new == "cycles"] <- "Cycles"
      }
      else {
        if (length(rightcols.new) == length(rightlabs))
          rightlabs.new <- rightlabs
        else if (length(rightcols.new) > length(rightlabs))
          stop("Too few labels defined in argument 'rightlabs'.")
        else {
          rightlabs.new <- rightcols.new
          #
          for (i in seq_along(rightcols.new)) {
            match1.i <- match(rightcols.new[i], rightcols)
            if (!is.na(rightlabs[match1.i]))
              rightlabs.new[i] <- rightlabs[match1.i]
            else {
              match2.i <- match(rightcols.new[i], colnames)
              if (!is.na(match2.i))
                rightlabs.new[i] <- labnames[match2.i]
              else if (rightcols.new[i] == "pval")
                rightlabs.new[i] <- "P-value"
              else if (rightcols.new[i] == "tau2")
                rightlabs.new[i] <- "Tau2"
              else if (rightcols.new[i] == "tau")
                rightlabs.new[i] <- "Tau"
              else if (rightcols.new[i] == "I2")
                rightlabs.new[i] <- "I2"
            }
          }
        }
      }
    }
    #
    if (length(leftcols.new) > 0) {
      if (missing.leftlabs) {
        leftlabs.new <- leftcols.new
        #
        if ((metacor | metaprop | metamean | metarate) &
            any(leftcols.new == "n"))
          leftlabs.new[leftlabs.new == "n"] <- "Total"
        #
        if (metamean & any(leftcols.new == "mean"))
          leftlabs.new[leftlabs.new == "mean"] <- "Mean"
        #
        if (metamean & any(leftcols.new == "sd"))
          leftlabs.new[leftlabs.new == "sd"] <- "SD"
        #
        if (metarate & any(leftcols.new == "time"))
          leftlabs.new[leftlabs.new == "time"] <- "Time"
        #
        if (any(leftcols.new == "pval"))
          leftlabs.new[leftlabs.new == "pval"] <- "P-value"
        #
        if (any(leftcols.new == "tau2"))
          leftlabs.new[leftlabs.new == "tau2"] <- "Tau2"
        #
        if (any(leftcols.new == "tau"))
          leftlabs.new[leftlabs.new == "tau"] <- "Tau"
        #
        if (any(leftcols.new == "I2"))
          leftlabs.new[leftlabs.new == "I2"] <- "I2"
        #
        if (three.level & any(leftcols.new == "cluster"))
          leftlabs.new[leftlabs.new == "cluster"] <- "Cluster"
        #
        if (n_of_1 & any(leftcols.new == "cycles"))
          leftlabs.new[leftlabs.new == "cycles"] <- "Cycles"
      }
      else {
        if (length(leftcols.new) == length(leftlabs))
          leftlabs.new <- leftlabs
        else if (length(leftcols.new) > length(leftlabs))
          stop("Too few labels defined in argument 'leftlabs'.")
        else {
          leftlabs.new <- leftcols.new
          #
          for (i in seq_along(leftcols.new)) {
            match1.i <- match(leftcols.new[i], leftcols)
            if (!is.na(leftlabs[match1.i]))
              leftlabs.new[i] <- leftlabs[match1.i]
            else {
              match2.i <- match(leftcols.new[i], colnames)
              if (!is.na(match2.i))
                leftlabs.new[i] <- labnames[match2.i]
              else if (leftcols.new[i] == "pval")
                leftlabs.new[i] <- "P-value"
              else if (leftcols.new[i] == "tau2")
                leftlabs.new[i] <- "Tau2"
              else if (leftcols.new[i] == "tau")
                leftlabs.new[i] <- "Tau"
              else if (leftcols.new[i] == "I2")
                leftlabs.new[i] <- "I2"
            }
          }
        }
      }
    }
  }
  #
  # Default set of columns if argument leftcols and / or
  # rightcols not specified
  #
  if (is.null(leftcols) && lsel) {
    #
    leftcols <- "studlab"
    #
    if (three.level)
      leftcols <- c(leftcols, "cluster")
    #
    if (n_of_1)
      leftcols <- c(leftcols, "cycles")
    #
    if (jama) {
      leftcols <- c(leftcols, "effect.ci")
    }
    else {
      if (metabin) {
        if (study.results) {
          if (bmj) {
            leftcols <- c(leftcols,
                          "event.n.e",
                          "event.n.c")
            label.e.attach <- "event.n.e"
            #
            if (is.null(bmj.text))
              label.e <- "No of events / total"
            else
              label.e <- bmj.text
          }
          else
            leftcols <- c(leftcols,
                          "event.e", "n.e",
                          "event.c", "n.c")
        }
        else {
          leftcols <- c(leftcols,
                        if (pooled.events) "event.e",
                        if (pooled.totals) "n.e",
                        if (pooled.events) "event.c",
                        if (pooled.totals) "n.c")
          if (pooled.events & !pooled.totals) {
            if (is.null(label.e.attach))
              label.e.attach <- "event.e"
            if (is.null(label.c.attach))
              label.c.attach <- "event.c"
          }
        }
      }
      #
      if (metacont) {
        if (study.results) {
          if (bmj) {
            leftcols <- c(leftcols,
                          "mean.sd.n.e",
                          "mean.sd.n.c")
            label.e.attach <- "mean.sd.n.e"
            #
            if (is.null(bmj.text))
              label.e <- "Mean (SD) / total"
            else
              label.e <- bmj.text
          }         
          else if (revman5)
            leftcols <- c(leftcols,
                          "mean.e", "sd.e", "n.e",
                          "mean.c", "sd.c", "n.c")
          else
            leftcols <- c(leftcols,
                          "n.e", "mean.e", "sd.e",
                          "n.c", "mean.c", "sd.c")
        }
        else if (pooled.totals) {
          leftcols <- c(leftcols, "n.e", "n.c")
          if (is.null(label.e.attach))
            label.e.attach <- "n.e"
          if (is.null(label.c.attach))
            label.c.attach <- "n.c"
        }
      }
      #
      if (metagen & study.results) {
        leftcols <- c(leftcols,
                      "TE", "seTE")
        if (!is.null(x$n.e)) {
          leftcols <- c(leftcols, "n.e")
          if (is.null(label.e.attach))
            label.e.attach <- "n.e"
        }
        if (!is.null(x$n.c)) {
          leftcols <- c(leftcols, "n.c")
          if (is.null(label.c.attach))
            label.c.attach <- "n.c"
        }
      }
      #
      if (metamean) {
        if (study.results) {
          if (revman5)
            leftcols <- c(leftcols,
                          "mean.e", "sd.e", "n.e")
          else
            leftcols <- c(leftcols,
                          "n.e", "mean.e", "sd.e")
        }
        else if (pooled.totals) {
          leftcols <- c(leftcols, "n.e")
          if (is.null(label.e.attach))
            label.e.attach <- "n.e"
        }
      }
      #
      if (metaprop) {
        if (study.results) {
          if (bmj) {
            leftcols <- c(leftcols,
                          "event.n.e")
            label.e.attach <- "event.n.e"
            #
            if (is.null(bmj.text))
              label.e <- "No of events / total"
            else
              label.e <- bmj.text
          
          }
          else
            leftcols <- c(leftcols, "event.e", "n.e")
        }
        else {
          leftcols <- c(leftcols,
                        if (pooled.events) "event.e",
                        if (pooled.totals) "n.e")
          if (pooled.events & !pooled.totals) {
            if (is.null(label.e.attach))
              label.e.attach <- "event.e"
          }
        }
      }
      #
      if (metarate) {
        if (study.results)
          leftcols <- c(leftcols,
                        "event.e", "time.e",
                        if (!is.null(x$n)) "n.e")
        else {
          leftcols <- c(leftcols,
                        if (pooled.events) "event.e",
                        if (pooled.times) "time.e")
          if (pooled.events & !pooled.times) {
            if (is.null(label.e.attach))
              label.e.attach <- "event.e"
          }
        }
      }
      #
      if (metacor) {
        if (study.results | pooled.totals)
          leftcols <- c(leftcols,
                        "n.e")
      }
      #
      if (metainc) {
        if (study.results)
          leftcols <- c(leftcols,
                        "event.e", "time.e",
                        "event.c", "time.c")
        else {
          leftcols <- c(leftcols,
                        if (pooled.events) "event.e",
                        if (pooled.times) "time.e",
                        if (pooled.events) "event.c",
                        if (pooled.times) "time.c")
          if (pooled.events & !pooled.times) {
            if (is.null(label.e.attach))
              label.e.attach <- "event.e"
            if (is.null(label.c.attach))
              label.c.attach <- "event.c"
          }
        }
      }
    }
    #
    # Add columns for RevMan 5 layout
    #
    if (revman5) {
      #
      if (overall & study.results & !any(x$method == "GLMM") & !metamerge) {
        if (common && !all(is.na(x$w.common)))
          leftcols <- c(leftcols, "w.common")
        if (random && !all(is.na(x$w.random)))
          leftcols <- c(leftcols, "w.random")
      }
      #
      leftcols <- c(leftcols, "effect.ci")
    }
    #
    # Add columns if risk of bias assessment is only information on
    # right side of the forest plot
    #
    if (!revman5.jama & rob.only) {
      #
      if (overall & study.results & !any(x$method == "GLMM") & !metamerge) {
        if (common)
          leftcols <- c(leftcols, "w.common")
        if (random)
          leftcols <- c(leftcols, "w.random")
      }
      #
      if (bmj)
        leftcols <- c(leftcols, "effect.ci")
      else
        leftcols <- c(leftcols, "effect", "ci")
    }
  }
  #
  leftcols[leftcols == "w.fixed"] <- "w.common"
  leftcols <- unique(leftcols)
  #
  if (is.null(rightcols) && rsel) {
    if (bmj)
      rightcols <- "effect.ci"
    else
      rightcols <- c("effect", "ci")
    #
    if (overall & study.results & !any(x$method == "GLMM") & !metamerge) {
      wcols <- c(if (common && !all(is.na(x$w.common))) "w.common",
                 if (random && !all(is.na(x$w.random))) "w.random")
      #
      if (bmj)
        rightcols <- c(wcols, rightcols)
      else
        rightcols <- c(rightcols, wcols)
    }
  }
  #
  if (RoB.available & rob.only &
      missing.leftcols & !revman5.jama) {
    if (overall & study.results & !any(x$method == "GLMM" & !metamerge)) {
      if (common)
        leftcols <- c(leftcols, "w.common")
      if (random)
        leftcols <- c(leftcols, "w.random")
    }
    #
    if (bmj)
      leftcols <- c(leftcols, "effect.ci")
    else
      leftcols <- c(leftcols, "effect", "ci")
  } 
  #
  rightcols[rightcols == "w.fixed"] <- "w.common"
  rightcols <- unique(rightcols)
  #
  if (rsel && !missing.rightlabs && length(rightlabs) > length(rightcols))
    stop("Too many labels defined in argument 'rightlabs': ",
         length(rightlabs), " label", if (length(rightlabs) > 1) "s",
         " for ", length(rightcols), " column",
         if (length(rightcols) > 1) "s",
         ".",
         call. = FALSE)
  #
  rightcols <- c(rightcols, rightcols.rob)
  #
  if (any(leftcols == "w.common") & any(rightcols == "w.common"))
    leftcols <- leftcols[!leftcols == "w.common"]
  #
  if (lsel && !missing.leftlabs && length(leftlabs) > length(leftcols))
    stop("Too many labels defined in argument 'leftlabs': ",
         length(leftlabs), " label", if (length(leftlabs) > 1) "s",
         " for ", length(leftcols), " column",
         if (length(leftcols) > 1) "s",
         ".",
         call. = FALSE)
  
  
  #
  #
  # (6) Select data for forest plot
  #
  #
  if (metacor) {
    x$n.e <- x$n
  }
  #
  if (metamean) {
    x$n.e <- x$n
    x$mean.e <- x$mean
    x$sd.e <- x$sd
    #
    if (!is.null(rightcols)) {
      if (any(rightcols == "n"))
        rightcols[rightcols == "n"] <- "n.e"
      if (any(rightcols == "mean"))
        rightcols[rightcols == "mean"] <- "mean.e"
      if (any(rightcols == "sd"))
        rightcols[rightcols == "sd"] <- "sd.e"
    }
    #
    if (!is.null(leftcols)) {
      if (any(leftcols == "n"))
        leftcols[leftcols == "n"] <- "n.e"
      if (any(leftcols == "mean"))
        leftcols[leftcols == "mean"] <- "mean.e"
      if (any(leftcols == "sd"))
        leftcols[leftcols == "sd"] <- "sd.e"
    }
  }
  #
  if (metaprop) {
    x$event.e <- x$event
    x$n.e <- x$n
    #
    if (!is.null(rightcols)) {
      if (any(rightcols == "n"))
        rightcols[rightcols == "n"] <- "n.e"
      if (any(rightcols == "event"))
        rightcols[rightcols == "event"] <- "event.e"
      if (any(rightcols == "event.n"))
        rightcols[rightcols == "event.n"] <- "event.n.e"
    }
    #
    if (!is.null(leftcols)) {
      if (any(leftcols == "n"))
        leftcols[leftcols == "n"] <- "n.e"
      if (any(leftcols == "event"))
        leftcols[leftcols == "event"] <- "event.e"
      if (any(leftcols == "event.n"))
        leftcols[leftcols == "event.n"] <- "event.n.e"
    }
  }
  #
  if (metarate) {
    x$event.e <- x$event
    x$time.e <- x$time
    x$n.e <- x$n
    #
    if (!is.null(rightcols)) {
      if (any(rightcols == "time"))
        rightcols[rightcols == "time"] <- "time.e"
      if (any(rightcols == "event"))
        rightcols[rightcols == "event"] <- "event.e"
      if (any(rightcols == "n"))
        rightcols[rightcols == "n"] <- "n.e"
    }
    #
    if (!is.null(leftcols)) {
      if (any(leftcols == "time"))
        leftcols[leftcols == "time"] <- "time.e"
      if (any(leftcols == "event"))
        leftcols[leftcols == "event"] <- "event.e"
      if (any(leftcols == "n"))
        leftcols[leftcols == "n"] <- "n.e"
    }
  }
  #
  # Calculate random effects weights in subgroup meta-analysis if
  # overall results are not shown (i.e., calculate weights based on
  # subgroup specific variance estimates)
  #
  if (random & by & !overall & !all(is.na(x$w.random)) &
      !(!is.null(x$tau.common) && x$tau.common)) {
    for (i in unique(subgroup)) {
      x$w.random[subgroup == i] <-
        1 / (x$seTE[subgroup == i]^2 + replaceNA(x$tau.w[i]^2, 0))
      x$w.random.w[i] <- sum(x$w.random[subgroup == i])
    }
  }
  #
  # Total number of studies to plot (*not* number of studies combined)
  #
  k.all <- length(x$TE)
  #
  if (allstudies)
    n.stud <- k.all # all studies
  else
    n.stud <- sum(!is.na(x$TE)) # number of studies with treatment estimates
  #
  if (length(type.study) == 1)
    type.study <- rep(type.study, k.all)
  else if (length(type.study) != k.all)
    stop("Argument 'type.study' must be a single character or of ",
         "same length as number of studies.")
  #
  if (!by)
    subgroup <- rep(1, k.all)
  #
  if (by & anyNA(subgroup))
    stop("Missing values in 'subgroup'")
  #
  if (allstudies)
    sel <- rep(TRUE, k.all)
  else
    sel <- !is.na(x$TE)
  #
  x$n.e <- x$n.e[sel]
  x$n.c <- x$n.c[sel]
  #
  x$event.e <- x$event.e[sel]
  x$event.c <- x$event.c[sel]
  #
  x$mean.e <- x$mean.e[sel]
  x$mean.c <- x$mean.c[sel]
  #
  x$sd.e <- x$sd.e[sel]
  x$sd.c <- x$sd.c[sel]
  #
  x$cor <- x$cor[sel]
  #
  x$time.e <- x$time.e[sel]
  x$time.c <- x$time.c[sel]
  #
  x$TE <- x$TE[sel]
  x$seTE <- x$seTE[sel]
  #
  x$lower <- x$lower[sel]
  x$upper <- x$upper[sel]
  #
  if (is.matrix(x$w.common))
    x$w.common <- x$w.common[sel, , drop = FALSE]
  else
    x$w.common <- x$w.common[sel]
  #
  if (is.matrix(x$w.random))
    x$w.random <- x$w.random[sel, , drop = FALSE]
  else
    x$w.random <- x$w.random[sel]
  #
  studlab <- studlab[sel]
  type.study <- type.study[sel]
  #
  x$n.harmonic.mean <- x$n.harmonic.mean[sel]
  x$t.harmonic.mean <- x$t.harmonic.mean[sel]
  #
  x$pval <- x$pval[sel]
  #
  x$cluster <- x$cluster[sel]
  #
  x$cycles <- x$cycles[sel]
  #
  subgroup <- subgroup[sel]
  sortvar <- sortvar[sel]
  #
  col.study <- col.study[sel]
  col.square <- col.square[sel]
  col.square.lines <- col.square.lines[sel]
  col.circle <- col.circle[sel]
  col.circle.lines <- col.circle.lines[sel]
  #
  col.inside <- col.inside[sel]
  #
  avail.exclude <- !is.null(x$exclude)
  if (avail.exclude)
    x$exclude <- x$exclude[sel]
  #
  if (sort | by) {
    if (sort.subgroup)
      subgroup.levels <- sort(x$subgroup.levels)
    else
      subgroup.levels <- x$subgroup.levels
    #
    subgroup.factor <- factor(subgroup, levels = subgroup.levels)
    o <- order(subgroup.factor, sortvar)
    #
    x$cluster <- x$cluster[o]
    x$cycles <- x$cycles[o]
    #
    x$n.e <- x$n.e[o]
    x$n.c <- x$n.c[o]
    #
    x$event.e <- x$event.e[o]
    x$event.c <- x$event.c[o]
    #
    x$mean.e <- x$mean.e[o]
    x$mean.c <- x$mean.c[o]
    #
    x$sd.e <- x$sd.e[o]
    x$sd.c <- x$sd.c[o]
    #
    x$cor <- x$cor[o]
    #
    x$time.e <- x$time.e[o]
    x$time.c <- x$time.c[o]
    #
    x$TE   <- x$TE[o]
    x$seTE <- x$seTE[o]
    #
    x$lower <- x$lower[o]
    x$upper <- x$upper[o]
    #
    if (is.matrix(x$w.common))
      x$w.common <- x$w.common[o, , drop = FALSE]
    else
      x$w.common <- x$w.common[o]
    #
    if (is.matrix(x$w.random))
      x$w.random <- x$w.random[o, , drop = FALSE]
    else
      x$w.random <- x$w.random[o]
    #
    studlab  <- studlab[o]
    type.study  <- type.study[o]
    #
    x$n.harmonic.mean <- x$n.harmonic.mean[o]
    x$t.harmonic.mean <- x$t.harmonic.mean[o]
    #
    x$pval <- x$pval[o]
    #
    subgroup <- subgroup[o]
    sortvar <- sortvar[o]
    #
    col.study <- col.study[o]
    col.square <- col.square[o]
    col.square.lines <- col.square.lines[o]
    col.circle <- col.circle[o]
    col.circle.lines <- col.circle.lines[o]
    #
    col.inside <- col.inside[o]
    #
    if (avail.exclude)
      x$exclude <- x$exclude[o]
    #
    if (newcols) {
      dataset1 <- dataset1[o, ]
      dataset2 <- dataset2[o, ]
      #
      if (RoB.available)
        rob <- rob[o, ]
    }
  }
  #
  if (bmj)
    studlab <- paste0(" ", studlab)
  #
  TE <- x$TE
  seTE <- x$seTE
  lowTE <- x$lower
  uppTE <- x$upper
  #
  TE.common <- x$TE.common
  lowTE.common <- x$lower.common
  uppTE.common <- x$upper.common
  #
  TE.random <- x$TE.random
  lowTE.random <- x$lower.random
  uppTE.random <- x$upper.random
  if (n.ran > 1 && length(TE.random) == 1)
    TE.random <- rep(TE.random, n.ran)
  #
  lowTE.predict <- x$lower.predict
  uppTE.predict <- x$upper.predict
  #
  if (LRT)
    Q <- x$Q.LRT
  else
    Q <- unlist(x$Q)
  #
  df.Q <- unlist(x$df.Q)
  pval.Q <- replaceNULL(unlist(x$pval.Q), pvalQ(Q, df.Q))
  #
  # Keep the first heterogeneity statistics
  #
  if (length(Q) > 1) {
    Q <- Q[1]
    df.Q <- df.Q[1]
    pval.Q <- pval.Q[1]
  }
  #
  tau2 <- x$tau2
  lower.tau2 <- x$lower.tau2
  upper.tau2 <- x$upper.tau2
  #
  if (is.null(tau2)) {
    tau2 <- x$tau^2
    lower.tau2 <- upper.tau2 <- NA
  }
  #
  if (length(tau2) > 1) {
    if (three.level) {
      tau2 <- sum(tau2)
      lower.tau2 <- NA
      upper.tau2 <- NA
    }
    else {
      tau2 <- tau2[1]
      lower.tau2 <- lower.tau2[1]
      upper.tau2 <- upper.tau2[1]
    }
  }
  #
  tau <- x$tau
  lower.tau <- x$lower.tau
  upper.tau <- x$upper.tau
  #
  if (length(tau) > 1) {
    if (three.level) {
      tau <- sqrt(sum(tau^2))
      lower.tau2 <- NA
      upper.tau2 <- NA
    }
    else {
      tau <- tau[1]
      lower.tau <- lower.tau[1]
      upper.tau <- upper.tau[1]
    }
  }
  #
  sign.lower.tau <- x$sign.lower.tau
  sign.upper.tau <- x$sign.lower.tau
  #
  if (is.null(lower.tau)) {
    lower.tau <- upper.tau <- NA
    sign.lower.tau <- sign.upper.tau <- ""
  }
  #
  I2 <- unlist(x$I2)
  lowI2 <- unlist(x$lower.I2)
  uppI2 <- unlist(x$upper.I2)
  #
  if (length(I2) > 1) {
    I2 <- I2[1]
    lowI2 <- lowI2[1]
    uppI2 <- uppI2[1]
  }
  #
  Rb <- unlist(x$Rb)
  lowRb <- unlist(x$lower.Rb)
  uppRb <- unlist(x$upper.Rb)
  #
  if (length(Rb) > 1) {
    Rb <- Rb[1]
    lowRb <- lowRb[1]
    uppRb <- uppRb[1]
  }
  #
  if (by) {
    Q.b.common <- unlist(x$Q.b.common)
    names(Q.b.common) <- colnames(lower.common.w)
    Q.w.common <- x$Q.w.common
    #
    Q.b.random <- unlist(x$Q.b.random)
    Q.w.random <- x$Q.w.random
    #
    Q.w <- x$Q.w
    #
    df.Q.w <- replaceNULL(x$df.Q.w, sum((k.w - 1)[!is.na(x$Q.w)]))
    df.Q.b <- replaceNULL(x$df.Q.b, (k - 1) - sum((k.w - 1)[!is.na(x$Q.w)]))
    df.Q.b.common <- replaceNULL(x$df.Q.b.common, df.Q.b)
    df.Q.b.random <- replaceNULL(x$df.Q.b.random, df.Q.b)    
    #
    pval.Q.b.common <-
      replaceNULL(x$pval.Q.b.common, pvalQ(Q.b.common, df.Q.b.common))
    pval.Q.w.common <-
      replaceNULL(x$pval.Q.w.common, pvalQ(Q.w.common, df.Q.w))
    #
    pval.Q.b.random <-
      unlist(replaceNULL(x$pval.Q.b.random, pvalQ(Q.b.random, df.Q.b.random)))
    pval.Q.w.random <-
      unlist(replaceNULL(x$pval.Q.w.random, pvalQ(Q.w.random, df.Q.w)))
    #
    Q.resid <- x$Q.resid
    df.Q.resid <- x$df.Q.resid
    pval.Q.resid <- x$pval.Q.resid
    #
    tau2.resid <- x$tau2.resid
    lower.tau2.resid <- x$lower.tau2.resid
    upper.tau2.resid <- x$upper.tau2.resid
    if (is.null(tau2.resid)) {
      tau2.resid <- x$tau.resid^2
      lower.tau2.resid <- upper.tau2.resid <- NA
    }
    tau.resid <- x$tau.resid
    lower.tau.resid <- x$lower.tau.resid
    upper.tau.resid <- x$upper.tau.resid
    sign.lower.tau.resid <- x$sign.lower.tau.resid
    sign.upper.tau.resid <- x$sign.lower.tau.resid
    if (is.null(lower.tau.resid)) {
      lower.tau.resid <- upper.tau.resid <- NA
      sign.lower.tau.resid <- sign.upper.tau.resid <- ""
    }
    #
    I2.resid <- x$I2.resid
    lowI2.resid <- x$lower.I2.resid
    uppI2.resid <- x$upper.I2.resid
  }
  else {
    Q.b.common <- Q.b.random <- NA
    df.Q.b <- df.Q.b.common <- df.Q.b.random <- NA
    pval.Q.b.common <- pval.Q.b.random <- NA
  }
  #
  hetstat.overall <- ""
  #
  if (overall.hetstat || is.character(hetstat)) {
    #
    hetstat.I2 <-
      paste0(hetseparator,
             formatN(100 * I2, digits.I2, "NA"), "%",
             if (print.I2.ci && !(is.na(lowI2) | is.na(uppI2)))
               pasteCI(100 * lowI2, 100 * uppI2,
                       digits.I2, big.mark,
                       text.NA = lab.NA, unit = "%"))
    #
    hetstat.tau2 <-
      paste0(formatPT(tau2, digits = digits.tau2, big.mark = big.mark,
                      lab = TRUE, labval = "", lab.NA = "NA"),
             if (print.tau2.ci && !(is.na(lower.tau2) | is.na(upper.tau2)))
               pasteCI(lower.tau2, upper.tau2, digits.tau2, big.mark,
                       sign.lower.tau, sign.upper.tau, lab.NA))
    #
    hetstat.tau <-
      paste0(formatPT(tau, digits = digits.tau, big.mark = big.mark,
                      lab = TRUE, labval = "", lab.NA = "NA"),
             if (print.tau.ci && !(is.na(lower.tau) | is.na(upper.tau)))
               pasteCI(lower.tau, upper.tau, digits.tau, big.mark,
                       sign.lower.tau, sign.upper.tau, lab.NA))
    #
    hetstat.Q <-
      paste0(hetseparator,
             formatN(Q, digits.Q, "NA", big.mark = big.mark),
             if (bmj.revman5) ", df",
             if (bmj.revman5) hetseparator,
             if (bmj.revman5) df.Q)
    #
    hetstat.pval.Q <-
      formatPT(pval.Q,
               lab = TRUE, labval = "",
               digits = digits.pval.Q,
               zero = zero.pval, JAMA = JAMA.pval,
               scientific = scientific.pval,
               lab.NA = "NA")
    #
    if (print.Rb)
      hetstat.Rb <-
        paste0(hetseparator,
               formatN(100 * Rb, digits.I2, "NA", big.mark = big.mark),
               "%",
               if (print.Rb.ci && !(is.na(lowRb) | is.na(uppRb)))
                 pasteCI(100 * lowRb, 100 * uppRb,
                         digits.I2, big.mark,
                         text.NA = lab.NA, unit = "%"))
    #
    # Remove superfluous spaces
    #
    while(grepl("  ", hetstat.I2))
      hetstat.I2 <- gsub("  ", " ", hetstat.I2)
    while(grepl("  ", hetstat.tau2))
      hetstat.tau2 <- gsub("  ", " ", hetstat.tau2)
    while(grepl("  ", hetstat.tau))
      hetstat.tau <- gsub("  ", " ", hetstat.tau)
    while(grepl("  ", hetstat.Q))
      hetstat.Q <- gsub("  ", " ", hetstat.Q)
    while(grepl("  ", hetstat.pval.Q))
      hetstat.pval.Q <- gsub("  ", " ", hetstat.pval.Q)
    if (print.Rb)
      while(grepl("  ", hetstat.Rb))
        hetstat.Rb <- gsub("  ", " ", hetstat.Rb)
    if (bmj) {
      hetstat.tau2 <- gsub(" = ", "=", hetstat.tau2)
      hetstat.tau <- gsub(" = ", "=", hetstat.tau)
      hetstat.Q <- gsub(" = ", "=", hetstat.Q)
      hetstat.pval.Q <- gsub(" = ", "=", hetstat.pval.Q)
      hetstat.I2 <- gsub(" = ", "=", hetstat.I2)
    }
    #
    if (bmj) {
      if (print.tau)
        hetstat.overall <-
          substitute(
            paste(hl,
                  tau, ht, "; ",
                  chi^2, hq,
                  ", P", hp, "; ",
                  I^2, hi),
            list(hl = hetlab,
                 ht = hetstat.tau,
                 hq = hetstat.Q,
                 hp = hetstat.pval.Q,
                 hi = hetstat.I2,
                 df = df.Q)
            )
      else if (print.tau2)
        hetstat.overall <-
          substitute(
            paste(hl,
                  tau^2, ht, "; ",
                  chi^2, hq,
                  ", P", hp, "; ",
                  I^2, hi),
            list(hl = hetlab,
                 ht = hetstat.tau2,
                 hq = hetstat.Q,
                 hp = hetstat.pval.Q,
                 hi = hetstat.I2,
                 df = df.Q)
            )
      else
        hetstat.overall <-
          substitute(
            paste(hl,
                  chi^2, hq,
                  ", P", hp, "; ",
                  I^2, hi),
            list(hl = hetlab,
                 hi = hetstat.I2,
                 hq = hetstat.Q,
                 hp = hetstat.pval.Q,
                 df = df.Q)
            )
    }
    else if (jama) {
      if (!missing.print.tau2 | !missing.print.tau) {
        if (print.tau)
          hetstat.overall <-
            substitute(
              paste(hl,
                    chi[df]^2, hq,
                    " (", italic(P), hp, "), ",
                    italic(I)^2, hi, ", ",
                    tau, ht),
              list(hl = hetlab,
                   df = df.Q,
                   hq = hetstat.Q,
                   hp = hetstat.pval.Q,
                   hi = hetstat.I2,
                   ht = hetstat.tau)
              )
        else if (print.tau2)
          hetstat.overall <-
            substitute(
              paste(hl,
                    chi[df]^2, hq,
                    " (", italic(P), hp, "), ",
                    italic(I)^2, hi, ", ",
                    tau^2, ht),
              list(hl = hetlab,
                   df = df.Q,
                   hq = hetstat.Q,
                   hp = hetstat.pval.Q,
                   hi = hetstat.I2,
                   ht = hetstat.tau2)
            )
        else
          hetstat.overall <-
            substitute(
              paste(hl,
                    chi[df]^2, hq,
                    " (", italic(P), hp, "), ",
                    italic(I)^2, hi),
              list(hl = hetlab,
                   df = df.Q,
                   hq = hetstat.Q,
                   hp = hetstat.pval.Q,
                   hi = hetstat.I2)
            )
      }
      else
        hetstat.overall <-
          substitute(
            paste(hl,
                  chi[df]^2, hq,
                  " (", italic(P), hp, "), ",
                  italic(I)^2, hi),
            list(hl = hetlab,
                 df = df.Q,
                 hq = hetstat.Q,
                 hp = hetstat.pval.Q,
                 hi = hetstat.I2)
          )
    }
    else if (revman5) {
      if (print.tau)
        hetstat.overall <-
          substitute(
            paste(hl,
                  "Tau", ht, "; ",
                  "Chi"^2, hq,
                  " (", P, hp, "); ",
                  I^2, hi),
            list(hl = hetlab,
                 ht = hetstat.tau,
                 hq = hetstat.Q,
                 hp = hetstat.pval.Q,
                 hi = hetstat.I2)
          )
      else
        hetstat.overall <-
          substitute(
            paste(hl,
                  "Tau"^2, ht, "; ",
                  "Chi"^2, hq,
                  " (", P, hp, "); ",
                  I^2, hi),
            list(hl = hetlab,
                 ht = hetstat.tau2,
                 hq = hetstat.Q,
                 hp = hetstat.pval.Q,
                 hi = hetstat.I2)
          )
    }
    else {
      #
      # One
      #
      if (print.I2 & !print.tau2.tau & !print.Q & !print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl, italic(I)^2, hi),
                     list(hl = hetlab, hi = hetstat.I2))
      else if (!print.I2 & print.tau2.tau & !print.Q & !print.pval.Q &
               !print.Rb)
        if (print.tau2)
          hetstat.overall <- substitute(paste(hl, tau^2, ht),
                                        list(hl = hetlab, ht = hetstat.tau2))
        else
          hetstat.overall <- substitute(paste(hl, tau, ht),
                                        list(hl = hetlab, ht = hetstat.tau))
      else if (!print.I2 & !print.tau2.tau & print.Q & !print.pval.Q &
               !print.Rb)
        hetstat.overall <-
          substitute(paste(hl, chi[df]^2, hq),
                     list(hl = hetlab, df = df.Q, hq = hetstat.Q))
      else if (!print.I2 & !print.tau2.tau & !print.Q & print.pval.Q &
               !print.Rb)
        hetstat.overall <-
          substitute(paste(hl, italic(p), hp),
                     list(hl = hetlab, hp = hetstat.pval.Q))
      else if (!print.I2 & !print.tau2.tau & !print.Q & !print.pval.Q &
               print.Rb)
        hetstat.overall <-
          substitute(paste(hl, italic(R)[italic(b)], hb),
                     list(hl = hetlab, hb = hetstat.Rb))
      #
      # Two
      #
      else if (print.I2 & print.tau2.tau & !print.Q & !print.pval.Q & !print.Rb)
        if (print.tau2)
          hetstat.overall <-
            substitute(paste(hl, italic(I)^2, hi,
                             ", ",
                             tau^2, ht),
                       list(hl = hetlab,
                            hi = hetstat.I2, ht = hetstat.tau2))
        else
          hetstat.overall <-
            substitute(paste(hl, italic(I)^2, hi,
                             ", ",
                             tau, ht),
                       list(hl = hetlab,
                            hi = hetstat.I2, ht = hetstat.tau))
      else if (print.I2 & !print.tau2.tau & print.Q & !print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl, italic(I)^2, hi,
                           ", ",
                           chi[df]^2, hq),
                     list(hl = hetlab, df = df.Q,
                          hi = hetstat.I2, hq = hetstat.Q))
      else if (print.I2 & !print.tau2.tau & !print.Q & print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl, italic(I)^2, hi,
                           ", ",
                           italic(p), hp),
                     list(hl = hetlab,
                          hi = hetstat.I2, hp = hetstat.pval.Q))
      else if (print.I2 & !print.tau2.tau & !print.Q & !print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl, italic(I)^2, hi,
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab,
                          hi = hetstat.I2, hb = hetstat.Rb))
      else if (!print.I2 & print.tau2.tau & print.Q & !print.pval.Q & !print.Rb)
        if (print.tau2)
          hetstat.overall <-
            substitute(paste(hl, tau^2, ht,
                             ", ",
                             chi[df]^2, hq),
                       list(hl = hetlab, df = df.Q,
                            ht = hetstat.tau2, hq = hetstat.Q))
        else
          hetstat.overall <-
            substitute(paste(hl, tau, ht,
                             ", ",
                             chi[df]^2, hq),
                       list(hl = hetlab, df = df.Q,
                            ht = hetstat.tau, hq = hetstat.Q))
      else if (!print.I2 & print.tau2.tau & !print.Q & print.pval.Q & !print.Rb)
        if (print.tau2)
          hetstat.overall <-
            substitute(paste(hl, tau^2, ht,
                             ", ",
                             italic(p), hp),
                       list(hl = hetlab,
                            ht = hetstat.tau2, hp = hetstat.pval.Q))
        else
          hetstat.overall <-
            substitute(paste(hl, tau, ht,
                             ", ",
                             italic(p), hp),
                       list(hl = hetlab,
                            ht = hetstat.tau, hp = hetstat.pval.Q))
      else if (!print.I2 & print.tau2.tau & !print.Q & !print.pval.Q & print.Rb)
        if (print.tau2)
          hetstat.overall <-
            substitute(paste(hl, tau^2, ht,
                             ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = hetlab,
                            ht = hetstat.tau2, hb = hetstat.Rb))
        else
          hetstat.overall <-
            substitute(paste(hl, tau, ht,
                             ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = hetlab,
                            ht = hetstat.tau, hb = hetstat.Rb))
      else if (!print.I2 & !print.tau2.tau & print.Q & print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl, chi[df]^2, hq,
                           " (",
                           italic(p), hp, ")"),
                     list(hl = hetlab, df = df.Q,
                          hq = hetstat.Q, hp = hetstat.pval.Q))
      else if (!print.I2 & !print.tau2.tau & print.Q & !print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl, chi[df]^2, hq,
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab, df = df.Q,
                          hq = hetstat.Q, hb = hetstat.Rb))
      else if (!print.I2 & !print.tau2.tau & !print.Q & print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl, italic(p), hp,
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab,
                          hp = hetstat.pval.Q, hb = hetstat.Rb))
      #
      # Three
      #
      else if (print.I2 & print.tau2.tau & print.Q & !print.pval.Q & !print.Rb)
        if (print.tau2)
          hetstat.overall <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau^2, ht, ", ",
                             chi[df]^2, hq),
                       list(hl = hetlab, df = df.Q,
                            hi = hetstat.I2, ht = hetstat.tau2,
                            hq = hetstat.Q))
        else
          hetstat.overall <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau, ht, ", ",
                             chi[df]^2, hq),
                       list(hl = hetlab, df = df.Q,
                            hi = hetstat.I2, ht = hetstat.tau,
                            hq = hetstat.Q))
      else if (print.I2 & print.tau2.tau & !print.Q & print.pval.Q & !print.Rb)
        if (print.tau2)
          hetstat.overall <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau^2, ht, ", ",
                             italic(p), hp),
                       list(hl = hetlab,
                            hi = hetstat.I2, ht = hetstat.tau2,
                            hp = hetstat.pval.Q))
        else
          hetstat.overall <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau, ht, ", ",
                             italic(p), hp),
                       list(hl = hetlab,
                            hi = hetstat.I2, ht = hetstat.tau,
                            hp = hetstat.pval.Q))

      else if (print.I2 & !print.tau2.tau & print.Q & print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           chi[df]^2, hq,
                           " (", italic(p), hp, ")"),
                     list(hl = hetlab, df = df.Q,
                          hi = hetstat.I2,
                          hq = hetstat.Q, hp = hetstat.pval.Q))
      else if (!print.I2 & print.tau2.tau & print.Q & print.pval.Q & !print.Rb)
        if (print.tau2)
          hetstat.overall <-
            substitute(paste(hl,
                             tau^2, ht, ", ",
                             chi[df]^2, hq,
                             " (", italic(p), hp, ")"),
                       list(hl = hetlab, df = df.Q,
                            ht = hetstat.tau2,
                            hq = hetstat.Q, hp = hetstat.pval.Q))
        else
          hetstat.overall <-
            substitute(paste(hl,
                             tau, ht, ", ",
                             chi[df]^2, hq,
                             " (", italic(p), hp, ")"),
                       list(hl = hetlab, df = df.Q,
                            ht = hetstat.tau,
                            hq = hetstat.Q, hp = hetstat.pval.Q))
      
      else if (print.I2 & print.tau2.tau & !print.Q & !print.pval.Q & print.Rb)
        if (print.tau2)
          hetstat.overall <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau^2, ht, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = hetlab,
                            hi = hetstat.I2, ht = hetstat.tau2,
                            hb = hetstat.Rb))
        else
          hetstat.overall <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau, ht, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = hetlab,
                            hi = hetstat.I2, ht = hetstat.tau,
                            hb = hetstat.Rb))          
      else if (print.I2 & !print.tau2.tau & print.Q & !print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           chi[df]^2, hq, ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab, df = df.Q,
                          hi = hetstat.I2,
                          hq = hetstat.Q,
                          hb = hetstat.Rb))
      else if (!print.I2 & print.tau2.tau & print.Q & !print.pval.Q & print.Rb)
        if (print.tau2)
          hetstat.overall <-
            substitute(paste(hl,
                             tau^2, ht, ", ",
                             chi[df]^2, hq, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = hetlab, df = df.Q,
                            ht = hetstat.tau2,
                            hq = hetstat.Q,
                            hb = hetstat.Rb))
        else
          hetstat.overall <-
            substitute(paste(hl,
                             tau, ht, ", ",
                             chi[df]^2, hq, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = hetlab, df = df.Q,
                            ht = hetstat.tau,
                            hq = hetstat.Q,
                            hb = hetstat.Rb))          
      else if (print.I2 & !print.tau2.tau & !print.Q & print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           italic(p), hp, ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab,
                          hi = hetstat.I2, hp = hetstat.pval.Q,
                          hb = hetstat.Rb))
      else if (!print.I2 & print.tau2.tau & !print.Q & print.pval.Q & print.Rb)
        if (print.tau2)
          hetstat.overall <-
            substitute(paste(hl,
                             tau^2, ht, ", ",
                             italic(p), hp, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = hetlab,
                            ht = hetstat.tau2,
                            hp = hetstat.pval.Q,
                            hb = hetstat.Rb))
        else
          hetstat.overall <-
            substitute(paste(hl,
                             tau, ht, ", ",
                             italic(p), hp, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = hetlab,
                            ht = hetstat.tau,
                            hp = hetstat.pval.Q,
                            hb = hetstat.Rb))
      else if (!print.I2 & !print.tau2.tau & print.Q & print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           chi[df]^2, hq,
                           " (", italic(p), hp, ")",
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab, df = df.Q,
                          hq = hetstat.Q, hp = hetstat.pval.Q,
                          hb = hetstat.Rb))      
      #
      # Four
      #
      if (print.I2 & print.tau2.tau & print.Q & print.pval.Q & !print.Rb)
        if (print.tau2)
          hetstat.overall <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau^2, ht, ", ",
                             chi[df]^2, hq,
                             " (", italic(p), hp, ")"),
                       list(hl = hetlab, df = df.Q,
                            hi = hetstat.I2, ht = hetstat.tau2,
                            hq = hetstat.Q, hp = hetstat.pval.Q))
        else
          hetstat.overall <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau, ht, ", ",
                             chi[df]^2, hq,
                             " (", italic(p), hp, ")"),
                       list(hl = hetlab, df = df.Q,
                            hi = hetstat.I2, ht = hetstat.tau,
                            hq = hetstat.Q, hp = hetstat.pval.Q))          
      else if (print.I2 & print.tau2.tau & print.Q & !print.pval.Q & print.Rb)
        if (print.tau2)
          hetstat.overall <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau^2, ht, ", ",
                             chi[df]^2, hq, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = hetlab, df = df.Q,
                            hi = hetstat.I2, ht = hetstat.tau2,
                            hq = hetstat.Q, hb = hetstat.Rb))
        else
          hetstat.overall <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau, ht, ", ",
                             chi[df]^2, hq, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = hetlab, df = df.Q,
                            hi = hetstat.I2, ht = hetstat.tau,
                            hq = hetstat.Q, hb = hetstat.Rb))
      else if (print.I2 & print.tau2.tau & !print.Q & print.pval.Q & print.Rb)
        if (print.tau2)
          hetstat.overall <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau^2, ht, ", ",
                             italic(p), hp, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = hetlab,
                            hi = hetstat.I2, ht = hetstat.tau2,
                            hp = hetstat.pval.Q, hb = hetstat.Rb))
        else
          hetstat.overall <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau, ht, ", ",
                             italic(p), hp, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = hetlab,
                            hi = hetstat.I2, ht = hetstat.tau,
                            hp = hetstat.pval.Q, hb = hetstat.Rb))
      else if (print.I2 & !print.tau2.tau & print.Q & print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           chi[df]^2, hq,
                           " (", italic(p), hp, ")",
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab, df = df.Q,
                          hi = hetstat.I2,
                          hq = hetstat.Q, hp = hetstat.pval.Q,
                          hb = hetstat.Rb))
      else if (!print.I2 & print.tau2.tau & print.Q & print.pval.Q & print.Rb)
        if (print.tau2)
          hetstat.overall <-
            substitute(paste(hl,
                             tau^2, ht, ", ",
                             chi[df]^2, hq,
                             " (", italic(p), hp, ")",
                             ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = hetlab, df = df.Q,
                            ht = hetstat.tau2,
                            hq = hetstat.Q, hp = hetstat.pval.Q,
                            hb = hetstat.Rb))
        else
          hetstat.overall <-
            substitute(paste(hl,
                             tau, ht, ", ",
                             chi[df]^2, hq,
                             " (", italic(p), hp, ")",
                             ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = hetlab, df = df.Q,
                            ht = hetstat.tau,
                            hq = hetstat.Q, hp = hetstat.pval.Q,
                            hb = hetstat.Rb))
      #
      # Five
      #
      else if (print.I2 & print.tau2.tau & print.Q & print.pval.Q & print.Rb)
        if (print.tau2)
          hetstat.overall <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau^2, ht, ", ",
                             chi[df]^2, hq,
                             " (", italic(p), hp, ")",
                             ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = hetlab, df = df.Q,
                            hi = hetstat.I2, ht = hetstat.tau2,
                            hq = hetstat.Q, hp = hetstat.pval.Q,
                            hb = hetstat.Rb))
        else
          hetstat.overall <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau, ht, ", ",
                             chi[df]^2, hq,
                             " (", italic(p), hp, ")",
                             ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = hetlab, df = df.Q,
                            hi = hetstat.I2, ht = hetstat.tau,
                            hq = hetstat.Q, hp = hetstat.pval.Q,
                            hb = hetstat.Rb))
    }
  }
  #
  # Line with residual heterogeneity
  #
  hetstat.resid <- ""
  #
  if (by && (length(tau2.resid) == 0 || is.na(tau2.resid)))
    print.tau2.tau.resid <- FALSE
  else
    print.tau2.tau.resid <- print.tau2.tau
  #
  if (resid.hetstat) {
    #
    hetstat.I2.resid <-
      paste0(hetseparator,
             formatN(100 * I2.resid, digits.I2, "NA"), "%",
             if (print.I2.ci && !(is.na(lowI2.resid) | is.na(uppI2.resid)))
               pasteCI(100 * lowI2.resid, 100 * uppI2.resid,
                       digits.I2, big.mark,
                       text.NA = lab.NA, unit = "%"))
    #
    hetstat.tau2.resid <-
      paste0(formatPT(tau2.resid, digits = digits.tau2, big.mark = big.mark,
                      lab = TRUE, labval = "", lab.NA = "NA"),
             if (print.tau2.ci &&
                 !(is.na(lower.tau2.resid) | is.na(upper.tau2.resid)))
               pasteCI(lower.tau2.resid, upper.tau2.resid, digits.tau2,
                       big.mark,
                       sign.lower.tau.resid, sign.upper.tau.resid, lab.NA))
    hetstat.tau.resid <-
      paste0(formatPT(tau.resid, digits = digits.tau, big.mark = big.mark,
                      lab = TRUE, labval = "", lab.NA = "NA"),
             if (print.tau.ci &&
                 !(is.na(lower.tau.resid) | is.na(upper.tau.resid)))
               pasteCI(lower.tau.resid, upper.tau.resid, digits.tau,
                       big.mark,
                       sign.lower.tau.resid, sign.upper.tau.resid, lab.NA))
    #
    hetstat.Q.resid <-
      paste0(hetseparator,
             formatN(Q.resid, digits.Q, "NA", big.mark = big.mark),
             if (bmj.revman5) ", df",
             if (bmj.revman5) hetseparator,
             if (bmj.revman5) df.Q.resid)
    #
    hetstat.pval.Q.resid <-
      formatPT(pval.Q.resid,
               lab = TRUE, labval = "",
               digits = digits.pval.Q,
               zero = zero.pval, JAMA = JAMA.pval,
               scientific = scientific.pval,
               lab.NA = "NA")
    #
    hetstat.Rb.resid <- ""
    print.Rb.resid <- FALSE
    #
    # Remove superfluous spaces
    #
    while(grepl("  ", hetstat.I2.resid))
      hetstat.I2.resid <- gsub("  ", " ", hetstat.I2.resid)
    while(grepl("  ", hetstat.tau2.resid))
      hetstat.tau2.resid <- gsub("  ", " ", hetstat.tau2.resid)
    while(grepl("  ", hetstat.tau.resid))
      hetstat.tau.resid <- gsub("  ", " ", hetstat.tau.resid)
    while(grepl("  ", hetstat.Q.resid))
      hetstat.Q.resid <- gsub("  ", " ", hetstat.Q.resid)
    while(grepl("  ", hetstat.pval.Q.resid))
      hetstat.pval.Q.resid <- gsub("  ", " ", hetstat.pval.Q.resid)
    while(grepl("  ", hetstat.Rb.resid))
      hetstat.Rb.resid <- gsub("  ", " ", hetstat.Rb.resid)
    #
    if (bmj) {
      if (print.tau)
        hetstat.resid <-
          substitute(
            paste(hl,
                  tau, ht, "; ",
                  chi^2, hq,
                  ", df=",
                  df,
                  ", P"=, hp,
                  "; ",
                  I^2, hi),
            list(hl = resid.hetlab,
                 ht = hetstat.tau.resid,
                 df = df.Q,
                 hq = hetstat.Q.resid,
                 hp = hetstat.pval.Q.resid,
                 hi = hetstat.I2.resid)
            )
      else
        hetstat.resid <-
          substitute(
            paste(hl,
                  tau^2, ht, "; ",
                  chi^2, hq,
                  ", df=",
                  df,
                  ", P"=, hp,
                  "; ",
                  I^2, hi),
              list(hl = resid.hetlab,
                   df = df.Q.resid,
                   hq = hetstat.Q.resid,
                   hp = hetstat.pval.Q.resid,
                   hi = hetstat.I2.resid)
          )
    }
    else if (jama) {
      if (!missing.print.tau2 | !missing.print.tau) {
        if (print.tau)
          hetstat.resid <-
            substitute(
              paste(hl,
                    chi[df]^2, hq,
                    " (",
                    italic(P), hp,
                    "), ",
                    italic(I)^2, hi,
                    tau, ht),
              list(hl = resid.hetlab,
                   df = df.Q.resid,
                   hq = hetstat.Q.resid,
                   hp = hetstat.pval.Q.resid,
                   hi = hetstat.I2.resid,
                   ht = hetstat.tau.resid)
            )
        else if (print.tau2)
          hetstat.resid <-
            substitute(
              paste(hl,
                    chi[df]^2, hq,
                    " (",
                    italic(P), hp,
                    "), ",
                    italic(I)^2, hi,
                    tau^2, ht),
              list(hl = resid.hetlab,
                   df = df.Q.resid,
                   hq = hetstat.Q.resid,
                   hp = hetstat.pval.Q.resid,
                   hi = hetstat.I2.resid,
                   ht = hetstat.tau2.resid)
            )
        else
          hetstat.resid <-
            substitute(
              paste(hl,
                    chi[df]^2, hq,
                    " (",
                    italic(P), hp,
                    "), ",
                    italic(I)^2, hi),
              list(hl = resid.hetlab,
                   df = df.Q.resid,
                   hq = hetstat.Q.resid,
                   hp = hetstat.pval.Q.resid,
                   hi = hetstat.I2.resid)
            )
      }
    }
    else if (revman5) {
      if (print.tau)
        hetstat.resid <-
          substitute(
            paste(hl,
                  "Tau", ht, "; ",
                  "Chi"^2, hq,
                  " (", P, hp, "); ",
                  I^2, hi),
            list(hl = resid.hetlab,
                 ht = hetstat.tau.resid,
                 hq = hetstat.Q.resid,
                 hp = hetstat.pval.Q.resid,
                 hi = hetstat.I2.resid)
            )
      else
        hetstat.resid <-
          substitute(
            paste(hl,
                  "Tau"^2, ht, "; ",
                  "Chi"^2, hq,
                  " (", P, hp, "); ",
                  I^2, hi),
            list(hl = resid.hetlab,
                 ht = hetstat.tau2.resid,
                 hq = hetstat.Q.resid,
                 hp = hetstat.pval.Q.resid,
                 hi = hetstat.I2.resid)
          )
    }
    else {
      #
      # One
      #
      if (print.I2 & !print.tau2.tau.resid & !print.Q & !print.pval.Q &
          !print.Rb.resid)
        hetstat.resid <-
          substitute(paste(hl, italic(I)^2, hi),
                     list(hl = resid.hetlab, hi = hetstat.I2.resid))
      else if (!print.I2 & print.tau2.tau.resid & !print.Q & !print.pval.Q &
               !print.Rb.resid)
        if (print.tau2)
          hetstat.resid <-
            substitute(paste(hl, tau^2, ht),
                       list(hl = resid.hetlab, ht = hetstat.tau2.resid))
        else
          hetstat.resid <-
            substitute(paste(hl, tau, ht),
                       list(hl = resid.hetlab, ht = hetstat.tau.resid))
      else if (!print.I2 & !print.tau2.tau.resid & print.Q & !print.pval.Q &
               !print.Rb.resid)
        hetstat.resid <-
          substitute(paste(hl, chi[df]^2, hq),
                     list(hl = resid.hetlab, df = df.Q.resid,
                          hq = hetstat.Q.resid))
      else if (!print.I2 & !print.tau2.tau.resid & !print.Q & print.pval.Q &
               !print.Rb.resid)
        hetstat.resid <-
          substitute(paste(hl, italic(p), hp),
                     list(hl = resid.hetlab, hp = hetstat.pval.Q.resid))
      else if (!print.I2 & !print.tau2.tau.resid & !print.Q & !print.pval.Q &
               print.Rb.resid)
        hetstat.resid <-
          substitute(paste(hl, italic(R)[italic(b)], hb),
                     list(hl = resid.hetlab, hb = hetstat.Rb.resid))
      #
      # Two
      #
      else if (print.I2 & print.tau2.tau.resid & !print.Q & !print.pval.Q &
               !print.Rb.resid)
        if (print.tau2)
          hetstat.resid <-
            substitute(paste(hl, italic(I)^2, hi,
                             ", ",
                             tau^2, ht),
                       list(hl = resid.hetlab,
                            hi = hetstat.I2.resid, ht = hetstat.tau2.resid))
        else
          hetstat.resid <-
            substitute(paste(hl, italic(I)^2, hi,
                             ", ",
                             tau, ht),
                       list(hl = resid.hetlab,
                            hi = hetstat.I2.resid, ht = hetstat.tau.resid))
      else if (print.I2 & !print.tau2.tau.resid & print.Q & !print.pval.Q &
               !print.Rb.resid)
        hetstat.resid <-
          substitute(paste(hl, italic(I)^2, hi,
                           ", ",
                           chi[df]^2, hq),
                     list(hl = resid.hetlab, df = df.Q.resid,
                          hi = hetstat.I2.resid, hq = hetstat.Q.resid))
      else if (print.I2 & !print.tau2.tau.resid & !print.Q & print.pval.Q &
               !print.Rb.resid)
        hetstat.resid <-
          substitute(paste(hl, italic(I)^2, hi,
                           ", ",
                           italic(p), hp),
                     list(hl = resid.hetlab,
                          hi = hetstat.I2.resid, hp = hetstat.pval.Q.resid))
      else if (print.I2 & !print.tau2.tau.resid & !print.Q & !print.pval.Q &
               print.Rb.resid)
        hetstat.resid <-
          substitute(paste(hl, italic(I)^2, hi,
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = resid.hetlab,
                          hi = hetstat.I2.resid, hb = hetstat.Rb.resid))
      else if (!print.I2 & print.tau2.tau.resid & print.Q & !print.pval.Q &
               !print.Rb.resid)
        if (print.tau2)
          hetstat.resid <-
            substitute(paste(hl, tau^2, ht,
                             ", ",
                             chi[df]^2, hq),
                       list(hl = resid.hetlab, df = df.Q.resid,
                            ht = hetstat.tau2.resid, hq = hetstat.Q.resid))
        else
          hetstat.resid <-
            substitute(paste(hl, tau, ht,
                             ", ",
                             chi[df]^2, hq),
                       list(hl = resid.hetlab, df = df.Q.resid,
                            ht = hetstat.tau.resid, hq = hetstat.Q.resid))
      else if (!print.I2 & print.tau2.tau.resid & !print.Q & print.pval.Q &
               !print.Rb.resid)
        if (print.tau2)
          hetstat.resid <-
            substitute(paste(hl, tau^2, ht,
                             ", ",
                             italic(p), hp),
                       list(hl = resid.hetlab,
                            ht = hetstat.tau2.resid, hp = hetstat.pval.Q.resid))
        else
          hetstat.resid <-
            substitute(paste(hl, tau, ht,
                             ", ",
                             italic(p), hp),
                       list(hl = resid.hetlab,
                            ht = hetstat.tau.resid, hp = hetstat.pval.Q.resid))
      else if (!print.I2 & print.tau2.tau.resid & !print.Q & !print.pval.Q &
               print.Rb.resid)
        if (print.tau2)
          hetstat.resid <-
            substitute(paste(hl, tau^2, ht,
                             ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = resid.hetlab,
                            ht = hetstat.tau2.resid, hb = hetstat.Rb.resid))
        else
          hetstat.resid <-
            substitute(paste(hl, tau, ht,
                             ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = resid.hetlab,
                            ht = hetstat.tau.resid, hb = hetstat.Rb.resid))
      else if (!print.I2 & !print.tau2.tau.resid & print.Q & print.pval.Q &
               !print.Rb.resid)
        hetstat.resid <-
          substitute(paste(hl, chi[df]^2, hq,
                           " (",
                           italic(p), hp, ")"),
                     list(hl = resid.hetlab, df = df.Q.resid,
                          hq = hetstat.Q.resid, hp = hetstat.pval.Q.resid))
      else if (!print.I2 & !print.tau2.tau.resid & print.Q & !print.pval.Q &
               print.Rb.resid)
        hetstat.resid <-
          substitute(paste(hl, chi[df]^2, hq,
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = resid.hetlab, df = df.Q.resid,
                          hq = hetstat.Q.resid, hb = hetstat.Rb.resid))
      else if (!print.I2 & !print.tau2.tau.resid & !print.Q & print.pval.Q &
               print.Rb.resid)
        hetstat.resid <-
          substitute(paste(hl, italic(p), hp,
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = resid.hetlab,
                          hp = hetstat.pval.Q.resid, hb = hetstat.Rb.resid))
      #
      # Three
      #
      else if (print.I2 & print.tau2.tau.resid & print.Q & !print.pval.Q &
               !print.Rb.resid)
        if (print.tau2)
          hetstat.resid <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau^2, ht, ", ",
                             chi[df]^2, hq),
                       list(hl = resid.hetlab, df = df.Q.resid,
                            hi = hetstat.I2.resid, ht = hetstat.tau2.resid,
                            hq = hetstat.Q.resid))
        else
          hetstat.resid <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau, ht, ", ",
                             chi[df]^2, hq),
                       list(hl = resid.hetlab, df = df.Q.resid,
                            hi = hetstat.I2.resid, ht = hetstat.tau.resid,
                            hq = hetstat.Q.resid))
      else if (print.I2 & print.tau2.tau.resid & !print.Q & print.pval.Q &
               !print.Rb.resid)
        if (print.tau2)
          hetstat.resid <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau^2, ht, ", ",
                             italic(p), hp),
                       list(hl = resid.hetlab,
                            hi = hetstat.I2.resid, ht = hetstat.tau2.resid,
                            hp = hetstat.pval.Q.resid))
        else
          hetstat.resid <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau, ht, ", ",
                             italic(p), hp),
                       list(hl = resid.hetlab,
                            hi = hetstat.I2.resid, ht = hetstat.tau.resid,
                            hp = hetstat.pval.Q.resid))
      else if (print.I2 & !print.tau2.tau.resid & print.Q & print.pval.Q &
               !print.Rb.resid)
        hetstat.resid <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           chi[df]^2, hq,
                           " (", italic(p), hp, ")"),
                     list(hl = resid.hetlab, df = df.Q.resid,
                          hi = hetstat.I2.resid,
                          hq = hetstat.Q.resid, hp = hetstat.pval.Q.resid))
      else if (!print.I2 & print.tau2.tau.resid & print.Q & print.pval.Q &
               !print.Rb.resid)
        if (print.tau2)
          hetstat.resid <-
            substitute(paste(hl,
                             tau^2, ht, ", ",
                             chi[df]^2, hq,
                             " (", italic(p), hp, ")"),
                       list(hl = resid.hetlab, df = df.Q.resid,
                            ht = hetstat.tau2.resid,
                            hq = hetstat.Q.resid, hp = hetstat.pval.Q.resid))
        else
          hetstat.resid <-
            substitute(paste(hl,
                             tau, ht, ", ",
                             chi[df]^2, hq,
                             " (", italic(p), hp, ")"),
                       list(hl = resid.hetlab, df = df.Q.resid,
                            ht = hetstat.tau.resid,
                            hq = hetstat.Q.resid, hp = hetstat.pval.Q.resid))
      else if (print.I2 & print.tau2.tau.resid & !print.Q & !print.pval.Q &
               print.Rb.resid)
        if (print.tau2)
          hetstat.resid <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau^2, ht, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = resid.hetlab,
                            hi = hetstat.I2.resid, ht = hetstat.tau2.resid,
                            hb = hetstat.Rb.resid))
        else
          hetstat.resid <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau, ht, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = resid.hetlab,
                            hi = hetstat.I2.resid, ht = hetstat.tau.resid,
                            hb = hetstat.Rb.resid))
      else if (print.I2 & !print.tau2.tau.resid & print.Q & !print.pval.Q &
               print.Rb.resid)
        hetstat.resid <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           chi[df]^2, hq, ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = resid.hetlab, df = df.Q.resid,
                          hi = hetstat.I2.resid,
                          hq = hetstat.Q.resid,
                          hb = hetstat.Rb.resid))
      else if (!print.I2 & print.tau2.tau.resid & print.Q & !print.pval.Q &
               print.Rb.resid)
        if (print.tau2)
          hetstat.resid <-
            substitute(paste(hl,
                             tau^2, ht, ", ",
                             chi[df]^2, hq, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = resid.hetlab, df = df.Q.resid,
                            ht = hetstat.tau2.resid,
                            hq = hetstat.Q.resid,
                            hb = hetstat.Rb.resid))
        else
          hetstat.resid <-
            substitute(paste(hl,
                             tau, ht, ", ",
                             chi[df]^2, hq, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = resid.hetlab, df = df.Q.resid,
                            ht = hetstat.tau.resid,
                            hq = hetstat.Q.resid,
                            hb = hetstat.Rb.resid))          
      else if (print.I2 & !print.tau2.tau.resid & !print.Q & print.pval.Q &
               print.Rb.resid)
        hetstat.resid <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           italic(p), hp, ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = resid.hetlab,
                          hi = hetstat.I2.resid, hp = hetstat.pval.Q.resid,
                          hb = hetstat.Rb.resid))
      else if (!print.I2 & print.tau2.tau.resid & !print.Q & print.pval.Q &
               print.Rb.resid)
        if (print.tau2)
          hetstat.resid <-
            substitute(paste(hl,
                             tau^2, ht, ", ",
                             italic(p), hp, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = resid.hetlab,
                            ht = hetstat.tau2.resid,
                            hp = hetstat.pval.Q.resid,
                            hb = hetstat.Rb.resid))
        else
          hetstat.resid <-
            substitute(paste(hl,
                             tau, ht, ", ",
                             italic(p), hp, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = resid.hetlab,
                            ht = hetstat.tau.resid,
                            hp = hetstat.pval.Q.resid,
                            hb = hetstat.Rb.resid))          
      else if (!print.I2 & !print.tau2.tau.resid & print.Q & print.pval.Q &
               print.Rb.resid)
        hetstat.resid <-
          substitute(paste(hl,
                           chi[df]^2, hq,
                           " (", italic(p), hp, ")",
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = resid.hetlab, df = df.Q.resid,
                          hq = hetstat.Q.resid, hp = hetstat.pval.Q.resid,
                          hb = hetstat.Rb.resid))      
      #
      # Four
      #
      if (print.I2 & print.tau2.tau.resid & print.Q & print.pval.Q &
          !print.Rb.resid)
        if (print.tau2)
          hetstat.resid <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau^2, ht, ", ",
                             chi[df]^2, hq,
                             " (", italic(p), hp, ")"),
                       list(hl = resid.hetlab, df = df.Q.resid,
                            hi = hetstat.I2.resid, ht = hetstat.tau2.resid,
                            hq = hetstat.Q.resid, hp = hetstat.pval.Q.resid))
        else
          hetstat.resid <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau, ht, ", ",
                             chi[df]^2, hq,
                             " (", italic(p), hp, ")"),
                       list(hl = resid.hetlab, df = df.Q.resid,
                            hi = hetstat.I2.resid, ht = hetstat.tau.resid,
                            hq = hetstat.Q.resid, hp = hetstat.pval.Q.resid))          
      else if (print.I2 & print.tau2.tau.resid & print.Q & !print.pval.Q &
               print.Rb.resid)
        if (print.tau2)
          hetstat.resid <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau^2, ht, ", ",
                             chi[df]^2, hq, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = resid.hetlab, df = df.Q.resid,
                            hi = hetstat.I2.resid, ht = hetstat.tau2.resid,
                            hq = hetstat.Q.resid, hb = hetstat.Rb.resid))
        else
          hetstat.resid <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau, ht, ", ",
                             chi[df]^2, hq, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = resid.hetlab, df = df.Q.resid,
                            hi = hetstat.I2.resid, ht = hetstat.tau.resid,
                            hq = hetstat.Q.resid, hb = hetstat.Rb.resid))          
      else if (print.I2 & print.tau2.tau.resid & !print.Q & print.pval.Q &
               print.Rb.resid)
        if (print.tau2)
          hetstat.resid <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau^2, ht, ", ",
                             italic(p), hp, ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = resid.hetlab,
                            hi = hetstat.I2.resid, ht = hetstat.tau2.resid,
                            hp = hetstat.pval.Q.resid, hb = hetstat.Rb.resid))
        else
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           tau, ht, ", ",
                           italic(p), hp, ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = resid.hetlab,
                          hi = hetstat.I2.resid, ht = hetstat.tau.resid,
                          hp = hetstat.pval.Q.resid, hb = hetstat.Rb.resid))
      else if (print.I2 & !print.tau2.tau.resid & print.Q & !print.pval.Q &
               print.Rb.resid)
        hetstat.resid <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           chi[df]^2, hq,
                           " (", italic(p), hp, ")",
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = resid.hetlab, df = df.Q.resid,
                          ht = hetstat.tau2.resid,
                          hq = hetstat.Q.resid, hp = hetstat.pval.Q.resid,
                          hb = hetstat.Rb.resid))
      else if (!print.I2 & print.tau2.tau.resid & print.Q & print.pval.Q &
               print.Rb.resid)
        if (print.tau2)
          hetstat.resid <-
            substitute(paste(hl,
                             tau^2, ht, ", ",
                             chi[df]^2, hq,
                             " (", italic(p), hp, ")",
                             ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = resid.hetlab, df = df.Q.resid,
                            ht = hetstat.tau2.resid,
                            hq = hetstat.Q.resid, hp = hetstat.pval.Q.resid,
                            hb = hetstat.Rb.resid))
        else
          hetstat.resid <-
            substitute(paste(hl,
                             tau, ht, ", ",
                             chi[df]^2, hq,
                             " (", italic(p), hp, ")",
                             ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = resid.hetlab, df = df.Q.resid,
                            ht = hetstat.tau.resid,
                            hq = hetstat.Q.resid, hp = hetstat.pval.Q.resid,
                            hb = hetstat.Rb.resid))          
      #
      # Five
      #
      else if (print.I2 & print.tau2.tau.resid & print.Q & print.pval.Q &
               print.Rb.resid)
        if (print.tau2)
          hetstat.resid <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau^2, ht, ", ",
                             chi[df]^2, hq,
                             " (", italic(p), hp, ")",
                             ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = resid.hetlab, df = df.Q.resid,
                            hi = hetstat.I2.resid, ht = hetstat.tau2.resid,
                            hq = hetstat.Q.resid, hp = hetstat.pval.Q.resid,
                            hb = hetstat.Rb.resid))
        else
          hetstat.resid <-
            substitute(paste(hl,
                             italic(I)^2, hi, ", ",
                             tau, ht, ", ",
                             chi[df]^2, hq,
                             " (", italic(p), hp, ")",
                             ", ",
                             italic(R)[italic(b)], hb),
                       list(hl = resid.hetlab, df = df.Q.resid,
                            hi = hetstat.I2.resid, ht = hetstat.tau.resid,
                            hq = hetstat.Q.resid, hp = hetstat.pval.Q.resid,
                            hb = hetstat.Rb.resid))          
    }
  }
  #
  # Label of test for overall effect
  #
  #
  pvals.overall <- formatPT(c(x$pval.common, x$pval.random),
                            lab = TRUE, labval = "",
                            digits = digits.pval,
                            zero = zero.pval, JAMA = JAMA.pval,
                            scientific = scientific.pval,
                            lab.NA = lab.NA)
  statistics.overall <- formatN(c(x$statistic.common, x$statistic.random),
                                digits.stat, lab.NA, big.mark = big.mark)
  #
  # Remove superfluous spaces
  #
  pvals.overall <- rmSpace(pvals.overall, end = TRUE)
  #
  while(any(grepl("  ", pvals.overall)))
    pvals.overall <- gsub("  ", " ", pvals.overall)
  while(any(grepl("  ", statistics.overall)))
    statistics.overall <- gsub("  ", " ", statistics.overall)
  #
  if (test.overall.common) {
    if (print.stat) {
      if (bmj)
        text.overall.common <-
          substitute(paste(tl,
                           Z, hetseparator, tt,
                           ", P", tp),
                     list(tl = label.test.overall.common,
                          hetseparator = hetseparator,
                          tt = statistics.overall[1],
                          tp = pvals.overall[1]))
      else if (revman5)
        text.overall.common <-
          substitute(paste(tl,
                           Z, hetseparator, tt,
                           " (P", tp, ")"),
                     list(tl = label.test.overall.common,
                          hetseparator = hetseparator,
                          tt = statistics.overall[1],
                          tp = pvals.overall[1]))
      else if (jama)
        text.overall.common <-
          substitute(paste(tl,
                           italic(z), hetseparator, tt,
                           " (", italic(P), tp, ")"),
                     list(tl = label.test.overall.common,
                          hetseparator = hetseparator,
                          tt = statistics.overall[1],
                          tp = pvals.overall[1]))
      else
        text.overall.common <-
          substitute(paste(tl,
                           italic(z), hetseparator, tt,
                           " (", italic(p), tp, ")"),
                     list(tl = label.test.overall.common,
                          hetseparator = hetseparator,
                          tt = statistics.overall[1],
                          tp = pvals.overall[1]))
    }
    else {
      if (bmj)
        text.overall.common <-
          substitute(paste(tl, " P", tp),
                     list(tl = label.test.overall.common,
                          tp = pvals.overall[1]))
      else if (revman5)
        text.overall.common <-
          substitute(paste(tl, " P", tp),
                     list(tl = label.test.overall.common,
                          tp = pvals.overall[1]))
      else if (jama)
        text.overall.common <-
          substitute(paste(tl, " ", italic(P), tp),
                     list(tl = label.test.overall.common,
                          tp = pvals.overall[1]))
      else
        text.overall.common <-
          substitute(paste(tl, " ", italic(p), tp),
                     list(tl = label.test.overall.common,
                          tp = pvals.overall[1]))     
    }
  }
  else
    text.overall.common <- ""
  #
  hakn.kero <- x$method.random.ci[1] %in% c("HK", "KR")
  #
  if (test.overall.random) {
    if (print.stat) {
      if (!hakn.kero) {
        if (bmj)
          text.overall.random <-
            substitute(paste(tl,
                             Z, hetseparator, tt,
                             ", P", tp),
                       list(tl = label.test.overall.random,
                            hetseparator = hetseparator,
                            tt = statistics.overall[2],
                            tp = pvals.overall[2]))
        else if (revman5)
          text.overall.random <-
            substitute(paste(tl,
                             Z, hetseparator, tt,
                             " (P", tp, ")"),
                       list(tl = label.test.overall.random,
                            hetseparator = hetseparator,
                            tt = statistics.overall[2],
                            tp = pvals.overall[2]))
        else if (jama)
          text.overall.random <-
            substitute(paste(tl,
                             italic(z), hetseparator, tt,
                             " (", italic(P), tp, ")"),
                       list(tl = label.test.overall.random,
                            hetseparator = hetseparator,
                            tt = statistics.overall[2],
                            tp = pvals.overall[2]))
        else
          text.overall.random <-
            substitute(paste(tl,
                             italic(z), hetseparator, tt,
                             " (", italic(p), tp, ")"),
                       list(tl = label.test.overall.random,
                            hetseparator = hetseparator,
                            tt = statistics.overall[2],
                            tp = pvals.overall[2]))
      }
      else {
        if (bmj)
          text.overall.random <-
            substitute(paste(tl,
                             t[df], hetseparator, tt,
                             ", P", tp),
                       list(tl = label.test.overall.random,
                            hetseparator = hetseparator,
                            tt = statistics.overall[2],
                            tp = pvals.overall[2],
                            df = round(x$df.random, 1)))
        else if (revman5)
          text.overall.random <-
            substitute(paste(tl,
                             t[df], hetseparator, tt,
                             " (P", tp, ")"),
                       list(tl = label.test.overall.random,
                            hetseparator = hetseparator,
                            tt = statistics.overall[2],
                            tp = pvals.overall[2],
                            df = round(x$df.random, 1)))
        else if (jama)
          text.overall.random <-
            substitute(paste(tl,
                             italic(t)[df], hetseparator, tt,
                             " (", italic(P), tp, ")"),
                       list(tl = label.test.overall.random,
                            hetseparator = hetseparator,
                            tt = statistics.overall[2],
                            tp = pvals.overall[2],
                            df = round(x$df.random, 1)))
        else
          text.overall.random <-
            substitute(paste(tl,
                             italic(t)[df], hetseparator, tt,
                             " (", italic(p), tp, ")"),
                       list(tl = label.test.overall.random,
                            hetseparator = hetseparator,
                            tt = statistics.overall[2],
                            tp = pvals.overall[2],
                            df = round(x$df.random, 1)))
      }
    }
    else {
      if (bmj)
        text.overall.random <-
          substitute(paste(tl, " P", tp),
                     list(tl = label.test.overall.random,
                          tp = pvals.overall[2]))
      else if (revman5)
        text.overall.random <-
          substitute(paste(tl, " P", tp),
                     list(tl = label.test.overall.random,
                          tp = pvals.overall[2]))
      else if (jama)
        text.overall.random <-
          substitute(paste(tl, " ", italic(P), tp),
                     list(tl = label.test.overall.random,
                          tp = pvals.overall[2]))
      else
        text.overall.random <-
          substitute(paste(tl, " ", italic(p), tp),
                     list(tl = label.test.overall.random,
                          tp = pvals.overall[2]))
    }
  }
  else
    text.overall.random <- ""
  #
  #
  # Label of test for subgroup differences
  #
  # if (by) {
  #   Q.bs <- c(Q.b.common, Q.b.random)
  #   pval.Q.bs <- c(pval.Q.b.common, pval.Q.b.random)
  # }
  # else {
  #   Q.bs <- NA
  #   pval.Q.bs <- NA
  # }
  #
  # hetstat.Q.bs <-
  #   paste0(hetseparator,
  #          gsub(" ", "", formatN(Q.bs, digits.Q, "NA", big.mark = big.mark)),
  #          if (!jama) ", df",
  #          if (!jama) hetseparator,
  #          if (!jama) c(df.Q.b.common,
  #                       if (length(df.Q.b.random) == 2)
  #                         paste(df.Q.b.random, collapse = ", ")
  #                       else
  #                         df.Q.b.random))
  #
  
  #
  if (by)
    hetstat.Q.b.common <- hetstat.pval.Q.b.common <-
      hetstat.Q.b.random <- hetstat.pval.Q.b.random <- NULL
  else
    hetstat.Q.b.common <- hetstat.pval.Q.b.common <-
      hetstat.Q.b.random <- hetstat.pval.Q.b.random <- NA
  #
  for (i in seq_along(Q.b.common)) {
    hetstat.i <-
      paste0(hetseparator,
             gsub(" ", "",
                  formatN(Q.b.common[i], digits.Q, "NA", big.mark = big.mark)),
             if (!jama) ", df",
             if (!jama) hetseparator,
             if (!jama) df.Q.b.common[i])
    #
    pval.i <-
      paste0(formatPT(pval.Q.b.common[i],
                      lab = TRUE, labval = "",
                      digits = digits.pval.Q,
                      zero = zero.pval, JAMA = JAMA.pval,
                      scientific = scientific.pval,
                      lab.NA = "NA"))
    #
    # Remove superfluous spaces
    #
    hetstat.i <- rmSpace(hetstat.i, end = TRUE)
    pval.i <- rmSpace(pval.i, end = TRUE)
    #
    while(any(grepl("  ", hetstat.i)))
      hetstat.i <- gsub("  ", " ", hetstat.i)
    while(any(grepl("  ", pval.i)))
      pval.i <- gsub("  ", " ", pval.i)
    #
    hetstat.Q.b.common <- c(hetstat.Q.b.common, hetstat.i)
    hetstat.pval.Q.b.common <- c(hetstat.pval.Q.b.common, pval.i)
  }
  #
  for (i in seq_along(Q.b.random)) {
    hetstat.i <-
      paste0(hetseparator,
             gsub(" ", "",
                  formatN(Q.b.random[i], digits.Q, "NA", big.mark = big.mark)),
             if (!jama) ", df",
             if (!jama) hetseparator,
             if (!jama) collapse(df.Q.b.random[[i]], quote = ""))
    #
    pval.i <-
      paste0(formatPT(pval.Q.b.random[i],
                      lab = TRUE, labval = "",
                      digits = digits.pval.Q,
                      zero = zero.pval, JAMA = JAMA.pval,
                      scientific = scientific.pval,
                      lab.NA = "NA"))
    #
    # Remove superfluous spaces
    #
    hetstat.i <- rmSpace(hetstat.i, end = TRUE)
    pval.i <- rmSpace(pval.i, end = TRUE)
    #
    while(any(grepl("  ", hetstat.i)))
      hetstat.i <- gsub("  ", " ", hetstat.i)
    while(any(grepl("  ", pval.i)))
      pval.i <- gsub("  ", " ", pval.i)
    #
    hetstat.Q.b.random <- c(hetstat.Q.b.random, hetstat.i)
    hetstat.pval.Q.b.random <- c(hetstat.pval.Q.b.random, pval.i)
  }
  #
  if (test.subgroup.common) {
    if (print.Q.subgroup) {
      if (bmj)
        text.subgroup.common <-
          substitute(paste(tl,
                           chi^2, tq,
                           ", P", tp),
                     list(tl = label.test.subgroup.common,
                          tq = hetstat.Q.b.common[1],
                          tp = hetstat.pval.Q.b.common[1]))
      else if (revman5)
        text.subgroup.common <-
          substitute(paste(tl,
                           "Chi"^2, tq,
                           " (P", tp, ")"),
                     list(tl = label.test.subgroup.common,
                          tq = hetstat.Q.b.common[1],
                          tp = hetstat.pval.Q.b.common[1]))
      else if (jama)
        text.subgroup.common <-
          substitute(paste(tl,
                           chi[df]^2, tq,
                           " (", italic(P), tp, ")"),
                     list(tl = label.test.subgroup.common,
                          tq = hetstat.Q.b.common[1],
                          tp = hetstat.pval.Q.b.common[1],
                          df = df.Q.b.common[1]))
      else
        text.subgroup.common <-
          substitute(paste(tl,
                           chi[df]^2, tq,
                           " (", italic(p), tp, ")"),
                     list(tl = label.test.subgroup.common,
                          tq = hetstat.Q.b.common[1],
                          tp = hetstat.pval.Q.b.common[1],
                          df = df.Q.b.common[1]))
    }
    else {
      if (bmj)
        text.subgroup.common <-
          substitute(paste(tl, " P", tp),
                     list(tl = label.test.subgroup.common,
                          tp = hetstat.pval.Q.b.common[1]))
      else if (revman5)
        text.subgroup.common <-
          substitute(paste(tl, " P", tp),
                     list(tl = label.test.subgroup.common,
                          tp = hetstat.pval.Q.b.common[1]))
      else if (jama)
        text.subgroup.common <-
          substitute(paste(tl, " ", italic(P), tp),
                     list(tl = label.test.subgroup.common,
                          tp = hetstat.pval.Q.b.common[1]))
      else
        text.subgroup.common <-
          substitute(paste(tl, " ", italic(p), tp),
                     list(tl = label.test.subgroup.common,
                          tp = hetstat.pval.Q.b.common[1]))
    }
  }
  else
    text.subgroup.common <- ""
  #
  if (test.subgroup.random) {
    if (print.Q.subgroup) {
      if (bmj) {
        if (length(df.Q.b.random[[1]]) == 2)
          text.subgroup.random <-
            substitute(paste(tl, "F", tq, ", P", tp),
                       list(tl = label.test.subgroup.random,
                            tq = hetstat.Q.b.random[1],
                            tp = hetstat.pval.Q.b.random[1]))
        else
          text.subgroup.random <-
            substitute(paste(tl, chi^2, tq, " P", tp),
                       list(tl = label.test.subgroup.random,
                            tq = hetstat.Q.b.random[1],
                            tp = hetstat.pval.Q.b.random[1]))
      }
      else if (revman5) {
        if (length(df.Q.b.random[[1]]) == 2)
          text.subgroup.random <-
            substitute(paste(tl,
                             "F", tq,
                             " (P", tp, ")"),
                       list(tl = label.test.subgroup.random,
                            tq = hetstat.Q.b.random[1],
                            tp = hetstat.pval.Q.b.random[1]))
        else
          text.subgroup.random <-
            substitute(paste(tl,
                             "Chi"^2, tq,
                             " (P", tp, ")"),
                       list(tl = label.test.subgroup.random,
                            tq = hetstat.Q.b.random[1],
                            tp = hetstat.pval.Q.b.random[1]))
      }
      else if (jama) {
        if (length(df.Q.b.random[[1]]) == 2)
          text.subgroup.random <-
            substitute(paste(tl,
                             F[df], tq,
                             " (", italic(P), tp, ")"),
                       list(tl = label.test.subgroup.random,
                            tq = hetstat.Q.b.random[1],
                            tp = hetstat.pval.Q.b.random[1],
                            df = paste(df.Q.b.random[[1]], collapse = ", ")))
        else
          text.subgroup.random <-
            substitute(paste(tl,
                             chi[df]^2, tq,
                             " (", italic(P), tp, ")"),
                       list(tl = label.test.subgroup.random,
                            tq = hetstat.Q.b.random[1],
                            tp = hetstat.pval.Q.b.random[1],
                            df = df.Q.b.random[[1]]))          
      }
      else {
        if (length(df.Q.b.random[[1]]) == 2)
          text.subgroup.random <-
            substitute(paste(tl,
                             F[df], tq,
                             " (", italic(p), tp, ")"),
                       list(tl = label.test.subgroup.random,
                            tq = hetstat.Q.b.random[1],
                            tp = hetstat.pval.Q.b.random[1],
                            df = paste(df.Q.b.random[[1]], collapse = ", ")))
        else
          text.subgroup.random <-
            substitute(paste(tl,
                             chi[df]^2, tq,
                             " (", italic(p), tp, ")"),
                       list(tl = label.test.subgroup.random,
                            tq = hetstat.Q.b.random[1],
                            tp = hetstat.pval.Q.b.random[1],
                            df = df.Q.b.random[[1]]))
      }
    }
    else {
      if (bmj)
        text.subgroup.random <-
          substitute(paste(tl, " P", tp),
                     list(tl = label.test.subgroup.random,
                          tp = hetstat.pval.Q.b.random[1]))
      else if (revman5)
        text.subgroup.random <-
          substitute(paste(tl, " P", tp),
                     list(tl = label.test.subgroup.random,
                          tp = hetstat.pval.Q.b.random[1]))
      else if (jama)
        text.subgroup.random <-
          substitute(paste(tl, " ", italic(P), tp),
                     list(tl = label.test.subgroup.random,
                          tp = hetstat.pval.Q.b.random[1]))
      else
        text.subgroup.random <-
          substitute(paste(tl, " ", italic(p), tp),
                     list(tl = label.test.subgroup.random,
                          tp = hetstat.pval.Q.b.random[1]))
    }
  }
  else
    text.subgroup.random <- ""
  
  
  #
  #
  # (7) Prepare data for subgroup analysis
  #
  #
  NAs <- rep(NA, n.com + n.ran + n.prd)
  NAs.com <- rep(NA, n.com)
  NAs.ran <- rep(NA, n.ran)
  NAs.prd <- rep(NA, n.prd)
  NAs.com1 <- rep(NA, n.com - 1)
  NAs.ran1 <- rep(NA, n.ran - 1)
  #
  blanks <- rep("", n.com + n.ran + n.prd)
  blanks.prd <- rep("", n.prd)
  #
  emp.com <- 1 + seq_len(n.com - 1)
  emp.ran <- n.com + 1 + seq_len(n.ran - 1)
  #
  all.com <- seq_len(n.com)
  all.ran <- n.com + seq_len(n.ran)
  all.prd <- n.com + n.ran + seq_len(n.prd)
  all.res <- c(all.com, all.ran, all.prd)
  #
  if (by) {
    NAs.by <- rep(NA, n.by)
    #
    k.w.hetstat <- if (metabind) x$k.w.orig else x$k.w
    o.w <- order(factor(x$subgroup.levels, levels = subgroup.levels))
    k.w.hetstat <- k.w.hetstat[o.w]
    if (hakn.kero)
      df.random.w <- x$df.random.w[o.w]
    #
    TE.common.w <- ordermat(TE.common.w, subgroup.levels)
    lower.common.w <- ordermat(lower.common.w, subgroup.levels)
    upper.common.w <- ordermat(upper.common.w, subgroup.levels)
    statistic.common.w <- ordermat(statistic.common.w, subgroup.levels)
    pval.common.w <- ordermat(pval.common.w, subgroup.levels)
    #
    TE.common.w <- list2vec(TE.common.w)
    lower.common.w <- list2vec(lower.common.w)
    upper.common.w <- list2vec(upper.common.w)
    statistic.common.w <- list2vec(statistic.common.w)
    pval.common.w <- list2vec(pval.common.w)
    #
    TE.random.w <- ordermat(TE.random.w, subgroup.levels)
    lower.random.w <- ordermat(lower.random.w, subgroup.levels)
    upper.random.w <- ordermat(upper.random.w, subgroup.levels)
    statistic.random.w <- ordermat(statistic.random.w, subgroup.levels)
    pval.random.w <- ordermat(pval.random.w, subgroup.levels)
    #
    TE.random.w <- list2vec(TE.random.w)
    lower.random.w <- list2vec(lower.random.w)
    upper.random.w <- list2vec(upper.random.w)
    statistic.random.w <- list2vec(statistic.random.w)
    pval.random.w <- list2vec(pval.random.w)
    #
    lower.predict.w <- ordermat(lower.predict.w, subgroup.levels)
    upper.predict.w <- ordermat(upper.predict.w, subgroup.levels)
    #
    lower.predict.w <- list2vec(lower.predict.w)
    upper.predict.w <- list2vec(upper.predict.w)
    #
    Q.w <- x$Q.w[o.w]
    pval.Q.w <- x$pval.Q.w[o.w]
    #
    I2.w <- x$I2.w[o.w]
    lowI2.w <- x$lower.I2.w[o.w]
    uppI2.w <- x$upper.I2.w[o.w]
    #
    Rb.w <- x$Rb.w[o.w]
    lowRb.w <- x$lower.Rb.w[o.w]
    uppRb.w <- x$upper.Rb.w[o.w]
    #
    if (is.list(x$tau2.w))
      x$tau2.w <- x$tau2.w[[1]]
    if (is.list(x$lower.tau2.w))
      x$lower.tau2.w <- x$lower.tau2.w[[1]]
    if (is.list(x$upper.tau2.w))
      x$upper.tau2.w <- x$upper.tau2.w[[1]]
    #
    if (is.list(x$tau.w))
      x$tau.w <- x$tau.w[[1]]
    if (is.list(x$lower.tau.w))
      x$lower.tau.w <- x$lower.tau.w[[1]]
    if (is.list(x$upper.tau.w))
      x$upper.tau.w <- x$upper.tau.w[[1]]
    #
    tau2.w <- x$tau2.w[o.w]
    lower.tau2.w <- x$lower.tau2.w[o.w]
    upper.tau2.w <- x$upper.tau2.w[o.w]
    #
    tau.w <- x$tau.w[o.w]
    lower.tau.w <- x$lower.tau.w[o.w]
    upper.tau.w <- x$upper.tau.w[o.w]
    #
    w.common.w <- collapsemat(x$w.common.w)
    if (is.list(w.common.w))
      w.common.w <- w.common.w[[1]]
    #
    w.random.w <- collapsemat(x$w.random.w)
    if (is.list(w.random.w))
      w.random.w <- w.random.w[[1]]
    #
    w.common.w <- collapsemat(w.common.w)[o.w]
    w.random.w <- collapsemat(w.random.w)[o.w]
    e.e.w <- if (metaprop | metarate) x$event.w[o.w] else x$event.e.w[o.w]
    t.e.w <- if (metainc | metarate) x$time.e.w[o.w] else NAs.by
    n.e.w <-
      if (metacor | metaprop | metamean | metarate)
        x$n.w[o.w]
      else
        x$n.e.w[o.w]
    e.c.w <- x$event.c.w[o.w]
    t.c.w <- if (metainc) x$time.c.w[o.w] else NAs.by
    n.c.w <- x$n.c.w[o.w]
    n.harmonic.mean.w <- x$n.harmonic.mean.w[o.w]
    t.harmonic.mean.w <- x$t.harmonic.mean.w[o.w]
    #
    k.all.w <- x$k.all.w[o.w]
    k.study.w <- x$k.study.w[o.w]
    k.w <- x$k.w[o.w]
    k.TE.w <- x$k.TE.w[o.w]
    #
    subgroup.logical <- subgroup.logical[o.w]
    subgroup.hetstat.logical <- subgroup.hetstat.logical[o.w]
    #
    common.subgroup.logical <- common.subgroup.logical[o.w]
    random.subgroup.logical <- random.subgroup.logical[o.w]
    prediction.subgroup.logical <- prediction.subgroup.logical[o.w]
    #
    # Do (not) drop subgroups without studies in subgroup
    # meta-analysis
    #
    if (allstudies)
      sel.w <- rep(TRUE, length(k.all.w))
    else
      sel.w <- k.w > 0
    #
    if (hakn.kero)
      df.random.w <- df.random.w[sel.w]
    #
    TE.common.w <- TE.common.w[repl(sel.w, n.com, n.by)]
    lower.common.w <- lower.common.w[repl(sel.w, n.com, n.by)]
    upper.common.w <- upper.common.w[repl(sel.w, n.com, n.by)]
    statistic.common.w <- statistic.common.w[repl(sel.w, n.com, n.by)]
    pval.common.w <- pval.common.w[repl(sel.w, n.com, n.by)]
    #
    TE.random.w <- TE.random.w[repl(sel.w, n.ran, n.by)]
    lower.random.w <- lower.random.w[repl(sel.w, n.ran, n.by)]
    upper.random.w <- upper.random.w[repl(sel.w, n.ran, n.by)]
    statistic.random.w <- statistic.random.w[repl(sel.w, n.ran, n.by)]
    pval.random.w <- pval.random.w[repl(sel.w, n.ran, n.by)]
    #
    lower.predict.w <- lower.predict.w[repl(sel.w, n.prd, n.by)]
    upper.predict.w <- upper.predict.w[repl(sel.w, n.prd, n.by)]
    #
    Q.w <- Q.w[sel.w]
    pval.Q.w <- pval.Q.w[sel.w]
    #
    I2.w    <- I2.w[sel.w]
    lowI2.w <- lowI2.w[sel.w]
    uppI2.w <- uppI2.w[sel.w]
    #
    Rb.w    <- Rb.w[sel.w]
    lowRb.w <- lowRb.w[sel.w]
    uppRb.w <- uppRb.w[sel.w]
    #
    tau2.w <- tau2.w[sel.w]
    lower.tau2.w <- lower.tau2.w[sel.w]
    upper.tau2.w <- upper.tau2.w[sel.w]
    #
    tau.w <- tau.w[sel.w]
    lower.tau.w <- lower.tau.w[sel.w]
    upper.tau.w <- upper.tau.w[sel.w]
    #
    w.common.w <- w.common.w[sel.w]
    w.random.w <- w.random.w[sel.w]
    e.e.w <- e.e.w[sel.w]
    t.e.w <- t.e.w[sel.w]
    n.e.w <- n.e.w[sel.w]
    e.c.w <- e.c.w[sel.w]
    t.c.w <- t.c.w[sel.w]
    n.c.w <- n.c.w[sel.w]
    n.harmonic.mean.w <- n.harmonic.mean.w[sel.w]
    t.harmonic.mean.w <- t.harmonic.mean.w[sel.w]
    #
    k.all.w <- k.all.w[sel.w]
    k.study.w <- k.study.w[sel.w]
    k.w <- k.w[sel.w]
    k.TE.w <- k.TE.w[sel.w]
    #
    subgroup.logical <- subgroup.logical[sel.w]
    subgroup.hetstat.logical <- subgroup.hetstat.logical[sel.w]
    #
    common.subgroup.logical <- common.subgroup.logical[sel.w]
    random.subgroup.logical <- random.subgroup.logical[sel.w]
    prediction.subgroup.logical <- prediction.subgroup.logical[sel.w]
    #
    subgroup.levels <- subgroup.levels[sel.w]
    #
    n.by <- length(subgroup.levels)
    n.com.w <- n.com * n.by
    n.ran.w <- n.ran * n.by
    n.prd.w <- n.prd * n.by
    n.stat.w <- 3 * n.by
    #
    NAs.by <- rep(NA, n.by)
    NAs.com.w <- rep(NA, n.com.w)
    NAs.ran.w <- rep(NA, n.ran.w)
    NAs.prd.w <- rep(NA, n.prd.w)
    NAs.stat.w <- rep(NA, n.stat.w)
    NAs.w <- rep(NA, n.com.w + n.ran.w + n.prd.w)
    #
    NAs.all <- c(NAs, NAs.w, NAs.stat.w)
    #
    blanks.com.w <- rep("", n.com.w)
    blanks.ran.w <- rep("", n.ran.w)
    blanks.prd.w <- rep("", n.prd.w)
    blanks.stat.w <- rep("", n.stat.w)
    #
    blanks.w <- rep("", n.com.w + n.ran.w + n.prd.w + n.stat.w)
    #
    idx.com.w <- 1 + (seq_len(n.by) * n.com) - n.com
    idx.ran.w <- 1 + (seq_len(n.by) * n.ran) - n.ran
    idx.prd.w <- 1 + (seq_len(n.by) * n.prd) - n.prd
    #
    notfirst.com.w <- seq_len(n.com.w)[!(seq_len(n.com.w) %in% idx.com.w)]
    notfirst.ran.w <- seq_len(n.ran.w)[!(seq_len(n.ran.w) %in% idx.ran.w)]
    notfirst.prd.w <- seq_len(n.prd.w)[!(seq_len(n.prd.w) %in% idx.prd.w)]
    #
    # Do not consider limits of confidence or prediction intervals in
    # subgroups to format confidence limits
    #
    sel.NA <- !repl(common.subgroup.logical, n.com, n.by)
    lower.common.w[sel.NA] <- NA
    upper.common.w[sel.NA] <- NA
    #
    sel.NA <- !repl(random.subgroup.logical, n.ran, n.by)
    lower.random.w[sel.NA] <- NA
    upper.random.w[sel.NA] <- NA
    #
    sel.NA <- !repl(prediction.subgroup.logical, n.prd, n.by)
    lower.predict.w[sel.NA] <- NA
    upper.predict.w[sel.NA] <- NA
    #
    first.com.w <- n.com + n.ran + n.prd + idx.com.w
    first.ran.w <- n.com + n.ran + n.prd + n.com.w + idx.ran.w
    first.prd.w <- n.com + n.ran + n.prd + n.com.w + n.ran.w + idx.prd.w
    #
    all.com.w <-
      c(seq_len(n.com),
        n.com + n.ran + n.prd + seq_len(n.com.w))
    all.ran.w <-
      c(n.com + seq_len(n.ran),
        n.com + n.ran + n.prd + n.com.w + seq_len(n.ran.w))
    all.prd.w <-
      c(all.prd, first.prd.w,
        n.com + n.ran + n.prd + n.com.w + n.ran.w + seq_len(n.prd.w))
    #
    all.het.w <- n.com + n.ran + n.prd + n.com.w + n.ran.w + n.prd.w +
      seq_len(n.by)
    all.efc.w <- n.com + n.ran + n.prd + n.com.w + n.ran.w + n.prd.w +
      n.by + seq_len(n.by)
    all.efr.w <- n.com + n.ran + n.prd + n.com.w + n.ran.w + n.prd.w +
      2 * n.by + seq_len(n.by)
    #
    all.stat.w <- c(all.het.w, all.efc.w, all.efr.w)
    #
    emp.w <- c(all.prd.w, all.stat.w)
    
    emp.com.w <- c(1 + seq_len(n.com - 1),
                   n.com + n.ran + n.prd + notfirst.com.w)
    emp.ran.w <- c(n.com + 1 + seq_len(n.ran - 1),
                   n.com + n.ran + n.prd + n.com.w + notfirst.ran.w)
    emp.prd.w <- c(all.prd,
                   first.prd.w,
                   n.com + n.ran + n.prd + n.com.w + n.ran.w + notfirst.prd.w)
    #
    emp.w <- c(emp.com.w, emp.ran.w, emp.prd.w,
               all.het.w, all.efc.w, all.efr.w)
    #
    sel.w <- c(first.com.w, first.ran.w, first.prd.w, all.stat.w)
    #
    if (common) {
      if (!overall) {
        i <- 0
        for (bylev.i in subgroup.levels) {
          i <- i + 1
          sel.i <- subgroup == bylev.i
          x$w.common[sel.i] <- x$w.common[sel.i] / w.common.w[i]
        }
        w.common.w.p <- ifelse(is.na(w.common.w), NA, 100)
      }
      else {
        if (!all(is.na(w.common.w)) && sum(w.common.w) > 0)
          w.common.w.p <-
            round(100 * w.common.w / sum(w.common.w, na.rm = TRUE),
                  digits.weight)
        else
          w.common.w.p <- w.common.w
      }
    }
    else {
      TE.common.w <- lower.common.w <- upper.common.w <- NAs.by
      #
      w.common.w.p <- NAs.by
      #
      if (!missing.text.common.w) {
        text.common.w <-
          deprecated(text.common.w, missing.text.common.w, args,
                     "text.common.w",
                     warn.deprecated)
        missing.text.common.w <- FALSE
      }
      if (missing.text.common.w)
        text.common.w <- rep("Overall", n.com.w)
    }
    #
    if (random) {
      if (!overall) {
        i <- 0
        for (bylev.i in subgroup.levels) {
          i <- i + 1
          sel.i <- subgroup == bylev.i
          x$w.random[sel.i] <- x$w.random[sel.i] / w.random.w[i]
        }
        w.random.w.p <- ifelse(is.na(w.random.w), NA, 100)
      }
      else {
        if (!all(is.na(w.random.w)) && sum(w.random.w) > 0)
          w.random.w.p <-
            round(100 * w.random.w / sum(w.random.w, na.rm = TRUE),
                  digits.weight)
        else
          w.random.w.p <- w.random.w
      }
    }
    else {
      TE.random.w <- lower.random.w <- upper.random.w <- NAs.by
      #
      w.random.w.p <- NAs.ran.w
      text.random.w <- blanks.ran.w
    }
    #
    hetstat.w <- vector("list", n.by)
    #
    if (is.character(hetstat) || hetstat) {
      #
      hetstat.I2.w <-
        paste0(hetseparator,
               round(100 * I2.w, digits.I2), "%",
               if (print.I2.ci)
                 ifelse(!(is.na(lowI2.w) | is.na(uppI2.w)),
                        pasteCI(100 * lowI2.w, 100 * uppI2.w,
                                digits.I2, big.mark,
                                text.NA = lab.NA, unit = "%"),
                        ""))
      #
      hetstat.tau2.w <-
        paste0(hetseparator,
               ifelse(is.na(tau.w), "NA",
               ifelse(tau2.w == 0, "0",
                      formatPT(tau2.w, digits = digits.tau2,
                               big.mark = big.mark, lab.NA = "NA"))))
      #
      if (print.tau2.ci)
        hetstat.tau2.w <-
        paste0(hetstat.tau2.w,
               ifelse(!(is.na(lower.tau2.w) | is.na(upper.tau2.w)),
                      pasteCI(lower.tau2.w, upper.tau2.w, digits.tau2, big.mark,
                              text.NA = lab.NA),
                      ""))
      #
      hetstat.tau.w <-
        paste0(hetseparator,
               ifelse(is.na(tau.w), "NA",
               ifelse(tau.w == 0, "0",
                      formatPT(tau.w, digits = digits.tau,
                               big.mark = big.mark, lab.NA = "NA"))))
      #
      if (print.tau.ci)
        hetstat.tau.w <-
        paste0(hetstat.tau.w,
               ifelse(!(is.na(lower.tau.w) | is.na(upper.tau.w)),
                      pasteCI(lower.tau.w, upper.tau.w, digits.tau, big.mark,
                              text.NA = lab.NA),
                      ""))
      #
      hetstat.Q.w <-
        paste0(hetseparator, round(Q.w, digits.Q),
               if (bmj.revman5)
                 paste0(", df",  hetseparator, k.w.hetstat - 1))
      #
      hetstat.pval.Q.w <-
        paste0(formatPT(replaceNULL(pval.Q.w, pvalQ(Q.w, k.w - 1)),
                        lab = TRUE, labval = "",
                        digits = digits.pval.Q,
                        zero = zero.pval, JAMA = JAMA.pval,
                        scientific = scientific.pval,
                        lab.NA = "NA"))
      #
      hetstat.Rb.w <-
        paste0(hetseparator,
               round(100 * Rb.w, digits.I2), "%",
               if (print.Rb.ci)
                 ifelse(!(is.na(lowRb.w) | is.na(uppRb.w)),
                        pasteCI(100 * lowRb.w, 100 * uppRb.w,
                                digits.I2, big.mark,
                                text.NA = lab.NA, unit = "%"),
                        ""))
      #
      # Remove superfluous spaces
      #
      hetstat.pval.Q.w <- rmSpace(hetstat.pval.Q.w, end = TRUE)
      #
      while(any(grepl("  ", hetstat.I2.w)))
        hetstat.I2.w <- gsub("  ", " ", hetstat.I2.w)
      while(any(grepl("  ", hetstat.tau2.w)))
        hetstat.tau2.w <- gsub("  ", " ", hetstat.tau2.w)
      while(any(grepl("  ", hetstat.tau.w)))
        hetstat.tau.w <- gsub("  ", " ", hetstat.tau.w)
      while(any(grepl("  ", hetstat.Q.w)))
        hetstat.Q.w <- gsub("  ", " ", hetstat.Q.w)
      while(any(grepl("  ", hetstat.pval.Q.w)))
        hetstat.pval.Q.w <- gsub("  ", " ", hetstat.pval.Q.w)
      while(any(grepl("  ", hetstat.Rb.w)))
        hetstat.Rb.w <- gsub("  ", " ", hetstat.Rb.w)
      #
      for (i in seq_len(n.by)) {
        if (bmj) {
          if (print.tau)
            hetstat.w[[i]] <-
              substitute(
                paste(hl,
                      tau, ht, "; ",
                      chi^2, hq,
                      ", P", hp,
                      "; ",
                      I^2, hi),
                list(hl = hetlab,
                     ht = hetstat.tau.w[i],
                     hq = hetstat.Q.w[i],
                     hp = hetstat.pval.Q.w[i],
                     hi = hetstat.I2.w[i])
                )
          else if (print.tau2)
            hetstat.w[[i]] <-
              substitute(
                paste(hl,
                      tau^2, ht, "; ",
                      chi^2, hq,
                      ", P", hp,
                      "; ",
                      I^2, hi),
                list(hl = hetlab,
                     ht = hetstat.tau2.w[i],
                     hq = hetstat.Q.w[i],
                     hp = hetstat.pval.Q.w[i],
                     hi = hetstat.I2.w[i])
              )
          else
            hetstat.w[[i]] <-
              substitute(
                paste(hl,
                      chi^2, hq,
                      ", P", hp,
                      "; ",
                      I^2, hi),
                list(hl = hetlab,
                     hq = hetstat.Q.w[i],
                     hp = hetstat.pval.Q.w[i],
                     hi = hetstat.I2.w[i])
              )
        }
        else if (jama) {
          if (!missing.print.tau2 | !missing.print.tau) {
            if (print.tau)
              hetstat.w[[i]] <-
                substitute(
                  paste(hl,
                        chi[df]^2, hq,
                        " (",
                        italic(P), hp,
                        "), ",
                        italic(I)^2, hi,
                        tau, ht),
                  list(hl = hetlab,
                       df = k.w.hetstat[i] - 1,
                       hq = hetstat.Q.w[i],
                       hp = hetstat.pval.Q.w[i],
                       hi = hetstat.I2.w[i],
                       ht = hetstat.tau.w[i])
                  )
            else if (print.tau2)
              hetstat.w[[i]] <-
                substitute(
                  paste(hl,
                        chi[df]^2, hq,
                        " (",
                        italic(P), hp,
                        "), ",
                        italic(I)^2, hi,
                        tau^2, ht),
                  list(hl = hetlab,
                       df = k.w.hetstat[i] - 1,
                       hq = hetstat.Q.w[i],
                       hp = hetstat.pval.Q.w[i],
                       hi = hetstat.I2.w[i],
                       ht = hetstat.tau2.w[i])
                )
            else
              hetstat.w[[i]] <-
                substitute(
                  paste(hl,
                        chi[df]^2, hq,
                        " (",
                        italic(P), hp,
                        "), ",
                        italic(I)^2, hi),
                  list(hl = hetlab,
                       df = k.w.hetstat[i] - 1,
                       hq = hetstat.Q.w[i],
                       hp = hetstat.pval.Q.w[i],
                       hi = hetstat.I2.w[i])
                  )
          }
          else
            hetstat.w[[i]] <-
              substitute(
                paste(hl,
                      chi[df]^2, hq,
                      " (",
                      italic(P), hp,
                      "), ",
                      italic(I)^2, hi),
                list(hl = hetlab,
                     df = k.w.hetstat[i] - 1,
                     hq = hetstat.Q.w[i],
                     hp = hetstat.pval.Q.w[i],
                     hi = hetstat.I2.w[i])
              )
        }
        else if (revman5) {
          if (print.tau)
            hetstat.w[[i]] <-
              substitute(
                paste(hl,
                      "Tau", ht, "; ",
                      "Chi"^2, hq,
                      " (",
                      P, hp,
                      "); ",
                      I^2, hi),
                list(hl = hetlab,
                     ht = hetstat.tau.w[i],
                     hq = hetstat.Q.w[i],
                     hp = hetstat.pval.Q.w[i],
                     hi = hetstat.I2.w[i])
              )
          else
            hetstat.w[[i]] <-
              substitute(
                paste(hl,
                      "Tau"^2, ht, "; ",
                      "Chi"^2, hq,
                      " (",
                      P, hp,
                      "); ",
                      I^2, hi),
                list(hl = hetlab,
                     ht = hetstat.tau2.w[i],
                     hq = hetstat.Q.w[i],
                     hp = hetstat.pval.Q.w[i],
                     hi = hetstat.I2.w[i])
              )
        }
        else {
          #
          # One
          #
          if (print.I2 & !print.tau2.tau & !print.Q & !print.pval.Q & !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, italic(I)^2, hi),
                         list(hl = hetlab, hi = hetstat.I2.w[i]))
          else if (!print.I2 & print.tau2.tau & !print.Q & !print.pval.Q &
                   !print.Rb)
            if (print.tau2)
              hetstat.w[[i]] <-
                substitute(paste(hl, tau^2, ht),
                           list(hl = hetlab, ht = hetstat.tau2.w[i]))
            else
              hetstat.w[[i]] <-
                substitute(paste(hl, tau, ht),
                           list(hl = hetlab, ht = hetstat.tau.w[i]))              
          else if (!print.I2 & !print.tau2.tau & print.Q & !print.pval.Q &
                   !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, chi[df]^2, hq),
                         list(hl = hetlab, df = k.w.hetstat[i] - 1,
                              hq = hetstat.Q.w[i]))
          else if (!print.I2 & !print.tau2.tau & !print.Q & print.pval.Q &
                   !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, italic(p), hp),
                         list(hl = hetlab, hp = hetstat.pval.Q.w[i]))
          else if (!print.I2 & !print.tau2.tau & !print.Q & !print.pval.Q &
                   print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, italic(R)[italic(b)], hb),
                         list(hl = hetlab, hb = hetstat.Rb.w[i]))
          #
          # Two
          #
          else if (print.I2 & print.tau2.tau & !print.Q & !print.pval.Q &
                   !print.Rb)
            if (print.tau2)
              hetstat.w[[i]] <-
                substitute(paste(hl, italic(I)^2, hi,
                                 ", ",
                                 tau^2, ht),
                           list(hl = hetlab,
                                hi = hetstat.I2.w[i], ht = hetstat.tau2.w[i]))
            else
              hetstat.w[[i]] <-
                substitute(paste(hl, italic(I)^2, hi,
                                 ", ",
                                 tau, ht),
                           list(hl = hetlab,
                                hi = hetstat.I2.w[i], ht = hetstat.tau.w[i]))              
          else if (print.I2 & !print.tau2.tau & print.Q & !print.pval.Q &
                   !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, italic(I)^2, hi,
                               ", ",
                               chi[df]^2, hq),
                         list(hl = hetlab,
                              hi = hetstat.I2.w[i], hq = hetstat.Q.w[i]))
          else if (print.I2 & !print.tau2.tau & !print.Q & print.pval.Q &
                   !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, italic(I)^2, hi,
                               ", ",
                               italic(p), hp),
                         list(hl = hetlab,
                              hi = hetstat.I2.w[i], hp = hetstat.pval.Q.w[i]))
          else if (print.I2 & !print.tau2.tau & !print.Q & !print.pval.Q &
                   print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, italic(I)^2, hi,
                               ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab,
                              hi = hetstat.I2.w[i], hb = hetstat.Rb.w[i]))
          else if (!print.I2 & print.tau2.tau & print.Q & !print.pval.Q &
                   !print.Rb)
            if (print.tau2)
              hetstat.w[[i]] <-
                substitute(paste(hl, tau^2, ht,
                                 ", ",
                                 chi[df]^2, hq),
                           list(hl = hetlab, df = k.w.hetstat[i] - 1,
                                ht = hetstat.tau2.w[i], hq = hetstat.Q.w[i]))
            else
              hetstat.w[[i]] <-
                substitute(paste(hl, tau, ht,
                                 ", ",
                                 chi[df]^2, hq),
                           list(hl = hetlab, df = k.w.hetstat[i] - 1,
                                ht = hetstat.tau.w[i], hq = hetstat.Q.w[i]))              
          else if (!print.I2 & print.tau2.tau & !print.Q & print.pval.Q &
                   !print.Rb)
            if (print.tau2)
              hetstat.w[[i]] <-
                substitute(paste(hl, tau^2, ht,
                                 ", ",
                                 italic(p), hp),
                           list(hl = hetlab,
                                ht = hetstat.tau2.w[i],
                                hp = hetstat.pval.Q.w[i]))
            else
              hetstat.w[[i]] <-
                substitute(paste(hl, tau, ht,
                                 ", ",
                                 italic(p), hp),
                           list(hl = hetlab,
                                ht = hetstat.tau.w[i],
                                hp = hetstat.pval.Q.w[i]))             
          else if (!print.I2 & print.tau2.tau & !print.Q & !print.pval.Q &
                   print.Rb)
            if (print.tau2)
              hetstat.w[[i]] <-
                substitute(paste(hl, tau^2, ht,
                                 ", ",
                                 italic(R)[italic(b)], hb),
                           list(hl = hetlab,
                                ht = hetstat.tau2.w[i], hb = hetstat.Rb.w[i]))
            else
              hetstat.w[[i]] <-
                substitute(paste(hl, tau, ht,
                                 ", ",
                                 italic(R)[italic(b)], hb),
                           list(hl = hetlab,
                                ht = hetstat.tau.w[i], hb = hetstat.Rb.w[i]))
          else if (!print.I2 & !print.tau2.tau & print.Q & print.pval.Q &
                   !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, chi[df]^2, hq,
                               " (",
                               italic(p), hp, ")"),
                         list(hl = hetlab, df = k.w.hetstat[i] - 1,
                              hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i]))
          else if (!print.I2 & !print.tau2.tau & print.Q & !print.pval.Q &
                   print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, chi[df]^2, hq,
                               ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab, df = k.w.hetstat[i] - 1,
                              hq = hetstat.Q.w[i], hb = hetstat.Rb.w[i]))
          else if (!print.I2 & !print.tau2.tau & !print.Q & print.pval.Q &
                   print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, italic(p), hp,
                               ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab,
                              hp = hetstat.pval.Q.w[i], hb = hetstat.Rb.w[i]))
          #
          # Three
          #
          else if (print.I2 & print.tau2.tau & print.Q & !print.pval.Q &
                   !print.Rb)
            if (print.tau2)
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 italic(I)^2, hi, ", ",
                                 tau^2, ht, ", ",
                                 chi[df]^2, hq),
                           list(hl = hetlab, df = k.w.hetstat[i] - 1,
                                hi = hetstat.I2.w[i], ht = hetstat.tau2.w[i],
                                hq = hetstat.Q.w[i]))
            else
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 italic(I)^2, hi, ", ",
                                 tau, ht, ", ",
                                 chi[df]^2, hq),
                           list(hl = hetlab, df = k.w.hetstat[i] - 1,
                                hi = hetstat.I2.w[i], ht = hetstat.tau.w[i],
                                hq = hetstat.Q.w[i]))              
          else if (print.I2 & print.tau2.tau & !print.Q & print.pval.Q &
                   !print.Rb)
            if (print.tau2)
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 italic(I)^2, hi, ", ",
                                 tau^2, ht, ", ",
                                 italic(p), hp),
                           list(hl = hetlab,
                                hi = hetstat.I2.w[i], ht = hetstat.tau2.w[i],
                                hp = hetstat.pval.Q.w[i]))
            else   
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 italic(I)^2, hi, ", ",
                                 tau, ht, ", ",
                                 italic(p), hp),
                           list(hl = hetlab,
                                hi = hetstat.I2.w[i], ht = hetstat.tau.w[i],
                                hp = hetstat.pval.Q.w[i]))
          else if (print.I2 & !print.tau2.tau & print.Q & print.pval.Q &
                   !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               italic(I)^2, hi, ", ",
                               chi[df]^2, hq,
                               " (", italic(p), hp, ")"),
                         list(hl = hetlab, df = k.w.hetstat[i] - 1,
                              hi = hetstat.I2.w[i],
                              hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i]))
          else if (!print.I2 & print.tau2.tau & print.Q & print.pval.Q &
                   !print.Rb)
            if (print.tau2)
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 tau^2, ht, ", ",
                                 chi[df]^2, hq,
                                 " (", italic(p), hp, ")"),
                           list(hl = hetlab, df = k.w.hetstat[i] - 1,
                                ht = hetstat.tau2.w[i],
                                hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i]))
            else
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 tau, ht, ", ",
                                 chi[df]^2, hq,
                                 " (", italic(p), hp, ")"),
                           list(hl = hetlab, df = k.w.hetstat[i] - 1,
                                ht = hetstat.tau.w[i],
                                hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i]))
          else if (print.I2 & print.tau2.tau & !print.Q & !print.pval.Q &
                   print.Rb)
            if (print.tau2)
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 italic(I)^2, hi, ", ",
                                 tau^2, ht, ", ",
                                 italic(R)[italic(b)], hb),
                           list(hl = hetlab,
                                hi = hetstat.I2.w[i], ht = hetstat.tau2.w[i],
                                hb = hetstat.Rb.w[i]))
            else
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 italic(I)^2, hi, ", ",
                                 tau, ht, ", ",
                                 italic(R)[italic(b)], hb),
                           list(hl = hetlab,
                                hi = hetstat.I2.w[i], ht = hetstat.tau.w[i],
                                hb = hetstat.Rb.w[i]))
          else if (print.I2 & !print.tau2.tau & print.Q & !print.pval.Q &
                   print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               italic(I)^2, hi, ", ",
                               chi[df]^2, hq, ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab, df = k.w.hetstat[i] - 1,
                              hi = hetstat.I2.w[i],
                              hq = hetstat.Q.w[i],
                              hb = hetstat.Rb.w[i]))
          else if (!print.I2 & print.tau2.tau & print.Q & !print.pval.Q &
                   print.Rb)
            if (print.tau2)
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 tau^2, ht, ", ",
                                 chi[df]^2, hq, ", ",
                                 italic(R)[italic(b)], hb),
                           list(hl = hetlab, df = k.w.hetstat[i] - 1,
                                ht = hetstat.tau2.w[i],
                                hq = hetstat.Q.w[i],
                                hb = hetstat.Rb.w[i]))
            else
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 tau, ht, ", ",
                                 chi[df]^2, hq, ", ",
                                 italic(R)[italic(b)], hb),
                           list(hl = hetlab, df = k.w.hetstat[i] - 1,
                                ht = hetstat.tau.w[i],
                                hq = hetstat.Q.w[i],
                                hb = hetstat.Rb.w[i]))
          else if (print.I2 & !print.tau2.tau & !print.Q & print.pval.Q &
                   print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               italic(I)^2, hi, ", ",
                               italic(p), hp, ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab,
                              hi = hetstat.I2.w[i], hp = hetstat.pval.Q.w[i],
                              hb = hetstat.Rb.w[i]))
          else if (!print.I2 & print.tau2.tau & !print.Q & print.pval.Q &
                   print.Rb)
            if (print.tau2)
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 tau^2, ht, ", ",
                                 italic(p), hp, ", ",
                                 italic(R)[italic(b)], hb),
                           list(hl = hetlab,
                                ht = hetstat.tau2.w[i],
                                hp = hetstat.pval.Q.w[i],
                                hb = hetstat.Rb.w[i]))
            else              
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 tau, ht, ", ",
                                 italic(p), hp, ", ",
                                 italic(R)[italic(b)], hb),
                           list(hl = hetlab,
                                ht = hetstat.tau.w[i],
                                hp = hetstat.pval.Q.w[i],
                                hb = hetstat.Rb.w[i]))
          else if (!print.I2 & !print.tau2.tau & print.Q & print.pval.Q &
                   print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               chi[df]^2, hq,
                               " (", italic(p), hp, ")",
                               ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab, df = k.w.hetstat[i] - 1,
                              hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i],
                              hb = hetstat.Rb.w[i]))      
          #
          # Four
          #
          if (print.I2 & print.tau2.tau & print.Q & print.pval.Q & !print.Rb)
            if (print.tau2)
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 italic(I)^2, hi, ", ",
                                 tau^2, ht, ", ",
                                 chi[df]^2, hq,
                                 " (", italic(p), hp, ")"),
                           list(hl = hetlab, df = k.w.hetstat[i] - 1,
                                hi = hetstat.I2.w[i], ht = hetstat.tau2.w[i],
                                hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i]))
            else
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 italic(I)^2, hi, ", ",
                                 tau, ht, ", ",
                                 chi[df]^2, hq,
                                 " (", italic(p), hp, ")"),
                           list(hl = hetlab, df = k.w.hetstat[i] - 1,
                                hi = hetstat.I2.w[i], ht = hetstat.tau.w[i],
                                hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i]))
          else if (print.I2 & print.tau2.tau & print.Q & !print.pval.Q &
                   print.Rb)
            if (print.tau2)
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 italic(I)^2, hi, ", ",
                                 tau^2, ht, ", ",
                                 chi[df]^2, hq, ", ",
                                 italic(R)[italic(b)], hb),
                           list(hl = hetlab, df = k.w.hetstat[i] - 1,
                                hi = hetstat.I2.w[i], ht = hetstat.tau2.w[i],
                                hq = hetstat.Q.w[i], hb = hetstat.Rb.w[i]))
            else
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 italic(I)^2, hi, ", ",
                                 tau, ht, ", ",
                                 chi[df]^2, hq, ", ",
                                 italic(R)[italic(b)], hb),
                           list(hl = hetlab, df = k.w.hetstat[i] - 1,
                                hi = hetstat.I2.w[i], ht = hetstat.tau.w[i],
                                hq = hetstat.Q.w[i], hb = hetstat.Rb.w[i]))
          else if (print.I2 & print.tau2.tau & !print.Q & print.pval.Q &
                   print.Rb)
            if (print.tau2)
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 italic(I)^2, hi, ", ",
                                 tau^2, ht, ", ",
                                 italic(p), hp, ", ",
                                 italic(R)[italic(b)], hb),
                           list(hl = hetlab,
                                hi = hetstat.I2.w[i], ht = hetstat.tau2.w[i],
                                hp = hetstat.pval.Q.w[i], hb = hetstat.Rb.w[i]))
            else
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 italic(I)^2, hi, ", ",
                                 tau, ht, ", ",
                                 italic(p), hp, ", ",
                                 italic(R)[italic(b)], hb),
                           list(hl = hetlab,
                                hi = hetstat.I2.w[i], ht = hetstat.tau.w[i],
                                hp = hetstat.pval.Q.w[i], hb = hetstat.Rb.w[i]))
          else if (print.I2 & !print.tau2.tau & print.Q & print.pval.Q &
                   print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               italic(I)^2, hi, ", ",
                               chi[df]^2, hq,
                               " (", italic(p), hp, ")",
                               ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab,
                              hi = hetstat.I2.w[i], df = k.w.hetstat[i] - 1,
                              hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i],
                              hb = hetstat.Rb.w[i]))
          else if (!print.I2 & print.tau2.tau & print.Q & print.pval.Q &
                   print.Rb)
            if (print.tau2)
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 tau^2, ht, ", ",
                                 chi[df]^2, hq,
                                 " (", italic(p), hp, ")",
                                 ", ",
                                 italic(R)[italic(b)], hb),
                           list(hl = hetlab, df = k.w.hetstat[i] - 1,
                                ht = hetstat.tau2.w[i],
                                hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i],
                                hb = hetstat.Rb.w[i]))
            else
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 tau, ht, ", ",
                                 chi[df]^2, hq,
                                 " (", italic(p), hp, ")",
                                 ", ",
                                 italic(R)[italic(b)], hb),
                           list(hl = hetlab, df = k.w.hetstat[i] - 1,
                                ht = hetstat.tau.w[i],
                                hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i],
                                hb = hetstat.Rb.w[i]))
          #
          # Five
          #
          else if (print.I2 & print.tau2.tau & print.Q & print.pval.Q &
                   print.Rb)
            if (print.tau2)
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 italic(I)^2, hi, ", ",
                                 tau^2, ht, ", ",
                                 chi[df]^2, hq,
                                 " (", italic(p), hp, ")",
                                 ", ",
                                 italic(R)[italic(b)], hb),
                           list(hl = hetlab, df = k.w.hetstat[i] - 1,
                                hi = hetstat.I2.w[i], ht = hetstat.tau2.w[i],
                                hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i],
                                hb = hetstat.Rb.w[i]))
            else
              hetstat.w[[i]] <-
                substitute(paste(hl,
                                 italic(I)^2, hi, ", ",
                                 tau, ht, ", ",
                                 chi[df]^2, hq,
                                 " (", italic(p), hp, ")",
                                 ", ",
                                 italic(R)[italic(b)], hb),
                           list(hl = hetlab, df = k.w.hetstat[i] - 1,
                                hi = hetstat.I2.w[i], ht = hetstat.tau.w[i],
                                hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i],
                                hb = hetstat.Rb.w[i]))
        }
        if (hetstat.pooled != "study" &
            !is.logical(text.subgroup.nohet) & k.w.hetstat[i] < 2)
          hetstat.w[[i]] <- paste0(hetlab, text.subgroup.nohet)
      }
    }
    #
    TE.w <- c(TE.common.w, TE.random.w, NAs.prd.w, NAs.stat.w)
    lowTE.w <-
      c(lower.common.w, lower.random.w, lower.predict.w, NAs.stat.w)
    uppTE.w <-
      c(upper.common.w, upper.random.w, upper.predict.w, NAs.stat.w)
    #
    n.harmonic.mean.w <-
      c(rep(n.harmonic.mean.w, n.com + n.ran + n.prd), NAs.stat.w)
    t.harmonic.mean.w <-
      c(rep(t.harmonic.mean.w, n.com + n.ran + n.prd), NAs.stat.w)
    #
    # Label of test for effect in subgroups
    #
    pvals.effect.w <-
      formatPT(c(pval.common.w, pval.random.w),
               lab = TRUE, labval = "",
               digits = digits.pval,
               zero = zero.pval, JAMA = JAMA.pval,
               scientific = scientific.pval,
               lab.NA = lab.NA)
    statistics.effect.w <-
      formatN(c(statistic.common.w, statistic.random.w),
              digits.stat, lab.NA, big.mark = big.mark)
    #
    # Remove superfluous spaces
    #
    pvals.effect.w <- rmSpace(pvals.effect.w, end = TRUE)
    #
    while(any(grepl("  ", pvals.effect.w)))
      pvals.effect.w <- gsub("  ", " ", pvals.effect.w)
    while(any(grepl("  ", statistics.effect.w)))
      statistics.effect.w <- gsub("  ", " ", statistics.effect.w)
    #
    if (any(test.effect.subgroup.common.logical)) {
      text.effect.subgroup.common <- vector("list", n.by)
      for (i in seq_len(n.by)) {
        if (print.stat) {
          if (bmj)
            text.effect.subgroup.common[[i]] <-
              substitute(paste(tl,
                               Z, hetseparator, tt,
                               ", P", tp),
                         list(tl = label.test.effect.subgroup.common,
                              hetseparator = hetseparator,
                              tt = statistics.effect.w[i],
                              tp = rmSpace(pvals.effect.w[i], end = TRUE)))
          else if (revman5)
            text.effect.subgroup.common[[i]] <-
              substitute(paste(tl,
                               Z, hetseparator, tt,
                               " (P", tp, ")"),
                         list(tl = label.test.effect.subgroup.common,
                              hetseparator = hetseparator,
                              tt = statistics.effect.w[i],
                              tp = rmSpace(pvals.effect.w[i], end = TRUE)))
          else if (jama)
            text.effect.subgroup.common[[i]] <-
              substitute(paste(tl,
                               italic(z), hetseparator, tt,
                               " (", italic(P), tp, ")"),
                         list(tl = label.test.effect.subgroup.common,
                              hetseparator = hetseparator,
                              tt = statistics.effect.w[i],
                              tp = rmSpace(pvals.effect.w[i], end = TRUE)))
          else
            text.effect.subgroup.common[[i]] <-
              substitute(paste(tl,
                               italic(z), hetseparator, tt,
                               " (", italic(p), tp, ")"),
                         list(tl = label.test.effect.subgroup.common,
                              hetseparator = hetseparator,
                              tt = statistics.effect.w[i],
                              tp = rmSpace(pvals.effect.w[i], end = TRUE)))
        }
        else {
          if (bmj)
            text.effect.subgroup.common[[i]] <-
              substitute(paste(tl,
                               " P", tp),
                         list(tl = label.test.effect.subgroup.common,
                              tp = rmSpace(pvals.effect.w[i], end = TRUE)))
          else if (revman5)
            text.effect.subgroup.common[[i]] <-
              substitute(paste(tl,
                               " P", tp),
                         list(tl = label.test.effect.subgroup.common,
                              tp = rmSpace(pvals.effect.w[i], end = TRUE)))
          else if (jama)
            text.effect.subgroup.common[[i]] <-
              substitute(paste(tl,
                               " ", italic(P), tp),
                         list(tl = label.test.effect.subgroup.common,
                              tp = rmSpace(pvals.effect.w[i], end = TRUE)))
          else
            text.effect.subgroup.common[[i]] <-
              substitute(paste(tl,
                               " ", italic(p), tp),
                         list(tl = label.test.effect.subgroup.common,
                              tp = rmSpace(pvals.effect.w[i], end = TRUE)))
        }
      }
    }
    else {
      text.effect.subgroup.common <- vector("list", n.by)
      for (i in seq_len(n.by))
        text.effect.subgroup.common[[i]] <- ""
    }
    #
    if (any(test.effect.subgroup.random.logical)) {
      text.effect.subgroup.random <- vector("list", n.by)
      for (i in seq_len(n.by)) {
        if (print.stat) {
          if (!hakn.kero) {
            if (bmj)
              text.effect.subgroup.random[[i]] <-
                substitute(paste(tl,
                                 Z, hetseparator, tt,
                                 ", P", tp),
                           list(tl = label.test.effect.subgroup.random,
                                hetseparator = hetseparator,
                                tt = statistics.effect.w[n.by + i],
                                tp = rmSpace(pvals.effect.w[n.by + i],
                                             end = TRUE)))
            else if (revman5)
              text.effect.subgroup.random[[i]] <-
                substitute(paste(tl,
                                 Z, hetseparator, tt,
                                 " (P", tp, ")"),
                           list(tl = label.test.effect.subgroup.random,
                                hetseparator = hetseparator,
                                tt = statistics.effect.w[n.by + i],
                                tp = rmSpace(pvals.effect.w[n.by + i],
                                             end = TRUE)))
            else if (jama)
              text.effect.subgroup.random[[i]] <-
                substitute(paste(tl,
                                 italic(z), hetseparator, tt,
                                 " (", italic(P), tp, ")"),
                           list(tl = label.test.effect.subgroup.random,
                                hetseparator = hetseparator,
                                tt = statistics.effect.w[n.by + i],
                                tp = rmSpace(pvals.effect.w[n.by + i],
                                             end = TRUE)))
            else
              text.effect.subgroup.random[[i]] <-
                substitute(paste(tl,
                                 italic(z), hetseparator, tt,
                                 " (", italic(p), tp, ")"),
                           list(tl = label.test.effect.subgroup.random,
                                hetseparator = hetseparator,
                                tt = statistics.effect.w[n.by + i],
                                tp = rmSpace(pvals.effect.w[n.by + i],
                                             end = TRUE)))
          }
          else {
            if (bmj)
              text.effect.subgroup.random[[i]] <-
                substitute(paste(tl,
                                 t[df], hetseparator, tt,
                                 ", P", tp),
                           list(tl = label.test.effect.subgroup.random,
                                hetseparator = hetseparator,
                                tt = statistics.effect.w[n.by + i],
                                tp = rmSpace(pvals.effect.w[n.by + i],
                                             end = TRUE),
                                df = df.random.w[i]))
            else if (revman5)
              text.effect.subgroup.random[[i]] <-
                substitute(paste(tl,
                                 t[df], hetseparator, tt,
                                 " (P", tp, ")"),
                           list(tl = label.test.effect.subgroup.random,
                                hetseparator = hetseparator,
                                tt = statistics.effect.w[n.by + i],
                                tp = rmSpace(pvals.effect.w[n.by + i],
                                             end = TRUE),
                                df = df.random.w[i]))
            else if (jama)
              text.effect.subgroup.random[[i]] <-
                substitute(paste(tl,
                                 italic(t)[df], hetseparator, tt,
                                 " (", italic(P), tp, ")"),
                           list(tl = label.test.effect.subgroup.random,
                                hetseparator = hetseparator,
                                tt = statistics.effect.w[n.by + i],
                                tp = rmSpace(pvals.effect.w[n.by + i],
                                             end = TRUE),
                                df = df.random.w[i]))
            else
              text.effect.subgroup.random[[i]] <-
                substitute(paste(tl,
                                 italic(t)[df], hetseparator, tt,
                                 " (", italic(p), tp, ")"),
                           list(tl = label.test.effect.subgroup.random,
                                hetseparator = hetseparator,
                                tt = statistics.effect.w[n.by + i],
                                tp = rmSpace(pvals.effect.w[n.by + i],
                                             end = TRUE),
                                df = df.random.w[i]))
          }
        }
        else {
          if (bmj)
            text.effect.subgroup.random[[i]] <-
              substitute(paste(tl,
                               " P", tp),
                         list(tl = label.test.effect.subgroup.random,
                              tp = rmSpace(pvals.effect.w[n.by + i],
                                           end = TRUE)))
          else if (revman5)
            text.effect.subgroup.random[[i]] <-
              substitute(paste(tl,
                               " P", tp),
                         list(tl = label.test.effect.subgroup.random,
                              tp = rmSpace(pvals.effect.w[n.by + i],
                                           end = TRUE)))
          else if (jama)
            text.effect.subgroup.random[[i]] <-
              substitute(paste(tl,
                               " ", italic(P)),
                         list(tl = label.test.effect.subgroup.random,
                              tp = rmSpace(pvals.effect.w[n.by + i],
                                           end = TRUE)))
          else
            text.effect.subgroup.random[[i]] <-
              substitute(paste(tl,
                               " ", italic(p), tp),
                         list(tl = label.test.effect.subgroup.random,
                              tp = rmSpace(pvals.effect.w[n.by + i],
                                           end = TRUE)))
        }
      }
    }
    else {
      text.effect.subgroup.random <- vector("list", n.by)
      for (i in seq_len(n.by))
        text.effect.subgroup.random[[i]] <- ""
    }
  }
  #
  x$pscale <- pscale
  x$irscale <- irscale
  #
  if (!is.null(x$sd.n_of_1))
    x$sd.n_of_1 <-
    formatPT(x$sd.n_of_1, digits = digits.sd, big.mark = big.mark)
  #
  text.details <- ""
  #
  if (details) {
    if (is.null(x$.text.details.methods)) {
      if (K.all == 1) {
        text.details <-
          catmeth(x,
                  common, random, prediction, overall, overall.hetstat,
                  #
                  func.transf = x$func.transf,
                  backtransf = backtransf, func.backtransf = fbt,
                  #
                  big.mark = big.mark, digits = digits,
                  digits.tau = digits.tau,
                  text.tau = gs("text.tau"), text.tau2 = gs("text.tau2"),
                  #
                  print.tau2 = FALSE,
                  #
                  forest = TRUE)
      }
      else {
        text.details <-
          catmeth(x,
                  common, random, prediction, overall, overall.hetstat,
                  #
                  func.transf = x$func.transf,
                  backtransf = backtransf, func.backtransf = fbt,
                  #
                  big.mark = big.mark, digits = digits,
                  digits.tau = digits.tau,
                  text.tau = gs("text.tau"), text.tau2 = gs("text.tau2"),
                  #
                  print.tau2 = print.tau2, print.tau2.ci = print.tau2.ci,
                  print.tau = print.tau, print.tau.ci = print.tau.ci,
                  #
                  print.I2 = print.I2 &
                    (overall.hetstat | (by && any(subgroup.hetstat.logical))),
                  #
                  forest = TRUE)
      }
    }
    else {
      text.details <- unlist(strsplit(x$.text.details.methods, "\n"))
      text.details <- text.details[text.details != ""]
    }
    #
    text.details <- unlist(strsplit(text.details, "\n"))
    text.details <- text.details[text.details != ""]
    #
    td <- vector("list", length(text.details))
    #
    for (i in seq_along(text.details))
      td[[i]] <- text.details[i]
    #
    is.tau <- grepl(" for tau", td, fixed = TRUE)
    is.tau.c <- grepl("assuming common", td, fixed = TRUE)
    is.tau <- ifelse(is.tau & is.tau.c, FALSE, is.tau)
    #
    any.tau <- sum(is.tau)
    any.tau.c <- sum(is.tau.c)
    #
    if (any.tau | any.tau.c) {
      id.tau <- seq_along(is.tau)[is.tau]
      id.tau.c <- seq_along(is.tau.c)[is.tau.c]
      #
      if (print.tau) {
        if (any.tau.c)
          for (i in id.tau.c)
            td[[i]] <-
              gsub(" for tau (assuming common tau in subgroups)", " for ",
                   td[[i]], fixed = TRUE)
        #
        if (any.tau)
          for (i in id.tau)
            td[[i]] <-
              gsub(" for tau", " for ", td[[i]], fixed = TRUE)
      }
      else {
        if (any.tau.c)
          for (i in id.tau.c)
            td[[i]] <-
              gsub(" for tau^2 (assuming common tau^2 in subgroups)", " for ",
                   td[[i]], fixed = TRUE)
        #
        if (any.tau)
          for (i in id.tau)
            td[[i]] <-
              gsub(" for tau^2", " for ", td[[i]], fixed = TRUE)
      }
      #
      if (print.tau) {
        if (revman5) {
          if (any.tau.c)
            for (i in id.tau.c)
              td[[i]] <-
                substitute(paste(txt1, txt2),
                           list(txt1 = td[[i]],
                                txt2 =
                                  "Tau (assuming common Tau in subgroups)"))
          #
          if (any.tau)
            for (i in id.tau)
              td[[i]] <-
                substitute(paste(txt1, txt2),
                           list(txt1 = td[[i]], txt2 = "Tau"))
        }
        else {
          if (any.tau.c)
            for (i in id.tau.c)
              td[[i]] <-
                substitute(paste(txt1, tau, txt2, tau, txt3),
                           list(txt1 = td[[i]],
                                txt2 = " (assuming common ",
                                txt3 = " in subgroups)"))
          #
          if (any.tau)
            for (i in id.tau)
              td[[i]] <-
                substitute(paste(txt, tau), list(txt = td[[i]]))
        }
      }
      else {
        if (revman5) {
          if (any.tau.c)
            for (i in id.tau.c)
              td[[i]] <-
                substitute(paste(txt1, txt2^2, txt3, txt4^2, txt5),
                           list(txt1 = td[[i]],
                                txt2 = "Tau",
                                txt3 = " (assuming common ",
                                txt4 = "Tau",
                                txt5 = " in subgroups)"))
          #
          if (any.tau)
            for (i in id.tau)
              td[[i]] <-
                substitute(paste(txt1, txt2^2),
                           list(txt1 = td[[i]], txt2 = "Tau"))
        }
        else {
          if (any.tau.c)
            for (i in id.tau.c)
              td[[i]] <-
                substitute(paste(txt1, tau^2, txt2, tau^2, txt3),
                           list(txt1 = td[[i]],
                                txt2 = " (assuming common ",
                                txt3 = " in subgroups)"))
          #
          if (any.tau)
            for (i in id.tau)
              td[[i]] <-
                substitute(paste(txt, tau^2), list(txt = td[[i]]))
        }
      }
    }
    #
    tau.ci.line <-
      grepl(" of tau^2 and tau", text.details, fixed = TRUE) |
      grepl(" of tau", text.details, fixed = TRUE) |
      grepl(" of tau^2", text.details, fixed = TRUE)
    #
    if (any(tau.ci.line)) {
      id <- seq_along(tau.ci.line)[tau.ci.line]
      #
      for (i in id) {
        td[[i]] <- gsub(" of tau^2 and tau", " of ", td[[i]], fixed = TRUE)
        td[[i]] <- gsub(" of tau", " of ", td[[i]], fixed = TRUE)
        td[[i]] <- gsub(" of tau^2", " of ", td[[i]], fixed = TRUE)
      }
      #
      if (print.tau.ci) {
        if (revman5)
          for (i in id)
            td[[i]] <-
              substitute(paste(txt1, txt2),
                         list(txt1 = td[[i]], txt2 = "Tau"))
        else
          for (i in id)
            td[[i]] <- substitute(paste(txt1, tau), list(txt1 = td[[i]]))
      }
      else {
        if (revman5)
          for (i in id)
            td[[i]] <-
              substitute(paste(txt1, txt2^2),
                         list(txt1 = td[[i]], txt2 = "Tau"))
        else
          for (i in id)
            td[[i]] <- substitute(paste(txt1, tau^2), list(txt1 = td[[i]]))
      }
    }
    #
    method.I2.line <- grepl("Calculation of I^2", text.details, fixed = TRUE)
    #
    if (any(method.I2.line)) {
      id <- seq_along(method.I2.line)[method.I2.line]
      #
      if (revman5) {
        for (i in id) {
          with.Q <-
            grepl("Calculation of I^2 based on Q", td[[i]], fixed = TRUE)
          #
          td[[i]] <- gsub("I^2 based on Q", "", td[[i]], fixed = TRUE)
          td[[i]] <- gsub("I^2 based on tau^2", "", td[[i]], fixed = TRUE)
          #
          if (with.Q)
            td[[i]] <-
            substitute(paste(txt1, txt2^2, txt3^2),
                       list(txt1 = td[[i]], txt2 = "I",
                            txt3 = " based on Chi"))
          else
            td[[i]] <-
            substitute(paste(txt1, txt2^2, txt3, txt4^2),
                       list(txt1 = td[[i]], txt2 = "I",
                            txt3 = " based on ", txt4 = "Tau"))
        }
      }
      else {
        for (i in id) {
          with.Q <-
            grepl("Calculation of I^2 based on Q", td[[i]], fixed = TRUE)
          #
          td[[i]] <- gsub("I^2 based on Q", "", td[[i]], fixed = TRUE)
          td[[i]] <- gsub("I^2 based on tau^2", "", td[[i]], fixed = TRUE)
          #
          if (with.Q)
            td[[i]] <-
            substitute(paste(txt1, italic(I)^2, txt2, Q),
                       list(txt1 = td[[i]], txt2 = " based on "))
          else
            td[[i]] <-
            substitute(paste(txt1, italic(I)^2, txt3, tau^2),
                       list(txt1 = td[[i]], txt2 = "I",
                            txt3 = " based on "))
        }
      }
    }
    #
    if (length(td) == 0)
      text.details <- ""
    else
      text.details <- td
  }
  
  
  #
  #
  # (8) Backtransform data
  #
  #
  TE.orig <- TE
  #
  if (backtransf) {
    #
    # Freeman-Tukey Arcsin transformation
    #
    if (metabind) {
      if (sm == "IRFT") {
        npft <- x$t.harmonic.mean
        npft.ma <- x$t.harmonic.mean.ma
      }
      else {
        npft <- x$n.harmonic.mean
        npft.ma <- x$n.harmonic.mean.ma
      }
    }
    else {
      if (sm == "IRFT") {
        npft <- x$time
        npft.ma <- 1 / mean(1 / x$time)
      }
      else {
        npft <- x$n
        npft.ma <- 1 / mean(1 / x$n)
      }
    }
    #
    # Individual study results
    #
    if (metaprop)
      TE <- x$event.e / x$n.e
    else if (metarate)
      TE <- x$event.e / x$time.e
    else if (!log.xaxis)
      TE <- backtransf(TE, sm, npft, npft, fbt, abt)
    #
    if (!log.xaxis) {
      lowTE <- backtransf(lowTE, sm, npft, npft, fbt, abt)
      uppTE <- backtransf(uppTE, sm, npft, npft, fbt, abt)
    }
    #
    # Results of meta-analysis
    #
    if (!log.xaxis) {
      TE.common    <- backtransf(TE.common, sm, npft.ma, npft.ma, fbt, abt)
      lowTE.common <- backtransf(lowTE.common, sm, npft.ma, npft.ma, fbt, abt)
      uppTE.common <- backtransf(uppTE.common, sm, npft.ma, npft.ma, fbt, abt)
      #
      TE.random <- backtransf(TE.random, sm, npft.ma, npft.ma, fbt, abt)
      lowTE.random <- backtransf(lowTE.random, sm, npft.ma, npft.ma, fbt, abt)
      uppTE.random <- backtransf(uppTE.random, sm, npft.ma, npft.ma, fbt, abt)
      #
      lowTE.predict <-
        backtransf(lowTE.predict, sm, npft.ma, npft.ma, fbt, abt)
      uppTE.predict <-
        backtransf(uppTE.predict, sm, npft.ma, npft.ma, fbt, abt)
      #
      if (by) {
        if (sm == "IRFT")
          npft.w <- t.harmonic.mean.w
        else
          npft.w <- n.harmonic.mean.w
        #
        TE.w    <- backtransf(TE.w, sm, npft.w, npft.w, fbt, abt)
        lowTE.w <- backtransf(lowTE.w, sm, npft.w, npft.w, fbt, abt)
        uppTE.w <- backtransf(uppTE.w, sm, npft.w, npft.w, fbt, abt)
      }
    }
    #
    # Apply argument 'pscale' or 'irscale'
    #
    TE <- scale * TE
    lowTE <- scale * lowTE
    uppTE <- scale * uppTE
    #
    TE.common    <- scale * TE.common
    lowTE.common <- scale * lowTE.common
    uppTE.common <- scale * uppTE.common
    #
    TE.random    <- scale * TE.random
    lowTE.random <- scale * lowTE.random
    uppTE.random <- scale * uppTE.random
    #
    lowTE.predict <- scale * lowTE.predict
    uppTE.predict <- scale * uppTE.predict
    #
    if (by) {
      TE.w    <- scale * TE.w
      lowTE.w <- scale * lowTE.w
      uppTE.w <- scale * uppTE.w
    }
    #
    # Switch lower and upper limit for VE if results have been
    # backtransformed
    #
    if (sm == "VE") {
      tmp.l <- lowTE
      lowTE <- uppTE
      uppTE <- tmp.l
      #
      tmp.l <- lowTE.common
      lowTE.common <- uppTE.common
      uppTE.common <- tmp.l
      #
      tmp.l <- lowTE.random
      lowTE.random <- uppTE.random
      uppTE.random <- tmp.l
      #
      tmp.l <- lowTE.predict
      lowTE.predict <- uppTE.predict
      uppTE.predict <- tmp.l
      #
      if (by) {
        tmp.l <- lowTE.w
        lowTE.w <- uppTE.w
        uppTE.w <- tmp.l
      }
    }
  }
  #
  # Set study confidence intervals to NA if they should not be printed
  #
  if (!study.results) {
    lowTE[!is.na(lowTE)] <- NA
    uppTE[!is.na(uppTE)] <- NA
  }
  #
  # Exclude study results from forest plot
  #
  TE.exclude <- TE
  lowTE.exclude <- lowTE
  uppTE.exclude <- uppTE
  #
  if (avail.exclude) {
    TE.exclude[x$exclude] <- NA
    lowTE.exclude[x$exclude] <- NA
    uppTE.exclude[x$exclude] <- NA
  }
  #
  if (!common) {
    TE.common    <- NAs.com
    lowTE.common <- NAs.com
    uppTE.common <- NAs.com
    #
    if (by) {
      TE.common.w    <- NAs.com.w
      lowTE.common.w <- NAs.com.w
      uppTE.common.w <- NAs.com.w
    }
  }
  #
  method.random.ci <- x$method.random.ci
  #
  if (!random) {
    TE.random    <- NAs.ran
    lowTE.random <- NAs.ran
    uppTE.random <- NAs.ran
    #
    if (by) {
      TE.random.w    <- NAs.ran.w
      lowTE.random.w <- NAs.ran.w
      uppTE.random.w <- NAs.ran.w
    }
  }
  #
  if (!prediction) {
    lowTE.predict <- NAs.prd
    uppTE.predict <- NAs.prd
  }
  #
  if (by & !overall)
    w.common.p <- round(100 * x$w.common, digits.weight)
  else {
    if (!all(is.na(x$w.common)) && sum(x$w.common) > 0) {
      if (is.matrix(x$w.common))
        w.common.p <-
          round(apply(x$w.common, 2, calcPercent), digits.weight)
      else
        w.common.p <-
          round(calcPercent(x$w.common), digits.weight)
    }
    else
      w.common.p <- x$w.common
  }
  #
  if (by & !overall)
    w.random.p <- round(100 * x$w.random, digits.weight)
  else {
    if (!all(is.na(x$w.random)) && sum(x$w.random) > 0) {
      if (is.matrix(x$w.random))
        w.random.p <-
          round(apply(x$w.random, 2, calcPercent), digits.weight)
      else
        w.random.p <-
          round(calcPercent(x$w.random), digits.weight)
    }
    else
      w.random.p <- x$w.random
  }
  #
  if (is.matrix(w.common.p))
    w.common.p <- as.vector(w.common.p[, 1])
  if (is.matrix(w.random.p))
    w.random.p <- as.vector(w.random.p[, 1])
  
  
  #
  #
  # (9) Determine column labels
  #
  #
  labs <- list()
  #
  if (lsel) {
    if (missing.leftlabs || length(leftcols) != length(leftlabs)) {
      for (i in seq_along(leftcols)) {
        j <- match(leftcols[i], colnames)
        if (!is.na(j))
          labs[[paste0("lab.", leftcols[i])]] <- labnames[j]
      }
    }
    else if (length(leftcols) == length(leftlabs)) {
      for (i in seq_along(leftcols)) {
        j <- match(leftcols[i], colnames)
        if (!is.na(leftlabs[i]))
          labs[[paste0("lab.", leftcols[i])]] <- leftlabs[i]
        else
          if (!is.na(j))
            labs[[paste0("lab.", leftcols[i])]] <- labnames[j]
      }
    }
  }
  #
  if (missing.rightlabs || length(rightcols) != length(rightlabs)) {
    for (i in seq_along(rightcols)) {
      j <- match(rightcols[i], colnames)
      if (!is.na(j))
        labs[[paste0("lab.", rightcols[i])]] <- labnames[j]
    }
  }
  else if (length(rightcols) == length(rightlabs)) {
    for (i in seq_along(rightcols)) {
      j <- match(rightcols[i], colnames)
      if (!is.na(rightlabs[i]))
        labs[[paste0("lab.", rightcols[i])]] <- rightlabs[i]
      else
        if (!is.na(j))
          labs[[paste0("lab.", rightcols[i])]] <- labnames[j]
    }
  }
  #
  if (!slab)
    labs[["lab.studlab"]] <- ""
  #
  # Check for "%" in weight labels
  #
  w.common.percent <- TRUE
  if (!is.null(labs[["lab.w.common"]])) {
    if (grepl("%", labs[["lab.w.common"]]))
      w.common.percent <- FALSE
  }
  #
  w.random.percent <- TRUE
  if (!is.null(labs[["lab.w.random"]])) {
    if (grepl("%", labs[["lab.w.random"]]))
      w.random.percent <- FALSE
  }
  # "studlab", "TE", "seTE",
  # "cluster", "cycles",
  # "n.e", "n.c",
  # "event.e", "event.c",
  # "mean.e", "mean.c",
  # "sd.e", "sd.c",
  # "cor",
  # "time.e", "time.c"
  #
  # Check for "\n" in label of column 'studlab'
  #
  clines <- twolines(labs[["lab.studlab"]], "studlab")
  #
  if (clines$newline) {
    newline.studlab <- TRUE
    labs[["lab.studlab"]] <- clines$bottom
    add.studlab <- clines$top
    longer.studlab <- clines$longer
  }
  else {
    newline.studlab <- FALSE
    longer.studlab <- labs[["lab.studlab"]]
  }
  #
  # Check for "\n" in label of column 'TE'
  #
  clines <- twolines(labs[["lab.TE"]], "TE")
  #
  if (clines$newline) {
    newline.TE <- TRUE
    labs[["lab.TE"]] <- clines$bottom
    add.TE <- clines$top
    longer.TE <- clines$longer
  }
  else {
    newline.TE <- FALSE
    longer.TE <- labs[["lab.TE"]]
  }
  #
  # Check for "\n" in label of column 'seTE'
  #
  clines <- twolines(labs[["lab.seTE"]], "seTE")
  #
  if (clines$newline) {
    newline.seTE <- TRUE
    labs[["lab.seTE"]] <- clines$bottom
    add.seTE <- clines$top
    longer.seTE <- clines$longer
  }
  else {
    newline.seTE <- FALSE
    longer.seTE <- labs[["lab.seTE"]]
  }
  #
  # Check for "\n" in label of column 'cluster'
  #
  clines <- twolines(labs[["lab.cluster"]], "cluster")
  #
  if (clines$newline) {
    newline.cluster <- TRUE
    labs[["lab.cluster"]] <- clines$bottom
    add.cluster <- clines$top
    longer.cluster <- clines$longer
  }
  else {
    newline.cluster <- FALSE
    longer.cluster <- labs[["lab.cluster"]]
  }
  #
  # Check for "\n" in label of column 'cycles'
  #
  clines <- twolines(labs[["lab.cycles"]], "cycles")
  #
  if (clines$newline) {
    newline.cycles <- TRUE
    labs[["lab.cycles"]] <- clines$bottom
    add.cycles <- clines$top
    longer.cycles <- clines$longer
  }
  else {
    newline.cycles <- FALSE
    longer.cycles <- labs[["lab.cycles"]]
  }
  #
  # Check for "\n" in label of column 'n.e'
  #
  clines <- twolines(labs[["lab.n.e"]], "n.e")
  #
  if (clines$newline) {
    newline.n.e <- TRUE
    labs[["lab.n.e"]] <- clines$bottom
    add.n.e <- clines$top
    longer.n.e <- clines$longer
  }
  else {
    newline.n.e <- FALSE
    longer.n.e <- labs[["lab.n.e"]]
  }
  #
  # Check for "\n" in label of column 'n.c'
  #
  clines <- twolines(labs[["lab.n.c"]], "n.c")
  #
  if (clines$newline) {
    newline.n.c <- TRUE
    labs[["lab.n.c"]] <- clines$bottom
    add.n.c <- clines$top
    longer.n.c <- clines$longer
  }
  else {
    newline.n.c <- FALSE
    longer.n.c <- labs[["lab.n.c"]]
  }
  #
  # Check for "\n" in label of column 'event.e'
  #
  clines <- twolines(labs[["lab.event.e"]], "event.e")
  #
  if (clines$newline) {
    newline.event.e <- TRUE
    labs[["lab.event.e"]] <- clines$bottom
    add.event.e <- clines$top
    longer.event.e <- clines$longer
  }
  else {
    newline.event.e <- FALSE
    longer.event.e <- labs[["lab.event.e"]]
  }
  #
  # Check for "\n" in label of column 'event.c'
  #
  clines <- twolines(labs[["lab.event.c"]], "event.c")
  #
  if (clines$newline) {
    newline.event.c <- TRUE
    labs[["lab.event.c"]] <- clines$bottom
    add.event.c <- clines$top
    longer.event.c <- clines$longer
  }
  else {
    newline.event.c <- FALSE
    longer.event.c <- labs[["lab.event.c"]]
  }
  #
  # Check for "\n" in label of column 'mean.e'
  #
  clines <- twolines(labs[["lab.mean.e"]], "mean.e")
  #
  if (clines$newline) {
    newline.mean.e <- TRUE
    labs[["lab.mean.e"]] <- clines$bottom
    add.mean.e <- clines$top
    longer.mean.e <- clines$longer
  }
  else {
    newline.mean.e <- FALSE
    longer.mean.e <- labs[["lab.mean.e"]]
  }
  #
  # Check for "\n" in label of column 'mean.c'
  #
  clines <- twolines(labs[["lab.mean.c"]], "mean.c")
  #
  if (clines$newline) {
    newline.mean.c <- TRUE
    labs[["lab.mean.c"]] <- clines$bottom
    add.mean.c <- clines$top
    longer.mean.c <- clines$longer
  }
  else {
    newline.mean.c <- FALSE
    longer.mean.c <- labs[["lab.mean.c"]]
  }
  #
  # Check for "\n" in label of column 'sd.e'
  #
  clines <- twolines(labs[["lab.sd.e"]], "sd.e")
  #
  if (clines$newline) {
    newline.sd.e <- TRUE
    labs[["lab.sd.e"]] <- clines$bottom
    add.sd.e <- clines$top
    longer.sd.e <- clines$longer
  }
  else {
    newline.sd.e <- FALSE
    longer.sd.e <- labs[["lab.sd.e"]]
  }
  #
  # Check for "\n" in label of column 'sd.c'
  #
  clines <- twolines(labs[["lab.sd.c"]], "sd.c")
  #
  if (clines$newline) {
    newline.sd.c <- TRUE
    labs[["lab.sd.c"]] <- clines$bottom
    add.sd.c <- clines$top
    longer.sd.c <- clines$longer
  }
  else {
    newline.sd.c <- FALSE
    longer.sd.c <- labs[["lab.sd.c"]]
  }
  #
  # Check for "\n" in label of column 'cor'
  #
  clines <- twolines(labs[["lab.cor"]], "cor")
  #
  if (clines$newline) {
    newline.cor <- TRUE
    labs[["lab.cor"]] <- clines$bottom
    add.cor <- clines$top
    longer.cor <- clines$longer
  }
  else {
    newline.cor <- FALSE
    longer.cor <- labs[["lab.cor"]]
  }
  #
  # Check for "\n" in label of column 'time.e'
  #
  clines <- twolines(labs[["lab.time.e"]], "time.e")
  #
  if (clines$newline) {
    newline.time.e <- TRUE
    labs[["lab.time.e"]] <- clines$bottom
    add.time.e <- clines$top
    longer.time.e <- clines$longer
  }
  else {
    newline.time.e <- FALSE
    longer.time.e <- labs[["lab.time.e"]]
  }
  #
  # Check for "\n" in label of column 'time.c'
  #
  clines <- twolines(labs[["lab.time.c"]], "time.c")
  #
  if (clines$newline) {
    newline.time.c <- TRUE
    labs[["lab.time.c"]] <- clines$bottom
    add.time.c <- clines$top
    longer.time.c <- clines$longer
  }
  else {
    newline.time.c <- FALSE
    longer.time.c <- labs[["lab.time.c"]]
  }
  #
  # Check for "\n" in label of column 'effect'
  #
  clines <- twolines(labs[["lab.effect"]], "effect")
  #
  if (clines$newline) {
    newline.effect <- TRUE
    labs[["lab.effect"]] <- clines$bottom
    add.effect <- clines$top
    longer.effect <- clines$longer
  }
  else {
    newline.effect <- FALSE
    longer.effect <- labs[["lab.effect"]]
  }
  #
  # Check for "\n" in label of column 'ci'
  #
  clines <- twolines(labs[["lab.ci"]], "ci")
  #
  if (clines$newline) {
    newline.ci <- TRUE
    labs[["lab.ci"]] <- clines$bottom
    add.ci <- clines$top
    longer.ci <- clines$longer
  }
  else {
    newline.ci <- FALSE
    longer.ci <- labs[["lab.ci"]]
  }
  #
  # Check for "\n" in label of column 'effect.ci'
  #
  clines <- twolines(labs[["lab.effect.ci"]], "effect.ci")
  #
  if (clines$newline) {
    newline.effect.ci <- TRUE
    labs[["lab.effect.ci"]] <- clines$bottom
    add.effect.ci <- clines$top
    longer.effect.ci <- clines$longer
  }
  else {
    newline.effect.ci <- FALSE
    longer.effect.ci <- labs[["lab.effect.ci"]]
  }
  #
  # Check for "\n" in label of column 'event.n.e'
  #
  clines <- twolines(labs[["lab.event.n.e"]], "event.n.e")
  #
  if (clines$newline) {
    newline.event.n.e <- TRUE
    labs[["lab.event.n.e"]] <- clines$bottom
    add.event.n.e <- clines$top
    longer.event.n.e <- clines$longer
  }
  else {
    newline.event.n.e <- FALSE
    longer.event.n.e <- labs[["lab.event.n.e"]]
  }
  #
  # Check for "\n" in label of column 'event.n.c'
  #
  clines <- twolines(labs[["lab.event.n.c"]], "event.n.c")
  #
  if (clines$newline) {
    newline.event.n.c <- TRUE
    labs[["lab.event.n.c"]] <- clines$bottom
    add.event.n.c <- clines$top
    longer.event.n.c <- clines$longer
  }
  else {
    newline.event.n.c <- FALSE
    longer.event.n.c <- labs[["lab.event.n.c"]]
  }
  #
  # Check for "\n" in label of column 'mean.sd.n.e'
  #
  clines <- twolines(labs[["lab.mean.sd.n.e"]], "mean.sd.n.e")
  #
  if (clines$newline) {
    newline.mean.sd.n.e <- TRUE
    labs[["lab.mean.sd.n.e"]] <- clines$bottom
    add.mean.sd.n.e <- clines$top
    longer.mean.sd.n.e <- clines$longer
  }
  else {
    newline.mean.sd.n.e <- FALSE
    longer.mean.sd.n.e <- labs[["lab.mean.sd.n.e"]]
  }
  #
  # Check for "\n" in label of column 'mean.sd.n.c'
  #
  clines <- twolines(labs[["lab.mean.sd.n.c"]], "mean.sd.n.c")
  #
  if (clines$newline) {
    newline.mean.sd.n.c <- TRUE
    labs[["lab.mean.sd.n.c"]] <- clines$bottom
    add.mean.sd.n.c <- clines$top
    longer.mean.sd.n.c <- clines$longer
  }
  else {
    newline.mean.sd.n.c <- FALSE
    longer.mean.sd.n.c <- labs[["lab.mean.sd.n.c"]]
  }
  #
  # Check for "\n" in label of column 'w.common'
  #
  clines <- twolines(labs[["lab.w.common"]], "w.common")
  #
  if (clines$newline) {
    newline.w.common <- TRUE
    labs[["lab.w.common"]] <- clines$bottom
    add.w.common <- clines$top
    longer.w.common <- clines$longer
  }
  else {
    newline.w.common <- FALSE
    longer.w.common <- labs[["lab.w.common"]]
  }
  #
  # Check for "\n" in label of column 'w.random'
  #
  clines <- twolines(labs[["lab.w.random"]], "w.random")
  #
  if (clines$newline) {
    newline.w.random <- TRUE
    labs[["lab.w.random"]] <- clines$bottom
    add.w.random <- clines$top
    longer.w.random <- clines$longer
  }
  else {
    newline.w.random <- FALSE
    longer.w.random <- labs[["lab.w.random"]]
  }
  #
  # Check for "\n" in label of column 'pval'
  #
  clines <- twolines(labs[["lab.pval"]], "pval")
  #
  if (clines$newline) {
    newline.pval <- TRUE
    labs[["lab.pval"]] <- clines$bottom
    add.pval <- clines$top
    longer.pval <- clines$longer
  }
  else {
    newline.pval <- FALSE
    longer.pval <- labs[["lab.pval"]]
  }
  #
  # Check for "\n" in label of column 'tau2'
  #
  clines <- twolines(labs[["lab.tau2"]], "tau2")
  #
  if (clines$newline) {
    newline.tau2 <- TRUE
    labs[["lab.tau2"]] <- clines$bottom
    add.tau2 <- clines$top
    longer.tau2 <- clines$longer
  }
  else {
    newline.tau2 <- FALSE
    longer.tau2 <- labs[["lab.tau2"]]
  }
  #
  # Check for "\n" in label of column 'tau'
  #
  clines <- twolines(labs[["lab.tau"]], "tau")
  #
  if (clines$newline) {
    newline.tau <- TRUE
    labs[["lab.tau"]] <- clines$bottom
    add.tau <- clines$top
    longer.tau <- clines$longer
  }
  else {
    newline.tau <- FALSE
    longer.tau <- labs[["lab.tau"]]
  }
  #
  # Check for "\n" in label of column 'I2'
  #
  clines <- twolines(labs[["lab.I2"]], "I2")
  #
  if (clines$newline) {
    newline.I2 <- TRUE
    labs[["lab.I2"]] <- clines$bottom
    add.I2 <- clines$top
    longer.I2 <- clines$longer
  }
  else {
    newline.I2 <- FALSE
    longer.I2 <- labs[["lab.I2"]]
  }
  #
  # Check for "\n" in argument 'smlab'
  #
  clines <- twolines(smlab, arg = TRUE)
  #
  if (clines$newline) {
    smlab1 <- clines$top
    smlab2 <- clines$bottom
    #
    newline.smlab <- TRUE
  }
  else {
    smlab1 <- smlab
    smlab2 <- ""
    #
    newline.smlab <- FALSE
  }
  #
  # Check for "\n" in argument 'xlab'
  #
  clines <- twolines(xlab, arg = TRUE)
  #
  if (clines$newline) {
    newline.xlab <- TRUE
    xlab <- clines$top
    xlab.add <- clines$bottom
  }
  else {
    xlab.add <- ""
    newline.xlab <- FALSE
  }
  #
  # Check for "\n" in argument 'label.left'
  #
  clines <- twolines(label.left, arg = TRUE)
  #
  if (clines$newline) {
    ll1 <- clines$top
    ll2 <- clines$bottom
    #
    newline.ll <- TRUE
  }
  else {
    ll1 <- label.left
    ll2 <- ""
    #
    newline.ll <- FALSE
  }
  #
  # Check for "\n" in argument 'label.right'
  #
  clines <- twolines(label.right, arg = TRUE)
  #
  if (clines$newline) {
    lr1 <- clines$top
    lr2 <- clines$bottom
    #
    newline.lr <- TRUE
  }
  else {
    lr1 <- label.right
    lr2 <- ""
    #
    newline.lr <- FALSE
  }
  #
  # Check for "\n" in additional columns
  #
  newline.addcol.left  <- FALSE
  newline.addcol.right <- FALSE
  #
  if (newcols) {
    if (length(leftcols.new) > 0) {
      for (i in seq_along(leftcols.new)) {
        # Check for "\n" in label of new column
        clines <- twolines(leftlabs.new[i], leftcols.new[i])
        newline.addcol.left <- c(newline.addcol.left, clines$newline)
      }
      newline.addcol.right <- sum(newline.addcol.right) > 0
    }
    if (length(leftcols.new) > 0) {
      for (i in seq_along(leftcols.new)) {
        # Check for "\n" in label of new column
        clines <- twolines(leftlabs.new[i], leftcols.new[i])
        newline.addcol.left <- c(newline.addcol.left, clines$newline)
      }
      newline.addcol.left <- sum(newline.addcol.left) > 0
    }
  }
  #
  newline <- newline.studlab | newline.effect |
    newline.ci | newline.effect.ci |
    newline.w.common | newline.w.random | newline.TE | newline.seTE |
    newline.cluster | newline.cycles |
    newline.n.e | newline.n.c | newline.event.e | newline.event.c |
    newline.mean.e | newline.mean.c | newline.sd.e | newline.sd.c |
    newline.cor | newline.time.e | newline.time.c |
    newline.pval | newline.tau2 | newline.tau | newline.I2 |
    newline.smlab | newline.addcol.left | newline.addcol.right
  #
  newline.all <- newline | (!newline & (newline.ll | newline.lr) & !addrow)
  
  
  #
  #
  # (10) Define columns in forest plot as well as x- and y-limits
  #
  #
  type.common <-
    deprecated(type.common, missing.type.common,
               args, "type.common",
               warn.deprecated)
  #
  type.subgroup.common <-
    deprecated(type.subgroup.common, missing.type.subgroup.common,
               args, "type.subgroup.common",
               warn.deprecated)
  #
  col.inside.common <-
    deprecated(col.inside.common, missing.col.inside.common,
               args, "col.inside.common",
               warn.deprecated)
  #
  # Set / check correct length of variables
  #
  col.diamond.common <-
    setlength(col.diamond.common, n.com, "number of common effect estimates")
  #
  col.diamond.random <-
    setlength(col.diamond.random, n.ran, "number of random effect estimates")
  #
  col.predict <-
    setlength(col.predict, n.prd, "number of prediction intervals")
  #
  col.diamond.lines.common <-
    setlength(col.diamond.lines.common, n.com,
              "number of common effect estimates")
  #
  col.diamond.lines.random <-
    setlength(col.diamond.lines.random, n.ran,
              "number of random effect estimates")
  #
  col.predict.lines <-
    setlength(col.predict.lines, n.prd, "number of prediction intervals")
  #
  col.inside.common <-
    setlength(col.inside.common, n.com, "number of common effect estimates")
  #
  col.inside.random <-
    setlength(col.inside.random, n.ran, "number of random effect estimates")
  #
  if (by) {
    #
    subgroup.name <-
      bylabel(subgroup.name, subgroup.levels, print.subgroup.name,
              sep.subgroup, big.mark = big.mark)
    #
    wrong.common <- FALSE
    wrong.random <- FALSE
    wrong.predict <- FALSE
    #
    if (n.by > 1) {
      if (length(text.common.w) == 1)
        text.common.w <- rep(text.common.w, n.com.w)
      else if (length(text.common.w) == n.com)
        text.common.w <- rep(text.common.w, n.by)
      else if (length(text.common.w) == n.by)
        text.common.w <- rep(text.common.w, n.com)
      else
        wrong.common <- TRUE
      #
      if (length(text.random.w) == 1)
        text.random.w <- rep(text.random.w, n.ran.w)
      else if (length(text.random.w) == n.ran)
        text.random.w <- rep(text.random.w, n.by)
      else if (length(text.random.w) == n.by)
        text.random.w <- rep(text.random.w, n.ran)
      else
        wrong.random <- TRUE
      #
      if (length(text.predict.w) == 1)
        text.predict.w <- rep(text.predict.w, n.prd.w)
      else if (length(text.predict.w) == n.prd)
        text.predict.w <- rep(text.predict.w, n.by)
      else if (length(text.predict.w) == n.by)
        text.predict.w <- rep(text.predict.w, n.prd)
      else
        wrong.predict <- TRUE
    }
    else if (n.by == 1) {
      if (length(text.common.w) != n.com)
        wrong.common <- TRUE
      if (length(text.random.w) != n.ran)
        wrong.random <- TRUE
      if (length(text.predict.w) != n.prd)
        wrong.predict <- TRUE
    }
    #
    if (wrong.common)
      stop("Argument 'text.common.w' must be a single character string",
           if (n.by > 1 | n.com > 1)
             paste0(" or of length\n  -", n.by, " (number of subgroups) or",
                    " \n  -", n.com, " (number of common effect analyses)")
           else
             ".")
    #
    if (wrong.random)
      stop("Argument 'text.random.w' must be a single character string",
           if (n.by > 1 | n.ran > 1)
             paste0(" or of length\n  -", n.by, " (number of subgroups) or",
                    " \n  -", n.ran, " (number of random effects analyses)")
           else
             ".")
    #
    if (wrong.predict)
      stop("Argument 'text.predict.w' must be a single character string",
           if (n.by > 1 | n.prd > 1)
             paste0(" or of length\n  -", n.by, " (number of subgroups) or",
                    " \n  -", n.prd, " (number of prediction intervals)")
           else
             ".")
    #
    if (hetstat.pooled == "common") {
      text.common <- hetstat.overall
      text.common.w <- hetstat.w
      hetstat <- FALSE
    }
    else if (hetstat.pooled == "random") {
      text.random <- hetstat.overall
      text.random.w <- hetstat.w
      hetstat <- FALSE
    }
    else if (hetstat.pooled == "study") {
      studlab <- hetstat.w
      hetstat <- FALSE
    }
    #
    modlabs <- c(text.common, text.random, text.predict,
                 hetstat.overall, hetstat.resid,
                 text.overall.common, text.overall.random,
                 text.subgroup.common, text.subgroup.random,
                 text.addline1, text.addline2,
                 if (details) text.details else NULL,
                 if (RoB.legend) text.rob else NULL,
                 subgroup.name, text.common.w, text.random.w, text.predict.w,
                 hetstat.w,
                 unlist(text.effect.subgroup.common),
                 unlist(text.effect.subgroup.random),
                 studlab)
    #
    TEs    <- c(TE.common, TE.random, NAs.prd, TE.w, TE)
    lowTEs <- c(lowTE.common, lowTE.random, lowTE.predict, lowTE.w, lowTE)
    uppTEs <- c(uppTE.common, uppTE.random, uppTE.predict, uppTE.w, uppTE)
    #
    TEs.exclude <- c(TE.common, TE.random, NAs.prd, TE.w, TE.exclude)
    lowTEs.exclude <- c(lowTE.common, lowTE.random, lowTE.predict, lowTE.w,
                        lowTE.exclude)
    uppTEs.exclude <- c(uppTE.common, uppTE.random, uppTE.predict, uppTE.w,
                        uppTE.exclude)
    #
    TEs.study <- c(blanks, blanks.w,
                   formatN(TE.orig, digits.TE, lab.NA, big.mark = big.mark,
                           monospaced = monospaced))
    #
    seTEs.study <- c(blanks, blanks.w,
                     formatN(seTE, digits.se, lab.NA, big.mark = big.mark,
                             monospaced = monospaced))
    #
    w.commons <- c(NAs.com, NAs.ran, NAs.prd,
                   NAs.com.w, NAs.ran.w, NAs.prd.w,
                   NAs.stat.w,
                   w.common.p)
    w.randoms <- c(NAs.com, NAs.ran, NAs.prd,
                   NAs.com.w, NAs.ran.w, NAs.prd.w,
                   NAs.stat.w,
                   w.random.p)
    #
    format.w.commons <- c(rep(100, n.com), NAs.ran, NAs.prd,
                          repl(w.common.w.p, n.com, n.by),
                          NAs.ran.w, NAs.prd.w, NAs.stat.w,
                          w.common.p)
    format.w.randoms <- c(NAs.com, rep(100, n.ran), NAs.prd,
                          NAs.com.w,
                          repl(w.random.w.p, n.ran, n.by),
                          NAs.prd.w, NAs.stat.w,
                          w.random.p)
    #
    Wc.format <-
      formatN(format.w.commons, digits.weight, text.NA = lab.NA.weight,
              monospaced = monospaced)
    Wr.format <-
      formatN(format.w.randoms, digits.weight, text.NA = lab.NA.weight,
              monospaced = monospaced)
    #
    sel.common <- Wc.format == lab.NA.weight
    sel.random <- Wr.format == lab.NA.weight
    #
    sel.common[first.ran.w] <- TRUE
    sel.random[first.com.w] <- TRUE
    #
    type.pooled <- c(rep(type.common, n.com),
                     rep(type.random, n.ran),
                     rep("predict", n.prd),
                     rep(type.subgroup.common, n.com.w),
                     rep(type.subgroup.random, n.ran.w),
                     rep("predict", n.prd.w),
                     blanks.stat.w)
    #
    col.diamond.pooled <-
      c(col.diamond.common, col.diamond.random, col.predict,
        rep(col.diamond.common, n.by),
        rep(col.diamond.random, n.by),
        rep(col.predict, n.by),
        NAs.stat.w)
    #
    col.diamond.lines.pooled <-
      c(col.diamond.lines.common, col.diamond.lines.random, col.predict.lines,
        rep(col.diamond.lines.common, n.by),
        rep(col.diamond.lines.random, n.by),
        rep(col.predict.lines, n.by),
        NAs.stat.w)
    #
    col.inside.pooled <-
      c(col.inside.common, col.inside.random, blanks.prd,
        rep(col.inside.common, n.by),
        rep(col.inside.random, n.by),
        rep("", n.prd.w),
        NAs.stat.w)
  }
  else {
    #
    if (hetstat.pooled == "common")
      text.common <- hetstat.overall
    else if (hetstat.pooled == "random")
      text.random <- hetstat.overall
    #
    modlabs <- c(text.common, text.random, text.predict,
                 hetstat.overall, hetstat.resid,
                 text.overall.common, text.overall.random,
                 text.subgroup.common, text.subgroup.random,
                 text.addline1, text.addline2,
                 if (details) text.details else NULL,
                 if (RoB.legend) text.rob else NULL,
                 studlab)
    #
    TEs    <- c(TE.common, TE.random, NAs.prd, TE)
    lowTEs <- c(lowTE.common, lowTE.random, lowTE.predict, lowTE)
    uppTEs <- c(uppTE.common, uppTE.random, uppTE.predict, uppTE)
    #
    TEs.exclude    <- c(TE.common, TE.random, NAs.prd, TE.exclude)
    lowTEs.exclude <- c(lowTE.common, lowTE.random,
                        lowTE.predict, lowTE.exclude)
    uppTEs.exclude <- c(uppTE.common, uppTE.random,
                        uppTE.predict, uppTE.exclude)
    #
    TEs.study <- c(blanks,
                   formatN(TE.orig, digits.TE, lab.NA, big.mark = big.mark))
    seTEs.study <- c(blanks,
                     formatN(seTE, digits.se, lab.NA, big.mark = big.mark))
    #
    w.commons <- c(NAs, w.common.p)
    w.randoms <- c(NAs, w.random.p)
    #
    format.w.common <- formatN(c(100, w.common.p), digits.weight,
                               monospaced = monospaced)
    format.w.random <- formatN(c(100, w.random.p), digits.weight,
                               monospaced = monospaced)
    #
    Wc.format <- c(rep(format.w.common[1], n.com),
                   rep(lab.NA.weight, n.ran),
                   blanks.prd,
                   format.w.common[-1])
    Wr.format <- c(rep(lab.NA.weight, n.com),
                   rep(format.w.random[1], n.ran),
                   blanks.prd,
                   format.w.random[-1])
    #
    sel.common <- Wc.format == lab.NA.weight
    sel.random <- Wr.format == lab.NA.weight
    #
    type.pooled <- c(rep(type.common, n.com),
                     rep(type.random, n.ran),
                     rep("predict", n.prd))
    #
    col.diamond.pooled <-
      c(col.diamond.common, col.diamond.random, col.predict)
    #
    col.diamond.lines.pooled <-
      c(col.diamond.lines.common, col.diamond.lines.random, col.predict.lines)
    #
    col.inside.pooled <-
      c(col.inside.common, col.inside.random, blanks.prd)
  }
  #
  # Treatment effect and confidence interval
  #
  if (backtransf & log.xaxis) {
    effect.format <-
      formatN(scale * exp(TEs / scale),
              digits, lab.NA.effect, big.mark = big.mark)
    ci.format <-
      ifelse(is.na(lowTEs) | is.na(uppTEs), lab.NA.effect,
             formatCI(formatN(scale * exp(lowTEs / scale),
                              digits = digits, text.NA = lab.NA,
                              big.mark = big.mark),
                      formatN(scale * exp(uppTEs / scale),
                              digits = digits, text.NA = lab.NA,
                              big.mark = big.mark)))
  }
  else {
    effect.format <- formatN(TEs, digits, lab.NA.effect, big.mark = big.mark,
                             monospaced = monospaced)
    ci.format <-
      ifelse(is.na(lowTEs) | is.na(uppTEs), lab.NA.effect,
             formatCI(formatN(lowTEs, digits = digits, text.NA = lab.NA,
                              big.mark = big.mark,
                              monospaced = monospaced),
                      formatN(uppTEs, digits = digits, text.NA = lab.NA,
                              big.mark = big.mark,
                              monospaced = monospaced)))
  }
  #
  effect.ci.format <- paste0(effect.format,
                             ifelse(is.na(TEs), "", " "),
                             ci.format)
  #
  # No treatment effect for prediction interval
  #
  effect.format[all.prd] <- ""
  #
  # Only print prediction interval if requested
  #
  if (!prediction) {
    ci.format[all.prd] <- ""
    effect.ci.format[all.prd] <- ""
  }
  #
  # No treatment effect and confidence interval in statistics lines
  #
  if (by) {
    effect.format[all.stat.w] <- ""
    ci.format[all.stat.w] <- ""
    effect.ci.format[all.stat.w] <- ""
  }
  #
  # Weights of common and random effects model
  #
  Wc.format <- paste0(Wc.format, if (w.common.percent & !bmj) "%")
  Wr.format <- paste0(Wr.format, if (w.random.percent & !bmj) "%")
  #
  Wc.format[Wc.format == "%"] <- ""
  Wr.format[Wr.format == "%"] <- ""
  #
  Wc.format[sel.common] <- lab.NA.weight
  Wr.format[sel.random] <- lab.NA.weight
  #
  # Treatment estimate and its standard error
  #
  TE.format <- TEs.study
  seTE.format <- seTEs.study
  #
  # Number of patients, events, and person times
  #
  if (!is.null(x$n.e.pooled))
    sum.n.e <- x$n.e.pooled
  else {
    if (avail.exclude)
      sum.n.e <- sum(x$n.e[!x$exclude], na.rm = TRUE)
    else
      sum.n.e <- sum(x$n.e, na.rm = TRUE)
  }
  #
  if (!is.null(x$n.c.pooled))
    sum.n.c <- x$n.c.pooled
  else {
    if (avail.exclude)
      sum.n.c <- sum(x$n.c[!x$exclude], na.rm = TRUE)
    else
      sum.n.c <- sum(x$n.c, na.rm = TRUE)
  }
  #
  if (!is.null(x$event.e.pooled))
    sum.e.e <- x$event.e.pooled
  else {
    if (avail.exclude)
      sum.e.e <- sum(x$event.e[!x$exclude], na.rm = TRUE)
    else
      sum.e.e <- sum(x$event.e, na.rm = TRUE)
  }
  #
  if (!is.null(x$event.c.pooled))
    sum.e.c <- x$event.c.pooled
  else {
    if (avail.exclude)
      sum.e.c <- sum(x$event.c[!x$exclude], na.rm = TRUE)
    else
      sum.e.c <- sum(x$event.c, na.rm = TRUE)
  }
  #
  if (!is.null(x$time.e.pooled))
    sum.t.e <- x$time.e.pooled
  else {
    if (avail.exclude)
      sum.t.e <- sum(x$time.e[!x$exclude], na.rm = TRUE)
    else
      sum.t.e <- sum(x$time.e, na.rm = TRUE)
  }
  #
  if (!is.null(x$time.c.pooled))
    sum.t.c <- x$time.c.pooled
  else {
    if (avail.exclude)
      sum.t.c <- sum(x$time.c[!x$exclude], na.rm = TRUE)
    else
      sum.t.c <- sum(x$time.c, na.rm = TRUE)
  }
  #
  if (is.character(x$cluster))
    as.character.cluster <- TRUE
  else if (is.factor(x$cluster)) {
    if (anyNA(suppressWarnings(as.numeric(as.character(x$cluster)))))
      as.character.cluster <- TRUE
    else
      as.character.cluster <- FALSE
    #
    x$cluster <- as.character(x$cluster)
  }
  else
    as.character.cluster <- FALSE
  #
  as.character.cycles <- FALSE
  #
  if (by) {
    if (pooled.totals) {
      Ne <- c(c(sum.n.e, NAs.com1),
              c(sum.n.e, NAs.ran1),
              NAs.prd,
              repl(n.e.w, n.com, n.by),
              repl(n.e.w, n.ran, n.by),
              NAs.prd.w, NAs.stat.w,
              x$n.e)
      Nc <- c(c(sum.n.c, NAs.com1),
              c(sum.n.c, NAs.ran1),
              NAs.prd,
              repl(n.c.w, n.com, n.by),
              repl(n.c.w, n.ran, n.by),
              NAs.prd.w, NAs.stat.w,
              x$n.c)
    }
    else {
      Ne <- c(NAs.all, x$n.e)
      Nc <- c(NAs.all, x$n.c)
    }
    if (pooled.events) {
      Ee <- c(c(sum.e.e, NAs.com1),
              c(sum.e.e, NAs.ran1),
              NAs.prd,
              repl(e.e.w, n.com, n.by),
              repl(e.e.w, n.ran, n.by),
              NAs.prd.w, NAs.stat.w,
              x$event.e)
      Ec <- c(c(sum.e.c, NAs.com1),
              c(sum.e.c, NAs.ran1),
              NAs.prd,
              repl(e.c.w, n.com, n.by),
              repl(e.c.w, n.ran, n.by),
              NAs.prd.w, NAs.stat.w,
              x$event.c)
    }
    else {
      Ee <- c(NAs.all, x$event.e)
      Ec <- c(NAs.all, x$event.c)
    }
    if (pooled.times) {
      Te <- c(c(sum.t.e, NAs.com1),
              c(sum.t.e, NAs.ran1),
              NAs.prd,
              repl(t.e.w, n.com, n.by),
              repl(t.e.w, n.ran, n.by),
              NAs.prd.w, NAs.stat.w,
              x$time.e)
      Tc <- c(c(sum.t.c, NAs.com1),
              c(sum.t.c, NAs.ran1),
              NAs.prd,
              repl(t.c.w, n.com, n.by),
              repl(t.c.w, n.ran, n.by),
              NAs.prd.w, NAs.stat.w,
              x$time.c)
    }
    else {
      Te <- c(NAs.all, x$time.e)
      Tc <- c(NAs.all, x$time.c)
    }
  }
  else {
    if (pooled.totals) {
      Ne <- c(c(sum.n.e, NAs.com1),
              c(sum.n.e, NAs.ran1),
              NAs.prd,
              x$n.e)
      Nc <- c(c(sum.n.c, NAs.com1),
              c(sum.n.c, NAs.ran1),
              NAs.prd,
              x$n.c)
    }
    else {
      Ne <- c(NAs,
              x$n.e)
      Nc <- c(NAs,
              x$n.c)
    }
    if (pooled.events) {
      Ee <- c(c(sum.e.e, NAs.com1),
              c(sum.e.e, NAs.ran1),
              NAs.prd,
              x$event.e)
      Ec <- c(c(sum.e.c, NAs.com1),
              c(sum.e.c, NAs.ran1),
              NAs.prd,
              x$event.c)
    }
    else {
      Ee <- c(NAs,
              x$event.e)
      Ec <- c(NAs,
              x$event.c)
    }
    if (pooled.times) {
      Te <- c(c(sum.t.e, NAs.com1),
              c(sum.t.e, NAs.ran1),
              NAs.prd,
              x$time.e)
      Tc <- c(c(sum.t.c, NAs.com1),
              c(sum.t.c, NAs.ran1),
              NAs.prd,
              x$time.c)
    }
    else {
      Te <- c(NAs,
              x$time.e)
      Tc <- c(NAs,
              x$time.c)
    }
  }
  #
  if (by) {
    bla <- rep("", length(Ne))
    scom.w <- sran.w <- sprd.w <- bla
    scom.w[first.com.w] <- "yes"
    sran.w[first.ran.w] <- "yes"
    sprd.w[first.prd.w] <- "yes"
    #
    ncom.w <- nran.w <- nprd.w <- notfirst <- bla
    ncom.w[emp.com.w] <- "yes"
    nran.w[emp.ran.w] <- "yes"
    nprd.w[emp.prd.w] <- "yes"
    notfirst[emp.w] <- "yes"
    #
    if (FALSE)
      print(data.frame(TEs = round(TEs, 2),
                       lowTEs = round(lowTEs, 2),
                       uppTEs = round(uppTEs, 2)
                       ))
    if (FALSE)
      print(data.frame(row = c(paste0("C", seq_len(n.com)),
                               paste0("R", seq_len(n.ran)),
                               paste0("P", seq_len(n.prd)),
                               paste0("C.w", seq_len(n.com.w)),
                               paste0("R.w", seq_len(n.ran.w)),
                               paste0("P.w", seq_len(n.prd.w)),
                               paste0("het.w", seq_len(n.by)),
                               paste0("efc.w", seq_len(n.by)),
                               paste0("efr.w", seq_len(n.by)),
                               paste0("TE", seq_len(k.all))),
                       Ne, notfirst))
  }
  Ne.format <-
    formatN(Ne, digits = digits.n, text.NA = lab.NA, big.mark = big.mark,
            monospaced = monospaced)
  Nc.format <-
    formatN(Nc, digits = digits.n, text.NA = lab.NA, big.mark = big.mark,
            monospaced = monospaced)
  Ee.format <-
    formatN(Ee, digits = digits.event, text.NA = lab.NA, big.mark = big.mark,
            monospaced = monospaced)
  Ec.format <-
    formatN(Ec, digits = digits.event, text.NA = lab.NA, big.mark = big.mark,
            monospaced = monospaced)
  #
  if (all(is_wholenumber(Te), na.rm = TRUE) & missing.digits.time)
    Te.format <-
      formatN(Te, digits = 0, text.NA = lab.NA, big.mark = big.mark)
  else
    Te.format <- formatN(Te, digits.time, lab.NA, big.mark = big.mark)
  #
  if (all(is_wholenumber(Tc), na.rm = TRUE) & missing.digits.time)
    Tc.format <-
      formatN(Tc, digits = 0, text.NA = lab.NA, big.mark = big.mark)
  else
    Tc.format <- formatN(Tc, digits.time, lab.NA, big.mark = big.mark)
  #
  # Print nothing in lines with prediction interval
  #
  Ne.format[all.prd] <- Nc.format[all.prd] <- ""
  Ee.format[all.prd] <- Ec.format[all.prd] <- ""
  Te.format[all.prd] <- Tc.format[all.prd] <- ""
  Wc.format[all.prd] <- Wr.format[all.prd] <- ""
  Wc.format[emp.com] <- Wr.format[emp.com] <- ""
  Wc.format[emp.ran] <- Wr.format[emp.ran] <- ""
  #
  if (by) {
    #
    # Print nothing in selected subgroup lines
    #
    Ne.format[emp.w] <- Nc.format[emp.w] <- ""
    Ee.format[emp.w] <- Ec.format[emp.w] <- ""
    Te.format[emp.w] <- Tc.format[emp.w] <- ""
    Wc.format[emp.w] <- Wr.format[emp.w] <- ""
  }
  #
  if (common.random) {
    #
    # Print nothing in lines with results for random effects model
    #
    Ne.format[all.ran] <- Nc.format[all.ran] <- ""
    Ee.format[all.ran] <- Ec.format[all.ran] <- ""
    Te.format[all.ran] <- Tc.format[all.ran] <- ""
    #
    if (by) {
      Ne.format[all.ran.w] <- Nc.format[all.ran.w] <- ""
      Ee.format[all.ran.w] <- Ec.format[all.ran.w] <- ""
      Te.format[all.ran.w] <- Tc.format[all.ran.w] <- ""
    }
  }
  #
  # Print nothing in lines for second, third etc. common effect or
  # random effects model
  #
  Ne.format[emp.com] <- Nc.format[emp.com] <- ""
  Ee.format[emp.com] <- Ec.format[emp.com] <- ""
  Te.format[emp.com] <- Tc.format[emp.com] <- ""
  #
  Ne.format[emp.ran] <- Nc.format[emp.ran] <- ""
  Ee.format[emp.ran] <- Ec.format[emp.ran] <- ""
  Te.format[emp.ran] <- Tc.format[emp.ran] <- ""
  #
  if (by) {
    Ne.format[emp.com.w] <- Nc.format[emp.com.w] <- ""
    Ee.format[emp.com.w] <- Ec.format[emp.com.w] <- ""
    Te.format[emp.com.w] <- Tc.format[emp.com.w] <- ""
    #
    Ne.format[emp.ran.w] <- Nc.format[emp.ran.w] <- ""
    Ee.format[emp.ran.w] <- Ec.format[emp.ran.w] <- ""
    Te.format[emp.ran.w] <- Tc.format[emp.ran.w] <- ""
  }
  #
  # Only print samples sizes if pooled.totals is TRUE
  #
  if (!pooled.totals) {
    Ne.format[all.com] <- Nc.format[all.com] <- ""
    Ne.format[all.ran] <- Nc.format[all.ran] <- ""
    #
    if (by) {
      Ne.format[all.com.w] <- Nc.format[all.com.w] <- ""
      Ne.format[all.ran.w] <- Nc.format[all.ran.w] <- ""
    }
  }
  #
  # Only print total number of events if pooled.events is TRUE
  #
  if (!pooled.events) {
    Ee.format[all.com] <- Ec.format[all.com] <- ""
    Ee.format[all.ran] <- Ec.format[all.ran] <- ""
    #
    if (by) {
      Ee.format[all.com.w] <- Ec.format[all.com.w] <- ""
      Ee.format[all.ran.w] <- Ec.format[all.ran.w] <- ""
    }
  }
  #
  # Only print total person times if pooled.times is TRUE
  #
  if (!pooled.times) {
    Te.format[all.com] <- Tc.format[all.com] <- ""
    Te.format[all.ran] <- Tc.format[all.ran] <- ""
    #
    if (by) {
      Te.format[all.com.w] <- Tc.format[all.com.w] <- ""
      Te.format[all.ran.w] <- Tc.format[all.ran.w] <- ""
    }
  }
  #
  # Mean and standard deviation
  #
  if (by) {
    Me <- c(NAs.all, x$mean.e)
    Mc <- c(NAs.all, x$mean.c)
    Se <- c(NAs.all, x$sd.e)
    Sc <- c(NAs.all, x$sd.c)
  }
  else {
    Me <- c(NAs, x$mean.e)
    Mc <- c(NAs, x$mean.c)
    Se <- c(NAs, x$sd.e)
    Sc <- c(NAs, x$sd.c)
  }
  #
  digits.R <- options()$digits
  #
  if (is.null(digits.mean)) {
    if (all(is_wholenumber(Me), na.rm = TRUE))
      Me.format <-
        formatN(Me, digits = 0, text.NA = lab.NA, big.mark = big.mark,
                monospaced = monospaced)
    else
      Me.format <-
        formatN(Me, digits = digits.R, text.NA = lab.NA, big.mark = big.mark,
                monospaced = monospaced)
    #
    if (all(is_wholenumber(Mc), na.rm = TRUE))
      Mc.format <-
        formatN(Mc, digits = 0, text.NA = lab.NA, big.mark = big.mark,
                monospaced = monospaced)
    else
      Mc.format <-
        formatN(Mc, digits = digits.R, text.NA = lab.NA, big.mark = big.mark,
                monospaced = monospaced)
  }
  else {
    Me.format <- formatN(Me, digits.mean, lab.NA, big.mark = big.mark,
                         monospaced = monospaced)
    Mc.format <- formatN(Mc, digits.mean, lab.NA, big.mark = big.mark,
                         monospaced = monospaced)
  }
  #
  if (is.null(digits.sd)) {
    if (all(is_wholenumber(Se), na.rm = TRUE))
      Se.format <-
        formatN(Se, digits = 0, text.NA = lab.NA, big.mark = big.mark,
                monospaced = monospaced)
    else
      Se.format <-
        formatN(Se, digits = digits.R, text.NA = lab.NA, big.mark = big.mark,
                monospaced = monospaced)
    #
    if (all(is_wholenumber(Sc), na.rm = TRUE))
      Sc.format <-
        formatN(Sc, digits = 0, text.NA = lab.NA, big.mark = big.mark,
                monospaced = monospaced)
    else
      Sc.format <-
        formatN(Sc, digits = digits.R, text.NA = lab.NA, big.mark = big.mark,
                monospaced = monospaced)
  }
  else {
    Se.format <- formatN(Se, digits.sd, lab.NA, big.mark = big.mark,
                         monospaced = monospaced)
    Sc.format <- formatN(Sc, digits.sd, lab.NA, big.mark = big.mark,
                         monospaced = monospaced)
  }
  #
  # Print nothing for lines with summary results
  #
  Me.format[all.res] <- Mc.format[all.res] <-
    Se.format[all.res] <- Sc.format[all.res] <- ""
  #
  if (by) {
    Me.format[sel.w] <- Mc.format[sel.w] <- ""
    Se.format[sel.w] <- Sc.format[sel.w] <- ""
  }
  #
  if (ev.n.bin) {
    Ee.bmj <- Ee.format
    Ec.bmj <- Ec.format
    #
    Ne.bmj <- Ne.format
    Nc.bmj <- Nc.format
    #
    Ne.bmj <- ifelse(Ee.bmj == "", "", Ne.bmj)
    Nc.bmj <- ifelse(Ec.bmj == "", "", Nc.bmj)
    #
    while(all(substring(Ee.bmj[Ee.bmj != ""], 1, 1) == " "))
      Ee.bmj[Ee.bmj != ""] <- substring(Ee.bmj[Ee.bmj != ""], 2)
    #
    while(all(substring(Ec.bmj[Ec.bmj != ""], 1, 1) == " "))
      Ec.bmj[Ec.bmj != ""] <- substring(Ec.bmj[Ec.bmj != ""], 2)
    #
    while(all(substring(Ne.bmj[Ne.bmj != ""], 1, 1) == " "))
      Ne.bmj[Ne.bmj != ""] <- substring(Ne.bmj[Ne.bmj != ""], 2)
    #
    while(all(substring(Nc.bmj[Nc.bmj != ""], 1, 1) == " "))
      Nc.bmj[Nc.bmj != ""] <- substring(Nc.bmj[Nc.bmj != ""], 2)
    #
    Ne.bmj <- rmSpace(Ne.bmj, end = FALSE)
    if (monospaced) {
      nchar.Ne.bmj <- nchar(Ne.bmj)
      Ne.bmj[nchar.Ne.bmj > 0] <-
        str_pad(Ne.bmj[nchar.Ne.bmj > 0], width = max(nchar.Ne.bmj),
                side = "right")
    }
    #
    Nc.bmj <- rmSpace(Nc.bmj, end = FALSE)
    if (monospaced) {
      nchar.Nc.bmj <- nchar(Nc.bmj)
      Nc.bmj[nchar.Nc.bmj > 0] <-
        str_pad(Nc.bmj[nchar.Nc.bmj > 0], width = max(nchar.Nc.bmj),
                side = "right")
    }
    #
    event.n.e.format <-
      ifelse(Ee.bmj != "",
             paste0(Ee.bmj, bmj.sep, Ne.bmj),
             "")
    event.n.c.format <-
      ifelse(Ec.format != "",
             paste0(Ec.bmj, bmj.sep, Nc.bmj),
             "")
    #
    if (pooled.totals & !pooled.events) {
      event.n.e.format[Ne.format != "" & Ee.format == ""] <-
        Ne.format[Ne.format != "" & Ee.format == ""]
      event.n.c.format[Nc.format != "" & Ec.format == ""] <-
        Nc.format[Nc.format != "" & Ec.format == ""]
    }
  }
  #
  if (m.s.n.cont) {
    Ne.bmj <- Ne.format
    Nc.bmj <- Nc.format
    #
    while(all(substring(Ne.bmj[Ne.bmj != ""], 1, 1) == " "))
      Ne.bmj[Ne.bmj != ""] <- substring(Ne.bmj[Ne.bmj != ""], 2)
    #
    while(all(substring(Nc.bmj[Nc.bmj != ""], 1, 1) == " "))
      Nc.bmj[Nc.bmj != ""] <- substring(Nc.bmj[Nc.bmj != ""], 2)
    #
    Ne.bmj <- rmSpace(Ne.bmj, end = FALSE)
    if (monospaced) {
      nchar.Ne.bmj <- nchar(Ne.bmj)
      Ne.bmj[nchar.Ne.bmj > 0] <-
        str_pad(Ne.bmj[nchar.Ne.bmj > 0], width = max(nchar.Ne.bmj),
                side = "right")
    }
    #
    Nc.bmj <- rmSpace(Nc.bmj, end = FALSE)
    if (monospaced) {
      nchar.Nc.bmj <- nchar(Nc.bmj)
      Nc.bmj[nchar.Nc.bmj > 0] <-
        str_pad(Nc.bmj[nchar.Nc.bmj > 0], width = max(nchar.Nc.bmj),
                side = "right")
    }
    #
    mean.sd.n.e.format <-
      ifelse(Me.format != "",
             paste0(Me.format, " (", Se.format, ")", bmj.sep, Ne.bmj),
             "")
    mean.sd.n.c.format <-
      ifelse(Mc.format != "",
             paste0(Mc.format, " (", Sc.format, ")", bmj.sep, Nc.bmj),
             "")
    #
    sel.e <- Ne.format != "" & mean.sd.n.e.format == ""
    mean.sd.n.e.format[sel.e] <- Ne.format[sel.e]
    #
    sel.c <- Nc.format != "" & mean.sd.n.c.format == ""
    mean.sd.n.c.format[sel.c] <- Nc.format[sel.c]
  }
  #
  if (ev.n.prop) {
    Ee.bmj <- Ee.format
    #
    Ne.bmj <- Ne.format
    #
    Ne.bmj <- ifelse(Ee.bmj == "", "", Ne.bmj)
    #
    while(all(substring(Ee.bmj[Ee.bmj != ""], 1, 1) == " "))
      Ee.bmj[Ee.bmj != ""] <- substring(Ee.bmj[Ee.bmj != ""], 2)
    #
    while(all(substring(Ne.bmj[Ne.bmj != ""], 1, 1) == " "))
      Ne.bmj[Ne.bmj != ""] <- substring(Ne.bmj[Ne.bmj != ""], 2)
    #
    Ne.bmj <- rmSpace(Ne.bmj, end = FALSE)
    if (monospaced) {
      nchar.Ne.bmj <- nchar(Ne.bmj)
      Ne.bmj[nchar.Ne.bmj > 0] <-
        str_pad(Ne.bmj[nchar.Ne.bmj > 0], width = max(nchar.Ne.bmj),
                side = "right")
    }
    #
    event.n.e.format <-
      ifelse(Ee.bmj != "",
             paste0(Ee.bmj, bmj.sep, Ne.bmj),
             "")
    #
    if (pooled.totals & !pooled.events) {
      event.n.e.format[Ne.format != "" & Ee.format == ""] <-
        Ne.format[Ne.format != "" & Ee.format == ""]
    }
  }
  #
  # Correlation
  #
  if (by)
    cor <- c(NAs.all, x$cor)
  else
    cor <- c(NAs, x$cor)
  #
  if (is.null(digits.cor))
    cor.format <- formatN(cor, digits = digits.R, text.NA = lab.NA)
  else
    cor.format <- formatN(cor, digits.cor, lab.NA, big.mark = big.mark)
  #
  # Print nothing for lines with summary results
  #
  cor.format[all.res] <- ""
  #
  if (by)
    cor.format[sel.w] <- ""
  #
  # P-value of effect
  #
  
  pval.format <-
    formatPT(c(x$pval.common, x$pval.random, NAs.prd,
               if (is.numeric(x$pval)) x$pval else rep(NA, length(x$pval))),
             digits = digits.pval,
             big.mark = big.mark,
             lab = FALSE, labval = "",
             zero = zero.pval, JAMA = JAMA.pval,
             scientific = scientific.pval,
             lab.NA = lab.NA)
  pval.format[all.prd] <- ""
  #
  # Cluster variable
  #
  if (by)
    cluster.format <- c(NAs.all, x$cluster)
  else
    cluster.format <- c(NAs, x$cluster)
  #
  # Print nothing for lines with summary results
  #
  cluster.format[all.res] <- ""
  #
  if (by)
    cluster.format[sel.w] <- ""
  #
  # N-of-1 variable
  #
  if (by)
    cycles.format <- c(NAs.all, x$cycles)
  else
    cycles.format <- c(NAs, x$cycles)
  #
  # Print nothing for lines with summary results
  #
  cycles.format[all.res] <- ""
  #
  if (by)
    cycles.format[sel.w] <- ""
  #
  #
  # y-axis:
  #
  #
  if ((!(metaprop | metacor) &
       (any(rightcols %in% c("n.e", "n.c")) |
        any(leftcols  %in% c("n.e", "n.c")))) |
      (metainc &
       (any(rightcols %in% c("time.e", "time.c")) |
        any(leftcols  %in% c("time.e", "time.c")))
      ) |
      (metacont &
       (any(rightcols %in% c("sd.e", "sd.c")) |
        any(leftcols  %in% c("sd.e", "sd.c")))
      ) |
      (metamean &
       (any(rightcols %in% c("sd.e")) |
        any(leftcols  %in% c("sd.e")))
      ) |
      (!is.null(label.e.attach) & !is.null(label.e)) |
      (!is.null(label.c.attach) & !is.null(label.c)) |
      RoB.available |
      newline.all
      ) {
    yHead <- 2
    yHeadadd <- 1
  }
  else {
    yHead <- 1
    yHeadadd <- NA
  }
  #
  if (!by) {
    N <- n.stud
    if (study.results)
      yTE <- 1:N
    else
      yTE <- rep(NA, N)
  }
  else {
    #
    j <- 1
    k <- 0
    yBylab <- NAs.by
    yTE <- rep(NA, n.stud)
    yTE.common.w <- rep(yBylab, n.com)
    yTE.random.w <- rep(yBylab, n.ran)
    yTE.predict.w <- rep(yBylab, n.prd)
    yTE.hetstat.w <- yBylab
    yTE.effect.common.w <- yBylab
    yTE.effect.random.w <- yBylab
    #
    seq.com.w <- seq_len(n.com)
    seq.ran.w <- seq_len(n.ran)
    seq.prd.w <- seq_len(n.prd)
    #
    n.by.i <- 0
    #
    for (i in seq_len(n.by)) {
      #
      if (allstudies)
        k.i <- k.all.w[i]
      else
        k.i <- k.TE.w[i]
      #
      k <- k + k.i
      #
      if (print.subgroup.labels) {
        yBylab[i] <- j
        j <- j + 1
      }
      #
      if (study.results) {
        yTE[(k - k.i + 1):k] <- j:(j + k.i - 1)
        j <- j + k.i
      }
      else
        yTE[(k - k.i + 1):k] <- NA
      #
      # Common effect model
      #
      if (common.subgroup.logical[i] & subgroup.logical[i]) {
        yTE.common.w[n.by.i * n.com + seq.com.w] <- j + seq.com.w - 1
        j <- j + max(seq.com.w)
      }
      #
      # Random effects model
      #
      if (random.subgroup.logical[i] & subgroup.logical[i]) {
        yTE.random.w[n.by.i * n.ran + seq.ran.w] <- j + seq.ran.w - 1
        j <- j + max(seq.ran.w)
      }
      #
      # Only pooled totals
      #
      # cat("pooled.totals\n")
      # print(pooled.totals)
      # cat("subgroup.logical[i]\n")
      # print(subgroup.logical[i])
      # cat("common.subgroup.logical[i]\n")
      # print(common.subgroup.logical[i])
      # cat("random.subgroup.logical[i]\n")
      # print(random.subgroup.logical[i])
      #
      if (pooled.totals & subgroup.logical[i] &
          !(common.subgroup.logical[i] | random.subgroup.logical[i])) {
        yTE.common.w[n.by.i * n.com + seq.com.w] <- j + seq.com.w - 1
        j <- j + max(seq.com.w)
      }
      #
      # Prediction interval in subgroups
      #
      if (prediction.subgroup.logical[i] & subgroup.logical[i]) {
        yTE.predict.w[n.by.i * n.prd + seq.prd.w] <-
          j + seq.prd.w - 1
        j <- j + max(seq.prd.w)
      }
      #
      # Heterogeneity statistics
      #
      if (subgroup.hetstat.logical[i]) {
        yTE.hetstat.w[i] <- j
        j <- j + 1
      }
      else
        yTE.hetstat.w[i] <- NA
      #
      # Test for effect in subgroup (common effect)
      #
      if (test.effect.subgroup.common.logical[i]) {
        yTE.effect.common.w[i] <- j
        j <- j + 1
      }
      else
        yTE.effect.common.w[i] <- NA
      #
      # Test for effect in subgroup (random effects)
      #
      if (test.effect.subgroup.random.logical[i]) {
        yTE.effect.random.w[i] <- j
        j <- j + 1
      }
      else
        yTE.effect.random.w[i] <- NA
      #
      y.w.i <- c(yTE.common.w[i], yTE.random.w[i], yTE.predict.w[i],
                 yTE.hetstat.w[i],
                 yTE.effect.common.w[i],
                 yTE.effect.random.w[i])
      #
      if (!study.results & !any(!is.na(y.w.i))) {
        yBylab[i] <- NA
        if (print.subgroup.labels)
          j <- j - 1
      }
      #
      if (!is.na(yBylab[i]) & addrow.subgroups)
        j <- j + 1
      #
      n.by.i <- n.by.i + 1
    }
    #
    if (!addrow.subgroups)
      j <- j + 1
    #
    yTE.w <- c(yTE.common.w, yTE.random.w, yTE.predict.w, yTE.hetstat.w,
               yTE.effect.common.w, yTE.effect.random.w)
  }
  #
  #
  # x-axis:
  #
  #
  if (avail.xlim && is.numeric(xlim[1]))
    if (log.xaxis)
      xlim <- log(xlim)
  #
  if (is.null(xlim)) {
    sel.low <- is.finite(lowTEs)
    sel.upp <- is.finite(uppTEs)
    #
    if (all(!sel.low))
      minTE <- -0.5
    else
      minTE <- min(lowTEs[sel.low], na.rm = TRUE)
    if (all(!sel.upp))
      maxTE <- 0.5
    else
      maxTE <- max(uppTEs[sel.upp], na.rm = TRUE)
    #
    xlim <- c(minTE, maxTE)
    #
    if (!is.na(ref) && ref < xlim[1])
      xlim[1] <- ref
    if (!is.na(ref) && ref > xlim[2])
      xlim[2] <- ref
    #
    if (!is.na(cid.below.null) && cid.below.null < xlim[1])
      xlim[1] <- cid.below.null
    if (!is.na(cid.below.null) && cid.below.null > xlim[2])
      xlim[2] <- cid.below.null
    #
    if (!is.na(cid.above.null) && cid.above.null < xlim[1])
      xlim[1] <- cid.above.null
    if (!is.na(cid.above.null) && cid.above.null > xlim[2])
      xlim[2] <- cid.above.null
  }
  #
  symmetric <- FALSE
  #
  if (!is.null(xlim) && is.character(xlim[1])) {
    #
    xlim <- setchar(xlim, "symmetric",
                    paste0("should be a numeric vector (min, max) or ",
                           "the character string \"symmetric\""))
    symmetric <- TRUE
    #
    sel.low <- is.finite(lowTEs)
    sel.upp <- is.finite(uppTEs)
    #
    if (all(!sel.low))
      minTE <- -0.5
    else
      minTE <- min(lowTEs[sel.low], na.rm = TRUE)
    if (all(!sel.upp))
      maxTE <- 0.5
    else
      maxTE <- max(uppTEs[sel.upp], na.rm = TRUE)
    #
    if (minTE < 0 & maxTE < 0)
      xlim <- c(minTE, -minTE)
    else if (minTE > 0 & maxTE > 0)
      xlim <- c(-maxTE, maxTE)
    else
      xlim <- c(-max(abs(c(minTE, maxTE))), max(abs(c(minTE, maxTE))))
  }
  #
  if (!is.na(ref) &&
      round(xlim[2] - ref, 6) == round(ref - xlim[1], 6))
    symmetric <- TRUE
  #
  if (by) {
    if (all(is.na(c(yTE, yTE.w))))
      max.yTE <- 0
    else
      max.yTE <- max(c(yTE, yTE.w), na.rm = TRUE)
  }
  else {
    if (all(is.na(yTE)))
      max.yTE <- 0
    else
      max.yTE <- max(yTE, na.rm = TRUE)
  }
  #
  yNext <- max.yTE + ifelse(max.yTE == 0 | !addrow.overall, 1, 2)
  #
  if (missing.xlab.pos)
    xlab.pos <- mean(xlim)
  #
  if (missing.smlab.pos)
    smlab.pos <- mean(xlim)
  #
  yTE.common  <- rep(NA, n.com)
  yTE.random <- rep(NA, n.ran)
  yPredict <- rep(NA, n.prd)
  yHetstat <- NA
  yResidHetstat <- NA
  yOverall.common  <- NA
  yOverall.random <- NA
  ySubgroup.common  <- NA
  ySubgroup.random <- NA
  yText.addline1 <- NA
  yText.addline2 <- NA
  #
  if (details)
    yText.details <- rep_len(NA, length(text.details))
  else
    yText.details <- NULL
  #
  if (RoB.legend)
    yText.rob <- rep_len(NA, length(text.rob))
  else
    yText.rob <- NULL
  #
  seq.com <- seq_len(n.com) - 1
  seq.ran <- seq_len(n.ran) - 1
  seq.prd <- seq_len(n.prd) - 1
  #
  if (common & random & overall) {
    yTE.common <- yNext + seq.com
    yNext <- yNext + 1 + max(seq.com)
    #
    yTE.random <- yNext + seq.ran
    yNext <- yNext + 1 + max(seq.ran)
  }
  #
  else if (common & !random & overall) {
    yTE.common <- yNext + seq.com
    yNext <- yNext + 1 + max(seq.com)
  }
  #
  else if (!common & random & overall) {
    yTE.random <- yNext + seq.ran
    yNext <- yNext + 1 + max(seq.ran)
  }
  #
  else if (!common & !random & pooled.totals & overall) {
    yTE.common <- yNext + seq.com
    yNext <- yNext + 1 + max(seq.com)
    #
    if (missing.text.common)
      text.common <- "Overall"
  }
  #
  if (prediction & overall) {
    yPredict <- yNext + seq.prd
    yNext <- yNext + 1 + max(seq.prd)
  }
  #
  yNext <- yNext + addrows.below.overall
  #
  if (overall.hetstat) {
    yHetstat <- yNext
    yNext <- yNext + 1
  }
  #
  if (resid.hetstat) {
    yResidHetstat <- yNext
    yNext <- yNext + 1
  }
  #
  if (test.overall.common) {
    yOverall.common <- yNext
    yNext <- yNext + 1
  }
  #
  if (test.overall.random) {
    yOverall.random <- yNext
    yNext <- yNext + 1
  }
  #
  if (test.subgroup.common) {
    ySubgroup.common <- yNext
    yNext <- yNext + 1
  }
  #
  if (test.subgroup.random) {
    ySubgroup.random <- yNext
    yNext <- yNext + 1
  }
  #
  if (!missing.text.addline1) {
    yText.addline1 <- yNext
    yNext <- yNext + 1
  }
  #
  if (!missing.text.addline2) {
    yText.addline2 <- yNext
    yNext <- yNext + 1
  }
  #
  if (details) {
    yNext <- yNext + 1
    yText.details <- seq(yNext, yNext + length(yText.details) - 1)
    yNext <- yNext + length(yText.details)
  }
  #
  if (RoB.legend) {
    yNext <- yNext + 1
    yText.rob <- seq(yNext, yNext + length(yText.rob) - 1)
  }
  #
  if (!common & !pooled.totals) text.common <- ""
  if (!random) text.random <- ""
  if (!prediction) text.predict <- ""
  #
  yTE <- yHead + yTE + addrow
  #
  yTE.common <- yHead + yTE.common + addrow
  yTE.random <- yHead + yTE.random + addrow
  yPredict   <- yHead + yPredict   + addrow
  #
  yHetstat <- yHead + yHetstat + addrow
  yResidHetstat <- yHead + yResidHetstat + addrow
  yOverall.common <- yHead + yOverall.common + addrow
  yOverall.random <- yHead + yOverall.random + addrow
  ySubgroup.common <- yHead + ySubgroup.common + addrow
  ySubgroup.random <- yHead + ySubgroup.random + addrow
  yText.addline1 <- yHead + yText.addline1 + addrow
  yText.addline2 <- yHead + yText.addline2 + addrow
  if (details)
    yText.details <- yHead + yText.details + addrow
  if (RoB.legend)
    yText.rob <- yHead + yText.rob + addrow
  #
  yStats <- c(yHetstat,
              yResidHetstat,
              yOverall.common, yOverall.random,
              ySubgroup.common, ySubgroup.random,
              yText.addline1, yText.addline2)
  #
  yStatsDetails <- c(yHetstat,
                     yResidHetstat,
                     yOverall.common, yOverall.random,
                     ySubgroup.common, ySubgroup.random,
                     yText.addline1, yText.addline2,
                     yText.details, yText.rob)
  #
  if (by) {
    yBylab <- yHead + yBylab + addrow
    yTE.w  <- yHead + yTE.w + addrow
    #
    yLab <- c(yHead,
              yTE.common, yTE.random, yPredict,
              yStatsDetails,
              yBylab, yTE.w,
              yTE)
    #
    yS <- c(yHead, yTE.common, yTE.random, yPredict, yTE.w, yTE)
    yS.noma <- c(yHead,
                 rep_len(NA, length(yTE.common)),
                 rep_len(NA, length(yTE.random)),
                 rep_len(NA, length(yPredict)),
                 rep_len(NA, length(yTE.w)),
                 yTE)
  }
  else {
    yLab <- c(yHead, yTE.common, yTE.random, yPredict,
              yStatsDetails,
              yTE)
    #
    yS <- c(yHead, yTE.common, yTE.random, yPredict, yTE)
    yS.noma <- c(yHead,
                 rep_len(NA, length(yTE.common)),
                 rep_len(NA, length(yTE.random)),
                 rep_len(NA, length(yPredict)),
                 yTE)
  }
  
  
  #
  #
  # (11) Format columns in forest plot
  #
  #
  col.studlab <- list(labels =
                        lapply(as.list(c(labs[["lab.studlab"]], modlabs)),
                               tg,
                               xpos = xpos.s, just = just.s,
                               fs = fs.study.labels,
                               ff = ff.study.labels,
                               fontfamily = fontfamily),
                      rows = yLab
                      )
  # Study label:
  col.studlab$labels[[1]] <- tg(labs[["lab.studlab"]], xpos.s,
                                just.s, fs.head, ff.head, fontfamily)
  # Common effect estimate:
  strt <- 1
  for (i in seq_len(n.com)) {
    col.studlab$labels[[strt + i]] <-
      tg(text.common[i], xpos.s, just.s,
         fs.common.labels, ff.common.labels, fontfamily)
  }
  # Random effects estimate:
  strt <- strt + i
  for (i in seq_len(n.ran)) {
    col.studlab$labels[[strt + i]] <-
      tg(text.random[i], xpos.s, just.s,
         fs.random.labels, ff.random.labels, fontfamily)
  }
  # Prediction interval:
  strt <- strt + i
  for (i in seq_len(n.prd)) {
    col.studlab$labels[[strt + i]] <-
      tg(text.predict[i], xpos.s, just.s,
         fs.predict.labels, ff.predict.labels, fontfamily)
  }
  # Heterogeneity statistics:
  strt <- strt + i
  col.studlab$labels[[strt + 1]] <-
    tg(hetstat.overall, xpos.s, just.s, fs.hetstat, ff.hetstat, fontfamily)
  # Statistic for residual heterogeneity:
  col.studlab$labels[[strt + 2]] <-
    tg(hetstat.resid, xpos.s, just.s, fs.hetstat, ff.hetstat, fontfamily)
  # Test for overall effect (common effect model):
  col.studlab$labels[[strt + 3]] <-
    tg(text.overall.common, xpos.s, just.s, fs.test.overall, ff.test.overall,
       fontfamily)
  # Test for overall effect (random effects model):
  col.studlab$labels[[strt + 4]] <-
    tg(text.overall.random, xpos.s, just.s, fs.test.overall, ff.test.overall,
       fontfamily)
  # Test for subgroup differences (common effect model):
  col.studlab$labels[[strt + 5]] <-
    tg(text.subgroup.common, xpos.s, just.s, fs.test.subgroup, ff.test.subgroup,
       fontfamily)
  # Test for subgroup differences (random effects model):
  col.studlab$labels[[strt + 6]] <-
    tg(text.subgroup.random, xpos.s, just.s, fs.test.subgroup, ff.test.subgroup,
       fontfamily)
  # First additional line:
  col.studlab$labels[[strt + 7]] <-
    tg(text.addline1, xpos.s, just.s, fs.addline, ff.addline, fontfamily)
  # Second additional line:
  col.studlab$labels[[strt + 8]] <-
    tg(text.addline2, xpos.s, just.s, fs.addline, ff.addline, fontfamily)
  # Details on meta-analysis methods:
  strt <- strt + 8
  i <- 0
  if (details) {
    for (i in seq_len(length(text.details))) {
      col.studlab$labels[[strt + i]] <-
        tg(text.details[[i]], xpos.s, just.s,
           if (i == 1) fs.heading else fs.details,
           if (i == 1) ff.heading else ff.details,
           fontfamily)
    }
  }
  # Risk of bias:
  strt <- strt + i
  i <- 0
  if (RoB.legend) {
    for (i in seq_len(length(text.rob))) {
      col.studlab$labels[[strt + i]] <-
        tg(text.rob[i], xpos.s, just.s,
           if (i == 1) fs.heading else fs.rob,
           if (i == 1) ff.heading else ff.rob,
           fontfamily)
    }
  }
  #
  if (by) {
    # Subgroup labels:
    strt <- strt + i
    for (i in seq_len(n.by)) {
      col.studlab$labels[[strt + i]] <-
        tg(subgroup.name[i], xpos.s, just.s,
           fs.head, ff.head, fontfamily, col.subgroup)
    }
    # Common effect estimates:
    strt <- strt + i
    for (i in seq_len(n.com.w)) {
      col.studlab$labels[[strt + i]] <-
        tg(text.common.w[i], xpos.s, just.s,
           fs.common.labels, ff.common.labels, fontfamily, col.subgroup)
    }
    # Random effects estimates:
    strt <- strt + i
    for (i in seq_len(n.ran.w)) {
      col.studlab$labels[[strt + i]] <-
        tg(text.random.w[i], xpos.s, just.s,
           fs.random.labels, ff.random.labels, fontfamily, col.subgroup)
    }
    # Prediction interval:
    strt <- strt + i
    for (i in seq_len(n.prd.w)) {
      col.studlab$labels[[strt + i]] <-
        tg(text.predict.w[i], xpos.s, just.s,
           fs.predict.labels, ff.predict.labels, fontfamily, col.subgroup)
    }
    # Heterogeneity statistics:
    strt <- strt + i
    for (i in seq_len(n.by)) {
      col.studlab$labels[[strt + i]] <-
        tg(hetstat.w[[i]], xpos.s, just.s,
           fs.hetstat, ff.hetstat, fontfamily, col.subgroup)
    }
    # Test for effect in subgroup (common effect model):
    strt <- strt + i
    for (i in seq_len(n.by)) {
      col.studlab$labels[[strt + i]] <-
        tg(text.effect.subgroup.common[[i]], xpos.s, just.s,
           fs.test.effect.subgroup, ff.test.effect.subgroup,
           fontfamily, col.subgroup)
    }
    # Test for effect in subgroup (random effects model):
    strt <- strt + i
    for (i in seq_len(n.by)) {
      col.studlab$labels[[strt + i]] <-
        tg(text.effect.subgroup.random[[i]], xpos.s, just.s,
           fs.test.effect.subgroup, ff.test.effect.subgroup,
           fontfamily, col.subgroup)
    }
  }
  #
  fcs <- list(fs.study = fs.study, ff.study = ff.study,
              fs.heading = fs.head, ff.heading = ff.head,
              fs.common = fs.common, ff.common = ff.common,
              fs.random = fs.random, ff.random = ff.random,
              fs.predict = fs.predict, ff.predict = ff.predict,
              by = by, n.by = n.by, col.subgroup = col.subgroup)
  #
  col.effect <- formatcol(labs[["lab.effect"]], effect.format,
                          yS, just.c, fcs, fontfamily,
                          n.com, n.ran, n.prd)
  #
  col.ci <- formatcol(labs[["lab.ci"]], ci.format, yS, just.c, fcs, fontfamily,
                      n.com, n.ran, n.prd)
  #
  col.effect.ci <-
    formatcol(labs[["lab.effect.ci"]], effect.ci.format, yS,
              if (bmj.revman5 & missing.just) "center" else just.c,
              fcs, fontfamily,
              n.com, n.ran, n.prd)
  #
  if (ev.n.bin) {
    col.event.n.e <-
      formatcol(labs[["lab.event.n.e"]], event.n.e.format, yS,
                if (bmj.revman5 & missing.just) "center" else just.c,
                fcs, fontfamily,
                n.com, n.ran, n.prd)
    #
    col.event.n.c <-
      formatcol(labs[["lab.event.n.c"]], event.n.c.format, yS,
                if (bmj.revman5 & missing.just) "center" else just.c,
                fcs, fontfamily,
                n.com, n.ran, n.prd)
  }
  #
  if (m.s.n.cont) {
    col.mean.sd.n.e <-
      formatcol(labs[["lab.mean.sd.n.e"]], mean.sd.n.e.format, yS,
                if (bmj.revman5 & missing.just) "center" else just.c,
                fcs, fontfamily,
                n.com, n.ran, n.prd)
    #
    col.mean.sd.n.c <-
      formatcol(labs[["lab.mean.sd.n.c"]], mean.sd.n.c.format, yS,
                if (bmj.revman5 & missing.just) "center" else just.c,
                fcs, fontfamily,
                n.com, n.ran, n.prd)
  }
  #
  if (ev.n.prop) {
    col.event.n.e <-
      formatcol(labs[["lab.event.n.e"]], event.n.e.format, yS,
                if (bmj.revman5 & missing.just) "center" else just.c,
                fcs, fontfamily,
                n.com, n.ran, n.prd)
  }
  #
  col.w.common  <- formatcol(labs[["lab.w.common"]], Wc.format, yS,
                             just.c, fcs, fontfamily,
                             n.com, n.ran, n.prd)
  col.w.random <- formatcol(labs[["lab.w.random"]], Wr.format, yS,
                            just.c, fcs, fontfamily,
                            n.com, n.ran, n.prd)
  #
  col.TE <- formatcol(labs[["lab.TE"]], TE.format, yS,
                      just.c, fcs, fontfamily,
                      n.com, n.ran, n.prd)
  col.seTE <- formatcol(labs[["lab.seTE"]], seTE.format, yS,
                        just.c, fcs, fontfamily,
                        n.com, n.ran, n.prd)
  #
  col.cluster <-
    formatcol(labs[["lab.cluster"]], cluster.format, yS,
              if (as.character.cluster) "left" else just.c, fcs, fontfamily,
              n.com, n.ran, n.prd)
  #
  col.cycles <-
    formatcol(labs[["lab.cycles"]], cycles.format, yS,
              if (as.character.cycles) "left" else just.c, fcs, fontfamily,
              n.com, n.ran, n.prd)
  #
  col.n.e <- formatcol(labs[["lab.n.e"]], Ne.format, yS, just.c, fcs,
                       fontfamily,
                       n.com, n.ran, n.prd)
  col.n.c <- formatcol(labs[["lab.n.c"]], Nc.format, yS, just.c, fcs,
                       fontfamily,
                       n.com, n.ran, n.prd)
  #
  col.event.e <- formatcol(labs[["lab.event.e"]], Ee.format, yS, just.c, fcs,
                           fontfamily,
                           n.com, n.ran, n.prd)
  col.event.c <- formatcol(labs[["lab.event.c"]], Ec.format, yS, just.c, fcs,
                           fontfamily,
                           n.com, n.ran, n.prd)
  #
  col.mean.e <- formatcol(labs[["lab.mean.e"]], Me.format, yS, just.c, fcs,
                          fontfamily,
                          n.com, n.ran, n.prd)
  col.mean.c <- formatcol(labs[["lab.mean.c"]], Mc.format, yS, just.c, fcs,
                          fontfamily,
                          n.com, n.ran, n.prd)
  #
  col.sd.e <- formatcol(labs[["lab.sd.e"]], Se.format, yS, just.c, fcs,
                        fontfamily,
                        n.com, n.ran, n.prd)
  col.sd.c <- formatcol(labs[["lab.sd.c"]], Sc.format, yS, just.c, fcs,
                        fontfamily,
                        n.com, n.ran, n.prd)
  #
  col.cor <- formatcol(labs[["lab.cor"]], cor.format, yS, just.c, fcs,
                       fontfamily,
                       n.com, n.ran, n.prd)
  #
  col.time.e <- formatcol(labs[["lab.time.e"]], Te.format, yS, just.c, fcs,
                          fontfamily,
                          n.com, n.ran, n.prd)
  col.time.c <- formatcol(labs[["lab.time.c"]], Tc.format, yS, just.c, fcs,
                          fontfamily,
                          n.com, n.ran, n.prd)
  #
  col.pval <- formatcol(labs[["lab.pval"]], pval.format, yS, just.c, fcs,
                        fontfamily,
                        n.com, n.ran, n.prd)
  #
  col.effect.calc <- formatcol(longer.effect, effect.format, yS, just.c, fcs,
                               fontfamily,
                               n.com, n.ran, n.prd)
  #
  col.ci.calc <- formatcol(longer.ci, ci.format, yS, just.c, fcs,
                           fontfamily,
                           n.com, n.ran, n.prd)
  #
  col.effect.ci.calc <- formatcol(longer.effect.ci, effect.ci.format, yS,
                                  just.c, fcs, fontfamily,
                                  n.com, n.ran, n.prd)
  #
  if (ev.n.bin) {
    col.event.n.e.calc <- formatcol(longer.event.n.e, event.n.e.format, yS,
                                    just.c, fcs, fontfamily,
                                    n.com, n.ran, n.prd)
    #
    col.event.n.c.calc <- formatcol(longer.event.n.c, event.n.c.format, yS,
                                    just.c, fcs, fontfamily,
                                    n.com, n.ran, n.prd)
  }
  #
  if (m.s.n.cont) {
    col.mean.sd.n.e.calc <-
      formatcol(longer.mean.sd.n.e, mean.sd.n.e.format, yS,
                just.c, fcs, fontfamily,
                n.com, n.ran, n.prd)
    #
    col.mean.sd.n.c.calc <-
      formatcol(longer.mean.sd.n.c, mean.sd.n.c.format, yS,
                just.c, fcs, fontfamily,
                n.com, n.ran, n.prd)
  }
  #
  if (ev.n.prop) {
    col.event.n.e.calc <- formatcol(longer.event.n.e, event.n.e.format, yS,
                                    just.c, fcs, fontfamily,
                                    n.com, n.ran, n.prd)
  }
  #
  col.w.common.calc  <- formatcol(longer.w.common, Wc.format, yS,
                                  just.c, fcs, fontfamily,
                                  n.com, n.ran, n.prd)
  col.w.random.calc <- formatcol(longer.w.random, Wr.format, yS,
                                 just.c, fcs, fontfamily,
                                 n.com, n.ran, n.prd)
  #
  col.TE.calc <- formatcol(longer.TE, TE.format, yS, just.c, fcs, fontfamily,
                           n.com, n.ran, n.prd)
  col.seTE.calc <- formatcol(longer.seTE, seTE.format, yS, just.c, fcs,
                             fontfamily,
                             n.com, n.ran, n.prd)
  #
  col.cluster.calc <- formatcol(longer.cluster, cluster.format,
                                yS, just.c, fcs, fontfamily,
                                n.com, n.ran, n.prd)
  #
  col.cycles.calc <- formatcol(longer.cycles, cycles.format,
                               yS, just.c, fcs, fontfamily,
                               n.com, n.ran, n.prd)
  #
  col.n.e.calc <- formatcol(longer.n.e, Ne.format, yS, just.c, fcs, fontfamily,
                            n.com, n.ran, n.prd)
  col.n.c.calc <- formatcol(longer.n.c, Nc.format, yS, just.c, fcs, fontfamily,
                            n.com, n.ran, n.prd)
  #
  col.event.e.calc <- formatcol(longer.event.e, Ee.format, yS,
                                just.c, fcs, fontfamily,
                                n.com, n.ran, n.prd)
  col.event.c.calc <- formatcol(longer.event.c, Ec.format, yS, just.c, fcs,
                                fontfamily,
                                n.com, n.ran, n.prd)
  #
  col.mean.e.calc <- formatcol(longer.mean.e, Me.format, yS, just.c, fcs,
                               fontfamily,
                               n.com, n.ran, n.prd)
  col.mean.c.calc <- formatcol(longer.mean.c, Mc.format, yS, just.c, fcs,
                               fontfamily,
                               n.com, n.ran, n.prd)
  #
  col.sd.e.calc <- formatcol(longer.sd.e, Se.format, yS, just.c, fcs,
                             fontfamily,
                             n.com, n.ran, n.prd)
  col.sd.c.calc <- formatcol(longer.sd.c, Sc.format, yS, just.c, fcs,
                             fontfamily,
                             n.com, n.ran, n.prd)
  #
  col.cor.calc <- formatcol(longer.cor, cor.format, yS, just.c, fcs,
                            fontfamily,
                            n.com, n.ran, n.prd)
  #
  col.time.e.calc <- formatcol(longer.time.e, Te.format, yS, just.c, fcs,
                               fontfamily,
                               n.com, n.ran, n.prd)
  col.time.c.calc <- formatcol(longer.time.c, Tc.format, yS, just.c, fcs,
                               fontfamily,
                               n.com, n.ran, n.prd)
  #
  col.pval.calc <- formatcol(longer.pval, pval.format, yS, just.c, fcs,
                             fontfamily,
                             n.com, n.ran, n.prd)
  #
  col.forest <- list(eff = TEs.exclude,
                     low = lowTEs.exclude,
                     upp = uppTEs.exclude,
                     rows = yS[-1],
                     #
                     # "square"  - normal CI with square
                     # "circle"  - normal CI with circle
                     # "squarediamond" - normal CI with squarediamond
                     # "diamond" - meta-analysis diamond
                     # "predict" - prediction interval
                     #
                     type = c(type.pooled, type.study),
                     #
                     col = c(col.diamond.lines.pooled, col.study),
                     col.square = c(col.diamond.pooled, col.square),
                     col.square.lines =
                       c(col.diamond.lines.pooled, col.square.lines),
                     col.inside = c(col.inside.pooled, col.inside),
                     col.circle = c(col.diamond.pooled, col.circle),
                     col.circle.lines =
                       c(col.diamond.lines.pooled, col.circle.lines),
                     #
                     col.diamond = c(col.diamond.pooled, col.square),
                     col.diamond.lines =
                       c(col.diamond.lines.pooled, col.square.lines),
                     #
                     lwd = lwd,
                     lwd.square = lwd.square, lwd.diamond = lwd.diamond,
                     #
                     arrow.type = arrow.type, arrow.length = arrow.length
                     )
  #
  # Sizes of squares
  #
  if (weight.study == "same") {
    information <- rep(0.9, length(TEs))
  }
  else {
    #
    if (weight.study == "common")
      information <- sqrt(w.commons)
    else if (weight.study == "random")
      information <- sqrt(w.randoms)
    # Square height equal to 1 for most precise study result
    if (!all(is.na(information))) {
      information <- information / max(information, na.rm = TRUE)
      if (length(unique(type.study)) == 1 &&
          unique(type.study) == "squarediamond")
        information[information < 0.35 * max(information, na.rm = TRUE)] <-
          0.35
    }
    else
      information <- rep(0.9, length(TEs))
    # Same / maximum polygon height for all meta-analytical results
    # (both overall and subgroup results)
    information[is.na(information)] <- 1
  }
  #
  col.forest$sizes <- information
  col.forest$sizes <- col.forest$sizes * squaresize
  #
  # Width of column 3
  col.forestwidth <- plotwidth
  #
  # Range on the x-axis for column 3
  col.forest$range <- xlim
  #
  cols <- list(col.studlab = col.studlab,
               col.effect = col.effect,
               col.ci = col.ci,
               col.effect.ci = col.effect.ci,
               col.w.common = col.w.common,
               col.w.random = col.w.random,
               col.TE = col.TE,
               col.seTE = col.seTE)
  #
  cols.new <- vector("list")
  #
  if (ev.n.bin) {
    cols$col.event.n.e <- col.event.n.e
    cols$col.event.n.c <- col.event.n.c
  }
  #
  if (m.s.n.cont) {
    cols$col.mean.sd.n.e <- col.mean.sd.n.e
    cols$col.mean.sd.n.c <- col.mean.sd.n.c
  }
  #
  if (ev.n.prop) {
    cols$col.event.n.e <- col.event.n.e
  }
  #
  cols.calc <- list(col.studlab = col.studlab,
                    col.effect = col.effect.calc,
                    col.ci = col.ci.calc,
                    col.effect.ci = col.effect.ci.calc,
                    col.w.common = col.w.common.calc,
                    col.w.random = col.w.random.calc,
                    col.TE = col.TE.calc,
                    col.seTE = col.seTE.calc)
  #
  if (ev.n.bin) {
    cols.calc$col.event.n.e <- col.event.n.e.calc
    cols.calc$col.event.n.c <- col.event.n.c.calc
  }
  #
  if (m.s.n.cont) {
    cols.calc$col.mean.sd.n.e <- col.mean.sd.n.e.calc
    cols.calc$col.mean.sd.n.c <- col.mean.sd.n.c.calc
  }
  #
  if (ev.n.prop) {
    cols.calc$col.event.n.e <- col.event.n.e.calc
  }
  #
  cols[["col.cluster"]] <- col.cluster
  #
  cols[["col.cycles"]] <- col.cycles
  #
  cols[["col.n.e"]] <- col.n.e
  cols[["col.n.c"]] <- col.n.c
  cols[["col.event.e"]] <- col.event.e
  cols[["col.event.c"]] <- col.event.c
  #
  cols[["col.mean.e"]] <- col.mean.e
  cols[["col.mean.c"]] <- col.mean.c
  cols[["col.sd.e"]] <- col.sd.e
  cols[["col.sd.c"]] <- col.sd.c
  #
  cols[["col.cor"]] <- col.cor
  #
  cols[["col.time.e"]] <- col.time.e
  cols[["col.time.c"]] <- col.time.c
  #
  cols[["col.pval"]] <- col.pval
  #
  # Calculate
  #
  cols.calc[["col.cluster"]] <- col.cluster.calc
  #
  cols.calc[["col.cycles"]] <- col.cycles.calc
  #
  cols.calc[["col.n.e"]] <- col.n.e.calc
  cols.calc[["col.n.c"]] <- col.n.c.calc
  cols.calc[["col.event.e"]] <- col.event.e.calc
  cols.calc[["col.event.c"]] <- col.event.c.calc
  #
  cols.calc[["col.mean.e"]] <- col.mean.e.calc
  cols.calc[["col.mean.c"]] <- col.mean.c.calc
  cols.calc[["col.sd.e"]] <- col.sd.e.calc
  cols.calc[["col.sd.c"]] <- col.sd.c.calc
  #
  cols.calc[["col.cor"]] <- col.cor.calc
  #
  cols.calc[["col.time.e"]] <- col.time.e.calc
  cols.calc[["col.time.c"]] <- col.time.c.calc
  #
  cols.calc[["col.pval"]] <- col.pval.calc
  #
  if (newcols) {
    #
    # Check just.addcols
    #
    if (length(leftcols.new) > 0)
      if (length(just.addcols.left) != 1) {
        if (length(just.addcols.left) != length(leftcols.new))
          stop("Length of argument 'just.addcols.left' must be one or ",
               "same as number of additional columms in argument 'leftcols'.")
      }
      else
        just.addcols.left <- rep(just.addcols.left, length(leftcols.new))
    #
    if (length(rightcols.new) > 0)
      if (length(just.addcols.right) != 1) {
        if (length(just.addcols.right) != length(rightcols.new))
          stop("Length of argument 'just.addcols.right' must be one or ",
               "same as number of additional columms in argument 'rightcols'.")
      }
      else
        just.addcols.right <- rep(just.addcols.right, length(rightcols.new))
    #
    # Check digits.addcols
    #
    if (length(leftcols.new) > 0)
      if (length(digits.addcols.left) != 1) {
        if (length(digits.addcols.left) != length(leftcols.new))
          stop("Length of argument 'digits.addcols.left' must be one or ",
               "same as number of additional columms in argument 'leftcols'.")
      }
      else
        digits.addcols.left <- rep(digits.addcols.left, length(leftcols.new))
    #
    if (length(rightcols.new) > 0)
      if (length(digits.addcols.right) != 1) {
        if (length(digits.addcols.right) != length(rightcols.new))
          stop("Length of argument 'digits.addcols.right' must be one or ",
               "same as number of additional columms in argument 'rightcols'.")
      }
      else
        digits.addcols.right <- rep(digits.addcols.right, length(rightcols.new))
    #
    if (by) {
      for (i in seq_along(rightcols.new)) {
        col.i <-
          newCol(rightcols.new[i], rightlabs.new[i],
                 rob, dataset1, dataset2, data.pooled,
                 n.com, n.ran, n.prd,
                 notavail.digits.addcols.right, lab.NA, big.mark,
                 zero.pval, JAMA.pval, scientific.pval,
                 digits.addcols.right[i],
                 digits.pval, digits.tau2, digits.tau, digits.I2,
                 all.prd,
                 n.com.w, n.ran.w, n.prd.w, n.stat.w)
        #
        cols.new[[col.i$colname]] <- col.i$format_var
        #
        cols[[col.i$colname]] <-
          formatcol(col.i$label, col.i$format_var,
                    if (!col.i$rob) yS else yS.noma,
                    if (rightcols.new[i] %in%
                        c("pval", "tau2", "tau", "I2", "Q", "pval.Q"))
                      just
                    else
                      just.addcols.right[i],
                    fcs, fontfamily,
                    n.com, n.ran, n.prd,
                    rightcols.new[i] %in% colnames(rob))
        #
        cols.calc[[col.i$colname]] <-
          formatcol(col.i$longer, col.i$format_var,
                    yS,
                    just.addcols.right[i],
                    fcs, fontfamily,
                    n.com, n.ran, n.prd)
      }
      #
      for (i in seq_along(leftcols.new)) {
        col.i <-
          newCol(leftcols.new[i], leftlabs.new[i],
                 rob, dataset1, dataset2, data.pooled,
                 n.com, n.ran, n.prd,
                 notavail.digits.addcols.left, lab.NA, big.mark,
                 zero.pval, JAMA.pval, scientific.pval,
                 digits.addcols.left[i],
                 digits.pval, digits.tau2, digits.tau, digits.I2,
                 all.prd,
                 n.com.w, n.ran.w, n.prd.w, n.stat.w)
        #
        cols.new[[col.i$colname]] <- col.i$format_var
        #
        cols[[col.i$colname]] <-
          formatcol(col.i$label, col.i$format_var,
                    if (!col.i$rob) yS else yS.noma,
                    if (leftcols.new[i] %in%
                        c("pval", "tau2", "tau", "I2", "Q", "pval.Q"))
                      just
                    else
                      just.addcols.left[i],
                    fcs, fontfamily,
                    n.com, n.ran, n.prd,
                    leftcols.new[i] %in% colnames(rob))
        #
        cols.calc[[col.i$colname]] <-
          formatcol(col.i$longer, col.i$format_var,
                    yS,
                    just.addcols.left[i],
                    fcs, fontfamily,
                    n.com, n.ran, n.prd)
      }
    }
    else {
      for (i in seq_along(rightcols.new)) {
        col.i <-
          newCol(rightcols.new[i], rightlabs.new[i],
                 rob, dataset1, dataset2, data.pooled,
                 n.com, n.ran, n.prd,
                 notavail.digits.addcols.right, lab.NA, big.mark,
                 zero.pval, JAMA.pval, scientific.pval,
                 digits.addcols.right[i],
                 digits.pval, digits.tau2, digits.tau, digits.I2,
                 all.prd)
        #
        cols.new[[col.i$colname]] <- col.i$format_var
        #
        cols[[col.i$colname]] <-
          formatcol(col.i$label, col.i$format_var,
                    if (!col.i$rob) yS else yS.noma,
                    if (rightcols.new[i] %in%
                        c("pval", "tau2", "tau", "I2", "Q", "pval.Q"))
                      just
                    else
                      just.addcols.right[i],
                    fcs, fontfamily,
                    n.com, n.ran, n.prd,
                    rightcols.new[i] %in% colnames(rob))
        #
        cols.calc[[col.i$colname]] <-
          formatcol(col.i$longer, col.i$format_var,
                    yS,
                    just.addcols.right[i],
                    fcs, fontfamily,
                    n.com, n.ran, n.prd)
      }
      #
      for (i in seq_along(leftcols.new)) {
        col.i <-
          newCol(leftcols.new[i], leftlabs.new[i],
                 rob, dataset1, dataset2, data.pooled,
                 n.com, n.ran, n.prd,
                 notavail.digits.addcols.left, lab.NA, big.mark,
                 zero.pval, JAMA.pval, scientific.pval,
                 digits.addcols.left[i],
                 digits.pval, digits.tau2, digits.tau, digits.I2,
                 all.prd)
        #
        cols.new[[col.i$colname]] <- col.i$format_var
        #
        cols[[col.i$colname]] <-
          formatcol(col.i$label, col.i$format_var,
                    if (!col.i$rob) yS else yS.noma,
                    if (leftcols.new[i] %in%
                        c("pval", "tau2", "tau", "I2", "Q", "pval.Q"))
                      just
                    else
                      just.addcols.left[i],
                    fcs, fontfamily,
                    n.com, n.ran, n.prd,
                    leftcols.new[i] %in% colnames(rob))
        #
        cols.calc[[col.i$colname]] <-
          formatcol(col.i$longer, col.i$format_var,
                    yS,
                    just.addcols.left[i],
                    fcs, fontfamily,
                    n.com, n.ran, n.prd)
      }
    }
  }
  #
  col.label.e <-
    tgl(label.e,
        if (bmj) bmj.xpos else xpos.label.e,
        if (bmj) "left" else just.label.e,
        fs.head, ff.head, fontfamily)
  #
  col.label.c <-
    tgl(label.c,
        if (bmj) bmj.xpos else xpos.label.c,
        if (bmj) "left" else just.label.c,
        fs.head, ff.head, fontfamily)
  #
  col.rob <- tgl(rob.text, rob.xpos, "left", fs.head, ff.head, fontfamily)
  #
  if (newline.studlab)
    col.add.studlab <- tgl(add.studlab, xpos.s, just.s, fs.head, ff.head,
                           fontfamily)
  #
  if (newline.effect)
    col.add.effect <- tgl(add.effect, xpos.c, just.c, fs.head, ff.head,
                          fontfamily)
  #
  if (newline.ci)
    col.add.ci <- tgl(add.ci, xpos.c, just.c, fs.head, ff.head, fontfamily)
  #
  if (newline.effect.ci)
    col.add.effect.ci <- tgl(add.effect.ci,
                             if (bmj.revman5) 0.5 else xpos.c,
                             if (bmj.revman5) "center" else just.c,
                             fs.head, ff.head, fontfamily)
  #
  if (ev.n.bin) {
    if (newline.event.n.e)
      col.add.event.n.e <- tgl(add.event.n.e,
                               if (bmj.revman5) 0.5 else xpos.c,
                               if (bmj.revman5) "center" else just.c,
                               fs.head, ff.head, fontfamily)
    #
    if (newline.event.n.c)
      col.add.event.n.c <- tgl(add.event.n.c,
                               if (bmj.revman5) 0.5 else xpos.c,
                               if (bmj.revman5) "center" else just.c,
                               fs.head, ff.head, fontfamily)
  }
  #
  if (m.s.n.cont) {
    if (newline.mean.sd.n.e)
      col.add.mean.sd.n.e <- tgl(add.mean.sd.n.e,
                               if (bmj.revman5) 0.5 else xpos.c,
                               if (bmj.revman5) "center" else just.c,
                               fs.head, ff.head, fontfamily)
    #
    if (newline.mean.sd.n.c)
      col.add.mean.sd.n.c <- tgl(add.mean.sd.n.c,
                               if (bmj.revman5) 0.5 else xpos.c,
                               if (bmj.revman5) "center" else just.c,
                               fs.head, ff.head, fontfamily)
  }
  #
  if (ev.n.prop) {
    if (newline.event.n.e)
      col.add.event.n.e <- tgl(add.event.n.e,
                               if (bmj.revman5) 0.5 else xpos.c,
                               if (bmj.revman5) "center" else just.c,
                               fs.head, ff.head, fontfamily)
  }
  #
  if (newline.w.common)
    col.add.w.common <- tgl(add.w.common, xpos.c, just.c, fs.head, ff.head,
                            fontfamily)
  #
  if (newline.w.random)
    col.add.w.random <- tgl(add.w.random, xpos.c, just.c, fs.head, ff.head,
                            fontfamily)
  #
  if (newline.TE)
    col.add.TE <- tgl(add.TE, xpos.c, just.c, fs.head, ff.head,
                      fontfamily)
  #
  if (newline.seTE)
    col.add.seTE <- tgl(add.seTE, xpos.c, just.c, fs.head, ff.head, fontfamily)
  #
  if (newline.cluster)
    col.add.cluster <-
      tgl(add.cluster,
          if (as.character.cluster) 0 else xpos.c,
          if (as.character.cluster) "left" else just.c,
          fs.head, ff.head, fontfamily)
  #
  if (newline.cycles)
    col.add.cycles <-
    tgl(add.cycles,
        if (as.character.cycles) 0 else xpos.c,
        if (as.character.cycles) "left" else just.c,
        fs.head, ff.head, fontfamily)
  #
  if (newline.n.e)
    col.add.n.e <- tgl(add.n.e, xpos.c, just.c, fs.head, ff.head, fontfamily)
  #
  if (newline.n.c)
    col.add.n.c <- tgl(add.n.c, xpos.c, just.c, fs.head, ff.head, fontfamily)
  #
  if (newline.event.e)
    col.add.event.e <- tgl(add.event.e, xpos.c, just.c, fs.head, ff.head,
                           fontfamily)
  #
  if (newline.event.c)
    col.add.event.c <- tgl(add.event.c, xpos.c, just.c, fs.head, ff.head,
                           fontfamily)
  #
  if (newline.mean.e)
    col.add.mean.e <- tgl(add.mean.e, xpos.c, just.c, fs.head, ff.head,
                          fontfamily)
  #
  if (newline.mean.c)
    col.add.mean.c <- tgl(add.mean.c, xpos.c, just.c, fs.head, ff.head,
                          fontfamily)
  #
  if (newline.sd.e)
    col.add.sd.e <- tgl(add.sd.e, xpos.c, just.c, fs.head, ff.head, fontfamily)
  #
  if (newline.sd.c)
    col.add.sd.c <- tgl(add.sd.c, xpos.c, just.c, fs.head, ff.head, fontfamily)
  #
  if (newline.cor)
    col.add.cor <- tgl(add.cor, xpos.c, just.c, fs.head, ff.head, fontfamily)
  #
  if (newline.time.e)
    col.add.time.e <- tgl(add.time.e, xpos.c, just.c, fs.head, ff.head,
                          fontfamily)
  #
  if (newline.time.c)
    col.add.time.c <- tgl(add.time.c, xpos.c, just.c, fs.head, ff.head,
                          fontfamily)
  #
  leftcols  <- paste0("col.", leftcols)
  rightcols <- paste0("col.", rightcols)
  
  
  #
  #
  # (12) Calculate width of columns in forest plot
  #
  #
  # Exclude lines with summary measures from calculation of column
  # width for study labels
  #
  del.lines <- NULL
  #
  if (!calcwidth.common)
    del.lines <- 1 + seq_len(n.com)
  #
  if (!calcwidth.random)
    del.lines <-
      c(del.lines, 1 + n.com + seq_len(n.ran))
  #
  if (!calcwidth.predict)
    del.lines <-
      c(del.lines, 1 + n.com + n.ran + seq_len(n.prd))
  #
  if (!calcwidth.hetstat)
    del.lines <-
      c(del.lines, 1 + n.com + n.ran + n.prd + seq(2))
  #
  if (!calcwidth.tests)
    del.lines <-
      c(del.lines, 1 + n.com + n.ran + n.prd + 2 + seq(4))
  #
  if (!calcwidth.addline)
    del.lines <-
      c(del.lines, 1 + n.com + n.ran + n.prd + 2 + 4 + seq(2))
  #
  if (details)
    del.lines <-
      c(del.lines,
        1 + n.com + n.ran + n.prd + 2 + 4 + 2 + seq_along(text.details))
  #
  if (RoB.legend)
    del.lines <-
      c(del.lines,
        1 + n.com + n.ran + n.prd + 2 + 4 + 2 + length(text.details) +
        seq_along(text.rob))
  #
  nd <- 1 + n.com + n.ran + n.prd + 2 + 4 + 2 +
    details * length(text.details) +
    RoB.legend * length(text.rob)
  #
  #nd <- n.com + n.ran + n.prd + 2 + 4 + 2
  #
  if (by) {
    if (!calcwidth.subgroup)
      del.lines <-
        c(del.lines, nd + seq_len(n.by))
    #
    if (!calcwidth.common)
      del.lines <-
        c(del.lines, nd + n.by + seq_len(n.by * n.com))
    #
    if (!calcwidth.random)
      del.lines <-
        c(del.lines,
          nd + n.by + n.by * n.com + seq_len(n.by * n.ran))
    #
    if (!calcwidth.predict)
      del.lines <-
        c(del.lines,
          nd + n.by + n.by * n.com + n.by * n.ran + seq_len(n.by * n.prd))
    #
    if (!calcwidth.hetstat)
      del.lines <-
        c(del.lines,
          nd + n.by + n.by * n.com + n.by * n.ran + n.by * n.prd +
          seq_len(n.by))
    # test for effect (CE)
    if (!calcwidth.tests)
      del.lines <-
        c(del.lines,
          nd + n.by + n.by * n.com + n.by * n.ran + n.by * n.prd +
          n.by + seq_len(n.by))
    # test for effect (RE)
    if (!calcwidth.tests)
      del.lines <-
        c(del.lines,
          nd + n.by + n.by * n.com + n.by * n.ran + n.by * n.prd +
          n.by + n.by + seq_len(n.by))
  }
  #
  if (lsel) {
    for (i in seq_along(leftcols)) {
      if (i == 1) {
        if (leftcols[[i]] == "col.studlab" & !is.null(del.lines))
          x1 <- unit.c(wcalc(cols.calc[[leftcols[i]]]$labels[-del.lines]))
        else
          x1 <- unit.c(wcalc(cols.calc[[leftcols[i]]]$labels))
      }
      else {
        if (leftcols[[i]] == "col.studlab" & !is.null(del.lines))
          x1 <- unit.c(x1,
                       colgap.left,
                       wcalc(cols.calc[[leftcols[i]]]$labels[-del.lines]))
        else
          x1 <- unit.c(x1,
                       if (leftcols[[i - 1]] == "col.studlab")
                         colgap.studlab
                       else colgap.left,
                       wcalc(cols.calc[[leftcols[i]]]$labels))
      }
    }
    #
    x1 <- unit.c(x1, colgap.forest.left, col.forestwidth)
  }
  else
    x1 <- unit.c(col.forestwidth)
  #
  rob.attach <- NULL
  #
  if (RoB.available) {
    colnames(rob) <- paste0("col.", colnames(rob))
    rob.attach <- colnames(rob)[1]
  }
  #
  if (rsel | RoB.available) {
    for (i in seq_along(rightcols)) {
      colgap.right.i <- colgap.right
      #
      if (RoB.available && i > 1 &&
          rightcols[i] %in% colnames(rob) &&
          rightcols[i - 1] %in% colnames(rob)) {
        if (rightcols[i] == "col.RoB.O")
          colgap.right.i <- colgap.rob.overall
        else
          colgap.right.i <- colgap.rob
      }
      #
      x1 <- unit.c(x1,
                   if (i == 1) colgap.forest.right else colgap.right.i,
                   if (rightcols[i] %in% colnames(rob))
                     wcalc(grid.text("aa",
                                     just = "center",
                                     gp = gpar(col = "transparent",
                                               fontsize = fs.head,
                                               fontface = ff.head,
                                               fontfamily = fontfamily)))
                   else
                     wcalc(cols.calc[[rightcols[i]]]$labels)
                   )
    }
  }
  
  
  #
  #
  # (13) Process arguments smlab, label.left and label.right
  #
  #
  if (by) {
    addline <- addrow * (!any(c(overall.hetstat,
                                test.overall.common, test.overall.random,
                                resid.hetstat,
                                test.subgroup.common, test.subgroup.random)))
    #
    nrow <- max(addline + c(yTE, yTE.common, yTE.random, yPredict,
                            yStatsDetails, yTE.w), na.rm = TRUE)
  }
  else {
    addline <- addrow * (!any(c(test.overall.common, test.overall.random,
                                overall.hetstat)))
    #
    nrow <- max(addline + c(yTE, yTE.common, yTE.random, yPredict,
                            yStatsDetails), na.rm = TRUE)
  }
  #
  # Determine minimal value on y-axis for lines of common / random
  # effect estimate or reference line
  #
  n.lines <- sum(!is.na(yStatsDetails)) + details + RoB.legend
  ymin.line <- max(addrow | addrow.overall,
                   n.lines + (n.lines > 0) * addrows.below.overall)
  #
  if (!by) {
    if (overall & (common | random | prediction)) {
      if (ymin.line == addrows.below.overall & !(!addrow.overall | !addrow))
        ymin.line <- ymin.line + 1
      #
      if ((!missing.text.addline1 | !missing.text.addline2) &
          (!missing.text.addline1 + !missing.text.addline2) == n.lines &
          !(!addrow.overall | !addrow))
        ymin.line <- ymin.line + 1
      #
      if (addrow.overall & !addrow & n.lines == 0)
        ymin.line <- ymin.line - 1
    }
    #
    if (!overall & addrow & n.lines == 0)
      ymin.line <- ymin.line - 1
  }
  else {
    if (!overall & !overall.hetstat & addrow &
        ((n.lines > !missing.text.addline1 + !missing.text.addline2) |
         n.lines == 0))
      ymin.line <- ymin.line - 1
  }
  #
  if (hetstat %in% c("common", "random") &
      (!missing.text.addline1 | !missing.text.addline2))
    ymin.line <- ymin.line + 1
  #
  ymin.common <- spacing * (ymin.line + prediction * n.prd + random * n.ran + 0.5)
  ymin.random <- spacing * (ymin.line + prediction * n.prd + 0.5)
  ymin.ref    <- spacing * (ymin.line + (!(overall | overall.hetstat) & addrow))
  #
  ymax <- spacing * (nrow - ifelse(is.na(yHeadadd), 1, 2) - 1 * addrow)
  #
  if (cid.pooled.only)
    ymax.ref <- spacing * (ymin.line + prediction * n.prd + random * n.ran +
                             common * n.com)
  else
    ymax.ref <- ymax
  #
  # Position on y-axis of left and right labels (at bottom of forest plot)
  #
  y.bottom.lr <- ymin.line - 2.5 + (!(overall | overall.hetstat) & addrow)
  #
  # Position on y-axis of label below x-axis
  #
  xlab.ypos <- y.bottom.lr - 1 * (print.label & bottom.lr) -
    1 * (print.label & bottom.lr & (newline.lr | newline.ll))
  #
  # Summary label at top of forest plot
  #
  smlab1 <- tgl(smlab1, unit(smlab.pos, "native"), "center",
                fs.smlab, ff.smlab,
                fontfamily, rows = 1 + (!is.na(yHeadadd) & !newline.smlab))
  #
  if (newline.smlab)
    smlab2 <- tgl(smlab2, unit(smlab.pos, "native"),
                  "center", fs.smlab, ff.smlab, fontfamily,
                  rows = 2)
  #
  # Left and right label on x-axis:
  #
  if (!bottom.lr & !is.na(ref)) {
    row1.lr <- if (!newline & (newline.ll | newline.lr) & !addrow)
                 1
               else if (!is.na(yHeadadd) & addrow)
                 2
               else if (is.na(yHeadadd))
                 1
               else
                 2
    #
    ll1 <- tgl(ll1, unit(ref - (xlim[2] - xlim[1]) / 30, "native"),
               "right", fs.lr, ff.lr, fontfamily, col.label.left,
               rows = row1.lr)
    #
    if (newline.ll)
      ll2 <- tgl(ll2, unit(ref - (xlim[2] - xlim[1]) / 30, "native"),
                 "right", fs.lr, ff.lr, fontfamily, col.label.left,
                 rows = row1.lr + 1)
    #
    lr1 <- tgl(lr1, unit(ref + (xlim[2] - xlim[1]) / 30, "native"),
               "left", fs.lr, ff.lr, fontfamily, col.label.right,
               rows = row1.lr)
    #
    if (newline.lr)
      lr2 <- tgl(lr2, unit(ref + (xlim[2] - xlim[1]) / 30, "native"),
                 "left", fs.lr, ff.lr, fontfamily, col.label.right,
                 rows = row1.lr + 1)
  }
  
  
  #
  #
  # (14) Generate forest plot
  #
  #
  if (!is.null(file) | !is.null(func.gr)) {
    if (is.null(func.gr)) {
      if (grepl("pdf$", tolower(file)))
        func.gr <- "pdf"
      else if (grepl("ps$", tolower(file)))
        func.gr <- "postscript"
      else if (grepl("svg$", tolower(file)))
        func.gr <- "svg"
      else if (grepl("bmp$", tolower(file)))
        func.gr <- "bmp"
      else if (grepl("jpg$", tolower(file)) |
               grepl("jpeg$", tolower(file)))
        func.gr <- "jpeg"
      else if (grepl("png$", tolower(file)))
        func.gr <- "png"
      else if (grepl("tif$", tolower(file)) |
               grepl("tiff$", tolower(file)))
        func.gr <- "tiff"
      else
        stop("Argument 'file' has unknown file extension; either provide ",
             "admissible file extension\n  (\".pdf\", \".ps\", \".svg\", ",
             "\".bmp\", \".jpg\", \".png\", \"tif\") or ",
             "\n  graphics function (argument 'func.gr').")
      #
      type.gr <- func.gr
    }
    else
      type.gr <- deparse(substitute(func.gr))
    #
    figheight <- gh(type.gr, rows.gr,
                    #
                    if (metabind) length(x$TE) else n.stud,
                    lowTE.common, lowTE.random, lowTE.predict,
                    x$subgroup, subgroup.levels,
                    lower.common.w, lower.random.w, lower.predict.w,
                    #
                    if (metabind) FALSE else common,
                    if (metabind) FALSE else random,
                    if (metabind) FALSE else overall,
                    if (metabind) FALSE else prediction,
                    if (metabind) FALSE else overall.hetstat,
                    study.results,
                    #
                    spacing,
                    #
                    xlab, xlab.add, label.right, label.left, bottom.lr,
                    #
                    prediction.subgroup, subgroup.hetstat,
                    if (metabind) FALSE else test.overall.common,
                    if (metabind) FALSE else test.overall.random,
                    if (metabind) FALSE else test.subgroup.common,
                    if (metabind) FALSE else test.subgroup.random,
                    #
                    text.addline1, text.addline2,
                    text.details, text.rob,
                    #
                    addrow, addrow.overall,
                    addrow.subgroups,
                    if (metabind) 0 else addrows.below.overall,
                    #
                    c(leftcols, rightcols), labs,
                    text.w.common, text.w.random)
    #
    args.gr.all <-
      c(list(file = file,
             height = figheight$total_height,
             width = width),
        args.gr)
    runNN(func.gr, args.gr.all)
  }
  else
    figheight <- gh("no_device_defined", rows.gr,
                    #
                    n.stud,
                    lowTE.common, lowTE.random, lowTE.predict,
                    x$subgroup, subgroup.levels,
                    lower.common.w, lower.random.w, lower.predict.w,
                    #
                    common, random, overall,
                    prediction, overall.hetstat,
                    study.results,
                    #
                    spacing,
                    #
                    xlab, xlab.add, label.right, label.left, bottom.lr,
                    #
                    prediction.subgroup, subgroup.hetstat,
                    test.overall.common, test.overall.random,
                    test.subgroup.common, test.subgroup.random,
                    #
                    text.addline1, text.addline2,
                    text.details, text.rob,
                    #
                    addrow, addrow.overall,
                    addrow.subgroups, addrows.below.overall,
                    #
                    c(leftcols, rightcols), labs,
                    text.w.common, text.w.random)
  #
  if (new)
    grid.newpage()
  #
  pushViewport(
    viewport(
      layout =
        grid.layout(
          nrow, length(x1), widths = x1,
          heights = unit(spacing, "lines"))))
  #
  # Left side of forest plot
  #
  j <- 1
  #
  if (lsel) {
    #
    # Add text for label.e and label.c (if position was specified by the user)
    #
    if (!is.na(yHeadadd)) {
      if (!is.null(label.e.attach)) {
        vars.e <- paste0("col.", label.e.attach)
        if (all(vars.e %in% leftcols)) {
          id.e <- seq_along(leftcols)[leftcols %in% vars.e]
          id.e <- 2 * (range(id.e) - 1) + 1
          id.e <- seq(min(id.e), max(id.e))
          #
          add.text(col.label.e, id.e)
        }
      }
      #
      if (!is.null(label.c.attach)) {
        vars.c <- paste0("col.", label.c.attach)
        if (all(vars.c %in% leftcols)) {
          id.c <- seq_along(leftcols)[leftcols %in% vars.c]
          id.c <- 2 * (range(id.c) - 1) + 1
          id.c <- seq(min(id.c), max(id.c))
          #
          add.text(col.label.c, id.c)
        }
      }
    }
    #
    for (i in seq_along(leftcols)) {
      add.text(cols[[leftcols[i]]], j)
      #
      if (!is.na(yHeadadd)) {
        if (is.null(label.e.attach)) {
          if (metabin) {
            if (leftcols[i] == "col.n.e" & just.label.e == "right")
              add.text(col.label.e, j)
            else if (leftcols[i] == "col.event.e" &
                     just.label.e %in% c("left", "center"))
              add.text(col.label.e, j)
          }
          else if (metacont) {
            if (leftcols[i] == "col.sd.e" & just.label.e == "right")
              add.text(col.label.e, j)
            else if (leftcols[i] == "col.mean.e" &
                     just.label.e %in% c("left", "center"))
              add.text(col.label.e, j)
          }
          else if (metainc) {
            if (leftcols[i] == "col.time.e" & just.label.e == "right")
              add.text(col.label.e, j)
            else if (leftcols[i] == "col.event.e" &
                     just.label.e %in% c("left", "center"))
              add.text(col.label.e, j)
          }
          else if (metamean) {
            if (revman5 & leftcols[i] == "col.n.e" &
                just.label.e == "right")
              add.text(col.label.e, j)
            else if (!revman5 & leftcols[i] == "col.sd.e" &
                     just.label.e == "right")
              add.text(col.label.e, j)
            else if (leftcols[i] == "col.sd.e" &
                     just.label.e %in% c("left", "center"))
              add.text(col.label.e, j)
          } 
        }
        #
        if (is.null(label.c.attach)) {
          if (metabin) {
            if (leftcols[i] == "col.n.c" & just.label.c == "right")
              add.text(col.label.c, j)
            else if (leftcols[i] == "col.event.c" &
                     just.label.c %in% c("left", "center"))
              add.text(col.label.c, j)
          }
          else if (metacont) {
            if (leftcols[i] == "col.sd.c" & just.label.c == "right")
              add.text(col.label.c, j)
            else if (leftcols[i] == "col.mean.c" &
                     just.label.c %in% c("left", "center"))
              add.text(col.label.c, j)
          }
          else if (metainc) {
            if (leftcols[i] == "col.time.c" & just.label.c == "right")
              add.text(col.label.c, j)
            else if (leftcols[i] == "col.event.c" &
                     just.label.c %in% c("left", "center"))
              add.text(col.label.c, j)
          }
        }
        #
        if (newline.studlab & leftcols[i] == "col.studlab")
          add.text(col.add.studlab, j)
        if (newline.effect & leftcols[i] == "col.effect")
          add.text(col.add.effect, j)
        if (newline.ci & leftcols[i] == "col.ci")
          add.text(col.add.ci, j)
        if (newline.effect.ci & leftcols[i] == "col.effect.ci")
          add.text(col.add.effect.ci, j)
        if (newline.event.n.e & leftcols[i] == "col.event.n.e")
          add.text(col.add.event.n.e, j)
        if (newline.event.n.c & leftcols[i] == "col.event.n.c")
          add.text(col.add.event.n.c, j)
        if (newline.mean.sd.n.e & leftcols[i] == "col.mean.sd.n.e")
          add.text(col.add.mean.sd.n.e, j)
        if (newline.mean.sd.n.c & leftcols[i] == "col.mean.sd.n.c")
          add.text(col.add.mean.sd.n.c, j)
        if (newline.w.common & leftcols[i] == "col.w.common")
          add.text(col.add.w.common, j)
        if (newline.w.random & leftcols[i] == "col.w.random")
          add.text(col.add.w.random, j)
        if (newline.TE & leftcols[i] == "col.TE")
          add.text(col.add.TE, j)
        if (newline.seTE & leftcols[i] == "col.seTE")
          add.text(col.add.seTE, j)
        if (newline.cluster & leftcols[i] == "col.cluster")
          add.text(col.add.cluster, j)
        if (newline.cycles & leftcols[i] == "col.cycles")
          add.text(col.add.cycles, j)
        if (newline.n.e & leftcols[i] == "col.n.e")
          add.text(col.add.n.e, j)
        if (newline.n.c & leftcols[i] == "col.n.c")
          add.text(col.add.n.c, j)
        if (newline.event.e & leftcols[i] == "col.event.e")
          add.text(col.add.event.e, j)
        if (newline.event.c & leftcols[i] == "col.event.c")
          add.text(col.add.event.c, j)
        if (newline.mean.e & leftcols[i] == "col.mean.e")
          add.text(col.add.mean.e, j)
        if (newline.mean.c & leftcols[i] == "col.mean.c")
          add.text(col.add.mean.c, j)
        if (newline.sd.e & leftcols[i] == "col.sd.e")
          add.text(col.add.sd.e, j)
        if (newline.sd.c & leftcols[i] == "col.sd.c")
          add.text(col.add.sd.c, j)
        if (newline.cor & leftcols[i] == "col.cor")
          add.text(col.add.cor, j)
        if (newline.time.e & leftcols[i] == "col.time.e")
          add.text(col.add.time.e, j)
        if (newline.time.c & leftcols[i] == "col.time.c")
          add.text(col.add.time.c, j)
        #
        # Add text in first line of forest plot for new columns
        #
        if (newcols)
          if (length(leftcols.new) > 0 &
              leftcols[i] %in% paste0("col.", leftcols.new)) {
            sel <- paste0("col.", leftcols.new) == leftcols[i]
            #
            # Check for "\n" in label of new column
            #
            clines <- twolines(leftlabs.new[sel], leftcols[i])
            #
            just.new <- just.addcols.left[sel]
            #
            if (just.new == "left")
              xpos.new <- 0
            else if (just.new == "center")
              xpos.new <- 0.5
            else if (just.new == "right")
              xpos.new <- 1
            #
            # Add first line
            #
            if (clines$newline)
              add.text(tgl(clines$top, xpos.new, just.new, fs.head, ff.head,
                           fontfamily), j)
          }
      }
      #
      j <- j + 2
    }
  }
  #
  # Produce forest plot
  #
  draw.lines(col.forest, j,
             ref, TE.common, unique(TE.random),
             overall, common, random, prediction,
             ymin.common, ymin.random, ymin.ref,
             ymax + 0.5 * header.line * addrow,
             ymax.ref + 0.5 * header.line * addrow,
             lwd, lty.common, lty.random, col.common, col.random,
             xlim[1], xlim[2],
             cid.below.null, cid.above.null, lty.cid, col.cid,
             fill.cid.below.null, fill.cid.above.null,
             fill,
             col.lines)
  #
  draw.axis(col.forest, j, yS, log.xaxis, at, label,
            fs.axis, ff.axis, fontfamily, lwd,
            xlim, avail.xlim,
            col.lines, col.label)
  #
  if (bottom.lr) {
    add.text(smlab1, j, xscale = col.forest$range)
    #
    if (newline.smlab)
      add.text(smlab2, j, xscale = col.forest$range)
  }
  #
  if (print.label) {
    if (!bottom.lr) {
      if (!is.na(ref)) {
        add.text(ll1, j, xscale = col.forest$range)
        #
        if (newline.ll)
          add.text(ll2, j, xscale = col.forest$range)
        #
        add.text(lr1, j, xscale = col.forest$range)
        #
        if (newline.lr)
          add.text(lr2, j, xscale = col.forest$range)
      }
    }
    else {
      add.label(ll1, j,
                if (bmj)
                  unit(xlim[1], "native")
                else
                  unit(ref - (xlim[2] - xlim[1]) / 30, "native"),
                unit(y.bottom.lr, "lines"),
                if (bmj) "left" else "right",
                fs.lr, ff.lr, col.label.left, fontfamily,
                xscale = col.forest$range)
      #
      if (newline.ll)
        add.label(ll2, j,
                  if (bmj)
                    unit(xlim[1], "native")
                  else
                    unit(ref - (xlim[2] - xlim[1]) / 30, "native"),
                  unit(y.bottom.lr - 1, "lines"),
                  if (bmj) "left" else "right",
                  fs.lr, ff.lr, col.label.left, fontfamily,
                  xscale = col.forest$range)
      #
      add.label(lr1, j,
                if (bmj)
                  unit(xlim[2], "native")
                else
                  unit(ref + (xlim[2] - xlim[1]) / 30, "native"),
                unit(y.bottom.lr, "lines"),
                if (bmj) "right" else "left",
                fs.lr, ff.lr, col.label.right, fontfamily,
                xscale = col.forest$range)
      #
      if (newline.lr)
        add.label(lr2, j,
                  if (bmj)
                    unit(xlim[2], "native")
                  else
                    unit(ref + (xlim[2] - xlim[1]) / 30, "native"),
                  unit(y.bottom.lr - 1, "lines"),
                  if (bmj) "right" else "left",
                  fs.lr, ff.lr, col.label.right, fontfamily,
                  xscale = col.forest$range)
    }
  }
  #
  add.xlab(col.forest, j, xlab, xlab.add, newline.xlab,
           xlab.pos, xlab.ypos, fs.xlab, ff.xlab,
           fontfamily)
  #
  draw.forest(col.forest, j)
  #
  j <- j + 2
  #
  #
  # Right side of forest plot
  #
  #
  if (rsel | RoB.available) {
    #
    # Add text for label.e and label.c (if position was specified by the user)
    #
    if (!is.na(yHeadadd)) {
      if (!is.null(label.e.attach)) {
        vars.e <- paste0("col.", label.e.attach)
        if (all(vars.e %in% rightcols)) {
          id.e <- seq_along(rightcols)[rightcols %in% vars.e]
          id.e <- 2 * (range(id.e) - 1)
          id.e <- seq(min(id.e), max(id.e))
          #
          add.text(col.label.e, j + id.e)
        }
      }
      #
      if (!is.null(label.c.attach)) {
        vars.c <- paste0("col.", label.c.attach)
        if (all(vars.c %in% rightcols)) {
          id.c <- seq_along(rightcols)[rightcols %in% vars.c]
          id.c <- 2 * (range(id.c) - 1)
          id.c <- seq(min(id.c), max(id.c))
          #
          add.text(col.label.c, j + id.c)
        }
      }
    }
    #
    i.rob <- 0
    #
    for (i in seq_along(rightcols)) {
      if (substring(rightcols[i], 1, 8) == "col.RoB.") {
        i.rob <- i.rob + 1
        #
        add.rob(cols[[rightcols[i]]], j, 0.85, fs.rob.symbols, ff.rob.symbols,
                fontfamily,
                rob[[rightcols[[i]]]],
                rob.categories[[i.rob]], rob.symbols[[i.rob]], rob.col[[i.rob]])
      }
      else
        add.text(cols[[rightcols[i]]], j)
      #
      if (!is.na(yHeadadd)) {
        if (is.null(label.e.attach)) {
          if (metabin) {
            if (rightcols[i] == "col.n.e" & just.label.e == "right")
              add.text(col.label.e, j)
            else if (rightcols[i] == "col.event.e" &
                     just.label.e %in% c("left", "center"))
              add.text(col.label.e, j)
          }
          else if (metacont) {
            if (rightcols[i] == "col.sd.e" & just.label.e == "right")
              add.text(col.label.e, j)
            else if (rightcols[i] == "col.mean.e" &
                     just.label.e %in% c("left", "center"))
              add.text(col.label.e, j)
          }
          else if (metainc) {
            if (rightcols[i] == "col.time.e" & just.label.e == "right")
              add.text(col.label.e, j)
            else if (rightcols[i] == "col.event.e" &
                     just.label.e %in% c("left", "center"))
              add.text(col.label.e, j)
          }
          else if (metamean) {
            if (revman5 & rightcols[i] == "col.n.e" & just.label.e == "right")
              add.text(col.label.e, j)
            else if (!revman5 & rightcols[i] == "col.sd.e" & just.label.e == "right")
              add.text(col.label.e, j)
            else if (rightcols[i] == "col.sd.e" & just.label.e %in% c("left", "center"))
              add.text(col.label.e, j)
          }
        }
        #
        if (is.null(label.c.attach)) {
          if (metabin) {
            if (rightcols[i] == "col.n.c" & just.label.c == "right")
              add.text(col.label.c, j)
            else if (rightcols[i] == "col.event.c" &
                     just.label.c %in% c("left", "center"))
              add.text(col.label.c, j)
          }
          else if (metacont) {
            if (rightcols[i] == "col.sd.c" & just.label.c == "right")
              add.text(col.label.c, j)
            else if (rightcols[i] == "col.mean.c" &
                     just.label.c %in% c("left", "center"))
              add.text(col.label.c, j)
          }
          else if (metainc) {
            if (rightcols[i] == "col.time.c" & just.label.c == "right")
              add.text(col.label.c, j)
            else if (rightcols[i] == "col.event.c" &
                     just.label.c %in% c("left", "center"))
              add.text(col.label.c, j)
          }
        }
        #
        if (!is.null(rob.attach)) {
          if (rightcols[i] == rob.attach)
            add.text(col.rob, j)
        }
        #
        if (newline.studlab & rightcols[i] == "col.studlab")
          add.text(col.add.studlab, j)
        if (newline.effect & rightcols[i] == "col.effect")
          add.text(col.add.effect, j)
        if (newline.ci & rightcols[i] == "col.ci")
          add.text(col.add.ci, j)
        if (newline.effect.ci & rightcols[i] == "col.effect.ci")
          add.text(col.add.effect.ci, j)
        if (newline.mean.sd.n.e & rightcols[i] == "col.mean.sd.n.e")
          add.text(col.add.mean.sd.n.e, j)
        if (newline.mean.sd.n.c & rightcols[i] == "col.mean.sd.n.c")
          add.text(col.add.mean.sd.n.c, j)
        if (newline.w.common & rightcols[i] == "col.w.common")
          add.text(col.add.w.common, j)
        if (newline.w.random & rightcols[i] == "col.w.random")
          add.text(col.add.w.random, j)
        if (newline.TE & rightcols[i] == "col.TE")
          add.text(col.add.TE, j)
        if (newline.seTE & rightcols[i] == "col.seTE")
          add.text(col.add.seTE, j)
        if (newline.cluster & rightcols[i] == "col.cluster")
          add.text(col.add.cluster, j)
        if (newline.cycles & rightcols[i] == "col.cycles")
          add.text(col.add.cycles, j)
        if (newline.n.e & rightcols[i] == "col.n.e")
          add.text(col.add.n.e, j)
        if (newline.n.c & rightcols[i] == "col.n.c")
          add.text(col.add.n.c, j)
        if (newline.event.e & rightcols[i] == "col.event.e")
          add.text(col.add.event.e, j)
        if (newline.event.c & rightcols[i] == "col.event.c")
          add.text(col.add.event.c, j)
        if (newline.mean.e & rightcols[i] == "col.mean.e")
          add.text(col.add.mean.e, j)
        if (newline.mean.c & rightcols[i] == "col.mean.c")
          add.text(col.add.mean.c, j)
        if (newline.sd.e & rightcols[i] == "col.sd.e")
          add.text(col.add.sd.e, j)
        if (newline.sd.c & rightcols[i] == "col.sd.c")
          add.text(col.add.sd.c, j)
        if (newline.cor & rightcols[i] == "col.cor")
          add.text(col.add.cor, j)
        if (newline.time.e & rightcols[i] == "col.time.e")
          add.text(col.add.time.e, j)
        if (newline.time.c & rightcols[i] == "col.time.c")
          add.text(col.add.time.c, j)
        #
        # Add text in first line of forest plot for new columns
        #
        if (newcols)
          if (length(rightcols.new) > 0 &
              rightcols[i] %in% paste0("col.", rightcols.new)) {
            sel <- paste0("col.", rightcols.new) == rightcols[i]
            #
            # Check for "\n" in label of new column
            #
            clines <- twolines(rightlabs.new[sel], rightcols[i])
            #
            just.new <- just.addcols.right[sel]
            #
            if (just.new == "left")
              xpos.new <- 0
            else if (just.new == "center")
              xpos.new <- 0.5
            else if (just.new == "right")
              xpos.new <- 1
            #
            # Add first line
            #
            if (clines$newline)
              add.text(tgl(clines$top, xpos.new, just.new,
                           fs.head, ff.head, fontfamily), j)
          }
      }
      #
      j <- j + 2
    }
  }
  #
  # Add header line
  #
  if (jama)
    hcols <- lsel * 2 * length(leftcols)
  else
    hcols <-
      lsel * 2 * length(leftcols) + 1 + rsel * 2 * length(rightcols)
  #
  if (ev.n.bin) {
    sel1 <- grep("col.event.n.e",
                 c(leftcols, if (all(rightcols != "col.")) rightcols))
    sel2 <- grep("col.event.n.c",
                 c(leftcols, if (all(rightcols != "col.")) rightcols))
    #
    if (length(sel1) > 0 & length(sel2) > 0) {
      if (sel1 > sel2) {
        sel3 <- sel2
        sel2 <- sel1
        sel1 <- sel2
      }
      #
      for (i in seq(2 * sel1 - 1, 2 * sel2 - 1)) {
        pushViewport(
          viewport(
            layout.pos.col = i,
            xscale = col.forest$range))
        #
        grid.lines(x = unit(0:1, "npc"),
                   y = unit(nrow - 1.5 + 0.5 * addrow, "lines"),
                   gp = gpar(lwd = lwd))
        #
        popViewport()
      }
    }
  }
  #
  if (m.s.n.cont) {
    sel1 <- grep("col.mean.sd.n.e",
                 c(leftcols, if (all(rightcols != "col.")) rightcols))
    sel2 <- grep("col.mean.sd.n.c",
                 c(leftcols, if (all(rightcols != "col.")) rightcols))
    #
    if (length(sel1) > 0 & length(sel2) > 0) {
      if (sel1 > sel2) {
        sel3 <- sel2
        sel2 <- sel1
        sel1 <- sel2
      }
      #
      for (i in seq(2 * sel1 - 1, 2 * sel2 - 1)) {
        pushViewport(
          viewport(
            layout.pos.col = i,
            xscale = col.forest$range))
        #
        grid.lines(x = unit(0:1, "npc"),
                   y = unit(nrow - 1.5 + 0.5 * addrow, "lines"),
                   gp = gpar(lwd = lwd))
        #
        popViewport()
      }
    }
  }
  #
  if (header.line) {
    if (header.line.pos == "both") {
      for (i in seq_len(hcols)) {
        pushViewport(
          viewport(
            layout.pos.col = i,
            xscale = col.forest$range))
        #
        grid.lines(x = unit(0:1, "npc"),
                   y = unit(nrow + 0.5 * addrow, "lines"),
                   gp = gpar(lwd = lwd, col = col.header.line))
        #
        popViewport()
      }
    }
    #
    for (i in seq_len(hcols)) {
      pushViewport(
        viewport(
          layout.pos.col = i,
          xscale = col.forest$range))
      #
      grid.lines(x = unit(0:1, "npc"),
                 y = unit(ymax + 0.5 * addrow, "lines"),
                 gp = gpar(lwd = lwd, col = col.header.line))
      #
      popViewport()
    }
  }
  #
  # Add JAMA lines
  #
  if (jama & header.line & !by) {
    for (i in seq_len(hcols)) {
      pushViewport(
        viewport(
          layout.pos.col = i,
          xscale = col.forest$range))
      #
      for (j in seq_len(k.all + n.com * common + n.ran * random +
                        n.prd * prediction))
        grid.lines(x = unit(0:1, "npc"),
                   y = unit(ymax + 0.5 * addrow - j, "lines"),
                   gp = gpar(lwd = 0.5 * lwd, col = col.jama.line))
      #
      popViewport()
    }
  }
  #
  popViewport()
  #
  if (dev.off)
    invisible(dev.off())

  
  res <- list(xlim = xlim, addrows.below.overall = addrows.below.overall,
              #
              colgap = colgap,
              colgap.left = colgap.left,
              colgap.right = colgap.right,
              colgap.studlab = colgap.studlab,
              colgap.forest = colgap.forest.left,
              colgap.forest.left = colgap.forest,
              colgap.forest.right = colgap.forest.right,
              #
              studlab = studlab,
              TE.format = TE.format,
              seTE.format = seTE.format,
              cluster.format = cluster.format,
              cycles.format = cycles.format,
              effect.format = effect.format,
              ci.format = ci.format,
              effect.ci.format = effect.ci.format)
  #
  if (ev.n.bin | ev.n.prop)
    res$effect.ci.format <- effect.ci.format
  #
  if (ev.n.bin)
    res$effect.ci.format <- effect.ci.format
  #
  if (length(cols.new) > 0) {
    for (i in names(cols.new))
      res[[i]] <- cols.new[[i]]
  }
  #
  res$figheight <- figheight
  #
  res$leftcols <- leftcols
  res$leftlabs <- leftlabs
  res$rightcols <- rightcols
  res$rightlabs <- rightlabs
  #
  invisible(res)
}





#' @rdname forest.meta
#' @method plot meta
#' @export
#'

plot.meta <- function(x, ...)
  forest(x, ...)





#' @rdname forest.meta
#' @export .forestArgs

.forestArgs <- function()
  formalArgs(forest.meta)
