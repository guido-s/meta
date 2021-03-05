#' Forest plot to display the result of a meta-analysis
#' 
#' @description
#' Draws a forest plot in the active graphics window (using grid
#' graphics system).
#' 
#' @aliases forest forest.meta
#' 
#' @param x An object of class \code{meta}.
#' @param sortvar An optional vector used to sort the individual
#'   studies (must be of same length as \code{x$TE}).
#' @param studlab A logical indicating whether study labels should be
#'   printed in the graph. A vector with study labels can also be
#'   provided (must be of same length as \code{x$TE} then).
#' @param layout A character string specifying the layout of the
#'   forest plot (see Details).
#' @param comb.fixed A logical indicating whether fixed effect
#'   estimate should be plotted.
#' @param comb.random A logical indicating whether random effects
#'   estimate should be plotted.
#' @param overall A logical indicating whether overall summaries
#'   should be plotted. This argument is useful in a meta-analysis
#'   with subgroups if summaries should only be plotted on group
#'   level.
#' @param text.fixed A character string used in the plot to label the
#'   pooled fixed effect estimate.
#' @param text.random A character string used in the plot to label the
#'   pooled random effects estimate.
#' @param lty.fixed Line type (pooled fixed effect estimate).
#' @param lty.random Line type (pooled random effects estimate).
#' @param col.fixed Line colour (pooled fixed effect estimate).
#' @param col.random Line colour (pooled random effects estimate).
#' @param text.w.fixed A character string used to label weights of
#'   fixed effect model.
#' @param text.w.random A character string used to label weights of
#'   random effects model.
#' @param prediction A logical indicating whether a prediction
#'   interval should be printed.
#' @param text.predict A character string used in the plot to label
#'   the prediction interval.
#' @param subgroup A logical indicating whether subgroup results
#'   should be shown in forest plot. This argument is useful in a
#'   meta-analysis with subgroups if summaries should not be plotted
#'   on group level.
#' @param print.subgroup.labels A logical indicating whether subgroup
#'   label should be printed.
#' @param bylab A character string with a label for the grouping
#'   variable.
#' @param print.byvar A logical indicating whether the name of the
#'   grouping variable should be printed in front of the group labels.
#' @param byseparator A character string defining the separator
#'   between label and levels of grouping variable.
#' @param text.fixed.w A character string to label the pooled fixed
#'   effect estimate within subgroups, or a character vector of same
#'   length as number of subgroups with corresponging labels.
#' @param text.random.w A character string to label the pooled random
#'   effect estimate within subgroups, or a character vector of same
#'   length as number of subgroups with corresponging labels.
#' @param text.predict.w A character string to label the prediction
#'   interval within subgroups, or a character vector of same length
#'   as number of subgroups with corresponging labels.
#' @param bysort A logical indicating whether groups should be ordered
#'   alphabetically.
#' @param pooled.totals A logical indicating whether total number of
#'   observations should be given in the figure.
#' @param pooled.events A logical indicating whether total number of
#'   events should be given in the figure.
#' @param pooled.times A logical indicating whether total person time
#'   at risk should be given in the figure.
#' @param study.results A logical indicating whether results for
#'   individual studies should be shown in the figure (useful to only
#'   plot subgroup results).
#' @param xlab A label for the x-axis.
#' @param xlab.pos A numeric specifying the center of the label on the
#'   x-axis.
#' @param smlab A label for the summary measurex (printed at top of
#'   figure).
#' @param smlab.pos A numeric specifying the center of the label for
#'   the summary measure.
#' @param xlim The x limits (min,max) of the plot, or the character
#'   "s" to produce symmetric forest plots.
#' @param allstudies A logical indicating whether studies with
#'   inestimable treatment effects should be plotted.
#' @param weight.study A character string indicating weighting used to
#'   determine size of squares or diamonds (argument
#'   \code{type.study}) to plot individual study results. One of
#'   missing, \code{"same"}, \code{"fixed"}, or \code{"random"}, can
#'   be abbreviated. Plot symbols have the same size for all studies
#'   or represent study weights from fixed effect or random effects
#'   model.
#' @param weight.subgroup A character string indicating weighting used
#'   to determine size of squares or diamonds (argument
#'   \code{type.subgroup}) to plot subgroup results. One of missing,
#'   \code{"same"}, or \code{"weight"}, can be abbreviated. Plot
#'   symbols have the same size for all subgroup results or represent
#'   subgroup weights from fixed effect or random effects model.
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
#' @param ref A numerical giving the reference value to be plotted as
#'   a line in the forest plot. No reference line is plotted if
#'   argument \code{ref} is equal to \code{NA}.
#' @param lower.equi A numerical giving the lower limit of equivalence
#'   to be plotted as a line in the forest plot. No line is plotted if
#'   argument \code{lower.equi} is equal to \code{NA}.
#' @param upper.equi A numerical giving the upper limit of equivalence
#'   to be plotted as a line in the forest plot. No line is plotted if
#'   argument \code{upper.equi} is equal to \code{NA}.
#' @param lty.equi Line type (limits of equivalence).
#' @param col.equi Line colour (limits of equivalence).
#' @param fill.equi Colour of area between limits of equivalence.
#' @param leftcols A character vector specifying (additional) columns
#'   to be plotted on the left side of the forest plot or a logical
#'   value (see Details).
#' @param rightcols A character vector specifying (additional) columns
#'   to be plotted on the right side of the forest plot or a logical
#'   value (see Details).
#' @param leftlabs A character vector specifying labels for
#'   (additional) columns on left side of the forest plot (see
#'   Details).
#' @param rightlabs A character vector specifying labels for
#'   (additional) columns on right side of the forest plot (see
#'   Details).
#' @param lab.e Label to be used for experimental group in table
#'   heading.
#' @param lab.c Label to be used for control group in table heading.
#' @param lab.e.attach.to.col A character specifying the column name
#'   where label \code{lab.e} should be attached to in table heading.
#' @param lab.c.attach.to.col A character specifying the column name
#'   where label \code{lab.c} should be attached to in table heading.
#' @param label.left Graph label on left side of forest plot.
#' @param label.right Graph label on right side of forest plot.
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
#' @param type.study A character string or vector specifying how to
#'   plot treatment effects and confidence intervals for individual
#'   studies (see Details).
#' @param type.fixed A character string specifying how to plot
#'   treatment effect and confidence interval for fixed effect
#'   meta-analysis (see Details).
#' @param type.random A character string specifying how to plot
#'   treatment effect and confidence interval for random effects
#'   meta-analysis (see Details).
#' @param type.subgroup A character string specifying how to plot
#'   treatment effect and confidence interval for subgroup results
#'   (see Details).
#' @param type.subgroup.fixed A character string specifying how to
#'   plot treatment effect and confidence interval for subgroup
#'   results (fixed effect model).
#' @param type.subgroup.random A character string specifying how to
#'   plot treatment effect and confidence interval for subgroup
#'   results (random effects model).
#' @param col.study The colour for individual study results and
#'   confidence limits.
#' @param col.inside The colour for individual study results and
#'   confidence limits if confidence limits are completely within
#'   squares.
#' @param col.square The colour for squares reflecting study's weight
#'   in the meta-analysis.
#' @param col.square.lines The colour for the outer lines of squares
#'   reflecting study's weight in the meta-analysis.
#' @param col.diamond The colour of diamonds representing the results
#'   for fixed effect and random effects models.
#' @param col.diamond.fixed The colour of diamonds for fixed effect
#'   estimates.
#' @param col.diamond.random The colour of diamonds for random effects
#'   estimates.
#' @param col.diamond.lines The colour of the outer lines of diamonds
#'   representing the results for fixed effect and random effects
#'   models.
#' @param col.diamond.lines.fixed The colour of the outer lines of
#'   diamond for fixed effect estimate.
#' @param col.diamond.lines.random The colour of the outer lines of
#'   diamond for random effects estimate.
#' @param col.inside.fixed The colour for result of fixed effect
#'   meta-analysis if confidence limit lies completely within square.
#' @param col.inside.random The colour for result of random effects
#'   meta-analysis if confidence limit lies completely within square.
#' @param col.predict Background colour of prediction interval.
#' @param col.predict.lines Colour of outer lines of prediction
#'   interval.
#' @param col.by The colour to print information on subgroups.
#' @param col.label.right The colour for label on right side of null
#'   effect.
#' @param col.label.left The colour for label on left side of null
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
#'   generalized linear mixed models.
#' @param test.overall A logical value indicating whether to print
#'   results of test for overall effect.
#' @param test.overall.fixed A logical value indicating whether to
#'   print results of test for overall effect (based on fixed effect
#'   model).
#' @param test.overall.random A logical value indicating whether to
#'   print results of test for overall effect (based on random effects
#'   model).
#' @param label.test.overall.fixed Label printed in front of results
#'   of test for overall effect (based on fixed effect model).
#' @param label.test.overall.random Label printed in front of results
#'   of test for overall effect (based on random effects model).
#' @param print.stat A logical value indicating whether z- or t-value
#'   for test of treatment effect should be printed.
#' @param test.subgroup A logical value indicating whether to print
#'   results of test for subgroup differences.
#' @param test.subgroup.fixed A logical value indicating whether to
#'   print results of test for subgroup differences (based on fixed
#'   effect model).
#' @param test.subgroup.random A logical value indicating whether to
#'   print results of test for subgroup differences (based on random
#'   effects model).
#' @param print.Q.subgroup A logical value indicating whether to print
#'   the value of the heterogeneity statistic Q (test for subgroup
#'   differences).
#' @param label.test.subgroup.fixed Label printed in front of results
#'   of test for subgroup differences (based on fixed effect model).
#' @param label.test.subgroup.random Label printed in front of results
#'   of test for subgroup differences (based on random effects model).
#' @param test.effect.subgroup A logical value indicating whether to
#'   print results of test for effect in subgroups.
#' @param test.effect.subgroup.fixed A logical value indicating
#'   whether to print results of test for effect in subgroups (based
#'   on fixed effect model).
#' @param test.effect.subgroup.random A logical value indicating
#'   whether to print results of test for effect in subgroups (based
#'   on random effects model).
#' @param label.test.effect.subgroup.fixed Label printed in front of
#'   results of test for effect in subgroups (based on fixed effect
#'   model).
#' @param label.test.effect.subgroup.random Label printed in front of
#'   results of test for effect in subgroups (based on random effects
#'   model).
#' @param text.addline1 Text for first additional line (below
#'   meta-analysis results).
#' @param text.addline2 Text for second additional line (below
#'   meta-analysis results).
#' @param fontsize The size of text (in points), see
#'   \code{\link{gpar}}.
#' @param fontfamily The font family, see \code{\link{gpar}}.
#' @param fs.heading The size of text for column headings, see
#'   \code{\link{gpar}}.
#' @param fs.fixed The size of text for results of fixed effect model,
#'   see \code{\link{gpar}}.
#' @param fs.random The size of text for results of random effects
#'   model, see \code{\link{gpar}}.
#' @param fs.predict The size of text for results of prediction
#'   interval, see \code{\link{gpar}}.
#' @param fs.fixed.labels The size of text for label of fixed effect
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
#' @param ff.heading The fontface for column headings, see
#'   \code{\link{gpar}}.
#' @param ff.fixed The fontface of text for results of fixed effect
#'   model, see \code{\link{gpar}}.
#' @param ff.random The fontface of text for results of random effects
#'   model, see \code{\link{gpar}}.
#' @param ff.predict The fontface of text for results of prediction
#'   interval, see \code{\link{gpar}}.
#' @param ff.fixed.labels The fontface of text for label of fixed
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
#' @param squaresize A numeric used to increase or decrease the size
#'   of squares in the forest plot.
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
#' @param calcwidth.pooled A logical indicating whether text for fixed
#'   effect and random effects model should be considered to calculate
#'   width of the column with study labels.
#' @param calcwidth.fixed A logical indicating whether text given in
#'   arguments \code{text.fixed} and \code{text.fixed.w} should be
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
#' @param spacing A numeric determining line spacing in a forest plot.
#' @param addrow A logical value indicating whether an empty row is
#'   printed above and below study results.
#' @param addrow.overall A logical value indicating whether an empty
#'   row is printed above overall meta-analysis results.
#' @param addrow.subgroups A logical value indicating whether an empty
#'   row is printed between results for subgroups.
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
#'   errors, see \code{print.default}.
#' @param digits.stat Minimal number of significant digits for z- or
#'   t-statistic for test of overall effect, see \code{print.default}.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance, see \code{print.default}.
#' @param digits.tau Minimal number of significant digits for square
#'   root of between-study variance, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall treatment effect, see \code{print.default}.
#' @param digits.pval.Q Minimal number of significant digits for
#'   p-value of heterogeneity test, see \code{print.default}.
#' @param digits.Q Minimal number of significant digits for
#'   heterogeneity statistic Q, see \code{print.default}.
#' @param digits.I2 Minimal number of significant digits for I-squared
#'   statistic, see \code{print.default}.
#' @param digits.weight Minimal number of significant digits for
#'   weights, see \code{print.default}.
#' @param digits.mean Minimal number of significant digits for means;
#'   only applies to \code{\link{metacont}} objects.
#' @param digits.sd Minimal number of significant digits for standard
#'   deviations; only applies to \code{\link{metacont}} objects.
#' @param digits.cor Minimal number of significant digits for
#'   correlations; only applies to \code{\link{metacor}} objects.
#' @param digits.time Minimal number of significant digits for times;
#'   only applies to \code{\link{metainc}} and \code{\link{metarate}}
#'   objects.
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
#' @param col.i Deprecated argument (replaced by \code{col.study}).
#' @param weight Deprecated argument (replaced by
#'   \code{weight.study}).
#' @param digits.zval Deprecated argument (replaced by
#'   \code{digits.stat}).
#' @param print.zval Deprecated argument (replaced by
#'   \code{print.stat}).
#' @param \dots Additional graphical arguments.
#' 
#' @details
#' A forest plot, also called confidence interval plot, is drawn in
#' the active graphics window. The forest functions in R package
#' \bold{meta} are based on the grid graphics system. In order to
#' print the forest plot, resize the graphics window and either use
#' \code{\link{dev.copy2eps}} or \code{\link{dev.copy2pdf}}. Another
#' possibility is to create a file using \code{\link{pdf}},
#' \code{\link{png}}, or \code{\link{svg}} and to specify the width and
#' height of the graphic (see Examples).
#' 
#' By default, treatment estimates and confidence intervals are
#' plotted in the following way:
#' \itemize{
#' \item For an individual study, a square with treatment estimate in
#'   the center and confidence interval as line extending either side
#'   of the square (\code{type.study = "square"})
#' \item For meta-analysis results, a diamond with treatment estimate
#'   in the center and right and left side corresponding to lower and
#'   upper confidence limits (\code{type.fixed = "diamond"},
#'   \code{type.random = "diamond"}, and \code{type.subgroup = "diamond"})
#' }
#' 
#' In a forest plot, size of the squares typically reflects the precision of
#' individual treatment estimates based either on the fixed effect
#' (\code{weight.study = "fixed"}) or random effects meta-analysis
#' (\code{weight.study = "random"}). Information from meta-analysis object
#' \code{x} is utilised if argument \code{weight.study} is missing. Weights
#' from the fixed effect model are used if argument \code{x$comb.fixed} is
#' \code{TRUE}; weights from the random effects model are used if argument
#' \code{x$comb.random} is \code{TRUE} and \code{x$comb.fixed} is \code{FALSE}.
#' The same square sizes are used if \code{weight.study = "same"}.
#' 
#' Arguments \code{text.fixed}, \code{text.random}, and
#' \code{text.predict} can be used to change the label to identify
#' overall results (fixed effect and random effects model as well as
#' prediction interval). By default the following text is printed:
#' \itemize{
#' \item "Fixed effect model" (argument \code{text.fixed})
#' \item "Random effects model" (\code{text.random})
#' \item "Prediction interval" (\code{text.predict})
#' }
#'
#' If confidence interval levels are different for individual studies,
#' meta-analysis, and prediction interval (arguments \code{level},
#' \code{level.comb}, \code{level.predict} in meta-analysis functions,
#' e.g., \code{\link{metabin}}), additional information is printed,
#' e.g., " (99\%-CI)" for a 99\% confidence interval in the
#' meta-analysis.
#' 
#' The following arguments can be used to print results for various
#' statistical tests:
#' \tabular{ll}{
#' \bold{Argument} \tab \bold{Statistical test} \cr
#' \code{test.overall.fixed} \tab Test for overall effect (fixed
#'   effect model) \cr
#' \code{test.overall.random} \tab Test for overall effect (random
#'   effects model) \cr
#' \code{test.effect.subgroup.fixed} \tab Test for effect in subgroup
#'   (FE model) \cr
#' \code{test.effect.subgroup.random} \tab Test for effect in subgroup
#'   (RE model) \cr
#' \code{test.subgroup.fixed} \tab Test for subgroup differences (FE
#'   model) \cr
#' \code{test.subgroup.random} \tab Test for subgroup differences (RE
#'   model)
#' }
#' 
#' By default, these arguments are \code{FALSE}. R function
#' \code{\link{settings.meta}} can be used to change this default for
#' the entire R session. For example, use the following command to
#' always print results of tests for an overall effect:
#' \code{settings.meta(test.overall = TRUE)}
#' 
#' The arguments \code{leftcols} and \code{rightcols} can be used to
#' specify columns which are plotted on the left and right side of the
#' forest plot, respectively. If argument \code{rightcols} is
#' \code{FALSE}, no columns will be plotted on the right side. By
#' default, i.e. if arguments \code{leftcols} and \code{rightcols} are
#' \code{NULL} and \code{layout = "meta"}, the following
#' \emph{\bold{columns}} will be printed \emph{\bold{on the right side
#' of the forest plot}}:
#' \tabular{ll}{
#' \bold{Meta-analysis results} \tab \bold{Value of argument
#'   rightcols} \cr
#' No summary \tab \code{c("effect", "ci")} \cr
#' Only fixed effect model \tab \code{c("effect", "ci", "w.fixed")}
#'   \cr
#' Only random effects model \tab \code{c("effect", "ci", "w.random")}
#'   \cr
#' Both models \tab \code{c("effect", "ci", "w.fixed", "w.random")}
#' }
#'
#' By default, estimated treatment effect and corresponding confidence
#' interval will be printed.  Depending on arguments \code{comb.fixed}
#' and \code{comb.random}, weights of the fixed effect and/or random
#' effects model will be given too. For an object of class
#' \code{\link{metacum}} or \code{\link{metainf}} only the estimated
#' treatment effect with confidence interval are plotted.
#' 
#' Depending on the class of the meta-analysis object (which is
#' defined by the R function used to generate the object) a different
#' set of \emph{\bold{columns}} is printed \emph{\bold{on the left
#' side of the forest plot}}:
#' \tabular{cl}{
#' \bold{Function} \tab \bold{Value of argument leftcols} \cr
#' \code{\link{metabin}} \tab \code{c("studlab", "event.e", "n.e",
#'   "event.c", "n.c")} \cr
#' \code{\link{metacont}} \tab \code{c("studlab", "n.e", "mean.e",
#'   "sd.e", "n.c", "mean.c", "sd.c")} \cr
#' \code{\link{metacor}} \tab \code{c("studlab", "n")} \cr
#' \code{\link{metagen}} \tab \code{c("studlab", "TE", "seTE")} \cr
#' \code{\link{metainc}} \tab \code{c("studlab", "event.e", "time.e",
#'   "event.c", "time.c")} \cr
#' \code{\link{metaprop}} \tab \code{c("studlab", "event", "n")} \cr
#' \code{\link{metarate}} \tab \code{c("studlab", "event", "time")}
#'   \cr
#' \code{\link{metacum}} \tab \code{"studlab"} \cr
#' \code{\link{metainf}} \tab \code{"studlab"}
#' }
#'
#' The arguments \code{leftlabs} and \code{rightlabs} can be used to
#' specify column headings which are plotted on left and right side of
#' the forest plot, respectively. For certain columns predefined
#' labels exist. If the arguments \code{leftlabs} and \code{rightlabs}
#' are \code{NULL}, the following default labels will be used:
#' \tabular{rcccccc}{
#' \bold{Column:} \tab \code{studlab} \tab \code{TE} \tab \code{seTE}
#'   \tab \code{n.e} \tab \code{n.c} \tab \code{n} \cr
#' \bold{Label:} \tab "Study" \tab "TE" \tab "seTE" \tab "Total" \tab
#'   "Total" \tab "Total" \cr
#' \cr
#' \bold{Column:} \tab \code{event.e} \tab \code{event.c} \tab
#'   \code{event} \tab \code{mean.e} \tab \code{mean.c} \tab \cr
#' \bold{Label:} \tab "Events" \tab "Events" \tab "Events" \tab "Mean"
#'   \tab "Mean" \tab \cr
#' \cr
#' \bold{Column:} \tab \code{sd.e} \tab \code{sd.c} \tab \code{time.e}
#'   \tab \code{time.c} \tab \code{effect} \tab \cr
#' \bold{Label:} \tab "SD" \tab "SD" \tab "Time" \tab "Time" \tab
#'   \code{x$sm} \tab \cr
#' \cr
#' \bold{Column:} \tab \code{ci} \tab \code{effect.ci} \tab
#'   \code{w.fixed} \tab \code{w.random} \tab \tab \cr
#' \bold{Label:} \tab \code{x$level}"\%-CI" \tab \emph{effect+ci} \tab
#'   "W(fixed)" \tab "W(random)" \tab \tab
#' }
#'
#' For additional columns, the column name will be used as a label. It
#' is possible to only provide labels for new columns (see
#' Examples). Otherwise the length of \code{leftlabs} and
#' \code{rightlabs} must be the same as the number of printed columns,
#' respectively. The value \code{NA} can be used to specify columns
#' which should use default labels (see Examples).
#' 
#' If argument \code{layout = "RevMan5"} (and arguments \code{leftcols} and
#' \code{rightcols} are \code{NULL}), the layout for forest plots used for
#' Cochrane reviews (which are generated with Review Manager 5,
#' \url{https://training.cochrane.org/online-learning/core-software-cochrane-reviews/revman})
#' is reproduced:
#' \enumerate{
#' \item All columns are printed on the left side of the forest plot
#'   (see arguments \code{leftcols} and \code{rightcols})
#' \item Tests for overall effect and subgroup differences are printed
#'   (\code{test.overall}, \code{test.effect.subgroup},
#'   \code{test.subgroup})
#' \item Diamonds representing meta-analysis results are printed in
#'   black (\code{diamond.fixed}, \code{diamond.random})
#' \item Colour of squares depends on the meta-analysis object
#'   (\code{col.square}, \code{col.square.lines})
#' \item Information on effect measure and meta-analysis method is
#'   printed above the forest plot (\code{smlab})
#' \item Label "Study or Subgroup" is printed for meta-analysis with
#'   subgroups (\code{leftlabs})
#' }
#'
#' If argument \code{layout = "JAMA"} (and arguments \code{leftcols} and
#' \code{rightcols} are \code{NULL}), instructions for authors of the
#' \emph{Journal of the American Medical Association}, see
#' \url{https://jamanetwork.com/journals/jama/pages/instructions-for-authors/},
#' are taken into account:
#' \enumerate{
#' \item Graph labels on right and left side are printed in bold font
#'   at top of forest plot (see arguments \code{bottom.lr} and
#'   \code{ff.lr})
#' \item Information on effect measure and level of confidence
#'   interval is printed at bottom of forest plot (\code{xlab})
#' \item Tests for overall effect are printed (\code{test.overall})
#' \item Diamonds representing meta-analysis results are printed in
#'   lightblue (\code{diamond.fixed}, \code{diamond.random})
#' \item Squares representing individual study results are printed in
#'   darkblue (\code{col.square}, \code{col.square.lines})
#' \item Between-study variance \eqn{\tau^2} is not printed
#' \item Empty rows are omitted (\code{addrow})
#' \item Label "Source" is printed instead of "Study" (\code{leftlabs})
#' \item P-values are printed without leading zeros (\code{zero.pval})
#' \item P-values are rounded to three digits (for 0.001 < p \eqn{\le}
#'   0.01) or two digits (p > 0.01) (\code{JAMA.pval})
#' }
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
#' 
#' If arguments \code{lab.e} and \code{lab.c} are \code{NULL},
#' "Experimental" and "Control" are used as labels for experimental
#' and control group, respectively.
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
#' 
#' A prediction interval for treatment effect of a new study (Higgins
#' et al., 2009) is given in the forest plot if arguments
#' \code{prediction} and \code{comb.random} are \code{TRUE}. For
#' graphical presentation of prediction intervals the approach by
#' Guddat et al. (2012) is used.
#' 
#' Argument \code{hetstat} can be a character string to specify where
#' to print heterogeneity information:
#' \itemize{
#' \item row with results for fixed effect model (\code{hetstat =
#' "fixed"}),
#' \item row with results for random effects model (\code{hetstat =
#' "random"}).
#' }
#' Otherwise, information on heterogeneity is printed in dedicated rows.
#' 
#' Note, in R package \bold{meta}, version 3.0-0 the following
#' arguments have been removed from R function forest.meta: byvar,
#' level, level.comb, level.predict. This functionality is now
#' provided by R function \code{\link{update.meta}} (or directly in R
#' functions, e.g., \code{\link{metabin}}, \code{\link{metacont}},
#' \code{\link{metagen}}, \code{\link{metacor}}, and
#' \code{\link{metaprop}}).
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}, \code{\link{metacont}},
#'   \code{\link{metagen}}, \code{\link{forest.metabind}},
#'   \code{\link{settings.meta}}
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
#'               data = Olkin1995, subset = c(41, 47, 51, 59),
#'               sm = "RR", method = "I",
#'               studlab = paste(author, year))
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
#' forest(m1, layout = "RevMan5", comb.fixed = FALSE,
#'        label.right = "Favours control", col.label.right = "red",
#'        label.left = "Favours experimental", col.label.left = "green",
#'        prediction = TRUE)
#' 
#' 
#' \dontrun{
#' # Create a PDF file forest-m1.pdf with the forest plot
#' #
#' pdf("forest-m1.pdf", width = 10, height = 3)
#' forest(m1)
#' dev.off()
#' 
#' # Sort studies by decreasing treatment effect within year subgroups
#' #
#' m2 <- update(m1, byvar = ifelse(year < 1987,
#'                                 "Before 1987", "1987 and later"),
#'              print.byvar = FALSE)
#' forest(m2, sortvar = -TE, comb.random = FALSE)
#' 
#' # Forest plot specifying argument xlim
#' #
#' forest(m1, xlim = c(0.01, 10))
#' 
#' # Print results of test for overall effect
#' #
#' forest(m1, test.overall.fixed = TRUE, test.overall.random = TRUE)
#' 
#' # Forest plot with 'classic' layout used in R package meta,
#' # version < 1.6-0
#' #
#' forest(m1, col.square = "black", hetstat = FALSE)
#' 
#' # Change set of columns printed on left side of forest plot
#' #
#' forest(m1, comb.random = FALSE, leftcols = "studlab")
#' 
#' # Do not print columns on right side of forest plot
#' #
#' forest(m1, rightcols = FALSE)
#' 
#' # Change study label to "Author"
#' #
#' forest(m1, comb.random = FALSE, leftlabs = c("Author", NA, NA, NA, NA))
#' 
#' # Just give effect estimate and 95% confidence interval on right
#' # side of forest plot (in one column)
#' #
#' forest(m1, rightcols = c("effect.ci"))
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
#'        leftcols = c("studlab", "n.e", "event.e", "n.c", "event.c"),
#'        lab.e.attach.to.col = "event.e",
#'        lab.c.attach.to.col = "event.c")
#' 
#' # Specify column labels only for variables 'year' and 'author'
#' # (and define digits for additional variables)
#' #
#' forest(m1,
#'        leftcols = c("studlab", "event.e", "n.e", "event.c", "n.c",
#'                     "author", "year"),
#'        leftlabs = c("Author", "Year of Publ"))
#' 
#' # Center text in all columns
#' #
#' forest(m1,
#'        leftcols = c("studlab", "event.e", "n.e", "event.c", "n.c",
#'                     "author", "year"),
#'        leftlabs = c("Author", "Year of Publ"), hetstat = FALSE,
#'        just = "center", just.addcols = "center", just.studlab = "center")
#' 
#' # Same result
#' #
#' forest(m1,
#'        leftcols = c("studlab", "event.e", "n.e", "event.c", "n.c",
#'                   "author", "year"),
#'        leftlabs = c("Author", "Year of Publ"), hetstat = FALSE,
#'        just = "c", just.addcols = "c", just.studlab = "c")
#' 
#' # Change some fontsizes and fontfaces
#' #
#' forest(m1,
#'        fs.study = 10, ff.study = "italic",
#'        fs.study.label = 11, ff.study.label = "bold",
#'        fs.axis = 5, ff.axis = "italic",
#'        ff.smlab = "bold.italic",
#'        ff.fixed = "plain", ff.hetstat = "plain")
#' 
#' # Change some colours
#' #
#' forest(m1,
#'        col.diamond = "green", col.diamond.lines = "red",
#'        col.study = c("green", "blue", "red", "orange"),
#'        col.square = "pink", col.square.lines = "black")
#' 
#' # Sort by weight in fixed effect model
#' #
#' forest(m1, sortvar = 1 / w.fixed, comb.random = FALSE)
#' 
#' # Sort by decreasing weight in fixed effect model
#' #
#' forest(m1, sortvar = -1 / w.fixed, comb.random = FALSE)
#' 
#' # Sort by size of treatment effect
#' #
#' forest(m1, sortvar = TE, comb.random = FALSE)
#' 
#' # Sort by size of treatment effect
#' #
#' forest(m1, sortvar = -TE, comb.random = FALSE)
#' 
#' # Sort by decreasing year of publication
#' #
#' forest(m1, sortvar = -year, comb.random = FALSE)
#' 
#' # Print results of test for subgroup differences (random effects
#' # model)
#' #
#' forest(m2,
#'        sortvar = -TE, comb.fixed = FALSE,
#'        test.subgroup.random = TRUE)
#' 
#' # Print only subgroup results
#' #
#' forest(m2, layout = "subgroup")
#' 
#' # Print only subgroup results (and consider text for heterogeneity
#' # measures in width of subgroup column)
#' #
#' forest(m2, layout = "subgroup", calcwidth.hetstat = TRUE)
#' }
#'
#' @method forest meta
#' @export
#' @export forest.meta


forest.meta <- function(x,
                        sortvar,
                        studlab = TRUE,
                        ##
                        layout = gs("layout"),
                        ##
                        comb.fixed = x$comb.fixed,
                        comb.random = x$comb.random,
                        overall = x$overall,
                        text.fixed = x$text.fixed,
                        text.random = x$text.random,
                        lty.fixed = 2, lty.random = 3,
                        col.fixed = "black", col.random = "black",
                        text.w.fixed = x$text.w.fixed,
                        text.w.random = x$text.w.random,
                        ##
                        prediction = x$prediction,
                        text.predict = x$text.predict,
                        ##
                        subgroup = TRUE,
                        print.subgroup.labels = TRUE,
                        bylab = x$bylab,
                        print.byvar = x$print.byvar,
                        byseparator = x$byseparator,
                        text.fixed.w = text.fixed,
                        text.random.w = text.random,
                        text.predict.w = text.predict,
                        bysort = FALSE,
                        ##
                        pooled.totals = comb.fixed | comb.random,
                        pooled.events = FALSE, pooled.times = FALSE,
                        ##
                        study.results = TRUE,
                        ##
                        xlab = "", xlab.pos,
                        smlab = NULL, smlab.pos, xlim = "symmetric",
                        ##
                        allstudies = TRUE,
                        weight.study,
                        weight.subgroup,
                        pscale = x$pscale,
                        irscale = x$irscale, irunit = x$irunit,
                        ##
                        ref =
                          ifelse(backtransf & is.relative.effect(x$sm), 1, 0),
                        ##
                        lower.equi = NA, upper.equi = NA,
                        lty.equi = 1, col.equi = "blue",
                        fill.equi = "transparent",
                        ##
                        leftcols = NULL, rightcols = NULL,
                        leftlabs = NULL, rightlabs = NULL,
                        ##
                        lab.e = x$label.e,
                        lab.c = x$label.c,
                        ##
                        lab.e.attach.to.col = NULL,
                        lab.c.attach.to.col = NULL,
                        ##
                        label.right = x$label.right,
                        label.left = x$label.left,
                        bottom.lr = TRUE,
                        ##
                        lab.NA = ".", lab.NA.effect = "", lab.NA.weight = "--",
                        ##
                        lwd = 1,
                        ##
                        at = NULL,
                        label = TRUE,
                        ##
                        type.study = "square",
                        type.fixed = "diamond",
                        type.random = type.fixed,
                        type.subgroup =
                          ifelse(study.results, "diamond", "square"),
                        type.subgroup.fixed = type.subgroup,
                        type.subgroup.random = type.subgroup,
                        ##
                        col.study = "black",
                        col.square = "gray",
                        col.square.lines = col.square,
                        col.inside = "white",
                        ##
                        col.diamond = "gray",
                        col.diamond.fixed = col.diamond,
                        col.diamond.random = col.diamond,
                        col.diamond.lines = "black",
                        col.diamond.lines.fixed = col.diamond.lines,
                        col.diamond.lines.random = col.diamond.lines,
                        ##
                        col.inside.fixed = col.inside,
                        col.inside.random = col.inside,
                        ##
                        col.predict = "red",
                        col.predict.lines = "black",
                        ##
                        col.by = "darkgray",
                        ##
                        col.label.right = "black",
                        col.label.left = "black",
                        ##
                        hetstat =
                          comb.fixed | comb.random | overall.hetstat,
                        overall.hetstat = x$overall.hetstat,
                        hetlab = "Heterogeneity: ",
                        resid.hetstat,
                        resid.hetlab = "Residual heterogeneity: ",
                        print.I2,
                        print.I2.ci = FALSE,
                        print.tau2,
                        print.tau2.ci = FALSE,
                        print.tau = FALSE,
                        print.tau.ci = FALSE,
                        print.Q = FALSE,
                        print.pval.Q,
                        print.Rb = FALSE,
                        print.Rb.ci = FALSE,
                        text.subgroup.nohet = "not applicable",
                        ##
                        LRT = FALSE,
                        ##
                        test.overall = gs("test.overall"),
                        test.overall.fixed =
                          comb.fixed & overall & test.overall,
                        test.overall.random =
                          comb.random & overall & test.overall,
                        label.test.overall.fixed,
                        label.test.overall.random,
                        ##
                        print.stat = TRUE,
                        ##
                        test.subgroup,
                        test.subgroup.fixed,
                        test.subgroup.random,
                        print.Q.subgroup = TRUE,
                        label.test.subgroup.fixed,
                        label.test.subgroup.random,
                        ##
                        test.effect.subgroup,
                        test.effect.subgroup.fixed,
                        test.effect.subgroup.random,
                        label.test.effect.subgroup.fixed,
                        label.test.effect.subgroup.random,
                        ##
                        text.addline1,
                        text.addline2,
                        ##
                        fontsize = 12,
                        fontfamily = NULL,
                        fs.heading = fontsize,
                        fs.fixed,
                        fs.random,
                        fs.predict,
                        fs.fixed.labels,
                        fs.random.labels,
                        fs.predict.labels,
                        fs.study = fontsize,
                        fs.study.labels = fs.study,
                        fs.hetstat,
                        fs.test.overall,
                        fs.test.subgroup,
                        fs.test.effect.subgroup,
                        fs.addline,
                        fs.axis = fontsize,
                        fs.smlab = fontsize,
                        fs.xlab = fontsize,
                        fs.lr = fontsize,
                        ##
                        ff.heading = "bold",
                        ff.fixed,
                        ff.random,
                        ff.predict,
                        ff.fixed.labels,
                        ff.random.labels,
                        ff.predict.labels,
                        ff.study = "plain",
                        ff.study.labels = ff.study,
                        ff.hetstat,
                        ff.test.overall,
                        ff.test.subgroup,
                        ff.test.effect.subgroup,
                        ff.addline,
                        ff.axis = "plain",
                        ff.smlab = "bold",
                        ff.xlab = "plain",
                        ff.lr = "plain",
                        ##
                        squaresize = 0.8 / spacing,
                        ##
                        plotwidth = if (layout == "JAMA") "8cm" else "6cm",
                        colgap = "2mm",
                        colgap.left = colgap,
                        colgap.right = colgap,
                        colgap.studlab = colgap.left,
                        colgap.forest = colgap,
                        colgap.forest.left = colgap.forest,
                        colgap.forest.right = colgap.forest,
                        ##
                        calcwidth.pooled =
                          (comb.fixed | comb.random) &
                          (overall | !is.null(x$byvar)),
                        calcwidth.fixed = calcwidth.pooled,
                        calcwidth.random = calcwidth.pooled,
                        calcwidth.predict = FALSE,
                        calcwidth.hetstat = FALSE,
                        calcwidth.tests  = FALSE,
                        calcwidth.subgroup = FALSE,
                        ##
                        just = if (layout == "JAMA") "left" else "right",
                        just.studlab = "left",
                        just.addcols = "center",
                        just.addcols.left = just.addcols,
                        just.addcols.right = just.addcols,
                        ##
                        spacing = 1,
                        addrow,
                        addrow.overall,
                        addrow.subgroups,
                        ##
                        new = TRUE,
                        ##
                        backtransf = x$backtransf,
                        digits = gs("digits.forest"),
                        digits.se = gs("digits.se"),
                        digits.stat = gs("digits.stat"),
                        digits.pval = max(gs("digits.pval") - 2, 2),
                        digits.pval.Q = max(gs("digits.pval.Q") - 2, 2),
                        digits.Q = gs("digits.Q"),
                        digits.tau2 = gs("digits.tau2"),
                        digits.tau = gs("digits.tau"),
                        digits.I2 = max(gs("digits.I2") - 1, 0),
                        digits.weight = gs("digits.weight"),
                        ##
                        digits.mean = digits,
                        digits.sd = digits.se,
                        digits.cor = digits,
                        digits.time = digits,
                        ##
                        digits.addcols = digits,
                        digits.addcols.right = digits.addcols,
                        digits.addcols.left = digits.addcols,
                        ##
                        scientific.pval = gs("scientific.pval"),
                        big.mark = gs("big.mark"),
                        zero.pval =
                          if (layout == "JAMA") FALSE else gs("zero.pval"),
                        JAMA.pval =
                          if (layout == "JAMA") TRUE else gs("JAMA.pval"),
                        ##
                        col.i = col.study,
                        weight = weight.study,
                        digits.zval = digits.stat,
                        print.zval = print.stat,
                        ...) {
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(x, "meta")
  x.name <- deparse(substitute(x))
  x <- updateversion(x)
  ##
  K.all <- length(x$TE)
  ##
  sm <- x$sm
  ##
  metabin <- inherits(x, "metabin")
  metacont <- inherits(x, "metacont")
  metacor <- inherits(x, "metacor")
  metagen <- inherits(x, "metagen")
  metainc <- inherits(x, "metainc")
  metaprop <- inherits(x, "metaprop")
  metarate <- inherits(x, "metarate")
  metabind <- inherits(x, "is.metabind")
  ##
  metainf.metacum <- inherits(x, "metainf") | inherits(x, "metacum")
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  mf <- match.call()
  error <- try(sortvar <- eval(mf[[match("sortvar", names(mf))]],
                               as.data.frame(x, stringsAsFactors = FALSE),
                               enclos = sys.frame(sys.parent())),
               silent = TRUE)
  if (class(error) == "try-error") {
    xd <- x$data
    sortvar <- eval(mf[[match("sortvar", names(mf))]],
                    xd, enclos = NULL)
    if (isCol(x$data, ".subset"))
      sortvar <- sortvar[x$data$.subset]
  }
  ##
  if (!is.null(sortvar) & metainf.metacum)
    warning("Argument 'sortvar' ignored for objects ",
            "created with metacum() or metainf().")
  sort <- !is.null(sortvar) & !metainf.metacum
  if (sort && (length(sortvar) != K.all))
    stop("Number of studies in object 'x' and argument 'sortvar' ",
         "have different length.")
  if (!sort)
    sortvar <- 1:K.all
  ##
  slab <- TRUE
  missing.studlab <- missing(studlab)
  ##
  if (!missing.studlab) {
    error <- try(studlab <- eval(mf[[match("studlab", names(mf))]],
                                 as.data.frame(x, stringsAsFactors = FALSE),
                                 enclos = sys.frame(sys.parent())),
                 silent = TRUE)
    if (class(error) == "try-error") {
      xd <- x$data
      studlab <- eval(mf[[match("studlab", names(mf))]],
                      xd, enclos = NULL)
      if (isCol(x$data, ".subset"))
        studlab <- studlab[x$data$.subset]
    }
  }
  ##
  if (length(studlab) == 1 & is.logical(studlab)) {
    if (studlab == FALSE) {
      studlab <- rep("", K.all)
      slab <- FALSE
    }
    else studlab <- x$studlab
  }
  ##
  if (length(studlab) != (K.all - 2 * (metainf.metacum & !missing.studlab)))
    stop("Number of studies in object 'x' and argument 'studlab' have ",
         "different length.")
  ##
  chklogical(comb.fixed)
  overall <- replaceNULL(overall, comb.fixed | comb.random)
  chklogical(comb.random)
  chklogical(overall)
  ##
  if (!is.null(lty.fixed))
    chknumeric(lty.fixed)
  if (!is.null(lty.random))
    chknumeric(lty.random)
  chkcolor(col.fixed)
  chkcolor(col.random)
  chklogical(prediction)
  chklogical(print.subgroup.labels)
  if (!is.null(print.byvar))
    chklogical(print.byvar)
  if (!is.null(byseparator))
    chkchar(byseparator)
  chklogical(bysort)
  chklogical(pooled.totals)
  chklogical(pooled.events)
  chklogical(pooled.times)
  chklogical(study.results)
  ## chknumeric(xlab.pos) ??
  ## chknumeric(smlab.pos) ??
  chklogical(allstudies)
  ##
  chklogical(backtransf)
  ##
  if (!is.null(pscale))
    chknumeric(pscale, length = 1)
  else
    pscale <- 1
  if (!is.null(irscale))
    chknumeric(irscale, length = 1)
  else
    irscale <- 1
  if (!is.null(irunit) && !is.na(irunit))
    chkchar(irunit)
  ##
  chknumeric(ref, length = 1)
  chknumeric(lower.equi, length = 1)
  chknumeric(upper.equi, length = 1)
  if (!is.na(lower.equi) && !is.na(upper.equi) && lower.equi > upper.equi)
    stop("Value for 'lower.equi' must be smaller than 'upper.equi'.")
  chknumeric(lty.equi)
  chkcolor(col.equi)
  ##
  layout <- setchar(layout, c("meta", "RevMan5", "JAMA", "subgroup"))
  if (layout == "subgroup" & is.null(x$byvar)) {
    warning("Argument 'layout' set to \"meta\" (default) as ",
            "no subgroup analysis was conducted.")
    layout <- "meta"
  }
  if (layout == "subgroup") {
    if (missing(type.subgroup))
      type.subgroup <- "square"
    if (missing(type.subgroup.fixed))
      type.subgroup.fixed <- "square"
    if (missing(type.subgroup.random))
      type.subgroup.random <- "square"
    ##
    if (missing(pooled.totals))
      pooled.totals <- FALSE
  }
  revman5 <- layout == "RevMan5"
  jama <- layout == "JAMA"
  revman5.jama <- revman5 | jama
  ##
  type.study <- setchar(type.study, c("square", "diamond", "predict"))
  type.fixed <- setchar(type.fixed, c("square", "diamond", "predict"))
  type.random <- setchar(type.random, c("square", "diamond", "predict"))
  type.subgroup <- setchar(type.subgroup, c("square", "diamond", "predict"))
  type.subgroup.fixed <- setchar(type.subgroup.fixed,
                                 c("square", "diamond", "predict"))
  type.subgroup.random <- setchar(type.subgroup.random,
                                  c("square", "diamond", "predict"))
  ##
  if (missing(weight.subgroup))
    weight.subgroup <- ifelse(type.subgroup == "square", "weight", "same")
  weight.subgroup <- setchar(weight.subgroup, c("weight", "same"))
  ##
  chklogical(bottom.lr)
  chkchar(lab.NA)
  chkchar(lab.NA.effect)
  chkchar(lab.NA.weight)
  if (!is.null(at))
    chknumeric(at)
  chkcolor(col.diamond)
  chkcolor(col.diamond.fixed)
  chkcolor(col.diamond.random)
  chkcolor(col.diamond.lines)
  chkcolor(col.diamond.lines.fixed)
  chkcolor(col.diamond.lines.random)
  chkcolor(col.inside.fixed)
  chkcolor(col.inside.random)
  chkcolor(col.predict)
  chkcolor(col.predict.lines)
  ##
  missing.hetstat <- missing(hetstat)
  missing.overall.hetstat <- missing(overall.hetstat)
  ##
  overall.hetstat <- replaceNULL(overall.hetstat, TRUE)
  chklogical(overall.hetstat)
  ##
  if (missing(print.I2))
    if (is.character(hetstat) || hetstat || overall.hetstat)
      print.I2 <- TRUE
    else
      print.I2 <- FALSE
  else
    chklogical(print.I2)
  ##
  chklogical(print.I2.ci)
  ##
  if (missing(print.tau2))
    if (is.character(hetstat) || hetstat || overall.hetstat)
      print.tau2 <- TRUE
    else
      print.tau2 <- FALSE
  else
    chklogical(print.tau2)
  ##
  chklogical(print.tau2.ci)
  chklogical(print.tau)
  chklogical(print.tau.ci)
  print.tau2.tau <- print.tau2 | print.tau
  if (print.tau2 & print.tau)
    print.tau2 <- FALSE
  if (print.tau2.ci & print.tau.ci)
    print.tau2.ci <- FALSE
  ##
  chklogical(print.Q)
  ##
  if (missing(print.pval.Q))
    if (is.character(hetstat) || hetstat || overall.hetstat)
      print.pval.Q <- TRUE
    else
      print.pval.Q <- FALSE
  else
    chklogical(print.pval.Q)
  ##
  chklogical(print.Rb)
  chklogical(print.Rb.ci)
  if (!is.logical(text.subgroup.nohet))
    chkchar(text.subgroup.nohet)
  else if (text.subgroup.nohet)
    text.subgroup.nohet <- "not applicable"
  chklogical(print.stat)
  ##
  hetstat.pooled <- ""
  if (is.character(hetstat)) {
    hetstat.pooled <- setchar(hetstat, c("fixed", "random", "study"))
    if (!metabind & hetstat.pooled == "study") {
      warning("Argument 'hetstat = \"study\"' ",
              "only considered for 'metabind' objects.")
      hetstat <- print.I2 | print.tau2.tau | print.Q | print.pval.Q | print.Rb
    }
  }
  else
    chklogical(hetstat)
  ##
  if (hetstat.pooled == "fixed") {
    comb.fixed <- TRUE
    overall <- TRUE
  }
  if (hetstat.pooled == "random") {
    comb.random <- TRUE
    overall <- TRUE
  }
  ##
  if (missing.overall.hetstat) {
    if (!missing.hetstat)
      if (is.character(hetstat))
        overall.hetstat <- FALSE
      else
        overall.hetstat <- hetstat
  }
  else
    chklogical(overall.hetstat)
  ##
  chklogical(LRT)
  if (LRT & x$method != "GLMM") {
    warning("Likelihood-Ratio test of heterogeneity only ",
            "available for generalized linear mixed models.")
    LRT <- FALSE
  }
  ##
  chkchar(hetlab)
  if (!missing(resid.hetstat))
    chklogical(resid.hetstat)
  else {
    if (overall && (is.character(hetstat) || hetstat) && !LRT &&
        !is.null(x$tau.common) && x$tau.common)
      resid.hetstat <- TRUE
    else
      resid.hetstat <- FALSE
  }
  chkchar(resid.hetlab)
  ##
  chklogical(test.overall.fixed)
  chklogical(test.overall.random)
  if (!missing(test.subgroup.fixed))
    chklogical(test.subgroup.fixed)
  if (!missing(test.subgroup.random))
    chklogical(test.subgroup.random)
  chklogical(print.Q.subgroup)
  if (!missing(test.effect.subgroup.fixed))
    chklogical(test.effect.subgroup.fixed)
  if (!missing(test.effect.subgroup.random))
    chklogical(test.effect.subgroup.random)
  chknumeric(fontsize, length = 1)
  chknumeric(fs.heading, length = 1)
  if (!missing(fs.fixed))
    chknumeric(fs.fixed, length = 1)
  if (!missing(fs.random))
    chknumeric(fs.random, length = 1)
  if (!missing(fs.predict))
    chknumeric(fs.predict, length = 1)
  if (!missing(fs.fixed.labels))
    chknumeric(fs.fixed.labels, length = 1)
  if (!missing(fs.random.labels))
    chknumeric(fs.random.labels, length = 1)
  if (!missing(fs.predict.labels))
    chknumeric(fs.predict.labels, length = 1)
  if (!missing(fs.study))
    chknumeric(fs.study, length = 1)
  if (!missing(fs.study.labels))
    chknumeric(fs.study.labels, length = 1)
  if (!missing(fs.hetstat))
    chknumeric(fs.hetstat, length = 1)
  if (!missing(fs.test.overall))
    chknumeric(fs.test.overall, length = 1)
  if (!missing(fs.test.subgroup))
    chknumeric(fs.test.subgroup, length = 1)
  if (!missing(fs.test.effect.subgroup))
    chknumeric(fs.test.effect.subgroup, length = 1)
  if (!missing(fs.addline))
    chknumeric(fs.addline, length = 1)
  chknumeric(fs.axis, length = 1)
  chknumeric(fs.smlab, length = 1)
  chknumeric(fs.xlab, length = 1)
  chknumeric(fs.lr, length = 1)
  chknumeric(squaresize, length = 1)
  chklogical(calcwidth.pooled)
  chklogical(calcwidth.fixed)
  chklogical(calcwidth.random)
  chklogical(calcwidth.predict)
  chklogical(calcwidth.hetstat)
  chklogical(calcwidth.tests)
  chklogical(calcwidth.subgroup)
  just.cols <- setchar(just, c("right", "center", "left"))
  just.studlab <- setchar(just.studlab, c("right", "center", "left"))
  just.addcols <- setchar(just.addcols, c("right", "center", "left"))
  just.addcols.left <- setchar(just.addcols.left, c("right", "center", "left"))
  just.addcols.right <- setchar(just.addcols.right, c("right", "center", "left"))
  ##
  if (missing(weight.study))
    weight.study <- ifelse(comb.random & !comb.fixed, "random", "fixed")
  weight.study <- setchar(weight.study, c("same", "fixed", "random"))
  ##
  chknumeric(spacing, length = 1)
  ##
  missing.text.addline1 <- missing(text.addline1)
  if (!missing.text.addline1)
    chkchar(text.addline1)
  else
    text.addline1 <- ""
  missing.text.addline2 <- missing(text.addline2)
  if (!missing.text.addline2)
    chkchar(text.addline2)
  else
    text.addline2 <- ""
  ##
  ## Check and set additional empty rows in forest plot
  ##
  if (!missing(addrow))
    chklogical(addrow)
  else
    addrow <- !revman5.jama
  if (!missing(addrow.overall))
    chklogical(addrow.overall)
  else
    addrow.overall <- !jama & overall & (comb.fixed | comb.random | prediction)
  if (!missing(addrow.subgroups))
    chklogical(addrow.subgroups)
  else
    addrow.subgroups <- !jama
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  chknumeric(digits.pval.Q, min = 1, length = 1)
  chknumeric(digits.Q, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  if (!missing(digits.mean))
    chknumeric(digits.mean, min = 0, length = 1)
  if (!missing(digits.sd))
    chknumeric(digits.sd, min = 0, length = 1)
  if (!missing(digits.cor))
    chknumeric(digits.cor, min = 0, length = 1)
  missing.digits.time <- missing(digits.time)
  if (!missing.digits.time)
    chknumeric(digits.time, min = 0, length = 1)
  missing.addcols.left <-
    missing(digits.addcols) & missing(digits.addcols.left)
  missing.addcols.right <-
    missing(digits.addcols) & missing(digits.addcols.right)
  if (!missing(digits.addcols))
    chknumeric(digits.addcols, min = 0)
  if (!missing(digits.addcols.right))
    chknumeric(digits.addcols.right, min = 0)
  if (!missing(digits.addcols.left))
    chknumeric(digits.addcols.left, min = 0)
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  ##
  cl <- paste0("update.meta() or ", class(x)[1], "()")
  addargs <- names(list(...))
  ##
  fun <- "forest.meta"
  ##
  warnarg("byvar", addargs, fun, cl)
  warnarg("level", addargs, fun, cl)
  warnarg("level.comb", addargs, fun, cl)
  warnarg("level.predict", addargs, fun, cl)
  ##
  ## Check for deprecated argument 'col.i'
  ##
  if (!missing(col.i))
    if (!missing(col.study))
      warning("Deprecated argument 'col.i' ignored as ",
              "argument 'col.study' is also provided.")
    else {
      warning("Deprecated argument 'col.i' has been replaced by ",
              "argument 'col.study'.")
      col.study <- col.i
    }
  ##
  ## Check for deprecated argument 'weight'
  ##
  if (!missing(weight))
    if (!missing(weight.study))
      warning("Deprecated argument 'weight' ignored as ",
              "argument 'weight.study' is also provided.")
    else {
      warning("Deprecated argument 'weight' has been replaced ",
              "by argument 'weight.study'.")
      weight.study <- weight
      weight.study <- setchar(weight.study, c("same", "fixed", "random"))
    }
  ##
  ## Check for deprecated argument 'digits.zval'
  ##
  if (!missing(digits.zval))
    if (!missing(digits.stat))
      warning("Deprecated argument 'digits.zval' ignored as ",
              "argument 'digits.stat' is also provided.")
    else {
      warning("Deprecated argument 'digits.zval' has been replaced by ",
              "argument 'digits.stat'.")
      digits.stat <- digits.zval
      chknumeric(digits.stat, min = 0, length = 1)
    }
  ##
  ## Check for deprecated argument 'print.zval'
  ##
  if (!missing(print.zval))
    if (!missing(print.stat))
      warning("Deprecated argument 'print.zval' ignored as ",
              "argument 'print.stat' is also provided.")
    else {
      warning("Deprecated argument 'print.zval' has been replaced by ",
              "argument 'print.stat'.")
      print.stat <- print.zval
      chklogical(print.stat)
    }
  ##
  ## Check for other deprecated arguments in '...'
  ##
  args  <- list(...)
  ## Check whether first argument is a list. In this case only use
  ## this list as input.
  if (length(args) > 0 && is.list(args[[1]]))
    args <- args[[1]]
  ##
  additional.arguments <- names(args)
  ##
  if (length(additional.arguments) > 0) {
    if ("labels" %in% additional.arguments)
      if (!missing(label))
        warning("Argument 'labels' ignored as both arguments ",
                "'label' and 'labels' are provided.")
      else
        label <- args[["labels"]]
    ##
    if (!is.na(charmatch("col.i.i", additional.arguments)))
      if (!missing(col.inside))
        warning("Deprecated argument 'col.i.inside.square' ignored as ",
                "argument 'col.inside' is also provided.")
      else {
        warning("Deprecated argument 'col.i.inside.square' has been ",
                "replaced by argument 'col.inside'.")
        col.inside <- args[[charmatch("col.i.i", additional.arguments)]]
      }
    ##
    if (!is.na(charmatch("col.diamond.f", additional.arguments)))
      if (!missing(col.diamond.lines.fixed))
        warning("Deprecated argument 'col.diamond.fixed.lines' ignored as ",
                "argument 'col.diamond.lines.fixed' is also provided.")
      else {
        warning("Deprecated argument 'col.diamond.fixed.lines' has been ",
                "replaced by argument 'col.diamond.lines.fixed'.")
        col.diamond.lines.fixed <- args[[charmatch("col.diamond.f",
                                                   additional.arguments)]]
      }
    ##
    if (!is.na(charmatch("col.diamond.r", additional.arguments)))
      if (!missing(col.diamond.lines.random))
        warning("Deprecated argument 'col.diamond.random.lines' ignored as ",
                "argument 'col.diamond.lines.random' is also provided.")
      else {
        warning("Deprecated argument 'col.diamond.random.lines' has been ",
                "replaced by argument 'col.diamond.lines.random'.")
        col.diamond.lines.random <- args[[charmatch("col.diamond.r",
                                                    additional.arguments)]]
      }
    ##
    if (!is.na(charmatch("adds", additional.arguments)))
      if (!missing(addrow))
        warning("Deprecated argument 'addspace' ignored as ",
                "argument 'addrow' is also provided.")
      else {
        warning("Deprecated argument 'addspace' has been replaced by ",
                "argument 'addrow'.")
        addrow <- args[[charmatch("adds", additional.arguments)]]
      }
    ##
    if (!is.na(charmatch("text.I2", additional.arguments)))
      warning("Argument 'text.I2' has been removed.")
    if (!is.na(charmatch("text.tau2", additional.arguments)))
      warning("Argument 'text.tau2' has been removed.")
  }
  ##
  ## Additional assignments
  ##
  if (jama) {
    if (missing(ff.fixed))
      ff.fixed <- "plain"
    if (missing(ff.random))
      ff.random <- ff.fixed
    if (missing(ff.predict))
      ff.predict <- ff.fixed
    if (missing(ff.fixed.labels))
      ff.fixed.labels <- ff.fixed
    if (missing(ff.random.labels))
      ff.random.labels <- ff.random
    if (missing(ff.predict.labels))
      ff.predict.labels <- ff.predict
    ##
    if (missing(fs.fixed))
      fs.fixed <- fontsize
    if (missing(fs.random))
      fs.random <- fs.fixed
    if (missing(fs.predict))
      fs.predict <- fs.fixed
    if (missing(fs.fixed.labels))
      fs.fixed.labels <- fs.fixed
    if (missing(fs.random.labels))
      fs.random.labels <- fs.random
    if (missing(fs.predict.labels))
      fs.predict.labels <- fs.predict
  }
  else {
    if (missing(ff.fixed))
      ff.fixed <- "bold"
    if (missing(ff.random))
      ff.random <- ff.fixed
    if (missing(ff.predict))
      ff.predict <- ff.fixed
    if (missing(ff.fixed.labels))
      ff.fixed.labels <- ff.fixed
    if (missing(ff.random.labels))
      ff.random.labels <- ff.random
    if (missing(ff.predict.labels))
      ff.predict.labels <- ff.predict
    ##
    if (missing(fs.fixed))
      fs.fixed <- fontsize
    if (missing(fs.random))
      fs.random <- fs.fixed
    if (missing(fs.predict))
      fs.predict <- fs.fixed
    if (missing(fs.fixed.labels))
      fs.fixed.labels <- fs.fixed
    if (missing(fs.random.labels))
      fs.random.labels <- fs.random
    if (missing(fs.predict.labels))
      fs.predict.labels <- fs.predict
  }
  hetseparator <- " = "
  ##
  if (revman5) {
    if (missing(ff.hetstat))
      ff.hetstat <- "plain"
    if (missing(ff.test.overall))
      ff.test.overall <- ff.hetstat
    if (missing(ff.test.subgroup))
      ff.test.subgroup <- ff.hetstat
    if (missing(ff.test.effect.subgroup))
      ff.test.effect.subgroup <- ff.hetstat
    if (missing(ff.addline))
      ff.addline <- ff.hetstat
    ##
    if (missing(fs.hetstat))
      fs.hetstat <- fontsize - 1
    if (missing(fs.test.overall))
      fs.test.overall <- fs.hetstat
    if (missing(fs.test.subgroup))
      fs.test.subgroup <- fs.hetstat
    if (missing(fs.test.effect.subgroup))
      fs.test.effect.subgroup <- fs.hetstat
    if (missing(fs.addline))
      fs.addline <- fs.hetstat
  }
  else if (jama) {
    if (missing(ff.hetstat))
      ff.hetstat <- "plain"
    if (missing(ff.test.overall))
      ff.test.overall <- ff.hetstat
    if (missing(ff.test.subgroup))
      ff.test.subgroup <- ff.hetstat
    if (missing(ff.test.effect.subgroup))
      ff.test.effect.subgroup <- ff.hetstat
    if (missing(ff.addline))
      ff.addline <- ff.hetstat
    ##
    if (missing(fs.hetstat))
      fs.hetstat <- fontsize - 1
    if (missing(fs.test.overall))
      fs.test.overall <- fs.hetstat
    if (missing(fs.test.subgroup))
      fs.test.subgroup <- fs.hetstat
    if (missing(fs.test.effect.subgroup))
      fs.test.effect.subgroup <- fs.hetstat
    if (missing(fs.addline))
      fs.addline <- fs.hetstat
  }
  else {
    if (missing(ff.hetstat))
      ff.hetstat <- "plain"
    if (missing(ff.test.overall))
      ff.test.overall <- ff.hetstat
    if (missing(ff.test.subgroup))
      ff.test.subgroup <- ff.hetstat
    if (missing(ff.test.effect.subgroup))
      ff.test.effect.subgroup <- ff.hetstat
    if (missing(ff.addline))
      ff.addline <- ff.hetstat
    ##
    if (missing(fs.hetstat))
      fs.hetstat <- fontsize - 1
    if (missing(fs.test.overall))
      fs.test.overall <- fs.hetstat
    if (missing(fs.test.subgroup))
      fs.test.subgroup <- fs.hetstat
    if (missing(fs.test.effect.subgroup))
      fs.test.effect.subgroup <- fs.hetstat
    if (missing(fs.addline))
      fs.addline <- fs.hetstat
  }
  
  
  ##
  ##
  ## (3) Check length of variables
  ##
  ##
  if (length(col.study) == 1)
    col.study <- rep(col.study, K.all)
  else
    chklength(col.study, K.all, fun)
  ##
  if (length(col.inside) == 1)
    col.inside <- rep(col.inside, K.all)
  else
    chklength(col.inside, K.all, fun)
  ##
  miss.col.square <- missing(col.square)
  miss.col.square.lines <- missing(col.square.lines)
  if (length(col.square) == 1)
    col.square <- rep(col.square, K.all)
  else
    chklength(col.square, K.all, fun)
  ##
  if (length(col.square.lines) == 1)
    col.square.lines <- rep(col.square.lines, K.all)
  else
    chklength(col.square.lines, K.all, fun)
  
  
  ##
  ##
  ## (4) Some assignments and additional checks
  ##
  ##
  prediction <- prediction & x$k >= 3
  ##
  level <- x$level
  level.comb <- x$level.comb
  level.predict <- x$level.predict
  ##
  chklevel(level)
  chklevel(level.comb)
  if (prediction)
    chklevel(level.predict)
  ##
  byvar <- x$byvar
  by <- !is.null(byvar)
  ##
  if (!by) {
    test.subgroup.fixed  <- FALSE
    test.subgroup.random <- FALSE
    test.effect.subgroup.fixed <- FALSE
    test.effect.subgroup.random <- FALSE
  }
  else {
    if (!missing(test.subgroup)) {
      if (missing(test.subgroup.fixed))
        test.subgroup.fixed <- comb.fixed & test.subgroup
      ##
      if (missing(test.subgroup.random))
        test.subgroup.random <- comb.random & test.subgroup
    }
    else {
      if (missing(test.subgroup.fixed))
        test.subgroup.fixed <- comb.fixed & gs("test.subgroup")
      ##
      if (missing(test.subgroup.random))
        test.subgroup.random <- comb.random & gs("test.subgroup")
    }
    ##
    if (!missing(test.effect.subgroup)) {
      if (missing(test.effect.subgroup.fixed))
        test.effect.subgroup.fixed <- comb.fixed & test.effect.subgroup
      ##
      if (missing(test.effect.subgroup.random))
        test.effect.subgroup.random <- comb.random & test.effect.subgroup
    }
    else {
      if (missing(test.effect.subgroup.fixed))
        test.effect.subgroup.fixed <- comb.fixed & gs("test.effect.subgroup")
      ##
      if (missing(test.effect.subgroup.random))
        test.effect.subgroup.random <- comb.random & gs("test.effect.subgroup")
    }
  }
  ##
  fixed.random <-
    (comb.fixed | test.subgroup.fixed | test.effect.subgroup.fixed) &
    (comb.random | test.subgroup.random | test.effect.subgroup.random)
  ##
  if (layout == "subgroup") {
    if (!missing(study.results) & study.results)
      warning("Argument 'study.results' set to FALSE as ",
              "argument 'layout' is \"subgroup\".")
    study.results <- FALSE
  }
  ##
  missing.text.fixed <- missing(text.fixed)
  if (missing.text.fixed | is.null(text.fixed)) {
    if (study.results & (x$level != x$level.comb | revman5)) {
      if (revman5.jama)
        text.fixed <- paste0("Total (",
                             if (fixed.random)
                               "fixed effect, ",
                             round(x$level.comb * 100), "% CI)")
      else if (!is.null(text.fixed))
        text.fixed <- paste0(text.fixed, " (",
                             round(x$level.comb * 100), "%-CI)")
      else
        text.fixed <- paste0("Fixed effect model (",
                             round(x$level.comb * 100), "%-CI)")
    }
    else {
      if (revman5.jama) {
        text.fixed <- "Total"
        if (fixed.random)
          text.fixed <- paste(text.fixed, "(fixed effect)")
      }
      else if (is.null(text.fixed))
        text.fixed <- "Fixed effect model"
    }
  }
  ##
  missing.text.random <- missing(text.random)
  if (missing.text.random | is.null(text.random)) {
    if (study.results & (x$level != x$level.comb | revman5)) {
      if (revman5.jama)
        text.random <- paste0("Total (",
                              if (fixed.random)
                                "random effects, ",
                              round(x$level.comb * 100), "% CI)")
      else if (!is.null(text.random))
        text.random <- paste0(text.random, " (",
                              round(x$level.comb * 100), "%-CI)")
      else
        text.random <- paste0("Random effects model (",
                              round(x$level.comb * 100), "%-CI)")
    }
    else {
      if (revman5.jama) {
        text.random <- "Total"
        if (fixed.random)
          text.random <- paste(text.random, "(random effects)")
      }
      else if (is.null(text.random))
        text.random <- "Random effects model"
    }
  }
  ##
  missing.text.predict <- missing(text.predict)
  if (missing.text.predict | is.null(text.predict)) {
    if (is.null(text.predict))
      text.predict <- "Prediction interval"
    if (!(length(x$level.predict) == 0) &&
        (study.results & (x$level != x$level.predict |
                          x$level.comb != x$level.predict)))
      text.predict <- paste0(text.predict, " (",
                             round(x$level.predict * 100), "%-PI)")
  }
  ##
  if (metainf.metacum) {
    overall.hetstat <- FALSE
    test.overall.fixed   <- FALSE
    test.overall.random  <- FALSE
    resid.hetstat <- FALSE
    test.subgroup.fixed  <- FALSE
    test.subgroup.random <- FALSE
    ##
    hetstat <- FALSE
    prediction <- FALSE
    test.effect.subgroup.fixed  <- FALSE
    test.effect.subgroup.random <- FALSE
  }
  ##
  if (is.null(x$null.effect) || is.na(x$null.effect)) {
    test.overall.fixed  <- FALSE
    test.overall.random <- FALSE
    test.effect.subgroup.fixed  <- FALSE
    test.effect.subgroup.random <- FALSE
  }
  ##
  if (!overall) {
    if (test.overall)
      test.overall <- FALSE
    if (test.overall.fixed)
      test.overall.fixed <- FALSE
    if (test.overall.random)
      test.overall.random <- FALSE
  }
  ##
  if (is.logical(leftcols)) {
    if (!leftcols) {
      text.fixed <- ""
      text.random <- ""
      text.predict <- ""
      leftcols <- "studlab"
      studlab <- rep("", K.all)
      slab <- FALSE
      hetstat <- FALSE
      overall.hetstat <- FALSE
      resid.hetstat <- FALSE
      colgap.left <- grid::unit(0, "mm")
      colgap.studlab <- grid::unit(0, "mm")
      colgap.forest.left <- grid::unit(0, "mm")
    }
    else
      leftcols <- NULL
  }
  ##
  notmiss.xlim <- !missing(xlim)
  ##
  if (just.studlab == "left")
    xpos.s <- 0
  else if (just.studlab == "center")
    xpos.s <- 0.5
  else if (just.studlab == "right")
    xpos.s <- 1
  ##
  if (just.cols == "left")
    xpos.c <- 0
  else if (just.cols == "center")
    xpos.c <- 0.5
  else if (just.cols == "right")
    xpos.c <- 1
  ##
  log.xaxis <- FALSE
  ##
  if (missing(ref) && (is.prop(sm) | is.rate(sm)))
    ref <- NA
  ##
  if (backtransf & is.relative.effect(sm)) {
    ref <- log(ref)
    lower.equi <- log(lower.equi)
    upper.equi <- log(upper.equi)
    log.xaxis <- TRUE
  }
  ##
  if (!backtransf & !missing(pscale) & pscale != 1 & !is.untransformed(sm)) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  ##
  if (!backtransf & !missing(irscale) & irscale != 1 & !is.untransformed(sm)) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  ##
  if (is.null(xlab))
    xlab <- xlab(sm, backtransf, newline = revman5.jama, revman5 = revman5,
                 big.mark = big.mark)
  ##
  smlab.null <- is.null(smlab)
  if (smlab.null)
    if (is.rate(sm))
      smlab <- xlab(sm, backtransf, irscale = irscale, irunit = irunit,
                    newline = !revman5.jama, revman5 = revman5,
                    big.mark = big.mark)
    else
      smlab <- xlab(sm, backtransf, pscale = pscale,
                    newline = !revman5.jama, revman5 = revman5,
                    big.mark = big.mark)
  ##
  if (is.null(label.right))
    label.right <- ""
  if (is.null(label.left))
    label.left <- ""
  ##
  print.label <- label.left != "" | label.right != "" & !is.na(ref)
  if (print.label & !bottom.lr) {
    if (!smlab.null)
      warning("Argument 'smlab' ignored as argument 'bottom.lr' is FALSE.")
    smlab <- ""
  }
  ##
  if (!by) {
    addrow.subgroups <- FALSE
    if (!missing(resid.hetstat) && resid.hetstat)
      warning("Information on residual heterogeneity only added to ",
              "forest plot of meta-analysis with subgroups",
              call. = FALSE)
    resid.hetstat <- FALSE
  }
  ##
  plotwidth <- setunit(plotwidth)
  colgap <- setunit(colgap)
  colgap.left <- setunit(colgap.left)
  colgap.right <- setunit(colgap.right)
  colgap.studlab <- setunit(colgap.studlab)
  colgap.forest <- setunit(colgap.forest)
  colgap.forest.left <- setunit(colgap.forest.left)
  colgap.forest.right <- setunit(colgap.forest.right)
  ##
  if (missing(label.test.overall.fixed))
    label.test.overall.fixed <-
      paste0("Test for overall effect",
             if (fixed.random) " (fixed effect)",
             ": ")
  if (missing(label.test.overall.random))
    label.test.overall.random <-
      paste0("Test for overall effect",
             if (fixed.random) " (random effects)",
             ": ")
  ##
  if (missing(label.test.subgroup.fixed))
    label.test.subgroup.fixed <-
      paste0("Test for subgroup differences",
             if (fixed.random) " (fixed effect)",
             ": ")
  if (missing(label.test.subgroup.random))
    label.test.subgroup.random <-
      paste0("Test for subgroup differences",
             if (fixed.random) " (random effects)",
             ": ")
  ##
  if (missing(label.test.effect.subgroup.fixed))
    label.test.effect.subgroup.fixed <-
      paste0(if (revman5.jama)
               "Test for overall effect"
             else
               "Test for effect in subgroup",
             if (fixed.random) " (fixed effect)",
             ": ")
  if (missing(label.test.effect.subgroup.random))
    label.test.effect.subgroup.random <-
      paste0(if (revman5.jama)
               "Test for overall effect"
             else
               "Test for effect in subgroup",
             if (fixed.random) " (random effects)",
             ": ")
  ##
  fs.head <- fs.heading
  ff.head <- ff.heading
  ##
  just.c <- just.cols
  just.s <- just.studlab
  
  
  ##
  ##
  ## (5) Determine columns on left and right side of forest plot
  ##
  ##
  ## Determine whether to print columns on right and / or left side
  ## of forest plot
  ##
  rsel <- !(is.logical(rightcols) && length(rightcols) == 1 && !rightcols)
  if (revman5.jama)
    rsel <- FALSE
  ##
  if (!rsel)
    rightcols <- NULL
  ##
  ## Check for duplicate columns
  ##
  if (length(c(rightcols, leftcols)) > 0 &&
      any(duplicated(c(rightcols, leftcols))))
    stop("Duplicate entries in 'leftcols' and 'rightcols'.")
  ##
  ## Predefined columns and labels
  ##
  sm.lab <- sm
  ##
  if (backtransf) {
    if (sm == "ZCOR")
      sm.lab <- "COR"
    else if (is.prop(sm)) {
      if (pscale == 1)
        sm.lab <- "Proportion"
      else
        sm.lab <- "Events"
    }
    else if (is.rate(sm)) {
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
  else 
    if (is.relative.effect(sm))
      sm.lab <- paste0("log", sm)
  ##
  colnames <- c("studlab", "TE", "seTE",
                "n.e", "n.c",
                "event.e", "event.c",
                "mean.e", "mean.c",
                "sd.e", "sd.c",
                "cor",
                "time.e", "time.c",
                "effect", "ci",
                "effect.ci",
                "w.fixed", "w.random")
  ##
  sel.studlab <- pmatch(layout, c("meta", "RevMan5", "JAMA", "subgroup"))
  lab.studlab <- c("Study", "Study", "Source", "Subgroup")[sel.studlab]
  if (revman5 & by)
    lab.studlab <- c("Study or\nSubgroup")
  ##
  if (revman5.jama)
    cisep <- " "
  else
    cisep <- "-"
  ##
  if (study.results)
    ci.lab <- paste0(100 * level, "%", cisep, "CI")
  else
    ci.lab <- paste0(100 * level.comb, "%", cisep, "CI")
  ##
  if (jama) {
    if (missing(ff.lr))
      ff.lr <- "bold"
    if (xlab == "")
      xlab <- paste0(smlab, " (", ci.lab, ")")
    if (miss.col.square)
      col.square <- rep("darkblue", K.all)
    if (miss.col.square.lines)
      col.square.lines <- rep("darkblue", K.all)
    if (missing(col.diamond.fixed))
      col.diamond.fixed <- "lightblue"
    if (missing(col.diamond.random))
      col.diamond.random <- "lightblue"
    ##
    smlab <- ""
    bottom.lr <- FALSE
  }
  else {
    if (revman5) {
      if (miss.col.square) {
        if (metacont)
          col.square <- rep("green", K.all)
        else if (metabin)
          col.square <- rep("blue", K.all)
        else
          col.square <- rep("red", K.all)
      }
      if (miss.col.square.lines) {
        if (metacont)
          col.square.lines <- rep("green", K.all)
        else if (metabin)
          col.square.lines <- rep("darkblue", K.all)
        else
          col.square.lines <- rep("red", K.all)
      }
      if (missing(col.diamond.fixed))
        col.diamond.fixed <- "black"
      if (missing(col.diamond.random))
        col.diamond.random <- "black"
      ##
      sel.method <- pmatch(x$method, c("Inverse", "MH", "Peto", "GLMM"))
      lab.method <- c("IV", "MH", "Peto", "GLMM")[sel.method]
      ##
      if (fixed.random)
        lab.model <- "Fixed + Random, "
      else if (comb.fixed)
        lab.model <- "Fixed, "
      else if (comb.random)
        lab.model <- "Random, "
      else
        lab.model <- ""
      ##
      if (smlab.null)
        smlab <- paste0(smlab, "\n", lab.method, ", ", lab.model, ci.lab)
    }
  }
  ##
  if (jama | gs("CIbracket") == "(")
    ci.lab.bracket <- paste0("(", ci.lab, ")")
  else if (gs("CIbracket") == "[")
    ci.lab.bracket <- paste0("[", ci.lab, "]")
  else if (gs("CIbracket") == "{")
    ci.lab.bracket <- paste0("{", ci.lab, "}")
  else if (gs("CIbracket") == "")
    ci.lab.bracket <- ci.lab
  ##
  if (!fixed.random) {
    text.w.fixed <- "Weight"
    text.w.random <- "Weight"
  }
  else {
    if (is.null(text.w.fixed))
      text.w.fixed <- paste0("Weight\n(", gs("text.w.fixed"), ")")
    else
      text.w.fixed <- paste0("Weight\n(", text.w.fixed, ")")
    ##
    if (is.null(text.w.random))
      text.w.random <- paste0("Weight\n(", gs("text.w.random"), ")")
    else
      text.w.random <- paste0("Weight\n(", text.w.random, ")")
  }
  ##
  labnames <- c(lab.studlab,
                "TE", if (revman5) "SE" else "seTE",
                "Total", "Total", "Events", "Events",
                "Mean", "Mean", "SD", "SD",
                "Cor",
                "Time", "Time",
                sm.lab,
                ci.lab,
                if (revman5 & smlab.null) smlab else paste(sm.lab, ci.lab.bracket),
                text.w.fixed,
                text.w.random)
  ##
  ## If any of the following list elements is NULL, these 'special'
  ## variable names are searched for in original data set (i.e., list
  ## element x$data)
  ##
  colnames.notNULL <- colnames
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
  ##
  ## Identify and process columns in addition to columns
  ## defined above in variables 'colnames' and 'labnames'
  ##
  colnames.new <- c(rightcols, leftcols)[!c(rightcols, leftcols) %in% colnames.notNULL]
  ##
  newcols <- length(colnames.new) > 0
  ##
  if (newcols) {
    dataset2 <- as.data.frame(x)
    ##
    if (is.null(x$data))
      dataset1 <- dataset2
    else
      dataset1 <- x$data
    ##
    if (!is.null(x$subset))
      dataset1 <- dataset1[x$subset, ]
    ##
    ## Check whether additional variables are
    ## part of meta-object
    ##
    for (i in colnames.new)
      if (length(dataset1[[i]]) == 0 & length(dataset2[[i]]) == 0)
        stop("Variable '", i, "' not available in '", x.name, "'.")
    ##
    rightcols.new <- rightcols[! rightcols %in% colnames.notNULL]
    leftcols.new  <- leftcols[! leftcols %in% colnames.notNULL]
    ##
    ## Determine label for new columns
    ## 1. Use column name as label if no label is given
    ##    argument right | left | labs
    ## 2. Otherwise use corresponding entry from
    ##    argument right | left | labs
    ##
    if (length(rightcols.new) > 0) {
      pos.rightcols.new <- match(rightcols.new, rightcols)
      ##
      rightlabs.new <- rightcols.new
      for (i in seq(along = rightcols.new)) {
        j <- match(rightcols.new[i], colnames)
        if (!is.na(j))
          rightlabs.new[i] <- labnames[j]
        else if (rightcols.new[i] == "pval")
          rightlabs.new[i] <- "P-value"
      }
      ##
      if (missing(rightlabs))
        rightlabs.new <- rightlabs.new
      else if (length(rightcols.new) == length(rightlabs))
        rightlabs.new <- rightlabs
      else if (max(pos.rightcols.new) <= length(rightlabs))
        rightlabs.new <- rightlabs[pos.rightcols.new]
      else if (max(pos.rightcols.new) > length(rightlabs))
        stop("Too few labels defined for argument 'rightcols'.")
      ##
      if ( (metacor | metaprop) & any(rightcols.new == "n"))
        rightlabs.new[rightlabs.new == "n"] <- "Total"
      ##
      if ( (metarate) & any(rightcols.new == "time"))
        rightlabs.new[rightlabs.new == "time"] <- "Time"
    }
    if (length(leftcols.new) > 0) {
      pos.leftcols.new <- match(leftcols.new, leftcols)
      ##
      leftlabs.new <- leftcols.new
      for (i in seq(along = leftcols.new)) {
        j <- match(leftcols.new[i], colnames)
        if (!is.na(j))
          leftlabs.new[i] <- labnames[j]
        else if (leftcols.new[i] == "pval")
          leftlabs.new[i] <- "P-value"
      }
      ##
      if (missing(leftlabs))
        leftlabs.new <- leftlabs.new
      else if (length(leftcols.new) == length(leftlabs))
        leftlabs.new <- leftlabs
      else if (max(pos.leftcols.new) <= length(leftlabs))
        leftlabs.new <- leftlabs[pos.leftcols.new]
      else if (max(pos.leftcols.new) > length(leftlabs))
        stop("Too few labels defined for argument 'leftcols'.")
      ##
      if ((metacor | metaprop) & any(leftcols.new == "n"))
        leftlabs.new[leftlabs.new == "n"] <- "Total"
      ##
      if ((metarate) & any(leftcols.new == "time"))
        leftlabs.new[leftlabs.new == "time"] <- "Time"
    }
  }
  ##
  ## Default set of columns if argument leftcols and / or
  ## rightcols not specified
  ##
  if (is.null(leftcols)) {
    ##
    leftcols <- "studlab"
    ##
    if (jama) {
      leftcols <- c(leftcols,
                    "effect.ci")
    }
    else {
      if (metabin) {
        if (study.results)
          leftcols <- c(leftcols,
                        "event.e", "n.e",
                        "event.c", "n.c")
        else {
          leftcols <- c(leftcols,
                        if (pooled.events) "event.e",
                        if (pooled.totals) "n.e",
                        if (pooled.events) "event.c",
                        if (pooled.totals) "n.c")
          if (pooled.events & !pooled.totals) {
            if (is.null(lab.e.attach.to.col))
              lab.e.attach.to.col <- "event.e"
            if (is.null(lab.c.attach.to.col))
              lab.c.attach.to.col <- "event.c"
          }
        }
      }
      ##
      if (metacont) {
        if (study.results) {
          if (revman5)
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
          if (is.null(lab.e.attach.to.col))
            lab.e.attach.to.col <- "n.e"
          if (is.null(lab.c.attach.to.col))
            lab.c.attach.to.col <- "n.c"
        }
      }
      ##
      if (metagen & study.results) {
        leftcols <- c(leftcols,
                      "TE", "seTE")
        if (!is.null(x$n.e)) {
          leftcols <- c(leftcols, "n.e")
          if (is.null(lab.e.attach.to.col))
            lab.e.attach.to.col <- "n.e"
        }
        if (!is.null(x$n.c)) {
          leftcols <- c(leftcols, "n.c")
          if (is.null(lab.c.attach.to.col))
            lab.c.attach.to.col <- "n.c"
        }
      }
      ##
      if (metaprop) {
        if (study.results)
          leftcols <- c(leftcols,
                        "event.e", "n.e")
        else {
          leftcols <- c(leftcols,
                        if (pooled.events) "event.e",
                        if (pooled.totals) "n.e")
          if (pooled.events & !pooled.totals) {
            if (is.null(lab.e.attach.to.col))
              lab.e.attach.to.col <- "event.e"
          }
        }
      }
      ##
      if (metarate) {
        if (study.results)
          leftcols <- c(leftcols,
                        "event.e", "time.e")
        else {
          leftcols <- c(leftcols,
                        if (pooled.events) "event.e",
                        if (pooled.times) "time.e")
          if (pooled.events & !pooled.times) {
            if (is.null(lab.e.attach.to.col))
              lab.e.attach.to.col <- "event.e"
          }
        }
      }
      ##
      if (metacor) {
        if (study.results | pooled.totals)
          leftcols <- c(leftcols,
                        "n.e")
      }
      ##
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
            if (is.null(lab.e.attach.to.col))
              lab.e.attach.to.col <- "event.e"
            if (is.null(lab.c.attach.to.col))
              lab.c.attach.to.col <- "event.c"
          }
        }
      }
    }
    ##
    ## Add columns for RevMan 5 layout
    ##
    if (revman5) {
      ##
      if (!metainf.metacum & overall & study.results & !x$method == "GLMM") {
        if (comb.fixed)
          leftcols <- c(leftcols, "w.fixed")
        if (comb.random)
          leftcols <- c(leftcols, "w.random")
      }
      ##
      leftcols <- c(leftcols, "effect.ci")
    }
  }
  ##
  if (is.null(rightcols) & rsel) {
    rightcols <- c("effect", "ci")
    ##
    if (!metainf.metacum & overall & study.results & !x$method == "GLMM") {
      if (comb.fixed)
        rightcols <- c(rightcols, "w.fixed")
      if (comb.random)
        rightcols <- c(rightcols, "w.random")
    }
  }
  
  
  ##
  ##
  ## (6) Select data for forest plot
  ##
  ##
  if (metacor) {
    x$n.e <- x$n
  }
  ##
  if (metaprop) {
    x$event.e <- x$event
    x$n.e <- x$n
    ##
    if (!is.null(rightcols)) {
      if (any(rightcols == "n"))
        rightcols[rightcols == "n"] <- "n.e"
      if (any(rightcols == "event"))
        rightcols[rightcols == "event"] <- "event.e"
    }
    ##
    if (!is.null(leftcols)) {
      if (any(leftcols == "n"))
        leftcols[leftcols == "n"] <- "n.e"
      if (any(leftcols == "event"))
        leftcols[leftcols == "event"] <- "event.e"
    }
  }
  if (metarate) {
    x$event.e <- x$event
    x$time.e <- x$time
    ##
    if (!is.null(rightcols)) {
      if (any(rightcols == "time"))
        rightcols[rightcols == "time"] <- "time.e"
      if (any(rightcols == "event"))
        rightcols[rightcols == "event"] <- "event.e"
    }
    ##
    if (!is.null(leftcols)) {
      if (any(leftcols == "time"))
        leftcols[leftcols == "time"] <- "time.e"
      if (any(leftcols == "event"))
        leftcols[leftcols == "event"] <- "event.e"
    }
  }
  ##
  if (metainf.metacum) {
    ##
    x$TE.fixed    <- rev(x$TE)[1]
    x$seTE.fixed  <- rev(x$seTE)[1]
    x$lower.fixed <- rev(x$lower)[1]
    x$upper.fixed <- rev(x$upper)[1]
    ##
    x$TE.random   <- rev(x$TE)[1]
    x$seTE.random <- rev(x$seTE)[1]
    x$lower.random <- rev(x$lower)[1]
    x$upper.random <- rev(x$upper)[1]
    ##
    x$n.harmonic.mean.ma <- rev(x$n.harmonic.mean)[1]
    x$t.harmonic.mean.ma <- rev(x$t.harmonic.mean)[1]
    ##
    x$w.all <- rev(x$w)[1]
    ##
    x$TE <- rev(rev(x$TE)[-(1:2)])
    x$seTE <- rev(rev(x$seTE)[-(1:2)])
    x$studlab <- rev(rev(x$studlab)[-(1:2)])
    ##
    x$lower <- rev(rev(x$lower)[-(1:2)])
    x$upper <- rev(rev(x$upper)[-(1:2)])
    ##
    x$w.fixed <- rev(rev(x$w)[-(1:2)])
    x$w.random <- rev(rev(x$w)[-(1:2)])
    ##
    x$n.harmonic.mean <- rev(rev(x$n.harmonic.mean)[-(1:2)])
    x$t.harmonic.mean <- rev(rev(x$t.harmonic.mean)[-(1:2)])
    ##
    if (overall & x$pooled == "fixed") {
      comb.fixed <- TRUE
      comb.random <- FALSE
      if (weight.study != "same")
        weight.study <- "fixed"
    }
    else if (overall & x$pooled == "random") {
      comb.fixed <- FALSE
      comb.random <- TRUE
      if (weight.study != "same")
        weight.study <- "random"
    }
  }
  ## Total number of studies to plot (*not* number of studies combined)
  k.all <- length(x$TE)
  ##
  if (allstudies)
    n.stud <- k.all # all studies
  else
    n.stud <- sum(!is.na(x$TE)) # number of studies with treatment estimates
  ##
  if (length(type.study) == 1)
    type.study <- rep(type.study, k.all)
  else if (length(type.study) != k.all)
    stop("Argument 'type.study' must be a single character or of ",
         "same length as number of studies.")
  ##
  if (!by)
    byvar <- rep(1, k.all)
  ##
  if (by & anyNA(byvar))
    stop("Missing values in 'byvar'")
  ##
  if (allstudies)
    sel <- 1:k.all
  else
    sel <- !is.na(x$TE)
  ##
  x$n.e <- x$n.e[sel]
  x$n.c <- x$n.c[sel]
  ##
  x$event.e <- x$event.e[sel]
  x$event.c <- x$event.c[sel]
  ##
  x$mean.e <- x$mean.e[sel]
  x$mean.c <- x$mean.c[sel]
  ##
  x$sd.e <- x$sd.e[sel]
  x$sd.c <- x$sd.c[sel]
  ##
  x$cor <- x$cor[sel]
  ##
  x$time.e <- x$time.e[sel]
  x$time.c <- x$time.c[sel]
  ##
  x$TE   <- x$TE[sel]
  x$seTE <- x$seTE[sel]
  ##
  x$lower <- x$lower[sel]
  x$upper <- x$upper[sel]
  ##
  x$w.fixed  <- x$w.fixed[sel]
  x$w.random <- x$w.random[sel]
  studlab  <- studlab[sel]
  type.study  <- type.study[sel]
  ##
  x$n.harmonic.mean <- x$n.harmonic.mean[sel]
  x$t.harmonic.mean <- x$t.harmonic.mean[sel]
  ##
  byvar   <- byvar[sel]
  sortvar <- sortvar[sel]
  ##
  col.study <- col.study[sel]
  col.square <- col.square[sel]
  col.square.lines <- col.square.lines[sel]
  ##
  col.inside <- col.inside[sel]
  ##
  null.exclude <- is.null(x$exclude)
  if (!null.exclude)
    exclude <- x$exclude[sel]
  ##
  if (sort | by) {
    if (bysort)
      bylevs <- sort(x$bylevs)
    else
      bylevs <- x$bylevs
    ##
    byvar.factor <- factor(byvar, levels = bylevs)
    o <- order(byvar.factor, sortvar)
    ##
    x$n.e <- x$n.e[o]
    x$n.c <- x$n.c[o]
    ##
    x$event.e <- x$event.e[o]
    x$event.c <- x$event.c[o]
    ##
    x$mean.e <- x$mean.e[o]
    x$mean.c <- x$mean.c[o]
    ##
    x$sd.e <- x$sd.e[o]
    x$sd.c <- x$sd.c[o]
    ##
    x$cor <- x$cor[o]
    ##
    x$time.e <- x$time.e[o]
    x$time.c <- x$time.c[o]
    ##
    x$TE   <- x$TE[o]
    x$seTE <- x$seTE[o]
    ##
    x$lower <- x$lower[o]
    x$upper <- x$upper[o]
    ##
    x$w.fixed  <- x$w.fixed[o]
    x$w.random <- x$w.random[o]
    studlab  <- studlab[o]
    type.study  <- type.study[o]
    ##
    x$n.harmonic.mean <- x$n.harmonic.mean[o]
    x$t.harmonic.mean <- x$t.harmonic.mean[o]
    ##
    byvar   <- byvar[o]
    sortvar <- sortvar[o]
    ##
    col.study <- col.study[o]
    col.square <- col.square[o]
    col.square.lines <- col.square.lines[o]
    ##
    col.inside <- col.inside[o]
    ##
    if (!null.exclude)
      exclude <- exclude[o]
    ##
    if (newcols) {
      dataset1 <- dataset1[o, ]
      dataset2 <- dataset2[o, ]
    }
  }
  ##
  if (by)
    n.by <- length(bylevs)
  else
    n.by <- 0
  ##
  if (metainf.metacum) {
    TE    <- x$TE
    seTE  <- x$seTE
    lowTE <- x$lower
    uppTE <- x$upper
    ##
    TE.fixed    <- x$TE.fixed
    lowTE.fixed <- x$lower.fixed
    uppTE.fixed <- x$upper.fixed
    ##
    TE.random    <- x$TE.random
    lowTE.random <- x$lower.random
    uppTE.random <- x$upper.random
    ##
    lowTE.predict <- NA
    uppTE.predict <- NA
    ##
    Q    <- NA
    df.Q <- NA
    pval.Q <- NA
    I2   <- NA
    tau2 <- NA
    lower.tau2 <- NA
    upper.tau2 <- NA
    tau <- NA
    lower.tau <- NA
    upper.tau <- NA
    sign.lower.tau <- ""
    sign.upper.tau <- ""
    lowI2 <- NA
    uppI2 <- NA
    Rb   <- NA
    lowRb <- NA
    uppRb <- NA
    ##
    Q.b.fixed  <- NA
    Q.b.random <- NA
    df.Q.b     <- NA
    pval.Q.b.fixed  <- NA
    pval.Q.b.random <- NA
  }
  else {
    TE <- x$TE
    seTE <- x$seTE
    lowTE <- x$lower
    uppTE <- x$upper
    ##
    if (metaprop & !backtransf) {
      ciTE <- ci(TE, seTE, level = level)
      lowTE <- ciTE$lower
      uppTE <- ciTE$upper
    }
    ##
    TE.fixed <- x$TE.fixed
    lowTE.fixed <- x$lower.fixed
    uppTE.fixed <- x$upper.fixed
    ##
    TE.random <- x$TE.random
    lowTE.random <- x$lower.random
    uppTE.random <- x$upper.random
    ##
    lowTE.predict <- x$lower.predict
    uppTE.predict <- x$upper.predict
    ##
    if (LRT)
      Q <- x$Q.LRT
    else
      Q <- x$Q
    df.Q <- x$df.Q
    pval.Q <- replaceNULL(x$pval.Q, pvalQ(Q, df.Q))
    ##
    tau2 <- x$tau2
    lower.tau2 <- x$lower.tau2
    upper.tau2 <- x$upper.tau2
    ##
    if (is.null(tau2)) {
      tau2 <- x$tau^2
      lower.tau2 <- upper.tau2 <- NA
    }
    ##
    if (length(tau2) > 1) {
      tau2 <- sum(tau2)
      lower.tau2 <- NA
      upper.tau2 <- NA
    }
    ##
    tau <- x$tau
    lower.tau <- x$lower.tau
    upper.tau <- x$upper.tau
    ##
    if (length(tau) > 1) {
      tau <- sqrt(sum(tau^2))
      lower.tau <- NA
      upper.tau <- NA
    }
    ##
    sign.lower.tau <- x$sign.lower.tau
    sign.upper.tau <- x$sign.lower.tau
    ##
    if (is.null(lower.tau)) {
      lower.tau <- upper.tau <- NA
      sign.lower.tau <- sign.upper.tau <- ""
    }
    ##
    I2 <- x$I2
    lowI2 <- x$lower.I2
    uppI2 <- x$upper.I2
    ##
    Rb <- x$Rb
    lowRb <- x$lower.Rb
    uppRb <- x$upper.Rb
    ##
    if (by) {
      Q.b.fixed <- x$Q.b.fixed
      Q.b.random <- x$Q.b.random
      df.Q.b <- x$df.Q.b
      pval.Q.b.fixed <-
        replaceNULL(x$pval.Q.b.fixed, pvalQ(Q.b.fixed, df.Q.b))
      pval.Q.b.random <-
        replaceNULL(x$pval.Q.b.random, pvalQ(Q.b.random, df.Q.b))
      ##
      Q.resid <- x$Q.resid
      df.Q.resid <- x$df.Q.resid
      pval.Q.resid <- x$pval.Q.resid
      ##
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
      ##
      I2.resid <- x$I2.resid
      lowI2.resid <- x$lower.I2.resid
      uppI2.resid <- x$upper.I2.resid
    }
    else {
      Q.b.fixed  <- NA
      Q.b.random <- NA
      df.Q.b     <- NA
    }
  }
  ##
  hetstat.overall <- ""
  ##
  if (overall.hetstat || is.character(hetstat)) {
    ##
    hetstat.I2 <-
      paste0(hetseparator,
             formatN(round(100 * I2, digits.I2),
                     digits.I2, "NA"), "%",
             if (print.I2.ci & !(is.na(lowI2) | is.na(uppI2)))
               pasteCI(100 * lowI2, 100 * uppI2,
                       digits.I2, big.mark,
                       sign.lower.tau, sign.upper.tau, lab.NA,
                       "%"))
    ##
    hetstat.tau2 <-
      paste0(formatPT(tau2, digits = digits.tau2, big.mark = big.mark,
                      lab = TRUE, labval = "", lab.NA = "NA"),
             if (print.tau2.ci & !(is.na(lower.tau2) | is.na(upper.tau2)))
               pasteCI(lower.tau2, upper.tau2, digits.tau2, big.mark,
                       sign.lower.tau, sign.upper.tau, lab.NA))
    ##
    hetstat.tau <-
      paste0(formatPT(tau, digits = digits.tau, big.mark = big.mark,
                      lab = TRUE, labval = "", lab.NA = "NA"),
             if (print.tau.ci & !(is.na(lower.tau) | is.na(upper.tau)))
               pasteCI(lower.tau, upper.tau, digits.tau, big.mark,
                       sign.lower.tau, sign.upper.tau, lab.NA))
    ##
    hetstat.Q <-
      paste0(hetseparator,
             formatN(round(Q, digits.Q), digits.Q, "NA", big.mark = big.mark),
             if (revman5) ", df",
             if (revman5) hetseparator,
             if (revman5) df.Q)
    ##
    hetstat.pval.Q <-
      formatPT(pval.Q,
               lab = TRUE, labval = "",
               digits = digits.pval.Q,
               zero = zero.pval, JAMA = JAMA.pval,
               scientific = scientific.pval,
               lab.NA = "NA")
    ##
    hetstat.Rb <-
      paste0(hetseparator,
             formatN(round(100 * Rb, digits.I2),
                     digits.I2, "NA", big.mark = big.mark),
             "%",
             if (print.Rb.ci & !(is.na(lowRb) | is.na(uppRb)))
               pasteCI(100 * lowRb, 100 * uppRb,
                       digits.I2, big.mark,
                       sign.lower.tau, sign.upper.tau, lab.NA,
                       "%"))
    ##
    ## Remove superfluous spaces
    ##
    while(grepl("  ", hetstat.I2))
      hetstat.I2 <- gsub("  ", " ", hetstat.I2)
    while(grepl("  ", hetstat.tau2))
      hetstat.tau2 <- gsub("  ", " ", hetstat.tau2)
    while(grepl("  ", hetstat.Q))
      hetstat.Q <- gsub("  ", " ", hetstat.Q)
    while(grepl("  ", hetstat.pval.Q))
      hetstat.pval.Q <- gsub("  ", " ", hetstat.pval.Q)
    while(grepl("  ", hetstat.Rb))
      hetstat.Rb <- gsub("  ", " ", hetstat.Rb)
    ##
    if (revman5)
      hetstat.overall <- substitute(paste(hl,
                                          "Tau"^2, ht, "; ",
                                          "Chi"^2, hq,
                                          " (",
                                          P, hp,
                                          "); ",
                                          I^2, hi),
                                    list(hl = hetlab,
                                         hi = hetstat.I2,
                                         ht = hetstat.tau2,
                                         hq = hetstat.Q,
                                         hp = hetstat.pval.Q))
    else if (jama)
      hetstat.overall <- substitute(paste(hl,
                                          chi[df]^2, hq,
                                          " (",
                                          italic(P), hp,
                                          "), ",
                                          italic(I)^2, hi
                                          ),
                                    list(hl = hetlab,
                                         hi = hetstat.I2,
                                         hq = hetstat.Q,
                                         hp = hetstat.pval.Q,
                                         df = df.Q))
    else {
      ##
      ## One
      ##
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
      ##
      ## Two
      ##
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
      ##
      ## Three
      ##
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
      ##
      ## Four
      ##
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
      ##
      ## Five
      ##
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
  ##
  ## Line with residual heterogeneity
  ##
  hetstat.resid <- ""
  ##
  if (by && (length(tau2.resid) == 0 || is.na(tau2.resid)))
    print.tau2.tau.resid <- FALSE
  else
    print.tau2.tau.resid <- print.tau2.tau
  ##
  if (resid.hetstat) {
    ##
    hetstat.I2.resid <-
      paste0(hetseparator,
             formatN(round(100 * I2.resid, digits.I2),
                     digits.I2, "NA"), "%",
             if (print.I2.ci & !(is.na(lowI2.resid) | is.na(uppI2.resid)))
               pasteCI(100 * lowI2.resid, 100 * uppI2.resid,
                       digits.I2, big.mark,
                       sign.lower.tau.resid, sign.upper.tau.resid, lab.NA,
                       "%"))
    ##
    hetstat.tau2.resid <-
      paste0(formatPT(tau2.resid, digits = digits.tau2, big.mark = big.mark,
                      lab = TRUE, labval = "", lab.NA = "NA"),
             if (print.tau2.ci &
                 !(is.na(lower.tau2.resid) | is.na(upper.tau2.resid)))
               pasteCI(lower.tau2.resid, upper.tau2.resid, digits.tau2,
                       big.mark,
                       sign.lower.tau.resid, sign.upper.tau.resid, lab.NA))
    hetstat.tau.resid <-
      paste0(formatPT(tau.resid, digits = digits.tau, big.mark = big.mark,
                      lab = TRUE, labval = "", lab.NA = "NA"),
             if (print.tau.ci &
                 !(is.na(lower.tau.resid) | is.na(upper.tau.resid)))
               pasteCI(lower.tau.resid, upper.tau.resid, digits.tau,
                       big.mark,
                       sign.lower.tau.resid, sign.upper.tau.resid, lab.NA))
    ##
    hetstat.Q.resid <-
      paste0(hetseparator,
             formatN(round(Q.resid, digits.Q), digits.Q, "NA",
                     big.mark = big.mark),
             if (revman5) ", df",
             if (revman5) hetseparator,
             if (revman5) df.Q.resid)
    ##
    hetstat.pval.Q.resid <-
      formatPT(pval.Q.resid,
               lab = TRUE, labval = "",
               digits = digits.pval.Q,
               zero = zero.pval, JAMA = JAMA.pval,
               scientific = scientific.pval,
               lab.NA = "NA")
    ##
    hetstat.Rb.resid <- ""
    print.Rb.resid <- FALSE
    ##
    ## Remove superfluous spaces
    ##
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
    ##
    if (revman5)
      hetstat.resid <- substitute(paste(hl,
                                          "Tau"^2, ht, "; ",
                                          "Chi"^2, hq,
                                          " (",
                                          P, hp,
                                          "); ",
                                          I^2, hi),
                                    list(hl = resid.hetlab,
                                         hi = hetstat.I2.resid,
                                         ht = hetstat.tau2.resid,
                                         hq = hetstat.Q.resid,
                                         hp = hetstat.pval.Q.resid))
    else if (jama)
      hetstat.resid <- substitute(paste(hl,
                                          chi[df]^2, hq,
                                          " (",
                                          italic(P), hp,
                                          "), ",
                                          italic(I)^2, hi
                                          ),
                                    list(hl = resid.hetlab,
                                         hi = hetstat.I2.resid,
                                         hq = hetstat.Q.resid,
                                         hp = hetstat.pval.Q.resid,
                                         df = df.Q.resid))
    else {
      ##
      ## One
      ##
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
                     list(hl = resid.hetlab, df = df.Q.resid, hq = hetstat.Q.resid))
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
      ##
      ## Two
      ##
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
      ##
      ## Three
      ##
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
      ##
      ## Four
      ##
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
      ##
      ## Five
      ##
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
  ##
  ## Label of test for overall effect
  ##
  ##
  if (test.overall.fixed | test.overall.random) {
    pvals.overall <- formatPT(c(x$pval.fixed, x$pval.random),
                              lab = TRUE, labval = "",
                              digits = digits.pval,
                              zero = zero.pval, JAMA = JAMA.pval,
                              scientific = scientific.pval,
                              lab.NA = "NA")
    statistics.overall <- formatN(round(c(x$statistic.fixed, x$statistic.random),
                                        digits = digits.stat),
                                  digits.stat, "NA", big.mark = big.mark)
    ##
    ## Remove superfluous spaces
    ##
    pvals.overall <- rmSpace(pvals.overall, end = TRUE)
    ##
    while(any(grepl("  ", pvals.overall)))
      pvals.overall <- gsub("  ", " ", pvals.overall)
    while(any(grepl("  ", statistics.overall)))
      statistics.overall <- gsub("  ", " ", statistics.overall)
  }
  ##
  if (test.overall.fixed) {
    if (print.stat) {
      if (revman5)
        text.overall.fixed  <- substitute(paste(tl,
                                                Z, hetseparator, tt,
                                                " (P", tp, ")"),
                                          list(tl = label.test.overall.fixed,
                                               hetseparator = hetseparator,
                                               tt = statistics.overall[1],
                                               tp = pvals.overall[1]))
      else if (jama)
        text.overall.fixed  <- substitute(paste(tl,
                                                italic(z), hetseparator, tt,
                                                " (", italic(P), tp, ")"),
                                          list(tl = label.test.overall.fixed,
                                               hetseparator = hetseparator,
                                               tt = statistics.overall[1],
                                               tp = pvals.overall[1]))
      else
        text.overall.fixed  <- substitute(paste(tl,
                                                italic(z), hetseparator, tt,
                                                " (", italic(p), tp, ")"),
                                          list(tl = label.test.overall.fixed,
                                               hetseparator = hetseparator,
                                               tt = statistics.overall[1],
                                               tp = pvals.overall[1]))
    }
    else {
      if (revman5)
        text.overall.fixed  <- substitute(paste(tl, " P", tp),
                                          list(tl = label.test.overall.fixed,
                                               tp = pvals.overall[1]))
      else if (jama)
        text.overall.fixed  <- substitute(paste(tl, " ", italic(P), tp),
                                          list(tl = label.test.overall.fixed,
                                               tp = pvals.overall[1]))
      else
        text.overall.fixed  <- substitute(paste(tl, " ", italic(p), tp),
                                          list(tl = label.test.overall.fixed,
                                               tp = pvals.overall[1]))     
    }
  }
  else
    text.overall.fixed <- ""
  ##
  if (test.overall.random) {
    if (print.stat) {
      if (!x$hakn) {
        if (revman5)
          text.overall.random  <- substitute(paste(tl,
                                                   Z, hetseparator, tt,
                                                   " (P", tp, ")"),
                                             list(tl = label.test.overall.random,
                                                  hetseparator = hetseparator,
                                                  tt = statistics.overall[2],
                                                  tp = pvals.overall[2]))
        else if (jama)
          text.overall.random  <- substitute(paste(tl,
                                                   italic(z), hetseparator, tt,
                                                   " (", italic(P), tp, ")"),
                                             list(tl = label.test.overall.random,
                                                  hetseparator = hetseparator,
                                                  tt = statistics.overall[2],
                                                  tp = pvals.overall[2]))
        else
          text.overall.random  <- substitute(paste(tl,
                                                   italic(z), hetseparator, tt,
                                                   " (", italic(p), tp, ")"),
                                             list(tl = label.test.overall.random,
                                                  hetseparator = hetseparator,
                                                  tt = statistics.overall[2],
                                                  tp = pvals.overall[2]))
      }
      else {
        if (revman5)
          text.overall.random  <- substitute(paste(tl,
                                                   t[df], hetseparator, tt,
                                                   " (P", tp, ")"),
                                             list(tl = label.test.overall.random,
                                                  hetseparator = hetseparator,
                                                  tt = statistics.overall[2],
                                                  tp = pvals.overall[2],
                                                  df = round(x$df.hakn, 1)))
        else if (jama)
          text.overall.random  <- substitute(paste(tl,
                                                   italic(t)[df], hetseparator, tt,
                                                   " (", italic(P), tp, ")"),
                                             list(tl = label.test.overall.random,
                                                  hetseparator = hetseparator,
                                                  tt = statistics.overall[2],
                                                  tp = pvals.overall[2],
                                                  df = round(x$df.hakn, 1)))
        else
          text.overall.random  <- substitute(paste(tl,
                                                   italic(t)[df], hetseparator, tt,
                                                   " (", italic(p), tp, ")"),
                                             list(tl = label.test.overall.random,
                                                  hetseparator = hetseparator,
                                                  tt = statistics.overall[2],
                                                  tp = pvals.overall[2],
                                                  df = round(x$df.hakn, 1)))
      }
    }
    else {
      if (revman5)
        text.overall.random  <- substitute(paste(tl, " P", tp),
                                           list(tl = label.test.overall.random,
                                                tp = pvals.overall[2]))
      else if (jama)
        text.overall.random  <- substitute(paste(tl, " ", italic(P), tp),
                                           list(tl = label.test.overall.random,
                                                tp = pvals.overall[2]))
      else
        text.overall.random  <- substitute(paste(tl, " ", italic(p), tp),
                                           list(tl = label.test.overall.random,
                                                tp = pvals.overall[2]))
    }
  }
  else
    text.overall.random <- ""
  ##
  ##
  ## Label of test for subgroup differences
  ##
  if (by) {
    Q.bs <- c(Q.b.fixed, Q.b.random)
    pval.Q.bs <- c(pval.Q.b.fixed, pval.Q.b.random)
  }
  else {
    Q.bs <- NA
    pval.Q.bs <- NA
  }
  ##
  hetstat.Q.bs <-
    paste0(hetseparator,
           gsub(" ", "", formatN(round(Q.bs, digits.Q), digits.Q, "NA",
                                 big.mark = big.mark)),
           if (!jama) ", df",
           if (!jama) hetseparator,
           if (!jama) df.Q.b)
  hetstat.Q.bs <- rmSpace(hetstat.Q.bs, end = TRUE)
  ##
  hetstat.pval.Q.bs <-
    paste0(formatPT(pval.Q.bs,
                    lab = TRUE, labval = "",
                    digits = digits.pval.Q,
                    zero = zero.pval, JAMA = JAMA.pval,
                    scientific = scientific.pval,
                    lab.NA = "NA"))
  ##
  ## Remove superfluous spaces
  ##
  hetstat.pval.Q.bs <- rmSpace(hetstat.pval.Q.bs, end = TRUE)
  ##
  while(any(grepl("  ", hetstat.pval.Q.bs)))
    hetstat.pval.Q.bs <- gsub("  ", " ", hetstat.pval.Q.bs)
  while(any(grepl("  ", hetstat.Q.bs)))
    hetstat.Q.bs <- gsub("  ", " ", hetstat.Q.bs)
  ##
  if (test.subgroup.fixed) {
    if (print.Q.subgroup) {
      if (revman5)
        text.subgroup.fixed <-
          substitute(paste(tl,
                           "Chi"^2, tq,
                           " (P", tp, ")"),
                     list(tl = label.test.subgroup.fixed,
                          tq = hetstat.Q.bs[1],
                          tp = hetstat.pval.Q.bs[1]))
      else if (jama)
        text.subgroup.fixed <-
          substitute(paste(tl,
                           chi[df]^2, tq,
                           " (", italic(P), tp, ")"),
                     list(tl = label.test.subgroup.fixed,
                          tq = hetstat.Q.bs[1],
                          tp = hetstat.pval.Q.bs[1],
                          df = df.Q.b))
      else
        text.subgroup.fixed <-
          substitute(paste(tl,
                           chi[df]^2, tq,
                           " (", italic(p), tp, ")"),
                     list(tl = label.test.subgroup.fixed,
                          tq = hetstat.Q.bs[1],
                          tp = hetstat.pval.Q.bs[1],
                          df = df.Q.b))
    }
    else {
      if (revman5)
        text.subgroup.fixed <-
          substitute(paste(tl, " P", tp),
                     list(tl = label.test.subgroup.fixed,
                          tp = hetstat.pval.Q.bs[1]))
      else if (jama)
        text.subgroup.fixed <-
          substitute(paste(tl, " ", italic(P), tp),
                     list(tl = label.test.subgroup.fixed,
                          tp = hetstat.pval.Q.bs[1]))
      else
        text.subgroup.fixed <-
          substitute(paste(tl, " ", italic(p), tp),
                     list(tl = label.test.subgroup.fixed,
                          tp = hetstat.pval.Q.bs[1]))
    }
  }
  else
    text.subgroup.fixed <- ""
  
  
  if (test.subgroup.random) {
    if (print.Q.subgroup) {
      if (revman5)
        text.subgroup.random <-
          substitute(paste(tl,
                           "Chi"^2, tq,
                           " (P", tp, ")"),
                     list(tl = label.test.subgroup.random,
                          tq = hetstat.Q.bs[2],
                          tp = hetstat.pval.Q.bs[2]))
      else if (jama)
        text.subgroup.random <-
          substitute(paste(tl,
                           chi[df]^2, tq,
                           " (", italic(P), tp, ")"),
                     list(tl = label.test.subgroup.random,
                          tq = hetstat.Q.bs[2],
                          tp = hetstat.pval.Q.bs[2],
                          df = df.Q.b))
      else
        text.subgroup.random <-
          substitute(paste(tl,
                           chi[df]^2, tq,
                           " (", italic(p), tp, ")"),
                     list(tl = label.test.subgroup.random,
                          tq = hetstat.Q.bs[2],
                          tp = hetstat.pval.Q.bs[2],
                          df = df.Q.b))
    }
    else {
      if (revman5)
        text.subgroup.random <-
          substitute(paste(tl, " P", tp),
                     list(tl = label.test.subgroup.random,
                          tp = hetstat.pval.Q.bs[2]))
      else if (jama)
        text.subgroup.random <-
          substitute(paste(tl, " ", italic(P), tp),
                     list(tl = label.test.subgroup.random,
                          tp = hetstat.pval.Q.bs[2]))
      else
        text.subgroup.random <-
          substitute(paste(tl, " ", italic(p), tp),
                     list(tl = label.test.subgroup.random,
                          tp = hetstat.pval.Q.bs[2]))
    }
  }
  else
    text.subgroup.random <- ""
  
  
  ##
  ##
  ## (7) Prepare data for subgroup analysis
  ##
  ##
  if (by) {
    k.w.hetstat <- if (metabind) x$k.w.orig else x$k.w
    ##
    o <- order(factor(x$bylevs, levels = bylevs))
    k.w <- x$k.w[o]
    k.w.hetstat <- k.w.hetstat[o]
    ##
    TE.fixed.w <- x$TE.fixed.w[o]
    lower.fixed.w <- x$lower.fixed.w[o]
    upper.fixed.w <- x$upper.fixed.w[o]
    statistic.fixed.w <- x$statistic.fixed.w[o]
    pval.fixed.w <- x$pval.fixed.w[o]
    ##
    TE.random.w <- x$TE.random.w[o]
    lower.random.w <- x$lower.random.w[o]
    upper.random.w <- x$upper.random.w[o]
    statistic.random.w <- x$statistic.random.w[o]
    pval.random.w <- x$pval.random.w[o]
    ##
    lower.predict.w <- x$lower.predict.w[o]
    upper.predict.w <- x$upper.predict.w[o]
    ##
    Q.w     <- x$Q.w[o]
    I2.w    <- x$I2.w[o]
    lowI2.w <- x$lower.I2.w[o]
    uppI2.w <- x$upper.I2.w[o]
    Rb.w    <- x$Rb.w[o]
    lowRb.w <- x$lower.Rb.w[o]
    uppRb.w <- x$upper.Rb.w[o]
    ##
    tau2.w <- x$tau2.w[o]
    lower.tau2.w <- x$lower.tau2.w[o]
    upper.tau2.w <- x$upper.tau2.w[o]
    tau.w <- x$tau.w[o]
    lower.tau.w <- x$lower.tau.w[o]
    upper.tau.w <- x$upper.tau.w[o]
    sign.lower.tau.w <- x$sign.lower.tau.w[o]
    sign.upper.tau.w <- x$sign.upper.tau.w[o]
    ##
    w.fixed.w  <- x$w.fixed.w[o]
    w.random.w <- x$w.random.w[o]
    e.e.w <- if (metaprop | metarate) x$event.w[o] else x$event.e.w[o]
    t.e.w <- if (metainc | metarate) x$time.e.w[o] else rep(NA, n.by)
    n.e.w <- if (metacor | metaprop) x$n.w[o] else x$n.e.w[o]
    e.c.w <- x$event.c.w[o]
    t.c.w <- if (metainc) x$time.c.w[o] else rep(NA, n.by)
    n.c.w <- x$n.c.w[o]
    n.harmonic.mean.w <- x$n.harmonic.mean.w[o]
    t.harmonic.mean.w <- x$t.harmonic.mean.w[o]
    k.all.w <- x$k.all.w[o]
    ##
    if (allstudies)
      sel <- 1:length(k.w)
    else
      sel <- k.w > 0
    k.w <- k.w[sel]
    TE.fixed.w <- TE.fixed.w[sel]
    lower.fixed.w <- lower.fixed.w[sel]
    upper.fixed.w <- upper.fixed.w[sel]
    statistic.fixed.w <- statistic.fixed.w[sel]
    pval.fixed.w <- pval.fixed.w[sel]
    ##
    TE.random.w <- TE.random.w[sel]
    lower.random.w <- lower.random.w[sel]
    upper.random.w <- upper.random.w[sel]
    statistic.random.w <- statistic.random.w[sel]
    pval.random.w <- pval.random.w[sel]
    ##
    lower.predict.w <- lower.predict.w[sel]
    upper.predict.w <- upper.predict.w[sel]
    ##
    Q.w     <- Q.w[sel]
    I2.w    <- I2.w[sel]
    lowI2.w <- lowI2.w[sel]
    uppI2.w <- uppI2.w[sel]
    Rb.w    <- Rb.w[sel]
    lowRb.w <- lowRb.w[sel]
    uppRb.w <- uppRb.w[sel]
    ##
    tau2.w <- tau2.w[sel]
    lower.tau2.w <- lower.tau2.w[sel]
    upper.tau2.w <- upper.tau2.w[sel]
    tau.w <- tau.w[sel]
    lower.tau.w <- lower.tau.w[sel]
    upper.tau.w <- upper.tau.w[sel]
    sign.lower.tau.w <- sign.lower.tau.w[sel]
    sign.upper.tau.w <- sign.upper.tau.w[sel]
    ##
    w.fixed.w  <- w.fixed.w[sel]
    w.random.w <- w.random.w[sel]
    e.e.w <- e.e.w[sel]
    t.e.w <- t.e.w[sel]
    n.e.w <- n.e.w[sel]
    e.c.w <- e.c.w[sel]
    t.c.w <- t.c.w[sel]
    n.c.w <- n.c.w[sel]
    n.harmonic.mean.w <- n.harmonic.mean.w[sel]
    t.harmonic.mean.w <- t.harmonic.mean.w[sel]
    k.all.w <- k.all.w[sel]
    ##
    bylevs <- bylevs[sel]
    n.by <- length(bylevs)
    ##
    sel.by.fixed         <- 3 + 0 * n.by + 1:n.by
    sel.by.random        <- 3 + 1 * n.by + 1:n.by
    sel.by.predict       <- 3 + 2 * n.by + 1:n.by
    sel.by.het           <- 3 + 3 * n.by + 1:n.by
    sel.by.effect.fixed  <- 3 + 4 * n.by + 1:n.by
    sel.by.effect.random <- 3 + 5 * n.by + 1:n.by
    sel.by <- c(sel.by.fixed, sel.by.random, sel.by.predict, sel.by.het,
                sel.by.effect.fixed, sel.by.effect.random)
    sel.by.noNA <- c(sel.by.het, sel.by.effect.fixed, sel.by.effect.random)
    ##
    if (!metainf.metacum & comb.fixed) {
      if (!overall) {
        i <- 0
        for (bylev.i in bylevs) {
          i <- i + 1
          sel.i <- byvar == bylev.i
          x$w.fixed[sel.i] <- x$w.fixed[sel.i] / w.fixed.w[i]
        }
        w.fixed.w.p <- ifelse(is.na(w.fixed.w), NA, 100)
      }
      else {
        if (!all(is.na(w.fixed.w)) && sum(w.fixed.w) > 0)
          w.fixed.w.p <-
            round(100 * w.fixed.w / sum(w.fixed.w, na.rm = TRUE),
                  digits.weight)
        else
          w.fixed.w.p <- w.fixed.w
      }
    }
    else {
      TE.fixed.w <- lower.fixed.w <- upper.fixed.w <- rep(NA, n.by)
      ##
      w.fixed.w.p <- rep(NA, n.by)
      if (missing(text.fixed.w))
        text.fixed.w <- rep("Overall", n.by)
    }
    ##
    if (!metainf.metacum & comb.random) {
      if (!overall) {
        i <- 0
        for (bylev.i in bylevs) {
          i <- i + 1
          sel.i <- byvar == bylev.i
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
      TE.random.w <- lower.random.w <- upper.random.w <- rep(NA, n.by)
      ##
      w.random.w.p <- rep(NA, n.by)
      text.random.w <- rep("", n.by)
    }
    ##
    if (metainf.metacum) {
      lower.predict.w <- upper.predict.w <- rep(NA, n.by)
      text.predict.w <- rep("", n.by)
    }
    ##
    hetstat.w <- vector("list", n.by)
    ##
    if (is.character(hetstat) || hetstat) {
      ##
      hetstat.I2.w <-
        paste0(hetseparator,
               round(100 * I2.w, digits.I2), "%",
               if (print.I2.ci)
                 ifelse(!(is.na(lowI2.w) | is.na(uppI2.w)),
                        pasteCI(100 * lowI2.w, 100 * uppI2.w,
                                digits.I2, big.mark,
                                sign.lower.tau.w, sign.upper.tau.w, lab.NA,
                                "%"),
                        ""))
      ##
      hetstat.tau2.w <-
        paste0(hetseparator,
               ifelse(is.na(tau.w), "NA",
               ifelse(tau2.w == 0, "0",
                      formatPT(tau2.w, digits = digits.tau2,
                               big.mark = big.mark, lab.NA = "NA"))),
               if (print.tau2.ci)
                 ifelse(!(is.na(lower.tau2.w) | is.na(upper.tau2.w)),
                        pasteCI(lower.tau2.w, upper.tau2.w,
                                digits.tau2, big.mark,
                                sign.lower.tau.w, sign.upper.tau.w, lab.NA),
                        ""))
      ##
      hetstat.tau.w <-
        paste0(hetseparator,
               ifelse(is.na(tau.w), "NA",
               ifelse(tau.w == 0, "0",
                      formatPT(tau.w, digits = digits.tau,
                               big.mark = big.mark, lab.NA = "NA"))),
               if (print.tau.ci)
                 ifelse(!(is.na(lower.tau.w) | is.na(upper.tau.w)),
                        pasteCI(lower.tau.w, upper.tau.w, digits.tau, big.mark,
                                sign.lower.tau.w, sign.upper.tau.w, lab.NA),
                        ""))
      ##
      hetstat.Q.w <-
        paste0(hetseparator, round(Q.w, digits.Q),
               if (revman5)
                 paste0(", df",  hetseparator, k.w.hetstat - 1))
      ##
      hetstat.pval.Q.w <-
        paste0(formatPT(replaceNULL(x$pval.Q.w, pvalQ(x$Q.w, x$df.Q.w)),
                        lab = TRUE, labval = "",
                        digits = digits.pval.Q,
                        zero = zero.pval, JAMA = JAMA.pval,
                        scientific = scientific.pval,
                        lab.NA = "NA"))
      ##
      hetstat.Rb.w <-
        paste0(hetseparator,
               round(100 * Rb.w, digits.I2), "%",
               if (print.Rb.ci)
                 ifelse(!(is.na(lowRb.w) | is.na(uppRb.w)),
                        pasteCI(100 * lowRb.w, 100 * uppRb.w,
                                digits.I2, big.mark,
                                sign.lower.tau.w, sign.upper.tau.w, lab.NA,
                                "%"),
                        ""))
      ##
      ## Remove superfluous spaces
      ##
      hetstat.pval.Q.w <- rmSpace(hetstat.pval.Q.w, end = TRUE)
      ##
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
      ##
      for (i in 1:n.by) {
        if (revman5)
          hetstat.w[[i]] <- substitute(paste(hl,
                                             "Tau"^2, ht, "; ",
                                             "Chi"^2, hq,
                                             " (",
                                             P, hp,
                                             "); ",
                                             I^2, hi),
                                       list(hl = hetlab,
                                            hi = hetstat.I2.w[i],
                                            ht = hetstat.tau2.w[i],
                                            hq = hetstat.Q.w[i],
                                            hp = hetstat.pval.Q.w[i]))
        else if (jama)
          hetstat.w[[i]] <- substitute(paste(hl,
                                             chi[df]^2, hq,
                                             " (",
                                             italic(P), hp,
                                             "), ",
                                             italic(I)^2, hi
                                             ),
                                       list(hl = hetlab,
                                            hi = hetstat.I2.w[i],
                                            hq = hetstat.Q.w[i],
                                            hp = hetstat.pval.Q.w[i],
                                            df = k.w.hetstat[i] - 1))
        else {
          ##
          ## One
          ##
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
          ##
          ## Two
          ##
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
          ##
          ## Three
          ##
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
          ##
          ## Four
          ##
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
          ##
          ## Five
          ##
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
    ##
    TE.w <- c(TE.fixed.w, TE.random.w, rep(NA, 4 * n.by))
    lowTE.w <-
      c(lower.fixed.w, lower.random.w, lower.predict.w,
        rep(NA, 3 * n.by))
    uppTE.w <-
      c(upper.fixed.w, upper.random.w, upper.predict.w,
        rep(NA, 3 * n.by))
    n.harmonic.mean.w <-
      c(n.harmonic.mean.w, n.harmonic.mean.w, n.harmonic.mean.w,
        rep(NA, 3 * n.by))
    t.harmonic.mean.w <-
      c(t.harmonic.mean.w, t.harmonic.mean.w, t.harmonic.mean.w,
        rep(NA, 3 * n.by))
    weight.w.p <- c(w.fixed.w.p, w.random.w.p, rep(NA, 4 * n.by))
    ##
    ##test.fixed.w <- ""
    ##
    ## Label of test for effect in subgroups
    ##
    if (test.effect.subgroup.fixed | test.effect.subgroup.random) {
      pvals.effect.w <- formatPT(c(pval.fixed.w, pval.random.w),
                                 lab = TRUE, labval = "",
                                 digits = digits.pval,
                                 zero = zero.pval, JAMA = JAMA.pval,
                                 scientific = scientific.pval,
                                 lab.NA = "NA")
      statistics.effect.w <- formatN(round(c(statistic.fixed.w,
                                             statistic.random.w),
                                           digits = digits.stat),
                                     digits.stat, "NA", big.mark = big.mark)
      ##
      ## Remove superfluous spaces
      ##
      pvals.effect.w <- rmSpace(pvals.effect.w, end = TRUE)
      ##
      while(any(grepl("  ", pvals.effect.w)))
        pvals.effect.w <- gsub("  ", " ", pvals.effect.w)
      while(any(grepl("  ", statistics.effect.w)))
        statistics.effect.w <- gsub("  ", " ", statistics.effect.w)
    }
    ##
    if (test.effect.subgroup.fixed) {
      text.effect.subgroup.fixed <- vector("list", n.by)
      for (i in 1:n.by) {
        if (print.stat) {
          if (revman5)
            text.effect.subgroup.fixed[[i]] <-
              substitute(paste(tl,
                               Z, hetseparator, tt,
                               " (P", tp, ")"),
                         list(tl = label.test.effect.subgroup.fixed,
                              hetseparator = hetseparator,
                              tt = statistics.effect.w[i],
                              tp = rmSpace(pvals.effect.w[i], end = TRUE)))
          else if (jama)
            text.effect.subgroup.fixed[[i]] <-
              substitute(paste(tl,
                               italic(z), hetseparator, tt,
                               " (", italic(P), tp, ")"),
                         list(tl = label.test.effect.subgroup.fixed,
                              hetseparator = hetseparator,
                              tt = statistics.effect.w[i],
                              tp = rmSpace(pvals.effect.w[i], end = TRUE)))
          else
            text.effect.subgroup.fixed[[i]] <-
              substitute(paste(tl,
                               italic(z), hetseparator, tt,
                               " (", italic(p), tp, ")"),
                         list(tl = label.test.effect.subgroup.fixed,
                              hetseparator = hetseparator,
                              tt = statistics.effect.w[i],
                              tp = rmSpace(pvals.effect.w[i], end = TRUE)))
        }
        else {
          if (revman5)
            text.effect.subgroup.fixed[[i]] <-
              substitute(paste(tl,
                               " P", tp),
                         list(tl = label.test.effect.subgroup.fixed,
                              tp = rmSpace(pvals.effect.w[i], end = TRUE)))
          else if (jama)
            text.effect.subgroup.fixed[[i]] <-
              substitute(paste(tl,
                               " ", italic(P), tp),
                         list(tl = label.test.effect.subgroup.fixed,
                              tp = rmSpace(pvals.effect.w[i], end = TRUE)))
          else
            text.effect.subgroup.fixed[[i]] <-
              substitute(paste(tl,
                               " ", italic(p), tp),
                         list(tl = label.test.effect.subgroup.fixed,
                              tp = rmSpace(pvals.effect.w[i], end = TRUE)))
        }
      }
    }
    else {
      text.effect.subgroup.fixed <- vector("list", n.by)
      for (i in 1:n.by)
        text.effect.subgroup.fixed[[i]] <- ""
    }
    ##
    if (test.effect.subgroup.random) {
      text.effect.subgroup.random <- vector("list", n.by)
      for (i in 1:n.by) {
        if (print.stat) {
          if (!x$hakn) {
            if (revman5)
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
            if (revman5)
              text.effect.subgroup.random[[i]] <-
                substitute(paste(tl,
                                 t[df], hetseparator, tt,
                                 " (P", tp, ")"),
                           list(tl = label.test.effect.subgroup.random,
                                hetseparator = hetseparator,
                                tt = statistics.effect.w[n.by + i],
                                tp = rmSpace(pvals.effect.w[n.by + i],
                                             end = TRUE),
                                df = k.w.hetstat - 1))
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
                                df = k.w.hetstat - 1))
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
                                df = k.w.hetstat - 1))
          }
        }
        else {
          if (revman5)
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
      for (i in 1:n.by)
        text.effect.subgroup.random[[i]] <- ""
    }
  }
  
  
  ##
  ##
  ## (8) Backtransform data
  ##
  ##
  TE.orig <- TE
  ##
  if (backtransf) {
    ##
    ## Freeman-Tukey Arcsin transformation
    ##
    if (metainf.metacum) {
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
    ##
    ## Individual study results
    ##
    if (metaprop) {
      TE <- x$event.e / x$n.e
    }
    ## Relative effect measures will be back transformed later
    else if (!is.relative.effect(sm)) {
      TE <- backtransf(TE, sm, "mean", npft)
      lowTE <- backtransf(lowTE, sm, "lower", npft)
      uppTE <- backtransf(uppTE, sm, "upper", npft)
    }
    ##
    ## Results of meta-analysis
    ##
    if (!is.relative.effect(sm)) {
      TE.fixed    <- backtransf(TE.fixed, sm, "mean",
                                npft.ma, warn = comb.fixed)
      lowTE.fixed <- backtransf(lowTE.fixed, sm, "lower",
                                npft.ma, warn = comb.fixed)
      uppTE.fixed <- backtransf(uppTE.fixed, sm, "upper",
                                npft.ma, warn = comb.fixed)
      ##
      TE.random <- backtransf(TE.random, sm, "mean",
                              npft.ma, warn = comb.random)
      lowTE.random <- backtransf(lowTE.random, sm, "lower",
                                 npft.ma, warn = comb.random)
      uppTE.random <- backtransf(uppTE.random, sm, "upper",
                                 npft.ma, warn = comb.random)
      ##
      if (!metainf.metacum) {
        lowTE.predict <- backtransf(lowTE.predict, sm, "lower",
                                    npft.ma, warn = prediction)
        uppTE.predict <- backtransf(uppTE.predict, sm, "upper",
                                    npft.ma, warn = prediction)
      }
      ##
      if (by) {
        if (sm == "IRFT")
          npft.w <- t.harmonic.mean.w
        else
          npft.w <- n.harmonic.mean.w
        ##
        TE.w     <- backtransf(TE.w, sm, "mean", npft.w)
        lowTE.w  <- backtransf(lowTE.w, sm, "lower", npft.w)
        uppTE.w  <- backtransf(uppTE.w, sm, "upper", npft.w)
      }
    }
    ##
    ## Apply argument 'pscale' to proportions / risk differences
    ##
    if (is.prop(sm) | sm == "RD") {
      TE <- pscale * TE
      lowTE <- pscale * lowTE
      uppTE <- pscale * uppTE
      ##
      TE.fixed    <- pscale * TE.fixed
      lowTE.fixed <- pscale * lowTE.fixed
      uppTE.fixed <- pscale * uppTE.fixed
      ##
      TE.random    <- pscale * TE.random
      lowTE.random <- pscale * lowTE.random
      uppTE.random <- pscale * uppTE.random
      ##
      lowTE.predict <- pscale * lowTE.predict
      uppTE.predict <- pscale * uppTE.predict
      ##
      if (by) {
        TE.w    <- pscale * TE.w
        lowTE.w <- pscale * lowTE.w
        uppTE.w <- pscale * uppTE.w
      }
    }
  }
  ##
  ## Apply argument 'irscale' to rates / incidence rate differences
  ##
  if (is.rate(sm) | sm == "IRD") {
    TE <- irscale * TE
    lowTE <- irscale * lowTE
    uppTE <- irscale * uppTE
    ##
    TE.fixed    <- irscale * TE.fixed
    lowTE.fixed <- irscale * lowTE.fixed
    uppTE.fixed <- irscale * uppTE.fixed
    ##
    TE.random    <- irscale * TE.random
    lowTE.random <- irscale * lowTE.random
    uppTE.random <- irscale * uppTE.random
    ##
    lowTE.predict <- irscale * lowTE.predict
    uppTE.predict <- irscale * uppTE.predict
    ##
    if (by) {
      TE.w    <- irscale * TE.w
      lowTE.w <- irscale * lowTE.w
      uppTE.w <- irscale * uppTE.w
    }
  }
  ##
  ## Exclude study results from forest plot
  ##
  TE.exclude <- TE
  lowTE.exclude <- lowTE
  uppTE.exclude <- uppTE
  ##
  if (!null.exclude) {
    TE.exclude[exclude] <- NA
    lowTE.exclude[exclude] <- NA
    uppTE.exclude[exclude] <- NA
  }
  ##
  if (!comb.fixed) {
    TE.fixed    <- NA
    lowTE.fixed <- NA
    uppTE.fixed <- NA
  }
  ##
  if (!comb.random) {
    TE.random    <- NA
    lowTE.random <- NA
    uppTE.random <- NA
  }
  ##
  if (!prediction) {
    lowTE.predict <- NA
    uppTE.predict <- NA
  }
  ##
  if (!metainf.metacum) {
    if (by & !overall)
      w.fixed.p <- round(100 * x$w.fixed, digits.weight)
    else {
      if (!all(is.na(x$w.fixed)) && sum(x$w.fixed) > 0)
        w.fixed.p <-
          round(100 * x$w.fixed / sum(x$w.fixed, na.rm = TRUE),
                digits.weight)
      else
        w.fixed.p <- x$w.fixed
    }
    ##
    if (by & !overall)
      w.random.p <- round(100 * x$w.random, digits.weight)
    else {
      if (!all(is.na(x$w.random)) && sum(x$w.random) > 0)
        w.random.p <-
          round(100 * x$w.random / sum(x$w.random, na.rm = TRUE),
                digits.weight)
      else
        w.random.p <- x$w.random
    }
  }
  else {
    w.fixed.p  <- rep(NA, length(TE))
    w.random.p <- rep(NA, length(TE))
  }
  
  
  ##
  ##
  ## (9) Determine column labels
  ##
  ##
  labs <- list()
  ##
  if (missing(leftlabs) || length(leftcols) != length(leftlabs)) {
    for (i in seq(along = leftcols)) {
      j <- match(leftcols[i], colnames)
      if (!is.na(j))
        labs[[paste0("lab.", leftcols[i])]] <- labnames[j]
    }
  }
  else if (length(leftcols) == length(leftlabs)) {
    for (i in seq(along = leftcols)) {
      j <- match(leftcols[i], colnames)
      if (!is.na(leftlabs[i]))
        labs[[paste0("lab.", leftcols[i])]] <- leftlabs[i]
      else
        if (!is.na(j))
          labs[[paste0("lab.", leftcols[i])]] <- labnames[j]
    }
  }
  ##
  if (missing(rightlabs) || length(rightcols) != length(rightlabs)) {
    for (i in seq(along = rightcols)) {
      j <- match(rightcols[i], colnames)
      if (!is.na(j))
        labs[[paste0("lab.", rightcols[i])]] <- labnames[j]
    }
  }
  else if (length(rightcols) == length(rightlabs)) {
    for (i in seq(along = rightcols)) {
      j <- match(rightcols[i], colnames)
      if (!is.na(rightlabs[i]))
        labs[[paste0("lab.", rightcols[i])]] <- rightlabs[i]
      else
        if (!is.na(j))
          labs[[paste0("lab.", rightcols[i])]] <- labnames[j]
    }
  }
  ##
  if (!slab)
    labs[["lab.studlab"]] <- ""
  ##
  ## Check for "%" in weight labels
  ##
  w.fixed.percent <- TRUE
  if (!is.null(labs[["lab.w.fixed"]])) {
    if (grepl("%", labs[["lab.w.fixed"]]))
      w.fixed.percent <- FALSE
  }
  ##
  w.random.percent <- TRUE
  if (!is.null(labs[["lab.w.random"]])) {
    if (grepl("%", labs[["lab.w.random"]]))
      w.random.percent <- FALSE
  }
  ## "studlab", "TE", "seTE",
  ## "n.e", "n.c",
  ## "event.e", "event.c",
  ## "mean.e", "mean.c",
  ## "sd.e", "sd.c",
  ## "cor",
  ## "time.e", "time.c",
  ##
  ## Check for "\n" in label of column 'studlab'
  ##
  clines <- twolines(labs[["lab.studlab"]], "studlab")
  ##
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
  ##
  ## Check for "\n" in label of column 'effect'
  ##
  clines <- twolines(labs[["lab.effect"]], "effect")
  ##
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
  ##
  ## Check for "\n" in label of column 'ci'
  ##
  clines <- twolines(labs[["lab.ci"]], "ci")
  ##
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
  ##
  ## Check for "\n" in label of column 'effect.ci'
  ##
  clines <- twolines(labs[["lab.effect.ci"]], "effect.ci")
  ##
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
  ##
  ## Check for "\n" in label of column 'w.fixed'
  ##
  clines <- twolines(labs[["lab.w.fixed"]], "w.fixed")
  ##
  if (clines$newline) {
    newline.w.fixed <- TRUE
    labs[["lab.w.fixed"]] <- clines$bottom
    add.w.fixed <- clines$top
    longer.w.fixed <- clines$longer
  }
  else {
    newline.w.fixed <- FALSE
    longer.w.fixed <- labs[["lab.w.fixed"]]
  }
  ##
  ## Check for "\n" in label of column 'w.random'
  ##
  clines <- twolines(labs[["lab.w.random"]], "w.random")
  ##
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
  ##
  ## Check for "\n" in label of column 'TE'
  ##
  clines <- twolines(labs[["lab.TE"]], "TE")
  ##
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
  ##
  ## Check for "\n" in label of column 'seTE'
  ##
  clines <- twolines(labs[["lab.seTE"]], "seTE")
  ##
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
  ##
  ## Check for "\n" in label of column 'n.e'
  ##
  clines <- twolines(labs[["lab.n.e"]], "n.e")
  ##
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
  ##
  ## Check for "\n" in label of column 'n.c'
  ##
  clines <- twolines(labs[["lab.n.c"]], "n.c")
  ##
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
  ##
  ## Check for "\n" in label of column 'event.e'
  ##
  clines <- twolines(labs[["lab.event.e"]], "event.e")
  ##
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
  ##
  ## Check for "\n" in label of column 'event.c'
  ##
  clines <- twolines(labs[["lab.event.c"]], "event.c")
  ##
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
  ##
  ## Check for "\n" in label of column 'mean.e'
  ##
  clines <- twolines(labs[["lab.mean.e"]], "mean.e")
  ##
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
  ##
  ## Check for "\n" in label of column 'mean.c'
  ##
  clines <- twolines(labs[["lab.mean.c"]], "mean.c")
  ##
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
  ##
  ## Check for "\n" in label of column 'sd.e'
  ##
  clines <- twolines(labs[["lab.sd.e"]], "sd.e")
  ##
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
  ##
  ## Check for "\n" in label of column 'sd.c'
  ##
  clines <- twolines(labs[["lab.sd.c"]], "sd.c")
  ##
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
  ##
  ## Check for "\n" in label of column 'cor'
  ##
  clines <- twolines(labs[["lab.cor"]], "cor")
  ##
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
  ##
  ## Check for "\n" in label of column 'time.e'
  ##
  clines <- twolines(labs[["lab.time.e"]], "time.e")
  ##
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
  ##
  ## Check for "\n" in label of column 'time.c'
  ##
  clines <- twolines(labs[["lab.time.c"]], "time.c")
  ##
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
  ##
  ## Check for "\n" in argument 'smlab'
  ##
  clines <- twolines(smlab, arg = TRUE)
  ##
  if (clines$newline) {
    smlab1 <- clines$top
    smlab2 <- clines$bottom
    ##
    newline.smlab <- TRUE
  }
  else {
    smlab1 <- smlab
    smlab2 <- ""
    ##
    newline.smlab <- FALSE
  }
  ##
  ## Check for "\n" in argument 'xlab'
  ##
  clines <- twolines(xlab, arg = TRUE)
  ##
  if (clines$newline) {
    newline.xlab <- TRUE
    xlab <- clines$top
    xlab.add <- clines$bottom
  }
  else {
    xlab.add <- ""
    newline.xlab <- FALSE
  }
  ##
  ## Check for "\n" in argument 'label.left'
  ##
  clines <- twolines(label.left, arg = TRUE)
  ##
  if (clines$newline) {
    ll1 <- clines$top
    ll2 <- clines$bottom
    ##
    newline.ll <- TRUE
  }
  else {
    ll1 <- label.left
    ll2 <- ""
    ##
    newline.ll <- FALSE
  }
  ##
  ## Check for "\n" in argument 'label.right'
  ##
  clines <- twolines(label.right, arg = TRUE)
  ##
  if (clines$newline) {
    lr1 <- clines$top
    lr2 <- clines$bottom
    ##
    newline.lr <- TRUE
  }
  else {
    lr1 <- label.right
    lr2 <- ""
    ##
    newline.lr <- FALSE
  }
  ##
  ## Check for "\n" in additional columns
  ##
  newline.addcol.left  <- FALSE
  newline.addcol.right <- FALSE
  ##
  if (newcols) {
    if (length(leftcols.new) > 0) {
      for (i in seq(along = leftcols.new)) {
        ## Check for "\n" in label of new column
        clines <- twolines(leftlabs.new[i], leftcols.new[i])
        newline.addcol.left <- c(newline.addcol.left, clines$newline)
      }
      newline.addcol.right <- sum(newline.addcol.right) > 0
    }
    if (length(leftcols.new) > 0) {
      for (i in seq(along = leftcols.new)) {
        ## Check for "\n" in label of new column
        clines <- twolines(leftlabs.new[i], leftcols.new[i])
        newline.addcol.left <- c(newline.addcol.left, clines$newline)
      }
      newline.addcol.left <- sum(newline.addcol.left) > 0
    }
  }
  ##
  newline <- newline.studlab | newline.effect | newline.ci | newline.effect.ci |
    newline.w.fixed | newline.w.random | newline.TE | newline.seTE |
    newline.n.e | newline.n.c | newline.event.e | newline.event.c |
    newline.mean.e | newline.mean.c | newline.sd.e | newline.sd.c |
    newline.cor | newline.time.e | newline.time.c |
    newline.smlab | newline.addcol.left | newline.addcol.right
  ##
  newline.all <- newline | (!newline & (newline.ll | newline.lr) & !addrow)
  
  
  ##
  ##
  ## (10) Define columns in forest plot as well as x- and y-limits
  ##
  ##
  if (by) {
    ##
    bylab <- bylabel(bylab, bylevs, print.byvar, byseparator,
                     big.mark = big.mark)
    ##
    wrong.fixed <- FALSE
    wrong.random <- FALSE
    wrong.predict <- FALSE
    ##
    if (n.by > 1) {
      if (length(text.fixed.w) == 1)
        text.fixed.w <- rep(text.fixed.w, n.by)
      else if (length(text.fixed.w) != n.by)
        wrong.fixed <- TRUE
      ##
      if (length(text.random.w) == 1)
        text.random.w <- rep(text.random.w, n.by)
      else if (length(text.random.w) != n.by)
        wrong.random <- TRUE
      ##
      if (length(text.predict.w) == 1)
        text.predict.w <- rep(text.predict.w, n.by)
      else if (length(text.predict.w) != n.by)
        wrong.predict <- TRUE
    }
    else if (n.by == 1) {
      if (length(text.fixed.w) != 1)
        wrong.fixed <- TRUE
      if (length(text.random.w) != 1)
        wrong.random <- TRUE
      if (length(text.predict.w) != 1)
        wrong.predict <- TRUE
    }
    ##
    if (wrong.fixed)
      stop("Argument 'text.fixed.w' must be a single character string",
           if (n.by > 1)
             paste0(" or of length ", n.by, "\n  (number of subgroups).")
           else
             ".")
    ##
    if (wrong.random)
      stop("Argument 'text.random.w' must be a single character string",
           if (n.by > 1)
             paste0(" or of length ", n.by, "\n  (number of subgroups).")
           else
             ".")
    ##
    if (wrong.predict)
      stop("Argument 'text.predict.w' must be a single character string",
           if (n.by > 1)
             paste0(" or of length ", n.by, "\n  (number of subgroups).")
           else
             ".")
    ##
    if (hetstat.pooled == "fixed") {
      text.fixed <- hetstat.overall
      text.fixed.w <- hetstat.w
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
    ##
    modlabs <- c(text.fixed, text.random, text.predict,
                 hetstat.overall, hetstat.resid,
                 text.overall.fixed, text.overall.random,
                 text.subgroup.fixed, text.subgroup.random,
                 text.addline1, text.addline2,
                 bylab, text.fixed.w, text.random.w, text.predict.w, hetstat.w,
                 unlist(text.effect.subgroup.fixed),
                 unlist(text.effect.subgroup.random),
                 studlab)
    ##
    TEs    <- c(TE.fixed, TE.random, NA, TE.w, TE)
    lowTEs <- c(lowTE.fixed, lowTE.random, lowTE.predict, lowTE.w, lowTE)
    uppTEs <- c(uppTE.fixed, uppTE.random, uppTE.predict, uppTE.w, uppTE)
    ##
    TEs.exclude    <- c(TE.fixed, TE.random, NA, TE.w, TE.exclude)
    lowTEs.exclude <- c(lowTE.fixed, lowTE.random, lowTE.predict, lowTE.w,
                        lowTE.exclude)
    uppTEs.exclude <- c(uppTE.fixed, uppTE.random, uppTE.predict, uppTE.w,
                        uppTE.exclude)
    ##
    TEs.study <- c("", "", "", rep("", 6 * n.by),
                   formatN(round(TE.orig, digits), digits, lab.NA,
                           big.mark = big.mark))
    seTEs.study <- c("", "", "", rep("", 6 * n.by),
                     formatN(round(seTE, digits.se), digits.se, lab.NA,
                             big.mark = big.mark))
    ##
    if (weight.subgroup == "same") {
      w.fixeds  <- c(NA, NA, NA, rep(NA, length(weight.w.p)), w.fixed.p)
      w.randoms <- c(NA, NA, NA, rep(NA, length(weight.w.p)), w.random.p)
    }
    else {
      w.fixeds  <- c(NA, NA, NA, weight.w.p, w.fixed.p)
      w.randoms <- c(NA, NA, NA, weight.w.p, w.random.p)
    }
    ##
    format.w.fixed  <- formatN(c(100, weight.w.p, w.fixed.p), digits.weight)
    format.w.random <- formatN(c(100, weight.w.p, w.random.p), digits.weight)
    w.fixeds.text  <- c(format.w.fixed[1], lab.NA.weight, "",
                        format.w.fixed[-1])
    w.randoms.text <- c(lab.NA.weight, format.w.random[1], "",
                        format.w.random[-1])
    ##
    sel.fixed  <- w.fixeds.text == lab.NA.weight
    sel.random <- w.randoms.text == lab.NA.weight
    ##
    sel.fixed[sel.by.random] <- TRUE
    sel.random[sel.by.fixed] <- TRUE
    ##
    type.pooled <- c(type.fixed,
                     type.random,
                     "predict",
                     rep(type.subgroup.fixed, n.by),
                     rep(type.subgroup.random, n.by),
                     rep("predict", n.by),
                     rep("", 3 * n.by))
    col.diamond.pooled <-
      c(col.diamond.fixed, col.diamond.random, col.predict,
        rep(col.diamond.fixed, n.by),
        rep(col.diamond.random, n.by),
        rep(col.predict, n.by),
        rep(NA, 3 * n.by))
    col.diamond.lines.pooled <-
      c(col.diamond.lines.fixed, col.diamond.lines.random,
        col.predict.lines,
        rep(col.diamond.lines.fixed, n.by),
        rep(col.diamond.lines.random, n.by),
        rep(col.predict.lines, n.by),
        rep(NA, 3 * n.by))
    col.inside.pooled <-
      c(col.inside.fixed, col.inside.random, "",
        rep(col.inside.fixed, n.by),
        rep(col.inside.random, n.by),
        rep("", n.by),
        rep(NA, 3 * n.by))
  }
  else {
    ##
    if (hetstat.pooled == "fixed")
      text.fixed <- hetstat.overall
    else if (hetstat.pooled == "random")
      text.random <- hetstat.overall
    ##
    modlabs <- c(text.fixed, text.random, text.predict,
                 hetstat.overall, hetstat.resid,
                 text.overall.fixed, text.overall.random,
                 text.subgroup.fixed, text.subgroup.random,
                 text.addline1, text.addline2,
                 studlab)
    ##
    TEs    <- c(TE.fixed, TE.random, NA, TE)
    lowTEs <- c(lowTE.fixed, lowTE.random, lowTE.predict, lowTE)
    uppTEs <- c(uppTE.fixed, uppTE.random, uppTE.predict, uppTE)
    ##
    TEs.exclude    <- c(TE.fixed, TE.random, NA, TE.exclude)
    lowTEs.exclude <- c(lowTE.fixed, lowTE.random,
                        lowTE.predict, lowTE.exclude)
    uppTEs.exclude <- c(uppTE.fixed, uppTE.random,
                        uppTE.predict, uppTE.exclude)
    ##
    TEs.study <- c("", "", "",
                   formatN(round(TE.orig, digits), digits, lab.NA,
                           big.mark = big.mark))
    seTEs.study <- c("", "", "",
                     formatN(round(seTE, digits.se), digits.se, lab.NA,
                             big.mark = big.mark))
    ##
    w.fixeds  <- c(NA, NA, NA, w.fixed.p)
    w.randoms <- c(NA, NA, NA, w.random.p)
    ##
    format.w.fixed  <- formatN(c(100, w.fixed.p), digits.weight)
    format.w.random <- formatN(c(100, w.random.p), digits.weight)
    w.fixeds.text  <- c(format.w.fixed[1], lab.NA.weight, "",
                        format.w.fixed[-1])
    w.randoms.text <- c(lab.NA.weight, format.w.random[1], "",
                        format.w.random[-1])
    ##
    sel.fixed <- w.fixeds.text == lab.NA.weight
    sel.random <- w.randoms.text == lab.NA.weight
    ##
    type.pooled <- c(type.fixed, type.random, "predict")
    col.diamond.pooled <-
      c(col.diamond.fixed, col.diamond.random, col.predict)
    col.diamond.lines.pooled <-
      c(col.diamond.lines.fixed, col.diamond.lines.random,
        col.predict.lines)
    col.inside.pooled <- c(col.inside.fixed, col.inside.random, "")
  }
  ##
  ## Treatment effect and confidence interval
  ##
  if (backtransf & is.relative.effect(sm)) {
    effect.format <- formatN(round(exp(TEs), digits), digits,
                             lab.NA.effect,
                             big.mark = big.mark)
    ci.format <- ifelse(is.na(lowTEs) | is.na(uppTEs), lab.NA.effect,
                        formatCI(format(round(exp(lowTEs), digits),
                                        scientific = FALSE,
                                        big.mark = big.mark),
                                 format(round(exp(uppTEs), digits),
                                        scientific = FALSE,
                                        big.mark = big.mark)))
  }
  else {
    effect.format <- formatN(round(TEs, digits), digits, lab.NA.effect,
                             big.mark = big.mark)
    ci.format <- ifelse(is.na(lowTEs) | is.na(uppTEs), lab.NA.effect,
                        formatCI(format(round(lowTEs, digits),
                                        scientific = FALSE,
                                        big.mark = big.mark),
                                 format(round(uppTEs, digits),
                                        scientific = FALSE,
                                        big.mark = big.mark)))
  }
  effect.ci.format <- paste(effect.format, ci.format)
  ##
  ## No treatment effect for prediction interval
  ##
  effect.format[3] <- ""
  ##
  ## Only print prediction interval if requested
  ##
  if (!prediction) {
    ci.format[3] <- ""
    effect.ci.format[3] <- ""
  }
  ##
  ## No treatment effect and confidence interval in heterogeneity line
  ##
  if (by) {
    effect.format[sel.by.noNA] <- ""
    ci.format[sel.by.noNA] <- ""
    effect.ci.format[sel.by.noNA] <- ""
  }
  ##
  ## Weights of fixed and random effects model
  ##
  w.fixed.format  <- paste0(w.fixeds.text, if (w.fixed.percent) "%")
  w.random.format <- paste0(w.randoms.text, if (w.random.percent) "%")
  ##
  w.fixed.format[w.fixed.format == "%"] <- ""
  w.random.format[w.random.format == "%"] <- ""
  ##
  w.fixed.format[sel.fixed] <- lab.NA.weight
  if (by)
    w.fixed.format[sel.by.noNA] <- ""
  w.random.format[sel.random] <- lab.NA.weight
  if (by)
    w.random.format[sel.by.noNA] <- ""
  ##
  ## Treatment estimate and its standard error
  ##
  TE.format <- TEs.study
  seTE.format <- seTEs.study
  ##
  ## Number of patients, events, and person times
  ##
  if (!is.null(x$n.e.pooled))
    sum.n.e <- x$n.e.pooled
  else
    sum.n.e <- sum(x$n.e, na.rm = TRUE)
  ##
  if (!is.null(x$n.c.pooled))
    sum.n.c <- x$n.c.pooled
  else
    sum.n.c <- sum(x$n.c, na.rm = TRUE)
  ##
  if (!is.null(x$event.e.pooled))
    sum.e.e <- x$event.e.pooled
  else
    sum.e.e <- sum(x$event.e, na.rm = TRUE)
  ##
  if (!is.null(x$event.c.pooled))
    sum.e.c <- x$event.c.pooled
  else
    sum.e.c <- sum(x$event.c, na.rm = TRUE)
  ##
  if (!is.null(x$time.e.pooled))
    sum.t.e <- x$time.e.pooled
  else
    sum.t.e <- sum(x$time.e, na.rm = TRUE)
  ##
  if (!is.null(x$time.c.pooled))
    sum.t.c <- x$time.c.pooled
  else
    sum.t.c <- sum(x$time.c, na.rm = TRUE)
  ##
  if (by) {
    if (pooled.totals) {
      Ne <- c(sum.n.e, sum.n.e, NA, n.e.w, n.e.w, rep(NA, 4 * n.by), x$n.e)
      Nc <- c(sum.n.c, sum.n.c, NA, n.c.w, n.c.w, rep(NA, 4 * n.by), x$n.c)
    }
    else {
      Ne <- c(NA, NA, NA, rep(NA, 6 * n.by), x$n.e)
      Nc <- c(NA, NA, NA, rep(NA, 6 * n.by), x$n.c)
    }
    if (pooled.events) {
      Ee <- c(sum.e.e, sum.e.e, NA, e.e.w, e.e.w, rep(NA, 4 * n.by), x$event.e)
      Ec <- c(sum.e.c, sum.e.c, NA, e.c.w, e.c.w, rep(NA, 4 * n.by), x$event.c)
    }
    else {
      Ee <- c(NA, NA, NA, rep(NA, 6 * n.by), x$event.e)
      Ec <- c(NA, NA, NA, rep(NA, 6 * n.by), x$event.c)
    }
    if (pooled.times) {
      Te <- c(sum.t.e, sum.t.e, NA, t.e.w, t.e.w, rep(NA, 4 * n.by), x$time.e)
      Tc <- c(sum.t.c, sum.t.c, NA, t.c.w, t.c.w, rep(NA, 4 * n.by), x$time.c)
    }
    else {
      Te <- c(NA, NA, NA, rep(NA, 6 * n.by), x$time.e)
      Tc <- c(NA, NA, NA, rep(NA, 6 * n.by), x$time.c)
    }
  }
  else {
    if (pooled.totals) {
      Ne <- c(sum.n.e, sum.n.e, NA, x$n.e)
      Nc <- c(sum.n.c, sum.n.c, NA, x$n.c)
    }
    else {
      Ne <- c(NA, NA, NA, x$n.e)
      Nc <- c(NA, NA, NA, x$n.c)
    }
    if (pooled.events) {
      Ee <- c(sum.e.e, sum.e.e, NA, x$event.e)
      Ec <- c(sum.e.c, sum.e.c, NA, x$event.c)
    }
    else {
      Ee <- c(NA, NA, NA, x$event.e)
      Ec <- c(NA, NA, NA, x$event.c)
    }
    if (pooled.times) {
      Te <- c(sum.t.e, sum.t.e, NA, x$time.e)
      Tc <- c(sum.t.c, sum.t.c, NA, x$time.c)
    }
    else {
      Te <- c(NA, NA, NA, x$time.e)
      Tc <- c(NA, NA, NA, x$time.c)
    }
  }
  ##
  Ne.format <- ifelse(is.na(Ne), lab.NA, format(Ne, scientific = FALSE,
                                                big.mark = big.mark))
  Nc.format <- ifelse(is.na(Nc), lab.NA, format(Nc, scientific = FALSE,
                                                big.mark = big.mark))
  Ee.format <- ifelse(is.na(Ee), lab.NA, format(Ee, scientific = FALSE,
                                                big.mark = big.mark))
  Ec.format <- ifelse(is.na(Ec), lab.NA, format(Ec, scientific = FALSE,
                                                big.mark = big.mark))
  ##
  if (all(is.wholenumber(Te), na.rm = TRUE) & missing.digits.time)
    Te.format <- ifelse(is.na(Te), lab.NA, format(Te, scientific = FALSE,
                                                  big.mark = big.mark))
  else
    Te.format <- formatN(round(Te, digits.time), digits.time, lab.NA,
                         big.mark = big.mark)
  ##
  if (all(is.wholenumber(Tc), na.rm = TRUE) & missing.digits.time)
    Tc.format <- ifelse(is.na(Tc), lab.NA, format(Tc, scientific = FALSE,
                                                  big.mark = big.mark))
  else
    Tc.format <- formatN(round(Tc, digits.time), digits.time, lab.NA,
                         big.mark = big.mark)
  ##
  ## Print nothing in line with prediction interval
  ##
  Ne.format[3] <- Nc.format[3] <- Ee.format[3] <- Ec.format[3] <-
    Te.format[3] <- Tc.format[3] <- ""
  ##
  if (by) {
    ##
    ## Print nothing in lines with heterogeneity results for subgroups
    ##
    Ne.format[sel.by.noNA] <- Nc.format[sel.by.noNA] <- ""
    Ee.format[sel.by.noNA] <- Ec.format[sel.by.noNA] <- ""
    Te.format[sel.by.noNA] <- Tc.format[sel.by.noNA] <- ""
  }
  ##
  if (fixed.random) {
    ##
    ## Print nothing in lines with results for random effects model
    ##
    Ne.format[2] <- Nc.format[2] <- ""
    Ee.format[2] <- Ec.format[2] <- ""
    Te.format[2] <- Tc.format[2] <- ""
    ##
    if (by) {
      Ne.format[sel.by.random] <- Nc.format[sel.by.random] <- ""
      Ee.format[sel.by.random] <- Ec.format[sel.by.random] <- ""
      Te.format[sel.by.random] <- Tc.format[sel.by.random] <- ""
    }
  }
  ##
  ## Only print total number of events if pooled.events is TRUE
  ##
  if (!pooled.events) {
    Ee.format[1:2] <- Ec.format[1:2] <- ""
    ##
    if (by) {
      Ee.format[sel.by.fixed]  <- Ec.format[sel.by.fixed]  <- ""
      Ee.format[sel.by.random] <- Ec.format[sel.by.random] <- ""
    }
  }
  ##
  ## Only print total person times if pooled.times is TRUE
  ##
  if (!pooled.times) {
    Te.format[1:2] <- Tc.format[1:2] <- ""
    ##
    if (by) {
      Te.format[sel.by.fixed]  <- Tc.format[sel.by.fixed]  <- ""
      Te.format[sel.by.random] <- Tc.format[sel.by.random] <- ""
    }
  }
  ##
  ## Mean and standard deviation
  ##
  if (by) {
    Me <- c(NA, NA, NA, rep(NA, 6 * n.by), x$mean.e)
    Mc <- c(NA, NA, NA, rep(NA, 6 * n.by), x$mean.c)
    Se <- c(NA, NA, NA, rep(NA, 6 * n.by), x$sd.e)
    Sc <- c(NA, NA, NA, rep(NA, 6 * n.by), x$sd.c)
  }
  else {
    Me <- c(NA, NA, NA, x$mean.e)
    Mc <- c(NA, NA, NA, x$mean.c)
    Se <- c(NA, NA, NA, x$sd.e)
    Sc <- c(NA, NA, NA, x$sd.c)
  }
  ##
  if (is.null(digits.mean)) {
    Me.format <- ifelse(is.na(Me), lab.NA, format(Me, scientific = FALSE,
                                                  big.mark = big.mark))
    Mc.format <- ifelse(is.na(Mc), lab.NA, format(Mc, scientific = FALSE,
                                                  big.mark = big.mark))
  }
  else {
    Me.format <- formatN(round(Me, digits.mean), digits.mean, lab.NA,
                         big.mark = big.mark)
    Mc.format <- formatN(round(Mc, digits.mean), digits.mean, lab.NA,
                         big.mark = big.mark)
  }
  if (is.null(digits.sd)) {
    Se.format <- ifelse(is.na(Se), lab.NA, format(Se, scientific = FALSE,
                                                  big.mark = big.mark))
    Sc.format <- ifelse(is.na(Sc), lab.NA, format(Sc, scientific = FALSE,
                                                  big.mark = big.mark))
  }
  else {
    Se.format <- formatN(round(Se, digits.sd), digits.sd, lab.NA,
                         big.mark = big.mark)
    Sc.format <- formatN(round(Sc, digits.sd), digits.sd, lab.NA,
                         big.mark = big.mark)
  }
  ##
  ## Print nothing for lines with summary results
  ##
  Me.format[1:3] <- Mc.format[1:3] <- Se.format[1:3] <- Sc.format[1:3] <- ""
  ##
  if (by) {
    Me.format[sel.by] <- Mc.format[sel.by] <- ""
    Se.format[sel.by] <- Sc.format[sel.by] <- ""
  }
  ##
  ## Correlation
  ##
  if (by)
    cor <- c(NA, NA, NA, rep(NA, 6 * n.by), x$cor)
  else
    cor <- c(NA, NA, NA, x$cor)
  ##
  if (is.null(digits.cor))
    cor.format <- ifelse(is.na(cor), lab.NA, format(cor, scientific = FALSE,
                                                    big.mark = big.mark))
  else
    cor.format <- formatN(round(cor, digits.cor), digits.cor, lab.NA,
                          big.mark = big.mark)
  ##
  ## Print nothing for lines with summary results
  ##
  cor.format[1:3] <- ""
  ##
  if (by)
    cor.format[sel.by] <- ""
  ##
  ##
  ## y-axis:
  ##
  ##
  if ((!(metaprop | metacor) &
       (any(rightcols %in% c("n.e", "n.c")) |
        any(leftcols  %in% c("n.e", "n.c")))
  ) |
  (metainc &
   (any(rightcols %in% c("time.e", "time.c")) |
    any(leftcols  %in% c("time.e", "time.c")))
  ) |
  (metacont &
   (any(rightcols %in% c("sd.e", "sd.c")) |
    any(leftcols  %in% c("sd.e", "sd.c")))
  ) |
  (!is.null(lab.e.attach.to.col) & !is.null(lab.e)) |
  (!is.null(lab.c.attach.to.col) & !is.null(lab.c)) |
  newline.all
  ) {
    yHead <- 2
    yHeadadd <- 1
  }
  else {
    yHead <- 1
    yHeadadd <- NA
  }
  ##
  if (!by) {
    N <- n.stud
    if (study.results)
      yTE <- 1:N
    else
      yTE <- rep(NA, N)
  }
  else {
    ##
    j <- 1
    k <- 0
    yBylab <- rep(NA, n.by)
    yTE <- rep(NA, n.stud)
    yTE.w.fixed <- yBylab
    yTE.w.random <- yBylab
    yTE.w.predict <- yBylab
    yTE.w.hetstat <- yBylab
    yTE.w.effect.fixed <- yBylab
    yTE.w.effect.random <- yBylab
    ##
    for (i in 1:n.by) {
      ##
      if (allstudies)
        k.i <- k.all.w[i]
      else
        k.i <- k.w[i]
      ##
      k <- k + k.i
      ##
      if (print.subgroup.labels) {
        yBylab[i] <- j
        j <- j + 1
      }
      ##
      if (study.results) {
        yTE[(k - k.i + 1):k] <- j:(j + k.i - 1)
        j <- j + k.i
      }
      else
        yTE[(k - k.i + 1):k] <- NA
      ##
      ## Fixed effect model
      ##
      if (comb.fixed & subgroup) {
        yTE.w.fixed[i] <- j
        j <- j + 1
      }
      else
        yTE.w.fixed[i] <- NA
      ##
      ## Random effect model
      ##
      if (comb.random & subgroup) {
        yTE.w.random[i] <- j
        j <- j + 1
      }
      else
        yTE.w.random[i] <- NA
      ##
      ## Only pooled totals
      ##
      if (pooled.totals & subgroup &
          !(comb.fixed | comb.random)) {
        yTE.w.fixed[i] <- j
        j <- j + 1
      }
      ##
      ## Prediction interval
      ##
      if (prediction & subgroup) {
        yTE.w.predict[i] <- j
        j <- j + 1
      }
      else
        yTE.w.predict[i] <- NA
      ##
      ## Heterogeneity statistics
      ##
      if (is.character(hetstat) || hetstat) {
        yTE.w.hetstat[i] <- j
        j <- j + 1
      }
      else
        yTE.w.hetstat[i] <- NA
      ##
      ## Test for effect in subgroup (fixed effect)
      ##
      if (test.effect.subgroup.fixed) {
        yTE.w.effect.fixed[i] <- j
        j <- j + 1
      }
      else
        yTE.w.effect.fixed[i] <- NA
      ##
      ## Test for effect in subgroup (random effects)
      ##
      if (test.effect.subgroup.random) {
        yTE.w.effect.random[i] <- j
        j <- j + 1
      }
      else
        yTE.w.effect.random[i] <- NA
      ##
      if (addrow.subgroups)
        j <- j + 1
    }
    ##
    if (!addrow.subgroups)
      j <- j + 1
    ##
    yTE.w <- c(yTE.w.fixed, yTE.w.random, yTE.w.predict, yTE.w.hetstat,
               yTE.w.effect.fixed, yTE.w.effect.random)
  }
  ##
  ##
  ## x-axis:
  ##
  ##
  if (notmiss.xlim && is.numeric(xlim[1]))
    if (is.relative.effect(sm))
      xlim <- log(xlim)
  ##
  if (is.null(xlim)) {
    if (metaprop) {
      xlim <- c(min(c(lowTE, lowTE.predict), na.rm = TRUE),
                max(c(uppTE, uppTE.predict), na.rm = TRUE))
      ##
      if (!is.na(ref) && ref < xlim[1])
        xlim[1] <- ref
      if (!is.na(ref) && ref > xlim[2])
        xlim[2] <- ref
      ##
      if (!is.na(lower.equi) && lower.equi < xlim[1])
        xlim[1] <- lower.equi
      if (!is.na(lower.equi) && lower.equi > xlim[2])
        xlim[2] <- lower.equi
      ##
      if (!is.na(upper.equi) && upper.equi < xlim[1])
        xlim[1] <- upper.equi
      if (!is.na(upper.equi) && upper.equi > xlim[2])
        xlim[2] <- upper.equi
    }
    else {
      sel.low <- is.finite(lowTE)
      sel.upp <- is.finite(uppTE)
      ##
      if (all(!sel.low))
        minTE <- -0.5
      else
        minTE <- min(c(lowTE[sel.low], lowTE.predict), na.rm = TRUE)
      if (all(!sel.upp))
        maxTE <- 0.5
      else
        maxTE <- max(c(uppTE[sel.upp], uppTE.predict), na.rm = TRUE)
      ##
      xlim <- c(minTE, maxTE)
      ##
      if (!is.na(ref) && ref < xlim[1])
        xlim[1] <- ref
      if (!is.na(ref) && ref > xlim[2])
        xlim[2] <- ref
      ##
      if (!is.na(lower.equi) && lower.equi < xlim[1])
        xlim[1] <- lower.equi
      if (!is.na(lower.equi) && lower.equi > xlim[2])
        xlim[2] <- lower.equi
      ##
      if (!is.na(upper.equi) && upper.equi < xlim[1])
        xlim[1] <- upper.equi
      if (!is.na(upper.equi) && upper.equi > xlim[2])
        xlim[2] <- upper.equi
    }
  }
  ##
  symmetric <- FALSE
  ##
  if (!is.null(xlim) && is.character(xlim[1])) {
    ##
    xlim <- setchar(xlim, "symmetric",
                    paste0("should be a numeric vector (min, max) or ",
                           "the character string \"symmetric\""))
    symmetric <- TRUE
    ##
    if (metaprop | metarate) {
      xlim <- c(min(c(lowTE, lowTE.predict), na.rm = TRUE),
                max(c(uppTE, uppTE.predict), na.rm = TRUE))
    }
    else {
      sel.low <- is.finite(lowTE)
      sel.upp <- is.finite(uppTE)
      ##
      if (all(!sel.low))
        minTE <- -0.5
      else
        minTE <- min(c(lowTE[sel.low], lowTE.predict), na.rm = TRUE)
      if (all(!sel.upp))
        maxTE <- 0.5
      else
        maxTE <- max(c(uppTE[sel.upp], uppTE.predict), na.rm = TRUE)
      ##
      if (minTE < 0 & maxTE < 0)
        xlim <- c(minTE, -minTE)
      else if (minTE > 0 & maxTE > 0)
        xlim <- c(-maxTE, maxTE)
      else
        xlim <- c(-max(abs(c(minTE, maxTE))), max(abs(c(minTE, maxTE))))
    }
  }
  ##
  if (!is.na(ref) &&
      round(xlim[2] - ref, 6) == round(ref - xlim[1], 6))
    symmetric <- TRUE
  ##
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
  ##
  yNext <- max.yTE + ifelse(max.yTE == 0 | !addrow.overall, 1, 2)
  ##
  if (missing(xlab.pos))
    xlab.pos <- mean(xlim)
  ##
  if (missing(smlab.pos))
    smlab.pos <- mean(xlim)
  ##
  yTE.fixed  <- NA
  yTE.random <- NA
  yPredict   <- NA
  yHetstat <- NA
  yResidHetstat <- NA
  yOverall.fixed  <- NA
  yOverall.random <- NA
  ySubgroup.fixed  <- NA
  ySubgroup.random <- NA
  yText.addline1 <- NA
  yText.addline2 <- NA
  ##
  if (comb.fixed & comb.random & overall) {
    yTE.fixed  <- yNext
    yTE.random <- yNext + 1
    yNext      <- yNext + 2
  }
  ##
  else if (comb.fixed & !comb.random & overall) {
    yTE.fixed <- yNext
    yNext     <- yNext + 1
  }
  ##
  else if (!comb.fixed & comb.random & overall) {
    yTE.random <- yNext
    yNext      <- yNext + 1
  }
  ##
  else if (!comb.fixed & !comb.random & pooled.totals & overall) {
    yTE.fixed  <- yNext
    yNext      <- yNext + 1
    if (missing(text.fixed))
      text.fixed <- "Overall"
  }
  ##
  if (prediction & overall) {
    yPredict <- yNext
    yNext    <- yNext + 1
  }
  ##
  if (overall.hetstat) {
    yHetstat <- yNext
    yNext    <- yNext + 1
  }
  ##
  if (resid.hetstat) {
    yResidHetstat <- yNext
    yNext  <- yNext + 1
  }
  ##
  if (test.overall.fixed) {
    yOverall.fixed <- yNext
    yNext          <- yNext + 1
  }
  ##
  if (test.overall.random) {
    yOverall.random <- yNext
    yNext           <- yNext + 1
  }
  ##
  if (test.subgroup.fixed) {
    ySubgroup.fixed <- yNext
    yNext           <- yNext + 1
  }
  ##
  if (test.subgroup.random) {
    ySubgroup.random <- yNext
    yNext            <- yNext + 1
  }
  ##
  if (!missing.text.addline1) {
    yText.addline1 <- yNext
    yNext      <- yNext + 1
  }
  ##
  if (!missing.text.addline2)
    yText.addline2 <- yNext
  ##
  if (!comb.fixed & !pooled.totals) text.fixed <- ""
  if (!comb.random) text.random <- ""
  if (!prediction) text.predict <- ""
  ##
  yTE <- yHead + yTE + addrow
  ##
  yTE.fixed  <- yHead + yTE.fixed + addrow
  yTE.random <- yHead + yTE.random + addrow
  yPredict   <- yHead + yPredict + addrow
  ##
  yHetstat <- yHead + yHetstat + addrow
  yResidHetstat <- yHead + yResidHetstat + addrow
  yOverall.fixed  <- yHead + yOverall.fixed + addrow
  yOverall.random <- yHead + yOverall.random + addrow
  ySubgroup.fixed  <- yHead + ySubgroup.fixed + addrow
  ySubgroup.random <- yHead + ySubgroup.random + addrow
  yText.addline1 <- yHead + yText.addline1 + addrow
  yText.addline2 <- yHead + yText.addline2 + addrow
  ##
  yStats <- c(yHetstat,
              yResidHetstat,
              yOverall.fixed, yOverall.random,
              ySubgroup.fixed, ySubgroup.random,
              yText.addline1, yText.addline2)
  ##
  if (by) {
    yBylab <- yHead + yBylab + addrow
    yTE.w  <- yHead + yTE.w + addrow
    ##
    yLab <- c(yHead,
              yTE.fixed, yTE.random, yPredict,
              yStats,
              yBylab, yTE.w,
              yTE)
    ##
    yS <- c(yHead, yTE.fixed, yTE.random, yPredict, yTE.w, yTE)
  }
  else {
    yLab <- c(yHead, yTE.fixed, yTE.random, yPredict,
              yStats,
              yTE)
    ##
    yS <- c(yHead, yTE.fixed, yTE.random, yPredict, yTE)
  }
  
  
  ##
  ##
  ## (11) Format columns in forest plot
  ##
  ##
  col.studlab <- list(labels = 
                        lapply(as.list(c(labs[["lab.studlab"]], modlabs)),
                               tg,
                               xpos = xpos.s, just = just.s,
                               fs = fs.study.labels,
                               ff = ff.study.labels,
                               fontfamily = fontfamily),
                      rows = yLab
                      )
  ## Study label:
  col.studlab$labels[[1]] <- tg(labs[["lab.studlab"]], xpos.s,
                                just.s, fs.head, ff.head, fontfamily)
  ## Fixed effect estimate:
  col.studlab$labels[[2]] <- tg(text.fixed, xpos.s, just.s,
                                fs.fixed.labels, ff.fixed.labels, fontfamily)
  ## Random effects estimate:
  col.studlab$labels[[3]] <- tg(text.random, xpos.s, just.s,
                                fs.random.labels, ff.random.labels, fontfamily)
  ## Prediction interval:
  col.studlab$labels[[4]] <- tg(text.predict, xpos.s, just.s,
                                fs.predict.labels, ff.predict.labels,
                                fontfamily)
  ## Heterogeneity statistics:
  col.studlab$labels[[5]] <- tg(hetstat.overall, xpos.s, just.s,
                                fs.hetstat, ff.hetstat, fontfamily)
  ## Statistic for residual heterogeneity:
  col.studlab$labels[[6]] <- tg(hetstat.resid, xpos.s, just.s,
                                fs.hetstat, ff.hetstat, fontfamily)
  ## Test for overall effect (fixed effect model):
  col.studlab$labels[[7]] <- tg(text.overall.fixed, xpos.s, just.s,
                                fs.test.overall, ff.test.overall, fontfamily)
  ## Test for overall effect (random effects model):
  col.studlab$labels[[8]] <- tg(text.overall.random, xpos.s, just.s,
                                fs.test.overall, ff.test.overall, fontfamily)
  ## Test for subgroup differences (fixed effect model):
  col.studlab$labels[[9]] <- tg(text.subgroup.fixed, xpos.s, just.s,
                                fs.test.subgroup, ff.test.subgroup, fontfamily)
  ## Test for subgroup differences (random effects model):
  col.studlab$labels[[10]] <- tg(text.subgroup.random, xpos.s, just.s,
                                 fs.test.subgroup, ff.test.subgroup, fontfamily)
  ## First additional line:
  col.studlab$labels[[11]] <- tg(text.addline1, xpos.s, just.s,
                                 fs.addline, ff.addline, fontfamily)
  ## Second additional line:
  col.studlab$labels[[12]] <- tg(text.addline2, xpos.s, just.s,
                                 fs.addline, ff.addline, fontfamily)
  ##
  n.summaries <- 12
  ##
  if (by) {
    for (i in 1:n.by) {
      ## Subgroup labels:
      col.studlab$labels[[n.summaries + i]] <-
        tg(bylab[i], xpos.s, just.s,
           fs.head, ff.head, fontfamily, col.by)
      ## Fixed effect estimates:
      col.studlab$labels[[n.summaries + 1 * n.by + i]] <-
        tg(text.fixed.w[[i]], xpos.s, just.s,
           fs.fixed.labels, ff.fixed.labels, fontfamily, col.by)
      ## Random effects estimates:
      col.studlab$labels[[n.summaries + 2 * n.by + i]] <-
        tg(text.random.w[[i]], xpos.s, just.s,
           fs.random.labels, ff.random.labels, fontfamily, col.by)
      ## Prediction interval:
      col.studlab$labels[[n.summaries + 3 * n.by + i]] <-
        tg(text.predict.w[[i]], xpos.s, just.s,
           fs.predict.labels, ff.predict.labels, fontfamily, col.by)
      ## Heterogeneity statistics:
      col.studlab$labels[[n.summaries + 4 * n.by + i]] <-
        tg(hetstat.w[[i]], xpos.s, just.s,
           fs.hetstat, ff.hetstat, fontfamily, col.by)
      ## Test for effect in subgroup (fixed effect model):
      col.studlab$labels[[n.summaries + 5 * n.by + i]] <-
        tg(text.effect.subgroup.fixed[[i]], xpos.s, just.s,
           fs.test.effect.subgroup, ff.test.effect.subgroup,
           fontfamily, col.by)
      ## Test for effect in subgroup (random effects model):
      col.studlab$labels[[n.summaries + 6 * n.by + i]] <-
        tg(text.effect.subgroup.random[[i]], xpos.s, just.s,
           fs.test.effect.subgroup, ff.test.effect.subgroup,
           fontfamily, col.by)
    }
  }
  ##
  fcs <- list(fs.study = fs.study, ff.study = ff.study,
              fs.heading = fs.head, ff.heading = ff.head,
              fs.fixed = fs.fixed, ff.fixed = ff.fixed,
              fs.random = fs.random, ff.random = ff.random,
              fs.predict = fs.predict, ff.predict = ff.predict,
              by = by, n.by = n.by, col.by = col.by)
  ##
  col.effect <- formatcol(labs[["lab.effect"]], effect.format,
                          yS, just.c, fcs, fontfamily)
  ##
  col.ci <- formatcol(labs[["lab.ci"]], ci.format, yS, just.c, fcs, fontfamily)
  ##
  col.effect.ci <-
    formatcol(labs[["lab.effect.ci"]], effect.ci.format, yS,
              if (revman5) "center" else just.c, fcs, fontfamily)
  ##
  col.w.fixed  <- formatcol(labs[["lab.w.fixed"]], w.fixed.format, yS,
                            just.c, fcs, fontfamily)
  col.w.random <- formatcol(labs[["lab.w.random"]], w.random.format, yS,
                            just.c, fcs, fontfamily)
  ##
  col.TE <- formatcol(labs[["lab.TE"]], TE.format, yS, just.c, fcs, fontfamily)
  col.seTE <- formatcol(labs[["lab.seTE"]], seTE.format, yS, just.c, fcs,
                        fontfamily)
  ##
  col.n.e <- formatcol(labs[["lab.n.e"]], Ne.format, yS, just.c, fcs,
                       fontfamily)
  col.n.c <- formatcol(labs[["lab.n.c"]], Nc.format, yS, just.c, fcs,
                       fontfamily)
  ##
  col.event.e <- formatcol(labs[["lab.event.e"]], Ee.format, yS, just.c, fcs,
                           fontfamily)
  col.event.c <- formatcol(labs[["lab.event.c"]], Ec.format, yS, just.c, fcs,
                           fontfamily)
  ##
  col.mean.e <- formatcol(labs[["lab.mean.e"]], Me.format, yS, just.c, fcs,
                          fontfamily)
  col.mean.c <- formatcol(labs[["lab.mean.c"]], Mc.format, yS, just.c, fcs,
                          fontfamily)
  ##
  col.sd.e <- formatcol(labs[["lab.sd.e"]], Se.format, yS, just.c, fcs,
                        fontfamily)
  col.sd.c <- formatcol(labs[["lab.sd.c"]], Sc.format, yS, just.c, fcs,
                        fontfamily)
  ##
  col.cor <- formatcol(labs[["lab.cor"]], cor.format, yS, just.c, fcs,
                       fontfamily)
  ##
  col.time.e <- formatcol(labs[["lab.time.e"]], Te.format, yS, just.c, fcs,
                          fontfamily)
  col.time.c <- formatcol(labs[["lab.time.c"]], Tc.format, yS, just.c, fcs,
                          fontfamily)
  ##
  ##
  ##
  col.effect.calc <- formatcol(longer.effect, effect.format, yS, just.c, fcs,
                               fontfamily)
  ##
  col.ci.calc <- formatcol(longer.ci, ci.format, yS, just.c, fcs,
                           fontfamily)
  ##
  col.effect.ci.calc <- formatcol(longer.effect.ci, effect.ci.format, yS,
                                  just.c, fcs, fontfamily)
  ##
  col.w.fixed.calc  <- formatcol(longer.w.fixed, w.fixed.format, yS,
                                 just.c, fcs, fontfamily)
  col.w.random.calc <- formatcol(longer.w.random, w.random.format, yS,
                                 just.c, fcs, fontfamily)
  ##
  col.TE.calc <- formatcol(longer.TE, TE.format, yS, just.c, fcs, fontfamily)
  col.seTE.calc <- formatcol(longer.seTE, seTE.format, yS, just.c, fcs,
                             fontfamily)
  ##
  col.n.e.calc <- formatcol(longer.n.e, Ne.format, yS, just.c, fcs, fontfamily)
  col.n.c.calc <- formatcol(longer.n.c, Nc.format, yS, just.c, fcs, fontfamily)
  ##
  col.event.e.calc <- formatcol(longer.event.e, Ee.format, yS,
                                just.c, fcs, fontfamily)
  col.event.c.calc <- formatcol(longer.event.c, Ec.format, yS, just.c, fcs,
                                fontfamily)
  ##
  col.mean.e.calc <- formatcol(longer.mean.e, Me.format, yS, just.c, fcs,
                               fontfamily)
  col.mean.c.calc <- formatcol(longer.mean.c, Mc.format, yS, just.c, fcs,
                               fontfamily)
  ##
  col.sd.e.calc <- formatcol(longer.sd.e, Se.format, yS, just.c, fcs,
                             fontfamily)
  col.sd.c.calc <- formatcol(longer.sd.c, Sc.format, yS, just.c, fcs,
                             fontfamily)
  ##
  col.cor.calc <- formatcol(longer.cor, cor.format, yS, just.c, fcs,
                            fontfamily)
  ##
  col.time.e.calc <- formatcol(longer.time.e, Te.format, yS, just.c, fcs,
                               fontfamily)
  col.time.c.calc <- formatcol(longer.time.c, Tc.format, yS, just.c, fcs,
                               fontfamily)
  ##
  ##
  ##
  col.forest <- list(eff = TEs.exclude,
                     low = lowTEs.exclude,
                     upp = uppTEs.exclude,
                     rows = yS[-1],
                     ##
                     ## "square"  - normal confidence interval
                     ## "diamond" - meta-analysis diamond
                     ## "predict" - prediction interval
                     ##
                     type = c(type.pooled, type.study),
                     ##
                     col = c(col.diamond.lines.pooled, col.study),
                     col.square = c(col.diamond.pooled, col.square),
                     col.square.lines = c(col.diamond.lines.pooled, col.square.lines),
                     col.inside = c(col.inside.pooled, col.inside),
                     ##
                     col.diamond = c(col.diamond.pooled, col.square),
                     col.diamond.lines = c(col.diamond.lines.pooled, col.square.lines),
                     ##
                     lwd = lwd
                     )
  ##
  ## Sizes of squares
  ##
  if (weight.study == "same") {
    information <- rep(0.9, length(TEs))
  }
  else {
    ##
    if (weight.study == "fixed")
      information <- sqrt(w.fixeds)
    else if (weight.study == "random")
      information <- sqrt(w.randoms)
    ## Square height equal to 1 for most precise study result
    if (!all(is.na(information)))
      information <- information / max(information, na.rm = TRUE)
    else
      information <- rep(0.9, length(TEs))
    ## Same / maximum polygon height for all meta-analytical results
    ## (both overall and subgroup results)
    information[is.na(information)] <- 1
  }
  ##
  col.forest$sizes <- information
  col.forest$sizes <- col.forest$sizes * squaresize
  ##
  ## Width of column 3
  col.forestwidth <- plotwidth
  ##
  ## Range on the x-axis for column 3
  col.forest$range <- xlim
  ##
  cols <- list(col.studlab = col.studlab,
               col.effect = col.effect,
               col.ci = col.ci,
               col.effect.ci = col.effect.ci,
               col.w.fixed = col.w.fixed,
               col.w.random = col.w.random,
               col.TE = col.TE,
               col.seTE = col.seTE)
  ##
  cols.calc <- list(col.studlab = col.studlab,
                    col.effect = col.effect.calc,
                    col.ci = col.ci.calc,
                    col.effect.ci = col.effect.ci.calc,
                    col.w.fixed = col.w.fixed.calc,
                    col.w.random = col.w.random.calc,
                    col.TE = col.TE.calc,
                    col.seTE = col.seTE.calc)
  ##
  cols[["col.n.e"]] <- col.n.e
  cols[["col.n.c"]] <- col.n.c
  cols[["col.event.e"]] <- col.event.e
  cols[["col.event.c"]] <- col.event.c
  ##
  cols[["col.mean.e"]] <- col.mean.e
  cols[["col.mean.c"]] <- col.mean.c
  cols[["col.sd.e"]] <- col.sd.e
  cols[["col.sd.c"]] <- col.sd.c
  ##
  cols[["col.cor"]] <- col.cor
  ##
  cols[["col.time.e"]] <- col.time.e
  cols[["col.time.c"]] <- col.time.c
  ##
  cols.calc[["col.n.e"]] <- col.n.e.calc
  cols.calc[["col.n.c"]] <- col.n.c.calc
  cols.calc[["col.event.e"]] <- col.event.e.calc
  cols.calc[["col.event.c"]] <- col.event.c.calc
  ##
  cols.calc[["col.mean.e"]] <- col.mean.e.calc
  cols.calc[["col.mean.c"]] <- col.mean.c.calc
  cols.calc[["col.sd.e"]] <- col.sd.e.calc
  cols.calc[["col.sd.c"]] <- col.sd.c.calc
  ##
  cols.calc[["col.cor"]] <- col.cor.calc
  ##
  cols.calc[["col.time.e"]] <- col.time.e.calc
  cols.calc[["col.time.c"]] <- col.time.c.calc
  ##
  if (newcols) {
    ##
    ## Check just.addcols
    ##
    if (length(leftcols.new) > 0)
      if (length(just.addcols.left) != 1) {
        if (length(just.addcols.left) != length(leftcols.new))
          stop("Length of argument 'just.addcols.left' must be one or ",
               "same as number of additional columms in argument 'leftcols'.")
      }
      else
        just.addcols.left <- rep(just.addcols.left, length(leftcols.new))
    ##
    if (length(rightcols.new) > 0)
      if (length(just.addcols.right) != 1) {
        if (length(just.addcols.right) != length(rightcols.new))
          stop("Length of argument 'just.addcols.right' must be one or ",
               "same as number of additional columms in argument 'rightcols'.")
      }
      else
        just.addcols.right <- rep(just.addcols.right, length(rightcols.new))
    ##
    ## Check digits.addcols
    ##
    if (length(leftcols.new) > 0)
      if (length(digits.addcols.left) != 1) {
        if (length(digits.addcols.left) != length(leftcols.new))
          stop("Length of argument 'digits.addcols.left' must be one or ",
               "same as number of additional columms in argument 'leftcols'.")
      }
      else
        digits.addcols.left <- rep(digits.addcols.left, length(leftcols.new))
    ##
    if (length(rightcols.new) > 0)
      if (length(digits.addcols.right) != 1) {
        if (length(digits.addcols.right) != length(rightcols.new))
          stop("Length of argument 'digits.addcols.right' must be one or ",
               "same as number of additional columms in argument 'rightcols'.")
      }
      else
        digits.addcols.right <- rep(digits.addcols.right, length(rightcols.new))
    ##
    if (by) {
      for (i in seq(along = rightcols.new)) {
        tname <- paste0("col.", rightcols.new[i])
        if (length(dataset1[[rightcols.new[i]]]) != 0)
          tmp.r <- dataset1[[rightcols.new[i]]]
        else if (length(dataset2[[rightcols.new[i]]]) != 0)
          tmp.r <- dataset2[[rightcols.new[i]]]
        ##
        if (!is.character(tmp.r)) {
          if (is.factor(tmp.r))
            tmp.r <- as.character(tmp.r)
          else if (all(is.wholenumber(tmp.r), na.rm = TRUE))
            tmp.r <- ifelse(is.na(tmp.r), lab.NA,
                            format(tmp.r, scientific = FALSE,
                                   big.mark = big.mark))
          else if (is.numeric(tmp.r)) {
            if (rightcols.new[i] == "pval")
              tmp.r <- formatPT(tmp.r, digits = digits.pval,
                                big.mark = big.mark)
            else
              tmp.r <- formatN(tmp.r, digits = digits.addcols.right[i],
                               text.NA = "", big.mark = big.mark)
          }
        }
        ##
        tmp.r <- ifelse(is.na(tmp.r), lab.NA, tmp.r)
        ##
        ## Check for "\n" in label of new column
        ##
        clines <- twolines(rightlabs.new[i], rightcols.new[i])
        ##
        if (clines$newline) {
          lab.new <- clines$bottom
          longer.new <- clines$longer
        }
        else
          lab.new <- longer.new <- rightlabs.new[i]
        cols[[tname]] <- formatcol(lab.new,
                                   c("", "", "", rep("", length(TE.w)), tmp.r),
                                   yS,
                                   if (rightcols.new[i] == "pval") just
                                   else just.addcols.right[i],
                                   fcs, fontfamily)
        cols.calc[[tname]] <- formatcol(longer.new,
                                        c("", "", "", rep("", length(TE.w)), tmp.r),
                                        yS,
                                        just.addcols.right[i],
                                        fcs, fontfamily)
      }
      for (i in seq(along = leftcols.new)) {
        tname <- paste0("col.", leftcols.new[i])
        if (length(dataset1[[leftcols.new[i]]]) != 0)
          tmp.l <- dataset1[[leftcols.new[i]]]        
        else if (length(dataset2[[leftcols.new[i]]]) != 0)
          tmp.l <- dataset2[[leftcols.new[i]]]
        ##
        if (!is.character(tmp.l)) {
          if (is.factor(tmp.l))
            tmp.l <- as.character(tmp.l)
          else if (all(is.wholenumber(tmp.l), na.rm = TRUE))
            tmp.l <- ifelse(is.na(tmp.l), lab.NA,
                            format(tmp.l, scientific = FALSE,
                                   big.mark = big.mark))
          else if (is.numeric(tmp.l)) {
            if (leftcols.new[i] == "pval")
              tmp.l <- formatPT(tmp.l, digits = digits.pval,
                                big.mark = big.mark)
            else
              tmp.l <- formatN(tmp.l, digits = digits.addcols.left[i],
                               text.NA = "", big.mark = big.mark)
          }
        }
        ##
        tmp.l <- ifelse(is.na(tmp.l), lab.NA, tmp.l)
        ##
        ## Check for "\n" in label of new column
        ##
        clines <- twolines(leftlabs.new[i], leftcols.new[i])
        ##
        if (clines$newline) {
          lab.new <- clines$bottom
          longer.new <- clines$longer
        }
        else
          lab.new <- longer.new <- leftlabs.new[i]
        ##
        cols[[tname]] <- formatcol(lab.new,
                                   c("", "", "",
                                     rep("", length(TE.w)), tmp.l),
                                   yS,
                                   if (leftcols.new[i] == "pval") just
                                   else just.addcols.left[i],
                                   fcs, fontfamily)
        ##
        cols.calc[[tname]] <- formatcol(longer.new,
                                        c("", "", "",
                                          rep("", length(TE.w)), tmp.l),
                                        yS,
                                        just.addcols.left[i],
                                        fcs, fontfamily)
      }
    }
    else {
      for (i in seq(along = rightcols.new)) {
        tname <- paste0("col.", rightcols.new[i])
        if (length(dataset1[[rightcols.new[i]]]) != 0)
          tmp.r <- dataset1[[rightcols.new[i]]]
        else if (length(dataset2[[rightcols.new[i]]]) != 0)
          tmp.r <- dataset2[[rightcols.new[i]]]
        ##
        if (!is.character(tmp.r)) {
          if (is.factor(tmp.r))
            tmp.r <- as.character(tmp.r)
          else if (all(is.wholenumber(tmp.r), na.rm = TRUE))
            tmp.r <- ifelse(is.na(tmp.r), lab.NA,
                            format(tmp.r, scientific = FALSE,
                                   big.mark = big.mark))
          else if (is.numeric(tmp.r)) {
            if (rightcols.new[i] == "pval")
              tmp.r <- formatPT(tmp.r, digits = digits.pval,
                                big.mark = big.mark,
                                lab = FALSE, labval = "",
                                zero = zero.pval, JAMA = JAMA.pval,
                                scientific = scientific.pval,
                                lab.NA = "NA")
            else
              tmp.r <- formatN(tmp.r, digits = digits.addcols.right[i],
                               text.NA = "", big.mark = big.mark)
          }
        }
        ##
        tmp.r <- ifelse(is.na(tmp.r), "", tmp.r)
        ##
        ## Check for "\n" in label of new column
        ##
        clines <- twolines(rightlabs.new[i], rightcols.new[i])
        ##
        if (clines$newline) {
          lab.new <- clines$bottom
          longer.new <- clines$longer
        }
        else
          lab.new <- longer.new <- rightlabs.new[i]
        ##
        cols[[tname]] <- formatcol(lab.new,
                                   c("", "", "", tmp.r),
                                   yS,
                                   if (rightcols.new[i] == "pval") just
                                   else just.addcols.right[i],
                                   fcs, fontfamily)
        cols.calc[[tname]] <- formatcol(longer.new,
                                        c("", "", "", tmp.r),
                                        yS,
                                        just.addcols.right[i],
                                        fcs, fontfamily)
      }
      for (i in seq(along = leftcols.new)) {
        tname <- paste0("col.", leftcols.new[i])
        if (length(dataset1[[leftcols.new[i]]]) != 0)
          tmp.l <- dataset1[[leftcols.new[i]]]        
        else if (length(dataset2[[leftcols.new[i]]]) != 0)
          tmp.l <- dataset2[[leftcols.new[i]]]
        ##
        if (!is.character(tmp.l)) {
          if (is.factor(tmp.l))
            tmp.l <- as.character(tmp.l)
          else if (all(is.wholenumber(tmp.l), na.rm = TRUE))
            tmp.l <- ifelse(is.na(tmp.l), lab.NA,
                            format(tmp.l, scientific = FALSE,
                                   big.mark = big.mark))
          else if (is.numeric(tmp.l)) {
            if (leftcols.new[i] == "pval")
              tmp.l <- formatPT(tmp.l, digits = digits.pval,
                                big.mark = big.mark,
                                lab = FALSE, labval = "",
                                zero = zero.pval, JAMA = JAMA.pval,
                                scientific = scientific.pval,
                                lab.NA = "NA")
            else
              tmp.l <- formatN(tmp.l, digits = digits.addcols.left[i],
                               text.NA = "", big.mark = big.mark)
          }
        }
        ##
        tmp.l <- ifelse(is.na(tmp.l), "", tmp.l)
        ##
        ## Check for "\n" in label of new column
        ##
        clines <- twolines(leftlabs.new[i], leftcols.new[i])
        ##
        if (clines$newline) {
          lab.new <- clines$bottom
          longer.new <- clines$longer
        }
        else
          lab.new <- longer.new <- leftlabs.new[i]
        ##
        cols[[tname]] <- formatcol(lab.new,
                                   c("", "", "", tmp.l),
                                   yS,
                                   if (leftcols.new[i] == "pval") just
                                   else just.addcols.left[i],
                                   fcs, fontfamily)
        cols.calc[[tname]] <- formatcol(longer.new,
                                        c("", "", "", tmp.l),
                                        yS,
                                        just.addcols.left[i],
                                        fcs, fontfamily)
      }
    }
  }
  ##
  col.lab.e <- tgl(lab.e, xpos.c, just.c, fs.head, ff.head, fontfamily)
  ##
  col.lab.c <- tgl(lab.c, xpos.c, just.c, fs.head, ff.head, fontfamily)
  ##
  ##
  ##
  if (newline.studlab)
    col.add.studlab <- tgl(add.studlab, xpos.s, just.s, fs.head, ff.head,
                           fontfamily)
  ##
  if (newline.effect)
    col.add.effect <- tgl(add.effect, xpos.c, just.c, fs.head, ff.head,
                          fontfamily)
  ##
  if (newline.ci)
    col.add.ci <- tgl(add.ci, xpos.c, just.c, fs.head, ff.head, fontfamily)
  ##
  if (newline.effect.ci)
    col.add.effect.ci <- tgl(add.effect.ci,
                             if (revman5) 0.5 else xpos.c,
                             if (revman5) "center" else just.c,
                             fs.head, ff.head, fontfamily)
  ##
  if (newline.w.fixed)
    col.add.w.fixed <- tgl(add.w.fixed, xpos.c, just.c, fs.head, ff.head,
                           fontfamily)
  ##
  if (newline.w.random)
    col.add.w.random <- tgl(add.w.random, xpos.c, just.c, fs.head, ff.head,
                            fontfamily)
  ##
  if (newline.TE)
    col.add.TE <- tgl(add.TE, xpos.c, just.c, fs.head, ff.head,
                      fontfamily)
  ##
  if (newline.seTE)
    col.add.seTE <- tgl(add.seTE, xpos.c, just.c, fs.head, ff.head, fontfamily)
  ##
  if (newline.n.e)
    col.add.n.e <- tgl(add.n.e, xpos.c, just.c, fs.head, ff.head, fontfamily)
  ##
  if (newline.n.c)
    col.add.n.c <- tgl(add.n.c, xpos.c, just.c, fs.head, ff.head, fontfamily)
  ##
  if (newline.event.e)
    col.add.event.e <- tgl(add.event.e, xpos.c, just.c, fs.head, ff.head,
                           fontfamily)
  ##
  if (newline.event.c)
    col.add.event.c <- tgl(add.event.c, xpos.c, just.c, fs.head, ff.head,
                           fontfamily)
  ##
  if (newline.mean.e)
    col.add.mean.e <- tgl(add.mean.e, xpos.c, just.c, fs.head, ff.head,
                          fontfamily)
  ##
  if (newline.mean.c)
    col.add.mean.c <- tgl(add.mean.c, xpos.c, just.c, fs.head, ff.head,
                          fontfamily)
  ##
  if (newline.sd.e)
    col.add.sd.e <- tgl(add.sd.e, xpos.c, just.c, fs.head, ff.head, fontfamily)
  ##
  if (newline.sd.c)
    col.add.sd.c <- tgl(add.sd.c, xpos.c, just.c, fs.head, ff.head, fontfamily)
  ##
  if (newline.cor)
    col.add.cor <- tgl(add.cor, xpos.c, just.c, fs.head, ff.head, fontfamily)
  ##
  if (newline.time.e)
    col.add.time.e <- tgl(add.time.e, xpos.c, just.c, fs.head, ff.head,
                          fontfamily)
  ##
  if (newline.time.c)
    col.add.time.c <- tgl(add.time.c, xpos.c, just.c, fs.head, ff.head,
                          fontfamily)
  ##
  leftcols  <- paste0("col.", leftcols)
  rightcols <- paste0("col.", rightcols)
  
  
  ##
  ##
  ## (12) Calculate width of columns in forest plot
  ##
  ##
  ## Exclude lines with summary measures from calculation of column
  ## width for study labels
  ##
  if (by) {
    del.lines <-
      c(if (!calcwidth.fixed)              # FE
          2,
        if (!calcwidth.random)             # RE
          3,
        if (!calcwidth.predict)            # PI
          4,
        if (!calcwidth.tests)              # tests
          5:n.summaries,
        ##
        ## Subgroups
        ##
        if (!calcwidth.subgroup)
          n.summaries + 0 * n.by + 1:n.by, # Labels
        if (!calcwidth.fixed)              # FE
          n.summaries + 1 * n.by + 1:n.by,
        if (!calcwidth.random)             # RE
          n.summaries + 2 * n.by + 1:n.by,
        if (!calcwidth.predict)            # PI
          n.summaries + 3 * n.by + 1:n.by,
        if (!calcwidth.hetstat)            # heterogeneity statistic
          n.summaries + 4 * n.by + 1:n.by,
        if (!calcwidth.tests)              # test for effect (FE)
          n.summaries + 5 * n.by + 1:n.by,
        if (!calcwidth.tests)              # test for effect (RE)
          n.summaries + 6 * n.by + 1:n.by
        )
  }
  else
    del.lines <- c(if (!calcwidth.fixed)   # FE
                     2,
                   if (!calcwidth.random)  # RE
                     3,
                   if (!calcwidth.predict) # PI
                     4,
                   if (!calcwidth.tests)   # tests
                     5:n.summaries)
  ##
  for (i in seq(along = leftcols)) {
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
  ##
  x1 <- unit.c(x1, colgap.forest.left, col.forestwidth)
  ##
  if (rsel) {
    for (i in seq(along = rightcols)) {
      x1 <- unit.c(x1,
                   if (i == 1) colgap.forest.right else colgap.right,
                   wcalc(cols.calc[[rightcols[i]]]$labels))
    }
  }
  
  
  ##
  ##
  ## (13) Process arguments smlab, label.left and label.right
  ##
  ##
  if (by) {
    addline <- addrow * (!any(c(overall.hetstat,
                                test.overall.fixed, test.overall.random,
                                resid.hetstat,
                                test.subgroup.fixed, test.subgroup.random)))
    ##
    nrow <- max(addline + c(yTE, yTE.fixed, yTE.random, yPredict,
                            yStats, yTE.w), na.rm = TRUE)
  }
  else {
    addline <- addrow * (!any(c(test.overall.fixed, test.overall.random,
                                overall.hetstat)))
    ##
    nrow <- max(addline + c(yTE, yTE.fixed, yTE.random, yPredict,
                            yStats), na.rm = TRUE)
  }
  ##
  ymin.line <- overall.hetstat + test.overall.fixed + test.overall.random +
    resid.hetstat + test.subgroup.fixed + test.subgroup.random +
    (1 - missing.text.addline1) + (1 - missing.text.addline2)
  ##
  ymin.line <- ymin.line + (overall & ymin.line == 0 &
                            !(!addrow.overall | !addrow))
  if (hetstat %in% c("fixed", "random") &
      (!missing.text.addline1 | !missing.text.addline2)) {
    ymin.line <- ymin.line + 1
  }
  ##
  ymin.fixed  <- spacing * (ymin.line + prediction + comb.random + 0.5)
  ymin.random <- spacing * (ymin.line + prediction + 0.5)
  ymin.ref    <- spacing * (ymin.line + (!(overall | overall.hetstat) & addrow))
  ##
  ymax <- spacing * (nrow - ifelse(is.na(yHeadadd), 1, 2) - 1 * addrow)
  ##
  ## Position on y-axis of left and right labels (at bottom of forest plot)
  ##
  y.bottom.lr <- ymin.line - 2.5 + (!(overall | overall.hetstat) & addrow)
  ##
  ## Position on y-axis of label below x-axis
  ##
  xlab.ypos <- y.bottom.lr - 1 * (print.label & bottom.lr) -
    1 * (print.label & bottom.lr & (newline.lr | newline.ll))
  ##
  ## Summary label at top of forest plot
  ##
  smlab1 <- tgl(smlab1, unit(smlab.pos, "native"), "center", fs.smlab, ff.smlab,
                fontfamily, rows = 1 + (!is.na(yHeadadd) & !newline.smlab))
  ##
  if (newline.smlab)
    smlab2 <- tgl(smlab2, unit(smlab.pos, "native"),
                  "center", fs.smlab, ff.smlab, fontfamily,
                  rows = 2)
  ##
  ## Left and right label on x-axis:
  ##
  if (!bottom.lr & !is.na(ref)) {
    row1.lr <- if (!newline & (newline.ll | newline.lr) & !addrow)
                 1
               else if (!is.na(yHeadadd) & addrow)
                 2
               else if (is.na(yHeadadd))
                 1
               else
                 2
    ##
    ll1 <- tgl(ll1, unit(ref - (xlim[2] - xlim[1]) / 30, "native"),
               "right", fs.lr, ff.lr, fontfamily, col.label.left,
               rows = row1.lr)
    ##
    if (newline.ll)
      ll2 <- tgl(ll2, unit(ref - (xlim[2] - xlim[1]) / 30, "native"),
                 "right", fs.lr, ff.lr, fontfamily, col.label.left,
                 rows = row1.lr + 1)
    ##
    lr1 <- tgl(lr1, unit(ref + (xlim[2] - xlim[1]) / 30, "native"),
               "left", fs.lr, ff.lr, fontfamily, col.label.right,
               rows = row1.lr)
    ##
    if (newline.lr)
      lr2 <- tgl(lr2, unit(ref + (xlim[2] - xlim[1]) / 30, "native"),
                 "left", fs.lr, ff.lr, fontfamily, col.label.right,
                 rows = row1.lr + 1)
  }
  
  
  ##
  ##
  ## (14) Generate forest plot
  ##
  ##
  if (new)
    grid.newpage()
  ##
  pushViewport(viewport(layout = grid.layout(
                          nrow,
                          length(x1),
                          widths = x1,
                          heights = unit(spacing, "lines"))))
  ##
  ## Left side of forest plot
  ##
  j <- 1
  ##
  for (i in seq(along = leftcols)) {
    add.text(cols[[leftcols[i]]], j)
    ##
    if (!is.na(yHeadadd)) {
      if (!is.null(lab.e.attach.to.col)) {
        if (leftcols[i] == paste0("col.", lab.e.attach.to.col))
          add.text(col.lab.e, j)
      }
      else if (metabin) {
        if (leftcols[i] == "col.n.e" & just.c == "right")
          add.text(col.lab.e, j)
        else if (leftcols[i] == "col.event.e" & just.c %in% c("left", "center"))
          add.text(col.lab.e, j)
      }
      else if (metacont) {
        if (leftcols[i] == "col.sd.e" & just.c == "right")
          add.text(col.lab.e, j)
        else if (leftcols[i] == "col.mean.e" & just.c %in% c("left", "center"))
          add.text(col.lab.e, j)
      }
      else if (metainc) {
        if (leftcols[i] == "col.time.e" & just.c == "right")
          add.text(col.lab.e, j)
        else if (leftcols[i] == "col.event.e" & just.c %in% c("left", "center"))
          add.text(col.lab.e, j)
      }
      ##
      if (!is.null(lab.c.attach.to.col)) {
        if (leftcols[i] == paste0("col.", lab.c.attach.to.col))
          add.text(col.lab.c, j)
      }
      else if (metabin) {
        if (leftcols[i] == "col.n.c" & just.c == "right")
          add.text(col.lab.c, j)
        else if (leftcols[i] == "col.event.c" & just.c %in% c("left", "center"))
          add.text(col.lab.c, j)
      }
      else if (metacont) {
        if (leftcols[i] == "col.sd.c" & just.c == "right")
          add.text(col.lab.c, j)
        else if (leftcols[i] == "col.mean.c" & just.c %in% c("left", "center"))
          add.text(col.lab.c, j)
      }
      else if (metainc) {
        if (leftcols[i] == "col.time.c" & just.c == "right")
          add.text(col.lab.c, j)
        else if (leftcols[i] == "col.event.c" & just.c %in% c("left", "center"))
          add.text(col.lab.c, j)
      }
      ##
      if (newline.studlab & leftcols[i] == "col.studlab")
        add.text(col.add.studlab, j)
      if (newline.effect & leftcols[i] == "col.effect")
        add.text(col.add.effect, j)
      if (newline.ci & leftcols[i] == "col.ci")
        add.text(col.add.ci, j)
      if (newline.effect.ci & leftcols[i] == "col.effect.ci")
        add.text(col.add.effect.ci, j)
      if (newline.w.fixed & leftcols[i] == "col.w.fixed")
        add.text(col.add.w.fixed, j)
      if (newline.w.random & leftcols[i] == "col.w.random")
        add.text(col.add.w.random, j)
      if (newline.TE & leftcols[i] == "col.TE")
        add.text(col.add.TE, j)
      if (newline.seTE & leftcols[i] == "col.seTE")
        add.text(col.add.seTE, j)
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
      ##
      ## Add text in first line of forest plot for new columns
      ##
      if (newcols)
        if (length(leftcols.new) > 0 &
            leftcols[i] %in% paste0("col.", leftcols.new)) {
          sel <- paste0("col.", leftcols.new) == leftcols[i]
          ##
          ## Check for "\n" in label of new column
          ##
          clines <- twolines(leftlabs.new[sel], leftcols[i])
          ##
          just.new <- just.addcols.left[sel]
          ##
          if (just.new == "left")
            xpos.new <- 0
          else if (just.new == "center")
            xpos.new <- 0.5
          else if (just.new == "right")
            xpos.new <- 1
          ##
          ## Add first line
          ##
          if (clines$newline)
            add.text(tgl(clines$top, xpos.new, just.new, fs.head, ff.head,
                         fontfamily), j)
        }
    }
    ##
    j <- j + 2
  }
  ##
  ## Produce forest plot
  ##
  draw.lines(col.forest, j,
             ref, TE.fixed, TE.random,
             overall, comb.fixed, comb.random, prediction,
             ymin.fixed, ymin.random, ymin.ref, ymax,
             lwd, lty.fixed, lty.random, col.fixed, col.random,
             xlim[1], xlim[2],
             lower.equi, upper.equi, lty.equi, col.equi, fill.equi)
  ##
  draw.axis(col.forest, j, yS, log.xaxis, at, label,
            fs.axis, ff.axis, fontfamily, lwd,
            xlim, notmiss.xlim)
  ##
  if (bottom.lr) {
    add.text(smlab1, j, xscale = col.forest$range)
    ##
    if (newline.smlab)
      add.text(smlab2, j, xscale = col.forest$range)
  }
  ##
  if (print.label) {
    if (!bottom.lr) {
      add.text(ll1, j, xscale = col.forest$range)
      ##
      if (newline.ll)
        add.text(ll2, j, xscale = col.forest$range)
      ##
      add.text(lr1, j, xscale = col.forest$range)
      ##
      if (newline.lr)
        add.text(lr2, j, xscale = col.forest$range)
    }
    else {
      add.label(ll1, j,
                unit(ref - (xlim[2] - xlim[1]) / 30, "native"),
                unit(y.bottom.lr, "lines"),
                "right",
                fs.lr, ff.lr, col.label.left, fontfamily,
                xscale = col.forest$range)
      ##
      if (newline.ll)
        add.label(ll2, j,
                  unit(ref - (xlim[2] - xlim[1]) / 30, "native"),
                  unit(y.bottom.lr - 1, "lines"),
                  "right",
                  fs.lr, ff.lr, col.label.left, fontfamily,
                  xscale = col.forest$range)
      ##
      add.label(lr1, j,
                unit(ref + (xlim[2] - xlim[1]) / 30, "native"),
                unit(y.bottom.lr, "lines"),
                "left",
                fs.lr, ff.lr, col.label.right, fontfamily,
                xscale = col.forest$range)
      ##
      if (newline.lr)
        add.label(lr2, j,
                  unit(ref + (xlim[2] - xlim[1]) / 30, "native"),
                  unit(y.bottom.lr - 1, "lines"),
                  "left",
                  fs.lr, ff.lr, col.label.right, fontfamily,
                  xscale = col.forest$range)
    }
  }
  ##
  add.xlab(col.forest, j, xlab, xlab.add, newline.xlab,
           xlab.pos, xlab.ypos, fs.xlab, ff.xlab,
           fontfamily)
  ##
  draw.forest(col.forest, j)
  ##
  j <- j + 2
  ##
  ##
  ## Right side of forest plot
  ##
  ##
  if (rsel) {
    for (i in seq(along = rightcols)) {
      add.text(cols[[rightcols[i]]], j)
      ##
      if (!is.na(yHeadadd)) {
        if (!is.null(lab.e.attach.to.col)) {
          if (rightcols[i] == paste0("col.", lab.e.attach.to.col))
            add.text(col.lab.e, j)
        }
        else if (metabin) {
          if (rightcols[i] == "col.n.e" & just.c == "right")
            add.text(col.lab.e, j)
          else if (rightcols[i] == "col.event.e" &
                   just.c %in% c("left", "center"))
            add.text(col.lab.e, j)
        }
        else if (metacont) {
          if (rightcols[i] == "col.sd.e" & just.c == "right")
            add.text(col.lab.e, j)
          else if (rightcols[i] == "col.mean.e" &
                   just.c %in% c("left", "center"))
            add.text(col.lab.e, j)
        }
        else if (metainc) {
          if (rightcols[i] == "col.time.e" & just.c == "right")
            add.text(col.lab.e, j)
          else if (rightcols[i] == "col.event.e" &
                   just.c %in% c("left", "center"))
            add.text(col.lab.e, j)
        }
        ##
        if (!is.null(lab.c.attach.to.col)) {
          if (rightcols[i] == paste0("col.", lab.c.attach.to.col))
            add.text(col.lab.c, j)
        }
        else if (metabin) {
          if (rightcols[i] == "col.n.c" & just.c == "right")
            add.text(col.lab.c, j)
          else if (rightcols[i] == "col.event.c" &
                   just.c %in% c("left", "center"))
            add.text(col.lab.c, j)
        }
        else if (metacont) {
          if (rightcols[i] == "col.sd.c" & just.c == "right")
            add.text(col.lab.c, j)
          else if (rightcols[i] == "col.mean.c" &
                   just.c %in% c("left", "center"))
            add.text(col.lab.c, j)
        }
        else if (metainc) {
          if (rightcols[i] == "col.time.c" & just.c == "right")
            add.text(col.lab.c, j)
          else if (rightcols[i] == "col.event.c" &
                   just.c %in% c("left", "center"))
            add.text(col.lab.c, j)
        }
        ##
        if (newline.studlab & rightcols[i] == "col.studlab")
          add.text(col.add.studlab, j)
        if (newline.effect & rightcols[i] == "col.effect")
          add.text(col.add.effect, j)
        if (newline.ci & rightcols[i] == "col.ci")
          add.text(col.add.ci, j)
        if (newline.effect.ci & rightcols[i] == "col.effect.ci")
          add.text(col.add.effect.ci, j)
        if (newline.w.fixed & rightcols[i] == "col.w.fixed")
          add.text(col.add.w.fixed, j)
        if (newline.w.random & rightcols[i] == "col.w.random")
          add.text(col.add.w.random, j)
        if (newline.TE & rightcols[i] == "col.TE")
          add.text(col.add.TE, j)
        if (newline.seTE & rightcols[i] == "col.seTE")
          add.text(col.add.seTE, j)
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
        ##
        ## Add text in first line of forest plot for new columns
        ##
        if (newcols)
          if (length(rightcols.new) > 0 &
              rightcols[i] %in% paste0("col.", rightcols.new)) {
            sel <- paste0("col.", rightcols.new) == rightcols[i]
            ##
            ## Check for "\n" in label of new column
            ##
            clines <- twolines(rightlabs.new[sel], rightcols[i])
            ##
            just.new <- just.addcols.right[sel]
            ##
            if (just.new == "left")
              xpos.new <- 0
            else if (just.new == "center")
              xpos.new <- 0.5
            else if (just.new == "right")
              xpos.new <- 1
            ##
            ## Add first line
            ##
            if (clines$newline)
              add.text(tgl(clines$top, xpos.new, just.new,
                           fs.head, ff.head, fontfamily), j)
          }
      }
      ##
      j <- j + 2
    }
  }
  ##
  popViewport()
  
  
  invisible(NULL)
}
