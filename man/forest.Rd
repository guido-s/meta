\name{forest}

\alias{forest}
\alias{forest.meta}
\alias{forest.metabind}

\title{Forest plot to display the result of a meta-analysis}

\description{
  Draws a forest plot in the active graphics window (using grid
  graphics system).
}

\usage{
forest(x, ...)

\method{forest}{meta}(x, sortvar, studlab=TRUE,
       layout=gs("layout"),
       comb.fixed=x$comb.fixed, comb.random=x$comb.random,
       overall=TRUE,
       text.fixed=NULL,
       text.random=NULL,
       lty.fixed=2, lty.random=3, col.fixed="black", col.random="black",
       prediction=x$prediction,
       text.predict=NULL,
       subgroup=TRUE,
       print.subgroup.labels=TRUE,
       bylab=x$bylab, print.byvar=x$print.byvar,
       byseparator=x$byseparator,
       text.fixed.w=text.fixed, text.random.w=text.random, bysort=FALSE, 
       pooled.totals=comb.fixed|comb.random, pooled.events=FALSE,
       pooled.times=FALSE, study.results=TRUE,
       xlab="", xlab.pos,
       smlab=NULL, smlab.pos, xlim="symmetric",
       allstudies=TRUE,
       weight.study, weight.subgroup,
       pscale=x$pscale, irscale=x$irscale, irunit=x$irunit,
       ref=ifelse(backtransf & is.relative.effect(x$sm), 1, 0),
       leftcols=NULL, rightcols=NULL,
       leftlabs=NULL, rightlabs=NULL,
       lab.e=x$label.e, lab.c=x$label.c,
       lab.e.attach.to.col=NULL, lab.c.attach.to.col=NULL,
       label.right=x$label.right, label.left=x$label.left, bottom.lr=TRUE,
       lab.NA=".",
       lab.NA.effect="",
       lwd=1,
       at=NULL, label=TRUE,
       type.study="square", type.fixed="diamond", type.random=type.fixed,
       type.subgroup=ifelse(study.results, "diamond", "square"),
       col.study="black",
       col.square="gray", col.square.lines=col.square,
       col.inside="white",
       col.diamond="gray",
       col.diamond.fixed=col.diamond, col.diamond.random=col.diamond,
       col.diamond.lines="black",
       col.diamond.lines.fixed=col.diamond.lines,
       col.diamond.lines.random=col.diamond.lines,
       col.inside.fixed=col.inside,
       col.inside.random=col.inside,
       col.predict="red", col.predict.lines="black",
       col.by="darkgray",
       col.label.right="black", col.label.left="black",
       hetstat = print.I2 | print.tau2 | print.Q | print.pval.Q | print.Rb,
       overall.hetstat = overall & hetstat,
       hetlab = "Heterogeneity: ",
       print.I2 = comb.fixed | comb.random,
       print.I2.ci = FALSE,
       print.tau2 = comb.fixed | comb.random,
       print.Q = FALSE,
       print.pval.Q = comb.fixed | comb.random,
       print.Rb = FALSE,
       print.Rb.ci = FALSE,
       text.subgroup.nohet = "not applicable",
       ##
       test.overall=gs("test.overall"),
       test.overall.fixed=comb.fixed&overall&test.overall,
       test.overall.random=comb.random&overall&test.overall,
       label.test.overall.fixed, label.test.overall.random,
       print.zval=TRUE,
       ##
       test.subgroup,
       test.subgroup.fixed, test.subgroup.random,
       print.Q.subgroup=TRUE,
       label.test.subgroup.fixed, label.test.subgroup.random,
       ##
       test.effect.subgroup,
       test.effect.subgroup.fixed, test.effect.subgroup.random,
       label.test.effect.subgroup.fixed,
       label.test.effect.subgroup.random,
       ##
       fontsize=12,
       fs.heading = fontsize,
       fs.fixed, fs.random, fs.predict,
       fs.fixed.labels, fs.random.labels, fs.predict.labels,
       fs.study = fontsize, fs.study.labels = fs.study,
       fs.hetstat, fs.test.overall,
       fs.test.subgroup, fs.test.effect.subgroup,
       fs.axis = fontsize, fs.smlab = fontsize, fs.xlab = fontsize,
       fs.lr = fontsize,
       ff.heading = "bold",
       ff.fixed, ff.random, ff.predict,
       ff.fixed.labels, ff.random.labels, ff.predict.labels,
       ff.study = "plain", ff.study.labels = ff.study,
       ff.hetstat, ff.test.overall,
       ff.test.subgroup, ff.test.effect.subgroup,
       ff.axis = "plain", ff.smlab = "bold", ff.xlab = "plain",
       ff.lr = "plain",
       #
       squaresize=0.8 / spacing,
       #
       plotwidth = if (layout != "JAMA") "6cm" else "8cm",
       colgap = "2mm",
       colgap.left = colgap, colgap.right = colgap,
       colgap.studlab = colgap.left, colgap.forest = colgap,
       colgap.forest.left = colgap.forest,
       colgap.forest.right = colgap.forest,
       #
       calcwidth.pooled=TRUE,
       calcwidth.fixed=calcwidth.pooled,
       calcwidth.random=calcwidth.pooled,
       calcwidth.predict=FALSE,
       calcwidth.hetstat=FALSE, calcwidth.tests=FALSE,
       #
       just=if (layout != "JAMA") "right" else "left",
       just.studlab="left", just.addcols="center",
       just.addcols.left=just.addcols, just.addcols.right=just.addcols,
       #
       spacing = 1,
       addrow, addrow.overall, addrow.subgroups,
       #
       new=TRUE,
       #
       backtransf=x$backtransf,
       digits=gs("digits.forest"), digits.se=gs("digits.se"),
       digits.zval=gs("digits.zval"),
       digits.pval=max(gs("digits.pval")-2, 2),
       digits.pval.Q=max(gs("digits.pval.Q")-2, 2),
       digits.Q=gs("digits.Q"),
       digits.tau2=gs("digits.tau2"),
       digits.I2=max(gs("digits.I2")-1, 0),
       digits.weight=gs("digits.weight"),
       #
       digits.mean = digits, digits.sd = digits.se,
       digits.cor = digits, digits.time = digits,
       #
       scientific.pval = gs("scientific.pval"), big.mark = gs("big.mark"),
       #
       col.i=col.study, weight=weight.study,
       ...)

\method{forest}{metabind}(x,
       leftcols, leftlabs,
       rightcols=c("effect", "ci"), rightlabs,
       #
       overall=FALSE, subgroup=FALSE, overall.hetstat=FALSE,
       #
       lab.NA="",
       #
       digits=gs("digits.forest"),
       digits.se=gs("digits.se"),
       digits.zval=gs("digits.zval"),
       digits.pval=max(gs("digits.pval") - 2, 2),
       digits.pval.Q=max(gs("digits.pval.Q") - 2, 2),
       digits.Q=gs("digits.Q"),
       digits.tau2=gs("digits.tau2"),
       digits.I2=max(gs("digits.I2") - 1, 0),
       #
       scientific.pval=gs("scientific.pval"),
       big.mark=gs("big.mark"),
       #
       smlab,
       ...)			    
}

\arguments{
  \item{x}{An object of class \code{meta} or \code{metabind}.}
  \item{sortvar}{An optional vector used to sort the individual studies
    (must be of same length as \code{x$TE}).}
  \item{studlab}{A logical indicating whether study labels should be
    printed in the graph. A vector with study labels can also be
    provided (must be of same length as \code{x$TE} then).}
  \item{layout}{A character string specifying the layout of the forest
    plot (see Details).}
  \item{comb.fixed}{A logical indicating whether fixed effect estimate
  should be plotted.}
  \item{comb.random}{A logical indicating whether random effects estimate
  should be plotted.}
  \item{overall}{A logical indicating whether overall summaries should
    be plotted. This argument is useful in a meta-analysis with
    subgroups if summaries should only be plotted on group level.}
  \item{text.fixed}{A character string used in the plot to label the
    pooled fixed effect estimate.}
  \item{text.random}{A character string used in the plot to label the
    pooled random effects estimate.}
  \item{lty.fixed}{Line type (pooled fixed effect estimate).}
  \item{lty.random}{Line type (pooled random effects estimate).}
  \item{col.fixed}{Line colour (pooled fixed effect estimate).}
  \item{col.random}{Line colour (pooled random effects estimate).}
  \item{prediction}{A logical indicating whether a prediction interval
    should be printed.}
  \item{text.predict}{A character string used in the plot to label the
    prediction interval.}
  \item{subgroup}{A logical indicating whether subgroup results should
    be shown in forest plot. This argument is useful in a
    meta-analysis with subgroups if summaries should not be plotted on
    group level.}
  \item{print.subgroup.labels}{A logical indicating whether subgroup
    label should be printed.}
  \item{bylab}{A character string with a label for the grouping variable.}
  \item{print.byvar}{A logical indicating whether the name of the grouping
    variable should be printed in front of the group labels.}
  \item{byseparator}{A character string defining the separator between
    label and levels of grouping variable.}
  \item{text.fixed.w}{A character string to label the pooled fixed
    effect estimate within subgroups, or a character vector of same
    length as number of subgroups with corresponging labels.}
  \item{text.random.w}{A character string to label the pooled random
    effect estimate within subgroups, or a character vector of same
    length as number of subgroups with corresponging labels.}
  \item{bysort}{A logical indicating whether groups should be ordered
    alphabetically.}
  \item{pooled.totals}{A logical indicating whether total number of
    observations should be given in the figure.}
  \item{pooled.events}{A logical indicating whether total number of
    events should be given in the figure.}  
  \item{pooled.times}{A logical indicating whether total person time
    at risk should be given in the figure.}
  \item{study.results}{A logical indicating whether results for
    individual studies should be shown in the figure (useful to only
    plot subgroup results).}
  \item{xlab}{A label for the x-axis.}
  \item{xlab.pos}{A numeric specifying the center of the label
    on the x-axis.}
  \item{smlab}{A label for the summary measurex (printed at top of figure).}
  \item{smlab.pos}{A numeric specifying the center of the label
    for the summary measure.}
  \item{xlim}{The x limits (min,max) of the plot, or the character "s"
    to produce symmetric forest plots.}
  \item{allstudies}{A logical indicating whether studies with inestimable
    treatment effects should be plotted.}    
  \item{weight.study}{A character string indicating weighting used to
    determine size of squares or diamonds (argument \code{type.study})
    to plot individual study results. One of missing, \code{"same"},
    \code{"fixed"}, or \code{"random"}, can be abbreviated. Plot
    symbols have the same size for all studies or represent study
    weights from fixed effect or random effects model.}
  \item{weight.subgroup}{A character string indicating weighting used
    to determine size of squares or diamonds (argument
    \code{type.subgroup}) to plot subgroup results. One of missing,
    \code{"same"}, or \code{"weight"}, can be abbreviated. Plot
    symbols have the same size for all subgroup results or represent
    subgroup weights from fixed effect or random effects model.}
  \item{pscale}{A numeric giving scaling factor for printing of single
    event probabilities or risk differences, i.e. if argument
    \code{sm} is equal to \code{"PLOGIT"}, \code{"PLN"},
    \code{"PRAW"}, \code{"PAS"}, \code{"PFT"}, or \code{"RD"}.}
  \item{irscale}{A numeric defining a scaling factor for printing of
    single incidence rates or incidence rate differences, i.e. if
    argument \code{sm} is equal to \code{"IR"}, \code{"IRLN"},
    \code{"IRS"}, \code{"IRFT"}, or \code{"IRD"}.}
  \item{irunit}{A character specifying the time unit used to calculate
    rates, e.g. person-years.}
  \item{ref}{A numerical giving the reference value to be plotted as a
    line in the forest plot. No reference line is plotted if argument
    \code{ref} is equal to \code{NA}.}
  \item{leftcols}{A character vector specifying (additional) columns
  to be plotted on the left side of the forest plot or a logical value
  (see Details).}
  \item{rightcols}{A character vector specifying (additional) columns
  to be plotted on the right side of the forest plot or a logical
  value (see Details).}
  \item{leftlabs}{A character vector specifying labels for (additional)
  columns on left side of the forest plot (see Details).}
  \item{rightlabs}{A character vector specifying labels for (additional)
  columns on right side of the forest plot (see Details).}
  \item{lab.e}{Label to be used for experimental group in table
    heading.}
  \item{lab.c}{Label to be used for control group in table
    heading.}
  \item{lab.e.attach.to.col}{A character specifying the column name where
    label \code{lab.e} should be attached to in table heading.}
  \item{lab.c.attach.to.col}{A character specifying the column name where
    label \code{lab.c} should be attached to in table heading.}
  \item{label.left}{Graph label on left side of forest plot.}
  \item{label.right}{Graph label on right side of forest plot.}
  \item{bottom.lr}{A logical indicating whether labels on right
    and left side should be printed at bottom or top of forest plot.}
  \item{lab.NA}{A character string to label missing values.}
  \item{lab.NA.effect}{A character string to label missing values in
    individual treatment estimates and confidence intervals.}
  \item{lwd}{The line width, see \code{\link{par}}.}
  \item{at}{The points at which tick-marks are to be drawn, see
  \code{grid.xaxis}.}
  \item{label}{A logical value indicating whether to draw the
    labels on the tick marks, or an expression or character vector which
    specify the labels to use. See \code{\link{grid.xaxis}}.}
  \item{type.study}{A character string or vector specifying how to
    plot treatment effects and confidence intervals for individual
    studies (see Details).}
  \item{type.fixed}{A character string specifying how to plot
    treatment effect and confidence interval for fixed effect
    meta-analysis (see Details).}
  \item{type.random}{A character string specifying how to plot
    treatment effect and confidence interval for random effects
    meta-analysis (see Details).}
  \item{type.subgroup}{A character string specifying how to plot
    treatment effect and confidence interval for subgroup results (see
    Details).}
  \item{col.study}{The colour for individual study results and
    confidence limits.}
  \item{col.inside}{The colour for individual study results and
    confidence limits if confidence limits are completely within
    squares.}
  \item{col.square}{The colour for squares reflecting study's weight in
    the meta-analysis.}
  \item{col.square.lines}{The colour for the outer lines of squares
    reflecting study's weight in the meta-analysis.}
  \item{col.diamond}{The colour of diamonds representing the results for
    fixed effect and random effects models.}
  \item{col.diamond.fixed}{The colour of diamonds for fixed effect estimates.}
  \item{col.diamond.random}{The colour of diamonds for random effects estimates.}
  \item{col.diamond.lines}{The colour of the outer lines of diamonds
    representing the results for fixed effect and random effects models.}
  \item{col.diamond.lines.fixed}{The colour of the outer lines of
    diamond for fixed effect estimate.}
  \item{col.diamond.lines.random}{The colour of the outer lines of
    diamond for random effects estimate.}
  \item{col.inside.fixed}{The colour for result of fixed effect
    meta-analysis if confidence limit lies completely within square.}
  \item{col.inside.random}{The colour for result of random effects
    meta-analysis if confidence limit lies completely within square.}  
  \item{col.predict}{Background colour of prediction interval.}
  \item{col.predict.lines}{Colour of outer lines of prediction interval.}
  \item{col.by}{The colour to print information on subgroups.}
  \item{col.label.right}{The colour for label on right side of null
    effect.}
  \item{col.label.left}{The colour for label on left side of null
    effect.}
  \item{hetstat}{A logical value indicating whether to print results for
    heterogeneity measures at all.}
  \item{overall.hetstat}{A logical value indicating whether to print
    results for heterogeneity measures for overall treatment
    comparisons. This argument is useful in a meta-analysis with
    subgroups if heterogeneity statistics should only be printed on
    subgroup level.}
  \item{hetlab}{Label printed in front of results for heterogeneity
    measures.}
  \item{print.I2}{A logical value indicating whether to print the value
    of the I-squared statistic.}
  \item{print.I2.ci}{A logical value indicating whether to print the
    confidence interval of the I-squared statistic.}
  \item{print.tau2}{A logical value indicating whether to print the value
    of the between-study variance tau-squared.}
  \item{print.Q}{A logical value indicating whether to print the value
    of the heterogeneity statistic Q.}
  \item{print.pval.Q}{A logical value indicating whether to print the
    p-value of the heterogeneity statistic Q.}
  \item{print.Rb}{A logical value indicating whether to print the value
    of the I-squared statistic.}
  \item{print.Rb.ci}{A logical value indicating whether to print the
    confidence interval of the I-squared statistic.}
  \item{text.subgroup.nohet}{A logical value or character string which
    is printed to indicate subgroups with less than two studies
    contributing to meta-analysis (and thus without heterogeneity). If
    FALSE, heterogeneity statistics are printed (with NAs).}
  \item{test.overall}{A logical value indicating whether to print
    results of test for overall effect.}
  \item{test.overall.fixed}{A logical value indicating whether to
    print results of test for overall effect (based on fixed effect
    model).}
  \item{test.overall.random}{A logical value indicating whether to
    print results of test for overall effect (based on random effects
    model).}
  \item{label.test.overall.fixed}{Label printed in front of results of
    test for overall effect (based on fixed effect model).}
  \item{label.test.overall.random}{Label printed in front of results
    of test for overall effect (based on random effects model).}
  \item{print.zval}{A logical value indicating whether z-value for
    test of treatment effect should be printed.}
  \item{test.subgroup}{A logical value indicating whether to print
    results of test for subgroup differences.}
  \item{test.subgroup.fixed}{A logical value indicating whether to
    print results of test for subgroup differences (based on fixed
    effect model).}
  \item{test.subgroup.random}{A logical value indicating whether to
    print results of test for subgroup differences (based on random
    effects model).}
  \item{print.Q.subgroup}{A logical value indicating whether to print
    the value of the heterogeneity statistic Q (test for subgroup
    differences).}
  \item{label.test.subgroup.fixed}{Label printed in front of results
    of test for subgroup differences (based on fixed effect model).}
  \item{label.test.subgroup.random}{Label printed in front of results
    of test for subgroup differences (based on random effects model).}
  \item{test.effect.subgroup}{A logical value indicating whether to
    print results of test for effect in subgroups.}
  \item{test.effect.subgroup.fixed}{A logical value indicating whether
    to print results of test for effect in subgroups (based on fixed
    effect model).}
  \item{test.effect.subgroup.random}{A logical value indicating
    whether to print results of test for effect in subgroups (based on
    random effects model).}
  \item{label.test.effect.subgroup.fixed}{Label printed in front of
    results of test for effect in subgroups (based on fixed effect
    model).}
  \item{label.test.effect.subgroup.random}{Label printed in front of
    results of test for effect in subgroups (based on random effects
    model).}
  \item{fontsize}{The size of text (in points), see \code{\link{gpar}}.}
  \item{fs.heading}{The size of text for column headings, see \code{\link{gpar}}.}
  \item{fs.fixed}{The size of text for results of fixed effect model, see \code{\link{gpar}}.}
  \item{fs.random}{The size of text for results of random effects model, see \code{\link{gpar}}.}
  \item{fs.predict}{The size of text for results of prediction interval, see \code{\link{gpar}}.}
  \item{fs.fixed.labels}{The size of text for label of fixed effect model, see \code{\link{gpar}}.}
  \item{fs.random.labels}{The size of text for label of random effects model, see \code{\link{gpar}}.}
  \item{fs.predict.labels}{The size of text for label of prediction interval, see \code{\link{gpar}}.}
  \item{fs.study}{The size of text for results of individual studies, see \code{\link{gpar}}.}
  \item{fs.study.labels}{The size of text for labels of individual studies, see \code{\link{gpar}}.}
  \item{fs.hetstat}{The size of text for heterogeneity measures, see \code{\link{gpar}}.}
  \item{fs.test.overall}{The size of text of test for overall effect, see \code{\link{gpar}}.}
  \item{fs.test.subgroup}{The size of text of test of subgroup differences, see \code{\link{gpar}}.}
  \item{fs.test.effect.subgroup}{The size of text of test of effect in subgroups, see \code{\link{gpar}}.}
  \item{fs.axis}{The size of text on x-axis, see \code{\link{gpar}}.}
  \item{fs.smlab}{The size of text of label for summary measure, see \code{\link{gpar}}.}
  \item{fs.xlab}{The size of text of label on x-axis, see \code{\link{gpar}}.}
  \item{fs.lr}{The size of text of label on left and right side of forest plot, see \code{\link{gpar}}.}
  \item{ff.heading}{The fontface for column headings, see \code{\link{gpar}}.}
  \item{ff.fixed}{The fontface of text for results of fixed effect model, see \code{\link{gpar}}.}
  \item{ff.random}{The fontface of text for results of random effects model, see \code{\link{gpar}}.}
  \item{ff.predict}{The fontface of text for results of prediction interval, see \code{\link{gpar}}.}
  \item{ff.fixed.labels}{The fontface of text for label of fixed effect model, see \code{\link{gpar}}.}
  \item{ff.random.labels}{The fontface of text for label of random effects model, see \code{\link{gpar}}.}
  \item{ff.predict.labels}{The fontface of text for label of prediction interval, see \code{\link{gpar}}.}
  \item{ff.study}{The fontface of text for results of individual studies, see \code{\link{gpar}}.}
  \item{ff.study.labels}{The fontface of text for labels of individual studies, see \code{\link{gpar}}.}
  \item{ff.hetstat}{The fontface of text for heterogeneity measures, see \code{\link{gpar}}.}
  \item{ff.test.overall}{The fontface of text of test for overall effect, see \code{\link{gpar}}.}
  \item{ff.test.subgroup}{The fontface of text for test of subgroup differences, see \code{\link{gpar}}.}
  \item{ff.test.effect.subgroup}{The fontface of text for test of effect in subgroups, see \code{\link{gpar}}.}
  \item{ff.axis}{The fontface of text on x-axis, see \code{\link{gpar}}.}
  \item{ff.smlab}{The fontface of text of label for summary measure, see \code{\link{gpar}}.}
  \item{ff.xlab}{The fontface of text of label on x-axis, see \code{\link{gpar}}.}
  \item{ff.lr}{The fontface of text of label on left and right side of forest plot, see \code{\link{gpar}}.}
  \item{squaresize}{A numeric used to increase or decrease the size of
    squares in the forest plot.}
  \item{plotwidth}{Either a character string, e.g., "8cm", "60mm", or
    "3inch", or a \code{\link[grid]{unit}} object specifying width of
    the forest plot.}
  \item{colgap}{Either a character string or a
    \code{\link[grid]{unit}} object specifying gap between columns
    printed on left and right side of forest plot.}
  \item{colgap.left}{Either a character string or a
    \code{\link[grid]{unit}} object specifying gap between columns
    printed on left side of forest plot.}
  \item{colgap.right}{Either a character string or a
    \code{\link[grid]{unit}} object specifying gap between columns
    printed on right side of forest plot.}
  \item{colgap.studlab}{Either a character string or a
    \code{\link[grid]{unit}} object specifying gap between column with
    study labels and subsequent column.}
  \item{colgap.forest}{Either a character string or a
    \code{\link[grid]{unit}} object specifying gap between column
    adjacent to forest plot and the forest plot.}
  \item{colgap.forest.left}{Either a character string or a
    \code{\link[grid]{unit}} object specifying gap between column on
    the left side of forest plot and the forest plot.}
  \item{colgap.forest.right}{Either a character string or a
    \code{\link[grid]{unit}} object specifying gap between column on
    the right side of forest plot and the forest plot.}
  \item{calcwidth.pooled}{A logical indicating whether text for fixed
    effect and random effects model should be considered to calculate
    width of column with study labels, see next two arguments. }
  \item{calcwidth.fixed}{A logical indicating whether text given in
    arguments \code{text.fixed} and \code{text.fixed.w} should be
    considered to calculate width of column with study labels.}
  \item{calcwidth.random}{A logical indicating whether text given in
    arguments \code{text.random} and \code{text.random.w} should be
    considered to calculate width of column with study labels.}
  \item{calcwidth.predict}{A logical indicating whether text given in
    argument \code{text.predict} should be considered to calculate
    width of column with study labels.}
  \item{calcwidth.hetstat}{A logical indicating whether text for
    heterogeneity statistics should be considered to calculate width
    of column with study labels.}
  \item{calcwidth.tests}{A logical indicating whether text for tests
    of overall effect or subgroup differences should be considered to
    calculate width of column with study labels.}
  \item{just}{Justification of text in all columns but columns
    with study labels and additional variables (possible values:
    "left", "right", "center").}
  \item{just.studlab}{Justification of text for study labels (possible
    values: "left", "right", "center").}
  \item{just.addcols}{Justification of text for additional columns (possible
    values: "left", "right", "center").}
  \item{just.addcols.left}{Justification of text for additional
    columns on left side of forest plot (possible values: "left",
    "right", "center"). Can be of same length as number of additional
    columns on left side of forest plot.}
  \item{just.addcols.right}{Justification of text for additional
    columns on right side of forest plot (possible values: "left",
    "right", "center"). Can be of same length as number of additional
    columns on right side of forest plot.}
  \item{spacing}{A numeric determining line spacing in a forest plot.}
  \item{addrow}{A logical value indicating whether an empty row is
    printed above and below study results.}
  \item{addrow.overall}{A logical value indicating whether an empty
    row is printed above overall meta-analysis results.}
  \item{addrow.subgroups}{A logical value indicating whether an empty
    row is printed between results for subgroups.}    
  \item{new}{A logical value indicating whether a new figure should
    be printed in an existing graphics window.}
  \item{backtransf}{A logical indicating whether results should be
    back transformed in forest plots. If \code{backtransf=TRUE},
    results for \code{sm="OR"} are presented as odds ratios rather
    than log odds ratios and results for \code{sm="ZCOR"} are
    presented as correlations rather than Fisher's z transformed
    correlations, for example.}
  \item{digits}{Minimal number of significant digits for treatment
    effects, see \code{print.default}.}
  \item{digits.se}{Minimal number of significant digits for standard
    errors, see \code{print.default}.}
  \item{digits.zval}{Minimal number of significant digits for z- or
    t-statistic for test of overall effect, see \code{print.default}.}
  \item{digits.tau2}{Minimal number of significant digits for
    between-study variance, see \code{print.default}.}
  \item{digits.pval}{Minimal number of significant digits for p-value
    of overall treatment effect, see \code{print.default}.}
  \item{digits.pval.Q}{Minimal number of significant digits for
    p-value of heterogeneity test, see \code{print.default}.}
  \item{digits.Q}{Minimal number of significant digits for
    heterogeneity statistic Q, see \code{print.default}.}
  \item{digits.I2}{Minimal number of significant digits for I-squared
    statistic, see \code{print.default}.}
  \item{digits.weight}{Minimal number of significant digits for
    weights, see \code{print.default}.}
       digits.cor=NULL, digits.time=NULL,
  \item{digits.mean}{Minimal number of significant digits for means;
    only applies to \code{\link{metacont}} objects.}
  \item{digits.sd}{Minimal number of significant digits for standard
    deviations; only applies to \code{\link{metacont}} objects.}
  \item{digits.cor}{Minimal number of significant digits for
    correlations; only applies to \code{\link{metacor}} objects.}
  \item{digits.time}{Minimal number of significant digits for times;
    only applies to \code{\link{metainc}} and \code{\link{metarate}}
    objects.}
  \item{scientific.pval}{A logical specifying whether p-values should
    be printed in scientific notation, e.g., 1.2345e-01 instead of
    0.12345.}
  \item{big.mark}{A character used as thousands separator.} 
  \item{col.i}{Deprecated argument (replaced by \code{col.study}).}
  \item{weight}{Deprecated argument (replaced by \code{weight.study}).}
  \item{\dots}{Additional graphical arguments (ignored at the moment).}
}

\details{
  A forest plot, also called confidence interval plot, is drawn in the
  active graphics window. The \code{forest} function is based on the
  grid graphics system. In order to print the forest plot, (i) resize
  the graphics window, (ii) either use \code{\link{dev.copy2eps}} or
  \code{\link{dev.copy2pdf}}.
  
  By default, treatment estimates and confidence intervals are plotted
  in the following way:
  %%
  \itemize{
   \item For an individual study, a square with treatment estimate in
     the center and confidence interval as line extending either side
     of the square (\code{type.study="square"})
   \item For meta-analysis results, a diamond with treatment estimate
     in the center and right and left side corresponding to lower and
     upper confidence limits (\code{type.fixed="diamond"},
     \code{type.random="diamond"}, and
     \code{type.subgroup="diamond"})
  }
  
  In a forest plot, size of the squares typically reflects the
  precision of individual treatment estimates based either on the
  fixed effect (\code{weight.study="fixed"}) or random effects
  meta-analysis (\code{weight.study="random"}). Information from
  meta-analysis object \code{x} is utilised if argument
  \code{weight.study} is missing. Weights from the fixed effect model
  are used if argument \code{x$comb.fixed} is \code{TRUE}; weights
  from the random effects model are used if argument
  \code{x$comb.random} is \code{TRUE} and \code{x$comb.fixed} is
  \code{FALSE}. The same square sizes are used if
  \code{weight.study="same"}.

  Arguments \code{text.fixed}, \code{text.random}, and
  \code{text.predict} can be used to change the label to identify
  overall results (fixed effect and random effects model as well as
  prediction interval). By default the following text is printed:
  %%
  \itemize{
   \item "Fixed effect model" (argument \code{text.fixed})
   \item "Random effects model" (\code{text.random})
   \item "Prediction interval" (\code{text.predict})
  }
  %%
  If confidence interval levels are different for individual
  studies, meta-analysis, and prediciton interval (arguments
  \code{level}, \code{level.comb}, \code{level.predict} in
  meta-analysis functions, e.g. \code{\link{metabin}}), additional
  information is printed, e.g. " (99\%-CI)" for a 99\% confidence
  interval in the meta-analysis.
  
  The following arguments can be used to print results for various statistical tests:
  \tabular{ll}{
   \bold{Argument} \tab \bold{Statistical test} \cr
   \code{test.overall.fixed} \tab Test for overall effect (fixed effect model) \cr
   \code{test.overall.random} \tab Test for overall effect (random effects model) \cr
   \code{test.effect.subgroup.fixed} \tab Test for effect in subgroup (FE model) \cr
   \code{test.effect.subgroup.random} \tab Test for effect in subgroup (RE model) \cr
   \code{test.subgroup.fixed} \tab Test for subgroup differences (FE model) \cr
   \code{test.subgroup.random} \tab Test for subgroup differences (RE model)
  }
  By default, these arguments are \code{FALSE}. R function
  \code{\link{settings.meta}} can be used to change this default for
  the entire R session. For example, use the following command to
  always print results of tests for an overall effect:
  
  \code{settings.meta(test.overall=TRUE)}
  
  The arguments \code{leftcols} and \code{rightcols} can be used to
  specify columns which are plotted on the left and right side of the
  forest plot, respectively. If argument \code{rightcols} is
  \code{FALSE}, no columns will be plotted on the right side. By
  default, i.e. if arguments \code{leftcols} and \code{rightcols} are
  \code{NULL} and \code{layout="meta"}, the following
  \emph{\bold{columns}} will be printed \emph{\bold{on the right side
  of the forest plot}}:
  \tabular{ll}{
   \bold{Meta-analysis results} \tab \bold{Value of argument rightcols} \cr
   No summary \tab \code{c("effect", "ci")} \cr
   Only fixed effect model \tab \code{c("effect", "ci", "w.fixed")} \cr
   Only random effects model \tab \code{c("effect", "ci", "w.random")} \cr
   Both models \tab \code{c("effect", "ci", "w.fixed", "w.random")}
  }
  By default, estimated treatment effect and corresponding confidence
  interval will be printed. Depending on arguments \code{comb.fixed}
  and \code{comb.random}, weights of the fixed effect and/or random
  effects model will be given too. For an object of class
  \code{\link{metacum}} or \code{\link{metainf}} only the estimated
  treatment effect with confidence interval are plotted.
  
  Depending on the class of the meta-analysis object (which is defined
  by the R function used to generate the object) a different set of
  \emph{\bold{columns}} is printed \emph{\bold{on the left side of the
  forest plot}}:
  \tabular{cl}{
   \bold{Function} \tab \bold{Value of argument leftcols} \cr
   \code{\link{metabin}} \tab \code{c("studlab", "event.e", "n.e",} \cr
   \tab \code{ "event.c", "n.c")} \cr
   \code{\link{metacont}} \tab \code{c("studlab", "n.e", "mean.e", "sd.e",} \cr
   \tab  \code{ "n.c", "mean.c", "sd.c")} \cr
   \code{\link{metacor}} \tab \code{c("studlab", "n")} \cr
   \code{\link{metagen}} \tab \code{c("studlab", "TE", "seTE")} \cr
   \code{\link{metainc}} \tab \code{c("studlab", "event.e", "time.e",} \cr
   \tab  \code{ "event.c", "time.c")} \cr
   \code{\link{metaprop}} \tab \code{c("studlab", "event", "n")} \cr
   \code{\link{metarate}} \tab \code{c("studlab", "event", "time")} \cr
   \code{\link{metacum}} \tab \code{"studlab"} \cr
   \code{\link{metainf}} \tab \code{"studlab"}
  }
  The arguments \code{leftlabs} and \code{rightlabs} can be used to
  specify column headings which are plotted on left and right side of
  the forest plot, respectively. For certain columns predefined labels
  exist. If the arguments \code{leftlabs} and \code{rightlabs} are
  \code{NULL}, the following default labels will be used:
  \tabular{rcccccc}{
  \bold{Column:} \tab \code{studlab} \tab \code{TE} \tab \code{seTE} \tab \code{n.e} \tab \code{n.c} \tab \code{n} \cr
  \bold{Label:} \tab "Study" \tab "TE" \tab "seTE"  \tab "Total" \tab "Total" \tab "Total" \cr\cr
  \bold{Column:} \tab \code{event.e} \tab \code{event.c} \tab \code{event} \tab \code{mean.e} \tab \code{mean.c} \tab \cr
  \bold{Label:} \tab "Events" \tab "Events"  \tab "Events" \tab "Mean" \tab "Mean" \tab \cr\cr
  \bold{Column:} \tab \code{sd.e} \tab \code{sd.c} \tab \code{time.e} \tab \code{time.c} \tab \code{effect} \tab \cr
  \bold{Label:} \tab "SD" \tab "SD" \tab "Time" \tab "Time" \tab \code{x$sm} \tab \cr\cr
  \bold{Column:} \tab \code{ci} \tab \code{effect.ci} \tab \code{w.fixed} \tab \code{w.random} \tab \tab \cr
  \bold{Label:} \tab \code{x$level}"\%-CI" \tab \emph{effect+ci} \tab "W(fixed)" \tab "W(random)" \tab \tab
  }
  For additional columns, the column name will be used as label. It is
  possible to only provide labels for new columns (see
  Examples). Otherwise the length of \code{leftlabs} and
  \code{rightlabs} must be the same as the number of printed columns,
  respectively. The value \code{NA} can be used to specify columns
  which should use default labels (see Examples).

  If argument \code{layout="RevMan5"} (and arguments \code{leftcols}
  and \code{rightcols} are \code{NULL}), the layout for forest plots
  used for Cochrane reviews (which are generated with Review Manager
  5,
  \url{http://community.cochrane.org/tools/review-production-tools/revman-5})
  is reproduced:
  %%
  \enumerate{
   \item All columns are printed on the left side of the forest plot
     (see arguments \code{leftcols} and \code{rightcols})
   \item Tests for overall effect and subgroup differences are printed
     (\code{test.overall}, \code{test.effect.subgroup},
     \code{test.subgroup})
   \item Diamonds representing meta-analysis results are printed in
     black (\code{diamond.fixed}, \code{diamond.random})
   \item Color of squares depends on the meta-analysis object
     (\code{col.square}, \code{col.square.lines})
   \item Information on effect measure and meta-analysis method is
     printed above the forest plot (\code{smlab})
    \item Label "Study or Subgroup" is printed for meta-analysis with
     subgroups (\code{leftlabs})
  }
  If argument \code{layout="JAMA"} (and arguments \code{leftcols} and
  \code{rightcols} are \code{NULL}), instructions for authors of the
  \emph{Journal of the American Medical Association}, see
  \url{http://jamanetwork.com/journals/jama/pages/instructions-for-authors},
  are taken into account:
  %%
  \enumerate{
   \item Graph labels on right and left side are printed in bold font
     at top of forest plot (see arguments \code{bottom.lr} and
     \code{ff.lr})
   \item Information on effect measure and level of confidence
     interval is printed at bottom of forest plot (\code{xlab})
   \item Tests for overall effect are printed (\code{test.overall})
   \item Diamonds representing meta-analysis results are printed in
     lightblue (\code{diamond.fixed}, \code{diamond.random})
   \item Squares representing individual study results are printed in
     darkblue (\code{col.square}, \code{col.square.lines})
   \item Between-study variance \eqn{\tau^2} is not printed
   \item Empty rows are omitted (\code{addrow})
   \item Label "Source" is printed instead of "Study"
     (\code{leftlabs})
  }
  The following changes are conducted if argument
  \code{layout="subgroup"} (and arguments \code{leftcols} and
  \code{rightcols} are \code{NULL}) and a subgroup analysis was
  conducted:
  %%
  \enumerate{
   \item Individual study results are omitted (see argument
     \code{study.results})
   \item Total number of observations is not printed
     (\code{pooled.totals})
   \item Label "Subgroup" is printed instead of "Study"
     (\code{leftlabs})
  }
  
  If arguments \code{lab.e} and \code{lab.c} are \code{NULL},
  "Experimental" and "Control" are used as labels for experimental and
  control group, respectively.
  
  Argument \code{pscale} can be used to rescale single proportions or
  risk differences, e.g. \code{pscale=1000} means that proportions are
  expressed as events per 1000 observations. This is useful in
  situations with (very) low event probabilities.
  
  Argument \code{irscale} can be used to rescale single rates or rate
  differences, e.g. \code{irscale=1000} means that rates are expressed
  as events per 1000 time units, e.g. person-years. This is useful in
  situations with (very) low rates. Argument \code{irunit} can be used
  to specify the time unit used in individual studies (default:
  "person-years"). This information is printed in summaries and forest
  plots if argument \code{irscale} is not equal to 1.
      
  A prediction interval for treatment effect of a new study (Higgins
  et al., 2009) is given in the forest plot if arguments
  \code{prediction} and \code{comb.random} are \code{TRUE}. For
  graphical presentation of prediction intervals the approach by
  Guddat et al. (2012) is used.
  
  Note, in R package \bold{meta}, version 3.0-0 the following
  arguments have been removed from R function forest.meta: byvar,
  level, level.comb, level.predict. This functionality is now provided
  by R function \code{\link{update.meta}} (or directly in R functions,
  e.g., \code{\link{metabin}}, \code{\link{metacont}},
  \code{\link{metagen}}, \code{\link{metacor}}, and
  \code{\link{metaprop}}).
}

\references{
  Guddat C, Grouven U, Bender R, Skipka G (2012), A note on the
  graphical presentation of prediction intervals in random-effects
  meta-analyses. \emph{Systematic Reviews}, \bold{1}, 34.

  Higgins JPT, Thompson SG, Spiegelhalter DJ (2009),
  A re-evaluation of random-effects meta-analysis.
  \emph{Journal of the Royal Statistical Society: Series A},
  \bold{172}, 137-159.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metagen}}, \code{\link{metabind}}, \code{\link{settings.meta}}}

\examples{
data(Olkin95)
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Olkin95, subset=c(41,47,51,59),
                 sm="RR", method="I",
                 studlab=paste(author, year))


\dontrun{
#
# Do standard (symmetric) forest plot
#
forest(meta1)
}


#
# Layout of forest plot similar to Review Mananager 5
# (see http://community.cochrane.org/tools/review-production-tools/revman-5)
#
# Furthermore, add labels on both sides of forest plot and prediction
# interval
#
forest(meta1, layout="RevMan5", comb.fixed=FALSE,
       label.right="Favours control", col.label.right="red",
       label.left="Favours experimental", col.label.left="green",
       prediction=TRUE)


\dontrun{
#
# Sort studies by decreasing treatment effect within year subgroups
#
meta2 <- update(meta1, byvar=ifelse(year < 1987,
                                    "Before 1987", "1987 and later"),
                print.byvar=FALSE)
forest(meta2,
       sortvar=-TE, comb.random=FALSE)


#
# Forest plot specifying argument xlim
#
forest(meta1, xlim=c(0.01, 10))


#
# Print results of test for overall effect
#
forest(meta1, test.overall.fixed=TRUE, test.overall.random=TRUE)


#
# Forest plot with 'classic' layout used in
# R package meta, version < 1.6-0
#
forest(meta1, col.square="black", hetstat=FALSE)


#
# Change set of columns printed on left side
# of forest plot
#
forest(meta1, comb.random=FALSE,
       leftcols="studlab")


#
# Do not print columns on right side of forest plot
#
forest(meta1, rightcols=FALSE)


#
# Change study label to "Author"
#
forest(meta1, comb.random=FALSE,
       leftlabs=c("Author", NA, NA, NA, NA))


#
# Just give effect estimate and 95% confidence interval
# on right side of forest plot (in one column)
#
forest(meta1, rightcols=c("effect.ci"))


#
# Just give effect estimate and 95% confidence interval
# on right side of forest plot
#
forest(meta1, rightcols=c("effect", "ci"))


#
# 1. Change order of columns on left side
# 2. Attach labels to columns 'event.e' and 'event.c'
#    instead of columns 'n.e' and 'n.c'
#
forest(meta1,
       leftcols=c("studlab", "n.e", "event.e", "n.c", "event.c"),
       lab.e.attach.to.col="event.e",
       lab.c.attach.to.col="event.c")


#
# Specify column labels only for newly created variables
# 'year' and 'author' (which are part of dataset Olkin95)
#
forest(meta1,
       leftcols=c("studlab", "event.e", "n.e", "event.c", "n.c",
                  "author", "year"),
       leftlabs=c("Author", "Year of Publ"))


#
# Center text in all columns
#
forest(meta1,
       leftcols=c("studlab", "event.e", "n.e", "event.c", "n.c",
                  "author", "year"),
       leftlabs=c("Author", "Year of Publ"), hetstat=FALSE,
       just="center", just.addcols="center", just.studlab="center")

#
# Same result
#
forest(meta1,
       leftcols=c("studlab", "event.e", "n.e", "event.c", "n.c",
                  "author", "year"),
       leftlabs=c("Author", "Year of Publ"), hetstat=FALSE,
       just="c", just.addcols="c", just.studlab="c")


#
# Change some fontsizes and fontfaces
#
forest(meta1,
       fs.study=10, ff.study="italic",
       fs.study.label=11, ff.study.label="bold",
       fs.axis=5, ff.axis="italic",
       ff.smlab="bold.italic",
       ff.fixed="plain", ff.hetstat="plain")


#
# Change some colours
#
forest(meta1,
       col.diamond="green", col.diamond.lines="red",
       col.study=c("green", "blue", "red", "orange"),
       col.square="pink", col.square.lines="black")


#
# Sort by weight in fixed effect model
#
forest(meta1, sortvar=1/w.fixed, comb.random=FALSE)


#
# Sort by decreasing weight in fixed effect model
#
forest(meta1, sortvar=-1/w.fixed, comb.random=FALSE)


#
# Sort by size of treatment effect
#
forest(meta1, sortvar=TE, comb.random=FALSE)


#
# Sort by size of treatment effect
#
forest(meta1, sortvar=-TE, comb.random=FALSE)


#
# Sort by decreasing year of publication
#
forest(meta1, sortvar=-year, comb.random=FALSE)


#
# Print results of test for subgroup differences (random effects model)
#
forest(meta2,
       sortvar=-TE, comb.fixed=FALSE,
       test.subgroup.random=TRUE)

#
# Print only subgroup results
#
forest(meta2, layout="subgroup")

#
# Print only subgroup results
# (and consider text for heterogeneity measures in width of subgroup
#  column)
#
forest(meta2, layout="subgroup", calcwidth.hetstat=TRUE)
}
}

\keyword{hplot}
