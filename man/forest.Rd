\name{forest}

\alias{forest}
\alias{forest.meta}

\title{Forest plot to display the result of a meta-analysis}

\description{
  Draws a forest plot in the active graphics window (using grid
  graphics system).
}

\usage{
forest(x, ...)

\method{forest}{meta}(x, sortvar, studlab=TRUE,
       comb.fixed=x$comb.fixed, comb.random=x$comb.random,
       overall=TRUE,
       text.fixed=if (x$level!=x$level.comb) paste("Fixed effect model (",
                    round(x$level.comb*100), "\%-CI)", sep="")
                  else "Fixed effect model",
       text.random=if (x$level!=x$level.comb) paste("Random effects model (",
                    round(x$level.comb*100), "\%-CI)", sep="")
                  else "Random effects model",
       lty.fixed=2, lty.random=3,
       prediction=x$prediction,
       text.predict=if (!(length(x$level.predict)==0) &&
                        x$level!=x$level.predict)
                     paste("Prediction interval (",
                           round(x$level.predict*100), "\%)",
                     sep="") else "Prediction interval",
       bylab=x$bylab, print.byvar=x$print.byvar,
       text.fixed.w=text.fixed, text.random.w=text.random,
       bysort=FALSE,
       pooled.totals=comb.fixed|comb.random, pooled.events=FALSE,
       xlab="", xlab.pos=ref,
       smlab=NULL, smlab.pos=ref, xlim="symmetric",
       allstudies=TRUE,
       weight,
       pscale=1,
       ref=ifelse(backtransf & is.relative.effect(x$sm), 1, 0),
       layout="meta",
       leftcols=NULL, rightcols=NULL,
       leftlabs=NULL, rightlabs=NULL,
       lab.e=x$label.e, lab.c=x$label.c,
       lab.e.attach.to.col=NULL, lab.c.attach.to.col=NULL,
       label.right=x$label.right, label.left=x$label.left,
       lab.NA=".",
       lab.NA.effect="",
       lwd=1,
       at=NULL, label=TRUE,
       col.i="black", col.i.inside.square="white",
       col.square="gray", col.square.lines=col.square,
       col.diamond="gray",
       col.diamond.fixed=col.diamond, col.diamond.random=col.diamond,
       col.diamond.lines="black",
       col.diamond.fixed.lines=col.diamond.lines,
       col.diamond.random.lines=col.diamond.lines,
       col.predict="red", col.predict.lines="black",
       col.by="darkgray",
       print.I2=comb.fixed|comb.random, print.I2.ci=FALSE,
       print.tau2=comb.fixed|comb.random,
       print.Q=FALSE, print.pval.Q=comb.fixed|comb.random,
       hetstat=print.I2|print.tau2|print.Q|print.pval.Q,
       overall.hetstat=overall&hetstat,
       hetlab="Heterogeneity: ",
       text.I2 = "I-squared",
       text.tau2 = "tau-squared",
       test.overall=.settings$test.overall,
       test.overall.fixed=comb.fixed&overall&test.overall,
       test.overall.random=comb.random&overall&test.overall,
       label.test.overall.fixed=paste("Test for overall effect",
         if (comb.fixed & comb.random) " (fixed effect)", ": ", sep=""),
       label.test.overall.random=paste("Test for overall effect",
         if (comb.fixed & comb.random) " (random effects)", ": ", sep=""),
       test.subgroup=.settings$test.subgroup,
       test.subgroup.fixed=if (missing(test.subgroup)) FALSE else test.subgroup,
       test.subgroup.random=if (missing(test.subgroup)) !is.null(x$byvar)&
         comb.random&test.subgroup else test.subgroup,
       print.Q.subgroup=print.Q,
       label.test.subgroup.fixed="Test for subgroup differences (fixed effect): ",
       label.test.subgroup.random=paste("Test for subgroup differences",
        if (test.subgroup.fixed | comb.fixed) " (random effects)", ": ", sep=""),
       fontsize=12,
       fs.heading=fontsize,
       fs.fixed=fontsize, fs.random=fs.fixed, fs.predict=fs.fixed, fs.study=fontsize,
       fs.fixed.labels=fs.fixed, fs.random.labels=fs.random, fs.predict.labels=fs.predict,
       fs.study.labels=fs.study, fs.hetstat=fontsize-2,
       fs.test.overall=fs.hetstat, fs.test.subgroup=fs.hetstat, fs.axis=fontsize,
       fs.smlab=fontsize, fs.xlab=fontsize, fs.lr=fontsize,
       ff.heading="bold",
       ff.fixed="bold", ff.random=ff.fixed, ff.predict=ff.fixed, ff.study="plain",
       ff.fixed.labels=ff.fixed, ff.random.labels=ff.random, ff.predict.labels=ff.predict,
       ff.study.labels=ff.study, ff.hetstat="bold.italic",
       ff.test.overall=ff.hetstat, ff.test.subgroup=ff.hetstat, ff.axis="plain",
       ff.smlab="bold", ff.xlab="plain", ff.lr="plain",
       #
       squaresize=0.8,
       #
       plotwidth=grid::unit(6, "cm"),
       colgap=grid::unit(2, "mm"),
       colgap.left=colgap, colgap.right=colgap,
       colgap.forest=colgap,
       colgap.forest.left=colgap.forest, colgap.forest.right=colgap.forest,
       #
       just="right", just.studlab="left", just.addcols="center",
       #
       addspace=TRUE,
       #
       new=TRUE,
       #
       backtransf=x$backtransf,
       digits=2, digits.se=4, digits.tau2=4, digits.pval=4, digits.pval.Q=digits.pval,
       digits.Q=1, digits.I2=1, ...)
}

\arguments{
  \item{x}{An object of class \code{meta}.}
  \item{sortvar}{An optional vector used to sort the individual studies
    (must be of same length as \code{x$TE}).}
  \item{studlab}{A logical indicating whether study labels should be
    printed in the graph. A vector with study labels can also be
    provided (must be of same length as \code{x$TE} then).}
  \item{level}{The level used to calculate confidence intervals for
    individual studies.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{comb.fixed}{A logical indicating whether fixed effect estimate
  should be plotted.}
  \item{comb.random}{A logical indicating whether random effects estimate
  should be plotted.}
  \item{overall}{A logical indicating whether overall summaries should
    be plotted. This argument is useful in combination with the
    argument \code{byvar} if summaries should only be plotted on group
    level.}
  \item{text.fixed}{A character string used in the plot to label the
    pooled fixed effect estimate.}
  \item{text.random}{A character string used in the plot to label the
    pooled random effects estimate.}
  \item{lty.fixed}{Line type (pooled fixed effect estimate).}
  \item{lty.random}{Line type (pooled random effects estimate).}
  \item{prediction}{A logical indicating whether a prediction interval
    should be printed.}
  \item{level.predict}{The level used to calculate prediction interval
    for a new study.}
  \item{text.predict}{A character string used in the plot to label the
    prediction interval.}
  \item{bylab}{A character string with a label for the grouping variable.}
  \item{print.byvar}{A logical indicating whether the name of the grouping
    variable should be printed in front of the group labels.}
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
  \item{weight}{A character string indicating which type of plotting
    symbols is to be used for individual treatment estimates. One of
    missing (see Details), \code{"same"}, \code{"fixed"}, or
    \code{"random"}, can be abbreviated. Plot symbols have the same size
    for all studies or represent study weights from fixed effect or
    random effects model.}
  \item{pscale}{A numeric giving scaling factor for probabilities for
    objects of class \code{metaprop}.}
  \item{ref}{A numerical giving the reference value to be plotted as a
    line in the forest plot. No reference line is plotted if argument
    \code{ref} is equal to \code{NA}.}
  \item{layout}{A character string specifying the layout of the forest
    plot (see Details).}
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
  \item{lab.NA}{A character string to label missing values.}
  \item{lab.NA.effect}{A character string to label missing values in
    individual treatment estimates and confidence intervals.}
  \item{lwd}{The line width, see \code{\link{par}}.}
  \item{at}{The points at which tick-marks are to be drawn, see
  \code{grid.xaxis}.}
  \item{label}{A logical value indicating whether to draw the
    labels on the tick marks, or an expression or character vector which
    specify the labels to use. See \code{\link{grid.xaxis}}.}
  \item{col.i}{The colour for individual study results and confidence limits.}
  \item{col.i.inside.square}{The colour for individual study results and
    confidence limits if confidence limits are completely within squares.}
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
  \item{col.diamond.fixed.lines}{The colour of the outer lines of diamond
    for fixed effect estimate.}
  \item{col.diamond.random.lines}{The colour of the outer lines of diamond
    for random effects estimate.}
  \item{col.predict}{Background colour of prediction interval.}
  \item{col.predict.lines}{Colour of outer lines of prediction interval.}
  \item{col.by}{The colour to print information on subgroups.}
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
  \item{hetstat}{A logical value indicating whether to print results for
    heterogeneity measures at all.}
  \item{overall.hetstat}{A logical value indicating whether to print
    results for heterogeneity measures for overall treatment
    comparisons. This argument is useful in combination with subgroup
    analyses (list object \code{byvar}) if heterogeneity statistics
    should only be printed on group level.}
  \item{hetlab}{Label printed in front of results for heterogeneity
    measures.}
  \item{text.I2}{Label for I-squared statistic.}
  \item{text.tau2}{Label for between-study variance tau-squared.}
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
  \item{test.subgroup}{A logical value indicating whether to print
    results of test for subgroup differences.}
  \item{test.subgroup.fixed}{A logical value indicating whether to
    print results of test for subgroup differences (based on fixed
    effect model).}
  \item{test.subgroup.random}{A logical value indicating whether to
    print results of test for subgroup differences (based on random
    effects model).}
  \item{label.test.subgroup.fixed}{Label printed in front of results
    of test for subgroup differences (based on fixed effect model).}
  \item{label.test.subgroup.random}{Label printed in front of results
    of test for subgroup differences (based on random effects model).}
  \item{print.Q.subgroup}{A logical value indicating whether to print
    the value of the heterogeneity statistic Q (test for subgroup
    differences).}
  \item{fontsize}{The size of text (in points), see \code{\link{gpar}}.}
  \item{fs.heading}{The size of text for column headings, see \code{\link{gpar}}.}
  \item{fs.fixed}{The size of text for results of fixed effect model, see \code{\link{gpar}}.}
  \item{fs.random}{The size of text for results of random effects model, see \code{\link{gpar}}.}
  \item{fs.predict}{The size of text for results of prediction interval, see \code{\link{gpar}}.}
  \item{fs.study}{The size of text for results of individual studies, see \code{\link{gpar}}.}
  \item{fs.fixed.labels}{The size of text for label of fixed effect model, see \code{\link{gpar}}.}
  \item{fs.random.labels}{The size of text for label of random effects model, see \code{\link{gpar}}.}
  \item{fs.predict.labels}{The size of text for label of prediction interval, see \code{\link{gpar}}.}
  \item{fs.study.labels}{The size of text for labels of individual studies, see \code{\link{gpar}}.}
  \item{fs.hetstat}{The size of text for heterogeneity measures, see \code{\link{gpar}}.}
  \item{fs.test.overall}{The size of text of test for overall effect, see \code{\link{gpar}}.}
  \item{fs.test.subgroup}{The size of text of test of subgroup differences, see \code{\link{gpar}}.}
  \item{fs.axis}{The size of text on x-axis, see \code{\link{gpar}}.}
  \item{fs.smlab}{The size of text of label for summary measure, see \code{\link{gpar}}.}
  \item{fs.xlab}{The size of text of label on x-axis, see \code{\link{gpar}}.}
  \item{fs.lr}{The size of text of label on left and right side of forest plot, see \code{\link{gpar}}.}
  \item{ff.heading}{The fontface for column headings, see \code{\link{gpar}}.}
  \item{ff.fixed}{The fontface of text for results of fixed effect model, see \code{\link{gpar}}.}
  \item{ff.random}{The fontface of text for results of random effects model, see \code{\link{gpar}}.}
  \item{ff.predict}{The fontface of text for results of prediction interval, see \code{\link{gpar}}.}
  \item{ff.study}{The fontface of text for results of individual studies, see \code{\link{gpar}}.}
  \item{ff.fixed.labels}{The fontface of text for label of fixed effect model, see \code{\link{gpar}}.}
  \item{ff.random.labels}{The fontface of text for label of random effects model, see \code{\link{gpar}}.}
  \item{ff.predict.labels}{The fontface of text for label of prediction interval, see \code{\link{gpar}}.}
  \item{ff.study.labels}{The fontface of text for labels of individual studies, see \code{\link{gpar}}.}
  \item{ff.hetstat}{The fontface of text for heterogeneity measures, see \code{\link{gpar}}.}
  \item{ff.test.overall}{The fontface of text of test for overall effect, see \code{\link{gpar}}.}
  \item{ff.test.subgroup}{The fontface of text for test of subgroup differences, see \code{\link{gpar}}.}
  \item{ff.axis}{The fontface of text on x-axis, see \code{\link{gpar}}.}
  \item{ff.smlab}{The fontface of text of label for summary measure, see \code{\link{gpar}}.}
  \item{ff.xlab}{The fontface of text of label on x-axis, see \code{\link{gpar}}.}
  \item{ff.lr}{The fontface of text of label on left and right side of forest plot, see \code{\link{gpar}}.}
  \item{squaresize}{A numeric used to increase or decrease the size of
    squares in the forest plot.}
  \item{plotwidth}{A unit object specifying width of the forest plot.}
  \item{colgap}{A unit object specifying gap between columns printed on
    left and right side of forest plot.}
  \item{colgap.left}{A unit object specifying gap between columns printed on
    left side of forest plot.}
  \item{colgap.right}{A unit object specifying gap between columns printed on
    right side of forest plot.}
  \item{colgap.forest}{A unit object specifying gap between column
    adjacent to forest plot and the forest plot.}
  \item{colgap.forest.left}{A unit object specifying gap between column
    on the left side of forest plot and the forest plot.}
  \item{colgap.forest.right}{A unit object specifying gap between column
    on the right side of forest plot and the forest plot.}
  \item{just}{Justification of text in all columns but columns
    with study labels and additional variables (possible values:
    "left", "right", "center").}
  \item{just.studlab}{Justification of text for study labels (possible
    values: "left", "right", "center").}
  \item{just.addcols}{Justification of text for additional columns (possible
    values: "left", "right", "center").}
  \item{addspace}{A logical value indicating whether additional space
    (i.e. a blank row) is printed above and below study results.}
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
 \item{\dots}{Additional graphical arguments (ignored at the moment).}
}

\details{
  A forest plot, also called confidence interval plot, is drawn in the
  active graphics window. Subgroup analyses are conducted and
  displayed in the plot if \code{byvar} is not missing.
  
  The \code{forest} function is based on the grid graphics
  system. In order to print the forest plot, (i) resize the graphics
  window, (ii) either use \code{\link{dev.copy2eps}} or
  \code{\link{dev.copy2pdf}}.

  Result of a test for overall effect or subgroup differences are
  printed if argument \code{test.overall.fixed},
  \code{test.overall.random}, \code{test.subgroup.fixed}, or
  \code{test.subgroup.random} is \code{TRUE}; by default, these
  arguments are \code{FALSE}. R function \code{\link{settings.meta}}
  can be used to change this default for the entire R session by
  setting argument \code{test.overall=TRUE} or
  \code{test.subgroup=TRUE}.
  
  Information from object \code{x} is utilised if argument \code{weight}
  is missing. Weights from the fixed effect model are used
  (\code{weight="fixed"}) if argument \code{x$comb.fixed} is
  \code{TRUE}; weights from the random effects model are used
  (\code{weight="random"}) if argument \code{x$comb.random}
  is \code{TRUE} and \code{x$comb.fixed} is \code{FALSE}.

  The arguments \code{leftcols} and \code{rightcols} can be used to
  specify columns which are plotted on the left and right side of the
  forest plot, respectively. If argument \code{rightcols} is
  \code{FALSE}, no columns will be plotted on the right side. By
  default, i.e. if arguments \code{leftcols} and \code{rightcols} are
  \code{NULL} and \code{layout="meta"}, the following default columns
  will be plotted.

  \tabular{cl}{
   Additional settings \tab Default for argument \code{rightcols} \cr
   --- \tab \code{c("effect", "ci")} \cr
   \code{comb.fixed=TRUE} & \code{comb.random=FALSE} \tab \code{c("effect", "ci", "w.fixed")} \cr
   \code{comb.fixed=FALSE} & \code{comb.random=TRUE} \tab \code{c("effect", "ci", "w.random")} \cr
   \code{comb.fixed=TRUE} & \code{comb.random=TRUE} \tab \code{c("effect", "ci", "w.fixed", "w.random")}
  }
  By default, estimated treatment effect and corresponding confidence
  interval will be printed. Depending on arguments \code{comb.fixed}
  and \code{comb.random}, weights of the fixed effect and/or random
  effects model will be given too. For an object of class
  \code{metacum} or \code{metainf} only the estimated treatment effect
  with confidence interval are plotted.
  
  Depending on the class of the meta-analysis object a different set
  of columns is printed on the left side of a forest plot.
  \tabular{cl}{
   Object type \tab Default for argument \code{leftcols} \cr
   metabin \tab \code{c("studlab", "event.e", "n.e", "event.c", "n.c")} \cr
   metacont \tab \code{c("studlab", "n.e", "mean.e",
    "sd.e", "n.c", "mean.c", "sd.c")} \cr
   metacor \tab \code{c("studlab", "n")} \cr
   metagen \tab \code{c("studlab", "TE", "seTE")} \cr
   metainc \tab \code{c("studlab", "event.e",
     "time.e", "event.c", "time.c")} \cr
   metaprop \tab \code{c("studlab", "event", "n")} \cr
   metacum \tab \code{c("studlab")} \cr
   metainf \tab \code{c("studlab")}
  }
  The arguments \code{leftlabs} and \code{rightlabs} can be used to
  specify column headings which are plotted on left and right side of
  the forest plot, respectively. For certain columns predefined labels
  exist. If the arguments \code{leftlabs} and \code{rightlabs} are
  \code{NULL}, the following default labels will be used:
  \tabular{ccccccc}{
  Column \tab \code{studlab} \tab \code{TE} \tab \code{seTE} \tab \code{n.e} \tab \code{n.c} \tab \code{n} \cr
  Label \tab "Study" \tab "TE" \tab "seTE"  \tab "Total" \tab "Total" \tab "Total" \cr\cr
  Column \tab \code{event.e} \tab \code{event.c} \tab \code{event} \tab \code{mean.e} \tab \code{mean.c} \tab \cr
  Label \tab "Events" \tab "Events"  \tab "Events" \tab "Mean" \tab "Mean" \tab \cr\cr
  Column \tab \code{sd.e} \tab \code{sd.c} \tab \code{time.e} \tab \code{time.c} \tab \tab \cr
  Label \tab "SD" \tab "SD" \tab "Time" \tab "Time" \tab \tab \cr\cr
  Column \tab \code{effect} \tab \code{ci} \tab \code{w.fixed} \tab \code{w.random} \tab \tab \cr
  Label \tab \code{x$sm} \tab \code{x$level}-CI \tab "W(fixed)" \tab "W(random)" \tab \tab
  }
  For additional columns the column name will be used
  as label. It is possible to only provide labels for new columns (see
  Examples). Otherwise the length of \code{leftlabs} and
  \code{rightlabs} must be the same as the number of printed
  columns. The value \code{NA} can be used to specify columns using
  the default labels (see Example).

  If argument \code{layout="revman5"} (and arguments \code{leftcols}
  and \code{rightcols} are \code{NULL}) all columns mentioned above
  will be printed on the left side of the forest plot.

  If arguments \code{lab.e} and \code{lab.c} are \code{NULL},
  "Experimental" and "Control" are used as labels for experimental and
  control group, respectively.
  
  The arguments \code{pscale} can be used to rescale proportions
  for objects of class \code{metaprop}, e.g. \code{pscale=100} means
  that proportions are expressed per 100 observations. This is useful in
  situations with (very) low proportions. For \code{pscale=100}, column
  heading and x-axis label are changed to "Prop (in \%)" and
  "Proportion (in \%)", respectively.

  A prediction interval for treatment effect of a new study (Higgins
  et al., 2009) is given in the forest plot if arguments
  \code{prediction} and \code{comb.random} are \code{TRUE}. For
  graphical presentation of prediction intervals the approach by
  Guddat et al. (2012) is used.
  
  Note, in R package meta, version 3.0-0 the following arguments have
  been removed from R function forest.meta: byvar, level, level.comb,
  level.predict. This functionality is now provided by R function
  \code{\link{update.meta}} (or directly in R functions
  \code{\link{metabin}}, \code{\link{metacont}},
  \code{\link{metagen}}, \code{\link{metacor}}, and
  \code{\link{metaprop}}).
}

\references{
  Guddat C, Grouven U, Bender R, Skipka G 2012), A note on the
  graphical presentation of prediction intervals in random-effects
  meta-analyses. \emph{Systematic Reviews}, \bold{1}, 34.

  Higgins JPT, Thompson SG, Spiegelhalter DJ (2009),
  A re-evaluation of random-effects meta-analysis.
  \emph{Journal of the Royal Statistical Society: Series A},
  \bold{172}, 137-159.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metagen}}}

\examples{
data(Olkin95)
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Olkin95, subset=c(41,47,51,59),
                 sm="RR", method="I",
                 studlab=paste(author, year))


#
# Do (symmetric) forest plot
#
forest(meta1)


#
# Forest plot specifying argument xlim
#
forest(meta1, xlim=c(0.01, 10))


#
# Add prediction interval to forest plot
#
forest(meta1, prediction=TRUE)


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
# Layout of forest plot similar to RevMan 5
#
forest(meta1, layout="revman5", comb.random=FALSE)


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
\dontrun{
#
# Same result
#
forest(meta1,
       leftcols=c("studlab", "event.e", "n.e", "event.c", "n.c",
                  "author", "year"),
       leftlabs=c("Author", "Year of Publ"), hetstat=FALSE,
       just="c", just.addcols="c", just.studlab="c")
}

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
       col.i=c("green", "blue", "red", "orange"),
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
# Sort studies by decreasing treatment effect within year subgroups
#
meta2 <- update(meta1, byvar=ifelse(year < 1987, "Before 1987", "1987 and later"),
                print.byvar=FALSE)
forest(meta2,
       sortvar=-TE, comb.random=FALSE)

#
# Print results of test for subgroup differences (random effects model)
#
forest(meta2,
       sortvar=-TE, comb.fixed=FALSE,
       test.subgroup.random=TRUE)
}

\keyword{hplot}
