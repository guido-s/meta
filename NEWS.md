## meta, version 7.0-0 (2024-01-11)

### Major changes

* Vignettes added
  - *vignette("meta-workflow")* (workflow)
  - *vignette("meta-tutorial")* (R commands from [Balduzzi et al.,
  2019](https://scholar.google.com/scholar?q=balduzzi+schwarzer+2019))

* R package **metadat** added to Depends (to access meta-analysis datasets)

* R package **robvis** added to Suggests (for risk of bias assessment)

* New functions rob(), barplot.rob() and traffic_light() for risk of
  bias assessment (RoB)

* New function read.cdir() to import Cochrane data package from
  Cochrane review of interventions

* New function blup.meta() to calculate best linear unbiased predictors (BLUPs)

* New functions estimates(), estimates.meta(), and estimates.blup.meta() to
  extract meta-analysis results

* Lower and upper confidence interval limits of individual study
  results stored as transformed limits for meta-analysis with single
  proportions or rates (in previous versions of **meta**, CI limits were
  back transformed for individual studies and not back transformed for
  meta-analysis results)

* Changes for forest plots:
  - forest plot can be directly saved to a file using common graphics
    device drivers (height of file is determined automatically)
  - BMJ layout implemented (layout = "BMJ")
  - details on meta-analysis methods can be shown in plot
  - risk of bias assessment automatically added for meta-analyses with
    RoB assessment
  - point estimates can be plotted as circles or diamonds instead of
    squares
  - default settings for columns on left or right side of forest plot
    can be defined in settings.meta()
  - common effect and random effects confidence intervals and
    prediction intervals are truncated if lower / upper limit is outside the
    limits of the x-axis

* New general setting "BMJ", i.e., R command *settings.meta("BMJ")*,
  to print results according to BMJ style and formating checklist,
  see, for example, [ BMJ
  Medicine](https://bmjmedicine.bmj.com/bmjmedicine/wp-content/uploads/sites/66/2023/06/BMJMED-style-formatting-checklist-for-original-research-pre-acceptance-1.pdf)

* R function metabind() can return both common effect and random
  effects results as well as prediction intervals

* Internal functions for (back) transformations made visible /
  accessible to the user; see help("meta-transf")

* Within-cluster correlation can be specified for three-level model
  (by default, rho = 0)

* Use approximate formulae for Hedges' g and Cohen's d for RevMan 5
  settings; see help page of settings.meta()

* Argument 'pscale' or 'irscale' could be used in principle with any
  effect measure which is useful for user-specified summary measure in
  metagen()

* New R function plot.meta() which calls forest.meta() internally

### User-visible changes

* metabin(), metacont(), metacor(), metainc(), metamean(), metaprop(),
  metarate(), update.meta():
  - new argument 'rho' to specify within-cluster correlation in
    three-level model

* print.meta():
  - show group sample sizes in printouts
  - print confidence intervals based on t- and normal distribution for
    metacont() or metamean() objects with a single study and argument
    'method.ci = "t"'
  - new argument 'print.Q' to suppress printing of heterogeneity
    statistic Q and test of heterogeneity
  - default for argument 'details.methods' can be defined using
    settings.meta()
  - do not print degrees of freedom for Hartung-Knapp or Kenward-Rogers
    intervals if argument \code{overall = FALSE}

* gs():
  - first argument can be a character vector instead of a character string
    to get the default setting of several arguments (which would be returned
    as a list)
  - new argument 'unname' to return named arguments (if unname = FALSE)

* metamerge():
  - can be used with object created with copas() or limitmeta() from R
    package **metasens** as only input (which adds the respective
    results to the standard meta-analysis object)
  - arguments 'label1' and 'label2' can be used to provide defaults to
    amend labels for common effect or random effects model, prediction
    intervals and subgroups
  - new arguments 'label1.common', 'label2.common', 'label1.random',
    'label2.random', 'label1.predict', 'label2.predict',
    'label1.subgroup', and 'label2.subgroup'

* metaadd():
  - argument 'data' can be a meta-analysis object created with R
    package **meta**
  - new arguments 'method.common', 'method.random', 'method.tau',
    'method.random.ci', and 'method.predict'

* metabind():
  - new arguments 'common', 'random' and 'prediction' replacing
    argument 'pooled'

* forest.meta():
  - new arguments 'file', 'width', 'rows.gr', 'func.gr', 'args.gr',
    and 'dev.off' to directly store a forest plot in a file
  - new arguments 'rob', 'rob.col', 'rob.symbols', 'rob.attach',
    'rob.xpos', 'rob.legend', 'fs.rob', 'fs.rob.symbols', 'ff.rob',
    'ff.rob.symbols', 'colgap.rob' and 'just.rob' for risk of bias
    assessment
  - new arguments 'details', 'fs.details' and 'ff.details' to add
    details on meta-analytical methods
  - point estimates can be plotted as circles instead of squares or
    diamonds (arguments 'type.study', 'type.common', 'type.random',
    'type.subgroup', 'type.subgroup.common', 'type.subgroup.random')
  - new arguments 'col.circle' and 'col.circle.lines' to define the
    colour of circles
  - new argument 'digits.n' and 'digits.event' to specify the number
    of significant digits for sample sizes and number of events
  - justification of results for effect + confidence interval can be
    specified by argument 'just' if argument 'layout = "RevMan5"'

* settings.meta():
  - new argument 'overall.hetstat' to specify whether to show
    information on between-study heterogeneity
  - new argument 'width' to specify width of graphics device
  - new arguments 'print.tau2', 'print.tau2.ci', 'print.tau' and
    'print.tau.ci' to specify whether to show (confidence intervals
    for) tau^2 or tau in printouts
  - new arguments 'leftcols', 'rightcols', 'leftlabs', 'rightlabs',
    'label.e.attach' and 'label.c.attach' to changes defaults for
    corresponding arguments in forest.meta()

* bubble.meta():
  - new arguments 'pscale' and 'irscale' added

* metaprop(), metarate():
  - list elements 'lower' and 'upper' contain the transformed lower
    and upper confidence interval limits for individual studies

* R functions for transformation:
  - transf(), cor2z(), p2asin(), p2logit(), VE2logVR()

* R functions for back transformation:
  - backtransf(), asin2ir(), asin2p(), logit2p(), logVR2VE(), z2cor()

* labels.meta():
  - R function is not exported

### Bug fixes

* funnel.meta():
  - automatically calculated limits on x-axis were to narrow for some
    settings

* metamean():
  - argument 'null.effect' was ignored to calculate the test statistic
    and p-value for individual studies (list elements 'statistic' and
    'pval')

* metabin():
  - use continuity correction if sm = "VE"
  
* metareg():
  - error if input to argument 'formula' was the name of an R function 

* print.summary.meta():
  - print correct backtransformed subgroup results for *metabind* objects
    with metaprop objects with Freeman-Tukey transformation as input

### Internal changes

* New internal function gh() to determine height of graphics file

* New internal function smlab() to determine the label for the summary measure

* Use of vcalc() from R package **metafor** to calculate the
  variance-covariance matrix in three-level model with within-cluster
  correlation not equal to 0

* funnel.meta(): list element of meta-analysis object can be directly
  specified in arguments 'text', 'col' and 'bg', e.g., argument 'text
  = studlab' to use study labels instead of plotting symbols

* Input to chkchar() can be a numeric vector

* settings.meta(): set defaults for arguments 'forest.tau2,' 'forest.tau',
  'forest.I2', 'forest.Q', 'forest.pval.Q', and 'forest.Rb' for BMJ, JAMA
   and RevMan5 layout


## meta, version 6.5-0 (2023-06-06)

### Major changes

* In R function metamerge(), user can decide whether to keep or ignore
  information from second meta-analysis on study weights and
  heterogeneity statistics

* New function metaadd() to add pooled results from external analysis
  to meta-analysis object

* Function update.meta() considers arguments 'method.mean' and 'method.sd'

* Variables with group specific information can be merged into a
  single variable in longarm()

* Additional thresholds can be specified to plot vertical lines in
  forest plots, e.g., to mark large, moderate and small effects

* Baujat plot can be used to evaluate influence of studies on random
  effects estimate

* Seed can be specified in meta-analysis functions to calculate
  reproducible bootstrap prediction intervals

### User-visible changes

* metamerge():
  - new arguments 'common1', 'random1', 'prediction1', 'common2',
    'random2', 'prediction2' to specify whether to keep common effect
    results, random effects results or prediction interval from first
    or second meta-analysis
  - new arguments 'keep', 'keep.Q', 'keep.I2' and 'keep.w' to
    determine whether additional information from second meta-analysis
    should be kept
  - new arguments 'common', 'random', 'prediction', 'overall' and
    'overall.hetstat' to specify which results to print
  - new arguments 'hetlabel1', 'hetlabel2', 'text.common1',
    'text.common2', 'text.random1', 'text.random2', 'text.predict1'
    and 'text.predict2' to label results from first or second
    meta-analysis

* longarm():
  - new arguments 'id1' and 'id2' to specify last character(s) of
    variable names with group specific information

### Bug fixes

* metabias():
  - do not conduct test for funnel plot asymmetry for three-level
    model (the test did not consider the cluster structure)

* forest.meta():
  - header line was concealed by equivalence region
  - error if argument 'resid.hetstat = TRUE' was used for subgroup
    meta-analysis without common between-study variance estimate in
    subgroups (argument 'tau.common = FALSE' in meta-analysis
    functions)

* read.rm5():
  - fix bug for error message *"In gsub("\x80", "EUR", txt) : unable
    to translate '<80>' to a wide string"* due to change in default
    settings in R function gsub()

### User-visible changes

* metabin(), metacont(), metacor(), metacr(), metainc(), metamean(),
  metaprop(), metarate(), update.meta():
  - new arguments 'seed.predict' and 'seed.predict.subgroup'
  - print an error message if bootstrap prediction interval is
    requested for three-level model

* trimfill.default():
  - new argument 'seed.predict'

* baujat.meta():
  - new argument 'pooled'

* forest.meta():
  - arguments 'lower.equi' and 'upper.equi' can be numeric vectors
  - new arguments 'fill.lower.equi' and 'fill.upper.equi' to specify
    fill colour(s) for lower or upper limits

* print.metabias():
  - do not print the intercept and its standard error (nuisance parameter)

### Internal changes

* full rewrite of function catmeth() to print details on meta-analysis
  methods

* internal functions is.wholenumber() etc. renamed to
  *is_wholenumber()* due to new CRAN policy regarding generic
  functions

* metabias():
  - intercept and its standard error are no longer part of list element
    'estimate' (nuisance parameter)


## meta, version 6.2-1 (2023-02-28)

### Bug fixes

* forest.meta():
  - correct order of cluster variable for three-level meta-analysis
    with subgroups or use of argument 'sortvar'

* metabin(), metacont(), metacor(), metainc(), metamean(), metaprop(),
  metarate():
  - recognise argument 'subset' for three-level meta-analysis

* *Ad hoc* variance correction for Hartung-Knapp method was not used
  for argument 'adhoc.hakn = "IQWiG6"' (bug was introduced in
  **meta**, version 6.0-0)

### User-visible changes

* metacr():
  - new argument 'Q.Cochrane'

* metareg():
  - using the regression formula as first unnamed argument will result
    in an error not a warning

### Internal changes

* metagen(): check argument 'sm' for known summary measures in lower
  case

* setchar(): new arguments 'return.NULL' and 'nchar.equal'

* New branch 'release' on GitHub starting with **meta**, version 6.2-1


## meta, version 6.2-0 (2023-02-14)

### Major changes

* New function trimfill.rm5() to conduct trim-and-fill method for all
  or selected meta-analyses of a Cochrane review

### User-visible changes

* nnt.meta(), print.nnt.meta():
  - NNTs for hazard ratios can be calculated following [Altman &
    Andersen (1999)](https://doi.org/10.1136/bmj.319.7223.1492)
  - print sensible confidence interval for NNTs if treatment effect is
    non-significant [(Altman,
    1998)](https://doi.org/10.1136/bmj.317.7168.1309)

### Internal changes

* Use generic functions metacum() for cumulative meta-analysis,
  metainf() for leave-one-out meta-analysis and metareg() for
  meta-regression

* New internal function chksuitable() to check for suitable classes


## meta, version 6.1-0 (2022-12-20)

### Major changes

* Meta-analysis of *Vaccine Efficacy* or *Vaccine Effectiveness*
  implemented

* For the generic inverse variance method,
  - untransformed values can be provided for treatment estimates and
    confidence limits, see argument 'transf'
  - original confidence limits for individual studies kept if
    arguments 'lower' and 'upper' are not missing
  - standard error set to missing if lower and upper confidence limits
    are identical (resulted in an error)

* Methods by [McGrath et.,
  (2020)](https://doi.org/10.1177/0962280219889080) and [Cai et.,
  (2021)](https://doi.org/10.1177/09622802211047348) implemented to
  approximate means and standard deviations from median and related
  statistics

* R package
  [**estmeansd**](https://cran.r-project.org/package=estmeansd)
  added to suggested packages to provided methods by [McGrath et.,
  (2020)](https://doi.org/10.1177/0962280219889080) and [Cai et.,
  (2021)](https://doi.org/10.1177/09622802211047348)

* Print header line in forest plots with JAMA or RevMan5 layout

### Bug fixes

* forest.meta():
  - no error for subgroup meta-analysis conducted with metarate()
    using argument 'n' to specify the sample size
  - same square sizes for a meta-analysis with or without subgroups
  - no error for Revman5 and JAMA layout in meta-analyses with more
    than one random effects method

* metarate():
  - calculate number of observations in subgroups if argument 'n' is
    provided

* metabind(), forest.metabind():
  - use correct study / method labels if no subgroup results are present

* settings.meta():
  - argument 'addrows.below.overall' not 'addrow.below.overall'

* metagen():
  - calculate correct prediction interval limits for three-level
    models
  - number of studies in meta-analysis equal to number of non-missing
    estimates and standard errors (only number of non-missing standard
    errors was considered)

### User-visible changes

* All meta-analysis functions:
  - print prediction interval(s) if argument 'method.predict' is not
    missing

* metabin, metagen(), metainc():
  - argument 'sm = "VE"' can be used for meta-analysis of vaccine
    efficacy or vaccine effectiveness

* metagen(), settings.meta():
  - new argument 'transf'

* forest.meta():
  - new argument 'header.line' to add header line
  - new argument 'digits.TE' to specify number of digits for
    transformed treatment estimates (list element 'TE')
  - use more informative column labels for 'TE' and 'seTE'
  
* metabin(), metainc(), metaprop() and metarate():
  - for GLMMs, stop with error if argument 'adhoc.hakn.ci' or
    'adhoc.hakn.pi' is unequal to ""

* nnt.meta(), nnt.default():
  - new argument 'small.values' to specify whether small treatment
    effects indicate a beneficial or harmful effect

* settings.meta():
  - new argument 'digits.TE.forest' to set default for argument
    'digits.TE' in forest.meta()

* Print blank space before negative upper confidence interval limit if
  separator is equal to "-"

* New help page *meta-sm* summarising available summary measures

* Help page of nnt() updated

* Change maintainer's email address

### Internal changes

* New internal functions transf(), cor2z(), p2asin(), logVR2VE() and
  VE2logVR()

* chknumeric():
  - new argument 'integer' to check for integer values

* List element 'df.Q.b.random' with degrees of freedom for test of
  subgroup differences under random effects model is a list instead of
  a vector if more than one random effects method was used (argument
  'method.random.ci')

* Several changes for meta-analysis using generalised linear mixed or
  three-level models


## meta, version 6.0-0 (2022-09-17)

### Major changes

* Meta-analysis object can contain results of several common effect or
  random effects methods, e.g., random effects meta-analysis with or
  without Hartung-Knapp method

* Kenward-Roger method implemented to estimate confidence or
  prediction interval [(Partlett & Riley,
  2017)](https://doi.org/10.1002/sim.7140)

* Bootstrap approach implemented to calculate prediction interval
  [(Nagashima et al., 2019)](https://doi.org/10.1177/0962280218773520)

* Rewrite of function metamerge() to merge pooled results of two
  meta-analyses into a single meta-analysis object

* Defaults for appearance of forest plots can be defined for the R
  session

* R package
  [**pimeta**](https://cran.r-project.org/package=pimeta)
  added to suggested packages in order to calculate bootstrap approach
  for prediction interval

* New argument 'method.random.ci' replaces argument 'hakn' to select
  method to calculate confidence interval for random effects estimate

* Major update of help pages:
  - help page for *meta-package* revised; content with details in
    meta-analysis functions moved to this help page
  - new help page *meta-object* describing content of meta-analysis
    functions; corresponding content moved from individual help pages

### Bug fixes

* forest.meta():
  - do not print label for subgroups with no information to print,
    e.g., if argument 'study.results = FALSE', for subgroups with only
    one or no study contributing to pooled estimate in the subgroup
  - do not show empty row before label on x-axis (argument 'xlab') if
    argument 'label.left' or 'label.right' is provided for
    meta-analysis without reference value (argument 'ref'), e.g.,
    meta-analysis of single means or proportions

* metareg():
  - use Paule-Mandel estimator if used in meta-analysis (instead of
    REML estimator)

### User-visible changes

* metabin(), metacont(), metacor(), metagen(), metainc(), metamean(),
  metaprop() and metarate():
  - argument 'hakn' replaced by 'method.random.ci'
  - argument 'adhoc.hakn' replaced by 'adhoc.hakn.ci'
  - new arguments 'method.predict' and 'adhoc.hakn.pi'

* settings.meta():
  - several new arguments added to define defaults for forest plots;
    see printout of command settings.meta(print = TRUE)

* forest.meta():
  - argument 'col.by' has been renamed to 'col.subgroup',
  - argument 'bysort' has been renamed to 'sort.subgroup'

* trimfill.meta():
  - arguments 'level', 'level.ma', 'method.random.ci', 'adhoc.hakn',
    'method.tau', 'method.tau.ci', 'level.predict', and
    'method.predict' removed


## meta, version 5.5-0 (2022-07-11)

### Major changes

* Use term 'common effect model' instead of 'fixed effect model' in
  the documentation and argument 'common' instead of 'fixed' to (not)
  show results for common effect model

* Three-level model implemented in all meta-analysis functions

* For continuity corrections, new argument 'method.incr' replaces
  arguments 'allincr' and 'addincr' for meta-analysis with binary
  outcome or incidence rates

* Exact Poisson confidence limits can be calculated for individual
  studies in meta-analysis of single rates

* Show information on statistical significance and between-study
  heterogeneity in forest plots of cumulative or leave-one-out
  meta-analysis

* Calculate Cochran's Q directly in **meta** for classic inverse
  variance meta-analysis (instead of taking it from **metafor**
  package)

* By default, do not print warnings for deprecated arguments; this can
  be changed with command 'settings.meta(warn.deprecated = TRUE)'


### Bug fixes

* Use correct standard error for Cox and Snell's method in smd2or()
  and or2smd()

* Three-level model did not work if variable from dataset was
  provided as input to argument 'id' in metacont()

* Argument 'tau.common = TRUE' was ignored in subgroup analysis of
  three-level model in metacont()

* Argument 'level' was ignored in the calculation of confidence limits
  for individual studies in metacont() and metamean() if argument
  'method.ci = "t"'

* Show correct studies in forest plot with subgroups and missing
  treatment effects if argument 'allstudies = FALSE'

* Show points in bubble plot of meta-regression with GLMM

### User-visible changes

* For three-level models,
  - argument 'id' has been renamed to 'cluster'
  - cluster variable is shown in forest plots

* New arguments 'common' and 'cluster' in functions metabin(),
  metacont(), metacor(), metagen(), metainc(), metamean(), metaprop()
  and metarate()

* New function subset.longarm() to select subset of a longarm object

* New argument 'method.ci' in function metarate()

* New argument 'method.ci.rate' in function settings.meta()

* New argument 'method.incr' in functions metabin(), metainc(),
  metaprop() and metarate()

* print.summary.meta():
  - for a single study and metabin() with method = "MH", sm = "RR" and
    RR.Cochrane = FALSE, print results using a continuity correction
    for sample sizes of 1x incr (individual study) and 2x incr
    (meta-analysis of single study)

### Internal changes

* forest.meta():
  - use meta:::formatN() instead of format() for formatting
  - print study label "1" instead of "" for a single study

* metarate():
  - list elements 'lower' and 'upper' contain untransformed confidence
    limits for individual studies

* New internal function update_needed() to check whether update of
  meta object is needed

* metabin(), metacont(), metacor(), metagen(), metainc(), metamean(),
  metaprop() and metarate():
  - new list element 'k.TE' with number of estimable effects


## meta, version 5.2-0 (2022-02-04)

### Major changes

* Forest plot for meta-analysis with subgroups:
  - more flexible printing of subgroup results
  - by default, do not show subgroup results (pooled estimates and
    information on heterogeneity) for subgroups consisting of a single
    study

* Prediction intervals in subgroups can be shown independently of
  prediction interval for overall meta-analysis in printouts and
  forest plots

* Bubble plot shows relative treatment effects on original scale
  instead of log scale and reference line is shown

* Trim and fill, limit meta-analysis and Copas selection model objects
  can be used in function metabind()

* New function longarm() to transform data from pairwise comparisons
  to long arm-based format

* New auxiliary function labels.meta() to create study labels for
  forest plots in JAMA or Lancet layout

* Printing of spaces in confidence intervals can be suppressed

* Help page of forest.meta() updated

### Bug fixes

* Use correct standard error to calculate prediction interval if
  Hartung-Knapp method was used

* In forest plots, show correct degrees of freedom for test of effect
  in subgroups for Hartung-Knapp method

* In update.meta(), consider input for arguments 'pscale', 'irscale'
  and 'irunit' for meta-analysis objects created with metagen()

### User-visible changes

* forest.meta():
  - new argument 'subgroup.hetstat'
  - arguments 'subgroup', 'subgroup.hetstat', 'prediction.subgroup',
    'test.effect.subgroup', 'test.effect.subgroup.fixed' and
    'test.effect.subgroup.random' can be a logical vector of same
    length as number of subgroups
  - arguments 'lab.e', 'lab.c', 'lab.e.attach.to.col' and
    'lab.c.attach.to.col' renamed to 'label.e', 'label.c',
    'label.e.attach' and 'label.c.attach'

* forest.meta(), metabin(), metacont(), metacor(), metacr(),
  metagen(), metainc(), metamean(), metaprop(), metarate(),
  print.meta(), update.meta():
  - new argument 'prediction.subgroup'

* metamerge():
  - first argument can be of class 'limitmeta' or 'copas'

* bubble.metareg():
  - new argument 'backtransf' to (not) back transform relative
    treatment effects on y-axis
  - new arguments 'ref', 'col.ref', 'lty.ref' and 'lwd.ref' for
    reference line

* settings.meta():
  - arguments 'print', 'reset' and 'setting' can be used like any
    other setting; for example, it is possible to fully reset the
    settings and switch to the RevMan 5 settings
  - R commands 'settings.meta("print")' and 'settings.meta()' produce
    the same printout
  - new global setting 'prediction.subgroup' for prediction intervals
    in subgroups
  - new global settings 'CIlower.blank' and 'CIupper.blank'

* cilayout():
  - new arguments 'lower.blank' and 'upper.blank' to suppress printing
    of spaces in confidence intervals
  - additional checks for length of arguments

### Internal changes

* metagen():
  - new list elements 'seTE.hakn' and 'seTE.hakn.adhoc' (with standard
    error for Hartung-Knapp method) and 'seTE.classic' for classic
    random effects inverse variance method

* forest.meta():
  - new code to assign missing column labels

* Internal function formatCI() considers values for 'lower.blank' and
  'upper.blank' in cilayout()

* New internal function catch() to catch value for an argument


## meta, version 5.1-1 (2021-12-02)

### Major changes

* For meta-analysis of single proportions,
  - export p-value of exact binomial test for individual studies if
    Clopper-Pearson method (method.ci = "CP") is used to calculate
    confidence intervals for individual studies
  - do not export p-value for individual studies if argument
    'method.ci' is not equal to "CP" or "NAsm" (normal approximation
    based on summary measure)

### Bug fixes

* Meta-analysis of continuous outcomes using Hedges' g or Cohen's d as
  summary measure resulted in [inestimable SMDs in individual
  studies](https://github.com/guido-s/meta/issues/42) if the total
  sample size was larger than 343 and argument 'exact.smd' was TRUE
  (default)

* Forest plot creation for meta-analysis of single means with
  subgroups resulted in an
  [error](https://github.com/guido-s/meta/issues/41)

### Internal changes

* New internal function ciClopperPearson() to calculate confidence
  limits and p-value for exact binomial method

* Exported list elements changed for internal functions
  ciAgrestiCoull(), ciSimpleAsymptotic() and ciWilsonScore()


## meta, version 5.1-0 (2021-11-17)

### Major changes

* By default, use exact formulae in estimation of the standardised
  mean difference (Hedges' g, Cohen's d) and its standard error
  [(White & Thomas, 2005)](https://doi.org/10.1191/1740774505cn081oa)

### Bug fixes

* Use of metagen() with argument 'id' (three-level model) does not
  result in an error if all estimates come from a single study

### Internal changes

* Fix errors due to extended checks of arguments equal to NULL in R
  package **metafor**, version 3.1 or above


## meta, version 5.0-1 (2021-10-20)

### Major changes

* For backward compatibility, use Q statistic based on Mantel-Haenszel
  estimate (argument 'Q.Cochrane') by default to calculate
  DerSimonian-Laird estimator of the between-study variance

### Bug fixes

* For small sample sizes, use correct entry from Table 2 in [Wan et.
  (2014)](https://doi.org/10.1186/1471-2288-14-135) to approximate
  standard deviation from median and related statistics


## meta, version 5.0-0 (2021-10-11)

### Major changes

* Behaviour of print.meta() and print.summary.meta() switched (to be
  in line with other print and summary functions in R)

* New default settings:
  - Restricted maximum likelihood (REML) instead of DerSimonian-Laird
    estimator used as default to estimate between-study heterogeneity
	(argument 'method.tau')
  - Do not use Q statistic based on Mantel-Haenszel estimate to
    calculate DerSimonian-Laird estimator of the between-study
    variance (argument 'Q.Cochrane')
  - Print 'Common effect model' instead of 'Fixed effect model'

* Default settings of **meta**, version 4 or lower, can be used with
  command *settings.meta("meta4")* - this does not change the new
  behaviour of print.meta() and print.summary.meta()

* Renamed arguments:
  - 'fixed' instead of 'comb.fixed'
  - 'random' instead of 'comb.random'
  - 'level.ma' instead of 'level.comb'
  - 'subgroup' instead of 'byvar'
  - 'subgroup.name' instead of 'bylab'
  - 'print.subgroup.name' instead of 'print.byvar'
  - 'sep.subgroup' instead of 'byseparator'
  - 'nchar.subgroup' instead of 'bylab.nchar'

### Internal changes

* Function gs() can be used to access internal settings

* Store internal auxiliary functions in files meta-aux.R to
  meta-xlab.R


## meta, version 4.19-2 (2021-09-29)

### Bug fixes

* Forest plots of meta-analyses assuming a common between-study
  heterogeneity variance in subgroups resulted in an error
  (bug was introduced in **meta**, version 4.16-0)

* For GLMMs, export Wald-type Q statistic for residual heterogeneity
  instead of missing value


## meta, version 4.19-1 (2021-09-14)

### Bug fixes

* metagen():
  - set random effects weights equal to zero for estimates with
    standard errors equal to NA (to fix error bubble.metareg)

* metareg():
  - for three-level model, use 'test = "t"' instead of 'test = "knha"'
    in internal call of rma.mv()

### User-visible changes

* summary.meta():
  - print tau2 and tau for subgroups with single study if argument
    'tau.common = TRUE'

* bubble.metareg():
  - show regression lines for a single categorical covariate


## meta, version 4.19-0 (2021-08-05)

### Major changes

* Subgroup analysis for three-level model fully implemented

* New default for forest plots to show results of test for subgroup
  differences in meta-analyses with subgroups

* Calculation of weights for three-level random effects model using
  weights.rma.mv() with argument 'type = "rowsum"' from R package
  **metafor**

* Print study label provided by argument 'studlab' for meta-analysis
  with a single study

* Total number of observations and events printed in summaries (if
  available)

### Bug fixes

* metagen():
  - treatment estimates for three-level models with subgroups were not
    based on common between-study variance despite argument
    'tau.common = TRUE'

* metareg():
  - use rma.mv() from R package **metafor** for three-level models
    instead of rma.uni()

### User-visible changes

* metabin(), metacont(), metacor(), metacr(), metagen(), metagen(),
  metainc(), metamean(), metaprop(), metarate():
  - new argument 'test.subgroup' to print results of test for subgroup
    differences

* print.meta():
  - for three-level models, column with grouping information added to
    study details

* metagen():
  - default for estimation of between-study variance has changed for
    three-level models with subgroups, i.e., tau2 is allowed to be
    different in subgroups by default

### Internal changes

* metagen():
  - new variable '.idx' with running index in meta-analysis dataset
    (list element 'data')
  - new logical list element 'three.level' indicating whether
    three-level model was used


## meta, version 4.18-2 (2021-06-11)

### Bug fixes

* For argument 'adhoc.hakn = "ci"', directly compare width of
  confidence intervals of Hartung-Knapp method and classic random
  effects meta-analysis


## meta, version 4.18-1 (2021-05-11)

### Major changes

* Calculate correct upper limit for confidence intervals of I2 and H2
  in very homogeneous meta-analyses (i.e., if Q < k - 1)

### Bug fixes

* forest.meta():
  - correct order of p-values for homogeneity tests within subgroups
    if argument 'bysort = TRUE'

* calcH():
  - set H = 1 in calculation of confidence interval for H if H < 1
    (i.e., if Q < k - 1)

* metabias():
  - bug fix for linear regression tests using **metafor**, version
    2.5-86

* metabind():
  - bug fix for a single meta-analysis object

### Internal changes

* metabias.bias():
  - argument '...' passed on to rma.uni()

* metagen():
  - set list element 'df.hakn' to NA instead of NULL if condition met
    for argument 'adhoc.hakn = "ci"'


## meta, version 4.18-0 (2021-03-05)

### Major changes

* Prediction intervals for subgroups implemented

### Bug fixes

* metacont():
  - use correct variance formula for Glass' delta

* metainc():
  - update command resulted in an error *Arguments 'event.e' and 'n.e'
    must have the same length* for meta-analysis with subgroups (due
    to list elements 'n.e.w' and 'n.c.w' which were interpreted as
    'n.e' and 'n.c' containing missing values instead of being NULL)

* print.meta():
  - use of argument 'details = TRUE' resulted in an error in
    meta-analyses with duplicated study labels

* Consider argument 'adhoc.hakn' to calculate confidence intervals in
  random effects subgroup meta-analyses

### User-visible changes

* print.meta():
  - column with information on subgroups added to details if argument
    'details = TRUE'

* forest.meta():
  - new argument 'text.predict.w' to label the prediction interval in
    subgroups
  - arguments 'text.fixed.w' and 'text.random.w' checked for correct
    length

* *Ad hoc* variance correction for Hartung-Knapp method not available
  for GLMMs

### Internal changes

* metacont():
  - get rid of warnings 'Unknown or uninitialised column' if argument
    'subset' is used

* subgroup():
  - calculate prediction intervals for subgroups


## meta, version 4.17-0 (2021-02-11)

### Major changes

* Tests of funnel plot asymmetry:
  - tests by [Macaskill et
    al. (2001)](https://doi.org/10.1002/sim.698) and [Pustejovsky &
    Rodgers (2019)](https://doi.org/10.1002/jrsm.1332) added
  - use regtest() from R package **metafor** internally for linear
    regression tests
  - new print layout providing more details

* New dataset Pagliaro1992 for meta-analysis on prevention of first
  bleeding in cirrhosis [(Pagliaro et
  al., 1992)](https://doi.org/10.7326/0003-4819-117-1-59)

### Bug fixes

* update.meta():
  - do not switch to three-level model if method.tau = "ML"

### User-visible changes

* metabias():
  - use name of first author to select test for funnel plot asymmetry
    instead of "rank", "linreg", "mm", "count", and "score" (can be
    abbreviated; old names are still recognised)

* print.metabias():
  - new arguments 'digits.stat', 'digits.se', 'digits.pval',
    'scientific.pval', 'big.mark', 'zero.pval', 'JAMA.pval'

### Internal changes

* linregcore():
  - complete rewrite using rma.uni() and regtest() from R package
    **metafor**


## meta, version 4.16-2 (2021-01-27)

### Bug fixes

* drapery():
  - use correct limits on y-axis for argument 'type = "zvalue"'

### User-visible changes

* funnel.meta():
  - inverse of square root of sample size can be plotted on y-axis
    (argument 'yaxis = "invsqrtsize"')

* forest.meta():
  - consider input for argument 'hetstat' to print heterogeneity
    statistics for overall results (see argument 'overall.hetstat')

* metabin(), metacont(), metacor(), metagen(), metagen(), metainc(),
  metamean(), metaprop(), metarate():
  - studies with missing values for subgroup variable (argument
    'byvar') can be excluded from meta-analysis using argument
    'subset'

### Internal changes

* funnel.meta():
  - try to derive sample sizes from list elements 'n.e' or 'n.c' if
    argument 'yaxis = "size"'


## meta, version 4.16-1 (2021-01-19)

### Bug fixes

* For argument 'adhoc.hakn = "ci"', use correct query to determine
  whether confidence interval of Hartung-Knapp method is smaller than
  classic random effects meta-analysis ([Hybrid method 2 in Jackson et
  al., 2017](https://doi.org/10.1002/sim.7411))


## meta, version 4.16-0 (2021-01-18)

### Major changes

* Three-level meta-analysis models can be fitted for generic and
  continuous outcomes ([Van den Noortgate et.,
  2013](https://doi.org/10.3758/s13428-012-0261-6)) by calling
  rma.mv() from R package **metafor** internally

* Measures I2 and H for residual heterogeneity are based on Q
  statistic for residual heterogeneity (instead of taken directly from
  **metafor** package)

* Additional *ad hoc* method implemented if confidence interval of
  Hartung-Knapp method is smaller than classic random effects
  meta-analysis ([Hybrid method 2 in Jackson et al.,
  2017](https://doi.org/10.1002/sim.7411))

* For funnel plot of a diagnostic test accuracy meta-analysis, use
  *effective sample size* ([Deeks et.,
  2005](https://doi.org/10.1016/j.jclinepi.2005.01.016)) by default on
  the y-axis

* New function metamerge() to merge pooled results of two
  meta-analyses into a single meta-analysis object

### Bug fixes

* metabin():
  - Mantel-Haenszel method of risk differences did not use continuity
    correction in case of studies with a zero cell count (argument
    'MH.exact = FALSE')

* metabin(), metainc(), metaprop(), metarate():
  - for GLMMs, confidence limits for classic random effects
    meta-analysis were calculated instead of confidence limits for
    Hartung-Knapp if argument 'hakn = TRUE'

* metabin(), metainc(), metaprop(), metarate():
  - works for GLMMs with zero events or number of events equal to
    number of patients in all studies

* forest.meta():
  - print results for test of subgroup effect in correct order if
    argument 'bysort = TRUE'

* read.rm5():
  - list elements 'method' and 'sm' had been encoded as a factor
    instead of character under R-versions below 4.0 which resulted in
    an error using metacr()
   

### User-visible changes

* Do not print empty confidence intervals for heterogeneity statistics

* metacont(), metagen(), update.meta():
  - new argument 'id' to specify which estimates belong to the same
    study (or laboratory) in order to use three-level model

* metabind():
  - argument '...' can be a single list of meta-analysis objects
  - meta-analyses can use different methods, e.g., different
    estimators of the between-study variance
  
* All meta-analysis functions:
  - argument 'adhoc.hakn = "iqwig6"' instead of 'adhoc.hakn = "ci"'
    uses the *ad hoc* method for Hartung-Knapp method described in
    General Methods 6.0 (IQWiG, 2020)
  - argument 'adhoc.hakn = "ci"' uses the *ad hoc* method described in
    Jackson et al. (2017)

* forest.meta():
  - column heading "Mean" instead of "MLN" for meta-analysis object
    created with metamean() with arguments 'sm = "MLN"' and
    'backtransf = TRUE'
  - study labels specified by argument 'studlab' tried to catch from
    meta-analysis object
  - do not print statistic for residual heterogeneity if argument
    'tau.common = FALSE' was used to conduct subgroup meta-analysis

* metainc():
  - square root transformed incidence rate difference added as new
    summary measure (sm = "IRSD")

* New arguments 'text.fixed', 'text.random', 'text.predict',
  'text.w.fixed' and 'text.w,random' in meta-analysis functions

* settings.meta():
  - new general setting "geneexpr" to print scientific p-values and
    not calculate confidence interval for between-study heterogeneity
    variance tau2
  - argument 'method.tau.ci' can be specified as a global setting
  - text for fixed effect and random effects model as well as
    prediction interval can be specified (arguments 'text.fixed',
    'text.random', 'text.predict', 'text.w.fixed', 'text.w.randon')

* print.meta(), print.summary.meta():
  - do not print information on continuity correction for exact
    Mantel-Haenszel method with single study

* metareg() can be used in loops to provide argument 'formula'

* New auxiliary function JAMAlabels() to create study labels in JAMA
  layout

### Internal changes

* Calculate measures of residual heterogeneity in hetcalc()


## meta, version 4.15-1 (2020-09-30)

### Bug fixes

* metacr():
  - set summary measure to "OR" for Peto odds ratio


## meta, version 4.15-0 (2020-09-29)

### Major changes

* Deeks' linear regression test for funnel plot asymmetry of funnel
  plots of diagnostic test accuracy studies implemented ([Deeks et.,
  2005](https://doi.org/10.1016/j.jclinepi.2005.01.016))

* *Effective sample size* ([Deeks et.,
  2005](https://doi.org/10.1016/j.jclinepi.2005.01.016)) can be used
  on y-axis of funnel plot
  
* Discard infinite estimates and standard errors from calculation of
  heterogeneity measures

* Diagnostic odds ratio (sm = "DOR") added as new effect measure in
  metabin() and metagen()

### User-visible changes

* forest.meta(), forest.metabind():
  - arguments 'digits.zval' and 'print.zval' renamed to 'digits.stat'
    and 'print.stat'

* print.summary.meta(), settings.meta():
  - argument 'digits.zval' renamed to 'digits.stat'
  
* metacr():
  - do not print a warning for inverse variance meta-analysis with
    binary outcome

* Help page for tests of funnel plot asymmetry updated

* Help pages for metabin() and metainc() updated


## meta, version 4.14-0 (2020-09-09)

### Major changes

* Median and related statistics can be used in meta-analysis with
  continuous outcomes to approximate means and standard deviations
  ([Wan et., 2014](https://doi.org/10.1186/1471-2288-14-135); [Luo
  et al., 2018](https://doi.org/10.1177/0962280216669183); [Shi et
  al., 2020](https://doi.org/10.1002/jrsm.1429))

* RevMan 5 analysis datasets can be imported directly using the
  RM5-file

* R package **xml2** added to Imports (RM5-files are in XML-format)

* Confidence intervals for individual studies can be based on quantile
  of t-distribution (only implemented for mean differences and raw
  untransformed means at the moment)

* For the generic inverse variance method,
  - methods by [Luo et
    al. (2018)](https://doi.org/10.1177/0962280216669183) implemented
    to estimate mean from sample size, median and other statistics
  - method by [Shi et al. (2020)](https://doi.org/10.1002/jrsm.1429)
    implemented to estimate the standard deviation from sample size,
    median, interquartile range and range

### Bug fixes

* forest.meta():
  - show all studies with estimable treatment effects if argument
    'allstudies = FALSE'

* metabind():
  - works with meta-analysis objects created with metacor()
  - calculate correct p-value for heterogeneity test if input are
    subgroup analyses of the same dataset
  - calculate correct p-value for within-subgroup heterogeneity test
    if input are subgroup analyses of the same dataset

* metacum():
  - works with Hartung-Knapp method

* metagen():
  - list element 'seTE' contained standard deviation instead of
    standard error for method by [Wan
    et. (2014)](https://doi.org/10.1186/1471-2288-14-135) to estimate
    mean and its standard error from median and other statistics

### User-visible changes

* read.rm5():
  - direct import of RM5-file possible
  - new argument 'debug' for debug messages while importing RM5-files
    directly

* metacr():
  - overall results not shown if this was specified in the Cochrane
    review (only applies to imported RM5-files)

* metagen(), metacont(), metamean():
  - new argument 'method.mean' to choose method to estimate mean from
    sample size, median and other statistics
  - new argument 'method.sd' to choose method to estimate standard
    deviation from sample size, median, interquartile range and range
  - new argument 'method.ci' to choose method for confidence intervals
    of individual studies (only applies to mean differences and raw
    untransformed means at the moment)

* metacont():
  - new arguments to estimate mean and standard deviation from median
    and related statistics:
	'median.e', 'q1.e', 'q3.e', 'min.e', 'max.e', 'median.c', 'q1.c',
    'q3.c', 'min.c', 'max.c', 'method.mean', 'method.sd',
    'approx.mean.e', 'approx.mean.c', 'approx.sd.e', 'approx.sd.c'

* metamean():
  - new arguments to estimate mean and standard deviation from median
    and related statistics:
	'median', 'q1', 'q3', 'min', 'max', 'method.mean', 'method.sd',
    'approx.mean', 'approx.sd'
  
* forest():
  - by default, show number of participants in forest plot if this
    information is available for meta-analysis objects created with
    metagen()
  - automatically format p-values for individual studies if added to
    forest plot using argument 'leftcols' or 'rightcols'

* Datasets renamed from Fleiss93, Fleiss93cont and Olkin95 to
  Fleiss1993bin, Fleiss1993cont and Olkin1995
 
* More sensible variable names in datasets Fleiss1993bin,
  Fleiss1993cont and Olkin1995

### Internal changes

* Previous R function read.rm5() for CSV-files renamed to
  read.rm5.csv()

* New auxiliary functions extract_outcomes(), oct2txt() and
  read.rm5.rm5() to import RevMan 5 analysis datasets

* ci():
  - list element 'z' renamed to 'statistic' as calculations can also
    be based on the t-distribution; list element 'z' is still part of
    the output for backward compatibility, however, could be removed
    in a future update

* metagen():
  - list elements 'zval', 'zval.fixed' and 'zval.random' renamed to
    'statistic', 'statistic.fixed' and 'statistic.random'; list
    elements 'zval', 'zval.fixed' and 'zval.random' are still part of
    the output for backward compatibility, however, could be removed
    in a future update

* Internal functions TE.seTE.iqr.range(), TE.seTE.iqr() and
  TE.seTE.range() renamed to mean.sd.iqr.range(), mean.sd.iqr() and
  mean.sd.range()

* mean.sd.iqr.range():
  - new arguments 'method.mean' and 'method.sd'

* mean.sd.iqr(), mean.sd.range():
  - new argument 'method.mean'

* chkchar(), chkcolor(), chklevel(), chknumeric():
  - argument 'single' renamed to 'length' (which can be used to test
    for a specific vector length instead whether it is a single value)
	(argument 'single' is still available for backward compatibility,
     however, will be removed in a future update)


## meta, version 4.13-0 (2020-07-02)

### Major changes

* Rely on generic functions from R package **metafor**, e.g., to
  produce forest or funnel plots (since R version 4.0.0 generic
  functions from an R package do not consider corresponding functions
  from another R package which can result in errors if R packages
  **meta** and **metafor** are both loaded)

* R function funnel.default() removed from **meta** due to conflict
  with **metafor**


## meta, version 4.12-0 (2020-05-04)

### Major changes

* Sample size method for meta-analysis of binary data with the odds
  ratio as summary measure implemented ([Bakbergenuly et al.,
  2020](https://www.doi.org/10.1002/jrsm.1404))

* *Ad hoc* variance correction for Hartung-Knapp method in the case of
  very homogeneous study results implemented ([Knapp and Hartung,
  2003](https://www.doi.org/10.1002/sim.1482); [IQWiG, General
  Methods: Draft of Version
  6.0](https://www.iqwig.de/en/about-us/methods/methods-paper/))

* Default settings according to recommendations in [General Methods of
  the Institute for Quality and Efficiency in Health Care (IQWIG),
  Germany](https://www.iqwig.de/en/about-us/methods/methods-paper/)
  
* Do not use predict.rma() from **metafor** package to calculate
  prediction intervals for generalised linear mixed models
  
### User-visible changes

* drapery():
  - study IDs or study labels can be printed at the top of the drapery
    plot to identify individual studies
  - more flexible plots, e.g., colours can be specified for individual
    studies based on p-value of treatment effect
  - possible value for argument 'type' renamed from "cvalue" to
    "zvalue" as drapery plots show test statistics, not critical
    values

* funnel.meta(), funnel.default():
  - argument 'log' is considered for relative summary measures, e.g.,
    odds or risk ratio

* metaprop():
  - can be used with non-integer number of events and sample sizes
  
* metabias.meta(), metabias.default():
  - third component of list element 'estimate' renamed from "slope" to
    "intercept" for linear regression tests
  
* settings.meta():
  - new possible general settings: "iqwig5" and "iqwig6", respectively

* Use Markdown for NEWS


## meta, version 4.11-0 (2020-02-20)

### Major changes

* New arguments 'overall' and 'overall.hetstat' in meta-analysis
  functions to control printing of overall meta-analysis results
  (useful to only show subgroup results)

* For GLMMs, use Wald-type Q statistic to calculate I2 of residual
  heterogeneity in meta-analysis with subgroups (instead of
  likelihood-ratio Q statistic)
  
### Bug fixes

* For GLMMs with subgroups, conduct the correct test for subgroup
  differences (bug was introduced in **meta**, version 4.9-7)

* summary.meta():
  - export the correct harmonic mean for fixed effect and random
    effects model (part of list elements 'fixed' and 'random')

* metabind():
  - do not produce an error if argument 'warn' or 'prediction' is not
    unique in meta-analyses
  
### User-visible changes

* forest.meta():
  - possible to print results for test of an overall effect or
    subgroup differences even if meta-analysis results are not shown
  - new defaults for arguments 'overall' and 'overall.hetstat' (which
    are now considered from meta-analysis objects)

* print.summary.meta():
  - for meta-analysis with subgroups, print information on Q and I^2
    with fixed effect results and information on tau and tau^2 with
    random effects results (previously, information on Q, I^2, tau,
    and tau^2 was reported twice)

### Internal changes

* do not calculate confidence limits for tau2 and tau in intermediate
  calculations of other quantities (i.e., use argument 'method.tau.ci
  = ""')


## meta, version 4.10-0 (2020-01-29)

### Major changes

* New function drapery() to generate a drapery plot which is based on
  p-value curves
  
### Bug fixes

* funnel.meta():
  - print contours in contour-enhanced funnel plots at correct
    position for relative effect measures (bug was introduced in
    **meta**, version 4.9-8)

### User-visible changes

* update.meta():
  - do not print a warning concerning argument 'Q.Cochrane' if
    argument 'sm = "ASD"' for meta-analysis objects created with
    metabin()

* print.summary.meta():
  - do not print z- and p-values if test for an overall effect was not
    conducted; see argument 'null.effect' in metamean(), metaprop(),
    and metarate()


## meta, version 4.9-9 (2019-12-19)

### Bug fixes

* forest.meta():
  - printing an additional column on the right side of the forest plot
    does not result in an error (bug was introduced in **meta**,
    version 4.9-8)

### User-visible changes

* labbe():
  - new argument 'pos.studlab'
  - argument checks implemented

* baujat(), bubble():
  - argument 'pos' renamed to 'pos.studlab'
  - argument checks implemented


## meta, version 4.9-8 (2019-12-16)

### Major changes

* Confidence intervals for the between-study variance tau2 and its
  square root tau are calculated

* Print tau as well as confidence intervals for tau2 and tau in
  outputs

* Square root of between-study variance can be printed in forest plots
  instead of between-study variance tau2; in addition, the confidence
  interval for tau2 or tau can be printed

* Use R package **metafor** to estimate between-study variance tau2
  for DerSimonian-Laird and Paule-Mandel method (which has been
  already used for all other methods to estimate tau2)

* For Mantel-Haenszel (MH) method, report results as MH method
  (instead of inverse variance, IV) for meta-analysis of binary
  outcome with a single study (results are identical for MH and IV
  method in this situation)

* Number of studies printed without digits in forest plots for R
  objects created with metabind()

* P-values can be printed according to [JAMA reporting
  standards](https://jamanetwork.com/journals/jama)

* In subgroup analyses, print the group labels instead of levels if
  the grouping variable is a factor

* In funnel plot, print funnel around random effects (instead of fixed
  effect) estimate if only random effects meta-analysis is conducted;
  only show funnel if either fixed effect or random effects
  meta-analysis was conducted

* New preferred citation of R package **meta**: [Balduzzi et
  al. (2019)](https://scholar.google.com/scholar?q=balduzzi+schwarzer+2019)

### User-visible changes

* print.summary.meta(), forest.meta():
  - new argument 'JAMA.pval' to print p-values according to JAMA
    reporting standards

* print.summary.meta():
  - new argument 'zero.pval' to remove leading zeros from p-values
  - print information on estimation of between-study variance even if
    only results for fixed effect model is shown
  - print information if Mantel-Haenszel estimate is used to calculate
    Q and tau2 (implemented similar to RevMan 5)
  - global setting for 'text.tau2' as defined in settings.meta() is
    considered in details of meta-analytical method

* print.meta():
  - do not print (missing) weights for GLMMs

* update.meta():
  - by default, do not print warnings (argument 'warn')
  - add information on variable defining subgroups (argument 'byvar')
    to meta-analysis dataset

* Command 'settings.meta("JAMA")' will change the settings for
  arguments 'zero.pval' and 'JAMA.pval'

* Help page with description of R package updated

* Major update of other help pages:
  - metacont(), metacor(), and metamean()

### Internal changes

* Function paulemandel() removed as R package **metafor** is used to
  estimate the between-study variance

* formatPT():
  - new argument 'JAMA'

* List elements 'C' and 'C.w' (scaling factor to estimate common
  between-study variance) removed from meta-analysis objects

* Import confint.rma.uni() from **metafor** to calculate confidence
  intervals for tau2 and tau

* New internal function pasteCI() to print formatted CIs

* New internal function is.wholenumber() to check for whole numbers


## meta, version 4.9-7 (2019-09-27)

### Major changes

* Subgroup analysis using argument 'byvar' possible for generalised
  linear mixed models (GLMMs)
    
### Bug fixes

* metaprop():
  - no error if argument 'tau.common' is TRUE for GLMM

* metabin(), metainc(), metarate():
  - consider argument 'control' in subgroup analysis

### User-visible changes

* Major update of help pages:
  - metabin(), metagen(), metainc(), metaprop(), metarate()


## meta, version 4.9-6 (2019-08-06)

### Major changes

* New functions to calculate the number needed to treat from the
  results of a meta-analysis

* Equivalence limits can be added to forest plots

* Font family can be specified in forest plots

* Print Wald-type test of heterogeneity for generalised linear mixed
  models (problem fixed in R package **metafor**, version 2.1-0)
    
### Bug fixes

* forest.meta():
  - (always) print correct length for reference line
  - (always) print label on x-axis at the correct vertical position
  - (always) print graph labels on the left and right side of the
    forest plot at the correct vertical position
  - no error if additional numeric variable is added to the right side
    of the forest plot (argument 'rightcols')

* summary.meta():
  - consider argument 'bylab'

* metaprop():
  - allow values 0 and 1 for argument 'null.effect'

### User-visible changes

* forest.meta():
  - new arguments 'lower.equi', 'upper.equi', 'lty.equi', 'col.e' and
    'fill.equi' to add equivalence limits
  - new argument 'fontfamily' to specify the font family

* forest.metabind():
  - information on heterogeneity printed for each meta-analysis

### Internal changes

* ciAgrestiCoull():
  - set lower confidence limit to 0 for negative values
  - set upper confidence limit to 1 for values above 1

* subgroup meta-analyses return new list element 'pval.Q.w'


## meta, version 4.9-5 (2019-04-11)

### Major changes

* For the generic inverse variance method, treatment estimates and
  standard errors of individual studies can be derived from
  - p-value or confidence limits
  - sample size, median, interquartile range and / or range (Wan et
    al. (2014), BMC Med Res Meth, 14, 135)

* New functions for the conversion of effect measures:
  - smd2or() - from standardised mean difference to log odds ratio
  - or2smd() - from log odds ratio to standardised mean difference

* Harbord test for funnel plot asymmetry implemented for risk ratio as
  effect measure

* Generalised linear mixed model is the new default method for
  meta-analysis of single proportions using the logit transformation

* R packages **metafor** and **lme4** moved from Suggests to Imports

* Suppress printing of Wald-type test of heterogeneity for generalised
  linear mixed models (problem in R function rma.glmm() from R package
  **metafor**, version 2.0-0)
 
* Use **roxygen2** for development of R package **meta**

### User-visible changes

* metagen():
  - new arguments 'pval', 'df', 'lower', 'upper', 'level.ci',
    'median', 'q1', 'q3', 'min', 'max', 'approx.TE', 'approx.seTE' to
    approximate treatment estimates and / or standard errors from
    other information

* forest.meta():
  - printing of leading zeros in p-values can be suppressed (new
    argument 'zero.pval')
  - rounding of values for additional numerical columns possible (new
    arguments 'digits.addcols', 'digits.addcols.left', and
    'digits.addcols.right')
  - argument 'big.mark' is considered for additional columns
  - new arguments 'type.subgroup.fixed', 'type.subgroup.random', and
    'lab.NA.weight'

* settings.meta(), gs():
  - argument names can be abbreviated

* Major update of help pages of metagen() and metaprop()

### Bug fixes

* metacum(), metainf():
  - consider argument 'method' for meta-analysis objects created with
    metaprop() or metarate()

* forest.meta():
  - argument 'studlab' can be used with objects created with metacum()
    or metainf()

* subgroup():
  - return subgroup sample sizes for objects created with metagen()

### Internal changes

* New internal functions TE.seTE.ci(), TE.seTE.iqr(),
  TE.seTE.iqr.range(), TE.seTE.range(), and seTE.ci.pval() to
  approximate treatment estimates or standard errors from other
  information

* setchar():
  - new argument 'stop.at.error'

* metagen():
  - list element 'data' contains the dataset of the meta-analysis
    object (i.e., list element 'data') instead of the whole
    meta-analysis object


## meta, version 4.9-4 (2019-01-02)

### Major changes

* Information on residual heterogeneity in meta-analyses with
  subgroups shown in printouts and forest plots

### User-visible changes

* forest.meta():
  - new arguments 'resid.hetstat' and 'resid.hetlab' to control
    printing of information on residual heterogeneity in meta-analyses
    with subgroups

### Bug fixes

* forest.meta():
  - works in meta-analyses with subgroups if argument 'allstudies =
    FALSE'


## meta, version 4.9-3 (2018-11-29)

### Major changes

* New argument 'control' in meta-analysis functions which is passed on
  to R function rma.uni() or rma.glmm() from R package **metafor** to
  control the iterative process to estimate the between-study variance
  tau^2

### User-visible changes

* metabin(), metacont(), metacor(), metagen(), metainc(), metamean(),
  metaprop(), metarate(), update.meta():
  - new argument 'control' (see major changes)

* forest.meta():
  - new argument 'calcwidth.subgroup'

### Bug fixes

* bubble.metareg():
  - ignore missing values in covariate to calculate limits on x-axis
  - works if dataset used to create meta-analysis object is a tibble
    instead of a data frame

### Internal changes

* metabind():
  - argument 'tau.common' only considered for subgroup analyses

* hetcalc():
  - argument 'control' passed on to R function rma.uni() from R
    package **metafor**

* metacum(), metainf(), subgroup():
  - argument 'control' from meta-analysis objects considered


## meta, version 4.9-2 (2018-06-06)

### Major changes

* All p-values of Q statistics are list elements of meta-analysis
  objects

### Bug fixes

* metareg():
  - consider argument 'intercept = FALSE' if argument 'formula' has
    been provided
    
### Internal changes

* New internal function replaceNULL()


## meta, version 4.9-1 (2018-03-21)

### Major changes

* Subgroup results consider the exclusion of individual studies (bug
  fix)

* For generalised linear mixed models, between-study variance set to
  NA if only a single study is considered in meta-analysis

### Bug fixes

* metamean():
  - use of argument 'byvar' for subgroup analyses possible

* metacor(), metamean(), metaprop(), metarate():
  - use as input to metabind() possible

* Internal function subgroup():
  - consider argument 'exclude' in subgroup analyses

* Internal function bylevs():
  - drop unused levels if subgroup variable is a factor variable

### User-visible changes

* print.summary.meta():
  - print information on Generalised Linear Mixed Model (GLMM) for
    metarate() objects
  - print information on increments added to calculate confidence
    intervals for individual studies (for metarate() with GLMM)

* funnel.meta():
  - new arguments 'ref.triangle', 'lty.ref', 'lwd.ref', 'col.ref', and
    'lty.ref.triangle' to add reference value (null effect) and
    corresponding confidence intervals to the funnel plot

* metabin():
  - new argument 'pscale' to change printout of risk differences

* metainc():
  - new arguments 'irscale' and 'irunit' to change printout of
    incidence rate differences

* forest.meta(), print.meta(), print.summary.meta(), update.meta():
  - consider arguments 'pscale', 'irscale', and 'irunit' for
    meta-analysis objects created with metabin() and metainc()

* print.meta():
  - new argument 'irunit'

### Internal changes

* metaprop():
  - for random effects model, rma.glmm() from package **metafor** is
    called internally with argument 'method = "FE"' if only a single
    study is available

* metareg():
  - for generalised linear mixed models, fallback to fixed effect
    model if number of studies is too small for random effects
    meta-regression

* asin2ir():
  - back transformation could result in (very small) negative zero
    values due to imprecisions (-1e-19); these values are set to zero
    now

* subgroup():
  - code for metamean() added

* chkchar():
  - new argument 'nchar' to test the length of character string(s)

* New internal function is.untransformed() to check for effect
  measures without (back) transformation


## meta, version 4.9-0 (2017-12-06)

### Major changes

* New function metamean() to conduct meta-analysis of single means

* New function metabind() to combine meta-analysis objects, e.g. to
  generate a forest plot with results of several subgroup analyses

* Subgroup analysis implemented for generalised linear mixed models
  (GLMMs) with and without assumption of common between-study variance
  (arguments 'byvar' and 'tau.common')

* Axis direction can be reversed for x-axis in forest plots

* Source code version of **meta** can be installed without
  compilation, i.e., without use of Rtools on Windows or 'Command-line
  tools for Xcode' on macOS

* Rank test for funnel plot asymmetry uses cor() from R package
  **stats** instead of internal C routine (negligibly slower, however,
  no need for compilation of source installs)

* Thousands separator can be used in printouts and forest plots for
  large numbers

* P-values equal to 0 are actually printed as "0" instead of "<
  0.0001"

### User-visible changes

* forest.meta(), print.meta(), print.summary.meta():
  - new argument 'big.mark' to specify character printed as thousands
    separator, e.g., big.mark = "," will result in printing of "1,000"
    instead of "1000"

* forest.meta():
  - sensible forest plot generated if first value in argument 'xlim'
    is larger than second value, e.g. xlim = c(10, -10)
  - separator between label and levels of grouping variable (argument
    'byseparator') is considered from meta-analysis object
  - for relative summary measures, e.g., odds ratio and risk ratio,
    labels on x-axis are not rounded to two digits (which resulted in
    the value 0 for a tick-mark at 0.001)
  - bug fix: lines for treatment effect in fixed effect and random
    effects model start in center of diamond if argument 'hetstat =
    FALSE'
  - bug fix: argument 'type.study' will be sorted according to
    arguments 'sortvar'

* metaprop():
  - arguments 'byvar' and 'tau.common' can be used for GLMMs

* Help page with overview of R functions in R package **meta** updated

### Internal changes

* New internal functions:
  - is.log.effect() to check for treatment effects combined on log
    scale
  - is.mean() to check whether summary measure refers to meta-analysis
    of single means

* Renamed internal functions:
  - formatCI() instead of p.ci()
  - formatN() instead of format.NA()
  - formatPT() instead of format.p()

* Removed R functions:
  - format.tau() as functionality is now provided by formatPT()
  - C program kenscore.c as cor() from R package **stats** is used
    instead to calculate Kendall's tau

* Deprecated functions: format.NA(), format.p(), p.ci()

* Check whether argument 'sm' is NULL in meta-analysis functions

* subgroup(): extended for GLMMs

* formatPT():
  - zero p-values are printed as "0" instead of "< 0.001"
  - NaNs are handled like NAs

* bylabel(), catmeth(), formatPT(), formatN(), xlab():
  - new argument 'big.mark' (see above)


## meta, version 4.8-4 (2017-08-11)

### User-visible changes

* forest.meta():
  - new arguments 'col.fixed' and 'col.random' to change colour of
    fixed effect and random effects lines

### Bug fixes

* bubble.metareg():
  - works if covariate in metareg() is not part of dataset used to
    generate meta-analysis object

* forest.meta():
  - lines for treatment effect in fixed effect and random effects
    model always start in center of diamond
    
* metacum(), metainf():
  - argument 'model.glmm' considered for metabin() and metainc()
    objects

* print.summary.meta():
  - print transformed null effect for meta-analysis of single
    correlations, proportions, or rates if argument 'backtransf =
    FALSE', i.e., for metacor(), metaprop(), and metarate() objects
    
* trimfill.meta():
  - argument 'null.effect' is considered to calculate p-value for
    fixed effect and random effects model for metacor(), metaprop(),
    and metarate() objects

### Internal changes

* New internal functions is.cor(), is.prop() and is.rate() to check
  whether summary measure refers to meta-analysis of correlations,
  proportions, or rates

* metabias.default(), radial.default(), trimfill.default():
  - call metagen() internally to create meta-analysis object
  - call metabias.meta(), radial.meta(), or trimfill.meta() internally
    to conduct analysis


## meta, version 4.8-3 (2017-07-21)

### Major changes

* Similar to RevMan 5, individual studies can be excluded from
  meta-analysis, however, will be shown in printouts and forest plots

* In forest plots, line spacing can be determined by the user

### User-visible changes

* metabin(), metacor(), metacont(), metagen(), metainc(), metaprop(),
  metarate():
  - new argument 'exclude' to exclude studies from meta-analysis

* forest.meta():
  - new argument 'spacing' to determine line spacing
  - bug fix for for meta-analysis with standardized mean difference
    (sm = "SMD") and argument 'layout = "RevMan5"'

* R function ci() can be used with vectors or matrices of treatment
  estimates and standard errors and a single value for argument 'df',
  i.e., degrees of freedom (which is used in R package **netmeta** to
  calculate prediction intervals for network meta-analysis estimates)

* metacum(), metainf():
  - argument 'null.effect' considered internally for objects generated
    with metacor(), metagen(), metaprop() and metarate()

### Internal changes
 
* baujat.meta(), metabias.meta(), metacum(), metainf(), forest.meta(),
  funnel.meta(), metareg(), print.meta(), radial.meta(),
  trimfill.meta(), update.meta():
  - changes to deal with excluded studies


## meta, version 4.8-2 (2017-05-24)

### Major changes
 
* Calculate confidence interval for I2 in a meta-analysis with two
  studies if the heterogeneity statistic Q is larger than 2

* P-values can be printed in scientific notation

* In forest plots, printing of z-values can be disabled and labels for
  tests can be changed by user

### User-visible changes

* forest.meta():
  - new argument 'print.zval' to print (default) or not print z-value
    for test of treatment effect
  - new argument 'print.Q.subgroup' to print (default) or not print
    Chi-squared statistic for test of subgroup differences
  - bug fix: print first line above second line if argument 'xlab'
    consists of two lines (bug was introduced in **meta**, version
    4.8-0)
  - labels of additional columns are printed in correct line if label
    consists of two lines
  - new argument 'scientific.pval' to print p-values in scientific
    notation, e.g., 1.2345e-01 instead of 0.12345
  - arguments 'label.test.overall.fixed', 'label.test.overall.random',
    'label.test.subgroup.fixed', 'label.test.subgroup.random',
    'label.test.effect.subgroup.fixed',
    'label.test.effect.subgroup.random' work as expected
  - new argument 'text.subgroup.nohet' to enable the user to change
    the text "not applicable" in the line with heterogeneity
    statistics for a subgroup with less than two studies contributing
    to the meta-analysis
  - forest plot without any study contributing to meta-analysis can be
    generated without an error, e.g., meta-analysis with binary
    outcome, sm="OR", and all event numbers equal to zero
    
* print.meta() and print.summary.meta():
  - new argument 'scientific.pval' to print p-values in scientific
    notation, e.g., "1.2345e-01" instead of "0.12345"
  - new arguments 'print.pval' and 'print.pval.Q' to specify number of
    significant digits for p-values

* R command 'help(meta)' can be used to show brief overview of R
  package **meta**

* Substantially decrease number of automatically run examples for
  forest.meta() as CRAN only allows a run time below 10 seconds for
  examples provided on a help page
      
### Internal changes

* new internal function pvalQ() to calculate p-value from
  heterogeneity tests

* calcH():
  - Calculate confidence interval for H in a meta-analysis with two
    studies if the heterogeneity statistic Q is larger than 2 (this
    confidence interval is used in isquared() to calculate a
    confidence interval for I2)

* hetcalc():
  - heterogeneity statistic Q set to 0 for a single study contributing
    to the meta-analysis (sometimes in this case Q was set to a value
    below 1e-30)

* subgroup():
  - list element 'df.Q.b' set to 0 if number of studies in
    meta-analysis is 0

* format.p():
  - new argument 'lab.NA' to change value printed for NAs

* forest.meta() and print.summary.meta():
  - use internal function pvalQ() instead of dedicated R code


## meta, version 4.8-1 (2017-03-17)

### User-visible changes

* metacum(), metainf():
  - bug fix for meta-analysis objects without continuity correction,
    i.e., metacont(), metacor(), metagen() (bug was introduced in
    **meta**, version 4.8-0)
  - bug fix for metarate() objects due to improper use of metaprop()
    internally


## meta, version 4.8-0 (2017-03-12)

### Major changes
 
* Continuity correction can be specified for each individual study in
  meta-analysis with proportions or incidence rates

### User-visible changes

* metabin(), metainc(), metaprop(), metarate():
  - argument 'incr' can be of same length as number of studies in
    meta-analysis

* metaprop():
  - bug fix in studies with missing information for events or sample
    size and argument 'method.ci = "CP"'
  - bug fix to calculate test for an overall effect

* forest.meta():
  - bug fix to print summary label (argument 'smlab') above forest
    plot if argument 'fontsize' is unequal to 12
  - by default, label on x-axis and text on top of forest plot are
    printed in center of forest plot (arguments 'xlab.pos',
    'smlab.pos')

* print.summary.meta():
  - print number of studies for fixed effect meta-analysis using
    Mantel-Haenszel method if different from number of studies in
    random effects model (only if summary measure is "RD" or "IRD" and
    at least one study has zero events)

* metainc():
  - bug fix to consider argument 'incr' for incidence rate difference

### Internal changes

* act on NOTE in CRAN checks with R version, 3.4.0, to register and
  declare native C routine 'kenscore'

* metabin(), metainc():
  - new list element 'k.MH' with number of studies in meta-analysis
    using Mantel-Haenszel method

* forest.meta():
  - auxiliary R functions removed from R code
  - cleaning / shortening of R code

* new auxiliary R functions used in forest.meta():
  - add.label(), add.text(), add.xlab(), draw.axis(),
    draw.ci.square(), draw.ci.diamond(), draw.ci.predict(),
    draw.forest(), draw.lines(), formatcol(), removeNULL(), tg(),
    tgl(), twolines(), wcalc()
      
* hetcalc(), calcH():
  - set heterogeneity statistics tau2, H and I2 to NA if only a single
    study contributes to meta-analysis

* updateversion():
  - use R function update.meta() if version of **meta** used to create
    R object is below 3.2


## meta, version 4.7-1 (2017-02-13)

### Major changes

* Null hypothesis for test of an overall effect can be specified for
  metacor(), metagen(), metaprop(), and metarate(); for all other
  meta-analysis functions implicit a null effect of zero is assumed
  (for relative effect measures, e.g., odds ratio and hazard ratio,
  the null effect is defined on the log scale)

* User can choose whether to print the following heterogeneity
  quantities: I^2, H, Rb (by default, heterogeneity measure Rb is not
  printed and thus revoking a change in **meta**, 4.7-0)

* In forest plots with subgroups, study weights are summed up to 100
  percent within each subgroup if no overall estimates are requested,
  i.e., argument 'overall = FALSE' (like before, by default, weights
  are not printed if argument 'overall = FALSE' and have to be
  explicitely requested using argument 'leftcols' or 'rightcols')
  
### User-visible changes

* forest.meta():
  - print line with heterogeneity statistics directly below individual
    study results if pooled effects are not shown in forest plot
    (overall = FALSE)
  - print right and left labels (arguments 'label.left',
    'label.right') in correct line if arguments 'overall = FALSE' and
    'addrow = FALSE'
  - bug fix: do not stop with an error if 'comb.fixed = FALSE',
    'comb.random = FALSE', and 'overall.hetstat = TRUE'

* ci(), metacor(), metagen(), metaprop(), metarate():
  - new argument 'null.effect' to specify null hypothesis for test of
    an overall effect, e.g., null.effect = 0.5 in metaprop() to test
    whether the overall proportion is equal to 0.5

* metagen():
  - Hartung-Knapp method only used for at least two studies in
    meta-analysis

* print.meta():
  - print covariate with subgroup information for each study, if
    subgroup analysis is conducted (argument 'byvar')

* print.summary.meta():
  - new arguments 'print.I2', 'print.H' and 'print.Rb' to specify
    heterogeneity measures shown in output
  - new arguments 'text.tau2', 'text.I2' and 'text.Rb' to change text
    printed to identify respective heterogeneity measure
  - only print information on double zero studies if argument
    'allstudies = TRUE'
  - print results for (empty) subgroup in meta-analysis with two
    studies and one subgroup with missing treatment estimate

* settings.meta():
  - new arguments 'print.I2', 'print.H', 'print.Rb', 'text.tau2',
    'text.I2' and 'text.Rb' to modify printing of heterogeneity
    measures

### Internal changes
    
* summary.meta():
  - bug fix renaming list element 'ircale' renamed to 'irscale'
  - list element 'within' removed which has not been used since
    **meta**, version 1.1-4


## meta, version 4.7-0 (2016-12-16)

### Major changes

* Forest plots:
  - forest plots with RevMan 5 and JAMA layout
  - use of mathematical symbols for I^2, tau^2, etc.
  - individual study results can be omitted from forest plot
    (especially useful to only print subgroup results)
  - labels can be printed at top of forest plot

* Measure of between-study heterogeneity added:
  - R_b ([Crippa et al. (2016)](https://www.doi.org/10.1002/sim.6980))

* Default settings of meta-analysis methods specified via gs() instead
  of extracting elements of list .settings (which makes output of
  args() easier to read, e.g., args(metabin))

* Version of suggested R package **metafor** must be at least 1.9-9
  (due to change in arguments of rma.uni() and rma.glmm())
  
### User-visible changes

* forest.meta():
  - argument 'layout':
      - new layout "JAMA" to produce forest plots with [JAMA
        style](https://jamanetwork.com/journals/jama/pages/instructions-for-authors/)
      - RevMan 5 layout extended
  - arguments can be specified without using grid::unit():
    'plotwidth', 'colgap', 'colgap.left', 'colgap.right',
    'colgap.studlab', 'colgap.forest', 'colgap.forest.left',
    'colgap.forest.right'
  - new argument 'study.results' to print (default) or omit individual
    study results from forest plot
  - new argument 'bottom.lr' to change position of labels on left and
    right side of forest plot
  - new arguments 'col.label.right' and 'col.label.left' to change
    colour of labels on left and right side of forest plot
  - argument 'weight' renamed to 'weight.study' and new argument
    'weight.subgroup' added to specify whether plotted subgroup
    results should be of same or different size
  - new arguments 'print.Rb', 'print.Rb.ci' and 'Rb.text' for
    heterogeneity measure Rb
  - new arguments to control printing: 'digits.cor', 'digits.mean',
    'digits.sd', 'digits.time', 'digits.zval'
  - new argument 'print.subgroup.labels' to print (default) or omit
    rows with subgroup label from forest plot
  - new argument 'type.subgroup' to change plotting of subgroup
    results
  - argument 'addspace' renamed to 'addrow'
  - new argument 'addrow.subgroups' to add a blank line between
    subgroup results
  - new argument 'addrow.overall' to add a blank before meta-analysis
    results
  - new argument 'blanks' to enhance printing of test statistics,
    heterogeneity measures, and p-values
  - new argument 'colgap.studlab' to specify space between column with
    study labels and subsequent column
  - new arguments to change width of column with study labels (these
    arguments are especially useful if only study labels are printed
    on left side of forest plot):
      - 'calcwidth.fixed' (consider text for fixed effect model)
      - 'calcwidth.random' (consider text for random effects model)
      - 'calcwidth.hetstat' (consider text for heterogeneity measures)
      - 'calcwidth.tests' (consider text for tests of effect or
        subgroup differences)
  - new column "effect.ci" with estimated treatment effect and
    confidence interval in one column
  - unnecessary arguments removed: 'text.I2', 'text.tau2'
    
* metabin(), metacont(), metacor(), metacr(), metacum(), metagen(),
  metainc(), metainf(), metaprop(), metarate(), trimfill.default(),
  trimfill.meta():
  - new measure of between-study heterogeneity implemented (list
    elements 'Rb', 'lower.Rb', 'upper.Rb')
  
* summary.meta():
  - new measure of between-study heterogeneity added (list element
    'Rb.w')
  
* print.meta(), print.summary.meta():
  - print heterogeneity measure Rb
  
* metabias.meta(), metabias.default():
  - checks for arguments implemented
      
* New function gs() to get default settings

* forest.meta(), metabin(), metacont(), metacor(), metacr(),
  metagen(), metainc(), metaprop(), metarate(), print.meta(),
  print.summary.meta():
  - use gs() to define defaults for arguments in meta-analysis
    functions, e.g. gs("hakn") instead of .settings$hakn

* metareg():
  - stop with an error if version of **metafor** package is below
    1.9-9

* metabin(), metainc(), metaprop(), metarate():
  - for GLMMs, stop with an error if version of **metafor** package is
    below 1.9-9
    
* metabin():
  - bug fix, do not stop with an error if no double zero events are
    present in a dataset with at least one study with NA event counts

* metareg():
  - bug fix, use of covariate 'x' does not result in an error

* settings.meta():
  - general settings for RevMan 5 and JAMA implemented
  - function can be used to change the layout of confidence intervals
    using arguments 'CIbracket' and 'CIseparator' which can also be
    set using cilayout()

* Several help pages updated, especially
  - forest.meta(), settings.meta(), meta-package

### Internal changes
    
* metabin(), metainc(), metaprop(), metarate(), metareg():
  - use argument 'test' instead of 'knha' and 'tdist' for calls of
    rma.uni() and rma.glmm(); change in R package **metafor**, version
    1.9-9
    
* subgroup():
  - new measure Rb of between-study heterogeneity implemented
    
* is.installed.package():
  - new check of version number of R package
  - use requireNamespace() instead of installed.packages()

* format.p():
  - for small p-values, print "p < 0.01" or "p < 0.001" instead of "p
    < 0.0001" if digits.pval is 2 or 3, respectively
  - new argument 'zero' to print ".001" instead of "0.001", etc

* meta-internal():
  - set defaults for arguments 'smrate' and 'layout'


## meta, version 4.6-0 (2016-10-12)

### Major changes

* New function metarate() to conduct meta-analysis of single incidence
  rates

* Peters' test for funnel plot asymmetry implemented for
  meta-analysis of single proportions

* Meta-analysis of ratio of means added to metacont()

* Justification of additional columns in forest plot can be
  specified individually for each additional column

* Justification of additional columns in forest plot can be
  specified individually for each additional column

* Calculation of Freeman-Tukey double arcsine transformation and
  back transformation slightly changed in meta-analysis of single
  proportions

* By default, do not print a warning if back transformation for
  metaprop() and metarate() objects results in values below 0 or
  above 1 (only for proportions); note, respective values are set
  to 0 or 1

### User-visible changes
    
* Help page with brief overview of **meta** package added
    
* Preferred citation of **meta** package in publications changed; see
  output of command 'citation("meta")'

* forest.meta(), metagen(), print.meta(), print.summary.meta(),
  summary.meta(), trimfill.default(), trimfill.meta(), update.meta():
  - new arguments 'irscale' and 'irunit' for meta-analysis objects
    created with metarate()

* settings.meta():
  - new arguments 'smrate' for meta-analysis objects created with
    metarate()

* funnel.meta(), funnel.default():
  - new argument 'pos.studlab' to change position of study labels

* forest.meta():
  - new arguments 'just.addcols.left' and 'just.addcols.right' to
    specify justification of additional columns on left and right side
    of forest plot

* metacont():
  - meta-analysis for ratio of means implemented (argument 'sm =
    "ROM"')
  - new argument 'backtransf' for ratio of means (argument 'sm =
    "ROM"')

* metaprop():
  - change in Freeman-Tukey double arcsine transformation only visible
    in printouts if argument 'backtransf = FALSE' or if list elements
    'TE', 'TE.fixed', and 'TE.random' (as well as confidence
    intervals) are extracted from a metaprop object

* print.summary.meta():
  - bug fix in subgroup() to print correct results for subgroup
    analyses of metaprop objects with argument 'sm = "PFT"'

* print.meta(), print.summary.meta():
  - new argument 'warn.backtransf' to specify whether a warning should
    be printed if backtransformed proportions and rates are below 0
    and back transformed proportions are above 1

* Help pages updated:
  - forest.meta(), metabias.meta(), metabin(), metacont(), metacor(),
    metagen(), metainc(), metainf(), metaprop(), print.meta(),
    print.summary.meta(), summary.meta(), trimfill.default(),
    trimfill.meta(), update.meta()

### Internal changes

* New function asin2ir() to back transform arcsine transformed
  incidence rates

* backtransf(), catmeth(), metacum(), metainf(), subgroup(), xlab():
  - extension to handle meta-analysis objects created with metarate()

* metaprop(), asin2p():
  - calculation of Freeman-Tukey double arcsine transformation changed
    to get similar estimates as arcsine transformation, i.e. multiply
    values by 0.5

* subgroup():
  - bux fix in calculation of harmonic mean of sample sizes for
    metaprop() objects with argument 'sm = "PFT"' and event times for
    metarate() objects with argument 'sm = "IRFT"'


## meta, version 4.5-0 (2016-08-17)

### Major changes

* New features in forest plots:
  - printing of columns on left side of forest plot can be omitted
  - total person time can be printed
  - text for fixed effect and random effects model can be omitted from
    calculation of width for study labels
  - plot type for confidence intervals (square or diamond) can be
    specified for each study as well as fixed effect and random
    effects estimate
  - printing of test for treatment effect in subgroups possible

* New function weights.meta() to calculate absolute and percentage
  weights in meta-analysis

* New argument 'byseparator' to define the separator between label and
  subgroup levels which is printed in meta-analysis summaries and
  forest plots - considered in all R functions dealing with
  meta-analysis and subgroups

* Argument 'pscale' - a scaling factor for printing of single event
  probabilities - considered in all R functions for single
  proportions; before this update, argument 'pscale' was only
  available in forest.meta()

### User-visible changes

* forest.meta():
  - argument 'ref' considered for metaprop() objects
  - argument 'leftcols = FALSE' omits printing of columns on left side
    of forest plot
  - new argument 'pooled.times' to print total person time
  - new argument 'calcwidth.pooled' to include or exclude text from
    pooled estimates to determine width of study labels
  - new argument names (old names can still be used at the moment,
    however, will result in an informative warning message):
    - 'col.i'                    -> 'col.study'
    - 'col.i.inside.square'      -> 'col.inside'
    - 'col.diamond.fixed.lines'  -> 'col.diamond.lines.fixed'
    - 'col.diamond.random.lines' -> 'col.diamond.lines.random'
  - new arguments:
    - 'type.study', 'type.fixed', 'type.random' to use squares or
      diamonds to plot treatment effects and confidence intervals
    - 'col.inside.fixed', 'col.inside.random' with information on
      colour to print confidence interval inside square
    - 'test.effect.subgroup', 'test.effect.subgroup.fixed',
      'test.effect.subgroup.random',
      'label.test.effect.subgroup.fixed',
      'label.test.effect.subgroup.random', 'fs.test.effect.subgroup',
      'ff.test.effect.subgroup' to print results for test of treatment
      effect in subgroups
  - bug fix to get correct length for reference line and lines for
    fixed effect and random effects estimate if argument 'test.overall
    = TRUE'
  - bug fix to consider arguments 'lab.e.attach.to.col' and
    'lab.c.attach.to.col' for metagen() objects

* metabin(), metacont(), metacor(), metagen(), metainc(), metaprop(),
  forest.meta(), print.summary.meta(), summary.meta(), update.meta(),
  settings.meta():
  - new argument 'byseparator'

* metagen(), metaprop(), print.meta(), print.summary.meta(),
  summary.meta(), trimfill.meta(), trimfill.default(), update.meta():
  - new argument 'pscale'

* labbe.metabin(), labbe.default():
  - transformed event probabilites can be plotted, e.g., log odds
    event probabilities for odds ratio as summary measure
  - line for null effect added by default; see arguments 'nulleffect',
    'lwd.nulleffect', 'col.nulleffect'

* metabin(), metainc(), metaprop():
  - use predict.rma() from **metafor** package to calculate prediction
    interval for GLMM method
  - print note for GLMM method that continuity correction is only used
    to calculate individual study results
    
* Help pages updated:
  - labbe.metabin(), labbe.default(), forest.meta(), metabin(),
    metacont(), metacor(), metagen(), metainc(), metaprop(),
    print.meta(), print.summary.meta(), summary.meta(),
    trimfill.meta(), trimfill.default(), update.meta()

### Internal changes

* New function bylabel() to print subgroup labels

* update.meta():
  - do not consider columns 'n.e' and 'n.c' as sample sizes for
    metagen() or metainc() object if not used in original call

* catmeth(), xlab():
  - new argument 'pscale'

* catmeth():
  - for GLMM, print information that continuity correction is only
    used to calculate individual study results

* metacum(), metainf():
  - list element 'pscale' added

* metabin():
  - list elements 'incr.e' and 'incr.c' contain zeros for Peto method
  - print warning that no continuity correction is used for Peto
    method if any of the following arguments is used: 'incr',
    'allincr', 'addincr', 'allstudies'

* metacr():
  - keep dataset used to conduct meta-analysis in list element 'data'

* paulemandel():
  - bug fix for error if used with a single study

* settings.meta():
  - bug fix to work as expected for argument 'method.tau'


## meta, version 4.4-1 (2016-06-20)

### User-visible changes

* metareg(), update.meta():
  - bug fix for error if used with metaprop() object and argument
    'method = "GLMM"'


## meta, version 4.4-0 (2016-05-13)

### Major changes

* Generalised linear mixed models (GLMMs) implemented by internal call
  of rma.glmm() from R package **metafor** by Wolfgang Viechtbauer

* R packages **lme4**, **numDeriv**, and **BiasedUrn** added to
  suggested packages which are required by rma.glmm()

* Print layout (especially number of printed digits) slightly modified
  which impacts output from print.meta(), print.summary.meta(), and
  forest.meta()

* New arguments to change number of digits in printouts and forest
  plots

### User-visible changes

* metabin(), metainc(), metaprop():
  - extension for meta-analysis based on GLMM; see argument 'method'
    and 'model.glmm' (not used in metaprop())
  - new argument '...' to provide additional arguments to rma.glmm()
  - some arguments can be used for other meta-analysis methods than
    inverse variance method: 'method.tau', 'hakn', 'tau.common',
    'TE.tau, 'tau.preset'

* metabin():
  - do not print warning that inverse variance instead of
    Mantel-Haenszel method is used for analysis of a single study
  - print warning if continuity correction (arguments 'incr',
    'allincr', 'addincr', 'allstudies') is used with arcsine
    difference, Peto method, or GLMM
  - check whether R package **BiasedUrn** is installed for conditional
    hypergeometric-normal GLMM (method = "GLMM", model.glmm = "CM.EL")

* forest.meta():
  - extension to plot meta-analysis based on GLMM
  - argument 'labels' can be used instead of argument 'label' to
    change labels on x-axis

* funnel.meta():
  - print default labels on y-axis with capital first letter

* metareg() and update.meta():
  - extension for meta-analysis based on GLMM

* print.meta():
  - new arguments to control printing: 'digits.se', 'digits.zval',
    'digits.Q', 'digits.tau2', 'digits.H', 'digits.I2', 'digits.prop',
    'digits.weight'
  - argument '...' passed on to internal call of print.summary.meta()

* print.summary.meta():
  - new arguments to control printing: 'digits.zval', 'digits.Q',
    'digits.tau2', 'digits.H', 'digits.I2'
  - print "--" for missing z-value instead of "NA"
  - only print confidence interval for H and I2 if lower and upper
    limits are not NA
  - print Wald-type and Likelihood-Ratio heterogeneity test for GLMMs

* settings.meta():
  - new arguments: 'model.glmm', 'digits', 'digits.se', 'digits.zval',
    'digits.Q', 'digits.tau2', 'digits.H', 'digits.I2', 'digits.prop',
    'digits.weight', 'digits.pval', 'digits.pval.Q'
  - check whether R package **metafor** is installed for specific
    values of argument 'method.tau'
  - check whether R packages required for GLMMs are available (if
    method = "GLMM"): **metafor**, **lme4**, **numDeriv**

* Help pages updated:
    metabin(), metainc(), metaprop(), metareg(), forest(),
    print.meta(), print.summary.meta(), settings.meta(), update.meta()

### Internal changes

* New function format.NA() to print other text than "NA" for missing
  values
  
* metagen():
  - only call paulemandel() if heterogeneity statistic Q is larger
      equal than number of studies minus 1
      (otherwise between-study heterogeneity tau2 is set equal to 0)

* metabin(), metainc(), metaprop(), summary.meta():
  - new list elements 'model.glmm', '.glmm.fixed', '.glmm.random',
    'version.metafor'

* metabin(), summary.meta():
  - new list element 'doublezeros' for odds ratio or risk
    ratio as summary measure

* Set defaults for arguments 'model.glmm', 'digits', 'digits.se',
  'digits.zval', 'digits.Q', 'digits.tau2', 'digits.H', 'digits.I2',
  'digits.prop', 'digits.weight', 'digits.pval', 'and digits.pval.Q'
       
* paulemandel():
  - more sensible warning if maximum number of iterations is reached
  - maximum number of iterations increased from 25 to 100

* format.p():
  - print trailing zeros

* catmeth():
  - print information for GLMMs
  - print information whether studies with double zeros are included
    in meta-analysis

* is.installed.package():
  - new arguments for more flexible error and warning messages:
    'func', 'argument', 'value', 'chksettings'


## meta, version 4.3-2 (2015-12-02)

* metacont():
  - bug fix to calculate correct treatment estimates for individual
    studies for Glass's delta

* metaprop():
  - print correct error message if number of events is larger than
    number of observations


## meta, version 4.3-1 (2015-11-13)

* forest.meta():
  - new arguments 'digits.se', 'digits.tau2', 'digits.pval',
    'digits.pval.Q', 'digits.Q', 'digits.I2' to control printing of
    standard errors, p-values, tau2 and heterogeneity statistics
  - new arguments 'test.overall' and 'test.subgroup' controlling
    whether information on test for overall effect and heterogeneity
    should be printed

* Internal function paulemandel():
  - bug fix to give studies with missing treatment effect and standard
    error zero weight in random effects meta-analysis
  - do not stop estimation algorithm if estimated tau2 is negative

* settings.meta():
  - bug fix for error if used with an unassigned argument

* format.p(), format.tau():
  - new argument 'digits' to round p-values and tau2 values

* chkchar(), chkclass(), chklength(), chklevel(), chklogical(),
  chkmiss(), chknull(), chknumeric(), setchar():
  - new argument 'name' to change name of checked argument in printout

* Help page of forest.meta() updated


## meta, version 4.3-0 (2015-07-02)

* metabin(), metainc(), and metaprop():
  - allow missing values in numbers of events or patients
    (corresponding studies get zero weight in meta-analysis)

* forest.meta():
  - print information on test for overall effect (arguments
    'test.overall.fixed' and 'test.overall.random')
  - print information on test for subgroup differences in
    meta-analysis with subgroups (arguments 'test.subgroup.fixed' and
    'test.subgroup.random')
  - new argument 'layout' to change layout of forest plot
  - argument 'lab.NA' considered for all columns in forest plot, e.g.,
    numbers of events and patients for metabin()
  - new argument 'lab.NA.effect' to label NAs in individual treatment
    estimates and confidence intervals
  - bug fix for error if random effects estimate is missing

* metareg():
  - additional arguments implemented ('hakn', 'level.comb',
    'intercept')
  - argument '...' is no longer ignored but passed on to rma.uni(),
    e.g., to control the iterative estimation process
  - bug fix to conduct fixed effect meta-regression (argument
    'method.tau = "FE"')

* metabin():
  - use inverse variance instead of Mantel-Haenszel method if only a
    single study has a non-missing treatment estimate or standard
    error

* settings.meta():
  - code added for new arguments in forest.meta() to print information
    on tests

* Help pages of metareg() and forest.meta() and link to RevMan webpage
  updated


## meta, version 4.2-0 (2015-05-08)

* Copyright changed (new names for Institute and Medical Center)

* metacont():
  - new argument 'exact.smd' to implement exact formulae for Hedges' g
    and Cohen's d (White and Thomas (2005; Hedges, 1981)
  - use formula from Borenstein et al. (2009) to calculate standard
    error for Cohen's d

* forest.meta():
  - bug fix to appropriately sort additional columns provided in
    arguments 'leftcols' and 'rightcols' if argument 'sortvar' is not
    missing
  - new argument 'print.I2.ci' to print confidence intervals for I2

* forest.meta(), print.meta, print.summary.meta():
  - prediction interval can be printed if random effects estimate is
    not shown

* settings.meta(), catmeth(), update.meta():
  - code added for new argument 'exact.smd' in metacont()

* ci(), kentau():
  - calculate p-values without floating point number representation
    problems, e.g., the command ci(9, 1) does not result in a p-value
    of 0 but 2.257177e-19

* Several help pages updated to reflect changes in metacont() and
  RevMan 5 reference


## meta, version 4.1-0 (2015-02-04)

* Title of R package changed

* metacont():
  - new argument 'method.smd' to implement Cohen's d (argument
    'method.smd = "Cohen"') and Glass' delta ('method.smd = "Glass"')
    as additional effect measures for the standardised mean difference
    ('sm = "SMD"')
  - new argument 'sd.glass' to choose the denominator for Glass' delta

* update.meta():
  - new arguments 'method.smd' and 'sd.glass' added

* summary.meta():
  - information for new arguments 'method.smd' and 'sd.glass' added to
    summary.meta object

* settings.meta():
  - code added for new arguments 'method.smd' and 'sd.glass' in
    metacont()

* forest.meta():
  - bug fix for staggered point estimates in metaprop() object with
    subgroups

* metagen():
  - bug fix to give studies with missing treatment effect but
    available standard error zero weight in meta-analysis

* paulemandel():
  - only consider studies without missing treatment effect and
    standard error in calculation of between-study variance

* chklevel():
  - print meaningful error message if confidence limit is outside the
    range [0, 1]

* catmeth():
  - print information on method to estimate standardised mean
    difference in metacont()

* Help pages updated for metacont() and update.meta()


## meta, version 4.0-3 (2015-01-06)

* metabin():
  - bug fix for error in printing of results for Mantel-Haenszel or
    Peto method if any study has zero events in both groups


## meta, version 4.0-2 (2014-12-06)

* metabin():
  - bug fix for error if Peto method is used
  - argument 'sm = "ASD"' for arcsine difference instead of 'sm =
    "AS"' (abbreviations 'sm = "AS"' and 'sm = "A"' can still be used)

* metabin(), metacont(), metacor(), metagen(), metainc(), and
  metaprop():
  - weights 'w.random.w' are calculated from random effects
    meta-analysis ignoring subgroup membership; internal function
    subgroup() changed accordingly
  - argument 'tau.common = TRUE' if argument 'tau.preset' is not NULL
    in subgroup analyses


## meta, version 4.0-1 (2014-11-19)

* forest.meta():
  - bug fix for meta-analyses with subgroups if additional columns
    were provided in argument 'leftcols' or 'rightcols'


## meta, version 4.0-0 (2014-11-19)

### Major revision

This update has been declared as major revision as R code to conduct
subgroup analyses has been moved from summary.meta() and forest.meta()
to metabin(), metacont(), metacor(), metagen(), metainc(), and
metaprop(). Accordingly, an R object generated with these functions
contains all results from subgroup analyses.

In the case of subgroups, the overall treatment effect in fixed effect
and random effects meta-analysis ignores subgroup membership. See
Borenstein et al. (2011), Introduction to Meta-Analysis, Wiley,
Chapter 19, "Obtaining an overall effect in the presence of subgroups,
Option 3.

Furthermore, several checks of function arguments have been
implemented in version 4.0-0 of meta.

### Details

* Function addvar() removed from R package **meta** as functionality
  is provided by forest.meta()

* forest.meta():
  - new meaning for argument 'just' which determines the justification
    of all columns but study labels (argument 'just.studlab') and
    columns added to the forest plot (argument 'just.addcols')
  - new argument 'just.addcols' to change justification of text in
    additional columns
  - new arguments 'text.I2' and 'text.tau2'
  - for metaprop objects, values "n" and "event" handled as standard
    columns in argument 'rightcols' and 'leftcols', i.e. justification
    is determined by argument 'just.cols'
  - subgroup results printed with the same polygon height as overall
    results, i.e. percentage weight is not considered to determine
    polygon height for subgroups

* bubble.metareg():
  - bug fix for meta-regression without intercept
  - bug fix for error in meta-regression using specific effect
    measure, e.g. 'sm = "RR"', "OR", or "HR"

* New internal R functions:
  - subgroup(), hetcalc()
  - updateversion()
  - bylevs(), byvarname()
  - chkchar(), chkclass(), chklength(), chklevel(), chklogical(),
    chkmetafor() chkmiss(), chknull(), chknumeric()
  - int2num(), npn()
  - setchar(), setstudlab()

* format.p(), format.tau(), catmeth(), print.summary.meta():
  - consider settings for option 'OutDec' (character used as decimal
    point in output conversions), e.g., options(OutDec = ",") will
    print "1,0" instead of "1.0"

* print.meta(), print.summary.meta():
  - print 'p-value' instead of 'p.value'

* print.summary.meta():
  - remove code for R objects created with version 2.0-0 or lower of
    **meta**

* Several help pages updated


## meta, version 3.8-0 (2014-09-12)

* forest.meta(), funnel.default(), funnel.meta(), metabin(), metacor,
  metacr(), metagen(), metainc(), metaprop(), print.meta(),
  print.summary.meta, summary.meta(), trimfill.default(),
  trimfill.meta():
  - new argument 'backtransf' indicating whether effect measures
    should be back transformed

* print.meta(), print.summary.meta():
  - argument 'logscale' replaced by 'backtransf'

* print.summary.meta(), forest.meta():
  - print prediction interval for Freeman-Tukey double arcsine
    transformation (sm = "PFT")

* forest.meta():
  - consider prediction interval to calculate limits on x-axis if
    argument 'prediction = TRUE'

* bubble.metareg():
  - new argument 'regline' indicating whether regression line should
    be added to plot

* settings.meta():
  - new argument 'print' to print listing of all settings as function
    call without arguments does not print settings any longer
  - list with previous settings can be provided as sole input

* New functions:
  - backtransf() to control back transformation of effect measures
  - is.relative.effect() to check for relative effect measures

* File DESCRIPTION:
  - R package **grid** defined as Imports instead of Depends

* Help pages updated to reflect changes in version 3.8-0


## meta, version 3.7-1 (2014-07-29)

* forest.meta():
  - bug fix to correctly sort lower and upper confidence interval
    limits if argument 'sortvar' is used (bug was introduced in
    **meta**, version 3.7-0)
  - argument 'sortvar' works without reference to meta-analysis
    object, e.g., command forest(meta1, sortvar = TE) can be used
    instead of forest(meta1, sortvar = m1$TE)

* Help page of forest.meta():
  - examples using argument 'sortvar' added


## meta, version 3.7-0 (2014-07-11)

* metaprop():
  - new argument 'method.ci' to implement various methods to calculate
    confidence intervals for individual studies (default:
    Clopper-Pearson method which is also called 'exact' binomial
    method)
  - list elements 'zval.fixed', 'pval.fixed', 'zval.random' and
    'pval.random' set to NA

* New internal functions:
  - ciWilsonScore()      used in metaprop()
  - ciAgrestiCoull()     used in metaprop()
  - ciSimpleAsymptotic() used in metaprop()
  - estimate.missing()   used in trimfill.meta() and trimfill.default()

* metacont():
  - new argument 'pooledvar' to conduct meta-analysis of mean
    differences based on pooled variance for individual studies

* update.meta():
  - function can be used to upgrade R objects created with older
    versions of **meta**, i.e. all versions between 0.5 and 3.6-0
  - extended to objects of the following classes:
	- trimfill()
	- metacum()
	- metainf()
  - new arguments:
	- 'method.ci' for metaprop() objects
	- 'pooledvar' for metacont() objects
	- 'left', 'ma.fixed', 'type' and 'n.iter.max' for trimfill()
      objects
  - new list element 'call.object' with call used to generate
    meta-analysis object

* as.data.frame.meta(), baujat.meta(), forest.meta(), funnel.meta(),
  labbe.metabin(), metacum(), metainf(), print.meta(), summary.meta,
  trimfill.meta():
  - call update.meta() to update meta-analysis objects created with
    **meta**, version < 3.7

* metabin(), metacont(), metacor(), metagen(), metainc(), metaprop(),
  trimfill.default(), trimfill.meta():
  - new list elements 'lower', 'upper', 'zval' and 'pval' with
    confidence limits, z- and p-values for individual studies

* print.meta(), print.summary.meta():
  - print information on method used for confidence intervals of
    individual studies

* metacum(), metainf():
  - add code for metainc() objects
  - new list element 'call' with function call
  - consider argument 'pooledvar' for metacont() objects

* metabin(), metacont(), metacor(), metagen(), metainc(), metaprop():
  - study labels will only be converted to characters for factor
    variables

* Help pages
  - updated to reflect changes in version 3.7-0
  - argument 'tau.preset' correctly described as the _square-root_ of
    the between-study variance


## meta, version 3.6.0 (2013-05-27)

* New functions:
  - baujat(), baujat.meta() for Baujat plot to explore heterogeneity
    in meta-analysis
  - bubble(), bubble.metareg() for bubble plot to display the result
    of a meta-regression

* metareg():
  - class 'metareg' added
  - new list element '.meta' with meta-analysis object used in
    function call

* update.meta():
  - argument 'studlab' fully functional (bug was introduced in
    **meta**, version 3.2-0)

* print.meta():
  - print study label for a single study in meta-analysis if argument
    'details = TRUE'; data.frame()) instead of cbind() used internally

* New internal function is.installed.package() replaces
  is.installed.metafor()

* Help pages datasets amlodipine and cisapride:
  - execute examples for Hartung-Knapp method

* Help pages merged:
  - forest(), forest.meta()
  - funnel(), funnel.meta()
  - labbe, labbe.metabin()
  - metabias(), metabias.meta()
  - trimfill(), trimfill.meta()


## meta, version 3.5-1 (2014-05-14)

* metabin():
  - inverse variance method used instead of Mantel-Haenszel method if
    argument 'tau.common = TRUE'

* metareg():
  - tilde sign not necessary in argument 'formula' to make this
    function more user friendly

* forest.meta():
  - print common tau2 for subgroups if argument 'tau.common = TRUE' in
    meta-analysis object

* metagen():
  - arguments 'n.e' and 'n.c' can be part of the dataset provided in
    argument 'data'
  - DerSimonian-Laird method used instead of Paule-Mandel method if
    argument 'tau.common = TRUE'

* metacor(), metainc(), and metaprop():
  - store value of arguments 'title', 'complab', and 'outclab' in
    meta-analysis object

* Some help pages (slightly) updated


## meta, version 3.5-0 (2014-04-19)

* New R function settings.meta() to define and print default settings
  for meta-analyses in R package **meta**

* metagen():
  - Hartung and Knapp method added; previously rma.uni() from R
    package **metafor** was called for this method
  - Paule-Mandel method to estimate between-study variance implemented
    using new internal function paulemandel() which is based on
    mpaule.default() from R package **metRology** by S.L.R. Ellison
    <s.ellison at lgc.co.uk> (Author of mpaule.default() is S. Cowen
    <simon.cowen at lgc.co.uk> with amendments by S.L.R. Ellison)

* metacont():
  - studies with missing treatment estimate get zero weight in
    meta-analysis

* metabin(), metacont(), metacor(), metacr(), metagen(), metainc(),
  metaprop():
  - default values changed according to settings.meta()

* metareg():
  - use argument 'method.tau = "REML"' if this argument is equal to
    "PM" for meta-analysis object

* Several help pages updated


## meta, version 3.2-1 (2014-03-26)

* forest.meta():
  - bug fix to show correct confidence limits for individual studies
    if argument 'level' is not equal to the default 0.95. (bug was
    introduced in **meta**, version 3.0-0)


## meta, version 3.2-0 (2014-03-12)

* metabin(), metacont(), metacor(), metagen(), metainc(), metaprop():
  - heterogeneity statistics I2 and H added to R object
  - column names changed in list element 'data'; columns starting with
    a "." used internally in update.meta()
  - string "byvar" is used as default label for grouping
    variable if argument 'bylab' is not provided

* metareg():
  - variable '.byvar' used instead of 'byvar' to reflect change in
    list element 'data'

* update.meta():
  - arguments 'byvar' and 'subset' fully functional
  - variables '.TE', ... used internally instead of TE, ... to reflect
    change in list element 'data'

* trimfill.default(), trimfill.meta():
  - heterogeneity statistics I2 and H added to R object

* metagen():
  - bug fix to correctly calculate weights (list elements 'w.fixed'
    and 'w.random') if any standard error is missing or zero for the
    Hartung-Knapp method (argument 'hakn = TRUE') or the DerSimonian
    Laird method is not used (argument 'method.tau' not equal to "DL")

* summary.meta():
  - subgroup analysis implemented for metainc() objects

* forest.meta():
  - groups will not be sorted automatically in alphabetical order (new
    argument 'bysort'). Use argument 'bysort = FALSE' to get the old
    behaviour of forest.meta()

* forest.meta(), summary.meta():
  - only (re)calculate heterogeneity statistics (Q, tau2, I2) for R
    objects generated with older versions of R package **meta**

* catmeth():
  - new argument 'tau.preset' to print information if between-study
    variance was pre-specified

* print.meta(), print.summary.meta():
  - argument 'tau.preset' used in catmeth()

* New internally used functions isquared() and calcH()

* Some help pages updated


## meta, version 3.1-2 (2013-12-01)

* forest.meta():
  - bug fix for error in meta-analyses with subgroups using any but
    metaprop() (bug was introduced in **meta**, version 3.1-1)


## meta, version 3.1-1 (2013-11-19)

* forest.meta():
  - bug fix to show random effects estimate in metaprop() objects with
    subgroups using argument 'sm = "PFT"'


## meta, version 3.1-0 (2013-11-11)

* New R function metainc() for meta-analysis of incidence rates

* Continuity correction:
  - metabin() and metaprop() do no longer print a warning in case of
    studies with a zero cell frequency
  - instead information on continuity correction is given under
    "Details on meta-analytical method" if a corresponding
    meta-analysis object is printed

* forest.meta(), funnel.default(), funnel.meta(), print.meta(),
  print.summary.meta(), update.meta(), catmeth(), xlab():
  - properly handle R objects of class "metainc"

* metaprop():
  - use correct variable names for 'event' and 'n' in list element
    'data' if metaprop() is called without argument 'data'

* metabin():
  - inverse variance method (argument 'sm = "Inverse"') is used
    automatically if argument 'tau.common = TRUE'
  - bug fix for error if argument 'tau.common = TRUE' and 'method =
    "MH"'

* catmeth():
  - print information on continuity correction for objects of class
    "metabin", "metaprop", and "metainc"

* summary.meta():
  - fixed effect and random effects estimates and confidence intervals
    are only (re)calculated for R objects created with **meta**,
    version < 2 if argument 'level.comb' has not been used

* trimfill.meta(), trimfill.default():
  - new list elements 'lower.fixed', 'upper.fixed', 'zval.fixed',
    'pval.fixed', 'lower.random', 'upper.random', 'zval.random',
    'pval.random' added to trimfill() object (bug was introduced in
    **meta**, version 2.0-0)

* New datasets smoking and lungcancer as examples for metainc()


## meta, version 3.0-1 (2013-09-17)

### Major revision

This update has been declared as major revision as the user interface
changed by dropping some arguments:
  - print.meta(), forest.meta(), summary.meta(): 'level',
    'level.prediction'
  - print.meta(), forest.meta(), metainf(), metacum(): 'level.comb'
  - in forest.meta(), summary.meta(): 'byvar'

This functionality is now provided by update.meta().

### Details

* New function update.meta() to update an existing meta-analysis
  object created with metabin(), metacont(), metagen(), metaprop(), or
  metacor()

* New function cilayout() to change layout of confidence intervals

* Deprecated function plot.meta() removed

* metabin(), metacont(), metagen(), metaprop(), metacor():
  - code cleaning for new R function update.meta()

* metabin(), metacont(), metagen(), metaprop(), metacor(),
  summary.meta():
  - new list components:
	- 'data' with original data used in function call
	- 'subset' with information on subset used in meta-analysis

* metareg():
  - argument 'data' renamed to 'x'
  - first two arguments interchanged (which is now in line with other
    R functions from R package **meta**)
  - information on grouping variable (list element 'byvar') is
    utilised if argument 'formula' is missing
  - any column from original dataset (list element 'data') can be
    used in meta-regression

* trimfill.meta():
  - new defaults for arguments 'comb.fixed' and 'comb.random' (by
    default only random effects estimate calculated)
  - arguments 'sm' and 'studlab' removed
  - new list elements (depending on class of meta-analysis object
    used in function call):
    - 'event.e', 'event.c', 'event' with number of events
    - 'n.e', 'n.c', 'n' with number of observations
    - 'mean.e', 'mean.c', 'sd.e', 'sd.c' with means and standard
      deviations
    - 'cor' with correlation
    - 'class.x' with class of meta-analysis object used in function
      call

* trimfill.default():
  - only calculate random effects estimate by default

* metacr():
  - new list elements 'event.e', 'n.e', 'event.c' and 'n.c' for Peto
    O-E method

* metaprop():
  - new list element 'incr.event'

* forest.meta():
  - bug fix for error if (i) the effect measure is equal to RR, OR, or
    HR and (ii) argument 'label' is not a logical value

* print.summary.meta():
  - print "0" instead of "< 0.0001" if between-study heterogeneity is
    zero

* New function format.tau() to print "0" instead of "< 0.0001" if tau2
  is zero

* p.ci():
  - new arguments 'bracket.left', 'separator' and 'bracket.right' for
    more flexible confidence interval layouts

* Several help pages updated


## meta, version 2.5-1 (2013-08-09)

* forest.meta():
  - new argument 'just.studlab' to change justification of study
    labels

* print.meta():
  - print correct information on method to calculate approximate
    confidence interval for metaprop() with a single study

* trimfill.meta():
  - new list elements 'title', 'complab', 'outclab', 'label.e',
    'label.c', 'label.left' and 'label.right'


## meta, version 2.5-0 (2013-07-25)

* metacr():
  - new arguments 'prediction' and 'level.predict' for prediction
    interval for a new study
  - new argument 'tau.common' for common tau2 across subgroups
  - new arguments 'level' and 'level.comb' for confidence interval of
    single studies or meta-analysis

* trimfill.meta(), trimfill.default():
  - new arguments 'prediction' and 'level.predict'

* forest.meta():
  - heterogeneity statistics are only shown if results for fixed
    effect or random effects model are plotted

* metagen(), metabin(), metacont(), metaprop(), metacor():
  - list elements 'comb.fixed', 'comb.random', and 'prediction' are
    set to FALSE for a single study

* print.meta(), print.summary.meta():
  - new argument 'logscale' to print results for summary measures
    "RR", "OR", "HR", or "PLN" on logarithmic scale

* Several help pages updated


## meta, version 2.4-1 (2013-06-20)

* metaprop():
  - bug fix for error in forest.meta() (bug was introduced in
    **meta**, version 2.4-0)
  - new list elements 'incr', 'allincr' and 'addincr' added (bug was
    introduced in **meta**, version 1.5-0)

* print.meta():
  - new arguments 'prediction' and 'level.predict' to print prediction
    interval for a new study

* forest.meta(), print.summary.meta():
  - only print warnings in internal call of asin2p() if result for
    fixed effect, random effects model or prediction interval are
    printed

* asin2p():
  - new argument 'warn' to only print warnings for meta-analysis results

* Example to generate forest plot added to help pages of metabin(),
  metacont(), metacor(), metacr(), metagen(), metaprop()


## meta, version 2.4-0 (2013-06-17)

* metagen(), metabin(), metacont(), metaprop(), metacor():
  - new arguments 'prediction' and 'level.predict' (prediction
    interval for a new study)
  - new argument 'tau.common' (common tau2 across subgroups)
  - help pages updated accordingly

* metaprop():
  - new default summary measure (sm = "PLOGIT")
  - deprecated argument 'freeman.tukey' removed

* summary.meta():
  - new arguments 'prediction' and 'level.predict'
  - list element 'tau.common' from meta-analysis object considered
  - correct values for list elements 'incr', 'allincr', and 'addincr'
    used in calculations for metaprop() objects

* forest.meta():
  - new arguments for prediction interval: 'prediction',
    'level.predict', 'text.predict', 'col.predict',
    'col.predict.lines', 'fs.predict', 'fs.predict.labels',
    'ff.predict', 'ff.predict.labels"
  - correct values for list elements 'incr', 'allincr', and 'addincr'
    used in calculations for metaprop() objects
  - information on confidence limit printed for pooled estimates if CI
    level is different from CI level for individual studies

* print.summary.meta():
  - new argument 'prediction'
  - new list element 'tau.common'

* catmeth():
  - print information on use of common tau2 across subgroups


## meta, version 2.3-0 (2013-05-12)

* forest.meta():
  - results for fixed effect and random effects models only
    (re)calculated for meta-analysis objects created with **meta**,
    version < 2

* metabin():
  - bug fix for error if argument 'sm = "RR"' and 'allstudies = TRUE'
    in meta-analysis with zero events in both groups


## meta, version 2.2-1 (2013-03-20)

* forest.meta():
  - new argument 'lab.NA' to label missing values (in older version of
    R package **meta** the fixed label 'NA' was used)
  - arguments 'colgap.forest.left' and 'colgap.forest.right'
    considered instead of only 'colgap.forest'

* labbe.metabin(), labbe.default():
  - bug fix for error if any event probability is equal to NA

* format.p():
  - bug fix for error if first argument contains any NAs


## meta, version 2.2-0 (2013-03-12)

* metabin():
  - studies with all events in both groups will be included in
    meta-analysis by default (in older **meta** versions these studies
    were only included if argument 'allstudies = TRUE')
  - standard error is positive for studies with all events in both
    groups (see Hartung & Knapp (2001), Stat Med, equation (18))

* forest.meta():
  - values provided by argument 'xlim' will be used as x-axis label
    for relative effect measures like risk ratio or odds ratio
  - default values for arguments 'smlab.pos' and 'xlab.pos' changed to
    always fall within plotting range


## meta, version 2.1-4 (2012-11-29)

* forest.meta(), metacum(), metainf(), print.meta():
  - correct back transformation of Freeman-Tukey Double arcsine
    transformation for metacum() and metainf() objects

* asin2p():
  - values outside the admissible range are set to the boundary values
    [0, 1]; a warning is printed in this casea

* Help pages:
  - new argument 'n.harmonic.mean' documented for metacum() and
    metainf()


## meta, version 2.1-3 (2012-11-20)

* forest.meta():
  - bug fix for metacum() or metainf() object with Freeman-Tukey
    double arcsine transformation (error message: 'Error in if
    (col$range[1] <= TE.fixed & TE.fixed <= col$range[2]) ...')


## meta, version 2.1-2 (2012-10-25)

* forest.meta(), metacum(), metainf(), summary.meta():
  - bug fix for error if argument 'method.tau' is different from the
    default (error message: 'Error in sqrt(W) %*% X : non-conformable
    arguments')

* forest.meta():
  - argument 'byvar' uses corresponding list element from
    meta-analysis object as default


## meta, version 2.1-1 (2012-08-12)

* summary.meta():
  - list element 'k0' added to trim-and-fill object

* print.summary.meta():
  - print number of added studies for trim-and-fill method


## meta, version 2.1-0 (2012-05-18)

* trimfill.meta(), trimfill.default():
  - new arguments 'hakn' and 'method.tau'

* metacum(), metainf():
  - add class "trimfill" for trim-and-fill objects

* catmeth(), print.meta(), print.summary.meta():
  - print information on trim-and-fill method

* metabias.meta(), funnel.meta():
  - print error message if used with metacum() or metainf() object

* funnel.meta():
  - use different plot symbols (argument 'pch') for trimfill() object

* .onLoad():
  - version nummer of **meta** is printed when library is loaded

* Help pages:
  - arguments 'hakn' and 'method.tau' documented in trimfill.meta()
    and trimfill.default()
  - changed default for argument 'pch' in funnel.meta() documented


## meta, version 2.0-4 (2012-05-03)

* metaprop():
  - variance estimate for log- and logit-transformation is based on
    Pettigrew et al. (1986)

* Reference Pettigrew et al. (1986) added to help pages for metabin()
  and metaprop()


## meta, version 2.0-2 (2012-04-17)

* metaprop():
  - warning is printed if any sample size is smaller than 10 for
    Freeman-Tukey double arcsine transformation


## meta, version 2.0-1 (2012-04-04)

* metabin(), metacont(), metacor(), metagen(), metaprop():
  - arguments 'subset' and 'byvar' can be of different length


## meta, version 2.0-0 (2012-03-20)

### Major revision

R package **meta** linked to R package **metafor** to provide
additional statistical methods, e.g. meta-regression and other
estimates for tau2 (REML, ...)

### Details

* New functions:
  - metareg()              meta-regression
  - metabias()             generic method for metabias()
  - metabias.default()     generic method for metabias()
  - metabias.meta()        generic method for metabias()
  - metabias.rm5()         generic method for metabias()
  - print.rm5()            generic method for rm5-object
  - print.summary.rm5()    generic method for rm5-object
  - summary.rm5()          generic method for rm5-object
  - catmeth()              internal function
  - crtitle()              internal function
  - hypergeometric()       internal function
  - is.installed.metafor() internal function
  - kentau()               internal function
  - linregcore()           internal function
  - p2logit()              internal function

* metabin(), metacont(), metagen():
  - new arguments 'label.left' and 'label.right' to add label on left
    or right side of forest plot

* metabin(), metacont(), metacor(), metagen(), metaprop():
  - new arguments:
	- 'hakn' (Hartung-Knapp method)
	- 'method.tau' (estimation method for tau2)
	- 'tau.preset' (fixed value for tau)
	- 'TE.tau' (pre-specified treatment effect to estimate tau)
	- 'method.bias' (test for funnel plot asymmetry used in metabias)
	- 'warn' (print warning messages)
  - new list elements in meta-analysis object:
	- 'se.tau2' with standard error of tau2
	- 'hakn' for Hartung-Knapp method
	- 'method.tau' with information on estimation method for tau2
	- 'tau.preset' for fixed tau value
	- 'TE.tau' for pre-specified treatment effect to estimate tau
	- 'method.bias' for test of funnel plot asymmetry used in
      metabias()
  - argument 'warn = FALSE' suppresses additional warning messages

* metabin():
  - studies are excluded from meta-analysis if (event.e > n.e |
    event.c > n.c) or (n.e <= 0 | n.c <= 0) or (event.e < 0 | event.c
    < 0)

* metacum(), metainf():
  - return NULL if function is used with a single study
  - arguments 'hakn', 'method.tau', 'tau.preset', 'method.bias',
    'label.left', 'label.right' are considered from meta-analysis
    object
  - argument 'level' removed

* metaprop():
  - correct variance 1 / (n + 0.5) instead of 1 / (n + 1) is used for
    the Freeman-Tukey double arcsine transformation (argument 'sm =
    "PFT"')

* asin2p():
  - completely rewritten as back transformation of Freeman-Tukey
    transformed proportions was inaccurate
  - back transformation of Freeman-Tukey proportions according to
    Miller (1978) - see help page of metaprop()

* print.metabias():
  - print a warning if number of studies is too small to conduct a
    test for funnel plot asymmetry

* print.summary.meta():
  - new argument 'bylab.nchar'
  - print test for subgroup differences for both fixed effect and
    random effects model
  - invisible(NULL) returned for metacum() and metainf() objects

* metacr():
  - new arguments:
	- 'sm' (summary measure)
	- 'method' (pooling method)
	- 'comb.fixed' (fixed effect model)
	- 'comb.random' (random effects model)
	- 'swap.events' (only for binary data)
	- 'method.tau' (estimation method for between-study variance)
	- 'hakn' (Hartung-Knapp adjustment)
	- 'title' (Title of Cochrane review)
	- 'complab' (Comparison label)
	- 'outclab' (Outcome label)
	- 'warn' (print warning messages)
  - removed argument:
	- 'smother' (summary measure)
  - return NULL if no data is available for selection of arguments
      'comp.no' and 'outcome.no'

* read.rm5():
  - changed substantially for reading of RevMan 5.1 files

* summary.meta():
  - arguments 'hakn', 'method.tau', 'tau.preset', 'method.bias' are
    considered from meta-analysis object
  - argument 'warn = FALSE' suppresses additional warning messages

* forest.meta():
  - treatment effect and 95% confidence interval is printed in the
    correct order for objects of class "metaprop" if argument 'sort'
    or 'order' is used
  - symmetric forest plot by default (argument xlim = "s")
  - new arguments:
	- 'smlab', 'smlab.pos', 'fs.smlab', 'fflab' for label of summary
      measure at top of figure
	- 'label.right', 'label.left', 'fs.lr', 'ff.lr' for label on right
      and left side below the x-axis
  - 'overall.hetstat' to show heterogeneity information for overall
    effect

* funnel.default(), funnel.meta():
  - arguments 'col.fixed' and 'col.random' are recognised

* metabias.default(), metabias.meta():
  - new argument 'k.min' to only conduct test for funnel plot
    asymmetry if number of studies in meta-analysis is larger or equal
    to 'k.min'
  - new argument '...' (ignored at the moment)

* trimfill.default(), trimfill.meta():
  - return 'invisible(NULL)' if number of studies is smaller than 3

* New datasets: amlodipine, cisapride

* File FLEISS93.MTV moved from directory data to directory extdata

* Several help pages updated

* Some new help pages added


## meta, version 1.6-1 (2010-10-28)

* forest.meta():
  - number of events is printed in the correct order for objects of
    class "metaprop" if argument 'sort' or 'order' is used
  - transformed proportions are printed for individual studies in
    column 'TE' if metagen() is used with argument 'sm' equal to
    either "PLN", "PLOGIT", "PAS", or "PFT"

* as.data.frame.meta():
  - function works for meta-analyses with six studies which previously
    resulted in an error message 'Error: evaluation nested too deeply:
    infinite recursion ...'
  - new argument '...'

* addvar():
  - option stringsAsFactors = FALSE added to internal call of
    as.data.frame.meta()
  - additional checks for existence of columns 'by.x' and 'by.y'
  - additional checks for situations with duplicate entries for
    columns 'by.x' and 'by.y' added

* print.meta():
  - back transformed proportions are printed for individual studies if
    metagen() is used with argument 'sm' equal to either "PLN",
    "PLOGIT", "PAS", or "PFT"

* Examples in help pages (slightly) updated:
  - read.mtv(), read.rm5(), metacr()


## meta, version 1.6-0 (2010-06-21)

* forest.meta():
  - for subgroup analyses (i.e. groups defined by argument 'byvar'),
    result for both fixed effect and random effects model are printed
    (in older versions of the **meta** package only results for either
    fixed effect or random effects model could be printed)
  - new arguments 'text.fixed.w' and 'text.random.w' to specify label
    for estimates within subgroups
  - new arguments to change colour of several parts of the plot:
    'col.i.inside.square', 'col.square', 'col.square.lines',
    'col.diamond', 'col.diamond.fixed', 'col.diamond.random',
    'col.diamond.lines', 'col.diamond.fixed.lines',
    'col.diamond.random.lines'
  - new arguments to print information on heterogeneity measures:
    'print.I2', 'print.tau2', 'print.Q', 'print.pval.Q', 'hetstat',
    'hetlab'
  - new arguments to change fontsize and fontface of several parts of the plot:
    'fs.heading', 'fs.fixed', 'fs.random', 'fs.study',
    'fs.fixed.labels', 'fs.random.labels', 'fs.study.labels',
    'fs.hetstat', 'fs.axis', 'fs.xlab', 'ff.heading', 'ff.fixed',
    'ff.random', 'ff.study', 'ff.fixed.labels', 'ff.random.labels',
    'ff.study.labels', 'ff.hetstat', 'ff.axis', 'ff.xlab'
  - new arguments to change gap between columns:
    'colgap.left', 'colgap.right=colgap', 'colgap.forest',
    'colgap.forest.left', 'colgap.forest.right'
  - new argument 'just' to change justification of text for additional
    columns
  - new argument 'addspace' to print a blank line at top and bottom of
    study results
  - argument 'squaresize' supersedes argument 'boxsize'
  - new argument 'new' indicating whether a new figure should be
    printed in an existing graphics window (internally, grid.newpage()
    is used if argument 'new = TRUE')
  - no line is printed for the fixed effect or random effects model if
    argument 'lty.fixed = NULL' or 'lty.random = NULL'
  - symmetric forest plots can be produced by setting argument 'xlim =
    "s"'

* print.summary.meta():
  - for subgroup analyses (i.e. groups defined by argument 'byvar'),
    result for test of heterogeneity printed separately for fixed
    effect and random effects model

* metabin(), summary.meta(), print.summary.meta():
  - new argument 'print.CMH' indicating whether
    Cochran-Mantel-Haenszel test for overall effect should be printed
    (default 'print.CMH = FALSE')

* Help pages updated:
  - forest.meta(), metabin(), print.summary.meta(), summary.meta()


## meta, version 1.5-0 (2010-05-07)

* Version jump to 1.5-0 as several changes have been implemented

* New functions:
  - metacor()        meta-analysis of correlations
  - forest()         generic method for forest plots
  - forest.meta()    generic method for forest plots
  - radial()         generic method for radial plots
  - radial.default() generic method for radial plots
  - radial.meta()    generic method for radial plots
  - asin2p()         internal function
  - logit2p()        internal function
  - xlab()           internal function
  - z2cor()          internal function

* forest.meta():
  - new arguments 'pooled.totals' and 'pooled.events' to specify
    whether total number of observations and events should be
    displayed in the plot
  - new argument 'pscale' to rescale proportions for objects of class
    "metaprop"
  - arguments 'label' and 'xlim' are recognised for other effect
    measures than RR, OR, and HR
  - arguments 'rightlabs' and 'leftlabs' accept NAs for columns using
    default labels
  - significant digits are printed uniformly
  - correct sum of percentage weight is printed for random effects
    model in forest plots with subgroups
  - x limits (min,max) of the plot are defined by the width of
    confidence intervals instead of (0,1) for objects of class "metaprop"

* metaprop():
  - implementation of additional transformations: log transformation,
    logit transformation, raw, i.e. untransformed, proportions
  - new argument 'sm' to choose summary measure (i.e. transformation)
  - use of argument 'freeman.tukey' is deprecated (replaced by
    argument 'sm')

* funnel(), funnel.meta(), labbe(), labbe.meta():
  - argument 'y' removed

* trimfill(), trimfill.meta():
  - argument 'seTE' removed

* summary.meta():
  - new list elements 'H.w', 'I2.w', 'Q.b.fixed' and 'Q.b.random' for
    heterogeneity statistics within subgroups

* forest.meta(), metacum(), metainf(), print.meta(),
  print.summary.meta(), summary.meta():
  - extension for meta-analysis of correlations

* plot.meta():
  - print warning that function was replaced by forest.meta()

* New list element 'version' with information on version number of
  **meta** package used to create an object; applies only to object
  creating functions, e.g. metabin() and metabias()

* Several help pages updated

* Use file NEWS instead of ChangeLog to document changes


## meta, version 1.1-8 (2010-01-12)

* summary.meta(), print.summary.meta():
  - test for subgroup differences is not calculated and printed for
    meta-analyses using the Mantel-Haenszel method for binary data


## meta, version 1.1-7 (2010-01-11)

* metabin(), metacont(), metagen(), metaprop():
  - sensible default value is used for argument 'bylab' if argument
   'byvar' is not missing


## meta, version 1.1-6 (2010-01-11)

* forest():
  - additional columns are printed in the correct order if argument
   'sort' or 'order' is used


## meta, version 1.1-5 (2009-12-21)

* forest():
  - new argument 'digits' specifying minimal number of significant
    digits for treatment estimate and its confidence interval


## meta, version 1.1-4 (2009-11-04)

* summary.meta():
  - results for subgroups (if byvar != NULL) are calculated for both
    fixed effect and random effects model:
   * list 'within' no longer returned by summary.meta()
   * lists 'within.fixed' and 'within.random' returned by
     summary.meta()
  - variable name of subgroups is printed correctly
  - check whether input is an object of class "summary.meta"

* print.summary.meta():
  - a warning is printed if both 'comb.fixed' and 'comb.random' are
    TRUE and results for subgroups are supposed to be printed

* Help pages of print.summary.meta() and forest() updated:
  - detailed information on printing and plotting of subgroup
   results if both comb.fixed and comb.random are TRUE

* Help page of metagen() updated:
  - new example with meta-analysis of survival data


## meta, version 1.1-3 (2009-10-30)

* Generic method for trim-and-fill method: trimfill(),
  trimfill.default(), trimfill.meta()


## meta, version 1.1-2 (2009-10-09)

* L'Abbe plot implemented: labbe(), labbe.default(), labbe.metabin()

* Generic method for funnel plots: funnel(), funnel.default(),
  funnel.meta()

* funnel.meta(), funnel.default():
  - contour-enhanced funnel plots can be produced (new arguments
    'contour.levels', 'col.contour', 'ref')
  - study labels can be printed on funnel plot (new arguments
    'studlab', 'cex.studlab')
  - line type, width and colour can be changed for fixed effect
    treatment effect (new arguments 'lty.fixed', 'lwd.fixed',
    'col.fixed')
  - random effects treatment effect can be plotted (new arguments
    'comb.random', 'lty.random', 'lwd.random', 'col.random')
  - new default values for some arguments:
	* 'pch = 21' (previously: 'pch = 1')
	* 'comb.fixed = x$comb.fixed'
  - background colour of points in funnel plot can be changed (new
    argument 'bg')

* forest():
  - new default values for arguments 'lab.e' and 'lab.c':
	* x$label.e and x$label.c (if these values are NULL the old
      default values "Experimental" and "Control" are used)

* metabin(), metacont(), metagen():
  - arguments 'label.e' and 'label.c' added

* metacr():
  - use arguments 'label.e' and 'label.c' in calls to metabin(),
    metacont(), metagen()


## meta, version 1.0-6 (2009-08-31)

* First major version released on CRAN (no ChangeLog available)

* Older versions (without ChangeLog):

  - meta, version 0.9-19 (2009-04-09)

  - meta, version 0.9-18 (2009-03-03)

  - meta, version 0.9-17 (2009-01-28)

  - meta, version 0.9-15 (2008-12-08)

  - meta, version 0.8-2 (2007-02-13)

  - meta, version 0.81 (2006-11-29)

  - meta, version 0.8 (2006-11-24)

  - First CRAN release (2006-02-08)

  - meta, version 0.5-0 (2005-02-23)
