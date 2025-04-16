#' Forest plot to display the result of a meta-analysis
#' 
#' @description
#' Draws a forest plot in the active graphics window (using grid
#' graphics system).
#' 
#' @aliases forest.metabind
#' 
#' @param x An object of class \code{\link{metabind}}.
#' @param leftcols A character vector specifying (additional) columns
#'   to be plotted on the left side of the forest plot or a logical
#'   value (see Details).
#' @param leftlabs A character vector specifying labels for
#'   (additional) columns on left side of the forest plot (see
#'   Details).
#' @param rightcols A character vector specifying (additional) columns
#'   to be plotted on the right side of the forest plot or a logical
#'   value (see Details).
#' @param rightlabs A character vector specifying labels for
#'   (additional) columns on right side of the forest plot (see
#'   Details).
#' @param common A logical indicating whether common effect estimates
#'   should be plotted.
#' @param random A logical indicating whether random effects estimates
#'   should be plotted.
#' @param overall A logical indicating whether overall summaries
#'   should be plotted. This argument is useful in a meta-analysis
#'   with subgroups if summaries should only be plotted on group
#'   level.
#' @param subgroup A logical indicating whether subgroup results
#'   should be shown in forest plot. This argument is useful in a
#'   meta-analysis with subgroups if summaries should not be plotted
#'   on group level.
#' @param hetstat Either a logical value indicating whether to print
#'   results for heterogeneity measures at all or a character string
#'   (see Details).
#' @param overall.hetstat A logical value indicating whether to print
#'   heterogeneity measures for overall treatment comparisons. This
#'   argument is useful in a meta-analysis with subgroups if
#'   heterogeneity statistics should only be printed on subgroup
#'   level.
#' @param prediction A logical indicating whether prediction
#'   interval(s) should be printed.
#' @param type A character string or vector specifying how to
#'   plot estimates.
#' @param type.common A single character string specifying how to
#'   plot common effect estimates.
#' @param type.random A single character string specifying how to
#'   plot random effects estimates.
#' @param type.predict A single character string specifying how to
#'   plot prediction intervals.
#' @param lab.NA A character string to label missing values.
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
#' @param col.diamond.common The colour of diamonds for common effect
#'   estimates.
#' @param col.diamond.random The colour of diamonds for random effects
#'   estimates.
#' @param col.diamond.lines The colour of the outer lines of diamonds
#'   representing the results for common effect and random effects
#'   models.
#' @param col.diamond.lines.common The colour of the outer lines of
#'   diamond for common effect estimates.
#' @param col.diamond.lines.random The colour of the outer lines of
#'   diamond for random effects estimates.
#' @param col.predict Background colour of prediction intervals.
#' @param col.predict.lines Colour of outer lines of prediction
#'   intervals.
#' @param digits Minimal number of significant digits for treatment
#'   effects, see \code{print.default}.
#' @param digits.se Minimal number of significant digits for standard
#'   errors, see \code{print.default}.
#' @param digits.stat Minimal number of significant digits for z- or
#'   t-statistic for test of overall effect, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall treatment effect, see \code{print.default}.
#' @param digits.pval.Q Minimal number of significant digits for
#'   p-value of heterogeneity test, see \code{print.default}.
#' @param digits.Q Minimal number of significant digits for
#'   heterogeneity statistic Q, see \code{print.default}.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance, see \code{print.default}.
#' @param digits.tau Minimal number of significant digits for square
#'   root of between-study variance, see \code{print.default}.
#' @param digits.I2 Minimal number of significant digits for I-squared
#'   statistic, see \code{print.default}.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param big.mark A character used as thousands separator.
#' @param print.subgroup.labels A logical indicating whether subgroup
#'   label should be printed.
#' @param addrow.subgroups A logical value indicating whether an empty
#'   row is printed between results for subgroups.
#' @param smlab A label for the summary measurex (printed at top of
#'   figure).
#' @param calcwidth.pooled A logical indicating whether text for
#'   common effect and random effects model should be considered to
#'   calculate width of the column with study labels.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional graphical arguments (passed on to
#'   \code{\link{forest.meta}}).
#' 
#' @details
#' A forest plot, also called confidence interval plot, is drawn in
#' the active graphics window. The forest functions in R package
#' \bold{meta} are based on the grid graphics system. In order to
#' print the forest plot, resize the graphics window and either use
#' \code{\link{dev.copy2eps}} or \code{\link{dev.copy2pdf}}. Another
#' possibility is to create a file using \code{\link{pdf}},
#' \code{\link{png}}, or \code{\link{svg}} and to specify the width
#' and height of the graphic (see \code{\link{forest.meta}} examples).
#' 
#' The arguments \code{leftcols} and \code{rightcols} can be used to
#' specify columns which are plotted on the left and right side of the
#' forest plot, respectively.
#' 
#' The arguments \code{leftlabs} and \code{rightlabs} can be used to
#' specify column headings which are plotted on left and right side of
#' the forest plot, respectively. For certain columns predefined
#' labels exist. For other columns, the column name will be used as a
#' label. It is possible to only provide labels for new columns (see
#' \code{\link{forest.meta}} examples). Otherwise the length of
#' \code{leftlabs} and \code{rightlabs} must be the same as the number
#' of printed columns, respectively. The value \code{NA} can be used
#' to specify columns which should use default labels.
#' 
#' Argument \code{hetstat} can be a character string to specify where
#' to print heterogeneity information:
#' \itemize{
#' \item row with results for common effect model (\code{hetstat =
#' "common"}),
#' \item row with results for random effects model (\code{hetstat =
#' "random"}),
#' \item rows with 'study' information (\code{hetstat = "study"}).
#' }
#' Otherwise, information on heterogeneity is printed in dedicated rows.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{forest.meta}}, \code{\link{metabin}},
#'   \code{\link{metacont}}, \code{\link{metagen}},
#'   \code{\link{metabind}}, \code{\link{settings.meta}}
#' 
#' @keywords hplot
#' 
#' @examples
#' data(Fleiss1993cont)
#' 
#' # Add some (fictitious) grouping variables:
#' #
#' Fleiss1993cont$age <- c(55, 65, 55, 65, 55)
#' Fleiss1993cont$region <- c("Europe", "Europe", "Asia", "Asia", "Europe")
#' 
#' m1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, sm = "SMD")
#'
#' # Conduct two subgroup analyses
#' #
#' mu1 <- update(m1, subgroup = age, subgroup.name = "Age group")
#' mu2 <- update(m1, subgroup = region, subgroup.name = "Region")
#'
#' # Combine subgroup meta-analyses and show forest plot with subgroup
#' # results
#' #
#' mb1 <- metabind(mu1, mu2)
#' mb1
#' forest(mb1)
#'
#' @method forest metabind
#' @export

forest.metabind <- function(x,
                            leftcols, leftlabs,
                            rightcols = c("effect", "ci"), rightlabs,
                            ##
                            common = x$common,
                            random = x$random,
                            overall = x$overall,
                            subgroup = FALSE,
                            hetstat = FALSE,
                            overall.hetstat = x$overall.hetstat,
                            prediction = x$prediction,
                            ##
                            lab.NA = "",
                            ##
                            col.square = gs("col.square"),
                            col.square.lines = col.square,
                            col.circle = gs("col.circle"),
                            col.circle.lines = col.circle,
                            ##
                            col.diamond = gs("col.diamond"),
                            col.diamond.common = col.diamond,
                            col.diamond.random = col.diamond,
                            col.diamond.lines = gs("col.diamond.lines"),
                            col.diamond.lines.common = col.diamond.lines,
                            col.diamond.lines.random = col.diamond.lines,
                            ##
                            col.predict = gs("col.predict"),
                            col.predict.lines = gs("col.predict.lines"),
                            #
                            type = NULL,
                            type.common = NULL,
                            type.random = NULL,
                            type.predict = NULL,
                            ##
                            digits = gs("digits.forest"),
                            digits.se = gs("digits.se"),
                            digits.stat = gs("digits.stat"),
                            digits.pval = max(gs("digits.pval") - 2, 2),
                            digits.pval.Q = max(gs("digits.pval.Q") - 2, 2),
                            digits.Q = gs("digits.Q"),
                            digits.tau2 = gs("digits.tau2"),
                            digits.tau = gs("digits.tau"),
                            digits.I2 = max(gs("digits.I2") - 1, 0),
                            ##
                            scientific.pval = gs("scientific.pval"),
                            big.mark = gs("big.mark"),
                            ##
                            print.subgroup.labels = x$with.subgroups,
                            addrow.subgroups = print.subgroup.labels,
                            ##
                            smlab,
                            calcwidth.pooled = overall,
                            ##
                            warn.deprecated = gs("warn.deprecated"),
                            ##
                            ...) {
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  chkclass(x, "metabind")
  x <- updateversion(x)
  ##
  missing.overall <- missing(overall)
  missing.overall.hetstat <- missing(overall.hetstat)
  ##
  overall <- replaceNULL(overall, FALSE)
  overall.hetstat <- replaceNULL(overall.hetstat, FALSE)
  ##
  chklogical(common)
  chklogical(random)
  chklogical(overall)
  chklogical(subgroup)
  chklogical(overall.hetstat)
  ##
  chklogical(print.subgroup.labels)
  chklogical(addrow.subgroups)
  ##
  chkchar(lab.NA)
  #
  if (missing(type))
    type <- x$data$type
  else {
    if (length(type) == 1 & nrow(x$data) > 1)
      type <- rep(type, nrow(x$data))
    else if (length(type) != nrow(x$data))
      stop("Argument 'type' must be a character vector of length 1 or ",
           nrow(x$data), ".",
           call. = FALSE)
  }
  #
  if (!is.null(type.common)) {
    if (length(type.common) != 1)
      stop("Argument 'type.common' must be of length 1.",
           call. = FALSE)
    ##
    type.common <-
      setchar(type.common, c("square", "diamond", "predict", "circle"))
  }
  ##
  if (!is.null(type.random)) {
    if (length(type.random) != 1)
      stop("Argument 'type.random' must be of length 1.",
           call. = FALSE)
    ##
    type.random <-
      setchar(type.random, c("square", "diamond", "predict", "circle"))
  }
  ##
  if (!is.null(type.predict)) {
    if (length(type.predict) != 1)
      stop("Argument 'type.predict' must be of length 1.",
           call. = FALSE)
    ##
    type.predict <-
      setchar(type.predict, c("square", "diamond", "predict", "circle"))
  }
  ##
  chkcolor(col.square, length = 1)
  chkcolor(col.square.lines, length = 1)
  chkcolor(col.circle, length = 1)
  chkcolor(col.circle.lines, length = 1)
  ##
  chkcolor(col.diamond, length = 1)
  chkcolor(col.diamond.common, length = 1)
  chkcolor(col.diamond.random, length = 1)
  chkcolor(col.diamond.lines, length = 1)
  chkcolor(col.diamond.lines.common, length = 1)
  chkcolor(col.diamond.lines.random, length = 1)
  ##
  chkcolor(col.predict, length = 1)
  chkcolor(col.predict.lines, length = 1)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  chknumeric(digits.pval.Q, min = 1, length = 1)
  chknumeric(digits.Q, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  ##
  chklogical(scientific.pval)
  chklogical(calcwidth.pooled)
  ##
  addargs <- names(list(...))
  ##
  idx <- charmatch(tolower(addargs), "hetstat", nomatch = NA)
  if (any(!is.na(idx)) && length(idx) > 0)
    if (list(...)[["hetstat"]])
      stop("Argument 'hetstat' must be FALSE for metabind objects.",
           call. = TRUE)
  ##
  for (i in addargs) {
    if (!is.null(setchar(i, "type.study", stop.at.error = FALSE)))
      stop("Use argument 'type' as argument 'type.study' is set internally.",
           call. = TRUE)
    if (!is.null(setchar(i, "fixed", stop.at.error = FALSE)))
      stop("Argument 'fixed' cannot be used with metabind objects.",
           call. = TRUE)
  }
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  digits.stat <-
    deprecated(digits.stat, missing(digits.stat), args, "digits.zval",
               warn.deprecated)
  digits.stat <- replaceNULL(digits.stat, gs("digits.stat"))
  chknumeric(digits.stat, min = 0, length = 1)
  
  
  x$k.w.orig <- x$k.w
  
  x$k.w <- x$k.study.w <- x$k.all.w <- x$k.TE.w <-
    as.vector(table(x$data$subgroup)[unique(x$data$subgroup)])
  ##as.vector(table(x$data$name)[unique(x$data$name)])
  
  
  missing.leftcols <- missing(leftcols)
  ##
  if (missing.leftcols) {
    leftcols <- c("studlab", "k")
    if (x$with.subgroups)
      leftcols <- c(leftcols, "pval.Q.b")
    else {
      leftcols <- c(leftcols, "tau2")
    }
  }
  ##
  if (!is.logical(rightcols)) {
    if ("pval.Q.b" %in% rightcols & "pval.Q.b" %in% leftcols)
      leftcols <- leftcols[leftcols != "pval.Q.b"]
    ##
    if ("k" %in% rightcols & "k" %in% leftcols)
      leftcols <- leftcols[leftcols != "k"]
    ##
    if ("tau2" %in% rightcols & "tau2" %in% leftcols)
      leftcols <- leftcols[leftcols != "tau2"]
    ##
    if ("Q" %in% rightcols & "Q" %in% leftcols)
      leftcols <- leftcols[leftcols != "Q"]
    ##
    if ("pval.Q" %in% rightcols & "pval.Q" %in% leftcols)
      leftcols <- leftcols[leftcols != "pval.Q"]
  }
  ##
  label.studlab <- ifelse(x$with.subgroups, "Subgroup", "Method")
  ##
  newlab <- function(col, lab, varname, value) {
    sel <- col == varname
    if (!any(sel) || !is.na(lab[sel]))
      return(lab)
    else
      lab[sel] <- value
    lab
  }
  ##
  if (missing(leftlabs))
    leftlabs <- rep(NA, length(leftcols))
  if (any(is.na(leftlabs))) {
    leftlabs <-
      newlab(leftcols, leftlabs, "studlab", label.studlab)
    leftlabs <-
      newlab(leftcols, leftlabs, "k", "Number of\nStudies")
    leftlabs <-
      newlab(leftcols, leftlabs, "tau2", "Between-study\nvariance")
    leftlabs <-
      newlab(leftcols, leftlabs, "tau", "Between-study\nSD")
    leftlabs <-
      newlab(leftcols, leftlabs, "pval.Q.b", "Interaction\nP-value")
    leftlabs <-
      newlab(leftcols, leftlabs, "pval.Q", "P-value")
  }
  ##
  if (missing(rightlabs))
    rightlabs <- rep(NA, length(rightcols))
  if (any(is.na(rightlabs))) {
    rightlabs <-
      newlab(rightcols, rightlabs, "studlab", label.studlab)
    rightlabs <-
      newlab(rightcols, rightlabs, "k", "Number of\nStudies")
    rightlabs <-
      newlab(rightcols, rightlabs, "tau2", "Between-study\nvariance")
    rightlabs <-
      newlab(rightcols, rightlabs, "tau", "Between-study\nSD")
    rightlabs <-
      newlab(rightcols, rightlabs, "pval.Q.b", "Interaction\nP-value")
    rightlabs <-
      newlab(rightcols, rightlabs, "pval.Q", "P-value")
  }
  
  
  ##
  ## Round and round ...
  ##
  if (x$with.subgroups) {
    x$data$Q.b <- round(x$data$Q.b, digits.Q)
    x$data$pval.Q.b <- formatPT(x$data$pval.Q.b,
                                lab = FALSE,
                                digits = digits.pval.Q,
                                scientific = scientific.pval,
                                lab.NA = lab.NA)
    x$data$pval.Q.b <- rmSpace(x$data$pval.Q.b)
  }
  else {
    x$data$Q <- round(x$data$Q, digits.Q)
    x$Q <- round(x$Q, digits.Q)
    x$data$pval.Q <-
      rmSpace(formatPT(x$data$pval.Q,
                       lab = FALSE,
                       digits = digits.pval.Q,
                       scientific = scientific.pval,
                       lab.NA = lab.NA))
  }
  ##
  x$data$tau2 <-
    rmSpace(formatPT(x$data$tau2, digits = digits.tau2,
                     big.mark = big.mark,
                     lab = FALSE, lab.NA = lab.NA))
  ##
  x$data$tau <-
    rmSpace(formatPT(x$data$tau, digits = digits.tau,
                     big.mark = big.mark,
                     lab = FALSE, lab.NA = lab.NA))
  ##
  I2.na <- is.na(x$data$I2)
  x$data$I2 <-
    formatN(round(100 * x$data$I2, digits.I2), digits.I2, "")
  x$data$lower.I2 <-
    formatN(round(100 * x$data$lower.I2, digits.I2), digits.I2, "")
  x$data$upper.I2 <-
    formatN(round(100 * x$data$upper.I2, digits.I2), digits.I2, "")
  ##
  x$data$I2 <-
    rmSpace(paste0(x$data$I2, ifelse(I2.na, "", "%")))
  x$data$lower.I2 <-
    rmSpace(paste0(x$data$lower.I2, ifelse(I2.na, "", "%")))
  x$data$upper.I2 <-
    rmSpace(paste0(x$data$upper.I2, ifelse(I2.na, "", "%")))
  ##
  x$data$k <- as.character(x$data$k)
  #
  model <- x$data$model
  ##
  if (!is.null(type.common))
    type[model == "common"] <- type.common
  ##
  if (!is.null(type.random))
    type[model == "random"] <- type.random
  ##
  if (!is.null(type.predict))
    type[model == "predict"] <- type.predict
  ##
  cols.square <- vector("character", length(type))
  cols.square.lines <- vector("character", length(type))
  ##
  for (i in seq_along(type)) {
    if (type[i] == "diamond") {
      if (model[i] == "common") {
        cols.square[i] <- col.diamond.common
        cols.square.lines[i] <- col.diamond.lines.common
      }
      else if (model[i] == "random") {
        cols.square[i] <- col.diamond.random
        cols.square.lines[i] <- col.diamond.lines.random
      }
      else {
        col.predict[i] <- col.diamond.random
        cols.square.lines[i] <- col.diamond.lines.random
      }
    }
    else if (type[i] == "predict") {
      cols.square[i] <- col.predict
      cols.square.lines[i] <- col.predict.lines
    }
    else if (type[i] == "circle") {
      cols.square[i] <- col.circle
      cols.square.lines[i] <- col.circle.lines
    }
    else {
      cols.square[i] <- gs("col.square")
      cols.square.lines[i] <- col.square.lines
    }
  }
  ##
  if (missing(smlab))
    smlab <- xlab_meta(x$sm, x$backtransf)
  ##
  if (!x$samedata) {
    overall <- FALSE
    overall.hetstat <- FALSE
  }
  
  
  x.forest <- x
  ##
  if (x$with.subgroups) {
    m.forest <- metagen(x.forest$TE,
                        lower = x.forest$lower,
                        upper = x.forest$upper,
                        subgroup = x.forest$subgroup,
                        studlab = x.forest$studlab,
                        data = x.forest,
                        common = common,
                        random = random,
                        prediction = prediction,
                        overall = overall,
                        print.subgroup.name = FALSE)
    #
    class(m.forest) <- c(class(m.forest), "is.metabind")
    m.forest$sm <- x.forest$sm
    #
    m.forest$TE[m.forest$data$type == "predict"] <- NA
    m.forest$statistic <- x.forest$statistic
    m.forest$pval <- x.forest$pval
    #
    m.forest$n.harmonic.mean <- x.forest$n.harmonic.mean
    m.forest$t.harmonic.mean <- x.forest$t.harmonic.mean
    m.forest$n.harmonic.mean.ma <- x.forest$n.harmonic.mean.ma
    m.forest$t.harmonic.mean.ma <- x.forest$t.harmonic.mean.ma
    m.forest$n.harmonic.mean.w <- x.forest$n.harmonic.mean.w
    m.forest$t.harmonic.mean.w <- x.forest$t.harmonic.mean.w
    #
    if (x$samedata) {
      m.forest$method <- x.forest$method
      m.forest$method.random <- x.forest$method.random
      m.forest$method.predict <- x.forest$method.predict
      m.forest$method.tau <- x.forest$method.tau
      #
      m.forest$TE.common <- x.forest$TE.common
      m.forest$seTE.common <- x.forest$seTE.common
      m.forest$lower.common <- x.forest$lower.common
      m.forest$upper.common <- x.forest$upper.common
      m.forest$statistic.common <- x.forest$statistic.common
      m.forest$pval.common <- x.forest$pval.common
      m.forest$zval.common <- x.forest$zval.common
      m.forest$text.common <- x.forest$text.common
      ##
      m.forest$w.common <-
        ifelse(is.na(m.forest$w.common), m.forest$w.common, NA)
      ##
      m.forest$TE.random <- x.forest$TE.random
      m.forest$seTE.random <- x.forest$seTE.random
      m.forest$lower.random <- x.forest$lower.random
      m.forest$upper.random <- x.forest$upper.random
      m.forest$statistic.random <- x.forest$statistic.random
      m.forest$pval.random <- x.forest$pval.random
      m.forest$zval.random <- x.forest$zval.random
      m.forest$df.random <- x.forest$df.random
      m.forest$text.random <- x.forest$text.random
      ##
      m.forest$w.random <-
        ifelse(is.na(m.forest$w.random), m.forest$w.random, NA)
      ##
      m.forest$lower.predict <- x.forest$lower.predict
      m.forest$upper.predict <- x.forest$upper.predict
      m.forest$text.predict <-
        replaceNULL(x.forest$text.predict, gs("text.predict"))
    }
    ##
    res <-
      forest(m.forest,
             leftcols = leftcols, leftlabs = leftlabs,
             rightcols = rightcols, rightlabs = rightlabs,
             print.subgroup.labels = print.subgroup.labels,
             addrow.subgroups = addrow.subgroups,
             ##
             hetstat = hetstat, overall.hetstat = overall.hetstat,
             calcwidth.pooled = calcwidth.pooled,
             lab.NA = lab.NA, smlab = smlab,
             ##
             type.study = type, col.square = cols.square,
             col.square.lines = cols.square.lines,
             ##
             pooled.totals = FALSE,
             common.subgroup = FALSE, random.subgroup = FALSE,
             prediction.subgroup = FALSE,
             ##
             col.predict = col.predict,
             ##
             weight.study = "same",
             ##
             digits = digits,
             digits.se = digits.se,
             digits.stat = digits.stat,
             digits.pval = digits.pval,
             digits.pval.Q = digits.pval.Q,
             digits.Q = digits.Q,
             digits.tau2 = digits.tau2,
             digits.tau = digits.tau,
             digits.I2 = digits.I2,
             ##
             scientific.pval = scientific.pval,
             big.mark = big.mark,
             ##
             test.subgroup = FALSE,
             details = FALSE,
             ...)
  }
  else {
    if (isCol(x.forest$data, "subgroup"))
      m.forest <- metagen(x.forest$data$TE,
                          lower = x.forest$data$lower,
                          upper = x.forest$data$upper,
                          studlab = x.forest$data$studlab,
                          subgroup = x.forest$data$subgroup,
                          print.subgroup.name = FALSE,
                          data = x.forest,
                          common = FALSE,
                          random = FALSE,
                          prediction = FALSE,
                          overall = FALSE,
                          method.tau = "DL", method.tau.ci = "")
    else
      m.forest <- metagen(x.forest$data$TE,
                          lower = x.forest$data$lower,
                          upper = x.forest$data$upper,
                          studlab = x.forest$data$studlab,
                          data = x.forest,
                          common = common,
                          random = random,
                          prediction = prediction,
                          overall = overall)
    #
    class(m.forest) <- c(class(m.forest), "is.metabind")
    m.forest$sm <- x.forest$sm
    m.forest$with.subgroups <- TRUE
    #
    m.forest$n.harmonic.mean <- x.forest$n.harmonic.mean
    m.forest$t.harmonic.mean <- x.forest$t.harmonic.mean
    #
    m.forest$TE[m.forest$data$type == "predict"] <- NA
    ##
    res <-
      forest(m.forest,
             leftcols = leftcols, leftlabs = leftlabs,
             rightcols = rightcols, rightlabs = rightlabs,
             ##
             hetstat = hetstat, overall.hetstat = overall.hetstat,
             calcwidth.pooled = calcwidth.pooled,
             lab.NA = lab.NA, smlab = smlab,
             ##
             type.study = type, col.square = cols.square,
             col.square.lines = cols.square.lines,
             ##
             pooled.totals = FALSE,
             ##
             col.predict = col.predict,
             ##
             weight.study = "same",
             ##
             digits = digits,
             digits.se = digits.se,
             digits.stat = digits.stat,
             digits.pval = digits.pval,
             digits.pval.Q = digits.pval.Q,
             digits.Q = digits.Q,
             digits.tau2 = digits.tau2,
             digits.tau = digits.tau,
             digits.I2 = digits.I2,
             ##
             scientific.pval = scientific.pval,
             big.mark = big.mark,
             ##
             test.subgroup = FALSE,
             details = FALSE,
             ...)
  }
  ##
  return(invisible(res))
}
