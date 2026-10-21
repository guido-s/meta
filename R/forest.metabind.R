#' Forest plot to display the result of a meta-analysis
#' 
#' @description
#' Draws a forest plot in the active graphics window (using grid
#' graphics system).
#' 
#' @aliases forest.metabind
#' 
#' @param x An object of class \code{\link{metabind}}.
#' @param layout A character string specifying the layout of the
#'   forest plot (see \code{\link{forest.meta}}).
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
#' @param col A named list specifying colours used in the forest
#'   plot; see \code{\link{forest-formatting}}.
#' @param print A named list specifying which information to print in the forest
#'   plot; see \code{\link{forest-formatting}}.
#' @param type A named list specifying graphical types for estimates, common
#'   effect estimates, random effect estimates, and prediction intervals;
#'   see \code{\link{forest-formatting}}.
#' @param digits A named list specifying digits used for displayed results; a
#'   scalar value is also accepted. See \code{\link{forest-formatting}}.
#' @param leftcols A character vector specifying columns to be printed on the
#'   left side of the forest plot or a logical value; see
#'   \code{\link{forest.meta}}.
#' @param rightcols A character vector specifying columns to be printed on the
#'   right side of the forest plot or a logical value; see
#'   \code{\link{forest.meta}}.
#' @param leftlabs A character vector specifying labels for columns on the left
#'   side of the forest plot (see Details).
#' @param rightlabs A character vector specifying labels for columns on the
#'   right side of the forest plot (see Details).
#' @param lab.NA A character string to label missing values.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param big.mark A character used as thousands separator.
#' @param addrow.subgroups A logical value indicating whether an empty
#'   row is printed between results for subgroups.
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
#' Arguments \code{leftcols} and \code{rightcols} can be used to specify
#' columns printed on the left or right side, respectively. Arguments
#' \code{leftlabs} and \code{rightlabs} can be used to provide headings for
#' these columns. Use \code{NA} for a column heading to retain its default
#' label.
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
#' ma <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, sm = "SMD")
#'
#' # Conduct two subgroup analyses
#' #
#' ma1 <- update(ma, subgroup = age, subgroup.name = "Age group")
#' ma2 <- update(ma, subgroup = region, subgroup.name = "Region")
#'
#' # Combine subgroup meta-analyses and show forest plot with subgroup
#' # results
#' #
#' mb <- metabind(ma1, ma2)
#' mb
#' forest(mb)
#'
#' @method forest metabind
#' @export

forest.metabind <- function(x,
                            #
                            layout = gs("layout"),
                            #
                            common = x$common,
                            random = x$random,
                            overall = x$overall,
                            subgroup = FALSE,
                            hetstat = FALSE,
                            overall.hetstat = x$overall.hetstat,
                            prediction = x$prediction,
                            #
                            col = NULL,
                            print = NULL,
                            type = NULL,
                            digits = NULL,
                            #
                            leftcols = NULL, leftlabs = NULL,
                            rightcols = NULL, rightlabs = NULL,
                            #
                            lab.NA = "",
                            #
                            scientific.pval = gs("scientific.pval"),
                            big.mark = gs("big.mark"),
                            #
                            addrow.subgroups = NULL,
                            #
                            ...) {
  
  #
  #
  # (1) Check and set arguments
  #
  #

  chkclass(x, "metabind")
  x <- updateversion(x)
  #
  x$k.w.orig <- x$k.w
  #
  x$k.w <- x$k.study.w <- x$k.all.w <- x$k.TE.w <-
    as.vector(table(x$data$subgroup)[unique(x$data$subgroup)])
  #
  layout <- setchar(layout, gs("layout4forest"), length = 1)
  #
  # Drop internally used arguments from '...'
  #
  drop_args <- NULL
  #
  replace_args <- list()
  #
  args <- drop_from_dots(list(...), drop_args, replace_args)
  #
  if (!x$with.subgroups) {
    overall <- FALSE
    overall.hetstat <- FALSE
  }
  #
  overall <- replaceNULL(overall, FALSE)
  overall.hetstat <- replaceNULL(overall.hetstat, FALSE)
  #
  chklogical(common)
  chklogical(random)
  chklogical(overall)
  chklogical(subgroup)
  chklogical(overall.hetstat)
  #
  chkchar(lab.NA)
  #
  # Types
  #
  type <- as_named_list(type)
  type.list <-
    merge_legacy(type, args, type_defaults(class = "metabind"), "type")
  #
  types <- type.list$res
  args <- type.list$args
  miss.type <- type.list$miss
  null.type <- type.list$null
  #
  types <- type.list$res
  #
  type <- types$type
  type.common <- types$common
  type.random <- types$random
  type.predict <- types$predict
  #
  miss.type.common <- miss.type$common
  miss.type.random <- miss.type$random
  miss.type.predict <- miss.type$predict
  #
  if (is.null(type))
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
  valid_types <- c("square", "diamond", "predict", "circle")
  #
  type <- setchar(type, valid_types)
  #
  if (!miss.type.common)
    type.common <- setchar(type.common, valid_types, length = 1)
  #
  if (!miss.type.random)
    type.random <- setchar(type.random, valid_types, length = 1)
  #
  if (!miss.type.predict)
    type.predict <- setchar(type.predict, valid_types, length = 1)
  #
  # Define colours
  #
  col <- as_named_list(col)
  col.list <- merge_legacy(col, args, col_defaults(x, "metabind"), "col")
  #
  col <- col.list$res
  args <- col.list$args
  #
  cols <- col.list$res
  #
  col <- cols$col
  col.square <- cols$square
  col.square.lines <- cols$square.lines
  col.circle <- cols$circle
  col.circle.lines <- cols$circle.lines
  col.diamond <- cols$diamond
  col.diamond.common <- cols$diamond.common
  col.diamond.random <- cols$diamond.random
  col.diamond.lines <- cols$diamond.lines
  col.diamond.lines.common <- cols$diamond.lines.common
  col.diamond.lines.random <- cols$diamond.lines.random
  col.predict <- cols$predict
  col.predict.lines <- cols$predict.lines
  #
  # Colour schemes for layouts:
  # - colors[1] - vertical line for common effect or random effects model
  # - colors[2] - diamond for meta-analysis results
  # - colors[3] - outer lines of diamonds
  # - colors[4] - square  colour
  # - colors[5] - outer lines of squares
  # - colors[6] - line color (axes, reference line, header lines)
  # - colors[7] - prediction interval
  # - colors[8] - outer lines of prediction intervals
  # - colors[9] - subgroups
  # - colors[10] - color within squares or circles
  # - colors[11] - color of circles and outer lines of circles
  #
  colors <- layout_colors(layout)
  #
  if (is.null(col.square))
    col.square <- colors[4]
  if (is.null(col.square.lines))
    col.square.lines <- col.square
  if (is.null(col.circle))
    col.circle <- colors[4]
  if (is.null(col.circle.lines))
    col.circle.lines <- col.circle
  if (is.null(col.diamond))
    col.diamond <- colors[2]
  if (is.null(col.diamond.common))
    col.diamond.common <- col.diamond
  if (is.null(col.diamond.random))
    col.diamond.random <- col.diamond
  if (is.null(col.diamond.lines))
    col.diamond.lines <- colors[3]
  if (is.null(col.diamond.lines.common))
    col.diamond.lines.common <- col.diamond.lines
  if (is.null(col.diamond.lines.random))
    col.diamond.lines.random <- col.diamond.lines
  if (is.null(col.predict))
    col.predict <- colors[7]
  if (is.null(col.predict.lines))
    col.predict.lines <- colors[8]
  #
  # Print information
  #
  print.list <-
    merge_legacy(print, args, print_defaults(x, "metabind"), "print")
  #
  print <- print.list$res
  args <- print.list$args
  #
  print.k <- print$k
  print.pval <- print$pval
  print.stat <- print$stat
  print.I2 <- print$I2
  print.tau2 <- print$tau2
  print.tau <- print$tau
  print.Q <- print$Q
  print.pval.Q <- print$pval.Q
  print.pval.Q.b <- print$pval.Q.b
  print.subgroup.labels <- print$subgroup.labels
  #
  print.pval.Q.b <- replaceNULL(print.pval.Q.b, x$with.subgroups)
  print.tau2 <- replaceNULL(print.tau2, !x$with.subgroups && random)
  print.subgroup.labels <- replaceNULL(print.subgroup.labels, x$with.subgroups)
  #
  chklogical(print.k)
  chklogical(print.pval)
  chklogical(print.stat)
  chklogical(print.I2)
  chklogical(print.tau2)
  chklogical(print.tau)
  chklogical(print.Q)
  chklogical(print.pval.Q)
  chklogical(print.pval.Q.b)
  chklogical(print.subgroup.labels)
  #
  # Digits
  #
  digits <- as_named_list(digits)
  digits.list <-
    merge_legacy(digits, args, digits_defaults("metabind"), "digits")
  #
  digitss <- digits.list$res
  args <- digits.list$args
  #
  digits <- digitss$digits
  digits.se <- digitss$se
  digits.stat <- digitss$stat
  digits.pval <- digitss$pval
  digits.pval.Q <- digitss$pval.Q
  digits.Q <- digitss$Q
  digits.tau2 <- digitss$tau2
  digits.tau <- digitss$tau
  digits.I2 <- digitss$I2
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  chknumeric(digits.pval.Q, min = 1, length = 1)
  chknumeric(digits.Q, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  #
  chklogical(scientific.pval)
  #
  addrow.subgroups <- replaceNULL(addrow.subgroups, print.subgroup.labels)
  chklogical(addrow.subgroups)
  
  #
  #
  # (2) Columns and labels shown in forest plot
  #
  #
  
  # Columns and labels on the left side
  #
  if (is.null(leftcols)) {
    leftcols <- "studlab"
    #
    if (print.k)
      leftcols <- c(leftcols, "k")
    #
    if (print.pval.Q.b)
      leftcols <- c(leftcols, "pval.Q.b")
    #
    if (print.tau2)
      leftcols <- c(leftcols, "tau2")
    #
    if (print.tau)
      leftcols <- c(leftcols, "tau")
    #
    if (print.I2)
      leftcols <- c(leftcols, "I2")
    #
    if (print.Q)
      leftcols <- c(leftcols, "Q")
    #
    if (print.pval.Q)
      leftcols <- c(leftcols, "pval.Q")
    #
    if (layout != "meta8") {
      if (layout == "RevMan5" || layout == "JAMA")
        leftcols <- c(leftcols, "effect.ci")
      else if (layout != "BMJ")
        leftcols <- c(leftcols, "effect", "ci")
      #
      if (print.stat && !(layout == "BMJ" || layout == "meta8"))
        leftcols <- c(leftcols, "statistic")
      #
      if (print.pval && !(layout == "BMJ" || layout == "meta8"))
        leftcols <- c(leftcols, "pval")
    }
  }
  #
  # Argument leftlabs
  #
  if (is.null(leftlabs))
    leftlabs <- rep_len(NA, length(leftcols))
  else {
    if (length(leftlabs) != length(leftcols))
      stop("Argument 'leftlabs' must provide ", length(leftcols),
           " labels (leftcols: ",
           paste(leftcols, collapse = ", "), ")",
           call. = FALSE)
  }
  #
  label.studlab <- ifelse(x$with.subgroups, "Subgroup", "Method")
  #
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

  # Columns and labels on the right side
  #
  if (is.null(rightcols)) {
    if (layout == "BMJ") {
      rightcols <- c("effect.ci")
      #
      if (print.stat)
        rightcols <- c(rightcols, "statistic")
      #
      if (print.pval)
        rightcols <- c(rightcols, "pval")
    }
    else if (layout == "meta8") {
      rightcols <- c("effect", "ci")
      #
      if (print.stat)
        rightcols <- c(rightcols, "statistic")
      #
      if (print.pval)
        rightcols <- c(rightcols, "pval")
      #
      if (print.tau)
        rightcols <- c(rightcols, "tau")
      #
      if (print.I2)
        rightcols <- c(rightcols, "I2")
    }
    else
      rightcols <- FALSE
  }
  #
  # Argument rightlabs
  #
  if (is.null(rightlabs))
    rightlabs <- rep_len(NA, length(rightcols))
  else {
    if (length(rightlabs) != length(rightcols))
      stop("Argument 'rightlabs' must provide ", length(rightcols),
           " labels (rightcols: ",
           paste(rightcols, collapse = ", "), ")",
           call. = FALSE)
  }
  #
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
  

  #
  #
  # (3) Round and format results
  #
  #

  if (x$with.subgroups) {
    x$data$Q.b <- round(x$data$Q.b, digits.Q)
    x$data$pval.Q.b <- formatPT(x$data$pval.Q.b,
                                lab = FALSE,
                                digits = digits.pval.Q,
                                scientific = scientific.pval,
                                lab.NA = lab.NA, trim = "left")
  }
  else {
    x$data$Q <- formatN(x$data$Q, digits = digits.Q, lab.NA = lab.NA,
                        big.mark = big.mark)
    x$Q <- round(x$Q, digits.Q)
    x$data$pval.Q <-
      formatPT(x$data$pval.Q,
               lab = FALSE,
               digits = digits.pval.Q,
               scientific = scientific.pval,
               lab.NA = lab.NA, trim = "left")
  }
  #
  x$data$tau2 <-
    formatPT(x$data$tau2, digits = digits.tau2,
             big.mark = big.mark,
             lab = FALSE, lab.NA = lab.NA, trim = "left")
  #
  x$data$tau <-
    formatPT(x$data$tau, digits = digits.tau,
             big.mark = big.mark,
             lab = FALSE, lab.NA = lab.NA, trim = "left")
  #
  x$data$I2 <-
    percent_label(x$data$I2, "", digits.I2, formatN, lab.NA = lab.NA)
  #
  x$data$k <- as.character(x$data$k)
  #
  model <- x$data$model
  #
  if (!is.null(type.common))
    type[model == "common"] <- type.common
  #
  if (!is.null(type.random))
    type[model == "random"] <- type.random
  #
  if (!is.null(type.predict))
    type[model == "predict"] <- type.predict
  #
  cols.square <- vector("character", length(type))
  cols.square.lines <- vector("character", length(type))
  #
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
      cols.square[i] <- col.square
      cols.square.lines[i] <- col.square.lines
    }
  }
  #
  if (!x$samedata) {
    overall <- FALSE
    overall.hetstat <- FALSE
  }

  #
  #
  #
  # (4) Create forest plot
  #
  #
  
  x.forest <- x
  #
  if (x$with.subgroups) {
    ma <- metagen(x.forest$TE, x.forest$seTE,
                  lower = x.forest$lower,
                  upper = x.forest$upper,
                  subgroup = x.forest$subgroup,
                  studlab = x.forest$studlab,
                  data = x.forest,
                  common = common,
                  random = random,
                  prediction = prediction,
                  overall = overall,
                  print.subgroup.name = FALSE,
                  transf = TRUE)
    #
    class(ma) <- c(class(ma), "is.metabind")
    ma$classes <- x$classes
    ma$MH.exact <- x$list$meth$MH.exact
    #
    ma$sm <- x.forest$sm
    #
    ma$TE[ma$data$type == "predict"] <- NA
    ma$statistic <- x.forest$statistic
    ma$pval <- x.forest$pval
    #
    ma$n.harmonic.mean <- x.forest$n.harmonic.mean
    ma$t.harmonic.mean <- x.forest$t.harmonic.mean
    ma$n.harmonic.mean.ma <- x.forest$n.harmonic.mean.ma
    ma$t.harmonic.mean.ma <- x.forest$t.harmonic.mean.ma
    ma$n.harmonic.mean.w <- x.forest$n.harmonic.mean.w
    ma$t.harmonic.mean.w <- x.forest$t.harmonic.mean.w
    #
    if (x$samedata) {
      ma$method <- x.forest$method
      ma$method.random <- x.forest$method.random
      ma$method.predict <- x.forest$method.predict
      ma$method.tau <- x.forest$method.tau
      #
      ma$TE.common <- x.forest$TE.common
      ma$seTE.common <- x.forest$seTE.common
      ma$lower.common <- x.forest$lower.common
      ma$upper.common <- x.forest$upper.common
      ma$statistic.common <- x.forest$statistic.common
      ma$pval.common <- x.forest$pval.common
      ma$zval.common <- x.forest$zval.common
      ma$text.common <- x.forest$text.common
      #
      ma$w.common <- ifelse(is.na(ma$w.common), ma$w.common, NA)
      #
      ma$TE.random <- x.forest$TE.random
      ma$seTE.random <- x.forest$seTE.random
      ma$lower.random <- x.forest$lower.random
      ma$upper.random <- x.forest$upper.random
      ma$statistic.random <- x.forest$statistic.random
      ma$pval.random <- x.forest$pval.random
      ma$zval.random <- x.forest$zval.random
      ma$df.random <- x.forest$df.random
      ma$text.random <- x.forest$text.random
      #
      ma$w.random <- ifelse(is.na(ma$w.random), ma$w.random, NA)
      #
      ma$lower.predict <- x.forest$lower.predict
      ma$upper.predict <- x.forest$upper.predict
      ma$text.predict <- replaceNULL(x.forest$text.predict, gs("text.predict"))
      #
      ma$backtransf <- x.forest$backtransf
      #
      ma$null.effect <- x.forest$null.effect
    }
    #
    forest_args <-
      list(
        x = ma,
        #
        layout = layout,
        #
        leftcols = leftcols, leftlabs = leftlabs,
        rightcols = rightcols, rightlabs = rightlabs,
        #
        print.subgroup.labels = print.subgroup.labels,
        addrow.subgroups = addrow.subgroups,
        #
        hetstat = hetstat, overall.hetstat = overall.hetstat,
        lab.NA = lab.NA,
        #
        type = list(study = type),
        col.square = cols.square,
        col.square.lines = cols.square.lines,
        #
        pooled.totals = FALSE,
        common.subgroup = FALSE, random.subgroup = FALSE,
        prediction.subgroup = FALSE,
        #
        col.predict = col.predict,
        #
        weight.study = "same",
        #
        digits = list(digits = digits,
                      se = digits.se,
                      stat = digits.stat,
                      pval = digits.pval,
                      pval.Q = digits.pval.Q,
                      Q = digits.Q,
                      tau2 = digits.tau2,
                      tau = digits.tau,
                      I2 = digits.I2),
        #
        scientific.pval = scientific.pval,
        big.mark = big.mark,
        #
        test.subgroup = FALSE)
  }
  else {
    if (isCol(x.forest$data, "subgroup"))
      ma <- metagen(x.forest$data$TE, x.forest$data$seTE,
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
                          method.tau = "DL", method.tau.ci = "",
                          transf = TRUE)
    else
      ma <- metagen(x.forest$data$TE, x.forest$data$seTE,
                    lower = x.forest$data$lower,
                    upper = x.forest$data$upper,
                    studlab = x.forest$data$studlab,
                    data = x.forest,
                    common = common,
                    random = random,
                    prediction = prediction,
                    overall = overall,
                    transf = TRUE)
    #
    class(ma) <- c(class(ma), "is.metabind")
    ma$classes <- x$classes
    ma$MH.exact <- x$list$meth$MH.exact
    #
    ma$sm <- x.forest$sm
    ma$with.subgroups <- TRUE
    #
    ma$n.harmonic.mean <- x.forest$n.harmonic.mean
    ma$t.harmonic.mean <- x.forest$t.harmonic.mean
    #
    ma$TE[ma$data$type == "predict"] <- NA
    #
    forest_args <-
      list(
        x = ma,
        #
        layout = layout,
        #
        leftcols = leftcols, leftlabs = leftlabs,
        rightcols = rightcols, rightlabs = rightlabs,
        #
        hetstat = hetstat, overall.hetstat = overall.hetstat,
        lab.NA = lab.NA,
        #
        type.study = type, col.square = cols.square,
        col.square.lines = cols.square.lines,
        #
        pooled.totals = FALSE,
        #
        col.predict = col.predict,
        #
        weight.study = "same",
        #
        digits = list(digits = digits,
                      se = digits.se,
                      stat = digits.stat,
                      pval = digits.pval,
                      pval.Q = digits.pval.Q,
                      Q = digits.Q,
                      tau2 = digits.tau2,
                      tau = digits.tau,
                      I2 = digits.I2),
        #
        scientific.pval = scientific.pval,
        big.mark = big.mark,
        #
        test.subgroup = FALSE)
  }
  #
  text.details <-
    catmeth(ma,
            common, random, ma$prediction, overall || print.pval, TRUE,
            #
            func.transf = x$func.transf,
            backtransf = x$backtransf, func.backtransf = x$func.backtransf,
            #
            big.mark = big.mark, digits = digits,
            digits.tau = digits.tau,
            text.tau = gs("text.tau"), text.tau2 = gs("text.tau2"),
            #
            print.tau2 = print.tau2, print.tau2.ci = FALSE,
            print.tau = print.tau, print.tau.ci = FALSE,
            #
            print.I2 = print.I2, text.I2 = gs("text.I2"),
            #
            print.df = TRUE, prediction.subgroup = FALSE,
            #
            forest = TRUE)
  #
  forest_args$x$.text.details.methods <- text.details
  #
  res <- do.call("forest.meta", c(forest_args, args))
  #
  return(invisible(res))
}
