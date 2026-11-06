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
#' @param prediction A logical indicating whether prediction intervals
#'   should be plotted.
#' @param overall A logical indicating whether overall summaries
#'   should be plotted. This argument is useful in a meta-analysis
#'   with subgroups if summaries should only be plotted on group
#'   level.
#' @param studlab A  variable containing labels for individual results.
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
#' @param details A logical specifying whether details on statistical
#'   methods should be printed.
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
                            common = !is.null(x$common),
                            random = !is.null(x$random),
                            prediction = !is.null(x$predict),
                            #
                            overall = common || random,
                            #
                            studlab = NULL,
                            subgroup = NULL,
                            hetstat = FALSE,
                            overall.hetstat = x$overall.hetstat,
                            #
                            col = NULL,
                            print = NULL,
                            type = NULL,
                            digits = NULL,
                            #
                            leftcols = NULL, leftlabs = NULL,
                            rightcols = NULL, rightlabs = NULL,
                            #
                            lab.NA = gs("lab.NA"),
                            #
                            scientific.pval = gs("scientific.pval"),
                            big.mark = gs("big.mark"),
                            #
                            addrow.subgroups = NULL,
                            #
                            details = FALSE,
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
  with.subgroups <- !is.null(x$subgroup)
  #
  if (!with.subgroups) {
    overall <- FALSE
    overall.hetstat <- FALSE
  }
  #
  overall <- replaceNULL(overall, FALSE)
  overall.hetstat <- replaceNULL(overall.hetstat, FALSE)
  #
  missing.common <- missing(common)
  missing.random <- missing(random)
  missing.prediction <- missing(prediction)
  #
  chklogical(common)
  chklogical(random)
  chklogical(prediction)
  #
  if (is.null(x$common)) {
    if (!missing.common && common)
      warning("Argument 'common' set to FALSE as ",
              "no common effect results are available.",
              call. = FALSE)
    #
    common <- FALSE
  }
  #
  if (is.null(x$random)) {
    if (!missing.random && random)
      warning("Argument 'random' set to FALSE as ",
              "no random effects results are available.",
              call. = FALSE)
    #
    random <- FALSE
  }
  #
  if (is.null(x$predict)) {
    if (!missing.prediction && prediction)
      warning("Argument 'prediction' set to FALSE as ",
              "no predict intervals are available.",
              call. = FALSE)
    #
    prediction <- FALSE
  }
  #
  chklogical(overall)
  chklogical(overall.hetstat)
  #
  if (!common && !random && !prediction) {
    warning("Nothing to show in forest plot.", call. = FALSE)
    return(invisible(NULL))
  }
  #
  if (!common)
    x$common <- NULL
  #
  if (!random)
    x$random <- NULL
  #
  if (!prediction)
    x$predict <- NULL
  #
  # Data set used in forest plot
  #
  fdat <- as.data.frame(x)
  #
  fdat$type <- ifelse(fdat$model == "Prediction interval", "predict", "square")
  #
  # Set certain information to missing
  #
  if (common && random) {
    sel.c <- fdat$model == "Common effect model"
    #
    fdat$tau2[sel.c] <- NA
    fdat$tau[sel.c] <- NA
    fdat$I2[sel.c] <- NA
  }
  #
  if (prediction && (common || random)) {
    sel.p <- fdat$model == "Prediction interval"
    #
    fdat$tau2[sel.p] <- NA
    fdat$tau[sel.p] <- NA
    fdat$I2[sel.p] <- NA
  }
  #
  # Act upon different summary measures
  #
  backtransf <- x$backtransf
  func.backtransf <- NULL
  #
  fdat$TE.orig <- fdat$TE
  fdat$seTE.orig <- fdat$seTE
  fdat$lower.orig <- fdat$lower
  fdat$upper.orig <- fdat$upper
  #
  pscale <- irscale <- irunit <- NULL
  #
  if (x$same.sm)
    label.sm <- x$sm
  else {
    if (x$sm.family == "relative"){
      label.sm <- "Effect"
      func.backtransf <- "exp"
    }
    else if (x$sm.family == "correlation") {
      label.sm <- "Correlation"
      backtransf <- FALSE
      #
      for (i in seq_len(nrow(fdat))) {
        fdat$TE[i] <- backtransf(fdat$TE[i], fdat$sm[i])
        fdat$lower[i] <- backtransf(fdat$lower[i], fdat$sm[i])
        fdat$upper[i] <- backtransf(fdat$upper[i], fdat$sm[i])
      }
    }
    else if (x$sm.family == "mean") {
      label.sm <- "Mean"
      backtransf <- FALSE
      #
      for (i in seq_len(nrow(fdat))) {
        fdat$TE[i] <- backtransf(fdat$TE[i], fdat$sm[i])
        fdat$lower[i] <- backtransf(fdat$lower[i], fdat$sm[i])
        fdat$upper[i] <- backtransf(fdat$upper[i], fdat$sm[i])
      }
    }
    else if (x$sm.family == "proportion") {
      pscale <- makeunique(x$pscale, 1)
      #
      if (pscale == 1)
        label.sm <- "Proportion"
      else
        label.sm <- "Events"
      #
      backtransf <- FALSE
      #
      for (i in seq_len(nrow(fdat))) {
        fdat$TE[i] <- pscale *
          backtransf(fdat$TE[i], fdat$sm[i], fdat$n.harmonic.mean[i])
        fdat$lower[i] <- pscale *
          backtransf(fdat$lower[i], fdat$sm[i], fdat$n.harmonic.mean[i])
        fdat$upper[i] <- pscale *
          backtransf(fdat$upper[i], fdat$sm[i], fdat$n.harmonic.mean[i])
      }
    }
    else if (x$sm.family == "rate") {
      irscale <- makeunique(x$irscale, 1)
      irunit <- makeunique(x$irunit, "units")
      #
      if (irscale == 1)
        label.sm <- "Rates"
      else
        label.sm <- "Events"
      #
      backtransf <- FALSE
      #
      for (i in seq_len(nrow(fdat))) {
        fdat$TE[i] <- irscale *
          backtransf(fdat$TE[i], fdat$sm[i], fdat$n.harmonic.mean[i])
        fdat$lower[i] <- irscale *
          backtransf(fdat$lower[i], fdat$sm[i], fdat$n.harmonic.mean[i])
        fdat$upper[i] <- irscale *
          backtransf(fdat$upper[i], fdat$sm[i], fdat$n.harmonic.mean[i])
      }
    }
    else
      label.sm <- "Effect"
  }
  #
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  studlab.model <-
    replaceNULL(mc$studlab == "model", FALSE, chklength = TRUE) ||
    (replaceNULL(mc$studlab == "text", FALSE, chklength = TRUE) &&
       all(fdat$text == fdat$model))
  #
  subgroup.model <-
    replaceNULL(mc$subgroup == "model", FALSE, chklength = TRUE) ||
    (replaceNULL(mc$subgroup == "text", FALSE, chklength = TRUE) &&
       all(fdat$text == fdat$model))
  #
  studlab <- catch("studlab", mc, fdat, sfsp)
  #
  if (!is.null(studlab)) {
    if (length(studlab) == 1)
      studlab <- rep(studlab, nrow(fdat))
    #
    chklength(studlab, nrow(fdat),
              text =
                paste0("Argument 'studlab' must be of length ",
                       nrow(fdat), (", i.e., number of results.")))
  }
  #
  subgroup <- catch("subgroup", mc, fdat, sfsp)
  #
  if (!is.null(subgroup)) {
    if (length(subgroup) == 1)
      subgroup <- rep(subgroup, nrow(fdat))
    #
    chklength(subgroup, nrow(fdat),
              text =
                paste0("Argument 'subgroup' must be of length ",
                       nrow(fdat), (", i.e., number of results.")))
  }
  
  
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
    type <- fdat$type
  else {
    if (length(type) == 1 & nrow(fdat) > 1)
      type <- rep(type, nrow(fdat))
    else if (length(type) != nrow(fdat))
      stop("Argument 'type' must be a character vector of length 1 or ",
           nrow(fdat), ".",
           call. = FALSE)
    #
    fdat$type <- type
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
  print.studlab <- print$studlab
  print.model <-
    replaceNULL(print$model,
                common + random > 1 && !(studlab.model || subgroup.model))
  print.sm <- replaceNULL(print$sm, !x$same.sm)
  print.k <- print$k
  print.pval <- print$pval
  print.stat <- print$stat
  print.null.effect <- print$null.effect
  print.I2 <- print$I2
  print.tau2 <- print$tau2
  print.tau2.ci <- print$tau2.ci
  print.tau <- print$tau
  print.tau.ci <- print$tau.ci
  print.Q <- print$Q
  print.pval.Q <- print$pval.Q
  print.pval.Q.b <- print$pval.Q.b
  print.subgroup.name <- print$subgroup.name
  #
  print.pval.Q.b <- replaceNULL(print.pval.Q.b, with.subgroups)
  print.tau2 <- replaceNULL(print.tau2, !with.subgroups && random)
  print.tau <- replaceNULL(print.tau, !with.subgroups && random)
  print.subgroup.name <- replaceNULL(print.subgroup.name, with.subgroups)
  #
  chklogical(print.studlab)
  chklogical(print.model)
  chklogical(print.sm)
  chklogical(print.k)
  chklogical(print.pval)
  chklogical(print.stat)
  chklogical(print.null.effect)
  chklogical(print.I2)
  chklogical(print.tau2)
  chklogical(print.tau2.ci)
  chklogical(print.tau)
  chklogical(print.tau.ci)
  chklogical(print.Q)
  chklogical(print.pval.Q)
  chklogical(print.pval.Q.b)
  chklogical(print.subgroup.name)
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
  addrow.subgroups <- replaceNULL(addrow.subgroups, print.subgroup.name)
  chklogical(addrow.subgroups)
  #
  chklogical(details)
  
  #
  #
  # (2) Columns and labels shown in forest plot
  #
  #
  
  # Columns on the left side
  #
  if (is.null(leftcols)) {
    if (print.studlab)
      leftcols <- "studlab"
    else
      NULL
    #
    if (print.model)
      leftcols <- c(leftcols, "model")
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
        leftcols <- c(leftcols, if (print.sm) "sm", "effect.ci")
      else if (layout != "BMJ")
        leftcols <- c(leftcols, if (print.sm) "sm", "effect", "ci")
      #
      if (print.stat && !(layout == "BMJ" || layout == "meta8"))
        leftcols <- c(leftcols, "statistic")
      #
      if (print.pval && !(layout == "BMJ" || layout == "meta8"))
        leftcols <- c(leftcols, "pval")
      #
      if (print.null.effect && !(layout == "BMJ" || layout == "meta8"))
        leftcols <- c(leftcols, "null.effect")
    }
  }
  #  
  # Columns on the right side
  #
  if (is.null(rightcols)) {
    if (layout == "BMJ") {
      rightcols <- c("effect.ci", if (print.sm) "sm")
      #
      if (print.stat)
        rightcols <- c(rightcols, "statistic")
      #
      if (print.pval)
        rightcols <- c(rightcols, "pval")
      #
      if (print.null.effect)
        rightcols <- c(rightcols, "null.effect")
    }
    else if (layout == "meta8") {
      rightcols <- c("effect", "ci", if (print.sm) "sm")
      #
      if (print.stat)
        rightcols <- c(rightcols, "statistic")
      #
      if (print.pval)
        rightcols <- c(rightcols, "pval")
      #
      if (print.null.effect)
        rightcols <- c(rightcols, "null.effect")
    }
    else
      rightcols <- FALSE
  }
  
  # Arguments leftlabs and rightlabs
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
  label.studlab <- ifelse(with.subgroups, "Subgroup", "Method")
  #
  cols <- c("studlab", "model",
            "k",
            "sm",
            "tau2", "tau",
            "pval.Q.b", "null.effect")
  labs <- c(label.studlab, "Model",
            if (with.subgroups) "Number of\nStudies" else NA,
            "Effect\nmeasure",
            if (with.subgroups) "Between-study\nvariance" else NA,
            if (with.subgroups) "Between-study\nSD" else NA,
            "Interaction\nP-value", "Null\neffect")
  #
  if (any(is.na(leftlabs))) {
    for (i in seq_along(cols))
      leftlabs <- newlab(leftcols, leftlabs, cols[i], labs[i])
  }
  #
  if (any(is.na(rightlabs))) {
    for (i in seq_along(cols))
      rightlabs <- newlab(rightcols, rightlabs, cols[i], labs[i])
  }
  
  #
  #
  # (3) Round and format results
  #
  #
  
  if (with.subgroups) {
    fdat$Q.b <- round(fdat$Q.b, digits.Q)
    fdat$pval.Q.b <- formatPT(fdat$pval.Q.b,
                              lab = FALSE,
                              digits = digits.pval.Q,
                              scientific = scientific.pval,
                              lab.NA = lab.NA, trim = "left")
    #
    fdat$Q.b.random <- round(fdat$Q.b.random, digits.Q)
    fdat$pval.Q.b.random <- formatPT(fdat$pval.Q.b.random,
                                     lab = FALSE,
                                     digits = digits.pval.Q,
                                     scientific = scientific.pval,
                                     lab.NA = lab.NA, trim = "left")
  }
  else {
    fdat$Q <- formatN(fdat$Q, digits = digits.Q, lab.NA = lab.NA,
                      big.mark = big.mark)
    fdat$pval.Q <-
      formatPT(fdat$pval.Q,
               lab = FALSE,
               digits = digits.pval.Q,
               scientific = scientific.pval,
               lab.NA = lab.NA, trim = "left")
  }
  #
  fdat$tau2 <-
    formatPT(fdat$tau2, digits = digits.tau2,
             big.mark = big.mark,
             lab = FALSE, lab.NA = lab.NA, trim = "left")
  #
  fdat$tau <-
    formatPT(fdat$tau, digits = digits.tau,
             big.mark = big.mark,
             lab = FALSE, lab.NA = lab.NA, trim = "left")
  #
  fdat$I2 <-
    percent_label(fdat$I2, "", digits.I2, formatN, lab.NA = lab.NA)
  #
  fdat$k <- as.character(fdat$k)
  #
  model <- fdat$model
  #
  if (!is.null(type.common))
    type[model == "Common effect model"] <- type.common
  #
  if (!is.null(type.random))
    type[model == "Random effects model"] <- type.random
  #
  if (!is.null(type.predict))
    type[model == "Prediction interval"] <- type.predict
  #
  cols.square <- vector("character", length(type))
  cols.square.lines <- vector("character", length(type))
  #
  for (i in seq_along(type)) {
    if (type[i] == "diamond") {
      if (model[i] == "Common effect model") {
        cols.square[i] <- col.diamond.common
        cols.square.lines[i] <- col.diamond.lines.common
      }
      else if (model[i] == "Random effects model") {
        cols.square[i] <- col.diamond.random
        cols.square.lines[i] <- col.diamond.lines.random
      }
      else {
        cols.square[i] <- col.diamond.random
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
  #
  #
  # (4) Create forest plot
  #
  #
  
  if (!with.subgroups) {
    if (is.null(studlab))
      studlab <- fdat$name
    #
    ma <- metagen(fdat$TE, fdat$seTE,
                  lower = fdat$lower, upper = fdat$upper,
                  studlab = studlab,
                  #
                  data = fdat,
                  #
                  sm = label.sm,
                  #
                  common = FALSE,
                  random = FALSE,
                  prediction = FALSE,
                  overall = FALSE,
                  #
                  transf = TRUE,
                  #
                  backtransf = backtransf,
                  func.backtransf = func.backtransf,
                  #
                  subgroup = subgroup, 
                  print.subgroup.name = print.subgroup.name,
                  #
                  label.left = makeunique(fdat$label.left, ""),
                  label.right = makeunique(fdat$label.right, ""),
                  #
                  method.tau = "DL", method.tau.ci = "")
    #
    ma$statistic <- fdat$statistic
    ma$pval <- fdat$pval
    #
    # Set certain information to missing
    #
    ma$TE[ma$data$model == "Prediction interval"] <- NA
    #
    forest_args_specific <- list(x = ma)
  }
  #
  if (with.subgroups) {
    if (is.null(studlab))
      studlab <- fdat$subgroup.levels
    #
    ma <- metagen(fdat$TE.w, fdat$seTE.w,
                  lower = fdat$lower.w,
                  upper = fdat$upper.w,
                  subgroup = fdat$subgroup.name,
                  studlab = studlab,
                  data = fdat,
                  common = common,
                  random = random,
                  prediction = prediction,
                  overall = overall,
                  backtransf = x$backtransf,
                  print.subgroup.name = FALSE,
                  transf = TRUE)
    #
    class(ma) <- c(class(ma), "is.metabind")
    ma$classes <- x$classes
    ma$MH.exact <- x$MH.exact
    #
    ma$sm <- fdat$sm
    #
    ma$TE[ma$data$model == "Prediction interval"] <- NA
    ma$statistic <- fdat$statistic.w
    ma$pval <- fdat$pval.w
    #
    ma$n.harmonic.mean <- fdat$n.harmonic.mean.w
    ma$t.harmonic.mean <- fdat$t.harmonic.mean.w
    ma$n.harmonic.mean.ma <- fdat$n.harmonic.mean
    ma$t.harmonic.mean.ma <- fdat$t.harmonic.mean
    #
    if (x$same.data) {
      ma$method <- fdat$method
      ma$method.random <- fdat$method.random
      ma$method.predict <- fdat$method.predict
      ma$method.tau <- fdat$method.tau
      #
      ma$sm <- makeunique(fdat$sm)
      #
      ma$TE.common <- fdat$TE.common
      ma$seTE.common <- fdat$seTE.common
      ma$lower.common <- fdat$lower.common
      ma$upper.common <- fdat$upper.common
      ma$statistic.common <- fdat$statistic.common
      ma$pval.common <- fdat$pval.common
      ma$zval.common <- fdat$zval.common
      ma$text.common <- fdat$text.common
      #
      ma$w.common <- ifelse(is.na(ma$w.common), ma$w.common, NA)
      #
      ma$TE.random <- fdat$TE.random
      ma$seTE.random <- fdat$seTE.random
      ma$lower.random <- fdat$lower.random
      ma$upper.random <- fdat$upper.random
      ma$statistic.random <- fdat$statistic.random
      ma$pval.random <- fdat$pval.random
      ma$zval.random <- fdat$zval.random
      ma$df.random <- fdat$df.random
      ma$text.random <- fdat$text.random
      #
      ma$w.random <- ifelse(is.na(ma$w.random), ma$w.random, NA)
      #
      ma$lower.predict <- fdat$lower.predict
      ma$upper.predict <- fdat$upper.predict
      ma$text.predict <- replaceNULL(fdat$text.predict, gs("text.predict"))
      #
      ma$null.effect <- makeunique(fdat$null.effect)
    }
    #
    forest_args_specific <-
      list(x = ma,
           #
           addrow.subgroups = addrow.subgroups,
           #
           common.subgroup = FALSE, random.subgroup = FALSE,
           prediction.subgroup = FALSE)
  }
  #
  forest_args_general <-
    list(
      layout = layout,
      #
      leftcols = leftcols, leftlabs = leftlabs,
      rightcols = rightcols, rightlabs = rightlabs,
      #
      hetstat = hetstat, overall.hetstat = overall.hetstat,
      lab.NA = lab.NA,
      #
      pooled.totals = FALSE,
      #
      weight.study = "same",
      #
      col = list(square = cols.square,
                 square.lines = cols.square.lines,
                 circle = col.circle,
                 circle.lines = col.circle.lines,
                 diamond = col.diamond,
                 diamond.common = col.diamond.common,
                 diamond.random = col.diamond.random,
                 diamond.lines = col.diamond.lines,
                 diamond.lines.common = col.diamond.lines.common,
                 diamond.lines.random = col.diamond.lines.random,
                 predict = col.predict,
                 predict.lines = col.predict.lines),
      #
      print = list(I2 = print.I2,
                   tau2 = print.tau2,
                   tau2.ci = print.tau2.ci,
                   tau = print.tau,
                   tau.ci = print.tau.ci,
                   Q = print.Q,
                   pval.Q = print.pval.Q,
                   subgroup.name = print.subgroup.name),
      #
      type = list(study = type),
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
      test.subgroup = FALSE,
      #
      details = details)
  #
  forest_args <- c(forest_args_general, forest_args_specific)
  #
  fdat$class[fdat$class== "metamiss"] <- "metabin"
  classes <- unique(fdat$class)
  #
  text.details <- catmeth_metabind(x, common, random, prediction,
                                   print.tau2 = print.tau2,
                                   print.tau2.ci = print.tau2.ci,
                                   print.tau = print.tau,
                                   print.tau.ci = print.tau.ci,
                                   backtransf = x$backtransf,
                                   digits = digits,
                                   pscale = pscale,
                                   irscale = irscale, irunit = irunit,
                                   forest = TRUE)
  #
  forest_args$x$.text.details.methods <- text.details
  #
  res <- do.call("forest.meta", c(forest_args, args))
  #
  res$data <- fdat
  #
  return(invisible(res))
}
