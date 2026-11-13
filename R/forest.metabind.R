#' Forest plot to display the results of several meta-analyses
#'
#' @description
#' Draws a forest plot combining common effect estimates, random effects
#' estimates, and prediction intervals from a metabind object.
#'
#' @param x An object of class \code{\link{metabind}}.
#' @param layout Layout of the forest plot; see \code{\link{forest.meta}}.
#' @param common,random,prediction Logical values selecting common effect
#'   estimates, random effects estimates, and prediction intervals.
#' @param overall A logical value indicating whether overall summaries should
#'   be plotted.
#' @param studlab A variable containing labels for individual results,
#'   evaluated in the combined results data frame or supplied as an external
#'   vector. Available with or without subgroup results.
#' @param subgroup For meta-analyses without subgroup results, a variable
#'   defining display groups, evaluated in the combined results data frame or
#'   supplied as an external vector. For example, \code{subgroup = type} groups
#'   results by common effect model, random effects model, or prediction
#'   interval. Unavailable when the included meta-analyses contain subgroup
#'   results.
#' @param overall.hetstat Ignored when meta-analyses have no subgroup results.
#'   Otherwise, a logical value indicating whether to display a shared overall
#'   heterogeneity row. The default \code{NULL} determines this automatically
#'   from agreement of the overall results and heterogeneity.
#' @param hetstat Ignored when meta-analyses have no subgroup results; use
#'   \code{print} to display heterogeneity information in columns.
#'   Otherwise, a logical value controlling separate heterogeneity rows,
#'   defaulting to \code{TRUE}. Setting \code{FALSE} suppresses all such rows,
#'   including the shared overall heterogeneity row.
#' @param backtransf A logical indicating whether results should be
#'   back-transformed. Different summary measures require back-transformation.
#' @param test A named list controlling statistical test rows.
#' @param label A named list controlling column and test labels. Logical values
#'   or character vectors control the x-axis tick labels instead.
#'   See Details for defaults specific to this function and
#'   \code{\link{meta-formatting}} for information on the elements.
#' @param col,print,type,text,digits Named lists controlling colours, displayed
#'   information, graphical types, model text. Scalar values are also accepted
#'   for \code{type} and \code{digits}.
#'   See Details for defaults specific to this function and
#'   \code{\link{meta-formatting}} for information on the elements.
#' @param leftcols,rightcols Columns printed to the left or right of the
#'   forest plot; see \code{\link{forest.meta}}. Explicit column selections
#'   take precedence over automatic column selection.
#' @param leftlabs,rightlabs Column headings. Use \code{NA} to retain a
#'   default heading.
#' @param lab.NA Label for missing values.
#' @param big.mark Character used as thousands separator.
#' @param pscale A numeric giving scaling factor for printing of
#'   single event probabilities or risk differences
#'   (see \code{\link{metaprop}}).
#' @param irscale A numeric defining a scaling factor for printing of
#'   single incidence rates or incidence rate differences
#'   (see \code{\link{metarate}}).
#' @param irunit A character specifying the time unit used to
#'   calculate rates, e.g., person-years (see \code{\link{metarate}}).
#' @param scientific.pval,zero.pval,JAMA.pval Logical values controlling
#'   the formatting of P values (see \code{\link{forest.meta}}).
#' @param details A logical specifying whether method details are printed.
#'   Defaults to \code{TRUE} when the included meta-analyses contain subgroup
#'   results, and \code{FALSE} otherwise.
#' @param \dots Additional arguments passed to \code{\link{forest.meta}}.
#'
#' @details
#' When the included meta-analyses contain subgroup results, the forest plot
#' groups these results by meta-analysis. Overall summaries are shared when
#' the datasets and selected overall results agree; otherwise they are shown
#' separately for each meta-analysis. Shared heterogeneity information
#' additionally requires agreement of the heterogeneity values and settings.
#'
#' Confidence and prediction intervals share a column. Heterogeneity information
#' is shown once per included meta-analysis, with tau and tau-squared
#' information also shown on its prediction interval rows when random effects
#' results are selected. Tau and tau-squared information is not shown on common
#' effect rows.
#'
#' Back-transformations are applied separately using the settings stored for
#' each meta-analysis, allowing results with different transformations
#' to be displayed together.
#'
#' \subsection{Formatting defaults}{
#' Unless listed below, defaults for supported formatting elements are
#' described in \code{\link{meta-formatting}}.
#'
#' \tabular{ll}{
#' \bold{Element} \tab \bold{Default} \cr
#' \code{print$stat}, \code{print$pval} \tab Automatic. \cr
#' \code{print$I2}, \code{print$pval.Q} \tab See below. \cr
#' \code{print$I2.ci}, \code{print$Q} \tab \code{FALSE}. \cr
#' \code{print$subgroup.name} \tab \code{FALSE}. \cr
#' \code{type$common}, \code{type$random} \tab \code{NULL} (automatic). \cr
#' \code{digits$pval} \tab \code{max(gs("digits.pval") - 2, 2)} \cr
#' \code{digits$pval.Q} \tab \code{max(gs("digits.pval.Q") - 2, 2)} \cr
#' }
#' Fill and border colours are determined automatically from the layout.
#' They can be overridden using the following elements of the named list
#' \code{col}:
#' \itemize{
#' \item \code{square}, \code{square.lines}: square fill and border;
#' \item \code{circle}, \code{circle.lines}: circle fill and border;
#' \item \code{diamond}, \code{diamond.lines}: diamond fill and border;
#' \item \code{diamond.common}, \code{diamond.lines.common}: common effect
#'   diamond fill and border;
#' \item \code{diamond.random}, \code{diamond.lines.random}: random effects
#'   diamond fill and border;
#' \item \code{predict}, \code{predict.lines}: prediction interval fill
#'   and border.
#' }
#' 
#' The named list \code{print} controls which columns are shown.
#' Columns \code{func}, \code{sm}, \code{method}, and \code{method.interval}
#' are selected automatically when their values differ between the displayed
#' results. Heterogeneity-method columns \code{method.tau},
#' \code{method.tau.ci}, and \code{method.I2} are selected similarly when
#' the corresponding heterogeneity statistics are shown.
#' Statistics and P values are shown when available, except in subgroup plots.
#' 
#' Without subgroup results, \code{print$I2} defaults to \code{FALSE},
#' and \code{print$pval.Q} follows \code{print$Q}.
#' With subgroup results, I-squared is shown by default. The global setting
#' \code{settings.meta(print.I2 = FALSE)} suppresses it unless
#' \code{print$I2} is set explicitly. Element \code{print$pval.Q}
#' defaults to \code{TRUE} for subgroup plots.
#' }
#'
#' @return
#' Invisibly returns a list with information about the forest plot,
#' including the selected columns (\code{leftcols}, \code{rightcols}),
#' their headings (\code{leftlabs}, \code{rightlabs}), formatted results,
#' and graphical settings (see \code{\link{forest.meta}}).
#' Element \code{data} contains the dataset used to create the forest plot.
#'
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @seealso \code{\link{metabind}}, \code{\link{print.metabind}},
#'   \code{\link{forest.meta}}
#'
#' @keywords hplot
#'
#' @examples
#' ma <- metagen(1:3, 3:1, sm = "md")
#' #
#' mb <- metabind(ma, update(ma, method.random.ci = "HK"))
#' #
#' forest(mb)
#' #
#' forest(mb, subgroup = type)
#'
#' @method forest metabind
#' @export

forest.metabind <- function(x,
                            layout = gs("layout"),
                            #
                            common = !is.null(x$common),
                            random = !is.null(x$random),
                            prediction = !is.null(x$predict),
                            overall = common || random || prediction,
                            #
                            studlab = NULL, subgroup = NULL,
                            #
                            overall.hetstat = NULL,
                            hetstat = !is.null(x$subgroup),
                            #
                            backtransf = x$backtransf,
                            #
                            test = NULL,
                            label,
                            col = NULL, print = NULL, type = NULL,
                            text = NULL, digits = NULL,
                            #
                            leftcols = NULL, rightcols = NULL,
                            leftlabs = NULL, rightlabs = NULL,
                            #
                            lab.NA =
                              if (is.null(x$subgroup)) gs("lab.NA") else "",
                            #
                            big.mark = gs("big.mark"),
                            #
                            pscale = NULL, irscale = NULL, irunit = NULL,
                            #
                            scientific.pval = gs("scientific.pval"),
                            zero.pval = gs("zero.pval"),
                            JAMA.pval = gs("JAMA.pval"),
                            #
                            details = FALSE, ...) {
  
  #
  #
  # (1) Check and set arguments
  #
  #
  
  chkclass(x, "metabind")
  x <- updateversion(x)
  #
  with.subgroups <- !is.null(x$subgroup)
  #
  layout <- setchar(layout, gs("layout4forest"), length = 1)
  args <- list(...)
  missing.backtransf <- missing(backtransf)
  #
  chklogical(common)
  chklogical(random)
  chklogical(prediction)
  chklogical(hetstat)
  #
  missing.overall.hetstat <- missing(overall.hetstat)
  chklogical(overall.hetstat, NULL.ok = TRUE)
  #
  chklogical(backtransf)
  #
  selections <- list(common = common, random = random, predict = prediction)
  supplied <- c(common = !missing(common), random = !missing(random),
                predict = !missing(prediction))
  for (i in names(selections)) {
    if (is.null(x[[i]])) {
      if (supplied[[i]] && selections[[i]])
        warning("Argument '", if (i == "predict") "prediction" else i,
                "' set to FALSE as no corresponding results are available.",
                call. = FALSE)
      selections[[i]] <- FALSE
    }
    if (!selections[[i]])
      x[[i]] <- NULL
  }
  common <- selections$common
  random <- selections$random
  prediction <- selections$predict
  #
  if (!(common || random || prediction)) {
    warning("Nothing to show in forest plot.", call. = FALSE)
    return(invisible(NULL))
  }
  #
  chklogical(overall)
  #
  if (!with.subgroups) {
    if (!missing.overall.hetstat && !is.null(overall.hetstat))
      warning("Input to argument 'overall.hetstat' ignored for meta-analyses ",
              "without subgroups.",
              call. = FALSE)
    #
    if (hetstat)
      warning("Input to argument 'hetstat' ignored for meta-analyses ",
              "without subgroups; use 'print' to display ",
              "heterogeneity information in columns.",
              call. = FALSE)
    #
    if (!missing(overall) && overall)
      warning("Input to argument 'overall' ignored for meta-analyses ",
              "without subgroups.",
              call. = FALSE)
    #
    overall <- overall.hetstat <- hetstat <- FALSE
  }
  #
  chkchar(lab.NA)
  chkchar(big.mark, length = 1)
  #
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  chklogical(details)
  
  #
  #
  # (2) Process named lists
  #
  #
  
  # Colours used in forest plot
  #
  col <- as_named_list(col)
  col.list <- merge_legacy(col, args, col_defaults(x, "metabind"))
  col <- col.list$res
  args <- col.list$args
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
  col <- replaceNULL_list(col, list(
    square = colors[4],
    circle = colors[4],
    diamond = colors[2],
    diamond.lines = colors[3],
    predict = colors[7],
    predict.lines = colors[8]))
  #
  col <- replaceNULL_list(col, list(
    square.lines = col[["square"]],
    circle.lines = col[["circle"]],
    diamond.common = col[["diamond"]],
    diamond.random = col[["diamond"]],
    diamond.lines.common = col[["diamond.lines"]],
    diamond.lines.random = col[["diamond.lines"]]))
  #
  # Columns shown in forest plot
  #
  print.defaults <- print_defaults(x, "metabind")
  #
  if (with.subgroups)
    print.defaults["I2"] <- list(gs("print.I2"))
  #
  print.list <- merge_legacy(print, args, print.defaults)
  print <- print.list$res
  args <- print.list$args
  #
  chklogical_list(print,
                  c("type", "func", "sm",
                    "backtransf", "args.backtransf",
                    "method", "method.interval", "method.tau",
                    "pval", "stat",
                    "level.ma", "level.predict", "level.hetstat",
                    "tau2", "tau", "I2",
                    "Q", "pval.Q", "Q.b", "pval.Q.b",
                    "subgroup.name"))
  #
  print.defaults <- list(
    tau2 = random,
    tau = random,
    Q.b = with.subgroups,
    pval.Q.b = with.subgroups,
    subgroup.name = FALSE)
  #
  if (with.subgroups)
    print.defaults <- c(print.defaults, list(I2 = TRUE, pval.Q = TRUE))
  else {
    print.defaults[["I2"]] <- FALSE
    print.defaults[["pval.Q"]] <- print[["Q"]]
  }
  #
  print <- replaceNULL_list(print, print.defaults)
  #
  if (!with.subgroups)
    print$pval.Q <- print[["pval.Q"]] && print[["Q"]]
  #
  print$tau2.ci <- print[["tau2"]] && print[["tau2.ci"]]
  print$tau.ci <- print[["tau"]] && print[["tau.ci"]]
  print$I2.ci <- print[["I2"]] && print[["I2.ci"]]
  print$H.ci <- print[["H"]] && print[["H.ci"]]
  #
  # Digits used to round results
  #
  digits <- as_named_list(digits)
  digits.list <- merge_legacy(digits, args, digits_defaults("metabind"))
  digits <- digits.list$res
  args <- digits.list$args
  #
  chknumeric_list(digits, min = list(pval = 1, pval.Q = 1))
  #
  # Statistical tests
  #
  test.defaults <- test_defaults(x, common, random, TRUE)
  test.defaults$subgroup <- FALSE
  test.defaults$subgroup.common <- FALSE
  test.defaults$subgroup.random <- FALSE
  test.defaults["effect.subgroup"] <- list(NULL)
  #
  test.list <- merge_legacy(
    test, args, test.defaults,
    inherit = list(overall.common = "overall", overall.random = "overall",
                   subgroup.common = "subgroup", subgroup.random = "subgroup",
                   effect.subgroup.common = "effect.subgroup",
                   effect.subgroup.random = "effect.subgroup"))
  test <- test.list$res
  args <- test.list$args
  #
  test <- replaceNULL_list(test, list(overall = TRUE, subgroup = FALSE))
  test <- replaceNULL_list(
    test,
    list(overall.common = test[["overall"]],
         overall.random = test[["overall"]],
         subgroup.common = test[["subgroup"]],
         subgroup.random = test[["subgroup"]],
         effect.subgroup.common = test[["overall.common"]],
         effect.subgroup.random = test[["overall.random"]]))
  #
  chklogical_list(test, NULL.ok = "effect.subgroup",
                  length = list(effect.subgroup = 0,
                                effect.subgroup.common = 0,
                                effect.subgroup.random = 0))
  #
  test <- replaceNULL_list(
    test,
    list(effect.subgroup = any(test[["effect.subgroup.common"]]) ||
           any(test[["effect.subgroup.random"]])))
  #
  if (any(unlist(test[c("subgroup", "subgroup.common", "subgroup.random")])))
    warning("Between-component tests are not available for metabind objects.",
            call. = FALSE)
  #
  test$subgroup <- test$subgroup.common <- test$subgroup.random <- FALSE
  #
  # Column and test labels
  #
  missing.label <- missing(label)
  if (missing.label)
    label <- TRUE
  label <- deprecated(label, missing.label, args, "labels",
                      replaceNULL(args[["warn.deprecated"]], gs("warn.deprecated")),
                      name = "label")
  args$labels <- NULL
  if (is.list(label))
    label.axis <- TRUE
  else {
    label.axis <- label
    label <- NULL
  }
  label.list <- merge_legacy(label, args, label_defaults(x))
  label <- label.list$res
  args <- label.list$args
  chkchar(label[["pval"]], length = 1, name = "label.pval")
  #
  # Plot types
  #
  type <- as_named_list(type)
  type.list <- merge_legacy(type, args, type_defaults(class = "metabind"))
  types <- type.list$res
  args <- type.list$args
  #
  # Model and subgroup text
  #
  text.list <- merge_legacy(
    text, args, text_defaults(x),
    inherit = list(common.w = "common", random.w = "random", predict.w = "predict"))
  text <- text.list$res
  args <- text.list$args
  
  #
  #
  # (3) General labels and transformation settings
  #
  #
  
  ids <- unique(c(x$common$id, x$random$id, x$predict$id))
  meta <- x$overall[x$overall$id %in% ids, , drop = FALSE]
  sm <- unique(meta$sm)
  same.sm <- length(sm) == 1
  #
  pscale <- replaceNULL(pscale, makeunique(meta$pscale, 1))
  irscale <- replaceNULL(irscale, makeunique(meta$irscale, 1))
  irunit <- replaceNULL(irunit, makeunique(meta$irunit, "units"))
  #
  chknumeric_strict(pscale, min = 1, length = 1)
  chknumeric_strict(irscale, min = 1, length = 1)
  chkchar(irunit, length = 1)
  #
  if (!same.sm && !backtransf) {
    if (!missing.backtransf)
      warning("Argument 'backtransf' set to TRUE as ",
              "meta-analyses have different summary measures.", call. = FALSE)
    backtransf <- TRUE
  }
  if (same.sm) {
    fbt <- makeunique(meta$func.backtransf, NULL)
    sm.lab <- smlab(sm, backtransf = backtransf, pscale = pscale,
                    irscale = irscale, func.backtransf = fbt, forest = TRUE)
    #
    if (is_relative_effect(sm) || identical(fbt, "exp"))
      lab.TE <- paste0(gs("log.prefix"), sm)
    else
      lab.TE <- sm
  }
  else {
    fbt <- NULL
    sm.lab <- switch(x$sm.family, relative = "Effect",
                     correlation = "Correlation", mean = "Mean",
                     proportion = if (pscale == 1) "Proportion" else "Events",
                     rate = if (irscale == 1) "Rate" else "Events", "Effect")
    lab.TE <- "TE"
  }
  if (sm.lab == "")
    sm.lab <- "Effect"
  if (lab.TE == "")
    lab.TE <- "TE"
  
  #
  #
  # (4) Prepare selected result data before formatting
  #
  #
  
  dat <- as.data.frame(x, subgroup = FALSE, drop_all_NA = FALSE)
  #
  if (with.subgroups) {
    if (!missing(subgroup))
      stop("Argument 'subgroup' is unavailable for subgroup meta-analyses.",
           call. = FALSE)
    #
    prepared <- prepare_metabind_subgroups(x)
    #
    dat <- prepared$dat %>%
      rename_with(~ sub("\\.pooled$", "", .x), ends_with(".pooled"))
    #
    subdat <- prepared$subdat
    #
    # Parent heterogeneity belongs to summary rows, not source subgroup rows
    #
    for (i in c("Q", "df.Q", "pval.Q", "I2", "lower.I2", "upper.I2",
                "H", "lower.H", "upper.H"))
      subdat[[i]] <- NA_real_
    #
    for (i in c("tau2", "lower.tau2", "upper.tau2", "se.tau2",
                "tau", "lower.tau", "upper.tau"))
      subdat[[i]] <- rep(list(NA_real_), nrow(subdat))
    #
    plotdat <- subdat
    #
    if (prepared$shared.overall && !prepared$shared.heterogeneity)
      hetstat <- FALSE
    #
    overall.hetstat <- prepared$shared.overall &&
      prepared$shared.heterogeneity &&
      hetstat &&
      replaceNULL(overall.hetstat, TRUE)
  }
  else
    plotdat <- dat
  #
  shared.overall <- with.subgroups && prepared$shared.overall
  omit.test.cols <- with.subgroups && (shared.overall || !prepared$same.data)
  if (with.subgroups && missing(details))
    details <- TRUE
  #
  for (i in c("TE", "seTE", "statistic", "pval", "df", "null.effect")) {
    if (is.null(plotdat[[i]]))
      plotdat[[i]] <- rep(NA_real_, nrow(plotdat))
  }
  #
  if (is.null(plotdat$func.backtransf))
    plotdat$func.backtransf <- rep(NA_character_, nrow(plotdat))
  #
  if (is.null(plotdat$args.backtransf))
    plotdat$args.backtransf <- rep(list(NULL), nrow(plotdat))
  #
  studlab <- catch(enquo(studlab), plotdat)
  #
  if (with.subgroups) {
    sort.subgroup <- replaceNULL(args$sort.subgroup, gs("sort.subgroup"))
    chklogical(sort.subgroup)
    args$sort.subgroup <- NULL
    subgroup.levels <- as.character(unique(plotdat$id))
    if (sort.subgroup)
      subgroup.levels <- sort(subgroup.levels)
    subgroup <- factor(plotdat$id, levels = subgroup.levels)
  }
  else
    subgroup <- catch(enquo(subgroup), plotdat)
  #
  if (is.null(studlab)) {
    if (with.subgroups && sum(c(common, random, prediction)) > 1)
      studlab <- paste0(plotdat$subgroup.levels, " (",
                        c("common", "random", "prediction")[
                          match(plotdat$type,
                                c("Common effect model", "Random effects model",
                                  "Prediction interval"))], ")")
    else if (with.subgroups)
      studlab <- as.character(plotdat$subgroup.levels)
    else
      studlab <- plotdat$name
  }
  else if (length(studlab) == 1)
    studlab <- rep(studlab, nrow(plotdat))
  #
  chklength(studlab, nrow(plotdat), name = "studlab")
  #
  if (!is.null(subgroup)) {
    if (length(subgroup) == 1)
      subgroup <- rep(subgroup, nrow(plotdat))
    chklength(subgroup, nrow(plotdat), name = "subgroup")
  }
  #
  if (is.null(subgroup))
    group <- rep(1L, nrow(plotdat))
  else if (with.subgroups)
    group <- as.integer(subgroup)
  else
    group <- match(subgroup, unique(subgroup))
  # Format source and summary rows together to keep scales and headings aligned
  if (with.subgroups) {
    n.display <- nrow(plotdat)
    plotdat <- bind_rows(plotdat, dat)
    studlab <- c(studlab, dat$text)
    if (prepared$shared.overall)
      group <- c(group, rep(1L, nrow(dat)))
    else
      group <- c(group, match(dat$id, levels(subgroup)))
  }
  #
  varies <- function(value, rows = rep(TRUE, nrow(plotdat))) {
    any(vapply(split(value[rows], group[rows]),
               function(z) length(unique(z[!is.na(z)])) > 1,
               logical(1)))
  }
  #
  is.predict <- plotdat$type == "Prediction interval"
  is.common <- plotdat$type == "Common effect model"
  #
  priority <-
    match(plotdat$type,
          c("Random effects model", "Common effect model",
            "Prediction interval"))
  #
  order.het <- order(priority, seq_len(nrow(plotdat)))
  #
  rows.het <- rep(FALSE, nrow(plotdat))
  rows.het[order.het[!duplicated(plotdat$id[order.het])]] <- TRUE
  #
  rows.tau <- rows.het |
    (is.predict & plotdat$id %in% plotdat$id[plotdat$type == "Random effects model"])
  if (with.subgroups) {
    summary.rows <- n.display + seq_len(nrow(dat))
    rows.het[] <- rows.tau[] <- FALSE
    order.het <- summary.rows[order(priority[summary.rows])]
    rows.het[order.het[!duplicated(plotdat$id[order.het])]] <- TRUE
    rows.tau <- rows.het | (seq_len(nrow(plotdat)) %in% summary.rows &
                              is.predict & plotdat$id %in%
                              plotdat$id[summary.rows[plotdat$type[summary.rows] == "Random effects model"]])
  }
  rows.tau[is.common] <- FALSE
  #
  if (shared.overall && !prepared$shared.heterogeneity)
    rows.het[] <- rows.tau[] <- FALSE
  #
  print$type <- replaceNULL(print$type,
                            !with.subgroups && varies(plotdat$type) &&
                              !identical(studlab, plotdat$type))
  #
  for (i in c("func", "sm", "method", "method.interval")) {
    if (shared.overall && i %in% c("method", "method.interval"))
      different <- any(vapply(split(plotdat[[i]], plotdat$type),
                              function(z) length(unique(z[!is.na(z)])) > 1,
                              logical(1)))
    else
      different <- varies(plotdat[[i]])
    #
    print[[i]] <- replaceNULL(print[[i]], different)
  }
  #
  has.pval <- any(!is.na(plotdat$pval))
  print$pval <- replaceNULL(print$pval, has.pval && !with.subgroups) && has.pval
  print$stat <- replaceNULL(print$stat, has.pval && !with.subgroups) &&
    any(!is.na(plotdat$statistic))
  #
  print$backtransf <- replaceNULL(
    print$backtransf,
    any(!is.na(plotdat$func.backtransf) & nzchar(plotdat$func.backtransf)))
  print$args.backtransf <- print$backtransf &&
    replaceNULL(print$args.backtransf, any(lengths(plotdat$args.backtransf) > 0))
  #
  print$subgroup.name <- replaceNULL(print$subgroup.name, FALSE)
  #
  levels.ci <- unique(round(plotdat$level.ma[!is.predict], 4))
  levels.pi <- unique(round(plotdat$level.predict[is.predict], 4))
  #
  generic.interval <- length(levels.ci) > 1 || length(levels.pi) > 1
  #
  if (length(levels.ci) > 0 && length(levels.pi) > 0) {
    if (generic.interval)
      ci.lab <- "CI/PI"
    else if (levels.ci == levels.pi)
      ci.lab <- percent_label(levels.ci, type = "CI/PI")
    else
      ci.lab <- paste0(percent_label(levels.ci), "/",
                       percent_label(levels.pi, type = "PI"))
  }
  else if (length(levels.pi) > 0) {
    if (generic.interval)
      ci.lab <- "PI"
    else
      ci.lab <- percent_label(levels.pi, type = "PI")
  }
  else {
    if (generic.interval)
      ci.lab <- "CI"
    else
      ci.lab <- percent_label(levels.ci)
  }
  #
  plotdat$level <- ifelse(is.predict, plotdat$level.predict, plotdat$level.ma)
  #
  show.level.ci <- generic.interval && replaceNULL(print$level.ma, TRUE)
  show.level.pi <- generic.interval && replaceNULL(print$level.predict, TRUE)
  show.level <- (any(!is.predict) && show.level.ci) ||
    (any(is.predict) && show.level.pi)
  #
  is.t <-
    plotdat$method.interval[!is.predict] %in% c("HK", "KR", "CR0", "CR1", "CR2")
  #
  if (all(is.t))
    stat.lab <- "t"
  else if (all(!is.t))
    stat.lab <- "z"
  else
    stat.lab <- "z|t"
  #
  show.df <- !omit.test.cols && any(!is.na(plotdat$df) & is.finite(plotdat$df))
  #
  plotdat$TE.orig <- plotdat$TE
  plotdat$seTE.orig <- plotdat$seTE
  plotdat$lower.orig <- plotdat$lower
  plotdat$upper.orig <- plotdat$upper
  #
  row_setting <- function(column, i, arguments = FALSE) {
    if (is.null(column))
      return(NULL)
    #
    value <- column[[i]]
    #
    if (is.atomic(value) && length(value) == 1 && is.na(value))
      return(NULL)
    #
    # Flat named argument columns must retain their argument names.
    #
    if (arguments && !is.null(value) && !is.list(value))
      column[i]
    else
      value
  }
  #
  # Resolve common forward-transformation metadata
  #
  transf.settings <- lapply(seq_len(nrow(meta)), function(i) {
    func <- row_setting(meta$func.transf, i)
    #
    if (identical(func, ""))
      func <- NULL
    #
    args <- row_setting(meta$args.transf, i, arguments = TRUE)
    #
    if (is.null(func) || length(args) == 0)
      args <- NULL
    #
    list(func = func, args = args)
  })
  #
  same.transf <- all(vapply(
    transf.settings,
    function(setting) identical(setting, transf.settings[[1]]),
    logical(1)))
  #
  common.transf <- if (same.transf) transf.settings[[1]] else NULL
  #
  # Only exponential transformations can share the log-axis shortcut
  #
  uses.exp <-
    vapply(seq_len(nrow(plotdat)),
           function(i) {
             fun <- row_setting(plotdat$func.backtransf, i)
             if (is.null(fun))
               is_relative_effect(plotdat$sm[i])
             else
               identical(fun, "exp") || identical(fun, exp)
           },
           logical(1))
  #
  needs.scaling <-
    ((is_prop(plotdat$sm) | plotdat$sm == "RD") & pscale != 1) |
    ((is_rate(plotdat$sm) | plotdat$sm == "IRD") & irscale != 1)
  #
  relative.axis <- backtransf && all(uses.exp) && !any(needs.scaling)
  #
  if (backtransf && !relative.axis) {
    for (i in seq_len(nrow(plotdat))) {
      for (j in c("TE", "lower", "upper")) {
        plotdat[[j]][i] <-
          backtransf(plotdat[[j]][i], plotdat$sm[i],
                     n = plotdat$n.harmonic.mean[[i]],
                     time = plotdat$t.harmonic.mean[[i]],
                     func = row_setting(plotdat$func.backtransf, i),
                     args = row_setting(plotdat$args.backtransf, i,
                                        arguments = TRUE))
        if (is_prop(plotdat$sm[i]) || plotdat$sm[i] == "RD")
          plotdat[[j]][i] <- pscale * plotdat[[j]][i]
        else if (is_rate(plotdat$sm[i]) || plotdat$sm[i] == "IRD")
          plotdat[[j]][i] <- irscale * plotdat[[j]][i]
      }
      #
      # Switch lower and upper limits
      #
      if (!is.na(plotdat$lower[i]) && !is.na(plotdat$upper[i]) &&
          plotdat$lower[i] > plotdat$upper[i]) {
        tmp.l <- plotdat$lower[i]
        plotdat$lower[i] <- plotdat$upper[i]
        plotdat$upper[i] <- tmp.l
      }
    }
  }
  #
  ref.original <- !backtransf && same.sm && is_single(sm)
  #
  ref.original <- !backtransf && same.sm && is_single(sm)
  #
  if (ref.original && !same.transf)
    stop("Cannot display results on a common analysis scale because ",
         "the selected meta-analyses use different forward-transformation ",
         "functions or arguments. Use 'backtransf = TRUE' or select ",
         "meta-analyses with the same transformation.",
         call. = FALSE)
  #
  null.ref <- plotdat$null.effect
  null.plot <- plotdat$null.effect
  #
  for (i in seq_len(nrow(plotdat))) {
    if (is.na(null.plot[i]))
      next
    #
    original <-
      plotdat$func[i] %in% c("metacor", "metamean", "metaprop", "metarate")
    #
    if (ref.original && !original)
      null.ref[i] <-
        backtransf(null.ref[i], plotdat$sm[i],
                   n = plotdat$n.harmonic.mean[[i]],
                   time = plotdat$t.harmonic.mean[[i]],
                   func = row_setting(plotdat$func.backtransf, i),
                   args = row_setting(plotdat$args.backtransf, i,
                                      arguments = TRUE))
    #
    if (!backtransf && original)
      null.plot[i] <-
        transf(null.plot[i], plotdat$sm[i],
               n = plotdat$n.harmonic.mean[[i]],
               time = plotdat$t.harmonic.mean[[i]],
               func = row_setting(plotdat$func.transf, i),
               args = row_setting(plotdat$args.transf, i, arguments = TRUE))
    else if (backtransf && !original)
      null.plot[i] <-
        backtransf(null.plot[i], plotdat$sm[i],
                   n = plotdat$n.harmonic.mean[[i]],
                   time = plotdat$t.harmonic.mean[[i]],
                   func = row_setting(plotdat$func.backtransf, i),
                   args = row_setting(plotdat$args.backtransf, i,
                                      arguments = TRUE))
    #
    if (backtransf && (is_prop(plotdat$sm[i]) || plotdat$sm[i] == "RD"))
      null.plot[i] <- pscale * null.plot[i]
    else if (backtransf && (is_rate(plotdat$sm[i]) || plotdat$sm[i] == "IRD"))
      null.plot[i] <- irscale * null.plot[i]
  }
  
  #
  #
  # (5) Format integrated result and heterogeneity columns
  #
  #
  
  if (with.subgroups)
    summary.numeric <- plotdat[summary.rows, , drop = FALSE]
  #
  plotdat$null.effect <- formatN(null.plot, digits = digits$digits,
                                 lab.NA = lab.NA, big.mark = big.mark)
  plotdat$func.backtransf[is.na(plotdat$func.backtransf) |
                            plotdat$func.backtransf == ""] <- lab.NA
  plotdat$args.backtransf <- vapply(plotdat$args.backtransf, format_args,
                                    character(1), lab.NA = lab.NA)
  plotdat$k <- as.character(plotdat$k)
  plotdat$level <- percent_label(plotdat$level, type = "")
  plotdat$level[(!is.predict & !show.level.ci) |
                  (is.predict & !show.level.pi)] <- ""
  plotdat$df <-
    formatN(plotdat$df, digits = 2, lab.NA = lab.NA, big.mark = big.mark)
  #
  if (all(is.na(plotdat$seTE)))
    print$se <- FALSE
  #
  show.adhoc <- print$method.interval &&
    any(plotdat$method.interval %in% c("HK", "HK-PR")) &&
    any(!is.na(plotdat$adhoc.hakn) & plotdat$adhoc.hakn != "")
  #
  plotdat$adhoc.hakn[is.na(plotdat$adhoc.hakn)] <- ""
  #
  het.levels <- unique(round(plotdat$level.hetstat[rows.het], 4))
  #
  if (length(het.levels) == 1)
    het.lab <- percent_label(het.levels)
  else
    het.lab <- "CI"
  #
  show.het.level <-
    replaceNULL(print$level.hetstat, length(het.levels) > 1) &&
    (print$tau2.ci || print$tau.ci || print$I2.ci || print$H.ci)
  #
  plotdat$level.hetstat <- percent_label(plotdat$level.hetstat, type = "")
  #
  for (i in c("tau2", "tau", "I2", "H")) {
    lower <- plotdat[[paste0("lower.", i)]]
    upper <- plotdat[[paste0("upper.", i)]]
    #
    if (i %in% c("I2", "H"))
      prec <- digits$I2
    else
      prec <- digits[[i]]
    #
    if (i %in% c("tau2", "tau"))
      plotdat[[paste0(i, ".ci")]] <-
      format_metabind_variance(lower, upper = upper, digits = prec,
                               lab.NA = lab.NA, big.mark = big.mark)
    else
      plotdat[[paste0(i, ".ci")]] <-
      formatCI(lower, upper, digits = prec, lab.NA = lab.NA, na.lab = TRUE,
               scale = if (i == "I2") 100 else 1,
               suffix = if (i == "I2") "%" else "", big.mark = big.mark)
    #
    if (i == "I2") {
      plotdat[[i]] <-
        percent_label(plotdat[[i]], "", prec, formatN, lab.NA = lab.NA)
    }
    else if (i == "H")
      plotdat[[i]] <- formatN(plotdat[[i]], digits = prec, lab.NA = lab.NA,
                              big.mark = big.mark)
    else
      plotdat[[i]] <- format_metabind_variance(
        plotdat[[i]], digits = prec, lab.NA = lab.NA, big.mark = big.mark)
  }
  #
  plotdat$Q <- formatN(plotdat$Q, digits = digits$Q, lab.NA = lab.NA,
                       big.mark = big.mark)
  #
  plotdat$pval.Q <- formatPT(plotdat$pval.Q, digits = digits$pval.Q,
                             lab.NA = lab.NA, lab = FALSE, trim = "left",
                             scientific = scientific.pval, zero = zero.pval,
                             JAMA = JAMA.pval, big.mark = big.mark)
  #
  subgroup.stat.lab <- "Q"
  #
  if (with.subgroups) {
    for (i in c("Q.b", "df.Q.b", "pval.Q.b")) {
      if (is.null(plotdat[[i]]))
        plotdat[[i]] <- rep(NA_real_, nrow(plotdat))
    }
    #
    if (is.null(plotdat[["test.Q.b"]]))
      plotdat$test.Q.b <- rep(NA_character_, nrow(plotdat))
    #
    first.level <- plotdat$subgroup.levels[match(plotdat$id, plotdat$id)]
    #
    keep <- !is.predict & !is.na(plotdat$subgroup.levels) &
      plotdat$subgroup.levels == first.level
    #
    subgroup.stat.lab <- subgroup_test_heading(plotdat[["test.Q.b"]][keep])
    #
    plotdat$Q.b <- formatN(plotdat$Q.b, digits = digits$Q, lab.NA = lab.NA,
                           big.mark = big.mark)
    #
    plotdat$df.Q.b <- vapply(plotdat$df.Q.b, format_df, character(1),
                             lab.NA = lab.NA, big.mark = big.mark)
    #
    plotdat$pval.Q.b <- formatPT(plotdat$pval.Q.b, digits = digits$pval.Q,
                                 lab.NA = lab.NA, lab = FALSE, trim = "left",
                                 scientific = scientific.pval, zero = zero.pval,
                                 JAMA = JAMA.pval, big.mark = big.mark)
    #
    for (i in c("Q.b", "df.Q.b", "pval.Q.b", "test.Q.b"))
      plotdat[[i]][!keep] <- lab.NA
  }
  #
  het.cols <- c("tau2", "tau", "I2", "H", "tau2.ci", "tau.ci", "I2.ci", "H.ci",
                "Q", "df.Q", "pval.Q", "level.hetstat")
  #
  for (i in het.cols) {
    if (i %in% c("tau2", "tau", "tau2.ci", "tau.ci"))
      keep <- rows.tau
    else
      keep <- rows.het
    #
    plotdat[[i]][!keep] <- lab.NA
  }
  #
  for (i in c("method.tau", "method.tau.ci", "method.I2")) {
    if (is.null(plotdat[[i]]))
      plotdat[[i]] <- NA_character_
    #
    plotdat[[i]][is.na(plotdat[[i]])] <- lab.NA
    if (i == "method.tau")
      plotdat[[i]][plotdat[[i]] == ""] <- lab.NA
    #
    if (i %in% c("method.tau", "method.tau.ci"))
      keep <- rows.tau
    else
      keep <- rows.het
    #
    plotdat[[i]][!keep] <- lab.NA
  }
  #
  print.het.method <-
    c(method.tau = (print$tau2 || print$tau),
      method.tau.ci = (print$tau2.ci || print$tau.ci),
      method.I2 = print$I2)
  #
  for (i in names(print.het.method)) {
    print.het.method[[i]] <- print.het.method[[i]] &&
      replaceNULL(print[[i]],
                  replaceNULL(print.list$res$method,
                              varies(plotdat[[i]], rows.het)))
  }
  #
  if (with.subgroups) {
    summary.display <- plotdat[summary.rows, , drop = FALSE]
    plotdat <- plotdat[seq_len(n.display), , drop = FALSE]
    studlab <- studlab[seq_len(n.display)]
    is.predict <- is.predict[seq_len(n.display)]
    is.common <- is.common[seq_len(n.display)]
  }
  #
  type.study <- types$type
  #
  if (is.null(type.study))
    type.study <- ifelse(is.predict, "predict", "square")
  else if (length(type.study) == 1)
    type.study <- rep(type.study, nrow(plotdat))
  #
  chklength(type.study, nrow(plotdat), name = "type")
  #
  valid.types <- c("square", "diamond", "predict", "circle")
  type.study <- setchar(type.study, valid.types)
  #
  for (i in c("common", "random", "predict")) {
    if (!is.null(types[[i]])) {
      types[[i]] <- setchar(types[[i]], valid.types, length = 1)
      sel <- plotdat$type == switch(i, common = "Common effect model",
                                    random = "Random effects model",
                                    predict = "Prediction interval")
      type.study[sel] <- types[[i]]
    }
  }
  #
  plotdat$.type.study <- type.study
  col.square <- col.square.lines <- rep("", nrow(plotdat))
  #
  for (i in seq_len(nrow(plotdat))) {
    if (type.study[i] == "diamond") {
      if (is.common[i])
        type.i <- "common"
      else
        type.i <- "random"
      #
      col.square[i] <- col[[paste0("diamond.", type.i)]]
      col.square.lines[i] <- col[[paste0("diamond.lines.", type.i)]]
    }
    else {
      if (type.study[i] == "predict")
        type.i <- "predict"
      else if (type.study[i] == "circle")
        type.i <- "circle"
      else
        type.i <- "square"
      #
      col.square[i] <- col[[type.i]]
      col.square.lines[i] <- col[[paste0(type.i, ".lines")]]
    }
  }
  #
  col$square <- col.square
  col$square.lines <- col.square.lines
  
  #
  #
  # (6) Select columns and their headings
  #
  #
  
  info.cols <- c(if (print$studlab) "studlab",
                 if (print$type) "type",
                 if (print$func) "func",
                 if (print$method) "method",
                 if (print$method.interval) "method.interval",
                 if (show.adhoc) "adhoc.hakn",
                 if (print$k) "k")
  #
  for (i in c("tau2", "tau")) {
    if (print[[i]])
      info.cols <- c(info.cols, i)
    if (print[[paste0(i, ".ci")]])
      info.cols <- c(info.cols, paste0(i, ".ci"))
  }
  #
  info.cols <- c(info.cols,
                 if (print.het.method[["method.tau"]]) "method.tau",
                 if (print.het.method[["method.tau.ci"]]) "method.tau.ci")
  #
  for (i in c("I2", "H")) {
    if (print[[i]])
      info.cols <- c(info.cols, i)
    if (print[[paste0(i, ".ci")]])
      info.cols <- c(info.cols, paste0(i, ".ci"))
    if (i == "I2" && print.het.method[["method.I2"]])
      info.cols <- c(info.cols, "method.I2")
  }
  #
  if (print$Q)
    info.cols <- c(info.cols, "Q", "df.Q")
  #
  if (print$pval.Q)
    info.cols <- c(info.cols, "pval.Q")
  #
  info.cols <- c(info.cols, if (show.het.level) "level.hetstat")
  # Shared heterogeneity is printed in the native overall heterogeneity row.
  if (with.subgroups && prepared$shared.overall)
    info.cols <- setdiff(info.cols, c(het.cols, names(print.het.method)))
  #
  if (with.subgroups && (common || random)) {
    info.cols <-
      c(info.cols,
        if (print$Q.b) "Q.b",
        if (print$Q.b) "df.Q.b",
        if (print$pval.Q.b) "pval.Q.b")
  }
  #
  transformation.cols <- c(if (print$backtransf) "func.backtransf",
                           if (print$args.backtransf) "args.backtransf")
  measure.cols <- c(if (print$sm) "sm", if (show.level) "level")
  if (layout %in% c("RevMan5", "JAMA", "BMJ"))
    effect.cols <- "effect.ci"
  else
    effect.cols <- c("effect", "ci")
  test.cols <- c(if (print$se) "seTE", if (print$stat) "statistic",
                 if (show.df) "df", if (print$pval) "pval",
                 if (print$null.effect) "null.effect")
  #
  if (is.null(leftcols))
    leftcols <- c(info.cols, if (!layout %in% c("BMJ", "meta8"))
      c(measure.cols, transformation.cols, effect.cols, test.cols))
  #
  if (is.null(rightcols)) {
    if (layout %in% c("BMJ", "meta8"))
      rightcols <- c(measure.cols, effect.cols, transformation.cols, test.cols)
    else
      rightcols <- FALSE
  }
  #
  if (with.subgroups)
    studlab.heading <- "Subgroup"
  else
    studlab.heading <- "Meta-Analysis"
  #
  labels <- c(studlab = studlab.heading,
              type = "Type",
              func = "Function",
              func.backtransf = "Back-\ntransformation",
              args.backtransf = "Arguments",
              method = "Method",
              method.interval = "Interval\nmethod",
              adhoc.hakn = "Ad hoc\ncorrection",
              k = "Number of\nStudies",
              sm = "Effect\nmeasure",
              level = "Level",
              effect = sm.lab,
              ci = ci.lab,
              effect.ci = paste0(sm.lab, "\n(", ci.lab, ")"),
              seTE = paste0("SE(", lab.TE, ")"),
              statistic = stat.lab,
              pval = label[["pval"]],
              df = "df",
              null.effect = "Null\neffect",
              tau2 = NA,
              tau = NA,
              I2 = NA,
              H = NA,
              Q = NA,
              df.Q = "df",
              pval.Q = paste(label[["pval"]], "Q-test", sep = "\n"),
              Q.b = subgroup.stat.lab,
              test.Q.b = "Test",
              df.Q.b = "df",
              pval.Q.b = paste("Interaction", label[["pval"]], sep = "\n"),
              tau2.ci = het.lab,
              tau.ci = het.lab,
              I2.ci = het.lab,
              H.ci = het.lab,
              method.tau = "Tau estimator",
              method.tau.ci = "Tau CI\nmethod",
              method.I2 = "I2 method",
              level.hetstat = "Het.\nlevel")
  #
  if (is.null(leftlabs))
    leftlabs <- rep(NA_character_, length(leftcols))
  else
    chklength(leftlabs, length(leftcols), name = "leftlabs")
  #
  if (is.null(rightlabs))
    rightlabs <- rep(NA_character_, length(rightcols))
  else
    chklength(rightlabs, length(rightcols), name = "rightlabs")
  #
  for (i in names(labels)) {
    leftlabs <- newlab(leftcols, leftlabs, i, labels[[i]])
    rightlabs <- newlab(rightcols, rightlabs, i, labels[[i]])
  }
  
  #
  #
  # (7) Create a display-only forest skeleton
  #
  #
  
  se.plot <- plotdat$seTE
  se.plot[is.na(se.plot) | se.plot <= 0] <- 1
  #
  if (relative.axis && same.sm)
    sm.plot <- sm
  else if (relative.axis)
    sm.plot <- "Effect"
  else if (!backtransf && same.sm)
    sm.plot <- sm
  else
    sm.plot <- ""
  #
  ma <- metagen(plotdat$TE, se.plot,
                lower = plotdat$lower, upper = plotdat$upper,
                studlab = studlab, data = plotdat,
                sm = sm.plot,
                common = FALSE, random = FALSE, prediction = FALSE,
                overall = FALSE, transf = TRUE,
                backtransf = relative.axis,
                func.backtransf =
                  if (relative.axis && !is_relative_effect(sm.plot))
                    "exp" else NULL,
                subgroup = subgroup, print.subgroup.name = print$subgroup.name,
                label.left = makeunique(plotdat$label.left, ""),
                label.right = makeunique(plotdat$label.right, ""),
                method.tau = "DL", method.tau.ci = "", warn = FALSE)
  #
  ma$TE[is.predict] <- NA
  ma$seTE <- plotdat$seTE
  #
  ma$statistic <- plotdat$statistic
  ma$pval <- plotdat$pval
  ma$null.effect <- NA_real_
  #
  if (ref.original) {
    ma$func.transf <- common.transf$func
    ma$args.transf <- common.transf$args
  }
  #
  # Define PAS and IRS as fallbacks - calculated in forest.meta()
  #
  if (ref.original && sm == "PFT")
    ma$n.harmonic.mean <- makeunique(plotdat$n.harmonic.mean, val = NULL)
  #
  if (ref.original && sm == "IRFT")
    ma$t.harmonic.mean <- makeunique(plotdat$t.harmonic.mean, val = NULL)
  #
  if (!"ref" %in% names(args)) {
    if (ref.original)
      args$ref <- makeunique(null.ref)
    else
      args$ref <- makeunique(null.plot)
  }
  #
  if (sm.plot == "" && !"smlab" %in% names(args))
    args$smlab <- sm.lab
  #
  forest_args_general <-
    list(
      layout = layout,
      #
      leftcols = leftcols,
      leftlabs = leftlabs,
      rightcols = rightcols,
      rightlabs = rightlabs,
      #
      hetstat = hetstat,
      overall.hetstat = overall.hetstat,
      #
      lab.NA = lab.NA,
      #
      pooled.totals = FALSE,
      weight.study = "same",
      #
      col = col,
      #
      print = print[c("I2", "I2.ci", "tau2", "tau2.ci", "tau", "tau.ci",
                      "H", "Q", "pval.Q", "subgroup.name")],
      #
      type = list(study = type.study),
      #
      digits = digits,
      #
      scientific.pval = scientific.pval,
      zero.pval = zero.pval,
      JAMA.pval = JAMA.pval,
      big.mark = big.mark,
      #
      details = details)
  #
  forest_args_specific <- list(x = ma)
  label.defaults <- text.defaults <- list()
  test.plot <- modifyList(
    test,
    list(overall = FALSE, overall.common = FALSE, overall.random = FALSE,
         subgroup = FALSE, subgroup.common = FALSE, subgroup.random = FALSE,
         effect.subgroup = FALSE, effect.subgroup.common = FALSE,
         effect.subgroup.random = FALSE))
  #
  if (with.subgroups) {
    mapped <- map_metabind_subgroups(ma, x, summary.numeric, summary.display,
                                     prepared$shared.overall,
                                     common, random, prediction, overall,
                                     shared.heterogeneity = prepared$shared.heterogeneity)
    #
    forest_args_specific <- c(list(x = mapped$ma), mapped$args)
    label.defaults <- mapped$args$label
    text.defaults <- mapped$args$text
    #
    forest_args_specific[c("label", "text", "test.overall", "test.subgroup",
                           "test.effect.subgroup",
                           "test.effect.subgroup.common",
                           "test.effect.subgroup.random")] <- NULL
    #
    forest_args_specific$sort.subgroup <- FALSE
    #
    forest_args_specific$subgroup.hetstat <- !shared.overall && hetstat
    #
    if (shared.overall) {
      has.null <- !is.na(mapped$ma$null.effect)
      test.plot$overall <- test[["overall"]] && has.null
      test.plot$overall.common <- common && test[["overall.common"]] && overall && has.null
      test.plot$overall.random <- random && test[["overall.random"]] && overall && has.null
      test.common <- test.plot[["overall.common"]]
      test.random <- test.plot[["overall.random"]]
    }
    else {
      for (i in c("effect.subgroup", "effect.subgroup.common",
                  "effect.subgroup.random")) {
        if (length(test[[i]]) != 1) {
          chklength(test[[i]], length(ma$subgroup.levels),
                    name = paste0("test.", i),
                    internal = FALSE)
          test[[i]] <- test[[i]][match(as.character(ma$subgroup.levels),
                                       unique(as.character(plotdat$id)))]
        }
      }
      #
      test.plot$effect.subgroup.common <-
        mapped$args$test.effect.subgroup.common &
        test[["effect.subgroup.common"]] & common
      #
      test.plot$effect.subgroup.random <-
        mapped$args$test.effect.subgroup.random &
        test[["effect.subgroup.random"]] & random
      #
      test.common <- any(test.plot[["effect.subgroup.common"]])
      test.random <- any(test.plot[["effect.subgroup.random"]])
      #
      test.plot$effect.subgroup <- test.common || test.random
    }
    #
    if (prepared$shared.overall) {
      forest_args_general$hetstat <- hetstat
      forest_args_general$overall.hetstat <- overall.hetstat
    }
  }
  else if (!is.null(subgroup)) {
    forest_args_specific <-
      c(forest_args_specific,
        list(common.subgroup = FALSE, random.subgroup = FALSE,
             prediction.subgroup = FALSE))
  }
  #
  if (!with.subgroups && is.null(subgroup))
    test.plot[c("subgroup", "subgroup.common", "subgroup.random",
                "effect.subgroup",
                "effect.subgroup.common", "effect.subgroup.random")] <- NULL
  forest_args_general$test <- test.plot
  #
  label.supplied <-
    names(label.list$miss)[!vapply(label.list$miss, isTRUE, logical(1))]
  #
  label.plot <-
    modifyList(label.defaults, label[label.supplied], keep.null = TRUE)
  #
  if (isTRUE(label.axis))
    forest_args_general$label <- label.plot
  else {
    forest_args_general$label <- label.axis
    args <- c(args, setNames(label.plot, paste0("label.", names(label.plot))))
  }
  #
  text.supplied <-
    names(text.list$miss)[!vapply(text.list$miss, isTRUE, logical(1))]
  #
  for (i in c("common", "random", "predict")) {
    if (i %in% text.supplied)
      text.supplied <- union(text.supplied, paste0(i, ".w"))
  }
  forest_args_general$text <-
    modifyList(text.defaults, text[text.supplied], keep.null = TRUE)
  #
  forest_args_general[names(forest_args_specific)] <- NULL
  forest_args <- c(forest_args_general, forest_args_specific)
  #
  if (details) {
    details.print <- print
    if (with.subgroups && (test.common || test.random))
      details.print$pval <- TRUE
    #
    forest_args$x$.text.details.methods <-
      catmeth_metabind(x,
                       #
                       common, random, prediction,
                       #
                       backtransf, details.print, digits,
                       #
                       pscale, irscale, irunit,
                       #
                       big.mark = big.mark,
                       #
                       forest = TRUE)
  }
  #
  res <- do.call("forest.meta", c(forest_args, args))
  #
  res$data <- plotdat
  #
  invisible(res)
}
