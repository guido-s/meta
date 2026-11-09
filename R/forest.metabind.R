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
#' @param overall A logical indicating whether overall summaries should
#'   be plotted. No additional pooled estimates are calculated for the
#'   component meta-analysis results.
#' @param studlab A variable containing labels for individual results,
#'   evaluated in the combined result data and then the calling environment.
#' @param subgroup A variable defining display groups, evaluated in the
#'   combined result data and then the calling environment. For example,
#'   \code{subgroup = type} groups results by result type.
#' @param hetstat,overall.hetstat Arguments passed to
#'   \code{\link{forest.meta}} for heterogeneity text. Heterogeneity text
#'   from the synthetic plotting object is disabled; use \code{print} to
#'   display the component meta-analyses' heterogeneity information.
#' @param backtransf A logical indicating whether results should be
#'   back-transformed. Different summary measures require back-transformation.
#' @param col,print,type,digits Named lists controlling colours, displayed
#'   information, graphical types, and rounding; see
#'   \code{\link{meta-formatting}}. Scalar values are also accepted for
#'   \code{type} and \code{digits}.
#' @param leftcols,rightcols Columns printed to the left or right of the
#'   forest plot; see \code{\link{forest.meta}}. Explicit column selections
#'   take precedence over automatic column selection.
#' @param leftlabs,rightlabs Column headings. Use \code{NA} to retain a
#'   default heading.
#' @param lab.NA Label for missing values.
#' @param scientific.pval,zero.pval,JAMA.pval Logical values controlling
#'   the formatting of P values.
#' @param big.mark Character used as thousands separator.
#' @param addrow.subgroups A logical indicating whether an empty row is
#'   printed between display groups.
#' @param details A logical specifying whether method details are printed.
#' @param \dots Additional arguments passed to \code{\link{forest.meta}}.
#'
#' @details
#' Automatic class, summary-measure, estimation-method, and interval-method
#' columns are shown when their values differ within a display group, or
#' across the selected results if no groups are specified. Set the
#' corresponding element of \code{print} to \code{TRUE} or \code{FALSE}
#' to override this behaviour. Element \code{method.tau} (or legacy argument
#' \code{print.method.tau}) controls the heterogeneity estimator column
#' independently of \code{print$method}. If P values are available, statistics and
#' P values are shown by default.
#'
#' Confidence and prediction intervals share a column. Its heading gives
#' their levels when each interval kind has one level, for example,
#' \code{95\% CI/PI} or \code{95\% CI/99\% PI}. Otherwise a generic
#' heading and a level column are used.
#'
#' Heterogeneity information is shown once per component meta-analysis:
#' on its first random effects row, falling back to its first common effect
#' row and then its first prediction interval row.
#'
#' @return Invisibly returns the result of \code{forest.meta()}, with the
#'   prepared result data in element \code{data}.
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' @seealso \code{\link{metabind}}, \code{\link{print.metabind}},
#'   \code{\link{forest.meta}}
#' @keywords hplot
#' @examples
#' ma <- metagen(1:3, 3:1, sm = "RR")
#' mb <- metabind(ma, update(ma, method.random.ci = "HK"))
#' forest(mb)
#' forest(mb, subgroup = type)
#'
#' @method forest metabind
#' @export

forest.metabind <- function(x,
                            layout = gs("layout"),
                            common = !is.null(x$common),
                            random = !is.null(x$random),
                            prediction = !is.null(x$predict),
                            overall = common || random,
                            studlab = NULL, subgroup = NULL,
                            hetstat = FALSE,
                            overall.hetstat = x$overall.hetstat,
                            backtransf = x$backtransf,
                            col = NULL, print = NULL, type = NULL,
                            digits = NULL,
                            leftcols = NULL, leftlabs = NULL,
                            rightcols = NULL, rightlabs = NULL,
                            lab.NA = gs("lab.NA"),
                            scientific.pval = gs("scientific.pval"),
                            zero.pval = gs("zero.pval"),
                            JAMA.pval = gs("JAMA.pval"),
                            big.mark = gs("big.mark"),
                            addrow.subgroups = NULL,
                            details = FALSE, ...) {
  
  #
  #
  # (1) Check and set arguments
  #
  #
  
  chkclass(x, "metabind")
  x <- updateversion(x)
  if (!is.null(x$subgroup))
    return(invisible(NULL))
  layout <- setchar(layout, gs("layout4forest"), length = 1)
  args <- list(...)
  missing.backtransf <- missing(backtransf)
  #
  chklogical(common)
  chklogical(random)
  chklogical(prediction)
  chklogical(backtransf)
  selections <- list(common = common, random = random, predict = prediction)
  supplied <- c(common = !missing(common), random = !missing(random),
                predict = !missing(prediction))
  for (nm in names(selections)) {
    if (is.null(x[[nm]])) {
      if (supplied[[nm]] && selections[[nm]])
        warning("Argument '", if (nm == "predict") "prediction" else nm,
                "' set to FALSE as no corresponding results are available.",
                call. = FALSE)
      selections[[nm]] <- FALSE
    }
    if (!selections[[nm]])
      x[[nm]] <- NULL
  }
  common <- selections$common
  random <- selections$random
  prediction <- selections$predict
  if (!(common || random || prediction)) {
    warning("Nothing to show in forest plot.", call. = FALSE)
    return(invisible(NULL))
  }
  # Component results are display rows, not observations to pool again.
  overall <- overall.hetstat <- FALSE
  if (!identical(hetstat, FALSE))
    warning("Argument 'hetstat' ignored; use 'print' to display ",
            "component heterogeneity information.", call. = FALSE)
  hetstat <- FALSE
  chkchar(lab.NA)
  chkchar(big.mark, length = 1)
  chklogical(scientific.pval)
  chklogical(zero.pval)
  chklogical(JAMA.pval)
  chklogical(details)
  #
  print.list <- merge_legacy(print, args, print_defaults(x, "metabind"))
  print <- print.list$res
  args <- print.list$args
  nullable <- c("type", "class", "sm", "method", "method.interval", "method.tau",
                "pval", "stat", "level.ma", "level.predict", "level.hetstat",
                "tau2", "tau", "pval.Q", "pval.Q.b", "subgroup.name")
  for (nm in names(print))
    chklogical(print[[nm]], name = paste0("print.", nm),
               NULL.ok = nm %in% nullable)
  print$tau2 <- replaceNULL(print$tau2, random)
  print$tau <- replaceNULL(print$tau, random)
  print$pval.Q <- replaceNULL(print$pval.Q, print$Q) && print$Q
  print$pval.Q.b <- replaceNULL(print$pval.Q.b, FALSE)
  print$subgroup.name <- replaceNULL(print$subgroup.name, FALSE)
  print$tau2.ci <- print$tau2 && print$tau2.ci
  print$tau.ci <- print$tau && print$tau.ci
  print$I2.ci <- print$I2 && print$I2.ci
  print$H.ci <- print$H && print$H.ci
  #
  digits <- as_named_list(digits)
  digits.list <- merge_legacy(digits, args, digits_defaults("metabind"))
  digitss <- digits.list$res
  args <- digits.list$args
  for (nm in names(digitss))
    chknumeric(digitss[[nm]], min = if (nm %in% c("pval", "pval.Q")) 1 else 0,
               length = 1, name = paste0("digits.", nm))
  digits <- digitss$digits
  #
  type <- as_named_list(type)
  type.list <- merge_legacy(type, args, type_defaults(class = "metabind"))
  types <- type.list$res
  args <- type.list$args
  col <- as_named_list(col)
  col.list <- merge_legacy(col, args, col_defaults(x, "metabind"))
  cols <- col.list$res
  args <- col.list$args
  colors <- layout_colors(layout)
  cols$square <- replaceNULL(cols$square, colors[4])
  cols$square.lines <- replaceNULL(cols$square.lines, cols$square)
  cols$circle <- replaceNULL(cols$circle, colors[4])
  cols$circle.lines <- replaceNULL(cols$circle.lines, cols$circle)
  cols$diamond <- replaceNULL(cols$diamond, colors[2])
  cols$diamond.common <- replaceNULL(cols$diamond.common, cols$diamond)
  cols$diamond.random <- replaceNULL(cols$diamond.random, cols$diamond)
  cols$diamond.lines <- replaceNULL(cols$diamond.lines, colors[3])
  cols$diamond.lines.common <-
    replaceNULL(cols$diamond.lines.common, cols$diamond.lines)
  cols$diamond.lines.random <-
    replaceNULL(cols$diamond.lines.random, cols$diamond.lines)
  cols$predict <- replaceNULL(cols$predict, colors[7])
  cols$predict.lines <- replaceNULL(cols$predict.lines, colors[8])
  
  #
  #
  # (2) General labels and transformation settings
  #
  #
  
  ids <- unique(c(x$common$id, x$random$id, x$predict$id))
  meta <- x$overall[x$overall$id %in% ids, , drop = FALSE]
  sm <- unique(meta$sm)
  same.sm <- length(sm) == 1
  pscale <- makeunique(meta$pscale, 1)
  irscale <- makeunique(meta$irscale, 1)
  irunit <- makeunique(meta$irunit, "units")
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
    lab.TE <- if (is_relative_effect(sm) || identical(fbt, "exp"))
      paste0(gs("log.prefix"), sm) else sm
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
  # (3) Prepare selected result data before formatting
  #
  #
  
  dat <- as.data.frame(x, drop_all_NA = FALSE)
  for (nm in c("TE", "seTE", "statistic", "pval", "df", "null.effect")) {
    if (is.null(dat[[nm]]))
      dat[[nm]] <- rep(NA_real_, nrow(dat))
  }
  #
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  studlab <- catch("studlab", mc, dat, sfsp)
  subgroup <- catch("subgroup", mc, dat, sfsp)
  #
  if (is.null(studlab))
    studlab <- dat$name
  else if (length(studlab) == 1)
    studlab <- rep(studlab, nrow(dat))
  chklength(studlab, nrow(dat), name = "studlab")
  if (!is.null(subgroup)) {
    if (length(subgroup) == 1)
      subgroup <- rep(subgroup, nrow(dat))
    chklength(subgroup, nrow(dat), name = "subgroup")
  }
  group <- if (is.null(subgroup)) rep(1L, nrow(dat)) else
    match(subgroup, unique(subgroup))
  varies <- function(value, rows = rep(TRUE, nrow(dat))) {
    any(vapply(split(value[rows], group[rows]),
               function(z) length(unique(z[!is.na(z)])) > 1,
               logical(1)))
  }
  is.predict <- dat$type == "Prediction interval"
  is.common <- dat$type == "Common effect model"
  priority <- match(dat$type, c("Random effects model", "Common effect model",
                               "Prediction interval"))
  order.het <- order(priority, seq_len(nrow(dat)))
  rows.het <- rep(FALSE, nrow(dat))
  rows.het[order.het[!duplicated(dat$id[order.het])]] <- TRUE
  #
  print$type <- replaceNULL(print$type,
                            varies(dat$type) && !identical(studlab, dat$type))
  for (nm in c("class", "sm", "method", "method.interval"))
    print[[nm]] <- replaceNULL(print[[nm]], varies(dat[[nm]]))
  has.pval <- any(!is.na(dat$pval))
  print$pval <- replaceNULL(print$pval, has.pval) && has.pval
  print$stat <- replaceNULL(print$stat, has.pval) && any(!is.na(dat$statistic))
  print$subgroup.name <- replaceNULL(print$subgroup.name, FALSE)
  addrow.subgroups <- replaceNULL(addrow.subgroups, !is.null(subgroup))
  chklogical(addrow.subgroups)
  #
  levels.ci <- unique(round(dat$level.ma[!is.predict], 4))
  levels.pi <- unique(round(dat$level.predict[is.predict], 4))
  generic.interval <- length(levels.ci) > 1 || length(levels.pi) > 1
  if (length(levels.ci) > 0 && length(levels.pi) > 0) {
    ci.lab <-
      if (generic.interval) "CI/PI"
      else if (levels.ci == levels.pi) percent_label(levels.ci, type = "CI/PI")
      else paste0(percent_label(levels.ci), "/",
                  percent_label(levels.pi, type = "PI"))
  }
  else if (length(levels.pi) > 0)
    ci.lab <- if (generic.interval) "PI" else
      percent_label(levels.pi, type = "PI")
  else
    ci.lab <- if (generic.interval) "CI" else percent_label(levels.ci)
  dat$level <- ifelse(is.predict, dat$level.predict, dat$level.ma)
  show.level.ci <- generic.interval && replaceNULL(print$level.ma, TRUE)
  show.level.pi <- generic.interval && replaceNULL(print$level.predict, TRUE)
  show.level <- (any(!is.predict) && show.level.ci) ||
    (any(is.predict) && show.level.pi)
  #
  is.t <- dat$method.interval[!is.predict] %in% c("HK", "KR", "CR0", "CR1", "CR2")
  stat.lab <- if (all(is.t)) "t" else if (all(!is.t)) "z" else "z|t"
  show.df <- any(!is.na(dat$df) & is.finite(dat$df))
  #
  # Keep relative effects on the log axis; other families use natural coordinates.
  relative.axis <- backtransf &&
    (all(is_relative_effect(sm)) || (same.sm && identical(fbt, "exp")))
  dat$TE.orig <- dat$TE
  dat$seTE.orig <- dat$seTE
  dat$lower.orig <- dat$lower
  dat$upper.orig <- dat$upper
  row_setting <- function(column, i, arguments = FALSE) {
    if (is.null(column))
      return(NULL)
    value <- column[[i]]
    if (is.atomic(value) && length(value) == 1 && is.na(value))
      return(NULL)
    # Flat named argument columns must retain their argument names.
    if (arguments && !is.null(value) && !is.list(value))
      column[i]
    else
      value
  }
  if (backtransf && !relative.axis) {
    for (i in seq_len(nrow(dat))) {
      for (nm in c("TE", "lower", "upper")) {
        dat[[nm]][i] <- backtransf(
          dat[[nm]][i], dat$sm[i],
          n = dat$n.harmonic.mean[[i]], time = dat$t.harmonic.mean[[i]],
          func = row_setting(dat$func.backtransf, i),
          args = row_setting(dat$args.backtransf, i, arguments = TRUE))
        if (is_prop(dat$sm[i]) || dat$sm[i] == "RD")
          dat[[nm]][i] <- pscale * dat[[nm]][i]
        else if (is_rate(dat$sm[i]) || dat$sm[i] == "IRD")
          dat[[nm]][i] <- irscale * dat[[nm]][i]
      }
    }
  }
  null.plot <- dat$null.effect
  for (i in seq_len(nrow(dat))) {
    if (is.na(null.plot[i]))
      next
    original <- dat$class[i] %in% c("metacor", "metamean", "metaprop", "metarate")
    if (!backtransf && original)
      null.plot[i] <- transf(null.plot[i], dat$sm[i],
                             n = dat$n.harmonic.mean[[i]],
                             time = dat$t.harmonic.mean[[i]],
                             func = row_setting(dat$func.transf, i),
                             args = row_setting(dat$args.transf, i, arguments = TRUE))
    else if (backtransf && !original)
      null.plot[i] <- backtransf(null.plot[i], dat$sm[i],
                                 n = dat$n.harmonic.mean[[i]],
                                 time = dat$t.harmonic.mean[[i]],
                                 func = row_setting(dat$func.backtransf, i),
                                 args = row_setting(dat$args.backtransf, i, arguments = TRUE))
    if (backtransf && (is_prop(dat$sm[i]) || dat$sm[i] == "RD"))
      null.plot[i] <- pscale * null.plot[i]
    else if (backtransf && (is_rate(dat$sm[i]) || dat$sm[i] == "IRD"))
      null.plot[i] <- irscale * null.plot[i]
  }
  
  #
  #
  # (4) Format integrated result and heterogeneity columns
  #
  #
  
  dat$null.effect <- formatN(null.plot, digits = digits, lab.NA = lab.NA,
                             big.mark = big.mark)
  dat$k <- as.character(dat$k)
  dat$level <- percent_label(dat$level, type = "")
  dat$level[(!is.predict & !show.level.ci) | (is.predict & !show.level.pi)] <- ""
  dat$df <- formatN(dat$df, digits = 2, lab.NA = lab.NA, big.mark = big.mark)
  if (all(is.na(dat$seTE)))
    print$se <- FALSE
  show.adhoc <- print$method.interval &&
    any(dat$method.interval %in% c("HK", "HK-PR")) &&
    any(!is.na(dat$adhoc.hakn) & dat$adhoc.hakn != "")
  dat$adhoc.hakn[is.na(dat$adhoc.hakn)] <- ""
  het.levels <- unique(round(dat$level.hetstat[rows.het], 4))
  het.lab <- if (length(het.levels) == 1) percent_label(het.levels) else "CI"
  show.het.level <- replaceNULL(print$level.hetstat, length(het.levels) > 1) &&
    (print$tau2.ci || print$tau.ci || print$I2.ci || print$H.ci)
  dat$level.hetstat <- percent_label(dat$level.hetstat, type = "")
  for (nm in c("tau2", "tau", "I2", "H")) {
    lower <- dat[[paste0("lower.", nm)]]
    upper <- dat[[paste0("upper.", nm)]]
    prec <- if (nm %in% c("I2", "H")) digitss$I2 else digitss[[nm]]
    dat[[paste0(nm, ".ci")]] <- formatCI(
      lower, upper, digits = prec, lab.NA = lab.NA, na.lab = TRUE,
      scale = if (nm == "I2") 100 else 1,
      suffix = if (nm == "I2") "%" else "", big.mark = big.mark)
    if (nm == "I2")
      dat[[nm]] <- percent_label(dat[[nm]], "", prec, formatN, lab.NA = lab.NA)
    else if (nm == "H")
      dat[[nm]] <- formatN(dat[[nm]], digits = prec, lab.NA = lab.NA,
                           big.mark = big.mark)
    else
      dat[[nm]] <- formatPT(dat[[nm]], digits = prec, lab.NA = lab.NA,
                            big.mark = big.mark, lab = FALSE, trim = "left")
  }
  dat$Q <- formatN(dat$Q, digits = digitss$Q, lab.NA = lab.NA,
                    big.mark = big.mark)
  dat$pval.Q <- formatPT(dat$pval.Q, digits = digitss$pval.Q,
                         lab.NA = lab.NA, lab = FALSE, trim = "left",
                         scientific = scientific.pval, zero = zero.pval,
                         JAMA = JAMA.pval, big.mark = big.mark)
  het.cols <- c("tau2", "tau", "I2", "H", "tau2.ci", "tau.ci", "I2.ci", "H.ci",
                "Q", "df.Q", "pval.Q", "level.hetstat")
  for (nm in het.cols)
    dat[[nm]][!rows.het] <- lab.NA
  #
  for (nm in c("method.tau", "method.tau.ci", "method.I2")) {
    if (is.null(dat[[nm]]))
      dat[[nm]] <- NA_character_
    dat[[nm]][is.na(dat[[nm]])] <- lab.NA
    dat[[nm]][!rows.het] <- lab.NA
  }
  print.het.method <- c(
    method.tau = (print$tau2 || print$tau),
    method.tau.ci = (print$tau2.ci || print$tau.ci),
    method.I2 = print$I2)
  for (nm in names(print.het.method))
    print.het.method[[nm]] <- print.het.method[[nm]] &&
      replaceNULL(print[[nm]],
                  replaceNULL(print.list$res$method, varies(dat[[nm]], rows.het)))
  #
  type.study <- types$type
  if (is.null(type.study))
    type.study <- ifelse(is.predict, "predict", "square")
  else if (length(type.study) == 1)
    type.study <- rep(type.study, nrow(dat))
  chklength(type.study, nrow(dat), name = "type")
  valid.types <- c("square", "diamond", "predict", "circle")
  type.study <- setchar(type.study, valid.types)
  for (nm in c("common", "random", "predict")) {
    if (!is.null(types[[nm]])) {
      types[[nm]] <- setchar(types[[nm]], valid.types, length = 1)
      sel <- dat$type == switch(nm, common = "Common effect model",
                                random = "Random effects model",
                                predict = "Prediction interval")
      type.study[sel] <- types[[nm]]
    }
  }
  dat$.type.study <- type.study
  cols.square <- cols.square.lines <- rep("", nrow(dat))
  for (i in seq_len(nrow(dat))) {
    if (type.study[i] == "diamond") {
      nm <- if (is.common[i]) "common" else "random"
      cols.square[i] <- cols[[paste0("diamond.", nm)]]
      cols.square.lines[i] <- cols[[paste0("diamond.lines.", nm)]]
    }
    else {
      nm <- if (type.study[i] == "predict") "predict" else
        if (type.study[i] == "circle") "circle" else "square"
      cols.square[i] <- cols[[nm]]
      cols.square.lines[i] <- cols[[paste0(nm, ".lines")]]
    }
  }
  
  #
  #
  # (5) Select columns and their headings
  #
  #
  
  info.cols <- c(if (print$studlab) "studlab", if (print$type) "type",
                 if (print$class) "class", if (print$method) "method",
                 if (print$method.interval) "method.interval",
                 if (show.adhoc) "adhoc.hakn", if (print$k) "k")
  #
  for (nm in c("tau2", "tau", "I2", "H")) {
    if (print[[nm]])
      info.cols <- c(info.cols, nm)
    if (print[[paste0(nm, ".ci")]])
      info.cols <- c(info.cols, paste0(nm, ".ci"))
  }
  if (print$Q)
    info.cols <- c(info.cols, "Q", "df.Q")
  #
  if (print$pval.Q)
    info.cols <- c(info.cols, "pval.Q")
  info.cols <- c(info.cols, names(print.het.method)[print.het.method],
                 if (show.het.level) "level.hetstat")
  result.cols <- c(if (print$sm) "sm", if (show.level) "level",
                   if (layout %in% c("RevMan5", "JAMA", "BMJ")) "effect.ci"
                   else c("effect", "ci"),
                   if (print$se) "seTE", if (print$stat) "statistic",
                   if (show.df) "df",
                   if (print$pval) "pval", if (print$null.effect) "null.effect")
  if (is.null(leftcols))
    leftcols <- c(info.cols, if (!layout %in% c("BMJ", "meta8")) result.cols)
  if (is.null(rightcols))
    rightcols <- if (layout %in% c("BMJ", "meta8")) result.cols else FALSE
  #
  labels <- c(studlab = "Meta-Analysis", type = "Type", class = "Function",
              method = "Method", method.interval = "Interval\nmethod",
              adhoc.hakn = "Ad hoc\ncorrection", k = "Number of\nStudies",
              sm = "Effect\nmeasure", level = "Level", effect = sm.lab,
              ci = ci.lab, effect.ci = paste0(sm.lab, "\n(", ci.lab, ")"),
              seTE = paste0("SE(", lab.TE, ")"), statistic = stat.lab,
              pval = "P value", df = "df", null.effect = "Null\neffect",
              tau2 = NA, tau = NA, I2 = NA,
              H = NA, Q = NA, df.Q = "df", pval.Q = "P value\nQ-test",
              tau2.ci = het.lab, tau.ci = het.lab, I2.ci = het.lab, H.ci = het.lab,
              method.tau = "Tau estimator", method.tau.ci = "Tau CI\nmethod",
              method.I2 = "I2 method", level.hetstat = "Het.\nlevel")
  if (is.null(leftlabs))
    leftlabs <- rep(NA_character_, length(leftcols))
  else
    chklength(leftlabs, length(leftcols), name = "leftlabs")
  if (is.null(rightlabs))
    rightlabs <- rep(NA_character_, length(rightcols))
  else
    chklength(rightlabs, length(rightcols), name = "rightlabs")
  for (nm in names(labels)) {
    leftlabs <- newlab(leftcols, leftlabs, nm, labels[[nm]])
    rightlabs <- newlab(rightcols, rightlabs, nm, labels[[nm]])
  }
  
  #
  #
  # (6) Create a display-only forest skeleton
  #
  #
  
  se.plot <- dat$seTE
  se.plot[is.na(se.plot) | se.plot <= 0] <- 1
  ma <- metagen(dat$TE, se.plot, lower = dat$lower, upper = dat$upper,
                studlab = studlab, data = dat,
                sm = if (relative.axis && same.sm) sm else sm.lab,
                common = FALSE, random = FALSE, prediction = FALSE,
                overall = FALSE, transf = TRUE,
                backtransf = relative.axis,
                func.backtransf =
                  if (relative.axis && (!same.sm || !is_relative_effect(sm)))
                    "exp" else NULL,
                subgroup = subgroup, print.subgroup.name = print$subgroup.name,
                label.left = makeunique(dat$label.left, ""),
                label.right = makeunique(dat$label.right, ""),
                method.tau = "DL", method.tau.ci = "", warn = FALSE)
  ma$TE[is.predict] <- NA
  #
  ma$statistic <- dat$statistic
  ma$pval <- dat$pval
  ma$null.effect <- NA_real_
  #
  if (!"ref" %in% names(args))
    args$ref <- makeunique(null.plot)
  #
  forest_args_general <- list(
    layout = layout,
    leftcols = leftcols, leftlabs = leftlabs,
    rightcols = rightcols, rightlabs = rightlabs,
    hetstat = hetstat, overall.hetstat = overall.hetstat,
    lab.NA = lab.NA, pooled.totals = FALSE, weight.study = "same",
    addrow.subgroups = addrow.subgroups,
    col = list(square = cols.square, square.lines = cols.square.lines,
               circle = cols$circle, circle.lines = cols$circle.lines,
               diamond = cols$diamond, diamond.common = cols$diamond.common,
               diamond.random = cols$diamond.random,
               diamond.lines = cols$diamond.lines,
               diamond.lines.common = cols$diamond.lines.common,
               diamond.lines.random = cols$diamond.lines.random,
               predict = cols$predict, predict.lines = cols$predict.lines),
    print = list(I2 = print$I2, I2.ci = print$I2.ci,
                 tau2 = print$tau2, tau2.ci = print$tau2.ci,
                 tau = print$tau, tau.ci = print$tau.ci,
                 H = print$H, Q = print$Q, pval.Q = print$pval.Q,
                 subgroup.name = print$subgroup.name),
    type = list(study = type.study), digits = digitss,
    scientific.pval = scientific.pval, zero.pval = zero.pval,
    JAMA.pval = JAMA.pval, big.mark = big.mark,
    test.subgroup = FALSE, test.overall = FALSE, details = details)
  #
  forest_args_specific <- list(x = ma)
  #
  if (!is.null(subgroup))
    forest_args_specific <- c(
      forest_args_specific,
      list(common.subgroup = FALSE, random.subgroup = FALSE,
           prediction.subgroup = FALSE))
  #
  forest_args <- c(forest_args_general, forest_args_specific)
  #
  if (details)
    forest_args$x$.text.details.methods <- catmeth_metabind(
      x, common, random, prediction,
      print.tau2 = print$tau2, print.tau2.ci = print$tau2.ci,
      print.tau = print$tau, print.tau.ci = print$tau.ci,
      print.I2 = print$I2, print.pval = print$pval,
      backtransf = backtransf, digits = digits,
      digits.tau2 = digitss$tau2, digits.tau = digitss$tau,
      pscale = pscale, irscale = irscale, irunit = irunit,
      big.mark = big.mark, forest = TRUE)
  #
  res <- do.call("forest.meta", c(forest_args, args))
  #
  res$data <- dat
  #
  invisible(res)
}
