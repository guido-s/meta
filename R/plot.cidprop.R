#' Plot density of prediction distribution highlighting areas of clinically
#' important benefit or harm
#'
#' @description
#' Plot density of prediction distribution highlighting areas of clinically
#' important benefit or harm
#'
#' @param x An object of class \code{cidprop}.
#' @param cid A numeric value or vector specifying clinically important
#'   differences (CID) / decision thresholds used to calculate expected
#'   proportions of clinically important benefit or harm (see Details).
#' @param cid.below.null A numeric value or vector specifying CID limits below
#'   the null effect (see Details).
#' @param cid.above.null A numeric value or vector specifying CID limits above
#'   the null effect (see Details).
#' @param label.cid A character string or vector specifying labels for
#'   clinically important differences. Must be of same length as argument
#'   \code{cid}.
#' @param label.cid.below.null A character string or vector specifying labels
#'   for clinically important differences below the null effect. Must be of
#'   same length as argument \code{cid.below.null} (or \code{cid}).
#' @param label.cid.above.null A character string or vector specifying labels
#'   for clinically important differences above the null effect. Must be of
#'   same length as argument \code{cid.above.null} (or \code{cid}).
#' @param small.values A character string specifying whether small
#'   treatment effects indicate a beneficial (\code{"desirable"}) or
#'   harmful (\code{"undesirable"}), can be abbreviated.
#' @param fill.cid.below.null Background colour(s) for CID areas below null
#'   effect.
#' @param fill.cid.above.null Background colour(s) for CID areas above null
#'   effect.
#' @param fill Background colour for area between decision thresholds.
#' @param legend A logical indicating whether to print a legend with
#'   expected proportions of beneficial, harmful, or not important effects.
#' @param studies A logical indicating whether to print estimates of individual
#'   studies.
#' @param random A logical indicating whether to show diamond of the random
#'   effects meta-analysis.
#' @param col.diamond The colour of the diamond representing the results
#'   for the random effects model.
#' @param col.diamond.lines The colour of the outer lines of the diamond
#'   representing the results of the random effects model.
#' @param prediction A logical indicating whether to show the prediction
#'   interval.
#' @param col.predict The colour of the prediction interval.
#' @param col.predict.lines The colour of the outer lines of the prediction
#'   interval.
#' @param big.mark A character used as thousands separator.
#' @param digits.cid Minimal number of significant digits for
#'   decision thresholds, see \code{\link{print.default}}.
#' @param digits.percent Minimal number of significant digits for
#'   expected proportions, printed as percentages, see
#'   \code{\link{print.default}}.
#' @param digits.xaxis Minimal number of significant digits for
#'   labels on x-axis, see \code{\link{print.default}}.
#' @param xlab Label on x-axis.
#' @param ylab Label on y-axis.
#' @param xlim Limits for x-axis.
#' @param ylim Limits for y-axis.
#' @param labels.x Predefined labels for tick marks on x-axis.
#' @param \dots Additional arguments (ignored)
#' 
#' @details
#' Arguments \code{cid}, \code{cid.below.null}, \code{cid.above.null},
#' \code{label.cid}, \code{label.cid.below.null}, \code{label.cid.above.null},
#' and \code{small.values} are identical to the main arguments of R function
#' \code{\link{cidprop}} which is called internally if any of these values has
#' been provided by the user.
#' 
#' R packages \bold{ggpubr} and \bold{gridExtra} must be installed in order to
#' add a legend to the plot with the CIDs, expected proportions of clinically
#' benefit or harm, and the area colours (due to using R functions
#' \code{ggarrange} and \code{tableGrob}). The data and colours shown in the
#' legend are stored in the attribute 'data.cid' of the returned ggplot object
#' (see Examples).
#' 
#' UTF-8 code for the less than or equal and greater than or equal signs are
#' used in the legend. Accordingly, graphic devices with full UTF-8 support are
#' required to save graphics, for example, \code{\link[grDevices]{cairo_pdf}}
#' instead of \code{\link[grDevices]{pdf}} from R package \bold{grDevices}.
#' 
#' @return
#' A ggplot object with additional class 'plot.cidprop'.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @seealso \code{\link{cidprop}}
#'
#' @examples
#' oldset <- settings.meta(digits.cid = 0)
#' 
#' m <- metagen(1:10 - 3, 1:10, sm = "MD")
#' 
#' pp1 <- cidprop(m, cid = 2)
#' pp1
#' plot(pp1, xlim = c(-4, 4))
#' \donttest{
#' pp2 <- cidprop(m, cid.below.null = 0.5, cid.above.null = 2)
#' pp2
#' plot(pp2, xlim = c(-4, 4))
#' 
#' pp3 <- cidprop(m, cid.below.null = 0.5, cid.above.null = 2,
#'   small.values = "u")
#' pp3
#' plot(pp3, xlim = c(-4, 4))
#' 
#' pp4 <- cidprop(m, cid = 1:2, label.cid = c("moderate", "large"))
#' pp4
#' plot(pp4, xlim = c(-4, 4))
#' 
#' pp5 <- cidprop(m, cid.below.null = -1.5, cid.above.null = 1:2,
#'   label.cid.below.null = "large",
#'   label.cid.above.null = c("moderate", "large"))
#' pp5
#' plpp5 <- plot(pp5, xlim = c(-4, 4))
#' plpp5
#' # Information on CIDs and colours
#' attr(plpp5, "data.cid")
#' }
#' \dontrun{
#' # R packages 'ggpubr' and 'gridExtra' must be available
#' if (requireNamespace("ggpubr", quietly = TRUE) &
#'     requireNamespace("gridExtra", quietly = TRUE)) {
#'   plot(pp1, xlim = c(-4, 4), legend = TRUE)
#' }
#' }
#' 
#' settings.meta(oldset)
#' 
#' @method plot cidprop
#' @export

plot.cidprop <- function(x,
                         cid = NULL,
                         cid.below.null = x$cid.below.null,
                         cid.above.null = x$cid.above.null,
                         #
                         label.cid = "",
                         label.cid.below.null = x$label.cid.below.null,
                         label.cid.above.null = x$label.cid.above.null,
                         #
                         small.values = x$small.values,
                         #
                         fill.cid.below.null = NULL,
                         fill.cid.above.null = NULL,
                         fill = "white",
                         #
                         legend = FALSE,
                         #
                         studies = TRUE,
                         #
                         random = TRUE,
                         col.diamond = gs("col.diamond"),
                         col.diamond.lines = gs("col.diamond.lines"),
                         #
                         prediction = TRUE,
                         col.predict = gs("col.predict"),
                         col.predict.lines = gs("col.predict.lines"),
                         #
                         big.mark = gs("big.mark"),
                         digits.cid = gs("digits.cid"), digits.percent = 1,
                         #
                         digits.xaxis = gs("digits.forest"),
                         #
                         xlab = NULL, ylab = NULL,
                         xlim = NULL, ylim = NULL,
                         labels.x = NULL,
                         #
                         ...) {
  
  #
  #
  # (1) Check arguments
  #
  #
  
  chkclass(x, "cidprop")
  #
  missing.xlim <- missing(xlim)
  #
  chklogical(legend)
  chklogical(studies)
  chklogical(random)
  chklogical(prediction)
  #
  chknumeric(digits.cid, min = 0, length = 1)
  chknumeric(digits.percent, min = 0, length = 1)
  chknumeric(digits.xaxis, min = 0, length = 1)
  #
  if (legend) {
    if (!is_installed_package("ggpubr", argument = "legend", value = TRUE,
                              stop = FALSE)) {
      warning("Please install R package 'ggpubr' in order to print a legend.",
      call. = FALSE)
      #
      legend <- FALSE
    }
  }
  
  #
  #
  # (2) Re-run cidprop() if settings changed
  #
  #
  
  avail.cid <- !missing(cid) & !is.null(cid) & !all(is.na(cid))
  avail.cid.below.null <- !missing(cid.below.null) &
    !is.null(cid.below.null) & !all(is.na(cid.below.null))
  avail.cid.above.null <- !missing(cid.above.null) &
    !is.null(cid.above.null) & !all(is.na(cid.above.null))
  #
  avail.label.cid <- !missing(label.cid)
  avail.label.cid.below.null <- !missing(label.cid.below.null)
  avail.label.cid.above.null <- !missing(label.cid.above.null)
  #
  avail.small.values <- !missing(small.values)
  #
  # Rerun cidprop()
  #
  if (avail.cid | avail.cid.below.null | avail.cid.above.null |
      avail.label.cid |
      avail.label.cid.below.null | avail.label.cid.above.null |
      avail.small.values) {
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
        warning("Argument 'cid.below.null' ignored as ", 
                "argument 'cid' is provided.",
                call. = FALSE)
      else if (avail.cid.above.null)
        warning("Argument 'cid.above.null' ignored as ",
                "argument 'cid' is provided.",
                call. = FALSE)
      #
      x <- cidprop(x$x,
                   cid = cid, label.cid = label.cid,
                   small.values = small.values)
    }
    else {
      x <- cidprop(x$x,
                   cid.below.null = cid.below.null,
                   cid.above.null = cid.above.null,
                   label.cid.below.null = label.cid.below.null,
                   label.cid.above.null = label.cid.above.null,
                   small.values = small.values)
    }
  }
  
  #
  #
  # (3) Extract results from cidprop-object
  #
  #
    
  x.meta <- x$x
  #
  sm <- x.meta$sm
  smlab <- smlab(sm, x.meta$backtransf, x.meta$pscale, x.meta$irscale)
  #
  if (is.null(xlab))
    xlab <- xlab_meta(sm, backtransf = TRUE)
  #
  if (is.null(ylab))
    ylab <- "Density"
  #
  is_relative <- is_relative_effect(sm)
  #
  cid.below.null <- x$cid.below.null
  cid.above.null <- x$cid.above.null
  #
  if (is_relative) {
    cid.below.null.transf <- log(cid.below.null)
    cid.above.null.transf <- log(cid.above.null)
  }
  else {
    cid.below.null.transf <- cid.below.null
    cid.above.null.transf <- cid.above.null
  }
  #
  label.cid.below.null <- x$label.cid.below.null
  label.cid.above.null <- x$label.cid.above.null
  #
  svd <- x$small.values == "desirable"
  #
  prop.cid.below.null <- x$prop.cid.below.null
  prop.cid.above.null <- x$prop.cid.above.null
  prop.within.cid <- x$prop.within.cid
  #
  method.predict <- x.meta$method.predict[1]
  #
  TE <- x.meta$TE
  lower <- x.meta$lower
  upper <- x.meta$upper
  w.random <- x.meta$w.random
  #
  TE.random <- x.meta$TE.random[1]
  lower.random <- x.meta$lower.random[1]
  upper.random <- x.meta$upper.random[1]
  lower.predict <- x.meta$lower.predict[1]
  upper.predict <- x.meta$upper.predict[1]
  seTE.predict <- x.meta$seTE.predict[1]
  df.predict <- x.meta$df.predict[1]
  #
  avail.cid.below.null <- !all(is.na(cid.below.null))
  avail.cid.above.null <- !all(is.na(cid.above.null))
  
  #
  #
  # (4) Get / set colors for CID areas
  #
  #
  
  if (svd) {
    if (length(cid.below.null) == 1)
      fill.cid.below.null <- replaceNULL(fill.cid.below.null, "lightgreen")
    else
      fill.cid.below.null <-
        replaceNULL(fill.cid.below.null,
                    hcl.colors(length(cid.below.null) + 1, palette = "Greens",
                               alpha = 1)[seq_along(cid.below.null)])
    #
    if (length(cid.above.null) == 1)
      fill.cid.above.null <- replaceNULL(fill.cid.above.null, "lightpink1")
    else
      fill.cid.above.null <-
        replaceNULL(fill.cid.above.null,
                    hcl.colors(length(cid.above.null) + 1, palette = "Reds",
                               alpha = 1)[rev(seq_along(cid.above.null))])
  }
  else {
    if (length(cid.below.null) == 1)
      fill.cid.below.null <- replaceNULL(fill.cid.below.null, "lightpink1")
    else
      fill.cid.below.null <-
        replaceNULL(fill.cid.below.null,
                    hcl.colors(length(cid.below.null) + 1, palette = "Reds",
                               alpha = 1)[seq_along(cid.below.null)])
    #
    if (length(cid.above.null) == 1)
      fill.cid.above.null <- replaceNULL(fill.cid.above.null, "lightgreen")
    else
      fill.cid.above.null <-
        replaceNULL(fill.cid.above.null,
                    hcl.colors(length(cid.above.null) + 1, palette = "Greens",
                               alpha = 1)[rev(seq_along(cid.above.null))])
  }
  #
  if (length(cid.below.null) != length(fill.cid.below.null))
    stop("Different length for arguments 'cid.below.null' and ",
         "'fill.cid.below.null'.",
         call. = FALSE)
  #
  if (length(cid.above.null) != length(fill.cid.above.null))
    stop("Different length for arguments 'cid.above.null' and ",
         "'fill.cid.above.null'.",
         call. = FALSE)
  #
  if (length(fill) != 1)
    stop("Argument 'fill' must be a single colour",
         call. = FALSE)
  
  #
  #
  # (5) Data set with CID information
  #
  #
  
  dat.l <- dat.u <- dat.w <- NULL
  #
  if (avail.cid.below.null) {
    dat.l <-
      data.frame(Colour = fill.cid.below.null,
                 Threshold = cid.below.null, prop = prop.cid.below.null,
                 label = label.cid.below.null,
                 category =
                   if (svd) "Beneficial effect" else "Harmful effect",
                 sign = "\u2264 ")
    #
    max.cid.below.null <- max(cid.below.null, na.rm = TRUE)
  }
  #
  if (avail.cid.above.null) {
    dat.u <-
      data.frame(Colour = fill.cid.above.null,
                 Threshold = cid.above.null, prop = prop.cid.above.null,
                 label = label.cid.above.null,
                 category =
                   if (svd) "Harmful effect" else "Beneficial effect",
                 sign = "\u2265 ")
    #
    min.cid.above.null <- min(cid.above.null, na.rm = TRUE)
  }
  #
  if (prop.within.cid > 0) {
    dat.w <- data.frame(Colour = fill,
                        Threshold = NA,
                        prop = prop.within.cid,
                        label = "",
                        category = "Not important effect",
                        sign = "")
    #
    if (avail.cid.below.null & avail.cid.above.null) {
      within.cid <-
        formatN(c(max.cid.below.null, min.cid.above.null),
                digits = digits.cid, big.mark = big.mark)
      #
     within.cid <- paste(">", within.cid[1], "to", "<", within.cid[2])
    }
    else if (avail.cid.below.null) {
      within.cid <-
        formatN(max.cid.below.null, digits = digits.cid, big.mark = big.mark)
      #
      within.cid <- paste(">", within.cid)
    }
    else if (avail.cid.above.null) {
      within.cid <-
        formatN(min.cid.above.null, digits = digits.cid, big.mark = big.mark)
      #
      within.cid <- paste("<", within.cid)
    }
  }
  #
  dat.cid <- rbind(dat.l, dat.w, dat.u)
  #
  Threshold <- prop <- label <- category <- sign <- NULL
  #
  dat.cid %<>%
    mutate(Threshold = formatN(Threshold, digits = digits.cid,
                               big.mark = big.mark, text.NA = ""),
           Threshold = if_else(category != "Not important effect",
                               paste0(sign, Threshold), within.cid),
           prop = paste0(formatPT(100 * prop, digits = digits.percent), "%"),
           category =
             if_else(label == "", category, paste(category, label))) %>%
    column_to_rownames("category") %>%
    rename(Percent = prop) %>%
    select(-label, -sign)
  
  #
  #
  # (6) Data set for plot
  #
  #
  
  # Get rid of warnings
  #
  cid.category <- xval <- yval <- y <- NULL
  #
  vals <- c(lower, upper, lower.predict, upper.predict)
  #
  if (missing.xlim) {
    from <- min(vals, na.rm = TRUE)
    to <- max(vals, na.rm = TRUE)
  }
  else {
    chknumeric(xlim, length = 2)
    #
    from <- if (is_relative) log(xlim[1]) else xlim[1]
    to <- if (is_relative) log(xlim[2]) else xlim[2]
  }
  #
  by <- (to - from) / 1000
  #
  seq <- c(seq(from, to, by),
           TE.random + c(-1e-8, 0, 1e-8), 
           lower.predict + c(-1e-8, 0, 1e-8),
           upper.predict + c(-1e-8, 0, 1e-8))
  #
  if (avail.cid.below.null) {
    for (i in seq_along(cid.below.null))
      seq <- c(seq, cid.below.null.transf[i] + c(-1e-8, 0, 1e-8))
  }
  #
  if (avail.cid.above.null) {
    for (i in seq_along(cid.above.null))
      seq <- c(seq, cid.above.null.transf[i] + c(-1e-8, 0, 1e-8))
  }
  #
  seq <- unique(sort(seq))
  #
  # Calculate density
  #
  dens <- dt((seq - TE.random) / seTE.predict, df.predict)
  #
  dat <- tibble(xval = seq, yval = dens)
  #
  if (avail.cid.below.null & avail.cid.above.null) {
    dat$cid.category <-
      cut(dat$xval,
          breaks = unique(c(-Inf, cid.below.null.transf,
                            cid.above.null.transf, Inf)))
    #
    if (prop.within.cid == 0)
      fill.category <- c(fill.cid.below.null, fill.cid.above.null)
    else
      fill.category <- c(fill.cid.below.null, fill, fill.cid.above.null)
  }
  else if (avail.cid.below.null) {
    dat$cid.category <-
      cut(dat$xval, breaks = c(-Inf, cid.below.null.transf, Inf))
    #
    fill.category <- c(fill.cid.below.null, fill)
  }
  #
  else if (avail.cid.above.null) {
    dat$cid.category <-
      cut(dat$xval, breaks = c(-Inf, cid.above.null.transf, Inf))
    #
    fill.category <- c(fill, fill.cid.above.null)
  }
  #
  if (is_relative) {
    dat$xval <- exp(dat$xval)
    #
    seq <- exp(seq)
  }
  #
  # Only keep values within limits of x-axis
  #
  if (!missing.xlim)
    dat %<>% filter(xval >= xlim[1] & xval <= xlim[2])
  #
  # Data for diamond of random effects model
  #
  x.diamond <- c(TE.random, upper.random, TE.random, lower.random)
  dat.diamond <- data.frame(
    x = if (is_relative) exp(x.diamond) else x.diamond,
    y = -0.03 + c(0.005, 0, -0.005, 0)
  )
  #
  # Data for prediction interval
  #
  dat.predict <- data.frame(
    x1 = if (is_relative) exp(lower.predict) else lower.predict,
    x2 = if (is_relative) exp(upper.predict) else upper.predict,
    y1 = -0.045, y2 = -0.05)
  
  #
  #
  # (7) Create plot
  #
  #
  
  p <- ggplot(dat, aes(xval, yval)) +
    geom_line() +
    geom_area(fill = fill) +
    # Add diamond for random effects model
    geom_polygon(data = dat.diamond, aes(x, y), 
                 fill = col.diamond, color = col.diamond.lines) +
    # Add prediction interval
    annotate("rect",
             xmin = if (is_relative) exp(lower.predict) else lower.predict,
             xmax = if (is_relative) exp(upper.predict) else upper.predict,
             ymin = -0.05, ymax = -0.045,
             fill = col.predict, color = col.predict.lines) +
  xlab(xlab) + ylab(ylab)
  #
  if (avail.cid.below.null | avail.cid.above.null) {
    levs.category <- levels(dat$cid.category)
    #
    for (i in seq_along(levs.category)) {
      p <- p +
        geom_area(data = filter(dat, cid.category == levs.category[i]),
                  fill = fill.category[i])
    }
  }
  #
  # Add individual study results
  #
  if (studies) {
    dat.TE <-
      data.frame(TE = if (is_relative) exp(TE) else TE,
                 yval = 0, size = 10 * w.random / max(w.random))
    #
    if (!missing.xlim)
      dat.TE %<>% filter(TE >= xlim[1] & TE <= xlim[2])
    #
    p <- p +
      geom_point(data = dat.TE, aes(x = TE, y = yval),
                 shape = 1, size = dat.TE$size)
  }
  #
  if (!is.null(ylim)) {
    chknumeric(ylim, length = 2)
    #
    p <- p + ylim(ylim[1], ylim[2])
  }
  #
  # Add labels on x-axis
  #
  if (missing(labels.x)) {
    if (is_relative)
      p <- p +
        scale_x_continuous(
          trans = "log",
          labels = number_format(accuracy = 1 / 10^digits.xaxis))
    else
      p <- p +
        scale_x_continuous(
          labels = number_format(accuracy = 1 / 10^digits.xaxis))
  }
  else {
    if (is_relative)
      p <- p +
        scale_x_continuous(trans = "log", breaks = labels.x, labels = labels.x)
    else
      p <- p +
        scale_x_continuous(breaks = labels.x, labels = labels.x)
  }
  #
  # Add legend with expected proportions
  #
  if (legend) {
    dat.legend <- dat.cid %>% mutate(Colour = "")
    #
    gtab <- gridExtra::tableGrob(dat.legend)
    #
    for (i in seq_len(nrow(gtab$layout))) {
      cell <- gtab$layout[i, ]
      # First column in data set is equal to 'cell$l == 2'
      if (cell$l == 2 && cell$name %in% c("core-fg", "colhead-fg")) {
        #
        row_index <- cell$t
        # Coloured background for each cell in the target column
        if (row_index > 1)
          gtab$grobs[[i]] <-
            grobTree(
              rectGrob(gp =
                         gpar(fill = dat.cid$Colour[row_index - 1], col = NA)),
              gtab$grobs[[i]]
          )
      }
    }
    #
    p <- ggpubr::ggarrange(p, gtab, ncol = 1, nrow = 2, heights = c(3, 1))
  }
  #
  attr(p, "data.cid") <- dat.cid
  #
  class(p) <- c(class(p), "plot.cidprop")
  #
  p
}
