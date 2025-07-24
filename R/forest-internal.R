##
##
## Definition of auxiliary functions for forest plots
##
##


add.label <- function(x, column,
                      xpos, ypos, just, fs.lr, ff.lr, col,
                      fontfamily,
                      ...) {
  ##
  pushViewport(viewport(layout.pos.col = column, ...))
  ##
  grid.text(x, x = xpos, y = ypos, just = just,
            gp = gpar(fontsize = fs.lr, fontface = ff.lr, col = col,
                      fontfamily = fontfamily))
  ##
  popViewport()
  ##
  invisible(NULL)
}


add.text <- function(x, column, ...) {
  ##
  for (i in seq_len(length(x$rows))) {
    if (!is.na(x$rows[i])) {
      pushViewport(
        viewport(
          layout.pos.row = x$rows[i],
          layout.pos.col = column, ...))
      #
      grid.draw(x$labels[[i]])
      #
      popViewport()
    }
  }
  ##
  invisible(NULL)
}


add.xlab <- function(x, column, xlab, xlab.add, newline.xlab,
                     xpos, ypos, fs.xlab, ff.xlab,
                     fontfamily) {
  ##
  pushViewport(viewport(layout.pos.col = column, xscale = x$range))
  ##
  ## Label on x-axis:
  ##
  grid.text(xlab,
            x = unit(xpos, "native"),
            y = unit(ypos, "lines"),
            just = "center",
            gp = gpar(fontsize = fs.xlab, fontface = ff.xlab,
                      fontfamily = fontfamily))
  ##
  if (newline.xlab)
    grid.text(xlab.add,
              x = unit(xpos, "native"),
              y = unit(ypos - 1, "lines"),
              just = "center",
              gp = gpar(fontsize = fs.xlab, fontface = ff.xlab,
                        fontfamily = fontfamily))
  ##
  popViewport()
  ##
  invisible(NULL)
}


add.rob <- function(x, column, size, fs, ff, fontfamily,
                    rob, rob.levels, rob.symbols, rob.colour, ...) {
  ##
  if (is.null(rob.levels)) {
    if (is.factor(rob))
      rob.levels <- levels(rob)
    else
      rob.levels <- unique(rob)
  }
  ##
  n.levs <- length(rob.levels)
  ##
  if (is.null(rob.symbols)) {
    if (n.levs == 3)
      rob.symbols <- c("-", "?", "+")
    else
      rob.symbols <- rev(seq_len(n.levs))
  }
  else if (length(rob.symbols) == 1 && is.logical(rob.symbols)) {
    if (rob.symbols) {
      if (n.levs == 3)
        rob.symbols <- c("-", "?", "+")
      else
        rob.symbols <- rev(seq_len(n.levs))
    }
  }
  else if (length(rob.symbols) != n.levs)
    stop("Wrong number of RoB symbols (argument 'rob.symbols').",
         call. = FALSE)
  ##
  if (is.null(rob.colour)) {
    if (n.levs == 3)
      rob.colour <- c("red", "yellow", "green")
    else
      rob.colour <- rev(seq_len(n.levs))
  }
  ##
  if (!(length(rob.symbols) == 1 && is.logical(rob.symbols) && !rob.symbols))
    txt.rob <-
      as.character(
        factor(rob, levels = rob.levels, labels = rob.symbols))
  else
    txt.rob <- NULL
  ##
  col.rob <-
    as.character(
      factor(rob, levels = rob.levels,
             labels = rob.colour))
  ##
  j <- 0
  ##
  for (i in seq_len(length(x$rows))) {
    if (!is.na(x$rows[i])) {
      pushViewport(
        viewport(
          layout.pos.row = x$rows[i],
          layout.pos.col = column, ...))
      ##
      if (i == 1)
        grid.draw(x$labels[[1]])
      else {
        j <- j + 1
        ##
        grid.circle(x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                    r = unit(size / 2, "snpc"),
                    gp = gpar(fill = col.rob[j], col = col.rob[j]))
        ##
        if (!is.null(txt.rob) && !is.na(txt.rob[j]))
          grid.text(txt.rob[j],
                    x = unit(0.5, "npc"),
                    y = unit(if (txt.rob[j] %in% c("-", "+")) 0.58 else 0.5,
                             "npc"),
                    gp = gpar(fontsize = fs, fontface = ff,
                              fontfamily = fontfamily),
                    just = c("center", "center"))
      }
      ##
      popViewport()
    }
  }
  ##
  invisible(NULL)
}


draw.axis <- function(x, column, yS, log.xaxis, at, label,
                      fs.axis, ff.axis, fontfamily, lwd,
                      xlim, notmiss.xlim,
                      col.line, col.label) {
  ##
  ## Function to draw x-axis
  ##
  pushViewport(
    viewport(
      layout.pos.row = max(yS, na.rm = TRUE),
      layout.pos.col = column,
      xscale = x$range))
  ##
  ## x-axis:
  ##
  if (log.xaxis) {
    if (is.null(at)) {
      x1000 <- c(0.001, 0.1, 1,  10, 1000)
      x100  <- c(0.01 , 0.1, 1,  10, 100)
      x10   <- c(0.1  , 0.5, 1,   2, 10)
      x5    <- c(0.2  , 0.5, 1,   2, 5)
      x2    <- c(0.5  , 1, 2)
      x1.5  <- c(0.75 , 1, 1.5)
      x1.25 <- c(0.8  , 1, 1.25)
      x1    <- c(0.9  , 1, 1.1)
      ##
      min.x <- min(exp(x$range[1]), 1)
      max.x <- max(exp(x$range[2]), 1)
      ##
      if (all(x1000 >= min.x) &
          all(x1000 <= max.x))
        label <- x1000
      else if (all(x100 >= min.x) &
               all(x100 <= max.x))
        label <- x100
      else if (all(x10 >= min.x) &
               all(x10 <= max.x))
        label <- x10
      else if (all(x5 >= min.x) &
               all(x5 <= max.x))
        label <- x5
      else if (all(x2 >= min.x) &
               all(x2 <= max.x))
        label <- x2
      else if (all(x1.5 >= min.x) &
               all(x1.5 <= max.x))
        label <- x1.5
      else if (all(x1.25 >= min.x) &
               all(x1.25 <= max.x))
        label <- x1.25
      else if (all(x1 >= min.x) &
               all(x1 <= max.x))
        label <- x1
      else
        label <- 1
      ##
      if (notmiss.xlim && is.numeric(xlim[1])) {
        if (exp(min(xlim)) < min(label))
          label <- c(exp(min(xlim)), label)
        if (exp(max(xlim)) > max(label))
          label <- c(label, exp(max(xlim)))
      }
      at <- log(label)
    }
    else {
      if (length(label) == 1 && is.logical(label) && label)
        label <- at
      at <- log(at)
    }
    ## Print x-axis labels
    grid.xaxis(name = "xaxis1",
               at = at, label = label,
               gp = gpar(fontsize = fs.axis, fontface = ff.axis,
                         fontfamily = fontfamily, lwd = lwd,
                         col = col.label, tcl = -0.1))
    ## Print xaxis and tick marks (in different colour)
    grid.xaxis(name = "xaxis2",
               at = at, label = FALSE,
               gp = gpar(fontsize = fs.axis, fontface = ff.axis,
                         fontfamily = fontfamily, lwd = lwd,
                         col = col.line, tcl = -0.1))
  }
  else {
    if (is.null(at)) {
      ## Print x-axis labels
      grid.xaxis(name = "xaxis1",
        gp = gpar(fontsize = fs.axis, fontface = ff.axis,
                  fontfamily = fontfamily, lwd = lwd,
                  col = col.label, tcl = -0.1))
      ## Print xaxis and tick marks (in different colour)
      grid.xaxis(name = "xaxis2",
        label = FALSE,
        gp = gpar(fontsize = fs.axis, fontface = ff.axis,
                  fontfamily = fontfamily, lwd = lwd,
                  col = col.line, tcl = -0.1))
    }
    else if ((length(label) == 1 && is.logical(label) && label) |
          (length(label) >= 1 & !is.logical(label))) {
      ## Print x-axis labels
      grid.xaxis(name = "xaxis1",
        at = at, label = label,
        gp = gpar(fontsize = fs.axis, fontface = ff.axis,
                  fontfamily = fontfamily, lwd = lwd,
                  col = col.label, tcl = -0.1))
      ## Print xaxis and tick marks (in different colour)
      grid.xaxis(name = "xaxis2",
        at = at, label = FALSE,
        gp = gpar(fontsize = fs.axis, fontface = ff.axis,
                  fontfamily = fontfamily, lwd = lwd,
                  col = col.line, tcl = -0.1))
    }
    else {
      ## Print x-axis labels
      grid.xaxis(name = "xaxis1",
                 at = at,
                 gp = gpar(fontsize = fs.axis, fontface = ff.axis,
                           fontfamily = fontfamily, lwd = lwd,
                           col = col.label, tcl = -0.1))
      ## Print xaxis and tick marks (in different colour)
      grid.xaxis(name = "xaxis2",
        at = at, label = FALSE,
        gp = gpar(fontsize = fs.axis, fontface = ff.axis,
                  fontfamily = fontfamily, lwd = lwd,
                  col = col.line, tcl = -0.1))
    }
  }
  ##
  popViewport()
  ##
  invisible(NULL)
}


draw.ci.diamond <- function(TE, lower, upper,
                            size, min, max,
                            col.diamond, col.diamond.lines,
                            lwd) {
  ##
  if (min > max) {
    tmp <- min
    min <- max
    max <- tmp
  }
  ##
  if (!is.na(TE) &&
      ((min <= TE & TE <= max) |
       (min <= lower & lower <= max) |
       (min <= upper & upper <= max))
      ) {
    if (min <= lower & max >= upper) {
      grid.polygon(x = unit(c(lower, TE, upper, TE), "native"),
                   y = unit(0.5 + c(0, 0.4 * size, 0, -0.4 * size), "npc"),
                   gp = gpar(fill = col.diamond, col = col.diamond.lines,
                             lwd = lwd))
    }
    ##
    else {
      if (min > lower) {
        x.min <- min
        y.min1 <- 0.5 + -0.4 * size * (lower - min) / (lower - TE)
        y.min2 <- 0.5 +  0.4 * size * (lower - min) / (lower - TE)
      }
      else {
        x.min <- lower
        y.min1 <- y.min2 <- 0.5
      }
      ##
      if (max < upper) {
        x.max <- max
        y.max1 <- 0.5 +  0.4 * size * (upper - max) / (upper - TE)
        y.max2 <- 0.5 + -0.4 * size * (upper - max) / (upper - TE)
      }
      else {
        x.max <- upper
        y.max1 <- y.max2 <- 0.5
      }
      ##
      grid.polygon(x = unit(c(x.min, x.min, TE, x.max, x.max, TE, x.min),
                            "native"),
                   y = unit(c(y.min1, y.min2, 0.5 + 0.4 * size,
                              y.max1, y.max2, 0.5 - 0.4 * size,
                              y.min1), "npc"),
                   gp = gpar(fill = col.diamond, col = col.diamond.lines,
                             lwd = lwd))
    }
  }
  ##
  invisible(NULL)
}


draw.ci.predict <- function(lower.predict, upper.predict,
                            size, min, max,
                            col.predict, col.predict.lines) {
  ##
  if (min > max) {
    tmp <- min
    min <- max
    max <- tmp
  }
  ##
  if (!(is.na(lower.predict) | is.na(upper.predict))) {
    range <- max - min
    ## Plot prediction interval only within plotting range
    if (min > lower.predict)
      x.min <- min + range / 40
    else
      x.min <- lower.predict
    ##
    if (max < upper.predict)
      x.max <- max - range / 40
    else
      x.max <- upper.predict
    ##
    grid.polygon(x = unit(c(x.min, x.min, x.max, x.max), "native"),
                 y = unit(0.5 + size * c(-1, 1, 1, -1) / 10, "npc"),
                 gp = gpar(fill = col.predict, col = col.predict.lines))
    ##
    if (min > lower.predict)
      grid.lines(x = unit(c(min, min + 0.00001), "native"),
                 y = 0.5,
                 gp = gpar(col = col.predict.lines,
                           fill = col.predict),
                 arrow = arrow(ends = "first",
                               length = unit(0.5, "npc"),
                               type = "closed"))
    ##
    if (max < upper.predict)
      grid.lines(x = unit(c(max - 0.00001, max), "native"),
                 y = 0.5,
                 gp = gpar(col = col.predict.lines,
                           fill = col.predict),
                 arrow = arrow(ends = "last",
                               length = unit(0.5, "npc"),
                               type = "closed"))
  }
  ##
  invisible(NULL)
}


draw.ci <- function(TE, lower, upper,
                    size, min, max,
                    lwd,
                    col,
                    col.square, col.square.lines,
                    col.circle, col.circle.lines,
                    col.inside,
                    type,
                    lwd.square,
                    arrow.type, arrow.length) {
  ##
  if (min > max) {
    tmp <- min
    min <- max
    max <- tmp
  }
  ##
  if (!is.na(TE)) {
    ##
    ## Draw lines in colour "col.inside" if totally inside rect
    ##
    TElineCol <- col
    ##
    if ((!is.na(size) & !is.na(lower) & !is.na(upper)))
      if (size > 0 &&
          (convertX(unit(TE, "native") + unit(0.5 * size, "lines"),
                    "native", valueOnly = TRUE) > upper) &&
          (convertX(unit(TE, "native") - unit(0.5 * size, "lines"),
                    "native", valueOnly = TRUE) < lower))
        TElineCol <- col.inside
  }
  ##
  if (!is.na(TE) && (TE >= min & TE <= max)) {
    if (!is.na(size) && size > 0 && !is.na(lower) && !is.na(upper)) {
      if (type == "square") {
        grid.rect(x = unit(TE, "native"),
                  width = unit(size, "snpc"),
                  height = unit(size, "snpc"),
                  gp = gpar(fill = col.square, col = col.square.lines,
                            lwd = lwd.square))
        ##
        grid.lines(x = unit(c(TE, TE), "native"),
                   y = unit(c(0.4, 0.6), "npc"),
                   gp = gpar(col = TElineCol, lwd = lwd))
      }
    }
    else
      grid.lines(x = unit(c(TE, TE), "native"),
                 y = unit(c(0.4, 0.6), "npc"),
                 gp = gpar(col = TElineCol, lwd = lwd))
  }
  ##
  if (!is.na(TE)) {
    ##
    ## Draw lines in colour "col.inside" if totally inside rect
    ##
    if (!is.na(size)) {
      lineCol <- col
      ##
      if (!is.na(lower) & !is.na(upper))
        if (size > 0 &&
            (convertX(unit(TE, "native") + unit(0.5 * size, "lines"),
                      "native", valueOnly = TRUE) > upper) &&
            (convertX(unit(TE, "native") - unit(0.5 * size, "lines"),
                      "native", valueOnly = TRUE) < lower))
          lineCol <- col.inside
      ##
      ## Draw arrow if exceed col range
      ## convertX() used to convert between coordinate systems
      ##
      if (!is.na(lower) && !is.na(upper) &&
          (lower >= min & upper <= max))
        grid.lines(x = unit(c(lower, upper), "native"), y = 0.5,
                   gp = gpar(col = lineCol, lwd = lwd))
      ##
      if (!is.na(lower) && !is.na(upper) &&
          (lower < min & upper > max))
        grid.lines(x = unit(c(min, max), "native"), y = 0.5,
                   gp = gpar(col = lineCol, lwd = lwd))
      ##
      if (!is.na(lower) && !is.na(upper) &&
          (lower < min & (upper <= max & upper > min)))
        grid.lines(x = unit(c(min, upper), "native"), y = 0.5,
                   gp = gpar(col = lineCol, lwd = lwd))
      ##
      if (!is.na(lower) && !is.na(upper) &&
          ((lower >= min & lower < max) & upper > max))
        grid.lines(x = unit(c(lower, max), "native"), y = 0.5,
                   gp = gpar(col = lineCol, lwd = lwd))
      ##
      if (!is.na(lower) && lower < min)
        grid.lines(x = unit(c(min - 0.00001, min), "native"), y = 0.5,
                   gp = gpar(col = lineCol, lwd = lwd, fill = lineCol),
                   arrow = arrow(ends = "first",
                                 length = unit(arrow.length, "inches"),
                                 type = arrow.type))
      if (!is.na(upper) && upper > max)
        grid.lines(x = unit(c(max, max + 0.00001), "native"), y = 0.5,
                   gp = gpar(col = lineCol, lwd = lwd, fill = lineCol),
                   arrow = arrow(ends = "last",
                                 length = unit(arrow.length, "inches"),
                                 type = arrow.type))
    }
  }
  ##
  if (!is.na(TE) && (TE >= min & TE <= max)) {
    if (!is.na(size) && size > 0 && !is.na(lower) && !is.na(upper)) {
      if (type == "circle")
        grid.circle(x = unit(TE, "native"), y = unit(0.5, "npc"),
                    r = unit(size / 2, "snpc"),
                    gp = gpar(fill = col.circle, col = col.circle.lines))
      ##
      else if (type == "squarediamond") {
        xmin <- convertX(unit(TE, "native") - unit(0.5 * size, "lines"),
                         "native", valueOnly = TRUE)
        xmax <- convertX(unit(TE, "native") + unit(0.5 * size, "lines"),
                         "native", valueOnly = TRUE)
        ##
        grid.polygon(x = unit(c(xmin, TE, xmax, TE), "native"),
                     y = unit(0.5 + c(0, 0.5 * size, 0, -0.5 * size), "npc"),
                     gp = gpar(col = col.square.lines, fill = col.square,
                               lwd = lwd.square))
      }
    }
    else
      grid.lines(x = unit(c(TE, TE), "native"),
                 y = unit(c(0.4, 0.6), "npc"),
                 gp = gpar(col = TElineCol, lwd = lwd))
  }
  ##
  invisible(NULL)
}


draw.forest <- function(x, column) {
  ##
  ## Function to plot results for individual studies and summaries
  ##
  for (i in seq_len(length(x$rows))) {
    if (!is.na(x$rows[i])) {
      pushViewport(
        viewport(
          layout.pos.row = x$rows[i],
          layout.pos.col = column,
          xscale = x$range))
      ##
      if (x$type[i] %in% c("square", "circle", "squarediamond"))
        draw.ci(x$eff[i], x$low[i], x$upp[i],
                x$sizes[i], x$range[1], x$range[2],
                x$lwd,
                x$col[i],
                x$col.square[i], x$col.square.lines[i],
                x$col.circle[i], x$col.circle.lines[i],
                x$col.inside[i],
                type = x$type[i],
                x$lwd.square,
                x$arrow.type, x$arrow.length)
      ##
      else if (x$type[i] == "diamond")
        draw.ci.diamond(x$eff[i], x$low[i], x$upp[i],
                        x$sizes[i], x$range[1], x$range[2],
                        x$col.diamond[i], x$col.diamond.lines[i],
                        x$lwd.diamond)
      ##
      else if (x$type[i] == "predict")
        draw.ci.predict(x$low[i], x$upp[i],
                        x$sizes[i], x$range[1], x$range[2],
                        x$col.diamond[i], x$col.diamond.lines[i])
      ##
      popViewport()
    }
  }
  ##
  invisible(NULL)
}


draw.lines <- function(x, column,
                       ref, TE.common, TE.random,
                       overall, common, random, prediction,
                       ymin.common, ymin.random, ymin.ref, ymax, ymax.ref,
                       lwd, lty.common, lty.random, col.common, col.random,
                       xmin, xmax,
                       cid.below.null, cid.above.null,
                       lty.cid, col.cid,
                       fill.cid.below.null, fill.cid.above.null,
                       fill,
                       col.line) {
  ##
  if (xmin > xmax) {
    xmin <- x$range[2]
    xmax <- x$range[1]
  }
  else {
    xmin <- x$range[1]
    xmax <- x$range[2]
  }
  ##
  pushViewport(viewport(layout.pos.col = column, xscale = x$range))
  ##
  ## Add background colour for confidence interval plot
  ##
  if (!is.null(fill))
    grid.polygon(x = unit(c(xmin, xmax, xmax, xmin), "native"),
                 y = unit(c(ymin.ref, ymin.ref, ymax, ymax),
                          "lines"),
                 gp = gpar(col = fill, fill = fill))
  ##
  ## Add equivalence region
  ##
  if (is.na(ref) & any(!is.na(cid.below.null)) & any(!is.na(cid.above.null)))
    ref.equi <- min(cid.below.null, na.rm = TRUE) +
      0.5 * (max(cid.above.null, na.rm = TRUE) - min(cid.below.null, na.rm = TRUE))
  else
    ref.equi <- ref
  ##
  if (!is.na(ref.equi) && (xmin <= ref.equi & ref.equi <= xmax) &&
      any(!is.na(cid.below.null))) {
    ##
    n.cid.below.null <- sum(!is.na(cid.below.null))
    n.fill.cid.below.null <- length(fill.cid.below.null)
    ##
    if (n.cid.below.null < n.fill.cid.below.null) {
      firstline <- FALSE
      cid.below.null <- c(xmin, cid.below.null, ref.equi)
    }
    else {
      firstline <- TRUE
      cid.below.null <- c(cid.below.null, ref.equi)
    }
    ##
    n.lo <- length(cid.below.null[-1])
    ##
    for (i in seq_len(n.lo)) {
      if (!is.na(cid.below.null[i]) && !is.na(cid.below.null[i + 1]) &&
          ((xmin <= cid.below.null[i] & cid.below.null[i] <= cid.below.null[i + 1]) &
           (cid.below.null[i + 1] <= xmax))) {
        ##
        grid.polygon(x = unit(c(cid.below.null[i], cid.below.null[i + 1],
                                cid.below.null[i + 1], cid.below.null[i]), "native"),
                     y = unit(c(ymin.ref, ymin.ref, ymax.ref, ymax.ref),
                              "lines"),
                     gp = gpar(lwd = lwd, col = "transparent",
                               fill = fill.cid.below.null[i]))
      }
    }
    ##
    for (i in seq_len(n.lo)) {
      if (!(i == 1 & !firstline) &
          !is.na(cid.below.null[i]) && !is.na(cid.below.null[i + 1]) &&
          ((xmin <= cid.below.null[i] & cid.below.null[i] <= cid.below.null[i + 1]) &
           (cid.below.null[i + 1] <= xmax)))
        grid.lines(x = unit(cid.below.null[i], "native"),
                   y = unit(c(ymin.ref, ymax.ref), "lines"),
                   gp = gpar(lwd = lwd, col = col.cid, lty = lty.cid))
    }
  }
  ##
  if (!is.na(ref.equi) && (xmin <= ref.equi & ref.equi <= xmax) &&
      any(!is.na(cid.above.null))) {
    ##
    n.cid.above.null <- sum(!is.na(cid.above.null))
    n.fill.cid.above.null <- length(fill.cid.above.null)
    ##
    if (n.cid.above.null < n.fill.cid.above.null) {
      lastline <- FALSE
      cid.above.null <- c(ref.equi, cid.above.null, xmax)
    }
    else {
      lastline <- TRUE
      cid.above.null <- c(ref.equi, cid.above.null)
    }
    ##
    n.up <- length(cid.above.null[-1])
    ##
    for (i in seq_len(n.up)) {
      if (!is.na(cid.above.null[i]) && !is.na(cid.above.null[i + 1]) &&
          ((xmin <= cid.above.null[i] & cid.above.null[i] <= cid.above.null[i + 1]) &
           (cid.above.null[i + 1] <= xmax))) {
        ##
        grid.polygon(x = unit(c(cid.above.null[i], cid.above.null[i + 1],
                                cid.above.null[i + 1], cid.above.null[i]), "native"),
                     y = unit(c(ymin.ref, ymin.ref, ymax.ref, ymax.ref),
                              "lines"),
                     gp = gpar(lwd = lwd, col = "transparent",
                               fill = fill.cid.above.null[i]))
      }
    }
    ##
    for (i in seq_len(n.up)) {
      if ((i != n.up | (i == n.up & lastline)) &
          !is.na(cid.above.null[i + 1]) &&
          (xmin <= cid.above.null[i + 1] & cid.above.null[i + 1] <= xmax))
          grid.lines(x = unit(cid.above.null[i + 1], "native"),
                     y = unit(c(ymin.ref, ymax.ref), "lines"),
                     gp = gpar(lwd = lwd, col = col.cid, lty = lty.cid))
    }
  }
  ##
  ## Reference line:
  ##
  if (!is.na(ref) && (xmin <= ref & ref <= xmax))
    grid.lines(x = unit(ref, "native"),
               y = unit(c(ymin.ref, ymax), "lines"),
               gp = gpar(lwd = lwd, col = col.line))
  ##
  ## Line for common effect estimate(s):
  ##
  if (common & overall)
    for (i in seq_along(TE.common))
      if (!is.na(TE.common[i]))
        if (xmin <= TE.common[i] & TE.common[i] <= xmax)
          if (!is.null(lty.common))
            grid.lines(x = unit(TE.common[i], "native"),
                       y = unit(c(ymin.common + length(TE.common) - i, ymax),
                                "lines"),
                       gp = gpar(lty = lty.common, lwd = lwd, col = col.common))
  ##
  ## Line for random effects estimate(s):
  ##
  if (random & overall)
    for (i in seq_along(TE.random))
      if (!is.na(TE.random[i]))
        if (xmin <= TE.random[i] & TE.random[i] <= xmax)
          if (!is.null(lty.random) & !is.na(TE.random[i]))
            grid.lines(x = unit(TE.random[i], "native"),
                       y = unit(c(ymin.random + length(TE.random) - i, ymax),
                                "lines"),
                       gp = gpar(lty = lty.random,
                                 lwd = lwd, col = col.random))
  ##
  popViewport()
  ##
  invisible(NULL)
}


formatcol <- function(x, y, rows, just = "right", settings,
                      fontfamily,
                      n.com, n.ran, n.prd,
                      rob = FALSE) {
  ##
  if (just == "left")
    xpos <- 0
  if (just == "center")
    xpos <- 0.5
  if (just == "right")
    xpos <- 1
  ##
  res <- list(labels = 
                lapply(c(x, as.list(y)),
                       textGrob, x = xpos, just = just,
                       gp = gpar(
                         fontsize = settings$fs.study,
                         fontface = settings$ff.study,
                         fontfamily = fontfamily)
                       ),
              rows = rows,
              rob = rob)
  ##
  ## Study label:
  ##
  res$labels[[1]] <- textGrob(x,
                              x = xpos, just = just,
                              gp = gpar(
                                fontsize = settings$fs.heading,
                                fontface = settings$ff.heading,
                                fontfamily = fontfamily)
                              )
  ##
  ## Common effect estimate:
  ##
  strt <- j <- 1
  for (i in seq_len(n.com)) {
    res$labels[[strt + i]] <- textGrob(y[strt - 1 + i],
                                       x = xpos, just = just,
                                       gp = gpar(
                                         fontsize = settings$fs.common,
                                         fontface = settings$ff.common,
                                         fontfamily = fontfamily)
                                       ) 
    j <- j + 1
  }
  ##
  ## Random effects estimate:
  ##
  strt <- j
  for (i in seq_len(n.ran)) {
    res$labels[[strt + i]] <- textGrob(y[strt - 1 + i],
                                       x = xpos, just = just,
                                       gp = gpar(
                                         fontsize = settings$fs.random,
                                         fontface = settings$ff.random,
                                         fontfamily = fontfamily)
                                       )
    j <- j + 1
  }
  ##
  ## Prediction interval:
  ##
  strt <- j
  for (i in seq_len(n.prd)) {
    res$labels[[strt + i]] <- textGrob(y[strt - 1 + i],
                                       x = xpos, just = just,
                                       gp = gpar(
                                         fontsize = settings$fs.predict,
                                         fontface = settings$ff.predict,
                                         fontfamily = fontfamily)
                                       )
    j <- j + 1
  }
  ##
  if (settings$by) {
    n.by <- settings$n.by
    strt <- j
    ##
    ## Common effect estimates:
    ##
    for (i in seq_len(n.by * n.com)) {
      res$labels[[strt + i]] <-
        textGrob(y[strt - 1 + i],
                 x = xpos, just = just,
                 gp = 
                   gpar(
                     fontsize = settings$fs.common,
                     fontface = settings$ff.common,
                     fontfamily = fontfamily,
                     col = settings$col.subgroup)
                 )
      j <- j + 1
    }
    ##
    ## Random effects estimates:
    ##
    strt <- j
    for (i in seq_len(n.by * n.ran)) {
      res$labels[[strt + i]] <-
        textGrob(y[strt - 1 + i],
                 x = xpos, just = just,
                 gp = 
                   gpar(
                     fontsize = settings$fs.random,
                     fontface = settings$ff.random,
                     fontfamily = fontfamily,
                     col = settings$col.subgroup)
                 )
      j <- j + 1
    }
    ##
    ## Prediction interval:
    ##
    strt <- j
    for (i in seq_len(n.by * n.prd)) {
      res$labels[[strt + i]] <-
        textGrob(y[strt - 1 + i],
                 x = xpos, just = just,
                 gp = 
                   gpar(
                     fontsize = settings$fs.predict,
                     fontface = settings$ff.predict,
                     fontfamily = fontfamily,
                     col = settings$col.subgroup)
                 )
      j <- j + 1
    }
  }
  ##
  res
}


removeNULL <- function(x, names, varname) {
  if (is.null(x[[varname]]))
    res <- names[names != varname]
  else
    res <- names
  ##
  res
}


tg <- function(x, xpos, just, fs, ff, fontfamily, col) {
  if (missing(col))
    res <- textGrob(x,
                    x = xpos, just = just,
                    gp = gpar(fontsize = fs, fontface = ff,
                              fontfamily = fontfamily))
  else
    res <- textGrob(x,
                    x = xpos, just = just,
                    gp = gpar(fontsize = fs, fontface = ff,
                              fontfamily = fontfamily,
                              col = col))
  ##
  res
}


tgl <- function(x, xpos, just, fs, ff, fontfamily, rows = 1, col) {
  ##
  if (missing(col))
    res <- list(labels = list(tg(x, xpos, just, fs, ff, fontfamily)),
                rows = rows)
  else
    res <- list(labels = list(tg(x, xpos, just, fs, ff, fontfamily, col)),
                rows = rows)
  ##
  res
}


twolines <- function(x, xname = deparse(substitute(x)), arg = FALSE) {
  newline <- FALSE
  bottom <- longer <- x
  top <- NULL
  ##
  if (!is.null(bottom)) {
    if (grepl("\n", bottom)) {
      wsplit <- unlist(strsplit(bottom, "\n"))
      if (length(wsplit) == 1) {
        if (substring(bottom, 1, 1) == "\n") {
          top <- ""
          bottom <- wsplit[1]
          longer <- bottom
          newline <- TRUE
        }
        else if (substring(bottom, nchar(bottom), nchar(bottom)) == "\n") {
          top <- wsplit[1]
          bottom <- ""
          longer <- top
          newline <- TRUE
        }
      }
      else if (length(wsplit) != 2) {
        if (arg)
          stop("Maximum of two lines for argument '", xname, "'.",
               call. = FALSE)
        else
          stop("Maximum of two lines for label of column '", xname, "'.",
               call. = FALSE)
      }
      else {
        top <- wsplit[1]
        bottom <- wsplit[2]
        longer <- ifelse(nchar(top) > nchar(bottom), top, bottom)
        newline <- TRUE
      }
    }
  }
  ##
  list(newline = newline, bottom = bottom, top = top, longer = longer)
}


wcalc <- function(x)
  max(unit(rep(1, length(x)), "grobwidth", x))


collapsemat <- function(x) {
  if (is.list(x)) {
    for (i in rev(seq_len(length(x)))) {
      if (is.null(x[[i]]))
        x[[i]] <- NULL
    }
  }
  ##
  if (is.list(x) & length(x) == 1)
    x <- x[[1]]
  ##
  if (is.matrix(x)) {
    res <- as.vector(t(x))
    names(res) <- rep(rownames(x), rep(ncol(x), nrow(x)))
  }
  else
    res <- x
  ##
  res
}


ordermat <- function(x, levs) {
  o <- order(factor(names(x), levels = levs))
  x[o]
}


selmat <- function(x, levs) {
  o <- order(factor(names(x), levels = levs))
  x[o]
}


notallNA <- function(x)
  any(!is.na(x))


repl <- function(x, n1, n2)
  rep(x, rep(n1, n2))


gh <- function(type.gr, rows.gr,
               ##
               n.stud,
               lower.common, lower.random, lower.predict,
               subgroup, subgroup.levels,
               lower.common.w, lower.random.w, lower.predict.w,
               ##
               common, random, overall,
               prediction, overall.hetstat,
               study.results,
               ##
               spacing,
               ##
               xlab, xlab.add, label.right, label.left, bottom.lr,
               ##
               prediction.subgroup, subgroup.hetstat,
               test.overall.common, test.overall.random,
               test.subgroup.common, test.subgroup.random,
               ##
               text.addline1, text.addline2,
               text.details, text.rob,
               ##
               addrow, addrow.overall,
               addrow.subgroups, addrows.below.overall,
               ##
               cols, labs,
               text.w.common, text.w.random) {
  
  
  ##
  ## (1) Determine height per row
  ##
  if (grepl("bmp$", tolower(type.gr)) ||
      grepl("jpg$", tolower(type.gr)) ||
      grepl("jpeg$", tolower(type.gr)) ||
      grepl("png$", tolower(type.gr)) ||
      grepl("tif$", tolower(type.gr)) ||
      grepl("tiff$", tolower(type.gr)))
    height_per_row <- 480 / 33
  else
    height_per_row <-  7 / 35
  ##
  ## (2) Column labels
  ##
  labs <- unlist(labs)
  if (any(grepl("col.w.common", cols)))
    labs <- c(labs, text.w.common)
  if (any(grepl("col.w.random", cols)))
    labs <- c(labs, text.w.random)
  ##
  rows_column_labels <- 1 + 1L * any(grepl("\n", unlist(labs))) + 1L * addrow
  ##
  ## (2) Study results
  ##
  rows_studies <-
    if (!study.results)
      0
    else
      n.stud
  ##
  ## (3) Meta-analysis results
  ##
  if (overall)
    rows_overall <-
      1L * addrow.overall +
      1L * common * length(lower.common) +
      1L * random * length(lower.random) +
      1L * prediction * length(lower.predict) +
      1L * ((any(common) + any(random) + any(prediction)) == 1)
  else
    rows_overall <- 1L * addrow.overall
  ##
  ## (4) Text below meta-analysis results
  ##
  n.details <- sum(text.details != "")
  n.rob <- sum(text.rob != "")
  rows_below_overall_labels <- addrows.below.overall + overall.hetstat +
    test.overall.common + test.overall.random +
    test.subgroup.common + test.subgroup.random +
    1L * (text.addline1 != "") +
    1L * (text.addline2 != "") +
    n.details + 1L * (n.details > 0) +
    n.rob + 1L * (n.rob > 0) +
    1L * (n.details > 0 | n.rob > 0) +
    2L * (n.details == 0 & n.rob == 0)
  ## + 0.25
  ## 
  ## (5) Information below confidence interval plot
  ##
  rows_xlab <-
    if (xlab == "")
      0
    else
      2
  ##
  rows_xlab.add <-
    if (xlab.add == "")
      0
    else
      2
  ##
  rows_xlab <- rows_xlab + rows_xlab.add
  ##
  if (!bottom.lr) {
    rows_label.left <- 0
    rows_label.right <- 0
  }
  else {
    rows_label.left <-
      if (is.null(label.left) || label.left == "")
        0
      else if (!grepl("\n", label.left))
        2
      else
        4
    ##
    rows_label.right <-
      if (is.null(label.right) || label.right == "")
        0
      else if (!grepl("\n", label.right))
        2
      else
        4
  }
  ##
  rows_label <- max(c(rows_label.left, rows_label.right))
  ##
  rows_below_forest <- 2 + rows_xlab + rows_label
  ##
  ## (6) Subgroup results
  ##
  if (is.null(subgroup))
    rows_subgroups <- 0
  else {
    n.subgr <- length(subgroup.levels)
    ##
    rows_subgroups_common <- common * n.subgr
    rows_subgroups_random <- random * n.subgr
    ##
    rows_subgroups_predict <- 
      if (is.vector(prediction.subgroup))
        sum(prediction.subgroup)
      else
        prediction.subgroup * n.subgr
    ##
    rows_subgroups_hetstat <- 
      if (length(subgroup.hetstat) > 1)
        sum(subgroup.hetstat)
      else
        subgroup.hetstat * n.subgr
    ##
    if (is.matrix(lower.common.w))
      rows_subgroups_common <- rows_subgroups_common * nrow(lower.common.w)
    ##
    if (is.matrix(lower.random.w))
      rows_subgroups_random <- rows_subgroups_random * nrow(lower.random.w)
    ##
    if (is.matrix(lower.predict.w))
      rows_subgroups_predict <- rows_subgroups_predict * nrow(lower.predict.w)
    ##
    rows_subgroups <-
      n.subgr + # Labels
      addrow.subgroups * (n.subgr - 1) +
      rows_subgroups_common + rows_subgroups_random + 
      rows_subgroups_predict + rows_subgroups_hetstat
  }
  ##
  ## (7) Determine total height of graphics device
  ##
  total_rows <-
    rows_column_labels +
    rows_studies +
    rows_overall +
    max(rows_below_overall_labels, rows_below_forest) +
    rows_subgroups +
    rows.gr +
    1
  ##
  total_height <- height_per_row * spacing * total_rows
  
  res <- data.frame(total_height, total_rows, height_per_row, spacing)
  ##
  res
}


show_subgroup_results <- function(x, n, lower, upper) {
  if (length(x) == 1) {
    if (is.matrix(lower))
      return(x &
             apply(lower, 1, notallNA) &
             apply(upper, 1, notallNA))
    else
      return(rep(x & notallNA(lower) & notallNA(upper), n))
  }
  else {
    chklength(x, n,
              text = paste0("Length of argument '",
                            deparse(substitute(x)),
                            "' must be equal to 1 or number of subgroups."))
    return(x)
  }
}

newCol <- function(varname, label,
                   rob, data1, data2, datap,
                   n.com, n.ran, n.prd,
                   notavail, lab.NA, big.mark,
                   zero.pval, JAMA.pval, scientific.pval,
                   digits, digits.pval, digits.tau2, digits.tau, digits.I2,
                   sel.prd,
                   n.com.w = 0, n.ran.w = 0, n.prd.w = 0, n.stat.w = 0) {
  
  by <- sum(c(n.com.w, n.ran.w, n.prd.w, n.stat.w) > 0)
  #
  n.subgroup <- n.com.w + n.ran.w + n.prd.w + n.stat.w
  #
  if (length(rob[[varname]]) != 0)
    fvar <- rob[[varname]]
  else if (length(data1[[varname]]) != 0)
    fvar <- data1[[varname]]
  else if (length(data2[[varname]]) != 0)
    fvar <- data2[[varname]]
  else
    stop("Variable '", varname,
         "' not available in meta-analysis object.",
         call. = FALSE)
  #
  if (isCol(datap, varname)) {
    if (nrow(datap) == 1) {
      fvar <- c(rep(datap[[varname]], n.com),
                rep(datap[[varname]], n.ran),
                rep(NA, n.prd),
                if (by) rep(NA, n.subgroup) else NULL,
                fvar)
    }
    else if (nrow(datap) == n.com + n.ran)
      fvar <- c(datap[[varname]],
                rep(NA, n.prd),
                if (by) rep(NA, n.subgroup) else NULL,
                fvar)
    else if (nrow(datap) == n.com + n.ran + n.com.w + n.ran.w)
      fvar <- c(datap[[varname]][seq_len(n.com + n.ran)],
                rep(NA, n.prd),
                datap[[varname]][n.com + n.ran + seq_len(n.com.w + n.ran.w)],
                rep(NA, n.prd.w + n.stat.w),
                fvar)
    else
      stop("Wrong number of row in data set 'data.pooled' ",
           "(must be 1",
           if (by)", " else " or ", n.com + n.ran,,
           if (by) " or ", n.com.w + n.ran.w, ").",
           call. = FALSE)
  }
  #
  if (!is.character(fvar)) {
    if (is.factor(fvar))
      fvar <- as.character(fvar)
    else if (notavail & all(is_wholenumber(fvar), na.rm = TRUE))
      fvar <-
        formatN(fvar, digits = 0, text.NA = lab.NA, big.mark = big.mark)
    else if (is.numeric(fvar)) {
      if (varname == "pval")
        fvar <- formatPT(fvar, digits = digits.pval,
                          big.mark = big.mark,
                          lab = FALSE, labval = "",
                          zero = zero.pval, JAMA = JAMA.pval,
                          scientific = scientific.pval,
                          lab.NA = lab.NA)
      else if (varname == "tau2")
        fvar <- formatPT(fvar, digits = digits.tau2,
                          big.mark = big.mark,
                          lab = FALSE, labval = "",
                          lab.NA = lab.NA)
      else if (varname == "tau")
        fvar <- formatPT(fvar, digits = digits.tau,
                          big.mark = big.mark,
                          lab = FALSE, labval = "",
                          lab.NA = lab.NA)
      else if (varname == "I2") {
        sel.r <- !is.na(fvar)
        fvar[sel.r] <-
          paste0(formatN(100 * fvar[sel.r], digits.I2), "%")
        fvar[!sel.r] <- lab.NA
      }
      else
        fvar <-
          formatN(fvar, digits = digits, text.NA = lab.NA, big.mark = big.mark)
    }
  }
  #
  fvar <- ifelse(is.na(fvar), "", fvar)
  #
  # Print nothing for lines with prediction interval
  #
  if (isCol(datap, varname))
    fvar[sel.prd] <- ""
  else
    fvar <- c(rep("", n.com), rep("", n.ran), rep("", n.prd),
              if (by) rep("", n.subgroup),
              fvar)
  #
  # Print nothing for lines for subgroups
  #
  if (by) {
    is_NA <- fvar == lab.NA
    is_w <- rep(TRUE, length(fvar))
    is_w[seq_len(n.com + n.ran + n.prd)] <- FALSE
    #
    fvar[is_NA & is_w] <- ""
  }
  #
  # Check for "\n" in label of new column
  #
  clines <- twolines(label, varname)
  #
  if (clines$newline) {
    label <- clines$bottom
    longer <- clines$longer
  }
  else
    label <- longer <- label
  #
  res <- list(format_var = fvar, colname = paste0("col.", varname),
              label = label, longer = longer,
              rob = length(rob[[varname]]) != 0)
  #
  res
}
