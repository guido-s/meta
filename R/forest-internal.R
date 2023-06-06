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
  for (i in 1:length(x$rows)) {
    if (!is.na(x$rows[i])) {
      pushViewport(viewport(layout.pos.row = x$rows[i],
                            layout.pos.col = column, ...))
      grid.draw(x$labels[[i]])
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


draw.axis <- function(x, column, yS, log.xaxis, at, label,
                      fs.axis, ff.axis, fontfamily, lwd,
                      xlim, notmiss.xlim) {
  ##
  ## Function to draw x-axis
  ##
  pushViewport(viewport(layout.pos.row = max(yS, na.rm = TRUE),
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
    grid.xaxis(at = at, label = label,
               gp = gpar(fontsize = fs.axis, fontface = ff.axis,
                         fontfamily = fontfamily, lwd = lwd))
  }
  else {
    if (is.null(at))
      grid.xaxis(gp = gpar(fontsize = fs.axis, fontface = ff.axis,
                           fontfamily = fontfamily, lwd = lwd))
    else
      if ((length(label) == 1 && is.logical(label) && label) |
          (length(label) >= 1 & !is.logical(label)))
        grid.xaxis(at = at, label = label,
                   gp = gpar(fontsize = fs.axis, fontface = ff.axis,
                             fontfamily = fontfamily, lwd = lwd))
    else
      grid.xaxis(at = at,
                 gp = gpar(fontsize = fs.axis, fontface = ff.axis,
                           fontfamily = fontfamily, lwd = lwd))
  }
  ##
  popViewport()
  ##
  invisible(NULL)
}


draw.ci.diamond <- function(TE, lower, upper,
                            size, min, max,
                            col.diamond, col.diamond.lines) {
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
      )
    grid.polygon(x = unit(c(lower, TE, upper, TE), "native"),
                 y = unit(0.5 + c(0, 0.3 * size, 0, -0.3 * size), "npc"),
                 gp = gpar(fill = col.diamond, col = col.diamond.lines))
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
    ## Plot prediction interval only within plotting range
    if ((min <= lower.predict & lower.predict <= max) |
        (min <= upper.predict & upper.predict <= max))
      grid.polygon(x = unit(c(lower.predict, lower.predict,
                              upper.predict, upper.predict), "native"),
                   y = unit(0.5 + size * c(-1, 1, 1, -1) / 10, "npc"),
                   gp = gpar(fill = col.predict, col = col.predict.lines))
  }
  ##
  invisible(NULL)
}


draw.ci.square <- function(TE, lower, upper,
                           size, min, max,
                           lwd,
                           col, col.square,
                           col.square.lines, col.inside) {
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
      grid.rect(x = unit(TE, "native"),
                width = unit(size, "snpc"),
                height = unit(size, "snpc"),
                gp = gpar(fill = col.square, col = col.square.lines))
      grid.lines(x = unit(c(TE, TE), "native"),
                 y = unit(c(0.4, 0.6), "npc"),
                 gp = gpar(col = TElineCol, lwd = lwd))
    }
    else
      grid.lines(x = unit(c(TE, TE), "native"),
                 y = unit(c(0.4, 0.6), "npc"),
                 gp = gpar(col = TElineCol, lwd = lwd))
  }
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
                   gp = gpar(col = lineCol, lwd = lwd),
                   arrow = arrow(ends = "first", length = unit(0.05, "inches")))
      if (!is.na(upper) && upper > max)
        grid.lines(x = unit(c(max, max + 0.00001), "native"), y = 0.5,
                   gp = gpar(col = lineCol, lwd = lwd),
                   arrow = arrow(ends = "last", length = unit(0.05, "inches")))
    }
  }
  ##
  invisible(NULL)
}


draw.forest <- function(x, column) {
  ##
  ## Function to plot results for individual studies and summaries
  ##
  for (i in 1:length(x$rows)) {
    if (!is.na(x$rows[i])) {
      pushViewport(viewport(layout.pos.row = x$rows[i],
                            layout.pos.col = column,
                            xscale = x$range))
      ##
      if (x$type[i] == "square")
        draw.ci.square(x$eff[i], x$low[i], x$upp[i],
                       x$sizes[i], x$range[1], x$range[2],
                       x$lwd,
                       x$col[i], x$col.square[i],
                       x$col.square.lines[i], x$col.inside[i])
      ##
      else if (x$type[i] == "diamond")
        draw.ci.diamond(x$eff[i], x$low[i], x$upp[i],
                        x$sizes[i], x$range[1], x$range[2],
                        x$col.diamond[i], x$col.diamond.lines[i])
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
                       ymin.common, ymin.random, ymin.ref, ymax,
                       lwd, lty.common, lty.random, col.common, col.random,
                       xmin, xmax,
                       lower.equi, upper.equi,
                       lty.equi, col.equi,
                       fill.lower.equi, fill.upper.equi) {
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
  ## Add equivalence region(s)
  ##
  if (is.na(ref) & any(!is.na(lower.equi)) & any(!is.na(upper.equi)))
    ref.equi <- min(lower.equi, na.rm = TRUE) +
      0.5 * (max(upper.equi, na.rm = TRUE) - min(lower.equi, na.rm = TRUE))
  else
    ref.equi <- ref
  ##
  if (!is.na(ref.equi) && (xmin <= ref.equi & ref.equi <= xmax) &&
      any(!is.na(lower.equi))) {
    ##
    n.lower.equi <- sum(!is.na(lower.equi))
    n.fill.lower.equi <- length(fill.lower.equi)
    ##
    if (n.lower.equi < n.fill.lower.equi) {
      firstline <- FALSE
      lower.equi <- c(xmin, lower.equi, ref.equi)
    }
    else {
      firstline <- TRUE
      lower.equi <- c(lower.equi, ref.equi)
    }
    ##
    n.lo <- length(lower.equi[-1])
    ##
    for (i in seq_len(n.lo)) {
      if (!is.na(lower.equi[i]) && !is.na(lower.equi[i + 1]) &&
          ((xmin <= lower.equi[i] & lower.equi[i] <= lower.equi[i + 1]) &
           (lower.equi[i + 1] <= xmax))) {
        ##
        grid.polygon(x = unit(c(lower.equi[i], lower.equi[i + 1],
                                lower.equi[i + 1], lower.equi[i]), "native"),
                     y = unit(c(ymin.ref, ymin.ref, ymax, ymax),
                              "lines"),
                     gp = gpar(lwd = lwd, col = "transparent",
                               fill = fill.lower.equi[i]))
      }
    }
    ##
    for (i in seq_len(n.lo)) {
      if (!(i == 1 & !firstline) &
          !is.na(lower.equi[i]) && !is.na(lower.equi[i + 1]) &&
          ((xmin <= lower.equi[i] & lower.equi[i] <= lower.equi[i + 1]) &
           (lower.equi[i + 1] <= xmax)))
        grid.lines(x = unit(lower.equi[i], "native"),
                   y = unit(c(ymin.ref, ymax), "lines"),
                   gp = gpar(lwd = lwd, col = col.equi, lty = lty.equi))
    }
  }
  ##
  if (!is.na(ref.equi) && (xmin <= ref.equi & ref.equi <= xmax) &&
      any(!is.na(upper.equi))) {
    ##
    n.upper.equi <- sum(!is.na(upper.equi))
    n.fill.upper.equi <- length(fill.upper.equi)
    ##
    if (n.upper.equi < n.fill.upper.equi) {
      lastline <- FALSE
      upper.equi <- c(ref.equi, upper.equi, xmax)
    }
    else {
      lastline <- TRUE
      upper.equi <- c(ref.equi, upper.equi)
    }
    ##
    n.up <- length(upper.equi[-1])
    ##
    for (i in seq_len(n.up)) {
      if (!is.na(upper.equi[i]) && !is.na(upper.equi[i + 1]) &&
          ((xmin <= upper.equi[i] & upper.equi[i] <= upper.equi[i + 1]) &
           (upper.equi[i + 1] <= xmax))) {
        ##
        grid.polygon(x = unit(c(upper.equi[i], upper.equi[i + 1],
                                upper.equi[i + 1], upper.equi[i]), "native"),
                     y = unit(c(ymin.ref, ymin.ref, ymax, ymax),
                              "lines"),
                     gp = gpar(lwd = lwd, col = "transparent",
                               fill = fill.upper.equi[i]))
      }
    }
    ##
    for (i in seq_len(n.up)) {
      if ((i != n.up | (i == n.up & lastline)) &
          !is.na(upper.equi[i + 1]) &&
          (xmin <= upper.equi[i + 1] & upper.equi[i + 1] <= xmax))
          grid.lines(x = unit(upper.equi[i + 1], "native"),
                     y = unit(c(ymin.ref, ymax), "lines"),
                     gp = gpar(lwd = lwd, col = col.equi, lty = lty.equi))
    }
  }
  ##
  ## Reference line:
  ##
  if (!is.na(ref) && (xmin <= ref & ref <= xmax))
    grid.lines(x = unit(ref, "native"),
               y = unit(c(ymin.ref, ymax), "lines"),
               gp = gpar(lwd = lwd))
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
                      n.com, n.ran, n.prd) {
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
              rows = rows)
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
