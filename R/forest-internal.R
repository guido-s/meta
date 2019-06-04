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
      tval.min <- min(exp(x$range[1]), 1)
      tval.max <- max(exp(x$range[2]), 1)
      ##
      if (all(x1000 >= tval.min) &
          all(x1000 <= tval.max))
        label <- x1000
      else if (all(x100 >= tval.min) &
               all(x100 <= tval.max))
        label <- x100
      else if (all(x10 >= tval.min) &
               all(x10 <= tval.max))
        label <- x10
      else if (all(x5 >= tval.min) &
               all(x5 <= tval.max))
        label <- x5
      else if (all(x2 >= tval.min) &
               all(x2 <= tval.max))
        label <- x2
      else if (all(x1.5 >= tval.min) &
               all(x1.5 <= tval.max))
        label <- x1.5
      else if (all(x1.25 >= tval.min) &
               all(x1.25 <= tval.max))
        label <- x1.25
      else if (all(x1 >= tval.min) &
               all(x1 <= tval.max))
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
    if (!is.na(size)) {
      TElineCol <- if (size > 0 &&
                       (convertX(unit(TE, "native") + unit(0.5 * size, "lines"),
                                 "native", valueOnly = TRUE) > upper) &&
                       (convertX(unit(TE, "native") - unit(0.5 * size, "lines"),
                                 "native", valueOnly = TRUE) < lower))
                     col.inside
                   else
                     col
    }
  }
  ##
  if (!is.na(TE) && (TE >= min & TE <= max)) {
    if (!is.na(size) && size > 0) {
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
      lineCol <- if (size > 0 &&
                     (convertX(unit(TE, "native") + unit(0.5 * size, "lines"),
                               "native", valueOnly = TRUE) > upper) &&
                     (convertX(unit(TE, "native") - unit(0.5 * size, "lines"),
                               "native", valueOnly = TRUE) < lower))
                   col.inside
                 else
                   col
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
                       ref, TE.fixed, TE.random,
                       overall, comb.fixed, comb.random, prediction,
                       ymin.fixed, ymin.random, ymin.ref, ymax,
		       lwd, lty.fixed, lty.random, col.fixed, col.random,
                       min, max,
                       lower.equi, upper.equi,
                       lty.equi, col.equi, fill.equi) {
  ##
  if (min > max) {
    min <- x$range[2]
    max <- x$range[1]
  }
  else {
    min <- x$range[1]
    max <- x$range[2]
  }
  ##
  pushViewport(viewport(layout.pos.col = column, xscale = x$range))
  ##
  ## Add equivalence region
  ##
  if ((!is.na(lower.equi) && (min <= lower.equi & lower.equi <= max)) &
      (!is.na(upper.equi) && (min <= upper.equi & upper.equi <= max)))
    grid.polygon(x = unit(c(lower.equi, upper.equi,
                            upper.equi, lower.equi), "native"),
                 y = unit(c(ymin.ref, ymin.ref, ymax, ymax), "lines"),
                 gp = gpar(lwd = lwd, col = fill.equi, fill = fill.equi))
  ##
  if (!is.na(lower.equi) && (min <= lower.equi & lower.equi <= max))
    grid.lines(x = unit(lower.equi, "native"),
               y = unit(c(ymin.ref, ymax), "lines"),
               gp = gpar(lwd = lwd, col = col.equi, lty = lty.equi))
  ##
  if (!is.na(upper.equi) && (min <= upper.equi & upper.equi <= max))
    grid.lines(x = unit(upper.equi, "native"),
               y = unit(c(ymin.ref, ymax), "lines"),
               gp = gpar(lwd = lwd, col = col.equi, lty = lty.equi))
  ##
  ## Reference line:
  ##
  if (!is.na(ref) && (min <= ref & ref <= max))
    grid.lines(x = unit(ref, "native"),
               y = unit(c(ymin.ref, ymax), "lines"),
               gp = gpar(lwd = lwd))
  ##
  ## Line for fixed effect estimate:
  ##
  if (comb.fixed & overall & !is.na(TE.fixed))
    if (min <= TE.fixed & TE.fixed <= max)
      if (!is.null(lty.fixed))
        grid.lines(x = unit(TE.fixed, "native"),
                   y = unit(c(ymin.fixed, ymax), "lines"),
                   gp = gpar(lty = lty.fixed, lwd = lwd, col = col.fixed))
  ##
  ## Line for random effects estimate:
  ##
  if (comb.random & overall & !is.na(TE.random))
    if (min <= TE.random & TE.random <= max)
      if (!is.null(lty.random) & !is.na(TE.random))
        grid.lines(x = unit(TE.random, "native"),
                   y = unit(c(ymin.random, ymax), "lines"),
                   gp = gpar(lty = lty.random, lwd = lwd, col = col.random))
  ##
  popViewport()
  ##
  invisible(NULL)
}


formatcol <- function(x, y, rows, just = "right", settings,
                      fontfamily) {
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
  ## Fixed effect estimate:
  ##
  res$labels[[2]] <- textGrob(y[1],
                              x = xpos, just = just,
                              gp = gpar(
                                fontsize = settings$fs.fixed,
                                fontface = settings$ff.fixed,
                                fontfamily = fontfamily)
                              )
  ##
  ## Random effects estimate:
  ##
  res$labels[[3]] <- textGrob(y[2],
                              x = xpos, just = just,
                              gp = gpar(
                                fontsize = settings$fs.random,
                                fontface = settings$ff.random,
                                fontfamily = fontfamily)
                              )
  ##
  ## Prediction interval:
  ##
  res$labels[[4]] <- textGrob(y[3],
                              x = xpos, just = just,
                              gp = gpar(
                                fontsize = settings$fs.predict,
                                fontface = settings$ff.predict,
                                fontfamily = fontfamily)
                              )
  ##
  if (settings$by)
    for (i in 1:settings$n.by) {
      ##
      ## Fixed effect estimates:
      ##
      res$labels[[4 + i]] <- textGrob(y[3 + i],
                                      x = xpos, just = just,
                                      gp = 
                                        gpar(
                                          fontsize = settings$fs.fixed,
                                          fontface = settings$ff.fixed,
                                          fontfamily = fontfamily,
                                          col = settings$col.by)
                                      )
      ##
      ## Random effects estimates:
      ##
      res$labels[[4 + settings$n.by + i]] <- textGrob(y[3 + settings$n.by + i],
                                             x = xpos, just = just,
                                             gp = 
                                               gpar(
                                                 fontsize = settings$fs.random,
                                                 fontface = settings$ff.random,
                                                 fontfamily = fontfamily,
                                                 col = settings$col.by)
                                             )
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
