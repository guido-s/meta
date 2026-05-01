forest_meta_internal <- function(
    new, nrow, x1, spacing, 
    yHeadadd,
    #
    cols, cols.new, newcols, by,
    #
    lsel, rsel,
    leftcols, leftcols.new,
    leftlabs, leftlabs.new,
    rightcols, rightcols.new,
    rightlabs, rightlabs.new,
    #
    label.e.attach, label.c.attach,
    col.label.e, col.label.c,
    just.label.e, just.label.c, just.addcols.left, just.addcols.right,
    #
    RoB.available, rob, col.rob, rob.attach,
    rob.categories, rob.symbols, rob.col,
    fs.rob.symbols, ff.rob.symbols,
    #
    clines, jama, col.jama.line, bmj, revman5,
    #
    ev.n.bin, m.s.n.cont, ev.n.prop,
    #
    fontfamily, fs.head, ff.head,
    #
    col.forest, col.lines, col.label, col.label.left, col.label.right,
    col.common, col.random, col.cid,
    #
    col.add.TE, col.add.ci, col.add.cluster, col.add.cor,
    col.add.cycles, col.add.effect, col.add.effect.ci,
    col.add.event.c, col.add.event.e, col.add.event.n.c,
    col.add.event.n.e, col.add.mean.c, col.add.mean.e,
    col.add.mean.sd.n.c, col.add.mean.sd.n.e, col.add.n.c,
    col.add.n.e, col.add.sd.c, col.add.sd.e, col.add.seTE,
    col.add.studlab, col.add.time.c, col.add.time.e,
    col.add.w.common, col.add.w.random,
    #
    ref, TE.common, TE.random,
    overall, common, random, prediction,
    #
    k.all, n.com, n.ran, n.prd,
    #
    ymin.common, ymin.random, ymin.ref, ymax, ymax.ref,
    #
    header.line, header.line.pos, col.header.line,
    #
    lwd, lty.common, lty.random, xlim, avail.xlim,
    #
    cid.below.null, cid.above.null, lty.cid,
    fill.cid.below.null, fill.cid.above.null, fill,
    #
    xlab, xlab.add, newline.xlab, xlab.pos, xlab.ypos, fs.xlab, ff.xlab,
    #
    bottom.lr, smlab1, smlab2, newline.smlab, newline.lr,
    print.label, ll1, ll2, newline.ll, lr1, lr2,
    y.bottom.lr, fs.lr, ff.lr,
    #
    yS, log.xaxis, at, label, fs.axis, ff.axis,
    #
    newline.TE, newline.ci, newline.cluster, newline.cor,
    newline.cycles, newline.effect,newline.effect.ci,
    newline.event.c, newline.event.e, newline.event.n.c,
    newline.event.n.e, newline.mean.c, newline.mean.e,
    newline.mean.sd.n.c, newline.mean.sd.n.e, newline.n.c,
    newline.n.e, newline.sd.c, newline.sd.e, newline.seTE,
    newline.studlab, newline.time.c, newline.time.e,
    newline.w.common, newline.w.random,
    #
    metabin, metacont, metacor, metagen, metainc, metamean,
    metaprop, metarate, metabind, metacum, metainf, metamerge, meta,
    #
    addrow, addrows.below.overall,
    #
    colgap, colgap.left, colgap.right,
    colgap.studlab, colgap.forest, colgap.forest.left, colgap.forest.right,
    #
    studlab, TE.format, seTE.format, cluster.format, cycles.format,
    effect.format, ci.format, effect.ci.format) {
    
  if (new)
    grid.newpage()
  #
  pushViewport(
    viewport(
      layout =
        grid.layout(
          nrow, length(x1), widths = x1,
          heights = unit(spacing, "lines"))))
  #
  # Left side of forest plot
  #
  j <- 1
  #
  if (lsel) {
    #
    # Add text for label.e and label.c (if position was specified by the user)
    #
    if (!is.na(yHeadadd)) {
      if (!is.null(label.e.attach)) {
        vars.e <- paste0("col.", label.e.attach)
        if (all(vars.e %in% leftcols)) {
          id.e <- seq_along(leftcols)[leftcols %in% vars.e]
          id.e <- 2 * (range(id.e) - 1) + 1
          id.e <- seq(min(id.e), max(id.e))
          #
          add.text(col.label.e, id.e)
        }
      }
      #
      if (!is.null(label.c.attach)) {
        vars.c <- paste0("col.", label.c.attach)
        if (all(vars.c %in% leftcols)) {
          id.c <- seq_along(leftcols)[leftcols %in% vars.c]
          id.c <- 2 * (range(id.c) - 1) + 1
          id.c <- seq(min(id.c), max(id.c))
          #
          add.text(col.label.c, id.c)
        }
      }
    }
    #
    for (i in seq_along(leftcols)) {
      add.text(cols[[leftcols[i]]], j)
      #
      if (!is.na(yHeadadd)) {
        if (is.null(label.e.attach)) {
          if (metabin) {
            if (leftcols[i] == "col.n.e" & just.label.e == "right")
              add.text(col.label.e, j)
            else if (leftcols[i] == "col.event.e" &
                     just.label.e %in% c("left", "center"))
              add.text(col.label.e, j)
          }
          else if (metacont) {
            if (leftcols[i] == "col.sd.e" & just.label.e == "right")
              add.text(col.label.e, j)
            else if (leftcols[i] == "col.mean.e" &
                     just.label.e %in% c("left", "center"))
              add.text(col.label.e, j)
          }
          else if (metainc) {
            if (leftcols[i] == "col.time.e" & just.label.e == "right")
              add.text(col.label.e, j)
            else if (leftcols[i] == "col.event.e" &
                     just.label.e %in% c("left", "center"))
              add.text(col.label.e, j)
          }
          else if (metamean) {
            if (revman5 & leftcols[i] == "col.n.e" &
                just.label.e == "right")
              add.text(col.label.e, j)
            else if (!revman5 & leftcols[i] == "col.sd.e" &
                     just.label.e == "right")
              add.text(col.label.e, j)
            else if (leftcols[i] == "col.sd.e" &
                     just.label.e %in% c("left", "center"))
              add.text(col.label.e, j)
          } 
        }
        #
        if (is.null(label.c.attach)) {
          if (metabin) {
            if (leftcols[i] == "col.n.c" & just.label.c == "right")
              add.text(col.label.c, j)
            else if (leftcols[i] == "col.event.c" &
                     just.label.c %in% c("left", "center"))
              add.text(col.label.c, j)
          }
          else if (metacont) {
            if (leftcols[i] == "col.sd.c" & just.label.c == "right")
              add.text(col.label.c, j)
            else if (leftcols[i] == "col.mean.c" &
                     just.label.c %in% c("left", "center"))
              add.text(col.label.c, j)
          }
          else if (metainc) {
            if (leftcols[i] == "col.time.c" & just.label.c == "right")
              add.text(col.label.c, j)
            else if (leftcols[i] == "col.event.c" &
                     just.label.c %in% c("left", "center"))
              add.text(col.label.c, j)
          }
        }
        #
        if (newline.studlab & leftcols[i] == "col.studlab")
          add.text(col.add.studlab, j)
        if (newline.effect & leftcols[i] == "col.effect")
          add.text(col.add.effect, j)
        if (newline.ci & leftcols[i] == "col.ci")
          add.text(col.add.ci, j)
        if (newline.effect.ci & leftcols[i] == "col.effect.ci")
          add.text(col.add.effect.ci, j)
        if (newline.event.n.e & leftcols[i] == "col.event.n.e")
          add.text(col.add.event.n.e, j)
        if (newline.event.n.c & leftcols[i] == "col.event.n.c")
          add.text(col.add.event.n.c, j)
        if (newline.mean.sd.n.e & leftcols[i] == "col.mean.sd.n.e")
          add.text(col.add.mean.sd.n.e, j)
        if (newline.mean.sd.n.c & leftcols[i] == "col.mean.sd.n.c")
          add.text(col.add.mean.sd.n.c, j)
        if (newline.w.common & leftcols[i] == "col.w.common")
          add.text(col.add.w.common, j)
        if (newline.w.random & leftcols[i] == "col.w.random")
          add.text(col.add.w.random, j)
        if (newline.TE & leftcols[i] == "col.TE")
          add.text(col.add.TE, j)
        if (newline.seTE & leftcols[i] == "col.seTE")
          add.text(col.add.seTE, j)
        if (newline.cluster & leftcols[i] == "col.cluster")
          add.text(col.add.cluster, j)
        if (newline.cycles & leftcols[i] == "col.cycles")
          add.text(col.add.cycles, j)
        if (newline.n.e & leftcols[i] == "col.n.e")
          add.text(col.add.n.e, j)
        if (newline.n.c & leftcols[i] == "col.n.c")
          add.text(col.add.n.c, j)
        if (newline.event.e & leftcols[i] == "col.event.e")
          add.text(col.add.event.e, j)
        if (newline.event.c & leftcols[i] == "col.event.c")
          add.text(col.add.event.c, j)
        if (newline.mean.e & leftcols[i] == "col.mean.e")
          add.text(col.add.mean.e, j)
        if (newline.mean.c & leftcols[i] == "col.mean.c")
          add.text(col.add.mean.c, j)
        if (newline.sd.e & leftcols[i] == "col.sd.e")
          add.text(col.add.sd.e, j)
        if (newline.sd.c & leftcols[i] == "col.sd.c")
          add.text(col.add.sd.c, j)
        if (newline.cor & leftcols[i] == "col.cor")
          add.text(col.add.cor, j)
        if (newline.time.e & leftcols[i] == "col.time.e")
          add.text(col.add.time.e, j)
        if (newline.time.c & leftcols[i] == "col.time.c")
          add.text(col.add.time.c, j)
        #
        # Add text in first line of forest plot for new columns
        #
        if (newcols)
          if (length(leftcols.new) > 0 &
              leftcols[i] %in% paste0("col.", leftcols.new)) {
            sel <- paste0("col.", leftcols.new) == leftcols[i]
            #
            # Check for "\n" in label of new column
            #
            clines <- twolines(leftlabs.new[sel], leftcols[i])
            #
            just.new <- just.addcols.left[sel]
            #
            if (just.new == "left")
              xpos.new <- 0
            else if (just.new == "center")
              xpos.new <- 0.5
            else if (just.new == "right")
              xpos.new <- 1
            #
            # Add first line
            #
            if (clines$newline)
              add.text(tgl(clines$top, xpos.new, just.new, fs.head, ff.head,
                           fontfamily), j)
          }
      }
      #
      j <- j + 2
    }
  }
  #
  # Produce forest plot
  #
  draw.lines(col.forest, j,
             ref, TE.common, unique(TE.random),
             overall, common, random, prediction,
             ymin.common, ymin.random, ymin.ref,
             ymax + 0.5 * header.line * addrow,
             ymax.ref + 0.5 * header.line * addrow,
             lwd, lty.common, lty.random, col.common, col.random,
             xlim[1], xlim[2],
             cid.below.null, cid.above.null, lty.cid, col.cid,
             fill.cid.below.null, fill.cid.above.null,
             fill,
             col.lines)
  #
  draw.axis(col.forest, j, yS, log.xaxis, at, label,
            fs.axis, ff.axis, fontfamily, lwd,
            xlim, avail.xlim,
            col.lines, col.label)
  #
  if (bottom.lr) {
    add.text(smlab1, j, xscale = col.forest$range)
    #
    if (newline.smlab)
      add.text(smlab2, j, xscale = col.forest$range)
  }
  #
  if (print.label) {
    if (!bottom.lr) {
      if (!is.na(ref)) {
        add.text(ll1, j, xscale = col.forest$range)
        #
        if (newline.ll)
          add.text(ll2, j, xscale = col.forest$range)
        #
        add.text(lr1, j, xscale = col.forest$range)
        #
        if (newline.lr)
          add.text(lr2, j, xscale = col.forest$range)
      }
    }
    else {
      add.label(ll1, j,
                if (bmj)
                  unit(xlim[1], "native")
                else
                  unit(ref - (xlim[2] - xlim[1]) / 30, "native"),
                unit(y.bottom.lr, "lines"),
                if (bmj) "left" else "right",
                fs.lr, ff.lr, col.label.left, fontfamily,
                xscale = col.forest$range)
      #
      if (newline.ll)
        add.label(ll2, j,
                  if (bmj)
                    unit(xlim[1], "native")
                  else
                    unit(ref - (xlim[2] - xlim[1]) / 30, "native"),
                  unit(y.bottom.lr - 1, "lines"),
                  if (bmj) "left" else "right",
                  fs.lr, ff.lr, col.label.left, fontfamily,
                  xscale = col.forest$range)
      #
      add.label(lr1, j,
                if (bmj)
                  unit(xlim[2], "native")
                else
                  unit(ref + (xlim[2] - xlim[1]) / 30, "native"),
                unit(y.bottom.lr, "lines"),
                if (bmj) "right" else "left",
                fs.lr, ff.lr, col.label.right, fontfamily,
                xscale = col.forest$range)
      #
      if (newline.lr)
        add.label(lr2, j,
                  if (bmj)
                    unit(xlim[2], "native")
                  else
                    unit(ref + (xlim[2] - xlim[1]) / 30, "native"),
                  unit(y.bottom.lr - 1, "lines"),
                  if (bmj) "right" else "left",
                  fs.lr, ff.lr, col.label.right, fontfamily,
                  xscale = col.forest$range)
    }
  }
  #
  add.xlab(col.forest, j, xlab, xlab.add, newline.xlab,
           xlab.pos, xlab.ypos, fs.xlab, ff.xlab,
           fontfamily)
  #
  draw.forest(col.forest, j)
  #
  j <- j + 2
  #
  #
  # Right side of forest plot
  #
  #
  if (rsel | RoB.available) {
    #
    # Add text for label.e and label.c (if position was specified by the user)
    #
    if (!is.na(yHeadadd)) {
      if (!is.null(label.e.attach)) {
        vars.e <- paste0("col.", label.e.attach)
        if (all(vars.e %in% rightcols)) {
          id.e <- seq_along(rightcols)[rightcols %in% vars.e]
          id.e <- 2 * (range(id.e) - 1)
          id.e <- seq(min(id.e), max(id.e))
          #
          add.text(col.label.e, j + id.e)
        }
      }
      #
      if (!is.null(label.c.attach)) {
        vars.c <- paste0("col.", label.c.attach)
        if (all(vars.c %in% rightcols)) {
          id.c <- seq_along(rightcols)[rightcols %in% vars.c]
          id.c <- 2 * (range(id.c) - 1)
          id.c <- seq(min(id.c), max(id.c))
          #
          add.text(col.label.c, j + id.c)
        }
      }
    }
    #
    i.rob <- 0
    #
    for (i in seq_along(rightcols)) {
      if (substring(rightcols[i], 1, 8) == "col.RoB.") {
        i.rob <- i.rob + 1
        #
        add.rob(cols[[rightcols[i]]], j, 0.85, fs.rob.symbols, ff.rob.symbols,
                fontfamily,
                rob[[rightcols[[i]]]],
                rob.categories[[i.rob]], rob.symbols[[i.rob]], rob.col[[i.rob]])
      }
      else
        add.text(cols[[rightcols[i]]], j)
      #
      if (!is.na(yHeadadd)) {
        if (is.null(label.e.attach)) {
          if (metabin) {
            if (rightcols[i] == "col.n.e" & just.label.e == "right")
              add.text(col.label.e, j)
            else if (rightcols[i] == "col.event.e" &
                     just.label.e %in% c("left", "center"))
              add.text(col.label.e, j)
          }
          else if (metacont) {
            if (rightcols[i] == "col.sd.e" & just.label.e == "right")
              add.text(col.label.e, j)
            else if (rightcols[i] == "col.mean.e" &
                     just.label.e %in% c("left", "center"))
              add.text(col.label.e, j)
          }
          else if (metainc) {
            if (rightcols[i] == "col.time.e" & just.label.e == "right")
              add.text(col.label.e, j)
            else if (rightcols[i] == "col.event.e" &
                     just.label.e %in% c("left", "center"))
              add.text(col.label.e, j)
          }
          else if (metamean) {
            if (revman5 & rightcols[i] == "col.n.e" & just.label.e == "right")
              add.text(col.label.e, j)
            else if (!revman5 & rightcols[i] == "col.sd.e" &
                     just.label.e == "right")
              add.text(col.label.e, j)
            else if (rightcols[i] == "col.sd.e" &
                     just.label.e %in% c("left", "center"))
              add.text(col.label.e, j)
          }
        }
        #
        if (is.null(label.c.attach)) {
          if (metabin) {
            if (rightcols[i] == "col.n.c" & just.label.c == "right")
              add.text(col.label.c, j)
            else if (rightcols[i] == "col.event.c" &
                     just.label.c %in% c("left", "center"))
              add.text(col.label.c, j)
          }
          else if (metacont) {
            if (rightcols[i] == "col.sd.c" & just.label.c == "right")
              add.text(col.label.c, j)
            else if (rightcols[i] == "col.mean.c" &
                     just.label.c %in% c("left", "center"))
              add.text(col.label.c, j)
          }
          else if (metainc) {
            if (rightcols[i] == "col.time.c" & just.label.c == "right")
              add.text(col.label.c, j)
            else if (rightcols[i] == "col.event.c" &
                     just.label.c %in% c("left", "center"))
              add.text(col.label.c, j)
          }
        }
        #
        if (!is.null(rob.attach)) {
          if (rightcols[i] == rob.attach)
            add.text(col.rob, j)
        }
        #
        if (newline.studlab & rightcols[i] == "col.studlab")
          add.text(col.add.studlab, j)
        if (newline.effect & rightcols[i] == "col.effect")
          add.text(col.add.effect, j)
        if (newline.ci & rightcols[i] == "col.ci")
          add.text(col.add.ci, j)
        if (newline.effect.ci & rightcols[i] == "col.effect.ci")
          add.text(col.add.effect.ci, j)
        if (newline.mean.sd.n.e & rightcols[i] == "col.mean.sd.n.e")
          add.text(col.add.mean.sd.n.e, j)
        if (newline.mean.sd.n.c & rightcols[i] == "col.mean.sd.n.c")
          add.text(col.add.mean.sd.n.c, j)
        if (newline.w.common & rightcols[i] == "col.w.common")
          add.text(col.add.w.common, j)
        if (newline.w.random & rightcols[i] == "col.w.random")
          add.text(col.add.w.random, j)
        if (newline.TE & rightcols[i] == "col.TE")
          add.text(col.add.TE, j)
        if (newline.seTE & rightcols[i] == "col.seTE")
          add.text(col.add.seTE, j)
        if (newline.cluster & rightcols[i] == "col.cluster")
          add.text(col.add.cluster, j)
        if (newline.cycles & rightcols[i] == "col.cycles")
          add.text(col.add.cycles, j)
        if (newline.n.e & rightcols[i] == "col.n.e")
          add.text(col.add.n.e, j)
        if (newline.n.c & rightcols[i] == "col.n.c")
          add.text(col.add.n.c, j)
        if (newline.event.e & rightcols[i] == "col.event.e")
          add.text(col.add.event.e, j)
        if (newline.event.c & rightcols[i] == "col.event.c")
          add.text(col.add.event.c, j)
        if (newline.mean.e & rightcols[i] == "col.mean.e")
          add.text(col.add.mean.e, j)
        if (newline.mean.c & rightcols[i] == "col.mean.c")
          add.text(col.add.mean.c, j)
        if (newline.sd.e & rightcols[i] == "col.sd.e")
          add.text(col.add.sd.e, j)
        if (newline.sd.c & rightcols[i] == "col.sd.c")
          add.text(col.add.sd.c, j)
        if (newline.cor & rightcols[i] == "col.cor")
          add.text(col.add.cor, j)
        if (newline.time.e & rightcols[i] == "col.time.e")
          add.text(col.add.time.e, j)
        if (newline.time.c & rightcols[i] == "col.time.c")
          add.text(col.add.time.c, j)
        #
        # Add text in first line of forest plot for new columns
        #
        if (newcols)
          if (length(rightcols.new) > 0 &
              rightcols[i] %in% paste0("col.", rightcols.new)) {
            sel <- paste0("col.", rightcols.new) == rightcols[i]
            #
            # Check for "\n" in label of new column
            #
            clines <- twolines(rightlabs.new[sel], rightcols[i])
            #
            just.new <- just.addcols.right[sel]
            #
            if (just.new == "left")
              xpos.new <- 0
            else if (just.new == "center")
              xpos.new <- 0.5
            else if (just.new == "right")
              xpos.new <- 1
            #
            # Add first line
            #
            if (clines$newline)
              add.text(tgl(clines$top, xpos.new, just.new,
                           fs.head, ff.head, fontfamily), j)
          }
      }
      #
      j <- j + 2
    }
  }
  #
  # Add header line
  #
  if (jama)
    hcols <- lsel * 2 * length(leftcols)
  else
    hcols <- lsel * 2 * length(leftcols) + 1 + rsel * 2 * length(rightcols)
  #
  if (ev.n.bin) {
    sel1 <- grep("col.event.n.e",
                 c(leftcols, if (all(rightcols != "col.")) rightcols))
    sel2 <- grep("col.event.n.c",
                 c(leftcols, if (all(rightcols != "col.")) rightcols))
    #
    if (length(sel1) > 0 & length(sel2) > 0) {
      if (sel1 > sel2) {
        sel3 <- sel2
        sel2 <- sel1
        sel1 <- sel2
      }
      #
      for (i in seq(2 * sel1 - 1, 2 * sel2 - 1)) {
        pushViewport(
          viewport(
            layout.pos.col = i,
            xscale = col.forest$range))
        #
        grid.lines(x = unit(0:1, "npc"),
                   y = unit(nrow - 1.5 + 0.5 * addrow, "lines"),
                   gp = gpar(lwd = lwd))
        #
        popViewport()
      }
    }
  }
  #
  if (m.s.n.cont) {
    sel1 <- grep("col.mean.sd.n.e",
                 c(leftcols, if (all(rightcols != "col.")) rightcols))
    sel2 <- grep("col.mean.sd.n.c",
                 c(leftcols, if (all(rightcols != "col.")) rightcols))
    #
    if (length(sel1) > 0 & length(sel2) > 0) {
      if (sel1 > sel2) {
        sel3 <- sel2
        sel2 <- sel1
        sel1 <- sel2
      }
      #
      for (i in seq(2 * sel1 - 1, 2 * sel2 - 1)) {
        pushViewport(
          viewport(
            layout.pos.col = i,
            xscale = col.forest$range))
        #
        grid.lines(x = unit(0:1, "npc"),
                   y = unit(nrow - 1.5 + 0.5 * addrow, "lines"),
                   gp = gpar(lwd = lwd))
        #
        popViewport()
      }
    }
  }
  #
  if (header.line) {
    if (header.line.pos == "both") {
      for (i in seq_len(hcols)) {
        pushViewport(
          viewport(
            layout.pos.col = i,
            xscale = col.forest$range))
        #
        grid.lines(x = unit(0:1, "npc"),
                   y = unit(nrow + 0.5 * addrow, "lines"),
                   gp = gpar(lwd = lwd, col = col.header.line))
        #
        popViewport()
      }
    }
    #
    for (i in seq_len(hcols)) {
      pushViewport(
        viewport(
          layout.pos.col = i,
          xscale = col.forest$range))
      #
      grid.lines(x = unit(0:1, "npc"),
                 y = unit(ymax + 0.5 * addrow, "lines"),
                 gp = gpar(lwd = lwd, col = col.header.line))
      #
      popViewport()
    }
  }
  #
  # Add JAMA lines
  #
  if (jama & header.line & !by) {
    for (i in seq_len(hcols)) {
      pushViewport(
        viewport(
          layout.pos.col = i,
          xscale = col.forest$range))
      #
      for (j in seq_len(k.all + n.com * common + n.ran * random +
                        n.prd * prediction))
        grid.lines(x = unit(0:1, "npc"),
                   y = unit(ymax + 0.5 * addrow - j, "lines"),
                   gp = gpar(lwd = 0.5 * lwd, col = col.jama.line))
      #
      popViewport()
    }
  }
  #
  popViewport()
  
  
  res <- list(xlim = xlim, addrows.below.overall = addrows.below.overall,
              #
              colgap = colgap,
              colgap.left = colgap.left,
              colgap.right = colgap.right,
              colgap.studlab = colgap.studlab,
              colgap.forest = colgap.forest,
              colgap.forest.left = colgap.forest.left,
              colgap.forest.right = colgap.forest.right,
              #
              studlab = studlab,
              TE.format = TE.format,
              seTE.format = seTE.format,
              cluster.format = cluster.format,
              cycles.format = cycles.format,
              effect.format = effect.format,
              ci.format = ci.format,
              effect.ci.format = effect.ci.format)
  #
  if (ev.n.bin | ev.n.prop)
    res$effect.ci.format <- effect.ci.format
  #
  if (ev.n.bin)
    res$effect.ci.format <- effect.ci.format
  #
  if (length(cols.new) > 0) {
    for (i in names(cols.new))
      res[[i]] <- cols.new[[i]]
  }
  #
  res$leftcols <- leftcols
  res$leftlabs <- leftlabs
  res$rightcols <- rightcols
  res$rightlabs <- rightlabs
  #
  invisible(res)
}
