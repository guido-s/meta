forest.meta <- function(x,
                        sortvar,
                        studlab = TRUE,
                        ##
                        layout = gs("layout"),
                        ##
                        comb.fixed = x$comb.fixed,
                        comb.random = x$comb.random,
                        overall = TRUE,
                        text.fixed = NULL,
                        text.random = NULL,
                        lty.fixed = 2, lty.random = 3,
                        ##
                        prediction = x$prediction,
                        text.predict = NULL,
                        ##
                        print.subgroup.labels = TRUE,
                        bylab = x$bylab,
                        print.byvar = x$print.byvar,
                        byseparator = gs("byseparator"),
                        text.fixed.w = text.fixed,
                        text.random.w = text.random,
                        bysort = FALSE,
                        ##
                        pooled.totals = comb.fixed | comb.random,
                        pooled.events = FALSE, pooled.times = FALSE,
                        ##
                        study.results = TRUE,
                        ##
                        xlab = "", xlab.pos,
                        smlab = NULL, smlab.pos, xlim = "symmetric",
                        ##
                        allstudies = TRUE,
                        weight.study,
                        weight.subgroup,
                        pscale = x$pscale,
                        irscale = x$irscale, irunit = x$irunit,
                        ##
                        ref = ifelse(backtransf & is.relative.effect(x$sm), 1, 0),
                        ##
                        leftcols = NULL, rightcols = NULL,
                        leftlabs = NULL, rightlabs = NULL,
                        ##
                        lab.e = x$label.e,
                        lab.c = x$label.c,
                        ##
                        lab.e.attach.to.col = NULL,
                        lab.c.attach.to.col = NULL,
                        ##
                        label.right = x$label.right,
                        label.left = x$label.left,
                        bottom.lr = TRUE,
                        ##
                        lab.NA = ".", lab.NA.effect = "",
                        ##
                        lwd = 1,
                        ##
                        at = NULL,
                        label = TRUE,
                        ##
                        type.study = "square",
                        type.fixed = "diamond",
                        type.random = type.fixed,
                        type.subgroup = ifelse(study.results, "diamond", "square"),
                        ##
                        col.study = "black",
                        col.square = "gray",
                        col.square.lines = col.square,
                        col.inside = "white",
                        ##
                        col.diamond = "gray",
                        col.diamond.fixed = col.diamond,
                        col.diamond.random = col.diamond,
                        col.diamond.lines = "black",
                        col.diamond.lines.fixed = col.diamond.lines,
                        col.diamond.lines.random = col.diamond.lines,
                        ##
                        col.inside.fixed = col.inside,
                        col.inside.random = col.inside,
                        ##
                        col.predict = "red",
                        col.predict.lines = "black",
                        ##
                        col.by = "darkgray",
                        ##
                        col.label.right = "black",
                        col.label.left = "black",
                        ##
                        hetstat = print.I2 | print.tau2 | print.Q | print.pval.Q | print.Rb,
                        overall.hetstat = overall & hetstat,
                        hetlab = "Heterogeneity: ",
                        print.I2 = comb.fixed | comb.random,
                        print.I2.ci = FALSE,
                        print.tau2 = comb.fixed | comb.random,
                        print.Q = FALSE,
                        print.pval.Q = comb.fixed | comb.random,
                        print.Rb = FALSE,
                        print.Rb.ci = FALSE,
                        text.subgroup.nohet = "not applicable",
                        ##
                        test.overall = gs("test.overall"),
                        test.overall.fixed = comb.fixed & overall & test.overall,
                        test.overall.random = comb.random & overall & test.overall,
                        label.test.overall.fixed,
                        label.test.overall.random,
                        ##
                        print.zval = TRUE,
                        ##
                        test.subgroup,
                        test.subgroup.fixed,
                        test.subgroup.random,
                        print.Q.subgroup = TRUE,
                        label.test.subgroup.fixed,
                        label.test.subgroup.random,
                        ##
                        test.effect.subgroup,
                        test.effect.subgroup.fixed,
                        test.effect.subgroup.random,
                        label.test.effect.subgroup.fixed,
                        label.test.effect.subgroup.random,
                        ##
                        fontsize = 12,
                        fs.heading = fontsize,
                        fs.fixed,
                        fs.random,
                        fs.predict,
                        fs.fixed.labels,
                        fs.random.labels,
                        fs.predict.labels,
                        fs.study = fontsize,
                        fs.study.labels = fs.study,
                        fs.hetstat,
                        fs.test.overall,
                        fs.test.subgroup,
                        fs.test.effect.subgroup,
                        fs.axis = fontsize,
                        fs.smlab = fontsize,
                        fs.xlab = fontsize,
                        fs.lr = fontsize,
                        ##
                        ff.heading = "bold",
                        ff.fixed,
                        ff.random,
                        ff.predict,
                        ff.fixed.labels,
                        ff.random.labels,
                        ff.predict.labels,
                        ff.study = "plain",
                        ff.study.labels = ff.study,
                        ff.hetstat,
                        ff.test.overall,
                        ff.test.subgroup,
                        ff.test.effect.subgroup,
                        ff.axis = "plain",
                        ff.smlab = "bold",
                        ff.xlab = "plain",
                        ff.lr = "plain",
                        ##
                        squaresize = 0.8 / spacing,
                        ##
                        plotwidth = if (layout != "JAMA") "6cm" else "8cm",
                        colgap = "2mm",
                        colgap.left = colgap,
                        colgap.right = colgap,
                        colgap.studlab = colgap.left,
                        colgap.forest = colgap,
                        colgap.forest.left = colgap.forest,
                        colgap.forest.right = colgap.forest,
                        ##
                        calcwidth.pooled = TRUE,
                        calcwidth.fixed = calcwidth.pooled,
                        calcwidth.random = calcwidth.pooled,
                        calcwidth.predict = FALSE,
                        calcwidth.hetstat = FALSE,
                        calcwidth.tests  = FALSE,
                        ##
                        just = if (layout != "JAMA") "right" else "left",
                        just.studlab = "left",
                        just.addcols = "center",
                        just.addcols.left = just.addcols,
                        just.addcols.right = just.addcols,
                        ##
                        spacing = 1,
                        addrow,
                        addrow.overall,
                        addrow.subgroups,
                        ##
                        new = TRUE,
                        ##
                        backtransf = x$backtransf,
                        digits = gs("digits.forest"),
                        digits.se = gs("digits.se"),
                        digits.zval = gs("digits.zval"),
                        digits.pval = max(gs("digits.pval") - 2, 2),
                        digits.pval.Q = max(gs("digits.pval.Q") - 2, 2),
                        digits.Q = gs("digits.Q"),
                        digits.tau2 = gs("digits.tau2"),
                        digits.I2 = max(gs("digits.I2") - 1, 0),
                        digits.weight = gs("digits.weight"),
                        ##
                        digits.mean = NULL,
                        digits.sd = NULL,
                        digits.cor = NULL,
                        digits.time = NULL,
                        ##
                        scientific.pval = gs("scientific.pval"),
                        ##
                        col.i = col.study,
                        weight = weight.study,
                        ...) {
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(x, "meta")
  x.name <- deparse(substitute(x))
  x <- updateversion(x)
  ##
  K.all <- length(x$TE)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  mf <- match.call()
  error <- try(sortvar <- eval(mf[[match("sortvar", names(mf))]],
                               as.data.frame(x, stringsAsFactors = FALSE),
                               enclos = sys.frame(sys.parent())),
               silent = TRUE)
  if (class(error) == "try-error") {
    xd <- x$data
    sortvar <- eval(mf[[match("sortvar", names(mf))]],
                    xd, enclos = NULL)
    if (!is.null(x$data$.subset))
      sortvar <- sortvar[x$data$.subset]
  }
  sort <- !is.null(sortvar)
  if (sort && (length(sortvar) != K.all))
    stop("Number of studies in object 'x' and argument 'sortvar' have different length.")
  if (!sort)
    sortvar <- 1:K.all
  ##
  slab <- TRUE
  if (length(studlab) == 1 & is.logical(studlab))
    if (studlab == FALSE) {
      studlab <- rep("", K.all)
      slab <- FALSE
    }
    else studlab <- x$studlab
  ##
  if (length(studlab) != K.all)
    stop("Number of studies in object 'x' and argument 'studlab' have different length.")
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(overall)
  ##
  if (!is.null(lty.fixed))
    chknumeric(lty.fixed)
  if (!is.null(lty.random))
    chknumeric(lty.random)
  chklogical(prediction)
  chklogical(print.subgroup.labels)
  if (!is.null(print.byvar))
    chklogical(print.byvar)
  chkchar(byseparator)
  chklogical(bysort)
  chklogical(pooled.totals)
  chklogical(pooled.events)
  chklogical(pooled.times)
  chklogical(study.results)
  ## chknumeric(xlab.pos) ??
  ## chknumeric(smlab.pos) ??
  chklogical(allstudies)
  if (!is.null(pscale))
    chknumeric(pscale, single = TRUE)
  else
    pscale <- 1
  if (!is.null(irscale))
    chknumeric(irscale, single = TRUE)
  else
    irscale <- 1
  chknumeric(ref)
  ##
  layout <- setchar(layout, c("meta", "RevMan5", "JAMA", "subgroup"))
  if (layout == "subgroup" & is.null(x$byvar)) {
    warning("Argument 'layout' set to \"meta\" (default) as no subgroup analysis was conducted.")
    layout <- "meta"
  }
  if (layout == "subgroup") {
    if (missing(type.subgroup))
      type.subgroup <- "square"
    ##
    if (missing(pooled.totals))
      pooled.totals <- FALSE
  }
  revman5 <- layout == "RevMan5"
  jama <- layout == "JAMA"
  revman5.jama <- revman5 | jama
  ##
  type.study <- setchar(type.study, c("square", "diamond"))
  type.fixed <- setchar(type.fixed, c("square", "diamond"))
  type.random <- setchar(type.random, c("square", "diamond"))
  type.subgroup <- setchar(type.subgroup, c("square", "diamond"))
  ##
  if (missing(weight.subgroup))
    weight.subgroup <- ifelse(type.subgroup == "square",
                              "weight", "same")
  weight.subgroup <- setchar(weight.subgroup,
                             c("weight", "same"))
  ##
  chklogical(bottom.lr)
  chkchar(lab.NA)
  chkchar(lab.NA.effect)
  if (!is.null(at))
    chknumeric(at)
  chkchar(col.diamond)
  chkchar(col.diamond.fixed)
  chkchar(col.diamond.random)
  chkchar(col.diamond.lines)
  chkchar(col.diamond.lines.fixed)
  chkchar(col.diamond.lines.random)
  chkchar(col.inside.fixed)
  chkchar(col.inside.random)
  chkchar(col.predict)
  chkchar(col.predict.lines)
  chklogical(print.I2)
  chklogical(print.I2.ci)
  chklogical(print.tau2)
  chklogical(print.Q)
  chklogical(print.pval.Q)
  chklogical(print.Rb)
  chklogical(print.Rb.ci)
  if (!is.logical(text.subgroup.nohet))
    chkchar(text.subgroup.nohet)
  else if (text.subgroup.nohet)
    text.subgroup.nohet <- "not applicable"
  chklogical(print.zval)
  chklogical(hetstat)
  chklogical(overall.hetstat)
  chklogical(test.overall.fixed)
  chklogical(test.overall.random)
  if (!missing(test.subgroup.fixed))
    chklogical(test.subgroup.fixed)
  if (!missing(test.subgroup.random))
    chklogical(test.subgroup.random)
  chklogical(print.Q.subgroup)
  if (!missing(test.effect.subgroup.fixed))
    chklogical(test.effect.subgroup.fixed)
  if (!missing(test.effect.subgroup.random))
    chklogical(test.effect.subgroup.random)
  chknumeric(fontsize, single = TRUE)
  chknumeric(fs.heading, single = TRUE)
  if (!missing(fs.fixed))
    chknumeric(fs.fixed, single = TRUE)
  if (!missing(fs.random))
    chknumeric(fs.random, single = TRUE)
  if (!missing(fs.predict))
    chknumeric(fs.predict, single = TRUE)
  if (!missing(fs.fixed.labels))
    chknumeric(fs.fixed.labels, single = TRUE)
  if (!missing(fs.random.labels))
    chknumeric(fs.random.labels, single = TRUE)
  if (!missing(fs.predict.labels))
    chknumeric(fs.predict.labels, single = TRUE)
  if (!missing(fs.study))
    chknumeric(fs.study, single = TRUE)
  if (!missing(fs.study.labels))
    chknumeric(fs.study.labels, single = TRUE)
  if (!missing(fs.hetstat))
    chknumeric(fs.hetstat, single = TRUE)
  if (!missing(fs.test.overall))
    chknumeric(fs.test.overall, single = TRUE)
  if (!missing(fs.test.subgroup))
    chknumeric(fs.test.subgroup, single = TRUE)
  if (!missing(fs.test.effect.subgroup))
    chknumeric(fs.test.effect.subgroup, single = TRUE)
  chknumeric(fs.axis, single = TRUE)
  chknumeric(fs.smlab, single = TRUE)
  chknumeric(fs.xlab, single = TRUE)
  chknumeric(fs.lr, single = TRUE)
  chknumeric(squaresize, single = TRUE)
  chklogical(calcwidth.pooled)
  chklogical(calcwidth.fixed)
  chklogical(calcwidth.random)
  chklogical(calcwidth.predict)
  chklogical(calcwidth.hetstat)
  chklogical(calcwidth.tests)
  just.cols <- setchar(just, c("right", "center", "left"))
  just.studlab <- setchar(just.studlab, c("right", "center", "left"))
  just.addcols <- setchar(just.addcols, c("right", "center", "left"))
  just.addcols.left <- setchar(just.addcols.left, c("right", "center", "left"))
  just.addcols.right <- setchar(just.addcols.right, c("right", "center", "left"))
  ##
  if (missing(weight.study))
    weight.study <- ifelse(comb.random & !comb.fixed, "random", "fixed")
  weight.study <- setchar(weight.study, c("same", "fixed", "random"))
  ##
  chknumeric(spacing, single = TRUE)
  ##
  ## Check and set additional empty rows in forest plot
  ##
  if (!missing(addrow))
    chklogical(addrow)
  else
    addrow <- !revman5.jama
  if (!missing(addrow.overall))
    chklogical(addrow.overall)
  else
    addrow.overall <- !jama & overall & (comb.fixed | comb.random | prediction)
  if (!missing(addrow.subgroups))
    chklogical(addrow.subgroups)
  else
    addrow.subgroups <- !jama
  ##
  chknumeric(digits, min = 0, single = TRUE)
  chknumeric(digits.tau2, min = 0, single = TRUE)
  chknumeric(digits.zval, min = 0, single = TRUE)
  chknumeric(digits.pval, min = 1, single = TRUE)
  chknumeric(digits.pval.Q, min = 1, single = TRUE)
  chknumeric(digits.Q, min = 0, single = TRUE)
  chknumeric(digits.I2, min = 0, single = TRUE)
  chknumeric(digits.se, min = 0, single = TRUE)
  if (!is.null(digits.mean))
    chknumeric(digits.mean, min = 0, single = TRUE)
  if (!is.null(digits.sd))
    chknumeric(digits.sd, min = 0, single = TRUE)
  if (!is.null(digits.cor))
    chknumeric(digits.cor, min = 0, single = TRUE)
  if (!is.null(digits.time))
    chknumeric(digits.time, min = 0, single = TRUE)
  chklogical(scientific.pval)
  ##
  cl <- paste("update.meta() or ", class(x)[1], "()", sep = "")
  addargs <- names(list(...))
  ##
  fun <- "forest.meta"
  ##
  warnarg("byvar", addargs, fun, cl)
  warnarg("level", addargs, fun, cl)
  warnarg("level.comb", addargs, fun, cl)
  warnarg("level.predict", addargs, fun, cl)
  ##
  ## Check for deprecated argument 'col.i'
  ##
  if (!missing(col.i))
    if (!missing(col.study))
      warning("Deprecated argument 'col.i' ignored as argument 'col.study' is also provided.")
    else {
      warning("Deprecated argument 'col.i' has been replaced by argument 'col.study'.")
      col.study <- col.i
    }
  ##
  ## Check for deprecated argument 'weight'
  ##
  if (!missing(weight))
    if (!missing(weight.study))
      warning("Deprecated argument 'weight' ignored as argument 'weight.study' is also provided.")
    else {
      warning("Deprecated argument 'weight' has been replaced by argument 'weight.study'.")
      weight.study <- weight
    }
  ##
  ## Check for other deprecated arguments in '...'
  ##
  args  <- list(...)
  ## Check whether first argument is a list. In this case only use
  ## this list as input.
  if (length(args) > 0 && is.list(args[[1]]))
    args <- args[[1]]
  ##
  additional.arguments <- names(args)
  ##
  if (length(additional.arguments) > 0) {
    if ("labels" %in% additional.arguments)
      if (!missing(label))
        warning("Argument 'labels' ignored as both arguments 'label' and 'labels' are provided.")
      else
        label <- args[["labels"]]
    ##
    if (!is.na(charmatch("col.i.i", additional.arguments)))
      if (!missing(col.inside))
        warning("Deprecated argument 'col.i.inside.square' ignored as argument 'col.inside' is also provided.")
      else {
        warning("Deprecated argument 'col.i.inside.square' has been replaced by argument 'col.inside'.")
        col.inside <- args[[charmatch("col.i.i", additional.arguments)]]
      }
    ##
    if (!is.na(charmatch("col.diamond.f", additional.arguments)))
      if (!missing(col.diamond.lines.fixed))
        warning("Deprecated argument 'col.diamond.fixed.lines' ignored as argument 'col.diamond.lines.fixed' is also provided.")
      else {
        warning("Deprecated argument 'col.diamond.fixed.lines' has been replaced by argument 'col.diamond.lines.fixed'.")
        col.diamond.lines.fixed <- args[[charmatch("col.diamond.f", additional.arguments)]]
      }
    ##
    if (!is.na(charmatch("col.diamond.r", additional.arguments)))
      if (!missing(col.diamond.lines.random))
        warning("Deprecated argument 'col.diamond.random.lines' ignored as argument 'col.diamond.lines.random' is also provided.")
      else {
        warning("Deprecated argument 'col.diamond.random.lines' has been replaced by argument 'col.diamond.lines.random'.")
        col.diamond.lines.random <- args[[charmatch("col.diamond.r", additional.arguments)]]
      }
    ##
    if (!is.na(charmatch("adds", additional.arguments)))
      if (!missing(addrow))
        warning("Deprecated argument 'addspace' ignored as argument 'addrow' is also provided.")
      else {
        warning("Deprecated argument 'addspace' has been replaced by argument 'addrow'.")
        addrow <- args[[charmatch("adds", additional.arguments)]]
      }
    ##
    if (!is.na(charmatch("text.I2", additional.arguments)))
      warning("Argument 'text.I2' has been removed.")
    if (!is.na(charmatch("text.tau2", additional.arguments)))
      warning("Argument 'text.tau2' has been removed.")
  }
  ##
  ## Additional assignments
  ##
  if (jama) {
    if (missing(ff.fixed))
      ff.fixed <- "plain"
    if (missing(ff.random))
      ff.random <- ff.fixed
    if (missing(ff.predict))
      ff.predict <- ff.fixed
    if (missing(ff.fixed.labels))
      ff.fixed.labels <- ff.fixed
    if (missing(ff.random.labels))
      ff.random.labels <- ff.random
    if (missing(ff.predict.labels))
      ff.predict.labels <- ff.predict
    ##
    if (missing(fs.fixed))
      fs.fixed <- fontsize
    if (missing(fs.random))
      fs.random <- fs.fixed
    if (missing(fs.predict))
      fs.predict <- fs.fixed
    if (missing(fs.fixed.labels))
      fs.fixed.labels <- fs.fixed
    if (missing(fs.random.labels))
      fs.random.labels <- fs.random
    if (missing(fs.predict.labels))
      fs.predict.labels <- fs.predict
  }
  else {
    if (missing(ff.fixed))
      ff.fixed <- "bold"
    if (missing(ff.random))
      ff.random <- ff.fixed
    if (missing(ff.predict))
      ff.predict <- ff.fixed
    if (missing(ff.fixed.labels))
      ff.fixed.labels <- ff.fixed
    if (missing(ff.random.labels))
      ff.random.labels <- ff.random
    if (missing(ff.predict.labels))
      ff.predict.labels <- ff.predict
    ##
    if (missing(fs.fixed))
      fs.fixed <- fontsize
    if (missing(fs.random))
      fs.random <- fs.fixed
    if (missing(fs.predict))
      fs.predict <- fs.fixed
    if (missing(fs.fixed.labels))
      fs.fixed.labels <- fs.fixed
    if (missing(fs.random.labels))
      fs.random.labels <- fs.random
    if (missing(fs.predict.labels))
      fs.predict.labels <- fs.predict
  }
  hetseparator <- " = "
  ##
  if (revman5) {
    if (missing(ff.hetstat))
      ff.hetstat <- "plain"
    if (missing(ff.test.overall))
      ff.test.overall <- ff.hetstat
    if (missing(ff.test.subgroup))
      ff.test.subgroup <- ff.hetstat
    if (missing(ff.test.effect.subgroup))
      ff.test.effect.subgroup <- ff.hetstat
    ##
    if (missing(fs.hetstat))
      fs.hetstat <- fontsize - 1
    if (missing(fs.test.overall))
      fs.test.overall <- fs.hetstat
    if (missing(fs.test.subgroup))
      fs.test.subgroup <- fs.hetstat
    if (missing(fs.test.effect.subgroup))
      fs.test.effect.subgroup <- fs.hetstat
  }
  else if (jama) {
    if (missing(ff.hetstat))
      ff.hetstat <- "plain"
    if (missing(ff.test.overall))
      ff.test.overall <- ff.hetstat
    if (missing(ff.test.subgroup))
      ff.test.subgroup <- ff.hetstat
    if (missing(ff.test.effect.subgroup))
      ff.test.effect.subgroup <- ff.hetstat
    ##
    if (missing(fs.hetstat))
      fs.hetstat <- fontsize - 1
    if (missing(fs.test.overall))
      fs.test.overall <- fs.hetstat
    if (missing(fs.test.subgroup))
      fs.test.subgroup <- fs.hetstat
    if (missing(fs.test.effect.subgroup))
      fs.test.effect.subgroup <- fs.hetstat
  }
  else {
    if (missing(ff.hetstat))
      ff.hetstat <- "plain"
    if (missing(ff.test.overall))
      ff.test.overall <- ff.hetstat
    if (missing(ff.test.subgroup))
      ff.test.subgroup <- ff.hetstat
    if (missing(ff.test.effect.subgroup))
      ff.test.effect.subgroup <- ff.hetstat
    ##
    if (missing(fs.hetstat))
      fs.hetstat <- fontsize - 1
    if (missing(fs.test.overall))
      fs.test.overall <- fs.hetstat
    if (missing(fs.test.subgroup))
      fs.test.subgroup <- fs.hetstat
    if (missing(fs.test.effect.subgroup))
      fs.test.effect.subgroup <- fs.hetstat
  }
  
  
  ##
  ##
  ## (3) Check length of variables
  ##
  ##
  if (length(col.study) == 1)
    col.study <- rep(col.study, K.all)
  else
    chklength(col.study, K.all, fun)
  ##
  if (length(col.inside) == 1)
    col.inside <- rep(col.inside, K.all)
  else
    chklength(col.inside, K.all, fun)
  ##
  miss.col.square <- missing(col.square)
  miss.col.square.lines <- missing(col.square.lines)
  if (length(col.square) == 1)
    col.square <- rep(col.square, K.all)
  else
    chklength(col.square, K.all, fun)
  ##
  if (length(col.square.lines) == 1)
    col.square.lines <- rep(col.square.lines, K.all)
  else
    chklength(col.square.lines, K.all, fun)
  
  
  ##
  ##
  ## (4) Some assignments and additional checks
  ##
  ##
  metabin <- inherits(x, "metabin")
  metacont <- inherits(x, "metacont")
  metacor <- inherits(x, "metacor")
  metagen <- inherits(x, "metagen")
  metainc <- inherits(x, "metainc")
  metainf.metacum <- inherits(x, "metainf") | inherits(x, "metacum")
  metaprop <- inherits(x, "metaprop")
  metarate <- inherits(x, "metarate")
  ##
  fixed.random <- comb.fixed & comb.random
  ##
  if (layout == "subgroup") {
    if (!missing(study.results) & study.results)
      warning("Argument 'study.results' set to FALSE as argument 'layout' is \"subgroup\".")
    study.results <- FALSE
  }
  ##
  if (is.null(text.fixed)) {
    if (study.results & (x$level != x$level.comb | revman5)) {
      if (revman5.jama)
        text.fixed <- paste("Total (",
                            if (fixed.random)
                              "fixed effect, ",
                            round(x$level.comb * 100), "% CI)", sep = "")
      else
        text.fixed <- paste("Fixed effect model (",
                            round(x$level.comb * 100), "%-CI)", sep = "")
    }
    else {
      if (revman5.jama) {
        text.fixed <- "Total"
        if (fixed.random)
          text.fixed <- paste(text.fixed, "(fixed effect)")
      }
      else
        text.fixed <- "Fixed effect model"
    }
  }
  if (is.null(text.random)) {
    if (study.results & (x$level != x$level.comb | revman5)) {
      if (revman5.jama)
        text.random <- paste("Total (",
                            if (fixed.random)
                              "random effects, ",
                            round(x$level.comb * 100), "% CI)", sep = "")
      else
        text.random <- paste("Random effects model (",
                            round(x$level.comb * 100), "%-CI)", sep = "")
    }
    else {
      if (revman5.jama) {
        text.random <- "Total"
        if (fixed.random)
          text.random <- paste(text.random, "(random effects)")
      }
      else
        text.random <- "Random effects model"
    }
  }
  if (is.null(text.predict))
    if (!(length(x$level.predict) == 0) &&
        (study.results & (x$level != x$level.predict | x$level.comb != x$level.predict)))
      text.predict <- paste("Prediction interval (",
                            round(x$level.predict * 100), "%-PI)", sep = "")
    else
      text.predict <- "Prediction interval"
  ##
  if (metainf.metacum) {
    overall.hetstat <- FALSE
    hetstat <- FALSE
    prediction <- FALSE
    test.overall.fixed   <- FALSE
    test.overall.random  <- FALSE
    test.subgroup.fixed  <- FALSE
    test.subgroup.random <- FALSE
    test.effect.subgroup.fixed  <- FALSE
    test.effect.subgroup.random <- FALSE
  }
  ##
  if (metaprop | metarate) {
    test.overall.fixed  <- FALSE
    test.overall.random <- FALSE
    test.effect.subgroup.fixed  <- FALSE
    test.effect.subgroup.random <- FALSE
  }
  ##
  if (!overall) {
    if (test.overall)
      test.overall <- FALSE
    if (test.overall.fixed)
      test.overall.fixed <- FALSE
    if (test.overall.random)
      test.overall.random <- FALSE
  }
  ##
  prediction <- prediction & x$k >= 3
  ##
  byvar <- x$byvar
  level <- x$level
  level.comb <- x$level.comb
  level.predict <- x$level.predict
  ##
  chklevel(level)
  chklevel(level.comb)
  if (prediction)
    chklevel(level.predict)
  ##
  if (is.logical(leftcols)) {
    if (!leftcols) {
      text.fixed <- ""
      text.random <- ""
      text.predict <- ""
      leftcols <- "studlab"
      studlab <- rep("", K.all)
      slab <- FALSE
      hetstat <- FALSE
      overall.hetstat <- FALSE
      colgap.left <- grid::unit(0, "mm")
      colgap.studlab <- grid::unit(0, "mm")
      colgap.forest.left <- grid::unit(0, "mm")
    }
    else
      leftcols <- NULL
  }
  ##
  notmiss.xlim <- !missing(xlim)
  ##
  if (just.studlab == "left")
    xpos.s <- 0
  else if (just.studlab == "center")
    xpos.s <- 0.5
  else if (just.studlab == "right")
    xpos.s <- 1
  ##
  if (just.cols == "left")
    xpos.c <- 0
  else if (just.cols == "center")
    xpos.c <- 0.5
  else if (just.cols == "right")
    xpos.c <- 1
  ##
  sm <- x$sm
  ##
  log.xaxis <- FALSE
  ##
  if (missing(ref) && sm %in% c("PLOGIT", "PLN", "PRAW", "PAS", "PFT",
                                "IR", "IRLN", "IRS", "IRFT"))
    ref <- NA
  ##
  if (backtransf & is.relative.effect(sm)) {
    ref <- log(ref)
    log.xaxis <- TRUE
  }
  ##
  if (!backtransf & !missing(pscale) & pscale != 1) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  ##
  if (!backtransf & !missing(irscale) & irscale != 1) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  ##
  if (is.null(xlab))
    xlab <- xlab(sm, backtransf, newline = revman5.jama, revman5 = revman5)
  ##
  smlab.null <- is.null(smlab)
  if (smlab.null)
    if (sm %in% c("IR", "IRLN", "IRS", "IRFT"))
      smlab <- xlab(sm, backtransf, irscale = irscale, irunit = irunit,
                    newline = !revman5.jama, revman5 = revman5)
    else
      smlab <- xlab(sm, backtransf, pscale = pscale,
                    newline = !revman5.jama, revman5 = revman5)
  ##
  if (is.null(label.right))
    label.right <- ""
  if (is.null(label.left))
    label.left <- ""
  ##
  print.label <- label.left != "" | label.right != "" & !is.na(ref)
  if (print.label & !bottom.lr) {
    if (!smlab.null)
      warning("Argument 'smlab' ignored as argument 'bottom.lr' is FALSE.")
    smlab <- ""
  }
  ##
  by <- !is.null(byvar)
  if (!by)
    addrow.subgroups <- FALSE
  ##
  plotwidth <- setunit(plotwidth)
  colgap <- setunit(colgap)
  colgap.left <- setunit(colgap.left)
  colgap.right <- setunit(colgap.right)
  colgap.studlab <- setunit(colgap.studlab)
  colgap.forest <- setunit(colgap.forest)
  colgap.forest.left <- setunit(colgap.forest.left)
  colgap.forest.right <- setunit(colgap.forest.right)
  ##
  if (missing(label.test.overall.fixed))
    label.test.overall.fixed <-
      paste("Test for overall effect",
            if (fixed.random) " (fixed effect)",
            ": ", sep = "")
  if (missing(label.test.overall.random))
    label.test.overall.random <-
      paste("Test for overall effect",
            if (fixed.random) " (random effects)",
            ": ", sep = "")
  ##
  if (missing(label.test.subgroup.fixed))
    label.test.subgroup.fixed <-
      paste("Test for subgroup differences",
            if (fixed.random) " (fixed effect)",
            ": ", sep = "")
  if (missing(label.test.subgroup.random))
    label.test.subgroup.random <-
      paste("Test for subgroup differences",
            if (fixed.random) " (random effects)",
            ": ", sep = "")
  ##
  if (missing(label.test.effect.subgroup.fixed))
    label.test.effect.subgroup.fixed <-
      paste(if (revman5.jama) "Test for overall effect" else "Test for effect in subgroup",
            if (fixed.random) " (fixed effect)",
            ": ", sep = "")
  if (missing(label.test.effect.subgroup.random))
    label.test.effect.subgroup.random <-
      paste(if (revman5.jama) "Test for overall effect" else "Test for effect in subgroup",
            if (fixed.random) " (random effects)",
            ": ", sep = "")
  ##
  fs.head <- fs.heading
  ff.head <- ff.heading
  ##
  just.c <- just.cols
  just.s <- just.studlab
  
  
  ##
  ##
  ## (5) Determine columns on left and right side of forest plot
  ##
  ##
  ## Determine whether to print columns on right and / or left side
  ## of forest plot
  ##
  rsel <- !(is.logical(rightcols) && length(rightcols) == 1 && !rightcols)
  if (revman5.jama)
    rsel <- FALSE
  ##
  if (!rsel)
    rightcols <- NULL
  ##
  ## Check for duplicate columns
  ##
  if (length(c(rightcols, leftcols)) > 0 &&
      any(duplicated(c(rightcols, leftcols))))
    stop("Duplicate entries in 'leftcols' and 'rightcols'.")
  ##
  ## Predefined columns and labels
  ##
  sm.lab <- sm
  ##
  if (backtransf) {
    if (sm == "ZCOR")
      sm.lab <- "COR"
    else if (sm %in% c("PLOGIT", "PLN", "PRAW", "PAS", "PFT")) {
      if (pscale == 1)
        sm.lab <- "Proportion"
      else
        sm.lab <- "Events"
    }
    else if (sm %in% c("IR", "IRLN", "IRS", "IRFT")) {
      if (irscale == 1)
        sm.lab <- "Rate"
      else
        sm.lab <- "Events"
    }
    else if (sm == "proportion")
      sm.lab <- "Proportion"
  }
  else 
    if (is.relative.effect(sm))
      sm.lab <- paste("log", sm, sep = "")
  ##
  colnames <- c("studlab", "TE", "seTE",
                "n.e", "n.c",
                "event.e", "event.c",
                "mean.e", "mean.c",
                "sd.e", "sd.c",
                "cor",
                "time.e", "time.c",
                "effect", "ci",
                "effect.ci",
                "w.fixed", "w.random")
  ##
  sel.studlab <- pmatch(layout, c("meta", "RevMan5", "JAMA", "subgroup"))
  lab.studlab <- c("Study", "Study", "Source", "Subgroup")[sel.studlab]
  if (revman5 & by)
    lab.studlab <- c("Study or\nSubgroup")
  ##
  if (revman5.jama)
    cisep <- " "
  else
    cisep <- "-"
  ##
  if (study.results)
    ci.lab <- paste(100 * level, "%", cisep, "CI", sep = "")
  else
    ci.lab <- paste(100 * level.comb, "%", cisep, "CI", sep = "")
  ##
  if (jama) {
    if (missing(ff.lr))
      ff.lr <- "bold"
    if (xlab == "")
      xlab <- paste(smlab, " (", ci.lab, ")", sep = "")
    if (miss.col.square)
      col.square <- rep("darkblue", K.all)
    if (miss.col.square.lines)
      col.square.lines <- rep("darkblue", K.all)
    if (missing(col.diamond.fixed))
      col.diamond.fixed <- "lightblue"
    if (missing(col.diamond.random))
      col.diamond.random <- "lightblue"
    ##
    smlab <- ""
    bottom.lr <- FALSE
  }
  else {
    if (revman5) {
      if (miss.col.square) {
        if (metacont)
          col.square <- rep("green", K.all)
        else if (metabin)
          col.square <- rep("blue", K.all)
        else
          col.square <- rep("red", K.all)
      }
      if (miss.col.square.lines) {
        if (metacont)
          col.square.lines <- rep("green", K.all)
        else if (metabin)
          col.square.lines <- rep("darkblue", K.all)
        else
          col.square.lines <- rep("red", K.all)
      }
      if (missing(col.diamond.fixed))
      col.diamond.fixed <- "black"
      if (missing(col.diamond.random))
        col.diamond.random <- "black"
      ##
      sel.method <- pmatch(x$method, c("Inverse", "MH", "Peto", "GLMM"))
      lab.method <- c("IV", "MH", "Peto", "GLMM")[sel.method]
      ##
      if (fixed.random)
        lab.model <- "Fixed + Random, "
      else if (comb.fixed)
        lab.model <- "Fixed, "
      else if (comb.random)
        lab.model <- "Random, "
      else
        lab.model <- ""
      ##
      if (smlab.null)
        smlab <- paste(smlab, "\n", lab.method, ", ",
                       lab.model,
                       ci.lab, sep = "")
    }
  }
  ##
  if (jama | gs("CIbracket") == "(")
    ci.lab.bracket <- paste("(", ci.lab, ")", sep = "")
  else if (gs("CIbracket") == "[")
    ci.lab.bracket <- paste("[", ci.lab, "]", sep = "")
  else if (gs("CIbracket") == "{")
    ci.lab.bracket <- paste("{", ci.lab, "}", sep = "")
  else if (gs("CIbracket") == "")
    ci.lab.bracket <- ci.lab
  ##
  labnames <- c(lab.studlab,
                "TE", if (revman5) "SE" else "seTE",
                "Total", "Total", "Events", "Events",
                "Mean", "Mean", "SD", "SD",
                "Cor",
                "Time", "Time",
                sm.lab,
                ci.lab,
                if (revman5 & smlab.null) smlab else paste(sm.lab, ci.lab.bracket),
                if (fixed.random) "Weight\n(fixed)" else "Weight",
                if (fixed.random) "Weight\n(random)" else "Weight")
  ##
  ## If any of the following list elements is NULL, these 'special'
  ## variable names are searched for in original data set (i.e., list
  ## element x$data)
  ##
  colnames.notNULL <- colnames
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "n.e")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "n.c")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "event.e")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "event.c")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "mean.e")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "mean.c")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "sd.e")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "sd.c")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "cor")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "time.e")
  colnames.notNULL <- removeNULL(x, colnames.notNULL, "time.c")
  ##
  ## Identify and process columns in addition to columns
  ## defined above in variables 'colnames' and 'labnames'
  ##
  colnames.new <- c(rightcols, leftcols)[!c(rightcols, leftcols) %in% colnames.notNULL]
  ##
  newcols <- length(colnames.new) > 0
  ##
  if (newcols) {
    dataset2 <- as.data.frame(x)
    ##
    if (is.null(x$data))
      dataset1 <- dataset2
    else
      dataset1 <- x$data
    ##
    if (!is.null(x$subset))
      dataset1 <- dataset1[x$subset, ]
    ##
    ## Check whether additional variables are
    ## part of meta-object
    ##
    for (i in colnames.new)
      if (length(dataset1[[i]]) == 0 & length(dataset2[[i]]) == 0)
        stop("Variable '", i, "' not available in '", x.name, "'.")
    ##
    rightcols.new <- rightcols[! rightcols %in% colnames.notNULL]
    leftcols.new  <- leftcols[! leftcols %in% colnames.notNULL]
    ##
    ## Determine label for new columns
    ## 1. Use column name as label if no label is given
    ##    argument right | left | labs
    ## 2. Otherwise use corresponding entry from
    ##    argument right | left | labs
    ##
    if (length(rightcols.new) > 0) {
      pos.rightcols.new <- match(rightcols.new, rightcols)
      ##
      rightlabs.new <- rightcols.new
      for (i in seq(along = rightcols.new)) {
        j <- match(rightcols.new[i], colnames)
        if (!is.na(j))
          rightlabs.new[i] <- labnames[j]
      }
      ##
      if (missing(rightlabs))
        rightlabs.new <- rightlabs.new
      else if (length(rightcols.new) == length(rightlabs))
        rightlabs.new <- rightlabs
      else if (max(pos.rightcols.new) <= length(rightlabs))
        rightlabs.new <- rightlabs[pos.rightcols.new]
      else if (max(pos.rightcols.new) > length(rightlabs))
        stop("Too few labels defined for argument 'rightcols'.")
      ##
      if ( (metacor | metaprop) & any(rightcols.new == "n"))
        rightlabs.new[rightlabs.new == "n"] <- "Total"
      ##
      if ( (metarate) & any(rightcols.new == "time"))
        rightlabs.new[rightlabs.new == "time"] <- "Time"
    }
    if (length(leftcols.new) > 0) {
      pos.leftcols.new <- match(leftcols.new, leftcols)
      ##
      leftlabs.new <- leftcols.new
      for (i in seq(along = leftcols.new)) {
        j <- match(leftcols.new[i], colnames)
        if (!is.na(j))
          leftlabs.new[i] <- labnames[j]
      }
      ##
      if (missing(leftlabs))
        leftlabs.new <- leftlabs.new
      else if (length(leftcols.new) == length(leftlabs))
        leftlabs.new <- leftlabs
      else if (max(pos.leftcols.new) <= length(leftlabs))
        leftlabs.new <- leftlabs[pos.leftcols.new]
      else if (max(pos.leftcols.new) > length(leftlabs))
        stop("Too few labels defined for argument 'leftcols'.")
      ##
      if ((metacor | metaprop) & any(leftcols.new == "n"))
        leftlabs.new[leftlabs.new == "n"] <- "Total"
      ##
      if ((metarate) & any(leftcols.new == "time"))
        leftlabs.new[leftlabs.new == "time"] <- "Time"
    }
  }
  ##
  ## Default set of columns if argument leftcols and / or
  ## rightcols not specified
  ##
  if (is.null(leftcols)) {
    ##
    leftcols <- "studlab"
    ##
    if (jama) {
      leftcols <- c(leftcols,
                    "effect.ci")
    }
    else {
      if (metabin) {
        if (study.results)
          leftcols <- c(leftcols,
                        "event.e", "n.e",
                        "event.c", "n.c")
        else {
          leftcols <- c(leftcols,
                        if (pooled.events) "event.e",
                        if (pooled.totals) "n.e",
                        if (pooled.events) "event.c",
                        if (pooled.totals) "n.c")
          if (pooled.events & !pooled.totals) {
            if (is.null(lab.e.attach.to.col))
              lab.e.attach.to.col <- "event.e"
            if (is.null(lab.c.attach.to.col))
              lab.c.attach.to.col <- "event.c"
          }
        }
      }
      ##
      if (metacont) {
        if (study.results) {
          if (revman5)
            leftcols <- c(leftcols,
                          "mean.e", "sd.e", "n.e",
                          "mean.c", "sd.c", "n.c")
          else
            leftcols <- c(leftcols,
                          "n.e", "mean.e", "sd.e",
                          "n.c", "mean.c", "sd.c")
        }
        else if (pooled.totals) {
          leftcols <- c(leftcols, "n.e", "n.c")
          if (is.null(lab.e.attach.to.col))
            lab.e.attach.to.col <- "n.e"
          if (is.null(lab.c.attach.to.col))
            lab.c.attach.to.col <- "n.c"
        }
      }
      ##
      if (metagen & study.results)
        leftcols <- c(leftcols,
                      "TE", "seTE")
      ##
      if (metaprop) {
        if (study.results)
          leftcols <- c(leftcols,
                        "event.e", "n.e")
        else {
          leftcols <- c(leftcols,
                      if (pooled.events) "event.e",
                      if (pooled.totals) "n.e")
          if (pooled.events & !pooled.totals) {
            if (is.null(lab.e.attach.to.col))
              lab.e.attach.to.col <- "event.e"
          }
        }
      }
      ##
      if (metarate) {
        if (study.results)
          leftcols <- c(leftcols,
                        "event.e", "time.e")
        else {
          leftcols <- c(leftcols,
                        if (pooled.events) "event.e",
                        if (pooled.times) "time.e")
          if (pooled.events & !pooled.times) {
            if (is.null(lab.e.attach.to.col))
              lab.e.attach.to.col <- "event.e"
          }
        }
      }
      ##
      if (metacor) {
        if (study.results | pooled.totals)
          leftcols <- c(leftcols,
                        "n.e")
      }
      ##
      if (metainc) {
        if (study.results)
          leftcols <- c(leftcols,
                        "event.e", "time.e",
                        "event.c", "time.c")
        else {
          leftcols <- c(leftcols,
                        if (pooled.events) "event.e",
                        if (pooled.times) "time.e",
                        if (pooled.events) "event.c",
                        if (pooled.times) "time.c")
          if (pooled.events & !pooled.times) {
            if (is.null(lab.e.attach.to.col))
              lab.e.attach.to.col <- "event.e"
            if (is.null(lab.c.attach.to.col))
              lab.c.attach.to.col <- "event.c"
          }
        }
      }
    }
    ##
    ## Add columns for RevMan 5 layout
    ##
    if (revman5) {
      ##
      if (!metainf.metacum & overall & study.results & !x$method == "GLMM") {
        if (comb.fixed)
          leftcols <- c(leftcols, "w.fixed")
        if (comb.random)
          leftcols <- c(leftcols, "w.random")
      }
      ##
      leftcols <- c(leftcols, "effect.ci")
    }
  }
  ##
  if (is.null(rightcols) & rsel) {
    rightcols <- c("effect", "ci")
    ##
    if (!metainf.metacum & overall & study.results & !x$method == "GLMM") {
      if (comb.fixed)
        rightcols <- c(rightcols, "w.fixed")
      if (comb.random)
        rightcols <- c(rightcols, "w.random")
    }
  }
  
  
  ##
  ##
  ## (6) Select data for forest plot
  ##
  ##
  if (metacor) {
    x$n.e <- x$n
  }
  ##
  if (metaprop) {
    x$event.e <- x$event
    x$n.e <- x$n
    ##
    if (!is.null(rightcols)) {
      if (any(rightcols == "n"))
        rightcols[rightcols == "n"] <- "n.e"
      if (any(rightcols == "event"))
        rightcols[rightcols == "event"] <- "event.e"
    }
    ##
    if (!is.null(leftcols)) {
      if (any(leftcols == "n"))
        leftcols[leftcols == "n"] <- "n.e"
      if (any(leftcols == "event"))
        leftcols[leftcols == "event"] <- "event.e"
    }
  }
  if (metarate) {
    x$event.e <- x$event
    x$time.e <- x$time
    ##
    if (!is.null(rightcols)) {
      if (any(rightcols == "time"))
        rightcols[rightcols == "time"] <- "time.e"
      if (any(rightcols == "event"))
        rightcols[rightcols == "event"] <- "event.e"
    }
    ##
    if (!is.null(leftcols)) {
      if (any(leftcols == "time"))
        leftcols[leftcols == "time"] <- "time.e"
      if (any(leftcols == "event"))
        leftcols[leftcols == "event"] <- "event.e"
    }
  }
  ##
  if (metainf.metacum) {
    ##
    x$TE.fixed    <- rev(x$TE)[1]
    x$seTE.fixed  <- rev(x$seTE)[1]
    x$lower.fixed <- rev(x$lower)[1]
    x$upper.fixed <- rev(x$upper)[1]
    ##
    x$TE.random   <- rev(x$TE)[1]
    x$seTE.random <- rev(x$seTE)[1]
    x$lower.random <- rev(x$lower)[1]
    x$upper.random <- rev(x$upper)[1]
    ##
    x$n.harmonic.mean.ma <- rev(x$n.harmonic.mean)[1]
    x$t.harmonic.mean.ma <- rev(x$t.harmonic.mean)[1]
    ##
    x$w.all <- rev(x$w)[1]
    ##
    x$TE <- rev(rev(x$TE)[-(1:2)])
    x$seTE <- rev(rev(x$seTE)[-(1:2)])
    x$studlab <- rev(rev(x$studlab)[-(1:2)])
    ##
    x$lower <- rev(rev(x$lower)[-(1:2)])
    x$upper <- rev(rev(x$upper)[-(1:2)])
    ##
    x$w.fixed <- rev(rev(x$w)[-(1:2)])
    x$w.random <- rev(rev(x$w)[-(1:2)])
    ##
    x$n.harmonic.mean <- rev(rev(x$n.harmonic.mean)[-(1:2)])
    x$t.harmonic.mean <- rev(rev(x$t.harmonic.mean)[-(1:2)])
    ##
    if (overall & x$pooled == "fixed") {
      comb.fixed <- TRUE
      comb.random <- FALSE
      if (weight.study != "same")
        weight.study <- "fixed"
    }
    else if (overall & x$pooled == "random") {
      comb.fixed <- FALSE
      comb.random <- TRUE
      if (weight.study != "same")
        weight.study <- "random"
    }
  }
  ## Total number of studies to plot (*not* number of studies combined)
  k.all <- length(x$TE)
  ##
  if (allstudies)
    n.stud <- k.all # all studies
  else
    n.stud <- x$k   # number of studies combined in meta-analysis
  ##
  if (!by)
    byvar <- rep(1, k.all)
  ##
  if (by & any(is.na(byvar)))
    stop("Missing values in 'byvar'")
  ##
  if (allstudies)
    sel <- 1:length(x$TE)
  else
    sel <- !is.na(x$TE)
  ##
  if (n.stud != sum(sel > 0))
    warning("n.stud != sum(sel)")
  ##
  x$n.e <- x$n.e[sel]
  x$n.c <- x$n.c[sel]
  ##
  x$event.e <- x$event.e[sel]
  x$event.c <- x$event.c[sel]
  ##
  x$mean.e <- x$mean.e[sel]
  x$mean.c <- x$mean.c[sel]
  ##
  x$sd.e <- x$sd.e[sel]
  x$sd.c <- x$sd.c[sel]
  ##
  x$cor <- x$cor[sel]
  ##
  x$time.e <- x$time.e[sel]
  x$time.c <- x$time.c[sel]
  ##
  x$TE   <- x$TE[sel]
  x$seTE <- x$seTE[sel]
  ##
  x$lower <- x$lower[sel]
  x$upper <- x$upper[sel]
  ##
  x$w.fixed  <- x$w.fixed[sel]
  x$w.random <- x$w.random[sel]
  studlab  <- studlab[sel]
  ##
  x$n.harmonic.mean <- x$n.harmonic.mean[sel]
  x$t.harmonic.mean <- x$t.harmonic.mean[sel]
  ##
  byvar   <- byvar[sel]
  sortvar <- sortvar[sel]
  ##
  col.study <- col.study[sel]
  col.square <- col.square[sel]
  col.square.lines <- col.square.lines[sel]
  ##
  col.inside <- col.inside[sel]
  ##
  if (!is.null(x$exclude))
    exclude <- x$exclude[sel]
  ##
  if (sort | by) {
    if (bysort)
      bylevs <- sort(x$bylevs)
    else
      bylevs <- x$bylevs
    ##
    byvar.factor <- factor(byvar, levels = bylevs)
    o <- order(byvar.factor, sortvar)
    ##
    x$n.e <- x$n.e[o]
    x$n.c <- x$n.c[o]
    ##
    x$event.e <- x$event.e[o]
    x$event.c <- x$event.c[o]
    ##
    x$mean.e <- x$mean.e[o]
    x$mean.c <- x$mean.c[o]
    ##
    x$sd.e <- x$sd.e[o]
    x$sd.c <- x$sd.c[o]
    ##
    x$cor <- x$cor[o]
    ##
    x$time.e <- x$time.e[o]
    x$time.c <- x$time.c[o]
    ##
    x$TE   <- x$TE[o]
    x$seTE <- x$seTE[o]
    ##
    x$lower <- x$lower[o]
    x$upper <- x$upper[o]
    ##
    x$w.fixed  <- x$w.fixed[o]
    x$w.random <- x$w.random[o]
    studlab  <- studlab[o]
    ##
    x$n.harmonic.mean <- x$n.harmonic.mean[o]
    x$t.harmonic.mean <- x$t.harmonic.mean[o]
    ##
    byvar   <- byvar[o]
    sortvar <- sortvar[o]
    ##
    col.study <- col.study[o]
    col.square <- col.square[o]
    col.square.lines <- col.square.lines[o]
    ##
    col.inside <- col.inside[o]
    ##
    if (!is.null(x$exclude))
      exclude <- exclude[o]
    ##
    if (newcols) {
      dataset1 <- dataset1[o, ]
      dataset2 <- dataset2[o, ]
    }
  }
  ##
  if (by) {
    n.by <- length(bylevs)
    sel.by.fixed         <- 3 + 0 * n.by + 1:n.by
    sel.by.random        <- 3 + 1 * n.by + 1:n.by
    sel.by.het           <- 3 + 2 * n.by + 1:n.by
    sel.by.effect.fixed  <- 3 + 3 * n.by + 1:n.by
    sel.by.effect.random <- 3 + 4 * n.by + 1:n.by
    sel.by <- c(sel.by.fixed, sel.by.random, sel.by.het,
                sel.by.effect.fixed, sel.by.effect.random)
    sel.by.noNA <- c(sel.by.het, sel.by.effect.fixed, sel.by.effect.random)
  }
  else
    n.by <- 0
  ##
  if (metainf.metacum) {
    TE    <- x$TE
    seTE  <- x$seTE
    lowTE <- x$lower
    uppTE <- x$upper
    ##
    TE.fixed    <- x$TE.fixed
    lowTE.fixed <- x$lower.fixed
    uppTE.fixed <- x$upper.fixed
    ##
    TE.random    <- x$TE.random
    lowTE.random <- x$lower.random
    uppTE.random <- x$upper.random
    ##
    lowTE.predict <- NA
    uppTE.predict <- NA
    ##
    Q    <- NA
    df   <- NA
    I2   <- NA
    tau2 <- NA
    lowI2 <- NA
    uppI2 <- NA
    Rb   <- NA
    lowRb <- NA
    uppRb <- NA
    ##
    Q.b.fixed  <- NA
    Q.b.random <- NA
    df.Q.b     <- NA
  }
  else {
    TE <- x$TE
    seTE <- x$seTE
    lowTE <- x$lower
    uppTE <- x$upper
    ##
    if (metaprop & !backtransf) {
      ciTE <- ci(TE, seTE, level = level)
      lowTE <- ciTE$lower
      uppTE <- ciTE$upper
    }
    ##
    TE.fixed <- x$TE.fixed
    lowTE.fixed <- x$lower.fixed
    uppTE.fixed <- x$upper.fixed
    ##
    TE.random <- x$TE.random
    lowTE.random <- x$lower.random
    uppTE.random <- x$upper.random
    ##
    lowTE.predict <- x$lower.predict
    uppTE.predict <- x$upper.predict
    ##
    Q    <- x$Q
    df   <- x$df.Q
    tau2 <- x$tau^2
    ##
    I2 <- x$I2
    lowI2 <- x$lower.I2
    uppI2 <- x$upper.I2
    ##
    Rb <- x$Rb
    lowRb <- x$lower.Rb
    uppRb <- x$upper.Rb
    ##
    Q.b.fixed  <- x$Q.b.fixed
    Q.b.random <- x$Q.b.random
    df.Q.b     <- x$df.Q.b
  }
  ##
  hetstat.overall <- ""
  ##
  if (overall.hetstat) {
    ##
    hetstat.I2 <-
      paste(hetseparator,
            format.NA(round(100 * I2, digits.I2),
                      digits.I2, "NA"), "%",
            if (print.I2.ci & x$k > 2)
              paste(" ",
                    p.ci(paste(format.NA(round(100 * lowI2, digits.I2),
                                         digits.I2, lab.NA),
                               "%", sep = ""),
                         paste(format.NA(round(100 * uppI2, digits.I2),
                                         digits.I2, lab.NA),
                               "%", sep = "")),
                    sep = ""),
            sep = "")
    ##
    hetstat.tau2 <-
      paste(hetseparator,
            if (is.na(tau2)) "NA"
            else if (tau2 == 0) "0"
            else format.tau(tau2, digits = digits.tau2),
            sep = "")
    ##
    hetstat.Q <-
      paste(hetseparator,
            format.NA(round(Q, digits.Q), digits.Q, "NA"),
            if (revman5) ", df",
            if (revman5) hetseparator,
            if (revman5) df,
            sep = "")
    ##
    hetstat.pval.Q <-
      paste(format.p(pvalQ(Q, df),
                     lab = TRUE, labval = "",
                     digits = digits.pval.Q,
                     zero = if (jama) FALSE else TRUE,
                     scientific = scientific.pval,
                     lab.NA = "NA"),
            sep = "")
    ##
    hetstat.Rb <-
      paste(hetseparator,
            format.NA(round(100 * Rb, digits.I2),
                      digits.I2, "NA"),
                      "%",
            if (print.Rb.ci & x$k > 2)
              paste(" ",
                    p.ci(paste(round(100 * lowRb, digits.I2), "%", sep = ""),
                         paste(round(100 * uppRb, digits.I2), "%", sep = "")),
                    sep = ""),
            sep = "")
    ##
    ## Remove superfluous spaces
    ##
    while(grepl("  ", hetstat.I2))
      hetstat.I2 <- gsub("  ", " ", hetstat.I2)
    while(grepl("  ", hetstat.tau2))
      hetstat.tau2 <- gsub("  ", " ", hetstat.tau2)
    while(grepl("  ", hetstat.Q))
      hetstat.Q <- gsub("  ", " ", hetstat.Q)
    while(grepl("  ", hetstat.pval.Q))
      hetstat.pval.Q <- gsub("  ", " ", hetstat.pval.Q)
    while(grepl("  ", hetstat.Rb))
      hetstat.Rb <- gsub("  ", " ", hetstat.Rb)
    ##
    if (revman5)
      hetstat.overall <- substitute(paste(hl,
                                          "Tau"^2, ht, "; ",
                                          "Chi"^2, hq,
                                          " (",
                                          P, hp,
                                          "); ",
                                          I^2, hi),
                                    list(hl = hetlab,
                                         hi = hetstat.I2,
                                         ht = hetstat.tau2,
                                         hq = hetstat.Q,
                                         hp = hetstat.pval.Q))
    else if (jama)
      hetstat.overall <- substitute(paste(hl,
                                          chi[df]^2, hq,
                                          " (",
                                          italic(P), hp,
                                          "), ",
                                          italic(I)^2, hi
                                          ),
                                    list(hl = hetlab,
                                         hi = hetstat.I2,
                                         hq = hetstat.Q,
                                         hp = hetstat.pval.Q,
                                         df = df))
    else {
      ##
      ## One
      ##
      if (print.I2 & !print.tau2 & !print.Q & !print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl, italic(I)^2, hi),
                     list(hl = hetlab, hi = hetstat.I2))
      else if (!print.I2 & print.tau2 & !print.Q & !print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl, tau^2, ht),
                     list(hl = hetlab, ht = hetstat.tau2))
      else if (!print.I2 & !print.tau2 & print.Q & !print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl, chi[df]^2, hq),
                     list(hl = hetlab, df = df, hq = hetstat.Q))
      else if (!print.I2 & !print.tau2 & !print.Q & print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl, italic(p), hp),
                     list(hl = hetlab, hp = hetstat.pval.Q))
      else if (!print.I2 & !print.tau2 & !print.Q & !print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl, italic(R)[italic(b)], hb),
                     list(hl = hetlab, hb = hetstat.Rb))
      ##
      ## Two
      ##
      else if (print.I2 & print.tau2 & !print.Q & !print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl, italic(I)^2, hi,
                           ", ",
                           tau^2, ht),
                     list(hl = hetlab,
                          hi = hetstat.I2, ht = hetstat.tau2))
      else if (print.I2 & !print.tau2 & print.Q & !print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl, italic(I)^2, hi,
                           ", ",
                           chi[df]^2, hq),
                     list(hl = hetlab, df = df,
                          hi = hetstat.I2, hq = hetstat.Q))
      else if (print.I2 & !print.tau2 & !print.Q & print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl, italic(I)^2, hi,
                           ", ",
                           italic(p), hp),
                     list(hl = hetlab,
                          hi = hetstat.I2, hp = hetstat.pval.Q))
      else if (print.I2 & !print.tau2 & !print.Q & !print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl, italic(I)^2, hi,
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab,
                          hi = hetstat.I2, hb = hetstat.Rb))
      else if (!print.I2 & print.tau2 & print.Q & !print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl, tau^2, ht,
                           ", ",
                           chi[df]^2, hq),
                     list(hl = hetlab, df = df,
                          ht = hetstat.tau2, hq = hetstat.Q))
      else if (!print.I2 & print.tau2 & !print.Q & print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl, tau^2, ht,
                           ", ",
                           italic(p), hp),
                     list(hl = hetlab,
                          ht = hetstat.tau2, hp = hetstat.pval.Q))
      else if (!print.I2 & print.tau2 & !print.Q & !print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl, tau^2, ht,
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab,
                          ht = hetstat.tau2, hb = hetstat.Rb))
      else if (!print.I2 & !print.tau2 & print.Q & print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl, chi[df]^2, hq,
                           " (",
                           italic(p), hp, ")"),
                     list(hl = hetlab, df = df,
                          hq = hetstat.Q, hp = hetstat.pval.Q))
      else if (!print.I2 & !print.tau2 & print.Q & !print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl, chi[df]^2, hq,
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab, df = df,
                          hq = hetstat.Q, hb = hetstat.Rb))
      else if (!print.I2 & !print.tau2 & !print.Q & print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl, italic(p), hp,
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab,
                          hp = hetstat.pval.Q, hb = hetstat.Rb))
      ##
      ## Three
      ##
      else if (print.I2 & print.tau2 & print.Q & !print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           tau^2, ht, ", ",
                           chi[df]^2, hq),
                     list(hl = hetlab, df = df,
                          hi = hetstat.I2, ht = hetstat.tau2,
                          hq = hetstat.Q))
      else if (print.I2 & print.tau2 & !print.Q & print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           tau^2, ht, ", ",
                           italic(p), hp),
                     list(hl = hetlab,
                          hi = hetstat.I2, ht = hetstat.tau2,
                          hp = hetstat.pval.Q))
      else if (print.I2 & !print.tau2 & print.Q & print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           chi[df]^2, hq,
                           " (", italic(p), hp, ")"),
                     list(hl = hetlab, df = df,
                          hi = hetstat.I2,
                          hq = hetstat.Q, hp = hetstat.pval.Q))
      else if (!print.I2 & print.tau2 & print.Q & print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           tau^2, ht, ", ",
                           chi[df]^2, hq,
                           " (", italic(p), hp, ")"),
                     list(hl = hetlab, df = df,
                          ht = hetstat.tau2,
                          hq = hetstat.Q, hp = hetstat.pval.Q))
      else if (print.I2 & print.tau2 & !print.Q & !print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           tau^2, ht, ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab,
                          hi = hetstat.I2, ht = hetstat.tau2,
                          hb = hetstat.Rb))
      else if (print.I2 & !print.tau2 & print.Q & !print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           chi[df]^2, hq, ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab, df = df,
                          hi = hetstat.I2,
                          hq = hetstat.Q,
                          hb = hetstat.Rb))
      else if (!print.I2 & print.tau2 & print.Q & !print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           tau^2, ht, ", ",
                           chi[df]^2, hq, ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab, df = df,
                          ht = hetstat.tau2,
                          hq = hetstat.Q,
                          hb = hetstat.Rb))
      else if (print.I2 & !print.tau2 & !print.Q & print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           italic(p), hp, ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab,
                          hi = hetstat.I2, hp = hetstat.pval.Q,
                          hb = hetstat.Rb))
      else if (!print.I2 & print.tau2 & !print.Q & print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           tau^2, ht, ", ",
                           italic(p), hp, ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab,
                          ht = hetstat.tau2,
                          hp = hetstat.pval.Q,
                          hb = hetstat.Rb))
      else if (!print.I2 & !print.tau2 & print.Q & print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           chi[df]^2, hq,
                           " (", italic(p), hp, ")",
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab, df = df,
                          hq = hetstat.Q, hp = hetstat.pval.Q,
                          hb = hetstat.Rb))      
      ##
      ## Four
      ##
      if (print.I2 & print.tau2 & print.Q & print.pval.Q & !print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           tau^2, ht, ", ",
                           chi[df]^2, hq,
                           " (", italic(p), hp, ")"),
                     list(hl = hetlab, df = df,
                          hi = hetstat.I2, ht = hetstat.tau2,
                          hq = hetstat.Q, hp = hetstat.pval.Q))
      else if (print.I2 & print.tau2 & print.Q & !print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           tau^2, ht, ", ",
                           chi[df]^2, hq, ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab, df = df,
                          hi = hetstat.I2, ht = hetstat.tau2,
                          hq = hetstat.Q, hb = hetstat.Rb))
      else if (print.I2 & print.tau2 & !print.Q & print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           tau^2, ht, ", ",
                           italic(p), hp, ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab,
                          hi = hetstat.I2, ht = hetstat.tau2,
                          hp = hetstat.pval.Q, hb = hetstat.Rb))
      else if (print.I2 & !print.tau2 & print.Q & !print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           chi[df]^2, hq,
                           " (", italic(p), hp, ")",
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab, df = df,
                          ht = hetstat.tau2,
                          hq = hetstat.Q, hp = hetstat.pval.Q,
                          hb = hetstat.Rb))
      else if (!print.I2 & print.tau2 & print.Q & print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           tau^2, ht, ", ",
                           chi[df]^2, hq,
                           " (", italic(p), hp, ")",
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab, df = df,
                          ht = hetstat.tau2,
                          hq = hetstat.Q, hp = hetstat.pval.Q,
                          hb = hetstat.Rb))
      ##
      ## Five
      ##
      else if (print.I2 & print.tau2 & print.Q & print.pval.Q & print.Rb)
        hetstat.overall <-
          substitute(paste(hl,
                           italic(I)^2, hi, ", ",
                           tau^2, ht, ", ",
                           chi[df]^2, hq,
                           " (", italic(p), hp, ")",
                           ", ",
                           italic(R)[italic(b)], hb),
                     list(hl = hetlab, df = df,
                          hi = hetstat.I2, ht = hetstat.tau2,
                          hq = hetstat.Q, hp = hetstat.pval.Q,
                          hb = hetstat.Rb))
    }
  }
  ##
  ## Text of test for subgroup differences
  ##
  if (!by) {
    test.subgroup.fixed  <- FALSE
    test.subgroup.random <- FALSE
    test.effect.subgroup.fixed <- FALSE
    test.effect.subgroup.random <- FALSE
    Q.b.fixed  <- NA
    Q.b.random <- NA
    df.Q.b     <- NA
  }
  else {
    if (!missing(test.subgroup)) {
      if (missing(test.subgroup.fixed))
        test.subgroup.fixed <- comb.fixed & test.subgroup
      else
        test.subgroup.fixed <- comb.fixed & test.subgroup.fixed
      ##
      if (missing(test.subgroup.random))
        test.subgroup.random <- comb.random & test.subgroup
      else
        test.subgroup.random <- comb.random & test.subgroup.random
    }
    else {
      if (missing(test.subgroup.fixed))
        test.subgroup.fixed <- comb.fixed & gs("test.subgroup")
      else
        test.subgroup.fixed <- comb.fixed & test.subgroup.fixed
      ##
      if (missing(test.subgroup.random))
        test.subgroup.random <- comb.random & gs("test.subgroup")
      else
        test.subgroup.random <- comb.random & test.subgroup.random
    }
    ##
    if (!missing(test.effect.subgroup)) {
      if (missing(test.effect.subgroup.fixed))
        test.effect.subgroup.fixed <- comb.fixed & test.effect.subgroup
      else
        test.effect.subgroup.fixed <- comb.fixed & test.effect.subgroup.fixed
      ##
      if (missing(test.effect.subgroup.random))
        test.effect.subgroup.random <- comb.random & test.effect.subgroup
      else
        test.effect.subgroup.random <- comb.random & test.effect.subgroup.random
    }
    else {
      if (missing(test.effect.subgroup.fixed))
        test.effect.subgroup.fixed <- comb.fixed & gs("test.effect.subgroup")
      else
        test.effect.subgroup.fixed <- comb.fixed & test.effect.subgroup.fixed
      ##
      if (missing(test.effect.subgroup.random))
        test.effect.subgroup.random <- comb.random & gs("test.effect.subgroup")
      else
        test.effect.subgroup.random <- comb.random & test.effect.subgroup.random
    }
  }
  ##
  ## Label of test for overall effect
  ##
  ##
  if (test.overall.fixed | test.overall.random) {
    pvals.overall <- format.p(c(x$pval.fixed, x$pval.random),
                              lab = TRUE, labval = "",
                              digits = digits.pval,
                              zero = if (jama) FALSE else TRUE,
                              scientific = scientific.pval,
                              lab.NA = "NA")
    zvals.overall <- format.NA(round(c(x$zval.fixed, x$zval.random),
                                     digits = digits.zval),
                               digits.zval, "NA")
    ##
    ## Remove superfluous spaces
    ##
    pvals.overall <- rmSpace(pvals.overall, end = TRUE)
    ##
    while(any(grepl("  ", pvals.overall)))
      pvals.overall <- gsub("  ", " ", pvals.overall)
    while(any(grepl("  ", zvals.overall)))
      zvals.overall <- gsub("  ", " ", zvals.overall)
  }
  ##
  if (test.overall.fixed) {
    if (print.zval) {
      if (revman5)
        text.overall.fixed  <- substitute(paste(tl,
                                                Z, hetseparator, tt,
                                                " (P", tp, ")"),
                                          list(tl = label.test.overall.fixed,
                                               hetseparator = hetseparator,
                                               tt = zvals.overall[1],
                                               tp = pvals.overall[1]))
      else if (jama)
        text.overall.fixed  <- substitute(paste(tl,
                                                italic(z), hetseparator, tt,
                                                " (", italic(P), tp, ")"),
                                          list(tl = label.test.overall.fixed,
                                               hetseparator = hetseparator,
                                               tt = zvals.overall[1],
                                               tp = pvals.overall[1]))
      else
        text.overall.fixed  <- substitute(paste(tl,
                                                italic(z), hetseparator, tt,
                                                " (", italic(p), tp, ")"),
                                          list(tl = label.test.overall.fixed,
                                               hetseparator = hetseparator,
                                               tt = zvals.overall[1],
                                               tp = pvals.overall[1]))
    }
    else {
      if (revman5)
        text.overall.fixed  <- substitute(paste(tl, " P", tp),
                                          list(tl = label.test.overall.fixed,
                                               tp = pvals.overall[1]))
      else if (jama)
        text.overall.fixed  <- substitute(paste(tl, " ", italic(P), tp),
                                          list(tl = label.test.overall.fixed,
                                               tp = pvals.overall[1]))
      else
        text.overall.fixed  <- substitute(paste(tl, " ", italic(p), tp),
                                          list(tl = label.test.overall.fixed,
                                               tp = pvals.overall[1]))     
    }
  }
  else
    text.overall.fixed <- ""
  ##
  if (test.overall.random) {
    if (print.zval) {
      if (!x$hakn) {
        if (revman5)
          text.overall.random  <- substitute(paste(tl,
                                                   Z, hetseparator, tt,
                                                   " (P", tp, ")"),
                                             list(tl = label.test.overall.random,
                                                  hetseparator = hetseparator,
                                                  tt = zvals.overall[2],
                                                  tp = pvals.overall[2]))
        else if (jama)
          text.overall.random  <- substitute(paste(tl,
                                                   italic(z), hetseparator, tt,
                                                   " (", italic(P), tp, ")"),
                                             list(tl = label.test.overall.random,
                                                  hetseparator = hetseparator,
                                                  tt = zvals.overall[2],
                                                  tp = pvals.overall[2]))
        else
          text.overall.random  <- substitute(paste(tl,
                                                   italic(z), hetseparator, tt,
                                                   " (", italic(p), tp, ")"),
                                             list(tl = label.test.overall.random,
                                                  hetseparator = hetseparator,
                                                  tt = zvals.overall[2],
                                                  tp = pvals.overall[2]))
      }
      else {
        if (revman5)
          text.overall.random  <- substitute(paste(tl,
                                                   t[df], hetseparator, tt,
                                                   " (P", tp, ")"),
                                             list(tl = label.test.overall.random,
                                                  hetseparator = hetseparator,
                                                  tt = zvals.overall[2],
                                                  tp = pvals.overall[2],
                                                  df = df))
        else if (jama)
          text.overall.random  <- substitute(paste(tl,
                                                   italic(t)[df], hetseparator, tt,
                                                   " (", italic(P), tp, ")"),
                                             list(tl = label.test.overall.random,
                                                  hetseparator = hetseparator,
                                                  tt = zvals.overall[2],
                                                  tp = pvals.overall[2],
                                                  df = df))
      else
        text.overall.random  <- substitute(paste(tl,
                                                 italic(t)[df], hetseparator, tt,
                                                 " (", italic(p), tp, ")"),
                                           list(tl = label.test.overall.random,
                                                hetseparator = hetseparator,
                                                tt = zvals.overall[2],
                                                tp = pvals.overall[2],
                                                df = df))
      }
    }
    else {
      if (revman5)
        text.overall.random  <- substitute(paste(tl, " P", tp),
                                           list(tl = label.test.overall.random,
                                                tp = pvals.overall[2]))
      else if (jama)
        text.overall.random  <- substitute(paste(tl, " ", italic(P), tp),
                                           list(tl = label.test.overall.random,
                                                tp = pvals.overall[2]))
      else
        text.overall.random  <- substitute(paste(tl, " ", italic(p), tp),
                                           list(tl = label.test.overall.random,
                                                tp = pvals.overall[2]))
    }
  }
  else
    text.overall.random <- ""
  ##
  ##
  ## Label of test for subgroup differences
  ##
  Q.bs <- c(Q.b.fixed, Q.b.random)
  ##
  hetstat.Q.bs <-
    paste(hetseparator,
          gsub(" ", "", format.NA(round(Q.bs, digits.Q), digits.Q, "NA")),
          if (!jama) ", df",
          if (!jama) hetseparator,
          if (!jama) df.Q.b,
          sep = "")
  hetstat.Q.bs <- rmSpace(hetstat.Q.bs, end = TRUE)
  ##
  hetstat.pval.Q.bs <-
    paste(format.p(pvalQ(Q.bs, df.Q.b),
                   lab = TRUE, labval = "",
                   digits = digits.pval.Q,
                   zero = if (jama) FALSE else TRUE,
                   scientific = scientific.pval,
                   lab.NA = "NA"),
          sep = "")
  ##
  ## Remove superfluous spaces
  ##
  hetstat.pval.Q.bs <- rmSpace(hetstat.pval.Q.bs, end = TRUE)
  ##
  while(any(grepl("  ", hetstat.pval.Q.bs)))
    hetstat.pval.Q.bs <- gsub("  ", " ", hetstat.pval.Q.bs)
  while(any(grepl("  ", hetstat.Q.bs)))
    hetstat.Q.bs <- gsub("  ", " ", hetstat.Q.bs)
  ##
  if (test.subgroup.fixed) {
    if (print.Q.subgroup) {
      if (revman5)
        text.subgroup.fixed <-
          substitute(paste(tl,
                           "Chi"^2, tq,
                           " (P", tp, ")"),
                     list(tl = label.test.subgroup.fixed,
                          tq = hetstat.Q.bs[1],
                          tp = hetstat.pval.Q.bs[1]))
      else if (jama)
        text.subgroup.fixed <-
          substitute(paste(tl,
                           chi[df]^2, tq,
                           " (", italic(P), tp, ")"),
                     list(tl = label.test.subgroup.fixed,
                          tq = hetstat.Q.bs[1],
                          tp = hetstat.pval.Q.bs[1],
                          df = df.Q.b))
      else
        text.subgroup.fixed <-
          substitute(paste(tl,
                           chi[df]^2, tq,
                           " (", italic(p), tp, ")"),
                     list(tl = label.test.subgroup.fixed,
                          tq = hetstat.Q.bs[1],
                          tp = hetstat.pval.Q.bs[1],
                          df = df.Q.b))
    }
    else {
      if (revman5)
        text.subgroup.fixed <-
          substitute(paste(tl, " P", tp),
                     list(tl = label.test.subgroup.fixed,
                          tp = hetstat.pval.Q.bs[1]))
      else if (jama)
        text.subgroup.fixed <-
          substitute(paste(tl, " ", italic(P), tp),
                     list(tl = label.test.subgroup.fixed,
                          tp = hetstat.pval.Q.bs[1]))
      else
        text.subgroup.fixed <-
          substitute(paste(tl, " ", italic(p), tp),
                     list(tl = label.test.subgroup.fixed,
                          tp = hetstat.pval.Q.bs[1]))
    }
  }
  else
    text.subgroup.fixed <- ""
  
  
  if (test.subgroup.random) {
    if (print.Q.subgroup) {
      if (revman5)
        text.subgroup.random <-
          substitute(paste(tl,
                           "Chi"^2, tq,
                           " (P", tp, ")"),
                     list(tl = label.test.subgroup.random,
                          tq = hetstat.Q.bs[2],
                          tp = hetstat.pval.Q.bs[2]))
      else if (jama)
        text.subgroup.random <-
          substitute(paste(tl,
                           chi[df]^2, tq,
                           " (", italic(P), tp, ")"),
                     list(tl = label.test.subgroup.random,
                          tq = hetstat.Q.bs[2],
                          tp = hetstat.pval.Q.bs[2],
                          df = df.Q.b))
      else
        text.subgroup.random <-
          substitute(paste(tl,
                           chi[df]^2, tq,
                           " (", italic(p), tp, ")"),
                     list(tl = label.test.subgroup.random,
                          tq = hetstat.Q.bs[2],
                          tp = hetstat.pval.Q.bs[2],
                          df = df.Q.b))
    }
    else {
      if (revman5)
        text.subgroup.random <-
          substitute(paste(tl, " P", tp),
                     list(tl = label.test.subgroup.random,
                          tp = hetstat.pval.Q.bs[2]))
      else if (jama)
        text.subgroup.random <-
          substitute(paste(tl, " ", italic(P), tp),
                     list(tl = label.test.subgroup.random,
                          tp = hetstat.pval.Q.bs[2]))
      else
        text.subgroup.random <-
          substitute(paste(tl, " ", italic(p), tp),
                     list(tl = label.test.subgroup.random,
                          tp = hetstat.pval.Q.bs[2]))
    }
  }
  else
    text.subgroup.random <- ""
  
  
  ##
  ##
  ## (7) Prepare data for subgroup analysis
  ##
  ##
  if (by) {
    o <- order(factor(x$bylevs, levels = bylevs))
    k.w <- x$k.w[o]
    TE.fixed.w <- x$TE.fixed.w[o]
    lower.fixed.w <- x$lower.fixed.w[o]
    upper.fixed.w <- x$upper.fixed.w[o]
    pval.fixed.w <- x$pval.fixed.w[o]
    TE.random.w <- x$TE.random.w[o]
    lower.random.w <- x$lower.random.w[o]
    upper.random.w <- x$upper.random.w[o]
    pval.random.w <- x$pval.random.w[o]
    Q.w        <- x$Q.w[o]
    I2.w       <- x$I2.w[o]
    lowI2.w    <- x$lower.I2.w[o]
    uppI2.w    <- x$upper.I2.w[o]
    Rb.w       <- x$Rb.w[o]
    lowRb.w    <- x$lower.Rb.w[o]
    uppRb.w    <- x$upper.Rb.w[o]
    tau.w      <- x$tau.w[o]
    w.fixed.w  <- x$w.fixed.w[o]
    w.random.w <- x$w.random.w[o]
    e.e.w <- if (metaprop | metarate) x$event.w[o] else x$event.e.w[o]
    t.e.w <- if (metainc | metarate) x$time.e.w[o] else rep(NA, n.by)
    n.e.w <- if (metacor | metaprop) x$n.w[o] else x$n.e.w[o]
    e.c.w <- x$event.c.w[o]
    t.c.w <- if (metainc) x$time.c.w[o] else rep(NA, n.by)
    n.c.w <- x$n.c.w[o]
    n.harmonic.mean.w <- x$n.harmonic.mean.w[o]
    t.harmonic.mean.w <- x$t.harmonic.mean.w[o]
    k.all.w <- x$k.all.w[o]
    ##
    if (allstudies)
      sel <- 1:length(k.w)
    else
      sel <- k.w > 0
    k.w <- k.w[sel]
    TE.fixed.w <- TE.fixed.w[sel]
    lower.fixed.w <- lower.fixed.w[sel]
    upper.fixed.w <- upper.fixed.w[sel]
    pval.fixed.w <- pval.fixed.w[sel]
    TE.random.w <- TE.random.w[sel]
    lower.random.w <- lower.random.w[sel]
    upper.random.w <- upper.random.w[sel]
    pval.random.w <- pval.random.w[sel]
    Q.w        <- Q.w[sel]
    I2.w       <- I2.w[sel]
    lowI2.w    <- lowI2.w[sel]
    uppI2.w    <- uppI2.w[sel]
    Rb.w       <- Rb.w[sel]
    lowRb.w    <- lowRb.w[sel]
    uppRb.w    <- uppRb.w[sel]
    tau.w      <- tau.w[sel]
    w.fixed.w  <- w.fixed.w[sel]
    w.random.w <- w.random.w[sel]
    e.e.w <- e.e.w[sel]
    t.e.w <- t.e.w[sel]
    n.e.w <- n.e.w[sel]
    e.c.w <- e.c.w[sel]
    t.c.w <- t.c.w[sel]
    n.c.w <- n.c.w[sel]
    n.harmonic.mean.w <- n.harmonic.mean.w[sel]
    t.harmonic.mean.w <- t.harmonic.mean.w[sel]
    k.all.w <- k.all.w[sel]
    bylevs <- bylevs[sel]
    ##
    if (!metainf.metacum & comb.fixed) {
      if (!overall) {
        i <- 0
        for (bylev.i in bylevs) {
          i <- i + 1
          sel.i <- byvar == bylev.i
          x$w.fixed[sel.i] <- x$w.fixed[sel.i] / w.fixed.w[i]
        }
        w.fixed.w.p <- ifelse(is.na(w.fixed.w), NA, 100)
      }
      else {
        if (!all(is.na(w.fixed.w)) && sum(w.fixed.w) > 0)
          w.fixed.w.p <- round(100 * w.fixed.w / sum(w.fixed.w, na.rm = TRUE), digits.weight)
        else
          w.fixed.w.p <- w.fixed.w
      }
    }
    else {
      TE.fixed.w <- lower.fixed.w <- upper.fixed.w <- rep(NA, n.by)
      ##
      w.fixed.w.p <- rep(NA, n.by)
      if (missing(text.fixed.w))
        text.fixed.w <- rep("Overall", n.by)
    }
    ##
    if (!metainf.metacum & comb.random) {
      if (!overall) {
        i <- 0
        for (bylev.i in bylevs) {
          i <- i + 1
          sel.i <- byvar == bylev.i
          x$w.random[sel.i] <- x$w.random[sel.i] / w.random.w[i]
        }
        w.random.w.p <- ifelse(is.na(w.random.w), NA, 100)
      }
      else {
        if (!all(is.na(w.random.w)) && sum(w.random.w) > 0)
          w.random.w.p <- round(100 * w.random.w / sum(w.random.w, na.rm = TRUE), digits.weight)
        else
          w.random.w.p <- w.random.w
      }
    }
    else {
      TE.random.w <- lower.random.w <- upper.random.w <- rep(NA, n.by)
      ##
      w.random.w.p <- rep(NA, n.by)
      text.random.w <- rep("", n.by)
    }
    ##
    hetstat.w <- vector("list", n.by)
    for (i in 1:n.by)
      hetstat.w[[i]] <- ""
    ##
    if (hetstat) {
      ##
      hetstat.I2.w <-
        paste(hetseparator,
              round(100 * I2.w, digits.I2), "%",
              if (print.I2.ci)
                ifelse(k.w > 2,
                       paste(" ",
                             p.ci(paste(round(100 * lowI2.w, digits.I2), "%", sep = ""),
                                  paste(round(100 * uppI2.w, digits.I2), "%", sep = "")),
                             sep = ""),
                       ""),
              sep = "")
      ##
      hetstat.tau2.w <-
        paste(hetseparator,
              ifelse(is.na(tau.w), "NA",
                     ifelse(tau.w == 0, "0",
                            format.tau(tau.w^2, digits = digits.tau2))),
              sep = "")
      ##
      hetstat.Q.w <-
        paste(hetseparator,
              round(Q.w, digits.Q),
              if (revman5) ", df",
              if (revman5) hetseparator,
              if (revman5) k.w - 1,
              sep = "")
      ##
      hetstat.pval.Q.w <-
        paste(format.p(pvalQ(Q.w, k.w - 1),
                       lab = TRUE, labval = "",
                       digits = digits.pval.Q,
                       zero = if (jama) FALSE else TRUE,
                       scientific = scientific.pval,
                       lab.NA = "NA"),
              sep = "")
      ##
      hetstat.Rb.w <-
        paste(hetseparator,
              round(100 * Rb.w, digits.I2), "%",
              if (print.Rb.ci)
                ifelse(k.w > 2,
                       paste(" ",
                             p.ci(paste(round(100 * lowRb.w, digits.I2), "%", sep = ""),
                                  paste(round(100 * uppRb.w, digits.I2), "%", sep = "")),
                             sep = ""),
                       ""),
              sep = "")
      ##
      ## Remove superfluous spaces
      ##
      hetstat.pval.Q.w <- rmSpace(hetstat.pval.Q.w, end = TRUE)
      ##
      while(any(grepl("  ", hetstat.I2.w)))
        hetstat.I2.w <- gsub("  ", " ", hetstat.I2.w)
      while(any(grepl("  ", hetstat.tau2.w)))
        hetstat.tau2.w <- gsub("  ", " ", hetstat.tau2.w)
      while(any(grepl("  ", hetstat.Q.w)))
        hetstat.Q.w <- gsub("  ", " ", hetstat.Q.w)
      while(any(grepl("  ", hetstat.pval.Q.w)))
        hetstat.pval.Q.w <- gsub("  ", " ", hetstat.pval.Q.w)
      while(any(grepl("  ", hetstat.Rb.w)))
        hetstat.Rb.w <- gsub("  ", " ", hetstat.Rb.w)
      ##
      for (i in 1:n.by) {
        if (revman5)
          hetstat.w[[i]] <- substitute(paste(hl,
                                             "Tau"^2, ht, "; ",
                                             "Chi"^2, hq,
                                             " (",
                                             P, hp,
                                             "); ",
                                             I^2, hi),
                                       list(hl = hetlab,
                                            hi = hetstat.I2.w[i],
                                            ht = hetstat.tau2.w[i],
                                            hq = hetstat.Q.w[i],
                                            hp = hetstat.pval.Q.w[i]))
        else if (jama)
          hetstat.w[[i]] <- substitute(paste(hl,
                                             chi[df]^2, hq,
                                             " (",
                                             italic(P), hp,
                                             "), ",
                                             italic(I)^2, hi
                                             ),
                                       list(hl = hetlab,
                                            hi = hetstat.I2.w[i],
                                            hq = hetstat.Q.w[i],
                                            hp = hetstat.pval.Q.w[i],
                                            df = k.w[i] - 1))
        else {
          ##
          ## One
          ##
          if (print.I2 & !print.tau2 & !print.Q & !print.pval.Q & !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, italic(I)^2, hi),
                         list(hl = hetlab, hi = hetstat.I2.w[i]))
          else if (!print.I2 & print.tau2 & !print.Q & !print.pval.Q & !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, tau^2, ht),
                         list(hl = hetlab, ht = hetstat.tau2.w[i]))
          else if (!print.I2 & !print.tau2 & print.Q & !print.pval.Q & !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, chi[df]^2, hq),
                         list(hl = hetlab, df = k.w[i] - 1, hq = hetstat.Q.w[i]))
          else if (!print.I2 & !print.tau2 & !print.Q & print.pval.Q & !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, italic(p), hp),
                         list(hl = hetlab, hp = hetstat.pval.Q.w[i]))
          else if (!print.I2 & !print.tau2 & !print.Q & !print.pval.Q & print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, italic(R)[italic(b)], hb),
                         list(hl = hetlab, hetstat.Rb.w[i]))
          ##
          ## Two
          ##
          else if (print.I2 & print.tau2 & !print.Q & !print.pval.Q & !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, italic(I)^2, hi,
                               ", ",
                               tau^2, ht),
                         list(hl = hetlab,
                              hi = hetstat.I2.w[i], ht = hetstat.tau2.w[i]))
          else if (print.I2 & !print.tau2 & print.Q & !print.pval.Q & !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, italic(I)^2, hi,
                               ", ",
                               chi[df]^2, hq),
                         list(hl = hetlab,
                              hi = hetstat.I2.w[i], hq = hetstat.Q.w[i]))
          else if (print.I2 & !print.tau2 & !print.Q & print.pval.Q & !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, italic(I)^2, hi,
                               ", ",
                               italic(p), hp),
                         list(hl = hetlab,
                              hi = hetstat.I2.w[i], hp = hetstat.pval.Q.w[i]))
          else if (print.I2 & !print.tau2 & !print.Q & !print.pval.Q & print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, italic(I)^2, hi,
                               ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab,
                              hi = hetstat.I2.w[i], hetstat.Rb.w[i]))
          else if (!print.I2 & print.tau2 & print.Q & !print.pval.Q & !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, tau^2, ht,
                               ", ",
                               chi[df]^2, hq),
                         list(hl = hetlab, df = k.w[i] - 1,
                              ht = hetstat.tau2.w[i], hq = hetstat.Q.w[i]))
          else if (!print.I2 & print.tau2 & !print.Q & print.pval.Q & !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, tau^2, ht,
                               ", ",
                               italic(p), hp),
                         list(hl = hetlab,
                              ht = hetstat.tau2.w[i], hp = hetstat.pval.Q.w[i]))
          else if (!print.I2 & print.tau2 & !print.Q & !print.pval.Q & print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, tau^2, ht,
                               ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab,
                              ht = hetstat.tau2.w[i], hetstat.Rb.w[i]))
          else if (!print.I2 & !print.tau2 & print.Q & print.pval.Q & !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, chi[df]^2, hq,
                               " (",
                               italic(p), hp, ")"),
                         list(hl = hetlab, df = k.w[i] - 1,
                              hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i]))
          else if (!print.I2 & !print.tau2 & print.Q & !print.pval.Q & print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, chi[df]^2, hq,
                               ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab, df = k.w[i] - 1,
                              hq = hetstat.Q.w[i], hetstat.Rb.w[i]))
          else if (!print.I2 & !print.tau2 & !print.Q & print.pval.Q & print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl, italic(p), hp,
                               ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab,
                              hp = hetstat.pval.Q.w[i], hetstat.Rb.w[i]))
          ##
          ## Three
          ##
          else if (print.I2 & print.tau2 & print.Q & !print.pval.Q & !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               italic(I)^2, hi, ", ",
                               tau^2, ht, ", ",
                               chi[df]^2, hq),
                         list(hl = hetlab, df = k.w[i] - 1,
                              hi = hetstat.I2.w[i], ht = hetstat.tau2.w[i],
                              hq = hetstat.Q.w[i]))
          else if (print.I2 & print.tau2 & !print.Q & print.pval.Q & !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               italic(I)^2, hi, ", ",
                               tau^2, ht, ", ",
                               italic(p), hp),
                         list(hl = hetlab,
                              hi = hetstat.I2.w[i], ht = hetstat.tau2.w[i],
                              hp = hetstat.pval.Q.w[i]))
          else if (print.I2 & !print.tau2 & print.Q & print.pval.Q & !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               italic(I)^2, hi, ", ",
                               chi[df]^2, hq,
                               " (", italic(p), hp, ")"),
                         list(hl = hetlab, df = k.w[i] - 1,
                              hi = hetstat.I2.w[i],
                              hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i]))
          else if (!print.I2 & print.tau2 & print.Q & print.pval.Q & !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               tau^2, ht, ", ",
                               chi[df]^2, hq,
                               " (", italic(p), hp, ")"),
                         list(hl = hetlab, df = k.w[i] - 1,
                              ht = hetstat.tau2.w[i],
                              hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i]))
          else if (print.I2 & print.tau2 & !print.Q & !print.pval.Q & print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               italic(I)^2, hi, ", ",
                               tau^2, ht, ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab,
                              hi = hetstat.I2.w[i], ht = hetstat.tau2.w[i],
                              hetstat.Rb.w[i]))
          else if (print.I2 & !print.tau2 & print.Q & !print.pval.Q & print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               italic(I)^2, hi, ", ",
                               chi[df]^2, hq, ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab, df = k.w[i] - 1,
                              hi = hetstat.I2.w[i],
                              hq = hetstat.Q.w[i],
                              hetstat.Rb.w[i]))
          else if (!print.I2 & print.tau2 & print.Q & !print.pval.Q & print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               tau^2, ht, ", ",
                               chi[df]^2, hq, ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab, df = k.w[i] - 1,
                              ht = hetstat.tau2.w[i],
                              hq = hetstat.Q.w[i],
                              hetstat.Rb.w[i]))
          else if (print.I2 & !print.tau2 & !print.Q & print.pval.Q & print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               italic(I)^2, hi, ", ",
                               italic(p), hp, ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab,
                              hi = hetstat.I2.w[i], hp = hetstat.pval.Q.w[i],
                              hetstat.Rb.w[i]))
          else if (!print.I2 & print.tau2 & !print.Q & print.pval.Q & print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               tau^2, ht, ", ",
                               italic(p), hp, ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab,
                              ht = hetstat.tau2.w[i],
                              hp = hetstat.pval.Q.w[i],
                              hetstat.Rb.w[i]))
          else if (!print.I2 & !print.tau2 & print.Q & print.pval.Q & print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               chi[df]^2, hq,
                               " (", italic(p), hp, ")",
                               ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab, df = k.w[i] - 1,
                              hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i],
                              hetstat.Rb.w[i]))      
          ##
          ## Four
          ##
          if (print.I2 & print.tau2 & print.Q & print.pval.Q & !print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               italic(I)^2, hi, ", ",
                               tau^2, ht, ", ",
                               chi[df]^2, hq,
                               " (", italic(p), hp, ")"),
                         list(hl = hetlab, df = k.w[i] - 1,
                              hi = hetstat.I2.w[i], ht = hetstat.tau2.w[i],
                              hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i]))
          else if (print.I2 & print.tau2 & print.Q & !print.pval.Q & print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               italic(I)^2, hi, ", ",
                               tau^2, ht, ", ",
                               chi[df]^2, hq, ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab, df = k.w[i] - 1,
                              hi = hetstat.I2.w[i], ht = hetstat.tau2.w[i],
                              hq = hetstat.Q.w[i], hetstat.Rb.w[i]))
          else if (print.I2 & print.tau2 & !print.Q & print.pval.Q & print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               italic(I)^2, hi, ", ",
                               tau^2, ht, ", ",
                               italic(p), hp, ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab,
                              hi = hetstat.I2.w[i], ht = hetstat.tau2.w[i],
                              hp = hetstat.pval.Q.w[i], hetstat.Rb.w[i]))
          else if (print.I2 & !print.tau2 & print.Q & !print.pval.Q & print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               italic(I)^2, hi, ", ",
                               chi[df]^2, hq,
                               " (", italic(p), hp, ")",
                               ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab, df = k.w[i] - 1,
                              ht = hetstat.tau2.w[i],
                              hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i],
                              hetstat.Rb.w[i]))
          else if (!print.I2 & print.tau2 & print.Q & print.pval.Q & print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               tau^2, ht, ", ",
                               chi[df]^2, hq,
                               " (", italic(p), hp, ")",
                               ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab, df = k.w[i] - 1,
                              ht = hetstat.tau2.w[i],
                              hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i],
                              hetstat.Rb.w[i]))
          ##
          ## Five
          ##
          else if (print.I2 & print.tau2 & print.Q & print.pval.Q & print.Rb)
            hetstat.w[[i]] <-
              substitute(paste(hl,
                               italic(I)^2, hi, ", ",
                               tau^2, ht, ", ",
                               chi[df]^2, hq,
                               " (", italic(p), hp, ")",
                               ", ",
                               italic(R)[italic(b)], hb),
                         list(hl = hetlab, df = k.w[i] - 1,
                              hi = hetstat.I2.w[i], ht = hetstat.tau2.w[i],
                              hq = hetstat.Q.w[i], hp = hetstat.pval.Q.w[i],
                              hetstat.Rb.w[i]))
        }
        if (!is.logical(text.subgroup.nohet) & k.w[i] < 2)
          hetstat.w[[i]] <- paste(hetlab, text.subgroup.nohet, sep = "")
      }
    }
    ##
    TE.w <- c(TE.fixed.w, TE.random.w, rep(NA, 3 * n.by))
    lowTE.w <- c(lower.fixed.w, lower.random.w, rep(NA, 3 * n.by))
    uppTE.w <- c(upper.fixed.w, upper.random.w, rep(NA, 3 * n.by))
    n.harmonic.mean.w <- c(n.harmonic.mean.w, n.harmonic.mean.w, rep(NA, 3 * n.by))
    t.harmonic.mean.w <- c(t.harmonic.mean.w, t.harmonic.mean.w, rep(NA, 3 * n.by))
    weight.w.p <- c(w.fixed.w.p, w.random.w.p, rep(NA, 3 * n.by))
    ##
    test.fixed.w <- ""
    ##
    ## Label of test for effect in subgroups
    ##
    if (test.effect.subgroup.fixed | test.effect.subgroup.random) {
      pvals.effect.w <- format.p(c(x$pval.fixed.w, x$pval.random.w),
                                 lab = TRUE, labval = "",
                                 digits = digits.pval,
                                 zero = if (jama) FALSE else TRUE,
                                 scientific = scientific.pval,
                                 lab.NA = "NA")
      zvals.effect.w <- format.NA(round(c(x$zval.fixed.w, x$zval.random.w),
                                        digits = digits.zval),
                                  digits.zval, "NA")
      ##
      ## Remove superfluous spaces
      ##
      pvals.effect.w <- rmSpace(pvals.effect.w, end = TRUE)
      ##
      while(any(grepl("  ", pvals.effect.w)))
        pvals.effect.w <- gsub("  ", " ", pvals.effect.w)
      while(any(grepl("  ", zvals.effect.w)))
        zvals.effect.w <- gsub("  ", " ", zvals.effect.w)
    }
    ##
    if (test.effect.subgroup.fixed) {
      text.effect.subgroup.fixed <- vector("list", n.by)
      for (i in 1:n.by) {
        if (print.zval) {
          if (revman5)
            text.effect.subgroup.fixed[[i]]  <- substitute(paste(tl,
                                                                 Z, hetseparator, tt,
                                                                 " (P", tp, ")"),
                                                           list(tl = label.test.effect.subgroup.fixed,
                                                                hetseparator = hetseparator,
                                                                tt = zvals.effect.w[i],
                                                                tp = rmSpace(pvals.effect.w[i], end = TRUE)))
          else if (jama)
            text.effect.subgroup.fixed[[i]]  <- substitute(paste(tl,
                                                                 italic(z), hetseparator, tt,
                                                                 " (", italic(P), tp, ")"),
                                                           list(tl = label.test.effect.subgroup.fixed,
                                                                hetseparator = hetseparator,
                                                                tt = zvals.effect.w[i],
                                                                tp = rmSpace(pvals.effect.w[i], end = TRUE)))
          else
            text.effect.subgroup.fixed[[i]]  <- substitute(paste(tl,
                                                                 italic(z), hetseparator, tt,
                                                                 " (", italic(p), tp, ")"),
                                                           list(tl = label.test.effect.subgroup.fixed,
                                                                hetseparator = hetseparator,
                                                                tt = zvals.effect.w[i],
                                                                tp = rmSpace(pvals.effect.w[i], end = TRUE)))
        }
        else {
          if (revman5)
            text.effect.subgroup.fixed[[i]]  <- substitute(paste(tl,
                                                                 " P", tp),
                                                           list(tl = label.test.effect.subgroup.fixed,
                                                                tp = rmSpace(pvals.effect.w[i], end = TRUE)))
          else if (jama)
            text.effect.subgroup.fixed[[i]]  <- substitute(paste(tl,
                                                                 " ", italic(P), tp),
                                                           list(tl = label.test.effect.subgroup.fixed,
                                                                tp = rmSpace(pvals.effect.w[i], end = TRUE)))
          else
            text.effect.subgroup.fixed[[i]]  <- substitute(paste(tl,
                                                                 " ", italic(p), tp),
                                                           list(tl = label.test.effect.subgroup.fixed,
                                                                tp = rmSpace(pvals.effect.w[i], end = TRUE)))
        }
      }
    }
    else {
      text.effect.subgroup.fixed <- vector("list", n.by)
      for (i in 1:n.by)
        text.effect.subgroup.fixed[[i]] <- ""
    }
    ##
    if (test.effect.subgroup.random) {
      text.effect.subgroup.random <- vector("list", n.by)
      for (i in 1:n.by) {
        if (print.zval) {
          if (!x$hakn) {
            if (revman5)
              text.effect.subgroup.random[[i]]  <- substitute(paste(tl,
                                                                    Z, hetseparator, tt,
                                                                    " (P", tp, ")"),
                                                              list(tl = label.test.effect.subgroup.random,
                                                                   hetseparator = hetseparator,
                                                                   tt = zvals.effect.w[n.by + i],
                                                                   tp = rmSpace(pvals.effect.w[n.by + i], end = TRUE)))
            else if (jama)
              text.effect.subgroup.random[[i]]  <- substitute(paste(tl,
                                                                    italic(z), hetseparator, tt,
                                                                    " (", italic(P), tp, ")"),
                                                              list(tl = label.test.effect.subgroup.random,
                                                                   hetseparator = hetseparator,
                                                                   tt = zvals.effect.w[n.by + i],
                                                                   tp = rmSpace(pvals.effect.w[n.by + i], end = TRUE)))
            else
              text.effect.subgroup.random[[i]]  <- substitute(paste(tl,
                                                                    italic(z), hetseparator, tt,
                                                                    " (", italic(p), tp, ")"),
                                                              list(tl = label.test.effect.subgroup.random,
                                                                   hetseparator = hetseparator,
                                                                   tt = zvals.effect.w[n.by + i],
                                                                   tp = rmSpace(pvals.effect.w[n.by + i], end = TRUE)))
          }
          else {
            if (revman5)
              text.effect.subgroup.random[[i]]  <- substitute(paste(tl,
                                                                    t[df], hetseparator, tt,
                                                                    " (P", tp, ")"),
                                                              list(tl = label.test.effect.subgroup.random,
                                                                   hetseparator = hetseparator,
                                                                   tt = zvals.effect.w[n.by + i],
                                                                   tp = rmSpace(pvals.effect.w[n.by + i], end = TRUE),
                                                                   df = k.w - 1))
            else if (jama)
              text.effect.subgroup.random[[i]]  <- substitute(paste(tl,
                                                                    italic(t)[df], hetseparator, tt,
                                                                    " (", italic(P), tp, ")"),
                                                              list(tl = label.test.effect.subgroup.random,
                                                                   hetseparator = hetseparator,
                                                                   tt = zvals.effect.w[n.by + i],
                                                                   tp = rmSpace(pvals.effect.w[n.by + i], end = TRUE),
                                                                   df = k.w - 1))
            else
              text.effect.subgroup.random[[i]]  <- substitute(paste(tl,
                                                                    italic(t)[df], hetseparator, tt,
                                                                    " (", italic(p), tp, ")"),
                                                              list(tl = label.test.effect.subgroup.random,
                                                                   hetseparator = hetseparator,
                                                                   tt = zvals.effect.w[n.by + i],
                                                                   tp = rmSpace(pvals.effect.w[n.by + i], end = TRUE),
                                                                   df = k.w - 1))
          }
        }
        else {
          if (revman5)
            text.effect.subgroup.random[[i]]  <- substitute(paste(tl,
                                                                  " P", tp),
                                                            list(tl = label.test.effect.subgroup.random,
                                                                 tp = rmSpace(pvals.effect.w[n.by + i], end = TRUE)))
          else if (jama)
            text.effect.subgroup.random[[i]]  <- substitute(paste(tl,
                                                                  " ", italic(P)),
                                                            list(tl = label.test.effect.subgroup.random,
                                                                 tp = rmSpace(pvals.effect.w[n.by + i], end = TRUE)))
          else
            text.effect.subgroup.random[[i]]  <- substitute(paste(tl,
                                                                  " ", italic(p), tp),
                                                            list(tl = label.test.effect.subgroup.random,
                                                                 tp = rmSpace(pvals.effect.w[n.by + i], end = TRUE)))
        }
      }
    }
    else {
      text.effect.subgroup.random <- vector("list", n.by)
      for (i in 1:n.by)
        text.effect.subgroup.random[[i]] <- ""
    }
  }
  
  
  ##
  ##
  ## (8) Backtransform data
  ##
  ##
  TE.orig <- TE
  ##
  if (backtransf) {
    ##
    ## Freeman-Tukey Arcsin transformation
    ##
    if (metainf.metacum) {
      if (sm == "IRFT") {
        npft <- x$t.harmonic.mean
        npft.ma <- x$t.harmonic.mean.ma
      }
      else {
        npft <- x$n.harmonic.mean
        npft.ma <- x$n.harmonic.mean.ma
      }
    }
    else {
      if (sm == "IRFT") {
        npft <- x$time
        npft.ma <- 1 / mean(1 / x$time)
      }
      else {
        npft <- x$n
        npft.ma <- 1 / mean(1 / x$n)
      }
    }
    ##
    ## Individual study results
    ##
    if (metaprop) {
      TE <- x$event.e / x$n.e
    }
    ## Relative effect measures will be back transformed later
    else if (!is.relative.effect(sm)) {
      TE <- backtransf(TE, sm, "mean", npft)
      lowTE <- backtransf(lowTE, sm, "lower", npft)
      uppTE <- backtransf(uppTE, sm, "upper", npft)
    }
    ##
    ## Results of meta-analysis
    ##
    if (!is.relative.effect(sm)) {
      TE.fixed    <- backtransf(TE.fixed, sm, "mean",
                                npft.ma, warn = comb.fixed)
      lowTE.fixed <- backtransf(lowTE.fixed, sm, "lower",
                                npft.ma, warn = comb.fixed)
      uppTE.fixed <- backtransf(uppTE.fixed, sm, "upper",
                                npft.ma, warn = comb.fixed)
      ##
      TE.random <- backtransf(TE.random, sm, "mean",
                              npft.ma, warn = comb.random)
      lowTE.random <- backtransf(lowTE.random, sm, "lower",
                                 npft.ma, warn = comb.random)
      uppTE.random <- backtransf(uppTE.random, sm, "upper",
                                 npft.ma, warn = comb.random)
      ##
      if (!metainf.metacum) {
        lowTE.predict <- backtransf(lowTE.predict, sm, "lower",
                                    npft.ma, warn = prediction)
        uppTE.predict <- backtransf(uppTE.predict, sm, "upper",
                                    npft.ma, warn = prediction)
      }
      ##
      if (by) {
        if (sm == "IRFT")
          npft.w <- t.harmonic.mean.w
        else
          npft.w <- n.harmonic.mean.w
        ##
        TE.w     <- backtransf(TE.w, sm, "mean", npft.w)
        lowTE.w  <- backtransf(lowTE.w, sm, "lower", npft.w)
        uppTE.w  <- backtransf(uppTE.w, sm, "upper", npft.w)
      }
    }
    ##
    ## Apply argument 'pscale' to proportions
    ##
    if (sm %in% c("PLOGIT", "PLN", "PRAW", "PAS", "PFT")) {
      TE <- pscale * TE
      lowTE <- pscale * lowTE
      uppTE <- pscale * uppTE
      ##
      TE.fixed    <- pscale * TE.fixed
      lowTE.fixed <- pscale * lowTE.fixed
      uppTE.fixed <- pscale * uppTE.fixed
      ##
      TE.random    <- pscale * TE.random
      lowTE.random <- pscale * lowTE.random
      uppTE.random <- pscale * uppTE.random
      ##
      lowTE.predict <- pscale * lowTE.predict
      uppTE.predict <- pscale * uppTE.predict
      ##
      if (by) {
        TE.w    <- pscale * TE.w
        lowTE.w <- pscale * lowTE.w
        uppTE.w <- pscale * uppTE.w
      }
    }
  }
  ##
  ## Apply argument 'irscale' to rates
  ##
  if (sm %in% c("IR", "IRLN", "IRS", "IRFT")) {
    TE <- irscale * TE
    lowTE <- irscale * lowTE
    uppTE <- irscale * uppTE
    ##
    TE.fixed    <- irscale * TE.fixed
    lowTE.fixed <- irscale * lowTE.fixed
    uppTE.fixed <- irscale * uppTE.fixed
    ##
    TE.random    <- irscale * TE.random
    lowTE.random <- irscale * lowTE.random
    uppTE.random <- irscale * uppTE.random
    ##
    lowTE.predict <- irscale * lowTE.predict
    uppTE.predict <- irscale * uppTE.predict
    ##
    if (by) {
      TE.w    <- irscale * TE.w
      lowTE.w <- irscale * lowTE.w
      uppTE.w <- irscale * uppTE.w
    }
  }
  ##
  ## Exclude study results from forest plot
  ##
  TE.exclude <- TE
  lowTE.exclude <- lowTE
  uppTE.exclude <- uppTE
  ##
  if (!is.null(x$exclude)) {
    TE.exclude[exclude] <- NA
    lowTE.exclude[exclude] <- NA
    uppTE.exclude[exclude] <- NA
  }
  ##
  if (!comb.fixed) {
    TE.fixed    <- NA
    lowTE.fixed <- NA
    uppTE.fixed <- NA
  }
  ##
  if (!comb.random) {
    TE.random    <- NA
    lowTE.random <- NA
    uppTE.random <- NA
  }
  ##
  if (!prediction) {
    lowTE.predict <- NA
    uppTE.predict <- NA
  }
  ##
  if (!metainf.metacum) {
    if (by & !overall)
      w.fixed.p <- round(100 * x$w.fixed, digits.weight)
    else {
      if (!all(is.na(x$w.fixed)) && sum(x$w.fixed) > 0)
        w.fixed.p <- round(100 * x$w.fixed / sum(x$w.fixed, na.rm = TRUE), digits.weight)
      else
        w.fixed.p <- x$w.fixed
    }
    ##
    if (by & !overall)
      w.random.p <- round(100 * x$w.random, digits.weight)
    else {
      if (!all(is.na(x$w.random)) && sum(x$w.random) > 0)
        w.random.p <- round(100 * x$w.random / sum(x$w.random, na.rm = TRUE), digits.weight)
      else
        w.random.p <- x$w.random
    }
  }
  else {
    w.fixed.p  <- rep(NA, length(TE))
    w.random.p <- rep(NA, length(TE))
  }


  ##
  ##
  ## (9) Determine column labels
  ##
  ##
  labs <- list()
  ##
  if (missing(leftlabs) || length(leftcols) != length(leftlabs)) {
    for (i in seq(along = leftcols)) {
      j <- match(leftcols[i], colnames)
      if (!is.na(j))
        labs[[paste("lab.", leftcols[i], sep = "")]] <- labnames[j]
    }
  }
  else if (length(leftcols) == length(leftlabs)) {
    for (i in seq(along = leftcols)) {
      j <- match(leftcols[i], colnames)
      if (!is.na(leftlabs[i]))
        labs[[paste("lab.", leftcols[i], sep = "")]] <- leftlabs[i]
      else
        if (!is.na(j))
          labs[[paste("lab.", leftcols[i], sep = "")]] <- labnames[j]
    }
  }
  ##
  if (missing(rightlabs) || length(rightcols) != length(rightlabs)) {
    for (i in seq(along = rightcols)) {
      j <- match(rightcols[i], colnames)
      if (!is.na(j))
        labs[[paste("lab.", rightcols[i], sep = "")]] <- labnames[j]
    }
  }
  else if (length(rightcols) == length(rightlabs)) {
    for (i in seq(along = rightcols)) {
      j <- match(rightcols[i], colnames)
      if (!is.na(rightlabs[i]))
        labs[[paste("lab.", rightcols[i], sep = "")]] <- rightlabs[i]
      else
        if (!is.na(j))
          labs[[paste("lab.", rightcols[i], sep = "")]] <- labnames[j]
    }
  }
  ##
  if (!slab)
    labs[["lab.studlab"]] <- ""
  ##
  ## Check for "%" in weight labels
  ##
  w.fixed.percent <- TRUE
  if (!is.null(labs[["lab.w.fixed"]])) {
    if (grepl("%", labs[["lab.w.fixed"]]))
      w.fixed.percent <- FALSE
  }
  ##
  w.random.percent <- TRUE
  if (!is.null(labs[["lab.w.random"]])) {
    if (grepl("%", labs[["lab.w.random"]]))
      w.random.percent <- FALSE
  }
  ## "studlab", "TE", "seTE",
  ## "n.e", "n.c",
  ## "event.e", "event.c",
  ## "mean.e", "mean.c",
  ## "sd.e", "sd.c",
  ## "cor",
  ## "time.e", "time.c",
  ##
  ## Check for "\n" in column studlab
  ##
  clines <- twolines(labs[["lab.studlab"]], "studlab")
  ##
  if (clines$newline) {
    newline.studlab <- TRUE
    labs[["lab.studlab"]] <- clines$bottom
    add.studlab <- clines$top
    longer.studlab <- clines$longer
  }
  else {
    newline.studlab <- FALSE
    longer.studlab <- labs[["lab.studlab"]]
  }
  ##
  ## Check for "\n" in column effect
  ##
  clines <- twolines(labs[["lab.effect"]], "effect")
  ##
  if (clines$newline) {
    newline.effect <- TRUE
    labs[["lab.effect"]] <- clines$bottom
    add.effect <- clines$top
    longer.effect <- clines$longer
  }
  else {
    newline.effect <- FALSE
    longer.effect <- labs[["lab.effect"]]
  }
  ##
  ## Check for "\n" in column ci
  ##
  clines <- twolines(labs[["lab.ci"]], "ci")
  ##
  if (clines$newline) {
    newline.ci <- TRUE
    labs[["lab.ci"]] <- clines$bottom
    add.ci <- clines$top
    longer.ci <- clines$longer
  }
  else {
    newline.ci <- FALSE
    longer.ci <- labs[["lab.ci"]]
  }
  ##
  ## Check for "\n" in column effect.ci
  ##
  clines <- twolines(labs[["lab.effect.ci"]], "effect.ci")
  ##
  if (clines$newline) {
    newline.effect.ci <- TRUE
    labs[["lab.effect.ci"]] <- clines$bottom
    add.effect.ci <- clines$top
    longer.effect.ci <- clines$longer
  }
  else {
    newline.effect.ci <- FALSE
    longer.effect.ci <- labs[["lab.effect.ci"]]
  }
  ##
  ## Check for "\n" in column w.fixed
  ##
  clines <- twolines(labs[["lab.w.fixed"]], "w.fixed")
  ##
  if (clines$newline) {
    newline.w.fixed <- TRUE
    labs[["lab.w.fixed"]] <- clines$bottom
    add.w.fixed <- clines$top
    longer.w.fixed <- clines$longer
  }
  else {
    newline.w.fixed <- FALSE
    longer.w.fixed <- labs[["lab.w.fixed"]]
  }
  ##
  ## Check for "\n" in column w.random
  ##
  clines <- twolines(labs[["lab.w.random"]], "w.random")
  ##
  if (clines$newline) {
    newline.w.random <- TRUE
    labs[["lab.w.random"]] <- clines$bottom
    add.w.random <- clines$top
    longer.w.random <- clines$longer
  }
  else {
    newline.w.random <- FALSE
    longer.w.random <- labs[["lab.w.random"]]
  }
  ##
  ## Check for "\n" in column TE
  ##
  clines <- twolines(labs[["lab.TE"]], "TE")
  ##
  if (clines$newline) {
    newline.TE <- TRUE
    labs[["lab.TE"]] <- clines$bottom
    add.TE <- clines$top
    longer.TE <- clines$longer
  }
  else {
    newline.TE <- FALSE
    longer.TE <- labs[["lab.TE"]]
  }
  ##
  ## Check for "\n" in column seTE
  ##
  clines <- twolines(labs[["lab.seTE"]], "seTE")
  ##
  if (clines$newline) {
    newline.seTE <- TRUE
    labs[["lab.seTE"]] <- clines$bottom
    add.seTE <- clines$top
    longer.seTE <- clines$longer
  }
  else {
    newline.seTE <- FALSE
    longer.seTE <- labs[["lab.seTE"]]
  }
  ##
  ## Check for "\n" in column n.e
  ##
  clines <- twolines(labs[["lab.n.e"]], "n.e")
  ##
  if (clines$newline) {
    newline.n.e <- TRUE
    labs[["lab.n.e"]] <- clines$bottom
    add.n.e <- clines$top
    longer.n.e <- clines$longer
  }
  else {
    newline.n.e <- FALSE
    longer.n.e <- labs[["lab.n.e"]]
  }
  ##
  ## Check for "\n" in column n.c
  ##
  clines <- twolines(labs[["lab.n.c"]], "n.c")
  ##
  if (clines$newline) {
    newline.n.c <- TRUE
    labs[["lab.n.c"]] <- clines$bottom
    add.n.c <- clines$top
    longer.n.c <- clines$longer
  }
  else {
    newline.n.c <- FALSE
    longer.n.c <- labs[["lab.n.c"]]
  }
  ##
  ## Check for "\n" in column event.e
  ##
  clines <- twolines(labs[["lab.event.e"]], "event.e")
  ##
  if (clines$newline) {
    newline.event.e <- TRUE
    labs[["lab.event.e"]] <- clines$bottom
    add.event.e <- clines$top
    longer.event.e <- clines$longer
  }
  else {
    newline.event.e <- FALSE
    longer.event.e <- labs[["lab.event.e"]]
  }
  ##
  ## Check for "\n" in column event.c
  ##
  clines <- twolines(labs[["lab.event.c"]], "event.c")
  ##
  if (clines$newline) {
    newline.event.c <- TRUE
    labs[["lab.event.c"]] <- clines$bottom
    add.event.c <- clines$top
    longer.event.c <- clines$longer
  }
  else {
    newline.event.c <- FALSE
    longer.event.c <- labs[["lab.event.c"]]
  }
  ##
  ## Check for "\n" in column mean.e
  ##
  clines <- twolines(labs[["lab.mean.e"]], "mean.e")
  ##
  if (clines$newline) {
    newline.mean.e <- TRUE
    labs[["lab.mean.e"]] <- clines$bottom
    add.mean.e <- clines$top
    longer.mean.e <- clines$longer
  }
  else {
    newline.mean.e <- FALSE
    longer.mean.e <- labs[["lab.mean.e"]]
  }
  ##
  ## Check for "\n" in column mean.c
  ##
  clines <- twolines(labs[["lab.mean.c"]], "mean.c")
  ##
  if (clines$newline) {
    newline.mean.c <- TRUE
    labs[["lab.mean.c"]] <- clines$bottom
    add.mean.c <- clines$top
    longer.mean.c <- clines$longer
  }
  else {
    newline.mean.c <- FALSE
    longer.mean.c <- labs[["lab.mean.c"]]
  }
  ##
  ## Check for "\n" in column sd.e
  ##
  clines <- twolines(labs[["lab.sd.e"]], "sd.e")
  ##
  if (clines$newline) {
    newline.sd.e <- TRUE
    labs[["lab.sd.e"]] <- clines$bottom
    add.sd.e <- clines$top
    longer.sd.e <- clines$longer
  }
  else {
    newline.sd.e <- FALSE
    longer.sd.e <- labs[["lab.sd.e"]]
  }
  ##
  ## Check for "\n" in column sd.c
  ##
  clines <- twolines(labs[["lab.sd.c"]], "sd.c")
  ##
  if (clines$newline) {
    newline.sd.c <- TRUE
    labs[["lab.sd.c"]] <- clines$bottom
    add.sd.c <- clines$top
    longer.sd.c <- clines$longer
  }
  else {
    newline.sd.c <- FALSE
    longer.sd.c <- labs[["lab.sd.c"]]
  }
  ##
  ## Check for "\n" in column cor
  ##
  clines <- twolines(labs[["lab.cor"]], "cor")
  ##
  if (clines$newline) {
    newline.cor <- TRUE
    labs[["lab.cor"]] <- clines$bottom
    add.cor <- clines$top
    longer.cor <- clines$longer
  }
  else {
    newline.cor <- FALSE
    longer.cor <- labs[["lab.cor"]]
  }
  ##
  ## Check for "\n" in column time.e
  ##
  clines <- twolines(labs[["lab.time.e"]], "time.e")
  ##
  if (clines$newline) {
    newline.time.e <- TRUE
    labs[["lab.time.e"]] <- clines$bottom
    add.time.e <- clines$top
    longer.time.e <- clines$longer
  }
  else {
    newline.time.e <- FALSE
    longer.time.e <- labs[["lab.time.e"]]
  }
  ##
  ## Check for "\n" in column time.c
  ##
  clines <- twolines(labs[["lab.time.c"]], "time.c")
  ##
  if (clines$newline) {
    newline.time.c <- TRUE
    labs[["lab.time.c"]] <- clines$bottom
    add.time.c <- clines$top
    longer.time.c <- clines$longer
  }
  else {
    newline.time.c <- FALSE
    longer.time.c <- labs[["lab.time.c"]]
  }
  ##
  ## Check for "\n" in argument smlab
  ##
  clines <- twolines(smlab, arg = TRUE)
  ##
  if (clines$newline) {
    smlab1 <- clines$top
    smlab2 <- clines$bottom
    ##
    newline.smlab <- TRUE
  }
  else {
    smlab1 <- smlab
    smlab2 <- ""
    ##
    newline.smlab <- FALSE
  }
  ##
  ## Check for "\n" in argument label.left
  ##
  clines <- twolines(label.left, arg = TRUE)
  ##
  if (clines$newline) {
    ll1 <- clines$top
    ll2 <- clines$bottom
    ##
    newline.ll <- TRUE
  }
  else {
    ll1 <- label.left
    ll2 <- ""
    ##
    newline.ll <- FALSE
  }
  ##
  ## Check for "\n" in argument label.right
  ##
  clines <- twolines(label.right, arg = TRUE)
  ##
  if (clines$newline) {
    lr1 <- clines$top
    lr2 <- clines$bottom
    ##
    newline.lr <- TRUE
  }
  else {
    lr1 <- label.right
    lr2 <- ""
    ##
    newline.lr <- FALSE
  }
  ##
  ## Check for newlines in additional columns
  ##
  newline.addcol.left  <- FALSE
  newline.addcol.right <- FALSE
  ##
  if (newcols) {
    if (length(leftcols.new) > 0) {
      for (i in seq(along = leftcols.new)) {
        ## Check for "\n" in label of new column
        clines <- twolines(leftlabs.new[i], leftcols.new[i])
        newline.addcol.left <- c(newline.addcol.left, clines$newline)
      }
      newline.addcol.right <- sum(newline.addcol.right) > 0
    }
    if (length(leftcols.new) > 0) {
      for (i in seq(along = leftcols.new)) {
        ## Check for "\n" in label of new column
        clines <- twolines(leftlabs.new[i], leftcols.new[i])
        newline.addcol.left <- c(newline.addcol.left, clines$newline)
      }
      newline.addcol.left <- sum(newline.addcol.left) > 0
    }
  }
  ##
  newline <- newline.studlab | newline.effect | newline.ci | newline.effect.ci |
    newline.w.fixed | newline.w.random | newline.TE | newline.seTE |
    newline.n.e | newline.n.c | newline.event.e | newline.event.c |
    newline.mean.e | newline.mean.c | newline.sd.e | newline.sd.c |
    newline.cor | newline.time.e | newline.time.c |
    newline.smlab | newline.addcol.left | newline.addcol.right
  ##
  newline.all <- newline | (!newline & (newline.ll | newline.lr) & !addrow)
  
  
  ##
  ##
  ## (10) Define columns in forest plot as well as x- and y-limits
  ##
  ##
  if (by) {
    ##
    bylab <- bylabel(bylab, bylevs, print.byvar, byseparator)
    ##
    if (length(text.fixed.w) == 1 & n.by > 1)
      text.fixed.w <- rep(text.fixed.w, n.by)
    if (length(text.random.w) == 1 & n.by > 1)
      text.random.w <- rep(text.random.w, n.by)
    ##
    modlabs <- c(text.fixed, text.random, text.predict,
                 hetstat.overall,
                 text.overall.fixed, text.overall.random,
                 text.subgroup.fixed, text.subgroup.random,
                 bylab, text.fixed.w, text.random.w, hetstat.w,
                 unlist(text.effect.subgroup.fixed),
                 unlist(text.effect.subgroup.random),
                 studlab)
    ##
    TEs    <- c(TE.fixed, TE.random, NA, TE.w, TE)
    lowTEs <- c(lowTE.fixed, lowTE.random, lowTE.predict, lowTE.w, lowTE)
    uppTEs <- c(uppTE.fixed, uppTE.random, uppTE.predict, uppTE.w, uppTE)
    ##
    TEs.exclude    <- c(TE.fixed, TE.random, NA, TE.w, TE.exclude)
    lowTEs.exclude <- c(lowTE.fixed, lowTE.random, lowTE.predict, lowTE.w,
                        lowTE.exclude)
    uppTEs.exclude <- c(uppTE.fixed, uppTE.random, uppTE.predict, uppTE.w,
                        uppTE.exclude)
    ##
    TEs.study <- c("", "", "", rep("", 5 * n.by),
                   format.NA(round(TE.orig, digits), digits, lab.NA))
    seTEs.study <- c("", "", "", rep("", 5 * n.by),
                     format.NA(round(seTE, digits.se), digits.se, lab.NA))
    ##
    if (weight.subgroup == "same") {
      w.fixeds  <- c(NA, NA, NA, rep(NA, length(weight.w.p)), w.fixed.p)
      w.randoms <- c(NA, NA, NA, rep(NA, length(weight.w.p)), w.random.p)
    }
    else {
      w.fixeds  <- c(NA, NA, NA, weight.w.p, w.fixed.p)
      w.randoms <- c(NA, NA, NA, weight.w.p, w.random.p)
    }
    ##
    format.w.fixed  <- format.NA(c(100, weight.w.p, w.fixed.p), digits.weight)
    format.w.random <- format.NA(c(100, weight.w.p, w.random.p), digits.weight)
    w.fixeds.text  <- c(format.w.fixed[1], "--", "", format.w.fixed[-1])
    w.randoms.text <- c("--", format.w.random[1], "", format.w.random[-1])
    ##
    sel.fixed  <- w.fixeds.text == "--"
    sel.random <- w.randoms.text == "--"
    ##
    sel.fixed[sel.by.random] <- TRUE
    sel.random[sel.by.fixed] <- TRUE
    ##
    type.pooled <- c(type.fixed,
                     type.random,
                     "predict",
                     rep(type.subgroup, n.by),
                     rep(type.subgroup, n.by),
                     rep("", 3 * n.by))
    col.diamond.pooled <- c(col.diamond.fixed,
                            col.diamond.random,
                            col.predict,
                            rep(col.diamond.fixed, n.by),
                            rep(col.diamond.random, n.by),
                            rep(NA, 3 * n.by))
    col.diamond.lines.pooled <- c(col.diamond.lines.fixed,
                                  col.diamond.lines.random,
                                  col.predict.lines,
                                  rep(col.diamond.lines.fixed, n.by),
                                  rep(col.diamond.lines.random, n.by),
                                  rep(NA, 3 * n.by))
    col.inside.pooled <- c(col.inside.fixed,
                           col.inside.random,
                           "",
                           rep(col.inside.fixed, n.by),
                           rep(col.inside.random, n.by),
                           rep(NA, 3 * n.by))
  }
  else {
    modlabs <- c(text.fixed, text.random, text.predict,
                 hetstat.overall,
                 text.overall.fixed, text.overall.random,
                 text.subgroup.fixed, text.subgroup.random, studlab)
    ##
    TEs    <- c(TE.fixed, TE.random, NA, TE)
    lowTEs <- c(lowTE.fixed, lowTE.random, lowTE.predict, lowTE)
    uppTEs <- c(uppTE.fixed, uppTE.random, uppTE.predict, uppTE)
    ##
    TEs.exclude    <- c(TE.fixed, TE.random, NA, TE.exclude)
    lowTEs.exclude <- c(lowTE.fixed, lowTE.random, lowTE.predict, lowTE.exclude)
    uppTEs.exclude <- c(uppTE.fixed, uppTE.random, uppTE.predict, uppTE.exclude)
    ##
    TEs.study <- c("", "", "",
                   format.NA(round(TE.orig, digits), digits, lab.NA))
    seTEs.study <- c("", "", "",
                     format.NA(round(seTE, digits.se), digits.se, lab.NA))
    ##
    w.fixeds  <- c(NA, NA, NA, w.fixed.p)
    w.randoms <- c(NA, NA, NA, w.random.p)
    ##
    format.w.fixed  <- format.NA(c(100, w.fixed.p), digits.weight)
    format.w.random <- format.NA(c(100, w.random.p), digits.weight)
    w.fixeds.text  <- c(format.w.fixed[1], "--", "", format.w.fixed[-1])
    w.randoms.text <- c("--", format.w.random[1], "", format.w.random[-1])
    ##
    sel.fixed <- w.fixeds.text == "--"
    sel.random <- w.randoms.text == "--"
    ##
    type.pooled <- c(type.fixed, type.random, "predict")
    col.diamond.pooled <- c(col.diamond.fixed, col.diamond.random, col.predict)
    col.diamond.lines.pooled <- c(col.diamond.lines.fixed, col.diamond.lines.random,
                                  col.predict.lines)
    col.inside.pooled <- c(col.inside.fixed, col.inside.random, "")
  }
  ##
  ## Treatment effect and confidence interval
  ##
  if (backtransf & is.relative.effect(sm)) {
    effect.format <- format.NA(round(exp(TEs), digits), digits, lab.NA.effect)
    ci.format <- ifelse(is.na(lowTEs) | is.na(uppTEs), lab.NA.effect,
                        p.ci(format(round(exp(lowTEs), digits), scientific = FALSE),
                             format(round(exp(uppTEs), digits), scientific = FALSE)))
  }
  else {
    effect.format <- format.NA(round(TEs, digits), digits, lab.NA.effect)
    ci.format <- ifelse(is.na(lowTEs) | is.na(uppTEs), lab.NA.effect,
                        p.ci(format(round(lowTEs, digits), scientific = FALSE),
                             format(round(uppTEs, digits), scientific = FALSE)))
  }
  effect.ci.format <- paste(effect.format, ci.format)
  ##
  ## No treatment effect for prediction interval
  ##
  effect.format[3] <- ""
  ##
  ## Only print prediction interval if requested
  ##
  if (!prediction) {
    ci.format[3] <- ""
    effect.ci.format[3] <- ""
  }
  ##
  ## No treatment effect and confidence interval in heterogeneity line
  ##
  if (by) {
    effect.format[sel.by.noNA] <- ""
    ci.format[sel.by.noNA] <- ""
    effect.ci.format[sel.by.noNA] <- ""
  }
  ##
  ## Weights of fixed and random effects model
  ##
  w.fixed.format  <- paste(w.fixeds.text, if (w.fixed.percent) "%", sep = "")
  w.random.format <- paste(w.randoms.text, if (w.random.percent) "%", sep = "")
  ##
  w.fixed.format[w.fixed.format == "%"] <- ""
  w.random.format[w.random.format == "%"] <- ""
  ##
  w.fixed.format[sel.fixed] <- "--"
  if (by)
    w.fixed.format[sel.by.noNA] <- ""
  w.random.format[sel.random] <- "--"
  if (by)
    w.random.format[sel.by.noNA] <- ""
  ##
  ## Treatment estimate and its standard error
  ##
  TE.format <- TEs.study
  seTE.format <- seTEs.study
  ##
  ## Number of patients, events, and person times
  ##
  sum.n.e <- sum(x$n.e, na.rm = TRUE)
  sum.n.c <- sum(x$n.c, na.rm = TRUE)
  sum.e.e <- sum(x$event.e, na.rm = TRUE)
  sum.e.c <- sum(x$event.c, na.rm = TRUE)
  sum.t.e <- sum(x$time.e, na.rm = TRUE)
  sum.t.c <- sum(x$time.c, na.rm = TRUE)
  ##
  if (by) {
    if (pooled.totals) {
      Ne <- c(sum.n.e, sum.n.e, NA, n.e.w, n.e.w, rep(NA, 3 * n.by), x$n.e)
      Nc <- c(sum.n.c, sum.n.c, NA, n.c.w, n.c.w, rep(NA, 3 * n.by), x$n.c)
    }
    else {
      Ne <- c(NA, NA, NA, rep(NA, 5 * n.by), x$n.e)
      Nc <- c(NA, NA, NA, rep(NA, 5 * n.by), x$n.c)
    }
    if (pooled.events) {
      Ee <- c(sum.e.e, sum.e.e, NA, e.e.w, e.e.w, rep(NA, 3 * n.by), x$event.e)
      Ec <- c(sum.e.c, sum.e.c, NA, e.c.w, e.c.w, rep(NA, 3 * n.by), x$event.c)
    }
    else {
      Ee <- c(NA, NA, NA, rep(NA, 5 * n.by), x$event.e)
      Ec <- c(NA, NA, NA, rep(NA, 5 * n.by), x$event.c)
    }
    if (pooled.times) {
      Te <- c(sum.t.e, sum.t.e, NA, t.e.w, t.e.w, rep(NA, 3 * n.by), x$time.e)
      Tc <- c(sum.t.c, sum.t.c, NA, t.c.w, t.c.w, rep(NA, 3 * n.by), x$time.c)
    }
    else {
      Te <- c(NA, NA, NA, rep(NA, 5 * n.by), x$time.e)
      Tc <- c(NA, NA, NA, rep(NA, 5 * n.by), x$time.c)
    }
  }
  else {
    if (pooled.totals) {
      Ne <- c(sum.n.e, sum.n.e, NA, x$n.e)
      Nc <- c(sum.n.c, sum.n.c, NA, x$n.c)
    }
    else {
      Ne <- c(NA, NA, NA, x$n.e)
      Nc <- c(NA, NA, NA, x$n.c)
    }
    if (pooled.events) {
      Ee <- c(sum.e.e, sum.e.e, NA, x$event.e)
      Ec <- c(sum.e.c, sum.e.c, NA, x$event.c)
    }
    else {
      Ee <- c(NA, NA, NA, x$event.e)
      Ec <- c(NA, NA, NA, x$event.c)
    }
    if (pooled.times) {
      Te <- c(sum.t.e, sum.t.e, NA, x$time.e)
      Tc <- c(sum.t.c, sum.t.c, NA, x$time.c)
    }
    else {
      Te <- c(NA, NA, NA, x$time.e)
      Tc <- c(NA, NA, NA, x$time.c)
    }
  }
  ##
  Ne.format <- ifelse(is.na(Ne), lab.NA, format(Ne, scientific = FALSE))
  Nc.format <- ifelse(is.na(Nc), lab.NA, format(Nc, scientific = FALSE))
  Ee.format <- ifelse(is.na(Ee), lab.NA, format(Ee, scientific = FALSE))
  Ec.format <- ifelse(is.na(Ec), lab.NA, format(Ec, scientific = FALSE))
  ##
  if (is.null(digits.time)) {
    Te.format <- ifelse(is.na(Te), lab.NA, format(Te, scientific = FALSE))
    Tc.format <- ifelse(is.na(Tc), lab.NA, format(Tc, scientific = FALSE))
  }
  else {
    Te.format <- format.NA(round(Te, digits.mean), digits.mean, lab.NA)
    Tc.format <- format.NA(round(Tc, digits.mean), digits.mean, lab.NA)
  }
  ##
  ## Print nothing in line with prediction interval
  ##
  Ne.format[3] <- Nc.format[3] <- Ee.format[3] <- Ec.format[3] <-
    Te.format[3] <- Tc.format[3] <- ""
  ##
  if (by) {
    ##
    ## Print nothing in lines with heterogeneity results for subgroups
    ##
    Ne.format[sel.by.noNA] <- Nc.format[sel.by.noNA] <- ""
    Ee.format[sel.by.noNA] <- Ec.format[sel.by.noNA] <- ""
    Te.format[sel.by.noNA] <- Tc.format[sel.by.noNA] <- ""
  }
  ##
  if (fixed.random) {
    ##
    ## Print nothing in lines with results for random effects model
    ##
    Ne.format[2] <- Nc.format[2] <- ""
    Ee.format[2] <- Ec.format[2] <- ""
    Te.format[2] <- Tc.format[2] <- ""
    ##
    if (by) {
      Ne.format[sel.by.random] <- Nc.format[sel.by.random] <- ""
      Ee.format[sel.by.random] <- Ec.format[sel.by.random] <- ""
      Te.format[sel.by.random] <- Tc.format[sel.by.random] <- ""
    }
  }
  ##
  ## Only print total number of events if pooled.events is TRUE
  ##
  if (!pooled.events) {
    Ee.format[1:2] <- Ec.format[1:2] <- ""
    ##
    if (by) {
      Ee.format[sel.by.fixed]  <- Ec.format[sel.by.fixed]  <- ""
      Ee.format[sel.by.random] <- Ec.format[sel.by.random] <- ""
    }
  }
  ##
  ## Only print total person times if pooled.times is TRUE
  ##
  if (!pooled.times) {
    Te.format[1:2] <- Tc.format[1:2] <- ""
    ##
    if (by) {
      Te.format[sel.by.fixed]  <- Tc.format[sel.by.fixed]  <- ""
      Te.format[sel.by.random] <- Tc.format[sel.by.random] <- ""
    }
  }
  ##
  ## Mean and standard deviation
  ##
  if (by) {
    Me <- c(NA, NA, NA, rep(NA, 5 * n.by), x$mean.e)
    Mc <- c(NA, NA, NA, rep(NA, 5 * n.by), x$mean.c)
    Se <- c(NA, NA, NA, rep(NA, 5 * n.by), x$sd.e)
    Sc <- c(NA, NA, NA, rep(NA, 5 * n.by), x$sd.c)
  }
  else {
    Me <- c(NA, NA, NA, x$mean.e)
    Mc <- c(NA, NA, NA, x$mean.c)
    Se <- c(NA, NA, NA, x$sd.e)
    Sc <- c(NA, NA, NA, x$sd.c)
  }
  ##
  if (is.null(digits.mean)) {
    Me.format <- ifelse(is.na(Me), lab.NA, format(Me, scientific = FALSE))
    Mc.format <- ifelse(is.na(Mc), lab.NA, format(Mc, scientific = FALSE))
  }
  else {
    Me.format <- format.NA(round(Me, digits.mean), digits.mean, lab.NA)
    Mc.format <- format.NA(round(Mc, digits.mean), digits.mean, lab.NA)
  }
  if (is.null(digits.sd)) {
    Se.format <- ifelse(is.na(Se), lab.NA, format(Se, scientific = FALSE))
    Sc.format <- ifelse(is.na(Sc), lab.NA, format(Sc, scientific = FALSE))
  }
  else {
    Se.format <- format.NA(round(Se, digits.sd), digits.sd, lab.NA)
    Sc.format <- format.NA(round(Sc, digits.sd), digits.sd, lab.NA)
  }
  ##
  ## Print nothing for lines with summary results
  ##
  Me.format[1:3] <- Mc.format[1:3] <- Se.format[1:3] <- Sc.format[1:3] <- ""
  ##
  if (by) {
    Me.format[sel.by] <- Mc.format[sel.by] <- ""
    Se.format[sel.by] <- Sc.format[sel.by] <- ""
  }
  ##
  ## Correlation
  ##
  if (by)
    cor <- c(NA, NA, NA, rep(NA, 5 * n.by), x$cor)
  else
    cor <- c(NA, NA, NA, x$cor)
  ##
  if (is.null(digits.cor))
    cor.format <- ifelse(is.na(cor), lab.NA, format(cor, scientific = FALSE))
  else
    cor.format <- format.NA(round(cor, digits.cor), digits.cor, lab.NA)
  ##
  ## Print nothing for lines with summary results
  ##
  cor.format[1:3] <- ""
  ##
  if (by)
    cor.format[sel.by] <- ""
  ##
  ##
  ## y-axis:
  ##
  ##
  if ((!(metaprop | metacor) &
         (any(rightcols %in% c("n.e", "n.c")) |
            any(leftcols  %in% c("n.e", "n.c")))
       ) |
      (metainc &
         (any(rightcols %in% c("time.e", "time.c")) |
            any(leftcols  %in% c("time.e", "time.c")))
       ) |
      (metacont &
         (any(rightcols %in% c("sd.e", "sd.c")) |
            any(leftcols  %in% c("sd.e", "sd.c")))
       ) |
      (!is.null(lab.e.attach.to.col) & !is.null(lab.e)) |
      (!is.null(lab.c.attach.to.col) & !is.null(lab.c)) |
      newline.all
      ) {
    yHead <- 2
    yHeadadd <- 1
  }
  else {
    yHead <- 1
    yHeadadd <- NA
  }
  ##
  if (!by) {
    N <- n.stud
    if (study.results)
      yTE <- 1:N
    else
      yTE <- rep(NA, N)
  }
  else {
    ##
    j <- 1
    k <- 0
    yBylab <- rep(NA, n.by)
    yTE <- rep(NA, n.stud)
    yTE.w.fixed <- yBylab
    yTE.w.random <- yBylab
    yTE.w.hetstat <- yBylab
    yTE.w.effect.fixed <- yBylab
    yTE.w.effect.random <- yBylab
    ##
    for (i in 1:n.by) {
      ##
      k.i <- k.all.w[i]
      k <- k + k.i
      ##
      if (print.subgroup.labels) {
        yBylab[i] <- j
        j <- j + 1
      }
      ##
      if (study.results) {
        yTE[(k - k.i + 1):k] <- j:(j + k.i - 1)
        j <- j + k.i
      }
      else
        yTE[(k - k.i + 1):k] <- NA
      ##
      ## Fixed effect model
      ##
      if (comb.fixed) {
        yTE.w.fixed[i] <- j
        j <- j + 1
      }
      else
        yTE.w.fixed[i] <- NA
      ##
      ## Random effect model
      ##
      if (comb.random) {
        yTE.w.random[i] <- j
        j <- j + 1
      }
      else
        yTE.w.random[i] <- NA
      ##
      ## Only pooled totals
      ##
      if (pooled.totals & !(comb.fixed | comb.random)) {
        yTE.w.fixed[i] <- j
        j <- j + 1
      }
      ##
      ## Heterogeneity statistics
      ##
      if (hetstat) {
        yTE.w.hetstat[i] <- j
        j <- j + 1
      }
      else
        yTE.w.hetstat[i] <- NA
      ##
      ## Test for effect in subgroup (fixed effect)
      ##
      if (test.effect.subgroup.fixed) {
        yTE.w.effect.fixed[i] <- j
        j <- j + 1
      }
      else
        yTE.w.effect.fixed[i] <- NA
      ##
      ## Test for effect in subgroup (random effects)
      ##
      if (test.effect.subgroup.random) {
        yTE.w.effect.random[i] <- j
        j <- j + 1
      }
      else
        yTE.w.effect.random[i] <- NA
      ##
      if (addrow.subgroups)
        j <- j + 1
    }
    ##
    if (!addrow.subgroups)
      j <- j + 1
    ##
    yTE.w <- c(yTE.w.fixed, yTE.w.random, yTE.w.hetstat,
               yTE.w.effect.fixed, yTE.w.effect.random)
  }
  ##
  ##
  ## x-axis:
  ##
  ##
  if (notmiss.xlim && is.numeric(xlim[1]))
    if (is.relative.effect(sm))
      xlim <- log(xlim)
  ##
  if (is.null(xlim)) {
    if (metaprop) {
      xlim <- c(min(c(lowTE, lowTE.predict), na.rm = TRUE),
                max(c(uppTE, uppTE.predict), na.rm = TRUE))
      ##
      if (!is.na(ref) && ref < xlim[1])
        xlim[1] <- ref
      if (!is.na(ref) && ref > xlim[2])
        xlim[2] <- ref
    }
    else {
      sel.low <- is.finite(lowTE)
      sel.upp <- is.finite(uppTE)
      ##
      if (all(!sel.low))
        minTE <- -0.5
      else
        minTE <- min(c(lowTE[sel.low], lowTE.predict), na.rm = TRUE)
      if (all(!sel.upp))
        maxTE <- 0.5
      else
        maxTE <- max(c(uppTE[sel.upp], uppTE.predict), na.rm = TRUE)
      ##
      xlim <- c(minTE, maxTE)
      ##
      if (!is.na(ref) && ref < xlim[1])
        xlim[1] <- ref
      if (!is.na(ref) && ref > xlim[2])
        xlim[2] <- ref
    }
  }
  ##
  symmetric <- FALSE
  ##
  if (!is.null(xlim) && is.character(xlim[1])) {
    ##
    xlim <- setchar(xlim, "symmetric",
                    "should be a numeric vector (min, max) or the character string \"symmetric\"")
    symmetric <- TRUE
    ##
    if (metaprop | metarate) {
      xlim <- c(min(c(lowTE, lowTE.predict), na.rm = TRUE),
                max(c(uppTE, uppTE.predict), na.rm = TRUE))
    }
    else {
      sel.low <- is.finite(lowTE)
      sel.upp <- is.finite(uppTE)
      ##
      if (all(!sel.low))
        minTE <- -0.5
      else
        minTE <- min(c(lowTE[sel.low], lowTE.predict), na.rm = TRUE)
      if (all(!sel.upp))
        maxTE <- 0.5
      else
        maxTE <- max(c(uppTE[sel.upp], uppTE.predict), na.rm = TRUE)
      ##
      if (minTE < 0 & maxTE < 0)
        xlim <- c(minTE, -minTE)
      else if (minTE > 0 & maxTE > 0)
        xlim <- c(-maxTE, maxTE)
      else
        xlim <- c(-max(abs(c(minTE, maxTE))), max(abs(c(minTE, maxTE))))
    }
  }
  ##
  if (!is.na(ref) &&
      round(xlim[2] - ref, 6) == round(ref - xlim[1], 6))
    symmetric <- TRUE
  ##
  if (by) {
    if (all(is.na(c(yTE, yTE.w))))
      max.yTE <- 0
    else
      max.yTE <- max(c(yTE, yTE.w), na.rm = TRUE)
  }
  else {
    if (all(is.na(yTE)))
      max.yTE <- 0
    else
    max.yTE <- max(yTE, na.rm = TRUE)
  }
  ##
  yNext <- max.yTE + ifelse(max.yTE == 0 | !addrow.overall, 1, 2)
  ##
  if (missing(xlab.pos))
    xlab.pos <- mean(xlim)
  ##
  if (missing(smlab.pos))
    smlab.pos <- mean(xlim)
  ##
  yTE.fixed  <- NA
  yTE.random <- NA
  yPredict   <- NA
  yHetstat <- NA
  yOverall.fixed  <- NA
  yOverall.random <- NA
  ySubgroup.fixed  <- NA
  ySubgroup.random <- NA
  ##
  if (fixed.random & overall) {
    yTE.fixed  <- yNext
    yTE.random <- yNext + 1
    yNext      <- yNext + 2
  }
  ##
  else if (comb.fixed & !comb.random & overall) {
    yTE.fixed <- yNext
    yNext     <- yNext + 1
  }
  ##
  else if (!comb.fixed & comb.random & overall) {
    yTE.random <- yNext
    yNext      <- yNext + 1
  }
  ##
  else if (!comb.fixed & !comb.random & pooled.totals & overall) {
    yTE.fixed  <- yNext
    yNext      <- yNext + 1
    if (missing(text.fixed))
      text.fixed <- "Overall"
  }
  ##
  if (prediction) {
    yPredict <- yNext
    yNext    <- yNext + 1
  }
  ##
  if (overall.hetstat) {
    yHetstat <- yNext
    yNext    <- yNext + 1
  }
  ##
  if (test.overall.fixed) {
    yOverall.fixed <- yNext
    yNext          <- yNext + 1
  }
  ##
  if (test.overall.random) {
    yOverall.random <- yNext
    yNext           <- yNext + 1
  }
  ##
  if (test.subgroup.fixed) {
    ySubgroup.fixed <- yNext
    yNext           <- yNext + 1
  }
  ##
  if (test.subgroup.random) {
    ySubgroup.random <- yNext
    yNext            <- yNext + 1
  }
  ##
  if (!comb.fixed & !pooled.totals) text.fixed <- ""
  if (!comb.random) text.random <- ""
  if (!prediction) text.predict <- ""
  ##
  yTE <- yHead + yTE + addrow
  ##
  yTE.fixed  <- yHead + yTE.fixed + addrow
  yTE.random <- yHead + yTE.random + addrow
  yPredict   <- yHead + yPredict + addrow
  ##
  yHetstat <- yHead + yHetstat + addrow
  yOverall.fixed  <- yHead + yOverall.fixed + addrow
  yOverall.random <- yHead + yOverall.random + addrow
  ySubgroup.fixed  <- yHead + ySubgroup.fixed + addrow
  ySubgroup.random <- yHead + ySubgroup.random + addrow
  ##
  yStats <- c(yHetstat,
              yOverall.fixed, yOverall.random,
              ySubgroup.fixed, ySubgroup.random)
  ##
  if (by) {
    yBylab <- yHead + yBylab + addrow
    yTE.w  <- yHead + yTE.w + addrow
    ##
    yLab <- c(yHead,
              yTE.fixed, yTE.random, yPredict,
              yStats,
              yBylab, yTE.w,
              yTE)
    ##
    yS <- c(yHead, yTE.fixed, yTE.random, yPredict, yTE.w, yTE)
  }
  else {
    yLab <- c(yHead, yTE.fixed, yTE.random, yPredict,
              yStats,
              yTE)
    ##
    yS <- c(yHead, yTE.fixed, yTE.random, yPredict, yTE)
  }
  
  
  ##
  ##
  ## (11) Format columns in forest plot
  ##
  ##
  col.studlab <- list(labels = 
                        lapply(as.list(c(labs[["lab.studlab"]], modlabs)),
                               tg,
                               xpos = xpos.s, just = just.s,
                               fs = fs.study.labels,
                               ff = ff.study.labels),
                      rows = yLab
                      )
  ## Study label:
  col.studlab$labels[[1]] <- tg(labs[["lab.studlab"]], xpos.s,
                                just.s, fs.head, ff.head)
  ## Fixed effect estimate:
  col.studlab$labels[[2]] <- tg(text.fixed, xpos.s, just.s,
                                fs.fixed.labels, ff.fixed.labels)
  ## Random effects estimate:
  col.studlab$labels[[3]] <- tg(text.random, xpos.s, just.s,
                                fs.random.labels, ff.random.labels)
  ## Prediction interval:
  col.studlab$labels[[4]] <- tg(text.predict,xpos.s, just.s,
                                fs.predict.labels, ff.predict.labels)
  ## Heterogeneity statistics:
  col.studlab$labels[[5]] <- tg(hetstat.overall, xpos.s, just.s,
                                fs.hetstat, ff.hetstat)
  ## Test for overall effect (fixed effect model):
  col.studlab$labels[[6]] <- tg(text.overall.fixed, xpos.s, just.s,
                                fs.test.overall, ff.test.overall)
  ## Test for overall effect (random effects model):
  col.studlab$labels[[7]] <- tg(text.overall.random, xpos.s, just.s,
                                fs.test.overall, ff.test.overall)
  ## Test for subgroup differences (fixed effect model):
  col.studlab$labels[[8]] <- tg(text.subgroup.fixed, xpos.s, just.s,
                                fs.test.subgroup, ff.test.subgroup)
  ## Test for subgroup differences (random effects model):
  col.studlab$labels[[9]] <- tg(text.subgroup.random, xpos.s,
                                just.s,fs.test.subgroup, ff.test.subgroup)
  ##
  n.summaries <- 9
  ##
  if (by) {
    for (i in 1:n.by) {
      ## Subgroup labels:
      col.studlab$labels[[n.summaries + i]] <-
        tg(bylab[i], xpos.s, just.s,
           fs.head, ff.head, col.by)
      ## Fixed effect estimates:
      col.studlab$labels[[n.summaries + n.by + i]] <-
        tg(text.fixed.w[i], xpos.s, just.s,
           fs.fixed.labels, ff.fixed.labels, col.by)
      ## Random effects estimates:
      col.studlab$labels[[n.summaries + 2 * n.by + i]] <-
        tg(text.random.w[i], xpos.s, just.s,
           fs.random.labels, ff.random.labels, col.by)
      ## Heterogeneity statistics:
      col.studlab$labels[[n.summaries + 3 * n.by + i]] <-
        tg(hetstat.w[[i]], xpos.s, just.s,
           fs.hetstat, ff.hetstat, col.by)
      ## Test for effect in subgroup (fixed effect model):
      col.studlab$labels[[n.summaries + 4 * n.by + i]] <-
        tg(text.effect.subgroup.fixed[[i]], xpos.s, just.s,
           fs.test.effect.subgroup, ff.test.effect.subgroup, col.by)
      ## Test for effect in subgroup (random effects model):
      col.studlab$labels[[n.summaries + 5 * n.by + i]] <-
        tg(text.effect.subgroup.random[[i]], xpos.s, just.s,
           fs.test.effect.subgroup, ff.test.effect.subgroup, col.by)
    }
  }
  ##
  fcs <- list(fs.study = fs.study, ff.study = ff.study,
              fs.heading = fs.head, ff.heading = ff.head,
              fs.fixed = fs.fixed, ff.fixed = ff.fixed,
              fs.random = fs.random, ff.random = ff.random,
              fs.predict = fs.predict, ff.predict = ff.predict,
              by = by, n.by = n.by, col.by = col.by)
  ##
  col.effect <- formatcol(labs[["lab.effect"]], effect.format, yS, just.c, fcs)
  ##
  col.ci <- formatcol(labs[["lab.ci"]], ci.format, yS, just.c, fcs)
  ##
  col.effect.ci <- formatcol(labs[["lab.effect.ci"]], effect.ci.format, yS,
                             if (revman5) "center" else just.c, fcs)
  ##
  col.w.fixed  <- formatcol(labs[["lab.w.fixed"]], w.fixed.format, yS,
                            just.c, fcs)
  col.w.random <- formatcol(labs[["lab.w.random"]], w.random.format, yS,
                            just.c, fcs)
  ##
  col.TE <- formatcol(labs[["lab.TE"]], TE.format, yS, just.c, fcs)
  col.seTE <- formatcol(labs[["lab.seTE"]], seTE.format, yS, just.c, fcs)
  ##
  col.n.e <- formatcol(labs[["lab.n.e"]], Ne.format, yS, just.c, fcs)
  col.n.c <- formatcol(labs[["lab.n.c"]], Nc.format, yS, just.c, fcs)
  ##
  col.event.e <- formatcol(labs[["lab.event.e"]], Ee.format, yS, just.c, fcs)
  col.event.c <- formatcol(labs[["lab.event.c"]], Ec.format, yS, just.c, fcs)
  ##
  col.mean.e <- formatcol(labs[["lab.mean.e"]], Me.format, yS, just.c, fcs)
  col.mean.c <- formatcol(labs[["lab.mean.c"]], Mc.format, yS, just.c, fcs)
  ##
  col.sd.e <- formatcol(labs[["lab.sd.e"]], Se.format, yS, just.c, fcs)
  col.sd.c <- formatcol(labs[["lab.sd.c"]], Sc.format, yS, just.c, fcs)
  ##
  col.cor <- formatcol(labs[["lab.cor"]], cor.format, yS, just.c, fcs)
  ##
  col.time.e <- formatcol(labs[["lab.time.e"]], Te.format, yS, just.c, fcs)
  col.time.c <- formatcol(labs[["lab.time.c"]], Tc.format, yS, just.c, fcs)
  ##
  ##
  ##
  col.effect.calc <- formatcol(longer.effect, effect.format, yS, just.c, fcs)
  ##
  col.ci.calc <- formatcol(longer.ci, ci.format, yS, just.c, fcs)
  ##
  col.effect.ci.calc <- formatcol(longer.effect.ci, effect.ci.format, yS,
                                  just.c, fcs)
  ##
  col.w.fixed.calc  <- formatcol(longer.w.fixed, w.fixed.format, yS,
                                 just.c, fcs)
  col.w.random.calc <- formatcol(longer.w.random, w.random.format, yS,
                                 just.c, fcs)
  ##
  col.TE.calc <- formatcol(longer.TE, TE.format, yS, just.c, fcs)
  col.seTE.calc <- formatcol(longer.seTE, seTE.format, yS, just.c, fcs)
  ##
  col.n.e.calc <- formatcol(longer.n.e, Ne.format, yS, just.c, fcs)
  col.n.c.calc <- formatcol(longer.n.c, Nc.format, yS, just.c, fcs)
  ##
  col.event.e.calc <- formatcol(longer.event.e, Ee.format, yS,
                                just.c, fcs)
  col.event.c.calc <- formatcol(longer.event.c, Ec.format, yS, just.c, fcs)
  ##
  col.mean.e.calc <- formatcol(longer.mean.e, Me.format, yS, just.c, fcs)
  col.mean.c.calc <- formatcol(longer.mean.c, Mc.format, yS, just.c, fcs)
  ##
  col.sd.e.calc <- formatcol(longer.sd.e, Se.format, yS, just.c, fcs)
  col.sd.c.calc <- formatcol(longer.sd.c, Sc.format, yS, just.c, fcs)
  ##
  col.cor.calc <- formatcol(longer.cor, cor.format, yS, just.c, fcs)
  ##
  col.time.e.calc <- formatcol(longer.time.e, Te.format, yS, just.c, fcs)
  col.time.c.calc <- formatcol(longer.time.c, Tc.format, yS, just.c, fcs)
  ##
  ##
  ##
  if (length(type.study) == 1)
    type.study <- rep(type.study, length(TE))
  else if (length(type.study) != length(TE))
    stop("Argument 'type.study' must be a single character or of same length as number of studies.")
  ##
  col.forest <- list(eff = TEs.exclude,
                     low = lowTEs.exclude,
                     upp = uppTEs.exclude,
                     rows = yS[-1],
                     ##
                     ## "square" means normal confidence interval, "diamond" means meta-analysis diamond,
                     ## "predict" means prediction interval
                     ##
                     type = c(type.pooled, type.study),
                     ##
                     col = c(col.diamond.lines.pooled, col.study),
                     col.square = c(col.diamond.pooled, col.square),
                     col.square.lines = c(col.diamond.lines.pooled, col.square.lines),
                     col.inside = c(col.inside.pooled, col.inside),
                     ##
                     col.diamond = c(col.diamond.pooled, col.square),
                     col.diamond.lines = c(col.diamond.lines.pooled, col.square.lines),
                     ##
                     lwd = lwd
                     )
  ##
  ## Sizes of squares
  ##
  if (weight.study == "same") {
    information <- rep(0.9, length(TEs))
  }
  else {
    ##
    if (weight.study == "fixed")
      information <- sqrt(w.fixeds)
    else if (weight.study == "random")
      information <- sqrt(w.randoms)
    ## Square height equal to 1 for most precise study result
    if (!all(is.na(information)))
      information <- information / max(information, na.rm = TRUE)
    else
      information <- rep(0.9, length(TEs))
    ## Same / maximum polygon height for all meta-analytical results
    ## (both overall and subgroup results)
    information[is.na(information)] <- 1
  }
  ##
  col.forest$sizes <- information
  col.forest$sizes <- col.forest$sizes * squaresize
  ##
  ## Width of column 3
  col.forestwidth <- plotwidth
  ##
  ## Range on the x-axis for column 3
  col.forest$range <- xlim
  ##
  cols <- list(col.studlab = col.studlab,
               col.effect = col.effect,
               col.ci = col.ci,
               col.effect.ci = col.effect.ci,
               col.w.fixed = col.w.fixed,
               col.w.random = col.w.random,
               col.TE = col.TE,
               col.seTE = col.seTE)
  ##
  cols.calc <- list(col.studlab = col.studlab,
                    col.effect = col.effect.calc,
                    col.ci = col.ci.calc,
                    col.effect.ci = col.effect.ci.calc,
                    col.w.fixed = col.w.fixed.calc,
                    col.w.random = col.w.random.calc,
                    col.TE = col.TE.calc,
                    col.seTE = col.seTE.calc)
  ##
  cols[["col.n.e"]] <- col.n.e
  cols[["col.n.c"]] <- col.n.c
  cols[["col.event.e"]] <- col.event.e
  cols[["col.event.c"]] <- col.event.c
  ##
  cols[["col.mean.e"]] <- col.mean.e
  cols[["col.mean.c"]] <- col.mean.c
  cols[["col.sd.e"]] <- col.sd.e
  cols[["col.sd.c"]] <- col.sd.c
  ##
  cols[["col.cor"]] <- col.cor
  ##
  cols[["col.time.e"]] <- col.time.e
  cols[["col.time.c"]] <- col.time.c
  ##
  cols.calc[["col.n.e"]] <- col.n.e.calc
  cols.calc[["col.n.c"]] <- col.n.c.calc
  cols.calc[["col.event.e"]] <- col.event.e.calc
  cols.calc[["col.event.c"]] <- col.event.c.calc
  ##
  cols.calc[["col.mean.e"]] <- col.mean.e.calc
  cols.calc[["col.mean.c"]] <- col.mean.c.calc
  cols.calc[["col.sd.e"]] <- col.sd.e.calc
  cols.calc[["col.sd.c"]] <- col.sd.c.calc
  ##
  cols.calc[["col.cor"]] <- col.cor.calc
  ##
  cols.calc[["col.time.e"]] <- col.time.e.calc
  cols.calc[["col.time.c"]] <- col.time.c.calc
  ##
  if (newcols) {
    ##
    if (length(leftcols.new) > 0)
      if (length(just.addcols.left) != 1) {
        if (length(just.addcols.left) != length(leftcols.new))
            stop("Length of argument 'just.addcols.left' must be one or same as number of additional columms in argument 'leftcols'.")
      }
      else
        just.addcols.left <- rep(just.addcols.left, length(leftcols.new))
    ##
    if (length(rightcols.new) > 0)
      if (length(just.addcols.right) != 1) {
        if (length(just.addcols.right) != length(rightcols.new))
            stop("Length of argument 'just.addcols.right' must be one or same as number of additional columms in argument 'rightcols'.")
      }
      else
        just.addcols.right <- rep(just.addcols.right, length(rightcols.new))
    ##
    if (by) {
      for (i in seq(along = rightcols.new)) {
        tname <- paste("col.", rightcols.new[i], sep = "")
        if (length(dataset1[[rightcols.new[i]]]) != 0)
          tmp.r <- dataset1[[rightcols.new[i]]]
        else if (length(dataset2[[rightcols.new[i]]]) != 0)
          tmp.r <- dataset2[[rightcols.new[i]]]
        if (is.factor(tmp.r))
          tmp.r <- as.character(tmp.r)
        tmp.r <- ifelse(is.na(tmp.r), lab.NA, tmp.r)
        ##
        ## Check for "\n" in label of new column
        ##
        clines <- twolines(rightlabs.new[i], rightcols.new[i])
        ##
        if (clines$newline) {
          lab.new <- clines$bottom
          longer.new <- clines$longer
        }
        else
          lab.new <- longer.new <- rightlabs.new[i]
        cols[[tname]] <- formatcol(lab.new,
                                   c("", "", "", rep("", length(TE.w)), tmp.r),
                                   yS,
                                   just.addcols.right[i],
                                   fcs)
        cols.calc[[tname]] <- formatcol(longer.new,
                                        c("", "", "", rep("", length(TE.w)), tmp.r),
                                        yS,
                                        just.addcols.right[i],
                                        fcs)
      }
      for (i in seq(along = leftcols.new)) {
        tname <- paste("col.", leftcols.new[i], sep = "")
        if (length(dataset1[[leftcols.new[i]]]) != 0)
          tmp.l <- dataset1[[leftcols.new[i]]]        
        else if (length(dataset2[[leftcols.new[i]]]) != 0)
          tmp.l <- dataset2[[leftcols.new[i]]]
        if (is.factor(tmp.l))
          tmp.l <- as.character(tmp.l)
        tmp.l <- ifelse(is.na(tmp.l), lab.NA, tmp.l)
        ##
        ## Check for "\n" in label of new column
        ##
        clines <- twolines(leftlabs.new[i], leftcols.new[i])
        ##
        if (clines$newline) {
          lab.new <- clines$bottom
          longer.new <- clines$longer
        }
        else
          lab.new <- longer.new <- leftlabs.new[i]
        ##
        cols[[tname]] <- formatcol(lab.new,
                                   c("", "", "",
                                     rep("", length(TE.w)), tmp.l),
                                   yS,
                                   just.addcols.left[i],
                                   fcs)
        ##
        cols.calc[[tname]] <- formatcol(longer.new,
                                        c("", "", "",
                                          rep("", length(TE.w)), tmp.l),
                                        yS,
                                        just.addcols.left[i],
                                        fcs)
      }
    }
    else {
      for (i in seq(along = rightcols.new)) {
        tname <- paste("col.", rightcols.new[i], sep = "")
        if (length(dataset1[[rightcols.new[i]]]) != 0)
          tmp.r <- dataset1[[rightcols.new[i]]]
        else if (length(dataset2[[rightcols.new[i]]]) != 0)
          tmp.r <- dataset2[[rightcols.new[i]]]
        if (is.factor(tmp.r))
          tmp.r <- as.character(tmp.r)
        tmp.r <- ifelse(is.na(tmp.r), "", tmp.r)
        ##
        ## Check for "\n" in label of new column
        ##
        clines <- twolines(rightlabs.new[i], rightcols.new[i])
        ##
        if (clines$newline) {
          lab.new <- clines$bottom
          longer.new <- clines$longer
        }
        else
          lab.new <- longer.new <- rightlabs.new[i]
        ##
        cols[[tname]] <- formatcol(lab.new,
                                   c("", "", "", tmp.r),
                                   yS,
                                   just.addcols.right[i],
                                   fcs)
        cols.calc[[tname]] <- formatcol(longer.new,
                                        c("", "", "", tmp.r),
                                        yS,
                                        just.addcols.right[i],
                                        fcs)
      }
      for (i in seq(along = leftcols.new)) {
        tname <- paste("col.", leftcols.new[i], sep = "")
        if (length(dataset1[[leftcols.new[i]]]) != 0)
          tmp.l <- dataset1[[leftcols.new[i]]]        
        else if (length(dataset2[[leftcols.new[i]]]) != 0)
          tmp.l <- dataset2[[leftcols.new[i]]]
        if (is.factor(tmp.l))
          tmp.l <- as.character(tmp.l)
        tmp.l <- ifelse(is.na(tmp.l), "", tmp.l)
        ##
        ## Check for "\n" in label of new column
        ##
        clines <- twolines(leftlabs.new[i], leftcols.new[i])
        ##
        if (clines$newline) {
          lab.new <- clines$bottom
          longer.new <- clines$longer
        }
        else
          lab.new <- longer.new <- leftlabs.new[i]
        ##
        cols[[tname]] <- formatcol(lab.new,
                                   c("", "", "", tmp.l),
                                   yS,
                                   just.addcols.left[i],
                                   fcs)
        cols.calc[[tname]] <- formatcol(longer.new,
                                        c("", "", "", tmp.l),
                                        yS,
                                        just.addcols.left[i],
                                        fcs)
      }
    }
  }
  ##
  col.lab.e <- tgl(lab.e, xpos.c, just.c, fs.head, ff.head)
  ##
  col.lab.c <- tgl(lab.c, xpos.c, just.c, fs.head, ff.head)
  ##
  ##
  ##
  if (newline.studlab)
    col.add.studlab <- tgl(add.studlab, xpos.s, just.s, fs.head, ff.head)
  ##
  if (newline.effect)
    col.add.effect <- tgl(add.effect, xpos.c, just.c, fs.head, ff.head)
  ##
  if (newline.ci)
    col.add.ci <- tgl(add.ci, xpos.c, just.c, fs.head, ff.head)
  ##
  if (newline.effect.ci)
    col.add.effect.ci <- tgl(add.effect.ci,
                             if (revman5) 0.5 else xpos.c,
                             if (revman5) "center" else just.c,
                             fs.head, ff.head)
  ##
  if (newline.w.fixed)
    col.add.w.fixed <- tgl(add.w.fixed, xpos.c, just.c, fs.head, ff.head)
  ##
  if (newline.w.random)
    col.add.w.random <- tgl(add.w.random, xpos.c, just.c, fs.head, ff.head)
  ##
  if (newline.TE)
    col.add.TE <- tgl(add.TE, xpos.c, just.c, fs.head, ff.head)
  ##
  if (newline.seTE)
    col.add.seTE <- tgl(add.seTE, xpos.c, just.c, fs.head, ff.head)
  ##
  if (newline.n.e)
    col.add.n.e <- tgl(add.n.e, xpos.c, just.c, fs.head, ff.head)
  ##
  if (newline.n.c)
    col.add.n.c <- tgl(add.n.c, xpos.c, just.c, fs.head, ff.head)
  ##
  if (newline.event.e)
    col.add.event.e <- tgl(add.event.e, xpos.c, just.c, fs.head, ff.head)
  ##
  if (newline.event.c)
    col.add.event.c <- tgl(add.event.c, xpos.c, just.c, fs.head, ff.head)
  ##
  if (newline.mean.e)
    col.add.mean.e <- tgl(add.mean.e, xpos.c, just.c, fs.head, ff.head)
  ##
  if (newline.mean.c)
    col.add.mean.c <- tgl(add.mean.c, xpos.c, just.c, fs.head, ff.head)
  ##
  if (newline.sd.e)
    col.add.sd.e <- tgl(add.sd.e, xpos.c, just.c, fs.head, ff.head)
  ##
  if (newline.sd.c)
    col.add.sd.c <- tgl(add.sd.c, xpos.c, just.c, fs.head, ff.head)
  ##
  if (newline.cor)
    col.add.cor <- tgl(add.cor, xpos.c, just.c, fs.head, ff.head)
  ##
  if (newline.time.e)
    col.add.time.e <- tgl(add.time.e, xpos.c, just.c, fs.head, ff.head)
  ##
  if (newline.time.c)
    col.add.time.c <- tgl(add.time.c, xpos.c, just.c, fs.head, ff.head)
  ##
  leftcols  <- paste("col.", leftcols, sep = "")
  rightcols <- paste("col.", rightcols, sep = "")
  
  
  ##
  ##
  ## (12) Calculate width of columns in forest plot
  ##
  ##
  ## Exclude lines with summary measures from calculation of column
  ## width for study labels
  ##
  if (by) {
    del.lines <- c(if (!calcwidth.fixed)              # FE
                     2,
                   if (!calcwidth.random)             # RE
                     3,
                   if (!calcwidth.predict)            # PI
                     4,
                   if (!calcwidth.tests)              # tests
                     5:n.summaries,
                   ##
                   n.summaries + 0 * n.by + 1:n.by,   # subgroup labels
                   if (!calcwidth.fixed)              # FE in subgroups
                     n.summaries + 1 * n.by + 1:n.by,
                   if (!calcwidth.random)             # RE in subgroups
                     n.summaries + 2 * n.by + 1:n.by,
                   if (!calcwidth.hetstat)            # heterogeneity statistic in subgroups
                     n.summaries + 3 * n.by + 1:n.by,
                   if (!calcwidth.tests)              # test for effect in subgroup (fixed effect)
                     n.summaries + 4 * n.by + 1:n.by,
                   if (!calcwidth.tests)              # test for effect in subgroup (random effects)
                     n.summaries + 5 * n.by + 1:n.by
                   )
  }
  else
    del.lines <- c(if (!calcwidth.fixed)   # FE
                     2,
                   if (!calcwidth.random)  # RE
                     3,
                   if (!calcwidth.predict) # PI
                     4,
                   if (!calcwidth.tests)   # tests
                     5:n.summaries)
  ##
  for (i in seq(along = leftcols)) {
    if (i == 1) {
      if (leftcols[[i]] == "col.studlab" & !is.null(del.lines))
        x1 <- unit.c(wcalc(cols.calc[[leftcols[i]]]$labels[-del.lines]))
       else
        x1 <- unit.c(wcalc(cols.calc[[leftcols[i]]]$labels))
    }
    else {
      if (leftcols[[i]] == "col.studlab" & !is.null(del.lines))
        x1 <- unit.c(x1,
                     colgap.left,
                     wcalc(cols.calc[[leftcols[i]]]$labels[-del.lines]))
      else
        x1 <- unit.c(x1,
                     if (leftcols[[i - 1]] == "col.studlab")
                       colgap.studlab
                     else colgap.left,
                     wcalc(cols.calc[[leftcols[i]]]$labels))
    }
  }
  ##
  x1 <- unit.c(x1, colgap.forest.left, col.forestwidth)
  ##
  if (rsel) {
    for (i in seq(along = rightcols)) {
      x1 <- unit.c(x1,
                   if (i == 1) colgap.forest.right else colgap.right,
                   wcalc(cols.calc[[rightcols[i]]]$labels))
    }
  }
  
  
  ##
  ##
  ## (13) Process arguments smlab, label.left and label.right
  ##
  ##
  if (by) {
    addline <- addrow * (!any(c(test.overall.fixed, test.overall.random,
                                overall.hetstat,
                                test.subgroup.fixed, test.subgroup.random)))
    ##
    nrow <- max(addline + c(yTE, yTE.fixed, yTE.random, yPredict,
                            yStats, yTE.w), na.rm = TRUE)
  }
  else {
    addline <- addrow * (!any(c(test.overall.fixed, test.overall.random,
                                overall.hetstat)))
    ##
    nrow <- max(addline + c(yTE, yTE.fixed, yTE.random, yPredict,
                            yStats), na.rm = TRUE)
  }
  ##
  summary.lines <- hetstat + test.overall.fixed + test.overall.random +
    test.subgroup.fixed + test.subgroup.random
  ymin.line <- summary.lines
  ymax.line <- nrow - ifelse(is.na(yHeadadd), 1, 2)
  ##
  ## Summary label at top of forest plot
  ##
  smlab1 <- tgl(smlab1, unit(smlab.pos, "native"), "center", fs.smlab, ff.smlab,
                rows = 1 + (!is.na(yHeadadd) & !newline.smlab))
  ##
  if (newline.smlab)
    smlab2 <- tgl(smlab2, unit(smlab.pos, "native"), "center", fs.smlab, ff.smlab,
                  rows = 2)
  ##
  ## Left and right label on x-axis:
  ##
  if (!bottom.lr & !is.na(ref)) {
    row1.lr <- if (!newline & (newline.ll | newline.lr) & !addrow)
                 1
               else if (!is.na(yHeadadd) & addrow)
                 2
               else if (is.na(yHeadadd))
                 1
               else
                 2
    ##
    ll1 <- tgl(ll1, unit(ref - (xlim[2] - xlim[1]) / 30, "native"),
               "right", fs.lr, ff.lr, col.label.left,
               rows = row1.lr)
    ##
    if (newline.ll)
      ll2 <- tgl(ll2, unit(ref - (xlim[2] - xlim[1]) / 30, "native"),
                 "right", fs.lr, ff.lr, col.label.left,
                 rows = row1.lr + 1)
    ##
    lr1 <- tgl(lr1, unit(ref + (xlim[2] - xlim[1]) / 30, "native"),
               "left", fs.lr, ff.lr, col.label.right,
               rows = row1.lr)
    ##
    if (newline.lr)
      lr2 <- tgl(lr2, unit(ref + (xlim[2] - xlim[1]) / 30, "native"),
                 "left", fs.lr, ff.lr, col.label.right,
                 rows = row1.lr + 1)
  }
  
  
  ##
  ##
  ## (14) Generate forest plot
  ##
  ##
  if (new)
    grid.newpage()
  ##
  pushViewport(viewport(layout = grid.layout(
                          nrow,
                          length(x1),
                          widths = x1,
                          heights = unit(spacing, "lines"))))
  ##
  ## Left side of forest plot
  ##
  j <- 1
  ##
  for (i in seq(along = leftcols)) {
    add.text(cols[[leftcols[i]]], j)
    ##
    if (!is.na(yHeadadd)) {
      if (!is.null(lab.e.attach.to.col)) {
        if (leftcols[i] == paste("col.", lab.e.attach.to.col, sep = ""))
          add.text(col.lab.e, j)
      }
      else if (metabin) {
        if (leftcols[i] == "col.n.e" & just.c == "right")
          add.text(col.lab.e, j)
        else if (leftcols[i] == "col.event.e" & just.c %in% c("left", "center"))
          add.text(col.lab.e, j)
      }
      else if (metacont) {
        if (leftcols[i] == "col.sd.e" & just.c == "right")
          add.text(col.lab.e, j)
        else if (leftcols[i] == "col.mean.e" & just.c %in% c("left", "center"))
          add.text(col.lab.e, j)
      }
      else if (metainc) {
        if (leftcols[i] == "col.time.e" & just.c == "right")
          add.text(col.lab.e, j)
        else if (leftcols[i] == "col.event.e" & just.c %in% c("left", "center"))
          add.text(col.lab.e, j)
      }
      ##
      if (!is.null(lab.c.attach.to.col)) {
        if (leftcols[i] == paste("col.", lab.c.attach.to.col, sep = ""))
          add.text(col.lab.c, j)
      }
      else if (metabin) {
        if (leftcols[i] == "col.n.c" & just.c == "right")
          add.text(col.lab.c, j)
        else if (leftcols[i] == "col.event.c" & just.c %in% c("left", "center"))
          add.text(col.lab.c, j)
      }
      else if (metacont) {
        if (leftcols[i] == "col.sd.c" & just.c == "right")
          add.text(col.lab.c, j)
        else if (leftcols[i] == "col.mean.c" & just.c %in% c("left", "center"))
          add.text(col.lab.c, j)
      }
      else if (metainc) {
        if (leftcols[i] == "col.time.c" & just.c == "right")
          add.text(col.lab.c, j)
        else if (leftcols[i] == "col.event.c" & just.c %in% c("left", "center"))
          add.text(col.lab.c, j)
      }
      ##
      if (newline.studlab & leftcols[i] == "col.studlab")
        add.text(col.add.studlab, j)
      if (newline.effect & leftcols[i] == "col.effect")
        add.text(col.add.effect, j)
      if (newline.ci & leftcols[i] == "col.ci")
        add.text(col.add.ci, j)
      if (newline.effect.ci & leftcols[i] == "col.effect.ci")
        add.text(col.add.effect.ci, j)
      if (newline.w.fixed & leftcols[i] == "col.w.fixed")
        add.text(col.add.w.fixed, j)
      if (newline.w.random & leftcols[i] == "col.w.random")
        add.text(col.add.w.random, j)
      if (newline.TE & leftcols[i] == "col.TE")
        add.text(col.add.TE, j)
      if (newline.seTE & leftcols[i] == "col.seTE")
        add.text(col.add.seTE, j)
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
      ##
      ## Add text in first line of forest plot for new columns
      ##
      if (newcols)
        if (length(leftcols.new) > 0 &
            leftcols[i] %in% paste("col.", leftcols.new, sep = "")) {
          sel <- paste("col.", leftcols.new, sep = "") == leftcols[i]
          ##
          ## Check for "\n" in label of new column
          ##
          clines <- twolines(leftlabs.new[sel], leftcols[i])
          ##
          just.new <- just.addcols.left[sel]
          ##
          if (just.new == "left")
            xpos.new <- 0
          else if (just.new == "center")
            xpos.new <- 0.5
          else if (just.new == "right")
            xpos.new <- 1
          ##
          ## Add first line
          ##
          if (clines$newline)
            add.text(tgl(clines$top, xpos.new, just.new, fs.head, ff.head), j)
        }
    }
    ##
    j <- j + 2
  }
  ##
  ## Produce forest plot
  ##
  draw.lines(col.forest, j,
             ref, TE.fixed, TE.random,
             overall, comb.fixed, comb.random, prediction,
             lwd, lty.fixed, lty.random,
             spacing * ymin.line, spacing * ymax.line,
             spacing * (ymin.line + 0.5), spacing * ymax.line,
             addrow, print.label, bottom.lr)
  ##
  draw.axis(col.forest, j, yS, log.xaxis, at, label,
            fs.axis, ff.axis, lwd,
            xlim, notmiss.xlim)
  ##
  if (bottom.lr) {
    add.text(smlab1, j, xscale = col.forest$range)
    ##
    if (newline.smlab)
      add.text(smlab2, j, xscale = col.forest$range)
  }
  ##
  if (print.label) {
    if (!bottom.lr) {
      add.text(ll1, j, xscale = col.forest$range)
      ##
      if (newline.ll)
        add.text(ll2, j, xscale = col.forest$range)
      ##
      add.text(lr1, j, xscale = col.forest$range)
      ##
      if (newline.lr)
        add.text(lr2, j, xscale = col.forest$range)
    }
    else {
      add.label(ll1, j,
                unit(ref - (xlim[2] - xlim[1]) / 30, "native"),
                unit(ymin.line - 2.5 - (!addrow & !overall), "lines"),
                "right",
                fs.lr, ff.lr, col.label.left, xscale = col.forest$range)
      ##
      if (newline.ll)
        add.label(ll2, j,
                  unit(ref - (xlim[2] - xlim[1]) / 30, "native"),
                  unit(ymin.line - 2.5 - (!addrow & !overall) - 1, "lines"),
                  "right",
                  fs.lr, ff.lr, col.label.left, xscale = col.forest$range)
      ##
      add.label(lr1, j,
                unit(ref + (xlim[2] - xlim[1]) / 30, "native"),
                unit(ymin.line - 2.5 - (!addrow & !overall), "lines"),
                "left",
                fs.lr, ff.lr, col.label.right, xscale = col.forest$range)
      ##
      if (newline.lr)
        add.label(lr2, j,
                  unit(ref + (xlim[2] - xlim[1]) / 30, "native"),
                  unit(ymin.line - 2.5 - (!addrow & !overall) - 1, "lines"),
                  "left",
                  fs.lr, ff.lr, col.label.right, xscale = col.forest$range)
    }
  }
  ##
  add.xlab(col.forest, j, xlab, xlab.pos, fs.xlab, ff.xlab, overall, ymin.line,
           addrow, print.label, bottom.lr, newline.lr, newline.ll)
  ##
  draw.forest(col.forest, j)
  ##
  j <- j + 2
  ##
  ##
  ## Right side of forest plot
  ##
  ##
  if (rsel) {
    for (i in seq(along = rightcols)) {
      add.text(cols[[rightcols[i]]], j)
      ##
      if (!is.na(yHeadadd)) {
        if (!is.null(lab.e.attach.to.col)) {
          if (rightcols[i] == paste("col.", lab.e.attach.to.col, sep = ""))
            add.text(col.lab.e, j)
        }
        else if (metabin) {
          if (rightcols[i] == "col.n.e" & just.c == "right")
            add.text(col.lab.e, j)
          else if (rightcols[i] == "col.event.e" & just.c %in% c("left", "center"))
            add.text(col.lab.e, j)
        }
        else if (metacont) {
          if (rightcols[i] == "col.sd.e" & just.c == "right")
            add.text(col.lab.e, j)
          else if (rightcols[i] == "col.mean.e" & just.c %in% c("left", "center"))
            add.text(col.lab.e, j)
        }
        else if (metainc) {
          if (rightcols[i] == "col.time.e" & just.c == "right")
            add.text(col.lab.e, j)
          else if (rightcols[i] == "col.event.e" & just.c %in% c("left", "center"))
            add.text(col.lab.e, j)
        }
        ##
        if (!is.null(lab.c.attach.to.col)) {
          if (rightcols[i] == paste("col.", lab.c.attach.to.col, sep = ""))
            add.text(col.lab.c, j)
        }
        else if (metabin) {
          if (rightcols[i] == "col.n.c" & just.c == "right")
            add.text(col.lab.c, j)
          else if (rightcols[i] == "col.event.c" & just.c %in% c("left", "center"))
            add.text(col.lab.c, j)
        }
        else if (metacont) {
          if (rightcols[i] == "col.sd.c" & just.c == "right")
            add.text(col.lab.c, j)
          else if (rightcols[i] == "col.mean.c" & just.c %in% c("left", "center"))
            add.text(col.lab.c, j)
        }
        else if (metainc) {
          if (rightcols[i] == "col.time.c" & just.c == "right")
            add.text(col.lab.c, j)
          else if (rightcols[i] == "col.event.c" & just.c %in% c("left", "center"))
            add.text(col.lab.c, j)
        }
        ##
        if (newline.studlab & rightcols[i] == "col.studlab")
          add.text(col.add.studlab, j)
        if (newline.effect & rightcols[i] == "col.effect")
          add.text(col.add.effect, j)
        if (newline.ci & rightcols[i] == "col.ci")
          add.text(col.add.ci, j)
        if (newline.effect.ci & rightcols[i] == "col.effect.ci")
          add.text(col.add.effect.ci, j)
        if (newline.w.fixed & rightcols[i] == "col.w.fixed")
          add.text(col.add.w.fixed, j)
        if (newline.w.random & rightcols[i] == "col.w.random")
          add.text(col.add.w.random, j)
        if (newline.TE & rightcols[i] == "col.TE")
          add.text(col.add.TE, j)
        if (newline.seTE & rightcols[i] == "col.seTE")
          add.text(col.add.seTE, j)
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
        ##
        ## Add text in first line of forest plot for new columns
        ##
        if (newcols)
          if (length(rightcols.new) > 0 &
              rightcols[i] %in% paste("col.", rightcols.new, sep = "")) {
            sel <- paste("col.", rightcols.new, sep = "") == rightcols[i]
            ##
            ## Check for "\n" in label of new column
            ##
            clines <- twolines(rightlabs.new[sel], rightcols[i])
            ##
            just.new <- just.addcols.right[sel]
            ##
            if (just.new == "left")
              xpos.new <- 0
            else if (just.new == "center")
              xpos.new <- 0.5
            else if (just.new == "right")
              xpos.new <- 1
            ##
            ## Add first line
            ##
            if (clines$newline)
              add.text(tgl(clines$top, xpos.new, just.new,
                           fs.head, ff.head), j)
          }
      }
      ##
      j <- j + 2
    }
  }
  ##
  popViewport()


  invisible(NULL)
}
