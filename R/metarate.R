metarate <- function(event, time, studlab,
                     ##
                     data = NULL, subset = NULL, method = "Inverse",
                     ##
                     sm = gs("smrate"),
                     ##
                     incr = gs("incr"), allincr = gs("allincr"),
                     addincr = gs("addincr"),
                     ##
                     level = gs("level"), level.comb = gs("level.comb"),
                     comb.fixed = gs("comb.fixed"),
                     comb.random = gs("comb.random"),
                     ##
                     hakn = gs("hakn"),
                     method.tau =
                       ifelse(!is.na(charmatch(tolower(method), "glmm",
                                               nomatch = NA)),
                              "ML", gs("method.tau")),
                     tau.preset = NULL, TE.tau = NULL,
                     tau.common = gs("tau.common"),
                     ##
                     prediction = gs("prediction"),
                     level.predict = gs("level.predict"),
                     ##
                     null.effect = NA,
                     ##
                     method.bias = gs("method.bias"),
                     ##
                     backtransf = gs("backtransf"),
                     irscale = 1, irunit = "person-years",
                     title = gs("title"), complab = gs("complab"),
                     outclab = "",
                     ##
                     byvar, bylab, print.byvar = gs("print.byvar"),
                     byseparator = gs("byseparator"),
                     ##
                     keepdata = gs("keepdata"),
                     warn = gs("warn"),
                     ...
                     ) {
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  chklevel(level)
  chklevel(level.comb)
  chklogical(comb.fixed)
  chklogical(comb.random)
  ##
  chklogical(hakn)
  method.tau <- setchar(method.tau,
                        c("DL", "PM", "REML", "ML", "HS", "SJ", "HE", "EB"))
  chklogical(tau.common)
  ##
  chklogical(prediction)
  chklevel(level.predict)
  ##
  chknumeric(null.effect, single = TRUE)
  ##
  method.bias <- setchar(method.bias,
                         c("rank", "linreg", "mm", "count", "score", "peters"))
  ##
  chklogical(backtransf)
  chknumeric(irscale, single = TRUE)
  if (!backtransf & irscale != 1) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  ##
  chklogical(keepdata)
  ##
  ## Additional arguments / checks for metainc objects
  ##
  fun <- "metarate"
  ##
  method <- setchar(method, c("Inverse", "GLMM"))
  if (method == "GLMM") {
    is.installed.package("lme4", fun, "method", " = \"GLMM\"")
    is.installed.package("numDeriv", fun, "method", " = \"GLMM\"")
    is.installed.package("metafor", fun, "method", " = \"GLMM\"",
                         version = .settings$metafor)
  }
  ##
  sm <- setchar(sm, c("IR", "IRLN", "IRS", "IRFT"))
  chklogical(allincr)
  chklogical(addincr)
  chklogical(warn)
  chkmetafor(method.tau, fun)
  ##
  if (method == "GLMM" & sm != "IRLN")
    stop("Generalised linear mixed models only possible with argument 'sm = \"IRLN\"'.")
  ##
  if (method == "GLMM" & method.tau != "ML")
    stop("Generalised linear mixed models only possible with argument 'method.tau = \"ML\"'.")
  
  
  ##
  ##
  ## (2) Read data
  ##
  ##
  nulldata <- is.null(data)
  ##
  if (nulldata)
    data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Catch event, n from data:
  ##
  event <- eval(mf[[match("event", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  chknull(event)
  k.All <- length(event)
  ##
  time <- eval(mf[[match("time", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  chknull(time)
  ##
  ## Catch incr from data:
  ##
  if (!missing(incr))
    incr <- eval(mf[[match("incr", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknumeric(incr, min = 0)
  ##
  ## Catch studlab, byvar, subset from data:
  ##
  studlab <- eval(mf[[match("studlab", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  studlab <- setstudlab(studlab, k.All)
  ##
  byvar <- eval(mf[[match("byvar", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  missing.byvar <- is.null(byvar)
  if (method == "GLMM" & !missing.byvar) {
    warning("Argument 'byvar' not considered for GLMMs. Use metareg function for subgroup analysis of GLMM meta-analyses.")
    byvar <- NULL
    missing.byvar <- is.null(byvar)
  }
  ##
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  missing.subset <- is.null(subset)
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  chklength(time, k.All, fun)
  chklength(studlab, k.All, fun)
  ##
  if (length(incr) > 1)
    chklength(incr, k.All, fun)
  ##
  if (!missing.byvar)
    chklength(byvar, k.All, fun)
  ##
  ## Additional checks
  ##
  if (method == "GLMM") {
    if (tau.common) {
      if (warn)
        warning("Argument 'tau.common' not considered for GLMM.")
      tau.common <- FALSE
    }
    if (!is.null(TE.tau)) {
      if (warn)
        warning("Argument 'TE.tau' not considered for GLMM.")
      TE.tau <- NULL
    }
    ##
    if (!is.null(tau.preset)) {
      if (warn)
        warning("Argument 'tau.preset' not considered for GLMM.")
      tau.preset <- NULL
    }
  }
  if (missing.byvar & tau.common) {
    warning("Value for argument 'tau.common' set to FALSE as argument 'byvar' is missing.")
    tau.common <- FALSE
  }
  if (!missing.byvar & !tau.common & !is.null(tau.preset)) {
    warning("Argument 'tau.common' set to TRUE as argument tau.preset is not NULL.")
    tau.common <- TRUE
  }
  
  
  ##
  ##
  ## (4) Subset and subgroups
  ##
  ##
  if (!missing.subset)
    if ((is.logical(subset) & (sum(subset) > k.All)) ||
        (length(subset) > k.All))
      stop("Length of subset is larger than number of studies.")
  ##  
  if (!missing.byvar) {
    chkmiss(byvar)
    byvar.name <- byvarname(mf[[match("byvar", names(mf))]])
    bylab <- if (!missing(bylab) && !is.null(bylab)) bylab else byvar.name
  }
  
  
  ##
  ##
  ## (5) Store complete dataset in list object data
  ##     (if argument keepdata is TRUE)
  ##
  ##
  if (keepdata) {
    if (nulldata)
      data <- data.frame(.event = event)
    else
      data$.event <- event
    ##
    data$.time <- time
    data$.studlab <- studlab
    ##
    data$.incr <- incr
    ##
    if (!missing.byvar)
      data$.byvar <- byvar
    ##
    if (!missing.subset) {
      if (length(subset) == dim(data)[1])
        data$.subset <- subset
      else {
        data$.subset <- FALSE
        data$.subset[subset] <- TRUE
      }
    }
  }
  
  
  ##
  ##
  ## (6) Use subset for analysis
  ##
  ##
  if (!missing.subset) {
    event <- event[subset]
    time  <- time[subset]
    studlab <- studlab[subset]
    ##
    if (length(incr) > 1)
      incr <- incr[subset]
    ##
    if (!missing.byvar)
      byvar <- byvar[subset]
  }
  ##
  ## Determine total number of studies
  ##
  k.all <- length(event)
  ##
  if (k.all == 0)
    stop("No studies to combine in meta-analysis.")
  ##
  ## No meta-analysis for a single study
  ##
  if (k.all == 1) {
    comb.fixed  <- FALSE
    comb.random <- FALSE
    prediction  <- FALSE
  }
  ##
  ## Check variable values
  ##
  chknumeric(event, 0)
  chknumeric(time, 0, zero = TRUE)
  ##
  ## Recode integer as numeric:
  ##
  event <- int2num(event)
  time  <- int2num(time)
  
  
  ##
  ##
  ## (7) Calculate results for individual studies
  ##
  ##
  sel <- switch(sm,
                IR   = event == 0,
                IRLN = event == 0,
                IRS  = rep(FALSE, length(event)),
                IRFT = rep(FALSE, length(event)))
  ##
  sparse <- any(sel, na.rm = TRUE)
  ##
  if (method == "GLMM" & sparse)
    if ((!missing(incr) & any(incr != 0)) |
        (!missing(allincr) & allincr ) |
        (!missing(addincr) & addincr)
        )
      warning("Note, for method = \"GLMM\", continuity correction only used to calculate individual study results.")
  ##
  ## No need to add anything to cell counts for arcsine transformation
  ##
  if (addincr)
    incr.event <- if (length(incr) == 1) rep(incr, k.all) else incr
  else
    if (sparse)
      if (allincr)
        incr.event <- if (length(incr) == 1) rep(incr, k.all) else incr
      else
        incr.event <- incr * sel
    else
      incr.event <- rep(0, k.all)
  ##
  if (sm == "IR") {
    TE <- (event + incr.event) / time
    seTE <- sqrt(TE / time)
    transf.null.effect <- null.effect
  }
  else if (sm == "IRLN") {
    TE <- log((event + incr.event) / time)
    seTE <- sqrt(1 / (event + incr.event))
    transf.null.effect <- log(null.effect)
  }
  else if (sm == "IRS") {
    TE <- sqrt(event / time)
    seTE <- sqrt(1 / (4 * time))
    transf.null.effect <- sqrt(null.effect)
  }
  else if (sm == "IRFT") {
    TE <- 0.5 * (sqrt(event / time) + sqrt((event + 1) / time))
    seTE <- sqrt(1 / (4 * time))
    transf.null.effect <- sqrt(null.effect)
  }
  ##
  ## Calculate confidence intervals
  ##
  ci.study <- ci(TE, seTE, level = level)
  ##
  lower.study <- ci.study$lower
  upper.study <- ci.study$upper
  
  
  ##
  ##
  ## (8) Do meta-analysis
  ##
  ##
  if (method == "GLMM") {
    glmm.fixed <- metafor::rma.glmm(xi = event, ti = time,
                                    method = "FE",
                                    test = ifelse(hakn, "t", "z"),
                                    level = 100 * level.comb,
                                    measure = "IRLN",
                                    ...)
    ##
    TE.fixed   <- as.numeric(glmm.fixed$b)
    seTE.fixed <- as.numeric(glmm.fixed$se)
    ##
    w.fixed <- rep(NA, length(event))
  }
  ##
  m <- metagen(TE, seTE, studlab,
               ##
               sm = sm,
               level = level,
               level.comb = level.comb,
               comb.fixed = comb.fixed,
               comb.random = comb.random,
               ##
               hakn = hakn,
               method.tau = method.tau,
               tau.preset = tau.preset,
               TE.tau = TE.tau,
               tau.common = FALSE,
               ##
               prediction = prediction,
               level.predict = level.predict,
               ##
               null.effect = transf.null.effect,
               ##
               method.bias = method.bias,
               ##
               backtransf = backtransf,
               title = title, complab = complab, outclab = outclab,
               ##
               keepdata = FALSE,
               warn = warn)
  ##
  if (!missing.byvar & tau.common) {
    ## Estimate common tau-squared across subgroups
    hcc <- hetcalc(TE, seTE, method.tau,
                   TE.tau,
                   byvar)
  }
  
  
  ##
  ##
  ## (9) Generate R object
  ##
  ##
  res <- list(event = event, time = time,
              incr = if (length(unique(incr)) == 1) unique(incr) else incr,
              sparse = sparse,
              allincr = allincr, addincr = addincr,
              incr.event = incr.event)
  ##
  ## Add meta-analysis results
  ## (after removing unneeded list elements)
  ##
  m$n.e <- NULL
  m$n.c <- NULL
  m$label.e <- NULL
  m$label.c <- NULL
  m$label.left <- NULL
  m$label.right <- NULL
  ##
  res <- c(res, m)
  res$null.effect <- null.effect
  ##
  ## Add data
  ##
  ##
  if (method == "GLMM") {
    ##
    ci.f <- ci(TE.fixed, seTE.fixed, level = level.comb,
               null.effect = transf.null.effect)
    ##
    res$method <- "GLMM"
    ##
    res$TE.fixed <- TE.fixed
    res$seTE.fixed <- seTE.fixed
    res$w.fixed <- w.fixed
    res$lower.fixed <- ci.f$lower
    res$upper.fixed <- ci.f$upper
    res$zval.fixed <- ci.f$z
    res$pval.fixed <- ci.f$p
    ##
    glmm.random <- metafor::rma.glmm(xi = event, ti = time,
                                     method = method.tau,
                                     test = ifelse(hakn, "t", "z"),
                                     level = 100 * level.comb,
                                     measure = "IRLN",
                                     ...)
    ##
    TE.random   <- as.numeric(glmm.random$b)
    seTE.random <- as.numeric(glmm.random$se)
    ##
    ci.r <- ci(TE.random, seTE.random, level = level.comb,
               null.effect = transf.null.effect)
    ##
    res$w.random <- rep(NA, length(event))
    ##
    res$TE.random <- TE.random
    res$seTE.random <- seTE.random
    res$lower.random <- ci.r$lower
    res$upper.random <- ci.r$upper
    res$zval.random <- ci.r$z
    res$pval.random <- ci.r$p
    ##
    res$se.tau2 <- NA
    ci.p <- metafor::predict.rma(glmm.random, level = 100 * level.predict)
    res$seTE.predict <- NA
    res$lower.predict <- ci.p$cr.lb
    res$upper.predict <- ci.p$cr.ub
    ##
    res$Q <- glmm.random$QE.Wld
    res$df.Q <- glmm.random$QE.df
    res$Q.LRT <- glmm.random$QE.LRT
    ##
    res$tau <- sqrt(glmm.random$tau2)
    ##
    res$H <- sqrt(glmm.random$H2)
    res$lower.H <- NA
    res$upper.H <- NA
    ##
    res$I2 <- glmm.random$I2 / 100
    res$lower.I2 <- NA
    res$upper.I2 <- NA
    ##
    res$.glmm.fixed  <- glmm.fixed
    res$.glmm.random <- glmm.random
    res$version.metafor <- packageDescription("metafor")$Version
  }
  ##
  res$lower <- lower.study
  res$upper <- upper.study
  ##
  res$irscale <- irscale
  res$irunit   <- irunit
  ##
  res$call <- match.call()
  ##
  if (keepdata) {
    res$data <- data
    if (!missing.subset)
      res$subset <- subset
  }
  ##
  class(res) <- c(fun, "meta")
  ##
  ## Add results from subgroup analysis
  ##
  if (!missing.byvar) {
    res$byvar <- byvar
    res$bylab <- bylab
    res$print.byvar <- print.byvar
    res$byseparator <- byseparator
    res$tau.common <- tau.common
    ##
    if (!tau.common)
      res <- c(res, subgroup(res))
    else if (!is.null(tau.preset))
      res <- c(res, subgroup(res, tau.preset))
    else {
      res <- c(res, subgroup(res, hcc$tau))
      res$Q.w.random <- hcc$Q
      res$df.Q.w.random <- hcc$df.Q
    }
    ##
    res$event.e.w <- NULL
    res$event.c.w <- NULL
    res$n.e.w <- NULL
    res$n.c.w <- NULL
    res$time.e.w <- NULL
    res$time.c.w <- NULL
  }
  ##
  class(res) <- c(fun, "meta")

  
  res
}
