metainc <- function(event.e, time.e, event.c, time.c, studlab,
                    ##
                    data = NULL, subset = NULL, exclude = NULL,
                    method = "MH",
                    ##
                    sm = gs("sminc"),
                    ##
                    incr = gs("incr"), allincr = gs("allincr"),
                    addincr = gs("addincr"),
                    model.glmm = "UM.FS",
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
                    method.bias = gs("method.bias"),
                    ##
                    n.e = NULL, n.c = NULL,
                    ##
                    backtransf = gs("backtransf"),
                    irscale = 1, irunit = "person-years",
                    title = gs("title"), complab = gs("complab"),
                    outclab = "",
                    label.e = gs("label.e"), label.c = gs("label.c"),
                    label.left = gs("label.left"),
                    label.right = gs("label.right"),
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
  ## (1) Check arguments
  ##
  ##
  chknull(sm)
  sm <- setchar(sm, c("IRR", "IRD"))
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
  method.bias <- setchar(method.bias,
                         c("rank", "linreg", "mm", "count", "score", "peters"))
  ##
  chklogical(backtransf)
  ##
  chknumeric(irscale, single = TRUE)
  chkchar(irunit)
  ##
  chklogical(keepdata)
  ##
  ## Additional arguments / checks for metainc objects
  ##
  fun <- "metainc"
  ##
  if (sm != "IRD" & irscale != 1) {
    warning("Argument 'irscale' only considered for incidence rate differences.")
    irscale <- 1
  }
  ##
  method <- setchar(method, c("Inverse", "MH", "Cochran", "GLMM"))
  if (method == "GLMM") {
    is.installed.package("lme4", fun, "method", " = \"GLMM\"")
    is.installed.package("numDeriv", fun, "method", " = \"GLMM\"")
    is.installed.package("metafor", fun, "method", " = \"GLMM\"",
                         version = .settings$metafor)
  }
  ##
  chklogical(allincr)
  chklogical(addincr)
  model.glmm <- setchar(model.glmm, c("UM.FS", "UM.RS", "CM.EL"))
  chklogical(warn)
  chkmetafor(method.tau, fun)
  ##
  if (method == "GLMM" & sm != "IRR")
    stop("Generalised linear mixed models only possible with argument 'sm = \"IRR\"'.")
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
  ## Catch 'event.e', 'time.e', 'event.c', 'time.c', 'n.e', and 'n.c' from data:
  ##
  event.e <- eval(mf[[match("event.e", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  chknull(event.e)
  k.All <- length(event.e)
  ##
  time.e <- eval(mf[[match("time.e", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknull(time.e)
  ##
  event.c <- eval(mf[[match("event.c", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  chknull(event.c)
  ##
  time.c <- eval(mf[[match("time.c", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknull(time.c)
  ##
  n.e <- eval(mf[[match("n.e", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  n.c <- eval(mf[[match("n.c", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  ##
  ## Catch 'incr' from data:
  ##
  if (!missing(incr))
    incr <- eval(mf[[match("incr", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknumeric(incr, min = 0)
  ##
  ## Catch 'studlab', 'byvar', 'subset', and 'exclude' from data:
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
  exclude <- eval(mf[[match("exclude", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  missing.exclude <- is.null(exclude)
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  chklength(time.e, k.All, fun)
  chklength(event.c, k.All, fun)
  chklength(time.c, k.All, fun)
  chklength(studlab, k.All, fun)
  ##
  if (length(incr) > 1)
    chklength(incr, k.All, fun)
  ##
  if (!missing.byvar)
    chklength(byvar, k.All, fun)
  ##
  if (!is.null(n.e))
    chklength(n.e, k.All, fun)
  if (!is.null(n.c))
    chklength(n.c, k.All, fun)
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
  ## (4) Subset, exclude studies, and subgroups
  ##
  ##
  if (!missing.subset)
    if ((is.logical(subset) & (sum(subset) > k.All)) ||
        (length(subset) > k.All))
      stop("Length of subset is larger than number of studies.")
  ##
  if (!missing.exclude) {
    if ((is.logical(exclude) & (sum(exclude) > k.All)) ||
        (length(exclude) > k.All))
      stop("Length of argument 'exclude' is larger than number of studies.")
    ##
    exclude2 <- rep(FALSE, k.All)
    exclude2[exclude] <- TRUE
    exclude <- exclude2
  }
  else
    exclude <- rep(FALSE, k.All)
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
      data <- data.frame(.event.e = event.e)
    else
      data$.event.e <- event.e
    ##
    data$.time.e <- time.e
    data$.event.c <- event.c
    data$.time.c <- time.c
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
    ##
    if (!missing.exclude)
      data$.exclude <- exclude
    ##
    if (!is.null(n.e))
      data$.n.e <- n.e
    if (!is.null(n.e))
      data$.n.c <- n.c
  }  
  
  
  ##
  ##
  ## (6) Use subset for analysis
  ##
  ##
  if (!missing.subset) {
    event.e <- event.e[subset]
    time.e <- time.e[subset]
    event.c <- event.c[subset]
    time.c <- time.c[subset]
    studlab <- studlab[subset]
    ##
    exclude <- exclude[subset]
    ##
    if (length(incr) > 1)
      incr <- incr[subset]
    ##
    if (!missing.byvar)
      byvar <- byvar[subset]
    ##
    if (!is.null(n.e))
      n.e <- n.e[subset]
    if (!is.null(n.c))
      n.c <- n.c[subset]
  }
  ##
  ## Determine total number of studies
  ##
  k.all <- length(event.e)
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
  chknumeric(event.e, 0)
  chknumeric(time.e, 0, zero = TRUE)
  chknumeric(event.c, 0)
  chknumeric(time.c, zero = TRUE)
  ##
  ## Recode integer as numeric:
  ##
  event.e <- int2num(event.e)
  time.e  <- int2num(time.e)
  event.c <- int2num(event.c)
  time.c  <- int2num(time.c)
  
  
  ##
  ##
  ## (7) Calculate results for individual studies
  ##
  ##
  sel <- switch(sm,
                IRD = event.e == 0 | event.c == 0,
                IRR = event.e == 0 | event.c == 0)
  ##
  ## Sparse computation
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
  if (sm == "IRR") {
    TE <- log(((event.e + incr.event) / time.e) / ((event.c + incr.event) / time.c))
    seTE <- sqrt(1 / (event.e + incr.event) + 1 / (event.c + incr.event))
  }
  else if (sm == "IRD") {
    TE <- event.e / time.e - event.c / time.c
    seTE <- sqrt((event.e + incr.event) / time.e^2 + (event.c + incr.event) / time.c^2)
  }
  
  
  ##
  ##
  ## (8) Do meta-analysis
  ##
  ##
  k <- sum(!is.na(event.e[!exclude]) & !is.na(event.c[!exclude]) &
           !is.na(time.e[!exclude]) & !is.na(time.c[!exclude]))
  ##
  if (method == "MH") {
    ##
    ## Greenland, Robins (1985)
    ## 
    x.k <- event.e
    y.k <- event.c
    n.k <- time.e
    m.k <- time.c
    ##
    N.k <- n.k + m.k
    t.k <- x.k + y.k
    ##
    if (sm == "IRR") {
      D <- n.k * m.k * t.k / N.k^2
      R <- x.k * m.k / N.k
      S <- y.k * n.k / N.k
      ##
      D[exclude] <- R[exclude] <- S[exclude] <- 0
      ##
      w.fixed <- S
      TE.fixed <- log(sum(R, na.rm = TRUE) / sum(S, na.rm = TRUE))
      seTE.fixed <- sqrt(sum(D, na.rm = TRUE) / (sum(R, na.rm = TRUE) *
                                                   sum(S, na.rm = TRUE)))
    }
    else if (sm == "IRD") {
      L <- (x.k * m.k^2 + y.k * n.k^2) / N.k^2
      S <- n.k * m.k / N.k
      ##
      L[exclude] <- S[exclude] <- 0
      ##
      w.fixed <- S
      TE.fixed <- weighted.mean(TE, w.fixed, na.rm = TRUE)
      seTE.fixed <- sqrt(sum(L, na.rm = TRUE) / sum(S, na.rm = TRUE)^2)
    }
  }
  ##
  else if (method == "Cochran") {
    ##
    ## Smoking and Health - Report of the Advisory Committee to the
    ## Surgeon General of the Public Health Service,
    ## Chapter 8
    ## 
    if (sm == "IRR") {
      w.fixed <- event.c * time.e / time.c
      w.fixed[exclude] <- 0
      TE.fixed <- weighted.mean(TE, w.fixed)
      seTE.fixed <- sqrt(1 / sum(event.e) + 1 / sum(event.c))
    }
    else if (sm == "IRD") {
      warning("Cochran method only available for Incidence Rate Ratio (sm = \"IRR\")")
      return(NULL)
    }
  }
  else if (method == "GLMM") {
    glmm.fixed <- metafor::rma.glmm(x1i = event.e[!exclude],
                                    t1i = time.e[!exclude],
                                    x2i = event.c[!exclude],
                                    t2i = time.c[!exclude],
                                    method = "FE",
                                    test = ifelse(hakn, "t", "z"),
                                    level = 100 * level.comb,
                                    measure = "IRR", model = model.glmm,
                                    ...)
    ##
    TE.fixed   <- as.numeric(glmm.fixed$b)
    seTE.fixed <- as.numeric(glmm.fixed$se)
    ##
    w.fixed <- rep(NA, length(time.e))
  }
  ##
  m <- metagen(TE, seTE, studlab,
               exclude = if (missing.exclude) NULL else exclude,
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
               TE.tau = if (method == "Inverse") TE.tau else TE.fixed,
               tau.common = FALSE,
               ##
               prediction = prediction,
               level.predict = level.predict,
               ##
               method.bias = method.bias,
               ##
               backtransf = backtransf,
               title = title, complab = complab, outclab = outclab,
               label.e = label.e, label.c = label.c,
               label.left = label.left, label.right = label.right,
               ##
               keepdata = FALSE,
               warn = warn)
  ##
  if (!missing.byvar & tau.common) {
    ## Estimate common tau-squared across subgroups
    hcc <- hetcalc(TE, seTE, method.tau,
                   if (method == "Inverse") TE.tau else TE.fixed,
                   byvar)
  }
  
  
  ##
  ##
  ## (9) Generate R object
  ##
  ##
  res <- list(event.e = event.e, time.e = time.e,
              event.c = event.c, time.c = time.c,
              method = method,
              incr = if (length(unique(incr)) == 1) unique(incr) else incr,
              sparse = sparse,
              allincr = allincr, addincr = addincr,
              incr.event = incr.event,
              k.MH = if (method == "MH") sum(w.fixed > 0) else NA)
  ##
  ## Add meta-analysis results
  ## (after removing unneeded list elements)
  ##
  m$n.e <- NULL
  m$n.c <- NULL
  m$method <- NULL
  ##
  res <- c(res, m)
  ##
  ## Add data
  ##
  res$n.e <- n.e
  res$n.c <- n.c
  res$TE.tau <- TE.tau
  ##
  res$irscale <- irscale
  res$irunit  <- irunit
  ##
  res$call <- match.call()
  ##
  if (method %in% c("MH", "Cochran", "GLMM")) {
    ##
    ci.f <- ci(TE.fixed, seTE.fixed, level = level.comb)
    ##
    res$TE.fixed <- TE.fixed
    res$seTE.fixed <- seTE.fixed
    res$w.fixed <- w.fixed
    res$lower.fixed <- ci.f$lower
    res$upper.fixed <- ci.f$upper
    res$zval.fixed <- ci.f$z
    res$pval.fixed <- ci.f$p
  }
  ##
  if (method == "GLMM") {
    ##
    if (sum(!exclude) > 1)
      glmm.random <- metafor::rma.glmm(x1i = event.e[!exclude],
                                       t1i = time.e[!exclude],
                                       x2i = event.c[!exclude],
                                       t2i = time.c[!exclude],
                                       method = method.tau,
                                       test = ifelse(hakn, "t", "z"),
                                       level = 100 * level.comb,
                                       measure = "IRR", model = model.glmm,
                                       ...)
    else {
      ##
      ## Fallback to fixed effect model due to small number of studies
      ##
      glmm.random <- glmm.fixed
    }
    ##
    TE.random   <- as.numeric(glmm.random$b)
    seTE.random <- as.numeric(glmm.random$se)
    ##
    ci.r <- ci(TE.random, seTE.random, level = level.comb)
    ##
    res$w.random <- rep(NA, length(event.e))
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
    if (is.null(res$lower.predict))
      res$lower.predict <- NA
    if (is.null(res$upper.predict))
      res$upper.predict <- NA
    ##
    res$model.glmm <- model.glmm
    ##
    res$Q <- glmm.random$QE.Wld
    res$df.Q <- glmm.random$QE.df
    res$Q.LRT <- glmm.random$QE.LRT
    ##
    if (k > 1)
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
    res$event.w <- NULL
    res$n.w <- NULL
  }
  ##
  class(res) <- c(fun, "meta")
  
  
  res
}
