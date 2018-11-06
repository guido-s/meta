metaprop <- function(event, n, studlab,
                     ##
                     data = NULL, subset = NULL, exclude = NULL,
                     method = "Inverse",
                     ##
                     sm = gs("smprop"),
                     ##
                     incr = gs("incr"), allincr = gs("allincr"),
                     addincr = gs("addincr"),
                     method.ci = gs("method.ci"),
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
                     pscale = 1,
                     title = gs("title"), complab = gs("complab"),
                     outclab = "",
                     ##
                     byvar, bylab, print.byvar = gs("print.byvar"),
                     byseparator = gs("byseparator"),
                     ##
                     keepdata = gs("keepdata"),
                     warn = gs("warn"),
                     ##
                     control = NULL,
                     ...
                     ) {
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  chknull(sm)
  sm <- setchar(sm, .settings$sm4prop)
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
  if (!anyNA(null.effect) | length(null.effect) != 1)
    chklevel(null.effect, name = "null.effect")
  ##
  method.bias <- setchar(method.bias,
                         c("rank", "linreg", "mm", "count", "score", "peters"))
  ##
  chklogical(backtransf)
  ##
  chknumeric(pscale, single = TRUE)
  if (!backtransf & pscale != 1 & !is.untransformed(sm)) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  ##
  chklogical(keepdata)
  ##
  ## Additional arguments / checks for metainc objects
  ##
  fun <- "metaprop"
  ##
  method <- setchar(method, c("Inverse", "GLMM"))
  is.glmm <- method == "GLMM"
  if (is.glmm) {
    is.installed.package("lme4", fun, "method", " = \"GLMM\"")
    is.installed.package("numDeriv", fun, "method", " = \"GLMM\"")
    is.installed.package("metafor", fun, "method", " = \"GLMM\"",
                         version = .settings$metafor)
  }
  ##
  chklogical(allincr)
  chklogical(addincr)
  method.ci <- setchar(method.ci, .settings$ci4prop)
  chklogical(warn)
  chkmetafor(method.tau, fun)
  ##
  if (is.glmm & sm != "PLOGIT")
    stop("Generalised linear mixed models only possible with argument 'sm = \"PLOGIT\"'.")
  ##
  if (is.glmm & method.tau != "ML")
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
  ## Catch 'event' and 'n' from data:
  ##
  event <- eval(mf[[match("event", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  chknull(event)
  k.All <- length(event)
  ##
  n <- eval(mf[[match("n", names(mf))]],
            data, enclos = sys.frame(sys.parent()))
  chknull(n)
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
  chklength(n, k.All, fun)
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
  if (is.glmm) {
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
      data <- data.frame(.event = event)
    else
      data$.event <- event
    ##
    data$.n <- n
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
  }
  
  
  ##
  ##
  ## (6) Use subset for analysis
  ##
  ##
  if (!missing.subset) {
    event <- event[subset]
    n   <- n[subset]
    studlab <- studlab[subset]
    ##
    exclude <- exclude[subset]
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
  chknumeric(n, 0, zero = TRUE)
  ##
  if (any(n < 10, na.rm = TRUE) & sm == "PFT")
    warning("Sample size very small (below 10) in at least one study. Accordingly, back transformation for pooled effect may be misleading for Freeman-Tukey double arcsine transformation. Please look at results for other transformations (e.g. sm = 'PAS' or sm = 'PLOGIT'), too.")
  ##
  if (any(event > n, na.rm = TRUE))
    stop("Number of events must not be larger than number of observations")
  ##
  ## Recode integer as numeric:
  ##
  event <- int2num(event)
  n     <- int2num(n)
  
  
  ##
  ##
  ## (7) Calculate results for individual studies
  ##
  ##
  sel <- switch(sm,
                PFT = rep(FALSE, length(event)),
                PAS = rep(FALSE, length(event)),
                PRAW =   event == 0 | (n - event) == 0,
                PLN =    event == 0 | (n - event) == 0,
                PLOGIT = event == 0 | (n - event) == 0)
  ##
  sparse <- any(sel, na.rm = TRUE)
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
  if (sm == "PFT") {
    TE <- 0.5 * (asin(sqrt(event / (n + 1))) + asin(sqrt((event + 1) / (n + 1))))
    seTE <- sqrt(1 / (4 * n + 2))
    transf.null.effect <- asin(sqrt(null.effect))
  }
  else if (sm == "PAS") {
    TE <- asin(sqrt(event / n))
    seTE <- sqrt(0.25 * (1 / n))
    transf.null.effect <- asin(sqrt(null.effect))
  }
  else if (sm == "PRAW") {
    TE <- event / n
    seTE <- sqrt((event + incr.event) * (n - event + incr.event) /
                   (n + 2 * incr.event)^3)
    transf.null.effect <- null.effect
  }
  else if (sm == "PLN") {
    TE <- log((event + incr.event) / (n + incr.event))
    ## Hartung, Knapp (2001), p. 3880, formula (18):
    seTE <- ifelse(event == n,
                   sqrt(1 / event                - 1 / (n + incr.event)),
                   sqrt(1 / (event + incr.event) - 1 / (n + incr.event))
                   )
    transf.null.effect <- log(null.effect)
  }
  else if (sm == "PLOGIT") {
    TE <- log((event + incr.event) / (n - event + incr.event))
    seTE <- sqrt(1 / (event + incr.event) +
                 1 / ((n - event + incr.event)))
    transf.null.effect <- log(null.effect / (1 - null.effect))
  }
  ##
  ## Calculate confidence intervals
  ##
  NAs <- rep(NA, k.all)
  ##
  if (method.ci == "CP") {
    lower.study <- upper.study <- NAs
    for (i in 1:k.all) {
      if (!is.na(event[i] & !is.na(n[i]))) {
        cint <- binom.test(event[i], n[i], conf.level = level)
        ##
        lower.study[i] <- cint$conf.int[[1]]
        upper.study[i] <- cint$conf.int[[2]]
      }
      else {
        lower.study[i] <- NA
        upper.study[i] <- NA
      }
    }
  }
  ##
  else {
    if (method.ci == "WS")
      ci.study <- ciWilsonScore(event, n, level = level)
    ##
    else if (method.ci == "WSCC")
      ci.study <- ciWilsonScore(event, n, level = level, correct = TRUE)
    ##
    else if (method.ci == "AC")
      ci.study <- ciAgrestiCoull(event, n, level = level)
    ##
    else if (method.ci == "SA")
      ci.study <- ciSimpleAsymptotic(event, n, level = level)
    ##
    else if (method.ci == "SACC")
      ci.study <- ciSimpleAsymptotic(event, n, level = level, correct = TRUE)
    ##
    else if (method.ci == "NAsm")
      ci.study <- ci(TE, seTE, level = level)
    ##
    lower.study <- ci.study$lower
    upper.study <- ci.study$upper
  }
  ##  
  if (method.ci == "NAsm") {
    if (sm == "PLN") {
      lower.study <- exp(lower.study)
      upper.study <- exp(upper.study)
    }
    ##
    else if (sm == "PLOGIT") {
      lower.study <- logit2p(lower.study)
      upper.study <- logit2p(upper.study)
    }
    ##
    else if (sm == "PAS") {
      lower.study <- asin2p(lower.study, value = "lower", warn = FALSE)
      upper.study <- asin2p(upper.study, value = "upper", warn = FALSE)
    }
    ##
    else if (sm == "PFT") {
      lower.study <- asin2p(lower.study, n, value = "lower", warn = FALSE)
      upper.study <- asin2p(upper.study, n, value = "upper", warn = FALSE)
    }
    ##
    lower.study[lower.study < 0] <- 0
    upper.study[upper.study > 1] <- 1
  }
  
  
  ##
  ##
  ## (8) Do meta-analysis
  ##
  ##
  k <- sum(!is.na(event[!exclude]) & !is.na(n[!exclude]))
  ##
  if (is.glmm & k > 0) {
    glmm.fixed <- metafor::rma.glmm(xi = event[!exclude], ni = n[!exclude],
                                    method = "FE",
                                    test = ifelse(hakn, "t", "z"),
                                    level = 100 * level.comb,
                                    measure = "PLO",
                                    control = control,
                                    ...)
    ##
    TE.fixed   <- as.numeric(glmm.fixed$b)
    seTE.fixed <- as.numeric(glmm.fixed$se)
    ##
    w.fixed <- rep(NA, length(event))
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
               warn = warn,
               ##
               control = control)
  ##
  if (method != "GLMM" & !missing.byvar & tau.common) {
    ## Estimate common tau-squared across subgroups
    hcc <- hetcalc(TE, seTE, method.tau, TE.tau, byvar,
                   control = control)
  }
  
  
  ##
  ##
  ## (9) Generate R object
  ##
  ##
  res <- list(event = event, n = n,
              incr = if (length(unique(incr)) == 1) unique(incr) else incr,
              sparse = sparse,
              allincr = allincr, addincr = addincr,
              method.ci = method.ci,
              incr.event = incr.event)
  ##
  ## Add meta-analysis results
  ## (after removing unneeded list elements)
  ##
  m$n.e <- NULL
  m$n.c <- NULL
  m$label.e <- ""
  m$label.c <- ""
  m$label.left <- ""
  m$label.right <- ""
  ##
  res <- c(res, m)
  res$null.effect <- null.effect
  ##
  ## Add data
  ##
  if (is.glmm & k > 0) {
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
    if (sum(!exclude) > 1)
      glmm.random <- metafor::rma.glmm(xi = event[!exclude], ni = n[!exclude],
                                       method = method.tau,
                                       test = ifelse(hakn, "t", "z"),
                                       level = 100 * level.comb,
                                       measure = "PLO",
                                       control = control,
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
    if (is.null(res$lower.predict))
      res$lower.predict <- NA
    if (is.null(res$upper.predict))
      res$upper.predict <- NA
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
  res$lower <- lower.study
  res$upper <- upper.study
  ##
  res$pscale <- pscale
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
      if (is.glmm)
        res <- c(res, subgroup(res, NULL,
                               factor(res$byvar, bylevs(res$byvar)), ...))
      else {
        res <- c(res, subgroup(res, hcc$tau))
        res$Q.w.random <- hcc$Q
        res$df.Q.w.random <- hcc$df.Q
      }
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
