metabin <- function(event.e, n.e, event.c, n.c, studlab,
                    ##
                    data = NULL, subset = NULL,
                    ##
                    method = ifelse(tau.common, "Inverse", .settings$method),
                    sm = 
                      ifelse(!is.na(charmatch(tolower(method), c("peto", "glmm"),
                                              nomatch = NA)),
                             "OR", .settings$smbin),
                    incr = .settings$incr, allincr = .settings$allincr,
                    addincr = .settings$addincr, allstudies = .settings$allstudies,
                    MH.exact = .settings$MH.exact, RR.cochrane = .settings$RR.cochrane,
                    model.glmm = "UM.FS",
                    ##
                    level = .settings$level, level.comb = .settings$level.comb,
                    comb.fixed = .settings$comb.fixed,
                    comb.random = .settings$comb.random,
                    ##
                    hakn = .settings$hakn,
                    method.tau =
                      ifelse(!is.na(charmatch(tolower(method), "glmm",
                                              nomatch = NA)),
                             "ML", .settings$method.tau),
                    tau.preset = NULL, TE.tau = NULL,
                    tau.common = .settings$tau.common,
                    ##
                    prediction = .settings$prediction,
                    level.predict = .settings$level.predict,
                    ##
                    method.bias = ifelse(sm == "OR", "score", .settings$method.bias),
                    ##
                    backtransf = .settings$backtransf,
                    title = .settings$title, complab = .settings$complab,
                    outclab = "",
                    label.e = .settings$label.e, label.c = .settings$label.c,
                    label.left = .settings$label.left, label.right = .settings$label.right,
                    ##
                    byvar, bylab, print.byvar = .settings$print.byvar,
                    byseparator = .settings$byseparator,
                    ##
                    print.CMH = .settings$print.CMH,
                    ##
                    keepdata = .settings$keepdata,
                    warn = .settings$warn,
                    ...
                    ) {
  
  
  ##
  ##
  ## (1) Check arguments
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
  method.bias <- setchar(method.bias,
                         c("rank", "linreg", "mm", "count", "score", "peters"))
  ##
  chklogical(backtransf)
  chklogical(keepdata)
  ##
  ## Additional arguments / checks for metabin objects
  ##
  fun <- "metabin"
  ##
  sm <- setchar(sm, c("OR", "RD", "RR", "ASD"))
  ##
  method <- setchar(method, c("Inverse", "MH", "Peto", "GLMM"))
  if (method == "GLMM") {
    is.installed.package("lme4", fun, "method", " = \"GLMM\"")
    is.installed.package("numDeriv", fun, "method", " = \"GLMM\"")
    is.installed.package("metafor", fun, "method", " = \"GLMM\"")
  }
  ##
  chklogical(allincr)
  chklogical(addincr)
  chklogical(allstudies)
  chklogical(MH.exact)
  chklogical(RR.cochrane)
  ##
  model.glmm <- setchar(model.glmm, c("UM.FS", "UM.RS", "CM.EL", "CM.AL"))
  if (method == "GLMM" & model.glmm == "CM.EL")
    is.installed.package("BiasedUrn", fun, "model.glmm", " = \"CM.EL\"")
  ##
  chklogical(print.CMH)
  chklogical(warn)
  chkmetafor(method.tau, fun)
  ##
  if (sm == "ASD")
    method <- "Inverse"
  ##
  if (!is.numeric(incr))
    incr <- setchar(incr, "TACC",
                    "should be numeric or the character string \"TACC\"")
  ##
  if (method == "Peto" & sm != "OR")
    stop("Peto's method only possible with argument 'sm = \"OR\"'")
  ##
  if (method == "GLMM" & sm != "OR")
    stop("Generalised linear mixed models only possible with argument 'sm = \"OR\"'.")
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
  ## Catch event.e, n.e, event.c, n.c from data:
  ##
  event.e <- eval(mf[[match("event.e", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  chknull(event.e)
  k.All <- length(event.e)
  ##
  n.e <- eval(mf[[match("n.e", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  chknull(n.e)
  ##
  event.c <- eval(mf[[match("event.c", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  chknull(event.c)
  ##
  n.c <- eval(mf[[match("n.c", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  chknull(n.c)
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
  chklength(n.e, k.All, fun)
  chklength(event.c, k.All, fun)
  chklength(n.c, k.All, fun)
  chklength(studlab, k.All, fun)
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
      data <- data.frame(.event.e = event.e)
    else
      data$.event.e <- event.e
    ##
    data$.n.e <- n.e
    data$.event.c <- event.c
    data$.n.c <- n.c
    data$.studlab <- studlab
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
    event.e <- event.e[subset]
    n.e <- n.e[subset]
    event.c <- event.c[subset]
    n.c <- n.c[subset]
    studlab <- studlab[subset]
    ##
    if (!missing.byvar)
      byvar <- byvar[subset]
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
    ##
    if (method == "MH")
      method <- "Inverse"
  }
  ##
  ## Check variable values
  ##
  chknumeric(event.e)
  chknumeric(n.e)
  chknumeric(event.c)
  chknumeric(n.c)
  ##
  ## Recode integer as numeric:
  ##
  event.e <- int2num(event.e)
  n.e     <- int2num(n.e)
  event.c <- int2num(event.c)
  n.c     <- int2num(n.c)
  
  
  ##
  ##
  ## (7) Calculate results for individual studies
  ##
  ##
  ## Include non-informative studies?
  ## (i.e. studies with either zero or all events in both groups)
  ##
  if (sm == "RD" | sm == "ASD")
    incl <- rep(1, k.all)
  else {
    allevents <- event.c == n.c & event.e == n.e
    if (allstudies)
      incl <- rep(1, k.all)
    else {
      if (sm == "OR")
        incl <- ifelse((event.c == 0   & event.e == 0) |
                         (event.c == n.c & event.e == n.e), NA, 1)
      if (sm == "RR")
        incl <- ifelse((event.c == 0 & event.e == 0), NA, 1)
    }
  }
  ##
  ## Exclude studies from meta-analysis:
  ##
  sel1 <- event.e > n.e
  sel2 <- event.c > n.c
  if ((any(sel1, na.rm = TRUE)) & warn)
    warning("Studies with event.e > n.e get no weight in meta-analysis.")
  if ((any(sel2, na.rm = TRUE)) & warn)
    warning("Studies with event.c > n.c get no weight in meta-analysis.")
  incl[sel1 | sel2] <- NA
  ##
  sel3 <- n.e <= 0 | n.c <= 0
  if ((any(sel3, na.rm = TRUE)) & warn)
    warning("Studies with non-positive values for n.e and / or n.c get no weight in meta-analysis.")
  incl[sel3] <- NA
  ##
  sel4 <- event.e < 0 | event.c < 0
  if ((any(sel4, na.rm = TRUE)) & warn)
    warning("Studies with negative values for event.e and / or event.c get no weight in meta-analysis.")
  incl[sel4] <- NA
  ##
  ## Sparse computation
  ##
  sel <- switch(sm,
                OR = ((n.e - event.e) == 0 | event.e == 0 |
                        (n.c - event.c) == 0 | event.c == 0),
                RD = ((n.e - event.e) == 0 | event.e == 0 |
                        (n.c - event.c) == 0 | event.c == 0),
                RR = ((n.e - event.e) == 0 | event.e == 0 |
                        (n.c - event.c) == 0 | event.c == 0),
                ASD = rep(FALSE, length(event.e)))
  ##
  sel[is.na(incl)] <- FALSE
  ##
  sparse <- any(sel, na.rm = TRUE)
  ##
  ## Check for studies with zero cell frequencies in both groups
  ##
  doublezeros <- FALSE
  if (sparse & sm %in% c("RR", "OR") & !(method %in% c("Peto", "GLMM"))) {
    sel.doublezeros <- switch(sm,
                              OR = (event.e == 0   & event.c ==   0) |
                                (event.c == n.c & event.e == n.e),
                              RR = (event.c == 0 & event.e == 0))
    if (any(sel.doublezeros))
      doublezeros <- TRUE
  }
  ##
  ## No need to add anything to cell counts for
  ##  (i)  arcsine difference as summary measure
  ##  (ii) Peto method or GLMM
  ##
  if (sm == "ASD" | method %in% c("Peto", "GLMM")) {
    if ((!missing(incr) & incr != 0) |
        (!missing(allincr) & allincr ) |
        (!missing(addincr) & addincr) |
        (!missing(allstudies) & allstudies)
        )
      if (sm == "ASD") {
        if (sparse | addincr) {
          warning("Note, no continuity correction considered for arcsine difference (sm = \"ASD\").")
        }
      }
      else if (method == "Peto") {
        if (sparse | addincr)
          warning("Note, no continuity correction considered for method = \"Peto\".")
      }
      else if (method == "GLMM") {
        if (sparse | addincr)
          warning("Note, for method = \"GLMM\", continuity correction only used to calculate individual study results.")
      }
  }
  ##
  ## Define continuity correction
  ##
  if (addincr) {
    ##
    if (is.numeric(incr)) {
      incr.e <- rep(incr, k.all)
      incr.c <- rep(incr, k.all)
    }
    else {
      if (incr == "TACC") {
        ##
        ## Treatment arm continuity correction:
        ##
        incr.e <- n.e / (n.e + n.c)
        incr.c <- n.c / (n.e + n.c)
      }
    }
  }
  else {
    if (sparse) {
      if (allincr) {
        ##
        if (is.numeric(incr)) {
          incr.e <- rep(incr, k.all)
          incr.c <- rep(incr, k.all)
        }
        else {
          if (incr == "TACC") {
            ##
            ## Treatment arm continuity correction:
            ##
            incr.e <- n.e / (n.e + n.c)
            incr.c <- n.c / (n.e + n.c)
          }
        }
      }
      else {
        ##
        ## Bradburn, Deeks, Altman, Stata-procedure "metan":
        ## & SAS PROC FREQ (for method = "Inverse")
        ##
        if (is.numeric(incr)) {
          incr.e <- incr * sel
          incr.c <- incr * sel
        }
        else {
          if (incr == "TACC") {
            ##
            ## Treatment arm continuity correction:
            ##
            incr.e <- n.e / (n.e + n.c) * sel
            incr.c <- n.c / (n.e + n.c) * sel
          }
        }
      }
    }
    else {
      incr.e <- rep(0, k.all)
      incr.c <- rep(0, k.all)
    }
  }
  ##
  ## No continuity correction for Peto method
  ##
  if (method == "Peto") {
    incr.e <- rep(0, k.all)
    incr.c <- rep(0, k.all)
  }
  ##  
  n11 <- event.e * incl
  n21 <- event.c * incl
  n1. <- n.e * incl
  n2. <- n.c * incl
  ##
  n.. <- n1. + n2.
  n12 <- n1. - n11
  n22 <- n2. - n21
  n.1 <- n11 + n21
  n.2 <- n12 + n22
  ##
  Q.CMH <- (sum(n11 - n1. * n.1 / n.., na.rm = TRUE)^2 /
              sum(n1. * n2. * n.1 * n.2 / n..^3, na.rm = TRUE))
  ##
  ## Estimation of treatment effects in individual studies
  ##
  if (sm == "OR") {
    if (method %in% c("MH", "Inverse", "GLMM")) {
      ## 
      ## Cooper & Hedges (1994), p. 251-2
      ## 
      TE <- log(((n11 + incr.e) * (n22 + incr.c)) /
                  ((n12 + incr.e) * (n21 + incr.c)))
      seTE <- sqrt((1 / (n11 + incr.e) + 1 / (n12 + incr.e) +
                      1 / (n21 + incr.c) + 1 / (n22 + incr.c)))
    }
    else if (method == "Peto") {
      ## 
      ## Cooper & Hedges (1994), p. 252
      ## 
      O <- n11
      E <- n1. * n.1 / n..
      V <- n1. * n2. * n.1 * n.2 / ((n.. - 1) * n..^2)
      ##
      TE <- (O - E) / V
      seTE <- sqrt(1 / V)
    }
  }
  else if (sm == "RR") {
    ## 
    ## Cooper & Hedges (1994), p. 247-8
    ##
    if (!RR.cochrane) {
      TE <- log(((n11 + incr.e) / (n1. + incr.e)) /
                  ((n21 + incr.c) / (n2. + incr.c)))
      ##
      ## Hartung & Knapp (2001), Stat Med, equation (18)
      ##
      seTE <- sqrt((1 / (n11 + incr.e * (!allevents)) - 1 / (n1. + incr.e) +
                      1 / (n21 + incr.c * (!allevents)) - 1 / (n2. + incr.c)))
    }
    else {
      TE <- log(((n11 + incr.e) / (n1. + 2 * incr.e)) /
                  ((n21 + incr.c) / (n2. + 2 * incr.c)))
      seTE <- sqrt((1 / (n11 + incr.e) - 1 / (n1. + 2 * incr.e) +
                      1 / (n21 + incr.c) - 1 / (n2. + 2 * incr.c)))
    }
  }
  else if (sm == "RD") {
    ## 
    ## Cooper & Hedges (1994), p. 246-7
    ## 
    TE <- n11 / n1. - n21 / n2.
    seTE <- sqrt((n11 + incr.e) * (n12 + incr.e) / (n1. + 2 * incr.e)^3 +
                   (n21 + incr.c) * (n22 + incr.c) / (n2. + 2 * incr.c)^3)
  }
  else if (sm == "ASD") {
    ## 
    ## Ruecker et al. (2009)
    ## 
    TE <- asin(sqrt(n11 / n1.)) - asin(sqrt(n21 / n2.))
    seTE <- sqrt(0.25 * (1 / n1. + 1 / n2.))
  }
  
  
  ##
  ##
  ## (8) Do meta-analysis
  ##
  ##
  if (sum(!is.na(TE)) == 1 & k.all > 1 & method == "MH") {
    if (warn)
      warning("For a single study, inverse variance method used instead of Mantel-Haenszel method.")
    method <- "Inverse"
  }
  ##  
  if (method == "MH") {
    incr.e <- incr.e * (!MH.exact)
    incr.c <- incr.c * (!MH.exact)
    ##
    if (sm == "OR") {
      ## 
      ## Cooper & Hedges (1994), p. 253-5 (MH.exact == TRUE)
      ##
      ## Bradburn, Deeks, Altman, Stata-procedure "metan"
      ## und RevMan 3.1 (MH.exact == FALSE)
      ## 
      A <- (n11 + incr.e) * (n22 + incr.c) / (n.. + 2 * incr.e + 2 * incr.c)
      B <- (n11 + incr.e + n22 + incr.c) / (n.. + 2 * incr.e + 2 * incr.c)
      C <- (n12 + incr.e) * (n21 + incr.c) / (n.. + 2 * incr.e + 2 * incr.c)
      D <- (n12 + incr.e + n21 + incr.c) / (n.. + 2 * incr.e + 2 * incr.c)
      ##
      ## Cooper & Hedges (1994), p. 265-6
      ##
      w.fixed <- C
      TE.fixed <- log(sum(A, na.rm = TRUE) / sum(C, na.rm = TRUE))
      seTE.fixed <- sqrt((1 / (2 * sum(A, na.rm = TRUE)^2)  * 
                            (sum(A * B, na.rm = TRUE) +
                               exp(TE.fixed) * (sum(B * C, na.rm = TRUE) +
                                                  sum(A * D, na.rm = TRUE)) +
                                                    exp(TE.fixed)^2 * sum(C * D, na.rm = TRUE))))
    }
    else if (sm == "RR") {
      ##
      ## Greenland, Robins (1985) (MH.exact == TRUE)
      ##
      ## Bradburn, Deeks, Altman, Stata-procedure "metan"
      ## (MH.exact == FALSE)
      ##
      D <- ((n1. + 2 * incr.e) * (n2. + 2 * incr.c) * (n.1 + incr.e + incr.c) -
              (n11 + incr.e) * (n21 + incr.c) * (n.. + 2 * incr.e + 2 * incr.c)) /
                (n.. + 2 * incr.e + 2 * incr.c)^2
      R <- (n11 + incr.e) * (n2. + 2 * incr.c) / (n.. + 2 * incr.e + 2 * incr.c)
      S <- (n21 + incr.c) * (n1. + 2 * incr.e) / (n.. + 2 * incr.e + 2 * incr.c)
      ##
      w.fixed <- S
      TE.fixed <- log(sum(R, na.rm = TRUE) / sum(S, na.rm = TRUE))
      seTE.fixed <- sqrt(sum(D, na.rm = TRUE) / (sum(R, na.rm = TRUE) * 
                                                   sum(S, na.rm = TRUE)))
    }
    else if (sm == "RD") {
      ##
      ## Jon Deeks (1999) (MH.exact == TRUE)
      ##
      ## Bradburn, Deeks, Altman, Stata-procedure "metan"
      ## und RevMan 3.1 (MH.exact == FALSE)
      ## 
      R <- ((n11 + incr.e) * (n12 + incr.e) * (n2. + 2 * incr.c)^3 +
              (n21 + incr.c) * (n22 + incr.c) * (n1. + 2 * incr.e)^3) /
                ((n1. + 2 * incr.e) * (n2. + 2 * incr.c) * (n.. + 2 * incr.e + 2 * incr.c)^2)
      ##
      S <- n1. * n2. / n..
      ##
      w.fixed <- S
      TE.fixed <- weighted.mean(TE, w.fixed, na.rm = TRUE)
      seTE.fixed <- sqrt(sum(R, na.rm = TRUE) / sum(S, na.rm = TRUE)^2)
    }
    ##
    w.fixed[is.na(w.fixed)] <- 0
  }
  else if (method == "Peto") {
    w.fixed <- 1 / seTE^2
    TE.fixed   <- weighted.mean(TE, w.fixed, na.rm = TRUE)
    seTE.fixed <- sqrt(1 / sum(w.fixed, na.rm = TRUE))
    ##
    w.fixed[is.na(w.fixed)] <- 0
  }
  else if (method == "GLMM") {
    glmm.fixed <- metafor::rma.glmm(ai = event.e, n1i = n.e,
                                    ci = event.c, n2i = n.c,
                                    method = "FE",
                                    tdist = hakn, level = 100 * level.comb,
                                    measure = "OR", model = model.glmm,
                                    ...)
    ##
    TE.fixed   <- as.numeric(glmm.fixed$b)
    seTE.fixed <- as.numeric(glmm.fixed$se)
    ##
    w.fixed <- rep(NA, length(event.e))
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
  res <- list(event.e = event.e, n.e = n.e,
              event.c = event.c, n.c = n.c,
              method = method,
              incr = incr, sparse = sparse,
              allincr = allincr, addincr = addincr,
              allstudies = allstudies,
              doublezeros = doublezeros,
              MH.exact = MH.exact, RR.cochrane = RR.cochrane,
              Q.CMH = Q.CMH, print.CMH = print.CMH,
              incr.e = incr.e, incr.c = incr.c)
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
  res$TE.tau <- TE.tau
  res$call <- match.call()
  ##
  if (method %in% c("MH", "Peto", "GLMM")) {
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
    glmm.random <- metafor::rma.glmm(ai = event.e, n1i = n.e,
                                     ci = event.c, n2i = n.c,
                                     method = method.tau,
                                     tdist = hakn, level = 100 * level.comb,
                                     measure = "OR", model = model.glmm,
                                     ...)
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
    ##
    res$model.glmm <- model.glmm
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
    res$time.e.w <- NULL
    res$time.c.w <- NULL
  }
  ##
  class(res) <- c(fun, "meta")
  
  
  res
}
