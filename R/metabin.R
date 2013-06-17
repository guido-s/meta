metabin <- function(event.e, n.e, event.c, n.c, studlab,
                    data=NULL, subset=NULL, method="MH",
                    sm=
                    ifelse(!is.na(charmatch(method, c("Peto", "peto"),
                                            nomatch = NA)),
                           "OR", "RR"),
                    incr=0.5, allincr=FALSE, addincr=FALSE, allstudies=FALSE,
                    MH.exact=FALSE, RR.cochrane=FALSE,
                    level=0.95, level.comb=level,
                    comb.fixed=TRUE, comb.random=TRUE,
                    hakn=FALSE,
                    method.tau="DL", tau.preset=NULL, TE.tau=NULL,
                    tau.common=FALSE,
                    prediction=FALSE, level.predict=level,
                    method.bias=NULL,
                    title="", complab="", outclab="",
                    label.e="Experimental", label.c="Control",
                    label.left="", label.right="",
                    byvar, bylab, print.byvar=TRUE,
                    print.CMH=FALSE,
                    warn=TRUE
                    ){
  
  
  if (is.null(data)) data <- sys.frame(sys.parent())
  ##
  ## Catch event.e, n.e, event.c, n.e, studlab (possibly), byvar (possibly) from data:
  ##
  mf <- match.call()
  mf$data <- mf$subset <- mf$method <- mf$sm <- NULL
  mf$incr <- mf$allincr <- mf$addincr <- mf$allstudies <- NULL
  mf$MH.exact <- mf$RR.cochrane <- NULL
  mf$level <- mf$level.comb <- mf$level.predict <- mf$prediction <- mf$warn <- NULL
  mf$hakn <- mf$method.tau <- mf$tau.preset <- mf$TE.tau <- mf$method.bias <- NULL
  mf[[1]] <- as.name("data.frame")
  mf <- eval(mf, data)
  ##
  ## Catch subset (possibly) from data:
  ##
  mf2 <- match.call()
  mf2$event.e <- mf2$n.e <- NULL
  mf2$event.c <- mf2$n.c <- NULL
  mf2$studlab <- NULL 
  mf2$data <- mf2$method <- mf2$sm <- NULL
  mf2$incr <- mf2$allincr <- mf2$addincr <- mf2$allstudies <- NULL
  mf2$MH.exact <- mf2$RR.cochrane <- NULL
  mf2$level <- mf2$level.comb <- mf2$level.predict <- mf2$prediction <- mf2$warn <- NULL
  mf2$hakn <- mf2$method.tau <- mf2$tau.preset <- mf2$TE.tau <- mf2$method.bias <- NULL
  mf2$byvar <- NULL
  mf2[[1]] <- as.name("data.frame")
  ##
  mf2 <- eval(mf2, data)
  ##
  if (!is.null(mf2$subset))
    if ((is.logical(mf2$subset) & (sum(mf2$subset) > length(mf$event.e))) ||
        (length(mf2$subset) > length(mf$event.e)))
      stop("Length of subset is larger than number of studies.")
    else
      mf <- mf[mf2$subset,]
  ##
  event.e <- mf$event.e
  n.e     <- mf$n.e
  event.c <- mf$event.c
  n.c     <- mf$n.c
  ##
  missing.byvar <- missing(byvar)
  if (!missing.byvar){
    byvar.name <- deparse(substitute(byvar))
    byvar <- mf$byvar
  }
  ##
  if (!missing(studlab))
    studlab <- as.character(mf$studlab)
  else
    studlab <- row.names(mf)
  
  
  k.all <- length(event.e)
  ##
  if (k.all == 0) stop("No studies to combine in meta-analysis.")


  if (match(sm, c("OR", "RD", "RR", "AS"), nomatch=0) == 0)
    stop("possible summary measures are \"OR\", \"RD\", \"RR\", and \"AS\"")
  ##
  if (!(is.numeric(event.e) & is.numeric(n.e) &
        is.numeric(event.c) & is.numeric(n.c)))
    stop("Non-numeric value for event.e, n.e, event.c or n.c")
  ##
  if ((any(n.e <= 0 | n.c <= 0)) & warn)
    warning("Studies with non-positive values for n.e and/or n.c get no weight in meta-analysis")
  ##
  if ((any(event.e < 0 | event.c < 0)) & warn)
    warning("Studies with negative values for event.e and/or event.c get no weight in meta-analysis")
  ##
  if ((any(event.e > n.e)) & warn)
    warning("Studies with event.e > n.e get no weight in meta-analysis")
  ##
  if ((any(event.c > n.c)) & warn)
    warning("Studies with event.c > n.c get no weight in meta-analysis")
  ##
  if (length(studlab) != k.all)
    stop("Number of studies and labels are different")
  ##
  if (!is.numeric(incr)){
    ##
    iincr <- charmatch(tolower(incr),
                       c("tacc"), nomatch = NA)
    ##
    if(is.na(iincr))
      stop("incr should be numeric or the character string \"TACC\"")
    ##
    incr <- c("TACC")[iincr]
  }
  
  
  ##
  ## Check for levels of confidence interval
  ##
  if (!is.numeric(level) | length(level)!=1)
    stop("parameter 'level' must be a numeric of length 1")
  if (level <= 0 | level >= 1)
    stop("parameter 'level': no valid level for confidence interval")
  ##
  if (!is.numeric(level.comb) | length(level.comb)!=1)
    stop("parameter 'level.comb' must be a numeric of length 1")
  if (level.comb <= 0 | level.comb >= 1)
    stop("parameter 'level.comb': no valid level for confidence interval")
  ##
  if (!is.numeric(level.predict) | length(level.predict)!=1)
    stop("parameter 'level.predict' must be a numeric of length 1")
  if (level.predict <= 0 | level.predict >= 1)
    stop("parameter 'level.predict': no valid level for confidence interval")
  
  
  if (sm == "AS") method <- "Inverse"
  
  
  imeth <- charmatch(tolower(method),
                     c("inverse", "mh", "peto"), nomatch = NA)
  ##
  if(is.na(imeth))
    stop("method should be \"Inverse\", \"MH\", or \"Peto\"")
  ##
  method <- c("Inverse", "MH", "Peto")[imeth]
  ##
  if (method == "Peto" & sm != "OR")
    stop("Peto's method only possible with \"sm=OR\"")
  ##
  if (k.all == 1 & method == "MH"){
    if (warn)
      warning("For a single study, inverse variance method used instead of Mantel Haenszel method.")
    method <- "Inverse"
  }
  
  
  if (method == "MH" & !is.logical(MH.exact))
    stop("MH.exact must be of type 'logical'")
  
  
  ##
  ## Set default for method.bias (see Sterne et al., BMJ 2011;343:d4002):
  ## - "score" for OR as effect measure
  ## - "linreg" otherwise
  ##
  if (missing(method.bias) & sm=="OR")
    method.bias <- "score"
  else
    method.bias <- "linreg"
  
  
  ##
  ## Recode integer as numeric:
  ##
  if (is.integer(event.e)) event.e <- as.numeric(event.e)
  if (is.integer(n.e))     n.e     <- as.numeric(n.e)
  if (is.integer(event.c)) event.c <- as.numeric(event.c)
  if (is.integer(n.c))     n.c     <- as.numeric(n.c)
  
  
  ##
  ##
  ## Include non-informative studies?
  ## (i.e. studies with either zero or all events in both groups)
  ##
  ##
  if (sm == "RD" | sm == "AS")
    incl <- rep(1, k.all)
  else{
    allevents <- event.c==n.c & event.e==n.e
    if (allstudies)
      incl <- rep(1, k.all)
    else{
      if (sm == "OR")
        incl <- ifelse((event.c==0 & event.e==0) |
                       (event.c==n.c & event.e==n.e), NA, 1)
      if (sm == "RR")
        incl <- ifelse((event.c==0 & event.e==0), NA, 1)
    }
  }
  ##
  ## Exclude studies from meta-analysis:
  ##
  incl[event.e > n.e | event.c > n.c] <- NA
  incl[n.e <= 0      | n.c <= 0     ] <- NA
  incl[event.e < 0   | event.c < 0  ] <- NA
  
  ##
  ## k: effective number of studies
  ##
  k <- sum(!is.na(incl))
  
  
  ##
  ##
  ## Sparse computation
  ##
  ##
  sel <- switch(sm,
                OR=((n.e - event.e) == 0 | event.e == 0 |
                    (n.c - event.c) == 0 | event.c == 0),
                RD=((n.e - event.e) == 0 | event.e == 0 |
                    (n.c - event.c) == 0 | event.c == 0),
                RR=((n.e - event.e) == 0 | event.e == 0 |
                    (n.c - event.c) == 0 | event.c == 0),
                AS=rep(FALSE, length(event.e)))
  ##
  sel[is.na(incl)] <- FALSE
  ##
  sparse <- any(sel)
  ##
  ## No need to add anything to cell counts for
  ##  (i)  arcsine difference
  ##  (ii) Peto method
  ## as summary measure.
  ## Accordingly, no warning will be printed.
  ##
  warn2 <- !(sm == "AS" | method == "Peto")
  ##
  if (addincr){
    ##
    if (is.numeric(incr)){
      incr.e <- rep(incr, k.all)
      incr.c <- rep(incr, k.all)
      ##
      if (warn & warn2 & incr > 0)
        warning(paste("Increment", incr, "added to each cell frequency of all studies"))
    }
    else{
      if (incr=="TACC"){
        ##
        ## Treatment arm continuity correction:
        ##
        incr.e <- n.e/(n.e+n.c)
        incr.c <- n.c/(n.e+n.c)
        ##
        if (warn & warn2)
          warning("Treatment arm continuity correction applied to all studies")
      }
    }
    ##
  }
  else{
    if (sparse){
      if (allincr){
        ##
        if (is.numeric(incr)){
          incr.e <- rep(incr, k.all)
          incr.c <- rep(incr, k.all)
          ##
          if (warn & warn2 & incr > 0)
            warning(paste("Increment", incr, "added to each cell frequency of all studies"))
        }
        else{
          if (incr=="TACC"){
            ##
            ## Treatment arm continuity correction:
            ##
            incr.e <- n.e/(n.e+n.c)
            incr.c <- n.c/(n.e+n.c)
            ##
            if (warn & warn2)
              warning("Treatment arm continuity correction applied to all studies")
          }
        }
      }
      else{
        ##
        ## Bradburn, Deeks, Altman, Stata-procedure "metan":
        ## & SAS PROC FREQ (for method="Inverse")
        ##
       ##
        if (is.numeric(incr)){
          incr.e <- incr*sel
          incr.c <- incr*sel
          ##
           if (warn & warn2 & incr > 0)
            warning(paste("Increment", incr, "added to each cell in 2x2 tables with zero cell frequencies"))
        }
        else{
          if (incr=="TACC"){
            ##
            ## Treatment arm continuity correction:
            ##
            incr.e <- n.e/(n.e+n.c)*sel
            incr.c <- n.c/(n.e+n.c)*sel
            ##
            if (warn & warn2)
              warning("Treatment arm continuity correction applied to studies with zero cell frequencies")
          }
        }
      }
    }
    else{
      incr.e <- rep(0, k.all)
      incr.c <- rep(0, k.all)
    }
  }
  
  
  n11 <- event.e*incl
  n21 <- event.c*incl
  n1. <- n.e*incl
  n2. <- n.c*incl
  ##
  n.. <- n1. + n2.
  n12 <- n1. - n11
  n22 <- n2. - n21
  n.1 <- n11 + n21
  n.2 <- n12 + n22
  
  
  Q.CMH <- (sum(n11 - n1.*n.1/n.., na.rm=TRUE)^2/
            sum(n1.*n2.*n.1*n.2/n..^3, na.rm=TRUE))
  
  
  ##
  ##
  ## Estimation of treatment effects in individual studies
  ##
  ##
  if (sm == "OR"){
    if (method == "MH" || method == "Inverse"){
      ## 
      ## Cooper & Hedges (1994), p. 251-2
      ## 
      TE <- log(((n11+incr.e)*(n22+incr.c)) /
                ((n12+incr.e)*(n21+incr.c)))
      seTE <- sqrt((1/(n11+incr.e) + 1/(n12+incr.e) +
                    1/(n21+incr.c) + 1/(n22+incr.c)))
    }
    else if (method == "Peto"){
      ## 
      ## Cooper & Hedges (1994), p. 252
      ## 
      O <- n11
      E <- n1.*n.1/n..
      V <- n1.*n2.*n.1*n.2/((n..-1)*n..^2)
      ##
      TE <- (O-E)/V
      seTE <- sqrt(1/V)
    }
  }
  else if (sm == "RR"){
    ## 
    ## Cooper & Hedges (1994), p. 247-8
    ##
    if (!RR.cochrane){
      TE <- log(((n11+incr.e)/(n1.+incr.e))/
                ((n21+incr.c)/(n2.+incr.c)))
      ##
      ## Hartung & Knapp (2001), Stat Med, equation (18)
      ##
      seTE <- sqrt((1/(n11+incr.e*(!allevents)) - 1/(n1.+incr.e) +
                    1/(n21+incr.c*(!allevents)) - 1/(n2.+incr.c)))
      }
    else{
      TE <- log(((n11+incr.e)/(n1.+2*incr.e))/
                ((n21+incr.c)/(n2.+2*incr.c)))
      seTE <- sqrt((1/(n11+incr.e) - 1/(n1.+2*incr.e) +
                    1/(n21+incr.c) - 1/(n2.+2*incr.c)))
    }
  }
  else if (sm == "RD"){
    ## 
    ## Cooper & Hedges (1994), p. 246-7
    ## 
    TE <- n11/n1. - n21/n2.
    seTE <- sqrt((n11+incr.e)*(n12+incr.e)/(n1.+2*incr.e)^3 +
                 (n21+incr.c)*(n22+incr.c)/(n2.+2*incr.c)^3)
  }
  else if (sm == "AS"){
    ## 
    ## Ruecker et al. (2009)
    ## 
    TE <- asin(sqrt(n11/n1.)) - asin(sqrt(n21/n2.))
    seTE <- sqrt(0.25*(1/n1. + 1/n2.))
  }
  
  
  ##
  ##
  ## Calculate fixed effect and random effects estimate
  ## for inverse variance and Peto method
  ##
  ##
  if (method == "Inverse" || method== "Peto"){
    if (method == "Inverse"){
      ## Subgroup analysis with equal tau^2:
      if (!missing.byvar & tau.common){
        if (!is.null(tau.preset))
          warning("Value for argument 'tau.preset' not considered as argument 'tau.common=TRUE'")
        ##
        sm1 <- summary(metagen(TE, seTE, byvar=byvar,
                               method.tau=method.tau,
                               tau.common=tau.common))
        sQ.w <- sum(sm1$Q.w)
        sk.w <- sum(sm1$k.w-1)
        sC.w <- sum(sm1$C.w)
        ##
        if (round(sQ.w, digits=18)<=sk.w) tau2 <- 0
        else tau2 <- (sQ.w-sk.w)/sC.w
        tau.preset <- sqrt(tau2)
      }
      ##
      if (!is.null(tau.preset))
        m <- metagen(TE, seTE,
                     hakn=hakn, method.tau=method.tau,
                     tau.preset=tau.preset, TE.tau=TE.tau,
                     level=level,
                     level.comb=level.comb,
                     prediction=prediction,
                     level.predict=level.predict)
      else
        m <- metagen(TE, seTE,
                     hakn=hakn, method.tau=method.tau,
                     TE.tau=TE.tau,
                     level=level,
                     level.comb=level.comb,
                     prediction=prediction,
                     level.predict=level.predict)
    }
    ##
    if (method == "Peto"){
      if (method.tau!="DL"){
        if (warn)
          warning("DerSimonian-Laird method used to estimate between-study variance for Peto method")
        method.tau <- "DL"
      }
      ##
      if (hakn){
        if (warn)
          warning("Hartung-Knapp method not available for Peto method")
        hakn <- FALSE
      }
      ##
      if (!is.null(TE.tau)){
        if (warn)
          warning("Argument 'TE.tau' not considered for Peto method")
        TE.tau <- NULL
      }
      ##
      if (tau.common){
        if (warn)
          warning("Argument 'tau.common' not considered for Peto method")
        tau.common <- FALSE
      }
      ##
      if (!is.null(tau.preset)){
        if (warn)
          warning("Argument 'tau.preset' not considered for Peto method")
        tau.preset <- NULL
      }
      ##
      m <- metagen(TE, seTE,
                   level=level,
                   level.comb=level.comb,
                   prediction=prediction,
                   level.predict=level.predict)
    }
    ##
    TE.fixed <- m$TE.fixed
    seTE.fixed <- m$seTE.fixed
    w.fixed <- m$w.fixed
    zval.fixed <- m$zval.fixed
    pval.fixed <- m$pval.fixed
    lower.fixed <- m$lower.fixed
    upper.fixed <- m$upper.fixed
    ##
    TE.random <- m$TE.random
    seTE.random <- m$seTE.random
    w.random <- m$w.random
    zval.random <- m$zval.random
    pval.random <- m$pval.random
    lower.random <- m$lower.random
    upper.random <- m$upper.random
    ##
    df.hakn <- m$df.hakn
    ##
    Q <- m$Q
    tau2 <- m$tau^2
    se.tau2 <- m$se.tau2
    Cval <- m$C
  }
  
  
  ##
  ##
  ## Calculate fixed effect and random effects estimate
  ## for Mantel Haenszel method
  ##
  ##
  if (method == "MH"){
    ##
    if (method.tau!="DL"){
      if (warn)
        warning("DerSimonian-Laird method used to estimate between-study variance for Mantel Haenszel method")
      method.tau <- "DL"
    }
    ##
    if (hakn){
      if (warn)
        warning("Hartung-Knapp method not available for Mantel Haenszel method")
      hakn <- FALSE
    }
    ##
    if (!missing.byvar & tau.common)
      if (warn)
        warning("Argument 'tau.common' not considered for Mantel Haenszel method")
    ##
    if (!is.null(TE.tau)){
      if (warn)
        warning("Argument 'TE.tau' not considered for Mantel Haenszel method")
      TE.tau <- NULL
    }
    ##
    if (!is.null(tau.preset)){
      if (warn)
        warning("Argument 'tau.preset' not considered for Mantel Haenszel method")
      tau.preset <- NULL
    }
    ##
    incr.e <- incr.e*(!MH.exact)
    incr.c <- incr.c*(!MH.exact)
    ##
    if (sm == "OR"){
      ## 
      ## Cooper & Hedges (1994), p. 253-5 (MH.exact==TRUE)
      ##
      ## Bradburn, Deeks, Altman, Stata-procedure "metan"
      ## und RevMan 3.1 (MH.exact==FALSE)
      ## 
      A <- (n11+incr.e)*(n22+incr.c)/(n..+2*incr.e+2*incr.c)
      B <- (n11+incr.e + n22+incr.c)/(n..+2*incr.e+2*incr.c)
      C <- (n12+incr.e)*(n21+incr.c)/(n..+2*incr.e+2*incr.c)
      D <- (n12+incr.e + n21+incr.c)/(n..+2*incr.e+2*incr.c)
      ##
      ## Cooper & Hedges (1994), p. 265-6
      ##
      w.fixed <- C
      TE.fixed <- log(sum(A, na.rm=TRUE)/sum(C, na.rm=TRUE))
      seTE.fixed <- sqrt((1/(2*sum(A, na.rm=TRUE)^2) *
                          (sum(A*B, na.rm=TRUE) +
                           exp(TE.fixed)*(sum(B*C, na.rm=TRUE)+
                                          sum(A*D, na.rm=TRUE)) +
                           exp(TE.fixed)^2*sum(C*D, na.rm=TRUE))))
    }
    else if (sm =="RR"){
      ##
      ## Greenland, Robins (1985) (MH.exact==TRUE)
      ##
      ## Bradburn, Deeks, Altman, Stata-procedure "metan"
      ## (MH.exact==FALSE)
      ##
      D <- ((n1.+2*incr.e)*(n2.+2*incr.c)*(n.1+incr.e+incr.c) -
            (n11+incr.e)*(n21+incr.c)*(n..+2*incr.e+2*incr.c))/
              (n..+2*incr.e+2*incr.c)^2
      R <- (n11+incr.e)*(n2.+2*incr.c)/(n..+2*incr.e+2*incr.c)
      S <- (n21+incr.c)*(n1.+2*incr.e)/(n..+2*incr.e+2*incr.c)
      ##
      w.fixed <- S
      TE.fixed <- log(sum(R, na.rm=TRUE)/sum(S, na.rm=TRUE))
      seTE.fixed <- sqrt(sum(D, na.rm=TRUE)/(sum(R, na.rm=TRUE)*
                                  sum(S, na.rm=TRUE)))
    }
    else if (sm == "RD"){
      ##
      ## Jon Deeks (1999) (MH.exact==TRUE)
      ##
      ## Bradburn, Deeks, Altman, Stata-procedure "metan"
      ## und RevMan 3.1 (MH.exact==FALSE)
      ## 
      R <- ((n11+incr.e)*(n12+incr.e)*(n2.+2*incr.c)^3 +
            (n21+incr.c)*(n22+incr.c)*(n1.+2*incr.e)^3)/
              ((n1.+2*incr.e)*(n2.+2*incr.c)*(n..+2*incr.e+2*incr.c)^2)
      ##
      S <- n1.*n2./n..
      ##
      w.fixed <- S
      TE.fixed <- weighted.mean(TE, w.fixed, na.rm=TRUE)
      seTE.fixed <- sqrt(sum(R, na.rm=TRUE)/sum(S, na.rm=TRUE)^2)
    }
    ##
    ci.f <- ci(TE.fixed, seTE.fixed, level=level.comb)
    zval.fixed <- ci.f$z
    pval.fixed <- ci.f$p
    lower.fixed <- ci.f$lower
    upper.fixed <- ci.f$upper
    
    m <- metagen(TE, seTE, method.tau="DL", TE.tau=TE.fixed,
                 level=level,
                 level.comb=level.comb,
                 prediction=prediction,
                 level.predict=level.predict)
    ##
    TE.random <- m$TE.random
    seTE.random <- m$seTE.random
    w.random <- m$w.random
    ##
    zval.random <- m$zval.random
    pval.random <- m$pval.random
    lower.random <- m$lower.random
    upper.random <- m$upper.random
    ##
    Q <- m$Q
    se.tau2 <- m$se.tau2
    tau2 <- m$tau^2
    Cval <- m$C
  }
  
  
  ##
  ## TODO - Is this check still necessary (version of meta >= 2.0) ???
  ## Modify fixed effects estimate:
  ##
  if (is.nan(TE.fixed)){
    TE.fixed <- NA
    seTE.fixed <- NA
  }
  ##
  w.fixed[is.nan(w.fixed)] <- 0
  w.fixed[is.na(w.fixed)] <- 0
  
  
  if (m$k>=3){
    seTE.predict <- sqrt(m$seTE.random^2 + m$tau^2)
    ci.p <- ci(m$TE.random, seTE.predict, level.predict, m$k-2)
    p.lower <- ci.p$lower
    p.upper <- ci.p$upper
  }
  else{
    seTE.predict <- NA
    p.lower <- NA
    p.upper <- NA
  }
  
  
  if (!missing.byvar & tau.common)
    tau.preset <- NULL
  
  
  res <- list(event.e=event.e, n.e=n.e,
              event.c=event.c, n.c=n.c,
              studlab=studlab,
              TE=TE, seTE=seTE,
              w.fixed=w.fixed,
              w.random=w.random,
              ##
              TE.fixed=TE.fixed, seTE.fixed=seTE.fixed,
              lower.fixed=lower.fixed, upper.fixed=upper.fixed,
              zval.fixed=zval.fixed, pval.fixed=pval.fixed,
              TE.random=TE.random, seTE.random=seTE.random,
              lower.random=lower.random, upper.random=upper.random,
              zval.random=zval.random, pval.random=pval.random,
              ##
              seTE.predict=seTE.predict,
              lower.predict=p.lower, upper.predict=p.upper,
              level.predict=level.predict,
              ##
              k=k, Q=Q, tau=sqrt(tau2), se.tau2=se.tau2,
              C=Cval,
              Q.CMH=Q.CMH,
              sm=sm, method=method,
              sparse=sparse,
              incr=incr,
              allincr=allincr,
              addincr=addincr,
              allstudies=allstudies,
              MH.exact=MH.exact,
              RR.cochrane=RR.cochrane,
              incr.e=incr.e,
              incr.c=incr.c,
              level=level,
              level.comb=level.comb,
              comb.fixed=comb.fixed,
              comb.random=comb.random,
              hakn=hakn,
              df.hakn=if (hakn) df.hakn else NULL,
              method.tau=method.tau,
              tau.preset=tau.preset,
              TE.tau=if (!is.null(TE.tau) & method.tau=="DL") TE.tau else NULL,
              tau.common=tau.common,
              prediction=prediction,
              method.bias=method.bias,
              title=title,
              complab=complab,
              outclab=outclab,
              label.e=label.e,
              label.c=label.c,
              label.left=label.left,
              label.right=label.right,
              call=match.call(),
              warn=warn)
  ##
  if (!missing(byvar)){
    res$byvar <- byvar
    res$bylab <- if (!missing(bylab)) bylab else byvar.name
  }
  res$print.byvar <- print.byvar
  
  res$print.CMH <- print.CMH
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metabin", "meta")
  
  res
}
