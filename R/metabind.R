#' Combine and summarize meta-analysis objects
#' 
#' @description
#' This function can be used to combine meta-analysis objects and is,
#' for example, useful to summarize results of various meta-analysis
#' methods or to generate a forest plot with results of several
#' subgroup analyses.
#' 
#' @param ... Any number of meta-analysis objects or a single list
#'   with meta-analyses.
#' @param name An optional character vector providing descriptive
#'   names for the meta-analysis objects.
#' @param pooled A character string indicating whether results of a
#'   fixed effect or random effects model should be considered. Either
#'   \code{"fixed"} or \code{"random"}, can be abbreviated.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=TRUE} (default), results for \code{sm="OR"} are
#'   printed as odds ratios rather than log odds ratios, for example.
#' @param outclab Outcome label for all meta-analyis objects.
#' 
#' @details
#' This function can be used to combine any number of meta-analysis
#' objects which is useful, for example, to summarize results of
#' various meta-analysis methods or to generate a forest plot with
#' results of several subgroup analyses (see Examples).
#'
#' Individual study results are not retained with
#' \code{metabind}. This is possible using R function
#' \code{\link{metamerge}} which, however, can only be used to combine
#' results of two meta-analyses.
#' 
#' @return
#' An object of class \code{c("metabind", "meta")} with corresponding
#' \code{print}, \code{summary}, and \code{forest} functions. See
#' \code{\link{metagen}} for more information on list elements.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{metagen}}, \code{\link{forest.metabind}},
#'   \code{\link{metamerge}}
#' 
#' @examples
#' data(Fleiss1993cont)
#' 
#' # Add some (fictitious) grouping variables:
#' #
#' Fleiss1993cont$age <- c(55, 65, 55, 65, 55)
#' Fleiss1993cont$region <- c("Europe", "Europe", "Asia", "Asia", "Europe")
#' 
#' m1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'                data = Fleiss1993cont, sm = "MD")
#'
#' # Conduct two subgroup analyses
#' #
#' mu1 <- update(m1, byvar = age, bylab = "Age group")
#' mu2 <- update(m1, byvar = region, bylab = "Region")
#'
#' # Combine subgroup meta-analyses and show forest plot with subgroup
#' # results
#' #
#' mb1 <- metabind(mu1, mu2)
#' mb1
#' forest(mb1)
#'
#' # Use various estimation methods for between-study heterogeneity
#' # variance
#' #
#' m1.pm <- update(m1, method.tau = "PM")
#' m1.reml <- update(m1, method.tau = "REML")
#' m1.ml <- update(m1, method.tau = "ML")
#' m1.hs <- update(m1, method.tau = "HS")
#' m1.sj <- update(m1, method.tau = "SJ")
#' m1.he <- update(m1, method.tau = "HE")
#' m1.eb <- update(m1, method.tau = "EB")
#'
#' # Combine meta-analyses and show results
#' #
#' taus <- c("DerSimonian-Laird estimator",
#'           "Paule-Mandel estimator",
#'           "Restricted maximum-likelihood estimator",
#'           "Maximum-likelihood estimator",
#'           "Hunter-Schmidt estimator",
#'           "Sidik-Jonkman estimator",
#'           "Hedges estimator",
#'           "Empirical Bayes estimator")
#' #
#' m1.taus <- metabind(m1, m1.pm, m1.reml, m1.ml, m1.hs, m1.sj, m1.he, m1.eb,
#'                     name = taus, pooled = "random")
#' m1.taus
#' forest(m1.taus, print.I2 = FALSE, print.pval.Q = FALSE)
#' 
#' @export metabind


metabind <- function(..., name, pooled, backtransf, outclab) {
  
  
  if (!missing(pooled))
    pooled <- setchar(pooled, c("fixed", "random"))
  ##
  if (!missing(backtransf))
    chklogical(backtransf)


  args <- list(...)
  ##
  n.meta <- length(args)
  n.i <- seq_len(n.meta)
  ##
  if (length(args) == 1) {
    if (inherits(args[[1]], "meta.rm5")) {
      args <- args[[1]]
      if (missing(name))
        name <- unlist(lapply(args, "[[" , "outclab"))
    }
    else if (!is.list(args[[1]]))
      stop("All elements of argument '...' must be of ",
           "class 'meta'.",
           call. = FALSE)
    ##
    if (!inherits(args[[1]], "meta")) {
      n.meta <- length(args[[1]])
      n.i <- seq_len(n.meta)
      ##
      args2 <- list()
      for (i in n.i)
        args2[[i]] <- args[[1]][[i]]
      ##
      args <- args2
    }
  }
  ##  
  for (i in n.i) {
    if (!inherits(args[[i]], "meta"))
      stop("All elements of argument '...' must be of ",
           "class 'meta'.",
           call. = FALSE)
    ##
    if (inherits(args[[i]], "metabind"))
      stop("Elements of argument '...' may not be of ",
           "class 'metabind'.",
           call. = FALSE)
  }
  
  
  is.subgroup <- rep(FALSE, n.meta)
  ##
  for (i in n.i) {
    if (!inherits(args[[i]], "meta"))
      stop("All elements of argument '...' must be of class 'meta'.",
           call. = FALSE)
    ##
    if (!is.null(args[[i]]$byvar))
      is.subgroup[i] <- TRUE
  }
  
  
  print.warning1 <- FALSE
  print.warning2 <- FALSE
  print.warning3 <- FALSE


  ##
  ## Name of meta-analysis object
  ##
  if (missing(name)) {
    name <- rep("", n.meta)
    ##
    for (i in n.i)
      name[i] <- replaceNULL(args[[i]]$bylab)
    ##
    if (all(is.na(name)))
      name <- paste0("meta", n.i)
    else if (anyNA(name))
      name[is.na(name)] <- paste0("meta", n.i[is.na(name)])
    ##
    if (length(unique(name)) != length(name)) {
      warning1 <-
        paste0("Meta-analyses are labelled 'meta1' to 'meta", n.meta,
               "' as argument 'name' is missing and information ",
               "from meta-analysis objects is not unique.")
      print.warning1 <- TRUE
      ##
      name <- paste0("meta", n.i)
    }
  }
  else {
    if (length(name) != length(is.subgroup))
      stop("Number of meta-analyses and names provided in ",
           "argument 'name' differ.",
           call. = FALSE)
    ##
    if (length(unique(name)) != length(name)) {
      warning2 <-
        paste0("Meta-analyses are labelled 'meta1' to 'meta", n.meta,
               "' as values of argument 'name' are not all disparate.")
      print.warning2 <- FALSE
      ##
      name <- paste0("meta", n.i)
    }
  }
  
  
  for (i in n.i) {
    m.i <- args[[i]]
    ##
    meth.i <- data.frame(sm = m.i$sm,
                         method = m.i$method,
                         level = m.i$level.comb,
                         level.comb = m.i$level.comb,
                         level.predict = m.i$level.predict,
                         comb.fixed = m.i$comb.fixed,
                         comb.random = m.i$comb.random,
                         hakn = m.i$hakn,
                         method.tau = m.i$method.tau,
                         tau.preset = replaceNULL(m.i$tau.preset),
                         TE.tau = replaceNULL(m.i$TE.tau),
                         tau.common = m.i$tau.common,
                         prediction = m.i$prediction,
                         method.bias = "",
                         null.effect = m.i$null.effect,
                         ##
                         title = m.i$title,
                         complab = m.i$complab,
                         outclab = if (missing(outclab)) m.i$outclab else outclab,
                         label.e = m.i$label.e,
                         label.c = m.i$label.c,
                         label.left = m.i$label.left,
                         label.right = m.i$label.right,
                         ##
                         print.byvar = FALSE,
                         byseparator = "",
                         warn = replaceNULL(m.i$warn, FALSE),
                         ##
                         backtransf = m.i$backtransf,
                         pscale = m.i$pscale,
                         irscale = m.i$irscale,
                         irunit = replaceNULL(m.i$ir.unit),
                         ##
                         stringsAsFactors = FALSE)
    ##
    if (i == 1)
      meth <- meth.i
    else
      meth <- rbind(meth, meth.i)
  }
  ##
  ## Unify some settings
  ##
  if (missing(pooled)) {
    if (all(meth$comb.fixed) & all(!meth$comb.random))
      pooled <- "fixed"
    else if (all(!meth$comb.fixed) & all(meth$comb.random))
      pooled <- "random"
    else {
      if (any(meth$comb.fixed)) {
        warning3 <-
          paste("Note, results from random effects model extracted.",
                "Use argument pooled = \"fixed\" for results of",
                "fixed effect model.")
        print.warning3 <- TRUE
      }
      pooled <- "random"
    }
  }
  ##
  if (pooled == "fixed") {
    meth$comb.fixed <- TRUE
    meth$comb.random <- FALSE
  }
  else {
    meth$comb.fixed <- FALSE
    meth$comb.random <- TRUE
  }
  
  
  for (i in n.i) {
    m.i <- args[[i]]
    ##
    if (length(m.i$tau) > 1)
      if (pooled == "random") {
        m.i$tau <- m.i$tau[2]
        m.i$tau2 <- m.i$tau2[2]
        ##
        if (length(m.i$lower.tau) > 1) {
          m.i$lower.tau <- m.i$lower.tau[2]
          m.i$upper.tau <- m.i$upper.tau[2]
          m.i$lower.tau2 <- m.i$lower.tau2[2]
          m.i$upper.tau2 <- m.i$upper.tau2[2]
        }
      }
      else {
        m.i$tau <- m.i$tau[1]
        m.i$tau2 <- m.i$tau2[1]
        ##
        if (length(m.i$lower.tau) > 1) {
          m.i$lower.tau <- m.i$lower.tau[1]
          m.i$upper.tau <- m.i$upper.tau[1]
          m.i$lower.tau2 <- m.i$lower.tau2[1]
          m.i$upper.tau2 <- m.i$upper.tau2[1]
        }
      } 
    ##
    subgroup.i <- data.frame(TE.fixed.w = m.i$TE.fixed,
                             seTE.fixed.w = m.i$seTE.fixed,
                             lower.fixed.w = m.i$lower.fixed,
                             upper.fixed.w = m.i$upper.fixed,
                             statistic.fixed.w = m.i$statistic.fixed,
                             pval.fixed.w = m.i$pval.fixed,
                             w.fixed.w = 0, # sum(m.i$w.fixed),
                             ##
                             TE.random.w = m.i$TE.random,
                             seTE.random.w = m.i$seTE.random,
                             lower.random.w = m.i$lower.random,
                             upper.random.w = m.i$upper.random,
                             statistic.random.w = m.i$statistic.random,
                             pval.random.w = m.i$pval.random,
                             df.hakn.w = replaceNULL(m.i$df.hakn),
                             w.random.w = 0, # sum(m.i$w.random),
                             ##
                             n.harmonic.mean.w =
                               1 / mean(1 / replaceNULL(m.i$n)),
                             t.harmonic.mean.w =
                               1 / mean(1 / replaceNULL(m.i$time)),
                             ##
                             n.e.w = sum(replaceNULL(m.i$n.e)),
                             n.c.w = sum(replaceNULL(m.i$n.c)),
                             ##
                             k.w = m.i$k,
                             k.all.w = length(m.i$TE),
                             ##
                             Q.w = m.i$Q,
                             df.Q.w = m.i$df.Q,
                             pval.Q.w = m.i$pval.Q,
                             ##
                             tau2.w = m.i$tau2,
                             se.tau2.w = m.i$se.tau2,
                             lower.tau2.w = m.i$lower.tau2,
                             upper.tau2.w = m.i$upper.tau2,
                             tau.w = m.i$tau,
                             lower.tau.w = m.i$lower.tau,
                             upper.tau.w = m.i$upper.tau,
                             H.w = m.i$H,
                             lower.H.w = m.i$lower.H,
                             upper.H.w = m.i$upper.H,
                             I2.w = m.i$I2,
                             lower.I2.w = m.i$lower.I2,
                             upper.I2.w = m.i$upper.I2,
                             Rb.w = m.i$Rb,
                             lower.Rb.w = m.i$lower.Rb,
                             upper.Rb.w = m.i$upper.Rb,
                             ##
                             stringsAsFactors = FALSE)
    ##
    if (is.subgroup[i]) {
      ##
      Q.b.fixed.i <- m.i$Q.b.fixed
      Q.b.random.i <- m.i$Q.b.random
      df.Q.b.i <- m.i$df.Q.b
      pval.Q.b.fixed.i  <- m.i$pval.Q.b.fixed
      pval.Q.b.random.i <- m.i$pval.Q.b.random
      ##
      n.bylevs.i <- length(m.i$k.w) - 1
      ##
      if (n.bylevs.i > 0) {
        Q.b.fixed.i <- c(Q.b.fixed.i, rep(NA, n.bylevs.i))
        Q.b.random.i <- c(Q.b.random.i, rep(NA, n.bylevs.i))
        df.Q.b.i <- c(df.Q.b.i, rep(NA, n.bylevs.i))
        pval.Q.b.fixed.i <- c(pval.Q.b.fixed.i, rep(NA, n.bylevs.i))
        pval.Q.b.random.i <- c(pval.Q.b.random.i, rep(NA, n.bylevs.i))
      }
      ##
      data.i <- data.frame(name = name[i],
                           bylevs = m.i$bylevs,
                           ##
                           n.e = replaceNULL(m.i$n.e.w),
                           n.c = replaceNULL(m.i$n.c.w),
                           df.hakn = replaceNULL(m.i$df.hakn.w),
                           ##
                           k = m.i$k.w,
                           k.all = m.i$k.all.w,
                           Q = m.i$Q.w,
                           df.Q = m.i$k.w - 1,
                           pval.Q = pvalQ(m.i$Q.w, m.i$k.w - 1),
                           ##
                           tau2 = m.i$tau.w^2,
                           lower.tau2 = m.i$lower.tau2.w,
                           upper.tau2 = m.i$upper.tau2.w,
                           tau = m.i$tau.w,
                           lower.tau = m.i$lower.tau.w,
                           upper.tau = m.i$upper.tau.w,
                           H = m.i$H.w,
                           lower.H = m.i$lower.H.w,
                           upper.H = m.i$upper.H.w,
                           I2 = m.i$I2.w,
                           lower.I2 = m.i$lower.I2.w,
                           upper.I2 = m.i$upper.I2.w,
                           Rb = m.i$Rb.w,
                           lower.Rb = m.i$lower.Rb.w,
                           upper.Rb = m.i$upper.Rb.w,
                           ##
                           Q.b.fixed = Q.b.fixed.i,
                           Q.b.random = Q.b.random.i,
                           df.Q.b = df.Q.b.i,
                           pval.Q.b.fixed = pval.Q.b.fixed.i,
                           pval.Q.b.random = pval.Q.b.random.i,
                           ##
                           stringsAsFactors = FALSE)
    }
    else
      data.i <- data.frame(name = name[i],
                           bylevs = "overall",
                           ##
                           n.e = sum(replaceNULL(m.i$n.e)),
                           n.c = sum(replaceNULL(m.i$n.c)),
                           df.hakn = replaceNULL(m.i$df.hakn),
                           ##
                           k = m.i$k,
                           k.all = length(m.i$TE),
                           Q = m.i$Q,
                           df.Q = m.i$df.Q,
                           pval.Q = pvalQ(m.i$Q, m.i$df.Q),
                           ##
                           tau = m.i$tau,
                           lower.tau = m.i$lower.tau,
                           upper.tau = m.i$upper.tau,
                           tau2 = m.i$tau^2,
                           lower.tau2 = m.i$lower.tau2,
                           upper.tau2 = m.i$upper.tau2,
                           H = m.i$H,
                           lower.H = m.i$lower.H,
                           upper.H = m.i$upper.H,
                           I2 = m.i$I2,
                           lower.I2 = m.i$lower.I2,
                           upper.I2 = m.i$upper.I2,
                           Rb = m.i$Rb,
                           lower.Rb = m.i$lower.Rb,
                           upper.Rb = m.i$upper.Rb,
                           ##
                           Q.b.fixed = NA,
                           Q.b.random = NA,
                           df.Q.b = NA,
                           pval.Q.b.fixed = NA,
                           pval.Q.b.random = NA,
                           ##
                           stringsAsFactors = FALSE)
    ##
    overall.i <- data.frame(name = name[i],
                            ##
                            TE.fixed = m.i$TE.fixed,
                            seTE.fixed = m.i$seTE.fixed,
                            lower.fixed = m.i$lower.fixed,
                            upper.fixed = m.i$upper.fixed,
                            statistic.fixed = m.i$statistic.fixed,
                            pval.fixed = m.i$pval.fixed,
                            ##
                            TE.random = m.i$TE.random,
                            seTE.random = m.i$seTE.random,
                            lower.random = m.i$lower.random,
                            upper.random = m.i$upper.random,
                            statistic.random = m.i$statistic.random,
                            pval.random = m.i$pval.random,
                            df.hakn = replaceNULL(m.i$df.hakn),
                            ##
                            seTE.predict = m.i$seTE.predict,
                            lower.predict = m.i$lower.predict,
                            upper.predict = m.i$upper.predict,
                            ##
                            k = m.i$k,
                            k.all = length(m.i$TE),
                            ##
                            Q = m.i$Q,
                            df.Q = m.i$df.Q,
                            tau2 = m.i$tau2,
                            lower.tau2 = m.i$lower.tau2,
                            upper.tau2 = m.i$upper.tau2,
                            se.tau2 = replaceNULL(m.i$se.tau2),
                            tau = m.i$tau,
                            lower.tau = m.i$lower.tau,
                            upper.tau = m.i$upper.tau,
                            ##
                            H = m.i$H,
                            lower.H = m.i$lower.H,
                            upper.H = m.i$upper.H,
                            I2 = m.i$I2,
                            lower.I2 = m.i$lower.I2,
                            upper.I2 = m.i$upper.I2,
                            Rb = m.i$Rb,
                            lower.Rb = m.i$lower.Rb,
                            upper.Rb = m.i$upper.Rb,
                            ##
                            Q.w.fixed = NA,
                            Q.w.random = NA,
                            ##
                            Q.b.fixed = NA,
                            pval.Q.b.fixed = NA,
                            Q.b.random = NA,
                            df.Q.b = NA,
                            pval.Q.b.random = NA,
                            ##
                            stringsAsFactors = FALSE)
    ##
    if (i == 1) {
      data <- data.i
      overall <- overall.i
      subgroup <- subgroup.i
    }
    else {
      data <- rbind(data, data.i)
      overall <- rbind(overall, overall.i)
      subgroup <- rbind(subgroup, subgroup.i)
    }
  }
  
  
  ## Unify more settings
  ##
  if (missing(backtransf)) {
    if (any(meth$backtransf))
      meth$backtransf <- TRUE
  }
  else
    meth$backtransf <- backtransf
  ##
  if (any(meth$warn))
    meth$warn <- TRUE
  ##
  if (any(meth$prediction))
    meth$prediction <- TRUE
  ##  
  ## Only consider argument 'tau.common' from subgroup meta-analyses
  ##
  if (any(is.subgroup) & any(!is.subgroup)) {
    tau.common.uniq <- unique(meth$tau.common[is.subgroup])
    if (length(tau.common.uniq) == 1)
      meth$tau.common[!is.subgroup] <- tau.common.uniq
  }
  ##
  show.studies <- TRUE
  overall.hetstat <- TRUE
  ##
  if (length(unique(meth$method)) != 1) {
    meth$method <- ""
    show.studies <- FALSE
    overall.hetstat <- FALSE
  }
  ##
  if (length(unique(meth$hakn)) != 1) {
    meth$hakn <- FALSE
  }
  ##
  if (length(unique(meth$method.tau)) != 1) {
    meth$method.tau <- ""
    show.studies <- FALSE
    overall.hetstat <- FALSE
  }
  ##
  if (length(unique(meth$method.bias)) != 1)
    meth$method.bias <- ""
  ##
  if (length(unique(meth$title)) != 1)
    meth$title <- ""
  ##
  if (length(unique(meth$complab)) != 1)
    meth$complab <- ""
  ##
  if (length(unique(meth$outclab)) != 1)
    meth$outclab <- ""
  ##
  if (length(unique(meth$label.e)) != 1)
    meth$label.e <- ""
  ##
  if (length(unique(meth$label.c)) != 1)
    meth$label.c <- ""
  ##
  if (length(unique(meth$label.left)) != 1)
    meth$label.left <- ""
  ##
  if (length(unique(meth$label.right)) != 1)
    meth$label.right <- ""
  ##
  if (length(unique(meth$pscale)) != 1)
    meth$pscale <- min(meth$pscale)
  ##
  if (length(unique(meth$irscale)) != 1)
    meth$irscale <- min(meth$irscale)
  ##
  if (length(unique(meth$irunit)) != 1)
    meth$irunit <- NA
  
  
  ## Check whether settings are unique
  ##
  meth2 <- meth
  meth2$level <- NULL
  n.meth <- apply(meth2, 2,
                  function(x)
                    length(unique(x)))
  ##
  if (any(n.meth != 1))
    stop("Setting for the following argument",
         if (sum(n.meth != 1) > 1) "s",
         " must be the same for all meta-analyses",
         ": ",
         paste0(paste0("'", names(meth2)[n.meth != 1], "'"),
                collapse = " - "))
  
  
  for (i in n.i) {
    m.i <- args[[i]]
    ##
    study.i <- data.frame(studlab = replaceNULL(m.i$bylevs, "overall"),
                          stringsAsFactors = FALSE)
    ##
    if (is.subgroup[i]) {
      study.i$n.e <- replaceNULL(m.i$n.e.w)
      study.i$n.c <- replaceNULL(m.i$n.c.w)
      ##
      if (pooled == "fixed") {
        study.i$TE <- m.i$TE.fixed.w
        study.i$seTE <- m.i$seTE.fixed.w
        study.i$lower <- m.i$lower.fixed.w
        study.i$upper <- m.i$upper.fixed.w
        study.i$statistic <- m.i$statistic.fixed.w
        study.i$pval <- m.i$pval.fixed.w
        study.i$w.fixed <- m.i$w.fixed.w
        study.i$w.random <- 0
      }
      else {
        study.i$TE <- m.i$TE.random.w
        study.i$seTE <- m.i$seTE.random.w
        study.i$lower <- m.i$lower.random.w
        study.i$upper <- m.i$upper.random.w
        study.i$statistic <- m.i$statistic.random.w
        study.i$pval <- m.i$pval.random.w
        study.i$w.fixed <- 0
        study.i$w.random <- m.i$w.random.w
      }
    }
    else {
      study.i$n.e <- sum(replaceNULL(m.i$n.e.w))
      study.i$n.c <- sum(replaceNULL(m.i$n.c.w))
      ##
      if (pooled == "fixed") {
        study.i$TE <- m.i$TE.fixed
        study.i$seTE <- m.i$seTE.fixed
        study.i$lower <- m.i$lower.fixed
        study.i$upper <- m.i$upper.fixed
        study.i$statistic <- m.i$statistic.fixed
        study.i$pval <- m.i$pval.fixed
        study.i$w.fixed <- 1
        study.i$w.random <- 0
      }
      else {
        study.i$TE <- m.i$TE.random
        study.i$seTE <- m.i$seTE.random
        study.i$lower <- m.i$lower.random
        study.i$upper <- m.i$upper.random
        study.i$statistic <- m.i$statistic.random
        study.i$pval <- m.i$pval.random
        study.i$w.fixed <- 0
        study.i$w.random <- 1
      }
    }
    ##
    study.i$byvar <- name[i]
    ##
    if (i == 1)
      study <- study.i
    else
      study <- rbind(study, study.i)
  }
  
  
  if (length(unique(study$byvar)) == 1) {
    res <- c(as.list(study), as.list(meth[1, ]), as.list(overall))
    res$byvar <- NULL
  }
  else
    res <- c(as.list(study), as.list(meth[1, ]),
             as.list(overall), as.list(subgroup))
  ##
  ##
  res$data <- data
  ##
  res$call <- match.call()
  res$version <- packageDescription("meta")$Version


  makeunique <- function(x, val = NA) {
    if (length(unique(x)) == 1)
      res <- unique(x)
    else
      res <- val
    ##
    res
  }
  ##
  res$TE.fixed <- makeunique(res$TE.fixed)
  res$seTE.fixed <- makeunique(res$seTE.fixed)
  res$lower.fixed <- makeunique(res$lower.fixed)
  res$upper.fixed <- makeunique(res$upper.fixed)
  res$statistic.fixed <- makeunique(res$statistic.fixed)
  res$pval.fixed <- makeunique(res$pval.fixed)
  ##
  res$TE.random <- makeunique(res$TE.random)
  res$seTE.random <- makeunique(res$seTE.random)
  res$lower.random <- makeunique(res$lower.random)
  res$upper.random <- makeunique(res$upper.random)
  res$statistic.random <- makeunique(res$statistic.random)
  res$pval.random <- makeunique(res$pval.random)
  res$df.hakn <- makeunique(res$df.hakn)
  ##
  res$seTE.predict <- makeunique(res$seTE.predict)
  res$lower.predict <- makeunique(res$lower.predict)
  res$upper.predict <- makeunique(res$upper.predict)
  ##
  res$k <- makeunique(res$k)
  res$k.all <- makeunique(res$k.all)
  ##
  res$Q <- makeunique(res$Q)
  res$df.Q <- makeunique(res$df.Q, 0)
  res$pval.Q <- makeunique(makeunique(res$pval.Q, pvalQ(res$Q, res$df.Q)))
  res$tau2 <- makeunique(res$tau2)
  res$se.tau2 <- makeunique(res$se.tau2)
  res$tau <- makeunique(res$tau)
  res$lower.tau <- res$upper.tau <- NA
  res$lower.tau2 <- res$upper.tau2 <- NA
  res$method.tau.ci <- ""
  ##
  res$H <- makeunique(res$H)
  res$lower.H <- makeunique(res$lower.H)
  res$upper.H <- makeunique(res$upper.H)
  ##
  res$I2 <- makeunique(res$I2)
  res$lower.I2 <- makeunique(res$lower.I2)
  res$upper.I2 <- makeunique(res$upper.I2)
  ##
  res$Rb <- makeunique(res$Rb)
  res$lower.Rb <- makeunique(res$lower.Rb)
  res$upper.Rb <- makeunique(res$upper.Rb)
  ##
  res$Q.w.fixed <- makeunique(res$Q.w.fixed)
  res$Q.w.random <- makeunique(res$Q.w.random)
  res$df.Q.w <- makeunique(res$df.Q.w, 0)
  ##
  res$Q.b.fixed <- makeunique(res$Q.b.fixed)
  res$Q.b.random <- makeunique(res$Q.b.random)
  ##
  res$df.Q.b <- makeunique(res$df.Q.b, 0)
  res$pval.Q.b.fixed <-
    makeunique(makeunique(res$pval.Q.b.fixed,
                          pvalQ(res$Q.b.fixed, res$df.Q.b)))
  res$pval.Q.b.random <-
    makeunique(makeunique(res$pval.Q.b.random,
                          pvalQ(res$Q.b.random, res$df.Q.b)))
  ##
  res$show.studies <- show.studies
  res$overall.hetstat <- overall.hetstat
  
  
  res$is.subgroup <- is.subgroup


  if (!is.null(res$byvar)) {
    res$bylab <- "meta-analysis"
    res$bylevs <- unique(res$byvar)
    res$w.fixed <- rep(0, length(res$w.fixed))
    res$w.fixed.w <- rep(0, length(res$w.fixed.w))
    res$w.random <- rep(0, length(res$w.random))
    res$w.random.w <- rep(0, length(res$w.random.w))
    res$lower.predict.w <- rep(NA, length(res$w.random.w))
    res$upper.predict.w <- rep(NA, length(res$w.random.w))
  }
  
  
  if (is.na(res$tau.preset))
    res$tau.preset <- NULL


  class(res) <- c("metabind", "meta")


  ##
  ##
  ##
  if (print.warning1)
    warning(warning1)
  if (print.warning2)
    warning(warning2)
  if (print.warning3)
    warning(warning3)


  res
}
