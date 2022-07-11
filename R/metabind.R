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
#' @param pooled A character string or vector indicating whether
#'   results of a common effect or random effects model should be
#'   considered. Either \code{"common"} or \code{"random"}, can be
#'   abbreviated.
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
#'   data = Fleiss1993cont, sm = "MD")
#'
#' # Conduct two subgroup analyses
#' #
#' mu1 <- update(m1, subgroup = age, subgroup.name = "Age group")
#' mu2 <- update(m1, subgroup = region, subgroup.name = "Region")
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
#' m1.dl <- update(m1, method.tau = "DL")
#' m1.ml <- update(m1, method.tau = "ML")
#' m1.hs <- update(m1, method.tau = "HS")
#' m1.sj <- update(m1, method.tau = "SJ")
#' m1.he <- update(m1, method.tau = "HE")
#' m1.eb <- update(m1, method.tau = "EB")
#'
#' # Combine meta-analyses and show results
#' #
#' taus <- c("Restricted maximum-likelihood estimator",
#'   "Paule-Mandel estimator",
#'   "DerSimonian-Laird estimator",
#'   "Maximum-likelihood estimator",
#'   "Hunter-Schmidt estimator",
#'   "Sidik-Jonkman estimator",
#'   "Hedges estimator",
#'   "Empirical Bayes estimator")
#' #
#' m1.taus <- metabind(m1, m1.pm, m1.dl, m1.ml, m1.hs, m1.sj, m1.he, m1.eb,
#'   name = taus, pooled = "random")
#' m1.taus
#' forest(m1.taus, print.I2 = FALSE, print.pval.Q = FALSE)
#' 
#' @export metabind


metabind <- function(..., name = NULL, pooled = NULL,
                     backtransf = NULL, outclab = NULL) {
  
  
  missing.name <- missing(name)
  missing.pooled <- missing(pooled)
  missing.backtransf <- missing(backtransf)
  missing.outclab <- missing(outclab)
  ##
  args <- list(...)
  ##
  n.meta <- length(args)
  n.i <- seq_len(n.meta)
  is.limit <- is.copas <- is.trimfill <- rep(FALSE, n.meta)
  ##
  if (!missing(pooled)) {
    pooled <- setchar(pooled, c("common", "random", "fixed"))
    pooled[pooled == "fixed"] <- "common"
  }
  ##
  if (!missing.backtransf)
    chklogical(backtransf)  
  
  
  ##
  ## Act on single meta-analysis object in '...'
  ##
  if (n.meta == 1) {
    if (inherits(args[[1]], "meta.rm5")) {
      args <- args[[1]]
      if (missing.name) {
        name <- unlist(lapply(args, "[[" , "outclab"))
        missing.name <- FALSE
      }
    }
    else if (inherits(args[[1]], c("limitmeta", "copas")))
      return(metamerge(args[[1]]))
    else if (inherits(args[[1]], "meta"))
      return(args[[1]])
    else if (!is.list(args[[1]]))
      stop("All elements of argument '...' must be of class 'meta', ",
           "'limitmeta', or 'copas'.",
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
  ## Act on limitmeta and copas objects
  ##
  name.i <- rep(NA, n.meta)
  ##
  for (i in n.i) {
    if (inherits(args[[i]], "metabind"))
      stop("Elements of argument '...' may not be of class 'metabind'.",
           call. = FALSE)
    ##
    if (inherits(args[[i]], "meta")) {
      args[[i]] <- updateversion(args[[i]])
      if (missing.name) {
        if (inherits(args[[i]], "trimfill")) {
          is.trimfill[i] <- TRUE
          name.i[i] <- "trimfill"
        }
        else
          name.i[i] <- replaceNULL(args[[i]]$subgroup.name)
        if (is.na(name.i[i]))
          name.i[i] <- class(args[[i]])[1]
      }
    }
    else if (inherits(args[[i]], c("limitmeta", "copas"))) {
      if (missing.name)
        name.i[i] <- class(args[[i]])
      if (inherits(args[[i]], "limitmeta"))
        is.limit[i] <- TRUE
      else
        is.copas[i] <- TRUE
      ##
      args[[i]] <- metamerge(args[[i]])
      args[[i]]$common <- FALSE
    }
    else
      stop("All elements of argument '...' must be of class 'meta', ",
           "'limitmeta', or 'copas'.",
           call. = FALSE)
  }
  ##
  is.limit.copas <- is.limit | is.copas
  ##
  is.subgroup <- rep(FALSE, n.meta)
  ##
  for (i in n.i) {
    if (!is.null(args[[i]]$subgroup))
      is.subgroup[i] <- TRUE
  }
  ##
  if (missing.pooled || length(pooled) == 1)
    pooled <- rep(pooled, n.meta)
  else
    chklength(pooled, n.meta,
              text = paste("Length of argument 'pooled' differs from",
                           "number of meta-analyses."))
  
  
  ##
  ## Name of meta-analysis object
  ##
  if (missing.name) {
    name <- name.i
    ##
    if (all(is.na(name)))
      name <- paste0("meta", n.i)
    else if (anyNA(name))
      name[is.na(name)] <- paste0("meta", n.i[is.na(name)])
  }
  else {
    if (length(name) != length(is.subgroup))
      stop("Number of meta-analyses and names provided in ",
           "argument 'name' differ.",
           call. = FALSE)
  }
  ##
  ## Names for meta-analyses must be unique
  ##
  if (length(unique(name)) != length(name)) {
    for (i in n.i)
      if (name[i] %in% c("metabin", "metainc", "metaprop", "metarate") &
          !is.trimfill[i])
        name[i] <- paste(name[i], args[[i]]$method, sep = ".")
  }
  ##
  if (length(unique(name)) != length(name)) {
    if (missing.pooled) {
      for (i in n.i)
        if (inherits(args[[i]], "meta") & !is.copas[i])
          name[i] <- paste(name[i], args[[i]]$method.tau, sep = ".")
    }
    else {
      for (i in n.i)
        if (inherits(args[[i]], "meta") & pooled[i] == "random" &
            !is.copas[i])
          name[i] <- paste(name[i], args[[i]]$method.tau, sep = ".")
    }
  }
  ##
  if (length(unique(name)) != length(name))
    name <- paste0("meta", n.i)
  
  
  for (i in n.i) {
    m.i <- args[[i]]
    ##
    meth.i <- data.frame(sm = m.i$sm,
                         method = m.i$method,
                         level = m.i$level.ma,
                         level.ma = m.i$level.ma,
                         level.predict = m.i$level.predict,
                         common = m.i$common,
                         random = m.i$random,
                         hakn = m.i$hakn,
                         method.tau = m.i$method.tau,
                         tau.preset = replaceNULL(m.i$tau.preset),
                         TE.tau = replaceNULL(m.i$TE.tau),
                         tau.common = replaceNULL(m.i$tau.common, FALSE),
                         prediction = m.i$prediction,
                         prediction.subgroup =
                           replaceNULL(m.i$prediction.subgroup, FALSE),
                         method.bias = "",
                         null.effect = m.i$null.effect,
                         ##
                         title = m.i$title,
                         complab = m.i$complab,
                         outclab =
                           if (missing.outclab) m.i$outclab else outclab,
                         label.e = m.i$label.e,
                         label.c = m.i$label.c,
                         label.left = m.i$label.left,
                         label.right = m.i$label.right,
                         ##
                         print.subgroup.name = FALSE,
                         sep.subgroup = "",
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
  if (missing.pooled) {
    if (all(meth$common) & all(!meth$random))
      pooled <- rep("common", n.meta)
    else
      pooled <- rep("random", n.meta)
  }
  ##
  unique.pooled <- length(unique(pooled)) == 1
  ##
  if (all(pooled == "random")) {
    meth$common <- FALSE
    meth$random <- TRUE
  }
  else {
    meth$common <- TRUE
    meth$random <- FALSE
  }
  
  
  for (i in n.i) {
    m.i <- args[[i]]
    ##
    if (length(m.i$tau) > 1)
      if (pooled[i] == "random") {
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
    if (unique.pooled) {
      sel.r <- TRUE
      sel.f <- !sel.r & !is.limit.copas[i]
    }
    else {
      sel.r <- pooled[i] == "random"
      sel.f <- !sel.r & !is.limit.copas[i]
    }
    ##
    subgroup.i <- data.frame(
      TE.common.w = if (sel.f) m.i$TE.common else m.i$TE.random,
      seTE.common.w = if (sel.f) m.i$seTE.common else m.i$seTE.random,
      lower.common.w = if (sel.f) m.i$lower.common else m.i$lower.random,
      upper.common.w = if (sel.f) m.i$upper.common else m.i$upper.random,
      statistic.common.w =
        if (sel.f) m.i$statistic.common else m.i$statistic.random,
      pval.common.w = if (sel.f) m.i$pval.common else m.i$pval.random,
      w.common.w = 0, # sum(m.i$w.common),
      ##
      TE.random.w = if (!sel.r) m.i$TE.common else m.i$TE.random,
      seTE.random.w = if (!sel.r) m.i$seTE.common else m.i$seTE.random,
      lower.random.w = if (!sel.r) m.i$lower.common else m.i$lower.random,
      upper.random.w = if (!sel.r) m.i$upper.common else m.i$upper.random,
      statistic.random.w =
        if (!sel.r) m.i$statistic.common else m.i$statistic.random,
      pval.common.w = if (!sel.r) m.i$pval.common else m.i$pval.random,
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
      k.study.w = m.i$k.study,
      k.all.w = m.i$k.all,
      k.TE.w = m.i$k.TE,
      ##
      Q.w = if (!sel.r) NA else m.i$Q,
      df.Q.w = if (!sel.r) NA else m.i$df.Q,
      pval.Q.w = if (!sel.r) NA else m.i$pval.Q,
      ##
      tau2.w = if (!sel.r) NA else m.i$tau2,
      se.tau2.w = if (!sel.r) NA else m.i$se.tau2,
      lower.tau2.w = if (!sel.r) NA else m.i$lower.tau2,
      upper.tau2.w = if (!sel.r) NA else m.i$upper.tau2,
      tau.w = if (!sel.r) NA else m.i$tau,
      lower.tau.w = if (!sel.r) NA else m.i$lower.tau,
      upper.tau.w = if (!sel.r) NA else m.i$upper.tau,
      H.w = if (!sel.r) NA else m.i$H,
      lower.H.w = if (!sel.r) NA else m.i$lower.H,
      upper.H.w = if (!sel.r) NA else m.i$upper.H,
      I2.w = if (!sel.r) NA else m.i$I2,
      lower.I2.w = if (!sel.r) NA else m.i$lower.I2,
      upper.I2.w = if (!sel.r) NA else m.i$upper.I2,
      Rb.w = if (!sel.r) NA else m.i$Rb,
      lower.Rb.w = if (!sel.r) NA else m.i$lower.Rb,
      upper.Rb.w = if (!sel.r) NA else m.i$upper.Rb,
      ##
      stringsAsFactors = FALSE)
    ##
    if (is.subgroup[i]) {
      ##
      Q.b.common.i <- m.i$Q.b.common
      Q.b.random.i <- m.i$Q.b.random
      df.Q.b.i <- m.i$df.Q.b
      pval.Q.b.common.i  <- m.i$pval.Q.b.common
      pval.Q.b.random.i <- m.i$pval.Q.b.random
      ##
      n.bylevs.i <- length(m.i$k.w) - 1
      ##
      if (n.bylevs.i > 0) {
        Q.b.common.i <- c(Q.b.common.i, rep(NA, n.bylevs.i))
        Q.b.random.i <- c(Q.b.random.i, rep(NA, n.bylevs.i))
        df.Q.b.i <- c(df.Q.b.i, rep(NA, n.bylevs.i))
        pval.Q.b.common.i <- c(pval.Q.b.common.i, rep(NA, n.bylevs.i))
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
                           n.harmonic.mean = m.i$n.harmonic.mean.w,
                           t.harmonic.mean = m.i$t.harmonic.mean.w,
                           ##
                           k = m.i$k.w,
                           k.study = m.i$k.study.w,
                           k.all = m.i$k.all.w,
                           k.TE = m.i$k.TE.w,
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
                           Q.b.common = Q.b.common.i,
                           Q.b.random = Q.b.random.i,
                           df.Q.b = df.Q.b.i,
                           pval.Q.b.common = pval.Q.b.common.i,
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
                           n.harmonic.mean =
                             1 / mean(1 / replaceNULL(m.i$n)),
                           t.harmonic.mean =
                             1 / mean(1 / replaceNULL(m.i$time)),
                           ##
                           k = m.i$k,
                           k.study = m.i$k.study,
                           k.all = m.i$k.all,
                           k.TE = m.i$k.TE,
                           Q = m.i$Q,
                           df.Q = m.i$df.Q,
                           pval.Q = pvalQ(m.i$Q, m.i$df.Q),
                           ##
                           tau = if (!sel.r) NA else m.i$tau,
                           lower.tau = if (!sel.r) NA else m.i$lower.tau,
                           upper.tau = if (!sel.r) NA else m.i$upper.tau,
                           tau2 = if (!sel.r) NA else m.i$tau^2,
                           lower.tau2 = if (!sel.r) NA else m.i$lower.tau2,
                           upper.tau2 = if (!sel.r) NA else m.i$upper.tau2,
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
                           Q.b.common = NA,
                           Q.b.random = NA,
                           df.Q.b = NA,
                           pval.Q.b.common = NA,
                           pval.Q.b.random = NA,
                           ##
                           stringsAsFactors = FALSE)
    ##
    overall.i <- data.frame(name = name[i],
                            ##
                            TE.common = m.i$TE.common,
                            seTE.common = m.i$seTE.common,
                            lower.common = m.i$lower.common,
                            upper.common = m.i$upper.common,
                            statistic.common = m.i$statistic.common,
                            pval.common = m.i$pval.common,
                            ##
                            TE.random = m.i$TE.random,
                            seTE.random = m.i$seTE.random,
                            lower.random = m.i$lower.random,
                            upper.random = m.i$upper.random,
                            statistic.random = m.i$statistic.random,
                            pval.random = m.i$pval.random,
                            df.hakn = replaceNULL(m.i$df.hakn),
                            ##
                            n.harmonic.mean.ma =
                              1 / mean(1 / replaceNULL(m.i$n)),
                            t.harmonic.mean.ma =
                              1 / mean(1 / replaceNULL(m.i$time)),
                            ##
                            seTE.predict = m.i$seTE.predict,
                            lower.predict = m.i$lower.predict,
                            upper.predict = m.i$upper.predict,
                            ##
                            k = m.i$k,
                            k.study = m.i$k.study,
                            k.all = m.i$k.all,
                            k.TE = m.i$k.TE,
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
                            Q.w.common = NA,
                            Q.w.random = NA,
                            ##
                            Q.b.common = NA,
                            pval.Q.b.common = NA,
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
  if (missing.backtransf) {
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
  if (any(meth$prediction.subgroup))
    meth$prediction.subgroup <- TRUE
  else if (is.null(meth$prediction.subgroup) || anyNA(meth$prediction.subgroup))
    meth$prediction.subgroup <- FALSE
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
      if (pooled[i] == "common") {
        study.i$TE <- m.i$TE.common.w
        study.i$seTE <- m.i$seTE.common.w
        study.i$lower <- m.i$lower.common.w
        study.i$upper <- m.i$upper.common.w
        study.i$statistic <- m.i$statistic.common.w
        study.i$pval <- m.i$pval.common.w
        study.i$w.common <- m.i$w.common.w
        study.i$w.random <- 0
      }
      else {
        study.i$TE <- m.i$TE.random.w
        study.i$seTE <- m.i$seTE.random.w
        study.i$lower <- m.i$lower.random.w
        study.i$upper <- m.i$upper.random.w
        study.i$statistic <- m.i$statistic.random.w
        study.i$pval <- m.i$pval.random.w
        study.i$w.common <- 0
        study.i$w.random <- m.i$w.random.w
      }
    }
    else {
      study.i$n.e <- sum(replaceNULL(m.i$n.e.w))
      study.i$n.c <- sum(replaceNULL(m.i$n.c.w))
      ##
      if (pooled[i] == "common") {
        study.i$TE <- m.i$TE.common
        study.i$seTE <- m.i$seTE.common
        study.i$lower <- m.i$lower.common
        study.i$upper <- m.i$upper.common
        study.i$statistic <- m.i$statistic.common
        study.i$pval <- m.i$pval.common
        study.i$w.common <- 1
        study.i$w.random <- 0
      }
      else {
        study.i$TE <- m.i$TE.random
        study.i$seTE <- m.i$seTE.random
        study.i$lower <- m.i$lower.random
        study.i$upper <- m.i$upper.random
        study.i$statistic <- m.i$statistic.random
        study.i$pval <- m.i$pval.random
        study.i$w.common <- 0
        study.i$w.random <- 1
      }
    }
    ##
    study.i$subgroup <- name[i]
    ##
    if (i == 1)
      study <- study.i
    else
      study <- rbind(study, study.i)
  }
  
  
  if (length(unique(study$subgroup)) == 1) {
    res <- c(as.list(study), as.list(meth[1, ]), as.list(overall))
    res$subgroup <- NULL
  }
  else
    res <- c(as.list(study), as.list(meth[1, ]),
             as.list(overall), as.list(subgroup))
  ##
  ##
  res$data <- data
  res$n.harmonic.mean <- data$n.harmonic.mean
  res$t.harmonic.mean <- data$t.harmonic.mean
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
  res$TE.common <- makeunique(res$TE.common)
  res$seTE.common <- makeunique(res$seTE.common)
  res$lower.common <- makeunique(res$lower.common)
  res$upper.common <- makeunique(res$upper.common)
  res$statistic.common <- makeunique(res$statistic.common)
  res$pval.common <- makeunique(res$pval.common)
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
  res$k.study <- makeunique(res$k.study)
  res$k.all <- makeunique(res$k.all)
  res$k.TE <- makeunique(res$k.TE)
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
  res$n.harmonic.mean.ma <- makeunique(res$n.harmonic.mean.ma)
  res$t.harmonic.mean.ma <- makeunique(res$t.harmonic.mean.ma)
  ##
  res$Rb <- makeunique(res$Rb)
  res$lower.Rb <- makeunique(res$lower.Rb)
  res$upper.Rb <- makeunique(res$upper.Rb)
  ##
  res$Q.w.common <- makeunique(res$Q.w.common)
  res$Q.w.random <- makeunique(res$Q.w.random)
  res$df.Q.w <- makeunique(res$df.Q.w, 0)
  ##
  res$Q.b.common <- makeunique(res$Q.b.common)
  res$Q.b.random <- makeunique(res$Q.b.random)
  ##
  res$df.Q.b <- makeunique(res$df.Q.b, 0)
  res$pval.Q.b.common <-
    makeunique(makeunique(res$pval.Q.b.common,
                          pvalQ(res$Q.b.common, res$df.Q.b)))
  res$pval.Q.b.random <-
    makeunique(makeunique(res$pval.Q.b.random,
                          pvalQ(res$Q.b.random, res$df.Q.b)))
  ##
  res$show.studies <- show.studies
  res$overall.hetstat <- overall.hetstat
  
  
  res$is.subgroup <- is.subgroup
  
  
  if (!is.null(res$subgroup)) {
    res$subgroup.name <- "meta-analysis"
    res$bylevs <- unique(res$subgroup)
    res$w.common <- rep(0, length(res$w.common))
    res$w.common.w <- rep(0, length(res$w.common.w))
    res$w.random <- rep(0, length(res$w.random))
    res$w.random.w <- rep(0, length(res$w.random.w))
    res$lower.predict.w <- rep(NA, length(res$w.random.w))
    res$upper.predict.w <- rep(NA, length(res$w.random.w))
  }
  
  
  if (is.na(res$tau.preset))
    res$tau.preset <- NULL
  ##
  if (!unique.pooled) {
    res$overall <- FALSE
    res$overall.hetstat <- FALSE
  }
  ##
  res$pooled <- pooled
  res$is.limit.copas <- is.limit.copas
  ##
  ## Backward compatibility
  ##
  res$fixed <- res$common
  res$comb.fixed <- res$common
  res$comb.random <- res$random
  res$level.comb <- res$level.ma
  ##
  res$w.fixed <- res$w.common
  res$TE.fixed <- res$TE.common
  res$seTE.fixed <- res$seTE.common
  res$lower.fixed <- res$lower.common
  res$upper.fixed <- res$upper.common
  res$statistic.fixed <- res$statistic.common
  res$pval.fixed <- res$pval.common
  res$zval.fixed <- res$zval.common
  ##
  res$text.fixed <- res$text.common
  res$text.w.fixed <- res$text.w.common
  ##
  if (!is.null(res$subgroup)) {
    res$byvar <- res$subgroup
    res$bylab <- res$subgroup.name
    res$print.byvar <- res$print.subgroup.name
    res$byseparator <- res$sep.subgroup
    ##
    res$TE.fixed.w <- res$TE.common.w
    res$seTE.fixed.w <- res$seTE.common.w
    res$lower.fixed.w <- res$lower.common.w
    res$upper.fixed.w <- res$upper.common.w
    res$statistic.fixed.w <- res$statistic.common.w
    res$pval.fixed.w <- res$pval.common.w
    res$zval.fixed.w <- res$zval.common.w
    res$w.fixed.w <- res$w.common.w
    ##
    res$Q.w.fixed <- res$Q.w.common
    res$pval.Q.w.fixed <- res$pval.Q.w.common
    res$Q.b.fixed <- res$Q.b.common
    res$pval.Q.b.fixed <- res$pval.Q.b.common
  }
  
  
  class(res) <- c("metabind", "meta")
  
  
  res
}
