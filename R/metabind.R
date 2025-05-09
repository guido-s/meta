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
#' @param subgroup An optional variable to generate a forest plot with
#'   subgroups.
#' @param common A logical vector indicating whether results of common
#'   effect model should be considered.
#' @param random A logical vector indicating whether results of random
#'   effects model should be considered.
#' @param prediction A logical vector indicating whether results of
#'   prediction intervals should be considered.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=TRUE} (default), results for \code{sm="OR"} are
#'   printed as odds ratios rather than log odds ratios, for example.
#' @param outclab Outcome label for all meta-analyis objects.
#' @param pooled Deprecated argument (replaced by \code{common} and
#'   \code{random}.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' 
#' @details
#' This function can be used to combine any number of meta-analysis
#' objects which is useful, for example, to summarize results of
#' various meta-analysis methods or to generate a forest plot with
#' results of several subgroup analyses (see Examples).
#' 
#' Individual study results are not retained with \code{metabind} as
#' the function allows to combine meta-analyses from different data
#' sets (e.g., with randomised or observational studies). Individual
#' study results are retained with R function \code{\link{metamerge}}
#' which can be used to combine results of meta-analyses of the
#' same dataset.
#' 
#' @return
#' An object of class \code{c("metabind", "meta")} with corresponding
#' generic functions (see \code{\link{meta-object}}).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
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
#'   data = Fleiss1993cont, sm = "SMD")
#'
#' # Conduct two subgroup analyses
#' #
#' mu1 <- update(m1, subgroup = age, subgroup.name = "Age group")
#' mu2 <- update(m1, subgroup = region, subgroup.name = "Region")
#'
#' # Combine random effects subgroup meta-analyses and show forest
#' # plot with subgroup results
#' #
#' mb1 <- metabind(mu1, mu2, common = FALSE)
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
#'   name = taus, common = FALSE)
#' m1.taus
#' forest(m1.taus)
#' 
#' @export metabind


metabind <- function(..., subgroup = NULL,
                     name = NULL,
                     common = NULL, random = NULL, prediction = NULL,
                     backtransf = NULL, outclab = NULL,
                     pooled = NULL,
                     warn.deprecated = gs("warn.deprecated")) {
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  
  missing.subgroup <- missing(subgroup)
  subgroup.meta <- subgroup
  missing.name <- missing(name)
  missing.pooled <- missing(pooled)
  missing.backtransf <- missing(backtransf)
  ##
  chklogical(warn.deprecated)
  ##
  args <- list(...)
  ##
  n.meta <- length(args)
  seq.meta <- seq_len(n.meta)
  ##
  if (n.meta == 1) {
    if (inherits(args[[1]], "meta.rm5")) {
      args <- args[[1]]
      ##
      if (missing.name) {
        name <- unlist(lapply(args, "[[" , "outclab"))
        missing.name <- FALSE
      }
    }
    ##
    else if (inherits(args[[1]], "copas"))
      return(metamerge(args[[1]]))
    ##
    else if (inherits(args[[1]], "limitmeta"))
      return(metamerge(args[[1]]))
    ##
    else if (inherits(args[[1]], "netpairwise"))
      stop("Elements of argument '...' may not be of class 'netpairwise'.",
           call. = FALSE)
    ##
    else if (inherits(args[[1]], "meta"))
      return(args[[1]])
    ##
    else if (!is.list(args[[1]]))
      stop("All elements of argument '...' must be of class 'meta', ",
           "'limitmeta', or 'copas'.",
           call. = FALSE)
    else {
      n.meta <- length(args[[1]])
      seq.meta <- seq_len(n.meta)
      ##
      args2 <- list()
      for (i in seq.meta)
        args2[[i]] <- args[[1]][[i]]
      ##
      args <- args2
    }
  }
  ##
  if (!missing.pooled) {
    pooled <- setchar(pooled, c("common", "random", "fixed"))
    pooled[pooled == "fixed"] <- "common"
    ##
    pooled.common <- pooled == "common"
    pooled.random <- pooled == "random"
  }
  ##
  common <-
    deprecated2(common, missing(common), pooled.common, missing.pooled,
                warn.deprecated, "pooled")
  ##
  if (!is.null(common)) {
    chklogical(common[1], "common")
    ##
    if (length(common) == 1)
      common <- rep_len(common, n.meta)
    else
      chklength(common, n.meta,
                text = paste("Length of argument 'common' differs from",
                             "number of meta-analyses."))
  }
  ##
  random <-
    deprecated2(random, missing(random), pooled.random, missing.pooled,
                warn.deprecated, "pooled")
  ##
  if (!is.null(random)) {
    chklogical(random[1], "random")
    ##
    if (length(random) == 1)
      random <- rep_len(random, n.meta)
    else
      chklength(random, n.meta,
                text = paste("Length of argument 'random' differs from",
                             "number of meta-analyses."))
  }
  ##
  if (!missing(prediction)) {
    if (length(prediction) == 1)
      prediction <- rep_len(prediction, n.meta)
    ##
    if (!is.logical(prediction))
      stop("Argument 'prediction' must contain logical values.",
           call. = TRUE)
    ##
    chklength(prediction, n.meta,
              text = paste("Length of argument 'prediction' differs from",
                           "number of meta-analyses."))
  }
  ##
  if (!missing.backtransf)
    chklogical(backtransf)
  ##  
  ## Act on limitmeta and copas objects
  ##
  name.i <- vector("character", n.meta)
  ##
  is.limit <- is.copas <- is.trimfill <- is.subgroup <- rep(FALSE, n.meta)
  samedata <- rep(NA, n.meta)
  ##
  for (i in seq.meta) {
    if (inherits(args[[i]], "metabind"))
      stop("Elements of argument '...' may not be of class 'metabind'.",
           call. = FALSE)
    ##
    if (inherits(args[[i]], "netpairwise"))
      stop("Elements of argument '...' may not be of class 'netpairwise'.",
           call. = FALSE)
    ##
    if (inherits(args[[i]], "meta")) {
      if (inherits(args[[i]], "metamerge"))
        stop("Meta-analysis objects created with metamerge() cannot ",
             "be used in metabind().",
             call. = FALSE)
      ##
      if (!inherits(args[[i]], "metaadd"))
        args[[i]] <- updateversion(args[[i]])
      ##
      if (missing.name) {
        if (inherits(args[[i]], "trimfill")) {
          is.trimfill[i] <- TRUE
          ##
          name.i[i] <- "trimfill"
        }
        else
          name.i[i] <- replaceNULL(args[[i]]$subgroup.name, "")
        ##
        if (name.i[i] == "")
          name.i[i] <- class(args[[i]])[1]
      }
    }
    ##
    else if (inherits(args[[i]], "copas")) {
      args[[i]] <-
        metamerge(update(args[[i]]$x, common = FALSE, random = FALSE),
                  args[[i]], label2 = "copas")
      ##
      is.copas[i] <- TRUE
      ##
      if (missing.name)
        name.i[i] <- "copas"
    }
    ##
    else if (inherits(args[[i]], "limitmeta")) {
      args[[i]] <-
        metamerge(update(args[[i]]$x, common = FALSE, random = FALSE),
                  args[[i]], label2 = "limit")
      ##
      is.limit[i] <- TRUE
      ##
      if (missing.name)
        name.i[i] <- "limitmeta"
    }
    ##
    else
      stop("All elements of argument '...' must be of class 'meta', ",
           "'limitmeta', or 'copas'.",
           call. = FALSE)
    ##
    is.subgroup[i] <- !is.null(args[[i]]$subgroup)
    samedata[i] <- !samedata(args[[1]], args[[i]], stop = FALSE)
  }
  ##
  ## Meta-analyses must all contain subgroups or none
  ##
  if (length(unique(is.subgroup)) != 1)
    stop("All or none meta-analyses must contain subgroups.",
         call. = FALSE)
  ##
  with.subgroups <- any(is.subgroup)
  is.limit.copas <- is.limit | is.copas
  ##
  samedata[is.limit.copas] <- TRUE
  samedata <- all(samedata)
  
  
  ##
  ##
  ## (2) Extract meta-analytical methods
  ##
  ##
  
  meth.list <- vector("list", n.meta)
  ##
  for (i in seq.meta)
    meth.list[[i]] <-
      meta2meth(args[[i]], outclab)
  ##
  meth <- list()
  ##
  for (i in names(meth.list[[1]]))
    meth[[i]] <- condense(meth.list, i)
  ##
  meth$common <- replaceNULL(common, meth$common)
  meth$common[is.limit.copas] <- FALSE
  meth$random <- replaceNULL(random, meth$random)
  meth$prediction <- replaceNULL(prediction, meth$prediction)
  ##
  ## Use common effect estimate if no result is selected
  ##
  meth$common[!(meth$common | meth$random | meth$prediction)] <- TRUE
  ##
  ## Check whether settings are unique
  ##
  meth2 <- meth[c("sm", "level.ma", "level.predict", "null.effect")]
  n.meth <- lapply(meth2, function(x) length(unique(x)))
  ##
  if (any(n.meth != 1))
    stop("Setting for the following argument",
         if (sum(n.meth != 1) > 1) "s",
         " must be the same for all meta-analyses: ",
         paste0(paste0("'", names(meth2)[n.meth != 1], "'"),
                collapse = " - "))
  ##
  ## Unify settings
  ##
  meth$sm <- makeunique(meth$sm)
  meth$method <- unique(meth$method)
  meth$method.random <- unique(meth$method.random)
  meth$three.level <- makeunique(meth$three.level, FALSE)
  ##
  meth$level <- meth$level.ma <- makeunique(meth$level.ma)
  meth$level.predict <- makeunique(meth$level.predict)
  ##
  common <- meth$common
  random <- meth$random
  prediction <- meth$prediction
  ##
  meth$overall <- with.subgroups & samedata
  meth$overall.hetstat <- with.subgroups & samedata
  ##
  if (missing.pooled & all(meth$common) & all(!meth$random))
    meth$common <- pooled == "onlycommon"
  else
    meth$common <- any(meth$common)
  ##
  meth$random <- any(meth$random)
  ##
  meth$prediction <- any(meth$prediction)
  ##
  meth$method.common.ci <- unique(meth$method.common.ci)
  meth$method.random.ci <- unique(meth$method.random.ci)
  meth$adhoc.hakn.ci <- unique(meth$adhoc.hakn.ci)
  meth$method.tau <- unique(meth$method.tau)
  meth$tau.preset <- unique(meth$tau.preset)
  meth$TE.tau <- unique(meth$TE.tau)
  meth$tau.common <-
    if (!with.subgroups)
      FALSE
    else
      unique(meth$tau.common)
  #
  meth$method.I2 <- unique(meth$method.I2)
  #
  meth$prediction.subgroup <-
    if (!with.subgroups)
      FALSE
    else
      unique(meth$prediction.subgroup)
  meth$method.predict <- unique(meth$method.predict)
  meth$adhoc.hakn.pi <- unique(meth$adhoc.hakn.pi)
  ##
  meth$method.bias <- unique(meth$method.bias)
  meth$null.effect <- makeunique(meth$null.effect) 
  ##
  meth$title <- makeunique(meth$title, "")
  meth$complab <- makeunique(meth$complab, "")
  meth$outclab <- makeunique(meth$outclab, "")
  ##
  meth$label.e <- makeunique(meth$label.e, "")
  meth$label.c <- makeunique(meth$label.c, "")
  meth$label.left <- makeunique(meth$label.left, "")
  meth$label.right <- makeunique(meth$label.right, "")
  ##
  meth$print.subgroup.name <- makeunique(meth$print.subgroup.name, FALSE)
  meth$sep.subgroup <- makeunique(meth$sep.subgroup, "")
  meth$warn <- makeunique(meth$warn, FALSE)
  ##
  meth$backtransf <-
    if (missing.backtransf)
      any(meth$backtransf)
    else
      backtransf
  meth$pscale <- makeunique(meth$pscale, 1)
  meth$irscale <- makeunique(meth$irscale, 1)
  meth$irunit <- makeunique(meth$irunit, "")
  ##
  if (!(any(common) | any(random) | any(prediction))) {
    warning("No results to bind.", call. = FALSE)
    return(NULL)
  }
  ##
  ## Name of meta-analysis object
  ##
  if (missing.name) {
    name <- name.i
    name[name == ""] <- paste0("meta", seq.meta[name == ""])
  }
  else
    chklength(name, length(is.subgroup),
              text =
                paste("Number of meta-analyses and names provided in",
                      "argument 'name' differ."))
  ##
  ## Names for meta-analyses must be unique
  ##
  if (length(unique(name)) != length(name) & missing.subgroup) {
    for (i in seq.meta)
      if (name[i] %in% c("metabin", "metainc", "metaprop", "metarate") &
          !is.trimfill[i])
        name[i] <- paste(name[i], args[[i]]$method, sep = ".")
  }
  ##
  if (length(unique(name)) != length(name) & missing.subgroup) {
    for (i in seq.meta)
      if (random[i])
        name[i] <- paste(name[i], args[[i]]$method.tau, sep = ".")
  }
  ##
  if (length(unique(name)) != length(name) & missing.subgroup)
    name <- paste0("meta", seq.meta)
  ##
  ## Check if more than one common effect / random effects CI or PI is
  ## provided
  ##
  for (i in seq.meta) {
    if (length(args[[i]]$lower.common) > 1) {
      ith <- min(i, 4)
      stop("More than one result for common effect model provided in ",
           i, switch(ith, "st", "nd", "rd", "th"),
           " meta-analysis.",
           call. = FALSE)
    }
  }
  ##
  for (i in seq.meta) {
    if (length(args[[i]]$lower.random) > 1) {
      ith <- min(i, 4)
      stop("More than one result for random effects model provided in ",
           i, switch(ith, "st", "nd", "rd", "th"),
           " meta-analysis.",
           call. = FALSE)
    }
  }
  ##
  for (i in seq.meta) {
    if (length(args[[i]]$lower.predict) > 1) {
      ith <- min(i, 4)
      stop("More than one prediction interval provided in ",
           i, switch(ith, "st", "nd", "rd", "th"),
           " meta-analysis.",
           call. = FALSE)
    }
  }
  
  
  ## ##
  ## ## Show individual results
  ## ##
  show.studies <- TRUE
  overall.hetstat <- TRUE
  ## ##
  ## if (length(unique(meth$method)) != 1 |
  ##     length(unique(meth$method.tau)) != 1) {
  ##   show.studies <- FALSE
  ##   overall.hetstat <- FALSE
  ## }
  
  
  ##
  ##
  ## (3) Extract meta-analysis data
  ##
  ##
  
  meta.list <- vector("list", n.meta)
  ##
  if (with.subgroups)
    for (i in seq.meta)
      meta.list[[i]] <-
        subgr2meta(args[[i]], common[i], random[i], prediction[i], name[i])
  else {
    for (i in seq.meta)
      meta.list[[i]] <-
        overall2meta(args[[i]], common[i], random[i], prediction[i], name[i])
  }
  ##
  meta <- list()
  ##
  for (i in names(meta.list[[1]]))
    meta[[i]] <- condense(meta.list, i)
  
  
  ##
  ##
  ## (4) Extract meta-analysis results and store them in subgroups
  ##
  ##
  
  subgroup.list <- vector("list", n.meta)
  ##
  for (i in seq.meta)
    subgroup.list[[i]] <- overall2subgr(args[[i]])
  ##
  subgroup <- list()
  ##
  for (i in names(subgroup.list[[1]]))
    subgroup[[i]] <- condense(subgroup.list, i)
  ##
  if (!with.subgroups) {
    subgroup$subgroup.levels <- name
  }
  
  
  ##
  ##
  ## (5) Extract study data
  ##
  ##
  
  data.list <- vector("list", n.meta)
  ##
  if (with.subgroups) {
    for (i in seq.meta)
      data.list[[i]] <-
        subgr2data(args[[i]], common[i], random[i], prediction[i], name[i])
  }
  else {
    if (!is.null(subgroup.meta)) {
      if (length(subgroup.meta == 1))
        subgroup.meta <- rep_len(subgroup.meta, n.meta)
      else if (length(subgroup.meta) != n.meta)
        stop("Argument 'subgroup' must be of same length as ",
             "number of meta-analysis results.",
             call. = FALSE)
    }
    else
      subgroup.meta <- rep("", n.meta)
    #
    for (i in seq.meta)
      data.list[[i]] <-
        overall2data(args[[i]], common[i], random[i], prediction[i], name[i],
                     subgroup.meta[[i]])
  }
  #
  data <- do.call("rbind", data.list)
  #
  if (!with.subgroups & all(subgroup.meta == ""))
    data <- data[, names(data) != "subgroup"]
  
  
  ##
  ##
  ## (6) Generate meta-analysis object
  ##
  ##
  
  res <- c(meta, meth, subgroup)
  ##
  res$data <- data
  ##
  res$common.meta <- common
  res$random.meta <- random
  res$prediction.meta <- prediction
  ##
  if (!with.subgroups) {
    res$n.harmonic.mean <- res$n.harmonic.mean.ma
    res$t.harmonic.mean <- res$t.harmonic.mean.ma
  }
  else {
    res$n.harmonic.mean.ma <- unique(res$n.harmonic.mean.ma)
    res$t.harmonic.mean.ma <- unique(res$t.harmonic.mean.ma)
  }
  #
  res$call <- match.call()
  res$version <- packageDescription("meta")$Version
  ##
  res$show.studies <- show.studies
  ##
  res$with.subgroups <- with.subgroups
  res$samedata <- samedata
  ##
  res$method <- unique(res$method)
  res$method.random <- unique(res$method.random)
  res$method.common.ci <- unique(res$method.common.ci)
  res$method.random.ci <- unique(res$method.random.ci)
  res$three.level <- unique(res$three.level)
  res$adhoc.hakn.ci <- unique(res$adhoc.hakn.ci)
  res$tau.common <- unique(res$tau.common)
  ##
  if (with.subgroups & samedata) {
    sel.c <-
      !duplicated(
         data.frame(res$TE.common,
                    res$lower.common,
                    res$upper.common))
    ##
    res$TE.common <- res$TE.common[sel.c]
    res$seTE.common <- res$seTE.common[sel.c]
    res$lower.common <- res$lower.common[sel.c]
    res$upper.common <- res$upper.common[sel.c]
    res$statistic.common <- res$statistic.common[sel.c]
    res$pval.common <- res$pval.common[sel.c]
    ##
    res$text.common <-
      if (is.null(res$text.common))
        gs("text.common")
      else
        res$text.common <- res$text.common[sel.c]
    ##
    sel.r <-
      !duplicated(
         data.frame(res$TE.random,
                    res$lower.random,
                    res$upper.random))
    ##
    res$TE.random <- res$TE.random[sel.r]
    res$seTE.random <- res$seTE.random[sel.r]
    res$lower.random <- res$lower.random[sel.r]
    res$upper.random <- res$upper.random[sel.r]
    res$statistic.random <- res$statistic.random[sel.r]
    res$pval.random <- res$pval.random[sel.r]
    ##
    res$text.random <-
      if (is.null(res$text.random))
        gs("text.random")
      else
        res$text.random <- res$text.random[sel.r]
    ##
    res$df.random <- unique(res$df.random)
    res$df.hakn <- unique(res$df.hakn)
    res$df.kero <- unique(res$df.kero)
    ##
    sel.p <-
      !duplicated(
         data.frame(res$lower.predict,
                    res$upper.predict))
    ##
    res$lower.predict <- res$lower.predict[sel.p]
    res$upper.predict <- res$upper.predict[sel.p]
    ##
    res$text.predict <- 
      if (is.null(res$text.predict))
        gs("text.predict")
      else
        res$text.predict[sel.p]
    ##
    res$n.e <- unique(res$n.e)
    res$n.c <- unique(res$n.c)
    ##
    res$k <- unique(res$k)
    res$k.study <- unique(res$k.study)
    res$k.all <- unique(res$k.all)
    res$k.TE <- unique(res$k.TE)
    ##
    res$tau2 <- unique(res$tau2)
    res$tau <- unique(res$tau)
    ##
    res$detail.tau <- ""
    ##
    res$H <- unique(res$H)
    res$lower.H <- unique(res$lower.H)
    res$upper.H <- unique(res$upper.H)
    ##
    res$I2 <- unique(res$I2)
    res$lower.I2 <- unique(res$lower.I2)
    res$upper.I2 <- unique(res$upper.I2)
    ##
    res$Rb <- unique(res$Rb)
    res$lower.Rb <- unique(res$lower.Rb)
    res$upper.Rb <- unique(res$upper.Rb)
    ##
    res$subgroup.name <- "meta-analysis"
    res$w.common <- res$w.random <- rep(0, length(res$TE))
    res$w.common.w <- rep(0, length(res$TE.common.w))
    res$w.random.w <- rep(0, length(res$TE.random.w))
    res$lower.predict.w <- rep(NA, length(res$w.random.w))
    res$upper.predict.w <- rep(NA, length(res$w.random.w))
    ##
    res$Q.b.common <- unique(res$Q.b.common)
    res$Q.b.random <- unique(res$Q.b.random)
    res$df.Q.b <- unique(res$df.Q.b)
    res$pval.Q.b.common <- unique(res$pval.Q.b.common)
    res$pval.Q.b.random <- unique(res$pval.Q.b.random)
  }
  else {
    res$TE.common <- res$seTE.common <-
      res$lower.common <- res$upper.common <-
        res$statistic.common <- res$pval.common <- NA
    ##
    res$TE.random <- res$seTE.random <-
      res$lower.random <- res$upper.random <-
        res$statistic.random <- res$pval.random <- NA
    ##
    res$lower.predict <- res$upper.predict <- NA
    ##
    if (!with.subgroups) {
      res$studlab <- res$data$studlab
      res$TE <- res$data$TE
      res$seTE <- res$data$seTE
      res$lower <- res$data$lower
      res$upper <- res$data$upper
      res$statistic <- res$data$statistic
      res$pval <- res$data$pval
      res$zval <- res$data$statistic
      ##
      res$w.common <- res$w.random <- rep_len(NA, length(res$studlab))
      ##
      res$text.w.common <- res$text.w.random <- ""
    }
  }
  ##
  if (all(is.na(res$tau.preset)))
    res$tau.preset <- NULL
  ##
  res$Q.w.common <- NA
  res$Q.w.random <- NA
  res$df.Q.w <- NA
  res$pval.Q.w.common <- NA
  res$pval.Q.w.random <- NA
  ##
  res$is.limit.copas <- is.limit.copas
  #
  res$classes <- unique(as.vector(sapply(args, function(x) class(x))))
  res$classes <- res$classes[res$classes != "meta"]
  #
  ## Backward compatibility
  ##
  res <- backward(res)
  ##  
  class(res) <- c("metabind", "meta")
  
  res$list <- list(meth = meth, meta = meta, subgroup = subgroup, data = data)
  
  res
}
