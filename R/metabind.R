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
#' @param common A logical vector indicating whether results of common
#'   effect model should be considered.
#' @param random A logical vector indicating whether results of random
#'   effects model should be considered.
#' @param prediction A logical vector indicating whether results of
#'   prediction intervals should be considered.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If \code{backtransf=TRUE}
#'   (default), results for \code{sm="OR"} are printed as odds ratios rather
#'   than log odds ratios, for example.
#' @param outclab Outcome label for all meta-analysis objects.
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
#' @seealso \code{\link{forest.metabind}}, \code{\link{as.data.frame.metabind}}, 
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
#' ma <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, sm = "SMD")
#'
#' # Conduct two subgroup analyses
#' #
#' ma1 <- update(ma, subgroup = age, subgroup.name = "Age group")
#' ma2 <- update(ma, subgroup = region, subgroup.name = "Region")
#'
#' # Combine random effects subgroup meta-analyses and show forest
#' # plot with subgroup results
#' #
#' mb <- metabind(ma1, ma2, common = FALSE)
#' mb
#' forest(mb)
#'
#' # Use various estimation methods for between-study heterogeneity
#' # variance
#' #
#' ma.pm <- update(ma, method.tau = "PM")
#' ma.dl <- update(ma, method.tau = "DL")
#' ma.ml <- update(ma, method.tau = "ML")
#' ma.hs <- update(ma, method.tau = "HS")
#' ma.sj <- update(ma, method.tau = "SJ")
#' ma.he <- update(ma, method.tau = "HE")
#' ma.eb <- update(ma, method.tau = "EB")
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
#' ma.taus <- metabind(ma, ma.pm, ma.dl, ma.ml, ma.hs, ma.sj, ma.he, ma.eb,
#'   name = taus, common = FALSE)
#' ma.taus
#' forest(ma.taus)
#' 
#' @export metabind

metabind <- function(...,
                     name = NULL,
                     common = NULL, random = NULL, prediction = NULL,
                     backtransf = NULL, outclab = NULL) {
  
  
  #
  #
  # (1) Check and set arguments
  #
  #
  
  missing.name <- missing(name)
  missing.backtransf <- missing(backtransf)
  #
  args <- list(...)
  #
  # Extract and drop list element 'pooled'
  #
  pooled <- args$pooled
  args$pooled <- NULL
  #
  if (!is.null(pooled)) {
    pooled <- setchar(pooled, c("common", "random", "fixed"))
    #
    if (pooled == "fixed")
      pooled <- "common"
    #
    if (!is.null(common) || !is.null(random)) {
      warning("Deprecated argument 'pooled' ignored; input considered for ",
              "arguments 'common' and 'random'.", call. = FALSE)
    }
    else {
      common <- pooled == "common"
      random <- pooled == "random"
    }
  }
  #
  args <- drop_from_dots(args, "warn.deprecated", deprecated = TRUE)
  #
  args <-
    drop_from_dots(
      args, "subgroup", deprecated = TRUE,
      text.deprecated = "Use argument 'subgroup' in forest.metabind().")
  #
  n.meta <- length(args)
  seq.meta <- seq_len(n.meta)
  #
  if (n.meta == 1) {
    if (inherits(args[[1]], "meta.rm5")) {
      args <- args[[1]]
      #
      if (missing.name) {
        name <- unlist(lapply(args, "[[" , "outclab"))
        missing.name <- FALSE
      }
    }
    #
    else if (inherits(args[[1]], "copas"))
      return(metamerge(args[[1]]))
    #
    else if (inherits(args[[1]], "limitmeta"))
      return(metamerge(args[[1]]))
    #
    else if (inherits(args[[1]], "netpairwise"))
      stop("Elements of argument '...' may not be of class 'netpairwise'.",
           call. = FALSE)
    #
    else if (inherits(args[[1]], "meta"))
      return(args[[1]])
    #
    else if (!is.list(args[[1]]))
      stop("All elements of argument '...' must be of class 'meta', ",
           "'limitmeta', or 'copas'.",
           call. = FALSE)
    else {
      n.meta <- length(args[[1]])
      seq.meta <- seq_len(n.meta)
      #
      args2 <- list()
      for (i in seq.meta)
        args2[[i]] <- args[[1]][[i]]
      #
      args <- args2
    }
  }
  #  
  # Act on limitmeta and copas objects
  #
  name.i <- vector("character", n.meta)
  #
  is.limit <- is.copas <- is.trimfill <- is.subgroup <- rep(FALSE, n.meta)
  samedata <- rep(NA, n.meta)
  #
  for (i in seq.meta) {
    if (inherits(args[[i]], "metabind"))
      stop("Elements of argument '...' may not be of class 'metabind'.",
           call. = FALSE)
    #
    if (inherits(args[[i]], "metacum"))
      stop("Elements of argument '...' may not be of class 'metacum'.",
           call. = FALSE)
    #
    if (inherits(args[[i]], "metainf"))
      stop("Elements of argument '...' may not be of class 'metainf'.",
           call. = FALSE)
    #
    if (inherits(args[[i]], "netpairwise"))
      stop("Elements of argument '...' may not be of class 'netpairwise'.",
           call. = FALSE)
    #
    if (inherits(args[[i]], "meta")) {
      if (inherits(args[[i]], "metamerge"))
        stop("Meta-analysis objects created with metamerge() cannot ",
             "be used in metabind().",
             call. = FALSE)
      #
      if (!inherits(args[[i]], "metaadd"))
        args[[i]] <- updateversion(args[[i]])
      #
      if (missing.name) {
        if (inherits(args[[i]], "trimfill")) {
          is.trimfill[i] <- TRUE
          #
          name.i[i] <- "trimfill"
        }
        else
          name.i[i] <- replaceNULL(args[[i]]$subgroup.name, "")
        #
        if (name.i[i] == "")
          name.i[i] <- class(args[[i]])[1]
      }
    }
    #
    else if (inherits(args[[i]], "copas")) {
      args[[i]] <-
        metamerge(update(args[[i]]$x, common = FALSE, random = FALSE),
                  args[[i]], label2 = "copas")
      #
      is.copas[i] <- TRUE
      #
      if (missing.name)
        name.i[i] <- "copas"
    }
    #
    else if (inherits(args[[i]], "limitmeta")) {
      args[[i]] <-
        metamerge(update(args[[i]]$x, common = FALSE, random = FALSE),
                  args[[i]], label2 = "limit")
      #
      is.limit[i] <- TRUE
      #
      if (missing.name)
        name.i[i] <- "limitmeta"
    }
    #
    else
      stop("All elements of argument '...' must be of class 'meta', ",
           "'limitmeta', or 'copas'.",
           call. = FALSE)
    #
    is.subgroup[i] <- !is.null(args[[i]]$subgroup)
    samedata[i] <- !samedata(args[[1]], args[[i]], stop = FALSE)
  }
  #
  # Meta-analyses must all contain subgroups or none
  #
  if (length(unique(is.subgroup)) != 1)
    stop("All or none meta-analyses must contain subgroups.",
         call. = FALSE)
  #
  with.subgroups <- any(is.subgroup)
  is.limit.copas <- is.limit | is.copas
  #
  samedata[is.limit.copas] <- TRUE
  samedata <- all(samedata)
  #
  if (!is.null(common)) {
    chklogical(common[1], "common")
    #
    if (length(common) == 1)
      common <- rep_len(common, n.meta)
    else
      chklength(common, n.meta,
                text = paste("Length of argument 'common' differs from",
                             "number of meta-analyses."))
  }
  else {
    for (i in seq.meta)
      common <- c(common, replaceNULL(args[[i]]$common, FALSE))
  }
  #
  if (!is.null(random)) {
    chklogical(random[1], "random")
    #
    if (length(random) == 1)
      random <- rep_len(random, n.meta)
    else
      chklength(random, n.meta,
                text = paste("Length of argument 'random' differs from",
                             "number of meta-analyses."))
  }
  else {
    for (i in seq.meta)
      random <- c(random, replaceNULL(args[[i]]$random, FALSE))
  }
  #
  if (!missing(prediction)) {
    if (length(prediction) == 1)
      prediction <- rep_len(prediction, n.meta)
    #
    if (!is.logical(prediction))
      stop("Argument 'prediction' must contain logical values.",
           call. = TRUE)
    #
    chklength(prediction, n.meta,
              text = paste("Length of argument 'prediction' differs from",
                           "number of meta-analyses."))
  }
  else {
    for (i in seq.meta)
      prediction <- c(prediction, replaceNULL(args[[i]]$prediction, FALSE))
  }
  #
  for (i in seq.meta) {
    args[[i]]$common <- common[i]
    args[[i]]$random <- random[i]
    args[[i]]$prediction <- prediction[i]
  }
  #
  # Name of meta-analysis object
  #
  if (missing.name) {
    name <- name.i
    name[name == ""] <- paste0("meta", seq.meta[name == ""])
  }
  else
    chklength(name, length(is.subgroup),
              text =
                paste("Number of meta-analyses and names provided in",
                      "argument 'name' differ."))
  #
  # Names for meta-analyses must be unique
  #
  if (length(unique(name)) != length(name)) {
    for (i in seq.meta)
      if (name[i] %in% c("metabin", "metainc", "metaprop", "metarate") &
          !is.trimfill[i])
        name[i] <- paste(name[i], args[[i]]$method, sep = ".")
  }
  #
  if (length(unique(name)) != length(name)) {
    for (i in seq.meta)
      if (random[i])
        name[i] <- paste(name[i], args[[i]]$method.tau, sep = ".")
  }
  #
  if (length(unique(name)) != length(name))
    name <- paste0("meta", seq.meta)
  #
  # Back-transformation
  #
  if (missing.backtransf)
    backtransf <- any(condense(args, "backtransf"))
  else
    chklogical(backtransf)


  #
  #
  # (2) Extract meta-analytical methods
  #
  #

  if (with.subgroups) {
    # TODO
    warning("Not yet implemented.")
    return(invisible(metagen(1, 1)))
  }
  else {
    #
    # Extract overall results
    #
    res.overall <- res.common <- res.random <- res.predict <- NULL
    #
    for (i in seq.meta) {
      res.overall <-
        bind_rows(res.overall,
                  extrMA(i, args, name, "overall", outclab))
      #
      if (args[[i]]$common) {
        res.common <-
          bind_rows(res.common,
                    extrMA(i, args, name, "common", outclab))
      }
      #
      if (args[[i]]$random) {
        res.random <-
          bind_rows(res.random,
                    extrMA(i, args, name, "random", outclab))
      }
      #
      if (args[[i]]$prediction) {
        res.predict <-
          bind_rows(res.predict,
                    extrMA(i, args, name, "prediction", outclab))
      }
    }
    #
    # Check whether IDs are unique and consecutive
    #
    ids <- res.overall$id
    #
    if (is.null(ids) || anyDuplicated(ids) ||
        !identical(sort(unique(ids)), seq_len(nrow(res.overall)))) {
      stop(
        "Internal error in metabind(): invalid IDs in element 'overall'. ",
        "Please report this issue with a reproducible example and ",
        "session information at ",
        "https://github.com/guido-s/meta/issues.",
        call. = FALSE
      )
    }
    #
    if (!is.null(res.common))
      rownames(res.common) <- seq_len(nrow(res.common))
    #
    if (!is.null(res.random))
      rownames(res.random) <- seq_len(nrow(res.random))
    #
    if (!is.null(res.predict))
      rownames(res.predict) <- seq_len(nrow(res.predict))
    #
    # Cleanup
    #
    res.overall <- drop_allNAs(res.overall)
    #
    if (!is.null(res.common))
      res.common <- drop_allNAs(res.common)
    #
    if (!is.null(res.random)) {
      res.random <- drop_allNAs(res.random)
      #
      res.random$df.hakn[res.random$method.interval != "HK"] <- NA
      res.random$df.kero[res.random$method.interval != "KR"] <- NA
    }
    #
    if (!is.null(res.predict)) {
      res.predict <- drop_allNAs(res.predict)
      #
      res.predict$df.hakn[!grepl("HK", res.predict$method.interval)] <- NA
      res.predict$df.kero[!grepl("KR", res.predict$method.interval)] <- NA
    }
    #
    # Return results
    #
    res <- list(
      overall = res.overall,
      #
      common = res.common,
      random = res.random,
      predict = res.predict,
      #
      with.subgroups = with.subgroups,
      #
      backtransf = backtransf,
      #
      version = packageDescription("meta")$Version,
      call = match.call())
    #
    class(res) <- "metabind"
    #
    return(res)
  }
}
