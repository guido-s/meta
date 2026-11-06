#' Coerce to a data frame
#' 
#' @description
#' The \code{as.data.frame} method returns a data frame containing
#' information on meta-analysis results.
#' 
#' @param x An object of class \code{metabind}.
#' @param row.names \code{NULL} or a character vector giving the row
#'   names for the data frame.
#' @param optional Argument is ignored.
#' @param common A logical indicating whether to return results for the common
#'   effect model.
#' @param random A logical indicating whether to return results for the random
#'   effects model.
#' @param prediction A logical indicating whether to return prediction
#'   intervals.
#' @param subgroup A logical indicating whether to return subgroup results in
#'   addition to the overall results.
#' @param drop_all_NA A logical indicating whether to drop variables containing
#'   only missing values.
#' @param \dots Other arguments (ignored).
#' 
#' @return A data frame is returned by the function \code{as.data.frame}.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metabind}}
#' 
#' @examples
#' data(Fleiss1993cont)
#' #
#' ma1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, sm = "SMD", studlab = paste(study, year))
#' #
#' ma2 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, sm = "SMD", studlab = paste(study, year),
#'   method.smd = "glass")
#' 
#' mb12 <-
#'   metabind(ma1, ma2, name = c("Hedges' g", "Glass' delta"), common = FALSE)
#' 
#' as.data.frame(mb12)
#'
#' @method as.data.frame metabind
#' @export

as.data.frame.metabind <- function(x, row.names = NULL, optional = FALSE,
                                   common = !is.null(x$common),
                                   random = !is.null(x$random),
                                   prediction = !is.null(x$predict),
                                   subgroup = !is.null(x$subgroup),
                                   drop_all_NA = TRUE,
                                   ...) {
  
  chkclass(x, "metabind")
  x <- updateversion(x)
  #
  chklogical(common)
  chklogical(random)
  chklogical(prediction)
  chklogical(subgroup)
  chklogical(drop_all_NA)
  #
  common <- common && !is.null(x$common)
  random <- random && !is.null(x$random)
  prediction <- prediction && !is.null(x$prediction)
  subgroup <- subgroup && !is.null(x$subgroup)
  
  if (subgroup) {
    full.common <- full.random <- full.predict <- NULL
    #
    if (common)
      full.common <-
        left_join(x$subgroup$common, x$common,
                  by = "id", suffix = c(".w", ""))
    #
    if (random)
      full.random <-
        left_join(x$subgroup$random, x$random,
                  by = c("id", "method.interval", "adhoc.hakn"),
                  suffix = c(".w", ""))
    #
    if (prediction)
      full.predict <-
        left_join(x$subgroup$predict, x$predict,
                  by = c("id", "method.interval", "adhoc.hakn"),
                  suffix = c(".w", ""))
    #
    full.results <- bind_rows(full.common, full.random, full.predict)
    #
    full.results  %<>% relocate(model, method.interval, .after = id)
    #
    full.results  %<>% relocate(seTE:text, .after = method.interval)
    #
    if (common || random) {
      full.results %<>%
        relocate(method, .after = model) %>%
        relocate(TE, .before = seTE)
    }
    #
    if (random || prediction) {
      full.results %<>%
        relocate(adhoc.hakn, df, df.hakn, df.kero, .after = method.interval)
    }
  }
  else {
      full.results <- bind_rows(x$common, x$random, x$predict)
  }
  #
  res <- left_join(full.results, x$overall, by = "id") %>%
    relocate(name, .after = id)
  #
  if (!is.null(row.names)) {
    if (length(row.names) != nrow(res))
      stop("Argument 'row.names' must be of length ", nrow(res))
    row.names(res) <- row.names
  }
  else
    row.names(res) <- seq_len(nrow(res))
  #
  if (!is.null(res$method.tau))
    res$method.tau[res$model == "Common effect model"] <- NA
  #
  if (!is.null(res$method.tau.ci))
    res$method.tau.ci[res$model == "Common effect model"] <- NA
  #
  if (!is.null(res$tau.preset))
    res$tau.preset[res$model == "Common effect model"] <- NA
  #
  if (!is.null(res$three.level))
    res$three.level[res$model == "Common effect model"] <- NA
  #
  if (!is.null(res$rho))
    res$rho[res$model == "Common effect model"] <- NA
  #
  if (drop_all_NA)
    res <- drop_allNAs(res)
  #
  attr(res, "version") <- packageDescription("meta")$Version
  
  res
}
