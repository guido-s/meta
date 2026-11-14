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
#'   For subgroup results, \code{df.Q.b} is a list column containing
#'   the complete degrees of freedom for each subgroup-difference test.
#'   Element \code{test.Q.b} identifies chi-squared (\code{"Q"}) and
#'   F (\code{"F"}) tests; unavailable types are \code{NA}.
#'   Variance fields \code{tau2}, \code{lower.tau2}, \code{upper.tau2},
#'   \code{se.tau2}, \code{tau}, \code{lower.tau}, and \code{upper.tau}
#'   are list columns containing the complete numeric vectors, including
#'   length-one vectors for single-component analyses. Component order and
#'   names are preserved.
#'   For random effects and prediction results, column \code{interval.slot}
#'   numbers the intervals within each meta-analysis and result type.
#'   Repeated uses of the same method have different interval numbers.
#'   Together, \code{id} and \code{interval.slot} link each subgroup result
#'   to the corresponding overall result.
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
  # Get rid of warning 'no visible binding for global variable'
  type <- method.interval <- id <- seTE <- method <- TE <- adhoc.hakn <-
    df <- df.hakn <- df.kero <- name <- NULL
  
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
  prediction <- prediction && !is.null(x$predict)
  subgroup <- subgroup && !is.null(x$subgroup)
  
  if (subgroup) {
    full.common <- full.random <- full.predict <- NULL
    #
    for (model in c(if (random) "random", if (prediction) "predict")) {
      if (anyDuplicated(x[[model]][c("id", "interval.slot")]))
        stop("Interval slots must uniquely identify pooled results ",
             "within each meta-analysis and result family.", call. = FALSE)
    }
    #
    if (common)
      full.common <-
        left_join(x$subgroup$common, x$common,
                  by = "id", suffix = c(".w", ""))
    #
    if (random)
      full.random <-
        left_join(select(x$subgroup$random, -method.interval, -adhoc.hakn), x$random,
                  by = c("id", "interval.slot"),
                  suffix = c(".w", ""))
    #
    if (prediction)
      full.predict <-
        left_join(select(x$subgroup$predict, -method.interval, -adhoc.hakn), x$predict,
                  by = c("id", "interval.slot"),
                  suffix = c(".w", ""))
    #
    full.results <- bind_rows(full.common, full.random, full.predict) %>%
      #
      relocate(type, method.interval, .after = id) %>%
      #
      relocate(seTE:text, .after = method.interval)
    #
    if (common || random) {
      full.results %<>%
        relocate(method, .after = type) %>%
        relocate(TE, .before = seTE)
    }
    #
    if (random || prediction) {
      full.results %<>%
        relocate(adhoc.hakn, df, df.hakn, df.kero, .after = method.interval)
    }
  }
  else {
    full.results <- bind_rows(
      if (common) x$common,
      if (random) x$random,
      if (prediction) x$predict)
  }
  #
  if (subgroup)
    res <- left_join(full.results, x$overall, by = "id", suffix = c(".w", ""))
  else
    res <- left_join(full.results, x$overall, by = "id")
  #
  res <- res %>%
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
    res$method.tau[res$type == "Common effect model"] <- NA
  #
  if (!is.null(res$method.tau.ci))
    res$method.tau.ci[res$type == "Common effect model"] <- NA
  #
  if (!is.null(res$tau.preset))
    res$tau.preset[res$type == "Common effect model"] <- NA
  #
  if (!is.null(res$three.level))
    res$three.level[res$type == "Common effect model"] <- NA
  #
  if (!is.null(res$rho))
    res$rho[res$type == "Common effect model"] <- NA
  #
  if (drop_all_NA)
    res <- drop_allNAs(res)
  #
  attr(res, "version") <- packageDescription("meta")$Version
  
  res
}
