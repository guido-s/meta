#' Get default for a meta-analysis setting.
#' 
#' @description
#' Get default for a meta-analysis setting in R package \bold{meta}.
#' 
#' @param x A character string or vector with setting name(s).
#' @param unname A logical whether to remove names from attributes.
#' 
#' @details
#' This function can be used to get the default for a meta-analysis
#' setting defined using \code{\link{settings.meta}}.
#' 
#' This function is primarily used to define default settings in
#' meta-analysis functions, e.g., \code{\link{metabin}} or
#' \code{\link{metacont}}. A list of all arguments with current
#' settings is printed using the command
#' \code{settings.meta("print")}.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{settings.meta}}
#' 
#' @examples
#' # Get default setting for confidence interval of random effects
#' # model
#' #
#' gs("method.random.ci")
#' 
#' # Get default setting for summary measure in metabin()
#' #
#' gs("smbin")
#' 
#' @export gs


gs <- function(x = NULL, unname = NULL) {
  
  if (is.null(x))
    return(NULL)
  #
  chkchar(x)
  if (!is.null(unname))
    chklogical(unname)
  
  x.list <- vector("list", length(x))
  #
  for (i in seq_along(x)) {
    x.list[[i]] <-
    setchar(x[i],
            c(.settings$argslist, .settings$argslist.internal),
            "unmatched", name = x[i],
            stop.at.error = FALSE)
  }
  x.list <- compact(x.list)
  #
  if (is.null(x.list) | length(x.list) == 0)
    return(NULL)
  
  res <- vector("list", length(x.list))
  #
  for (i in seq_along(x.list))
    res[[i]] <- settings.meta(quietly = TRUE)[[x.list[[i]]]]
  #
  if (is.null(res) | length(res) == 0)
    return(NULL)
  #
  if (is.null(unname))
    unname <- length(res) <= 1
  #
  if (!unname)
    names(res) <- unlist(x.list)
  
  if (length(res) == 1)
    res <- unlist(res)
  
  res
}
