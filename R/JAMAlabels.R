#' Create study labels in JAMA layout (deprecated function)
#' 
#' @description
#' Deprecated function to create study labels in JAMA layout (for
#' forest plot). Replaced by \code{\link{labels.meta}}.
#' 
#' @details
#' This auxiliary function can be used to create study labels in JAMA
#' layout which can be added to a forest plot using argument
#' 'studlab'.
#' 
#' @param author A vector providing study authors.
#' @param year A vector providing year of publication.
#' @param citation A vector providing citation numbers.
#' @param data An optional data frame containing the study
#'   information.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{labels.meta}}, \code{\link{forest.meta}}
#' 
#' @examples
#' data(Fleiss1993bin)
#' 
#' refs <- 20 + 1:7
#' 
#' Fleiss1993bin$mylabs <-
#'   JAMAlabels(study, year, refs, data = Fleiss1993bin)
#'
#' m <- metabin(d.asp, n.asp, d.plac, n.plac, data = Fleiss1993bin,
#'   studlab = paste(study, year),
#'   sm = "OR", random = FALSE)
#' 
#' forest(m, studlab = mylabs, layout = "JAMA",
#'   fontfamily = "Times", fontsize = 10)
#' 
#' @export JAMAlabels


JAMAlabels <- function(author, year, citation, data = NULL) {
  
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  if (nulldata)
    data <- sfsp
  ##
  if (missing(author) | missing(year) | missing(citation))
    stop("Mandatory arguments: 'author', 'year', and 'citation'.",
         call. = FALSE)
  ##
  ## Catch 'author', 'year', and 'citation' from data:
  ##
  author <- catch("author", mc, data, sfsp)
  year <- catch("year", mc, data, sfsp)
  citation <- catch("citation", mc, data, sfsp)
  ##
  if (length(author) != length(year))
    stop("Arguments 'author' and 'year' must be of same length.",
         call. = FALSE)
  ##
  if (length(author) != length(citation))
    stop("Arguments 'author' and 'citation' must be of same length.",
         call. = FALSE)
  
  
  res <- NA
  for (i in seq(along = author))
    res[i] <-
      as.expression(substitute(study^ref~(year),
                               list(study = author[i],
                                    ref = citation[i],
                                    year = year[i])))
  ##
  res
}
