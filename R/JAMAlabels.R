#' Create study labels in JAMA layout (for forest plot)
#' 
#' @description
#' Create study labels in JAMA layout (for forest plot).
#' 
#' @details
#' This auxiliary function can be used to create study labels in JAMA
#' layout which can be added to a forest plot using argument
#' 'studlab'.
#' 
#' @param author A vector providing study authors.
#' @param year A vector providing year of publication.
#' @param citation A vector proving citation numbers.
#' @param data An optional data frame containing the study
#'   information.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{forest.meta}}
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
#'              studlab = paste(study, year),
#'              sm = "OR", comb.random = FALSE)
#' 
#' forest(m, studlab = mylabs, layout = "JAMA",
#'        fontfamily = "Times", fontsize = 10)
#' 
#' @export JAMAlabels


JAMAlabels <- function(author, year, citation, data = NULL) {
  
  nulldata <- is.null(data)
  ##
  if (nulldata)
    data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  if (missing(author) | missing(year) | missing(citation))
    stop("Mandatory arguments: 'author', 'year', and 'citation'.",
         call. = FALSE)
  ##
  ## Catch 'author', 'year', and 'citation' from data:
  ##
  author <- eval(mf[[match("author", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  ##
  year <- eval(mf[[match("year", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  ##
  citation <- eval(mf[[match("citation", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
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

