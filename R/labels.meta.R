#' Create study labels for forest plot
#' 
#' @description
#' Create study labels for forest plot.
#' 
#' @details
#' This auxiliary function can be used to create study labels in JAMA
#' or Lancet layout which can be added to a forest plot using argument
#' 'studlab'.
#' 
#' @param object An object of class \code{meta}.
#' @param author An optional vector providing study authors.
#' @param year An optional vector providing year of publication.
#' @param citation An optional vector providing citation numbers.
#' @param layout A character string specifying layout. Either "JAMA"
#'   or "Lancet".
#' @param data An optional data frame containing the study
#'   information.
#' @param \dots Additional arguments (ignored at the moment).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{forest.meta}}
#' 
#' @examples
#' data(Fleiss1993bin)
#' 
#' refs <- 20 + 1:7
#'
#' m <- metabin(d.asp, n.asp, d.plac, n.plac, data = Fleiss1993bin,
#'  studlab = study, sm = "OR", random = FALSE)
#' 
#' forest(m,
#'   studlab = labels(m, year = year, citation = refs, layout = "JAMA"),
#'   layout = "JAMA", fontfamily = "Times", fontsize = 10)
#' 
#' forest(m,
#'   studlab = labels(m, year = year, citation = refs, layout = "Lancet"))
#' 
#' @method labels meta
#' @export


labels.meta <- function(object,
                        author = object$studlab, year = "", citation = NULL,
                        layout = "JAMA",
                        data = object$data,
                        ...) {
  
  chkclass(object, "meta")
  layout <- setchar(layout, c("JAMA", "Lancet"))
  
  
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  ## Catch 'author', 'year', and 'citation' from data:
  ##
  if (!missing(author))
    author <- catch("author", mc, data, sfsp)
  ##
  if (!missing(year))
    year <- catch("year", mc, data, sfsp)
  else
    year <- rep_len("", length(author))
  ##
  if (!missing(citation))
    citation <- catch("citation", mc, data, sfsp)
  else
    citation <- seq_along(author)
  ##
  if (length(author) != length(year))
    stop("Arguments 'author' and 'year' must be of same length.",
         call. = FALSE)
  ##
  if (length(author) != length(citation))
    stop("Arguments 'author' and 'citation' must be of same length.",
         call. = FALSE)
  
  
  res <- rep_len(NA, length(author))
  ##
  if (layout == "JAMA") {
    for (i in seq(along = author)) {
      res[i] <-
        if (year[i] == "" | is.na(year[i]))
          as.expression(substitute(study^ref,
                                   list(study = author[i],
                                        ref = citation[i],
                                        year = year[i])))
        else
          as.expression(substitute(study^ref~(year),
                                   list(study = author[i],
                                        ref = citation[i],
                                        year = year[i])))
    }
  }
  else if (layout == "Lancet") {
    for (i in seq(along = author)) {
      res[i] <-
        if (year[i] == "" | is.na(year[i]))
          as.expression(substitute(study^ref,
                                   list(study = author[i],
                                        ref = citation[i],
                                        year = year[i])))
        else
          as.expression(substitute(study~(year)^ref,
                                   list(study = author[i],
                                        ref = citation[i],
                                        year = year[i])))
    }
  }
  ##
  res
}
