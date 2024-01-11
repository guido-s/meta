setcat <- function(x, labels) {
  if (is.null(labels))
    return(x)
  ##
  x <- setchar(x, labels, stop.at.error = FALSE, setNA = TRUE)
  factor(x, levels = labels)
}


setdom <- function(dm, tool, domains, domain.available) {

  n <- length(dm)
  ##
  n.names <- length(domains)
  n.domains <- sum(domain.available)
  ##
  n.tool <-
    if (tool %in% c("RoB1", "ROBINS-I", "ROBINS-E"))
      7
    else if (tool %in% c("RoB2-cluster", "RoB2-crossover"))
      6
    else if (tool == "RoB2")
      5
  
  if (n.names == n.domains)
    return(domains)
  else if (n.names == n - n.tool |
           n.names == n.domains - n.tool) {
    domains.tool <- dm[seq_len(n.tool)]
    domains.add <- vector(mode = "character")
    for (i in seq_len(n.names))
      domains.add[i] <- domains[i]
    return(c(domains.tool, domains.add))
  }
  else if (n.names == n.domains - n.tool) {
    domains.add <- 
    for (i in seq_len(n.names))
      domains[n.tool + i] <- domains[i]
    return(domains[domain.available])
  }
  else
    stop("Wrong number of domains provided for '", tool,
         "' (must be ", n.domains, n.tool, " or ", n.domains - n.tool, ").",
         call. = FALSE)
  
  invisible(NULL)
}


catleg <- function(x) {
  domains <- attr(x, "domains")
  ##
  vars <- names(x)
  vars <- vars[!(vars %in% c("Study", "Weight"))]
  vars[vars == "Overall"] <- "O"
  ##
  paste0("(", vars, ") ", domains)
}


catcat <- function(x) {
  cat <- attr(x, "categories")
  ##
  vars <- names(x)
  vars <- vars[!(vars %in% c("Study", "Weight"))]
  vars[vars == "Overall"] <- "O"
  ##
  unique.list <- function(x) {
    if (length(unique(sapply(x, length))) != 1)
      return(FALSE)
    ##
    for (i in seq_len(length(x[[1]])))
      if (length(unique(sapply(x, "[[", i))) != 1)
        return(FALSE)
    ##
    return(TRUE)
  }
  if (unique.list(cat)) {
    if (suppressWarnings(all(!is.na(as.numeric(cat[[1]])))))
      return(collapse(cat[[1]], collapse = ", ", quote = ""))
    else
      return(collapse(cat[[1]], collapse = ", "))
  }
  else
    return(paste0("(", vars, ") ", lapply(cat, collapse, collapse = ", ")))
}


catcol <- function(x) {
  col <- attr(x, "col")
  ##
  vars <- names(x)
  vars <- vars[!(vars %in% c("Study", "Weight"))]
  vars[vars == "Overall"] <- "O"
  ##
  unique.list <- function(x) {
    if (length(unique(sapply(x, length))) != 1)
      return(FALSE)
    ##
    for (i in seq_len(length(x[[1]])))
      if (length(unique(sapply(x, "[[", i))) != 1)
        return(FALSE)
    ##
    return(TRUE)
  }
  if (unique.list(col)) {
    if (suppressWarnings(all(!is.na(as.numeric(col[[1]])))))
      return(collapse(col[[1]], collapse = ", ", quote = ""))
    else
      return(collapse(col[[1]], collapse = ", "))
  }
  else
    return(paste0("(", vars, ") ",
                  lapply(col, collapse, collapse = ", ", quote = "")))
}


catsymb <- function(x) {
  symb <- attr(x, "symbols")
  ##
  vars <- names(x)
  vars <- vars[!(vars %in% c("Study", "Weight"))]
  vars[vars == "Overall"] <- "O"
  ##
  unique.list <- function(x) {
    if (length(unique(sapply(x, length))) != 1)
      return(FALSE)
    ##
    for (i in seq_len(length(x[[1]])))
      if (length(unique(sapply(x, "[[", i))) != 1)
        return(FALSE)
    ##
    return(TRUE)
  }
  if (unique.list(symb)) {
    if (suppressWarnings(all(!is.na(as.numeric(symb[[1]])))))
      return(collapse(symb[[1]], collapse = ", ", quote = ""))
    else
      return(collapse(symb[[1]], collapse = ", "))
  }
  else
    return(paste0("(", vars, ") ",
                  lapply(symb, collapse, collapse = ", ", quote = "")))
}


definecat <- function(avail, x, var, tool, warn) {

  if (!avail)
    return(NULL)
  
  varname <- deparse(substitute(x))
  ##
  robins <- tolower(substring(tool, 1, 6)) == "robins"
  rob <- tolower(substring(tool, 1, 3)) == "rob" & !robins
  
  if (is.null(x)) {
    if (tool == "RoB1")
      x <-
        c("Low risk of bias", "Unclear risk of bias", "High risk of bias")
    else if (rob)
      x <-
        c("Low risk of bias", "Some concerns", "High risk of bias")
    else if (robins)
      x <-
        c("Low risk", "Some concerns", "High risk", "Very high risk", "NI")
    else {
      if (warn)
        warning("Alpha-numeric order of available categories used as ",
                "argument '",  varname, "' is missing.",
                call. = FALSE)
      x <- unique(sort(var))
    }
  }
  ##
  x
}


definesymb <- function(avail, x, cat, tool, warn) {
  
  if (!avail)
    return(NULL)
  
  varname <- deparse(substitute(x))
  ##
  robins <- tolower(substring(tool, 1, 6)) == "robins"
  rob <- tolower(substring(tool, 1, 3)) == "rob" & !robins
  
  if (is.null(x)) {
    if (tool == "user-defined")
      x <- FALSE
    else {
      if (rob)
        x <- c("+", "?", "-")
      else if (robins)
        x <- FALSE
    }
  }
  else {
    chkchar(x, nchar = 1)
    ##
    if (length(x) == 1 && is.logical(x)) {
      if (x) {
        if (tool == "user-defined")
          x <- seq_along(cat)
        else if (rob)
          x <- c("+", "?", "-")
        else
          x <- FALSE
      }
    }
    else 
      chklength(x, length(cat),
                text =
                  paste0("Wrong number of RoB symbols ",
                         "(argument '", varname, "')."))
  }
  ##
  x
}


definecol <- function(avail, x, cat, tool, warn) {
  
  if (!avail)
    return(NULL)
  
  varname <- deparse(substitute(x))
  ##
  robins <- tolower(substring(tool, 1, 6)) == "robins"
  rob <- tolower(substring(tool, 1, 3)) == "rob" & !robins
  
  if (is.null(x)) {
    if (tool == "user-defined")
      x <- seq_len(length(cat))
    else {
      if (rob)
        x <- c("green", "yellow", "red")
      else if (robins)
        x <-  c("green", "yellow", "red", "darkred", "darkgrey")
    }
  }
  else
    chklength(x, length(cat),
              text =
                paste0("Wrong number of RoB colours ",
                       "(argument '", varname, "')."))
  ##
  x
}
