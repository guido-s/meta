setcat <- function(x, labels) {
  if (is.null(labels))
    return(x)
  ##
  x <- setchar(x, labels)
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

setleg <- function(x) {
  domains <- attr(x, "domains")
  ##
  vars <- names(x)
  vars <- vars[!(vars %in% c("Study", "Weight"))]
  vars[vars == "Overall"] <- "O"
  ##
  paste0("(", vars, ") ", domains)
}
