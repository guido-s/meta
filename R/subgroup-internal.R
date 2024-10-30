extrVec <- function(x, varname, levs, first = FALSE) {
  res <- sapply(x, "[", varname)
  ##
  if (max(sapply(res, length)) == 0)
    return(NULL)
  ##
  if (first)
    res <- lapply(res, first)
  ##
  if (length(res[[1]]) == 1) {
    res <- unlist(res)
    names(res) <- levs
    return(res)
  }
  ##
  print(res)
  stop("List element '", varname, "' does not contain single values.",
         call. = FALSE)
}

extrMat <- function(x, varname, clab, rlab) {
  res <- sapply(x, "[", varname)
  ##
  if (max(sapply(res, length)) == 0)
    return(NULL)
  ##
  maxlen <- max(sapply(res, length))
  res <- lapply(res, addNAs, max = maxlen)
  ##
  if (length(res[[1]]) == 1) {
    res <- unlist(res)
    names(res) <- clab
    return(res)
  }
  else if (length(res[[1]]) == length(rlab)) {
    res <- matrix(unlist(res), nrow = length(rlab))
    rownames(res) <- rlab
    colnames(res) <- clab
    return(t(res))
  }
  ##
  print(res)
  stop("Wrong number of labels for list element '", varname, "'.",
       call. = FALSE)
}

addNAs <- function(x, max) {
  if (length(x) == 0)
    return(rep(NA, max))
  else if (length(x) == 1 && length(x) < max)
    return(c(x, rep(NA, max - 1)))
  else
    return(x)
}

first <- function(x)
  x[1]

subsetVar <- function(x, subset) {
  if (is.null(x))
    return(x)
  else if (length(x) == 1)
    return(x)
  else
    return(x[subset])
}
