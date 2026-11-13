subgroup_result_matrix <- function(x, reference) {
  template <- list2mat(reference)
  #
  if (is.null(template))
    return(NULL)
  #
  if (is.list(reference)) {
    if (is.null(x))
      x <- rep(list(NULL), length(reference))
    #
    if (!is.list(x) || length(x) != length(reference))
      stop("Subgroup result blocks do not match their confidence intervals.",
           call. = FALSE)
    #
    for (i in seq_along(reference)) {
      if (is.null(x[[i]]) && !is.null(reference[[i]])) {
        block <- list2mat(reference[[i]])
        x[i] <- list(matrix(
          NA_real_, nrow = nrow(block), ncol = ncol(block),
          dimnames = dimnames(block)))
      }
    }
  }
  else if (is.null(x))
    return(matrix(NA_real_, nrow = nrow(template), ncol = ncol(template),
                  dimnames = dimnames(template)))
  #
  res <- list2mat(x, reference)
  #
  if (!is.null(res) &&
      nrow(res) == nrow(template) &&
      ncol(res) == 1L && ncol(template) > 1L)
    res <- matrix(rep(res, ncol(template)),
                    nrow = nrow(template), ncol = ncol(template))
  #
  if (!identical(dim(res), dim(template)))
    stop("Subgroup results have incompatible dimensions.",
         call. = FALSE)
  #
  dimnames(res) <- dimnames(template)
  res
}
