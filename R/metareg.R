metareg <- function(formula, data, method.tau=data$method.tau){
  
  if (!inherits(data, "meta"))
    stop("Argument 'data' must be an object of class \"meta\"")

  if (is.null(method.tau))
    method.tau <- "DL"
  
  ##
  ## Check whether R package metafor is installed
  ##
  is.installed.metafor()
  
  res <- metafor::rma.uni(yi=data$TE, sei=data$seTE,
                          data=data,
                          mods=formula, method=method.tau)
  
  res
}
