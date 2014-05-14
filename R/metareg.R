metareg <- function(x, formula,
                    method.tau=x$method.tau, ...){
  
  if ("data" %in% names(list(...))){
    warning("Please note, argument 'data' has been renamed to 'x' in version 3.0-0 of R package meta (see help page of R function metareg). No meta-regression conducted.")
    return(invisible(NULL))
  }
  
  if (is.call(x)){
    warning("Please note, first two arguments of R function metareg have been interchanged in version 3.0-0 of R package meta. No meta-regression conducted.")
    return(invisible(NULL))
  }
  
  if (!inherits(x, "meta"))
    stop("Argument 'x' must be an object of class \"meta\"")
  
  if (missing(formula))
    if (!is.null(x$data$.byvar))
      formula <- as.call(~.byvar)
    else{
      warning("No meta-regression conducted as argument 'formula' is missing and no information is provided on subgroup variable, i.e. list element 'byvar' in meta-analysis object 'x' (see help page of R function metareg).")
      return(invisible(NULL))
    }
  else{
    formula.text <- deparse(substitute(formula))
    formula.text <- gsub("~", "", formula.text)
    formula.text <- gsub("\\\"", "", formula.text)
    formula.text <- gsub("\\\'", "", formula.text)
    formula <- as.formula(paste("~", formula.text))
  }
  
  if (is.null(method.tau))
    method.tau <- "DL"
  ##
  if (method.tau=="PM"){
    warning("Meta-regresion method not available for method.tau=\"PM\". Using REML method instead (method.tau=\"REML\").")
    method.tau <- "REML"
  }
  
  ##
  ## Check whether R package metafor is installed
  ##
  is.installed.metafor()
  
  if (is.null(x$data)){
    warning("Necessary data not available. Please, recreate meta-analysis object without option 'keepdata=FALSE'.")
    return(invisible(NULL))
  }
  
  if (!is.null(x$subset))
    dataset <- x$data[x$subset,]
  else
    dataset <- x$data
  
  res <- metafor::rma.uni(yi=x$TE, sei=x$seTE,
                          data=dataset,
                          mods=formula, method=method.tau)
  
  res
}
