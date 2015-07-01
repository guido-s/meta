metareg <- function(x, formula,
                    method.tau=x$method.tau,
                    hakn=x$hakn,
                    level.comb=x$level.comb,
                    intercept=TRUE,
                    ...){
  
  if ("data" %in% names(list(...))){
    warning("Please note, argument 'data' has been renamed to 'x' in version 3.0-0 of R package meta (see help page of R function metareg). No meta-regression conducted.")
    return(invisible(NULL))
  }
  
  if (is.call(x)){
    warning("Please note, first two arguments of R function metareg have been interchanged in version 3.0-0 of R package meta. No meta-regression conducted.")
    return(invisible(NULL))
  }
  
  
  ##
  ##
  ## (1) Check for meta object
  ##
  ##
  chkclass(x, "meta")
  
  
  ##
  ## Check whether R package metafor is installed
  ##
  is.installed.package("metafor", "'metareg'")
  
  
  if (missing(formula))
    if (!is.null(x$data$.byvar))
      if (intercept)
        formula <- as.call(~.byvar)
      else
        formula <- as.call(~.byvar-1)
        
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
  method.tau <- setchar(method.tau,
                        c("DL", "PM", "REML", "ML", "HS", "SJ", "HE", "EB", "FE"))
  ##
  if (method.tau=="PM"){
    warning("Meta-regresion method not available for method.tau=\"PM\". Using REML method instead (method.tau=\"REML\").")
    method.tau <- "REML"
  }
  ##
  chklogical(hakn)
  ##
  chklevel(level.comb)
  chklogical(intercept)
  
  if (is.null(x$data)){
    warning("Necessary data not available. Please, recreate meta-analysis object without option 'keepdata=FALSE'.")
    return(invisible(NULL))
  }
  
  if (!is.null(x$subset))
    dataset <- x$data[x$subset,]
  else
    dataset <- x$data
  
  res <- metafor::rma.uni(yi=x$TE,
                          sei=x$seTE,
                          data=dataset,
                          mods=formula, method=method.tau,
                          knha=hakn, level=100*level.comb,
                          ...)

  res$.meta <- list(x=x,
                    formula=formula,
                    method.tau=method.tau,
                    hakn=hakn,
                    level.comb=level.comb,
                    intercept=intercept,
                    dots=list(...),
                    call=match.call(),
                    version=packageDescription("meta")$Version,
                    version.metafor=packageDescription("metafor")$Version)

  class(res) <- c("metareg", class(res))
  
  res
}
