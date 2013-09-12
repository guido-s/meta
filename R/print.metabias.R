print.metabias <- function(x, ...){
  
  if (!inherits(x, "metabias"))
    stop("Argument 'x' must be an object of class \"metabias\"")
  
  ## Print information for meta-analysis from Cochrane Review
  ##
  crtitle(x)
  
  class(x) <- "htest"
  ##
  if (length(x$p.value)!=0)
    print(x)
  else{
    ##
    ## Check whether number of studies is too small:
    ##
    if (length(x$k)!=0 & length(x$k.min!=0)){
      if (x$k <= x$k.min){
        if (x$k <= 2)
          warning("Number of studies (k=",  x$k,
                  ") too small to test for small study effects.")
        else
          warning("Number of studies (k=",  x$k,
                  ") too small to test for small study effects (k.min=",
                  x$k.min, "). Change parameter 'k.min' if appropriate.\n")
      }
    }
    ##
    ## Check whether meta-analysis has subgroups:
    ##
    if (length(x$subgroup)!=0)
      warning("No test for small study effects conducted for meta-analysis with subgroups.")
  }
  
  invisible(NULL)
}
