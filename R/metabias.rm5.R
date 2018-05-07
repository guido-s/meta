metabias.rm5 <- function(x, comp.no, outcome.no,
                         method.bias = "linreg",
                         method.bias.binary = method.bias,
                         method.bias.or = "score",
                         k.min = 10, ...) {
  
  
  ##
  ##
  ## (1) Check for rm5 object
  ##
  ##
  chkclass(x, "rm5")
  
  
  if (missing(comp.no))
    comp.no <- unique(x$comp.no)
  
  
  n <- 1
  ##
  for (i in comp.no) {
    if (missing(outcome.no))
      jj <- unique(x$outcome.no[x$comp.no == i])
    else
      jj <- outcome.no
    for (j in jj) {
      ##
      m1 <- metacr(x, i, j)
      ##
      if (inherits(m1, "metabin")) {
        if (m1$sm == "OR")
          mb1 <- metabias(m1, k.min = k.min, method.bias = method.bias.or)
        else 
          mb1 <- metabias(m1, k.min = k.min, method.bias = method.bias.binary)
      }
      else
        mb1 <- metabias(m1, k.min = k.min, method.bias = method.bias)
      ##
      if (!is.null(mb1$estimate)) {
        if (n > 1)
          cat("\n*****\n\n")
        print(mb1)
        ##
        n <- n + 1
      }
    }
  }
  
  invisible(NULL)
}
