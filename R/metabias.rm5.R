metabias.rm5 <- function(x, comp.no, outcome.no,
                         method.bias="linreg",
                         method.bias.binary=method.bias,
                         method.bias.or="score",
                         k.min=10, ...){
  
  
  ##
  ##
  ## (1) Check for rm5 object
  ##
  ##
  chkclass(x, "rm5")
  
  
  if (missing(comp.no))
    comp.no <- unique(x$comp.no)
  
  n <- 1
  for (i in comp.no){
    if (missing(outcome.no))
      jj <- unique(x$outcome.no[x$comp.no==i])
    else
      jj <- outcome.no
    for (j in jj){
      ##
      if (n>1)
        cat("\n*****\n\n")
      n <- n+1
      ##
      m1 <- metacr(x, i, j)
      ##
      if (inherits(m1, "metabin")){
        if (m1$sm=="OR")
            print(metabias(m1, k.min=k.min, method.bias=method.bias.or))
        else 
          print(metabias(m1, k.min=k.min, method.bias=method.bias.binary))
      }
      else
        print(metabias(m1, k.min=k.min, method.bias=method.bias))
    }
  }
  
  invisible(NULL)
}

