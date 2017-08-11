metabias.default <- function(x, seTE,
                             method.bias = "linreg",
                             plotit = FALSE, correct = FALSE,
                             k.min = 10, ...) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  k.All <- length(x)
  ##
  chknumeric(x)
  chknumeric(seTE)
  ##
  fun <- "metabias"
  chklength(seTE, k.All, fun)
  ##
  method.bias <- setchar(method.bias, c("rank", "linreg", "mm"))
  
  
  ##
  ##
  ## (2) Do meta-analysis
  ##
  ##
  m <- metagen(x, seTE)
  
  
  ##
  ##
  ## (3) Conduct test for funnel plot asymmetry
  ##
  ##
  res <- metabias(m, method.bias = method.bias,
                  plotit = plotit, correct = correct,
                  k.min = k.min, ...)
  
  
  res
}
