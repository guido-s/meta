format.NA <- function(x, digits = 2, text.NA = "--") {
  outdec <- options()$OutDec
  
  res <- format(ifelse(is.na(x),
                       text.NA,
                       formatC(x, decimal.mark = outdec,
                               format = "f", digits = digits)
                       )
                )
  ##
  res <-  rmSpace(res, end = TRUE)
  ##
  res
}
