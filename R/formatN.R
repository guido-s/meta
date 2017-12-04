formatN <- function(x, digits = 2, text.NA = "--", big.mark = "") {
  
  outdec <- options()$OutDec
  
  res <- format(ifelse(is.na(x),
                       text.NA,
                       formatC(x, decimal.mark = outdec,
                               format = "f", digits = digits,
                               big.mark = big.mark)
                       )
                )
  ##
  res <-  rmSpace(res, end = TRUE)
  ##
  res
}
