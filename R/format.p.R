format.p <- function(p, lab=FALSE, labval="p", noblanks=FALSE,
                     digits=4){
  if (is.null(p))
    return("")
  outdec <- options()$OutDec
  if (lab)
    res <- format(ifelse(is.na(p), paste(labval, "= --"),
                         ifelse(p < 0.0001,
                                paste(labval, " < 0", outdec, "0001", sep=""),
                                paste(paste(labval, "="),
                                      formatC(round(p, digits), decimal.mark=outdec)
                                      )
                                )
                         )
                  )
  else
    res <- format(ifelse(is.na(p), "      --",
                         ifelse(p < 0.0001,
                                paste("< 0", outdec, "0001", sep=""),
                                paste(" ", formatC(round(p, digits), decimal.mark=outdec)
                                      )
                                )
                         )
                  )
  ##
  if (noblanks)
    res <- gsub(" ", "", res)
  ##
  res
}
