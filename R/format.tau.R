format.tau <- function(tau, lab=FALSE, labval="tau", noblanks=FALSE, digits=4){
  outdec <- options()$OutDec
  if (lab)
    res <- format(ifelse(is.na(tau), paste(labval, "= --"),
                         ifelse(tau > 0 & tau < 0.0001,
                                paste(labval, " < 0", outdec, "0001", sep=""),
                                paste(paste(labval, "="),
                                      formatC(round(tau, digits), decimal.mark=outdec)
                                      )
                                )
                         )
                  )
  else
    res <- format(ifelse(is.na(tau), "--",
                         ifelse(tau > 0 & tau < 0.0001,
                                paste("< 0", outdec, "0001", sep=""),
                                paste(" ", formatC(round(tau, digits), decimal.mark=outdec))
                                )
                         )
                  )
  ##
  if (noblanks)
    res <- gsub(" ", "", res)
  ##
  res
}
