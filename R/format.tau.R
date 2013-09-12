format.tau <- function(tau, lab=FALSE, labval="tau", noblanks=FALSE){
  labval <- labval
  if (lab)
    res <- format(ifelse(is.na(tau), paste(labval, "= --"),
                         ifelse(tau > 0 & tau < 0.0001,
                                paste(labval, "< 0.0001"),
                                paste(paste(labval, "="),
                                      formatC(round(tau, 4))
                                      )
                                )
                         )
                  )
  else
    res <- format(ifelse(is.na(tau), "--",
                         ifelse(tau > 0 & tau < 0.0001, "< 0.0001",
                                paste(" ", formatC(round(tau, 4)))
                                )
                         )
                  )
  ##
  if (noblanks)
    res <- gsub(" ", "", res)
  ##
  res
}
