format.p <- function(p, lab=FALSE, labval="p", noblanks=FALSE){
  if (lab)
    res <- format(ifelse(is.na(p), paste(labval, "= --",),
                         ifelse(p < 0.0001, paste(labval, "< 0.0001"),
                                paste(paste(labval, "="),
                                      formatC(round(p, 4))))))
  else
    res <- format(ifelse(is.na(p), "      --",
                         ifelse(p < 0.0001, "< 0.0001",
                                paste(" ", formatC(round(p, 4))))))
  ##
  if (noblanks)
    res <- gsub(" ", "", res)
  ##
  res
}
