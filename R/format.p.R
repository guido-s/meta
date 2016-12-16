format.p <- function(p, lab = FALSE, labval = "p", noblanks = FALSE,
                     digits = 4, zero = TRUE) {
  if (is.null(p))
    return("")
  outdec <- options()$OutDec
  if (lab)
    res <- format(ifelse(is.na(p), paste(labval, "= --"),
                         ifelse(p < 1 / 10^digits,
                                paste(labval, " < 0", outdec,
                                      paste(rep("0", digits - 1), collapse = ""),
                                      "1", sep = ""),
                                paste(paste(labval, " = "),
                                      formatC(round(p, digits), decimal.mark = outdec,
                                              format = "f", digits = digits)
                                      )
                                )
                         )
                  )
  else
    res <- format(ifelse(is.na(p), "      --",
                         ifelse(p < 1 / 10^digits,
                                paste("< 0", outdec,
                                      paste(rep("0", digits - 1), collapse = ""),
                                      "1", sep = ""),
                                paste(" ", formatC(round(p, digits), decimal.mark = outdec,
                                                   format = "f", digits = digits)
                                      )
                                )
                         )
                  )
  ##
  if (noblanks)
    res <- gsub(" ", "", res)
  if (!zero)
    res <- gsub("0\\.", "\\.", res)
  ##
  res
}
