formatPT <- function(x, lab = FALSE, labval = "p", noblanks = FALSE,
                     digits = 4, zero = TRUE, scientific = FALSE,
                     lab.NA = "--", big.mark = "",
                     JAMA = FALSE) {
  
  if (is.null(x))
    return("")
  
  outdec <- options()$OutDec
  
  n.zeros <- digits - 1
  n.zeros[n.zeros < 0] <- 0
  
  if (!scientific) {
    if (lab) {
      if (!JAMA)
        res <- format(ifelse(is.na(x) | is.nan(x),
                             paste(labval, "=", lab.NA),
                      ifelse(x == 0,
                             paste(labval, "= 0"),
                      ifelse(x < 1 / 10^digits,
                             paste0(labval, " < 0", outdec,
                                    paste(rep("0",
                                             n.zeros), collapse = ""),
                                    "1"),
                             paste(paste(labval, "="),
                                   formatC(round(x, digits),
                                           decimal.mark = outdec,
                                           big.mark = big.mark,
                                           format = "f", digits = digits)
                                   )
                             )
                      )
                      )
                      )
      else
        res <- format(ifelse(is.na(x) | is.nan(x),
                             paste(labval, "=", lab.NA),
                      ifelse(x < 0.001,
                             paste0(labval, " < 0", outdec,
                                    paste(rep("0", 2), collapse = ""), "1"),
                      ifelse(x >= 0.001 & x < 0.01,
                             paste(paste(labval, "="),
                                   formatC(x,
                                           decimal.mark = outdec,
                                           big.mark = big.mark,
                                           format = "f", digits = 3)),
                      ifelse(x >= 0.01 & x <= 0.99,
                             paste(paste(labval, "="),
                                   formatC(x,
                                           decimal.mark = outdec,
                                           big.mark = big.mark,
                                           format = "f", digits = 2)),
                             paste(paste(labval, ">"),
                                   formatC(0.99,
                                           decimal.mark = outdec,
                                           big.mark = big.mark,
                                           format = "f", digits = 2)))
                      )
                      )
                      )
                      )
      
    }
    else {
      if (!JAMA)
        res <- format(ifelse(is.na(x) | is.nan(x),
                             lab.NA,
                      ifelse(x == 0,
                             0,
                      ifelse(x < 1 / 10^digits,
                             paste0("< 0", outdec,
                                    paste(rep("0", n.zeros), collapse = ""),
                                    "1"),
                             formatC(round(x, digits),
                                     decimal.mark = outdec,
                                     big.mark = big.mark,
                                     format = "f", digits = digits)
                             )
                      )
                      ),
                      justify = "right")
      else
        res <- format(ifelse(is.na(x) | is.nan(x),
                             lab.NA,
                      ifelse(x < 0.001,
                             paste0("< 0", outdec,
                                    paste(rep("0", 2), collapse = ""), "1"),
                      ifelse(x >= 0.001 & x < 0.01,
                             formatC(x,
                                     decimal.mark = outdec,
                                     big.mark = big.mark,
                                     format = "f", digits = 3),
                      ifelse(x >= 0.01 & x <= 0.99,
                             formatC(x,
                                     decimal.mark = outdec,
                                     big.mark = big.mark,
                                     format = "f", digits = 2),
                             paste(">",
                                   formatC(0.99,
                                           decimal.mark = outdec,
                                           big.mark = big.mark,
                                           format = "f", digits = 2)))
                      )
                      )
                      ),
                      justify = "right")
    }
  }
  else {
    if (lab)
      res <- format(ifelse(is.na(x) | is.nan(x),
                           paste(labval, "=", lab.NA),
                           paste(labval, "=",
                                 formatC(x, decimal.mark = outdec,
                                         big.mark = big.mark,
                                         format = "e", digits = digits)
                                 )
                           )
                    )
    else
      res <- formatC(x, decimal.mark = outdec,
                     big.mark = big.mark, format = "e", digits = digits)
  }
  ##
  if (noblanks)
    res <- gsub(" ", "", res)
  if (!zero)
    res <- gsub("0\\.", "\\.", res)
  ##
  ## Treat NaNs as NAs
  ##
  res[grep("NaN", res)] <- lab.NA
  
  res
}
