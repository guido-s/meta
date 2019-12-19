format.p <- function(x, lab = FALSE, labval = "p", noblanks = FALSE,
                     digits = 4, zero = TRUE, scientific = FALSE,
                     lab.NA = "--", big.mark = "") {
  
  warning("Use of function format.p() from R package meta is deprecated; ",
          "use instead formatPT().")
  
  if (is.null(x))
    return("")
  
  outdec <- options()$OutDec
  
  n.zeros <- digits - 1
  n.zeros[n.zeros < 0] <- 0
  
  if (!scientific) {
    if (lab)
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
