cathet <- function(k,
                   print.tau, text.tau2, tau, digits.tau2, big.mark,
                   print.H, H, lowH, uppH, digits.H,
                   print.I2, print.ci.I2, text.I2,
                   I2, lowI2, uppI2, digits.I2,
                   print.Rb, text.Rb, Rb, lowRb, uppRb
                   ) {
  
  
  pasteCI <- function(lower, upper, digits, big.mark, char = "")
    paste0(" ",
           formatCI(paste0(formatN(lower, digits, big.mark = big.mark), char),
                    paste0(formatN(upper, digits, big.mark = big.mark), char)))
  
  
  cat(
    paste0(
      if (print.tau)
        formatPT(tau^2,
                 lab = TRUE, labval = text.tau2,
                 digits = digits.tau2,
                 lab.NA = "NA",
                 big.mark = big.mark),
      ##
      if (print.H)
        paste0(if (print.tau)
                 "; ",
               "H = ",
               if (is.na(H))
                 "NA"
               else
                 formatN(H, digits.H, "NA", big.mark = big.mark),
               if (k > 2 & !(is.na(lowH) | is.na(uppH)))
                 pasteCI(lowH, uppH, digits.H, big.mark)
               ),
      ##
      if (print.I2)
        paste0(if (print.tau | print.H)
                 "; ",
               text.I2, " = ",
               if (is.na(I2))
                 "NA"
               else
                 paste0(formatN(I2, digits.I2), "%"),
               if (print.ci.I2)
                 pasteCI(lowI2, uppI2, digits.I2, big.mark, "%")
               ),
      ##
      if (print.Rb)
        paste0(if (print.tau | print.H | print.I2)
                 "; ",
               text.Rb, " = ",
               if (is.na(Rb))
                 "NA"
               else
                 paste0(formatN(Rb, digits.I2, big.mark = big.mark), "%"),
               if (k > 2 & !(is.na(lowRb) | is.na(uppRb)))
                 pasteCI(lowRb, uppRb, digits.I2, big.mark, "%")
               ),
      ##
      if (print.tau | print.H | print.I2 | print.Rb)
        "\n"
    )
  )
  
  
  invisible(NULL)
}
