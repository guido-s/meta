cathet <- function(k,
                   tau,
                   print.tau2, text.tau2, digits.tau2,
                   print.tau, text.tau, digits.tau,
                   I2, lowI2, uppI2, 
                   print.I2, print.ci.I2, text.I2, digits.I2,
                   H, lowH, uppH,
                   print.H, digits.H,
                   Rb, lowRb, uppRb,
                   print.Rb, text.Rb,
                   big.mark) {
  
  
  pasteCI <- function(lower, upper, digits, big.mark, char = "")
    paste0(" ",
           formatCI(paste0(formatN(lower, digits, big.mark = big.mark), char),
                    paste0(formatN(upper, digits, big.mark = big.mark), char)))
  
  
  cat(
    paste0(
      if (print.tau2)
        formatPT(tau^2,
                 lab = TRUE, labval = text.tau2,
                 digits = digits.tau2,
                 lab.NA = "NA",
                 big.mark = big.mark),
      ##
      if (print.tau)
        paste0(if (print.tau2)
                 "; ",
               formatPT(tau,
                        lab = TRUE, labval = text.tau,
                        digits = digits.tau,
                        lab.NA = "NA",
                        big.mark = big.mark)),
      ##
      if (print.I2)
        paste0(if (print.tau2 | print.tau)
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
      if (print.H)
        paste0(if (print.tau2 | print.tau | print.I2)
                 "; ",
               "H = ",
               if (is.na(H))
                 "NA"
               else
                 formatN(H, digits.H, "NA", big.mark = big.mark),
               if (!(is.na(lowH) | is.na(uppH)))
                 pasteCI(lowH, uppH, digits.H, big.mark)
               ),
      ##
      if (print.Rb)
        paste0(if (print.tau2 | print.tau | print.I2 | print.H)
                 ";\n",
               text.Rb, " = ",
               if (is.na(Rb))
                 "NA"
               else
                 paste0(formatN(Rb, digits.I2, big.mark = big.mark), "%"),
               if (!(is.na(lowRb) | is.na(uppRb)))
                 pasteCI(lowRb, uppRb, digits.I2, big.mark, "%")
               ),
      ##
      if (print.tau2 | print.tau | print.I2 | print.H | print.Rb)
        "\n"
    )
  )
  
  
  invisible(NULL)
}
