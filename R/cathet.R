cathet <- function(k,
                   tau2, lower.tau2, upper.tau2,
                   print.tau2, print.tau2.ci, text.tau2, digits.tau2,
                   tau, lower.tau, upper.tau,
                   print.tau, print.tau.ci, text.tau, digits.tau,
                   sign.lower.tau, sign.upper.tau,
                   I2, lowI2, uppI2, 
                   print.I2, print.I2.ci, text.I2, digits.I2,
                   H, lowH, uppH,
                   print.H, digits.H,
                   Rb, lowRb, uppRb,
                   print.Rb, text.Rb,
                   big.mark) {
  
  
  if (is.null(lower.tau2))
    lower.tau2 <- NA
  if (is.null(upper.tau2))
    upper.tau2 <- NA
  if (is.null(lower.tau))
    lower.tau <- NA
  if (is.null(upper.tau))
    upper.tau <- NA
  
  
  cat(
    paste0(
      " ",
      if (print.tau2)
        paste0(formatPT(tau^2,
                        lab = TRUE, labval = text.tau2,
                        digits = digits.tau2,
                        lab.NA = "NA",
                        big.mark = big.mark),
               if (print.tau2.ci)
                 pasteCI(lower.tau2, upper.tau2, digits.tau2, big.mark,
                         sign.lower.tau, sign.upper.tau)),
      ##
      if (print.tau)
        paste0(if (print.tau2)
                 "; ",
               formatPT(tau,
                        lab = TRUE, labval = text.tau,
                        digits = digits.tau,
                        lab.NA = "NA",
                        big.mark = big.mark),
               if (print.tau.ci)
                 pasteCI(lower.tau, upper.tau, digits.tau, big.mark,
                         sign.lower.tau, sign.upper.tau)),
      ##
      if (print.I2)
        paste0(if (print.tau2 | print.tau)
                 ";",
               if (print.tau2.ci | print.tau.ci)
                 "\n",
               if (print.tau2 | print.tau)
                 " ",
               text.I2, " = ",
               if (is.na(I2))
                 "NA"
               else
                 paste0(formatN(I2, digits.I2), "%"),
               if (print.I2.ci)
                 pasteCI(lowI2, uppI2, digits.I2, big.mark, unit = "%")
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
                 pasteCI(lowRb, uppRb, digits.I2, big.mark, unit = "%")
               ),
      ##
      if (print.tau2 | print.tau | print.I2 | print.H | print.Rb)
        "\n"
    )
  )
  
  
  invisible(NULL)
}
