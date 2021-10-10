## Auxiliary functions to format print output
##
## Package: meta
## Author: Guido Schwarzer <sc@imbi.uni-freiburg.de>
## License: GPL (>= 2)
##
bylabel <- function(subgroup.name, bylevs, print.subgroup.name,
                    sep.subgroup, big.mark = "") {
  if (print.subgroup.name) {
    if (length(subgroup.name) == 0 || subgroup.name == "")
      res <- format(bylevs, big.mark = big.mark)
    else
      res <- paste0(subgroup.name, sep.subgroup,
                    format(bylevs, big.mark = big.mark))
  }
  else
    res <- format(bylevs, big.mark = big.mark)
  ##
  res
}
crtitle <- function(x) {
  tl <- options()$width - 12
  newline <- FALSE
  ##  
  if (!is.null(x$title)) {
    if (x$title != "") {
      newline <- TRUE
      if (nchar(x$title) <= tl)
        cat(paste0("Review:     ", x$title, "\n"))
      else
        cat(paste0("Review:     ", substring(x$title, 1, tl - 4), " ...\n"))
    }
  }
  if (!is.null(x$complab)) {
    if (x$complab != "") {
      newline <- TRUE
      if (nchar(x$complab) <= tl)
        cat(paste0("Comparison: ", x$complab, "\n"))
      else
        cat(paste0("Comparison: ", substring(x$complab, 1, tl - 4), " ...\n"))
    }
  }
  if (!is.null(x$outclab)) {
    if (x$outclab != "") {
      newline <- TRUE
      if (nchar(x$outclab) <= tl)
        cat(paste0("Outcome:    ", x$outclab, "\n"))
      else
        cat(paste0("Outcome:    ", substring(x$outclab, 1, tl - 4), " ...\n"))
    }
  }
  ##
  if (newline)
    cat("\n")
}
format.NA <- function(x, digits = 2, text.NA = "--", big.mark = "") {
  
  warning("Use of function format.NA() from R package meta is deprecated; ",
          "use instead formatN().")
  
  outdec <- options()$OutDec
  
  res <- format(ifelse(is.na(x),
                       text.NA,
                       formatC(x, decimal.mark = outdec,
                               format = "f", digits = digits,
                               big.mark = big.mark)
                       )
                )
  ##
  res <-  rmSpace(res, end = TRUE)
  ##
  res
}
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
formatCI <- function(lower, upper, rmspace = TRUE,
                     bracket.left = gs("CIbracket"),
                     separator = gs("CIseparator"),
                     bracket.right,
                     justify.lower = "right",
                     justify.upper = justify.lower,
                     ...
                     ) {
  
  ## Change layout of CIs
  ##
  bracks <- c("[", "(", "{", "")
  ibracket <- charmatch(bracket.left,
                        bracks,
                        nomatch = NA)
  ##
  if (is.na(ibracket) | ibracket == 0)
    stop("No valid bracket type specified. ",
         "Admissible values: '[', '(', '{', '\"\"'")
  ##
  bracktype <- bracks[ibracket]
  ##
  if (bracktype == "[") {
    bracketLeft <- "["
    bracketRight <- "]"
  }
  else if (bracktype == "(") {
    bracketLeft <- "("
    bracketRight <- ")"
  }
  else if (bracktype == "{") {
    bracketLeft <- "{"
    bracketRight <- "}"
  }
  else if (bracktype == "") {
    bracketLeft <- ""
    bracketRight <- ""
  }
  ##
  if (missing(bracket.left))
    bracket.left <- bracketLeft
  ##
  if (missing(bracket.right))
    bracket.right <- bracketRight
  
  if (rmspace) {
    lower <- rmSpace(lower)
    upper <- rmSpace(upper)
  }
  ##
  res <- ifelse(lower != "NA" & upper != "NA",
                paste0(bracket.left,
                       format(lower, justify = justify.lower),
                       separator,
                       format(upper, justify = justify.upper),
                       bracket.right),
                "")
  ##
  res
}
formatN <- function(x, digits = 2, text.NA = "--", big.mark = "") {
  
  outdec <- options()$OutDec
  
  res <- format(ifelse(is.na(x),
                       text.NA,
                       formatC(x, decimal.mark = outdec,
                               format = "f", digits = digits,
                               big.mark = big.mark)
                       )
                )
  ##
  res <-  rmSpace(res, end = TRUE)
  ##
  res
}
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
p.ci <- function(lower, upper, rmspace = TRUE,
                 bracket.left = gs("CIbracket"),
                 separator = gs("CIseparator"),
                 bracket.right = "]",
                 justify.lower = "right",
                 justify.upper = justify.lower,
                 ...
                 ) {
  
  warning("Use of function p.ci() from R package meta is deprecated; ",
          "use instead formatCI().")
  
  ## Change layout of CIs
  ##
  ibracktype <- charmatch(gs("CIbracket"),
                          c("[", "(", "{", ""), nomatch = NA)
  if (is.na(ibracktype) | ibracktype == 0) {
    warning("No valid bracket type specified globally for R package meta: ",
            gs("CIbracket"),
            "\n  Using default bracket type: '['. See help page on ",
            "R command 'cilayout' for further information.")
    bracktype <- "["
  }
  else
    bracktype <- c("[", "(", "{", "")[ibracktype]
  ##
  if (bracktype == "[") {
    bracketLeft <- "["
    bracketRight <- "]"
  }
  else if (bracktype == "(") {
    bracketLeft <- "("
    bracketRight <- ")"
  }
  else if (bracktype == "{") {
    bracketLeft <- "{"
    bracketRight <- "}"
  }
  else if (bracktype == "") {
    bracketLeft <- ""
    bracketRight <- ""
  }
  ##
  if (missing(bracket.left))
    bracket.left <- bracketLeft
  ##
  if (missing(bracket.right))
    bracket.right <- bracketRight
  ##
  
  if (rmspace) {
    lower <- rmSpace(lower)
    upper <- rmSpace(upper)
  }
  ##
  res <- ifelse(lower != "NA" & upper != "NA",
                paste0(bracket.left,
                       format(lower, justify = justify.lower),
                       separator,
                       format(upper, justify = justify.upper),
                       bracket.right),
                "")
  ##
  res
}
pasteCI <- function(lower, upper, digits, big.mark,
                    sign.lower = "", sign.upper = "", text.NA = "NA",
                    unit = "")
  paste0(" ",
         formatCI(paste0(sign.lower,
                         formatN(lower, digits, big.mark = big.mark,
                                 text.NA = text.NA), unit),
                  paste0(sign.upper,
                         formatN(upper, digits, big.mark = big.mark,
                                 text.NA = text.NA), unit)))
rmSpace <- function(x, end = FALSE, pat = " ") {
  
  if (!end) {
    while (any(substring(x, 1, 1) == pat, na.rm = TRUE)) {
      sel <- substring(x, 1, 1) == pat
      x[sel] <- substring(x[sel], 2)
    }
  }
  else {
    last <- nchar(x)
    
    while (any(substring(x, last, last) == pat, na.rm = TRUE)) {
      sel <- substring(x, last, last) == pat
      x[sel] <- substring(x[sel], 1, last[sel] - 1)
      last <- nchar(x)
    }
  }
  
  x
}
