# Auxiliary functions to format print output
#
# Package: meta
# Author: Guido Schwarzer <guido.schwarzer@uniklinik-freiburg.de>
# License: GPL (>= 2)
#

bylabel <- function(subgroup.name, subgroup.levels, print.subgroup.name,
                    sep.subgroup, big.mark = "") {
  if (print.subgroup.name) {
    if (length(subgroup.name) == 0 || subgroup.name == "")
      res <- format(subgroup.levels, big.mark = big.mark)
    else
      res <- paste0(subgroup.name, sep.subgroup,
                    format(subgroup.levels, big.mark = big.mark))
  }
  else
    res <- format(subgroup.levels, big.mark = big.mark)
  #
  res
}

crtitle <- function(x, prefix = NULL, forest = FALSE) {
  avail <- function(value)
    !(is.null(value) || identical(value, ""))
  #
  title_line <- function(value, prefix = NULL, gap = 0) {
    if (!avail(value))
      return(character(0))
    #
    if (nchar(value) > tl)
      value <- paste0(substring(value, 1, max(1, tl - 4)), " ...")
    #
    paste0(prefix, if (!is.null(prefix)) strrep(" ", gap), value)
  }
  
  tl <- if (forest) 666 else getOption("width") - 12
  
  # Define prefixes and gaps
  #
  prefix <- replaceNULL(prefix, avail(x$complab) || avail(x$outclab))
  #
  if (!prefix) {
    prefixes <- rep(NULL, 3)
    gaps <- rep(0, 3)
  }
  else {
    prefixes <- c("Review:", "Comparison:", "Outcome:")
    #
    gaps <- if (forest)
      rep(1, 3)
    else {
      c(max(1, 5 * avail(x$complab), 2 * avail(x$outclab)),
        1,
        max(1, 4 * avail(x$complab)))
    }
  }
  
  # Title lines
  #
  lines <- c(title_line(x$title, prefixes[1], gaps[1]),
             title_line(x$complab, prefixes[2], gaps[2]),
             title_line(x$outclab, prefixes[3], gaps[3]))
  #
  text <- paste(lines, collapse = "\n")
  #
  if (forest)
    text
  else if (length(lines) > 0)
    cat(text, "\n\n", sep = "")
}

formatCI <- function(lower, upper,
                     bracket.left = gs("CIbracket"),
                     separator = gs("CIseparator"),
                     bracket.right,
                     justify.lower = "right",
                     justify.upper = justify.lower,
                     lower.blank = gs("CIlower.blank"),
                     upper.blank = gs("CIupper.blank"),
                     #
                     digits = NULL,
                     lab.NA = gs("lab.NA"),
                     big.mark = "",
                     scale = 1,
                     suffix = "",
                     na.lab = FALSE,
                     ...) {
  
  # Change layout of CIs
  #
  chkchar(bracket.left, length = 1)
  chkchar(separator, length = 1)
  #
  if (missing(bracket.left)) {
    bracktype <- setchar(bracket.left, c("[", "(", "{", ""))
    #
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
    else {
      bracketLeft <- ""
      bracketRight <- ""
    }
    #
    bracket.left <- bracketLeft
  }
  else {
    bracketRight <- switch(bracket.left,
                           "[" = "]",
                           "(" = ")",
                           "{" = "}",
                           "")
  }
  #
  if (!missing(bracket.right))
    chkchar(bracket.right, length = 1)
  else
    bracket.right <- bracketRight
  
  # Identify missing limits before formatting
  #
  missing.lower <- is.na(lower)
  missing.upper <- is.na(upper)
  #
  if (is.character(lower))
    missing.lower <- missing.lower | lower == "NA"
  #
  if (is.character(upper))
    missing.upper <- missing.upper | upper == "NA"
  #
  missing.ci <- missing.lower | missing.upper
  
  # Format limits
  #
  if (is.null(digits)) {
    format.lower <- format(lower, justify = justify.lower)
    format.upper <- format(upper, justify = justify.upper)
  }
  else {
    format.lower <-
      formatN(scale * lower, digits = digits, lab.NA = lab.NA,
              big.mark = big.mark)
    format.upper <-
      formatN(scale * upper, digits = digits, lab.NA = lab.NA,
              big.mark = big.mark)
    #
    format.lower <- paste0(format.lower, suffix)
    format.upper <- paste0(format.upper, suffix)
  }
  
  # Change layout of CIs
  #
  if (!lower.blank)
    format.lower <- rmSpace(format.lower)
  if (!upper.blank)
    format.upper <- rmSpace(format.upper)
  #
  if (separator == "-") {
    format.upper <-
      paste0(ifelse(substring(format.upper, 1, 1) == "-", " ", ""),
             format.upper)
  }
  
  res <- paste0(bracket.left, format.lower,
                separator,
                format.upper, bracket.right)
  #
  if (na.lab)
    res[missing.ci] <- lab.NA
  else
    res[missing.ci] <- ""
  #
  res
}

formatN <- function(x, digits = 2, lab.NA = "--", big.mark = "",
                    format.whole.numbers = TRUE,
                    monospaced = FALSE) {
  
  outdec <- options()$OutDec  

  if (!monospaced) {
    if (format.whole.numbers) {
      res <- format(ifelse(is.na(x),
                           lab.NA,
                           formatC(x, decimal.mark = outdec,
                                   format = "f", digits = digits,
                                   big.mark = big.mark)
                           )
                    )
    }
    else {
      res <- format(ifelse(is.na(x),
                           lab.NA,
                    ifelse(is_wholenumber(x),
                           x,
                           formatC(x, decimal.mark = outdec,
                                   format = "f", digits = digits,
                                   big.mark = big.mark)
                           )
                    )
                    )
    }
  }
  else {
    x <- round(x, digits)
    res <- ifelse(is.na(x),
                  lab.NA,
                  format(x, decimal.mark = outdec, big.mark = big.mark))
  }
  #
  res <- rmSpace(res, end = TRUE)
  #
  res
}

formatPT <- function(x, lab = FALSE, labval = "p", noblanks = FALSE,
                     digits = 4, zero = TRUE, scientific = FALSE,
                     lab.NA = "--", big.mark = "",
                     JAMA = FALSE, trim = NULL) {
  
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
  
  # Remove any whitespace
  #
  if (noblanks)
    res <- gsub(" ", "", res)
  if (!zero)
    res <- gsub("0\\.", "\\.", res)
  
  # Treat NaNs as NAs
  #
  res[grep("NaN", res)] <- lab.NA
  
  # Trim whitespace
  #
  if (!is.null(trim)) {
    res <- trimws(res, which = trim)
  }
  
  res
}

pasteCI <- function(lower, upper, digits, big.mark,
                    sign.lower = "", sign.upper = "", lab.NA = "NA",
                    unit = "")
  paste0(" ",
         formatCI(paste0(sign.lower,
                         formatN(lower, digits, big.mark = big.mark,
                                 lab.NA = lab.NA), unit),
                  paste0(sign.upper,
                         formatN(upper, digits, big.mark = big.mark,
                                 lab.NA = lab.NA), unit)))

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
