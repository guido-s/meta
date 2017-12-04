forest.metabind <- function(x,
                            leftcols,
                            leftlabs,
                            rightcols = c("effect", "ci"),
                            rightlabs,
                            ##
                            overall = FALSE,
                            subgroup = FALSE,
                            overall.hetstat = FALSE,
                            ##
                            lab.NA = "",
                            ##
                            digits = gs("digits.forest"),
                            digits.se = gs("digits.se"),
                            digits.zval = gs("digits.zval"),
                            digits.pval = max(gs("digits.pval") - 2, 2),
                            digits.pval.Q = max(gs("digits.pval.Q") - 2, 2),
                            digits.Q = gs("digits.Q"),
                            digits.tau2 = gs("digits.tau2"),
                            digits.I2 = max(gs("digits.I2") - 1, 0),
                            ##
                            scientific.pval = gs("scientific.pval"),
                            big.mark = gs("big.mark"),
                            ##
                            smlab,
                            ...) {
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  chkclass(x, "metabind")
  ##
  chklogical(overall)
  chklogical(subgroup)
  chklogical(overall.hetstat)
  ##
  chkchar(lab.NA)
  ##
  chknumeric(digits, min = 0, single = TRUE)
  chknumeric(digits.se, min = 0, single = TRUE)
  chknumeric(digits.zval, min = 0, single = TRUE)
  chknumeric(digits.pval, min = 1, single = TRUE)
  chknumeric(digits.pval.Q, min = 1, single = TRUE)
  chknumeric(digits.Q, min = 0, single = TRUE)
  chknumeric(digits.tau2, min = 0, single = TRUE)
  chknumeric(digits.I2, min = 0, single = TRUE)
  ##
  chklogical(scientific.pval)
  ##
  addargs <- names(list(...))
  ##
  idx <- charmatch(tolower(addargs), "hetstat", nomatch = NA)
  if (!is.na(idx) && length(idx) > 0)
    if (list(...)[["hetstat"]])
      stop("Argument 'hetstat' must be FALSE for metabind objects.")
  ##
  idx <- charmatch(tolower(addargs), "comb.fixed", nomatch = NA)
  if (!is.na(idx) && length(idx) > 0)
    stop("Argument 'comb.fixed' cannot be used with metabind objects.")
  ##
  idx <- charmatch(tolower(addargs), "comb.random", nomatch = NA)
  if (!is.na(idx) && length(idx) > 0)
    stop("Argument 'comb.random' cannot be used with metabind objects.")
  
  
  x$k.w <- x$k.all.w <- as.vector(table(x$data$name)[unique(x$data$name)])
  
  
  missing.leftcols <- missing(leftcols)
  ##
  if (missing.leftcols) {
    leftcols <- c("studlab", "k")
    if (any(x$is.subgroup))
      leftcols <- c(leftcols, "pval.Q.b")
  }
  ##
  if (!is.logical(rightcols)) {
    if ("pval.Q.b" %in% rightcols & "pval.Q.b" %in% leftcols)
      leftcols <- leftcols[leftcols != "pval.Q.b"]
    ##
    if ("k" %in% rightcols & "k" %in% leftcols)
      leftcols <- leftcols[leftcols != "k"]
  }
  ##
  label.studlab <- ifelse(any(x$is.subgroup), "Subgroup", "Meta-Analysis")
  ##
  if (missing(leftlabs)) {
    leftlabs <- rep(NA, length(leftcols))
    leftlabs[leftcols == "studlab"] <- label.studlab
    leftlabs[leftcols == "k"] <- "Number of\nStudies"
    leftlabs[leftcols == "Q.b"] <- "Q.b"
    leftlabs[leftcols == "tau2"] <- "Between-study\nvariance"
    leftlabs[leftcols == "pval.Q.b"] <- "Interaction\nP-value"
    leftlabs[leftcols == "I2"] <- "I2"
  }
  ##
  if (missing(rightlabs)) {
    rightlabs <- rep(NA, length(rightcols))
    rightlabs[rightcols == "studlab"] <- label.studlab
    rightlabs[rightcols == "k"] <- "Number of\nStudies"
    rightlabs[rightcols == "Q.b"] <- "Q.b"
    rightlabs[rightcols == "tau2"] <- "Between-study\nvariance"
    rightlabs[rightcols == "pval.Q.b"] <- "Interaction\nP-value"
    rightlabs[rightcols == "I2"] <- "I2"
  }


  ##
  ## Set test for interaction
  ##
  if (any(x$is.subgroup)) {
    if (x$comb.fixed) {
      x$data$Q.b <- x$data$Q.b.fixed
      x$data$pval.Q.b <- x$data$pval.Q.b.fixed
    }
    else {
      x$data$Q.b <- x$data$Q.b.random
      x$data$pval.Q.b <- x$data$pval.Q.b.random
    }
  }
  
  
  ##
  ## Round and round ...
  ##
  x$TE <- round(x$TE, digits)
  x$TE.fixed <- round(x$TE.fixed, digits)
  x$TE.random <- round(x$TE.random, digits)
  if (!is.null(x$byvar)) {
    x$TE.fixed.w <- round(x$TE.fixed.w, digits)
    x$TE.random.w <- round(x$TE.random.w, digits)
  }
  ##
  x$lower <- round(x$lower, digits)
  x$lower.fixed <- round(x$lower.fixed, digits)
  x$lower.random <- round(x$lower.random, digits)
  x$lower.predict <- round(x$lower.predict, digits)
  ##
  x$upper <- round(x$upper, digits)
  x$upper.fixed <- round(x$upper.fixed, digits)
  x$upper.random <- round(x$upper.random, digits)
  x$upper.predict <- round(x$upper.predict, digits)
  ##
  x$seTE <- round(x$seTE, digits.se)
  ##
  x$zval <- round(x$zval, digits.zval)
  ##
  if (any(x$is.subgroup)) {
    x$data$Q.b <- round(x$data$Q.b, digits.Q)
    x$data$pval.Q.b <- formatPT(x$data$pval.Q.b,
                                lab = FALSE,
                                digits = digits.pval.Q,
                                scientific = scientific.pval,
                                lab.NA = lab.NA)
    x$data$pval.Q.b <- rmSpace(x$data$pval.Q.b)
  }
  ##
  x$Q <- round(x$Q, digits.Q)
  ##
  x$data$tau2 <- formatPT(x$data$tau2, digits = digits.tau2,
                          big.mark = big.mark,
                          lab = FALSE, lab.NA = lab.NA)
  ##
  I2.na <- is.na(x$data$I2)
  x$data$I2 <- formatN(round(100 * x$data$I2, digits.I2),
                       digits.I2, "")
  x$data$lower.I2 <- formatN(round(100 * x$data$lower.I2,
                                   digits.I2),
                             digits.I2, "")
  x$data$upper.I2 <- formatN(round(100 * x$data$upper.I2,
                                   digits.I2),
                             digits.I2, "")
  ##
  x$data$I2 <- paste(x$data$I2,
                     ifelse(I2.na, "", "%"),
                     sep = "")
  x$data$lower.I2 <- paste(x$data$lower.I2,
                           ifelse(I2.na, "", "%"),
                           sep = "")
  x$data$upper.I2 <- paste(x$data$upper.I2,
                           ifelse(I2.na, "", "%"),
                           sep = "")
  ##
  x$data$tau2 <- rmSpace(x$data$tau2)
  ##
  x$data$I2 <- rmSpace(x$data$I2)
  x$data$lower.I2 <- rmSpace(x$data$lower.I2)
  x$data$upper.I2 <- rmSpace(x$data$upper.I2)
  
  
  if (missing(smlab))
    smlab <- paste(if (x$comb.fixed)
                     "Fixed Effect Model"
                   else
                     "Random Effects Model",
                   "\n(",
                   xlab(x$sm, x$backtransf),
                   ")", sep = "")
  
  
  class(x) <- "meta"
  
  
  forest(x,
         leftcols = leftcols, leftlabs = leftlabs,
         rightcols = rightcols, rightlabs = rightlabs,
         overall = overall, subgroup = subgroup, overall.hetstat = overall.hetstat,
         hetstat = FALSE,
         lab.NA = lab.NA, smlab = smlab,
         ##
         digits = digits,
         digits.se = digits.se,
         digits.zval = digits.zval,
         digits.pval = digits.pval,
         digits.pval.Q = digits.pval.Q,
         digits.Q = digits.Q,
         digits.tau2 = digits.tau2,
         digits.I2 = digits.I2,
         ##
         scientific.pval = scientific.pval,
         big.mark = big.mark,
         ...)
  
  
  invisible(NULL)
}
