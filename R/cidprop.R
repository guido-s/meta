#' Calculate expected proportion of comparable studies with clinically important
#' benefit or harm
#'
#' @description
#' Calculate expected proportion of comparable studies with clinically important
#' benefit or harm which is derived from the prediction interval.
#'
#' @aliases cidprop cidprop.meta
#'
#' @param x An object of class \code{meta}.
#' @param cid A numeric value or vector specifying clinically important
#'   differences (CID) / decision thresholds used to calculate expected
#'   proportions of clinically important benefit or harm (see Details).
#' @param cid.below.null A numeric value or vector specifying CID limits below
#'   the null effect (see Details).
#' @param cid.above.null A numeric value or vector specifying CID limits above
#'   the null effect (see Details).
#' @param label.cid A character string or vector specifying labels for
#'   clinically important differences. Must be of same length as argument
#'   \code{cid}.
#' @param label.cid.below.null A character string or vector specifying labels
#'   for clinically important differences below the null effect. Must be of
#'   same length as argument \code{cid.below.null} (or \code{cid}).
#' @param label.cid.above.null A character string or vector specifying labels
#'   for clinically important differences above the null effect. Must be of
#'   same length as argument \code{cid.above.null} (or \code{cid}).
#' @param small.values A character string specifying whether small
#'   treatment effects indicate a beneficial (\code{"desirable"}) or
#'   harmful (\code{"undesirable"}), can be abbreviated.
#' @param big.mark A character used as thousands separator.
#' @param digits.cid Minimal number of significant digits for
#'   CIDs / decision thresholds, see \code{print.default}.
#' @param digits.percent Minimal number of significant digits for
#'   expected proportions, printed as percentages, see \code{print.default}.
#' @param details.methods A logical specifying whether details on
#'   statistical methods should be printed.
#' @param \dots Additional arguments (ignored)
#' 
#' @details
#' Expected proportions of comparable studies with clinically important
#' benefit or harm are derived from the prediction interval in
#' the meta-analysis.
#' 
#' Clinically important benefit or harm can be defined using either argument
#' \code{cid} or \code{cid.below.null} and \code{cid.above.null}.
#' Input for the later arguments will be ignored if argument \code{cid} was
#' specified. In this case, the values of \code{cid.below.null} and
#' \code{cid.above.null} will be equal to
#' \itemize{
#' \item \code{cid} and \code{1 / cid} for ratio measures,
#' \item \code{cid} and \code{-cid} for difference measures.
#' }
#' 
#' Thresholds based on argument \code{cid} will always be symmetric. Asymmetric
#' thresholds can be defined using arguments \code{cid.below.null} and
#' \code{cid.above.null}.
#'
#' @return
#' A list with elements
#' \item{prop.cid.below.null}{Expected proportion of comparable studies below
#'   lower CID(s)}
#' \item{prop.cid.above.null}{Expected proportion of comparable studies above
#'   upper CID(s)}
#' \item{prop.within.cid}{Expected proportion of comparable studies between
#'   lower and upper CID(s)}
#' \item{cid, cid.below.null, cid.above.null, small.values, x}{As defined
#'   above}
#' \item{label.cid, label.cid.below.null, label.cid.above.null}{As defined
#'   above}
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @seealso \code{\link{plot.cidprop}}
#'
#' @examples
#' oldset <- settings.meta(digits.cid = 0)
#' 
#' m <- metagen(1:10 - 3, 1:10, sm = "MD")
#' #
#' pp1 <- cidprop(m, cid = 2)
#' pp1
#' #
#' pp2 <- cidprop(m, cid.below = 0.5, cid.above = 2)
#' pp2
#' #
#' pp3 <- cidprop(m, cid.below = 0.5, cid.above = 2, small.values = "u")
#' pp3
#' 
#' pp4 <- cidprop(m, cid = 1:2, label.cid = c("moderate", "large"))
#' pp4
#' #
#' pp5 <- cidprop(m, cid.below = -1.5, cid.above = 1:2,
#'   label.cid.below = "large", label.cid.above = c("moderate", "large"))
#' pp5
#' 
#' settings.meta(oldset)
#'
#' @rdname cidprop
#' @method cidprop meta
#' @export

cidprop.meta <- function(x,
                          cid = NULL,
                          cid.below.null = NULL, cid.above.null = NULL,
                          #
                          label.cid = "",
                          label.cid.below.null = NULL,
                          label.cid.above.null = NULL,
                          #
                          small.values = "desirable",
                          ...) {
  
  #
  #
  # (1) Check meta-analysis object and arguments
  #
  #
  
  chkclass(x, "meta")
  #
  sm <- x$sm
  is_relative <- is_relative_effect(sm)
  #
  missing.cid <- missing(cid)
  #
  if (is_prop(sm) | is_rate(sm) | is_mean(sm)) {
    if (is_prop(sm))
      outcome <- "proportions"
    else if (is_prop(sm))
      outcome <- "rates"
    else
      outcome <- "means"
    #
    ref <- replaceNULL(x$null.effect)
    #
    if (is.na(ref) & !missing.cid)
      stop("Argument 'cid' can only be used for meta-analysis of single ",
           outcome, " if argument 'null.effect' was specified.",
           call. = FALSE)
  }
  else if (is_relative)
    ref <- 1
  else
    ref <- 0
  #
  small.values <- setchar(small.values, c("desirable", "undesirable"))
  #
  # CID
  #
  avail.cid <- !missing.cid & !is.null(cid) & !all(is.na(cid))
  avail.cid.below.null <-
    !missing(cid.below.null) & !is.null(cid.below.null) &
    !all(is.na(cid.below.null))
  avail.cid.above.null <-
    !missing(cid.above.null) & !is.null(cid.above.null) &
    !all(is.na(cid.above.null))
  #
  if (!avail.cid & !avail.cid.below.null & !avail.cid.above.null)
    stop("At least one decision threshold (argument 'cid', ",
         "'cid.below.null', or 'cid.above.null') must be specified.",
         call. = FALSE)
  #
  if (avail.cid) {
    if (any(is.na(cid)))
      stop("Missing values not allows in argument 'cid'.",
           call. = FALSE)
    #
    if (avail.cid.below.null + avail.cid.above.null == 2)
      warning("Arguments 'cid.below.null' and 'cid.above.null' ignored as ",
              "argument 'cid' is provided.",
              call. = FALSE)
    else if (avail.cid.below.null)
      warning(
        "Argument 'cid.below.null' ignored as argument 'cid' is provided.",
        call. = FALSE)
    else if (avail.cid.above.null)
      warning(
        "Argument 'cid.above.null' ignored as argument 'cid' is provided.",
        call. = FALSE)
    #
    if (any(diff(cid) <= 0))
      stop("Values for argument 'cid' must be increasing.",
           call. = FALSE)
    #
    if (any(cid < ref) & any(cid > ref))
      stop("All values provided for argument 'cid' must be either ",
           "smaller or larger than reference value of ", ref, ".",
           call. = FALSE)
    #
    if (all(cid < ref)) {
      cid.below.null <- cid
      #
      if (is_relative)
        cid.above.null <- rev(1 / cid)
      else
        cid.above.null <- rev(-cid)
    }
    else {
      cid.above.null <- cid
      #
      if (is_relative)
        cid.below.null <- rev(1 / cid)
      else
        cid.below.null <- rev(-cid)
    }
    #
    avail.cid.below.null <- TRUE
    avail.cid.above.null <- TRUE
  }
  #
  if (avail.cid.below.null) {
    chknumeric(cid.below.null)
    #
    cid.below.null.transf <- cid.below.null
    #
    if (is_relative)
      cid.below.null.transf <- log(cid.below.null)
  }
  else {
    cid.below.null <- NA
    cid.below.null.transf <- NA
  }
  #
  if (avail.cid.above.null) {
    chknumeric(cid.above.null)
    #
    cid.above.null.transf <- cid.above.null
    #
    if (is_relative)
      cid.above.null.transf <- log(cid.above.null)
  }
  else {
    cid.above.null <- NA
    cid.above.null.transf <- NA
  }
  #
  # CID labels
  #
  avail.label.cid <- !missing(label.cid) & !is.null(label.cid)
  avail.label.cid.below.null <-
    !missing(label.cid.below.null) & !is.null(label.cid.below.null)
  avail.label.cid.above.null <-
    !missing(label.cid.above.null) & !is.null(label.cid.above.null)
  #
  if (avail.cid & avail.label.cid && length(cid) != length(label.cid))
    stop("Arguments 'cid' and 'label.cid' must be of same length.",
         call. = FALSE)
  #
  if (avail.cid.below.null & avail.label.cid.below.null &&
      length(cid.below.null) != length(label.cid.below.null))
    stop("Arguments 'cid.below.null' and 'label.cid.below.null' must be of same length.",
         call. = FALSE)
  #
  if (avail.cid.above.null & avail.label.cid.above.null &&
      length(cid.above.null) != length(label.cid.above.null))
    stop("Arguments 'cid.above.null' and 'label.cid.above.null' must be of same length.",
         call. = FALSE)
  #
  if (avail.label.cid) {
    if (avail.label.cid.below.null + avail.label.cid.above.null == 2)
      warning("Arguments 'label.cid.below.null' and 'label.cid.above.null' ignored as ",
              "argument 'label.cid' is provided.",
              call. = FALSE)
    else if (avail.label.cid.below.null)
      warning("Argument 'label.cid.below.null' ignored as argument 'label.cid' is provided.",
              call. = FALSE)
    else if (avail.label.cid.above.null)
      warning("Argument 'label.cid.above.null' ignored as argument 'label.cid' is provided.",
              call. = FALSE)
    #
    #
    if (all(cid < ref)) {
      label.cid.below.null <- label.cid
      label.cid.above.null <- rev(label.cid)
    }
    else {
      label.cid.below.null <- rev(label.cid)
      label.cid.above.null <- label.cid
    }
    #
    avail.label.cid.below.null <- TRUE
    avail.label.cid.above.null <- TRUE
  }
  #
  if (!avail.label.cid.below.null)
    label.cid.below.null <-
    if (length(cid.below.null) == 1) "" else rev(seq_along(cid.below.null))
  #
  if (!avail.label.cid.above.null)
    label.cid.above.null <-
    if (length(cid.above.null) == 1) "" else seq_along(cid.above.null)
  
  
  #
  #
  # (2) Only consider results of first random effects meta-analysis
  #
  #
  
  method.predict <- x$method.predict[1]
  #
  TE.random <- x$TE.random[1]
  lower.random <- x$lower.random[1]
  upper.random <- x$upper.random[1]
  lower.predict <- x$lower.predict[1]
  upper.predict <- x$upper.predict[1]
  seTE.predict <- x$seTE.predict[1]
  df.predict <- x$df.predict[1]
  
  
  #
  #
  # (3) Calculate expected proportions
  #
  #
  
  prop.cid.below.null <-
    pt((cid.below.null.transf - TE.random) / seTE.predict, df.predict)
  #
  prop.cid.above.null <-
    pt((cid.above.null.transf - TE.random) / seTE.predict, df.predict,
       lower.tail = FALSE)
  #
  prop.within.cid <-
    1 -
    replaceNA(max(prop.cid.below.null), 0) -
    replaceNA(max(prop.cid.above.null), 0)
  #
  if (is_zero(prop.within.cid))
    prop.within.cid <- 0
  
  
  #
  #
  # (4) Return cidprop object
  #
  #
  
  res <- list(prop.cid.below.null = prop.cid.below.null, 
              prop.cid.above.null = prop.cid.above.null,
              prop.within.cid = prop.within.cid,
              #
              cid = cid,
              cid.below.null = cid.below.null, cid.above.null = cid.above.null,
              #
              label.cid = label.cid,
              label.cid.below.null = label.cid.below.null,
              label.cid.above.null = label.cid.above.null,
              #
              small.values = small.values,
              #
              ref = ref,
              #
              x = x)
  #
  class(res) <- "cidprop"
  #
  res
}


#' @rdname cidprop
#' @export cidprop

cidprop <- function(x, ...)
  UseMethod("cidprop")


#' @rdname cidprop
#' @method print cidprop
#' @export

print.cidprop <- function(x,
                           digits.cid = gs("digits.cid"), digits.percent = 1,
                           big.mark = gs("big.mark"),
                           details.methods = gs("details"),
                           ...) {
  
  #
  #
  # (1) Check arguments
  #
  #
    
  chkclass(x, "cidprop")
  #
  chknumeric(digits.cid, min = 0, length = 1)
  chknumeric(digits.percent, min = 0, length = 1)
  chkchar(big.mark, length = 1)
  chklogical(details.methods)
  
  #
  #
  # (2) Extract information from cidprop object
  #
  #
    
  cid.below.null <- x$cid.below.null
  cid.above.null <- x$cid.above.null
  #
  avail.cid.below.null <- !all(is.na(cid.below.null))
  avail.cid.above.null <- !all(is.na(cid.above.null))
  #
  prop.cid.below.null <- x$prop.cid.below.null
  prop.cid.above.null <- x$prop.cid.above.null
  prop.within.cid <- x$prop.within.cid
  #
  label.cid.below.null <- x$label.cid.below.null
  label.cid.above.null <- x$label.cid.above.null
  #
  # If meta-analysis object is available
  #
  if (!is.null(x$x)) {
    if (!x$x$backtransf) {
      cid.below.null <- transf(cid.below.null, x$x$sm)
      cid.above.null <- transf(cid.above.null, x$x$sm)
    }
    #
    smlab <- smlab(x$x$sm, x$x$backtransf, x$x$pscale, x$x$irscale)
    #
    crtitle(x$x)
  }
  else
    smlab <- ""
  #
  svd <- x$small.values == "desirable"
  
  #
  #
  # (3) Data set with CID information
  #
  #
  
  dat.l <- dat.u <- dat.w <- NULL
  #
  if (avail.cid.below.null) {
    dat.l <-
      data.frame(Threshold = cid.below.null, prop = prop.cid.below.null,
                 label = label.cid.below.null,
                 category =
                   if (svd) "Beneficial effect" else "Harmful effect",
                 sign = "\u2264 ")
    #
    max.cid.below.null <- max(cid.below.null, na.rm = TRUE)
  }
  #
  if (avail.cid.above.null) {
    dat.u <-
      data.frame(Threshold = cid.above.null, prop = prop.cid.above.null,
                 label = label.cid.above.null,
                 category =
                   if (svd) "Harmful effect" else "Beneficial effect",
                 sign = "\u2265 ")
    #
    min.cid.above.null <- min(cid.above.null, na.rm = TRUE)
  }
  #
  if (prop.within.cid > 0) {
    dat.w <- data.frame(Threshold = NA,
                        prop = prop.within.cid,
                        label = "",
                        category = "Not important effect",
                        sign = "")
    #
    if (avail.cid.below.null & avail.cid.above.null) {
      within.cid <-
        formatN(c(max.cid.below.null, min.cid.above.null),
                digits = digits.cid, big.mark = big.mark)
      #
      within.cid <- paste(">", within.cid[1], "to", "<", within.cid[2])
    }
    else if (avail.cid.below.null) {
      within.cid <-
        formatN(max.cid.below.null, digits = digits.cid, big.mark = big.mark)
      #
      within.cid <- paste(">", within.cid)
    }
    else if (avail.cid.above.null) {
      within.cid <-
        formatN(min.cid.above.null, digits = digits.cid, big.mark = big.mark)
      #
      within.cid <- paste("<", within.cid)
    }
  }
  #
  dat.cid <- rbind(dat.l, dat.w, dat.u)
  #
  Threshold <- prop <- label <- category <- sign <- NULL
  #
  dat.cid %<>%
    mutate(Threshold = formatN(Threshold, digits = digits.cid,
                               big.mark = big.mark, text.NA = ""),
           Threshold = if_else(category != "Not important effect",
                               paste0(sign, Threshold), within.cid),
           prop = paste0(formatPT(100 * prop, digits = digits.percent), "%"),
           category =
             if_else(label == "", category, paste(category, label)),
           category = paste0(category, " ")) %>%
    column_to_rownames("category") %>%
    rename(Percent = prop) %>%
    select(-label, -sign)
  #
  print(dat.cid)
  #
  if (details.methods & !is.null(x$x)) {
    catmeth(x$x,
            FALSE, FALSE, TRUE, FALSE, TRUE,
            #
            func.transf = x$x$func.transf,
            backtransf = FALSE, func.backtransf = x$x$func.backtransf,
            #
            big.mark = "", digits = 2,
            digits.tau = gs("digits.tau"),
            text.tau = gs("text.tau"), text.tau2 = gs("text.tau2"),
            #
            print.tau2 = TRUE,
            print.tau2.ci = FALSE,
            print.tau = TRUE,
            print.tau.ci = FALSE,
            #
            print.I2 = FALSE, text.I2 = "",
            #
            print.df = TRUE, prediction.subgroup = FALSE)
  }
  #
  invisible(NULL)
}
