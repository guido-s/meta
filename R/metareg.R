#' Meta-regression
#' 
#' @description
#' Meta-regression for objects of class \code{meta}. This is a wrapper
#' function for the R function \code{\link[metafor]{rma.uni}} in the R
#' package \bold{metafor} (Viechtbauer 2010).
#' 
#' @details
#' This R function is a wrapper function for R function
#' \code{\link[metafor]{rma.uni}} in the R package \bold{metafor}
#' (Viechtbauer 2010).
#' 
#' Note, results are not back-transformed in printouts of
#' meta-analyses using summary measures with transformations, e.g.,
#' log risk ratios are printed instead of the risk ratio if argument
#' \code{sm = "RR"} and logit transformed proportions are printed if
#' argument \code{sm = "PLOGIT"}.
#' 
#' Argument '\dots{}' can be used to pass additional arguments to R
#' function \code{\link[metafor]{rma.uni}}. For example, argument
#' \code{control} to provide a list of control values for the
#' iterative estimation algorithm. See help page of R function
#' \code{\link[metafor]{rma.uni}} for more details.
#' 
#' @param x An object of class \code{meta}.
#' @param formula Either a character string or a formula object.
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance tau-squared. Either
#'   \code{"FE"}, \code{"DL"}, \code{"REML"}, \code{"ML"},
#'   \code{"HS"}, \code{"SJ"}, \code{"HE"}, or \code{"EB"}, can be
#'   abbreviated.
#' @param hakn A logical indicating whether the method by Hartung and
#'   Knapp should be used to adjust test statistics and confidence
#'   intervals.
#' @param level.comb The level used to calculate confidence intervals
#'   for parameter estimates in the meta-regression model.
#' @param intercept A logical indicating whether an intercept should
#'   be included in the meta-regression model.
#' @param \dots Additional arguments passed to R function
#'   \code{\link[metafor]{rma.uni}}.
#' 
#' @return
#' An object of class \code{c("metareg", "rma.uni","rma")}. Please
#' look at the help page of R function \code{\link[metafor]{rma.uni}}
#' for more details on the output from this function.
#' 
#' In addition, a list \code{.meta} is added to the output containing
#' the following components:
#' \item{x, formula, method.tau, hakn, level.comb, intercept}{As
#'   defined above.}
#' \item{dots}{Information provided in argument '\dots{}'.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package \bold{meta} used to create
#'   object.}
#' \item{version.metafor}{Version of R package \bold{metafor} used to
#'   create object.}
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{bubble}}, \code{\link{summary.meta}},
#'   \code{\link{metagen}}
#' 
#' @references
#' Viechtbauer W (2010):
#' Conducting Meta-Analyses in R with the Metafor Package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
#' 
#' @keywords models regression
#' 
#' @examples
#' data(Fleiss1993cont)
#' # Add some (fictitious) grouping variables:
#' Fleiss1993cont$age <- c(55, 65, 55, 65, 55)
#' Fleiss1993cont$region <- c("Europe", "Europe", "Asia", "Asia", "Europe")
#' 
#' m1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'                data = Fleiss1993cont, sm = "MD")
#' \dontrun{
#' # Warnings due to wrong ordering of arguments (order has changed
#' # with version 3.0-0 of R package meta)
#' #
#' metareg(~ region, m1)
#' metareg(~ region, data = m1)
#' 
#' # Warning as no information on covariate is available
#' #
#' metareg(m1)
#' }
#' 
#' # Do meta-regression for covariate region
#' #
#' mu2 <- update(m1, byvar = region, tau.common = TRUE, comb.fixed = FALSE)
#' metareg(mu2)
#' 
#' # Same result for
#' # - tau-squared
#' # - test of heterogeneity
#' # - test for subgroup differences
#' # (as argument 'tau.common' was used to create mu2)
#' #
#' mu2
#' metareg(mu2, intercept = FALSE)
#' metareg(m1, region)
#' 
#' # Different result for
#' # - tau-squared
#' # - test of heterogeneity
#' # - test for subgroup differences
#' # (as argument 'tau.common' is - by default - FALSE)
#' #
#' mu1 <- update(m1, byvar = region)
#' mu1
#' 
#' # Generate bubble plot
#' #
#' bubble(metareg(mu2))
#' 
#' # Do meta-regression with two covariates
#' #
#' metareg(mu1, region + age)
#' 
#' # Do same meta-regressions using formula notation
#' #
#' metareg(m1, ~ region)
#' metareg(mu1, ~ region + age)
#' 
#' # Do meta-regression using REML method and print intermediate
#' # results for iterative estimation algorithm; furthermore print
#' # results with three digits.
#' #
#' metareg(mu1, region, method.tau = "REML",
#'         control = list(verbose = TRUE), digits = 3)
#' 
#' # Use Hartung-Knapp method
#' #
#' mu3 <- update(mu2, hakn = TRUE)
#' mu3
#' metareg(mu3, intercept = FALSE)
#' 
#' @export metareg


metareg <- function(x, formula, method.tau = x$method.tau,
                    hakn = x$hakn, level.comb = x$level.comb,
                    intercept = TRUE, ...) {

  if ("data" %in% names(list(...))) {
    warning("Please note, argument 'data' has been renamed to 'x' ",
            "in version 3.0-0 of R package meta ",
            "(see help page of R function metareg). ",
            "No meta-regression conducted.")
    return(invisible(NULL))
  }

  if (is.call(x)) {
    warning("Please note, first two arguments of R function metareg ",
            "have been interchanged in version 3.0-0 of R package meta. ",
            "No meta-regression conducted.")
    return(invisible(NULL))
  }


  ##
  ##
  ## (1) Check for meta object
  ##
  ##
  chkclass(x, "meta")
  
  
  ##
  ## Assignments
  ##
  TE <- x$TE
  seTE <- x$seTE
  method <- x$method
  ##
  model.glmm <- x$model.glmm
  ##
  metabin <- inherits(x, "metabin")
  metainc <- inherits(x, "metainc")
  metaprop <- inherits(x, "metaprop")
  metarate <- inherits(x, "metarate")
  ##
  if (metabin) {
    event.e <- x$event.e
    n.e <- x$n.e
    event.c <- x$event.c
    n.c <- x$n.c
  }
  else if (metainc) {
    event.e <- x$event.e
    time.e <- x$time.e
    event.c <- x$event.c
    time.c <- x$time.c
  }
  else if (metaprop) {
    event <- x$event
    n <- x$n
  }
  else if (metarate) {
    event <- x$event
    time <- x$time
  }


  if (missing(formula)) {
    if (isCol(x$data, ".byvar"))
      if (intercept)
        formula <- as.call(~ .byvar)
      else
        formula <- as.call(~ .byvar - 1)
    else {
      warning("No meta-regression conducted as argument 'formula' ",
              "is missing and no information is provided on subgroup ",
              "variable, i.e. list element 'byvar' in meta-analysis object ",
              "'x' (see help page of R function metareg).")
      return(invisible(NULL))
    }
  }
  else {
    formula.text <- deparse(substitute(formula))
    formula.text <- gsub("~", "", formula.text)
    formula.text <- gsub("\\\"", "", formula.text)
    formula.text <- gsub("\\\'", "", formula.text)
    if (!intercept)
      formula.text <- paste0(formula.text, " - 1")
    formula <- as.formula(paste("~", formula.text))
  }


  if (is.null(method.tau))
    method.tau <- "DL"
  ##
  method.tau <- setchar(method.tau, c(.settings$meth4tau, "FE"))
  ##
  if (method.tau == "PM") {
    warning("Meta-regresion method not available for method.tau = \"PM\". ",
            "Using REML method instead (method.tau = \"REML\").")
    method.tau <- "REML"
  }
  ##
  chklogical(hakn)
  ##
  chklevel(level.comb)
  chklogical(intercept)

  if (is.null(x$data)) {
    warning("Necessary data not available. Please, recreate meta-analysis ",
            "object without option 'keepdata = FALSE'.")
    return(invisible(NULL))
  }


  ##
  ## Use subset of studies in meta-regression
  ##
  if (!is.null(x$subset))
    dataset <- x$data[x$subset, ]
  else
    dataset <- x$data


  ##
  ## Exclude studies from meta-regression
  ##
  if (!is.null(x$exclude)) {
    exclude <- dataset$.exclude
    dataset <- dataset[!dataset$.exclude, ]
  }
  else
    exclude <- rep(FALSE, nrow(dataset))


  ##
  ## Argument test in rma.uni() and rma.glmm()
  ##
  test <- ifelse(!hakn, "z",
                 ifelse(method != "GLMM", "knha", "t"))

  ##
  ## Covariate 'x' make problems without removing meta-analysis object x
  ##
  ..x <- x
  rm(x)
  ##
  warn.FE <- paste("Fallback to fixed effect model (argument",
                   "method.tau = \"FE\") due to small number of studies.")
  ##
  if (method != "GLMM")
    res <- rma.uni(yi = TE[!exclude], sei = seTE[!exclude],
                   data = dataset,
                   mods = formula, method = method.tau,
                   test = test, level = 100 * level.comb,
                   ...)
  else {
    if (metabin) {
      if (sum(!exclude) > 2)
        res <- rma.glmm(ai = event.e[!exclude], n1i = n.e[!exclude],
                        ci = event.c[!exclude], n2i = n.c[!exclude],
                        data = dataset,
                        mods = formula, method = method.tau,
                        test = test, level = 100 * level.comb,
                        measure = "OR", model = model.glmm,
                        ...)
      else {
        if (method.tau != "FE")
          warning(warn.FE)
        res <- rma.glmm(ai = event.e[!exclude], n1i = n.e[!exclude],
                        ci = event.c[!exclude], n2i = n.c[!exclude],
                        data = dataset,
                        mods = formula, method = "FE",
                        test = test, level = 100 * level.comb,
                        measure = "OR", model = model.glmm,
                        ...)
      }
    }
    else if (metainc) {
      if (sum(!exclude) > 2)
        res <- rma.glmm(x1i = event.e[!exclude], t1i = time.e[!exclude],
                        x2i = event.c[!exclude], t2i = time.c[!exclude],
                        data = dataset,
                        mods = formula, method = method.tau,
                        test = test, level = 100 * level.comb,
                        measure = "IRR", model = model.glmm,
                        ...)
      else {
        if (method.tau != "FE")
          warning(warn.FE)
        res <- rma.glmm(x1i = event.e[!exclude], t1i = time.e[!exclude],
                        x2i = event.c[!exclude], t2i = time.c[!exclude],
                        data = dataset,
                        mods = formula, method = "FE",
                        test = test, level = 100 * level.comb,
                        measure = "IRR", model = model.glmm,
                        ...)
      }
    }
    else if (metaprop) {
      if (sum(!exclude) > 2)
        res <- rma.glmm(xi = event[!exclude], ni = n[!exclude],
                        data = dataset,
                        mods = formula, method = method.tau,
                        test = test, level = 100 * level.comb,
                        measure = "PLO",
                        ...)
      else {
        if (method.tau != "FE")
          warning(warn.FE)
        res <- rma.glmm(xi = event[!exclude], ni = n[!exclude],
                        data = dataset,
                        mods = formula, method = "FE",
                        test = test, level = 100 * level.comb,
                        measure = "PLO",
                        ...)
      }
    }
    else if (metarate) {
      if (sum(!exclude) > 2)
        res <- rma.glmm(xi = event[!exclude], ti = time[!exclude],
                        data = dataset,
                        mods = formula, method = method.tau,
                        test = test, level = 100 * level.comb,
                        measure = "IRLN",
                        ...)
      else {
        if (method.tau != "FE")
          warning(warn.FE)
        res <- rma.glmm(xi = event[!exclude], ti = time[!exclude],
                        data = dataset,
                        mods = formula, method = "FE",
                        test = test, level = 100 * level.comb,
                        measure = "IRLN",
                        ...)
      }
    }
  }


  res$.meta <- list(x = ..x,
                    formula = formula,
                    method.tau = method.tau,
                    hakn = hakn,
                    level.comb = level.comb,
                    intercept = intercept,
                    dots = list(...),
                    call = match.call(),
                    version = packageDescription("meta")$Version,
                    version.metafor = packageDescription("metafor")$Version)

  class(res) <- c("metareg", class(res))

  res
}
