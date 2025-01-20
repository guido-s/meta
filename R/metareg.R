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
#' @param level.ma The level used to calculate confidence intervals
#'   for parameter estimates in the meta-regression model.
#' @param intercept A logical indicating whether an intercept should
#'   be included in the meta-regression model.
#' @param \dots Additional arguments passed to R function
#'   \code{\link[metafor]{rma.uni}}.
#' 
#' @return
#' An object of class \code{c("metareg", "rma.uni", "rma")}. Please
#' look at the help page of R function \code{\link[metafor]{rma.uni}}
#' for more details on the output from this function.
#' 
#' In addition, a list \code{.meta} is added to the output containing
#' the following components:
#' \item{x, formula, method.tau, hakn, level.ma, intercept}{As
#'   defined above.}
#' \item{dots}{Information provided in argument '\dots{}'.}
#' \item{call}{Function call.}
#' \item{version}{Version of R package \bold{meta} used to create
#'   object.}
#' \item{version.metafor}{Version of R package \bold{metafor} used to
#'   create object.}
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
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
#'   data = Fleiss1993cont, sm = "SMD")
#' \dontrun{
#' # Error due to wrong ordering of arguments (order has changed in
#' # R package meta, version 3.0-0)
#' #
#' try(metareg(~ region, m1))
#' try(metareg(~ region, data = m1))
#' 
#' # Warning as no information on covariate is available
#' #
#' metareg(m1)
#' }
#' 
#' # Do meta-regression for covariate region
#' #
#' mu2 <- update(m1, subgroup = region, tau.common = TRUE, common = FALSE)
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
#' mu1 <- update(m1, subgroup = region)
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
#'   control = list(verbose = TRUE), digits = 3)
#' 
#' # Use Hartung-Knapp method
#' #
#' mu3 <- update(mu2, method.random.ci = "HK")
#' mu3
#' metareg(mu3, intercept = FALSE)
#' 
#' @rdname metareg
#' @method metareg meta
#' @export


metareg.meta <- function(x, formula, method.tau = x$method.tau,
                         hakn = x$method.random.ci == "HK",
                         level.ma = x$level.ma,
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
  chksuitable(x, "Meta-regression", c("metamerge", "netpairwise"),
               check.mlm = FALSE)
  ##
  x <- updateversion(x)
  
  
  ##
  ## Assignments
  ##
  TE <- x$TE
  seTE <- x$seTE
  method <- x$method
  ##
  model.glmm <- x$model.glmm
  ##
  three.level <- !is.null(x$k.study) && x$k != x$k.study
  ##
  metabin <- inherits(x, "metabin")
  metainc <- inherits(x, "metainc")
  metaprop <- inherits(x, "metaprop")
  metarate <- inherits(x, "metarate")
  ##
  if (metabin) {
    if (x$method == "LRP")
      stop("Meta-regression not implemented for penalised logistic regression.",
           call. = FALSE)
    #
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
    if (isCol(x$data, ".subgroup"))
      if (intercept)
        formula <- as.call(~ .subgroup)
      else
        formula <- as.call(~ .subgroup - 1)
    else {
      warning("No meta-regression conducted as argument 'formula' ",
              "is missing and no information is provided on subgroup ",
              "variable, i.e. list element 'subgroup' in meta-analysis object ",
              "'x' (see help page of R function metareg).")
      return(invisible(NULL))
    }
  }
  else {
    as.char <-
      try(is.character(formula) || is.numeric(formula) ||
          is.function(formula), silent = TRUE)
    ##
    if (inherits(as.char, "try-error"))
      formula.text <- deparse(substitute(formula))
    else {
      if ((as.char & length(formula) != 1) ||
          is.function(formula))
        formula.text <- paste("~", deparse(substitute(formula)))
      else
        formula.text <- deparse(formula)
    }
    ##
    formula.text <- gsub("~", "", formula.text)
    formula.text <- gsub("\\\"", "", formula.text)
    formula.text <- gsub("\\\'", "", formula.text)
    ##
    nulldata <- is.null(x$data)
    ##
    for (i in as.vector(sapply(strsplit(formula.text, "+", fixed = TRUE),
                               trimws, which = "both"))) {
      if (".GlobalEnv" %in% find(i)) {
        if (nulldata) {
          if (is.null(x$data)) {
            x$data <- data.frame(first = .GlobalEnv[[i]])
            names(x$data) <- i
          }
          else
            x$data[[i]] <- .GlobalEnv[[i]]
        }
        else {
          if (isCol(x$data, i))
            warning("R object '", i, "' found in .GlobalEnv used in ",
                    "meta-regression instead of data from ",
                    "meta-analysis object.")
          x$data[[i]] <- .GlobalEnv[[i]]
        }
      }
    }
    ##
    if (!intercept)
      formula.text <- paste0(formula.text, " - 1")
    ##
    formula <- as.formula(paste("~", formula.text))
  }
  
  
  if (is.null(method.tau))
    method.tau <- "DL"
  ##
  method.tau <- setchar(method.tau, c(gs("meth4tau"), "FE"))
  ##
  chklogical(hakn)
  ##
  chklevel(level.ma)
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
  ## Argument test in rma.uni(), rma.glmm() and rma.mv()
  ##
  test <- ifelse(!hakn, "z",
                 ifelse(method == "GLMM" | three.level, "t", "knha"))

  ##
  ## Covariate 'x' makes problems without removing meta-analysis object x
  ##
  ..x <- x
  rm(x)
  ##
  warn.FE <- paste("Fallback to common effect model (argument",
                   "method.tau = \"FE\") due to small number of studies.")
  ##
  if (method != "GLMM") {
    ##
    ## Three-level model
    ##
    if (three.level) {
      ##
      dataset$.idx <- seq_len(nrow(dataset))
      ##
      res <-
        runNN(rma.mv,
              list(yi = TE[!exclude], V = seTE[!exclude]^2,
                   data = dataset,
                   mods = formula, method = method.tau,
                   random = as.call(~ 1 | .id / .idx),
                   test = test, level = 100 * level.ma,
                   ...))
    }
    else
      res <-
        runNN(rma.uni,
              list(yi = TE[!exclude], sei = seTE[!exclude],
                   data = dataset,
                   mods = formula, method = method.tau,
                   test = test, level = 100 * level.ma,
                   ...))
  }
  else {
    if (metabin) {
      if (sum(!exclude) > 2)
        res <-
          runNN(rma.glmm,
                list(ai = event.e[!exclude], n1i = n.e[!exclude],
                     ci = event.c[!exclude], n2i = n.c[!exclude],
                     data = dataset,
                     mods = formula, method = method.tau,
                     test = test, level = 100 * level.ma,
                     measure = "OR", model = model.glmm,
                     ...))
      else {
        if (method.tau != "FE")
          warning(warn.FE)
        res <-
          runNN(rma.glmm,
                list(ai = event.e[!exclude], n1i = n.e[!exclude],
                     ci = event.c[!exclude], n2i = n.c[!exclude],
                     data = dataset,
                     mods = formula, method = "FE",
                     test = test, level = 100 * level.ma,
                     measure = "OR", model = model.glmm,
                     ...))
      }
    }
    else if (metainc) {
      if (sum(!exclude) > 2)
        res <-
          runNN(rma.glmm,
                list(x1i = event.e[!exclude], t1i = time.e[!exclude],
                     x2i = event.c[!exclude], t2i = time.c[!exclude],
                     data = dataset,
                     mods = formula, method = method.tau,
                     test = test, level = 100 * level.ma,
                     measure = "IRR", model = model.glmm,
                     ...))
      else {
        if (method.tau != "FE")
          warning(warn.FE)
        res <-
          runNN(rma.glmm,
                list(x1i = event.e[!exclude], t1i = time.e[!exclude],
                     x2i = event.c[!exclude], t2i = time.c[!exclude],
                     data = dataset,
                     mods = formula, method = "FE",
                     test = test, level = 100 * level.ma,
                     measure = "IRR", model = model.glmm,
                     ...))
      }
    }
    else if (metaprop) {
      if (sum(!exclude) > 2)
        res <-
          runNN(rma.glmm,
                list(xi = event[!exclude], ni = n[!exclude],
                      data = dataset,
                     mods = formula, method = method.tau,
                     test = test, level = 100 * level.ma,
                     measure = "PLO",
                     ...))
      else {
        if (method.tau != "FE")
          warning(warn.FE)
        res <-
          runNN(rma.glmm,
                 list(xi = event[!exclude], ni = n[!exclude],
                      data = dataset,
                      mods = formula, method = "FE",
                      test = test, level = 100 * level.ma,
                      measure = "PLO",
                      ...))
      }
    }
    else if (metarate) {
      if (sum(!exclude) > 2)
        res <-
          runNN(rma.glmm,
                list(xi = event[!exclude], ti = time[!exclude],
                     data = dataset,
                     mods = formula, method = method.tau,
                     test = test, level = 100 * level.ma,
                     measure = "IRLN",
                     ...))
      else {
        if (method.tau != "FE")
          warning(warn.FE)
        res <-
          runNN(rma.glmm,
                list(xi = event[!exclude], ti = time[!exclude],
                     data = dataset,
                     mods = formula, method = "FE",
                     test = test, level = 100 * level.ma,
                     measure = "IRLN",
                     ...))
      }
    }
  }
  
  
  res$.meta <- list(x = ..x,
                    formula = formula,
                    method.tau = method.tau,
                    hakn = hakn,
                    level.ma = level.ma,
                    intercept = intercept,
                    dots = list(...),
                    call = match.call(),
                    version = packageDescription("meta")$Version,
                    version.metafor = packageDescription("metafor")$Version)
  
  class(res) <- c("metareg", class(res))
  res$call <- NULL
  
  res
}





#' @rdname metareg
#' @export metareg


metareg <- function(x, ...) 
  UseMethod("metareg")





#' @rdname metareg
#' @method metareg default
#' @export


metareg.default <- function(x, ...)
  stop("Meta-regression not available for an object of class '",
       class(x)[1], "'.",
       call. = FALSE)
