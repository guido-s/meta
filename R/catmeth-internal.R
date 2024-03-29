methtxt <- function(x, i, random, method) {
  ##
  if (x$method[i] == "MH")
    txt <- text_MH(x, i, random)
  else if (x$method[i] == "Peto")
    txt <- text_Peto(x, i, random)
  else if (x$method[i] == "Inverse")
    txt <- text_Inverse(x, i, random, method)
  else if (x$method[i] == "GLMM")
    txt <- text_GLMM(x, i, random, method)
  else if (x$method[i] == "Cochran")
    txt <- text_Cochran(x, i, random)
  else if (x$method[i] == "SSW")
    txt <- text_SSW(x, i, random)
  else if (x$method[i] != "") 
    txt <- paste("\n-", x$method[i])
  else
    txt <- ""
  ##
  txt
}


text_MH <- function(x, i, random) {
  meth.i <- x[i, , drop = FALSE]
  ##
  txt <- "\n- Mantel-Haenszel method"
  ##
  if (!is.null(meth.i$sparse)) {
    if ((meth.i$sparse | meth.i$method.incr == "all") &
        (!is.null(meth.i$MH.exact) && meth.i$MH.exact))
      txt <-
        paste0(txt,
               if (random)
                 ", without continuity correction)"
               else
                 " (without continuity correction)")
  }
  ##
  if (random)
    txt <- paste0(txt, " (", gs("text.w.common"), " effect model)")
  ##
  txt
}


text_Peto <- function(x, i, random) {
  meth.i <- x[i, , drop = FALSE]
  ##
  txt <- "\n- Peto method"
  ##
  if (meth.i$model == "common" &
      any(x$method[x$model == "random"] != "Peto"))
    txt <- paste0(txt, " (", gs("text.w.common"), " effect model)")
  else if (meth.i$model == "random" &
      any(x$method[x$model == "common"] != "Peto"))
    txt <- paste0(txt, " (", gs("text.w.random"), " effects model)")
  ##
  txt
}


text_Inverse <- function(x, i, random, method) {
  meth.i <- x[i, , drop = FALSE]
  ##
  txt <- "\n- Inverse variance method"
  ##
  if (meth.i$three.level)
    txt <-
      paste0(txt,
             " (three-level model",
             if (meth.i$rho != 0)
               paste0(", rho = ", meth.i$rho, ")")
             else
               ")")
  else {
    if (meth.i$model == "common" &
        any(x$method[x$model == "random"] != "Inverse"))
      txt <- paste0(txt, " (", gs("text.w.common"), " effect model)")
    else if (meth.i$model == "random" &
             any(x$method[x$model == "common"] != "Inverse"))
      txt <- paste0(txt, " (", gs("text.w.random"), " effects model)")
  }
  ##
  txt <-
    paste0(txt,
           if (method == "metacont" &&
               !is.null(meth.i$pooledvar) &&
               !is.na(meth.i$pooledvar) && meth.i$pooledvar)
             " (with pooled variance for individual studies)"
           else
             "")
  ##
  txt
}


text_GLMM <- function(x, i, random, method) {
  meth.i <- x[i, , drop = FALSE]
  ##
  if (method == "metabin") {
    txt <-
      if (meth.i$model.glmm == "UM.FS")
        "\n- Logistic regression model (fixed study effects)"
      else if (meth.i$model.glmm == "UM.RS")
        paste0(
          "\n- Mixed-effects logistic regression model ",
          "(random study effects)")
      else if (meth.i$model.glmm == "CM.EL")
        paste0(
          "\n- Generalised linear mixed model ",
          "(conditional Hypergeometric-Normal)")
      else if (meth.i$model.glmm == "CM.AL")
        "\n- Generalised linear mixed model (conditional Binomial-Normal)"
  }
  else if (method == "metainc") {
    txt <-
      if (meth.i$model.glmm == "UM.FS")
        "\n- Poisson regression model (fixed study effects)"
      else if (meth.i$model.glmm == "UM.RS")
        paste0(
          "\n- Mixed-effects Poisson regression model ",
          "(random study effects)")
      else if (meth.i$model.glmm == "CM.EL")
        paste0(
          "\n- Generalised linear mixed model ",
          "(conditional Poisson-Normal)")
  }
  else if (method == "metaprop")
    txt <-
      "\n- Random intercept logistic regression model"
  else if (method == "metarate")
    txt <-
      "\n- Random intercept Poisson regression model"
  ##
  if (meth.i$model == "common" &
      any(x$method[x$model == "random"] != "GLMM"))
    txt <- paste0(txt, " (", gs("text.w.common"), " effect model)")
  else if (meth.i$model == "random" &
      any(x$method[x$model == "common"] != "GLMM"))
    txt <- paste0(txt, " (", gs("text.w.random"), " effects model)")
  ##
  txt
}


text_Cochran <- function(x, i, random) {
  meth.i <- x[i, , drop = FALSE]
  ##
  txt <- "\n- Cochran method"
  ##
  if (meth.i$model == "common" &
      any(x$method[x$model == "random"] != "Peto"))
    txt <- paste0(txt, " (", gs("text.w.common"), " effect model)")
  ##
  txt
}


text_SSW <- function(x, i, random) {
  meth.i <- x[i, , drop = FALSE]
  ##
  txt <- "\n- Sample size method"
  ##
  if (meth.i$model == "common" &
      any(x$method[x$model == "random"] != "Peto"))
    txt <- paste0(txt, " (", gs("text.w.common"), " effect model)")
  else if (meth.i$model == "random" &
      any(x$method[x$model == "common"] != "Peto"))
    txt <- paste0(txt, " (", gs("text.w.random"), " effects model)")
  ##
  txt
}
