text_meth <- function(x, i, random, method) {
  ##
  if (x$method[i] == "MH" & x$model[i] == "common")
    txt <- text_MH(x, i, random)
  else if (x$method[i] == "Peto")
    txt <- text_Peto(x, i, random)
  else if (x$method[i] == "Inverse" |
           (x$method[i] == "MH" & x$model[i] == "random"))
    txt <- text_Inverse(x, i, random, method)
  else if (x$method[i] == "GLMM")
    txt <- text_GLMM(x, i, random, method)
  else if (x$method[i] == "LRP")
    txt <- text_LRP(x, i, random, method)
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
  txt <- paste0("\n- ",
                if (!is.null(meth.i$MH.exact) && meth.i$MH.exact)
                  "Exact ",
                "Mantel-Haenszel method")
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
      !any(x$method[x$model == "random"] == "Peto"))
    txt <- paste0(txt, " (", gs("text.w.common"), " effect model)")
  else if (meth.i$model == "random" &
      !any(x$method[x$model == "common"] == "Peto"))
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
    if (meth.i$model == "common" && any(x$model == "random") &&
        !any(x$method[x$model == "random"] == "Inverse"))
      txt <- paste0(txt, " (", gs("text.w.common"), " effect model)")
    else if (meth.i$model == "random" && any(x$model == "common") &&
             !any(x$method[x$model == "common"] == "Inverse"))
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
  if ("metabin" %in% method) {
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
  else if ("metainc" %in% method) {
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
  else if ("metaprop" %in% method)
    txt <-
      "\n- Random intercept logistic regression model"
  else if ("metarate" %in% method)
    txt <-
      "\n- Random intercept Poisson regression model"
  ##
  if (meth.i$model == "common" && any(x$model == "random") &&
      !any(x$method[x$model == "random"] == "GLMM"))
    txt <- paste0(txt, " (", gs("text.w.common"), " effect model)")
  else  if (meth.i$model == "random" && any(x$model == "common") &&
            !any(x$method[x$model == "common"] == "GLMM"))
    txt <- paste0(txt, " (", gs("text.w.random"), " effects model)")
  ##
  txt
}


text_LRP <- function(x, i, random, method) {
  meth.i <- x[i, , drop = FALSE]
  #
  txt <- "\n- Penalised logistic regression model"
  #
  if (any(x$method == "LRP" & x$model == "random"))
    txt <- paste0(txt, " (sqrt(phi) = ", round(sqrt(x$phi[i]), 4), ")")
  #
  txt
}


text_Cochran <- function(x, i, random) {
  meth.i <- x[i, , drop = FALSE]
  #
  txt <- "\n- Cochran method"
  ##
  if (random)
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
      !any(x$method[x$model == "random"] == "SSW"))
    txt <- paste0(txt, " (", gs("text.w.common"), " effect model)")
  else if (meth.i$model == "random" &
      !any(x$method[x$model == "common"] != "SSW"))
    txt <- paste0(txt, " (", gs("text.w.random"), " effects model)")
  ##
  txt
}
