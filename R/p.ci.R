p.ci <- function(lower, upper, rmspace=TRUE,
                 bracket.left="[",
                 separator="; ",
                 bracket.right="]",
                 justify.lower="right",
                 justify.upper=justify.lower,
                 ...
                 ){

  ## Change layout of CIs
  ##
  ibracktype <- charmatch(.settings$CIbracket,
                          c("[", "(", "{", ""), nomatch = NA)
  if (is.na(ibracktype) | ibracktype==0){
    warning("No valid bracket type specified globally for R package meta: ",
            .settings$CIbracket,
            "\n  Using default bracket type: '['. See help page on R command 'cilayout' for further information.")
    bracktype <- "["
  }
  else
    bracktype <- c("[", "(", "{", "")[ibracktype]
  ##
  if (bracktype=="["){
    bracketLeft <- "["
    bracketRight <- "]"
  }
  else if (bracktype=="("){
    bracketLeft <- "("
    bracketRight <- ")"
  }
  else if (bracktype=="{"){
    bracketLeft <- "{"
    bracketRight <- "}"
  }
  else if (bracktype==""){
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
  if (missing(separator))
    separator <- .settings$CIseparator
  
  if (rmspace){
    lower <- rmSpace(lower)
    upper <- rmSpace(upper)
  }
  ##
  ifelse (lower=="NA" & upper=="NA",
          "",
          paste(bracket.left,
                format(lower, justify=justify.lower),
                separator,
                format(upper, justify=justify.upper),
                bracket.right, sep=""))
}
