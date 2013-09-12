cilayout <- function(bracket="[", separator="; "){
  
  ibracket <- charmatch(bracket,
                        c("[", "(", "{", ""),
                        nomatch = NA)
  ##
  if (is.na(ibracket) | ibracket==0)
    stop("No valid bracket type specified. Admissible values: '[', '(', '{', '\"\"'")
  
  bracket <- c("[", "(", "{", "")[ibracket]
  
  setOption("CIbracket", bracket)
  setOption("CIseparator", separator)
  
  invisible(NULL)
}
