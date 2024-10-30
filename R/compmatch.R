compmatch <- function(x, split) {
  
  if (split %in% gs("special.characters"))
    split <- paste0("\\", split)
  
  res <- any(grepl(split, x))
  
  res
}
