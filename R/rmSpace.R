rmSpace <- function(x, end=FALSE, pat=" "){
  
  if ( !end ){
    while ( any(substring(x, 1, 1) == pat, na.rm=TRUE) ){
      sel <- substring(x, 1, 1) == pat
      x[sel] <- substring(x[sel], 2)
  }
  }
  else{
    last <- nchar(x)
    
    while ( any(substring(x, last, last) == pat, na.rm=TRUE) ){
      sel <- substring(x, last, last) == pat
      x[sel] <- substring(x[sel], 1, last[sel]-1)
      last <- nchar(x)
    }
  }
  
  x
}
