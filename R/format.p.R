format.p <- function(p, lab=FALSE){
  if (lab)
    format(ifelse(is.na(p), "p = --",
                  ifelse(p < 0.0001, "p < 0.0001",
                         paste("p =", formatC(round(p, 4))))))
  else
    format(ifelse(is.na(p), "      --",
                  ifelse(p < 0.0001, "< 0.0001",
                         paste(" ", formatC(round(p, 4))))))
}
