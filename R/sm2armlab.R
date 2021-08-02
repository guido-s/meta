sm2armlab <- function(sm) {
  
  sms <- c("SMD", "MD", "RD", "RR", "OR", "DOR", "ASD",
           "IRD", "IRR", "ROM", "HR")
  
  labs <- c("Std. Mean", "Mean", "Risk", "Risk", "Odds", "Odds", "Arcsine",
            "Rate", "Rate", "Mean", "Hazard")
  
  labs <- labs[pmatch(sm, sms)]
  labs[is.na(labs)] <- ""
  
  labs
}
