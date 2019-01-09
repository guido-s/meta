##
## (1) Make R packages available
##
library(devtools)
library(roxygen2)


##
## (2) Create documentation file(s)
##
document("../meta")


##
## (3) Build R package and PDF file with help pages
##
build("../meta")
build_manual("../meta")


##
## (4) Install R package
##
install("../meta")


##
## (5) Check R package
##
check("../meta")


##
## (6) Check R package (with dontrun examples)
##
check("../meta", run_dont_test = TRUE)
