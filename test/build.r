library(devtools)

setwd("~/git/")
build('ImbalancedPG/')
install.packages("ImbalancedPG_1.0.tar.gz", repos = NULL, type = "source")

require("RcppArmadillo")

require("ImbalancedPG")
