library(devtools)

setwd("~/git/ImbalancedPG/")
build('scalableDA')
install.packages("scalableDA_1.0.tar.gz", repos = NULL, type = "source")

require("RcppArmadillo")

require("scalableDA")
