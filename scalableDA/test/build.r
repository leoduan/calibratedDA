library("devtools")

setwd("~/git/calibratedDA/")
build('scalableDA')
install.packages("scalableDA_1.0.tar.gz", repos = NULL, type = "source")

require("RcppArmadillo")

require("scalableDA")
