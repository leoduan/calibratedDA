zinb <- read.csv("http://www.karlin.mff.cuni.cz/~pesta/prednasky/NMFM404/Data/fish.csv")
# zinb <- within(zinb, {
  # nofish <- factor(nofish)
  # livebait <- factor(livebait)
  # camper <- factor(camper)
# })


require("ImbalancedPG")

n<- nrow(zinb)
X<- zinb

y<- zinb[,8]
X<- cbind(1,zinb[,1:7])
B<- diag(1000,8,8)
b<- rep(0,8)

as.matrix((X))

zinb[1,]

fit<- ImbalancedPG::poisson_reg( y, as.matrix(X), b,B, r0ini =  2,c = 1,burnin = 1000,run = 1000)

ts.plot(fit$beta[,1])
ts.plot(fit$beta[,2])
ts.plot(fit$beta[,3])
ts.plot(fit$beta[,4])

acf(rowSums(fit$beta))
