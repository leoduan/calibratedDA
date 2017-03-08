load("./jags_zip_fit.Rda")

ts.plot(jags_zip_fit[[1]][,1])


fit3<- jags_zip_fit[[1]]
dim(fit3)



acfADA<- apply(fit3[,1:50], 2, getAcf)

df<- melt(acfADA)

df$X1<- df$X1-1
df$X2<- as.factor(df$X2)

names(df)<- c("Lag","ParIndex","ACF")

pdf("./poisson_zip_da_acf.pdf",3,3)
ggplot(data=df, aes(x=Lag, y=ACF,col= as.factor(ParIndex)))+ geom_line(size=.75)+   theme(legend.position="none")+  scale_colour_manual(values = rep("black",100))
dev.off()

pdf("./traceplot_poisson_zip_da.pdf",6,6)
par(mfrow=c(3,1))
ts.plot((fit3[,1]),xlab="Iteration",ylab=TeX('$\\beta_0$'))
ts.plot((fit3[,9]),xlab="Iteration",ylab=TeX('$\\beta_8$'))
ts.plot((fit3[,45]),xlab="Iteration",ylab=TeX('$\\beta_45$'))
dev.off()
