
require('R.matlab')

mp<-readMat("~/Downloads/convdat.mat")


dim(mp$Ntr)
dim(mp$Ytr)

dim(mp$Nts)
dim(mp$Yts)




data<- as.matrix(mp$Ytr)
keep<- colSums(data)>0
data<- data[,keep]

# data2<- as.matrix(mp$Ytr)[30001:59792,]
data2<- as.matrix(mp$Yts)
data2<- data2[,keep]


colSums(data)

cor_Y<- cor(data)


highCorIdx<- c(1:101)[rowSums(cor_Y==max(cor_Y[cor_Y<0.9],na.rm = T),na.rm = T)==1]

# cor_Y[highCorIdx[2],]

pick<- highCorIdx[2]
# pick<- 1

# pick<- 10

N<- nrow(data)

# N<- 10000

o<- order(cor_Y[pick,],decreasing = T)
feature<-   o<= 100
feature[pick]<-FALSE



X<- data[,feature]
X2<- data2[,feature]

y<- data[,pick]
y2<- data[,pick]

keep <- apply(X, 2, sd)>0
X<- X[,keep]

X2<- X2[,keep]



X<- log(X+1)
X<- cbind(1,X)

X2<- log(X2+1)
X2<- cbind(1,X2)

X<- cbind(X, log(mp$Ntr+1))
X2<- cbind(X2, log(mp$Nts+1))


p<- ncol(X)


B<- diag(1000,p,p)
b<- rep(0,p)

sum(y>0)/N

