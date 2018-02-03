library(alphashape3d)
library(mvtnorm)

Sigma1 <- diag(c(1,1,1))# Covariance matrix
Sigma2 <- diag(c(1,1,1))# Covariance matrix
mu1 <- c(0,0,0)
mu2<-c(0,2,2)

alpha<-100

n1<-5000
n2<-5000
x1<-rmvnorm(n1,mean=mu1,sigma=Sigma1) # Sample1
x2<-rmvnorm(n2,mean=mu2,sigma=Sigma2) # Sample2


a1<-ashape3d(x1,alpha=alpha) # ashape sample 1
a2<-ashape3d(x2,alpha=alpha) # ashape sample 2
plot(a1,trans=0.2)
plot(a2,clear=FALSE,col=c(4,4,4),trans=0.2)

x<-rbind(x1,x2) 
in1<-inashape3d(a1, indexAlpha = 1, x)
in2<-inashape3d(a2, indexAlpha = 1, x)

xin<-x[in1&in2,] # Points in intersection
a3<-ashape3d(xin,alpha=alpha)
plot(a3,clear=FALSE,col=c(3,3,3),trans=1) # Intersection in green
volume_ashape3d(a3) # Interseciton volume
