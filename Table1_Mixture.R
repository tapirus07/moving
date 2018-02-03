# ====================================================
# Example of estimation of the volume of HDR 
# ====================================================
set.seed(1234)
source("hdr3d.R") # Function that estimates the HDR
library(alphashape3d)

# Example 
# Standard Multivariate Normal
# =======================================================================
library(mvtnorm) # Generate multivariate normals

# Dataset (3)
#============
Sigma1 <- diag(c(1,1,1))# Covariance matrix
Sigma2 <- diag(c(1,1,1))# Covariance matrix
mu1 <- c(0,0,0)
mu2<-c(0,15,15)
p1 <- 0.5 # Probability of the first group

# Dataset (4)
#============
Sigma1 <- diag(c(1,1,1))# Covariance matrix
Sigma2 <- diag(c(1,3,1))# Covariance matrix
mu1 <- c(0,0,0)
mu2<-c(0,11,8)
p1 <- 0.5 # Probability of the first group


# Input parameters
# ===========================
lev<-c(0.5,0.95) # Level of the desired highest density region (it can be a vector)
grid<-20         # Number of grid points
h<-NULL          # h=NULL corresponds to the Normal optimal smoothing (see Input in hdr3d)

alpha<-3         # Value of alpha for the alpha-shape reconstruction
scale<-FALSE     # Set scale<-TRUE if the scales are very different in the three dimensions (for the alpha-shape computation)
plotHDR<-TRUE    # Do you want to plot the results?
aspect<-FALSE    # Do you want to adjusts the aspect ratio in the plot?

nsize<-c(30,50,500,1000)  # Sample sizes
B<-0                   # Size of the sample generated from the kernel density to compute the quantile for the HDR estimation
                       # B=0 -> no sample is generated (it works with the sample points)
nsamples<-100             # Number of samples


# Approximated volume MonteCarlo integration 
# ----------------------------------------------------------------------
na <- 500000
n1 <- rbinom(1,size=na,prob=p1)  ## how many from first distribution?
n2 <- na-n1
x1<-rmvnorm(n1,mean=mu1,sigma=Sigma1)
x2<-rmvnorm(n2,mean=mu2,sigma=Sigma2)
x<-rbind(x1,x2)
fnx<-p1*dmvnorm(x,mean=mu1,sigma=Sigma1)+(1-p1)*dmvnorm(x,mean=mu2,sigma=Sigma2)
qfn<-quantile(fnx,1-lev)
rg<-apply(x,2,range)
xmc<-cbind(runif(na,rg[1,1],rg[2,1]),runif(na,rg[1,2],rg[2,2]),runif(na,rg[1,3],rg[2,3]))
fnxmc<-p1*dmvnorm(xmc,mean=mu1,sigma=Sigma1)+(1-p1)*dmvnorm(xmc,mean=mu2,sigma=Sigma2)
mc_HDRvol<-numeric()
for (k in 1:length(lev)){
mc_HDRvol[k]<-mean(fnxmc>=qfn[k])*prod(diff(rg))
}

# Volume estimation
# ------------------
volumeHDR<-array(dim=c(length(nsize),nsamples,length(lev)),dimnames=list(NULL,NULL,c("level 0.95","level 0.5"))) 
for(ns in 1:length(nsize)){
n<-nsize[ns]

for (j in 1:nsamples){
n1 <- rbinom(1,size=n,prob=p1)  ## how many from first distribution?
n2 <- n-n1
x1<-rmvnorm(n1,mean=mu1,sigma=Sigma1)
x2<-rmvnorm(n2,mean=mu2,sigma=Sigma2)
x<-rbind(x1,x2)



# Estimation of the highest density region (HDR) in 3D
# ====================================================================================
hdrout<-hdr3d(x,h=h,grid=grid,lev=lev,B=B)

for (i in 1:length(lev)){
hdr<-hdrout$HDR[[i]]

# alpha-shape and volume computation
# ====================================================================================
if(scale){
hdrorig<-hdr
real_range<-diff(apply(hdr,2,range))
center<-apply(hdr,2,min)
hdr<-t(t(hdr)-center)
hdr<-hdr/matrix(apply(hdr,2,max),nrow=nrow(hdr),ncol=3,byrow=TRUE)
}
else{
real_range<-rep(1,3)
}
as3d<-ashape3d(hdr,alpha=alpha)
vas3d<-volume_ashape3d(as3d)
volumeHDR[ns,j,i]=vas3d*prod(real_range)
if(scale){
as3d$x<-hdrorig
}
}
}
}


# Results
# ===================================
Vol_mean<-apply(volumeHDR,c(1,3),mean)
rownames(Vol_mean)<-paste("n=",nsize,sep="")
colnames(Vol_mean)<-paste("Level ",lev,sep="")
Vol_sd<-apply(volumeHDR,c(1,3),sd)
rownames(Vol_sd)<-paste("n=",nsize,sep="")
colnames(Vol_sd)<-paste("Level ",lev,sep="")

mc_HDRvol
Vol_mean
Vol_sd
