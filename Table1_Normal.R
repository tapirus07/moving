# ====================================================
# Example of estimation of the volume of HDR 
# ====================================================
set.seed(1234)
source("hdr3d.R") # Function that estimates the HDR
library(alphashape3d)

# Example 
# Standard Multivariate Normal
# =======================================================================
library(mvtnorm)          # Generate multivariate normals
mu<-c(0,0,0)              # Mean vector
Sigma<- diag(c(1,1,1))    # Covariance matrix Dataset (1)
#Sigma<- diag(c(1,3,7))   # Covariance matrix Dataset (2)


# Input parameters
# ===========================
lev<-c(0.5,0.95) # Level of the desired highest density region (it can be a vector)
grid<-20         # Number of grid points
h<-NULL          # h=NULL corresponds to the Normal optimal smoothing (see Input in hdr3d)

alpha<-10        # Value of alpha for the alpha-shape reconstruction
scale<-FALSE     # Set scale<-TRUE if the scales are very different in the three dimensions (for the alpha-shape computation)
plotHDR<-TRUE    # Do you want to plot the results?
aspect<-FALSE    # Do you want to adjusts the aspect ratio in the plot?

nsize<-c(30,50,500)       # Sample sizes
B<-3*nsize                # Size of the sample generated from the kernel density to compute the quantile for the HDR estimation
                          # B=0 -> no sample is generated (it works with the sample points)
nsamples<-100             # Number of samples


# Exact volume 
# ------------------
d<-3
exact_HDRvol<-2*pi^(d/2)/(d*gamma(d/2))*qchisq(lev,d)^(d/2)*sqrt(det(Sigma))


# Volume estimation
# ------------------
volumeHDR<-array(dim=c(length(nsize),nsamples,length(lev)),dimnames=list(NULL,NULL,c("level 0.95","level 0.5"))) 
for(ns in 1:length(nsize)){
n<-nsize[ns]

for (j in 1:nsamples){
x<-rmvnorm(n,mean=mu,sigma=Sigma)


# Estimation of the highest density region (HDR) in 3D
# ====================================================================================
hdrout<-hdr3d(x,h=h,grid=grid,lev=lev,B=B[ns])

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


exact_HDRvol
Vol_mean
Vol_sd


