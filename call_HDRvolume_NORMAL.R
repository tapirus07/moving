# ====================================================
# Example of estimation of the volume of HDR 
# ====================================================

source("hdr3d.R") # Function that estimates the HDR
library(alphashape3d)

# Example 
# Standard Multivariate Normal
# =======================================================================
library(mvtnorm) # Generate multivariate normals
n<-2000
mu<-c(0,0,0) # Mean vector
Sigma<- diag(c(1,3,7))# Covariance matrix
x<-rmvnorm(n,mean=mu,sigma=Sigma)


# Input parameters
# ===========================
lev<-c(0.95,0.5) # Level of the desired highest density region (it can be a vector)
grid<-20         # Number of grid points
B<-5000          # Size of the sample generated from the kernel density to compute the quantile for the HDR estimation
h<-NULL          # h=NULL corresponds to the Normal optimal smoothing (see Input in hdr3d)

alpha<-10        # Value of alpha for the alpha-shape reconstruction
scale<-FALSE     # Set scale<-TRUE if the scales are very different in the three dimensions (for the alpha-shape computation)
plotHDR<-TRUE    # Do you want to plot the results?
aspect<-FALSE    # Do you want to adjusts the aspect ratio in the plot?


# Estimation of the highest density region (HDR) in 3D
# ====================================================================================
hdrout<-hdr3d(x,h=h,grid=grid,lev=lev,B=B)

volumeHDR<-numeric() 
as3d<-list()
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
as3d[[i]]<-ashape3d(hdr,alpha=alpha)
vas3d<-volume_ashape3d(as3d[[i]])
volumeHDR[i]=vas3d*prod(real_range)
if(scale){
as3d[[i]]$x<-hdrorig
}
}
names(volumeHDR)<-paste("Level",lev)


# Results
# ===================================
print(volumeHDR) # Estimated volume


# Plots
# ===================================
# The plot shows:
#   - Original sample points in black
#   - Point clouds for the HDR in differnt colors 
#   - alpha-shapes for the different levels

# If aspect=FALSE it shows the original scale
# If aspect=TRUE  it adjusts the aspect ratio

if(plotHDR){
open3d()
# Original sample
plot3d(x,xlab="x",ylab="y",zlab="z",aspect=aspect) #Relocations
# HDR point cloud and alpha-shape
for (i in 1:length(lev)){
plot3d(hdrout$HDR[[i]],col=(i+1),add=TRUE) 
plot(as3d[[i]],clear=FALSE,col=(i+1),transparency=0.4)
}
}


