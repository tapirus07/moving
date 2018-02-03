# ====================================================
# Example of estimation of the volume of HDR 
# ====================================================

source("hdr3d.R") # Function that estimates the HDR
library(alphashape3d)

# Opossum data
# =======================================================================
data<-read.table("opossum_locations.txt",T)
id<-as.vector(data[,4])
op<-as.matrix(data[,1:3])
ind<-unique(id)

# Input parameters
# ===========================
lev<-c(0.95,0.5) # Level of the desired highest density region (it can be a vector)
grid<-100         # Number of grid points
B<-10000         # Size of the sample generated from the kernel density to compute the quantile for the HDR estimation
h<-hs           # h=NULL corresponds to the Normal optimal smoothing (see Input in function hdr3d)
#h="href"
opt<-1           # The HDR regions are computed in 1=B points, 2=uniform sample, 3=grid.
canopy<-25       # If not NULL the HDR regions are restricted to the interval [0,canopy] in the z-axis.

alpha<-0.1       # Value of alpha for the alpha-shape reconstruction
scale<-TRUE      # Set scale<-TRUE if the scales are very different in the three dimensions (for the alpha-shape computation)
plotHDR<-TRUE    # Do you want to plot the results?
aspect<-TRUE     # Do you want to adjusts the aspect ratio in the plot?

# Estimation for each ID
# ===========================
volumeHDRopp<-list() # volumeHDRopp[[i]] contains the estimated volume of the HDR of animals i
hdrID<-list()
ashape=list()
for(k in 1:length(ind)){
x<-op[id==ind[k],] # Data set

# Estimation of the highest density region (HDR) in 3D
# ====================================================================================
hdrID[[k]]<-hdr3d(x,h=h,grid=grid,lev=lev,B=B,zmax=canopy,opt=opt)
hdrout<-hdrID[[k]]

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
volumeHDRopp[[k]]<-volumeHDR


# Plots
# ===================================
# The plot shows:
#   - Original sample points in black
#   - Point clouds for the HDR in different colors 
#   - alpha-shapes for the different levels
# If aspect=FALSE it shows the original scale
# If aspect=TRUE  it adjusts the aspect ratio

if(plotHDR){
open3d()
plot3d(x,xlab="x",ylab="y",zlab="z",aspect=aspect,main=ind[k]) #Relocations
for (i in 1:length(lev)){
plot3d(hdrout$HDR[[i]],col=(i+1),aspect=aspect,add=TRUE) # HDR point clouds 
plot(as3d[[i]],clear=FALSE,col=(i+1),transparency=0.4)   # alpha-shape
}
}
ashape[[k]]=as3d
}

# Volumes of HDR for each ID
# ============================
names(volumeHDRopp)<-ind
volumeHDRopp

# Bandwith matrices for each ID:  H=diag(h1^2,h2^2,h3^2)
# =======================================================
for (k in 1:length(ind)){
print(ind[k])
print(hdrID[[k]]$H)
}

#MEAN OF "HREF"
hs=apply(sqrt(do.call("rbind",lapply(hdrID,function(x){colSums(x$H)}))),2,median)


# Volumes ASHAPE OVERLAPING
# ============================
iso=1
alpha=0.1
over=matrix(0,length(ashape),length(ashape))
diag(over)=1
colnames(over)=rownames(over)=ind
for(i in 1:length(ashape)){
	for(j in c(1:length(ashape))[1:length(ashape)>i]){
	#GETTING HDR POINTS FROM TWO INDIVIDUALS
	x<-rbind(data.frame(hdrID[[i]]$HDR[[iso]]),data.frame(hdrID[[j]]$HDR[[2]]))
	#CHECKING FOR HDR POINTS INSIDE ASHAPE
	in1<-inashape3d(ashape[[i]][[iso]], indexAlpha = 1, as.matrix(x))
	in2<-inashape3d(ashape[[j]][[iso]], indexAlpha = 1, as.matrix(x))
	xin<-x[in1&in2,]
	if(nrow(xin)==0){over[i,j]=0} #IF ANY HDR POINTS FALL INSIDE ASHAPE
	else{
	#RESCALING xin POINTS FOR ASHAPE COMPUTATION
	real_range<-diff(apply(xin,2,range))
	center<-apply(xin,2,min)
	xin<-t(t(xin)-center)
	xin<-xin/matrix(apply(xin,2,max),nrow=nrow(xin),ncol=3,byrow=TRUE)
	a3<-ashape3d(xin,alpha=alpha)
	over[i,j]=volume_ashape3d(a3)*prod(real_range)/volumeHDRopp[[i]][iso] # Interseciton volume
	over[j,i]=volume_ashape3d(a3)*prod(real_range)/volumeHDRopp[[j]][iso]
	#FOR PLOTING
	a3$x=x[in1&in2,]
	open3d()
	plot3d(x)
	title3d(paste(i,j))
	plot(ashape[[i]][[iso]],add=T,clear=F,transparency=0.4,col=2)
	plot(ashape[[j]][[iso]],add=T,clear=F,transparency=0.4,col=4)
	plot(a3,clear=FALSE,col=c(3,3,3),trans=1) # Intersection in green
	}
	}
}
over50=over
over95=over


#BIDIMENSIONAL OVERLAP
library(adehabitatHR)
k=kernelUD(SpatialPointsDataFrame(op[,1:2],data=data.frame(id=id)))
hss=unlist(lapply(k,function(x){x@h$h}))
k=kernelUD(SpatialPointsDataFrame(op[,1:2],data=data.frame(id=id)),h=mean(hss),ex=2,same4all=T,grid=200)
pol95=getverticeshr(k,95)
pol50=getverticeshr(k,50)
plot(pol95,col=1:6)
plot(pol50,col=1:6)
overbi95=kerneloverlaphr(k,method="HR",percent=95)
overbi50=kerneloverlaphr(k,method="HR",percent=50)
overbi95=overbi95[ind,ind]
overbi50=overbi50[ind,ind]

par(mfrow=c(1,2))
diag(overbi95)=diag(over95)=diag(overbi50)=diag(over50)=NA
plot(as.vector(overbi95),as.vector(over95),pch=16,xlab="Standard Kernel2d Overlap(%)",ylab="Kernel3D Overlap (%)")
lines(c(0,1),c(0,1))
plot(as.vector(overbi50),as.vector(over50),pch=16,xlab="Standard Kernel2d Overlap(%)",ylab="Kernel3D Overlap (%)")
lines(c(0,1),c(0,1))

cor.test(as.vector(overbi95),as.vector(over95))
cor.test(as.vector(overbi50),as.vector(over50))


summary(as.vector(over95))
summary(as.vector(over50))

