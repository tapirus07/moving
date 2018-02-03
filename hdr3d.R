# ====================================================================================
# Function hdr3d: Estimation of the highest density region (HDR) in 3D
# ====================================================================================


# Input
# ============
# x --------- Data matrix n x 3
# h --------- Either a 3-dimensional vector with the bandwitdh parameters (h1,h2,h3) or: 
#             h = NULL ---- Normal optimal smoothing (default)
#             h = "href" ---- Same X-Y coordinates
#             h = "LSCV" ---- LSCV bandwidth (package ks)
#             h = "Plug-in" ---- Plug-in bandwidth (package ks)
#             h = "Plug-in diag" ---- Diagonal Plug-in bandwidth (package ks)
#             h = "SCV" ---- SCV bandwidth (package ks)
# grid ------ Number of grid points
# lev ------- Level of the desired highest density region (it can be a vector). By default 0.95
# B --------- Size of the sample generated from the kernel density. If B=0 (default), no sample is generated (it works with the sample points)
# opt ------- Net of points representing the HDR
#             opt=1 Original sample points (if B=0 or B<n) or sample generated from the kernel density (if B>n)
#             opt=2 Uniform sample on the range of x (by default with size nopt=50000)
#             opt=3 Regular grid
# nopt ------ When opt=2, sample size of the uniform sample on the range of x 
# zmax ------ If not NULL the HDR regions do not include the space above zmax and below 0
# ...  ------ Arguments to be passed (such as arguments to kde function)


# Output
# ============
# List with three components
#   hdr ----- List of HDR points clouds. hdr[[i]] contains the HDR point cloud for level lev[i]
#   H   ----- Bandwidth matrix
#   level --- Input level of the desired highest density region (it can be a vector)

hdr3d<-function(x,h=NULL,grid=30,lev=0.95,B=0,opt=1,nopt=5000,zmax=NULL,...){

library(ks)
library(mvtnorm)

n<-nrow(x)
###################################
##### Bandwidth parameters in 3D
###################################
{
if(!is.numeric(h)){
if(is.null(h)){
h=c(sd(x[,1])*(4/(5*n))^(1/7),sd(x[,2])*(4/(5*n))^(1/7),sd(x[,3])*(4/(5*n))^(1/7)) # Normal optimal smoothing
H<-diag(h^2)
}
else if (h=="href"){
h=c(rep(0.5*(sd(x[,1])+sd(x[,2]))*(4/(5*n))^(1/7),2),sd(x[,3])*(4/(5*n))^(1/7)) # Same scale in X-Y coordinates
H<-diag(h^2)
}
else if (h=="LSCV"){
H<-Hlscv(x)
}
else if (h=="Plug-in"){
H<-Hpi(x)
}
else if (h=="Plug-in diag"){
H<-Hpi.diag(x)
}
else if (h=="SCV"){
H<-Hscv(x)
}
}
else{
H<-diag(h^2)
}
}

###################################
##### Kernel density estimation
###################################
kest<-kde(x,H=H,gridsize=rep(grid,3),...) # 3D-kernel estimation
fn<-kest$estimate # fn(x), for x in the grid

###################################
##### Highest density region
###################################
if(n>B){ # For large sample size y=fn(x) can be used directly to estimate the quantile 
rkest<-x
fnx<-predict(kest,x=rkest) # fn(x), for the original x
qfn<-quantile(fnx,1-lev)
}
else{ # If the sample size is moderate, it is preferable to generate observations from fn
x.ind <- sample(1:n, size = B, replace = TRUE)
rkest <- x[x.ind,] + rmvnorm(B,sigma=H) # Realization from fn (Silverman, page 143)
fnx<-predict(kest,x=rkest) # fn(x), for the generated x
qfn<-quantile(fnx,1-lev)
}

# HDR from generated points from fn (or original sample points if n large)
if(opt==2){
rg<-apply(rkest,2,range)
if(!is.null(zmax)){
rkest<-cbind(runif(nopt,rg[1,1],rg[2,1]),runif(nopt,rg[1,2],rg[2,2]),runif(nopt,0,zmax))
}
else{
rkest<-cbind(runif(nopt,rg[1,1],rg[2,1]),runif(nopt,rg[1,2],rg[2,2]),runif(nopt,rg[1,3],rg[2,3]))
}
fnx<-predict(kest,x=rkest) # fn(x), for the generated x
}
if(opt==3){
nopt<-grid
rg<-apply(rkest,2,range)
if(!is.null(zmax)){
rkest<-as.matrix(expand.grid(seq(rg[1,1],rg[2,1],length=nopt),seq(rg[1,2],rg[2,2],length=nopt),seq(0,zmax,length=nopt)))
}
else{
rkest<-as.matrix(expand.grid(seq(rg[1,1],rg[2,1],length=nopt),seq(rg[1,2],rg[2,2],length=nopt),seq(rg[1,3],rg[2,3],length=nopt)))
}
fnx<-predict(kest,x=rkest) # fn(x), for the generated x
}

hdr<-list()
for (k in 1:length(qfn)){
zlev<-which(fnx>=qfn[k])
aux<-rkest[zlev,]
if(opt==3){aux<-aux+rmvnorm(nrow(aux),sigma=diag(10^-6*apply(aux,2,sd)))}
if(!is.null(zmax)){
aux<-aux[(aux[,3]>0&aux[,3]<zmax),]
}
hdr[[k]]<-aux
}
return(list("HDR"=hdr,"H"=H,"level"=lev))
}



