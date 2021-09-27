library(png)
library('spatialfil')
png('blur.png',width=1000,height=250)
par(mfrow=c(1,3),mar=c(1,0,1,0))
mat=readPNG("blue.png")[,,1]
nrmlz=function(mat)
{
  tot=sum(mat)
  mat/tot
}

b1=nrmlz(applyFilter(mat,convKernel(sigma = 1,k="gaussian")))
b2=nrmlz(applyFilter(mat,convKernel(sigma = 8,k="gaussian")))
b3=nrmlz(applyFilter(mat,convKernel(sigma = 20,k="gaussian")))

mid=dim(b1)[1]/2
b1.s=b1[(mid-30):(mid+30),(mid-200):(mid+200)]
b2.s=b2[(mid-30):(mid+30),(mid-200):(mid+200)]
b3.s=b3[(mid-30):(mid+30),(mid-200):(mid+200)]
d=dim(b1.s)

image(t(b1.s[rev(1:d[1]),]),axes=F,col=rgb(1,(0:256)/256,(0:256)/256))
image(t(b2.s[rev(1:d[1]),]),axes=F,col=rgb(1,(0:256)/256,(0:256)/256))
image(t(b3.s[rev(1:d[1]),]),axes=F,col=rgb(1,(0:256)/256,(0:256)/256))
dev.off()