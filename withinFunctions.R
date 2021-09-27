library(BayesFactor)
library(diagram)

drawWithinMods8=function(){
  bColVal=c("white","darkred","red")
  textCol=c("black","white","black")
  bColI=c(1,2,1,3,2,2,1,2)
  names=c("Y~A+S+A:S","Y~A+A:S","Y~A+S","Y~S+A:S","Y~A","Y~A:S","Y~S","Y~.")
  M=matrix(nrow=length(names),ncol=length(names),data=0)
  MM <- as.data.frame(M)
  MM[[1,4]] <- "F"
  MM[[3,7]] <- "1"
  MM[[1,7]] <- "2"
  plotmat(MM,pos=c(1,3,3,1),name=names,curve=0,box.type="round",box.size=.04,
          box.prop=1,box.col=bColVal[bColI],
          arr.length=0,txt.col=textCol[bColI])
  
}


drawTheta=function(){
  theta=seq(-60,160,1)
  clear=dnorm(theta,60,20)
  null=dnorm(theta,0,20)
  top=1.2*max(clear)
  plot(theta,clear,typ='l',ylim=c(0,top),axes=F,
       ylab="",xlab=expression(paste("Individual's True Effects ",theta," (ms)")))
  #lines(theta,null)
  axis(1)
  text(60,.9*top,"Y~A+S+A:S",adj=-.1)
  abline(v=60)
  text(62,.1*top,adj=-.2,"A")
  shape::Arrows(64,.45*top,80,.45*top,code=3)
  text(90,.45*top,adj=-.2,"A:S")
  lines(theta,null,lty=2)
  text(0,.9*top,"Y~S+A:S",adj=1.1)
}

drawBlurMods=function(){
  par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2,1,0))
  theta=seq(-60,160,1)
  clear=dgamma(theta,shape=2,scale=30)
  blur=dgamma(theta,shape=2,scale=10)
  top=1.5*max(blur)
  plot(theta,blur,typ='l',ylim=c(0,top),axes=F,
       ylab="",xlab=expression(paste("Individual Effects, ",theta)))
  lines(theta,clear)
  shape::Arrows(0,0,0,.9*top)
  axis(1)
  text(0,top,"Full Blur")
  text(60,.5*top,"Some Blur")
  text(100,.2*top,"Clear")
  mtext(side=3,adj=0,cex=1.2,"A",line=-1)
  
  theta=seq(-60,160,1)
  clear=dexGAUS(theta, mu = 0, sigma = 20, nu = 40)
  blur=dexGAUS(theta, mu = 0, sigma = 5, nu = 10)
  top=1.5*max(blur)
  plot(theta,blur,typ='l',ylim=c(0,top),axes=F,
       ylab="",xlab=expression(paste("Individual Effects, ",theta)))
  lines(theta,clear)
  shape::Arrows(0,0,0,.9*top)
  axis(1)
  text(0,top,"Full Blur")
  text(50,.5*top,"Some Blur")
  text(100,.1*top,"Clear")
  mtext(side=3,adj=0,cex=1.2,"B",line=-1)
  
  theta=seq(-60,160,1)
  clear=dnorm(theta,60,20)
  blur=dnorm(theta,30,20)
  null=dnorm(theta,0,20)
  top=1.2*max(blur)
  plot(theta,blur,typ='l',ylim=c(0,top),axes=F,
       ylab="",xlab=expression(paste("Individual Effects, ",theta)))
  lines(theta,clear)
  lines(theta,null)
  axis(1)
  text(-18,.7*top,"Full\n Blur",adj=1)
  text(30,.9*top,"Some Blur")
  text(75,.7*top,"Clear",adj=0)
  mtext(side=3,adj=0,cex=1.2,"C",line=-1)}


makeData=function(t.beta,t.gamma.s,t.s){
  t.mu=1000
  t.alpha=rnorm(I,0,200)
  t.beta=t.beta
  t.gamma=matrix(ncol=J,rnorm(I*2,0,t.gamma.s))
  sub=rep(1:I,each=J*L)
  cond=rep(rep(1:J,each=L),I)
  t.mean=t.mu+t.alpha[sub]+t.beta*(as.integer(cond)-3/2)+t.gamma[cbind(sub,cond)]
  y=rnorm(N,t.mean,t.s)
  dat=data.frame(sub,cond,y)
  return(dat)
}

Q=function(a){
  I=diag(a)
  J=matrix(1,nrow=a,ncol=a)
  S=I-J/a
  E=eigen(S)
  return((E$vectors[1:a,1:(a-1)]))
}

circMult=function(Xa,Xb){
  N=dim(Xa)[1]
  a=dim(Xa)[2]
  b=dim(Xb)[2]
  X=matrix(nrow=N,ncol=a*b)
  for (n in 1:N) X[n,]=as.vector(t(outer(Xa[n,],Xb[n,])))
  return(X)
}

bf8=function(dat,scale){
  I=max(dat$sub)
  J=max(dat$cond)
  N=length(dat$sub)
  XsubR=matrix(0,nrow=N,ncol=I)
  for (n in 1:N) XsubR[n,dat$sub[n]]=1
  XsubF=XsubR%*%Q(I)
  XcondR=matrix(0,nrow=N,ncol=J)
  for (n in 1:N) XcondR[n,dat$cond[n]]=1
  XcondF=XcondR%*%Q(J)
  
  # Subject Only
  X=cbind(XsubR)
  gmap=rep(0,I)
  bfSub=nWayAOV(dat$y,X,gmap,rscale=scale[1])
  
  # Additive Fixed
  X=cbind(XsubR,XcondF)
  gmap=rep(0:1,c(I,1))
  #chain=nWayAOV(y,X,gmap,rscale=c(1,1,1),posterior=TRUE)
  bfAdd=nWayAOV(dat$y,X,gmap,rscale=scale[1:2])
  
  # Random Interaction
  XsubcondR=circMult(XcondR,XsubR)
  X=cbind(XsubR,XcondR,XsubcondR)
  gmap=rep(0:2,c(I,2,I*J))
  bfRan=nWayAOV(dat$y,X,gmap,rscale=scale)
  X=cbind(XsubR,XsubcondR)
  gmap=rep(0:1,c(I,I*J))
  bfRan0=nWayAOV(dat$y,X,gmap,rscale=scale[c(1,3)])
  
  # Fixed 
  XsubcondF=circMult(XcondF,XsubF)
  X=cbind(XsubF,XcondF,XsubcondF)
  gmap=rep(0:2,c(I-1,1,(I-1)*(J-1)))
  bfFix=nWayAOV(dat$y,X,gmap,rscale=scale)
  X=cbind(XsubR,XsubcondF)
  gmap=rep(0:1,c(I,(I-1)*(J-1)))
  bfFix0=nWayAOV(dat$y,X,gmap,rscale=scale[c(1,3)])
  
  #CT
  XsubcondCT=circMult(XcondF,XsubR)
  X=cbind(XsubR,XcondF,XsubcondCT)
  gmap=rep(0:2,c(I,1,(I*(J-1))))
  bfCT=nWayAOV(dat$y,X,gmap,rscale=scale)
  X=cbind(XsubR,XsubcondCT)
  gmap=rep(0:1,c(I,I*(J-1)))
  bfCT0=nWayAOV(dat$y,X,gmap,rscale=scale[c(1,3)])
  myBF=c(bfSub$bf,bfAdd$bf,bfRan$bf,bfFix$bf,bfCT$bf,bfRan0$bf,bfFix0$bf,bfCT0$bf)
  names(myBF)=c("Y~S","Y~S+A","Y~S+A+A:S (ran)","Y~S+A+A:S (fix)","Y~S+A+A:S (CT)","Y~S+A:S (ran)","Y~S+A:S (fix)","Y~S+A:S (CT)")
  return(myBF)
}

myAOV=function(dat){
  dat$sub=as.factor(dat$sub)
  dat$cond=as.factor(dat$cond)
  g=summary(aov(y~cond+Error(sub/cond),data=dat))
  return(g)}  
