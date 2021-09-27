library(diagram)
bColVal=c('lightyellow','white','salmon1')
bDens=c(0,0,20)

drawAOV4=function(){
  names=c("Y~A+B+A:B","Y~A+A:B","Y~A+B","Y~B+A:B")
  M=matrix(nrow=length(names),ncol=length(names),data=0)
  MM <- as.data.frame(M)
  MM[[1,2]] <- "F[B]"
  MM[[1,3]] <- "F[AB]"
  MM[[1,4]] <- "F[A]"
  plotmat(MM,pos=c(1,3),name=names,curve=0,box.type="round",box.size=.04,
          box.prop=1,box.col=bColVal[1],arr.length=0)}


drawAOV8=function(){
  bColI=c(1,3,1,3,1,3,1,1)
  names=c("Y~A+B+A:B","Y~A+A:B","Y~A+B","Y~B+A:B","Y~A","Y~A:B","Y~B","Y~.")
  M=matrix(nrow=length(names),ncol=length(names),data=0)
  plotmat(M,pos=c(1,3,3,1),name=names,curve=0,box.type="round",
          box.size=.04,box.prop=1,box.col=bColVal[bColI],arr.length=0)
}


drawIceCream=function(){
  par(mfrow=c(2,2))
  s=rep(1:2,2)
  f=rep(1:2,each=2)
  data=c(5,2.5,2.5,5)
  dat=tapply(data,list(f,s),mean)
  matplot(dat,typ='l',ylab="Price Offered ($)",xlab="Fat Content",axes=F,col="black",lty=1,ylim=c(0,6))
  matpoints(dat,pch=21,cex=2,bg=c('darkred','white'),col='black')
  axis(2)
  axis(1,at=1:2,label=c("Low","High"))
  legend(1.35,6,c("High Sugar","Low Sugar"),pch=22,pt.bg=c("white","darkred"),pt.cex=1)
  mtext(side=3,adj=0,"A.")
  
  fat=seq(0,10,.1)
  sugar=seq(0,10,.1)
  
  myfun = function(x,y) dmvnorm(cbind(x,y),c(5,5),matrix(c(3,1.3,1.3,2),nrow=2))
  a=outer(fat,sugar,myfun)
  image(fat,sugar,a,xlab="Sugar Content",ylab="Fat Content",col=heat.colors(100,alpha=.5))
  contour(fat,sugar,a,add=T,labels="")
  points(c(4,6,4,6),c(4,4,6,6),pch=19)
  abline(h=c(4,6),lty=2)
  abline(v=c(4,6),lty=2)
  mtext(side=3,adj=0,"B.")
  
  slope=2.5/(a[sugar==4,fat==4]-a[sugar==4,fat==6])
  intercept=5-a[sugar==4,fat==4]*slope
  a2dollar=function(a) slope*a+intercept
  d=c(a[sugar==4,fat==4],
      a[sugar==6,fat==4],
      a[sugar==4,fat==7],
      a[sugar==6,fat==7])
  data=a2dollar(d)
  dat=tapply(data,list(f,s),mean)
  matplot(dat,typ='l',ylab="Price Offered ($)",xlab="Fat Content",axes=F,col="black",lty=1,ylim=c(0,6))
  matpoints(dat,pch=21,cex=2,bg=c('darkred','white'),col='black')
  axis(2)
  axis(1,at=1:2,label=c("Low","High"))
  mtext(side=3,adj=0,"C.")
  
  image(fat,sugar,a,xlab="Sugar Content",ylab="Fat Content",col=heat.colors(100,alpha=.5))
  contour(fat,sugar,a,add=T,labels="")
  points(c(4,7,4,7),c(4,4,6,6),pch=19)
  abline(h=c(4,6),lty=2)
  abline(v=c(4,7),lty=2)
  mtext(side=3,adj=0,"D.")
}


makeSmallEffect=function(){
  set.seed(1234)
  A=2
  B=2
  S=20
  evenF=.218
  t.a=c(-.5,.5)*evenF
  t.b=c(-.5,.5)*evenF
  t.ab=matrix(c(-.5,.5,.5,-.5)*evenF*2,nrow=2)
  a=as.factor(rep(1:A,each=S*B))
  b=as.factor(rep(rep(1:B,each=S),A))
  sub=as.factor(rep(1:(S*A*B)))
  t.su=rnorm(S*A*B)
  m=tapply(t.su,list(a,b),mean)
  t.s=t.su-m[cbind(a,b)]
  y=100*(t.s+t.a[a]+t.b[b]+t.ab[cbind(a,b)])+700
  return(data.frame(a,b,sub,y))
}

makeSSE=function(dat){
  y=dat$y
  A=2
  B=2
  S=20
  X0=rep(1,S*A*B)
  Xa=2*(as.integer(dat$a)-1)-1
  Xb=2*(as.integer(dat$b)-1)-1
  Xi=Xa*Xb
  SSE=1:8
  SSE[1]=sum(lm.fit(x=cbind(X0),y=y)$residual^2)
  SSE[2]=sum(lm.fit(x=cbind(X0,Xa),y=y)$residual^2)
  SSE[3]=sum(lm.fit(x=cbind(X0,Xb),y=y)$residual^2)
  SSE[4]=sum(lm.fit(x=cbind(X0,Xi),y=y)$residual^2)
  SSE[5]=sum(lm.fit(x=cbind(X0,Xa,Xb),y=y)$residual^2)
  SSE[6]=sum(lm.fit(x=cbind(X0,Xa,Xi),y=y)$residual^2)
  SSE[7]=sum(lm.fit(x=cbind(X0,Xb,Xi),y=y)$residual^2)
  SSE[8]=sum(lm.fit(x=cbind(X0,Xa,Xb,Xi),y=y)$residual^2)
  return(SSE)
}

makeAIC=function(dat){
  SSE=makeSSE(dat)
  k=c(1,2,2,2,3,3,3,4)
  n=length(dat$y)
  return(n*log(SSE/n)+2*k)}

makeBIC=function(dat){
  SSE=makeSSE(dat)
  k=c(1,2,2,2,3,3,3,4)
  n=length(dat$y)
  return(n*log(SSE/n)+k*log(n))}