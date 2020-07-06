#Análise de influência local para os modelos multivariados com regressao

#Pacotes e variaveis globais

#library(MASS)
#library(expm)
#library(moments)
#library(abind)
#library(mvtnorm)
#require(MCMCpack)
#require(doParallel)

#n=300
#p=2
#q=3
#p0=p*(p+1)/2
#beta=c(4,-3,5)
#rSigma=matrix(c(2.5,-1,-1,1.5),p,p) 
#Sigma=rSigma%*%rSigma
#lambda=c(2,-1)
#nu=3.5
#eta=0.25
#gamma=0.4

#P<-c(beta,vech(rSigma),lambda,nu) #stn ou ssn
#P<-c(beta,vech(rSigma),lambda,eta,gamma) #scn
#r<-length(P)

#Perturbacoes

#pertres<-function(y,w,D){
# n=nrow(y)
# p=ncol(y)
# Sy=sqrt(diag(cov(y)))
# yw<-matrix(0,n,p)
# for(i in 1:n){  
# yw[i,]<-y[i,]+w[i]*diag(D)%*%Sy
# } 
# return(yw)  
#} #perturbacao da resposta

#pertexp<-function(X,w,V){
#  p=nrow(X)
#  q=ncol(X[[n]])
#  n=length(X)/(nrow(X)*ncol(X))
#  Xw<-array(0,c(p,q,n))
#  E%*%Mx<-matrix(0,p,q)
#  E%*%Mx[V[1],V[2]]<-1
#  for(i in 1:n){  
#    Xw[,,i]<-X[[i]]+w[i]*E%*%Mx
#  } 
#  return(Xw)  
#} #perturbacao na variavel explicativa de ordem V


#Modelo normal

hesqmnr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=matrix(0,p,1)
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  d2qfdb2<-matrix(0,q,q)
  d2qfdbda<-matrix(0,q,p0)
  d2qfdbdh<-matrix(0,q,p)
  d2qfda2<-matrix(0,p0,p0)
  d2qfdhda<-matrix(0,p,p0)
  d2qfdh2<-matrix(0,p,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-1
    auxi<-pmax(t(h)%*%invB%*%e[i,],-37)
    t[i]<-0
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    d2qfdb2<-d2qfdb2-t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%X[[i]]
    d2qfdbdh<-d2qfdbdh-as.numeric(t[i]-Ai)*t(X[[i]])%*%invB+t(X[[i]])%*%invB%*%h%*%t(e[i,])%*%invB
    d2qfdh2<-d2qfdh2-invB%*%e[i,]%*%t(e[i,])%*%invB
    dAida<-matrix(0,p0,1)
    d2Aidhda<-matrix(0,p,p0) 
    d2didbda<-matrix(0,q,p0)
    d2Aidbda<-matrix(0,q,p0)
    d2ldda2<-matrix(0,p0,p0)
    d2dida2<-matrix(0,p0,p0)
    d2Aida2<-matrix(0,p0,p0)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      d2Aidhda[,j]=-invB%*%Bj%*%invB%*%e[i,]
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidbda[,j]=t(X[[i]])%*%invB%*%Bj%*%invB%*%h
      d2didbda[,j]=2*t(X[[i]])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]
      for(k in 1:p0){
        ind=rep(0,p0)
        ind[k]=1
        Bk=xpnd(ind)
        d2ldda2[j,k]=-sum(diag(invB%*%Bk%*%invB%*%Bj))
        d2Aida2[j,k]=-t(h)%*%invB%*%(Bk%*%invB%*%Bj+Bj%*%invB%*%Bk)%*%invB%*%e[i,]
        d2dida2[j,k]=t(e[i,])%*%invB%*%(Bk%*%invB%*%Bj%*%invB+Bj%*%invB%*%Bk%*%invB+Bj%*%invB%*%invB%*%Bk+Bk%*%invB%*%invB%*%Bj+invB%*%Bk%*%invB%*%Bj+invB%*%Bj%*%invB%*%Bk)%*%invB%*%e[i,]
      }
    }
    d2qfdhda<-d2qfdhda+as.numeric(t[i]-Ai)*d2Aidhda-invB%*%e[i,]%*%t(dAida)
    d2qfdbda<-d2qfdbda-1/2*as.numeric(u[i])*d2didbda+as.numeric(t[i]-Ai)*d2Aidbda+t(X[[i]])%*%invB%*%h%*%t(dAida) 
    d2qfda2<-d2qfda2-1/2*d2ldda2-1/2*as.numeric(u[i])*d2dida2+as.numeric(t[i]-Ai)*d2Aida2-dAida%*%t(dAida)
  }
  d2qfdb=cbind(d2qfdb2,d2qfdbda)
  d2qfda=cbind(t(d2qfdbda),d2qfda2)
  d2qfdh=cbind(t(d2qfdbdh),d2qfdhda,d2qfdh2,matrix(0,p,1))
  Hqf=rbind(d2qfdb,d2qfda)
  return(Hqf)
}

Dw0ponmnr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=matrix(0,p,1)
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-1
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-0
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%e[i,]-as.numeric(t[i])*t(X[[i]])%*%invB%*%h
    ddida<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    dldda<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dldda[j]=2*sum(diag(invB%*%Bj))
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      ddida[j]=-t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]
    }
    dqfdai=-1/2*dldda-1/2*as.numeric(u[i])*ddida+as.numeric(t[i]-Ai)*dAida
    dqfdhi=as.numeric(t[i]-Ai)*invB%*%e[i,]
    can=rep(0,n)
    can[i]=1
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai))
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta ponderacao de casos na normal multivariada com regressao

Dw0resmnr<-function(P,y,X,D)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=matrix(0,p,1)
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  Sy=sqrt(diag(cov(y)))
  dAidaw=t(h)%*%invB%*%diag(D)%*%Sy
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-1
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-0
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%diag(D)%*%Sy
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=-t(h)%*%invB%*%Bj%*%invB%*%diag(D)%*%Sy
      ddidaw[j]=-2*t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%diag(D)%*%Sy
    }
    dqfdai=-1/2*as.numeric(u[i])*ddidaw+as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=as.numeric(t[i]-Ai)*invB%*%diag(D)%*%Sy-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai))
    can=rep(0,n)
    can[i]=1
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao da resposta na normal multivariada com regressao

Dw0expmnr<-function(P,y,X,V)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=matrix(0,p,1)
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  mX<-X[[1]][1,2:(q/2)]
  for(i in 2:n){
    mX<-rbind(mX,X[[i]][1,2:(q/2)])
  }
  Mx=diag(c(0,sqrt(diag(cov(mX)))),q)
  E=matrix(0,p,q)
  E[V[1],(V[1]-1)*(q/2)+V[2]+1]<-1
  dAidaw=-t(h)%*%invB%*%E%*%Mx%*%b
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-1
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-0
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=as.numeric(u[i])*(t(E%*%Mx)%*%invB%*%invB%*%e[i,]-t(X[[i]])%*%invB%*%invB%*%E%*%Mx%*%b)-as.numeric(t[i,]-Ai)*t(E%*%Mx)%*%invB%*%h-t(X[[i]])%*%invB%*%h%*%t(h)%*%invB%*%E%*%Mx%*%b
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=t(h)%*%invB%*%Bj%*%invB%*%E%*%Mx%*%b
      ddidaw[j]=2*t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%E%*%Mx%*%b
    }
    dqfdai=-1/2*as.numeric(u[i])*ddidaw+as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=-as.numeric(t[i]-Ai)*invB%*%E%*%Mx%*%b-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai))
    can=rep(0,n)
    can[i]=1
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao nas explicativas na normal multivariada com regressao

Dw0escmnr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=matrix(0,p,1)
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-1
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-0
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%e[i,]-1/2*as.numeric(t[i]-Ai)*t(X[[i]])%*%invB%*%h
    dAidaw=1/2*t(h)%*%invB%*%e[i,]
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=-1/2*t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      ddidaw[j]=-t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]
    }
    dqfdai=-1/2*as.numeric(u[i])*ddidaw+as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=1/2*as.numeric(t[i]-Ai)*invB%*%e[i,]-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai))
    can=rep(0,n)
    can[i]=1
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
}


#Modelo skew normal

hesqmsnr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(length(P))])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  d2qfdb2<-matrix(0,q,q)
  d2qfdbda<-matrix(0,q,p0)
  d2qfdbdh<-matrix(0,q,p)
  d2qfda2<-matrix(0,p0,p0)
  d2qfdhda<-matrix(0,p,p0)
  d2qfdh2<-matrix(0,p,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-1
    auxi<-pmax(t(h)%*%invB%*%e[i,],-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    d2qfdb2<-d2qfdb2-t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%X[[i]]
    d2qfdbdh<-d2qfdbdh-as.numeric(t[i]-Ai)*t(X[[i]])%*%invB+t(X[[i]])%*%invB%*%h%*%t(e[i,])%*%invB
    d2qfdh2<-d2qfdh2-invB%*%e[i,]%*%t(e[i,])%*%invB
    dAida<-matrix(0,p0,1)
    d2Aidhda<-matrix(0,p,p0) 
    d2didbda<-matrix(0,q,p0)
    d2Aidbda<-matrix(0,q,p0)
    d2ldda2<-matrix(0,p0,p0)
    d2dida2<-matrix(0,p0,p0)
    d2Aida2<-matrix(0,p0,p0)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      d2Aidhda[,j]=-invB%*%Bj%*%invB%*%e[i,]
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidbda[,j]=t(X[[i]])%*%invB%*%Bj%*%invB%*%h
      d2didbda[,j]=2*t(X[[i]])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]
      for(k in 1:p0){
        ind=rep(0,p0)
        ind[k]=1
        Bk=xpnd(ind)
        d2ldda2[j,k]=-sum(diag(invB%*%Bk%*%invB%*%Bj))
        d2Aida2[j,k]=-t(h)%*%invB%*%(Bk%*%invB%*%Bj+Bj%*%invB%*%Bk)%*%invB%*%e[i,]
        d2dida2[j,k]=t(e[i,])%*%invB%*%(Bk%*%invB%*%Bj%*%invB+Bj%*%invB%*%Bk%*%invB+Bj%*%invB%*%invB%*%Bk+Bk%*%invB%*%invB%*%Bj+invB%*%Bk%*%invB%*%Bj+invB%*%Bj%*%invB%*%Bk)%*%invB%*%e[i,]
      }
    }
    d2qfdhda<-d2qfdhda+as.numeric(t[i]-Ai)*d2Aidhda-invB%*%e[i,]%*%t(dAida)
    d2qfdbda<-d2qfdbda-1/2*as.numeric(u[i])*d2didbda+as.numeric(t[i]-Ai)*d2Aidbda+t(X[[i]])%*%invB%*%h%*%t(dAida) 
    d2qfda2<-d2qfda2-1/2*d2ldda2-1/2*as.numeric(u[i])*d2dida2+as.numeric(t[i]-Ai)*d2Aida2-dAida%*%t(dAida)
  }
  d2qfdb=cbind(d2qfdb2,d2qfdbda,d2qfdbdh)
  d2qfda=cbind(t(d2qfdbda),d2qfda2,t(d2qfdhda))
  d2qfdh=cbind(t(d2qfdbdh),d2qfdhda,d2qfdh2)
  Hqf=rbind(d2qfdb,d2qfda,d2qfdh)
  return(Hqf)
}

Dw0ponmsnr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-1
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%e[i,]-as.numeric(t[i])*t(X[[i]])%*%invB%*%h
    ddida<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    dldda<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dldda[j]=2*sum(diag(invB%*%Bj))
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      ddida[j]=-t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]
    }
    dqfdai=-1/2*dldda-1/2*as.numeric(u[i])*ddida+as.numeric(t[i]-Ai)*dAida
    dqfdhi=as.numeric(t[i]-Ai)*invB%*%e[i,]
    can=rep(0,n)
    can[i]=1
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi))
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta ponderacao de casos na normal assimetrica multivariada com regressao

Dw0resmsnr<-function(P,y,X,D)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  Sy=sqrt(diag(cov(y)))
  dAidaw=t(h)%*%invB%*%diag(D)%*%Sy
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-1
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%diag(D)%*%Sy
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=-t(h)%*%invB%*%Bj%*%invB%*%diag(D)%*%Sy
      ddidaw[j]=-2*t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%diag(D)%*%Sy
    }
    dqfdai=-1/2*as.numeric(u[i])*ddidaw+as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=as.numeric(t[i]-Ai)*invB%*%diag(D)%*%Sy-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi))
    can=rep(0,n)
    can[i]=1
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao da resposta na normal assimetrica multivariada com regressao

Dw0expmsnr<-function(P,y,X,V)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  mX<-X[[1]][1,2:(q/2)]
  for(i in 2:n){
    mX<-rbind(mX,X[[i]][1,2:(q/2)])
  }
  Mx=diag(c(0,sqrt(diag(cov(mX)))),q)
  E=matrix(0,p,q)
  E[V[1],(V[1]-1)*(q/2)+V[2]+1]<-1
  dAidaw=-t(h)%*%invB%*%E%*%Mx%*%b
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-1
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=as.numeric(u[i])*(t(E%*%Mx)%*%invB%*%invB%*%e[i,]-t(X[[i]])%*%invB%*%invB%*%E%*%Mx%*%b)-as.numeric(t[i,]-Ai)*t(E%*%Mx)%*%invB%*%h-t(X[[i]])%*%invB%*%h%*%t(h)%*%invB%*%E%*%Mx%*%b
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=t(h)%*%invB%*%Bj%*%invB%*%E%*%Mx%*%b
      ddidaw[j]=2*t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%E%*%Mx%*%b
    }
    dqfdai=-1/2*as.numeric(u[i])*ddidaw+as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=-as.numeric(t[i]-Ai)*invB%*%E%*%Mx%*%b-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi))
    can=rep(0,n)
    can[i]=1
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao nas explicativas na normal assimetrica multivariada com regressao

Dw0escmsnr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-1
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%e[i,]-1/2*as.numeric(t[i]-Ai)*t(X[[i]])%*%invB%*%h
    dAidaw=1/2*t(h)%*%invB%*%e[i,]
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=-1/2*t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      ddidaw[j]=-t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]
    }
    dqfdai=-1/2*as.numeric(u[i])*ddidaw+as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=1/2*as.numeric(t[i]-Ai)*invB%*%e[i,]-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi))
    can=rep(0,n)
    can[i]=1
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao da escala na normal assimetrica multivariada com regressao

Dw0assmsnr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-1
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%h%*%t(h)%*%invB%*%e[i,]-as.numeric(t[i]-Ai)*t(X[[i]])%*%invB%*%h
    dAidaw=t(h)%*%invB%*%e[i,]
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
    }
    dqfdai=as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=as.numeric(t[i]-Ai)*invB%*%e[i,]-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi))
    can=rep(0,n)
    can[i]=1
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao da assimetria na normal assimetrica multivariada com regressao


#Modelo skew t

Dw0ponmstr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(r-1)])
  v=as.numeric(P[r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  lu<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-(v+p)/(v+d[i])
    lu[i]<-digamma((v+p)/2)-log((v+d[i])/2)
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%e[i,]-as.numeric(t[i])*t(X[[i]])%*%invB%*%h
    ddida<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    dldda<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dldda[j]=2*sum(diag(invB%*%Bj))
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      ddida[j]=-t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]
    }
    dqfdai=-1/2*dldda-1/2*as.numeric(u[i])*ddida+as.numeric(t[i]-Ai)*dAida
    dqfdhi=as.numeric(t[i]-Ai)*invB%*%e[i,]
    dqfdvi=1/2*(log(v/2)-digamma(v/2)+1-u[i]+lu[i])
    can=rep(0,n)
    can[i]=1
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi),as.matrix(dqfdvi))
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta ponderacao de casos na t assimetrica multivariada com regressao

#Dw0ponmstr(theta_st,y_st,X_st) #teste

Dw0resmstr<-function(P,y,X,D)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(r-1)])
  v=as.numeric(P[r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  Sy=sqrt(diag(cov(y)))
  dAidaw=t(h)%*%invB%*%diag(D)%*%Sy
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-(v+p)/(v+d[i])
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%diag(D)%*%Sy
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=-t(h)%*%invB%*%Bj%*%invB%*%diag(D)%*%Sy
      ddidaw[j]=-2*t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%diag(D)%*%Sy
    }
    dqfdai=-1/2*as.numeric(u[i])*ddidaw+as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=as.numeric(t[i]-Ai)*invB%*%diag(D)%*%Sy-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi),as.matrix(0))
    can=rep(0,n)
    can[i]=1
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao da resposta na t assimetrica multivariada com regressao

#Dw0resmstr(theta_st,y_st,X_st,c(1,1)) #teste

Dw0expmstr<-function(P,y,X,V)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(r-1)])
  v=as.numeric(P[r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  mX<-X[[1]][1,2:(q/2)]
  for(i in 2:n){
    mX<-rbind(mX,X[[i]][1,2:(q/2)])
  }
  Mx=diag(c(0,sqrt(diag(cov(mX)))),q)
  E=matrix(0,p,q)
  E[V[1],(V[1]-1)*(q/2)+V[2]+1]<-1
  dAidaw=-t(h)%*%invB%*%E%*%Mx%*%b
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-(v+p)/(v+d[i])
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=as.numeric(u[i])*(t(E%*%Mx)%*%invB%*%invB%*%e[i,]-t(X[[i]])%*%invB%*%invB%*%E%*%Mx%*%b)-as.numeric(t[i,]-Ai)*t(E%*%Mx)%*%invB%*%h-t(X[[i]])%*%invB%*%h%*%t(h)%*%invB%*%E%*%Mx%*%b
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=t(h)%*%invB%*%Bj%*%invB%*%E%*%Mx%*%b
      ddidaw[j]=2*t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%E%*%Mx%*%b
    }
    dqfdai=-1/2*as.numeric(u[i])*ddidaw+as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=-as.numeric(t[i]-Ai)*invB%*%E%*%Mx%*%b-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi),as.matrix(0))
    can=rep(0,n)
    can[i]=1
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao nas explicativas na t assimetrica multivariada com regressao

#Dw0expmstr(theta_st,y_st,X_st,c(1,2)) #teste

Dw0escmstr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(r-1)])
  v=as.numeric(P[r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-(v+p)/(v+d[i])
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%e[i,]-1/2*as.numeric(t[i]-Ai)*t(X[[i]])%*%invB%*%h
    dAidaw=1/2*t(h)%*%invB%*%e[i,]
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=-1/2*t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      ddidaw[j]=-t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]
    }
    dqfdai=-1/2*as.numeric(u[i])*ddidaw+as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=1/2*as.numeric(t[i]-Ai)*invB%*%e[i,]-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi),as.matrix(0))
    can=rep(0,n)
    can[i]=1
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao da escala na t assimetrica multivariada com regressao

#Dw0escmstr(theta_st,y_st,X_st) #teste

Dw0assmstr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(r-1)])
  v=as.numeric(P[r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-(v+p)/(v+d[i])
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%h%*%t(h)%*%invB%*%e[i,]-as.numeric(t[i]-Ai)*t(X[[i]])%*%invB%*%h
    dAidaw=t(h)%*%invB%*%e[i,]
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
    }
    dqfdai=as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=as.numeric(t[i]-Ai)*invB%*%e[i,]-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi),as.matrix(0))
    can=rep(0,n)
    can[i]=1
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao da assimetria na t assimetrica multivariada com regressao

#Dw0assmstr(theta_st,y_st,X_st) #teste

#Medidas de influencia local 

#w1=rep(0,n)
#w1[150]=-3
#yw_st<-pertres(y_st,w1,c(1,0)) 
#yw_st #perturbacao nas respostas
#w2=rep(0,n)
#w2[150]=6
#Xw_st<-pertexp(X_st,w2,c(1,3)) 
#Xw_st #perturbacao nas explicativas

#hist(yw_st[,1])
#plot(yw_st[,1]~Xw_st[1,2,])
#identify(Xw_st[1,2,],yw_st[,1],label=1:n)
#plot(yw_st[,1]~Xw_st[1,3,])
#identify(Xw_st[1,3,],yw_st[,1],label=1:n)
#hist(yw_st[,2])
#plot(yw_st[,2]~Xw_st[2,2,])
#identify(Xw_st[2,2,],yw_st[,2],label=1:n)
#plot(yw_st[,2]~Xw_st[2,3,])
#identify(Xw_st[2,3,],yw_st[,2],label=1:n)
#hist(Xw_st[1,2,])
#hist(Xw_st[2,2,])
#hist(Xw_st[1,3,])
#hist(Xw_st[2,3,]) #testes 

#thetaw_st<-EMDMSTRM(yw_st,Xw_st)$estmax
#thetaw_st
#theta_st<-EMDMSTRM(y_st,X_st)$estmax
#theta_st #comparacao das estimativas dos parametros com e sem perturbacao

#Dw0mstr=Dw0ponmstr(thetaw_st,yw_st,Xw_st) #ponderacao de casos
#Dw0mstr=Dw0resmstr(thetaw_st,yw_st,Xw_st,c(1,0)) #perturbacao na resposta escolhida
#Dw0mstr=Dw0expmstr(thetaw_st,yw_st,Xw_st,c(1,3)) #perturbacao na explicativa escolhida
#Dw0mstr=Dw0escmstr(thetaw_st,yw_st,Xw_st) #perturbacao na escala
#Dw0mstr=Dw0assmstr(thetaw_st,yw_st,Xw_st) #perturbacao na assimetria

#HQw0st=-2*t(Dw0mstr)%*%solve(hesqmstr(thetaw_st,yw_st,Xw_st))%*%Dw0mstr
#plot(eigen(HQw0st)$values) #teste de direcoes dadas pela matriz de curvatura

#avlwst=rep(0,n)
#for(k in 1:r){
#  avlwst[k]=eigen(HQw0st)$values[k]/sum(eigen(HQw0st)$values[1:r])
#}
#print(avlwst[1:r]) #autovalores positivos ponderados

#avtwst=eigen(HQw0st)$vectors[,1:r] 
#plot(avtwst[,1]) #autovetores associados

#M0st<-rep(0,n)
#for(j in 1:n){
# M0stj=0   
# for(k in 1:r){  
#  M0stj=M0stj+avlwst[k]*avtwst[j,k]^2  
# }
# M0st[j]=M0stj  
#}  
#mean(M0st)

#plot(M0st)
#for(c in 1:4){
#  abline(h=mean(M0st)+c*sd(M0st),col=c)
#}
#identify(1:n,M0st) #analise baseada em poon e poon as vezes nao detecta

#M0st<-rep(0,n)
#for(j in 1:n){
#  M0st[j]=avlwst[1]*avtwst[j,1]^2
#}  
#mean(M0st)

#plot(M0st)
#for(c in 1:4){
#  abline(h=mean(M0st)+c*sd(M0st),col=c)
#}
#identify(1:n,M0st) #analise baseada em cook em geral nao detecta


#Modelo skew slash

Dw0ponmssr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(r-1)])
  v=as.numeric(P[r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  lu<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-((2*v+p)/d[i])*(pgamma(1,v+1+p/2,scale=2/d[i])/pgamma(1,v+p/2,scale=2/d[i]))
    lu[i]<-integrate(fauxss,0,1,v,p,d[i],1,1)$val/pgamma(1,v+p/2,scale=2/d[i])*(d[i]/2)^(v+p/2)/gamma(v+p/2)
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%e[i,]-as.numeric(t[i])*t(X[[i]])%*%invB%*%h
    ddida<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    dldda<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dldda[j]<-2*sum(diag(invB%*%Bj))
      dAida[j]<--t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      ddida[j]<--t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]
    }
    dqfdai=-1/2*dldda-1/2*as.numeric(u[i])*ddida+as.numeric(t[i]-Ai)*dAida
    dqfdhi=as.numeric(t[i]-Ai)*invB%*%e[i,]
    dqfdvi=1/v+lu[i]
    dqfdti=c(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi),as.matrix(dqfdvi))
    can=rep(0,n)
    can[i]=1
    Dw0=Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta ponderacao de casos na slash assimetrica multivariada com regressao

#Dw0ponmssr(theta_ss,y_ss,X_ss) #teste

Dw0resmssr<-function(P,y,X,D)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(r-1)])
  v=as.numeric(P[r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  Sy=sqrt(diag(cov(y)))
  dAidaw=t(h)%*%invB%*%diag(D)%*%Sy
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-((2*v+p)/d[i])*(pgamma(1,v+1+p/2,scale=2/d[i])/pgamma(1,v+p/2,scale=2/d[i]))
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%diag(D)%*%Sy
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=-t(h)%*%invB%*%Bj%*%invB%*%diag(D)%*%Sy
      ddidaw[j]=-2*t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%diag(D)%*%Sy
    }
    dqfdai=-1/2*as.numeric(u[i])*ddidaw+as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=as.numeric(t[i]-Ai)*invB%*%diag(D)%*%Sy-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi),as.matrix(0))
    can=rep(0,n)
    can[i]=1
    Dw0=Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao resposta na slash assimetrica multivariada com regressao

#Dw0resmssr(theta_ss,y_ss,X_ss,c(0,1)) #teste

Dw0expmssr<-function(P,y,X,V)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(r-1)])
  v=as.numeric(P[r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  mX<-X[[1]][1,2:(q/2)]
  for(i in 2:n){
    mX<-rbind(mX,X[[i]][1,2:(q/2)])
  }
  Mx=diag(c(0,sqrt(diag(cov(mX)))),q)
  E=matrix(0,p,q)
  E[V[1],(V[1]-1)*(q/2)+V[2]+1]<-1
  dAidaw=-t(h)%*%invB%*%E%*%Mx%*%b
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-((2*v+p)/d[i])*(pgamma(1,v+1+p/2,scale=2/d[i])/pgamma(1,v+p/2,scale=2/d[i]))
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=as.numeric(u[i])*(t(E%*%Mx)%*%invB%*%invB%*%e[i,]-t(X[[i]])%*%invB%*%invB%*%E%*%Mx%*%b)-as.numeric(t[i,]-Ai)*t(E%*%Mx)%*%invB%*%h-t(X[[i]])%*%invB%*%h%*%t(h)%*%invB%*%E%*%Mx%*%b
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=t(h)%*%invB%*%Bj%*%invB%*%E%*%Mx%*%b
      ddidaw[j]=2*t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%E%*%Mx%*%b
    }
    dqfdai=-1/2*as.numeric(u[i])*ddidaw+as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=-as.numeric(t[i]-Ai)*invB%*%E%*%Mx%*%b-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi),as.matrix(0))
    can=rep(0,n)
    can[i]=1
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao nas explicativas na slash assimetrica multivariada com regressao

#Dw0expmssr(theta_ss,y_ss,X_ss,c(2,3)) #teste

Dw0escmssr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(r-1)])
  v=as.numeric(P[r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-((2*v+p)/d[i])*(pgamma(1,v+1+p/2,scale=2/d[i])/pgamma(1,v+p/2,scale=2/d[i]))
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%e[i,]-1/2*as.numeric(t[i]-Ai)*t(X[[i]])%*%invB%*%h
    dAidaw=1/2*t(h)%*%invB%*%e[i,]
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=-1/2*t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      ddidaw[j]=-t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]
    }
    dqfdai=-1/2*as.numeric(u[i])*ddidaw+as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=1/2*as.numeric(t[i]-Ai)*invB%*%e[i,]-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi),as.matrix(0))
    can=rep(0,n)
    can[i]=1
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao da escala na slash assimetrica multivariada com regressao

#Dw0escmssr(theta_ss,y_ss,X_ss) #teste

Dw0assmssr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(r-1)])
  v=as.numeric(P[r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-((2*v+p)/d[i])*(pgamma(1,v+1+p/2,scale=2/d[i])/pgamma(1,v+p/2,scale=2/d[i]))
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%h%*%t(h)%*%invB%*%e[i,]-as.numeric(t[i]-Ai)*t(X[[i]])%*%invB%*%h
    dAidaw=t(h)%*%invB%*%e[i,]
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      ddidaw[j]=-t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]
    }
    dqfdai=as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=as.numeric(t[i]-Ai)*invB%*%e[i,]-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi),as.matrix(0))
    can=rep(0,n)
    can[i]=1
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao da assimetria na slash assimetrica multivariada com regressao

#Dw0assmssr(theta_ss,y_ss,X_ss) #teste

#Medidas de influencia local 

#w1=rep(0,n)
#w1[150]=4
#yw_ss<-pertres(y_ss,w1,c(0,1)) 
#yw_ss #perturbacao nas respostas
#w2=rep(0,n)
#w2[150]=-3
#Xw_ss<-pertexp(X_ss,w2,c(2,3)) 
#Xw_ss #perturbacao nas explicativas

#hist(yw_ss[,1])
#plot(yw_ss[,1]~Xw_ss[1,2,])
#identify(Xw_ss[1,2,],yw_ss[,1],label=1:n)
#plot(yw_ss[,1]~Xw_ss[1,3,])
#identify(Xw_ss[1,3,],yw_ss[,1],label=1:n)
#hist(yw_ss[,2])
#plot(yw_ss[,2]~Xw_ss[2,2,])
#identify(Xw_ss[2,2,],yw_ss[,2],label=1:n)
#plot(yw_ss[,2]~Xw_ss[2,3,])
#identify(Xw_ss[2,3,],yw_ss[,2],label=1:n)
#hist(Xw_ss[1,2,])
#hist(Xw_ss[2,2,])
#hist(Xw_ss[1,3,])
#hist(Xw_ss[2,3,]) #testes 

#thetaw_ss<-EMDMSSRM(yw_ss,Xw_ss)$estmax
#thetaw_ss
#theta_ss<-EMDMSSRM(y_ss,X_ss)$estmax
#theta_ss #comparacao das estimativas dos parametros com e sem perturbacao

#Dw0mssr=Dw0ponmssr(thetaw_ss,yw_ss,Xw_ss) #ponderacao de casos
#Dw0mssr=Dw0resmssr(thetaw_ss,yw_ss,Xw_ss,c(0,1)) #perturbacao na resposta escolhida 
#Dw0mssr=Dw0expmssr(thetaw_ss,yw_ss,Xw_ss,c(2,3)) #perturbacao na explicativa escolhida 
#Dw0mssr=Dw0escmssr(thetaw_ss,yw_ss,Xw_ss) #perturbacao na escala
#Dw0mssr=Dw0assmssr(thetaw_ss,yw_ss,Xw_ss) #perturbacao na assimetria

#HQw0ss=-2*t(Dw0mssr)%*%solve(hesqmssr(thetaw_ss,yw_ss,Xw_ss))%*%Dw0mssr
#plot(eigen(HQw0ss)$values) #teste de direcoes dadas pela matriz de curvatura

#avlwss=rep(0,n)
#for(k in 1:r){
#  avlwss[k]=eigen(HQw0ss)$values[k]/sum(eigen(HQw0ss)$values[1:r])
#}
#print(avlwss[1:r]) #autovalores positivos ponderados

#avtwss=eigen(HQw0ss)$vectors[,1:r] 
#plot(avtwss[,1]) #autovetores associados

#M0ss<-rep(0,n)
#for(j in 1:n){
#  M0ssj=0   
#  for(k in 1:r){  
#    M0ssj=M0ssj+avlwss[k]*avtwss[j,k]^2  
#  }
#  M0ss[j]=M0ssj  
#}  
#mean(M0ss)

#plot(M0ss)
#for(c in 1:4){
#  abline(h=mean(M0ss)+c*sd(M0ss),col=c)
#}
#identify(1:n,M0ss) #analise baseada em poon e poon as vezes nao detecta

#M0ss<-rep(0,n)
#for(j in 1:n){
#  M0ss[j]=avlwss[1]*avtwss[j,1]^2
#}  
#mean(M0ss)

#plot(M0ss)
#for(c in 1:4){
#  abline(h=mean(M0ss)+c*sd(M0ss),col=c)
#}
#identify(1:n,M0ss) #analise baseada em cook em geral nao detecta


#Modelo skew contaminado

Dw0ponmscr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(r-2)])
  v=as.numeric(P[r-1])
  g=as.numeric(P[r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  v1<-matrix(0,n,1)
  v2<-matrix(0,n,1)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    v1[i]<-(v*g^(p/2)*exp(-g*d[i]/2))/(v*g^(p/2)*exp(-g*d[i]/2)+(1-v)*exp(-d[i]/2))
    v2[i]<-1-v1[i]
    u[i]<-g*v1[i]+v2[i]
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi<-t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%e[i,]-as.numeric(t[i])*t(X[[i]])%*%invB%*%h
    ddida<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    dldda<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dldda[j]<-2*sum(diag(invB%*%Bj))
      dAida[j]<--t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      ddida[j]<--t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]
    }
    dqfdai=-1/2*dldda-1/2*as.numeric(u[i])*ddida+as.numeric(t[i]-Ai)*dAida
    dqfdhi=as.numeric(t[i]-Ai)*invB%*%e[i,]
    dqfdvi=v1[i]/v-v2[i]/(1-v)
    dqfdgi=1/(2*g)*v1[i]*(p-g*d[i])
    dqfdti=c(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi),as.matrix(dqfdvi),as.matrix(dqfdgi))
    can=rep(0,n)
    can[i]=1
    Dw0=Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta ponderacao de casos na contaminada assimetrica multivariada com regressao

#Dw0mscr(theta_sc,y_sc,X_sc) #teste

Dw0resmscr<-function(P,y,X,D)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(r-2)])
  v=as.numeric(P[r-1])
  g=as.numeric(P[r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  v1<-matrix(0,n,1)
  v2<-matrix(0,n,1)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  Sy=sqrt(diag(cov(y)))
  dAidaw=t(h)%*%invB%*%diag(D)%*%Sy
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    v1[i]<-(v*g^(p/2)*exp(-g*d[i]/2))/(v*g^(p/2)*exp(-g*d[i]/2)+(1-v)*exp(-d[i]/2))
    v2[i]<-1-v1[i]
    u[i]<-g*v1[i]+v2[i]
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%diag(D)%*%Sy
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=-t(h)%*%invB%*%Bj%*%invB%*%diag(D)%*%Sy
      ddidaw[j]=-2*t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%diag(D)%*%Sy
    }
    dqfdai=-1/2*as.numeric(u[i])*ddidaw+as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=as.numeric(t[i]-Ai)*invB%*%Sy-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi),as.matrix(c(0,0)))
    can=rep(0,n)
    can[i]=1
    Dw0=Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao da resposta na contaminada assimetrica multivariada com regressao

#Dw0resmscr(theta_sc,y_sc,X_sc,c(1,1)) #teste

Dw0expmscr<-function(P,y,X,V)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(r-2)])
  v=as.numeric(P[r-1])
  g=as.numeric(P[r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  v1<-matrix(0,n,1)
  v2<-matrix(0,n,1)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  mX<-X[[1]][1,2:(q/2)]
  for(i in 2:n){
    mX<-rbind(mX,X[[i]][1,2:(q/2)])
  }
  Mx=diag(c(0,sqrt(diag(cov(mX)))),q)
  E=matrix(0,p,q)
  E[E[1],(E[1]-1)*(q/2)+E[2]+1]<-1
  dAidaw=-t(h)%*%invB%*%E%*%Mx%*%b
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    v1[i]<-(v*g^(p/2)*exp(-g*d[i]/2))/(v*g^(p/2)*exp(-g*d[i]/2)+(1-v)*exp(-d[i]/2))
    v2[i]<-1-v1[i]
    u[i]<-g*v1[i]+v2[i]
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=as.numeric(u[i])*(t(E%*%Mx)%*%invB%*%invB%*%e[i,]-t(X[[i]])%*%invB%*%invB%*%E%*%Mx%*%b)-as.numeric(t[i,]-Ai)*t(E%*%Mx)%*%invB%*%h-t(X[[i]])%*%invB%*%h%*%t(h)%*%invB%*%E%*%Mx%*%b
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=t(h)%*%invB%*%Bj%*%invB%*%E%*%Mx%*%b
      ddidaw[j]=2*t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%E%*%Mx%*%b
    }
    dqfdai=-1/2*as.numeric(u[i])*ddidaw+as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=-as.numeric(t[i]-Ai)*invB%*%E%*%Mx%*%b-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi),as.matrix(c(0,0)))
    can=rep(0,n)
    can[i]=1
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao nas explicativas na contaminada assimetrica multivariada com regressao

#Dw0expmscr(theta_sc,y_sc,X_sc,c(1,3)) #teste

Dw0escmscr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(r-2)])
  v=as.numeric(P[r-1])
  g=as.numeric(P[r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  v1<-matrix(0,n,1)
  v2<-matrix(0,n,1)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    v1[i]<-(v*g^(p/2)*exp(-g*d[i]/2))/(v*g^(p/2)*exp(-g*d[i]/2)+(1-v)*exp(-d[i]/2))
    v2[i]<-1-v1[i]
    u[i]<-g*v1[i]+v2[i]
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%e[i,]-1/2*as.numeric(t[i]-Ai)*t(X[[i]])%*%invB%*%h
    dAidaw=1/2*t(h)%*%invB%*%e[i,]
    ddidaw<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=-1/2*t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      ddidaw[j]=-t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]
    }
    dqfdai=-1/2*as.numeric(u[i])*ddidaw+as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=1/2*as.numeric(t[i]-Ai)*invB%*%e[i,]-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi),as.matrix(c(0,0)))
    can=rep(0,n)
    can[i]=1
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao da escala na contaminada assimetrica multivariada com regressao

#Dw0escmscr(theta_sc,y_sc,X_sc) #teste

Dw0assmscr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  r=length(P)
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(r-2)])
  v=as.numeric(P[r-1])
  g=as.numeric(P[r])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  v1<-matrix(0,n,1)
  v2<-matrix(0,n,1)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  Dw0<-matrix(0,r,n)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    v1[i]<-(v*g^(p/2)*exp(-g*d[i]/2))/(v*g^(p/2)*exp(-g*d[i]/2)+(1-v)*exp(-d[i]/2))
    v2[i]<-1-v1[i]
    u[i]<-g*v1[i]+v2[i]
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdbi=t(X[[i]])%*%invB%*%h%*%t(h)%*%invB%*%e[i,]-as.numeric(t[i]-Ai)*t(X[[i]])%*%invB%*%h
    dAidaw=t(h)%*%invB%*%e[i,]
    dAida<-matrix(0,p0,1)
    d2Aidaw<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidaw[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
    }
    dqfdai=as.numeric(t[i]-Ai)*d2Aidaw-as.numeric(dAidaw)*dAida
    dqfdhi=as.numeric(t[i]-Ai)*invB%*%e[i,]-as.numeric(dAidaw)*invB%*%e[i,]
    dqfdti=rbind(as.matrix(dqfdbi),as.matrix(dqfdai),as.matrix(dqfdhi),as.matrix(c(0,0)))
    can=rep(0,n)
    can[i]=1
    Dw0<-Dw0+dqfdti%*%t(can)
  }
  return(Dw0)
} #delta perturbacao da assimetria na contaminada assimetrica multivariada com regressao

#Dw0assmscr(theta_sc,y_sc,X_sc) #teste

#Medidas de influencia local 

#w1=rep(0,n)
#w1[c(100,200)]=c(-3,2)
#yw_sc<-pertres(y_sc,w1,c(1,1)) 
#yw_sc #perturbacao nas respostas
#w2=rep(0,n)
#w2[150]=-6
#Xw_sc<-pertexp(X_sc,w2,c(2,2)) 
#Xw_sc #perturbacao nas explicativas

#hist(yw_sc[,1])
#plot(yw_sc[,1]~Xw_sc[1,2,])
#identify(Xw_sc[1,2,],yw_sc[,1],label=1:n)
#plot(yw_sc[,1]~Xw_sc[1,3,])
#identify(Xw_sc[1,3,],yw_sc[,1],label=1:n)
#hist(yw_sc[,2])
#plot(yw_sc[,2]~Xw_sc[2,2,])
#identify(Xw_sc[2,2,],yw_sc[,2],label=1:n)
#plot(yw_sc[,2]~Xw_sc[2,3,])
#identify(Xw_sc[2,3,],yw_sc[,2],label=1:n)
#hist(Xw_sc[1,2,])
#hist(Xw_sc[2,2,])
#hist(Xw_sc[1,3,])
#hist(Xw_sc[2,3,]) #testes 

#thetaw_sc<-EMDMSCRM(yw_sc,Xw_sc)$estmax
#thetaw_sc
#theta_sc<-EMDMSCRM(y_sc,X_sc)$estmax
#theta_sc #comparacao das estimativas dos parametros com e sem perturbacao

#Dw0mscr=Dw0ponmscr(thetaw_sc,yw_sc,Xw_sc) #ponderacao de casos
#Dw0mscr=Dw0resmscr(thetaw_sc,yw_sc,Xw_sc,c(1,1)) #perturbacao na respostas escolhida
#Dw0mscr=Dw0expmscr(thetaw_sc,yw_sc,Xw_sc,c(2,2)) #perturbacao na explicativa escolhida 
#Dw0mscr=Dw0escmscr(thetaw_sc,yw_sc,Xw_sc) #perturbacao na escala
#Dw0mscr=Dw0assmscr(thetaw_sc,yw_sc,Xw_sc) #perturbacao na assimetria

#HQw0sc=-2*t(Dw0mscr)%*%solve(hesqmscr(thetaw_sc,yw_sc,Xw_sc))%*%Dw0mscr
#plot(eigen(HQw0sc)$values) #teste de direcoes dadas pela matriz de curvatura

#avlwsc=rep(0,n)
#for(k in 1:r){
# avlwsc[k]=eigen(HQw0sc)$values[k]/sum(eigen(HQw0sc)$values[1:r])
#}
#print(avlwsc[1:r]) #autovalores positivos ponderados

#avtwsc=eigen(HQw0sc)$vectors[,1:r] 
#plot(avtwsc[,1]) #autovetores associados

#M0sc<-rep(0,n)
#for(j in 1:n){
#  M0scj=0   
#  for(k in 1:r){  
#    M0scj=M0scj+avlwsc[k]*avtwsc[j,k]^2  
#  }
#  M0sc[j]=M0scj  
#}  
#mean(M0sc)

#plot(M0sc)
#for(c in 1:4){
#  abline(h=mean(M0sc)+c*sd(M0sc),col=c)
#}
#identify(1:n,M0sc) #analise baseada em poon e poon as vezes nao detecta

#M0sc<-rep(0,n)
#for(j in 1:n){
#  M0sc[j]=avlwsc[1]*avtwsc[j,1]^2
#}  
#mean(M0sc)

#plot(M0sc)
#for(c in 1:4){
#  abline(h=mean(M0sc)+c*sd(M0sc),col=c)
#}
#identify(1:n,M0sc) #analise baseada em cook em geral nao detecta
