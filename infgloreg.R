#Analise de influencia global para os modelos multivariados com regressao

#Pacotes e variaveis globais

#library(MASS)
#library(expm)
#library(moments)
#library(abind)
#library(mvtnorm)
#require(MCMCpack)
#require(doParallel)

#source(file="C:\\Users\\GRACILIANO\\Documents\\Topicos\\stmreg.R")
#source(file="C:\\Users\\GRACILIANO\\Documents\\Topicos\\ssmreg.R")
#source(file="C:\\Users\\GRACILIANO\\Documents\\Topicos\\scmreg.R")

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
#X<-run1(n)

#Modelo skew t

#y_st<-rmstr(n,beta,rSigma,lambda,nu,X_st)
#colMeans(y_st)

#Q function 

qmstr<-function(P,PH,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p*(p+1)/2)])
  invB=solve(B)
  h=as.matrix(P[(q+p*(p+1)/2+1):(length(P)-1)])
  v=as.numeric(P[length(P)])
  bH=as.matrix(PH[1:q])
  BH=xpnd(PH[(q+1):(q+p*(p+1)/2)])
  invBH=solve(BH)
  hH=as.matrix(PH[(q+p*(p+1)/2+1):(length(PH)-1)])
  vH=as.numeric(PH[length(PH)])
  e<-matrix(0,n,p)
  d<-matrix(0,n,1)
  eH<-matrix(0,n,p)
  dH<-matrix(0,n,1)
  uH<-matrix(0,n,1)
  luH<-matrix(0,n,1)
  tH<-matrix(0,n,1)
  qf1<-0
  qf2<-0
  for(i in 1:n){
    eH[i,]<-y[i,]-X[[i]]%*%bH
    dH[i]<-t(eH[i,])%*%invBH%*%invBH%*%eH[i,]
    uH[i]<-(vH+p)/(vH+dH[i])
    luH[i]<-digamma((vH+p)/2)-log((vH+dH[i])/2)
    auxHi<-pmax(t(hH)%*%invBH%*%eH[i,],-37)
    tH[i]<-auxHi+dnorm(auxHi)/pnorm(auxHi)
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    auxi<-pmax(t(h)%*%invB%*%e[i,],-37)
    qf1<-qf1-1/2*as.numeric(log(det(B)^2)+uH[i]*d[i]-2*tH[i]*auxi+auxi^2)
    qf2<-qf2+as.numeric(v/2*(log(v/2)-uH[i]+luH[i])-log(gamma(v/2)))
  }
  qf=qf1+qf2
  return(qf)
} #funcao q da t normal assimetrica multivariada com regressao

#theta_st=EMDMSTRM(y_st,X_st)$estmax
#print(theta_st)
#theta0_st=EMDMSTRP(y_st,X_st)$estmax
#print(theta0_st)
#qmstr(theta_st,theta_st,y_st,X_st)
#qmstr(theta0_st,theta0_st,y_st,X_st)

#passoe_st<-function(P,y,X)
#{
#  n=length(X)/(nrow(X)*ncol(X))
#  p=nrow(X)
#  if(missing(X)) {X<-array(1,c(p,1,n))}
#  q=ncol(X)
#  p0=p*(p+1)/2
#  b=as.matrix(P[1:q])
#  B=xpnd(P[(q+1):(q+p0)])
#  invB=solve(B)
#  h=as.matrix(P[(q+p0+1):(length(P)-1)])
#  v=as.numeric(P[length(P)])
#  e<-matrix(0,n,p)
#  d<-matrix(0,n,1)
#  u<-matrix(0,n,1)
#  lu<-matrix(0,n,1)
#  t<-matrix(0,n,1)
#  for(i in 1:n){
#    e[i,]<-y[i,]-X[[i]]%*%b
#    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
#    u[i]<-(v+p)/(v+d[i])
#    lu[i]<-digamma((v+p)/2)-log((v+d[i])/2)
#    auxi<-pmax(t(h)%*%invB%*%e[i,],-37)
#    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
#  }  
#  return(list(u=u,lu=lu,t=t))
#}

#obj_st=EMDMSTRM(y_st,X_st)
#theta_st=obj_st$estmax
#ec_st=passoe_st(theta_st,y_st,X_st)
#ec_st$u-obj_st$u
#ec_st$lu-obj_st$lu
#ec_st$t-obj_st$t

#qamstr<-function(P,u,lu,t,y,X)
#{
#  n=nrow(y)
#  p=ncol(y)
#  if(missing(X)) {X<-array(1,c(p,1,n))}
#  q=ncol(X)
#  b=as.matrix(P[1:q])
#  B=xpnd(P[(q+1):(q+p*(p+1)/2)])
#  invB=solve(B)
#  h=as.matrix(P[(q+p*(p+1)/2+1):(length(P)-1)])
#  v=as.numeric(P[length(P)])
#  e<-matrix(0,n,p)
#  d<-matrix(0,n,1)
#  qf1<-0
#  qf2<-0
#  for(i in 1:n){
#    e[i,]<-y[i,]-X[[i]]%*%b
#    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
#    auxi<-pmax(t(h)%*%invB%*%e[i,],-37)
#    qf1<-qf1-1/2*as.numeric(log(det(B)^2)+u[i]*d[i]-2*t[i]*auxi+auxi^2)
#    qf2<-qf2+as.numeric(v/2*(log(v/2)-u[i]+lu[i])-log(gamma(v/2)))
#  }
#  qf=qf1+qf2
#  return(qf)
#} #funcao q alternativa da t normal assimetrica multivariada com regressao

#qamstr(theta_st,ec_st$u,ec_st$lu,ec_st$t,y_st,X_st)
#qmstr(theta_st,theta_st,y_st,X_st) 
#qamstr(theta_st,obj_st$u,obj_st$lu,obj_st$t,y_st,X_st)
#qmstr(theta0_st,theta0_st,y_st,X_st)
#lmstr(theta_st,y_st,X_st)
#lmstr(theta0_st,y_st,X_st)
#lmstvr(theta_st[9],theta_st[1:3],xpnd(theta_st[4:6]),theta_st[7:8],y_st,X_st)
#lmstvr(theta0_st[9],theta0_st[1:3],xpnd(theta0_st[4:6]),theta0_st[7:8],y_st,X_st) #testes

gradqmstr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(length(P)-1)])
  v=as.numeric(P[length(P)])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  lu<-matrix(0,n,1)
  t<-matrix(0,n,1)
  dqfdb1<-matrix(0,q,1)
  dqfdb2<-matrix(0,q,1)
  dqfda<-matrix(0,p0,1)
  dqfdh<-matrix(0,p,1)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-(v+p)/(v+d[i])
    lu[i]<-digamma((v+p)/2)-log((v+d[i])/2)
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdb1<-dqfdb1+t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%e[i,]
    dqfdb2<-dqfdb2+as.numeric(t[i])*t(X[[i]])%*%invB%*%h
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
    dqfda<-dqfda-1/2*dldda-1/2*as.numeric(u[i])*ddida+as.numeric(t[i]-Ai)*dAida
    dqfdh<-dqfdh+as.numeric(t[i]-Ai)*invB%*%e[i,]
  }
  dqfdb=dqfdb1-dqfdb2
  dqfdv=n/2*(log(v/2)-digamma(v/2)+1-1/n*sum(u-lu))
  gqf=c(as.vector(dqfdb),as.vector(dqfda),as.vector(dqfdh),dqfdv)
  return(gqf)
} #gradiente da funcao q da t normal assimetrica multivariada com regressao

#gradqmstr(theta_st,y_st,X_st)
#gradqmstr(theta0_st,y_st,X_st)
#escmstr(theta_st,y_st,X_st)
#escmstr(theta0_st,y_st,X_st) #teste

hesqmstr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(length(P)-1)])
  v=as.numeric(P[length(P)])
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
    u[i]<-(v+p)/(v+d[i])
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
  d2qfdv2=n/4*(2/v-trigamma(v/2))
  d2qfdb=cbind(d2qfdb2,d2qfdbda,d2qfdbdh,matrix(0,q,1))
  d2qfda=cbind(t(d2qfdbda),d2qfda2,t(d2qfdhda),matrix(0,p0,1))
  d2qfdh=cbind(t(d2qfdbdh),d2qfdhda,d2qfdh2,matrix(0,p,1))
  d2qfdv=cbind(matrix(0,1,q+p0+p),matrix(d2qfdv2))
  Hqf=rbind(d2qfdb,d2qfda,d2qfdh,d2qfdv)
  return(Hqf)
} #hessiana da funcao q da t assimetrica multivariada com regressao

#hesqmstr(theta_st,y_st,X_st)
#eigen(hesqmstr(theta_st,y_st,X_st))$values
#mifmstr(theta_st,y_st,X_st)
#eigen(mifmstr(theta_st,y_st,X_st))$values
#hesqmstr(theta0_st,y_st,X_st)
#eigen(hesqmstr(theta0_st,y_st,X_st))$values
#mifmstr(theta0_st,y_st,X_st)
#eigen(mifmstr(theta0_st,y_st,X_st))$values #testes

#Deleção para tamanho de amostra n=300

#hist(y_st[,1],breaks=50)
#hist(y_st[,2],breaks=50)
#which.max(y_st[,1])
#which.min(y_st[,2])
#which.max(y_st[,2])
#plot(y_st[,1]~y_st[,2])
#identify(y_st[,1],y_st[,2],label=1:n)

#n=300
#X<-run1(n)
#y_st<-rmstr(n,beta,rSigma,lambda,nu,X)
#y_st
#theta_st<-EMDMSTRM(y_st,X)$estmax
#theta_st

#cl<-makeCluster(detectCores()-1)
#registerDoParallel(cl)
#Theta_st<-NULL
#Theta_st<-foreach(j=1:n,.combine=rbind)%dopar%{
#  library(MASS)
#  library(expm)
#  library(moments)
#  library(abind)
#  library(mvtnorm)
#  require(MCMCpack)
#  Theta_st<-EMDMSTRM(y_st[-j,],X_st[,,-j])$estmax
  #nTheta_st<-c(EMDMSTRM(y_st[-j,],X_st[,,-j])$estmax,j)
#}
#print(Theta_st)
#stopCluster(cl)

#Theta_st<-nTheta_st[,1:9]

#lTheta_st<-matrix(0,q+p0+p+1,n)
#qfTheta_st<-matrix(0,q+p0+p+1,n)
#for(i in 1:n){
  #lTheta_st[,i]<-theta_st-solve(mifmstr(theta_st,y_st,X))%*%escmstr(theta_st,y_st[-i,],X[,,-i])
#  qfTheta_st[,i]<-theta_st-solve(hesqmstr(theta_st,y_st,X_st))%*%gradqmstr(theta_st,y_st[-i,],X_st[,,-i])
#}
#lTheta_st

#theta_st
#apply(Theta_st,2,summary)
#apply(qfTheta_st,1,summary)
#which.min(Theta_st[-c(291,184,81),8])
#which.max(Theta_st[-c(6,184,291,166,27,177),7])
#EMDMSTRM(y_st[-6,],X_st[,,-6])$estmax
#EMDMSTRM(y_st[-214,],X_st[,,-214])$estmax

#Medidas classicas de influencia global

#QD0_st<-rep(0,n)
#QD_st<-rep(0,n)
#GD0_st<-rep(0,n)
#GD_st<-rep(0,n)
#for(i in 1:n){
#  QD0_st[i]<-2*(qmstr(theta_st,theta_st,y_st,X_st)-qmstr(Theta_st[i,],theta_st,y_st,X_st))
#  QD_st[i]<-2*(qmstr(theta_st,theta_st,y_st,X_st)-qmstr(qfTheta_st[,i],theta_st,y_st,X_st))
#  GD0_st[i]<--t(theta_st-Theta_st[i,])%*%hesqmstr(theta_st,y_st,X_st)%*%(theta_st-Theta_st[i,])
#  GD_st[i]<--t(gradqmstr(theta_st,y_st[-i,],X_st[,,-i]))%*%solve(hesqmstr(theta_st,y_st,X_st))%*%gradqmstr(theta_st,y_st[-i,],X_st[,,-i])
#}
#par(mfrow=c(1,2))
#plot(QD0_st~QD_st)
#plot(GD0_st~GD_st)
#plot(QD0_st~GD0_st)

#plot(QD_st)
#identify(1:n,QD_st)
#plot(GD_st)
#identify(1:n,GD_st)

#Ponto de corte:

#c=1,2,3,4 

#plot(QD0_st)
#for(c in 1:4){
#abline(h=mean(QD0_st)+c*sd(QD0_st),col=c)
#}
#identify(1:n,QD0_st)

#plot(GD0_st)
#for(c in 1:4){
#  abline(h=mean(GD0_st)+c*sd(GD0_st),col=c)
#}
#identify(1:n,GD0_st)

#plot(QD_st)
#for(c in 1:4){
#  abline(h=mean(QD_st)+c*sd(QD_st),col=c)
#}
#identify(1:n,QD_st)

#plot(GD_st)
#for(c in 1:4){
#  abline(h=mean(GD_st)+c*sd(GD_st),col=c)
#}
#identify(1:n,GD_st)
#GD_st[291]

#Distancia de Cook generalizada por grupo de parametros

#beta-alpha

#GDba0_st<-rep(0,n)
#GDba_st<-rep(0,n)
#for(i in 1:n){
#  GDba0_st[i]<--t(theta_st-Theta_st[i,])[1:(q+p0)]%*%hesqmstr(theta_st,y_st,X_st)[1:(q+p0),1:(q+p0)]%*%(theta_st-Theta_st[i,])[1:(q+p0)]
#  GDba_st[i]<--t(gradqmstr(theta_st,y_st[-i,],X_st[,,-i]))[1:(q+p0)]%*%solve(hesqmstr(theta_st,y_st,X_st))[1:(q+p0),1:(q+p0)]%*%gradqmstr(theta_st,y_st[-i,],X_st[,,-i])[1:(q+p0)]
#}
#which.max(GDba0_st)
#plot(GDba0_st~GDba_st)
#plot(GDba0_st)
#for(c in 1:4){
#  abline(h=mean(GDba0_st)+c*sd(GDba0_st),col=c)
#}
#identify(1:n,GDba0_st)

#lambda

#GDh0_st<-rep(0,n)
#GDh_st<-rep(0,n)
#for(i in 1:n){
#  GDh0_st[i]<--t(theta_st-Theta_st[i,])[(q+p0+1):(q+p0+p)]%*%hesqmstr(theta_st,y_st,X_st)[(q+p0+1):(q+p0+p),(q+p0+1):(q+p0+p)]%*%(theta_st-Theta_st[i,])[(q+p0+1):(q+p0+p)]
#  GDh_st[i]<--t(gradqmstr(theta_st,y_st[-i,],X_st[,,-i]))[(q+p0+1):(q+p0+p)]%*%solve(hesqmstr(theta_st,y_st,X_st))[(q+p0+1):(q+p0+p),(q+p0+1):(q+p0+p)]%*%gradqmstr(theta_st,y_st[-i,],X_st[,,-i])[(q+p0+1):(q+p0+p)]
#}
#which.max(GDh0_st)
#plot(GDh0_st~GDh_st)
#plot(GDh0_st)
#for(c in 1:4){
#  abline(h=mean(GDh0_st)+c*sd(GDh0_st),col=c)
#}
#identify(1:n,GDh_st)
#plot(GDh_st)
#for(c in 1:4){
#  abline(h=mean(GDh_st)+c*sd(GDh_st),col=c)
#}
#identify(1:n,GDh_st)

#tau

#GDt0_st<-rep(0,n)
#GDt_st<-rep(0,n)
#for(i in 1:n){
#  GDt0_st[i]<--t(theta_st-Theta_st[i,])[q+p0+p+1]%*%hesqmstr(theta_st,y_st,X_st)[q+p0+p+1,q+p0+p+1]%*%(theta_st-Theta_st[i,])[q+p0+p+1]
#  GDt_st[i]<--t(gradqmstr(theta_st,y_st[-i,],X_st[,,-i]))[q+p0+p+1]%*%solve(hesqmstr(theta_st,y_st,X_st))[q+p0+p+1,q+p0+p+1]%*%gradqmstr(theta_st,y_st[-i,],X_st[,,-i])[q+p0+p+1]
#}
#which.max(GDt0_st)
#plot(GDt0_st~GDt_st)
#plot(GDt0_st)
#for(c in 1:4){
#  abline(h=mean(GDt0_st)+c*sd(GDt0_st),col=c)
#}
#identify(1:n,GDt0_st)
#plot(GDt_st)
#plot(GDt_st)
#for(c in 1:4){
#  abline(h=mean(GDt_st)+c*sd(GDt_st),col=c)
#}
#identify(1:n,GDt_st)


#Mahalanobis x pesos

mahalanobis_st<-function(P,y,X)
{
  n=length(X)
  p=nrow(X[[n]])
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(length(P)-1)])
  v=as.numeric(P[length(P)])
  e<-matrix(0,n,p)
  d<-matrix(0,n,1)
  u<-matrix(0,n,1)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-(v+p)/(v+d[i])
  }  
  return(list(d=d,u=u))
}
#plot(mahalanobis_st(theta_st,y_st,X_st)$d,mahalanobis_st(theta_st,y_st,X_st)$u)
#identify(mahalanobis_st(theta_st,y_st,X_st)$d,mahalanobis_st(theta_st,y_st,X_st)$u,labels=1:n)
#mahalanobis_st(theta_st,y_st,X_st)$d[291]


#TRC e MRC

#par(mfrow=c(1,2))

#TRC_st<-rep(0,n)
#for(i in 1:n){
#  TRC_st[i]=sum(abs((theta_st-Theta_st[i,])/theta_st))
#}
#plot(TRC_st)
#for(c in 1:4){
#  abline(h=mean(TRC_st)+c*sd(TRC_st),col=c)
#}
#identify(1:n,TRC_st)

#MRC_st<-rep(0,n)
#for(i in 1:n){
#  MRC_st[i]=max(abs((theta_st-Theta_st[i,])/theta_st))
#}
#plot(MRC_st)
#for(c in 1:4){
#  abline(h=mean(MRC_st)+c*sd(MRC_st),col=c)
#}
#identify(1:n,MRC_st)

#Com valor esperado

#f0st<-function(u,v,h){
#  (v/2)^(v/2)/gamma(v/2)*(u^((v-3)/2)*exp(-u*v/2))/sqrt(u+as.numeric(t(h)%*%h))
#}

#f0st(0.5,nu,lambda)
#integrate(f0st,lower=0,upper=Inf,nu,lambda)$val #teste

#Ehesqmstr<-function(P,y,X)
#{
#  n=length(X)/(nrow(X)*ncol(X))
#  p=nrow(X)
#  if(missing(X)) {X<-array(1,c(p,1,n))}
#  q=ncol(X)
#  p0=p*(p+1)/2
#  b=as.matrix(P[1:q])
#  B=xpnd(P[(q+1):(q+p0)])
#  invB=solve(B)
#  h=as.matrix(P[(q+p0+1):(length(P)-1)])
#  v=as.numeric(P[length(P)])
#  Ee=sqrt(2/pi)*integrate(f0st,lower=0,upper=Inf,v,h)$val*B%*%h #esperanca residuos
#  EA=t(h)%*%invB%*%Ee
#  EeeT=(v/2)^(v/2)/gamma(v/2)*(2/v)^(v/2-1)*gamma(v/2-1)*B%*%B #segundo momento residuos
#  VZ=(v/2)^(v/2)/gamma(v/2)*(2/v)^(v/2-1)*gamma(v/2-1)*diag(p)-2/pi*(integrate(f0st,lower=0,upper=Inf,v,h)$val)^2*h%*%t(h) #variancia residuos padrao  
#  e<-matrix(0,n,p)
#  d<-matrix(0,n,1)
#  u<-matrix(0,n,1)
#  t<-matrix(0,n,1)
#  Ed2qfdb2<-matrix(0,q,q)
#  Ed2qfdbda<-matrix(0,q,p0)
#  Ed2qfdbdh<-matrix(0,q,p)
#  Ed2qfda2<-matrix(0,p0,p0)
#  Ed2qfdhda<-matrix(0,p,p0)
#  Ed2qfdh2<-matrix(0,p,p)
#  for(i in 1:n){
#    e[i,]<-y[i,]-X[[i]]%*%b
#    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
#    u[i]<-(v+p)/(v+d[i])
#    auxi<-pmax(t(h)%*%invB%*%e[i,],-37)
#    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    #Ai<-as.numeric(t(h)%*%invB%*%e[i,])
#    Ed2qfdb2<-Ed2qfdb2-t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%X[[i]]
#    Ed2qfdbdh<-Ed2qfdbdh-as.numeric(t[i]-EA)*t(X[[i]])%*%invB+t(X[[i]])%*%invB%*%h%*%t(Ee)%*%invB
#    Ed2qfdh2<-Ed2qfdh2-invB%*%EeeT%*%invB
#    EdAida<-matrix(0,p0,1)
#    Ed2Aidhda<-matrix(0,p,p0) 
#    Ed2didbda<-matrix(0,q,p0)
#    Ed2Aidbda<-matrix(0,q,p0)
#    Ed2ldda2<-matrix(0,p0,p0)
#    Ed2dida2<-matrix(0,p0,p0)
#    Ed2Aida2<-matrix(0,p0,p0)
#    EAd2Aidhda<-matrix(0,p,p0)
#    EdAidhdAidaT<-matrix(0,p,p0)
#    EAd2Aida2<-matrix(0,p0,p0) 
#    EdAidadAidaT<-matrix(0,p0,p0)
#    for(j in 1:p0){
#      ind=rep(0,p0)
#      ind[j]=1
#      Bj=xpnd(ind)
#      Ed2Aidhda[,j]=-invB%*%Bj%*%invB%*%Ee
#      EdAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%Ee
#      Ed2Aidbda[,j]=t(X[[i]])%*%invB%*%Bj%*%invB%*%h
#      Ed2didbda[,j]=2*t(X[[i]])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%Ee
#      EAd2Aidhda[,j]=-invB%*%Bj%*%invB%*%EeeT%*%invB%*%h 
#      EdAidhdAidaT[,j]=-invB%*%EeeT%*%invB%*%Bj%*%invB%*%h
#      for(k in 1:p0){
#        ind=rep(0,p0)
#        ind[k]=1
#        Bk=xpnd(ind)
#        Ed2ldda2[j,k]=-sum(diag(invB%*%Bk%*%invB%*%Bj))
#        Ed2Aida2[j,k]=-t(h)%*%invB%*%(Bk%*%invB%*%Bj+Bj%*%invB%*%Bk)%*%invB%*%Ee
#        Ed2dida2[j,k]=sum(diag((Bk%*%invB%*%Bj%*%invB+Bj%*%invB%*%Bk%*%invB+Bj%*%invB%*%invB%*%Bk+Bk%*%invB%*%invB%*%Bj+invB%*%Bk%*%invB%*%Bj+invB%*%Bj%*%invB%*%Bk)%*%VZ))+t(Ee)%*%invB%*%(Bk%*%invB%*%Bj%*%invB+Bj%*%invB%*%Bk%*%invB+Bj%*%invB%*%invB%*%Bk+Bk%*%invB%*%invB%*%Bj+invB%*%Bk%*%invB%*%Bj+invB%*%Bj%*%invB%*%Bk)%*%invB%*%Ee
#        EAd2Aida2[j,k]=t(h)%*%invB%*%(Bk%*%invB%*%Bj+Bj%*%invB%*%Bk)%*%invB%*%EeeT%*%invB%*%h
#        EdAidadAidaT[j,k]=t(h)%*%invB%*%Bj%*%invB%*%EeeT%*%invB%*%Bk%*%invB%*%h
#      }
#    }
#    Ed2qfdhda<-Ed2qfdhda+as.numeric(t[i])*Ed2Aidhda-EAd2Aidhda-EdAidhdAidaT 
#    Ed2qfdbda<-Ed2qfdbda-1/2*as.numeric(u[i])*Ed2didbda+as.numeric(t[i]-EA)*Ed2Aidbda+t(X[[i]])%*%invB%*%h%*%t(EdAida) 
#    Ed2qfda2<-Ed2qfda2-1/2*Ed2ldda2-1/2*as.numeric(u[i])*Ed2dida2+as.numeric(t[i])*Ed2Aida2-EAd2Aida2-EdAidadAidaT 
#  }
#  Ed2qfdv2=n/4*(2/v-trigamma(v/2))
#  Ed2qfdb=cbind(Ed2qfdb2,Ed2qfdbda,Ed2qfdbdh,matrix(0,q,1))
#  Ed2qfda=cbind(t(Ed2qfdbda),Ed2qfda2,t(Ed2qfdhda),matrix(0,p0,1))
#  Ed2qfdh=cbind(t(Ed2qfdbdh),Ed2qfdhda,Ed2qfdh2,matrix(0,p,1))
#  Ed2qfdv=cbind(matrix(0,1,q+p0+p),matrix(Ed2qfdv2))
#  EHqf=rbind(Ed2qfdb,Ed2qfda,Ed2qfdh,Ed2qfdv)
#  return(EHqf)
#} #Esperança da hessiana da funcao q da t normal assimetrica multivariada com regressao

#Ehesqmstr(theta_st,y_st,X_st)
#eigen(Ehesqmstr(theta_st,y_st,X_st))$values #teste

#EqfTheta_st<-matrix(0,q+p0+p+1,n)
#for(i in 1:n){
#  EqfTheta_st[,i]<-theta_st-solve(Ehesqmstr(theta_st,y_st,X_st))%*%gradqmstr(theta_st,y_st[-i,],X_st[,,-i])
#}
#colMeans(Theta_st-t(EqfTheta_st))
#apply(Theta_st-t(EqfTheta_st),2,min)


#Modelo skew slash

#y_ss<-rmssr(n,beta,rSigma,lambda,nu,X_ss)
#colMeans(y_ss)

#Q function 

qmssr<-function(P,PH,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p*(p+1)/2)])
  invB=solve(B)
  h=as.matrix(P[(q+p*(p+1)/2+1):(length(P)-1)])
  v=as.numeric(P[length(P)])
  bH=as.matrix(PH[1:q])
  BH=xpnd(PH[(q+1):(q+p*(p+1)/2)])
  invBH=solve(BH)
  hH=as.matrix(PH[(q+p*(p+1)/2+1):(length(PH)-1)])
  vH=as.numeric(PH[length(PH)])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  eH<-matrix(0,n,p)
  dH<-matrix(0,n,1)
  uH<-matrix(0,n,1)
  luH<-matrix(0,n,1)
  tH<-matrix(0,n,1)
  qf1<-0
  qf2<-0
  for(i in 1:n){
    eH[i,]<-y[i,]-X[[i]]%*%bH
    dH[i]<-t(eH[i,])%*%invBH%*%invBH%*%eH[i,]
    uH[i]<-(2*vH+p)/dH[i]*pgamma(1,vH+1+p/2,scale=2/dH[i])/pgamma(1,vH+p/2,scale=2/dH[i])
    luH[i]<-integrate(fauxss,0,1,vH,p,dH[i],1,1)$val/pgamma(1,vH+p/2,scale=2/dH[i])*(dH[i]/2)^(vH+p/2)/gamma(vH+p/2)
    auxHi<-pmax(t(hH)%*%invBH%*%eH[i,],-37)
    tH[i]<-auxHi+dnorm(auxHi)/pnorm(auxHi)
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    auxi<-pmax(t(h)%*%invB%*%e[i,],-37)
    qf1<-qf1-1/2*as.numeric(log(det(B)^2)+uH[i]*d[i]-2*tH[i]*auxi+auxi^2)
    qf2<-qf2+as.numeric(log(v)+(v-1)*luH[i])
  }
  qf=qf1+qf2
  return(qf)
} #funcao q da slash normal assimetrica multivariada

#theta_ss=EMDMSSRM(y_ss,X)$estmax
#print(theta_ss)
#qmssr(theta_ss,theta_ss,y_ss,X_ss)

gradqmssr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(length(P)-1)])
  v=as.numeric(P[length(P)])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  lu<-matrix(0,n,1)
  t<-matrix(0,n,1)
  dqfdb<-matrix(0,q,1)
  dqfdb1<-matrix(0,q,1)
  dqfdb2<-matrix(0,q,1)
  dqfda<-matrix(0,p0,1)
  dqfdh<-matrix(0,p,1)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-((2*v+p)/d[i])*(pgamma(1,v+1+p/2,scale=2/d[i])/pgamma(1,v+p/2,scale=2/d[i]))
    lu[i]<-integrate(fauxss,0,1,v,p,d[i],1,1)$val/pgamma(1,v+p/2,scale=2/d[i])*(d[i]/2)^(v+p/2)/gamma(v+p/2)
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdb1<-dqfdb1+t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%e[i,]
    dqfdb2<-dqfdb2+as.numeric(t[i])*t(X[[i]])%*%invB%*%h
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
    dqfda<-dqfda-1/2*dldda-1/2*as.numeric(u[i])*ddida+as.numeric(t[i]-Ai)*dAida
    dqfdh<-dqfdh+as.numeric(t[i]-Ai)*invB%*%e[i,]
  }
  dqfdb=dqfdb1-dqfdb2
  dqfdv=n/v+sum(lu)
  gqf=c(as.vector(dqfdb),as.vector(dqfda),as.vector(dqfdh),dqfdv)
  return(gqf)
} #gradiente da funcao q da slash normal assimetrica multivariada com regressao

#gradqmssr(theta_ss,y_ss,X_ss)
#escmssr(theta_ss,y_ss,X_ss) #teste

hesqmssr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(length(P)-1)])
  v=as.numeric(P[length(P)])
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
    u[i]<-(2*v+p)/d[i]*pgamma(1,v+1+p/2,scale=2/d[i])/pgamma(1,v+p/2,scale=2/d[i])
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
  d2qfdv2=-n/v^2
  d2qfdb=cbind(d2qfdb2,d2qfdbda,d2qfdbdh,matrix(0,q,1))
  d2qfda=cbind(t(d2qfdbda),d2qfda2,t(d2qfdhda),matrix(0,p0,1))
  d2qfdh=cbind(t(d2qfdbdh),d2qfdhda,d2qfdh2,matrix(0,p,1))
  d2qfdv=cbind(matrix(0,1,q+p0+p),matrix(d2qfdv2))
  Hqf=rbind(d2qfdb,d2qfda,d2qfdh,d2qfdv)
  return(Hqf)
} #hessiana da funcao q da slash normal assimetrica multivariada com regressao

#hesqmssr(theta_ss,y_ss,X_ss)
#eigen(hesqmssr(theta_ss,y_ss,X_ss))$values
#mifmssr(theta_ss,y_ss,X_ss)
#eigen(mifmssr(theta_ss,y_ss,X_ss))$values #teste

#Delecao para tamanho de amostra n=300

#hist(y_ss[,1],breaks=50)
#hist(y_ss[,2],breaks=50)
#plot(y_ss[,1]~y_ss[,2])
#which.max(y_ss[,1])
#which.min(y_ss[,2])
#which.min(y_ss[,1])
#which.max(y_ss[,2])

#n=300
#X<-run1(n)
#y_ss<-rmssr(n,beta,rSigma,lambda,nu,X_ss)
#y_ss
#theta_ss<-EMDMSSRM(y_ss,X_ss)$estmax
#theta_ss

#cl<-makeCluster(detectCores()-1)
#registerDoParallel(cl)
#Theta_ss<-NULL
#Theta_ss<-foreach(j=1:n,.combine=rbind)%dopar%{
#  library(MASS)
#  library(expm)
#  library(moments)
#  library(abind)
#  library(mvtnorm)
#  require(MCMCpack)
#  Theta_ss<-EMDMSSRM(y_ss[-j,],X_ss[,,-j])$estmax
#}
#print(Theta_ss)
#stopCluster(cl)

#lTheta_ss<-matrix(0,q+p0+p+1,n)
#qfTheta_ss<-matrix(0,q+p0+p+1,n)
#for(i in 1:n){
  #lTheta_ss[,i]<-theta_ss-solve(mifmssr(theta_ss,y_ss,X_ss))%*%escmssr(theta_ss,y_ss[-i,],X_ss[,,-i])
#  qfTheta_ss[,i]<-theta_ss-solve(hesqmssr(theta_ss,y_ss,X_ss))%*%gradqmssr(theta_ss,y_ss[-i,],X_ss[,,-i])
#}
#lTheta_ss

#theta_ss
#apply(Theta_ss,2,summary)
#apply(qfTheta_ss,1,summary)
#plot(Theta_ss[,9]~qfTheta_ss[9,])
#which.min(Theta_ss[,8])
#Theta_ss[214,]
#qfTheta_ss[,214]
#EMDMSSRM(y_ss[-214,],X_ss[,,-214])$estmax
#which.max(Theta_ss[-214,4])

#Medidas classicas de influencia global

#QD0_ss<-rep(0,n)
#QD_ss<-rep(0,n)
#GD0_ss<-rep(0,n)
#GD_ss<-rep(0,n)
#for(i in 1:n){
#  QD0_ss[i]<-2*(qmssr(theta_ss,theta_ss,y_ss,X_ss)-qmssr(Theta_ss[i,],theta_ss,y_ss,X_ss))
#  QD_ss[i]<-2*(qmssr(theta_ss,theta_ss,y_ss,X_ss)-qmssr(qfTheta_ss[,i],theta_ss,y_ss,X_ss))
#  GD0_ss[i]<--t(theta_ss-Theta_ss[i,])%*%hesqmssr(theta_ss,y_ss,X_ss)%*%(theta_ss-Theta_ss[i,])
#  GD_ss[i]<--t(gradqmssr(theta_ss,y_ss[-i,],X_ss[,,-i]))%*%solve(hesqmssr(theta_ss,y_ss,X_ss))%*%gradqmssr(theta_ss,y_ss[-i,],X_ss[,,-i])
#}
#par(mfrow=c(1,2))
#plot(QD0_ss~QD_ss)
#plot(GD0_ss~GD_ss)

#plot(QD_ss)
#identify(1:n,QD_ss)
#plot(QD0_ss)
#identify(1:n,QD0_ss)
#plot(GD_ss)
#identify(1:n,GD_ss)
#plot(GD0_ss)
#identify(1:n,GD0_ss)

#Ponto de corte:

#c=1,2,3,4 

#plot(QD0_ss)
#for(c in 1:4){
#  abline(h=mean(QD0_ss)+c*sd(QD0_ss),col=c)
#}
#identify(1:n,QD0_ss)

#plot(GD0_ss)
#for(c in 1:4){
#  abline(h=mean(GD0_ss)+c*sd(GD0_ss),col=c)
#}
#identify(1:n,GD0_ss)

#plot(QD_ss)
#for(c in 1:4){
#  abline(h=mean(QD_ss)+c*sd(QD_ss),col=c)
#}
#identify(1:n,QD_ss)

#plot(GD_ss)
#for(c in 1:4){
#  abline(h=mean(GD_ss)+c*sd(GD_ss),col=c)
#}
#identify(1:n,GD_ss)

#Distancia de Cook generalizada por grupo de parametros

#beta-alpha

#GDba0_ss<-rep(0,n)
#GDba_ss<-rep(0,n)
#for(i in 1:n){
#  GDba0_ss[i]<--t(theta_ss-Theta_ss[i,])[1:(q+p0)]%*%hesqmssr(theta_ss,y_ss,X_ss)[1:(q+p0),1:(q+p0)]%*%(theta_ss-Theta_ss[i,])[1:(q+p0)]
#  GDba_ss[i]<--t(gradqmssr(theta_ss,y_ss[-i,],X_ss[,,-i]))[1:(q+p0)]%*%solve(hesqmssr(theta_ss,y_ss,X_ss))[1:(q+p0),1:(q+p0)]%*%gradqmssr(theta_ss,y_ss[-i,],X_ss[,,-i])[1:(q+p0)]
#}
#which.max(GDba0_ss)
#plot(GDba0_ss~GDba_ss)
#plot(GDba0_ss)
#for(c in 1:4){
#  abline(h=mean(GDba0_ss)+c*sd(GDba0_ss),col=c)
#}
#identify(1:n,GDba0_ss)

#lambda

#GDh0_ss<-rep(0,n)
#GDh_ss<-rep(0,n)
#for(i in 1:n){
#  GDh0_ss[i]<--t(theta_ss-Theta_ss[i,])[(q+p0+1):(q+p0+p)]%*%hesqmssr(theta_ss,y_ss,X_ss)[(q+p0+1):(q+p0+p),(q+p0+1):(q+p0+p)]%*%(theta_ss-Theta_ss[i,])[(q+p0+1):(q+p0+p)]
#  GDh_ss[i]<--t(gradqmssr(theta_ss,y_ss[-i,],X_ss[,,-i]))[(q+p0+1):(q+p0+p)]%*%solve(hesqmssr(theta_ss,y_ss,X_ss))[(q+p0+1):(q+p0+p),(q+p0+1):(q+p0+p)]%*%gradqmssr(theta_ss,y_ss[-i,],X_ss[,,-i])[(q+p0+1):(q+p0+p)]
#}
#which.max(GDh0_ss[-214])
#plot(GDh0_ss~GDh_ss)
#plot(GDh0_ss)
#for(c in 1:4){
#  abline(h=mean(GDh0_ss)+c*sd(GDh0_ss),col=c)
#}
#identify(1:n,GDh0_ss)
#which.max(Theta_ss[,1])
#theta_ss

#tau

#GDt0_ss<-rep(0,n)
#GDt_ss<-rep(0,n)
#for(i in 1:n){
#  GDt0_ss[i]<--t(theta_ss-Theta_ss[i,])[q+p0+p+1]%*%hesqmssr(theta_ss,y_ss,X_ss)[q+p0+p+1,q+p0+p+1]%*%(theta_ss-Theta_ss[i,])[q+p0+p+1]
#  GDt_ss[i]<--t(gradqmssr(theta_ss,y_ss[-i,],X_ss[,,-i]))[q+p0+p+1]%*%solve(hesqmssr(theta_ss,y_ss,X_ss))[q+p0+p+1,q+p0+p+1]%*%gradqmssr(theta_ss,y_ss[-i,],X_ss[,,-i])[q+p0+p+1]
#}
#which.max(GDt0_ss[-214])
#plot(GDt0_ss~GDt_ss)
#plot(GDt0_ss)
#for(c in 1:4){
#  abline(h=mean(GDt0_ss)+c*sd(GDt0_ss),col=c)
#}
#identify(1:n,GDt0_ss)

#Mahalanobis x pesos

mahalanobis_ss<-function(P,y,X)
{
  n=length(X)
  p=nrow(X[[n]])
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(length(P)-1)])
  v=as.numeric(P[length(P)])
  e<-matrix(0,n,p)
  d<-matrix(0,n,1)
  u<-matrix(0,n,1)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-(2*v+p)/d[i]*pgamma(1,v+1+p/2,scale=2/d[i])/pgamma(1,v+p/2,scale=2/d[i])
  }  
  return(list(d=d,u=u))
}
#plot(mahalanobis_ss(theta_ss,y_ss,X_ss)$d,mahalanobis_ss(theta_ss,y_ss,X_ss)$u)
#identify(mahalanobis_ss(theta_ss,y_ss,X_ss)$d,mahalanobis_ss(theta_ss,y_ss,X_ss)$u,labels=1:n)
#which.max(mahalanobis_ss(theta_ss,y_ss,X_ss)$d)
#mahalanobis_ss(theta_ss,y_ss,X_ss)$d[165]

#TRC e MRC

#par(mfrow=c(1,2))

#TRC_ss<-rep(0,n)
#for(i in 1:n){
#  TRC_ss[i]=sum(abs((theta_ss-Theta_ss[i,])/theta_ss))
#}
#plot(TRC_ss)
#for(c in 1:4){
#  abline(h=mean(TRC_ss)+c*sd(TRC_ss),col=c)
#}
#identify(1:n,TRC_ss)

#MRC_ss<-rep(0,n)
#for(i in 1:n){
#  MRC_ss[i]=max(abs((theta_ss-Theta_ss[i,])/theta_ss))
#}
#plot(MRC_ss)
#for(c in 1:4){
#  abline(h=mean(MRC_ss)+c*sd(MRC_ss),col=c)
#}
#identify(1:n,MRC_ss)

#Com valor esperado

#f0ss<-function(u,v,h){
#  v*(u^(v-3/2)/sqrt(u+as.numeric(t(h)%*%h)))
#}

#f0ss(0.5,nu,lambda)
#integrate(f0ss,lower=0,upper=1,nu,lambda)$val #teste

#Ehesqmssr<-function(P,y,X)
#{
#  n=length(X)/(nrow(X)*ncol(X))
#  p=nrow(X)
#  if(missing(X)) {X<-array(1,c(p,1,n))}
#  q=ncol(X)
#  p0=p*(p+1)/2
#  b=as.matrix(P[1:q])
#  B=xpnd(P[(q+1):(q+p0)])
#  invB=solve(B)
#  h=as.matrix(P[(q+p0+1):(length(P)-1)])
#  v=as.numeric(P[length(P)])
#  Ee=sqrt(2/pi)*integrate(f0ss,lower=0,upper=1,v,h)$val*B%*%h #esperanca residuos
#  EA=t(h)%*%invB%*%Ee
#  EeeT=v/(v-1)*B%*%B #segundo momento residuos
#  VZ=v/(v-1)*diag(p)-2/pi*(integrate(f0ss,lower=0,upper=1,v,h)$val)^2*h%*%t(h) #variancia residuos padrao  
#  e<-matrix(0,n,p)
#  d<-matrix(0,n,1)
#  u<-matrix(0,n,1)
#  t<-matrix(0,n,1)
#  Ed2qfdb2<-matrix(0,q,q)
#  Ed2qfdbda<-matrix(0,q,p0)
#  Ed2qfdbdh<-matrix(0,q,p)
#  Ed2qfda2<-matrix(0,p0,p0)
#  Ed2qfdhda<-matrix(0,p,p0)
#  Ed2qfdh2<-matrix(0,p,p)
#  for(i in 1:n){
#    e[i,]<-y[i,]-X[[i]]%*%b
#    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
#    u[i]<-(2*v+p)/d[i]*pgamma(1,v+1+p/2,scale=2/d[i])/pgamma(1,v+p/2,scale=2/d[i])
#    auxi<-pmax(t(h)%*%invB%*%e[i,],-37)
#    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    #Ai<-as.numeric(t(h)%*%invB%*%e[i,])
#    Ed2qfdb2<-Ed2qfdb2-t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%X[[i]]
#    Ed2qfdbdh<-Ed2qfdbdh-as.numeric(t[i]-EA)*t(X[[i]])%*%invB+t(X[[i]])%*%invB%*%h%*%t(Ee)%*%invB
#    Ed2qfdh2<-Ed2qfdh2-invB%*%EeeT%*%invB
#    EdAida<-matrix(0,p0,1)
#    Ed2Aidhda<-matrix(0,p,p0) 
#    Ed2didbda<-matrix(0,q,p0)
#    Ed2Aidbda<-matrix(0,q,p0)
#    Ed2ldda2<-matrix(0,p0,p0)
#    Ed2dida2<-matrix(0,p0,p0)
#    Ed2Aida2<-matrix(0,p0,p0)
#    EAd2Aidhda<-matrix(0,p,p0)
#    EdAidhdAidaT<-matrix(0,p,p0)
#    EAd2Aida2<-matrix(0,p0,p0) 
#    EdAidadAidaT<-matrix(0,p0,p0)
#    for(j in 1:p0){
#      ind=rep(0,p0)
#      ind[j]=1
#      Bj=xpnd(ind)
#      Ed2Aidhda[,j]=-invB%*%Bj%*%invB%*%Ee
#      EdAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%Ee
#      Ed2Aidbda[,j]=t(X[[i]])%*%invB%*%Bj%*%invB%*%h
#      Ed2didbda[,j]=2*t(X[[i]])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%Ee
#      EAd2Aidhda[,j]=-invB%*%Bj%*%invB%*%EeeT%*%invB%*%h 
#      EdAidhdAidaT[,j]=-invB%*%EeeT%*%invB%*%Bj%*%invB%*%h
#      for(k in 1:p0){
#        ind=rep(0,p0)
#        ind[k]=1
#        Bk=xpnd(ind)
#        Ed2ldda2[j,k]=-sum(diag(invB%*%Bk%*%invB%*%Bj))
#        Ed2Aida2[j,k]=-t(h)%*%invB%*%(Bk%*%invB%*%Bj+Bj%*%invB%*%Bk)%*%invB%*%Ee
#        Ed2dida2[j,k]=sum(diag((Bk%*%invB%*%Bj%*%invB+Bj%*%invB%*%Bk%*%invB+Bj%*%invB%*%invB%*%Bk+Bk%*%invB%*%invB%*%Bj+invB%*%Bk%*%invB%*%Bj+invB%*%Bj%*%invB%*%Bk)%*%VZ))+t(Ee)%*%invB%*%(Bk%*%invB%*%Bj%*%invB+Bj%*%invB%*%Bk%*%invB+Bj%*%invB%*%invB%*%Bk+Bk%*%invB%*%invB%*%Bj+invB%*%Bk%*%invB%*%Bj+invB%*%Bj%*%invB%*%Bk)%*%invB%*%Ee
#        EAd2Aida2[j,k]=t(h)%*%invB%*%(Bk%*%invB%*%Bj+Bj%*%invB%*%Bk)%*%invB%*%EeeT%*%invB%*%h
#        EdAidadAidaT[j,k]=t(h)%*%invB%*%Bj%*%invB%*%EeeT%*%invB%*%Bk%*%invB%*%h
#      }
#    }
#    Ed2qfdhda<-Ed2qfdhda+as.numeric(t[i])*Ed2Aidhda-EAd2Aidhda-EdAidhdAidaT 
#    Ed2qfdbda<-Ed2qfdbda-1/2*as.numeric(u[i])*Ed2didbda+as.numeric(t[i]-EA)*Ed2Aidbda+t(X[[i]])%*%invB%*%h%*%t(EdAida) 
#    Ed2qfda2<-Ed2qfda2-1/2*Ed2ldda2-1/2*as.numeric(u[i])*Ed2dida2+as.numeric(t[i])*Ed2Aida2-EAd2Aida2-EdAidadAidaT 
#  }
#  Ed2qfdv2=-n/v^2
#  Ed2qfdb=cbind(Ed2qfdb2,Ed2qfdbda,Ed2qfdbdh,matrix(0,q,1))
#  Ed2qfda=cbind(t(Ed2qfdbda),Ed2qfda2,t(Ed2qfdhda),matrix(0,p0,1))
#  Ed2qfdh=cbind(t(Ed2qfdbdh),Ed2qfdhda,Ed2qfdh2,matrix(0,p,1))
#  Ed2qfdv=cbind(matrix(0,1,q+p0+p),matrix(Ed2qfdv2))
#  EHqf=rbind(Ed2qfdb,Ed2qfda,Ed2qfdh,Ed2qfdv)
#  return(EHqf)
#} #Esperança da hessiana da funcao q da slash normal assimetrica multivariada com regressao

#Ehesqmssr(theta_ss,y_ss,X_ss)
#eigen(Ehesqmssr(theta_ss,y_ss,X_ss))$values #teste

#EqfTheta_ss<-matrix(0,q+p0+p+1,n)
#for(i in 1:n){
# EqfTheta_ss[,i]<-theta_ss-solve(Ehesqmssr(theta_ss,y_ss,X_ss))%*%gradqmssr(theta_ss,y_ss[-i,],X_ss[,,-i])
#}
#colMeans(Theta_ss-t(EqfTheta_ss))
#apply(Theta_ss-t(EqfTheta_ss),2,median)


#Modelo skew contaminado

#y_sc<-rmscr(n,beta,rSigma,lambda,eta,gamma,X_sc)
#colMeans(y_sc)

#Q function 

qmscr<-function(P,PH,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p*(p+1)/2)])
  invB=solve(B)
  h=as.matrix(P[(q+p*(p+1)/2+1):(length(P)-2)])
  v=as.numeric(P[length(P)-1])
  g=as.numeric(P[length(P)])
  bH=as.matrix(PH[1:q])
  BH=xpnd(PH[(q+1):(q+p*(p+1)/2)])
  invBH=solve(BH)
  hH=as.matrix(PH[(q+p*(p+1)/2+1):(length(PH)-2)])
  vH=as.numeric(PH[length(PH)-1])
  gH=as.numeric(PH[length(PH)])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  eH<-matrix(0,n,p)
  dH<-matrix(0,n,1)
  v1H<-matrix(0,n,1)
  v2H<-matrix(0,n,1)
  tH<-matrix(0,n,1)
  qf1<-0
  qf2<-0
  for(i in 1:n){
    eH[i,]<-y[i,]-X[[i]]%*%bH
    dH[i]<-t(eH[i,])%*%invBH%*%invBH%*%eH[i,]
    v1H[i]<-(vH*gH^(p/2)*exp(-gH*dH[i]/2))/(vH*gH^(p/2)*exp(-gH*dH[i]/2)+(1-vH)*exp(-dH[i]/2))
    v2H[i]<-1-v1H[i]
    auxHi<-pmax(t(hH)%*%invBH%*%eH[i,],-37)
    tH[i]<-auxHi+dnorm(auxHi)/pnorm(auxHi)
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    auxi<-pmax(t(h)%*%invB%*%e[i,],-37)
    qf1<-qf1-1/2*as.numeric(log(det(B)^2)+v2H[i]*d[i]-2*tH[i]*auxi+auxi^2)
    qf2<-qf2+1/2*as.numeric(v1H[i]*(p*log(g)-g*d[i]+2*log(v))+2*v2H[i]*log(1-v))
  }
  qf=qf1+qf2
  return(qf)
} #funcao q da normal contaminada assimetrica multivariada

#theta_sc=EMDMSCRM(y_sc,X_sc)$estmax
#print(theta_sc)
#qmscr(theta_sc,theta_sc,y_sc,X_sc)

gradqmscr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(length(P)-2)])
  v=as.numeric(P[length(P)-1])
  g=as.numeric(P[length(P)])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  v1<-matrix(0,n,1)
  v2<-matrix(0,n,1)
  u<-matrix(0,n,1)
  t<-matrix(0,n,1)
  dqfdb1<-matrix(0,q,1)
  dqfdb2<-matrix(0,q,1)
  dqfda<-matrix(0,p0,1)
  dqfdh<-matrix(0,p,1)
  dqfdv<-0
  dqfdg<-0
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    v1[i]<-(v*g^(p/2)*exp(-g*d[i]/2))/(v*g^(p/2)*exp(-g*d[i]/2)+(1-v)*exp(-d[i]/2))
    v2[i]<-1-v1[i]
    u[i]<-g*v1[i]+v2[i]
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    Ai<-as.numeric(t(h)%*%invB%*%e[i,])
    dqfdb1<-dqfdb1+t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%e[i,]
    dqfdb2<-dqfdb2+as.numeric(t[i])*t(X[[i]])%*%invB%*%h
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
    dqfda<-dqfda-1/2*dldda-1/2*as.numeric(u[i])*ddida+as.numeric(t[i]-Ai)*dAida
    dqfdh<-dqfdh+as.numeric(t[i]-Ai)*invB%*%e[i,]
    dqfdv<-dqfdv+v1[i]/v-v2[i]/(1-v)
    dqfdg<-dqfdg+1/(2*g)*v1[i]*(p-g*d[i])
  }
  dqfdb=dqfdb1-dqfdb2
  gqf=c(as.vector(dqfdb),as.vector(dqfda),as.vector(dqfdh),dqfdv,dqfdg)
  return(gqf)
} #gradiente da q function da normal contaminada assimetrica multivariada com regressao

#gradqmscr(theta_sc,y_sc,X_sc) 
#escmscr(theta_sc,y_sc,X_sc) #teste

hesqmscr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(length(P)-2)])
  v=as.numeric(P[length(P)-1])
  g=as.numeric(P[length(P)])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  v1<-matrix(0,n,1)
  v2<-matrix(0,n,1)
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
    v1[i]<-(v*g^(p/2)*exp(-g*d[i]/2))/(v*g^(p/2)*exp(-g*d[i]/2)+(1-v)*exp(-d[i]/2))
    v2[i]<-1-v1[i]
    u[i]<-g*v1[i]+v2[i]
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
  d2qfdv2=-(sum(v1)/v^2+sum(v2)/(1-v)^2)
  d2qfdg2=-p*sum(v1)/(2*g^2)
  d2qfdb=cbind(d2qfdb2,d2qfdbda,d2qfdbdh,matrix(0,q,2))
  d2qfda=cbind(t(d2qfdbda),d2qfda2,t(d2qfdhda),matrix(0,p0,2))
  d2qfdh=cbind(t(d2qfdbdh),d2qfdhda,d2qfdh2,matrix(0,p,2))
  d2qfdv=cbind(matrix(0,1,q+p0+p),matrix(d2qfdv2),matrix(0,1,1))
  d2qfdg=cbind(matrix(0,1,q+p0+p),matrix(0,1,1),matrix(d2qfdg2))
  Hqf=rbind(d2qfdb,d2qfda,d2qfdh,d2qfdv,d2qfdg)
  return(Hqf)
} #hessiana da funcao q da normal contaminada assimetrica multivariada com regressao

#hesqmscr(theta_sc,y_sc,X_sc)
#eigen(hesqmscr(theta_sc,y_sc,X_sc))$values 
#mifmscr(theta_sc,y_sc,X_sc)
#eigen(mifmscr(theta_sc,y_sc,X_sc))$values #teste

#Delecao para tamanho de amostra n=300

#hist(y_sc[,1],breaks=50)
#hist(y_sc[,2],breaks=50)
#plot(y_sc[,1]~y_sc[,2])
#which.max(y_sc[,1])
#which.min(y_sc[,2])
#which.min(y_sc[,1])
#which.max(y_sc[,2])

#n=300
#X<-run1(n)
#y_sc<-rmscr(n,beta,rSigma,lambda,eta,gamma,X_sc)
#theta_sc<-EMDMSCRM(y_sc,X_sc)$estmax
#theta_sc

#cl<-makeCluster(detectCores()-1)
#registerDoParallel(cl)
#Theta_sc<-NULL
#Theta_sc<-foreach(j=1:n,.combine=rbind)%dopar%{
#  library(MASS)
#  library(expm)
#  library(moments)
#  library(abind)
#  library(mvtnorm)
#  require(MCMCpack)
#  Theta_sc<-EMDMSCRM(y_sc[-j,],X_sc[,,-j])$estmax
#}
#print(Theta_sc)
#stopCluster(cl)

#lTheta_sc<-matrix(0,q+p0+p+2,n)
#qfTheta_sc<-matrix(0,q+p0+p+2,n)
#for(i in 1:n){
  #lTheta_sc[,i]<-theta_sc-solve(mifmscr(theta_sc,y_sc,X_sc))%*%escmscr(theta_sc,y_sc[-i,],X_sc[,,-i])
#  qfTheta_sc[,i]<-theta_sc-solve(hesqmscr(theta_sc,y_sc,X_sc))%*%gradqmscr(theta_sc,y_sc[-i,],X_sc[,,-i])
#}
#lTheta_sc
#qfTheta_sc

#theta_sc
#apply(Theta_sc,2,summary)
#apply(qfTheta_sc,1,summary)
#plot(Theta_sc[,10]~qfTheta_sc[10,])
#which.max(Theta_sc[-162,7])
#EMDMSCRM(y_sc[-226,],X_sc[,,-226])$estmax

#Medidas classicas de influencia global

#QD0_sc<-rep(0,n)
#QD_sc<-rep(0,n)
#GD0_sc<-rep(0,n)
#GD_sc<-rep(0,n)
#for(i in 1:n){
#  QD0_sc[i]<-2*(qmscr(theta_sc,theta_sc,y_sc,X_sc)-qmscr(Theta_sc[i,],theta_sc,y_sc,X_sc))
#  QD_sc[i]<-2*(qmscr(theta_sc,theta_sc,y_sc,X_sc)-qmscr(qfTheta_sc[,i],theta_sc,y_sc,X_sc))
#  GD0_sc[i]<--t(theta_sc-Theta_sc[i,])%*%hesqmscr(theta_sc,y_sc,X_sc)%*%(theta_sc-Theta_sc[i,])
#  GD_sc[i]<--t(gradqmscr(theta_sc,y_sc[-i,],X_sc[,,-i]))%*%solve(hesqmscr(theta_sc,y_sc,X_sc))%*%gradqmscr(theta_sc,y_sc[-i,],X_sc[,,-i])
#}
#par(mfrow=c(1,2))
#plot(QD0_sc~QD_sc)
#plot(GD0_sc~GD_sc)

#plot(QD_sc)
#identify(1:n,QD_sc)
#plot(QD0_sc)
#identify(1:n,QD0_sc)
#plot(GD_sc)
#identify(1:n,GD_sc)
#plot(GD0_sc)
#identify(1:n,GD0_sc)

#Distancia de Cook generalizada por grupo de parametros

#beta-alpha

#GDba0_sc<-rep(0,n)
#GDba_sc<-rep(0,n)
#for(i in 1:n){
#  GDba0_sc[i]<--t(theta_sc-Theta_sc[i,])[1:(q+p0)]%*%hesqmscr(theta_sc,y_sc,X_sc)[1:(q+p0),1:(q+p0)]%*%(theta_sc-Theta_sc[i,])[1:(q+p0)]
#  GDba_sc[i]<--t(gradqmscr(theta_sc,y_sc[-i,],X_sc[,,-i]))[1:(q+p0)]%*%solve(hesqmscr(theta_sc,y_sc,X_sc))[1:(q+p0),1:(q+p0)]%*%gradqmscr(theta_sc,y_sc[-i,],X_sc[,,-i])[1:(q+p0)]
#}
#which.max(GDba0_sc)
#plot(GDba0_sc~GDba_sc)
#plot(GDba0_sc)
#for(c in 1:4){
#  abline(h=mean(GDba0_sc)+c*sd(GDba0_sc),col=c)
#}
#identify(1:n,GDba0_sc)

#lambda

#GDh0_sc<-rep(0,n)
#GDh_sc<-rep(0,n)
#for(i in 1:n){
#  GDh0_sc[i]<--t(theta_sc-Theta_sc[i,])[(q+p0+1):(q+p0+p)]%*%hesqmscr(theta_sc,y_sc,X_sc)[(q+p0+1):(q+p0+p),(q+p0+1):(q+p0+p)]%*%(theta_sc-Theta_sc[i,])[(q+p0+1):(q+p0+p)]
#  GDh_sc[i]<--t(gradqmscr(theta_sc,y_sc[-i,],X_sc[,,-i]))[(q+p0+1):(q+p0+p)]%*%solve(hesqmscr(theta_sc,y_sc,X_sc))[(q+p0+1):(q+p0+p),(q+p0+1):(q+p0+p)]%*%gradqmscr(theta_sc,y_sc[-i,],X_sc[,,-i])[(q+p0+1):(q+p0+p)]
#}
#which.max(GDh0_sc)
#plot(GDh0_sc~GDh_sc)
#plot(GDh0_sc)
#for(c in 1:4){
#  abline(h=mean(GDh0_sc)+c*sd(GDh0_sc),col=c)
#}
#identify(1:n,GDh0_sc)

#tau

#GDt0_sc<-rep(0,n)
#GDt_sc<-rep(0,n)
#for(i in 1:n){
#  GDt0_sc[i]<--t(theta_sc-Theta_sc[i,])[(q+p0+p+1):(q+p0+p+2)]%*%hesqmscr(theta_sc,y_sc,X_sc)[(q+p0+p+1):(q+p0+p+2),(q+p0+p+1):(q+p0+p+2)]%*%(theta_sc-Theta_sc[i,])[(q+p0+p+1):(q+p0+p+2)]
#  GDt_sc[i]<--t(gradqmscr(theta_sc,y_sc[-i,],X_sc[,,-i]))[(q+p0+p+1):(q+p0+p+2)]%*%solve(hesqmscr(theta_sc,y_sc,X_sc))[(q+p0+p+1):(q+p0+p+2),(q+p0+p+1):(q+p0+p+2)]%*%gradqmscr(theta_sc,y_sc[-i,],X_sc[,,-i])[(q+p0+p+1):(q+p0+p+2)]
#}
#which.max(GDt0_sc)
#plot(GDt0_sc~GDt_sc)
#plot(GDt0_sc)
#for(c in 1:4){
#  abline(h=mean(GDt0_sc)+c*sd(GDt0_sc),col=c)
#}
#identify(1:n,GDt0_sc)

#Mahalanobis x pesos

mahalanobis_sc<-function(P,y,X)
{
  n=length(X)
  p=nrow(X[[n]])
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(length(P)-2)])
  v=as.numeric(P[length(P)-1])
  g=as.numeric(P[length(P)])
  e<-matrix(0,n,p)
  d<-matrix(0,n,1)
  v1<-matrix(0,n,1)
  v2<-matrix(0,n,1)
  u<-matrix(0,n,1)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    v1[i]<-(v*g^(p/2)*exp(-g*d[i]/2))/(v*g^(p/2)*exp(-g*d[i]/2)+(1-v)*exp(-d[i]/2))
    v2[i]<-1-v1[i]
    u[i]<-g*v1[i]+v2[i]
  }  
  return(list(d=d,u=u))
}
#plot(mahalanobis_sc(theta_sc,y_sc,X_sc)$d,mahalanobis_sc(theta_sc,y_sc,X_sc)$u)
#identify(mahalanobis_sc(theta_sc,y_sc,X_sc)$d,mahalanobis_sc(theta_sc,y_sc,X_sc)$u,labels=1:n)

#which.max(mahalanobis_sc(theta_sc,y_sc,X_sc)$d)
#mahalanobis_sc(theta_sc,y_sc,X_sc)$d[226]

#Ponto de corte:

#c=1,2,3,4 

#plot(QD0_sc)
#for(c in 1:4){
#  abline(h=mean(QD0_sc)+c*sd(QD0_sc),col=c)
#}
#identify(1:n,QD0_sc)

#plot(GD0_sc)
#for(c in 1:4){
#  abline(h=mean(GD0_sc)+c*sd(GD0_sc),col=c)
#}
#identify(1:n,GD0_sc)

#plot(QD_sc)
#for(c in 1:4){
#  abline(h=mean(QD_sc)+c*sd(QD_sc),col=c)
#}
#identify(1:n,QD_sc)

#plot(GD_sc)
#for(c in 1:4){
#  abline(h=mean(GD_sc)+c*sd(GD_sc),col=c)
#}
#identify(1:n,GD_sc)

#TRC e MRC

#par(mfrow=c(1,2))

#TRC_sc<-rep(0,n)
#for(i in 1:n){
#  TRC_sc[i]=sum(abs((theta_sc-Theta_sc[i,])/theta_sc))
#}
#plot(TRC_sc)
#for(c in 1:4){
#  abline(h=mean(TRC_sc)+c*sd(TRC_sc),col=c)
#}
#identify(1:n,TRC_sc)

#MRC_sc<-rep(0,n)
#for(i in 1:n){
#  MRC_sc[i]=max(abs((theta_sc-Theta_sc[i,])/theta_sc))
#}
#plot(MRC_sc)
#for(c in 1:4){
#  abline(h=mean(MRC_sc)+c*sd(MRC_sc),col=c)
#}
#identify(1:n,MRC_sc)

#Com valor esperado

#f0sc<-function(u,h){
#  1/sqrt(u*(u+as.numeric(t(h)%*%h)))
#}

#eta*f0sc(gamma,lambda)+(1-eta)*f0sc(1,lambda) #teste

#Ehesqmscr<-function(P,y,X)
#{
#  n=length(X)/(nrow(X)*ncol(X))
#  p=nrow(X)
#  if(missing(X)) {X<-array(1,c(p,1,n))}
#  q=ncol(X)
#  p0=p*(p+1)/2
#  b=as.matrix(P[1:q])
#  B=xpnd(P[(q+1):(q+p0)])
#  invB=solve(B)
#  h=as.matrix(P[(q+p0+1):(length(P)-2)])
#  v=as.numeric(P[length(P)-1])
#  g=as.numeric(P[length(P)])
#  Ee=sqrt(2/pi)*(eta*f0sc(g,h)+(1-eta)*f0sc(1,h))*B%*%h #esperanca residuos
#  EA=t(h)%*%invB%*%Ee
#  EeeT=(v/g+1-v)*B%*%B #segundo momento residuos
#  VZ=(v/g+1-v)*diag(p)-2/pi*(eta*f0sc(g,h)+(1-eta)*f0sc(1,h))^2*h%*%t(h) #variancia residuos padrao  
#  e<-matrix(0,n,p)
#  d<-matrix(0,n,1)
#  v1<-matrix(0,n,1)
#  v2<-matrix(0,n,1)
#  u<-matrix(0,n,1)
#  t<-matrix(0,n,1)
#  Ed2qfdb2<-matrix(0,q,q)
#  Ed2qfdbda<-matrix(0,q,p0)
#  Ed2qfdbdh<-matrix(0,q,p)
#  Ed2qfda2<-matrix(0,p0,p0)
#  Ed2qfdhda<-matrix(0,p,p0)
#  Ed2qfdh2<-matrix(0,p,p)
#  for(i in 1:n){
#    e[i,]<-y[i,]-X[[i]]%*%b
#    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
#    v1[i]<-(v*g^(p/2)*exp(-g*d[i]/2))/(v*g^(p/2)*exp(-g*d[i]/2)+(1-v)*exp(-d[i]/2))
#    v2[i]<-1-v1[i]
#    u[i]<-g*v1[i]+v2[i]
#    auxi<-pmax(t(h)%*%invB%*%e[i,],-37)
#    t[i]<-auxi+dnorm(auxi)/pnorm(auxi)
    #Ai<-as.numeric(t(h)%*%invB%*%e[i,])
#    Ed2qfdb2<-Ed2qfdb2-t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%X[[i]]
#    Ed2qfdbdh<-Ed2qfdbdh-as.numeric(t[i]-EA)*t(X[[i]])%*%invB+t(X[[i]])%*%invB%*%h%*%t(Ee)%*%invB
#    Ed2qfdh2<-Ed2qfdh2-invB%*%EeeT%*%invB
#    EdAida<-matrix(0,p0,1)
#    Ed2Aidhda<-matrix(0,p,p0) 
#    Ed2didbda<-matrix(0,q,p0)
#    Ed2Aidbda<-matrix(0,q,p0)
#    Ed2ldda2<-matrix(0,p0,p0)
#    Ed2dida2<-matrix(0,p0,p0)
#    Ed2Aida2<-matrix(0,p0,p0)
#    EAd2Aidhda<-matrix(0,p,p0)
#    EdAidhdAidaT<-matrix(0,p,p0)
#    EAd2Aida2<-matrix(0,p0,p0) 
#    EdAidadAidaT<-matrix(0,p0,p0)
#    for(j in 1:p0){
#      ind=rep(0,p0)
#      ind[j]=1
#      Bj=xpnd(ind)
#      Ed2Aidhda[,j]=-invB%*%Bj%*%invB%*%Ee
#      EdAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%Ee
#      Ed2Aidbda[,j]=t(X[[i]])%*%invB%*%Bj%*%invB%*%h
#      Ed2didbda[,j]=2*t(X[[i]])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%Ee
#      EAd2Aidhda[,j]=-invB%*%Bj%*%invB%*%EeeT%*%invB%*%h 
#      EdAidhdAidaT[,j]=-invB%*%EeeT%*%invB%*%Bj%*%invB%*%h
#      for(k in 1:p0){
#        ind=rep(0,p0)
#        ind[k]=1
#        Bk=xpnd(ind)
#        Ed2ldda2[j,k]=-sum(diag(invB%*%Bk%*%invB%*%Bj))
#        Ed2Aida2[j,k]=-t(h)%*%invB%*%(Bk%*%invB%*%Bj+Bj%*%invB%*%Bk)%*%invB%*%Ee
#        Ed2dida2[j,k]=sum(diag((Bk%*%invB%*%Bj%*%invB+Bj%*%invB%*%Bk%*%invB+Bj%*%invB%*%invB%*%Bk+Bk%*%invB%*%invB%*%Bj+invB%*%Bk%*%invB%*%Bj+invB%*%Bj%*%invB%*%Bk)%*%VZ))+t(Ee)%*%invB%*%(Bk%*%invB%*%Bj%*%invB+Bj%*%invB%*%Bk%*%invB+Bj%*%invB%*%invB%*%Bk+Bk%*%invB%*%invB%*%Bj+invB%*%Bk%*%invB%*%Bj+invB%*%Bj%*%invB%*%Bk)%*%invB%*%Ee
#        EAd2Aida2[j,k]=t(h)%*%invB%*%(Bk%*%invB%*%Bj+Bj%*%invB%*%Bk)%*%invB%*%EeeT%*%invB%*%h
#        EdAidadAidaT[j,k]=t(h)%*%invB%*%Bj%*%invB%*%EeeT%*%invB%*%Bk%*%invB%*%h
#      }
#    }
#    Ed2qfdhda<-Ed2qfdhda+as.numeric(t[i])*Ed2Aidhda-EAd2Aidhda-EdAidhdAidaT 
#    Ed2qfdbda<-Ed2qfdbda-1/2*as.numeric(u[i])*Ed2didbda+as.numeric(t[i]-EA)*Ed2Aidbda+t(X[[i]])%*%invB%*%h%*%t(EdAida) 
#    Ed2qfda2<-Ed2qfda2-1/2*Ed2ldda2-1/2*as.numeric(u[i])*Ed2dida2+as.numeric(t[i])*Ed2Aida2-EAd2Aida2-EdAidadAidaT 
#  }
#  Ed2qfdv2=-(sum(v1)/v^2+sum(v2)/(1-v)^2)
#  Ed2qfdg2=-p*sum(v1)/(2*g^2)
#  Ed2qfdb=cbind(Ed2qfdb2,Ed2qfdbda,Ed2qfdbdh,matrix(0,q,2))
#  Ed2qfda=cbind(t(Ed2qfdbda),Ed2qfda2,t(Ed2qfdhda),matrix(0,p0,2))
#  Ed2qfdh=cbind(t(Ed2qfdbdh),Ed2qfdhda,Ed2qfdh2,matrix(0,p,2))
#  Ed2qfdv=cbind(matrix(0,1,q+p0+p),matrix(Ed2qfdv2),matrix(0,1,1))
#  Ed2qfdg=cbind(matrix(0,1,q+p0+p),matrix(0,1,1),matrix(Ed2qfdg2))
#  EHqf=rbind(Ed2qfdb,Ed2qfda,Ed2qfdh,Ed2qfdv,Ed2qfdg)
#  return(EHqf)
#} #Esperanca da hessiana da funcao q da normal contaminada assimetrica multivariada com regressao

#Ehesqmscr(theta_sc,y_sc,X_sc)
#eigen(Ehesqmscr(theta_sc,y_sc,X_sc))$values #teste

#EqfTheta_sc<-matrix(0,q+p0+p+2,n)
#for(i in 1:n){
#  EqfTheta_sc[,i]<-theta_sc-solve(Ehesqmscr(theta_sc,y_sc,X))%*%gradqmscr(theta_sc,y_sc[-i,],X[,,-i])
#}
#colMeans(Theta_sc-t(EqfTheta_sc))
#apply(Theta_sc-t(EqfTheta_sc),2,median)
