#Estimacao de maxima verossimilhanca via EM para um modelo slash assimetrico multivariado com regressao


# install.packages("TruncatedDistributions", repos="http://R-Forge.R-project.org") # pacote para rodar trunc dist. Tbm tem o pacote truncd dist

qtgamma <-function (p, shape, scale = 1, a = 0, b = Inf) 
{
    stopifnot(all(p >= 0 & p <= 1) & all(scale > 0) & all(shape > 
        0))
    Fa <- pgamma(a, shape, scale = scale)
    Fb <- pgamma(b, shape, scale = scale)
    pNew <- p * (Fb - Fa) + Fa
    x <- qgamma(pNew, shape, scale = scale)
    return(x)
}

dtgamma <-function (x, shape, scale = 1, a = 0, b = Inf) 
{
    stopifnot(all(shape > 0) & all(scale > 0))
    Fa <- pgamma(a, shape, scale = scale)
    Fb <- pgamma(b, shape, scale = scale)
    y <- dgamma(x, shape, scale = scale)
    inda <- which(x < a)
    indb <- which(x > b)
    if (length(inda) > 0) 
        y[inda] <- 0
    if (length(indb) > 0) 
        y[indb] <- 0
    return(y/(Fb - Fa))
}

#Geracao dos dados

#Rmssr<-function(n,b,B,h,v,X){
#  p=length(h)
#  if(missing(X)) {X<-array(1,c(p,1,n))}
#  y<-matrix(0,n,p)
#  for(i in 1:n){
#    u=rbeta(1,shape1=v,shape2=1) #fator de escala da Slash  
#    y[i,]=X[[i]]%*%as.matrix(b)+B%*%((u*(u+sum(h*h)))^(-1/2)*(h*abs(rnorm(1)))+solve(sqrtm(u*diag(p)+h%*%t(h)))%*%mvrnorm(1,rep(0,p),diag(p)))
#  } 
#  return(y) 
#}

rmssr<-function(n,b,B,h,v,X){
  p=length(h)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  y<-matrix(0,n,p)
  for(i in 1:n){
    u=rbeta(1,shape1=v,shape2=1) #fator de escala da Slash  
    y[i,]=X[[i]]%*%as.matrix(b)+B%*%((u*(u+sum(h*h)))^(-1/2)*(h*abs(rnorm(1)))+solve(sqrtm(u*diag(p)+h%*%t(h)))%*%mvrnorm(1,rep(0,p),diag(p)))
  } 
  return(y) 
}

#y<-rmssr(n,beta,rSigma,lambda,nu,X) 
#y<-rmssr(n,beta[1],rSigma,lambda,nu) #geracao independente de X
#y #dados observados resposta


#Funcoes dos dados incompletos

dmssr<-function(P,y,X)
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
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
  }
  aux=pmax(as.vector(t(h)%*%invB%*%t(e)),-37)
  f=2*v/(det(B)*(2*pi)^(p/2))*gamma(v+p/2)*pgamma(1,v+p/2,scale=2/d)/((d/2)^(v+p/2))*pnorm(aux)
  return(f)
} #densidade da slash normal assimetrica multivariada com regressao


lmssr<-function(P,y,X)
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
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
  }
  aux=pmax(as.vector(t(h)%*%invB%*%t(e)),-37)
  l=sum(log(2*v/(det(B)*(2*pi)^(p/2))*gamma(v+p/2)*pgamma(1,v+p/2,scale=2/d)/((d/2)^(v+p/2))*pnorm(aux)))
  return(l)
} #log verossimilhanca da slash normal assimetrica multivariada com regressao 

#lmssr(P,y,X) #teste
#lmssr(P0,y) #sem explicativas

lmssvr<-function(v,b,B,h,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  invB=solve(B)
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%as.matrix(b)
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
  }
  lv=sum(log(2*v*gamma(v+p/2)*pgamma(1,v+p/2,scale=2/d)/((d/2)^(v+p/2))))
  return(lv)
} #log verossimilhanca da slash normal assimetrica multivariada com regressao em nu

#lmssvr(nu,beta,Sigma,lambda,y,X) #teste
#lmssvr(nu,beta[1],Sigma,lambda,y) #sem explicativas


fauxss<-function(u,v,p,d,r,s){
 u^(v+p/2-r)*exp(-u*d/2)*log(u)^s
}

elogu<-function(u,nu,p,d){
 log(u)*dtgamma(u, shape=nu+p/2, scale=2/d, a=0, b=1)
}


#fauxss(0.5,nu,p,2,1,1)
#integrate(fauxss,lower=0,upper=1,nu,p,2,1,1)$val #teste

escmssr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b=as.matrix(P[1:q])
  p0=p*(p+1)/2
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p0+1):(length(P)-1)])
  v=as.numeric(P[length(P)])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  W<-matrix(0,n,1)
  u<-matrix(0,n,1)
  dldb1<-matrix(0,q,1)
  dldb2<-matrix(0,q,1)
  dlda<-matrix(0,p0,1)
  dldv0<-0
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-(2*v+p)/d[i]*pgamma(1,v+1+p/2,scale=2/d[i])/pgamma(1,v+p/2,scale=2/d[i])
    auxi<-pmax(t(h)%*%invB%*%e[i,],-37)
    W[i]<-dnorm(auxi)/pnorm(auxi)
    dldb1<-dldb1+as.numeric(u[i])*t(X[[i]])%*%invB%*%invB%*%e[i,]
    dldb2<-dldb2+as.numeric(W[i])*t(X[[i]])%*%invB%*%h
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
    dlda<-dlda-1/2*dldda-1/2*as.numeric(u[i])*ddida+as.numeric(W[i])*dAida
    dldv0<-dldv0+integrate(fauxss,lower=0,upper=1,v,p,d[i],1,1)$val/integrate(fauxss,lower=0,upper=1,v,p,d[i],1,0)$val
  }
  dldb=dldb1-dldb2
  dldv=dldv0+n/v
  dldh=invB%*%t(e)%*%W
  el=c(as.vector(dldb),as.vector(dlda),as.vector(dldh),dldv)
  return(el)
} #escore da slash normal assimetrica multivariada com regressao

#escmssr(P,y,X) #teste
#escmssr(P0,y) #sem explicativas

mifmssr<-function(P,y,X)
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
  vu<-matrix(0,n,1)
  W<-matrix(0,n,1)
  Wl<-matrix(0,n,1)
  d2ldb2<-matrix(0,q,q)
  d2ldbda<-matrix(0,q,p0)
  d2ldbdh<-matrix(0,q,p)
  d2lda2<-matrix(0,p0,p0)
  d2ldhda<-matrix(0,p,p0)
  d2ldh2<-matrix(0,p,p)
  d2ldbdv<-matrix(0,q,1)
  d2ldadv<-matrix(0,p0,1)
  d2ldv2<-matrix(0,1,1)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-(2*v+p)/d[i]*pgamma(1,v+1+p/2,scale=2/d[i])/pgamma(1,v+p/2,scale=2/d[i])
    vu[i]<-(2*v+p+2)*(2*v+p)/d[i]^2*pgamma(1,v+2+p/2,scale=2/d[i])/pgamma(1,v+p/2,scale=2/d[i])-((2*v+p)/d[i])^2*(pgamma(1,v+1+p/2,scale=2/d[i])/pgamma(1,v+p/2,scale=2/d[i]))^2
    auxi<-pmax(t(h)%*%invB%*%e[i,],-37)
    W[i]<-dnorm(auxi)/pnorm(auxi)
    Wl[i]<--W[i]*(auxi+W[i])
    d2ldb2<-d2ldb2-as.numeric(u[i])*t(X[[i]])%*%invB%*%invB%*%X[[i]]+as.numeric(vu[i])*t(X[[i]])%*%invB%*%invB%*%e[i,]%*%t(e[i,])%*%invB%*%invB%*%X[[i]]+as.numeric(Wl[i])*t(X[[i]])%*%invB%*%h%*%t(h)%*%invB%*%X[[i]]
    d2ldbdh<-d2ldbdh-as.numeric(W[i])*t(X[[i]])%*%invB-as.numeric(Wl[i])*t(X[[i]])%*%invB%*%h%*%t(e[i,])%*%invB
    d2ldh2<-d2ldh2+as.numeric(Wl[i])*invB%*%e[i,]%*%t(e[i,])%*%invB
    ddida<-matrix(0,p0,1)
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
      ddida[j]=-t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]
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
    d2ldhda<-d2ldhda+as.numeric(W[i])*d2Aidhda+as.numeric(Wl[i])*invB%*%e[i,]%*%t(dAida)
    d2ldbda<-d2ldbda-1/2*as.numeric(u[i])*d2didbda-1/2*as.numeric(vu[i])*t(X[[i]])%*%invB%*%invB%*%e[i,]%*%t(ddida)+as.numeric(W[i])*d2Aidbda-as.numeric(Wl[i])*t(X[[i]])%*%invB%*%h%*%t(dAida)
    d2lda2<-d2lda2-1/2*d2ldda2-1/2*as.numeric(u[i])*d2dida2+1/4*as.numeric(vu[i])*ddida%*%t(ddida)+as.numeric(W[i])*d2Aida2+as.numeric(Wl[i])*dAida%*%t(dAida)
    d2ldbdv<-d2ldbdv+(integrate(fauxss,lower=0,upper=1,v,p,d[i],0,1)$val/integrate(fauxss,lower=0,upper=1,v,p,d[i],1,0)$val-integrate(fauxss,lower=0,upper=1,v,p,d[i],1,1)$val*integrate(fauxss,lower=0,upper=1,v,p,d[i],0,0)$val/integrate(fauxss,lower=0,upper=1,v,p,d[i],1,0)$val^2)*t(X[[i]])%*%invB%*%invB%*%e[i,]
    d2ldadv<-d2ldadv-1/2*(integrate(fauxss,lower=0,upper=1,v,p,d[i],0,1)$val/integrate(fauxss,lower=0,upper=1,v,p,d[i],1,0)$val-integrate(fauxss,lower=0,upper=1,v,p,d[i],1,1)$val*integrate(fauxss,lower=0,upper=1,v,p,d[i],0,0)$val/integrate(fauxss,lower=0,upper=1,v,p,d[i],1,0)$val^2)*ddida
    d2ldv2<-d2ldv2+integrate(fauxss,lower=0,upper=1,v,p,d[i],1,2)$val/integrate(fauxss,lower=0,upper=1,v,p,d[i],1,0)$val-(integrate(fauxss,lower=0,upper=1,v,p,d[i],1,1)$val/integrate(fauxss,lower=0,upper=1,v,p,d[i],1,0)$val)^2-1/v^2
  }  
  d2ldb=cbind(d2ldb2,d2ldbda,d2ldbdh,d2ldbdv)
  d2lda=cbind(t(d2ldbda),d2lda2,t(d2ldhda),d2ldadv)
  d2ldh=cbind(t(d2ldbdh),d2ldhda,d2ldh2,matrix(0,p,1))
  d2ldv=cbind(t(d2ldbdv),t(d2ldadv),matrix(0,1,p),d2ldv2)
  Hl=rbind(d2ldb,d2lda,d2ldh,d2ldv)
  return(Hl)
} #matriz de informação da slash normal assimetrica multivariada com regressao

#mifmssr(P,y,X) #teste
#mifmssr(P0,y) #sem explicativas


#Algoritmos EM para o modelo slash normal assimetrico multivariado com regressao

#Padrao

EMDMSSRP<-function(y,X)
{
  #t0=Sys.time()
  aa=system.time({
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b0<-matrix(0,q,q)
  b1<-matrix(0,q,1)
  for(i in 1:n){
    b0<-b0+t(X[[i]])%*%X[[i]]
    b1<-b1+t(X[[i]])%*%y[i,] 
  }
  b<-solve(b0)%*%b1
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
  }
  S<-cov(e)
  invS<-solve(S)
  B<-sqrtm(S)
  invB<-solve(B)
  h<-as.matrix(skewness(e))
  v<-1.01
  P_0<-c(as.vector(b),vech(B),as.vector(h),v)
  crit<-1
  iter<-0
  while((crit>=10^(-6))&&(iter<=5000))
  { 
    iter<-iter+1 
    d<-matrix(0,n,1)
    e<-matrix(0,n,p)
    for(i in 1:n){
      e[i,]<-y[i,]-X[[i]]%*%b
      d[i]<-t(e[i,])%*%invS%*%e[i,]
    }
    u<-((2*v+p)/d)*(pgamma(1,v+1+p/2,scale=2/d)/pgamma(1,v+p/2,scale=2/d))
    aux<-pmax(as.vector(e%*%invB%*%h),-37)
    W<-as.matrix(dnorm(aux)/pnorm(aux))
    t<-e%*%invB%*%h+W #passo E
    S<-1/n*(t(e)%*%diag(as.vector(u),n)%*%e)
    invS<-solve(S)
    B<-sqrtm(S)
    invB<-solve(B)
    h<-B%*%solve(t(e)%*%e)%*%t(e)%*%t
    b0<-matrix(0,q,q)
    b1<-matrix(0,q,1)
    for(i in 1:n){
      b0<-b0+t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%X[[i]]
      b1<-b1+t(X[[i]])%*%(as.numeric(u[i])*(invB%*%invB)%*%y[i,]-as.numeric(t[i])*invB%*%h+invB%*%h%*%t(h)%*%invB%*%y[i,])
    }
    b<-solve(b0)%*%b1 #passo M
    V<-optim(v,lmssvr,gr=NULL,b,B,h,y,X,method="L-BFGS-B",lower=0.1,upper=40,control=list(fnscale=-1,maxit=50))
    v<-as.numeric(V$par) 
    P<-c(as.vector(b),vech(B),as.vector(h),v)
    crit<-sqrt(sum((P-P_0)^2))
    #crit<-sqrt(sum(escmssr(P,y,X)*escmssr(P,y,X)))
    #crit<-sqrt(sum((lmssr(P,y,X)-lmssr(P_0,y,X))^2))
    #print(crit)
    P_0<-P
  }# o algoritmo EM
  })
  #t1=Sys.time()
  estmax=P
  logvero=lmssr(P,y,X)
  escore=escmssr(P,y,X)
  informacao=mifmssr(P,y,X)   
  tempo=as.numeric(aa[3])#as.numeric(t1-t0)
  object.out<-list(estmax=estmax,logvero=logvero,escore=escore,informacao=informacao,iter=iter,tempo=tempo)
} 

#X<-run(n)
#y<-rmssr(n,beta,rSigma,lambda,nu,X)
#theta<-tryCatch(EMDMSSRP(y,X))
#y<-rmssr(n,beta[1],Sigma,lambda,nu)
#theta<-tryCatch(EMDMSSRP(y))
#print(theta)
#P
#lmssr(P,y,X)
#P0
#lmssr(P0,y)

#Modificado

EMDMSSRM<-function(y,X)
{
  #t0=Sys.time()
  aa=system.time({
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b0<-matrix(0,q,q)
  b1<-matrix(0,q,1)
  for(i in 1:n){
    b0<-b0+t(X[[i]])%*%X[[i]]
    b1<-b1+t(X[[i]])%*%y[i,]
  }
  b<-solve(b0)%*%b1
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
  }
  S<-cov(e)
  invS<-solve(S)
  B<-sqrtm(S)
  invB<-solve(B)
  h<-as.matrix(skewness(e))
  v<-1.01
  P_0<-c(as.vector(b),vech(B),as.vector(h),v)
  crit<-1
  iter<-0
  while ((crit>=10^(-6))&&(iter<=5000))
  { 
    iter<-iter+1 
    d<-matrix(0,n,1)
    e<-matrix(0,n,p)
    lu<-0
    for(i in 1:n){
      e[i,]<-y[i,]-X[[i]]%*%b
      d[i]<-t(e[i,])%*%invS%*%e[i,]
      #lu[i]<-integrate(fauxss,lower=0,upper=1,v,p,d[i],1,1)$val/integrate(fauxss,lower=0,upper=1,v,p,d[i],1,0)$val
      #lu[i]<-integrate(fauxss,lower=0,upper=1,v,p,d[i],1,1)$val/pgamma(1,v+p/2,scale=2/d[i])*(d[i]/2)^(v+p/2)/gamma(v+p/2)
      lu[i]<-integrate(elogu,lower=qtgamma(0.001,v+p/2,2/d[i],0,1),upper=qtgamma(0.999,v+p/2,2/d[i],0,1),v,p,d[i])$val
    }
    u<-((2*v+p)/d)*(pgamma(1,v+1+p/2,scale=2/d)/pgamma(1,v+p/2,scale=2/d))
    aux<-pmax(as.vector(e%*%invB%*%h),-37)
    W<-as.matrix(dnorm(aux)/pnorm(aux))
    t<-e%*%invB%*%h+W #passo E
    S<-1/n*(t(e)%*%diag(as.vector(u))%*%e)
    invS<-solve(S)
    B<-sqrtm(S)
    invB<-solve(B)
    h<-B%*%solve(t(e)%*%e)%*%t(e)%*%t
    b0<-matrix(0,q,q)
    b1<-matrix(0,q,1)
    for(i in 1:n){
      b0<-b0+t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%X[[i]]
      b1<-b1+t(X[[i]])%*%(as.numeric(u[i])*invS%*%y[i,]-as.numeric(t[i])*invB%*%h+invB%*%h%*%t(h)%*%invB%*%y[i,])  
    }
    C<-solve(b0)
    b<-C%*%b1
    v<--n/sum(lu)
    P<-c(as.vector(b),vech(B),as.vector(h),v)
    crit<-sqrt(sum((P-P_0)^2))
    #crit<-sqrt(sum(escmssr(P,y,X)*escmssr(P,y,X)))
    #crit<-sqrt(sum((lmssr(P,y,X)-lmssr(P_0,y,X))^2))
    #print(crit)
    P_0<-P
  }# o algoritmo EM
  })
  #t1=Sys.time()
  estmax=P
  logvero=lmssr(P,y,X)
  escore=0#escmssr(P,y,X)
  informacao=0#mifmssr(P,y,X)
  tempo=as.numeric(aa[3])#as.numeric(t1-t0)
  object.out<-list(estmax=estmax,logvero=logvero,escore=escore,informacao=informacao,iter=iter,tempo=tempo,u=u,C=C)
} 

#X<-run(n)
#y<-rmssr(n,beta,rSigma,lambda,nu,X)
#theta<-tryCatch(EMDMSSRM(y,X))
#y<-rmssr(n,beta[1],rSigma,lambda,nu)
#theta<-tryCatch(EMDMSSRM(y))
#print(theta)
#P
#lmssr(P,y,X)
#P0
#lmssr(P0,y)

#Tradicional

EMDMSSRT<-function(y,X)
{
  t0=Sys.time()
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b0<-matrix(0,q,q)
  b1<-matrix(0,q,1)
  for(i in 1:n){
    b0<-b0+t(X[[i]])%*%X[[i]]
    b1<-b1+t(X[[i]])%*%y[i,] 
  }
  b<-solve(b0)%*%b1
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
  }
  S<-cov(e)
  invS<-solve(S)
  B<-sqrtm(S)
  invB<-solve(B)
  h<-as.matrix(skewness(e))
  v<-1.01
  P_0<-c(as.vector(b),vech(B),as.vector(h),v)
  O<-optim(P_0,lmssr,gr=escmssr,y,X,method="L-BFGS-B",lower=abind(rep(-Inf,q),0,-Inf,0,rep(-Inf,p),0.001),upper=rep(Inf,q+p*(p+1)/2+p+1),control=list(fnscale=-1,maxit=5000,factr=1e3))
  P<-O$par
  estmax=P
  logvero=lmssr(P,y,X)
  escore=escmssr(P,y,X)
  informacao=mifmssr(P,y,X)
  iter=as.numeric(O$counts[1])
  t1=Sys.time()
  tempo=as.numeric(t1-t0)
  object.out<-list(estmax=estmax,logvero=logvero,escore=escore,informacao=informacao,iter=iter,tempo=tempo)
} 

#Envelopes skew slash normal

envelSS<-function(Y,X,beta,Sigma,nu){
  n=nrow(Y)
  p=ncol(Y)
  d2<-matrix(0,n,1)
  for(i in 1:n){
    d2[i]=t(Y[i,]-X[[i]]%*%beta)%*%solve(Sigma)%*%(Y[i,]-X[[i]]%*%beta)
  }
  d2s1=sort(d2)
  d2s=t(d2s1)
  xq2 <- qchisq(ppoints(n), p)
  Xsim<-matrix(0,200,n)
  for(i in 1:200){
    u0<-rep(0,n)
    aux<-matrix(0,n,p)
    for(j in 1:n){
     u0[j]=rbeta(1,shape1=nu,shape2=1)
     aux[j,] <- mvrnorm(1,rep(0, p),diag(p)/u0[j])}
     Xsim[i,]<-apply(aux*aux,1,sum)}
  Xsim2<-apply(Xsim,1,sort)
  d21<-matrix(0,n,1)
  d22<-matrix(0,n,1)
  for(i in 1:n){
    d21[i]  <- quantile(Xsim2[i,],0.025)
    d22[i]  <- quantile(Xsim2[i,],0.975)}
  d2med <-apply(Xsim2,1,mean)
  fy <- range(d2s,d21,d22)
  plot(xq2,d2s,xlab = expression(paste("Theoretical ",chi[p]^2, " quantiles")),
       ylab="Sample values and sibetalated envelope",pch=20,ylim=fy)
  par(new=T)
  plot(xq2,d21,type="l",ylim=fy,xlab="",ylab="")
  par(new=T)
  plot(xq2,d2med,type="l",ylim=fy,xlab="",ylab="")
  par(new=T)
  plot(xq2,d22,type="l",ylim=fy,xlab="",ylab="")
  print(sum((as.vector(d2s)<as.vector(d21))|(as.vector(d2s)>as.vector(d22))))
}

envel0SS <-function(Y,X,beta,Sigma,lambda,nu){
  n=nrow(Y)
  p=ncol(Y)
  q=ncol(X[[n]])
  d2<-matrix(0,n,1)
  for(i in 1:n){
    d2[i]=t(Y[i,]-X[[i]]%*%beta)%*%solve(Sigma)%*%(Y[i,]-X[[i]]%*%beta)
  }
  d2s1=sort(d2)
  d2s=t(d2s1)
  xq2 <- qchisq(ppoints(n), p)
  Xsim<-matrix(0,100,n)
  for(i in 1:100){
   Y0=rmssr(n,beta,sqrtm(Sigma),lambda,nu,X)
   Theta0=EMDMSSRM(Y0,X)$estmax
   d20<-matrix(0,n,1)
   for(j in 1:n){
    d20[j]=t(Y[j,]-X[[j]]%*%Theta0[1:q])%*%solve(xpnd(Theta0[(q+1):(q+p*(p+1)/2)])%*%xpnd(Theta0[(q+1):(q+p*(p+1)/2)]))%*%(Y[j,]-X[[j]]%*%Theta0[1:q])
   }
  Xsim[i,]<-d20}
  Xsim2<-apply(Xsim,1,sort)
  d21<-matrix(0,n,1)
  d22<-matrix(0,n,1)
  for(i in 1:n){
    d21[i]  <- quantile(Xsim2[i,],0.025)
    d22[i]  <- quantile(Xsim2[i,],0.975)}
  d2med <-apply(Xsim2,1,mean)

  fy <- range(d2s,d21,d22)
  plot(xq2,d2s,xlab = expression(paste("Theoretical ",P(theta), " quantiles")),
       ylab="Sample values and sibetalated envelope",pch=20,ylim=fy)
  par(new=T)
  plot(xq2,d21,type="l",ylim=fy,xlab="",ylab="")
  par(new=T)
  plot(xq2,d2med,type="l",ylim=fy,xlab="",ylab="")
  par(new=T)
  plot(xq2,d22,type="l",ylim=fy,xlab="",ylab="")
  print(sum((as.vector(d2s)<as.vector(d21))|(as.vector(d2s)>as.vector(d22))))
}


