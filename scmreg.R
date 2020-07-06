#Estimacao de maxima verossimilhanca via EM para um modelo normal contaminado assimetrico multivariado com regressao

#Geracao dos dados

rmscr<-function(n,b,B,h,tau,X){
  v=tau[1]
  g=tau[2]
  p=length(h)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  a<-rep(0,n)
  y<-matrix(0,n,p)
  for(i in 1:n){
    a[i]<-runif(1,0,1) #numero aleatorio 
    if(a[i]<v){
      y[i,]=X[[i]]%*%as.matrix(b)+1/sqrt(g)*B%*%((1+sum(h*h/g))^(-1/2)*(h/sqrt(g)*abs(rnorm(1)))+solve(sqrtm(diag(p)+h%*%t(h)/g))%*%mvrnorm(1,rep(0,p),diag(p)))
    }
    else
      y[i,]=X[[i]]%*%as.matrix(b)+B%*%((1+sum(h*h))^(-1/2)*(h*abs(rnorm(1)))+solve(sqrtm(diag(p)+h%*%t(h)))%*%mvrnorm(1,rep(0,p),diag(p)))
  } 
  return(y) 
}

#y<-rmscr(n,beta,rSigma,lambda,nu,gamma,X) 
#y<-rmscr(n,beta[1],rSigma,lambda,nu,gamma) #geracao independente de X
#y #dados observados resposta

#Funcoes dos dados incompletos

dmscr<-function(P,y,X)
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
  #d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    #d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
  }
  aux=pmax(as.vector(t(h)%*%invB%*%t(e)),-37)
  f=2*(v/det(1/sqrt(g)*B)*dmvnorm(sqrt(g)*e%*%invB)+(1-v)/det(B)*dmvnorm(e%*%invB))*pnorm(aux)
  return(f)
} #densidade da normal contaminada assimetrica multivariada com regressao

lmscr<-function(P,y,X)
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
  #d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    #d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
  }
  aux=pmax(as.vector(t(h)%*%invB%*%t(e)),-37)
  l=sum(log(2*(v/det(1/sqrt(g)*B)*dmvnorm(sqrt(g)*e%*%invB)+(1-v)/det(B)*dmvnorm(e%*%invB))*pnorm(aux)))
  return(l)
} #log-verossimilhanca da normal contaminada assimetrica multivariada com regressao

#lmscr(P,y,X) #teste
#lmscr(P0,y) #sem explicativas

lmscvgr<-function(PA,b,B,h,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  v=PA[1]
  g=PA[2]
  invB=solve(B)
  #d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%as.matrix(b)
    #d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
  }
  lvg=sum(log(v/det(1/sqrt(g)*B)*dmvnorm(sqrt(g)*e%*%invB)+(1-v)/det(B)*dmvnorm(e%*%invB)))
  return(lvg)
} #log-verossimilhanca da normal contaminada assimetrica multivariada com regressao em nu e gamma

#lmscvgr(c(nu,gamma),beta,rSigma,lambda,y,X) #teste
#lmscvgr(c(nu,gamma),beta[1],Sigma,lambda,y) #sem explicativas

escmscr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b=as.matrix(P[1:q])
  p0=p*(p+1)/2
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p*(p+1)/2+1):(length(P)-2)])
  v=as.numeric(P[length(P)-1])
  g=as.numeric(P[length(P)])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  W<-matrix(0,n,1)
  K<-matrix(0,n,1)
  u<-matrix(0,n,1)
  dldb1<-matrix(0,q,1)
  dldb2<-matrix(0,q,1)
  dlda<-matrix(0,p0,1)
  dldv<-0
  dldg<-0
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    K[i]<-sqrt(2*pi)*(v*g^(p/2)*dnorm(sqrt(g*d[i]))+(1-v)*dnorm(sqrt(d[i])))
    u[i]<-sqrt(2*pi)*(v*g^(p/2+1)*dnorm(sqrt(g*d[i]))+(1-v)*dnorm(sqrt(d[i])))/K[i]
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
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
    dldv<-dldv+sqrt(2*pi)*(g^(p/2)*dnorm(sqrt(g*d[i]))-dnorm(sqrt(d[i])))/K[i]
    dldg<-dldg+sqrt(2*pi)*((p-g*d[i])*v*g^(p/2-1)*dnorm(sqrt(g*d[i])))/(2*K[i])
  }
  dldb=dldb1-dldb2
  dldh=invB%*%t(e)%*%W
  el=c(as.vector(dldb),as.vector(dlda),as.vector(dldh),dldv,dldg)
  return(el)
} #escore da normal contaminada assimetrica multivariada com regressao

#escmscr(P,y,X) #teste
#escmscr(P0,y) #sem explicativas

mifmscr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p*(p+1)/2+1):(length(P)-2)])
  v=as.numeric(P[length(P)-1])
  g=as.numeric(P[length(P)])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  K<-matrix(0,n,1)
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
  d2ldbdg<-matrix(0,q,1)
  d2ldadg<-matrix(0,p0,1)
  d2ldv2<-0
  d2ldvdg<-0
  d2ldg2<-0
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    K[i]<-sqrt(2*pi)*(v*g^(p/2)*dnorm(sqrt(g*d[i]))+(1-v)*dnorm(sqrt(d[i])))
    u[i]<-sqrt(2*pi)*(v*g^(p/2+1)*dnorm(sqrt(g*d[i]))+(1-v)*dnorm(sqrt(d[i])))/K[i]
    vu[i]<-sqrt(2*pi)*(v*(1-v)*(1-g)^2*g^(p/2)*dnorm(sqrt((g+1)*d[i])))/K[i]^2
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
    d2ldbdv<-d2ldbdv-sqrt(2*pi)*((1-g)*g^(p/2)*dnorm(sqrt((g+1)*d[i])))/K[i]^2*t(X[[i]])%*%invB%*%invB%*%e[i,]
    d2ldadv<-d2ldadv+sqrt(2*pi)*((1-g)*g^(p/2)*dnorm(sqrt((g+1)*d[i])))/(2*K[i]^2)*ddida
    d2ldbdg<-d2ldbdg-sqrt(2*pi)*v*((1-v)*(p-g*(2+p+d[i])+g^2*d[i])*g^(p/2-1)*dnorm(sqrt((g+1)*d[i]))-2*v*g^p*dnorm(sqrt(2*g*d[i])))/(2*K[i]^2)*t(X[[i]])%*%invB%*%invB%*%e[i,]
    d2ldadg<-d2ldadg+sqrt(2*pi)*v*((1-v)*(p-g*(2+p+d[i])+g^2*d[i])*g^(p/2-1)*dnorm(sqrt((g+1)*d[i]))-2*v*g^p*dnorm(sqrt(2*g*d[i])))/(4*K[i]^2)*ddida
    d2ldv2<-d2ldv2-(sqrt(2*pi)*(g^(p/2)*dnorm(sqrt(g*d[i]))-dnorm(sqrt(d[i])))/K[i])^2
    d2ldvdg<-d2ldvdg+sqrt(2*pi)*(p-g*d[i])*g^(p/2-1)*dnorm(sqrt((g+1)*d[i]))/(2*K[i]^2)  
    d2ldg2<-d2ldg2-sqrt(2*pi)*v*(2*p*v*g^(p-2)*dnorm(sqrt(2*g*d[i]))-(1-v)*((p-g*d[i])^2-2*p)*g^(p/2-2)*dnorm(sqrt((g+1)*d[i])))/(4*K[i]^2)
  }  
  d2ldb=cbind(d2ldb2,d2ldbda,d2ldbdh,d2ldbdv,d2ldbdg)
  d2lda=cbind(t(d2ldbda),d2lda2,t(d2ldhda),d2ldadv,d2ldadg)
  d2ldh=cbind(t(d2ldbdh),d2ldhda,d2ldh2,matrix(0,p,2))
  d2ldv=cbind(t(d2ldbdv),t(d2ldadv),matrix(0,1,p),as.matrix(d2ldv2),as.matrix(d2ldvdg))
  d2ldg=cbind(t(d2ldbdg),t(d2ldadg),matrix(0,1,p),as.matrix(d2ldvdg),as.matrix(d2ldg2))
  Hl=rbind(d2ldb,d2lda,d2ldh,d2ldv,d2ldg)
  return(Hl)
} #matriz de informação da normal contaminada assimetrica multivariada com regressao

#mifmscr(P,y,X)-xpnd(vech(mifmscr(P,y,X)))
#mifmscr(P0,y) #sem explicativas


#Algoritmos EM para o modelo normal contaminado assimetrico multivariado com regressao

#Padrao

EMDMSCRP<-function(y,X)
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
  v<-0.5
  g<-0.5
  P_0<-c(as.vector(b),vech(B),as.vector(h),v,g)
  crit<-1
  iter<-0
  while((crit>=10^(-6))&&(iter<=5000))
  { 
    iter<-iter+1 
    u<-(g*v/det(1/sqrt(g)*B)*dmvnorm(sqrt(g)*e%*%invB)+(1-v)/det(B)*dmvnorm(e%*%invB))/(v/det(1/sqrt(g)*B)*dmvnorm(sqrt(g)*e%*%invB)+(1-v)/det(B)*dmvnorm(e%*%invB))
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
    b<-solve(b0)%*%b1
    #d<-matrix(0,n,1)
    e<-matrix(0,n,p)
    for(i in 1:n){
      e[i,]<-y[i,]-X[[i]]%*%b
      #d[i]<-t(e[i,])%*%invS%*%e[i,]
    }
    U<-optim(c(v,g),lmscvgr,gr=NULL,b,B,h,y,X,method="L-BFGS-B",lower=c(0.01,0.01),upper=c(0.99,0.99),control=list(fnscale=-1,maxit=50))
    v<-as.vector(U$par[1])
    g<-as.vector(U$par[2]) #passo M 
    P<-c(as.vector(b),vech(B),as.vector(h),v,g)
    crit<-sqrt(sum((P-P_0)^2))
    P_0<-P
  }# o algoritmo EM
  })
  #t1=Sys.time()
  estmax=P
  logvero=lmscr(P,y,X)
  escore=escmscr(P,y,X)
  informacao=mifmscr(P,y,X)
  tempo=as.numeric(aa[3])#as.numeric(t1-t0)
  object.out<-list(estmax=estmax,logvero=logvero,escore=escore,informacao=informacao,iter=iter,tempo=tempo)
} 

#X<-run(n)
#y<-rmscr(n,beta,rSigma,lambda,nu,gamma,X)  
#theta<-tryCatch(EMDMSCRP(y,X))
#y<-rmscreg(n,beta[1],rSigma,lambda,nu,gamma)
#theta<-tryCatch(EMDMSCRP(y)))
#print(theta)
#P
#lmscr(P,y,X)
#P0
#lmscr(P0,y)

#Modificado

EMDMSCRM<-function(y,X)
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
  B<-sqrtm(S)
  invB<-solve(B)
  h<-as.matrix(skewness(e))
  v<-0.5
  g<-0.5
  P_0<-c(as.vector(b),vech(B),as.vector(h),v,g)
  crit<-1
  iter<-0
  while((crit>=10^(-6))&&(iter<=5000))
  { 
    iter<-iter+1 
    v1<-v/det(1/sqrt(g)*B)*dmvnorm(sqrt(g)*e%*%invB)/(v/det(1/sqrt(g)*B)*dmvnorm(sqrt(g)*e%*%invB)+(1-v)/det(B)*dmvnorm(e%*%invB))
    v2<-(1-v)/det(B)*dmvnorm(e%*%invB)/(v/det(1/sqrt(g)*B)*dmvnorm(sqrt(g)*e%*%invB)+(1-v)/det(B)*dmvnorm(e%*%invB))
    u<-g*v1+v2
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
    d<-matrix(0,n,1)
    e<-matrix(0,n,p)
    for(i in 1:n){
      e[i,]<-y[i,]-X[[i]]%*%b
      d[i]<-t(e[i,])%*%invS%*%e[i,]
    }
    v<-sum(v1)/n
    g<-n*p*v/sum(v1*d) #passo M  
    P<-c(as.vector(b),vech(B),as.vector(h),v,g)
    crit<-sqrt(sum((P-P_0)^2))
    P_0<-P
  }# o algoritmo EM
  })
  #t1=Sys.time()
  estmax=P
  logvero=lmscr(P,y,X)
  escore=0#escmscr(P,y,X)
  informacao=0#mifmscr(P,y,X)
  tempo=as.numeric(aa[3])#as.numeric(t1-t0)
  object.out<-list(estmax=estmax,logvero=logvero,escore=escore,informacao=informacao,iter=iter,tempo=tempo,u=u,C=C)
} 

#X<-run(n)
#y<-rmscr(n,beta,rSigma,lambda,nu,gamma,X)
#theta<-tryCatch(EMDMSCRM(y,X))
#y<-rmscr(n,beta[1],rSigma,lambda,nu,gamma)
#theta<-tryCatch(EMDMSCRM(y))
#print(theta)
#P
#lmscr(P,y,X)
#P0
#lmscr(P0,y)

EMDMSCRT<-function(y,X)
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
  v<-0.5
  g<-0.5
  P_0<-c(as.vector(b),vech(B),as.vector(h),v,g)
  O<-optim(P_0,lmscr,gr=escmscr,y,X,method="L-BFGS-B",lower=abind(rep(-Inf,q),0,-Inf,0,rep(-Inf,p),0,0.001),upper=abind(rep(Inf,q+p*(p+1)/2+p),1,1),control=list(fnscale=-1,maxit=5000,factr=1e3))
  P<-O$par
  estmax=P
  logvero=lmscr(P,y,X)
  escore=escmscr(P,y,X)
  informacao=mifmscr(P,y,X)
  iter=as.numeric(O$counts[1])
  t1=Sys.time()
  tempo=as.numeric(t1-t0)
  object.out<-list(estmax=estmax,logvero=logvero,escore=escore,informacao=informacao,iter=iter,tempo=tempo)
} 

#Envelopes skew contaminated normal

envelSC<-function(Y,X,beta,Sigma,nu,gama){
  n=nrow(Y)
  p=ncol(Y)
  d2<-matrix(0,n,1)
  for(i in 1:n){
    d2[i]=t(Y[i,]-X[[i]]%*%beta)%*%solve(Sigma)%*%(Y[i,]-X[[i]]%*%beta)
  }
  d2s=sort(d2)
  d2s=t(d2s)
  xq2 <- qchisq(ppoints(n), p)
  Xsim<-matrix(0,200,n)
  for(i in 1:200){
    u1<-rchisq(n, p)/gama
    u2<-rchisq(n, p)
    u3<-runif(n)
    id<-(u3<nu)
    u2[id]<-u1[id]
    Xsim[i,]=u2}
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

envel0SC<-function(Y,X,beta,Sigma,lambda,nu,gama){
  n=nrow(Y)
  p=ncol(Y)
  q=ncol(X[[n]])
  d2<-matrix(0,n,1)
  for(i in 1:n){
    d2[i]=t(Y[i,]-X[[i]]%*%beta)%*%solve(Sigma)%*%(Y[i,]-X[[i]]%*%beta)
  }
  d2s=sort(d2)
  d2s=t(d2s)
  xq2 <- qchisq(ppoints(n), p)
  Xsim<-matrix(0,100,n)
  for(i in 1:100){
  Y0=rmscr(n,beta,sqrtm(Sigma),lambda,nu,gama,X)
  Theta0=EMDMSCRM(Y0,X)$estmax
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

