#Estimacao de maxima verossimilhanca via EM para um modelo normal assimetrico multivariado com regressao

#n=200
#p=2
#q=3
#beta=c(40,-30,15)
#rSigma=matrix(c(2.5,-1,-1,1.5),p,p) 
#Sigma=rSigma%*%rSigma
#lambda=c(2,-1)

#P<-c(beta,vech(rSigma),lambda)
#P0<-c(beta[1],vech(rSigma),lambda) #sem variavel explicativa

#Geracao dos dados

rmsnr<-function(n,b,B,h,X){
  p=length(h)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  y<-matrix(0,n,p)
  for(i in 1:n){
      y[i,]=X[[i]]%*%as.matrix(b)+B%*%((1+sum(h*h))^(-1/2)*(h*abs(rnorm(1)))+solve(sqrtm(diag(p)+h%*%t(h)))%*%mvrnorm(1,rep(0,p),diag(p)))
  } 
  return(y) 
}

#y<-rmsnr(n,beta,rSigma,lambda,X) 
#y<-rmscr(n,beta[1],rSigma,lambda) #geracao independente de X
#y #dados observados resposta

#Funcoes dos dados incompletos

dmsnr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p*(p+1)/2)])
  invB=solve(B)
  h=as.matrix(P[(q+p*(p+1)/2+1):length(P)])
  #d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    #d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
  }
  aux=pmax(as.vector(t(h)%*%invB%*%t(e)),-37)
  f=2/det(B)*dmvnorm(e%*%invB)*pnorm(aux)
  return(f)
} #densidade da normal assimetrica multivariada com regressao

lmsnr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p*(p+1)/2)])
  invB=solve(B)
  h=as.matrix(P[(q+p*(p+1)/2+1):length(P)])
  #d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    #d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
  }
  aux=pmax(as.vector(t(h)%*%invB%*%t(e)),-37)
  l=sum(log(2/det(B)*dmvnorm(e%*%invB)*pnorm(aux)))
  return(l)
} #log-verossimilhanca da normal assimetrica multivariada com regressao

#lmsnr(P,y,X) #teste
#lmscr(P0,y) #sem explicativas

escmsnr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b=as.matrix(P[1:q])
  p0=p*(p+1)/2
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p*(p+1)/2+1):length(P)])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  W<-matrix(0,n,1)
  u<-matrix(0,n,1)
  dldb1<-matrix(0,q,1)
  dldb2<-matrix(0,q,1)
  dlda<-matrix(0,p0,1)
  dldv<-0
  dldg<-0
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-1
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
   }
  dldb=dldb1-dldb2
  dldh=invB%*%t(e)%*%W
  el=c(as.vector(dldb),as.vector(dlda),as.vector(dldh))
  return(el)
} #escore da normal contaminada assimetrica multivariada com regressao

#escmsnr(P,y,X) #teste
#escmsnr(P0,y) #sem explicativas

mifmsnr<-function(P,y,X)
{
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p*(p+1)/2+1):length(P)])
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
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-1
    vu[i]<-0
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
  }  
  d2ldb=cbind(d2ldb2,d2ldbda,d2ldbdh)
  d2lda=cbind(t(d2ldbda),d2lda2,t(d2ldhda))
  d2ldh=cbind(t(d2ldbdh),d2ldhda,d2ldh2)
  Hl=rbind(d2ldb,d2lda,d2ldh)
  return(Hl)
} #matriz de informação da normal assimetrica multivariada com regressao

#mifmsnr(P,y,X) #teste
#mifmsnr(P0,y) #sem explicativas


#Algoritmos EM para o modelo normal assimetrico multivariado com regressao

#Padrao (normal)

EMDMSNRP<-function(y,X)
{
  t0=Sys.time()
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  array(1,c(p,1,n))
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
  h<-matrix(0,p,1)
  P_0<-c(as.vector(b),vech(B),as.vector(h))
  crit<-1
  iter<-0
  while(crit>=10^(-6))
  { 
    iter<-iter+1 
    u<-rep(1,n)
    aux<-pmax(as.vector(e%*%invB%*%h),-37)
    #W<-as.matrix(dnorm(aux)/pnorm(aux))
    t<-rep(0,n) #passo E
    S<-1/n*(t(e)%*%diag(as.vector(u),n)%*%e)
    invS<-solve(S)
    B<-sqrtm(S)
    invB<-solve(B)
    h<-matrix(0,p,1)
    b0<-matrix(0,q,q)
    b1<-matrix(0,q,1)
    for(i in 1:n){
      b0<-b0+t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%X[[i]]
      b1<-b1+t(X[[i]])%*%(as.numeric(u[i])*(invB%*%invB)%*%y[i,]-as.numeric(t[i])*invB%*%h+invB%*%h%*%t(h)%*%invB%*%y[i,])
    }
    C<-solve(b0)
    b<-C%*%b1
    #d<-matrix(0,n,1)
    e<-matrix(0,n,p)
    for(i in 1:n){
      e[i,]<-y[i,]-X[[i]]%*%b
      #d[i]<-t(e[i,])%*%invS%*%e[i,]
    } #passo M 
    P<-c(as.vector(b),vech(B),as.vector(h))
    crit<-sqrt(sum((P-P_0)^2))
    P_0<-P
  }# o algoritmo EM
  t1=Sys.time()
  estmax=P
  logvero=lmsnr(P,y,X)
  escore=escmsnr(P,y,X)
  informacao=mifmsnr(P,y,X)
  tempo=as.numeric(t1-t0)
  object.out<-list(estmax=estmax,logvero=logvero,escore=escore,informacao=informacao,iter=iter,tempo=tempo,u=u,C=C)
}

#Modificado (normal assimetrica)

EMDMSNRM<-function(y,X)
{
  t0=Sys.time()
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  #array(1,c(p,1,n))
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
  P_0<-c(as.vector(b),vech(B),as.vector(h))
  crit<-1
  iter<-0
  while(crit>=10^(-6))
  { 
    iter<-iter+1 
    u<-rep(1,n)
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
    C<-solve(b0)
    b<-C%*%b1
    #d<-matrix(0,n,1)
    e<-matrix(0,n,p)
    for(i in 1:n){
      e[i,]<-y[i,]-X[[i]]%*%b
      #d[i]<-t(e[i,])%*%invS%*%e[i,]
    } #passo M 
    P<-c(as.vector(b),vech(B),as.vector(h))
    crit<-sqrt(sum((P-P_0)^2))
    P_0<-P
  }# o algoritmo EM
  t1=Sys.time()
  estmax=P
  logvero=lmsnr(P,y,X)
  escore=0#escmsnr(P,y,X)
  informacao=0#mifmsnr(P,y,X)
  tempo=as.numeric(t1-t0)
  object.out<-list(estmax=estmax,logvero=logvero,escore=escore,informacao=informacao,iter=iter,tempo=tempo,u=u,C=C)
} 

#X<-run(n)
#y<-rmsnr(n,beta,rSigma,lambda,X)  
#theta<-tryCatch(EMDMSNRP(y,X))
#print(theta)
#y<-rmsnr(n,beta[1],rSigma,lambda)
#theta<-tryCatch(EMDMSNRP(y)))
#print(theta)
#P
#lmscr(P,y,X)
#P0
#lmscr(P0,y)

EMDMSNRT<-function(y,X)
{
  t0=Sys.time()
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  array(1,c(p,1,n))
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
  P_0<-c(as.vector(b),vech(B),as.vector(h))
  O<-optim(P_0,lmsnr,gr=escmsnr,y,X,method="L-BFGS-B",lower=abind(rep(-Inf,q),0,-Inf,0,rep(-Inf,p)),upper=rep(Inf,q+p*(p+1)/2+p),control=list(fnscale=-1,maxit=5000,factr=1e3))
  P<-O$par
  estmax=P
  logvero=lmsnr(P,y,X)
  escore=escmsnr(P,y,X)
  informacao=mifmsnr(P,y,X)
  iter=as.numeric(O$counts[1])
  t1=Sys.time()
  tempo=as.numeric(t1-t0)
  object.out<-list(estmax=estmax,logvero=logvero,escore=escore,informacao=informacao,iter=iter,tempo=tempo)
} 

#Envelopes skew normal 

envelSN<-function(Y,X,beta,Sigma){
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
  for(i in 1:200){Xsim[i,]<-rchisq(n, p)}
  Xsim2<-apply(Xsim,1,sort)# aplica nas linhas, neste caso, coloca cada coluna em ordem crescente
  d21<-matrix(0,n,1)
  d22<-matrix(0,n,1)
  for(i in 1:n){
    d21[i]  <- quantile(Xsim2[i,],0.025)
    d22[i]  <- quantile(Xsim2[i,],0.975)}
  d2med <-apply(Xsim2,1,mean) # tira media de cada coluna; tira medias dos vetores menores...ate os maiores
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

envel0SN<-function(Y,X,beta,Sigma,lambda){
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
   Y0=rmsnr(n,beta,sqrtm(Sigma),lambda,X)
   Theta0=EMDMSNRM(Y0,X)$estmax
   d20<-matrix(0,n,1)
   for(j in 1:n){
    d20[j]=t(Y[j,]-X[[j]]%*%Theta0[1:q])%*%solve(xpnd(Theta0[(q+1):(q+p*(p+1)/2)])%*%xpnd(Theta0[(q+1):(q+p*(p+1)/2)]))%*%(Y[j,]-X[[j]]%*%Theta0[1:q])
   }
  Xsim[i,]<-d20}
  Xsim2<-apply(Xsim,1,sort)# aplica nas linhas, neste caso, coloca cada coluna em ordem crescente
  d21<-matrix(0,n,1)
  d22<-matrix(0,n,1)
  for(i in 1:n){
    d21[i]  <- quantile(Xsim2[i,],0.025)
    d22[i]  <- quantile(Xsim2[i,],0.975)}
  d2med <-apply(Xsim2,1,mean) # tira media de cada coluna; tira medias dos vetores menores...ate os maiores
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

#Envelopes normal 

envelN<-function(Y,X,beta,Sigma){
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
  for(i in 1:200){Xsim[i,]<-rchisq(n, p)}
  Xsim2<-apply(Xsim,1,sort)# aplica nas linhas, neste caso, coloca cada coluna em ordem crescente
  d21<-matrix(0,n,1)
  d22<-matrix(0,n,1)
  for(i in 1:n){
    d21[i]  <- quantile(Xsim2[i,],0.025)
    d22[i]  <- quantile(Xsim2[i,],0.975)}
  d2med <-apply(Xsim2,1,mean) # tira media de cada coluna; tira medias dos vetores menores...ate os maiores
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

envel0N<-function(Y,X,beta,Sigma){
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
  Y0<-matrix(0,n,p)
  for(i in 1:100){
   for(k in 1:n){
    Y0[k,]=mvrnorm(1,X[[k]]%*%beta,Sigma)
   }
   Theta0=EMDMSNRP(Y0,X)$estmax
   d20<-matrix(0,n,1)
   for(j in 1:n){
    d20[j]=t(Y[j,]-X[[j]]%*%Theta0[1:q])%*%solve(xpnd(Theta0[(q+1):(q+p*(p+1)/2)])%*%xpnd(Theta0[(q+1):(q+p*(p+1)/2)]))%*%(Y[j,]-X[[j]]%*%Theta0[1:q])
   }
  Xsim[i,]<-d20}
  Xsim2<-apply(Xsim,1,sort)# aplica nas linhas, neste caso, coloca cada coluna em ordem crescente
  d21<-matrix(0,n,1)
  d22<-matrix(0,n,1)
  for(i in 1:n){
    d21[i]  <- quantile(Xsim2[i,],0.025)
    d22[i]  <- quantile(Xsim2[i,],0.975)}
  d2med <-apply(Xsim2,1,mean) # tira media de cada coluna; tira medias dos vetores menores...ate os maiores
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





# LS estimators
EMLS0<-function(y,X)
{
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
  return(b)
  }

EMLS1<-function(y,X)
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
  P_0<-c(as.vector(b),vech(B))
  crit<-1
  iter<-0
  while(crit>=10^(-6))
  { 
    iter<-iter+1 
    b0<-matrix(0,q,q)
    b1<-matrix(0,q,1)
    for(i in 1:n){
      b0<-b0+t(X[[i]])%*%invS%*%X[[i]]
      b1<-b1+t(X[[i]])%*%invS%*%y[i,]
    }
    C<-solve(b0)
    b<-C%*%b1
    #d<-matrix(0,n,1)
    e<-matrix(0,n,p)
    for(i in 1:n){
      e[i,]<-y[i,]-X[[i]]%*%b}
    S<-1/(n-q)*t(e)%*%e
    invS<-solve(S)
    B<-sqrtm(S)
    P<-c(as.vector(b),vech(B))
    crit<-sqrt(sum((P-P_0)^2))
    P_0<-P
  }# o algoritmo EM
  t1=Sys.time()
  estmax=P
  tempo=as.numeric(t1-t0)
  object.out<-list(estmax=estmax,iter=iter,tempo=tempo)
}


# Normal model
EMN<-function(y,X)
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
  P_0<-c(as.vector(b),vech(B))
  crit<-1
  iter<-0
  while(crit>=10^(-6))
  { 
    iter<-iter+1 
    b0<-matrix(0,q,q)
    b1<-matrix(0,q,1)
    for(i in 1:n){
      b0<-b0+t(X[[i]])%*%invS%*%X[[i]]
      b1<-b1+t(X[[i]])%*%invS%*%y[i,]
    }
    C<-solve(b0)
    b<-C%*%b1
    #d<-matrix(0,n,1)
    e<-matrix(0,n,p)
    for(i in 1:n){
      e[i,]<-y[i,]-X[[i]]%*%b}
    S<-1/n*t(e)%*%e
    invS<-solve(S)
    B<-sqrtm(S)
    P<-c(as.vector(b),vech(B))
    crit<-sqrt(sum((P-P_0)^2))
    P_0<-P
  }# o algoritmo EM
  t1=Sys.time()
  estmax=P
  tempo=as.numeric(t1-t0)
  object.out<-list(estmax=estmax,iter=iter,tempo=tempo)
}
