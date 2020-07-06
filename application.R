
source("snmreg.R")
source("stmreg.R")
source("ssmreg.R")
source("scmreg.R")

source("infgloreg.R")
source("inflocreg.R")


#Manipulacao dos dados 

dados=read.csv("AtlasBrasil_Consulta.csv")
fix(dados)
attach(dados)


esperanca_vida=Esperan√.a.de.vida.ao.nascer.2010
renda_per_capita=Renda.per.capita.2010


taxa_atividade=Taxa.de.atividade...18.anos.ou.mais.2010
taxa_fecundidade=Taxa.de.fecundidade.total.2010
razao_dependencia=Raz√.o.de.depend√™ncia.2010
taxa_envelhecimento=Taxa.de.envelhecimento.2010
taxa_analfabetismo=Taxa.de.analfabetismo...18.anos.ou.mais.2010 


municipios=Espacialidades
n=length(municipios)
p=2

y=cbind(esperanca_vida,renda_per_capita)
y #matriz de respostas

ThetaN0=EMDMSNP(y)
ThetaSN0=EMDMSNM(y)
ThetaSS0=EMDMSSM(y)
ThetaSC0=EMDMSCM(y)
ThetaST0=EMDMSTM(y)

q=5
X1<-list()
for(i in 1:n){
    X1[[i]]=matrix(0,p,p*q)
    for(j in 1:p){
     X1[[i]][j,((j-1)*q+1):(j*q)]=cbind(1,taxa_atividade[i],taxa_fecundidade[i],razao_dependencia[i],taxa_analfabetismo[i])
    }    
}
X1 #matriz de planejamento 1 com as 4 explicativas para cada resposta

ThetaN1=EMDMSNRP(y,X1)
ThetaSN1=EMDMSNRM(y,X1)
ThetaSS1=EMDMSSRM(y,X1)
ThetaSC1=EMDMSCRM(y,X1)
ThetaST1=EMDMSTRM(y,X1)

ThetaN1$estmax
ThetaSN1$estmax
ThetaSS1$estmax
ThetaSC1$estmax
ThetaST1$estmax

ThetaN1$estmax[1:(p*q+p*(p+1)/2)]/sqrt(diag(-solve(ThetaN1$informacao[c(1:(p*q+p*(p+1)/2)),c(1:(p*q+p*(p+1)/2))])))
ThetaSN1$estmax/sqrt(diag(-solve(ThetaSN1$informacao)))
ThetaSS1$estmax/sqrt(diag(-solve(ThetaSS1$informacao)))
ThetaSC1$estmax/sqrt(diag(-solve(ThetaSC1$informacao)))
ThetaST1$estmax/sqrt(diag(-solve(ThetaST1$informacao))) #retirar taxa_fecundidade

q1=4
q2=4

X2<-list()
for(i in 1:n){
    X2[[i]]=matrix(0,p,q1+q2)
    X2[[i]][1,1:q1]=cbind(1,taxa_atividade[i],razao_dependencia[i],taxa_analfabetismo[i])
    X2[[i]][2,(q1+1):(q1+q2)]=cbind(1,taxa_atividade[i],razao_dependencia[i],taxa_analfabetismo[i])  
}
X2 #matriz de planejamento 2 retirando a segunda explicativa

ThetaN2=EMDMSNRP(y,X2)
ThetaSN2=EMDMSNRM(y,X2)
ThetaSS2=EMDMSSRM(y,X2)
ThetaSC2=EMDMSCRM(y,X2)
ThetaST2=EMDMSTRM(y,X2)

ThetaN2$estmax
ThetaSN2$estmax
ThetaSS2$estmax
ThetaSC2$estmax
ThetaST2$estmax

ThetaN2$estmax[1:(q1+q2+p*(p+1)/2)]/sqrt(diag(-solve(ThetaN2$informacao[c(1:(q1+q2+p*(p+1)/2)),c(1:(q1+q2+p*(p+1)/2))])))
ThetaSN2$estmax/sqrt(diag(-solve(ThetaSN2$informacao)))
ThetaSS2$estmax/sqrt(diag(-solve(ThetaSS2$informacao)))
ThetaSC2$estmax/sqrt(diag(-solve(ThetaSC2$informacao)))
ThetaST2$estmax/sqrt(diag(-solve(ThetaST2$informacao)))


#Analise grafica das estimativas 

X=X2
vecy=c(y[,1],y[,2])
vecy
mX<-X[[1]][1,2:4]
for(i in 2:n){
 mX<-rbind(mX,X[[i]][1,2:4])
}
mX

#N

bn=ThetaN2$estmax[1:(q1+q2)]
an=ThetaN2$estmax[(q1+q2+1):(q1+q2+p*(p+1)/2)]
Bn=xpnd(an)
hn=ThetaN2$estmax[(q1+q2+p*(p+1)/2+1):(q1+q2+p*(p+1)/2+p)]
un=ThetaN2$u
Cn=ThetaN2$C
Wn<-list()
for(i in 1:n){
  Wn[[i]]=solve(Bn)%*%(un[i]*diag(p)+hn%*%t(hn))%*%solve(Bn)%*%X[[i]]
}
matWn<-Wn[[1]]
for(i in 2:n){
 matWn<-rbind(matWn,Wn[[i]])
}
Hn=matWn%*%Cn%*%t(matWn)
plot(diag(Hn)) 
abline(h=2*sum(diag(Hn))/n,col=2)
contn<-0
for(i in 1:(2*n)){
if(diag(Hn)[i]>2*sum(diag(Hn))/n){contn<-contn+1}
else{contn<-contn}
} 
contn #alavanca n

En=matrix(0,p,1)
Vn=Bn%*%Bn #momentos n

yn<-matrix(0,n,p)
for(i in 1:n){
 yn[i,]=X[[i]]%*%bn+En
}
en=y-yn

plot(en[,1])
abline(h=0)
summary(en[,1])
sd(en[,1])
skewness(en[,1]) 
boxplot(en[,1])  

plot(en[,2])
abline(h=0)
summary(en[,2])
sd(en[,2])
skewness(en[,2]) 
boxplot(en[,2]) #residuos n


#SN

bsn=ThetaSN2$estmax[1:(q1+q2)]
asn=ThetaSN2$estmax[(q1+q2+1):(q1+q2+p*(p+1)/2)]
Bsn=xpnd(asn)
hsn=ThetaSN2$estmax[(q1+q2+p*(p+1)/2+1):(q1+q2+p*(p+1)/2+p)]
usn=ThetaSN2$u
Csn=ThetaSN2$C
Wsn<-list()
for(i in 1:n){
  Wsn[[i]]=solve(Bsn)%*%(usn[i]*diag(p)+hsn%*%t(hsn))%*%solve(Bsn)%*%X[[i]]
}
matWsn<-Wsn[[1]]
for(i in 2:n){
 matWsn<-rbind(matWsn,Wsn[[i]])
}
Hsn=matWsn%*%Csn%*%t(matWsn)
plot(diag(Hsn)) 
abline(h=2*sum(diag(Hsn))/n,col=2)
contsn<-0
for(i in 1:(2*n)){
if(diag(Hsn)[i]>2*sum(diag(Hsn))/n){contsn<-contsn+1}
else{contsn<-contsn}
} 
contsn #alavanca sn

Esn=sqrt(2/pi)*Bsn%*%hsn/sqrt(1+sum(hsn*hsn))
Euhsn=1/sqrt(1+sum(hsn*hsn))
Vsn=Bsn%*%(diag(p)-2/pi*Euhsn^2*hsn%*%t(hsn))%*%Bsn #momentos sn

ysn<-matrix(0,n,p)
for(i in 1:n){
 ysn[i,]=X[[i]]%*%bsn+Esn
}
esn=y-ysn

plot(esn[,1])
abline(h=0)
summary(esn[,1])
sd(esn[,1])
skewness(esn[,1]) 
boxplot(esn[,1]) #ver boxplot assimetrico 

plot(esn[,2])
abline(h=0)
summary(esn[,2])
sd(esn[,2])
skewness(esn[,2]) 
boxplot(esn[,2]) #residuos sn


#SC

bsc=ThetaSC2$estmax[1:(q1+q2)]
asc=ThetaSC2$estmax[(q1+q2+1):(q1+q2+p*(p+1)/2)]
Bsc=xpnd(asc)
hsc=ThetaSC2$estmax[(q1+q2+p*(p+1)/2+1):(q1+q2+p*(p+1)/2+p)]
vsc=ThetaSC2$estmax[q1+q2+p*(p+1)/2+p+1]
gsc=ThetaSC2$estmax[q1+q2+p*(p+1)/2+p+2]
usc=ThetaSC2$u
Csc=ThetaSC2$C

Wsc<-list()
for(i in 1:n){
  Wsc[[i]]=solve(Bsc)%*%(usc[i]*diag(p)+hsc%*%t(hsc))%*%solve(Bsc)%*%X[[i]]
}
matWsc<-Wsc[[1]]
for(i in 2:n){
 matWsc<-rbind(matWsc,Wsc[[i]])
}
Hsc=matWsn%*%Csc%*%t(matWsc)
plot(diag(Hsc)) 
abline(h=2*sum(diag(Hsc))/n,col=2)
contsc<-0
for(i in 1:(2*n)){
if(diag(Hsc)[i]>2*sum(diag(Hsc))/n){contsc<-contsc+1}
else{contsc<-contsc}
} 
contsc #alavanca sc

Esc=sqrt(2/pi)*(vsc/sqrt(vsc*(vsc+sum(hsc*hsc)))+(1-vsc)/sqrt(1+sum(hsc*hsc)))*Bsc%*%hsc
Euhsc=vsc/sqrt(vsc*(vsc+sum(hsc*hsc)))+(1-vsc)/sqrt(1+sum(hsc*hsc))
Vsc=Bsc%*%((vsc/gsc+1-vsc)*diag(p)-2/pi*Euhsc^2*hsc%*%t(hsc))%*%Bsc #momentos sc

ysc<-matrix(0,n,p)
for(i in 1:n){
 ysc[i,]=X[[i]]%*%bsc+Esc
}
esc=y-ysc

plot(esc[,1])
abline(h=0)
summary(esc[,1])
sd(esc[,1])
skewness(esc[,1]) 
boxplot(esc[,1]) #ver boxplot assimetrico  

plot(esc[,2])
abline(h=0)
summary(esc[,2])
sd(esc[,2])
skewness(esc[,2]) 
boxplot(esc[,2]) #residuos sc


#SS

bss=ThetaSS2$estmax[1:(q1+q2)]
ass=ThetaSS2$estmax[(q1+q2+1):(q1+q2+p*(p+1)/2)]
Bss=xpnd(ass)
hss=ThetaSS2$estmax[(q1+q2+p*(p+1)/2+1):(q1+q2+p*(p+1)/2+p)]
vss=ThetaSS2$estmax[q1+q2+p*(p+1)/2+p+1]
uss=ThetaSS2$u
Css=ThetaSS2$C

Wss<-list()
for(i in 1:n){
  Wss[[i]]=solve(Bss)%*%(uss[i]*diag(p)+hss%*%t(hss))%*%solve(Bss)%*%X[[i]]
}
matWss<-Wss[[1]]
for(i in 2:n){
 matWss<-rbind(matWss,Wss[[i]])
}
Hss=matWss%*%Css%*%t(matWss)
plot(diag(Hss)) 
abline(h=2*sum(diag(Hss))/n,col=2)
contss<-0
for(i in 1:(2*n)){
if(diag(Hss)[i]>2*sum(diag(Hss))/n){contss<-contss+1}
else{contss<-contss}
} 
contss #alavanca ss

fss<-function(u){
fss=u^(vss-3/2)/sqrt(u+sum(hss*hss))
return(fss)
}
Euhss=vss*integrate(fss,0,1)$value
Ess=sqrt(2/pi)*Euhss*Bss%*%hss
Vss=Bss%*%((vss/(vss-1))*diag(p)-2/pi*Euhss^2*hss%*%t(hss))%*%Bss #momentos ss

yss<-matrix(0,n,p)
for(i in 1:n){
 yss[i,]=X[[i]]%*%bss+Ess
}
ess=y-yss

plot(ess[,1])
abline(h=0)
summary(ess[,1])
sd(ess[,1])
skewness(ess[,1]) 
boxplot(ess[,1]) 

plot(ess[,2])
abline(h=0)
summary(ess[,2])
sd(ess[,2])
skewness(ess[,2]) 
boxplot(ess[,2]) #residuos ss


#ST

bst=ThetaST2$estmax[1:(q1+q2)]
ast=ThetaST2$estmax[(q1+q2+1):(q1+q2+p*(p+1)/2)]
Bst=xpnd(ast)
hst=ThetaST2$estmax[(q1+q2+p*(p+1)/2+1):(q1+q2+p*(p+1)/2+p)]
vst=ThetaST2$estmax[q1+q2+p*(p+1)/2+p+1]
ust=ThetaST2$u
Cst=ThetaST2$C
Wst<-list()
for(i in 1:n){
  Wst[[i]]=solve(Bst)%*%(ust[i]*diag(p)+hst%*%t(hst))%*%solve(Bst)%*%X[[i]]
}
matWst<-Wst[[1]]
for(i in 2:n){
 matWst<-rbind(matWst,Wst[[i]])
}
Hst=matWst%*%Cst%*%t(matWst)
plot(diag(Hst)) 
abline(h=2*sum(diag(Hst))/n,col=2)
contst<-0
for(i in 1:(2*n)){
if(diag(Hst)[i]>2*sum(diag(Hst))/n){contst<-contst+1}
else{contst<-contst}
} 
contst #alavanca st

fst<-function(u){
fst=u^((vst-3)/2)*exp(-u*vst/2)/sqrt(u+sum(hst*hst))
return(fst)
}
Euhst=(vst/2)^(vst/2)/gamma(vst/2)*integrate(fst,0,Inf)$value
Est=sqrt(2/pi)*Euhst*Bst%*%hst
Vst=Bst%*%((vst/(vst-2))*diag(p)-2/pi*Euhst^2*hst%*%t(hst))%*%Bst #momentos st

yst<-matrix(0,n,p)
for(i in 1:n){
 yst[i,]=X[[i]]%*%bst+Est
}
est=y-yst

plot(est[,1])
abline(h=0)
summary(est[,1])
sd(est[,1])
skewness(est[,1]) 
boxplot(est[,1]) 

plot(est[,2])
abline(h=0)
summary(est[,2])
sd(est[,2])
skewness(est[,2]) 
boxplot(est[,2]) #residuos st



#Selecao do melhor modelo

z=seq(67.5,80,length=250)
w=seq(150,1800,length=250)
MX<-list()
MX[[1]]<-cbind(matrix(rbind(c(1,colMeans(mX)),rep(0,q1)),p,q1),matrix(rbind(rep(0,q2),c(1,colMeans(mX))),p,q2))

f=matrix(0,250,250)
for(i in 1:250){
 for(j in 1:250){
  #f[i,j]=dmsnr(ThetaN2$estmax,cbind(z[i],w[j]),MX)
  #f[i,j]=dmsnr(ThetaSN2$estmax,cbind(z[i],w[j]),MX)
  f[i,j]=dmssr(ThetaSS2$estmax,cbind(z[i],w[j]),MX)
  #f[i,j]=dmstr(ThetaST2$estmax,cbind(z[i],w[j]),MX)
  #f[i,j]=dmscr(ThetaSC2$estmax,cbind(z[i],w[j]),MX)
 }  
}
color<-colorRampPalette(c("black","white"))
#barplot(1:50,col=color(50))
contour(z,w,f,nlevels=1000,drawlabels=F,col=color(50)[25])
par(new=T)
plot(y,xlab='Y1',ylab='Y2',xlim=c(min(z),max(z)),ylim=c(min(w),max(w)),main='Slash normal assimÈtrica',pch=19)
 

MODELO=c('N','SN','SS','ST','SC')
LOGVERO=c(ThetaN2$logvero,ThetaSN2$logvero,ThetaSS2$logvero,ThetaST2$logvero,ThetaSC2$logvero)
AICn=-2*ThetaN2$logvero+2*length(ThetaN2$estmax[1:(q1+q2+p*(p+1)/2)])
AICsn=-2*ThetaSN2$logvero+2*length(ThetaSN2$estmax)
AICss=-2*ThetaSS2$logvero+2*length(ThetaSS2$estmax)
AICst=-2*ThetaST2$logvero+2*length(ThetaST2$estmax)
AICsc=-2*ThetaSC2$logvero+2*length(ThetaSC2$estmax)
AIC=c(AICn,AICsn,AICss,AICst,AICsc)
BICn=-2*ThetaN2$logvero+log(n)*length(ThetaN2$estmax[1:(q1+q2+p*(p+1)/2)])
BICsn=-2*ThetaSN2$logvero+log(n)*length(ThetaSN2$estmax)
BICss=-2*ThetaSS2$logvero+log(n)*length(ThetaSS2$estmax)
BICst=-2*ThetaST2$logvero+log(n)*length(ThetaST2$estmax)
BICsc=-2*ThetaSC2$logvero+log(n)*length(ThetaSC2$estmax)
BIC=c(BICn,BICsn,BICss,BICst,BICsc)
rbind(MODELO,LOGVERO,AIC,BIC)


#Interpretacao do modelo

taxa_atividade0=(taxa_atividade-mean(taxa_atividade))/sd(taxa_atividade)
razao_dependencia0=(razao_dependencia-mean(razao_dependencia))/sd(razao_dependencia)
taxa_analfabetismo0=(taxa_analfabetismo-mean(taxa_analfabetismo))/sd(taxa_analfabetismo)

X0<-list()
for(i in 1:n){
    X0[[i]]=matrix(0,p,q1+q2)
    X0[[i]][1,1:q1]=cbind(1,taxa_atividade0[i],razao_dependencia0[i],taxa_analfabetismo0[i])
    X0[[i]][2,(q1+1):(q1+q2)]=cbind(1,taxa_atividade0[i],razao_dependencia0[i],taxa_analfabetismo0[i])  
}
X0 #matriz de planejamento 3 com variaveis explicativas padronizadas

ThetaN0=EMDMSNRP(y,X0)
ThetaSN0=EMDMSNRM(y,X0)
ThetaSS0=EMDMSSRM(y,X0)
ThetaSC0=EMDMSCRM(y,X0)
ThetaST0=EMDMSTRM(y,X0)

sqrt(diag(solve(-ThetaN0$informacao)))
sqrt(diag(solve(-ThetaSN0$informacao)))
sqrt(diag(solve(-ThetaSS0$informacao)))
sqrt(diag(solve(-ThetaSC0$informacao)))
sqrt(diag(solve(-ThetaST0$informacao)))

ThetaN0$estmax
ThetaSN0$estmax
ThetaSS0$estmax
ThetaSC0$estmax
ThetaST0$estmax



#Influencia global


#Estimacoes com a retiradad de cada observacao

cl<-makeCluster(detectCores())
registerDoParallel(cl)

Theta_n<-NULL
Theta_n<-foreach(j=1:n,.combine=rbind)%dopar%{
  library(MASS)
  library(expm)
  library(moments)
  library(abind)
  library(mvtnorm)
  require(MCMCpack)
  Theta_n<-EMDMSNRP(y[-j,],X[-j])$estmax
}
Theta_sn<-NULL
Theta_sn<-foreach(j=1:n,.combine=rbind)%dopar%{
  library(MASS)
  library(expm)
  library(moments)
  library(abind)
  library(mvtnorm)
  require(MCMCpack)
  Theta_sn<-EMDMSNRM(y[-j,],X[-j])$estmax
}
Theta_st<-NULL
Theta_st<-foreach(j=1:n,.combine=rbind)%dopar%{
  library(MASS)
  library(expm)
  library(moments)
  library(abind)
  library(mvtnorm)
  require(MCMCpack)
  Theta_st<-EMDMSTRM(y[-j,],X[-j])$estmax
}
Theta_sc<-NULL
Theta_sc<-foreach(j=1:n,.combine=rbind)%dopar%{
  library(MASS)
  library(expm)
  library(moments)
  library(abind)
  library(mvtnorm)
  require(MCMCpack)
  Theta_sc<-EMDMSCRM(y[-j,],X[-j])$estmax
}
Theta_ss<-NULL
Theta_ss<-foreach(j=1:n,.combine=rbind)%dopar%{
  library(MASS)
  library(expm)
  library(moments)
  library(abind)
  library(mvtnorm)
  require(MCMCpack)
  Theta_ss<-EMDMSSRM(y[-j,],X[-j])$estmax
}
stopCluster(cl)


par(mfrow=c(1,3))


#N

maxgd=max(GD_n,GD_sn,GD_ss)
#ThetaN=ThetaN2$estmax
GD_n<-rep(0,n)
for(i in 1:n){
  GD_n[i]<--t(ThetaN-Thetan[i,])%*%hesqmnr(ThetaN,y,X)%*%(ThetaN-Thetan[i,])    
}
plot(GD_n, xlab='Index', ylab='GD', font.lab=2,ps=0.1,main='Normal model', ylim=c(0,maxgd)) 
identify(GD_n,n=4)
#DISTANCIA DE COOK GENERALIZADA N


#SN

#ThetaSN=ThetaSN2$estmax
GD_sn<-rep(0,n)
for(i in 1:n){
  GD_sn[i]<--t(ThetaSN-Thetasn[i,])%*%hesqmsnr(ThetaSN,y,X)%*%(ThetaSN-Thetasn[i,])
}
plot(GD_sn, xlab='Index', ylab='GD', font.lab=2,ps=0.1,main='SN model', ylim=c(0,maxgd)) 
identify(GD_sn,n=4)
#DISTANCIA DE COOK GENERALIZADA SN


#ST


#SS

#ThetaSS=ThetaSS2$estmax
GD_ss<-rep(0,n)
for(i in 1:n){
  GD_ss[i]<--t(ThetaSS-Thetass[i,])%*%hesqmssr(ThetaSS,y,X)%*%(ThetaSS-Thetass[i,])
}
plot(GD_ss,xlab='Index', ylab='GD', font.lab=2,ps=0.1,main='SSN model', ylim=c(0,maxgd))
identify(GD_ss,n=2)
#DISTANCIA DE COOK GENERALIZADA SS


#SC



#Influencia local

par(mfrow=c(1,3))


#N

#r=length(ThetaN)
ThetaN=ThetaN2$estmax[1:11]
Dw0mnr=Dw0ponmnr(ThetaN,y,X2) #ponderacao de casos
Dw0mnr=Dw0resmnr(ThetaN,y,X2,c(0,1)) #perturbacao na resposta
Dw0mnr=Dw0expmnr(ThetaN,y,X2,c(2,3)) #perturbacao na explicativa 
#Dw0mnr=Dw0escmnr(ThetaN,y,X) #perturbacao na escala

HQw0n=-2*t(Dw0mnr)%*%solve(hesqmnr(ThetaN,y,X2))%*%Dw0mnr
plot(Re(eigen(HQw0n)$values)) #teste de direcoes dadas pela matriz de curvatura

avlwn=rep(0,n)
for(k in 1:r){
 avlwn[k]=eigen(HQw0n)$values[k]/sum(eigen(HQw0n)$values[1:r])
}
print(Re(avlwn[1:r])) #autovalores positivos ponderados

avtwn=eigen(HQw0n)$vectors[,1:r] 
plot(Re(avtwn[,1])) #autovetores associados

M0n<-rep(0,n)
for(j in 1:n){
  M0nj=0   
  for(k in 1:r){  
    M0nj=M0nj+Re(avlwn)[k]*Re(avtwn)[j,k]^2  
  }
  M0n[j]=M0nj  
}  
mean(M0n)

max3=max(M0n,M0sn,M0ss) # rodar as 3 para pegar o maximo

plot(M0n,xlab='Index',ylab='M(0)',font.lab=2,ps=0.1,main='Normal model',ylim=c(0,max3))
abline(h=mean(M0n)+4*sd(M0n),col='gray48',lty=2,lwd=2)
identify(M0n,n=2) #MEDIDA M0 Normal


#SN

#r=length(ThetaSN)
ThetaSN=ThetaSN2$estmax
Dw0msnr=Dw0ponmsnr(ThetaSN,y,X2) #ponderacao de casos
Dw0msnr=Dw0resmsnr(ThetaSN,y,X2,c(0,1)) #perturbacao na resposta
Dw0msnr=Dw0expmsnr(ThetaSN,y,X2,c(2,3)) #perturbacao na explicativa 
#Dw0msnr=Dw0escmsnr(ThetaSN,y,X) #perturbacao na escala

HQw0sn=-2*t(Dw0msnr)%*%solve(hesqmsnr(ThetaSN,y,X2))%*%Dw0msnr
plot(Re(eigen(HQw0sn)$values)) #teste de direcoes dadas pela matriz de curvatura

avlwsn=rep(0,n)
for(k in 1:r){
 avlwsn[k]=eigen(HQw0sn)$values[k]/sum(eigen(HQw0sn)$values[1:r])
}
print(Re(avlwsn[1:r])) #autovalores positivos ponderados

avtwsn=eigen(HQw0sn)$vectors[,1:r] 
plot(Re(avtwsn[,1])) #autovetores associados

M0sn<-rep(0,n)
for(j in 1:n){
  M0snj=0   
  for(k in 1:r){  
    M0snj=M0snj+Re(avlwsn)[k]*Re(avtwsn)[j,k]^2  
  }
  M0sn[j]=M0snj  
}  
mean(M0sn)

plot(M0sn,xlab='Index',ylab='M(0)',font.lab=2,ps=0.1,main='SN model',ylim=c(0,max3))
abline(h=mean(M0sn)+4*sd(M0sn),col='gray48',lty=2,lwd=2)
identify(M0sn,n=2) #MEDIDA M0 SN


#SS

#r=length(ThetaSS)
ThetaSS=ThetaSS2$estmax
Dw0mssr=Dw0ponmssr(ThetaSS,y,X2) #ponderacao de casos
Dw0mssr=Dw0resmssr(ThetaSS,y,X2,c(0,1)) #perturbacao na resposta
Dw0mssr=Dw0expmssr(ThetaSS,y,X2,c(2,3)) #perturbacao na explicativa 
#Dw0mssr=Dw0escmssr(ThetaSS,y,X) #perturbacao na escala

HQw0ss=-2*t(Dw0mssr)%*%solve(hesqmssr(ThetaSS,y,X2))%*%Dw0mssr
plot(Re(eigen(HQw0ss)$values)) #teste de direcoes dadas pela matriz de curvatura

avlwss=rep(0,n)
for(k in 1:r){
 avlwss[k]=eigen(HQw0ss)$values[k]/sum(eigen(HQw0ss)$values[1:r])
}
print(Re(avlwss[1:r])) #autovalores positivos ponderados

avtwss=eigen(HQw0ss)$vectors[,1:r] 
plot(Re(avtwss[,1])) #autovetores associados

M0ss<-rep(0,n)
for(j in 1:n){
  M0ssj=0   
  for(k in 1:r){  
    M0ssj=M0ssj+Re(avlwss)[k]*Re(avtwss)[j,k]^2  
  }
  M0ss[j]=M0ssj  
}  
mean(M0ss)

plot(M0ss,xlab='Index',ylab='M(0)',font.lab=2,ps=0.1,main='SSN model',ylim=c(0,max3))
abline(h=mean(M0ss)+5*sd(M0ss),col='gray48',lty=2,lwd=2)
identify(M0ss,n=2) #MEDIDA M0 SLASH



