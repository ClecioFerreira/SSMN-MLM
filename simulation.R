source('stmreg.R')
source('snmreg.R')
source('scmreg.R')
source('ssmreg.R')

library(MASS)
library(expm)
library(moments)
library(abind)
library(mvtnorm)
require(MCMCpack)

################### Simulation
p=2
beta=c(40,-30,15)
rSigma=matrix(c(2.5,-1,-1,1.5),p,p) 
Sigma=rSigma%*%rSigma
lambda=c(2,-1)
nu=3.5

n=100


X<-list()
for (i in 1:n){
    aux<-matrix(1,2,3)
    aux[,2]<-runif(2,1,10)
    aux[,3]<-runif(2,-5,0)
    X[[i]]<- aux    
    }
mu<-matrix(0,n,p)
for (i in 1:n) mu[i,]=X[[i]]%*%beta

# StN model
erro=rmstr(n,rep(0,lb),rSigma,lambda,nu,X)
y=mu+erro
thetaECME=EMDMSTRP(y,X)# ECME
thetaECM=EMDMSTRM(y,X)# ECM

