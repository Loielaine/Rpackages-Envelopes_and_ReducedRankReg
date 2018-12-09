#R_reduced_rank_envelope
#liyixi
#2018.11.30
library(Renvlp)
library(expm)
library(readr)
setwd("/Users/Loielaine/Desktop/Good_Good_Study/2018Fall/BIOSTAT615_Statistic_Computing/project/Rpackages-Envelopes_and_ReducedRankReg/")
dyn.load("R_reduced_rank_envelope.so")
R_reduced_rank_envelope_given_u_d <- function(X,Y,u,d){
    X = as.matrix(X)
    Y = as.matrix(Y)
    u = as.integer(u)
    d = as.integer(d)
    Gamma = env(t(X),t(Y),u)$Gamma
    Gamma0 = env(t(X),t(Y),u)$Gamma0
    if(u==nrow(Y)){Gamma0 = matrix(0)}
    R_beta = rep(0.0, length=nrow(Y)*nrow(X))
    R_fitvalues = rep(0.0, length=nrow(Y)*ncol(Y))
    R_residuals = rep(0.0, length=nrow(Y)*ncol(Y))
    R_residualvariances = rep(0.0, length=nrow(Y)*nrow(Y))
   res = .C("R_reduced_rank_envelope",
            X,nrow(X),ncol(X), 
            Y,nrow(Y),ncol(Y),
            u,d,
            Gamma,nrow(Gamma),ncol(Gamma),
            Gamma0,nrow(Gamma0),ncol(Gamma0),
            R_beta = R_beta,
            R_fitvalues = R_fitvalues,
            R_residuals = R_residuals,
            R_residualvariances = R_residualvariances)
   return(list(beta=matrix(res$R_beta,nrow=nrow(Y),ncol=nrow(X)),
               fitvalues=matrix(res$R_fitvalues,nrow=nrow(Y),ncol=ncol(Y)),
               residuals=matrix(res$R_residuals,nrow=nrow(Y),ncol=ncol(Y)),
               residual_variance=matrix(res$R_residualvariances,nrow=nrow(Y),ncol=nrow(Y))))
}
#example 
Y<-as.matrix(read.table("Y.txt"))#Y file from data_simulation.cpp
X<-as.matrix(read.table("X.txt"))#X file from data_simulation.cpp
res <- R_reduced_rank_envelope_given_u_d(X,Y,20,10)
#returns
beta<-res$beta
fitvalues<-res$fitvalues
residuals<-res$residuals
residal_variance<-res$residual_variance


