#R_reduced_rank_envelope

library(Renvlp)
library(expm)
library(readr)

rrenv = function(X,Y,u=NULL,d=NULL,Beta=NULL){
  X = as.matrix(X)
  Y = as.matrix(Y)
  source("rrenv_choose_du.R") 
  dyn.load("R_reduced_rank_envelope.so")
  R_reduced_rank_envelope_given_u_d <- function(X,Y,u,d){
    X = as.matrix(X)
    Y = as.matrix(Y)
    u = as.integer(u)
    d = as.integer(d)
    p = dim(X)[1]
    r = dim(Y)[1]
    Gamma = env(t(X),t(Y),u)$Gamma
    Gamma0 = env(t(X),t(Y),u)$Gamma0
    if(u==nrow(Y)){Gamma0 = matrix(0)}
    R_beta = rep(0.0, length=nrow(Y)*nrow(X))
    R_fitvalues = rep(0.0, length=nrow(Y)*ncol(Y))
    R_estiduals = rep(0.0, length=nrow(Y)*ncol(Y))
    R_estidualvariances = rep(0.0, length=nrow(Y)*nrow(Y))
    est = .C("R_reduced_rank_envelope",
             X,nrow(X),ncol(X), 
             Y,nrow(Y),ncol(Y),
             u,d,
             Gamma,nrow(Gamma),ncol(Gamma),
             Gamma0,nrow(Gamma0),ncol(Gamma0),
             R_beta = R_beta,
             R_fitvalues = R_fitvalues,
             R_estiduals = R_estiduals,
             R_estidualvariances = R_estidualvariances)
    beta=matrix(est$R_beta,nrow=nrow(Y),ncol=nrow(X))
    fitvalues=matrix(est$R_fitvalues,nrow=nrow(Y),ncol=ncol(Y))
    estiduals=matrix(est$R_estiduals,nrow=nrow(Y),ncol=ncol(Y))
    estidual_variance=matrix(est$R_estidualvariances,nrow=nrow(Y),ncol=nrow(Y))
    if((u>d)&(u<=r)&(d<=min(p,r))){
      err.table = matrix(0, ncol=3, nrow=1)
      NR = (p+u-d)*d + ((r*(r+1))/2)
      Lud = 0
      Lud = sapply(1:dim(Y)[2],function(i){dmvnorm(Y[,i],fitvalues[,i],estidual_variance,log=T)})
      Lud = sum(unlist(Lud))
      err = 2*NR - 2*Lud
      err.table[1,1]=d; err.table[1,2]=u; err.table[1,3]=err
      err.table = as.data.frame(err.table)
      names(err.table) = c('d','u','AIC')
    return(list(beta=beta,fit_values=fitvalues,residuals=estiduals,
                  residual_variance=estidual_variance,error_table=err.table))
      }
    else{
      print("Error: Invalid u and d!")
    }
  }
  if((is.null(u))&(is.null(d))){
    choose_du =  rrenv.choose_du(X=X, Y=Y, Beta=Beta)
    u = choose_du$best_combination$u
    d = choose_du$best_combination$d
    print(choose_du)
  }
  else{
    u = as.integer(u)
    print(paste0("u=",u))
    d = as.integer(d)
    print(paste0("d=",d))
  }
  est =  R_reduced_rank_envelope_given_u_d(X=X, Y=Y,u=u,d=d)
  return(est)

}

#example 
Y<-as.matrix(read.table("Y.txt"))#Y file from data_simulation.cpp
X<-as.matrix(read.table("X.txt"))#X file from data_simulation.cpp
est <- rrenv(X,Y,20,10,NULL)
#returns
beta<-est$beta
fit_values<-est$fit_values
residuals<-est$residuals
residual_variance<-est$residual_variance
error_table<-est$error_table