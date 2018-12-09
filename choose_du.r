library(Renvlp)
library(expm)
library(readr)
library(data.table)
library(ggplot2)
library(plyr)
library(mvtnorm)

rrenv.choose_du = function(X,Y, Beta=NULL){
  
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
  
  
  p = dim(X)[1]
  r = dim(Y)[1]
  max.u = r
  max.d = min(p,r)
  c = max.u-max.d
  rows = ((max.u*(max.u+1))/2) - ((c*(c+1))/2) - 1
  err.table = matrix(0, ncol=3, nrow=rows)
  cnt = 0
  
  if(is.null(Beta)){
    for(u in 2:max.u){
      d=0
      while(d<max.d & d<u){
        d = d+1
        cnt = cnt+1
        est = R_reduced_rank_envelope_given_u_d(X,Y,u,d)
        est.Betas = est$beta
        est.Ys = est$fitvalues
        est.residual_variance = est$residual_variance
        NR = (p+u-d)*d + ((r*(r+1))/2)
        Lud = 0
#        for(i in 1:r){
#          L = dmvnorm(x = t(Y)[i,], mean=t(est.Ys)[i,], sigma=(est.residual_variance), log=T)
#          Lub = Lud + L
#        }
        Lud = sapply(1:dim(Y)[2],function(i){dmvnorm(Y[,i],est.Ys[,i],est.residual_variance,log=T)})
        Lud = sum(unlist(Lud))
        err = 2*NR - 2*Lud
        err.table[cnt,1]=d; err.table[cnt,2]=u; err.table[cnt,3]=err
      }
    }
    err.table = as.data.frame(err.table)
    names(err.table) = c('d','u','AIC')
    m <- ddply(err.table, .(d), transform, rescale = scales::rescale(AIC))
    heatmap <- ggplot(m, aes(d, u)) + theme_gray() + geom_tile(aes(fill = rescale), colour = "white") + scale_fill_gradient(low = "white",  high = "steelblue")
    return(list(error_table = err.table, best_combination = err.table[which.min(err.table[,3]),], heatmap=heatmap))
  } 
  
  
  else if(all(dim(Beta)==c(r,p))) {
    for(u in 2:max.u){
      d=0
      while(d<max.d & d<u){
        d = d+1
        cnt = cnt+1
        NR = (p+u-d)*d + ((r*(r+1))/2)
        est = R_reduced_rank_envelope_given_u_d(X,Y,u,d)
        est.Betas = est$beta
        err = norm(est.Betas-Beta,type = "F")
        err.table[cnt,1]=d; err.table[cnt,2]=u; err.table[cnt,3]=err
      }
    }
    err.table = as.data.frame(err.table)
    names(err.table) = c('d','u','Error_in_Beta')
    m <- ddply(err.table, .(d), transform, rescale = scales::rescale(Error_in_Beta))
    heatmap <- ggplot(m, aes(d, u)) + theme_gray() + geom_tile(aes(fill = rescale), colour = "white") + scale_fill_gradient(low = "white",  high = "steelblue")
    return(list(error_table = err.table, best_combination = err.table[which.min(err.table[,3]),], heatmap=heatmap))
  }
  else{
    print("Error: Invalid Beta Matrix!")}
}


X = as.matrix(read.table('/Users/shnlee/Desktop/615final/X.txt'))
Y = as.matrix(read.table('/Users/shnlee/Desktop/615final/Y.txt'))
Beta = as.matrix(read.table('/Users/shnlee/Desktop/615final/Beta.txt'))

a = rrenv.choose_du(X=X, Y=Y, Beta=Beta)
b = rrenv.choose_du(X=X, Y=Y, Beta=NULL)

a$error_table
a$best_combination
a$heatmap
