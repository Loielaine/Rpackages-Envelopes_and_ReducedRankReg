#R_reduced_rank_envelope
#check packages
packages<-c("Renvlp","expm","data.table","ggplot2","plyr","mvtnorm","rstudioapi","MASS")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
#change path
#script.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#setwd(script.dir)
#setwd(paste0(script.dir,"/../src/"))
setwd("/Users/Loielaine/Desktop/Good_Good_Study/2018Fall/BIOSTAT615_Statistic_Computing/project/Rpackage_rrenv/rrenv/src/")
#functions
rrenv.choose_du = function(X,Y, Beta=NULL){
  
  dyn.load("rrenv.so")
  rrenv.given_du <- function(X,Y,u,d){
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
    res = .C("rrenv",
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
        est = rrenv.given_du(X,Y,u,d)
        est.Betas = est$beta
        est.Ys = est$fitvalues
        est.residual_variance = est$residual_variance
        NR = (p+u-d)*d + ((r*(r+1))/2)
        Lud = 0
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
        est = rrenv.given_du(X,Y,u,d)
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
#main function
rrenv = function(X,Y,u=NULL,d=NULL,Beta=NULL){
  X = as.matrix(X)
  Y = as.matrix(Y)
  dyn.load("rrenv.so")
  rrenv.given_du <- function(X,Y,u,d){
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
    est = .C("rrenv",
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
      print(err.table)
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
  est =  rrenv.given_du(X=X, Y=Y,u=u,d=d)
  return(est)

}

#example 
# Generate Y
mu <- sample(0:3, 15, replace = T)      
Sigma = matrix(0.5 , nrow=15, ncol=15) + diag(1.5, nrow=15, ncol=15)
Y <- mvrnorm(5000, mu = mu, Sigma = Sigma)
# Calculate kernel density estimate
X <- matrix(NA, nrow=5000, ncol=15)
for(i in 1:15){
  X[,i] = rnorm(5000, mean = sample(1:3, 1), sd = sample(1:3, 1))
}
est <- rrenv(t(X),t(Y),u=10,d=5,Beta=NULL)
#returns
beta<-est$beta
fit_values<-est$fit_values
residuals<-est$residuals
residual_variance<-est$residual_variance
error_table<-est$error_table
