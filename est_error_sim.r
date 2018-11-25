#Sunghun Lee
#Evaluate and visualize our simulations
library(Renvlp)
library(expm)
library(readr)
library(data.table)
setwd("/Users/shnlee/Desktop/615final/simulation data") #This folder is uploaded on the githup repo

samples = c(160,270,450,740,1200,2000)
cases = list(1=c(d=1,u=10,p=10,r=20),2=c(d=4,u=5,p=6,r=20),3=c(d=5,u=10,p=15,r=20))



err = data.frame(CASE=rep(1:3,each=6), SIZE=rep(samples,3),OLS=rep(NA,6), 
                 ENV=rep(NA,3*6), RR=rep(NA,3*6), RRENV=rep(NA,3*6))
par(mfrow=c(1,3))
for(i in 1:3){
  for(j in samples){
    
    ## Read inputs
    txt = paste(i,"_",j, sep="")
    X = as.matrix(fread(paste('X',txt,'.txt',sep="")))
    Y = as.matrix(fread(paste('Y',txt,'.txt',sep="")))
    Beta = as.matrix(fread(paste('Betas',txt,'.txt',sep="")))
    
    ## Test r and p
    if(cases[[i]]["r"]==dim(Y)[1]&cases[[i]]["p"]==dim(X)[1]){
      
      
      ## Get Estimation errors - From Pedro's work
      red.rank.envelope<-function(X,Y,u,d){
        E<-env(t(X),t(Y),u)
        Gamma<-E$Gamma
        Gamma0<-E$Gamma0
        N<-dim(Y)[2]
        ifelse(u!=dim(Y)[1],Y.g<-t(Gamma)%*%as.matrix(Y),Y.g<-Y) 
        r.g<-dim(Y.g)[1]
        Sx<-X%*%t(X)/N
        Sxyg<-t(X%*%scale(t(Y.g), scale = FALSE))/N
        Syg<-Y.g%*%t(Y.g)/N
        Cyg.x<-sqrtm(solve(Syg))%*%Sxyg%*%sqrtm(solve(Sx))
        D.truncated<-diag(c(svd(Cyg.x)$d[1:d],rep(0,length(svd(Cyg.x)$d)-d)))
        Cyg.x.d<-svd(Cyg.x)$u%*%D.truncated%*%t(svd(Cyg.x)$v)
        beta.r.e<-Gamma%*%sqrtm(Syg)%*%Cyg.x.d%*%solve(sqrtm(Sx))
        return(beta.r.e)
      }
      
      OLS<-function(X,Y){
        r<-dim(Y)[1]
        p<-dim(X)[1]
        betas<-matrix(rep(0,r*p),ncol = p)
        for(i in 1:r)
          betas[i,]<-lm(Y[i,]~t(X))$coefficients[-1]
        return(betas)
      }
      
      # Get estimation errors for each model
      p = cases[[i]]['p']; d = cases[[i]]['d'] ; u = cases[[i]]['u'] ; r = cases[[i]]['r']
      ols = norm(OLS(X,Y)-Beta,type="F")
      env = norm(env(t(X),t(Y),u)$beta-Beta,type="F")
      rr = norm(red.rank.envelope(X,Y,dim(Y)[1],2)-Beta,type = "F")
      rrenv = norm(red.rank.envelope(X,Y,u,d)-Beta,type = "F")
      
      errs = c(ols,env,rr,rrenv)
      err[err$CASE==i&err$SIZE==j,3:6] =  errs
    }
  }
  
  # Visualize errors by each case
  d = err[err$CASE==i,-c(1,2)]
  matplot(d, type = c("b"), pch=1:4, col = 1:4, ylim=c(0,4.1),
          xlab="Sample size", ylab="Averaged estimation error", xaxt = "n",
          main=paste("(d,u,p,r) = (",cases[[i]][1],",",
                     cases[[i]][2],",",cases[[i]][3],",",cases[[i]][4],")", sep=""))
  axis(side=1,at=1:6,labels=c("160","270","450","740","1200", "2000"))
  legend("topright", legend = c("OLS","ENV", "RR", "RRENV"), col=1:4, pch=1:4)
}
