# R package: Envelopes and Reduced-rank Regression

A project developing a R-package which enables us to build Reduced-rank envelopes regression models

## Contributors

* __Pedro Orozco del Pino__ (PhD student in Biostatistics at University of Michigan, Ann Arbor)
* __Sunghun Lee__ (Masters student in Applied Statistics at University of Michigan, Ann Arbor)
* __Yixi Li__ (Masters student in Survery Methodology, Statistical Science at University of Michigan, Ann Arbor)

### Reference

This project is based on a research paper: 
R. Dennis Cook, Liliana Forzani, Xin Zhang; Envelopes and reduced-rank regression, Biometrika, Volume 102, Issue 2, 1 June 2015, Pages 439–456, https://doi.org/10.1093/biomet/asv001


### Description

In multivariate regression analysis there is the case in which the coefficients matrix is not full rank. In such cases the reduced-rank methodology can get a better estimate than OLS. Another improvement that can be done in multivariate analysis is when the estimation of the beta matrix can be done with less information without compromissing accuracy. For this, we look a subspace in which the projection the response is independent of the covariates. The envelope regression methodology estimates this subspace and also outperforms OLS. In this project we build a program that combines reduced-rank and envelope regression methodologies. Each of the methodologies don't compite with each other since they solve different challenges. It is because of this that the combine methodology outperforms OLS, reduced-rank and envelopes.

### Data simulation

We compared the performance (in terms of the fit of model) of our new Reduced-rank envelope package and existing R packages (for OLS, Reduced-rank regression and Envelope regression). As you can see from the below plot, the Reduced-rank envelope regression model with our new R package always showed the best performance in all the three cases. In addition, it is also impressive that Reduced-rank envelope was stable between different sample sizes from 160 to 2000. In particular, it is noteworthy that in a small sample size, its performance was superior to other methods.

![ScreenShot](https://github.com/shnlee-ds/Rpackages-Envelopes_and_ReducedRankReg/blob/master/simulation.png)

### Implementation of reduced-rank envelope regression
The pseudocode of the implementation of reduced-rank envelope regression is as follows:
#### Data: Input multivariate data (X,Y), with dimension of X is p\*N, and dimension of Y is r\*N
#### Result: Return the best selection of (u,d), estimates of beta, fitted values, residuals, variance of residuals, error table
- if (u,d) is not given, do:
  - select (u,d) by _rrenv_choose_du(X, Y, Beta=NULL)_
  - return the best selction of (u,d)
- if (u,d) is given, do:
  - check whether the input (u,d) is valid
- redeced-rank envelope regression:
  - envelope estimation:
    - estimate the envelope by _envlp_(X,Y,u,d)_ 
    - return Gamma, Gamma0
  - reduecd-rank regression:
    - estimate beta by _R_reduced_rank_envelope_given_u_d(X,Y,u,d, Gamma, Gamma0)_ in C++
    - return estimates of beta, fitted values, residuals, variance of residuals, error table

Reduced-rank enevelope regression is implemented in R and C++, and is warpped into a R package, _rrenv_.
```ruby
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
```
A example to use _rrenv_ in R is:
```ruby
#input X and Y
Y<-as.matrix(read.table("Y.txt"))
X<-as.matrix(read.table("X.txt"))
est <- rrenv(X,Y,20,10,NULL)
#returns
beta<-est$beta
fit_values<-est$fit_values
residuals<-est$residuals
residual_variance<-est$residual_variance
error_table<-est$error_table
```

### Selection of parameter d and u

When you have a data, a function _rrenv.choose_du(X, Y, Beta=NULL)_ in our package recommends the combination of parameter (d, u) for a reduced-rank envelope regression modeling based on the error calculation. When you input the matrices of X, Y and the true Beta matrix(the left case of below figure), the function calculates  ![](https://latex.codecogs.com/gif.latex?%5Cleft%20%5C%7C%20%5Chat%7B%5Cbeta%7D-%5Cbeta%20%5Cright%20%5C%7C_F), while it calculates ![](https://latex.codecogs.com/gif.latex?AIC%20%3D%202N_%7BRE%7D%20-%202L_%7Bu%2Cd%7D%20%3D%20%28p&plus;r-d%29d&plus;r%28r&plus;1%29/2%20-%202%5Chat%7BL%7D_%7Bu%2Cd%7D) when you only input the matrices of X and Y without the true Beta matrix (the right case of below figure).

The function _rrenv.choose_du(X, Y, Beta=NULL)_ returns not only the heatmap but also the best combination of (d,u) based on the calculated errors.

![ScreenShot](https://github.com/shnlee-ds/Rpackages-Envelopes_and_ReducedRankReg/blob/master/choose_du.png) 
