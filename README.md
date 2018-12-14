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

In multivariate regression analysis there is the case in which the coefficients matrix is not full rank. In such cases the reduced-rank methodology can get a better estimate than OLS. Another improvement that can be done in multivariate analysis is when the estimation of the beta matrix can be done with less information without compromissing accuracy. For this, we look a subspace in which the projection of the response is immaterial of the covariates. This means that the kernel projection of the response is independent of the projection given the covariates and also is independent of the covariates. The envelope regression methodology estimates this subspace and also outperforms OLS. Each of this methodologies is an improvement of OLS multivariate regression that could be combined. In this project we build a program that combines reduced-rank and envelope regression methodologies into an R package. Each of the methodologies don't compite with each other since they attain different challenges. The reduced-rank envelope regression combines strenghts of each of the methodologies and it doesn't amplifies weaknes, it is because of this that the combine methodology outperforms OLS, reduced-rank and envelopes.

### Data simulation

We compared the performance (in terms of the fit of model) of our new Reduced-rank envelope package and existing R packages (for OLS, Reduced-rank regression and Envelope regression). As you can see from the below plot, the Reduced-rank envelope regression model with our new R package always showed the best performance in all the three cases. In addition, it is also impressive that Reduced-rank envelope was stable between different sample sizes from 160 to 2000. In particular, it is noteworthy that in a small sample size, its performance was superior to other methods. The first escenario (left most panel) is design to be ideal for reduced rank regression because the rank of the betas is only 1 and the envelope was the same as the amounts of predictios. For this case the RRE outperforms all three methodologies by almost one magnitud in small sample sizes. The second escenario shows that when the rank of beta matrix is close the dimension of the envelope the RRE gives a modest improvement over envelope regression. Lastly we show a case in which the reduce rank and the envelope methodologies by themselves should be useful. RRE still outperforms all of them in all of the sample sizes. The takeaway from this is that RRE adapts when one of the methodologies is not ideal and when they are combines them adequatly. 
![ScreenShot](https://github.com/shnlee-ds/Rpackages-Envelopes_and_ReducedRankReg/blob/master/simulation.png)

### Implementation of reduced-rank envelope regression
The pseudocode of the implementation of reduced-rank envelope regression is as follows:
#### Data: Input multivariate data (X,Y), with dimension of X is p\*N, and dimension of Y is r\*N
#### Result: Return the best selection of (u,d), estimates of beta, fitted values, residuals, variance of residuals, error table
- if (u,d) is not given, do:
  - select (u,d) by _rrenv.choose_du(X, Y, Beta=NULL)_
  - return the best selction of (u,d)
- if (u,d) is given, do:
  - check whether the input (u,d) is valid
- redeced-rank envelope regression:
  - envelope estimation:
    - estimate the envelope by _envlp_(X,Y,u,d)_ 
    - return Gamma, Gamma0
  - reduecd-rank regression:
    - estimate beta by _rrenv_given_du(X,Y,u,d, Gamma, Gamma0)_ in C++
    - return estimates of beta, fitted values, residuals, variance of residuals, error table

Reduced-rank enevelope regression is implemented in R and C++, and is warpped into a R package, _rrenv_.

A example to use _rrenv_ in R is:
```ruby
library(MASS)
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
```

### Selection of parameter d and u

When you have a data, a function _rrenv.choose_du(X, Y, Beta=NULL)_ in our package recommends the combination of parameter (d, u) for a reduced-rank envelope regression modeling based on the error calculation. When you input the matrices of X, Y and the true Beta matrix(the left case of below figure), the function calculates  ![](https://latex.codecogs.com/gif.latex?%5Cleft%20%5C%7C%20%5Chat%7B%5Cbeta%7D-%5Cbeta%20%5Cright%20%5C%7C_F), while it calculates ![](https://latex.codecogs.com/gif.latex?AIC%20%3D%202N_%7BRE%7D%20-%202L_%7Bu%2Cd%7D%20%3D%20%28p&plus;r-d%29d&plus;r%28r&plus;1%29/2%20-%202%5Chat%7BL%7D_%7Bu%2Cd%7D) when you only input the matrices of X and Y without the true Beta matrix (the right case of below figure).

The function _rrenv.choose_du(X, Y, Beta=NULL)_ returns not only the heatmap but also the best combination of (d,u) based on the calculated errors.

![ScreenShot](https://github.com/shnlee-ds/Rpackages-Envelopes_and_ReducedRankReg/blob/master/choose_du.png) 
