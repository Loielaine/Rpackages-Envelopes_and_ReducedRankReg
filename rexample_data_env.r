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
