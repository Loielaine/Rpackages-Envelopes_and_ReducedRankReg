\name{rrenv-package}
\alias{rrenv-package}
\alias{rrenv}
\docType{package}
\title{
\packageTitle{rrenv}
}
\description{
\packageDescription{rrenv}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{rrenv}
\packageIndices{rrenv}
~~ An overview of how to use the package, including the most important functions ~~
This package is to build reduced-rank envelope regression models.
}
\author{
\packageAuthor{rrenv}
Pedro Orozco del Pino, Sunghun Lee, Yixi Li
Maintainer: \packageMaintainer{rrenv}
Yixi Li <liyixi@umich.edu>

}
\references{rrenv}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
~~ documentation directory ~~
R. Dennis Cook, Liliana Forzani, Xin Zhang; Envelopes and reduced-rank regression, Biometrika, Volume 102, Issue 2, 1 June 2015, Pages 439–456, https://doi.org/10.1093/biomet/asv001 
\keyword{rrenv}
reduced rank envelope regression
\seealso{
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
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
}
