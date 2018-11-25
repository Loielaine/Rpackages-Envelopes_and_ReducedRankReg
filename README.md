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
 
![ScreenShot](https://github.com/shnlee-ds/Rpackages-Envelopes_and_ReducedRankReg/blob/master/Data%20simulation.png)
