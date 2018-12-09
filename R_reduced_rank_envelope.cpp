//Project: Reduced_Rank_Envelope_Regression
//liyixi
//2018.11.15
//Reduced_rank_envelope




#include <stdio.h>
//#include <RcppEigen.h>
//#include <BH.h>
#include "ReadMatrix615.h"

// [[Rcpp:depends(RcppEigen)]]
// [[Rcpp:depends(BH)]]

using namespace std;
using namespace boost;
using namespace Eigen;

MatrixXd sqrtSingularvalues(MatrixXd D){
	MatrixXd sqrtD = MatrixXd::Zero(D.rows(),D.rows()); 
	ArrayXXd sqrtArrayD = D.array().sqrt();
  	int i=0;
  	while(i<D.rows()){
  		sqrtD(i,i) = sqrtArrayD(i,0);
  		i++;
  	}
  	return sqrtD;
}

class ReduceRankEnvelope{
public:
	MatrixXd X;
	MatrixXd Y;
	MatrixXd Gamma;
	MatrixXd Gamma0;
	MatrixXd Sx_sqrt;
	MatrixXd Syg_sqrt;
	MatrixXd Syg_sqrt2;
	MatrixXd Cancov_truncated;
	int u, d;
	MatrixXd Beta;
	MatrixXd FitValues;
	MatrixXd Residuals;
	MatrixXd ResidualVariance;
	ReduceRankEnvelope(MatrixXd dX, MatrixXd dY, int du, int dd, MatrixXd dGamma, MatrixXd dGamma0){
		X = dX;
		Y = dY;
		u = du;
		d = dd;
		Gamma = dGamma;
		Gamma0 = dGamma0;
	}

	void initialize(){
	cout << "---------Reduced-rank Envelope Regression--------" << endl;
	int p = X.rows(); //X: p*N
	int r = Y.rows(); //Y: r*N
	int N = Y.cols();

	//estimate envelope given u
	MatrixXd Y_gamma(u,N); 
	if(u==r){ 
		Y_gamma = Y;
	}
	else{
		Y_gamma = Gamma.transpose() * Y;
	}
	 // Y_gamma: u * N

	//estimate the reduced rank envelope estimator
	//variance and covariance
	MatrixXd Sx = (X * X.transpose())/ double(N);  //p*p
	MatrixXd Syg = (Y_gamma*Y_gamma.transpose())/ double(N); //u*u 
	MatrixXd Sxyg = (X * Y_gamma.transpose()) /double(N); //p*u

	//canonical covariance
	
	JacobiSVD<MatrixXd> svdx(Sx.inverse(),ComputeThinU | ComputeThinV);
	JacobiSVD<MatrixXd> svdy(Syg.inverse(),ComputeThinU | ComputeThinV);
	JacobiSVD<MatrixXd> svdyp(Syg,ComputeThinU | ComputeThinV);
	
	Sx_sqrt = svdx.matrixU() * sqrtSingularvalues(svdx.singularValues()) * svdx.matrixU().inverse();
	Syg_sqrt = svdy.matrixU() * sqrtSingularvalues(svdy.singularValues()) * svdy.matrixU().inverse();
	Syg_sqrt2 = svdyp.matrixU() * sqrtSingularvalues(svdyp.singularValues())* svdyp.matrixU().inverse();
	MatrixXd Cancov = Syg_sqrt * Sxyg.transpose() * Sx_sqrt; // u*u * u*p * p*p = u*p
	
	//truncate canonical covariance
	JacobiSVD<MatrixXd> svd(Cancov,ComputeThinU | ComputeThinV);
	MatrixXd truncated = MatrixXd::Zero(u,u);
  	int i=0;
  	while(i<d){
  		truncated(i,i) = svd.singularValues()(i,0);
  		i++;
  	}
  		
  	Cancov_truncated = svd.matrixU() * truncated * svd.matrixV().transpose(); // u*p, with rank = d
  	}

  	//estimate beta
  	MatrixXd beta(){
  		cout << "------------Beta Estimation-----------" << endl;
  		Beta = Gamma * Syg_sqrt2 * Cancov_truncated * Sx_sqrt; //r*u * u*u * u*p * p*p = r*p
  		return Beta;
  	}

  	//estimate fit values and residuals
  	MatrixXd fitValues(){
  		cout << "-------------Fitted Values------------" << endl;
  		FitValues = Beta * X;
  		return FitValues;
  	}

  	MatrixXd residualValues(){
  		cout << "---------------Residuals--------------" << endl;
  		Residuals = Y - FitValues;
  		return Residuals;
  	}

  	MatrixXd residualVariance(){
  		cout << "---------Variance of Residuals---------" << endl;
  		MatrixXd Sy = Y * Y.transpose() / double(Y.cols()); //r*r
  		if(u== Y.rows()){
  			ResidualVariance = Gamma * Syg_sqrt2 * (MatrixXd::Identity(u,u) - Cancov_truncated * Cancov_truncated.transpose()) * Syg_sqrt2 * Gamma.transpose();
  		}
  		else{
  			ResidualVariance = Gamma * Syg_sqrt2 * (MatrixXd::Identity(u,u) - Cancov_truncated * Cancov_truncated.transpose()) * Syg_sqrt2 * Gamma.transpose() + Gamma0 * Gamma0.transpose() *  Sy * Gamma0 * Gamma0.transpose();
  		}
  		return ResidualVariance;
  	}
  	
};

extern "C"{
void R_reduced_rank_envelope(double* R_X, int* nrow_X, int* ncol_X,
	                         double* R_Y, int* nrow_Y, int* ncol_Y,
	                         int* R_u, int* R_d,
	                         double* R_Gamma, int* nrow_G, int* ncol_G,
	                         double* R_Gamma0,int* nrow_G0, int* ncol_G0,
	                         double* R_beta, double* R_fitvalues, 
	                         double* R_residuals, double* R_residualvariances){
 	

    Map<MatrixXd> X(R_X,*nrow_X,*ncol_X);
    Map<MatrixXd> Y(R_Y,*nrow_Y,*ncol_Y);
    Map<MatrixXd> Gamma(R_Gamma,*nrow_G,*ncol_G);
    Map<MatrixXd> Gamma0(R_Gamma0,*nrow_G0,*ncol_G0);

//parameters
	int u = *R_u;
	int d = *R_d;

//reduced rank envelope
ReduceRankEnvelope rre(X, Y, u, d, Gamma, Gamma0); //input(TODO: use BIC to choose u,d)
rre.initialize();
MatrixXd Beta_hat = rre.beta();  //r*p
MatrixXd FitValues = rre.fitValues(); //r*N
MatrixXd Residuals = rre.residualValues(); //r*N
MatrixXd ResidualVariance = rre.residualVariance(); //r*r 

//output
for(int i=0;i<*nrow_Y;i++){
	for(int j=0; j<*nrow_X;j++){
		R_beta[*nrow_Y*j+i] = Beta_hat(i,j);
	}
}

for(int i=0;i<*nrow_Y;i++){
	for(int j=0; j<*ncol_Y;j++){
		R_fitvalues[*nrow_Y*j+i] = FitValues(i,j);
		R_residuals[*nrow_Y*j+i] = Residuals(i,j);
	}
}

for(int i=0;i<*nrow_Y;i++){
	for(int j=0; j<*nrow_Y;j++){
		R_residualvariances[*nrow_Y*j+i] = ResidualVariance(i,j);
	}
}


}
}