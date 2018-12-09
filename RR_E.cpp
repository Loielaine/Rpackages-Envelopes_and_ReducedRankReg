#include <RcppEigen.h>

//[[Rcpp::depends(RcppEigen)]]

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
//[[Rcpp::export]]
MatrixXd RrEnvCpp(MatrixXd X,
              MatrixXd Y,
              int u, int d,
              MatrixXd Gamma,
              MatrixXd Gamma0){
           
            int r = Y.rows();
            int N = Y.cols();
            
            MatrixXd Y_gamma(u,N); 
            if(u==r){ 
            Y_gamma = Y;
            }
            else{
            Y_gamma = Gamma.transpose() * Y;
            }
            //int r_gamma = Y_gamma.rows(); 
            
            MatrixXd Sx = (X * X.transpose())/ double(N);  
            
            MatrixXd Syg = (Y_gamma*Y_gamma.transpose())/ double(N);
            MatrixXd Sxyg = (X * Y_gamma.transpose()) /double(N); 
            
            JacobiSVD<MatrixXd> svdx(Sx.inverse(),ComputeThinU | ComputeThinV);
            JacobiSVD<MatrixXd> svdy(Syg.inverse(),ComputeThinU | ComputeThinV);
            JacobiSVD<MatrixXd> svdyp(Syg,ComputeThinU | ComputeThinV);
            
            MatrixXd Sx_sqrt = svdx.matrixU() * sqrtSingularvalues(svdx.singularValues()) * svdx.matrixU().inverse();
            MatrixXd Syg_sqrt = svdy.matrixU() * sqrtSingularvalues(svdy.singularValues()) * svdy.matrixU().inverse();
            MatrixXd Syg_sqrt2 = svdyp.matrixU() * sqrtSingularvalues(svdyp.singularValues())* svdyp.matrixU().inverse();
            MatrixXd Cancov = Syg_sqrt * Sxyg.transpose() * Sx_sqrt;
            
            JacobiSVD<MatrixXd> svd(Cancov,ComputeThinU | ComputeThinV);
            MatrixXd truncated = MatrixXd::Zero(u,u);
            int i=0;
            while(i<d){
            truncated(i,i) = svd.singularValues()(i,0);
            i++;
            }	
            MatrixXd Cancov_truncated = svd.matrixU() * truncated * svd.matrixV().transpose(); // u*p, with rank = d
            MatrixXd Beta = Gamma * Syg_sqrt2 * Cancov_truncated * Sx_sqrt; //r*u * u*u * u*p * p*p = r*p
           return Beta;
            
            }


