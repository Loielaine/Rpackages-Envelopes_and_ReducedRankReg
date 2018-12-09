//Sunghun Lee
//Data_simulation for choosing d and u


#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<cstdlib>
#include<climits>
#include<cmath>
#include<Eigen/Dense>
#include<boost/random/mersenne_twister.hpp>
#include<boost/random/uniform_real.hpp>
#include<boost/random/normal_distribution.hpp>

using namespace std;
using namespace boost;
using namespace Eigen;

int main(int argc, char *argv[]){

	int p,r,N;
	p = atoi(argv[1]);
	r = atoi(argv[2]);
	N = atoi(argv[3]);


	for(int u=1; u<r; u++){
		for(int d=1; d<=u; d++){

			cout << "u: " << u << " " << " d: " << d << endl;
	
		//random uniform distribution
			mt19937 rng;
			rng.seed(d+u);
			uniform_real<> uni(0 ,1);
			normal_distribution<> norm(0 ,1);
			double rn;

		//Generating
		//(1)Beta generating
			//Gamma
			MatrixXd Gamma(r,u);
			for(int i=0; i<r; i++){
				for(int j=0; j<u; j++){
					rn = uni(rng);
					Gamma(i,j) = rn;
				}	
			}

			JacobiSVD<MatrixXd> svd_gamma(Gamma,ComputeThinV | ComputeThinU);
			MatrixXd Gamma_orthogonal = svd_gamma.matrixU();

			FullPivLU<MatrixXd> lu_gamma(Gamma.transpose());
			MatrixXd Gamma_0 = lu_gamma.kernel();

			MatrixXd eta(u,d);
			for(int i=0; i<u; i++){
				for(int j=0; j<d; j++){
					rn = norm(rng);
					eta(i,j) = rn;
				}	
			}

			JacobiSVD<MatrixXd> svd_eta(eta,ComputeThinV | ComputeThinU);
			MatrixXd eta_orthogonal = svd_eta.matrixU();

			//B
			MatrixXd B(d,p);
			for(int i=0; i<d; i++){
				for(int j=0; j<p; j++){
					rn = uni(rng);
					B(i,j) = rn;
				}	
			}

			//Beta
			MatrixXd Beta = Gamma_orthogonal * eta_orthogonal * B; //r*u * u*d * d*p
			

			double f_norm = Beta.norm();
			
			//double f_norm = 0;
			for(int i=0; i<Beta.rows(); i++){
				for(int j=0; j<Beta.cols(); j++){
					Beta(i,j) = Beta(i,j)/f_norm;
				}	
			}

			
		//(2)Sigma generating
			//Omega
			MatrixXd Omega(u,u);
			for(int i=0; i<u; i++){
				for(int j=0; j<u; j++){
					Omega(i,j) = pow((-0.9),abs(i-j));
				}	
			}

			//Omega_0
			MatrixXd Omega_0(r-u,r-u);
			for(int i=0; i<(r-u); i++){
				for(int j=0; j<(r-u); j++){
					Omega_0(i,j) = 5 * pow((-0.5),abs(i-j));
				}	
			}

			MatrixXd Sigma = Gamma_orthogonal * Omega * Gamma_orthogonal.transpose() + Gamma_0 * Omega_0 * Gamma_0.transpose();

		//(3)X generating
			MatrixXd X(p,N);
			
			for(int i=0; i<p; i++){
				for(int j=0; j<N; j++){
					rn = norm(rng);
					X(i,j) = rn;
				}	
			}

			MatrixXd X_centered = X.rowwise() - X.colwise().mean();

		//(4)Y generating
			MatrixXd Z(r,N);
			for(int i=0; i<r; i++){
				for(int j=0; j<N; j++){
					rn = norm(rng);
					Z(i,j) = rn;
				}	
			}
			// decompose Sigma = A * A.transpose()
			LLT<MatrixXd> ll_sigma(Sigma); //LL^T Cholesky decomposition
			MatrixXd A = ll_sigma.matrixL();
			
			//error
			MatrixXd error(r,N);
			error = A * Z;
			//Y
			MatrixXd Y = Beta * X_centered + error;
			

			string sx = "_X_";
			string sy = "_Y_";
			string sb = "_Beta_";
			stringstream ssd;
			ssd << d;
			string sd = ssd.str();
			stringstream ssu;
			ssu << u;
			string su = ssu.str();
			stringstream sN;
			sN << N;
			string s0 = sN.str();
			string stxt = ".txt";
			string Xname = su + sx+ sd + stxt;
			string Yname = su + sy+ sd + stxt;
			string Bname = su + sb+ sd + stxt;

				//Write X, Y and Beta in three files 
			// cout<< result(1,1) << " " << filename.c_str() << endl;
			ofstream Xfile, Bfile, Yfile;
		        Xfile.open (Xname.c_str());
		        Xfile << X;
		        Xfile.close();
		        Yfile.open (Yname.c_str());
		        Yfile << Y;
		        Yfile.close();
		        Bfile.open (Bname.c_str());
		        Bfile << Beta;
		        Bfile.close();
	   	}
	}
   	return 0;
}
