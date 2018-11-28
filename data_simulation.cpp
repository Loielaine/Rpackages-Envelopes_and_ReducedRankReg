//Project: Reduced_Rank_Envelope_Regression
//liyixi
//2018.11.04
//Data_simulation


#include<iostream>
#include<fstream>
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

//parameters
	int p,r,u,d,N;
	p = atoi(argv[1]);
	r = atoi(argv[2]);
	u = atoi(argv[3]);
	d = atoi(argv[4]);
	N = atoi(argv[5]);

//random uniform distribution
	mt19937 rng;
	rng.seed(time(0));
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
	//cout << Gamma_orthogonal << endl;
	//cout << Gamma_orthogonal.rows() << "*" << Gamma_orthogonal.cols() << endl;
	//cout << Gamma_orthogonal.transpose() * Gamma_orthogonal << endl;

	FullPivLU<MatrixXd> lu_gamma(Gamma.transpose());
	MatrixXd Gamma_0 = lu_gamma.kernel();
	//cout << Gamma_0.rows() << "*" << Gamma_0.cols() << endl;
	

	//eta
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
	MatrixXd Sigma(r,r);
	//Omega_0
	if(r==u){
		Sigma = Gamma_orthogonal * Omega * Gamma_orthogonal.transpose();
	}
	else{
	MatrixXd Omega_0(r-u,r-u);
	for(int i=0; i<(r-u); i++){
		for(int j=0; j<(r-u); j++){
			Omega_0(i,j) = 5 * pow((-0.5),abs(i-j));
		}	
	}

	Sigma = Gamma_orthogonal * Omega * Gamma_orthogonal.transpose() + Gamma_0 * Omega_0 * Gamma_0.transpose(); //r*u * u*u * u*r 
	}
//(3)X generating
	//X~N(0,Ip) p*N
	MatrixXd X(p,N);
	
	for(int i=0; i<p; i++){
		for(int j=0; j<N; j++){
			rn = norm(rng);
			X(i,j) = rn;
		}	
	}

	MatrixXd X_centered = X.rowwise() - X.colwise().mean();

//(4)Y generating
	//error~N(0,Sigma) r*N
	MatrixXd Z(r,N); //Z~N(0,IrN) r*N
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
	
	//Write X, Y and Beta in three files 
	ofstream Xfile, Bfile, Yfile;
        Xfile.open ("X.txt");
        Xfile << X;
        Xfile.close();
        Yfile.open ("Y.txt");
        Yfile << Y;
        Yfile.close();
        Bfile.open ("Betas.txt");
        Bfile << Beta;
        Bfile.close();

}






