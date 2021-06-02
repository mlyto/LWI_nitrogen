#pragma once
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#include <vector>
#include "Eigen/Sparse"
#include "Eigen/Dense"

using std::vector;

typedef Eigen::RowVectorXd DenVec;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;
typedef Eigen::MatrixXd DenMat;  


class initial
{
public:
	initial();
	~initial();
	
	vector<SpMat> BoltzM_0;
	vector<SpMat> BoltzM_X;
	vector<SpMat> BoltzM_B;

	void cut();

private:

	double Rot_0;     
	double E_term;
	double Norm_0;
	
	vector<vector<T> > tripletBoltz_0;
	vector<vector<T> > tripletBoltz_X;
	vector<vector<T> > tripletBoltz_B;	

	DenVec Pb_X;
	DenVec Pb_B;

};


