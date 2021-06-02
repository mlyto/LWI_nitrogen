#pragma once
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#include <vector>
#include "Eigen/Sparse"

using std::vector;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

class legendre_matrices
{
public:
	legendre_matrices();
	~legendre_matrices();


	vector<SpMat> CR_X_sparse;
	vector<SpMat> CR_B_sparse;
	vector<SpMat> CR_0_sparse;
	vector<SpMat> CR_sparse;
	vector<SpMat> CR_theta_sparse;
	vector<SpMat> CT_sparse;	

	void cut();

private:

	double RX_JJ, RB_JJ, R0_JJ, R_JJ;
	double RX_2, RB_2, R0_2, R_2;
	double S_1;

	vector<vector<T> > tripletX;
	vector<vector<T> > tripletB;
	vector<vector<T> > triplet0;
	vector<vector<T> > tripletCT;
	vector<vector<T> > tripletCR;
	vector<vector<T> > tripletCR_theta;
};
