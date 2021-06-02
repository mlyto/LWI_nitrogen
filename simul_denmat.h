#pragma once
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#include <vector>
#include <complex>
#include "Eigen/Sparse"
#include "Eigen/Dense"
#include <fftw3.h>
#include "initial.h"
#include "mesh.h"
#include "legendre_matrices.h"
#include <cstdlib>

using std::size_t;
using std::vector;
using std::complex;

typedef Eigen::MatrixXcd DenMatComp;
typedef Eigen::RowVectorXd DenVec;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseMatrix<complex<double> >SpMatComp;
typedef Eigen::Triplet<double> T;


class simul_denmat
{
public:
	simul_denmat(mesh*, legendre_matrices*, initial*, double, double);
	~simul_denmat();

	vector<double> nonlin1;  
	vector<double> nonlin2; 
	vector<double> nonlin3; 
	vector<double> P_theta, dPdt_theta;
	vector<double> P_theta1, dPdt_theta1;
	vector<double> Pol_mu, P_Xt, P_Bt;
	
        void Cos2_FT();
	void initial_data(size_t, int, vector<double>&);
	void Runge_Kutta_pump(size_t);
	void Runge_Kutta_seed(vector<double>&, vector<double>&);
	void cut();
	void show_H();
	void X_B_population();
	void X_B_InTime();
        void show_polarization();
	
private:

	double Eng_X, Eng_B;
	double valH, delta_Eij;
	size_t nbT, nbt, nbT2;
	size_t i_del0;
	size_t i_Tp, ConvC, i_sigma;
	size_t start_ff;
	size_t pump_peak;
	double dTp, dt, dz, C_ref;
	double dPdt_mu, d2Pdt_mu, d2Pdt_theta;
	vector<double> Tscale;
	vector<double> time_p, node, Sp1, Sp2, Sp3;
	vector<double> d2Pdt_theta1;
	vector<double> Efield, Uot;
	vector<double> E_half, Uo_half;
	vector<double> Potential;
	vector<double> dPdt_mu_prev;
	int mult;
	double dT_2, dT_6;
	double dt_2, dt_6;
	double NormC;

	SpMat Hamiltonian_0, Hamiltonian_X, Hamiltonian_B;
	SpMat PermuteD;
	vector<T> tripletH_0, tripletH_X, tripletH_B;
	vector<T> tripletP, tripletOne;

        vector<SpMatComp> Rho_in;
	vector<SpMatComp> CR_o, CR_x, CR_b, CR_, CR_sum, CT_, CR_theta;
	vector<SpMatComp> Rho_N20, Rho_X0, Rho_B0;
	
	DenMatComp K1, K2, K3, K4;
	DenMatComp K1_o, K2_o, K3_o, K4_o;
	DenMatComp K1_x, K2_x, K3_x, K4_x;
	DenMatComp K1_b, K2_b, K3_b, K4_b;
	DenMatComp Rho_N2, Rho_X, Rho_B, Rho, Rho_Bs, Rho_theta;
	DenMatComp Rho_N2_asym, Aodd, Aeven;
	double odd_even;	
	DenMatComp AUX;

	SpMat ONE_tr;
	SpMatComp IDNT_tr;
	SpMatComp H0, H0_o, H0_x, H0_b;
	SpMatComp H_mat_0, H_mat_X, H_mat_B, H_mat;
	SpMatComp Permute;

	DenMatComp ExpHdt_0, ExpHdt_X, ExpHdt_B;

	vector<double> Cos2_T_N2, Cos2X_T, Cos2B_T, Cos2_T;
	vector<double> Cos2t;
	vector<double> Pol;
	vector<double> P_T, dP_theta_dT, d2P_theta_dT;

	DenVec P0_peak, PX_peak, PB_peak, PBs_peak;
	DenVec Pa_X;
	DenVec Pa_B;
	double P_X_after, P_B_after;
	

	vector<double> omeg;
	size_t size_tempo;
	fftw_complex *Cos2T;
	fftw_complex *COS2;

	void Hamiltonian0();
	void Permutation();
	void Exp_dHdt();
	void sym2asym();
};


