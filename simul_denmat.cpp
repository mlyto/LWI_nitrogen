///////////////////////////////////////////////////////////////////////////////////                    
//                                                                               //                    
//  M.Lytova, 2018-2019, Project: N2+ lasing in rotationaly excited medium       //                    
//                                                                               //                    
/////////////////////////////////////////////////////////////////////////////////// 
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "simul_denmat.h"
#include "Global.h"
#include "spline.h"

#define REAL 0
#define IMAG 1

const complex<double> ImagUnit(0.0, 1.0);

simul_denmat::simul_denmat(mesh* MyMesh, legendre_matrices* MyLeg, initial* MyDM, double Cref, double dz_max)
{
	Eng_X = omega_x / 2. - omX_x / 4. + omY_x / 8.;
	Eng_B = omega_b / 2. - omX_b / 4. + omY_b / 8. + Te_B;

	this->dTp = MyMesh->dTp;
	dT_2 = dTp * 0.5;
	dT_6 = dTp * 0.166666666666666666666667;

	this->dt = MyMesh->dt;
	dt_2 = dt * 0.5;
	dt_6 = dt * 0.166666666666666666666667;

	nbt = MyMesh->nb_t;
	this->i_Tp = MyMesh->i_Tp;
	this->ConvC = MyMesh->ConvC;
	this->i_del0 = MyMesh->i_del0;
	this->i_sigma = MyMesh->i_sigma;
	nbT = MyMesh->nb_T;
	nbT2 = 2 * nbT;              // need intermidiate points for R-K
	Tscale = MyMesh->Tsc;        // to represent data	
	Potential = MyMesh->Upt;

	C_ref = Cref;
	dz = dz_max;

// number of time iterations for the pump before the regime of field-free propagation starts
	start_ff = (int)(2 * Tau_on / AU_t / dTp + 0.49);
// number of time iterations corresponding to the pump peak
	pump_peak = (int)(Tau_on / AU_t / dTp + 0.49);	
	
	std::cerr << " dt = " << dt << ",  dTp = " << dTp << ",  i_Tp = " << i_Tp
		  << ",  ConvC = " << ConvC << std::endl;
	std::cerr << " i_del0 = " << i_del0 <<  ",    start_ff = " << start_ff 
		<< ",   i_sigma = " << i_sigma << std::endl;
	std::cerr << " pump_peak = " << pump_peak << std::endl;

	Rho_in.resize(Jmax01);
	for (size_t m = 0; m < Jmax01; m++)
	{
		Rho_in[m] = SpMatComp(Jmax_XB2, Jmax_XB2);
	}

	AUX = DenMatComp(Jmax_XB2, Jmax_XB2);

	Cos2_T_N2.resize(nbT);
	Cos2X_T.resize(nbT);
	Cos2B_T.resize(nbT);
	Cos2_T.resize(nbT);
	
	Pol.resize(nbt);
	Cos2t.resize(nbt);
	P_Xt.resize(nbt);
	P_Bt.resize(nbt);
	P_theta1.resize(nbt);
	dPdt_theta1.resize(nbt);
	d2Pdt_theta1.resize(nbt);
	P_theta.resize(nbt);
	dPdt_theta.resize(nbt);
	nonlin1.resize(nbt);
	nonlin2.resize(nbt);
	nonlin3.resize(nbt);
	
	Uot.resize(nbt);
	Uo_half.resize(nbt);
	Pol_mu.resize(nbt);
	dPdt_mu_prev.resize(nbt);

	dPdt_mu = 0.0;
	d2Pdt_mu = 0.0;
	d2Pdt_theta = 0.0;
	
	Permutation(); // auxiliary matrix
	Hamiltonian0(); // filling sparse matrices H0_0, H0_X, H0_B and H0
	
	for (size_t j = 0; j < Jmax_XB2; j++)                 
	{
	  if (j < Jmax1) {
	    tripletOne.push_back(T(j, j, alpha_tr_X));
	  }
	  else
	    tripletOne.push_back(T(j, j, alpha_tr_B));
	}

	ONE_tr = SpMat(Jmax_XB2, Jmax_XB2);
	ONE_tr.setFromTriplets(tripletOne.begin(), tripletOne.end());
	IDNT_tr = SpMatComp(Jmax_XB2, Jmax_XB2);
	IDNT_tr = ONE_tr.cast<complex<double>>(); 

	CR_.resize(Jmax01);
	CR_theta.resize(Jmax01);
	CR_sum.resize(Jmax01);
	CT_.resize(Jmax01);
	CR_o.resize(Jmax01);
	CR_x.resize(Jmax01);
	CR_b.resize(Jmax01);
	Rho_N20.resize(Jmax01);
	Rho_X0.resize(Jmax01);
	Rho_B0.resize(Jmax01);
	
	for (size_t m = 0; m < Jmax01; m++)
	  {
	    CR_[m] = (MyLeg->CR_sparse[m]).cast<complex<double>>();  // includes all delta_alpha_X and delta_alpha_B
	    CR_theta[m] = (MyLeg->CR_theta_sparse[m]).cast<complex<double>>();
	    CR_sum[m] = CR_[m] + IDNT_tr;
	    CT_[m] = (MyLeg->CT_sparse[m]).cast<complex<double>>();  // includes mu
	    CR_o[m] = (MyLeg->CR_0_sparse[m]).cast<complex<double>>(); //includes delta_alpha
	    CR_x[m] = (MyLeg->CR_X_sparse[m]).cast<complex<double>>(); //includes delta_alpha_X
	    CR_b[m] = (MyLeg->CR_B_sparse[m]).cast<complex<double>>(); //includes delta_alpha_B
	    Rho_X0[m] = (MyDM->BoltzM_X[m]).cast<complex<double>>(); 
	    Rho_B0[m] = (MyDM->BoltzM_B[m]).cast<complex<double>>();
	    Rho_N20[m] = (MyDM->BoltzM_0[m]).cast<complex<double>>();
	  //	  std::cerr << " m = " << m << std::endl;
	  //	  std::cerr << CR_theta[m];
	}

	P0_peak = DenVec(Jmax1);
	PX_peak = DenVec(Jmax1);
	PB_peak = DenVec(Jmax1);
	PBs_peak = DenVec(Jmax1);	
	Pa_X = DenVec(Jmax1);
	Pa_B = DenVec(Jmax1);
	
}

simul_denmat::~simul_denmat()
{

}

void simul_denmat::Hamiltonian0()
{
	for (size_t j = 0; j < Jmax1; j++)
	{
		valH = Bo * j*(j + 1) - Do * j*j*(j + 1)*(j + 1);
		tripletH_0.push_back(T(j, j, valH));
		valH = Bo_X * j*(j + 1) - Do_X * j*j*(j + 1)*(j + 1) + Eng_X;
		tripletH_X.push_back(T(j, j, valH));
		valH = Bo_B * j*(j + 1) - Do_B * j*j*(j + 1)*(j + 1) + Eng_B;
		tripletH_B.push_back(T(j, j, valH));
	}

	Hamiltonian_0 = SpMat(Jmax1, Jmax1);
	Hamiltonian_X = SpMat(Jmax1, Jmax1);
	Hamiltonian_B = SpMat(Jmax1, Jmax1);
	Hamiltonian_0.setFromTriplets(tripletH_0.begin(), tripletH_0.end());
	Hamiltonian_X.setFromTriplets(tripletH_X.begin(), tripletH_X.end());
	Hamiltonian_B.setFromTriplets(tripletH_B.begin(), tripletH_B.end());

	H0_o = SpMatComp(Jmax1, Jmax1);
	H0_x = SpMatComp(Jmax1, Jmax1);
	H0_b = SpMatComp(Jmax1, Jmax1);
	H0_o = Hamiltonian_0.cast<complex<double>>();
	H0_x = Hamiltonian_X.cast<complex<double>>();
	H0_b = Hamiltonian_B.cast<complex<double>>();
	H0 = SpMatComp(Jmax_XB2, Jmax_XB2);
	H0.rightCols(Jmax1) = H0_b;
	H0 = Permute * H0;
	H0.leftCols(Jmax1) = H0_x;
}

void simul_denmat::Exp_dHdt()
{
	ExpHdt_0 = DenMatComp(Jmax1, Jmax1);
	ExpHdt_X = DenMatComp(Jmax1, Jmax1);
	ExpHdt_B = DenMatComp(Jmax1, Jmax1);
	for (size_t j = 0; j < Jmax1; j++)
	{
		for (size_t i = 0; i <= j; i++)
		{
			delta_Eij = Bo * (j*(j + 1) - i * (i + 1)) - Do * (j*j*(j + 1)*(j + 1) -
				i * i*(i + 1)*(i + 1));
			ExpHdt_0(i, j) = cos(dTp*delta_Eij) + ImagUnit * sin(dTp*delta_Eij);
			ExpHdt_0(j, i) = conj(ExpHdt_0(i, j));

			delta_Eij = Bo_X * (j*(j + 1) - i * (i + 1)) - Do_X * (j*j*(j + 1)*(j + 1) -
				i * i*(i + 1)*(i + 1));
			ExpHdt_X(i, j) = cos(dTp*delta_Eij) + ImagUnit * sin(dTp*delta_Eij);
			ExpHdt_X(j, i) = conj(ExpHdt_X(i, j));

			delta_Eij = Bo_B * (j*(j + 1) - i * (i + 1)) - Do_B * (j*j*(j + 1)*(j + 1) -
				i * i * (i + 1)*(i + 1));
			ExpHdt_B(i, j) = cos(dTp*delta_Eij) + ImagUnit * sin(dTp*delta_Eij);
			ExpHdt_B(j,i) = conj(ExpHdt_B(i, j));
		}
	}
}

void simul_denmat::sym2asym()
{
  //  Aodd = Rho_N2;
  Aodd = Rho_theta;  
  for (size_t i = 0; i < Jmax1; i += 2)
    {
      for(size_t j = i; j < Jmax1; j++)
	{   Aodd(i,j) = 0.0;
	  Aodd(j,i) = 0.0;
	}
    }
  //    Aeven = Rho_N2 - Aodd;
    Aeven = Rho_theta - Aodd;
    odd_even = Aodd.trace().real()/Aeven.trace().real();
  if (odd_even > 1.0e-10)
    {
      Rho_N2_asym = odd_even*Aeven + 1./odd_even*Aodd;
    }
  else
    {
      Rho_N2_asym = Aeven;
    } 
}

void simul_denmat::Permutation()
{
	for (size_t j = 0; j < Jmax1; j++)
	{
		tripletP.push_back(T(j + Jmax1, j, 1));
	}

	PermuteD = SpMat(Jmax_XB2, Jmax_XB2);
	PermuteD.setFromTriplets(tripletP.begin(), tripletP.end());
	Permute = SpMatComp(Jmax_XB2, Jmax_XB2);
	Permute = PermuteD.cast<complex<double>>();
}

void simul_denmat::Runge_Kutta_pump(size_t probe)
{
        Exp_dHdt(); // filling dense matrix for field-free propagation
	
	for (size_t m = 0; m < Jmax01; m++)
	{
		std::cerr << " m=" << m << std::endl;

		Rho_N2 = Rho_N20[m];
		Rho_X = Rho_X0[m];
		Rho_B = Rho_B0[m]; 

		mult = 2;
		if (m == 0) mult = 1;
		Cos2_T_N2[0] += mult * (Rho_N2*CR_o[m]).trace().real();
		Cos2X_T[0] += mult * (Rho_X*CR_x[m]).trace().real();
		Cos2B_T[0] += mult * (Rho_B*CR_b[m]).trace().real();
		H_mat_0 = H0_o - CR_o[m] * Potential[0];   // H[0]     
		for (size_t i = 1; i < nbT; i++)
		{
			if (i <= start_ff) {
			  //==============  for neutral N2 - ground state =================
				K1_o = Rho_N2 * H_mat_0 * ImagUnit;
				K1_o += DenMatComp(K1_o.transpose()).conjugate();
				H_mat_0 = H0_o - CR_o[m] * Potential[2 * i - 1];  // H[i-1/2]
				K2_o = (Rho_N2 + K1_o * dT_2)*H_mat_0 * ImagUnit;
				K2_o += DenMatComp(K2_o.transpose()).conjugate();
				K3_o = (Rho_N2 + K2_o * dT_2)*H_mat_0 * ImagUnit;
				K3_o += DenMatComp(K3_o.transpose()).conjugate();
				H_mat_0 = H0_o - CR_o[m] * Potential[2 * i];   // H[i]
				K4_o = (Rho_N2 + K3_o * dTp)*H_mat_0 * ImagUnit;
				K4_o += DenMatComp(K4_o.transpose()).conjugate();
				Rho_N2 += (K1_o + 2. * K2_o + 2. * K3_o + K4_o) * dT_6;			
				if (i == pump_peak)  {
				  Rho_theta = CR_theta[m] * Rho_N2 * CR_theta[m];
				  NormC = Rho_theta.trace().real();
				  Rho_theta *= Rho_N2.trace().real()/NormC;
				  Rho_X = Rho_theta * P_X;
				  sym2asym();
				  Rho_B = Rho_N2_asym * (1.-P_X);
				  Rho_Bs = Rho_theta * (1.-P_X);
				  H_mat_X = H0_x - CR_x[m] * Potential[2 * i];   // H[0]
				  H_mat_B = H0_b - CR_b[m] * Potential[2 * i];
				  P0_peak += mult * (Rho_N2.diagonal()).real();
				  PX_peak += mult * (Rho_X.diagonal()).real();
				  PB_peak += mult * (Rho_B.diagonal()).real();
				  PBs_peak += mult * (Rho_Bs.diagonal()).real();				  
				}	       
				
				if (i > pump_peak)  {
			//==============  for level X =================
				  K1_x = Rho_X * H_mat_X * ImagUnit;
				  K1_x += DenMatComp(K1_x.transpose()).conjugate();
				  H_mat_X = H0_x - CR_x[m] * Potential[2 * i - 1];  // H[i-1/2]
				  K2_x = (Rho_X + K1_x * dT_2)*H_mat_X * ImagUnit;
				  K2_x += DenMatComp(K2_x.transpose()).conjugate();
				  K3_x = (Rho_X + K2_x * dT_2)*H_mat_X * ImagUnit;
				  K3_x += DenMatComp(K3_x.transpose()).conjugate();
				  H_mat_X = H0_x - CR_x[m] * Potential[2 * i];   // H[i]
				  K4_x = (Rho_X + K3_x * dTp)*H_mat_X * ImagUnit;
				  K4_x += DenMatComp(K4_x.transpose()).conjugate();
				  Rho_X += (K1_x + 2. * K2_x + 2. * K3_x + K4_x) * dT_6;
			//==============  for level B =================/	
				  K1_b = Rho_B * H_mat_B * ImagUnit;
				  K1_b += DenMatComp(K1_b.transpose()).conjugate();
				  H_mat_B = H0_b - CR_b[m] * Potential[2 * i - 1];  // H[i-1/2]
				  K2_b = (Rho_B + K1_b * dT_2)*H_mat_B * ImagUnit;
				  K2_b += DenMatComp(K2_b.transpose()).conjugate();
				  K3_b = (Rho_B + K2_b * dT_2)*H_mat_B * ImagUnit;
				  K3_b += DenMatComp(K3_b.transpose()).conjugate();
				  H_mat_B = H0_b - CR_b[m] * Potential[2 * i];   // H[i]
				  K4_b = (Rho_B + K3_b * dTp)*H_mat_B * ImagUnit;
				  K4_b += DenMatComp(K4_b.transpose()).conjugate();
				  Rho_B += (K1_b + 2. * K2_b + 2. * K3_b + K4_b) * dT_6;
				}
			//===============================================*/
			}
			else {
			        Rho_N2 = Rho_N2.cwiseProduct(ExpHdt_0);
			  	Rho_X = Rho_X.cwiseProduct(ExpHdt_X);
			  	Rho_B = Rho_B.cwiseProduct(ExpHdt_B);
			}
			//===========================================================		
			Cos2_T_N2[i] += mult * (Rho_N2*CR_o[m]).trace().real();   // BN: it's already multiplied by delta_alpha because of CR_o
			if (i >= pump_peak) {
			  Cos2X_T[i] += mult * (Rho_X*CR_x[m]).trace().real();      // BN: it's already multiplied by delta_alpha_X because of CR_x and P_X due to Rho_X
			  Cos2B_T[i] += mult * (Rho_B*CR_b[m]).trace().real();      // BN: it's already multiplied by delta_alpha_B because of CR_b and (1-P_X) due to Rho_B
			}
			if (i == pump_peak + i_del0 * (probe - 1) - seed_start * i_sigma) {      
			  AUX.topLeftCorner(Jmax1, Jmax1) = Rho_X;
			  AUX.bottomRightCorner(Jmax1, Jmax1) = Rho_B;
			  Rho_in[m] = AUX.sparseView();
			  Pa_X += mult * (Rho_X.diagonal()).real();
			  Pa_B += mult * (Rho_B.diagonal()).real();			  
			}
		}
	}

	P_X_after = Pa_X.sum(); // P_X after passing of the pump pulses
	P_B_after = Pa_B.sum();  // P_B after passing of the pump pulses

	std::cerr << "At peak: P_X0 =" << PX_peak.sum() << "  P_B0 =" << PB_peak.sum() <<
	  "  P_Bs0 =" << PBs_peak.sum() << std::endl;
	std::cerr << "After pumping: P_X = " << P_X_after << ", P_B = " << P_B_after << std::endl;

	P_T.resize(nbT);
	dP_theta_dT.resize(nbT);
	d2P_theta_dT.resize(nbT);
	
	for (size_t i = 0; i < nbT; i++)
	  {
	    Cos2_T[i] = Cos2X_T[i] + Cos2B_T[i];  
	    P_T[i] = ioniz  * Cos2_T[i] + (1. - ioniz) * Cos2_T_N2[i];
	  }
	
	for (size_t i = 1; i < nbT - 1; i++)
	  {
	    dP_theta_dT[i] = (P_T[i + 1] - P_T[i - 1]) / 2. / dTp;
	    d2P_theta_dT[i] = (P_T[i + 1] - 2. * P_T[i] + P_T[i - 1]) / dTp / dTp;
	  }
	
}

void simul_denmat::Runge_Kutta_seed(vector<double>& Efield_m, vector<double>& E_half_m)
{	
	Efield = Efield_m;
	E_half = E_half_m;

	for (size_t i = 0; i < nbt; i++)
	{
		Uot[i] = Efield[i] * Efield[i] / 2.;
		Uo_half[i] = E_half[i] * E_half[i] / 2.;
		Pol[i] = 0.0;
		Cos2t[i] = 0.0;
		P_Xt[i] = 0.0;
		P_Bt[i] = 0.0;
		
	}	
	
	P_Xt[0] = P_X_after;
	P_Bt[0] = P_B_after;
	
	for (size_t m = 0; m < Jmax01; m++)
	{
	        std::cerr << " m=" << m << std::endl;	  
		mult = 2;
		if (m == 0) mult = 1;
		Rho = Rho_in[m];
		Pol[0] += mult * (Rho*CT_[m]).trace().real();
		Cos2t[0] += mult * (Rho*CR_[m]).trace().real();
		H_mat = H0 - CR_sum[m] * Uot[0] - CT_[m] * Efield[0];   // H[0]

		for (size_t i = 1; i < nbt; i++)
		{			
			K1 = (Rho * H_mat) * ImagUnit;
			K1 += DenMatComp(K1.transpose()).conjugate();
			H_mat = H0 - CR_sum[m] * Uo_half[i] - CT_[m] * E_half[i];  // H[i-1/2]
			K2 = ((Rho + K1 * dt_2)*H_mat) * ImagUnit;
			K2 += DenMatComp(K2.transpose()).conjugate();
			K3 = ((Rho + K2 * dt_2)*H_mat) * ImagUnit;
			K3 += DenMatComp(K3.transpose()).conjugate();
			H_mat = H0 - CR_sum[m] * Uot[i] - CT_[m] * Efield[i];   // H[i]
			K4 = ((Rho + K3 * dt)*H_mat) * ImagUnit;
			K4 += DenMatComp(K4.transpose()).conjugate();
			Rho += (K1 + K2 * 2. + K3 * 2. + K4) * dt_6;
//=====================================================
			Pol[i] += mult * (Rho*CT_[m]).trace().real();
			Cos2t[i] += mult * ((Rho.topLeftCorner(Jmax1, Jmax1)*CR_x[m]).trace().real()+(Rho.bottomRightCorner(Jmax1, Jmax1)*CR_b[m]).trace().real()); // delta_alpha_X&B and P_X are included
			P_Xt[i] += mult * (Rho.topLeftCorner(Jmax1, Jmax1)).trace().real();
			P_Bt[i] += mult * (Rho.bottomRightCorner(Jmax1, Jmax1)).trace().real();
		}
	}

 	for (size_t i = 1; i < nbt - 1; i++)
	{
		dPdt_mu = (Pol[i + 1] - Pol[i - 1]) / 2. / dt;
		d2Pdt_mu = (Pol[i + 1] - 2. * Pol[i] + Pol[i - 1]) / dt / dt;
		Pol_mu[i] = Nmol * ioniz * PI * dz / cel *
						(- dz * d2Pdt_mu / C_ref - 3. * dPdt_mu + dPdt_mu_prev[i]);
		dPdt_mu_prev[i] = dPdt_mu;
		dPdt_theta[i] = ioniz * (Cos2t[i + 1] - Cos2t[i - 1]) / 2. / dt	+ dPdt_theta1[i];
		d2Pdt_theta = ioniz * (Cos2t[i + 1] - 2.* Cos2t[i] + Cos2t[i - 1]) / dt / dt + d2Pdt_theta1[i];
		P_theta[i] = ioniz * Cos2t[i] + P_theta1[i];
		nonlin1[i] = Nmol * PI*dz / cel * (dz / C_ref * d2Pdt_theta + 3.*dPdt_theta[i]);
		nonlin2[i] = Nmol * PI*C_ref / cel * (3.*((1. - ioniz) * alpha_tr + ioniz * (P_Xt[i] * alpha_tr_X + (1.-P_Xt[i]) * alpha_tr_B) + P_theta[i]) + 2.*CFL*dt*dPdt_theta[i]);
		nonlin3[i] = Nmol * 2.*PI*C_ref / cel * ((1. - ioniz) * alpha_tr + ioniz * (P_Xt[i] * alpha_tr_X + (1.-P_Xt[i]) * alpha_tr_B) + P_theta[i]);
	}
}


void simul_denmat::Cos2_FT()
{
	size_tempo = nbT;
	omeg.resize(size_tempo);

	Cos2T = (fftw_complex*)fftw_malloc(size_tempo * sizeof(fftw_complex));
	COS2 = (fftw_complex*)fftw_malloc(size_tempo * sizeof(fftw_complex));

	for (size_t i = 0; i < nbT; i++)
	{	
		Cos2T[i][REAL] = P_T[i] * (a0 - a1 * cos(2. * PI* i / size_tempo) + 
							a2 * cos(4. * PI* i / size_tempo)) - a3 * cos(6. * PI* i / size_tempo);
		Cos2T[i][IMAG] = 0.0;
	}
		
	fftw_plan pCos2;
	
	pCos2 = fftw_plan_dft_1d(size_tempo, Cos2T, COS2, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(pCos2);

	std::fstream out_harmonics_Cos2("Cos2_harmonics", std::ios::out);
	{
	  for (size_t i = 0; i < size_tempo; i++)
	    {
	      omeg[i] = 2 * PI*i / dTp / Bo / nbT;
	      out_harmonics_Cos2 << omeg[i] << "  " << (COS2[i][REAL] * COS2[i][REAL] +
							COS2[i][IMAG] * COS2[i][IMAG]) << "  " << std::endl;
	    }
	  out_harmonics_Cos2.close();
	}
	
	
	fftw_destroy_plan(pCos2);
	fftw_free(Cos2T);
	fftw_free(COS2);
}

void simul_denmat::initial_data(size_t probe, int i_start, vector<double>& time_m)
{	
	size_t j1 = 0;
	size_t j1_in = 0;
	size_t j2 = 0;
	size_t j2_in = 0;
	size_t ip_min = 0;
	size_t ip_max = 0;
	size_t start = (size_t)(abs(i_start));
	tk::spline Spline1, Spline2, Spline3;
	
	time_p = time_m;   // for spline	

	if (i_start < 0) {
		ip_min = 0;
		ip_max = i_Tp - start;
		for (size_t j = 0; j < ConvC*start; j++)
		{
		  P_theta1[j] = delta_alpha / 3.;
			dPdt_theta1[j] = 0.0;
			d2Pdt_theta1[j] = 0.0;
		}
		j1_in = ConvC * start;
		j2_in = j1_in + ConvC;
	}
	else {
		j1_in = 0;
		j2_in = ConvC;
		ip_min = start;
		ip_max = start + i_Tp;
	}

	node.resize(ip_max - ip_min);
	Sp1.resize(ip_max - ip_min);
	Sp2.resize(ip_max - ip_min);
	Sp3.resize(ip_max - ip_min);

	j1 = j1_in;
	j2 = j2_in;
	size_t ind = 0;
	size_t ind2 = 0;
	for (size_t ip = ip_min; ip < ip_max; ip++)
	{
		ind = ip - ip_min;
		ind2 = (size_t)((j1 + j2) / 2. + 0.5);
		node[ind] = time_p[ind2];
		Sp1[ind] = P_T[ip];
		Sp2[ind] = dP_theta_dT[ip];
		Sp3[ind] = d2P_theta_dT[ip];
		j2 = j2 + ConvC;
		j1 = j1 + ConvC;
	}
	Spline1.set_points(node, Sp1);
	Spline2.set_points(node, Sp2);
	Spline3.set_points(node, Sp3);

	for (size_t j = j1_in; j < nbt; j++)
	{
		P_theta1[j] = Spline1(time_p[j]);
		dPdt_theta1[j] = Spline2(time_p[j]);
		d2Pdt_theta1[j] = Spline3(time_p[j]);
	}

}

void simul_denmat::X_B_population()
{
     std::fstream outf1("XB_population", std::ios::out | std::ios::app);  // to add   ios::out|ios::app
     {
       outf1 <<  P0_peak;
       outf1 << std::endl;
       outf1 <<  PX_peak;
       outf1 << std::endl;
       outf1 <<  PB_peak;
       outf1 << std::endl;
       outf1 <<  PBs_peak;
       outf1 << std::endl;
       outf1 <<  Pa_X;
       outf1 << std::endl;
       outf1 <<  Pa_B;
       outf1 << std::endl;
       outf1.close();       
     }
}

void simul_denmat::show_polarization()
{
       std::fstream out_polarization("Polarization", std::ios::out);
       {
               for (size_t i = 0; i < nbt; i++)
		 {
		   out_polarization << time_p[i] << "  " << Pol[i]  << "  " << std::endl;
		 }
	       out_polarization.close();
       }
  
}

void simul_denmat::cut()
{
	std::fstream out_Cos2_T("Cos2_T", std::ios::out);
	{
		for (size_t i = 1; i < nbT - 1; i++)
		{
		  out_Cos2_T << Tscale[i] << "  " << Cos2X_T[i]/delta_alpha_X/P_X << "  " << Cos2B_T[i]/delta_alpha_B/(1.-P_X) << "  "
			     << Cos2_T_N2[i]/delta_alpha << "  " << P_T[i] << "  " << dP_theta_dT[i] << "  "<< d2P_theta_dT[i] << "   " << std::endl;
		}
		out_Cos2_T.close();
	}


	std::fstream out_Rho_in("Rho_in", std::ios::out);
	{
	
                for (size_t m = 0; m < Jmax01; m++)
		{
			out_Rho_in << "			m=" << m << std::endl;
			out_Rho_in << Rho_in[m];
		}
		out_Rho_in.close();
	}

}

void simul_denmat::X_B_InTime()
{
  std::fstream out_pX_pB("XB_seed", std::ios::out);
  {
    for (size_t i = 0; i < nbt - 1; i++)
      {
	out_pX_pB << i*dt << "  " << P_Xt[i] << "  " << P_Bt[i]  << std::endl;
      }
    out_pX_pB.close();
  }

}

void simul_denmat::show_H()
{
	std::fstream out_ham("HamiltonianS", std::ios::out);
	{
		out_ham << "Identity" << std::endl;
		out_ham << IDNT_tr;
		out_ham << "state X" << std::endl;
		out_ham << H0_x;
		out_ham << "state B" << std::endl;
		out_ham << H0_b;
		out_ham << "X+B" << std::endl;
		out_ham << H0;
		out_ham << "Matrix Ext(i*dH_X*dt)" << std::endl;
		out_ham << ExpHdt_X<< std::endl;
		out_ham << "Matrix Ext(i*dH_B*dt)" << std::endl;
		out_ham << ExpHdt_B<< std::endl;
				
		out_ham.close();
	}
}

