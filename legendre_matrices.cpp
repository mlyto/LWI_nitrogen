///////////////////////////////////////////////////////////////////////////////////                    
//                                                                               //                    
//  M.Lytova, 2018-2019, Project: N2+ lasing in rotationaly excited medium       //                    
//                                                                               //                    
/////////////////////////////////////////////////////////////////////////////////// 
#include "legendre_matrices.h"
#include "Global.h"
#include <fstream>
#include <iostream>
#include <cmath>

using std::size_t;

legendre_matrices::legendre_matrices()
{
	tripletX.resize(Jmax01);
	tripletB.resize(Jmax01);
	triplet0.resize(Jmax01);
	tripletCR.resize(Jmax01);
	tripletCR_theta.resize(Jmax01);
	tripletCT.resize(Jmax01);	
	
	CR_X_sparse.resize(Jmax01);
	CR_B_sparse.resize(Jmax01);
	CR_0_sparse.resize(Jmax01);
	CR_theta_sparse.resize(Jmax01);
	CR_sparse.resize(Jmax01);	
	CT_sparse.resize(Jmax01);

	for (size_t J = 0; J < Jmax1; J++)
	{
	  for (size_t M = 0; M <= Min(J, Jmax0); M++)
		{
		  R_JJ = (2.*(J*J-M*M)+2.*J-1.)/(2.*J-1.)/(2.*J+3.)+0.5;
		  RX_JJ =delta_alpha_X*(2.*(J*J-M*M)+2.*J-1.)/(2.*J-1.)/(2.*J+3.);
 		  RB_JJ =delta_alpha_B*(2.*(J*J-M*M)+2.*J-1.)/(2.*J-1.)/(2.*J+3.);
		  R0_JJ =delta_alpha*(2.*(J*J-M*M)+2.*J-1.)/(2.*J-1.)/(2.*J+3.);
		  tripletX[M].push_back(T(J,J,RX_JJ));
		  tripletB[M].push_back(T(J,J,RB_JJ));
		  triplet0[M].push_back(T(J,J,R0_JJ));
		  tripletCR[M].push_back(T(J,J,RX_JJ));
		  tripletCR[M].push_back(T(J+Jmax1,J+Jmax1,RB_JJ));
		  tripletCR_theta[M].push_back(T(J,J,R_JJ));
		  if ((J > 0)&&(M!=J)) {
				        S_1 = sqrt((J*J-M*M)/(2.*J+1.)/(2.*J-1.))*Mu_XB;  // matrices CT include mu
					tripletCT[M].push_back(T(J,J-1+Jmax1, S_1));
					tripletCT[M].push_back(T(J-1+Jmax1,J, S_1));
					tripletCT[M].push_back(T(J+Jmax1,J-1, S_1));
					tripletCT[M].push_back(T(J-1,J+Jmax1, S_1));
					if(J < Jmax1-1) {
					  R_2 = 1./(2.*J+1.)*sqrt(((J+1.)*(J+1.)-M*M)*(J*J-M*M)/(2.*J+3.)/(2.*J-1.));
					  RX_2 = delta_alpha_X/(2.*J+1.)*sqrt(((J+1.)*(J+1.)-M*M)*(J*J-M*M)/(2.*J+3.)/(2.*J-1.));
					  RB_2 = delta_alpha_B/(2.*J+1.)*sqrt(((J+1.)*(J+1.)-M*M)*(J*J-M*M)/(2.*J+3.)/(2.*J-1.));
  					  R0_2 = delta_alpha/(2.*J+1.)*sqrt(((J+1.)*(J+1.)-M*M)*(J*J-M*M)/(2.*J+3.)/(2.*J-1.));	
					  tripletX[M].push_back(T(J+1, J-1, RX_2));
					  tripletX[M].push_back(T(J-1, J+1, RX_2));
					  tripletB[M].push_back(T(J+1, J-1, RB_2));
					  tripletB[M].push_back(T(J-1, J+1, RB_2));
					  triplet0[M].push_back(T(J+1, J-1, R0_2));
					  triplet0[M].push_back(T(J-1, J+1, R0_2));
					  tripletCR_theta[M].push_back(T(J+1,J-1,R_2));
					  tripletCR_theta[M].push_back(T(J-1,J+1,R_2));					  
					  tripletCR[M].push_back(T(J+1,J-1,RX_2));
					  tripletCR[M].push_back(T(J-1,J+1,RX_2));					  
					  tripletCR[M].push_back(T(J+1+Jmax1,J-1+Jmax1,RB_2));
					  tripletCR[M].push_back(T(J-1+Jmax1,J+1+Jmax1,RB_2));
        				}	
			       }
		}
	}

	for (size_t m = 0; m < Jmax01; m++)
	{
		CR_X_sparse[m] = SpMat(Jmax1, Jmax1);
		CR_X_sparse[m].setFromTriplets(tripletX[m].begin(), tripletX[m].end());
		CR_B_sparse[m] = SpMat(Jmax1, Jmax1);
		CR_B_sparse[m].setFromTriplets(tripletB[m].begin(), tripletB[m].end());
		CR_0_sparse[m] = SpMat(Jmax1, Jmax1);
		CR_0_sparse[m].setFromTriplets(triplet0[m].begin(), triplet0[m].end());
		CT_sparse[m] = SpMat(Jmax_XB2, Jmax_XB2);
		CT_sparse[m].setFromTriplets(tripletCT[m].begin(), tripletCT[m].end());
		CR_sparse[m] = SpMat(Jmax_XB2, Jmax_XB2);
		CR_sparse[m].setFromTriplets(tripletCR[m].begin(), tripletCR[m].end());
		CR_theta_sparse[m] = SpMat(Jmax1, Jmax1);
		CR_theta_sparse[m].setFromTriplets(tripletCR_theta[m].begin(), tripletCR_theta[m].end());
	}
}

legendre_matrices::~legendre_matrices()
{

}


void legendre_matrices::cut()    // shows all sparse matrices
{
	std::fstream out_CR_X_sparse("mat_CR_X_sparse", std::ios::out);
	{
		for (size_t m = 0; m < Jmax01; m++)
		{
			out_CR_X_sparse << "  m=" << m << std::endl;
			out_CR_X_sparse <<  CR_X_sparse[m];			
		}
		out_CR_X_sparse.close();
	}

	std::fstream out_CR_B_sparse("mat_CR_B_sparse", std::ios::out);
	{
		for (size_t m = 0; m < Jmax01; m++)
		{
			out_CR_B_sparse << "  m=" << m << std::endl;
			out_CR_B_sparse <<  CR_B_sparse[m];
		}
		out_CR_B_sparse.close();
		
	}

	std::fstream out_CR_sparse("mat_CR_sparse", std::ios::out);
	{
		for (size_t m = 0; m < Jmax01; m++)
		{
			out_CR_sparse << "  m=" << m << std::endl;
			out_CR_sparse <<  CR_0_sparse[m];
		}
		out_CR_sparse.close();
	}

	std::fstream out_CT_sparse("mat_CT_sparse", std::ios::out);
	{
		for (size_t m = 0; m < Jmax01; m++)
		{
			out_CT_sparse << "  m=" << m << std::endl;
			out_CT_sparse <<  CT_sparse[m];
		}
		out_CT_sparse.close();
	}

	std::fstream out_CR_theta_sparse("mat_CR_theta", std::ios::out);
	{
		for (size_t m = 0; m < Jmax01; m++)
		{
			out_CR_theta_sparse << "  m=" << m << std::endl;
			out_CR_theta_sparse <<  CR_theta_sparse[m];
		}
		out_CR_theta_sparse.close();
	}
}
	
