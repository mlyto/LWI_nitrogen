///////////////////////////////////////////////////////////////////////////////////                    
//                                                                               //                    
//  M.Lytova, 2018-2019, Project: N2+ lasing in rotationaly excited medium       //                    
//                                                                               //                    
/////////////////////////////////////////////////////////////////////////////////// 
#include <fstream>
#include <iostream>
#include "initial.h"
#include "Global.h"

using std::size_t;

initial::initial()
{
  
        E_term = kB * Temperature / AU_En;
	
	BoltzM_0.resize(Jmax01);
	BoltzM_X.resize(Jmax01);
	BoltzM_B.resize(Jmax01);
	tripletBoltz_0.resize(Jmax01);
	tripletBoltz_X.resize(Jmax01);
	tripletBoltz_B.resize(Jmax01);

	Norm_0 = 0;
	for (size_t j = 0; j < Jmax01; j++)
	{
		Rot_0 = (2 - j % 2) * exp(-(Bo * j* (j + 1) - Do * j * j * (j + 1) * (j + 1)) / E_term);
		Norm_0 += (2 * j + 1) * Rot_0;

		for (size_t m = 0; m <= j; m++)
		  {
		    tripletBoltz_0[m].push_back(T(j, j, Rot_0));
		  }
	}

	DenMat ToTrace;
	double summa = 0.0;
	int mult = 0;
	
		for (size_t m = 0; m < Jmax01; m++)
	{
		BoltzM_0[m] = SpMat(Jmax1, Jmax1);
		BoltzM_X[m] = SpMat(Jmax1, Jmax1);
		BoltzM_B[m] = SpMat(Jmax1, Jmax1);
		BoltzM_0[m].setFromTriplets(tripletBoltz_0[m].begin(), tripletBoltz_0[m].end());
		BoltzM_0[m] *= 1. / Norm_0;
		ToTrace = BoltzM_0[m];
		mult = 2;
		if (m == 0) mult = 1;
		summa += mult * ToTrace.trace();
	}
		std::cerr << " summa = " << summa  << std::endl;	
			
}
	
initial::~initial()
{  }

void initial::cut()
{
	std::fstream initial_dens("initial_DM", std::ios::out);
	{
		for (size_t m = 0; m < Jmax01; m++)
		{
			initial_dens << "  m=" << m << std::endl;
			initial_dens << "Neutral N2" << std::endl;
			initial_dens << BoltzM_0[m];
			initial_dens << "state X" << std::endl;
			initial_dens << BoltzM_X[m];
			initial_dens << "state B" << std::endl;
			initial_dens << BoltzM_B[m];
			initial_dens << std::endl;
		}
		initial_dens.close();
	}
}

