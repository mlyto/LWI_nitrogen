///////////////////////////////////////////////////////////////////////////////////                    
//                                                                               //                    
//  M.Lytova, 2018-2019, Project: N2+ lasing in rotationaly excited medium       //                    
//                                                                               //                    
/////////////////////////////////////////////////////////////////////////////////// 
#include "Global.h"
#include "simul_maxwell.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#define REAL 0
#define IMAG 1

using std::cout;
using std::cerr;
using std::endl;

simul_maxwell::simul_maxwell(mesh* MyMesh)
{
	nbt = MyMesh->nb_t;
	this->dt = MyMesh->dt;
	this->i_sigma_dt = MyMesh->i_sigma*MyMesh->ConvC;
	
	C_ref = cel / 2. / PI / Nmol / ((1. - ioniz) * (alpha_tr + delta_alpha / 3.) + ioniz * (P_X * (alpha_tr_X + delta_alpha_X / 3.) + (1. - P_X) * (alpha_tr_B + delta_alpha_B / 3.)));	
	dz = CFL * C_ref * dt;
	
	nb_z = (int)(Zmax / AU / dz + 0.5);   // note: +0.5 as (int) truncates	

	dPdt_theta_prev.resize(nbt);
	E_half.resize(nbt);
	Pol_mu.resize(nbt);
	
// === want to draw the spectra in interval omega = [(omega - spread), (omega + spread)]
	nbt_zp = nbt;
	spread = 3.*2*sqrt(2.) / sigma * AU_t;
	s_min = (size_t)((omega - spread) * nbt_zp * dt / 2. / PI);
	s_max = (size_t)((omega + spread) * nbt_zp * dt / 2. / PI) + 1;
	s_mid = (size_t)((s_max + s_min) / 2);
	omeg2.resize(s_max - s_min);
	IntSp.resize(s_max - s_min);
	PolSp.resize(s_max - s_min);
	d_omeg2 = 2. * PI / dt / omega / (double)nbt_zp;
	for (size_t i = 0; i < s_max-s_min; i++)
	{
	    omeg2[i] = (i + s_min) * d_omeg2;

            if ((omeg2[i] < band_left)&&(omeg2[i]+d_omeg2 >= band_left)) P_edge = i;
	    if ((omeg2[i] < band_right)&&(omeg2[i]+d_omeg2 >= band_right)) R_edge = i;
	}
	cerr << " omeg2[0] = " << omeg2[0] << " , omeg2[s_max-s_min-1] = " << omeg2[s_max-s_min-1] << std::endl;
	cerr << " band_left = " << band_left << " band_right = " << band_right << std::endl;
	cerr << " C_ref=" << C_ref << ",  dz=" << dz << ", nb_z=" << nb_z << ", nbt_zp="<< nbt_zp <<endl;
	cerr << " Zmax = " << Zmax << ", actual Zmax = " << dz*nb_z*AU << endl;
	cerr << " s_min = " << s_min << ", s_max = " << s_max << ", s_mid = " << s_mid <<
	  ", d_omeg2 = " << d_omeg2 << ", P_edge = " << P_edge << ", R_edge = " << R_edge << endl;
}

simul_maxwell::~simul_maxwell()
{
	
}

void simul_maxwell::lax_wendroff(mesh* MyMesh, simul_denmat* MyDenmat)
{
	this->t_del = MyMesh->t_del;
	Efield_1 = Efield;
	time_p = MyMesh->time_probe;
	dPdt_1 = MyDenmat->dPdt_theta1;
	P_t_1 = MyDenmat->P_theta1;

	for (size_t k = 0; k < 1; k++) // nb_z; k++)
	{
		cerr << "	k = " << k << endl;
		Efield[0] = Efield_1[0];
		Efield[nbt - 1] = Efield_1[nbt - 1];
		E_half[1] = (Efield_1[1] + Efield_1[0]) / 2.;
		E_half[nbt - 1] = (Efield[nbt - 1] + Efield[nbt - 2]) / 2.; 
//===== Lagrange interp of 4th order at i - 1 / 2
		for (size_t i = 2; i < nbt - 1; i++)		
		E_half[i] = (-Efield_1[i + 1] + 9. * Efield_1[i] + 9. * Efield_1[i - 1] - Efield_1[i - 2]) / 16.;

		MyDenmat->Runge_Kutta_seed(Efield_1, E_half);

		//		if (k == nb_z-1)  {
		//		    MyDenmat->show_polarization();
		    //		  } 

		NL1 = MyDenmat->nonlin1;
		NL2 = MyDenmat->nonlin2;
		NL3 = MyDenmat->nonlin3;
		Pol_mu = MyDenmat->Pol_mu;
		P_Xt = MyDenmat->P_Xt;

		for (size_t i = 1; i < nbt - 1; i++)
		{	
			Efield[i] = Efield_1[i] * (1. - NL1[i]) + CFL * ((1. - NL2[i])*(Efield_1[i + 1] - Efield_1[i - 1]) +
				CFL * (1. - NL3[i])*(Efield_1[i + 1] - 2. * Efield_1[i] + Efield_1[i - 1])) / 2.
				+ Pol_mu[i] + Nmol * PI*dz / cel * dPdt_theta_prev[i];

			dEdt_1 = (Efield_1[i + 1] - Efield_1[i - 1]) / 2. / dt;
			dPdt_theta_prev[i] = dPdt_1[i] * Efield_1[i] + ((1. - ioniz) * alpha_tr + ioniz * P_Xt[i] * alpha_tr_X +
						ioniz * (1. - P_Xt[i]) * alpha_tr_B + P_t_1[i])*dEdt_1;
		}
		Efield_1 = Efield;		
		dPdt_1 = MyDenmat->dPdt_theta;
		P_t_1 = MyDenmat->P_theta;
	}
}

void simul_maxwell::Field_FT()
{
	size_tempo2 = nbt_zp;
	
	fftw_plan pField;
	fftw_plan pPolar;
	
	Field_window = (fftw_complex*)fftw_malloc(size_tempo2 * sizeof(fftw_complex));
	Field_spectra = (fftw_complex*)fftw_malloc(size_tempo2 * sizeof(fftw_complex));
	Polar_window = (fftw_complex*)fftw_malloc(size_tempo2 * sizeof(fftw_complex));
	Polar_spectra = (fftw_complex*)fftw_malloc(size_tempo2 * sizeof(fftw_complex));

	pField = fftw_plan_dft_1d(size_tempo2, Field_window, Field_spectra, FFTW_FORWARD, FFTW_ESTIMATE);
	pPolar = fftw_plan_dft_1d(size_tempo2, Polar_window, Polar_spectra, FFTW_FORWARD, FFTW_ESTIMATE);
	
	for (size_t i = 0; i < size_tempo2; i++)
	{
	     Field_window[i][REAL] = Efield[i];
	     Field_window[i][IMAG] = 0.0;
	     Polar_window[i][REAL] = Pol_mu[i];
	     Polar_window[i][IMAG] = 0.0;	     
	}

	fftw_execute(pField);
	fftw_execute(pPolar);
	
	size_t m = 0;
	for (size_t i = 0; i < s_max-s_min; i++)
	{
		m = i + s_min;
		IntSp[i] = (Field_spectra[m][REAL] * Field_spectra[m][REAL] +
			Field_spectra[m][IMAG] * Field_spectra[m][IMAG])*dt*dt;
		PolSp[i] = (Polar_spectra[m][REAL] * Polar_spectra[m][REAL] +
			    Polar_spectra[m][IMAG] * Polar_spectra[m][IMAG])*dt*dt;		
	}

	fftw_destroy_plan(pField);
	fftw_destroy_plan(pPolar);
	
	fftw_free(Field_window);
	fftw_free(Field_spectra);
	fftw_free(Polar_window);
	fftw_free(Polar_spectra);
}

void simul_maxwell::norm_area(mesh* MyMesh, size_t probe)
{
         
        MyMesh->seed(probe);
	Efield = MyMesh->Efield;
        Field_FT();
	IntSp0 = IntSp;

	integral0 = 0.0;

        for (size_t i = 1; i < s_max-s_min-1; i++)   // trapezoidal rule
        {
                integral0 += IntSp[i]*d_omeg2;		
        }

	integral0 += (IntSp[0] + IntSp[s_max-s_min-1])*d_omeg2 / 2.;
	
	integral_band = 0.0;

	for (size_t i = P_edge-1; i < R_edge; i++)   // trapezoidal rule
	  {
	    integral_band += IntSp[i]*d_omeg2;
	  }

	integral_band += (IntSp[P_edge] + IntSp[R_edge])*d_omeg2 / 2.;
	
	
	cerr << "integral0 = " << integral0 << ", integral_band = " << integral_band << endl;
}

void simul_maxwell::gain()
{
        Field_FT(); 
        integral = 0.0;
 
	for (size_t i = P_edge-1; i < R_edge; i++)  // trapezoidal rule
	{
	     integral += (IntSp[i] - IntSp0[i])*d_omeg2;
	}

             integral += (IntSp[P_edge] - IntSp0[P_edge] + IntSp[R_edge] - IntSp0[R_edge])*d_omeg2 / 2.;
}

void simul_maxwell::cut_spectra(char* number)
{
        char harmonics_name[19] = "elec_harmonics";

         strcat(harmonics_name, number);

	std::fstream out_harmonics_Efield(harmonics_name, std::ios::out);
	{
		for (size_t i = 0; i < s_max-s_min; i++)
		{
		  out_harmonics_Efield << omeg2[i] << "  " <<IntSp0[i] << "  " << IntSp[i] << "  " << PolSp[i]  << std::endl;
		}
		out_harmonics_Efield.close();
	}

}

void simul_maxwell::cut_efield(char* number)
{
        char efield_name[11] = "Efield";

        strcat(efield_name, number);

        std::fstream out_efield(efield_name, std::ios::out);
        {
                for (size_t i = 0; i < nbt; i++)
                {
                 out_efield << time_p[i] << "  " << Efield[i]  << std::endl;
                }
                out_efield.close();
        }

}



