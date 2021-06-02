#pragma once
#include <vector>
#include <string>
#include <fftw3.h>
#include "simul_denmat.h"

using std::size_t;
using std::vector;

class simul_maxwell
{
public:
	simul_maxwell(mesh*);
	~simul_maxwell();

	double C_ref;
	double dz;
		
	void lax_wendroff(mesh*, simul_denmat*);
	void Field_FT();
	void norm_area(mesh*, size_t);	
	void gain();
	void cut_spectra(char* );
        void cut_efield(char* );
	
private:

	size_t nb_z;
	size_t i_sigma_dt;
	int i_start;
	size_t nbt;
	size_t nbt_zp;
	double dt;
	double dEdt_1;
	double spread;
	
	vector<double> time_p;
	vector<double> Efield; // k+1
	vector<double> E_half;
	vector<double> Efield_1; // k
	vector<double> dPdt_theta_prev; 
	vector<double> NL1, NL2, NL3;
	vector<double> P_Xt;
	vector<double> P_t_1, dPdt_1;
	vector<double> Pol_mu;  

	double d_omeg2, integral, integral0, integral_band;
	double t_del;
	vector<double> omeg2, IntSp, IntSp0, PolSp;
	size_t size_tempo2, s_min, s_max, s_mid, P_edge, R_edge;
	fftw_complex *Field_window, *Polar_window;
	fftw_complex *Field_spectra, *Polar_spectra;

	vector<double>::iterator pos;
	size_t it, sizeIntens;
};

