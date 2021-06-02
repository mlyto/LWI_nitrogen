///////////////////////////////////////////////////////////////////////////////////                    
//                                                                               //                    
//  M.Lytova, 2018-2019, Project: N2+ lasing in rotationaly excited medium       //                    
//                                                                               //                    
/////////////////////////////////////////////////////////////////////////////////// 
#include "mesh.h"
#include "Global.h"
#include <cmath> 
#include <iostream>
#include <fstream> 

mesh::mesh() 
{  
// ==== SEED ===============

  E_bs = 80; // multiplyer to obtain better spectra of the seed field
  T_probe = E_bs * T_probe0;   // s, time scale for a probe
  nb_t = pow(2, 20); // 22 number of nodes
  dt = T_probe  / AU_t / (nb_t - 1);   // AU
  d_time = T_probe * 1e12 / (nb_t - 1);  // same in ps
  Eo = sqrt(seed_intensity) * AU_I_to_E; // AU
  Efield.resize(nb_t);
  time_probe.resize(nb_t);
  
// ==== PUMP ===============

  C_bs = 4;               // multiplyer to obtain better spectra of Cos2
  nb_T = 12500 * C_bs + 1;         // number of nodes 
  pump_probe = 50;        // duration_pump0/duration_probe0 (without multipliers)
  T_pump=pump_probe*T_probe0*C_bs;  // time scale for the pump, s 
  dTsc = T_pump * 1e12 / (nb_T - 1);  // ps
  dTp = T_pump / AU_t / (nb_T - 1) ;   // AU

  Tsc.resize(nb_T);
  for (size_t i = 0; i < nb_T; i++) Tsc[i] = i * dTsc; // for data output

  N_on = (int)(2 * 2 * Tau_on *1e12 / dTsc  + 0.49); // for pump we need intermid points
  Upump = AU_I_to_E*AU_I_to_E*pump_intensity / 4.;   // in AU
  
  Upt.resize(2 * nb_T);
  
  for (size_t i = 0; i <= N_on; i++)
    {
      Upt[i] = Upump * sin(PI*i / N_on)*sin(PI*i / N_on);
    } 

  //============================

  nb_seeds = (int)(T_probe0 * pump_probe / DelayInc + 0.49);
  i_del0 = (int)(DelayInc * 1e12 / dTsc + 0.49);
  ConvC = (int)(dTp / dt);
  i_Tp = (int)(T_probe * 1e12 / dTsc + 0.49);
  i_sigma = (int)(sigma * 1e12 / dTsc + 0.49);

  std::cerr << " E_bs = " << E_bs << ", nb_t = " << nb_t << ", nb_T = " << nb_T << ", nb_seeds = " << nb_seeds << std::endl;
  std::cerr << " Nmol = " << Nmol / AU / AU / AU << ", ioniz = " << ioniz << ", Temperature = " << Temperature << std::endl;
  std::cerr << " Ipump = " << pump_intensity << ",  Iseed = " << seed_intensity << std::endl;
  std::cerr << " Tau_on = " << Tau_on << ", sigma = " << sigma << ", lambda = " << lambda << std::endl;
  std::cerr << " Jmax0 = " << Jmax0 << ", Jmax = " << Jmax << std::endl;
  
}
mesh::~mesh()
{ 
	
}

void mesh::draw_pump()
{
	std::fstream outUp("pump_profile", std::ios::out);
	{
		for (size_t i = 0; i < nb_T; i++)
		outUp << Tsc[i] << "   " << Upt[2*i] << std::endl;
		outUp.close();
	}
}

void mesh::seed(size_t probe)
{
	t_del = DelayInc * (probe-1);
	to = t_del;   // max of the pump, s
	t_start = (to - seed_start * sigma)*1e12;  // start time for the probe to avoid BCs problems, ps
	i_start = (int)(t_start / dTsc+0.5);  // yes, can be <0 
	t_to = -seed_start * sigma / AU_t;	
	t_to_sigma = -seed_start;  //  -seed_start * sigma / sigma;

	std::cerr << " t_del = " << t_del*1e12 << ", ps" << std::endl;
	
	for (size_t i = 0; i < nb_t; i++)
	{		
		Efield[i] = Eo * exp(-4 * log(2.0)*(t_to_sigma*t_to_sigma))*cos(omega*t_to);
		t_to += dt;  // AU
		t_to_sigma += dt * (AU_t / sigma); 
		time_probe[i] = i * d_time + t_start;   // ps
	}

}
  
