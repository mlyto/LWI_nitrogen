#pragma once
#include <cstddef>

using std::size_t;

const double PI = 3.1415926;
const double cel = 137.035999;  // AU
const double AU = 5.291772109e-9; // cm
const double AU_t = 2.4188843e-17; // s
const double AU_En = 4.35974418e-18; // J
const double AU_I_to_E = 5.338027e-9; // sqrt(intensity)* -> AU
const double En_cm_to_AU = 4.5562878e-6; // * for cm^{-1} -> AU 
const double kB = 1.38064852e-23; // J/K
//===== gas parameters ===========================
const double Temperature = 298; //  K for all
const double pressure = 13332.2; // Pa = 100 torr
const double Nmol = 5.0e18 * AU * AU * AU; // pressure / kB / Temperature * (1e-2*AU)*(1e-2*AU)*(1e-2*AU); number density of N2
const double ioniz = 0.001; // degree of ionization
const size_t Jmax0 = 30; // in initial distribution
const size_t Jmax01 = Jmax0 + 1; 
const size_t Jmax = 40; // J=0,1,...,Jmax for N2 and N2+
const size_t Jmax1 = Jmax + 1;
const size_t Jmax_XB2 = 2 * Jmax1;
const double P_X = 0.3; // relative population for X-state
const double Mu_XB = -0.740; // parallel
//===== N2 and N2+ rotational constants ===========================
const double Te_B = 0.11600946;  // AU
const double Bo = 9.0651031e-6; // AU N2
const double Do = 2.624e-11; // AU N2
const double al_x = 8.57038e-8;
const double al_b = 1.09350e-7;
const double Bo_X = 8.8016541e-6-al_x/2.; // AU  N2+
const double Do_X = 2.779e-11; // AU
const double Bo_B = 9.4522919e-6-al_b/2.; // AU
const double Do_B = 2.811e-11; // AU
const double omega_x = 1.005573e-2; // AU
const double omega_b = 1.102549e-2; // AU
const double omX_x = 7.335623e-5; // AU
const double omX_b = 1.0561475e-4; // AU
const double omY_x = -1.822515e-7; // AU
const double omY_b = -2.4467264e-6; // AU
const double delta_alpha = 4.349; // this is from Maria's paper 4.266; // 6.74833449;  for N2 in AU, or 1 A^3 
const double alpha_tr = 9.252; //10.6623691;  for N2 in AU, or 1.58 A^3
const double delta_alpha_X = 9.695; // 6.343434;  for N2+ in AU, or 0.94 A^3 
const double alpha_tr_X = 8.509; // 10.055018;  for N2+ in AU, or 1.49 A^3
const double delta_alpha_B = -4.68;  // 6.343434; for N2+ in AU, or 0.94 A^3
const double alpha_tr_B = 6.582; // 10.055018; for N2+ in AU, or 1.49 A^3 
//====== pump and seed pulse parameters =========================
const double lambda = 391.; // seed wavelength, nm
const double omega = 45.55 / lambda;  // seed frequency AU <-> 397 nm
const double band_left = lambda / 391.5;
const double band_right = lambda / 388.35;
const double seed_intensity = 10.0e10; // W/cm^2
const double Zmax = 0.05; // propagation length, cm (in experiment 300 - 500 micrometers)
const double T_probe0 = 250.0e-15;  // time scale for a probe, s
const double sigma = 20.0e-15; // Full-Width-Half-Max, s   seed
const int seed_start = 5; // time for the seed starts at to-seed_start*sigma
const double DelayInc = 10.0e-15; // time delay increment, s
const double CFL = 0.99;
const double pump_intensity = 1.0e14; // W/cm^2
const double Tau_on = 23.0e-15; // pump pulse duration 2*Tau_on=50 fs
//============ windows function coefficients ==================
const double a0 = 0.3635819;
const double a1 = 0.4891775;
const double a2 = 0.1365995;
const double a3 = 0.0106411;

inline size_t Min(size_t a, const size_t b)
{
  return a < b ? a : b;
}

