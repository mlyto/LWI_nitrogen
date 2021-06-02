#pragma once
#include<vector>

using std::size_t;

class mesh
{

 public:
	 mesh();
	 ~mesh();
	 size_t nb_T;
	 size_t nb_t;
	 size_t nb_seeds;
	 size_t i_del0;
	 size_t i_Tp;
	 size_t ConvC;
	 int i_start;
	 size_t i_sigma;
	 double dt;
	 double dTp;
	 double Eo;
	 double t_del;
	 
	 std::vector<double> Upt;
	 std::vector<double> Efield;
	 std::vector<double> time_probe;
	 std::vector<double> Tsc;


	 void draw_pump();
	 void seed(size_t);
  
 private:
  
	 size_t E_bs, C_bs;
	 double d_time, dTsc;
	 double pump_probe;
	 double T_probe, T_pump;
	 double Upump;
	 double to, t_to, t_to_sigma, t_start;
	 size_t N_on;
};

