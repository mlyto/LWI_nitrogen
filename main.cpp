///////////////////////////////////////////////////////////////////////////////////
//                                                                               //
//  M.Lytova, 2018-2019, Project: N2+ lasing in rotationaly excited medium       //
//                                                                               //
///////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <stdio.h>
#include <iostream>
#include <chrono>
#include <fstream>
#include <string>
#include "Global.h"
#include "mesh.h"
#include "legendre_matrices.h"
#include "initial.h"
#include "simul_denmat.h"
#include "simul_maxwell.h"

using namespace std::chrono;
using std::cout;
using std::cerr;
using std::endl;


int main(int argc, char **argv) 
{
  
        mesh *my_mesh = new mesh();   // grids for the pump and seed pulses
	initial *init_DM = new initial();    // initial DM <- Boltzmann distribution 
	legendre_matrices *my_leg = new legendre_matrices();  // sparse Legendre matrices
	simul_maxwell *my_maxw = new simul_maxwell(my_mesh);
	simul_denmat *my_denmat = new simul_denmat(my_mesh, my_leg, init_DM, my_maxw->C_ref, my_maxw->dz);
	
	int ext_n = std::stoi(argv[1]);    // from submit file
  	size_t probe = (size_t)ext_n;
	my_denmat->Runge_Kutta_pump(probe);

	if (ext_n == 401) {
	  /*my_mesh->draw_pump();      
          my_leg->cut();                                                                                                                                           
          my_denmat->show_H(); */  
          init_DM->cut();        
	  my_denmat->cut();  
	  my_denmat->Cos2_FT();
	  my_denmat->X_B_population(); 
	}
	  
	
	my_maxw->norm_area(my_mesh, probe);  // initial spectral intensity 

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	my_denmat->initial_data(probe, my_mesh->i_start, my_mesh->time_probe);
	my_maxw->lax_wendroff(my_mesh, my_denmat);
	my_maxw->gain();
	       
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	auto duration = duration_cast<seconds>(t2 - t1).count();
	
	cerr << "one probe: " << duration << "  seconds" << endl;

	/*/	if ((ext_n == 1)||(ext_n == 50*(int)(ext_n/50)))
	my_maxw->cut_efield(argv[1]); */
	
	my_maxw->cut_spectra(argv[1]);

	if (ext_n == 401)
	my_denmat->X_B_InTime(); 
	
	delete my_denmat;
	delete my_maxw;
	delete init_DM;
	delete my_leg;
	delete my_mesh;
}  
