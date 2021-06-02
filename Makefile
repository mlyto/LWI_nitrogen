.PHONY : all clean

CXX 	 = g++
CXXFLAGS   = -Wall -O2 
LIBS    = -lfftw3
LDFLAGS   = $(LIBS)
CXXFLAGS  += -c -std=c++11  

SOURCES=main.cpp initial.cpp legendre_matrices.cpp mesh.cpp simul_denmat.cpp simul_maxwell.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXEC=azot

all: $(SOURCES) $(EXEC) 

$(EXEC): $(OBJECTS) 
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS)  $< -o $@

clean:
	rm -f $(EXEC) $(OBJECTS)
	rm -f mat_* Efield* elec_* Cos2* Intens* P_theta* pump_* Rho_* Ef_U_P Pol* XB* Hamilton* initial_DM current* std* slurm* *~
