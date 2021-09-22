/// Select one of the following for the tutorial:
//~ #include "../Tutorials/tutorialSimpleMD/simple.h"
//~ #include "../Tutorials/tutorialSimpleMD/liquid.h"
//~ #include "../Tutorials/tutorialHeat/heat.h"
#include "../Tutorials/tutorialChemotaxis/chemotaxis.h"

int main(int argc, char *argv[]) {
	
	MPI_Init(&argc, &argv);
	
	MPI_Comm mpicomm = MPI_COMM_WORLD;
	vamde::Parameters * param = new vamde::Parameters(mpicomm);
	
	init_input(param);
	
	runSimulation(param);
	
	MPI_Finalize();
    
	return 0;
}
