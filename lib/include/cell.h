// IMPETUS
// version 1.0
// Jul 2016
// Author: Vi Q. Ha




#ifndef VMD_CELLSPACE_H
#define VMD_CELLSPACE_H


#include "gridspace.h"
#include "defines.h"
#include <stdlib.h> 
namespace vamde {

class CellSpace : public GridSpace {
public:
	
	
	struct Cinit{
		double minimum_cell_size;
		double sigma;
		double mass;
		int particle_distribution[3];
		double position_shift[3];
		int velocity_distribution_type ;
		double orderly_velocity_distribution [3];
		double uniformly_randomly_velocity_distribution_lower[3];
		double uniformly_randomly_velocity_distribution_upper[3];
		
		char *output_directory_name;
		void readinput(char * inputfilename );
		Cinit();
	};
	
	Cinit input;

	CellSpace(char * inputfilename ,vamde::Parameters *_param);
	CellSpace(Cinit _cinit,vamde::Parameters *_param);
	void refreshCell();

	/// nc is the 3 components representing number of cells in the whole simulation
	int nc[DIM];	 // number of cells in simulation domain
	
	int n_total_cell;
	int n_local_cell;
	int n_nonlocal_cell;
	int * nonlocal_cell_list;
	
	int n_total_global_particles;
	Domain shell;
	
	Cell *grid ;
	void initCell();
	
	void deleteWholeCell( int c);
	
	void deleteBorders();
		
	void loopCells_applyForce();
	
	void move_particles_to_the_correct_cells();
	
	void setup();

	void setCommunication(
					int d,
					int *lowericstart, int *lowericstop,
					int *lowericstartreceive, int *lowericstopreceive,
					int *uppericstart, int *uppericstop,
					int *uppericstartreceive, int *uppericstopreceive) ;
	void transfer(
					int dd,
					int *ic_number,											/// vector with the number of cells in each direction
					int sendproc, 											/// the 'lower' processor number along this direction [d]
					int receiveproc, 
					int *sendstart, int *sendstop,					/// vector with the cells range to send to the lower processor
					int *receivestart, int *receivestop,
					bool isCopyNotMove,
					bool sendUpNotDown
					);
					
	void sendReceiveCell(
					int d,
					int *ic_number,											/// vector with the number of cells in each direction
					int lowerproc, 											/// the 'lower' processor number along this direction [d]
					int *lowericstart, int *lowericstop,					/// vector with the cells range to send to the lower processor
					int *lowericstartreceive, int *lowericstopreceive,		/// vector with the cells range to recieve from the lower processor
					int upperproc, 
					int *uppericstart, int *uppericstop,
					int *uppericstartreceive, int *uppericstopreceive,
					bool isCopyNotMove
					);
					
	void transfer(int dd,bool isCopyNotMove,bool sendUpNotDown);
	void applyBoundaryConditions(Particle * p);
	bool is_in_box(ParticleList *pl) ;
	bool particle_is_in_processor_space( Particle * p ) ;
	void addParticles(Particle *p) ;
	void moveParticlesToCell(ParticleList *pl,Cell * grid);
	ParticleList * createParticle(int i, double x, double y, double z) ;
	void initParticles();
	int get_n_local();
	int get_n_global();
	int reset_n_total();
	void reset_all_id();
	
	
class AtomEyesCfgSegmentWriter_CellSpace{
	
	public:
		double atomic_length;
        AtomEyesCfgSegmentWriter_CellSpace(vamde::CellSpace *_s) ;
		void set_atom_name( char * _atom_name);
		void print();

	private:
		int entry_count;
		int tick ;
		int rank ;
		std::ofstream out;
		double l[3];
		double lo[3];
		int cell_print_interval ;
		vamde::CellSpace * s;
		char *output_directory_name;
		char *atom_name;
		void print_part(Particle *i)  ;
		void printInfo();
		void printHeader();
		void printType();
		void printBody ();
		void loop();
};
	AtomEyesCfgSegmentWriter_CellSpace * cfgwriter;
};

} // end namespace
#endif

