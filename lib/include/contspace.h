// IMPETUS
// version 1.0
// Jul 2016
// Author: Vi Q. Ha




#ifndef VMD_CONTSPACE_H
#define VMD_CONTSPACE_H


#include "gridspace.h"
#include "defines.h"
#include "Continuum.h"


namespace vamde {

class ContSpace : public GridSpace , public FiniteDifference {
public:
	struct Cinit{
		int nx, ny, nz , ghost_layer;
		char *output_directory_name;
		void readinput(char * inputfilename ) ;
	};
	
	Cinit input;
	int n_real[3];
	int n_real_total;
	int n_real_total_global;
	int ghost_layer;
	double lo_coordinate[3];
	ContinuumNode * Globalboundary[3][2];

	ContSpace(Cinit _cinit,vamde::Parameters *_param);
	ContSpace(int _nx, int _ny, int _nz, int _ghost_layer,vamde::Parameters *_param);
	virtual void setNodePositions();
	void getGap();
	void setup() ;
	void setupCont();
	void ERROR_DIM(){
		std::cout << " ***** Error ***** \n Finite Difference continuum will only work with 3D \n Please change DIM to 3 in your parameters.\n ";
		std::cout << "Please contact Vi Q. Ha if you have any further difficulties.\n For more information, please check out the VAMDE offcial website:\n http://engr.uconn.edu/~vqh06001/vamde.html \n *****************\n";
	}
	
	void transfer(int dd, bool isCopyNotMove, bool sendUpNotDown );
	void setGlobalBoundary(int l, int m, ContinuumNode * _node);
	void setGlobalSourceValue(double val);
	void setGlobalBoundarySource(int l, int m);
	void setGlobalBoundarySink(int l, int m);
	void setGlobalBoundaryInsulated(int l, int m);
	void linkAllNodeGlobalBoundary();
	void linkNodeGlobalBoundary(int l, int m);
	void setAllGlobalBoundaryInsulated();
	void setAllGlobalBoundarySink();
	void setAllGlobalBoundarySource();
	
	
class AtomEyesCfgSegmentWriter_ContSpace{
	
	public:
		double atomic_length;
        AtomEyesCfgSegmentWriter_ContSpace(vamde::ContSpace * _cs);
		void set_atom_name( char * _atom_name);
		void print();

	private:
		
		vamde::ClusterInfo * space;
		vamde::ContSpace * cs;
		int cont_print_interval ;
		int entry_count;
		int tick ;
		int rank ;
		std::ofstream out;
		double l[3];
		char *output_directory_name;
		char *atom_name;
		void printInfo();
		void printHeader();
		void printType();
		void printBody ();
		void loop();
};

	AtomEyesCfgSegmentWriter_ContSpace * cfgwriter;
};

} // end namespace
#endif

