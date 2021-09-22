// IMPETUS
// version 1.0
// Jul 2016
// Author: Vi Q. Ha




#ifndef VMD_DEFINE_H
#define VMD_DEFINE_H


#include "defines.h"
#include <fstream>
namespace vamde {
	
class Clock
{
    private:
		double dt;
		
    public:
		
		int step;
		double t;
		Clock();
		Clock(double _dt) ;
		void set_dt(double _dt);
		
		void set_step(int _step);
		
		
		int getStep();
		double getTime();
		void advance();
};


class Shell{
public:
	real shell_lo[DIM];
	real shell_hi[DIM];
	real total_l[DIM];
};

class Domain{
public:
		
	real lo[DIM]; 		 		// lower coordinate of the world domain
	real hi[DIM]; 	 			// upper coordinate of the world domain
	real l[DIM]; 	 // size of simulation domain
		
	
	void printvars(){
		std::cout << "lo = " << lo[0]<< " " << lo[1] << " " << lo[2] << std::endl;
		std::cout << "hi = " << hi[0]<< " " << hi[1] << " " << hi[2] << std::endl;
		std::cout << "l = " << l[0]<< " " << l[1] << " " << l[2] << std::endl;
	}
};

class Interprocess {
public:
	/// world.np is the 3 components representing number of processors for the simulation

	int np[DIM];	 // 3 components number of subdomains 
	int np_total;	 // total number of subdomains
	
	int numprocs;	 // number of processes provided
	int myrank; 
	
	/// ip the 3 components index id from rank number
	int ip[DIM];	 // position of process in the process mesh
	/// determining the 6 neighor process rank number and store in ip_lower and ip_upper
	int ip_lower[DIM];	 // process number of the neighbor processes
	int ip_upper[DIM];
	
	void printvars(){
		std::cout << "np = " << np[0]<< " " << np[1] << " " << np[2] << std::endl;
		std::cout << "np_total = " << np_total << std::endl;
		std::cout << "numprocs = " << numprocs << std::endl;
		std::cout << "myrank = " << myrank << std::endl;
		std::cout << "ip = " << ip[0]<< " " << ip[1] << " " << ip[2] << std::endl;
		std::cout << "ip_lower = " << ip_lower[0]<< " " << ip_lower[1] << " " << ip_lower[2] << std::endl;
		std::cout << "ip_upper = " << ip_upper[0]<< " " << ip_upper[1] << " " << ip_upper[2] << std::endl;
	}
} ;

struct ClusterInfo{
	MPI_Comm mpicomm;
	Domain world;
	Domain home;
	Interprocess proc;
	ClusterInfo(MPI_Comm _mpicomm);

	void AssignSpace();

	void printvars(){
		std::cout << "world " << std::endl; 
		world . printvars();
		std::cout << "home " << std::endl;
		home . printvars();

		proc . printvars();
	}
} ;

struct Parameters {
	
	TransformationMatrix transformation_matrix;
	TransformationMatrix original_matrix;
	int n_net;
	int n_cell;
	int n_cont;
	
	double delta_t;
	int end_step;
	double desired_temperature;
	int cell_print_interval;
	int cont_print_interval;
	int net_print_interval;
	ClusterInfo * gridinfo;
	Clock * clock;
	int rank(){return gridinfo->proc.myrank;}
	double time(){return clock->getTime();}
	int step(){return clock->getStep();}
    Parameters(MPI_Comm mpicomm);
	void readinput(char * inputfilename ) ;
	void buildParams();
	void initMatrix();
	void viewinfo(){
		std::cout << "n_cell = " << n_cell << std::endl;
		std::cout << "n_cont = " << n_cont << std::endl;
		std::cout << "delta_t = " << delta_t << std::endl;
		std::cout << "end_step = " << end_step << std::endl;
		std::cout << "desired_temperature = " << desired_temperature << std::endl;
		std::cout << "cell_print_interval = " << cell_print_interval << std::endl;
		std::cout << "cont_print_interval = " << cont_print_interval << std::endl;
		gridinfo -> printvars();
		
	}
	
	
};


#include <time.h>
class TimeKeeper
{
    private:
	Parameters * param;
	public:
	clock_t t0,t1;
	
	TimeKeeper(Parameters * _param);
	void initTime();
	float getTime();
	void print_progress(int intv);
	void print_final();
};
}
#endif
