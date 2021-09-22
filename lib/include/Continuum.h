// IMPETUS
// version 1.0
// Jul 2016
// Author: Vi Q. Ha

#ifndef VAMDE_CONT
#define VAMDE_CONT

#include <iostream>

namespace vamde {
	
	class ContinuumNode {
	private:
	public:
		int index[3];
		double position[3];
		bool active;
		
		double gradient[3];
		double concentration;
		double c_v;
		double c_a;
		double conductivity;
		class ContinuumNode * next[3][2];
		double neighborConcentration[3][2];
		ContinuumNode(){concentration = 0;c_v = 0;c_a = 0; conductivity = 1; active = 1;}
	
	};
	
	class Continuum {
	private:
	public:
		int time;
		double dt;
		double alacrity;
		double diffusivity;
		double elasticity;
		int distributions[3];
		int n_total;
		double size[3];
		double lo_coordinate[3];
		double gap[3];
		int nx(){return distributions[0];}
		int ny(){return distributions[1];}
		int nz(){return distributions[2];}
		
		Continuum(int _nx, int _ny, int _nz);
		Continuum();
		~Continuum();
		ContinuumNode ***node;
		ContinuumNode heavy, sink, source;
		ContinuumNode * boundary[3][2];
		
		void setSize(double _dx, double _dy, double _dz);
		void setDistributions(int _nx, int _ny, int _nz);
		virtual void getGap();
		
		void defaultSetup(int _nx, int _ny, int _nz);
		void BuildContinuum();
		void DestroyContinuum();
		void setDefaultBoundary();
		void setAllBoundarySource();
		void setAllBoundarySink();
		void setBoundarySink(int l, int m);
		
		void setSourceValue(double _boundary_source);
		void setBoundarySource(int l, int m);
		
		void setAllBoundaryInsulated();
		void setBoundaryInsulated(int l, int m);
		
		void setIndex();
		virtual void setNodePositions();
		void linkAllNodes();
		
		void setNodeTo(ContinuumNode * _node, ContinuumNode * _set, bool _active);
		void setNodeTo(ContinuumNode *, ContinuumNode *);
		void setNodeToHeavy(ContinuumNode *);
		void setNodeToSink(ContinuumNode *);
		
		void setToHeavy(int, int , int);
		void setToSink(int, int , int);
		
		void print();
		void print2();
		virtual void printAtomEyeCfg();
		virtual void Iterate(){ std::cout << "Continuum has no Iteration commands" ;};
		
		void UpdateBoundaryConditions();
	};

	class FiniteDifference : public Continuum {
	private:
	public:
		FiniteDifference();
		FiniteDifference(int _nx, int _ny, int _nz) ;
		virtual void Iterate();
		void IterateParabolicPDE();
		/// 3 components starting point of a 'local' cell
		int node_start[3];	 // width of border neighborhood, corresponds to the first local index in the interior of the subdomain
		/// 3 components ending point of a 'local' cell
		int node_stop[3]; 	 // first local index in the upper border neighborhood
		
		void FD_DefaultSetup();
		void setDefaultIteration();
		
		void parabolicPDE();
		void hyperbolicPDE();
		void IntegrateConcentration();
		void collectNeighborConcentration();
	};

}

#endif
