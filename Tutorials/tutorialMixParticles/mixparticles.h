#include "partList.h"

vamde::Parameters *  init_input ( vamde::Parameters * param){
	
		param->readinput("Tutorials/tutorialMixParticles/input/param.input" );
		return param;
}
void initParticles(vamde::CellSpace ** cellpointers) {
	vamde::Parameters * param = cellpointers[0]->param;
	
	int npart[3];
	int nside = 20;
	npart[0] = nside;
	npart[1] = nside;
	npart[2] = nside;
	int N = npart[0]*npart[1]*npart[2];
	int w=0;
	double x,y,z;
		
	for (int i=0; i<npart[0]; i++) {
		for (int j=0; j<npart[1]; j++) {
			for (int k=0; k<npart[2]; k++) {
				x= param->gridinfo->world.lo[0] + 
					(0.5+i)*(param->gridinfo->world.l[0] / (double) npart[0] ) ;
				y= param->gridinfo->world.lo[1] + 
					(0.5+j)*(param->gridinfo->world.l[1] / (double) npart[1] ) ;
				z= param->gridinfo->world.lo[2] + 
					(0.5+k)*(param->gridinfo->world.l[2] / (double) npart[2] ) ;
				
				double p = (100* ( (double)rand()/ RAND_MAX) )  ;
				
				int c = 0;
				
				if (p < 10) c = 1;
				cellpointers[c]->createParticle(w,x,y,z);
				
				w++;
			}
		}
	}
}


class Spring3 : public getParticlesAndNeighborsMulti {
	public:
	
		real sigma ;
		real epsilon;
        Spring3(vamde::CellSpace ** _cell_pointer_array, int _n_cellspaces) : getParticlesAndNeighborsMulti( _cell_pointer_array,_n_cellspaces ) {
		
			 sigma = 1.0;
			 epsilon = 5.0;
		}
		
		void NeighborAction(Particle *i, Particle *j) {
			real r = 0.0;
			for (int d=0; d<DIM; d++)
				r += sqr(j->x[d] - i->x[d]);
			real dx = sqrt(r);
			if (dx <= 1.1) {
				real f =  epsilon * (dx - sigma) ;
				for (int d=0; d<DIM; d++){
					i->F[d] += f * (j->x[d] - i->x[d]) ;
				}
			}
		}
};

class LeapFrog3: public getParticlesAndNeighborsMulti {
	
	public:
        LeapFrog3(vamde::CellSpace ** _cell_pointer_array, int _n_cellspaces) : getParticlesAndNeighborsMulti( _cell_pointer_array,_n_cellspaces ) {
			PART = 0;
			dt = param-> delta_t;
		}
		void step(int _step){
			PART = _step; 
			loopLocalParticles();
			PART = 0;
		}
	private:
		int PART;
		double dt;
		void ParticleAction(Particle *i)  {
			for (int d=0; d<DIM; d++) {
				if (PART == 1) {
					i->v[d] = i->v[d] + ( 0.5 * dt * i->F[d] / i->m );
					i->x[d] = i->x[d] + ( dt * i->v[d] );
					/// orientation
					i->nv[d] = i->nv[d] + ( 0.5 * dt * i->na[d]);
					i->nx[d] = i->nx[d] + ( dt * i->nv[d] );
				} else {
					i->v[d] = i->v[d] + ( 0.5 * dt * i->F[d] / i->m );
					/// orientation
					i->nv[d] = i->nv[d] + ( 0.5 * dt * i->na[d]);
				}
			}
		}
};

class ZeroForce3 : public getParticlesAndNeighborsMulti {
	
	public:
        ZeroForce3(vamde::CellSpace ** _cell_pointer_array, int _n_cellspaces) : getParticlesAndNeighborsMulti( _cell_pointer_array,_n_cellspaces ) {
		}
		void ParticleAction(Particle *i)  {
			for (int d=0; d<DIM; d++) {
				i->F[d] =0;
				i->na[d] =0;
			}
		}
};


void runSimulation(vamde::Parameters * param) {
	vamde::TimeKeeper * tk;
	tk = new vamde::TimeKeeper(param);
	std::cout << "New Vamde " << std::endl;
	
	
	std::cout << "Simulation Start " << std::endl;
	/// parameters: 
	double delta_t = param-> delta_t;
	double end_step = param-> end_step;
	double & t = param->clock->t;
	int & step = param->clock->step;
	
	/// Initiate clock
	param->clock->set_step(0);
	
	/// Creating the dynamic space
	vamde::CellSpace ** s;
	s = new vamde::CellSpace  * [param->n_cell];
	s[0] = new vamde::CellSpace("Tutorials/tutorialMixParticles/input/cell0.input",param);
	s[1] = new vamde::CellSpace("Tutorials/tutorialMixParticles/input/cell1.input",param);
	
	/// Create particles
	initParticles(s);
	
	/// Reset the unique ID of all particles
	for (int n=0; n<param->n_cell; n++) {
		s[n]->refreshCell();
	}
	
	
	/// creating objects
	//~ LenardJones3 lenardjones(s,param->n_cell );
	//~ lenardjones.setRepulsive();
	
	Spring3 spring(s,param->n_cell );
	ZeroForce3 zeroforce(s,param->n_cell);
	LeapFrog3 leapfrog(s,param->n_cell);
	
	/// print Cfgs
	for (int n=0; n<param->n_cell; n++) {
		s[n]->deleteBorders();
		s[n]->cfgwriter->print();
	}
	/// shear.transform();
	
	while (step < end_step) {
		
		/// Print Runtime progress
		tk->print_progress(100);
		/// Advace the clock by one step
		param->clock->advance();
		/// LeapFrog step 1
		leapfrog.step(1);
		for (int n=0; n<param->n_cell; n++) {
			/// delete ghost particles
			s[n]->deleteBorders();
			/// move particles to new cells after integration before moving interprocessor
			s[n]->move_particles_to_the_correct_cells();
			/// move particles across processors.
			s[n]->moveParticles();
			/// copy new ghost particles
			s[n]->copyParticles();
		}
		/// print Cfgs
		for (int n=0; n<param->n_cell; n++) {
			s[n]->cfgwriter->print();
		}
		
		zeroforce.loopAllParticles();
		
		/// Apply Lenard Jones Pair potential
		//~ lenardjones.loopNeighbors();
		spring.loopNeighbors();

		/// LeapFrog step 2
		leapfrog.step(2);
		
		
	}
	/// Print Runtime progress
	tk->print_progress(100);
}


