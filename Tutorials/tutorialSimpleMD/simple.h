#include "partList.h"

/// input the location of param.input
vamde::Parameters *  init_input ( vamde::Parameters * param){
	
		param->readinput("Tutorials/tutorialSimpleMD/input/param.input" );
		return param;
}

/// Spring Potential using the pair iteration class getCellList
class Spring1 : public getCellList {
	public:
        Spring1(vamde::CellSpace *_s) : getCellList(_s) {}
		
		void action(Particle *i, Particle *j) {
			real sigma = 1.0;
			real epsilon = 5.0;
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

/// Leapfrog Algorithm using the single particle iteration class getPartList
class LeapFrog1: public getPartList {
	
	public:
        LeapFrog1(vamde::CellSpace *_s) : getPartList(_s) {
			PART = 0;
			dt =_s->param-> delta_t;
		}
		void step(int _step){
			PART = _step; 
			apply();
			PART = 0;
		}
	private:
		int PART;
		double dt;
		void action(Particle *i)  {
			for (int d=0; d<DIM; d++) {
				if (PART == 1) {
					i->v[d] = i->v[d] + ( 0.5 * dt * i->F[d] / i->m );
					i->x[d] = i->x[d] + ( dt * i->v[d] );
				} else {
					i->v[d] = i->v[d] + ( 0.5 * dt * i->F[d] / i->m );
				}
			}
		}
};

/// Set All forces to zero, using the single particle iteration class getPartList
class ZeroForce1 : public getPartList {
	
	public:
        ZeroForce1(vamde::CellSpace *_s) : getPartList(_s) {
		}
		void action(Particle *i)  {
			for (int d=0; d<DIM; d++) {
				i->F[d] =0;
				i->na[d] =0;
			}
		}
};



void runSimulation(vamde::Parameters * param0) {
	
	/// Initiating the Parameters object
	vamde::Parameters * param = new vamde::Parameters(MPI_COMM_WORLD);
	param->readinput("Tutorials/tutorialSimpleMD/input/param.input" );
	
	/// Creating the dynamic space
	vamde::CellSpace * s0;
	s0 = new vamde::CellSpace("Tutorials/tutorialSimpleMD/input/cell0.input",param);
	
	/// creating objects
	Spring1 spring(s0);
	ZeroForce1 zeroforce(s0);
	LeapFrog1 leapfrog(s0);
	
	/// print Cfgs
	s0->deleteBorders();
	s0->cfgwriter->print();
	while (param->step() < param->end_step) {
		
		/// Advace the clock by one step
		param->clock->advance();
		
		if ( (param->step() % 100 == 0) && (param->rank() == 0) )std::cout << "Step now: " << param->step() << std::endl;
		
		/// LeapFrog step 1
		leapfrog.step(1);

		/// delete ghost particles
		s0->deleteBorders();
		/// move particles to new cells after integration before moving interprocessor
		s0->move_particles_to_the_correct_cells(); 
		/// move particles across processors.
		s0->moveParticles();
		/// copy new ghost particles
		s0->copyParticles();
		
		/// Set All forces to zero
		zeroforce.apply();
		
		/// Apply Lenard Jones Pair potential
		spring.apply();
		
		/// LeapFrog step 2
		leapfrog.step(2);
		/// print Cfgs
		s0->cfgwriter->print();
		
	}
}


