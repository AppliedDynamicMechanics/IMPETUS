#include "field.h"
#include "partList.h"
vamde::Parameters *  init_input ( vamde::Parameters * param){
		std::cout << "Project chemotaxis" << std::endl;
		param->readinput("Tutorials/tutorialChemotaxis/input/param.input" );
		return param;
}

class LenardJones3 : public getParticlesAndNeighborsMulti {
	public:
		real range;
		bool isAttractive;
		real amplification;
        LenardJones3(vamde::CellSpace ** _cell_pointer_array, int _n_cellspaces) : getParticlesAndNeighborsMulti( _cell_pointer_array,_n_cellspaces ) {
			default_lj_init();
		}
		void default_lj_init(){
			
			setAttractive();
			setAmplification(1.0);
		}
		void setRepulsive(){
			isAttractive = false;
		}
		
		void setAttractive(){
			isAttractive = true;
		}
		
		void setAmplification(real _amp){
			amplification = _amp;
		}
		
		void NeighborAction(Particle *i, Particle *j) {
			//~ std::cout << i << " rax " << j << std::endl;
			real sigma0 = icell ->input.sigma;
			real sigma1 = jcell ->input.sigma;
			
			real sigma 	= ( sigma0 + sigma1 ) / 2;
			
			
			if (isAttractive) {
				range = sigma* 2.5 ;
			} else {
				range = sigma* 1.122462048309373;
			}
				
			real epsilon = 1.0;
			real r = 0.0;
			for (int d=0; d<DIM; d++){
				r += sqr(j->x[d] - i->x[d]); // squared distance r=rij
			}
			if (r < (range * range))
			{
				real s = sqr(sigma) / r;
				s = sqr(s) * s;
				real f = 24 * epsilon * s / r * (1 - 2 * s);
				f *= amplification;
				for (int d=0; d<DIM; d++){
					i->F[d] += f * (j->x[d] - i->x[d]);
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


class Viscosity3 : public getParticlesAndNeighborsMulti {
	
	public:
		double coefficient; 
        Viscosity3(vamde::CellSpace ** _cell_pointer_array, int _n_cellspaces) : getParticlesAndNeighborsMulti( _cell_pointer_array,_n_cellspaces ) {
			coefficient = 0.1; 
		}
		void ParticleAction(Particle *i)  {
			for (int d=0; d<DIM; d++){
				i->F[d] += - coefficient * i->v[d];
			}
		}
};


/// Migration function built using getPartList on the interactive field 
class Migration_tutorial : public getPartList {
	
	public:
	
		double threshold;
		double f;
        Migration_tutorial(vamde::CellSpace *_s, vamde::InteractiveField * _c) : getPartList(_s) , c(_c){
			threshold = 1 ;
			f =  0.25;
		}
	private:
		vamde::InteractiveField  * c;
		
		void action(Particle *i){
			/// Each particles obtain the concentration value of the field at it's position
			double concentration =  c -> getConcentration(i->x[0],i->x[1],i->x[2]);
			/// Each particles obtain the concentration gradient of the field at it's position 
			vec3 gradient = c -> getGradient(i->x[0],i->x[1],i->x[2]);
			
			double dr[3];
			double rr = 0;
			for (int d=0; d<DIM; d++) {
				dr[d] = gradient.r[d];
				rr += dr[d] * dr[d];
			}
			double x =sqrt(rr); 
			double drhat[3];
			for (int d=0; d<DIM; d++) {
				drhat[d] = dr[d] / x;
			}
			if (x > threshold) {
				for (int d=0; d<DIM; d++){
					/// Particle influence by a force along the unit vector of the concentration gradient at it's position
					i->F[d] += f * drhat[d];
				}
			}
		}
};

void runSimulation(vamde::Parameters * param)  {
	
	vamde::TimeKeeper * tk;
	tk = new vamde::TimeKeeper(param);
	
	/// parameters: 
	double delta_t = param-> delta_t;
	double end_step = param-> end_step;
	double & t = param->clock->t;
	int & step = param->clock->step;
	
	/// Initiate clock
	param->clock->set_step(0);
	
	/// Creating the dynamic space
	vamde::CellSpace ** s;
	std::cout << " param->n_cell " << param->n_cell << std::endl;
	s = new vamde::CellSpace  * [param->n_cell];
	s[0] = new vamde::CellSpace("Tutorials/tutorialChemotaxis/input/cell0.input",param);
	s[0]->cfgwriter->set_atom_name( "Cs");
	
	/// Reset the unique ID of all particles
	for (int n=0; n<param->n_cell; n++) {
		s[n]->refreshCell();
	}
	
	/// Creating the interactive fields 
	vamde::InteractiveField ** c;
	std::cout << " param->n_cont " << param->n_cont << std::endl;
	c = new vamde::InteractiveField  * [param->n_cont];
	
	vamde::InteractiveField::Cinit continit;
	continit.readinput("Tutorials/tutorialChemotaxis/input/cont0.input");
	c[0] = new vamde::InteractiveField(continit,param);
	c[0]->diffusivity = 0.1;
	
	double x_mid = (param->gridinfo->world.lo[0] + param->gridinfo->world.hi[0]) /2;
	double y_mid = (param->gridinfo->world.lo[1] + param->gridinfo->world.hi[1]) /2;
	double z_mid = (param->gridinfo->world.lo[2] + param->gridinfo->world.hi[2]) /2;

	c[0]->setConcentration(x_mid+0.4*x_mid,y_mid,z_mid, 10 );
	c[0]->setConcentration(x_mid-0.4*x_mid,y_mid,z_mid,-10 );
	c[0]->setAllGlobalBoundarySink();
	
	for (int n=0; n<param->n_cont; n++) {
		c[n]->cfgwriter->print();
	}
	
	/// creating objects
	LenardJones3 lenardjones(s,param->n_cell );
	ZeroForce3 zeroforce(s,param->n_cell);
	LeapFrog3 leapfrog(s,param->n_cell);
	Viscosity3 viscosity(s,param->n_cell);
	
	for (int n=0; n<param->n_cell; n++) {
		s[n]->deleteBorders();
		s[n]->cfgwriter->print();
	}
	
	while (step < end_step) {
		/// Print Runtime progress
		tk->print_progress(10);
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
		/// Set All forces to zero
		zeroforce.loopLocalParticles();

		/// Apply Lenard Jones Pair potential
		lenardjones.loopNeighbors();
		
		viscosity.loopLocalParticles();
		/// LeapFrog step 2
		leapfrog.step(2);

		/// print Cfgs
		for (int n=0; n<param->n_cell; n++) {
			s[n]->cfgwriter->print();
		}
		
		for (int n=0; n<param->n_cont; n++) {
			c[n]->copyParticles();
		}
		
		for (int n=0; n<param->n_cont; n++) {
			c[n]->Iterate();
		}
		
		/// continuum
		c[0]->cfgwriter->print();
		c[0]->setConcentration(x_mid+0.4*x_mid,y_mid,z_mid,  10 );
		c[0]->setConcentration(x_mid-0.4*x_mid,y_mid,z_mid, -10 );
		
		Migration_tutorial mg0(s[0], c[0]);
		mg0.apply();
	}
	
}

