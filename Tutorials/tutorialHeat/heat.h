#include "field.h"

/// input the location of param.input
vamde::Parameters *  init_input ( vamde::Parameters * param){
		std::cout << "Project heat tutorial" << std::endl;
		param->readinput("Tutorials/tutorialHeat/input/param.input" );
		return param;
}

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
	
	/// Creating the interactive fields 
	vamde::InteractiveField * c0;
	vamde::InteractiveField::Cinit continit;
	continit.readinput("Tutorials/tutorialHeat/input/cont0.input");
	c0 = new vamde::InteractiveField(continit,param);
	c0->diffusivity = 0.05;
	
	double x_mid = (param->gridinfo->world.lo[0] + param->gridinfo->world.hi[0]) /2;
	double y_mid = (param->gridinfo->world.lo[1] + param->gridinfo->world.hi[1]) /2;
	double z_mid = (param->gridinfo->world.lo[2] + param->gridinfo->world.hi[2]) /2;
	
	/// set concentration to field at value of 10
	c0->setConcentration(x_mid,y_mid,z_mid,10 );
	/// set global boundary condition to sink (constant value of zero)
	c0->setAllGlobalBoundarySink();
	/// print visualization files 
	c0->cfgwriter->print();
	
	while (step < end_step) {
		
		/// Print Runtime progress
		tk->print_progress(10);
		/// Advace the clock by one step
		param->clock->advance();
		/// LeapFrog step 1
		c0->copyParticles();
		/// Apply parabolic equation to the whole continuum
		c0->IterateParabolicPDE();
		/// print visualization files 
		c0->cfgwriter->print();

		
	}
}

