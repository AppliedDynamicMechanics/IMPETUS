#include "partList.h"

/// input the location of param.input
vamde::Parameters *  init_input ( vamde::Parameters * param){
	
		param->readinput("Tutorials/tutorialSimpleMD/param.input" );
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

/// Lennard-Jones Potential using the pair iteration class getCellList
class LennardJones1 : public getCellList {
	public:
		real range;
		bool isAttractive;
        LennardJones1(vamde::CellSpace *_s) : getCellList(_s) {
			default_lj_init();
		}
		void default_lj_init(){
			
			setAttractive();
		}
		void setRepulsive(){
			isAttractive = false;
		}
		
		void setAttractive(){
			isAttractive = true;
		}
		
		
		
		void action(Particle *i, Particle *j) {
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
				r += sqr(j->x[d] - i->x[d]);
			}
			if (r < (range * range))
			{
				
				real s = sqr(sigma) / r;
				s = sqr(s) * s;
				real f = 24 * epsilon * s / r * (1 - 2 * s);
				for (int d=0; d<DIM; d++){
					i->F[d] += f * (j->x[d] - i->x[d]);
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

/// Mean Square Displacement recording using the single particle iteration class getPartList
class MSD1 : public getPartList {

	public:
		double MSD_sum;
		int computingRank;
		int step_delay;
		
        MSD1(vamde::CellSpace *_s) : getPartList(_s)  {
			MSD_sum = 0;
			computingRank = 0;
			std::ofstream out;
			out.open("Tutorials/tutorialSimpleMD/output/MSD/files/MSD.out" , std::ofstream::out | std::ofstream::trunc);
			out.close();
			step_delay = 0;
		}
		void set_step_delay( int _step_delay){
			step_delay = _step_delay;
		}
		
		
		
			
		void print() {
			double step = s->param->step();
			if (step >= step_delay) {
				clear_sum();
				apply();
				comp_sum();
			}
		}
	private:
		void clear_sum(){
			MSD_sum = 0;
		}
		void comp_sum() {
			
			double MSD_GLOBAL_SUM = sum_mpi_double(s->param->gridinfo->mpicomm, MSD_sum,computingRank);
	
			if ( s->param->rank() == computingRank ) {
				char *buf = new char[100];
				sprintf (buf, "Tutorials/tutorialSimpleMD/output/MSD/files/MSD.out" );
				std::ofstream out;
				double timeNow = s->param->clock->getTime();
				out.open(buf,std::ios::app);
				
				double MSD_avg = MSD_GLOBAL_SUM/(s->n_total_global_particles);
				out  << timeNow << " " <<  MSD_avg << std::endl;
				out.close();
			}
		}
		void action(Particle *i)  {
			double dr_dot_dr = 0;
			double dr[DIM];
			for (int d=0; d<DIM; d++) {
				dr[d] = i->x[d] - i->x_original[d] + (s->param->gridinfo->world.l[d] * i->boundary_cross_count[d]);
				dr_dot_dr += dr[d] * dr[d];
			}
			MSD_sum += dr_dot_dr;
		}
};

/// Radial Distribution Function recording using the single particle iteration class getPartList
class RDF1 : public getCellList {
	
	public:

		double rangeRdf;
		/// the number of indexes of the array histRdf
		int sizeHistRdf;
		/// the number of steps the program collect before printing
		int limitRdf ;
		/// choose the rank you wish to do the calculations
		int computing_rank;

		RDF1(vamde::CellSpace *_s) : getCellList(_s) {
			rangeRdf = 4.0;
			/// the number of indexes of the array histRdf
			sizeHistRdf = 200;
			/// the number of steps the program collect before printing
			limitRdf = 100; // default
			countRdf = 0;
			deltaR= 0;
			/// choose the rank you wish to do the calculations
			computing_rank = 0;
			
			histRdf = new double[sizeHistRdf];
			histRdf_global = new double[sizeHistRdf];
			deltaR = rangeRdf / sizeHistRdf;
			
			writemany = false;
			output_directory_name = "output/RDF/";
		}
		void init(){
			delete histRdf;
			delete histRdf_global;
			histRdf = new double[sizeHistRdf];
			histRdf_global = new double[sizeHistRdf];
			deltaR = rangeRdf / sizeHistRdf;
		}
		void print(){
			//~ std::cout << "pring rdf \n" ;
			clearArray();
			apply();
			compute();
		}
		void set_output_directory(char * _output_directory_name)
		{
			output_directory_name = _output_directory_name;
		}
		void set_writemany(){
			writemany = true;
		}
		void set_limitRdf(int num){
			limitRdf = num;
		}
			
		void set_writeone(){
			writemany = false;
		}
	private:
		char *output_directory_name;
		bool writemany;
		int countRdf;
		real deltaR;
		/// array to store distances between particles
		double *histRdf;
		double *histRdf_global;
		
		void clearArray() {
		/// parallel version
			if (countRdf == 0) {
				for (int n = 0; n < sizeHistRdf; n ++) {
					histRdf[n] = 0.;
					histRdf_global[n] = 0.;
				}
			}
		}

		void action(Particle *i, Particle *j) {

			  double dr[DIM];
			  for (int d=0; d<DIM; d++) 
				dr[d] = i->x[d] - j->x[d];
			  
			  double rr = 0;
			  for (int d=0; d<DIM; d++) 
				rr += dr[d] * dr[d];
			  int n;
			  if (rr < (rangeRdf*rangeRdf)) {
				n = sqrt (rr) / deltaR;
				++ histRdf[n];
			  }
		}
		
		void compute () {
			real normFac;
			++ countRdf;
			if (countRdf == limitRdf) {
				get_global_sum_histRdf();
				double boxVolume = 1;
				for (int d=0; d<DIM; d++) 
					boxVolume *= cs->param->gridinfo->world.l[d];
				
				int nMol = cs->n_total_global_particles;
				normFac = boxVolume/ (2. * M_PI * (deltaR) *(deltaR) *(deltaR) *
				nMol*nMol * countRdf);
				for (int n = 0; n < sizeHistRdf; n ++)
					histRdf_global[n] *= normFac / (n - 0.5) / (n - 0.5) /2; 
				if (cs->param->rank() == computing_rank) write_data();
				countRdf = 0;
			}
		}
		
		void get_global_sum_histRdf() {
			MPI_Reduce(histRdf,histRdf_global,sizeHistRdf,MPI_DOUBLE,MPI_SUM,computing_rank,
              cs->param->gridinfo->mpicomm);
		}
		void write_data (){
			
			char *buf = new char[100];
			sprintf (buf, "%sgr_%d.out", output_directory_name, cs->param->step());
			std::ofstream out;
			if (writemany) out.open(buf);

			sprintf (buf, "%sgr_current.out", output_directory_name);
			std::ofstream out_new;
			out_new.open(buf);
			
			long double rb;
			int n;
			for (n = 0; n < sizeHistRdf; n ++) {
				rb = (n + 0.5) * rangeRdf / sizeHistRdf;
				if (writemany) out << rb << " " << histRdf_global[n] << std::endl;
				out_new << rb << " " << histRdf_global[n] << std::endl;
			}
			if (writemany) out.close();
			out_new.close();
			delete buf;
		}
		
};


void runSimulation(vamde::Parameters * param0) {
	
	/// Initiating the Parameters object
	vamde::Parameters * param = new vamde::Parameters(MPI_COMM_WORLD);
	param->readinput("Tutorials/tutorialSimpleMD/input/param.input" );
	
	
	vamde::TimeKeeper * tk;
	tk = new vamde::TimeKeeper(param);
	std::cout << "New Vamde " << std::endl;
	/// Creating the dynamic space
	vamde::CellSpace * s0;
	s0 = new vamde::CellSpace("Tutorials/tutorialSimpleMD/input/liquid.input",param);
	
	/// creating objects
	LennardJones1 lennardjones(s0);
	lennardjones.setAttractive();
	ZeroForce1 zeroforce(s0);
	LeapFrog1 leapfrog(s0);
	
	MSD1 msd(s0);
	RDF1 rdf(s0);
	rdf.set_output_directory("Tutorials/tutorialSimpleMD/output/RDF/files/");
	rdf.init();
	
	/// print Cfgs
	s0->deleteBorders();
	s0->cfgwriter->print();
	
	while (param->step() < param->end_step) {
		
		/// Print Runtime progress
		tk->print_progress(100);
		/// Advace the clock by one step
		param->clock->advance();
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
		lennardjones.apply();
		
		/// LeapFrog step 2
		leapfrog.step(2);
		
		/// print msd
		msd.print();
		
		/// print rdf
		rdf.print();
		/// print Cfgs
		s0->cfgwriter->print();
		
	}
}


