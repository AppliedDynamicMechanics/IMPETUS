// IMPETUS
// version 1.0
// Jul 2016
// Author: Vi Q. Ha


#ifndef VMD_NETWORK_H
#define VMD_NETWORK_H


#include "defines.h"
#include "parameters.h"

namespace vamde {
class Node {
    public:
		Particle * p;

} ;
typedef struct NodeLink {
	Node * n;
	struct NodeLink *next;
} NodeLink;

class Edge {
    public:
		int id;
		Particle  * A;
		Particle  * B;
		double weight;
		class Edge * next;
		int n_periodic[3];
		Edge() {
			id = -1;
			weight = 0;
			n_periodic[0] = 0;
			n_periodic[1] = 0;
			n_periodic[2] = 0;
		}
} ;
	
class Network
{
    private:
	protected:
		int n_nodes;
		bool is_periodic_boundary_condition;
    public:
		Cell nodes ;
		Edge ** edges;
		Parameters * param;
		int computing_rank;
		int reset_n_nodes();
		int get_n_nodes();
		Network(vamde::Parameters * _param) ;
		void init_master();
		virtual void init_nodes();
		virtual ParticleList * createParticleNode(int i, double x, double y, double z);
		void view_nodes();
		virtual void init_edges();
		void view_edges();
		void insertEdge(Edge **root_list, Edge *i);
		void applyPeriodicBoundaryCondition();
};

class MasterNetwork : public Network
{
	
	private:
		Particle * particlearray; 
		Node ** nodepointerarray; /// return this back to private maybe
	public:

	class LinkTable{
		public:
			int  ** links;
			int size;
			
			void view() {
				for (int i=0; i<size; i++) {
					std::cout << i << ":: " ;
					for (int j=0; j<6; j++) {
						std::cout << links[i][j] << " ";
					}
					std::cout << std::endl;
				}
			}
	};
	
	
	struct Cinit{
		double sigma;
		double mass;
		int n_pattern[2];
		double position_shift[3];
		int network_pattern_type ;
		int velocity_distribution_type ;
		double orderly_velocity_distribution [3];
		double uniformly_randomly_velocity_distribution_lower[3];
		double uniformly_randomly_velocity_distribution_upper[3];
		
		char *output_directory_name;
		void readinput(char * inputfilename );

	};
	
	Cinit input;
		Particle * get_p_by_id(int id);
		int rank();
		MasterNetwork(vamde::Parameters * _param);
		MasterNetwork(char * inputfilename ,vamde::Parameters *_param);
		void default_pn_init();
		int get_total_nodes();
		void wrap(Particle * p);
		void wrap_all();
		void scatter_all_to_slaves();
		void gather_local_to_master();
		void build_nodepointerarray();
		void update_nodes_after_gather(Particle * rbuf, int rsize);
		ParticleList * createParticleNode(int i, double x, double y, double z);
		void createAndInsertPairEdge(int u,Particle * a, Particle * b,int x, int y, int z);
		void createAndInsertPairEdge_byID(int u,int a, int b,int x, int y, int z);
		void refreshNetwork();
		void init_square();
		void init_hex(int n_x) ;
		bool is_particle_home(Particle * p, Parameters * param) ;


class AtomEyesCfgSegmentWriter_MasterNetwork{
	
	public:
		double atomic_length;
        AtomEyesCfgSegmentWriter_MasterNetwork(vamde::MasterNetwork *_s);
		void set_atom_name( char * _atom_name);
		void print();
		void printFrag();
		void printSingle();

	private:
		int entry_count;
		int tick ;
		int rank ;
		std::ofstream out;
		double l[3];
		double lo[3];
		int net_print_interval ;
		long double H[3][3];
		vamde::MasterNetwork * s;
		char *output_directory_name;
		char *atom_name;
		void print_part(Particle *i)  ;
		void printInfo();
		void printHeader();
		void printType();
		void printBody ();
		void loop();
};

	AtomEyesCfgSegmentWriter_MasterNetwork * cfgwriter;
	
};


} // end namespace
#endif
