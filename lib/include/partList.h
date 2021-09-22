// IMPETUS
// version 1.0
// Jul 2016
// Author: Vi Q. Ha




#ifndef VMD_LIST_H
#define VMD_LIST_H


#include "defines.h"
#include "cell.h"
#include "field.h"
#include "network.h"

class getPartList
{
    private:
		
	protected:
		
		vamde::CellSpace *s;
		Cell *grid;
		vamde::CellSpace *icell;
	
    public:
		void apply();
		getPartList(vamde::CellSpace *_s) ;
		virtual void action(Particle *i) {
			std::cout << "You need to redefined action() for the child class" << std::endl;
		}

};

class getLocalList /// NOTE: THIS MIGHT NOT WORK ANYMORE, NEED TO CHECK IC_START AND IC_STOP
{
    private:
		
	protected:
		
		vamde::CellSpace *s;
		Cell *grid;
		vamde::CellSpace *icell;
	
    public:
		void apply();
		getLocalList(vamde::CellSpace *_s) ;
		virtual void action(Particle *i) {
			std::cout << "You need to redefined action() for the child class" << std::endl;
		}

};

class getCellList
{
    private:
		
	protected:
		
		vamde::CellSpace *cs;
		Cell *grid;
		vamde::CellSpace *icell;
		vamde::CellSpace *jcell;
	
    public:
		void apply();
		void loopCells_applyAction();
		
		getCellList(vamde::CellSpace *_cs);
		virtual void action(Particle *i, Particle *j) {
			std::cout << "You need to redefined action() for the child class" << std::endl;
		}

};

class getParticlesAndNeighborsMulti {
	
	
	int pairTotal ;
	int pairCount ;
	protected:
	
	vamde::Parameters * param;
	int n_cellspaces;
	protected:
	vamde::CellSpace ** cell_pointer_array;
	double * cell_minsize_array;
	int * cell_order_list;
	double isigma;
	double jsigma;
	vamde::CellSpace *icell;
	vamde::CellSpace *jcell;
	void ExchangeSort();
	
	void loopTwoCells(vamde::CellSpace *s0, vamde::CellSpace *s1);
	
	public: 
	
	class ParticlePair {
		public:
			Particle  * i;
			Particle  * j;
	} ;
	ParticlePair * pairlist;
	
	getParticlesAndNeighborsMulti(vamde::CellSpace ** _cell_pointer_array, int _n_cellspaces);
	getParticlesAndNeighborsMulti(vamde::CellSpace * _s0, vamde::CellSpace *_s1);
	void default_init();
	
	void countPairsTwoCells(vamde::CellSpace *s0, vamde::CellSpace *s1) ;
	int countPairs();
	void buildPairTwoCells(vamde::CellSpace *s0, vamde::CellSpace *s1) ;
	void destroyPairList() ;
	void buildPairList() ;
	void loopPairList();
	void loopAllParticles() ;
	void loopLocalParticles() ;
	
	
	
	void loopIntraCell() ;
	void loopInterCell() ;
	void loopNeighbors();
	virtual void NeighborAction(Particle *i, Particle *j) {
		std::cout << "You need to redefined action() for the child class" << std::endl;
	}
	virtual void ParticleAction(Particle *i) {
		std::cout << "You need to redefined action() for the child class" << std::endl;
	}
	
};


class getGloballyUniquePairList
{
    private:
		
	protected:
		vamde::CellSpace *s;
		Cell *grid;
	
    public:
		void apply();
		void loopCells_applyAction();
		
		getGloballyUniquePairList(vamde::CellSpace *_s);
		virtual void action(Particle *i, Particle *j) {
			std::cout << "You need to redefined action() for the child class" << std::endl;
		}

};

#include "network.h"

class getNodesAndEdges
{
    private:
		
	protected:
		vamde::MasterNetwork * n;
		vamde::Parameters * param;
		vamde::Edge * e;
    public:
		getNodesAndEdges(vamde::MasterNetwork * _n):n(_n){
			param = n->param;
		}
		void loopNodes();
		
		void loopEdges();
		virtual void nodeAction(Particle *i) {
			std::cout << "You need to redefined action() for the child class" << std::endl;
		}
		virtual void edgeAction(Particle *i, Particle *j) {
			std::cout << "You need to redefined action() for the child class" << std::endl;
		}

};
class getNodesAndEdgesMulti
{
    private:
		
	protected:
		vamde::MasterNetwork ** net_pointer_array;
		vamde::Parameters * param;
		vamde::Edge * e;
		int n_networks;
    public:
		getNodesAndEdgesMulti(vamde::MasterNetwork ** _net_pointer_array, int _n_networks):net_pointer_array(_net_pointer_array),n_networks(_n_networks){
			param = net_pointer_array[0]->param;
		}
		void loopNodes();
		void loopEdges();
		virtual void NodeAction(Particle *i) {
			std::cout << "You need to redefined action() for the child class" << std::endl;
		}
		virtual void EdgeAction(Particle *i, Particle *j) {
			std::cout << "You need to redefined action() for the child class" << std::endl;
		}


};


class NetworkCellInteraction : public getParticlesAndNeighborsMulti,  public getNodesAndEdgesMulti  {
	private:
public:
		vamde::Parameters * param;
	NetworkCellInteraction (vamde::CellSpace ** _cell_pointer_array, int _n_cellspaces,vamde::MasterNetwork ** _net_pointer_array, int _n_networks):
														getParticlesAndNeighborsMulti( _cell_pointer_array, _n_cellspaces),
														getNodesAndEdgesMulti(_net_pointer_array, _n_networks)
	{
		param = getNodesAndEdgesMulti::param;
	}
	void loopLocalSingles(){
		loopLocalParticles();
		loopNodes();
	}
	void loopSingles(){
		loopLocalSingles();
	}
	void loopAllPairs(){
		loopNeighbors();
		loopParticlesWithNodes();
	}
	virtual void SingleAction(Particle *i){
		std::cout << "You need to redefined action() for the child class" << std::endl;
	}
	void ParticleAction(Particle *i){
		SingleAction(i);
	}
	void NodeAction(Particle *i){
		SingleAction(i);
	}
	
	virtual void PairAction(Particle *i, Particle *j){
		std::cout << "You need to redefined action() for the child class" << std::endl;
	}
	
	void NeighborAction(Particle *i, Particle *j){
		PairAction(i,j);
	}
	void EdgeAction(Particle *i, Particle *j){
		PairAction(i,j);
	}
	
	void loopParticlesWithNodes();
	
	void loopNetCell_Interactions( vamde::CellSpace * s , vamde::MasterNetwork * n );
	
	void loopNetCell_Interactions( vamde::MasterNetwork * n , vamde::CellSpace * s );
	void singleNetCell_Interactions_Bruteforce(  Particle * p , vamde::MasterNetwork * n , vamde::CellSpace * s );
	void singleNetCell_Interactions_ByCells(  Particle * p , vamde::MasterNetwork * n , vamde::CellSpace * s );
	
	
};

#endif
