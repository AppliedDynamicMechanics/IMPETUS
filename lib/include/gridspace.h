// IMPETUS
// version 1.0
// Jul 2016
// Author: Vi Q. Ha


#ifndef VMD_GRIDSPACE_H
#define VMD_GRIDSPACE_H

#include "defines.h"
#include "parameters.h"

namespace vamde {

	enum sendDir { down, up};
	enum syncType { copy, move };

typedef struct {
	int sendproc;
	int receiveproc;
	int sendstart[DIM]; 
	int sendstop[DIM];
	int receivestart[DIM]; 
	int receivestop[DIM];
} Transfer;

typedef struct {
	int lowericstart[DIM], lowericstop[DIM];
	int uppericstartreceive[DIM], uppericstopreceive[DIM];
	
	int uppericstart[DIM], uppericstop[DIM];
	int lowericstartreceive[DIM], lowericstopreceive[DIM];
		
		void print() {
			
			
			for (int d=0; d<DIM; d++) {
				std::cout << 
				" lowericstart[d] " <<
				lowericstart[d] <<
				" lowericstop[d] " <<
				lowericstop[d] <<
				" uppericstartreceive[d] " <<
				uppericstartreceive[d] <<
				" uppericstopreceive[d] " <<
				uppericstopreceive[d] <<
				
				" uppericstart[d] " <<
				uppericstart[d] <<
				" uppericstop[d] " <<
				uppericstop[d] <<
				" lowericstartreceive[d] " <<
				lowericstartreceive[d] <<
				" lowericstopreceive[d] " <<
				lowericstopreceive[d] <<
				std::endl;
				
				
			}
			
			
		}
	
} Communication;



class GridSpace {

public:

	
	Parameters * param;
	/// 3 components physical size of a cell
	real cell_size[DIM]; 	 // dimension of a cell
	/// 3 components starting point of a 'local' cell
	int ic_start[DIM];	 // width of border neighborhood, corresponds to the first local index in the interior of the subdomain
	/// 3 components ending point of a 'local' cell
	int ic_stop[DIM]; 	 // first local index in the upper border neighborhood
	/// 3 components number of local cells (I think) // Edited By Vi on 2015/09/25
	int ic_local[DIM];
	/// I ... don't even know anymore ... the number of cells locally and non-locally?
	int ic_number[DIM]; // number of cells in subdomain, including border neighborhood
	
	/// communicator
	Communication copy_comm[DIM];
	Communication move_comm[DIM];
	void setCopyCommunication();
	void setMoveCommunication();
	
	int rank();
	
	GridSpace(vamde::Parameters *_param);
	void initGrid();
	
	Transfer getTransfer(int dd, bool isCopyNotMove, bool sendUpNotDown) ;
	Transfer getTransferDOWN(int dd, Communication comm) ;
	Transfer getTransferUP(int dd, Communication comm) ;
	
	void copyParticles();
	void moveParticles();
	void sendReceiveCell(bool isCopyNotMove);
	virtual void setup();
	virtual void transfer(int dd,bool isCopyNotMove,bool sendUpNotDown);
	virtual void deleteBorders();
};


} // end namespace
#endif
