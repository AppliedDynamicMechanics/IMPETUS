// IMPETUS
// version 1.0
// Jul 2016
// Author: Vi Q. Ha


#ifndef VMD_FIELD_H
#define VMD_FIELD_H

#include "contspace.h"

namespace vamde {
	
class InteractiveField : public ContSpace{
public:
	struct nodeOctuplets{
		ContinuumNode * neighbor[2][2][2];
	};
	
	InteractiveField(Cinit _cinit,vamde::Parameters * _param);
	InteractiveField(int _nx, int _ny, int _nz, int _ghost_layer,vamde::Parameters * _param);
	bool is_at_home(double px, double py, double pz);
	bool is_at_home_or_haunted(double px, double py, double pz);

	double linearInterpolate( double in,  double x0, double x1, double y0, double y1) ;
	nodeOctuplets getContNeighbors(double px, double py, double pz);
	
	
	vec3 getGradient(double px, double py, double pz);
	vec3 getGradientNew(double px, double py, double pz);
	
	double getConcentration(double px, double py, double pz);
	void setConcentration(double px, double py, double pz, double secretion_value);
	void leaveConcentration(double px, double py, double pz, double secretion_value);

	void setConcentrationPulse(double px, double py, double pz, double pulse_value);
	void leaveConcentrationPulse(double px, double py, double pz, double pulse_value);
	void setC_v(double px, double py, double pz, double pulse_value);
	void leaveC_v(double px, double py, double pz, double pulse_value);
	void clearGhostC_v();
	void clearGhostConcentrations();
	void sendGhostC_vToNeighborNodes();
};

} // end namespace
#endif
