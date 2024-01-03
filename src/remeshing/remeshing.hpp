#ifndef INCLUDED_REMESHING
#define INCLUDED_REMESHING

#include "../../Utils.hpp"
#include "../neighbor/neighbor.hpp"
#include "../adaptation/adaptation.hpp"

class remeshing
{
	// **Internal variable for particle redistribution
	neighbor neighbor_step;			// Neighbor evaluation and base Cell List (*@param lifetime: throughout simulation)
	adaptation adapt_step;			// The adaptation step (*@param lifetime: one time calculation)
	GridNode *tempGrid;				// Temporaray grid for particle redistribution
	Particle *particleBase;         // Base particle for LSMPS calculation (@param lifetime: throughout simulation)
	bool adap_flag;					 // A flag indicate adaptation is performed

	// std::vector<bool> sign_par;  // sign particles. @param lifetime: singleton
	// std::vector<double> adtParSize;  // Real size of particle after adaptation
	// int initial_num;

	// The private method

	void redistribute_particles(Particle &_particle);
	void redistribute_active_particles(Particle &_particle);
	void set_particle_sign();

public:
	// Remeshing initialization and neighbor search
	void set_neighbor(Particle &_particle, const GridNode &_baseGrid);

	// Remeshing method (*adaptation + redistribution)
	void get_remeshing(Particle &_particle, GridNode &_baseGrid, 
					   const int _iteration, const std::vector<Body> &_bodyList);
};

#endif
