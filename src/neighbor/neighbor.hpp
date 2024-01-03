#ifndef INCLUDED_NEIGHBOR
#define INCLUDED_NEIGHBOR

#include "../../Utils.hpp"
#include "neighbor_utilities/base_grid.hpp"		// Utilities (base grid, ...)
#include "direct/direct_find.hpp"				// Direct find method
#include "link_list/link_list.hpp"				// Link list method
#include "spatial_hash/spatial_hash.hpp"		// Spatial hash method
#include "cell_list/cell_list.hpp"				// Cell list method
#include "inter_search/inter_search.hpp"		// Inter search method
#include "../grid_block/gridNodeNgh.hpp"		// Grid node method

/**
 *  @brief  Neighbor evaluation class consisted of method to generate neighbor data 
 *  ID list.
 * 
 *  @headerfile neighbor.hpp
*/
class neighbor
{
private:
	// The instances inside neighbor
	CellList cellListData;		// A cell list lifetime throughout the simulation

	// Utilities function
	void create_temp_grid(const Particle &_parEval, NghBaseGrid &_tempGrid);

public:
	// Main neighbor search manager
	void neigbor_search(Particle &_evalPar, const GridNode &_baseGrid); // Search neighbor toward itself
	void neigbor_search(Particle &_evalPar, Particle &_srcPar);         // Search neighbor toward source particle data
	
	// Obsolete function [NEED FURTHER Change]
	void cell_list_init(Particle& parEval);
	bool particle_adaptation(const Particle& parEval, Particle& baseParticle, std::vector<double>& PARsize);

	/* Integrated with parallel computing [?]
		> Direct search	: [X]
		> Link List		: [X]
		> Cell List		: [X]
		> Spatial hash	: Yes
		> Grid Node		: Yes
		> Inter search	: Yes
	*/

};

#endif

/** Two types of neighbor search
 *   [1] One-time use search lifetime	-> All variable only used for one time search
 * 		-> Link_list, spatial_hashing, inter_search
 * 		# No need particle storage, only use a spatial data
 * 
 *   [2] Throughout simulation lifetime -> The member data is saved throughout the simulation
 * 		-> Cell list, Grid node
 *		# Need particle data storage for data integration
*/
