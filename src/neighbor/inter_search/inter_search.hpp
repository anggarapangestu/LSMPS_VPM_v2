#ifndef INCLUDED_NEIGHBOR_INTER_SEARCH
#define INCLUDED_NEIGHBOR_INTER_SEARCH

#include "../../../Utils.hpp"
#include "../neighbor_utilities/base_grid.hpp"

/**
 *  @brief  A one-time neighbor search between two set data (source and target 
 *  point data). Featured with single resolution evenly distributed grid size 
 *  that support for 2D and 3D simulation. Follow the spatial hashing method.
 * 
 *  @headerfile inter_search.hpp
*/
class InterSearchNgh
{
private:
    /** Illustration of inter search
     *  # For example the target is distributed evently
     *     while the source is scattered
     *      
     *          Domain of the                
     *         source particle                   Domain of the 
     *        ___________________               target particle
     *       |      *          * |            ___________________
     *       |   *     *   *     |           | x    x    x    x  |
     *       |  *                |           |                   |
     *       |        *      *   |           | x    x    x    x  |
     *       |      *   *        |           |                   |
     *       |  *      *     *   |           | x    x    x    x  |
     *       |___________________|           |___________________|
     *      
     *                                        Grid grouped Domain
     *          Overlayed Domain            _______________________  
     *        ___________________          |_____|_____|_____|_    |
     *       |      *          * |         |_____|*____|_____*_|___|
     *       | x *  x  * x *  x  |         | x * |x  * x *  x| |   |
     *       |  *                |         |__*__|_____|_____|_|___|
     *       | x    x *  x   *x  |         | x   |x *  x   *x| |   |
     *       |      *   *        |         |_____|*___*|_____|_|___|
     *       | x*   x  * x   *x  |         | x*  |x  * x   *x| |   |
     *       |___________________|         |_____|_____|_____|_|___|
     *       
     *  # In this subroutine, the neighbor ID list is the list of 
     *     source particle ID that is a neighbor to the target particle.
     * 
     *  # For each particle [x], check the particle [*] that lies under
     *     the support domain of particle [x].
     * 
     *  # NOTE: This method doesn't seek any neighbor correlation between 
     *     particle [x]_i and [x]_j.
    */

    // Member list
    // ***********
    // Spatial hash parameter
    std::vector<std::vector<int>> pntIDListSrc; // List of source point ID inside each grid

public:
    // Neighbor search between 2 point data set
    
	void find_neighbor(const Particle &_sourcePar, const Particle &_targetPar, 
                       std::vector<std::vector<int>> &_nghList);
};

#endif