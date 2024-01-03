#ifndef INCLUDED_NEIGHBOR_SPATIAL_HASH
#define INCLUDED_NEIGHBOR_SPATIAL_HASH

#include "../../../global.hpp"
#include "../neighbor_utilities/base_grid.hpp"

/**
 *  @brief  Spatial hashing neighbor search class. A one-time neighbor search subroutine
 *  featured by spatial hash method. Input of data point coordinates and size properties.
 *  NOTE:
 *      [1] Efficient evaluation limited for domain with less then 2 resolution levels!
 *      [2] The input data structure is container of each properties ordered by ID!
 * 
 *  @headerfile link_list.hpp
*/
class SpatialHashNgh
{
private:
    /* Short procedure note:
        [1] Divide the domain into grids (one resolution only)
        [2] Group all point in the domain into its corresponding grid
        [3] At each grid collect neighboring grid (for point data neighbor evaluation)
        [4] Evaluate neighbor
    */
    /*/ Neighbor Feature\
         <+> Point group : 1 res grid\
         <+> Evaluation  : neighboring group\
         <+> Interaction : spatial hash
    /*/ 
    
    // Member list
    // ***********
    // Spatial hash parameter
    std::vector<std::vector<int>> pntIDList;    // List of point ID inside each grid
    std::vector<bool> evalFlag;                 // A flag denote current point already done evaluated

public:
    // Spatial Hashing utilities
    //3D: spatial hashing: Ical
    /*void hashing(const int np, const std::vector<double> &xp, const std::vector<double> &yp, 
                const std::vector<double> &zp, std::vector<int> &hx, std::vector<int> &hy,
                std::vector<int> &hz, std::vector<std::vector<int>> &hinside, double cell_size);
    */
    void find_neighbor(std::vector<std::vector<int>> &_nghIDList, NghBaseGrid &_baseGrid,
                       const std::vector<double> &_size,
                       const std::vector<double> &_xp,
                       const std::vector<double> &_yp,
                       const std::vector<double> &_zp);
    // void spatial_hashing(Particle &p, Cell &cell, const int np, const double cell_size);
    
};

#endif