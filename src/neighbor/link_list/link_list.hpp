#ifndef INCLUDED_NEIGHBOR_LINK_LIST
#define INCLUDED_NEIGHBOR_LINK_LIST

#include "../../../global.hpp"
#include "../neighbor_utilities/base_grid.hpp"

/**
 *  @brief  Link list neighbor search class. A one-time neighbor search subroutine
 *  featured by link list method. Input of data point coordinates and size properties.
 *  NOTE:
 *      [1] Efficient evaluation limited for domain with less then 2 resolution levels!
 *      [2] The input data structure is container of each properties ordered by ID!
 * 
 *  @headerfile link_list.hpp
*/
class LinkListNgh
{
private:
    /* Short procedure note:
        [1] Divide the domain into grids (one resolution only)
        [2] Group all point in the domain into its corresponding grid
        [3] At each grid collect neighboring grid (for point data neighbor evaluation)
        [4] Evaluate neighbor using link-list interpretation
    */
    /*/ Neighbor Feature
            <+> Point group : 1 res grid
            <+> Evaluation  : neighboring group
            <+> Interaction : link list
    /*/ 

    // Member list
    // ***********
    // Link-list parameter
    std::vector<int> nextID;        // Point the next point ID in the LL chain of current point
    std::vector<int> headID;        // The head of the LL chain at each grid (*hold the point ID data)
    std::vector<int> pntGridID;     // Pair of grid ID where it contains the pnt ID

public:
    void find_neighbor(std::vector<std::vector<int>> &_nghIDList, NghBaseGrid &_baseGrid,
                       const std::vector<double> &_size,
                       const std::vector<double> &_xp,
                       const std::vector<double> &_yp,
                       const std::vector<double> &_zp);
};

#endif