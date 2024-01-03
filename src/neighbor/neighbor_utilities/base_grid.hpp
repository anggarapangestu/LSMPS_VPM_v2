#ifndef INCLUDED_NEIGHBOR_BASE
#define INCLUDED_NEIGHBOR_BASE

#include "../../../global.hpp"

/**
 *  @brief  A grid data management for one-time neighbor search. Featured with single
 *  resolution evenly distributed grid size that support for 2D and 3D simulation.
 * 
 *  @headerfile base_grid.hpp
*/
class NghBaseGrid
{
private:
    /* GRID ILLUSTRATION
        Grid Illustration for 2D model (*the expansion to 3D also work)
         (One resolution grid)
         _____________________________             
        |  x  |x 16 | {17}| {18}|  x  | 3          Description:
        |_x___|_____|_____|x____|_____|              [x] : example of point data
        |  x  | {11}|  x  | {13}| x   | 2  [index]   {#} : Grid ID
        |_____|_____|_____|_____|_____|    [in y ] 
        | {5}x|   x | {7} |   x | {9} | 1  [basis] Terminology:
        |_____|_____|_____|_____|_____|              > index : position coordinate
        | {0} | {1} |  x x| {3} |x{4} | 0                       of the grid in each basis
        |_x___|_____|_x___|_____|_____|              > ID    : label of grid in a spesific
           0     1     2     3     4                            order (see illust.)
               [index in x basis]

        Example:
        * Grid {2} located at index [2,0] contains 3 point
        * Grid {11} located at index [1,2] contains 0 point
    */
    
    // Member list
    // ***********
    // Grid data
    int indexCount[DIM];    // Number of index at each basis
    double pivCoor[DIM];    // The pivot coordinate (reference coordinate -> global minimum location)
    int gridCount;          // Total number of grid
    double gridSize;        // The size of grid edge

    // Grid index transformation
    int flattenMul[DIM];    // Multiplier for index flatten [1, cx, cx*cy]; c: index count

public:
    // Method list
    // ***********
    // Access internal data
    double get_size();      // Access grid size
    int get_count();        // Access total grid count

    // Grid initialization transformation
    void initialize_grid(double maxCoor[DIM], double minCoor[DIM],
                         const std::vector<double> &_parSize);

    // Grid index transformation
    int get_grid_ID(double _pntCoor[DIM]) const;
    void get_grid_index(int _index[DIM], int ID) const;

    // Neighbor evaluation utilities
    void find_grid_ngh(std::vector<int> &_nghIDList, int _gridID) const;
};

#endif
