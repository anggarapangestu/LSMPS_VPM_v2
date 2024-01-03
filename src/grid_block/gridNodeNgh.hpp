#ifndef INCLUDED_NEIGHBOR_GRID_NODE
#define INCLUDED_NEIGHBOR_GRID_NODE

#include "gridNode.hpp"

class GridNodeNgh
{
public:
    void find_neighbor(std::vector<std::vector<int>> &_nghIDList,
                       const GridNode &_baseGrid,
                       const Particle &_evalPar);
};

#endif