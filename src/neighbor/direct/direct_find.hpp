#ifndef INCLUDED_NEIGHBOR_DIRECT_SEARCH
#define INCLUDED_NEIGHBOR_DIRECT_SEARCH

#include "../../../global.hpp"

class directFindNgh
{
public:
    void find_neighbor(std::vector<std::vector<int>> &_nghIDList,
                       const std::vector<double> &_size,
                       const std::vector<double> &_xp,
                       const std::vector<double> &_yp,
                       const std::vector<double> &_zp);
};

#endif