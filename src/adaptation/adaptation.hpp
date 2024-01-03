#ifndef PARTICLE_ADAPTATION
#define PARTICLE_ADAPTATION

#include "../grid_block/gridNodeAdapt.hpp"

class adaptation
{
private:
    // * Internal variables
public:
    bool get_adaptation(const Particle &_parEval, Particle *&_parBase, GridNode &_baseGrid);
};

#endif
