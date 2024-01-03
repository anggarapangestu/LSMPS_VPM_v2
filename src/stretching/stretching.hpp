#ifndef INCLUDED_STR
#define INCLUDED_STR

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

/**
 *  @brief A particle stretching subroutine. Perform the stretching to the particle
 *  only for the active region. Also update the vorticity value.
 *
 *  @headerfile stretching.hpp
 */
class stretching
{
public:
    // Stretching calculation
	void main_stretching(Particle &_particle);

    // Combination of diffusion and stretching only for 3D simulation
	void calc_diff_stretch(Particle &_particle);

};

#endif
