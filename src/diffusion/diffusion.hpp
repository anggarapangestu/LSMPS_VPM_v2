#ifndef INCLUDED_DIFFUSION
#define INCLUDED_DIFFUSION

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

/**
 *  @brief A particle diffusion subroutine. Perform diffusion to the particle
 *  only for the active region. Also update the vorticity value.
 *
 *  @headerfile diffusion.hpp
 */
class diffusion
{
	// Diffusion calculation
	void diffusion_2d(Particle &_particle);
	void diffusion_3d(Particle &_particle);
public:
	// Diffusion manager calculation
	void main_diffusion(Particle &_particle);
};

#endif
