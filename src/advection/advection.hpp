#ifndef INCLUDED_ADVECTION
#define INCLUDED_ADVECTION

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

/**
 *  @brief A particle advection subroutine. Perform advection to the particle
 *  using the given type. There are euler and runge kutta 2 type.
 *
 *  @headerfile advection.hpp
 */
class advection
{
	// The euler advection
	void advection_euler(Particle &_particle);
	
	// Trying using runge kutta 2nd order
	void advection_rk2(Particle &_particle);
public:
	// The public advection
	void main_advection(Particle &_particle);
};

#endif
