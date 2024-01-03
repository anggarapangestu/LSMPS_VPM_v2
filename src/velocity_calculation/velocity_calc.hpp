#ifndef INCLUDED_VELOCITY_POISSON
#define INCLUDED_VELOCITY_POISSON

#include "../../Utils.hpp"
#include "../FMM/treeCell.hpp"

/**
 *  @brief Velocity calculation class. Manage all type of velocity calculation.
 *  NOTE: Still limited with FMM calculation only. Direct biot savart is not included.
 *  
 *  @headerfile	velocity_calc.hpp
 */
class VelocityCalc
{
private:
	// Internal Variable
	TreeCell treeData;		// The basis of tree data of poisson solver [@param lifetime: throughout simulation]

	// Velocity calculation

	void velocity_old(Particle &_particle, const int _step);
	void velocity_fmm_2d(Particle &_particle, const int _step);
	void velocity_fmm_3d_fast(Particle &_particle, const int _step);
	void velocity_fmm_3d(Particle &_particle, const int _step);

public:
	// Velocity manager

	void get_velocity(Particle &_particle, const int _step);
};

#endif
