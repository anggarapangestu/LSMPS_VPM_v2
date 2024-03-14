#ifndef INCLUDED_VELOCITY_POISSON
#define INCLUDED_VELOCITY_POISSON

#include "../../Utils.hpp"
#include "../FMM/treeCell.hpp"
#include "../grid_block/gridNode.hpp"
#include "../LSMPS_poisson/poissonLSMPS.hpp"

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
	PoissonLSMPS poissonSolver;		// The poisson solver based on LSMPS [@param lifetime: throughout simulation]

	// Velocity calculation subroutine in 2D
	// =====================================

	void velocity_old(Particle &_particle, const int _step);
	void velocity_fmm_2d(Particle &_particle, const int _step);
	void velocity_LSMPS_poisson_2d(Particle &_particle, const int _step);
	void velocity_LSMPS_poisson_2d_vel(Particle &_particle, const int _step);
	void velocity_LSMPS_poisson_2d_mod(Particle &_particle, const int _step);
	void velocity_LSMPS_poisson_2d_mres(Particle &_particle, const GridNode &_baseGrid, const int _step);
	void velocity_LSMPS_poisson_2d_mres_2(Particle &_particle);

	// Velocity calculation subroutine in 3D
	// =====================================

	void velocity_fmm_3d_fast(Particle &_particle, const int _step);
	void velocity_fmm_3d(Particle &_particle, const int _step);

public:
	// Velocity manager

	void get_velocity(Particle &_particle, const GridNode &_baseGrid, const int _step);
	void get_velocity(Particle &_particle, const int _step);

	// Testing LSMPS Poisson
	void LSMPS_poisson_2d(Particle &_particle, const int _step);
};

#endif
