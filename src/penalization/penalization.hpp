#ifndef INCLUDED_PENALIZATION
#define INCLUDED_PENALIZATION

#include "../../Utils.hpp"

/**
 *  @brief Penalization calculation class.
 *  NOTE: Still only designed for 2D simulation.
 *  
 *  @headerfile	penalization.hpp
 */
class penalization
{
	// Internal variables
	// double lambda;				// Porosity parameter @life_time : throughout simulation
    std::vector<int> baseID;    // The list of ID pointing to original particle ID

	// Definition Method
	void lambda_def();			// Define the value of lambda

	// Calculation Function

	void no_slip(Particle &_tempParticle, 
				 Particle &_baseParticle, 
				 const std::vector<Body> &_bodyList, 
				 int _step) const;
	void no_slip_3d(Particle &_tempParticle,
				    Particle &_baseParticle,
				    const std::vector<Body> &_bodyList,
				    int _step) const;
	void no_slip_iterative(Particle &_tempParticle, 
						   Particle &_baseParticle, 
						   const std::vector<Body> &_bodyList, 
						   int _step) const;

public:
	// Main penalization manager

	void get_penalization(Particle &_particle, const std::vector<Body> &_bodyList, int _step);

	// Calculation Method
	
	void calculate_chi(Particle &_particle, const std::vector<Body> &_bodyList) const;  	// Calculate kai parameter
	double get_chi(const double _normalDist) const;	// Get the value of chi based on the normal distance to body surface

	// Penalization constructor
	penalization(){
		// Update lambda value on the construction
    	this->lambda_def();
	}
};

#endif
