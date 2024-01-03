#include "velocity_biot_savart.hpp"
#include <omp.h>
#define MIN_DISTANCE_FACTOR 1e-5  // The nearest source point still in consideration

void velocity_biot_savart::biotsavart_direct_2d(Particle &tarPar, Particle &srcPar)
{
	// Calculate biot-savart for each particle
	#pragma omp parallel for
	for (int i = 0; i < tarPar.num; i++){
		// Internal variables
		double dx, dy;
		double rij2, rij, sij, sij2;
		double gsig;

		// Calculate the biot-savart sum toward all particle in the domain
		for (int j = 0; j < srcPar.num; j++){
			// Position difference between target and source
			dx = tarPar.x[i] - srcPar.x[j];
			dy = tarPar.y[i] - srcPar.y[j];
			
			// Calculating the target distance toward source
			rij2 = dx*dx + dy*dy;
			rij = std::sqrt(rij2);

			if (rij > (MIN_DISTANCE_FACTOR*Pars::sigma)){
				// Calculating the regulation function
				sij2 = (tarPar.s[i] * tarPar.s[i]) + (srcPar.s[j] * srcPar.s[j]);
				sij = 0.5 * std::sqrt(sij2);	// Average size by averaging the area
				
				// Calculating the regulation function (Effect of attempting continous to discrete)
				regul_func_2d(rij, sij, gsig);

				// Sum the biot-savart velocity
				tarPar.u[i] -= (dy * srcPar.gz[j]) * gsig / rij2;
				tarPar.v[i] += (dx * srcPar.gz[j]) * gsig / rij2;
			}
		}
	}
}
