#include "geometry.hpp"

/**
 *  @brief This code assigns the body translational velocity at each time step.
 * 	Body velocity value is given in the "global.cpp" that are manually inputed.
 *  
 *  @param	_body	The body for translational velocity to be updated.
 * 
 *  @headerfile geometry.hpp
 */
void geometry::define_vel_variation(Body &b)
{
	printf("<+> Define body velocity at each iteration ... \n");
	double _time = omp_get_wtime();

	// Reserve the body velocity translation container
	b.uT.clear();
	b.uT.resize(Pars::max_iter);
	if (DIM > 1){
		b.vT.clear();
		b.vT.resize(Pars::max_iter);
	}
	if (DIM > 2){
		b.wT.clear();
		b.wT.resize(Pars::max_iter);
	}
	
	// Set up the velocity in x-direction
	if (Pars::opt_body_motion == 1){
		// [OPTION 1]: Sudden move and stop at 1.0s [Square curve]
		for (int _it = 0; _it < Pars::max_iter; _it++)
		{
			// Set the real time at this iteration
			double time = _it * Pars::dt;
			// Update velocity
			if (time <= 1.0e0) b.uT[_it] = Pars::ubody;		// Constantly move
			else /*time>=1.0*/ b.uT[_it] = 0.0e0;			// Stop at time >= 1.0 s
		}
	}
	else if (Pars::opt_body_motion == 2){
		// [OPTION 2]: Smoothly start and stop [Square curve with slightly curve at each edge]
		double dT = 0.01e0;		// Time increment
		double t1 = 0.1e0;		// Time set 1
		double t2 = 1.1e0;		// Time set 2
		double t_mid = 0.5 * (t1 + t2);	// Time in the middle between 1 and 2
		double vel;

		for (int _it = 0; _it < Pars::max_iter; _it++)
		{
			// Set the real time at this iteration
			double time = _it * Pars::dt;

			// Calculate the velocity value
			if (time <= t1)
				vel = Pars::ubody * (0.5 * std::exp((time - t1) / dT));
			else if ((time > t1) && (time <= t_mid))
				vel = Pars::ubody * (1.0 - (0.5 * std::exp(-(time - t1) / dT)));
			else if ((time > t_mid) && (time <= t2))
				vel = Pars::ubody * (1.0 - (0.5 * std::exp((time - t2) / dT)));
			else /*(time > t2)*/
				vel = Pars::ubody * (0.5 * std::exp(-(time - t2) / dT));
			
			// Assign the body velocity value
			b.uT[_it] = vel;
		}
	}
	else if (Pars::opt_body_motion == 3){
		// [OPTION 3]: Constant velocity
		for (int _it = 0; _it < Pars::max_iter; _it++) b.uT[_it] = Pars::ubody;
	}
	
	// Set up the velocity in y-direction (constant velocity)
	if (DIM > 1) 
	for (int _it = 0; _it < Pars::max_iter; _it++) b.vT[_it] = Pars::vbody;

	// Set up the velocity in z-direction (constant velocity)
	if (DIM > 2) 
	for (int _it = 0; _it < Pars::max_iter; _it++) b.wT[_it] = Pars::wbody;

	// Simulation time counter
	_time = omp_get_wtime() - _time;
	printf("<-> Set body velocity at each iteration\n");
	printf("    computational time:                [%f s]\n", _time);

}
