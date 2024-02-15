#include "remeshing.hpp"
#include "../geometry/geometry.hpp"

/**
 *  @brief  Particle redistribution manager. Evaluate the particle adaptation and 
 *  particle redistribution.
 *         
 *  @param  _particle  Particle data container.
 *  @param  _baseGrid  Grid node data for adaptation tools.
 *  @param  _iteration The current iteration.
 *  @param	_bodyList  Body data container.
*/
void remeshing::get_remeshing(Particle &_par, GridNode &_gridNode, const int _step, 
							  const std::vector<Body> &_bodyList)
{
	/*
	Inside Method:
	> Perform Distribution Adaptation
	  * Procedure (Use the one from cell list or grid Node)
	    - <Check whether the Adaptation is needed in the current state>
		- Evaluate the adapted particle (create a list)
		- Perform the split (available) and merging (still not available)
		- Update the cell List
		- Update the neighbor list
		- Update the particle base
	  * Needs:
	    - Determine the update region (Evaluate each particle vorticity value)
	
	*/

	/** Illustration of data update between adaptation and redistribution
     *      PARTICLE DATA     |      ADAPTATION        |      REDISTRIBUTION
     *  ----------------------|------------------------|-------------------------
     *   * Coordinate         | - Adapted              | - Interpolated to base data
     *                        | - Store in the base    |       
     *   * Size               | - Set to old base      | - Update to new base
     *   * Num                | - New data             | - New data
     *   * Neigbor list       | - Connect to old base  | - Evaluate new
     *   * Velocity           | - Resize with new num  | - Resize with 0
     *   * Vorticity          | - Resize with new num  | - Interpolated to new position 
     *   * Vortex strength    | - Resize with new num  | - Recalculated from vorticity
     *   * Pressure           | [-]                    | [-]
     *   * Particle flag      | - All flags is         | [-]
	 *      > chi             |   re-evaluated from    |                               
	 *      > isActive        |   geometry             | - Re-evaluated from vorticity
	 *      > All body flag   |                        |                               
    */
	
	// Computational time manager
	double _time;
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick;
        std::chrono::duration<double> span;
    #endif

	// this->initial_num = this->particleBase->num;
	// this->adtParSize = this->particleBase->s;

	// ============================================================= //
	// ========== Particle Distribution Adaptation Update ========== //
	// ============================================================= //
    // Initial set for the adaptation flag
	this->adap_flag = false;

	if (Pars::flag_adaptive_dist && (_step%Pars::adapt_inv == 0)/* && (_step != 0)*/){
		// Evaluating and Performing particle adaptation
		printf("\nPerform Particle Adaptation ...\n");
		#if (TIMER_PAR == 0)
			// Timer using super clock (chrono)
			tick = std::chrono::system_clock::now();
		#elif (TIMER_PAR == 1)
			// Timer using paralel package
			_time = omp_get_wtime();
		#endif


		// // OLD METHOD [Using cell list]
		// // The current method will evaluate and peforming particle adaptation
		// this->adap_flag = this->neighbor_step.particle_adaptation(_par, *this->particleBase, this->adtParSize);

		// Perform the particle adaptation
		this->adap_flag = this->adapt_step.get_adaptation(_par, this->particleBase, _gridNode);

		// Update the particle flag data
		if (this->adap_flag == false){
			MESSAGE_LOG << "Adaptation checked! No adaptation!\n";
		}
		else{
			printf("Evaluate the near body part ...\n");
			geometry geom_tool;
			geom_tool.eval_near_body(*this->particleBase, _bodyList);    // Calculate the near body object (for penalization)
			geom_tool.eval_body_flag(*this->particleBase, _bodyList);    // Calculate isSurface (for adaptation) and also chi (for penalization) and active flag
		}
		
		// Particle adaptation summary time display
		#if (TIMER_PAR == 0)
			// Timer using super clock (chrono)
			span = std::chrono::system_clock::now() - tick;
			_time = span.count();
		#elif (TIMER_PAR == 1)
			// Timer using paralel package
			_time = omp_get_wtime() - _time;
		#endif
		printf("<-> Particle adaptation comp. time:    [%f s]\n", _time);
	}
	
	// ============================================= //
	// ========== Particle Redistribution ========== //
	// ============================================= //
	// Update the current grid
	this->tempGrid = &_gridNode;

	// Evaluating and Performing particle redistribution
	printf("\nPerform Particle Redistribution ...\n");
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime();
    #endif

	// Perform the particle redistribution
	if (Pars::flag_fast_remesh){
		// [Type 2] Redistribute near active region only [Turns out not very effective]
		this->redistribute_active_particles(_par);
	}
	else{
		// [Type 1] Redistribute all particle
		this->redistribute_particles(_par);
	}
	
	// Particle redistribution summary time display
	printf("<+> Number of particle after remesh     : %8d\n", _par.num);
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        span = std::chrono::system_clock::now() - tick;
        _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<-> Particle redistribution calculation \n");
	printf("    comp. time:                        [%f s]\n", _time);

	return;
}


/**
 *  @brief  Particle neighbor search initialization. In this method the particle neighbor ID
 *  toward itself is evaluated. The data of current particle are now completed, assign all data 
 *  to the base particle for further particle redistribution.
 *         
 *  @param  _particle  Particle data container for neighbor evaluation.
 *  @param  _baseGrid  Grid node data container used for particle neighbor evaluation.
*/
void remeshing::set_neighbor(Particle &_par, const GridNode &_gridNode){
	/* Procedure
		[1] Evaluate the particle neighbor toward the data itself
		[2] Assign the particle data to Base Particle
	*/
	printf("\n+-------------- Neighbor Evaluation --------------+\n");

	// Initialize the cell List only for neighbor type 2 [CELL LIST only]
	if (Pars::opt_neighbor == 2){this->neighbor_step.cell_list_init(_par);}

	// [1] Evaluate the neighbor list
	this->neighbor_step.neigbor_search(_par, _gridNode);	
	
	// [2] Assign the particle base
	this->particleBase = new Particle;
	*(this->particleBase) = _par;
	
	printf("+-------------------------------------------------+\n");
}