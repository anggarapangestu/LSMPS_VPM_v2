#include "velocity_calc.hpp"
#include "velocity_biot_savart.hpp"
#include "../FMM/fmm2D.hpp"
#include "../FMM/fmm3D.hpp"
#include "../LSMPS/LSMPSa.hpp"
#include "../LSMPS/LSMPSb.hpp"
#include "../grid_block/generateGrid.hpp"
#include "../grid_block/gridNodeNgh.hpp"


#include "../neighbor/neighbor.hpp"
#include "../save_data/save_data.hpp"

// A debugger to check the neighbor evaluation after adaptation and redistribution
void save_ngh_all(const Particle &_myPar, int ID);      // For DEBUGGING

/** The type of 2D simulation FMM calculation
 *   [1] : Old fashion code
 *   [2] : New code using tree code built in unordered map of cell pointer
 *   [3] : Poisson solver using LSMPS on Stream function
 *   [4] : Poisson solver using LSMPS on Velocity
 *   [5] : Multiresolution poisson solver using LSMPS
 *   [6] : Multiresolution poisson solver using LSMPS (Modification)
*/
#define VELOCITY_2D_TYPE 3

/** The type of 3D simulation FMM calculation
 *   [1] : Refactor of old code (tree data built in vector, turns out to be 1.5 times faster)
 *   [2] : New code using tree code built in unordered map of cell pointer
*/
#define VELOCITY_3D_TYPE 1

/** Flag for activation of boundary evaluation
 *   0:= Only take the equation inside domain, 
 *   1:= Take boundary condition into account
*/
#define BOUNDARY_ACTIVE 1

#define ALL_PARTICLE true	// Code construction is still on going

/**
 *  @brief Velocity calculation manager.
 *  NOTE: Works well both for 2D and 3D simulation.
 *  
 *  @param	_particle  The particle for velocity evaluation.
 *  @param  _step  The current simulation iteration.
 */
void VelocityCalc::get_velocity(Particle &p, const int step)
{
	// Print log to console
	printf("\nCalculate Velocity ...\n");
	printf("<+> Calculate the rotational velocity term\n");

	// Computational time manager
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

	// **Calculate the rotational part
	if (DIM == 2){
		// Resize the velocity of particle
		p.u.clear(); p.u.resize(p.num,0.0e0);
		p.v.clear(); p.v.resize(p.num,0.0e0);

		switch (VELOCITY_2D_TYPE){
		case 1:	// Old fashing FMM code of velocity calculation
			printf("%s<+> Type 1: Old based velocity calculation %s\n", FONT_CYAN, FONT_RESET);
			this->velocity_old(p, step);
			break;
		case 2:	// FMM based velocity calculation
			printf("%s<+> Type 2: Unordered map tree cell pointer %s\n", FONT_CYAN, FONT_RESET);
			this->velocity_fmm_2d(p, step);
			break;
		case 3:	// LSMPS based poisson solver velocity calculation
			printf("%s<+> Type 3: LSMPS poisson solver %s\n", FONT_CYAN, FONT_RESET);
			this->velocity_LSMPS_poisson_2d(p, step);
			break;
		case 4:
			break;
		default:
			break;
		}
	}

	else if (DIM == 3){
		// Resize the velocity of particle
		p.u.clear(); p.u.resize(p.num,0.0e0);
		p.v.clear(); p.v.resize(p.num,0.0e0);
		p.w.clear(); p.w.resize(p.num,0.0e0);

		switch (VELOCITY_3D_TYPE){
		case 1:	// FMM based velocity calculation using basic tree cell
			printf("%s<+> Type 1: Vector tree cell data %s\n", FONT_CYAN, FONT_RESET);
			this->velocity_fmm_3d_fast(p, step);
			break;
		case 2:	// FMM based velocity calculation using pointer tree cell
			printf("%s<+> Type 2: Unordered map tree cell pointer %s\n", FONT_CYAN, FONT_RESET);
			this->velocity_fmm_3d(p, step);
			break;
		default:
			break;
		}
	}
	
	
	// **Calculate the irrotational part
	// Helmholtz decomposition (u = u_rot + u_irrot)
	printf("<+> Adding the irrotational velocity term \n");
	printf("    [helmholtz backward decomposition]\n");
	#pragma omp parallel for
	for (int i = 0; i < p.num; i++){
		p.u[i] += Pars::u_inf;
		p.v[i] += Pars::v_inf;
	}

	if (DIM == 3){
	#pragma omp parallel for
	for (int i = 0; i < p.num; i++){
		p.w[i] += Pars::w_inf;
	}
	}

	// Display computational time
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<-> Velocity calculation comp. time:   [%f s]\n", _time);

	return;
}

/**
 *  @brief Velocity calculation manager for multiresolution.
 *  NOTE: Works well both for 2D and 3D simulation.
 *  
 *  @param	_particle  The particle for velocity evaluation.
 *  @param	_baseGrid  The node grid of particle data.
 *  @param  _step  The current simulation iteration.
 */
void VelocityCalc::get_velocity(Particle &p, const GridNode &bGrd, const int step)
{
	// Print log to console
	printf("\nCalculate Velocity ...\n");
	printf("<+> Calculate the rotational velocity term\n");

	// Computational time manager
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

	// **Calculate the rotational part
	if (DIM == 2){
		// Resize the velocity of particle
		p.u.clear(); p.u.resize(p.num,0.0e0);
		p.v.clear(); p.v.resize(p.num,0.0e0);

		switch (VELOCITY_2D_TYPE){
		case 1:	// Old fashing FMM code of velocity calculation
			printf("%s<+> Type 1: Old based velocity calculation %s\n", FONT_CYAN, FONT_RESET);
			this->velocity_old(p, step);
			break;
		case 2:	// FMM based velocity calculation
			printf("%s<+> Type 2: Unordered map tree cell pointer %s\n", FONT_CYAN, FONT_RESET);
			this->velocity_fmm_2d(p, step);
			break;
		case 3:	// LSMPS based poisson solver velocity calculation
			printf("%s<+> Type 3: LSMPS poisson on stream function%s\n", FONT_CYAN, FONT_RESET);
			this->velocity_LSMPS_poisson_2d(p, step);
			// this->velocity_LSMPS_poisson_2d_mod(p, step);
			break;
		case 4:	// LSMPS based poisson solver velocity calculation
			printf("%s<+> Type 4: LSMPS poisson on velocity %s\n", FONT_CYAN, FONT_RESET);
			this->velocity_LSMPS_poisson_2d_vel(p, step);
			break;
		case 5:	// LSMPS based poisson solver velocity calculation w/ multiresolution modification
			printf("%s<+> Type 5: LSMPS poisson solver multires %s\n", FONT_CYAN, FONT_RESET);
			this->velocity_LSMPS_poisson_2d_mres(p, bGrd, step);
			break;
		case 6:	// A simple test for fast debugger
			printf("%s<+> Type 6: LSMPS poisson multires debug %s\n", FONT_CYAN, FONT_RESET);
			this->velocity_LSMPS_poisson_2d_mres_2(p);
			break;
		default:
			break;
		}
	}

	else if (DIM == 3){
		// Resize the velocity of particle
		p.u.clear(); p.u.resize(p.num,0.0e0);
		p.v.clear(); p.v.resize(p.num,0.0e0);
		p.w.clear(); p.w.resize(p.num,0.0e0);

		switch (VELOCITY_3D_TYPE){
		case 1:	// FMM based velocity calculation using basic tree cell
			printf("%s<+> Type 1: Vector tree cell data %s\n", FONT_CYAN, FONT_RESET);
			this->velocity_fmm_3d_fast(p, step);
			break;
		case 2:	// FMM based velocity calculation using pointer tree cell
			printf("%s<+> Type 2: Unordered map tree cell pointer %s\n", FONT_CYAN, FONT_RESET);
			this->velocity_fmm_3d(p, step);
			break;
		default:
			break;
		}
	}
	
	
	// **Calculate the irrotational part
	// Helmholtz decomposition (u = u_rot + u_irrot)
	printf("<+> Adding the irrotational velocity term \n");
	printf("    [helmholtz backward decomposition]\n");
	#pragma omp parallel for
	for (int i = 0; i < p.num; i++){
		p.u[i] += Pars::u_inf;
		p.v[i] += Pars::v_inf;
	}

	if (DIM == 3){
	#pragma omp parallel for
	for (int i = 0; i < p.num; i++){
		p.w[i] += Pars::w_inf;
	}
	}

	// Display computational time
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<-> Velocity calculation comp. time:   [%f s]\n", _time);

	return;
}

// =====================================================
// +-------------- ADDITIONAL CALCULATOR --------------+
// =====================================================
// Testing LSMPS Poisson
void VelocityCalc::LSMPS_poisson_2d(Particle &p, const int _step){
	/**
	 * DATA:
	 *  - Coordinate p.x, p.y, p.z
	 *  - Size		 p.s
	 *  - Neighbor   p.neighbor
	 *  - Source     p.vorticity
	 *  - Boundary   p.isBoundary, p.boundaryVal 
	 *  - Analytics  p.chi ?
	 *  - Result     p.gz  ?
	*/

	save_data save_manager;
	std::ofstream _write;
	double T1, T2;
	
	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		double _time = omp_get_wtime();
	#endif
	
    // PROCEDURE 1: Data initialization
    // ********************************************************************
	std::vector<std::vector<double>> particle_POS (DIM);	// Particle position
	std::vector<double> RHS = p.vorticity;	// The particle source value container (Poisson)
	// std::vector<double> RHS(p.num, 0.0);	// The particle source value container (Laplace)
	std::vector<double> psi;

	// Assign the particle data into container
	basis_loop(d) particle_POS[d].resize(p.num, 0.0e0);
	#pragma omp parallel for
	for (int i = 0; i < p.num; i++){
		particle_POS[0][i] = p.x[i];
		particle_POS[1][i] = p.y[i];
		
		// Assign the value of right hand side
		#if (BOUNDARY_ACTIVE == 1)
			// if (p.isBoundary[i] == true) RHS[i] = p.chi[i];
			if (p.boundaryLoc[i] != 0) RHS[i] = p.boundaryVal[i];
		#endif
	}
	
	// PROCEDURE 2: Generate the global matrix
    // ********************************************************************
	// Generate the global matrix 
	//  * At initial of simulation step 0
	//  * At initial of resume step (resume_step + 1)
	//  * When adaptation is done
	// if ((p.isAdapt == true) || 
	// 	(Pars::opt_start_state == 0 && _step == 0) ||
	// 	(Pars::opt_start_state == 1 && _step == Pars::resume_step+1))
	// {
		// #if (TIMER_PAR == 0)
		// 	// Timer using super clock (chrono)
		// 	std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
		// #elif (TIMER_PAR == 1)
		// 	// Timer using paralel package
		// 	double _time = omp_get_wtime();
		// #endif

		// Generate the global matrix
		printf("Generate the global matrix ... \n");
		poissonSolver.create_global_matrix(particle_POS, p.neighbor, p.s, p.boundaryLoc /*p.isBoundary*/);

		#if (TIMER_PAR == 0)
			// Timer using super clock (chrono)
			std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
			double _time = span.count();
		#elif (TIMER_PAR == 1)
			// Timer using paralel package
			_time = omp_get_wtime() - _time;
		#endif
		// printf("<-> Global matrix computation time  :  [%f s]\n", _time);
		printf("<+> Initialize global matrix     :  %f s\n", _time);
	// }

	// Update timer container
	T1 = _time;

	// PROCEDURE 3: Solve the poisson equation
    // ********************************************************************
	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		tick = std::chrono::system_clock::now();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		// double _time = omp_get_wtime();
		_time = omp_get_wtime();
	#endif

	// Solve the global matrix
	printf("\nSolve the potential of poisson equation ... \n");
	poissonSolver.solve(RHS);
	
	// Get the stream function calculated from poisson equation
	p.gz = poissonSolver.get_phi();

	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
		double _time = span.count();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime() - _time;
	#endif
	// printf("<-> Poisson solver computation time :  [%f s]\n", _time);
	printf("<+> Step 1: Solve global matrix  :  %f s\n", _time);

	// Update timer container
	T2 = _time;

	_write.open(save_manager.get_log_directory(), std::ofstream::out | std::ofstream::app);
	_write << "Matrix Generation Time             : "  << wR << T1 << "\n";
	_write << "Solving Matrix Time                : "  << wR << T2 << "\n";
	_write << "Total computational Time           : "  << wR << T1+T2 << "\n";
	_write.close();
	
	return;
}

// =====================================================
// +------------- VELOCITY CALCULATOR 2D --------------+
// =====================================================
// #pragma region VELOCITY_2D

/**
 *  @brief Old velocity calculation. <!> Please check to any reference <!>.
 *  NOTE: This code supposed to work well, but the content may not robust.
 *  
 *  @param	_particle  The particle for velocity evaluation.
 *  @param  _step  The current simulation iteration.
*/
void VelocityCalc::velocity_old(Particle &p, const int step){
	// Internal variables
	std::vector<int> _index;	// The original index container
	Particle _particle;		  	// store current active particles box (to increase robustness)
	Particle _particleDense;  	// store particle with core size as smallest from current active particles
	
	// Initialize the 
	_particle.num = 0;
	_particleDense.num = 0;

	
	// TODO: store current active particles
	for (int i = 0; i < p.num; i++){
		//if (p.isActive[i] == true){
			_particle.x.push_back(p.x[i]);
			_particle.y.push_back(p.y[i]);
			_particle.s.push_back(p.s[i]);
			_particle.gz.push_back(p.gz[i]);
			_particle.u.push_back(0.0e0);
			_particle.v.push_back(0.0e0);
			_index.push_back(i);

			_particle.num++;
			_particleDense.x.push_back(p.x[i]);
			_particleDense.y.push_back(p.y[i]);
			_particleDense.s.push_back(p.s[i]);
			_particleDense.gz.push_back(p.gz[i]);
			_particleDense.u.push_back(p.u[i]);
			_particleDense.v.push_back(p.v[i]);
			_particleDense.num++;
		//}
	}

	// The tools for FMM calculation
	velocity_biot_savart FMM;
	
	// Define the parameter for FMM calculation
	int _iCutoff = Pars::icutoff;
	int _nS = Pars::n_s;
	int _nInter = 1;
	int _ndp = Pars::ndp;
	
	FMM.biotsavart_fmm_2d( _particle, _particleDense,  _iCutoff, _nS,  _nInter,  _ndp);
	//FMM.biotsavart_direct_2d(_particle, _particleDense);

	// Update base particle velocity
	for (int i = 0; i < _particle.num; i++){
		// Alias to the original index
		const int &ori_ID = _index[i];
		p.u[ori_ID] = _particle.u[i];
		p.v[ori_ID] = _particle.v[i];
	}

	return;
}

/**
 *  @brief FMM based velocity calculation.
 *  NOTE: Works well only for 2D simulation.
 *  
 *  @param	_particle  The particle for velocity evaluation.
 *  @param  _step  The current simulation iteration.
*/
void VelocityCalc::velocity_fmm_2d(Particle &p, const int step){
	/* Procedure:
       1. Prepare the data for velocity calculation (by FMM)
       2. Create the tree cell data
       3. Calculate the velocity by FMM
       4. Update the particle data
    */

    // PROCEDURE 1: Data initialization
    // ********************************************************************
	std::vector<std::vector<double>> particle_POS	// Particle position
	(p.num, std::vector<double>(DIM,0.0));
	std::vector<double> particle_SRC(p.num, 0.0);	// The particle source value container
	std::vector<bool> particle_mark(p.num, false);	// The particle active mark container
	
	// Assign the particle data into container
	for (int i = 0; i < p.num; i++){
		// Particle position
		particle_POS[i][0] = p.x[i];
		particle_POS[i][1] = p.y[i];

		// Particle source
		particle_SRC[i] = - p.vorticity[i] * std::pow(p.s[i], 2.0) / (2 * M_PI);  // Vortex strength
		
		// Particle mark
		particle_mark[i] = p.isActive[i];
	}

	
	// PROCEDURE 2: Create tree cell data
    // ********************************************************************
	// Initialize the tree
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif
	
	// Initialization Tree Map
	if (step == 0){
		printf("<+> Initialize tree ... \n");
		this->treeData.initializeTree(particle_POS, particle_mark);
	}else{
		printf("<+> Update tree ... \n");
		this->treeData.updateTree(particle_POS, particle_mark);
	}

	// [!] Additional information
	printf("    >> Leaf min level %2d\n", this->treeData.min_level);
	printf("    >> Leaf max level %2d\n", this->treeData.max_level);

	// // [DEBUGGING LINE]
    // this->treeData.saveLeafTree(treeData, std::to_string(step));
    // this->treeData.saveTree(treeData, std::to_string(step));
	
	// // A debug line to investigate the mozaic pattern problem source
	// std::vector<int> _group1 = {809,207,50,11};
	// std::vector<int> _group2 = {194,51,12};
	// std::vector<int> _group3 = {13600,3386};
	// std::vector<int> _list1,_list2,_list3,_list4;
	// for (const int &ID : _group1){
	// 	this->treeData.intList(ID,_list1,_list3);
	// 	this->treeData.extList(ID,_list2,_list4);
	// 	this->treeData.saveSelTree(treeData, std::to_string(ID)+"_L1_", _list1);
	// 	this->treeData.saveSelTree(treeData, std::to_string(ID)+"_L2_", _list2);
	// 	this->treeData.saveSelTree(treeData, std::to_string(ID)+"_L3_", _list3);
	// 	this->treeData.saveSelTree(treeData, std::to_string(ID)+"_L4_", _list4);
	// }
	// ERROR_LOG << "DATA 1 has just wrote\n";
	// for (const int &ID : _group2){
	// 	this->treeData.intList(ID,_list1,_list3);
	// 	this->treeData.extList(ID,_list2,_list4);
	// 	this->treeData.saveSelTree(treeData, std::to_string(ID)+"_L1_", _list1);
	// 	this->treeData.saveSelTree(treeData, std::to_string(ID)+"_L2_", _list2);
	// 	this->treeData.saveSelTree(treeData, std::to_string(ID)+"_L3_", _list3);
	// 	this->treeData.saveSelTree(treeData, std::to_string(ID)+"_L4_", _list4);
	// }
	// ERROR_LOG << "DATA 2 has just wrote\n";
	// for (const int &ID : _group3){
	// 	this->treeData.intList(ID,_list1,_list3);
	// 	this->treeData.extList(ID,_list2,_list4);
	// 	this->treeData.saveSelTree(treeData, std::to_string(ID)+"_L1_", _list1);
	// 	this->treeData.saveSelTree(treeData, std::to_string(ID)+"_L2_", _list2);
	// 	this->treeData.saveSelTree(treeData, std::to_string(ID)+"_L3_", _list3);
	// 	this->treeData.saveSelTree(treeData, std::to_string(ID)+"_L4_", _list4);
	// }
	// ERROR_LOG << "DATA 3 has just wrote\n";
	// throw std::exception();

	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<+> Tree finished in       : %f s\n", _time);
	
	
	// PROCEDURE 3: Calculate the velocity (by FMM)
    // ********************************************************************
	// The FMM operator class
	fmm2D FMM_tool;
	
	// Calculate the FMM accelerated method
	FMM_tool.calcField(this->treeData, particle_POS, particle_mark, particle_SRC);
	
	// Get the final result
	std::vector<double> Ex = FMM_tool.get_Field_x();
	std::vector<double> Ey = FMM_tool.get_Field_y();
	
	// PROCEDURE 4: Update the particle velocity
    // ********************************************************************
	// Update the velocity
	for (int i = 0; i < p.num; i++){
		p.u[i] = Ey[i];
		p.v[i] = -Ex[i];
	}
	
	return;
}

/**
 *  @brief The velocity solution using LSMPS poisson solver on the stream function.
 *  
 *  @param	_particle  The particle for velocity evaluation.
 *  @param  _step  The current simulation iteration.
*/
void VelocityCalc::velocity_LSMPS_poisson_2d(Particle &p, const int _step){
	/* Procedure:
       1. Prepare the data for velocity calculation (by FMM)
       2. Generate the global matrix
       3. Solve the global matrix
       4. Update the particle data
    */

    // PROCEDURE 1: Data initialization
    // ********************************************************************
	std::vector<std::vector<double>> particle_POS (DIM);	// Particle position
	// std::vector<double> RHS = p.vorticity;	// The particle source value container
	std::vector<double> RHS(p.num, 0.0);		// The particle source value container
	std::vector<double> psi;

	// Assign the particle data into container
	basis_loop(d) particle_POS[d].resize(p.num, 0.0e0);
	#pragma omp parallel for
	for (int i = 0; i < p.num; i++){
		particle_POS[0][i] = p.x[i];
		particle_POS[1][i] = p.y[i];
		
		// Assign the value of right hand side
		#if (BOUNDARY_ACTIVE == 0)
			RHS[i] = -p.vorticity[i];
		#elif (BOUNDARY_ACTIVE == 1)
			// if (p.isBoundary[i] == true)	RHS[i] = 0; // p.gz[i];  // p.boundaryVal[i];
			if (p.boundaryLoc[i] != 0)
				RHS[i] = 0.0; // p.boundaryVal[i];
			else
				RHS[i] = -p.vorticity[i];
		#endif
	}
	
	// PROCEDURE 2: Generate the global matrix
    // ********************************************************************
	// Generate the global matrix 
	//  * At initial of simulation step 0
	//  * At initial of resume step (resume_step + 1)
	//  * When adaptation is done
	if ((p.isAdapt == true) || 
		(Pars::opt_start_state == 0 && _step == 0) ||
		(Pars::opt_start_state == 1 && _step == Pars::resume_step+1))
	{
		#if (TIMER_PAR == 0)
			// Timer using super clock (chrono)
			std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
		#elif (TIMER_PAR == 1)
			// Timer using paralel package
			double _time = omp_get_wtime();
		#endif

		// Generate the global matrix
		printf("Generate the global matrix ... \n");
		poissonSolver.create_global_matrix(particle_POS, p.neighbor, p.s, p.boundaryLoc /*p.isBoundary*/);

		#if (TIMER_PAR == 0)
			// Timer using super clock (chrono)
			std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
			double _time = span.count();
		#elif (TIMER_PAR == 1)
			// Timer using paralel package
			_time = omp_get_wtime() - _time;
		#endif
		// printf("<-> Global matrix computation time  :  [%f s]\n", _time);
		printf("<+> Initialize global matrix     :  %f s\n", _time);
	}


	// PROCEDURE 3: Solve the poisson equation
    // ********************************************************************
	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		double _time = omp_get_wtime();
	#endif

	// Solve the global matrix
	printf("\nSolve the potential of poisson equation ... \n");
	poissonSolver.solve(RHS);
	
	// Get the stream function calculated from poisson equation
	psi = poissonSolver.get_phi();

	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
		double _time = span.count();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime() - _time;
	#endif
	// printf("<-> Poisson solver computation time :  [%f s]\n", _time);
	printf("<+> Step 1: Solve global matrix  :  %f s\n", _time);
	
	
	// PROCEDURE 4: Solve the poisson equation
    // ********************************************************************
	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		tick = std::chrono::system_clock::now();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime();
	#endif
	
	// Calculate the differential of stream function
	printf("Calculate the rotational part of velocity ... \n");
	LSMPSa LSMPSpsi;
	LSMPSpsi.set_LSMPS(p.x, p.y, p.s, psi, p.neighbor);
	std::vector<double> dPsidx = LSMPSpsi.get_ddx();
	std::vector<double> dPsidy = LSMPSpsi.get_ddy();

	// Reserve the data for phi
	p.P = std::vector<double>(p.num);

	// Update the velocity
	for (int i = 0; i < p.num; i++){
		p.u[i] = dPsidy[i];
		p.v[i] = -dPsidx[i];
		p.P[i] = psi[i];
	}

	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		span = std::chrono::system_clock::now() - tick;
		_time = span.count();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime() - _time;
	#endif
	// printf("<-> Velocity calculation comp. time :  [%f s]\n", _time);
	printf("<+> Step 2: Compute velocity     :  %f s\n", _time);

	return;
}

/**
 *  @brief The velocity solution using LSMPS poisson solver on the stream function.
 *  
 *  @param	_particle  The particle for velocity evaluation.
 *  @param  _step  The current simulation iteration.
*/
void VelocityCalc::velocity_LSMPS_poisson_2d_vel(Particle &p, const int _step){
	/* Procedure:
       1. Prepare the data for velocity calculation (by FMM)
       2. Generate the global matrix
       3. Solve the global matrix
       4. Update the particle data
    */
   	
	// PROCEDURE 0: Calculate the poisson source
    // ********************************************************************
	// Here the poisson equation for 2D
	//    lap(u) = -dw/dy
	//    lap(v) = dw/dx
	// Calculate the vorticity differential
	LSMPSa lsmps_tool;
	lsmps_tool.set_LSMPS(p.x, p.y, p.s, p.vorticity, p.neighbor);
	std::vector<double> dwdx = lsmps_tool.get_ddx();
	std::vector<double> dwdy = lsmps_tool.get_ddy();


    // PROCEDURE 1: Data initialization
    // ********************************************************************
	std::vector<std::vector<double>> particle_POS (DIM);	// Particle position
	// std::vector<double> RHS = p.vorticity;	// The particle source value container
	std::vector<double> RHS_v(p.num, 0.0);		// The particle source value container
	std::vector<double> RHS_u(p.num, 0.0);		// The particle source value container

	// Note on boundary condition (All dirichlet)
	//  > Farfield      : U = 0, V = 0
	//  > Wall          : U = 0, V = 0
	//  > Inlet/Outlet  : U = U_stream, V = 0

	// Assign the particle data into container
	basis_loop(d) particle_POS[d].resize(p.num, 0.0e0);
	#pragma omp parallel for
	for (int i = 0; i < p.num; i++){
		particle_POS[0][i] = p.x[i];
		particle_POS[1][i] = p.y[i];
		
		// Assign the value of right hand side
		#if (BOUNDARY_ACTIVE == 0)
			RHS_u[i] = -dwdy[i];	// lap(u) = -dw/dy
			RHS_v[i] = dwdx[i];		// lap(v) = dw/dx
		#elif (BOUNDARY_ACTIVE == 1)
			// Set the RHS for u poisson and v poisson equation
			if (p.boundaryLoc[i] == 0){
				// Inside domain
				RHS_u[i] = -dwdy[i];	// lap(u) = -dw/dy
				RHS_v[i] = dwdx[i];		// lap(v) = dw/dx
			}
			else{
				// Domain boundary 
				RHS_u[i] = p.vortx[i];  // 0.0
				RHS_v[i] = p.vorty[i];  // 0.0
			}
		#endif
	}
	
	// PROCEDURE 2: Generate the global matrix
    // ********************************************************************
	// Generate the global matrix 
	//  * At initial of simulation step 0
	//  * At initial of resume step (resume_step + 1)
	//  * When adaptation is done
	if ((p.isAdapt == true) || 
		(Pars::opt_start_state == 0 && _step == 0) ||
		(Pars::opt_start_state == 1 && _step == Pars::resume_step+1))
	{
		#if (TIMER_PAR == 0)
			// Timer using super clock (chrono)
			std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
		#elif (TIMER_PAR == 1)
			// Timer using paralel package
			double _time = omp_get_wtime();
		#endif

		// Generate the global matrix
		printf("Generate the global matrix ... \n");
		poissonSolver.create_global_matrix(particle_POS, p.neighbor, p.s, p.boundaryLoc /*p.isBoundary*/);

		#if (TIMER_PAR == 0)
			// Timer using super clock (chrono)
			std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
			double _time = span.count();
		#elif (TIMER_PAR == 1)
			// Timer using paralel package
			_time = omp_get_wtime() - _time;
		#endif
		// printf("<-> Global matrix computation time  :  [%f s]\n", _time);
		printf("<+> Initialize global matrix     :  %f s\n", _time);
	}


	// PROCEDURE 3: Solve the poisson equation
    // ********************************************************************
	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		double _time = omp_get_wtime();
	#endif

	// Solve the global matrix
	printf("\nSolve the x velocity of poisson equation ... \n");
	poissonSolver.solve(RHS_u);		// Solve the global matrix
	p.u = poissonSolver.get_phi();	// Get the velocity x

	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
		double _time = span.count();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime() - _time;
	#endif
	printf("<+> Velocity x comp. time        :  %f s\n", _time);


	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		tick = std::chrono::system_clock::now();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime();
	#endif

	// Solve the global matrix
	printf("\nSolve the y velocity of poisson equation ... \n");
	poissonSolver.solve(RHS_v);			// Solve the global matrix
	p.v = poissonSolver.get_phi();	// Get the velocity y

	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
		double _time = span.count();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime() - _time;
	#endif
	printf("<+> Velocity y comp. time        :  %f s\n", _time);
	
	return;
}

/**
 *  @brief The velocity solution using LSMPS poisson solver.
 *  NOTE: <?> Need to be check. Change the order of particle data point.
 *  This solver make a new order of particle before doing the solver.
 *  
 *  @param	_particle  The particle for velocity evaluation.
 *  @param  _step  The current simulation iteration.
*/
void VelocityCalc::velocity_LSMPS_poisson_2d_mod(Particle &p, const int _step){
	/* Procedure:
       1. Prepare the data for velocity calculation (by FMM)
       2. Generate the global matrix
       3. Solve the global matrix
       4. Update the particle data
    */

    // PROCEDURE 1: Data initialization
    // ********************************************************************
	std::vector<std::vector<double>> particle_POS (DIM);	// Particle position
	// std::vector<double> RHS = p.vorticity;	// The particle source value container
	std::vector<double> RHS(p.num, 0.0);		// The particle source value container
	std::vector<double> psi;
	
	Particle orderedPar;						// A new particle container for the ordered data
	orderedPar.num = 0;

	std::vector<int> indexBoundary;
	std::vector<int> indexNew(p.num);	// The index of old particle inside the new ordered container [original -> ordered]

	// [!] Rearrange the data order
	// Collect all particle inside the 
	for (int i = 0; i < p.num; i++){
		// if (p.isBoundary[i] == true){
		if (p.boundaryLoc[i] != 0){
			indexBoundary.push_back(i);
		}
	}
	
	// Separate the data by [Inside, Boundary]
	// Assign the Inside data
	for (int i = 0; i < p.num; i++){
		// if (p.isBoundary[i] == false){
		if (p.boundaryLoc[i] == 0){
			// Get the indexing
			indexNew[i] = orderedPar.num;

			// Update the data
			orderedPar.x.push_back(p.x[i]);
			orderedPar.y.push_back(p.y[i]);
			orderedPar.s.push_back(p.s[i]);
			orderedPar.vorticity.push_back(p.vorticity[i]);
			// orderedPar.isBoundary.push_back(p.isBoundary[i]);
			orderedPar.boundaryLoc.push_back(p.boundaryLoc[i]);
			orderedPar.boundaryVal.push_back(p.boundaryVal[i]);
			orderedPar.num++;
		}
	}

	// Assign the Boundary data
	for (const int &ID : indexBoundary){
		// Get the indexing
		indexNew[ID] = orderedPar.num;

		// Update the data
		orderedPar.x.push_back(p.x[ID]);
		orderedPar.y.push_back(p.y[ID]);
		orderedPar.s.push_back(p.s[ID]);
		orderedPar.vorticity.push_back(p.vorticity[ID]);
		// orderedPar.isBoundary.push_back(p.isBoundary[ID]);
		orderedPar.boundaryLoc.push_back(p.boundaryLoc[ID]);
		orderedPar.boundaryVal.push_back(p.boundaryVal[ID]);
		orderedPar.num++;
	}

	// Update the neighbor data
	orderedPar.neighbor.resize(p.num);
	for (int i = 0; i < p.num; i++){
		// Aliasing the ID
		const int &oriID = i;
		const int &newID = indexNew[oriID];

		// Check the neighbor and put into the new ID
		for (size_t j = 0; j < p.neighbor[oriID].size(); j++){
			// Aliasing the neighbor ID
			const int &oriNghID = p.neighbor[oriID][j];
			const int &newNghID = indexNew[oriNghID];

			// Assign the neighbor
			orderedPar.neighbor[newID].push_back(newNghID);
		}
	}


	// // // // [DEBUGGING TOOLS]
	// // // // =============================== START ===============================
	// save_data save_tool;
	// orderedPar.chi.resize(orderedPar.num, 0.0e0);
	// orderedPar.u.resize(orderedPar.num, 0.0e0);
	// orderedPar.v.resize(orderedPar.num, 0.0e0);
	// orderedPar.gx.resize(orderedPar.num, 0.0e0);
	// orderedPar.gy.resize(orderedPar.num, 0.0e0);
	// orderedPar.gz.resize(orderedPar.num, 0.0e0);
	// orderedPar.isActive.resize(orderedPar.num, false);
	// std::cout << "LOOKS GOOD up here\n";
	// save_tool.save_par_state(orderedPar, "Ordered_Data", 0);
	// throw std::exception();
	// // // // ================================ END ================================

	// Assign the particle data into container
	basis_loop(d) particle_POS[d].resize(p.num, 0.0e0);
	#pragma omp parallel for
	for (int i = 0; i < p.num; i++){
		particle_POS[0][i] = orderedPar.x[i];
		particle_POS[1][i] = orderedPar.y[i];
		RHS[i] = -orderedPar.vorticity[i];
	}

	// Assign the particle data into container
	basis_loop(d) particle_POS[d].resize(p.num, 0.0e0);
	#pragma omp parallel for
	for (int i = 0; i < p.num; i++){
		particle_POS[0][i] = orderedPar.x[i];
		particle_POS[1][i] = orderedPar.y[i];
		RHS[i] = -orderedPar.vorticity[i];
	}
	
	// PROCEDURE 2: Generate the global matrix
    // ********************************************************************
	// Generate the global matrix 
	//  * At initial of simulation step 0
	//  * At initial of resume step (resume_step + 1)
	//  * When adaptation is done
	if ((p.isAdapt == true) || 
		(Pars::opt_start_state == 0 && _step == 0) ||
		(Pars::opt_start_state == 1 && _step == Pars::resume_step+1))
	{
		#if (TIMER_PAR == 0)
			// Timer using super clock (chrono)
			std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
		#elif (TIMER_PAR == 1)
			// Timer using paralel package
			double _time = omp_get_wtime();
		#endif

		// Generate the global matrix
		printf("Generate the global matrix ... \n");
		poissonSolver.create_global_matrix(particle_POS, orderedPar.neighbor, orderedPar.s, orderedPar.boundaryLoc /*orderedPar.isBoundary*/);

		#if (TIMER_PAR == 0)
			// Timer using super clock (chrono)
			std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
			double _time = span.count();
		#elif (TIMER_PAR == 1)
			// Timer using paralel package
			_time = omp_get_wtime() - _time;
		#endif
		// printf("<-> Global matrix computation time  :  [%f s]\n", _time);
		printf("<+> Initialize global matrix     :  %f s\n", _time);
	}


	// PROCEDURE 3: Solve the poisson equation
    // ********************************************************************
	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		double _time = omp_get_wtime();
	#endif

	// Solve the global matrix
	printf("\nSolve the potential of poisson equation ... \n");
	poissonSolver.solve(RHS);
	
	// Get the stream function calculated from poisson equation
	psi = poissonSolver.get_phi();

	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
		double _time = span.count();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime() - _time;
	#endif
	// printf("<-> Poisson solver computation time :  [%f s]\n", _time);
	printf("<+> Step 1: Solve global matrix  :  %f s\n", _time);
	
	
	// PROCEDURE 4: Solve the poisson equation
    // ********************************************************************
	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		tick = std::chrono::system_clock::now();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime();
	#endif

	// Reorder the psi
	std::vector<double> orderedPsi = psi;
	for (int ID = 0; ID < p.num; ID++){
		// Create alias to ID
		const int &newID = indexNew[ID];
		psi[ID] = orderedPsi[newID];
	}
	
	// Calculate the differential of stream function
	printf("Calculate the rotational part of velocity ... \n");
	LSMPSa LSMPSpsi;
	LSMPSpsi.set_LSMPS(p.x, p.y, p.s, psi, p.neighbor);
	std::vector<double> dPsidx = LSMPSpsi.get_ddx();
	std::vector<double> dPsidy = LSMPSpsi.get_ddy();

	// Update the velocity
	for (int i = 0; i < p.num; i++){
		p.u[i] = dPsidy[i];
		p.v[i] = -dPsidx[i];
		// p.v[i] = psi[i];
	}

	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		span = std::chrono::system_clock::now() - tick;
		_time = span.count();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime() - _time;
	#endif
	// printf("<-> Velocity calculation comp. time :  [%f s]\n", _time);
	printf("<+> Step 2: Compute velocity     :  %f s\n", _time);

	return;
}

/**
 *  @brief The velocity solution using LSMPS poisson solver.
 *  NOTE: <?> Need to be check. (Using patch only for 2 resolution)
 *  
 *  @param	_particle  The particle for velocity evaluation.
*/
void VelocityCalc::velocity_LSMPS_poisson_2d_mres_2(Particle &p){
	/* Procedure:
       1. Prepare the data for velocity calculation (by FMM)
       2. Generate the global matrix
       3. Solve the global matrix
       4. Update the particle data
    */

    // PROCEDURE 0: Data collection
    // ********************************************************************
	// Create a base particle only for 

	// ====== Copy of parameter from initialization ====== //
    // Set the particle division type
    int type = 1;               // 0:= split in x direction; 1:= split in y direction
    bool up_region = true;      // 0:= divide the lower coordinate region; 1:= divide the higher coordinate region
    double split_loc = -0.5;    // The location of splitting

	// The additional data
	Particle PatchLow, PatchHigh;

	// Generate the particle inside patch region
	{
		// A reference to each particle
		Particle *Patch1, *Patch2;
		if (up_region){
			Patch1 = &PatchHigh;
			Patch2 = &PatchLow;
		}else{
			Patch1 = &PatchLow;
			Patch2 = &PatchHigh;
		}
		
		// The type two : splitted domain
        std::cout << "|A| A splitted domain distribution\n";
        
        // Number of particle in x and y direction
        int _nx, _ny;       // Number of particle in x and y direction
        int _nparticle;     // Number of all particle
        
        // Simulation domain pivot coordinate
        double pivX, pivY;  // Pivot location of the most left (x) and the most bottom (y)
		double _x, _y;	// Temporary data for x and y coordinate
        double size;	// Temporary data for size

        // [1] Generate the first region (low coordinate)
        // **********************************************
        size = up_region ? (Pars::sigma) : (2*Pars::sigma);
        if (type == 0){
            // Splitted in x direction
            // Update the number of particle in x and y direction
            _nx = 4;
            _ny = std::ceil(Pars::lydom / size);
            _nparticle = _nx * _ny;
            
            // Simulation domain pivot coordinate
            pivX = split_loc - (_nx*size);
            pivY = - (_ny*size) / 2.0;
        }
        else if (type == 1){
            // Splitted in y direction
            // Update the number of particle in x and y direction
            _nx = std::ceil(Pars::lxdom / size);
            _ny = 4;
            _nparticle = _nx * _ny;
            
            // Simulation domain pivot coordinate
            pivX = -Pars::xdom;
            pivY = split_loc - (_ny*size);
        }

        // Adjust the size of vector into the calculated particle number
        Patch1->num = _nparticle;
        Patch1->x.resize(Patch1->num, 0.0e0);
        Patch1->y.resize(Patch1->num, 0.0e0);
        Patch1->s.resize(Patch1->num, up_region ? (2*Pars::sigma) : (Pars::sigma));
        // Patch1->level.resize(Patch1->num, up_region ? Pars::max_level - 1 : Pars::max_level);

        // Generate the particle position
        int count = 0;
        for (int i = 0; i < _nx; i++){
            for (int j = 0; j < _ny; j++){
                _x = pivX + (i+0.5)*size;
                _y = pivY + (j+0.5)*size + Pars::mp_shift;   // Additional of unbalance: Pars::mp_shift
                Patch1->x[count] = _x;
                Patch1->y[count] = _y;
                count++;
            }
        }
        std::cout << "|X| Done generating the first (1st) part\n";


        // [1] Generate the second region (high coordinate)
        // ************************************************
        size = up_region ? 2*Pars::sigma : Pars::sigma;
        if (type == 0){
            // Splitted in x direction
            // Update the number of particle in x and y direction
            _nx = 4;
            _ny = std::ceil(Pars::lydom / size);
            _nparticle = _nx * _ny;
            
            // Simulation domain pivot coordinate
            pivX = split_loc;
            pivY = - (_ny*size) / 2.0;
        }
        else if (type == 1){
            // Splitted in y direction
            // Update the number of particle in x and y direction
            _nx = std::ceil(Pars::lxdom / size);
            _ny = 4;
            _nparticle = _nx * _ny;
            
            // Simulation domain pivot coordinate
            pivX = -Pars::xdom;
            pivY = split_loc;
        }

        // Adjust the size of vector into the calculated particle number
        Patch2->num = _nparticle;
        Patch2->x.resize(Patch2->num, 0.0e0);
        Patch2->y.resize(Patch2->num, 0.0e0);
        Patch2->s.resize(Patch2->num, up_region ? (Pars::sigma) : (2*Pars::sigma));
        // Patch2->level.resize(Patch2->num, up_region ? Pars::max_level : Pars::max_level - 1);

        // Generate the particle position
		count = 0;
        for (int i = 0; i < _nx; i++){
            for (int j = 0; j < _ny; j++){
                _x = pivX + (i+0.5)*size;
                _y = pivY + (j+0.5)*size + Pars::mp_shift;   // Additional of unbalance: Pars::mp_shift
                Patch2->x[count] = _x;
                Patch2->y[count] = _y;
                count++;
            }
        }

        std::cout << "|X| Done generating the second (2nd) part\n";
	}

	// Interpolate the data of each particle
	neighbor ngh_tool;
	ngh_tool.neigbor_search(PatchLow, p);
	ngh_tool.neigbor_search(PatchHigh, p);

	// Interpolation
	LSMPSb LSMPScalc;
	LSMPScalc.set_LSMPS_2D(PatchHigh.x, PatchHigh.y, PatchHigh.s,
						   p.x, p.y, p.s, p.vorticity, PatchHigh.neighbor);
	PatchHigh.vorticity = LSMPScalc.get_d0();

	LSMPScalc.set_LSMPS_2D(PatchLow.x, PatchLow.y, PatchLow.s,
						   p.x, p.y, p.s, p.vorticity, PatchLow.neighbor);
	PatchLow.vorticity = LSMPScalc.get_d0();

	// Seperate and combine the data
	Particle dataLow, dataHigh;
	Particle dataCombine, baseCombine;

	// Update the neighbor of the combined data (4 times)
	dataCombine = p;
	dataCombine.num += PatchLow.num;
	for (int i = 0; i < PatchLow.num; i++){
		// Collect the data of low resolution
		dataCombine.x.push_back(PatchLow.x[i]);
		dataCombine.y.push_back(PatchLow.y[i]);
		dataCombine.s.push_back(Pars::sigma * 2);
		dataCombine.vorticity.push_back(PatchLow.vorticity[i]);
	}
	
	dataCombine.num += PatchHigh.num;
	for (int i = 0; i < PatchHigh.num; i++){
		// Collect the data of high resolution
		dataCombine.x.push_back(PatchHigh.x[i]);
		dataCombine.y.push_back(PatchHigh.y[i]);
		dataCombine.s.push_back(Pars::sigma);
		dataCombine.vorticity.push_back(PatchHigh.vorticity[i]);
	}

	// Split the original data from size
	std::vector<int> indexConvertLow, indexConvertHigh;
	
	// High resolution
	indexConvertHigh.clear();
	dataHigh.num = 0;
	double size = Pars::sigma;
	for (int i = 0; i < p.num; i++){
		double ratio = p.s[i] / size;
		if (ratio > 0.99 && ratio < 1.01){
			// Into the criteria
			indexConvertHigh.push_back(i);
			dataHigh.x.push_back(p.x[i]);
			dataHigh.y.push_back(p.y[i]);
			dataHigh.s.push_back(p.s[i]);
			dataHigh.num++;
		}
	}

	// Low resolution
	indexConvertLow.clear();
	dataLow.num = 0;
	size = Pars::sigma * 2;
	for (int i = 0; i < p.num; i++){
		double ratio = p.s[i] / size;
		if (ratio > 0.99 && ratio < 1.01){
			// Into the criteria
			indexConvertLow.push_back(i);
			dataLow.x.push_back(p.x[i]);
			dataLow.y.push_back(p.y[i]);
			dataLow.s.push_back(p.s[i]);
			dataLow.num++;
		}
	}

	// Update the neighbor data
	// Update the original data of high resolution
	baseCombine = dataCombine;
	size = Pars::sigma;
	for (int i = 0; i < baseCombine.num; i++){
		double ratio = baseCombine.s[i] / size;
		if (ratio < 0.99 || ratio > 1.01){
			// Outside the criteria
			baseCombine.x[i] = Pars::lxdom;
			baseCombine.y[i] = Pars::lxdom;
		}
	}
	ngh_tool.neigbor_search(dataHigh, baseCombine);
	// ngh_tool.neigbor_search(PatchLow, baseCombine);

	// // Update the original data of low resolution
	baseCombine = dataCombine;
	size = Pars::sigma * 2;
	for (int i = 0; i < baseCombine.num; i++){
		double ratio = baseCombine.s[i] / size;
		if (ratio < 0.99 || ratio > 1.01){
			// Outside the criteria
			baseCombine.x[i] = Pars::lxdom;
			baseCombine.y[i] = Pars::lxdom;
		}
	}
	ngh_tool.neigbor_search(dataLow, baseCombine);
	// ngh_tool.neigbor_search(PatchHigh, baseCombine);

	// Change the neighbor of the high and low resolution patch
	PatchHigh.s = std::vector<double>(PatchHigh.num, Pars::sigma);
	PatchLow.s = std::vector<double>(PatchLow.num, 2*Pars::sigma);
	ngh_tool.neigbor_search(PatchHigh, p);
	ngh_tool.neigbor_search(PatchLow, p);

	// Update the final data
	dataCombine.neighbor.clear();
	dataCombine.neighbor.resize(p.num);
	
	for (int i = 0; i < dataLow.num; i++){
		int id = indexConvertLow[i];
		dataCombine.neighbor[id] = dataLow.neighbor[i];
	}
	for (int i = 0; i < dataHigh.num; i++){
		int id = indexConvertHigh[i];
		dataCombine.neighbor[id] = dataHigh.neighbor[i];
	}

	for (int i = 0; i < PatchLow.num; i++){
		dataCombine.neighbor.push_back(PatchLow.neighbor[i]);
		// dataCombine.s[p.num + i] = Pars::sigma;
	}
	for (int i = 0; i < PatchHigh.num; i++){
		dataCombine.neighbor.push_back(PatchHigh.neighbor[i]);
		// dataCombine.s[p.num + PatchLow.num + i] = Pars::sigma * 2;
	}


	// // // // [DEBUGGING TOOLS]
	// // // // =============================== START ===============================
	// save_data save_tool;
	// int number = 100;
	// int interval = dataCombine.num / number;
    
    // for (int i = 0; i < number; i++){
	// 	int ID = interval * i;
    //     MESSAGE_LOG << "List of random ID " << std::abs(ID) << "\n";
    //     save_ngh_all(dataCombine, std::abs(ID));
    // }

	// throw std::exception();
	// // // // ================================ END ================================


	// // // [DEBUGGING TOOLS]
	// // // =============================== START ===============================
	save_data save_tool;
	save_tool.save_par_state(p, "All", 0);

	PatchHigh.chi.resize(PatchHigh.num, 0.0e0);
	PatchHigh.u.resize(PatchHigh.num, 0.0e0);
	PatchHigh.v.resize(PatchHigh.num, 0.0e0);
	PatchHigh.gz.resize(PatchHigh.num, 0.0e0);
	PatchHigh.isActive.resize(PatchHigh.num, false);
	save_tool.save_par_state(PatchHigh, "PatchHigh", 0);
	
	PatchLow.chi.resize(PatchLow.num, 0.0e0);
	PatchLow.u.resize(PatchLow.num, 0.0e0);
	PatchLow.v.resize(PatchLow.num, 0.0e0);
	PatchLow.gz.resize(PatchLow.num, 0.0e0);
	PatchLow.isActive.resize(PatchLow.num, false);
	save_tool.save_par_state(PatchLow, "PatchLow", 0);

	dataCombine.chi.resize(dataCombine.num, 0.0e0);
	dataCombine.u.resize(dataCombine.num, 0.0e0);
	dataCombine.v.resize(dataCombine.num, 0.0e0);
	dataCombine.gz.resize(dataCombine.num, 0.0e0);
	dataCombine.isActive.resize(dataCombine.num, false);
	save_tool.save_par_state(dataCombine, "Data_All", 0);
	// throw std::exception();
	// // // ================================ END ================================


	// PROCEDURE 1: Data initialization
    // ********************************************************************
	std::vector<std::vector<double>> particle_POS (DIM);	// Particle position
	// std::vector<double> RHS = p.vorticity;	// The particle source value container
	std::vector<double> RHS(dataCombine.num, 0.0);		// The particle source value container
	std::vector<double> psi;

	// Assign the particle data into container
	basis_loop(d) particle_POS[d].resize(dataCombine.num, 0.0e0);
	#pragma omp parallel for
	for (int i = 0; i < dataCombine.num; i++){
		particle_POS[0][i] = dataCombine.x[i];
		particle_POS[1][i] = dataCombine.y[i];
		RHS[i] = -dataCombine.vorticity[i];
	}
	
	// PROCEDURE 2: Generate the global matrix
    // ********************************************************************
	// Generate the global matrix 
	//  * At initial of simulation step 0
	//  * At initial of resume step (resume_step + 1)
	//  * When adaptation is done
	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		double _time = omp_get_wtime();
	#endif

	// Generate the global matrix
	printf("Generate the global matrix ... \n");
	poissonSolver.create_global_matrix(particle_POS, dataCombine.neighbor, dataCombine.s, dataCombine.boundaryLoc /*dataCombine.isBoundary*/);

	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
		double _time = span.count();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime() - _time;
	#endif
	// printf("<-> Global matrix computation time  :  [%f s]\n", _time);
	printf("<+> Initialize global matrix     :  %f s\n", _time);


	// PROCEDURE 3: Solve the poisson equation
    // ********************************************************************
	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime();
	#endif

	// Solve the global matrix
	printf("\nSolve the potential of poisson equation ... \n");
	poissonSolver.solve(RHS);
	
	// Get the stream function calculated from poisson equation
	psi = poissonSolver.get_phi();

	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
		double _time = span.count();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime() - _time;
	#endif
	// printf("<-> Poisson solver computation time :  [%f s]\n", _time);
	printf("<+> Step 1: Solve global matrix  :  %f s\n", _time);
	
	
	// PROCEDURE 4: Solve the poisson equation
    // ********************************************************************
	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		tick = std::chrono::system_clock::now();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime();
	#endif
	
	// Calculate the differential of stream function
	printf("Calculate the rotational part of velocity ... \n");
	LSMPSa LSMPSpsi;
	LSMPSpsi.set_LSMPS(p.x, p.y, p.s, psi, p.neighbor);
	std::vector<double> dPsidx = LSMPSpsi.get_ddx();
	std::vector<double> dPsidy = LSMPSpsi.get_ddy();

	// Update the velocity
	for (int i = 0; i < p.num; i++){
		p.u[i] = dPsidy[i];
		// p.v[i] = -dPsidx[i];
		p.v[i] = psi[i];
	}

	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		span = std::chrono::system_clock::now() - tick;
		_time = span.count();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime() - _time;
	#endif
	// printf("<-> Velocity calculation comp. time :  [%f s]\n", _time);
	printf("<+> Step 2: Compute velocity     :  %f s\n", _time);

	return;
}


/**
 *  @brief The velocity solution using LSMPS poisson solver.
 *  NOTE: <?> Need to be check. The serious MRES calculation.
 *  
 *  @param	_particle  The particle for velocity evaluation.
 *  @param	_baseGrid  The node grid of particle data.
 *  @param  _step  The current simulation iteration.
 */
void VelocityCalc::velocity_LSMPS_poisson_2d_mres(Particle &basePar, const GridNode &baseGrd, const int step){
	/**
	 *  NOTE: There is a moderate modification for multiresolution scheme
	 *  PS: This modification is an attempt to apply the LSMPS poisson equation 
	 * 		for multiresolution data
	 * 
	 *  The procedure is broken down as follow, before get into the procedure 
	 *   I want to discuss about the resolution interface.
	 * 
	 *  Supposed we have 3 resolution level in which the resolution only
	 *   jump within one level
	 *     __ __ __ __   ____  ____  ____   _________  _________ 
	 *    |__|__|__|__|-|    ||    ||    |-|         ||         |
	 *    |__|__|__|__|-|____||____||____|-|         ||         |
	 *    |__|__|__|__|-|    ||    ||    |-|         ||         |
	 *    |__|__|__|__|-|____||____||____|-|_________||_________|
	 *    |__|__|__|__|-|    ||    ||    |-|         ||         |
	 *    |__|__|__|__|-|____||____||____|-|         ||         |
	 *    |__|__|__|__|-|    ||    ||    |-|         ||         |
	 *    |__|__|__|__|-|____||____||____|-|_________||_________|
	 *         Res 1   |       Res 2      |        Res 3
	 *                 |              Interface 2 -> Contact region between
	 *                 |                              resolution 2 and 3
	 *            Interface 1    
	 *              -> Contact region between resolution 1 and 2
	 *   
	 * 
	 *  At each multiresolution interface, both the low and high resolution region are extended 
	 *   in one resolution size. The extension part both from low and high resolution will
	 *   overlap with the original distribution. See below: (Only the Interface 2)
	 * 
	 *                 Interface 2
	 *     ____  _________ | _________  _________      NOTE:
	 *    |    ||    |    |-|    |    ||         |       > Both the extension patch overlap
	 *    |____||____|____|-|____|    ||         |          the original region.
	 *    |    ||    |    |-|    |    ||         |       > The extension only add one resolution
	 *    |____||____|____|-|____|____||_________|          therefore, the extension of finer resolution
	 *    |    ||    |    |-|    |    ||         |          only takes half size of the coarse resolution.
	 *    |____||____|____|-|____|    ||         |          Otherwise, the coarse resolution takes twice
	 *    |    ||    |    |-|    |    ||         |          the size of finer resolution region.
	 *    |____||____|____|-|____|____||_________|
	 *           \___ ___/   \__/   
	 *               v         \---> The extension of Res 2 
	 *     The extension of           in the Res 3 region
	 *    Res 3 in Res 2 region
	 * 
	 *  IDEA:
	 *  By implementing the additional interface patch of the above scheme the number of particle
	 *   is quite dense near the interface. However, the Poisson solver modification is not about
	 *   adding particle density on interface, instead it creates additional particle to maintain
	 *   the particle spacing is (relatively) uniform for the particle that lies near the interface.
	 * 
	*/

	/**
	 * PROCEDURE:
	 *  1. Data initialization
	 *     a. Collect all nodes that near the resolution interface
	 *     b. Generate new patch particle
	 *     c. Interpolate the new patch particle
	 *     d. Update the neighbor connection
	 *     e. Re-arrange the Position, RHS, Size, and Neighbor container
	 *  2. Solve the global matrix [DIRECT]
	 *  3. Update particle data [Use the index mapping]
	 *     a. Calculate velocity from stream function [LSMPS A]
	 * 
	 * NOTE:
	 *  > The additional patch of particle data is saved to maintain data sequence in order to
	 *     reuse the constructed global matrix.
	 *  > Component used in calculation:
	 *     - Coordinate      : @lifetime until next distribution adaptation
	 *     - Neighbor        : @lifetime until next distribution adaptation
	 *     - Size            : @lifetime until next distribution adaptation
	 *     - RHS (Vorticity) : @lifetime one iteration
	*/

    // PROCEDURE 1: Data initialization
    // ********************************************************************
	/** Particle Data container
	 *    [1] Base particle (original data)
	 *    [2] Patched particle (ordered per node group)
	 * 
	 *  Group Data container
	 *    [1] Near Interface particle list (nIFaceID) -> For neighbor update
	 *    [2] Near Interface Node list
	 *   _________________		
	 *  |                 |		
	 *  |                 |		
	 *  |                 |
	 *  |_________________|
	 * 
	*/

	if (true)	// ((basePar.isAdapt == true) || 
	// 	(Pars::opt_start_state == 0 && step == 0) ||
	// 	(Pars::opt_start_state == 1 && step == Pars::resume_step+1))
	{
		// [COLLECTION 1]
		// Collect the near interface node
		std::unordered_map<int, bool> nIFaceFlag;	// The flag for face node ID
		poissonSolver.nIFaceNodeID.clear();

		// Find the new near interface node
		std::vector<bool> _tempFlag(basePar.num, false);// The flag for face node ID
		
		// Check each particle whether it have neighbor with different size or not
		#pragma omp parallel for
		for (int ID = 0; ID < basePar.num; ID++){
			// Alias to current particle size
			const double &currSize = basePar.s[ID];

			// Iterate each neighbor on current particle
			for (const int &nghID : basePar.neighbor[ID]){
				// Alias to neighbor particle size
				const double &nghSize = basePar.s[nghID];

				// Compare the size between particle
				double ratio = nghSize / currSize;
				if (ratio < 0.95 || ratio > 1.05){
					// The current neighbor have different size
					// TODO: Collect the node [Hold and put]
					_tempFlag[ID] = true;
					break;
				}
			}
		}

		// TODO: Collect the node
		for (int ID = 0; ID < basePar.num; ID++){
			if(_tempFlag[ID] == true){
				// Get the node ID occupies the current particle
				const int &_nodeID = basePar.nodeID[ID];
				if(nIFaceFlag.count(_nodeID) == 0){
					poissonSolver.nIFaceNodeID.push_back(_nodeID);
					nIFaceFlag.insert({_nodeID, true});
				}
			}
		}

		// Collect near interface particle
		// Reserve the other data
		poissonSolver.nIFaceParID.clear();
		poissonSolver.nIFaceParFlag.clear();
		poissonSolver.nIFaceParFlag.resize(basePar.num, false);

		// Collect all Reserve the other data
		for (const int &_nodeID : poissonSolver.nIFaceNodeID){
			// Exception check
			if(baseGrd.nodeMap.count(_nodeID) == 0){
				ERROR_LOG << "Node is not existed!\n";
				throw std::runtime_error("Particle hold a non existed node!!");
			}

			// Alias to the current node
			const Node *const &currNode = baseGrd.nodeMap.at(_nodeID);

			for (const int &ID : currNode->parList){
				poissonSolver.nIFaceParID.push_back(ID);
				poissonSolver.nIFaceParFlag[ID] = true;
			}
		}



		// [COLLECTION 2]
		// Collect the patched node
		std::unordered_map<int, bool> nPatchFlag;	// The flag for patched node ID
		poissonSolver.nPatchNodeID.clear();

		for (const int &_nodeID : poissonSolver.nIFaceNodeID){
			// Exception check
			if(baseGrd.nodeMap.count(_nodeID) == 0){
				ERROR_LOG << "Node is not existed!\n";
				throw std::runtime_error("Particle hold a non existed node!!");
			}

			// Find the neighbor node
			std::vector<int> nghNodeIDList;
			baseGrd.findNghLvl(nghNodeIDList, baseGrd.nodeMap.at(_nodeID));

			// Evaluate each neighbor
			for (const int &_nghID : nghNodeIDList){
				// // Exception of evaluating current node
				// if (_nghID == _nodeID) continue;

				// Check node existence
				if (baseGrd.nodeMap.count(_nghID) == 0){
					// Current node is not existed
					// Collect the data
					if (nPatchFlag.count(_nghID) == 0){
						// TYPE 1 : Finer resolution extension in coarse region
						poissonSolver.nPatchNodeID.push_back(_nghID);
						nPatchFlag.insert({_nghID, true});
					}
					continue;
				}

				// Check if current node not a leaf
				if (baseGrd.nodeMap.at(_nghID)->isLeaf != true){
					// Current node is not a leaf node
					// Collect the data
					if (nPatchFlag.count(_nghID) == 0){
						// TYPE 2 : Coarse resolution extension in finer region
						poissonSolver.nPatchNodeID.push_back(_nghID);
						nPatchFlag.insert({_nghID, true});
					}
					continue;
				}
			}
		}


		// [GENERATE PATCH DATA]
		// Reserve all patch particle data (additional resolution region)
		poissonSolver.patchPar.num = 0;
		poissonSolver.patchPar.x.clear();
		poissonSolver.patchPar.y.clear();
		poissonSolver.patchPar.s.clear();
		poissonSolver.patchPar.nodeID.clear();
		
		poissonSolver.nIFaceNodeParMap.clear();

		// The number of particle inside a node
		const int nodeParNum = Pars::intPow(baseGrd.baseParNum, DIM);
		int div[DIM] = {1, baseGrd.baseParNum};

		// Generate particle from node
		for (const int &_nodeID : poissonSolver.nPatchNodeID){
			// Internal variable
			double _parSize;
			double pivCoor[DIM];

			// Check whether the node is existed on the grid or not
			if (baseGrd.nodeMap.count(_nodeID) == 0){
				// Type 1 generation
				// Get the necessary data
				int _nodeLvl = baseGrd.getLevel(_nodeID);
				double _nodeSize = baseGrd.gridSize / Pars::intPow(2, _nodeLvl);
				_parSize = _nodeSize / baseGrd.baseParNum;

				// Obtain the pivot location by node index
				int _nodeIndex[DIM];
				baseGrd.ID2Index(_nodeIndex, _nodeID, _nodeLvl);

				// Update the pivot coordinate
				basis_loop(d) pivCoor[d] = baseGrd.pivotCoor[d] + _nodeIndex[d]*_nodeSize;
			}else{
				// Type 2 generation
				// Alias to the current node
				Node *const &currNode = baseGrd.nodeMap.at(_nodeID);

				// The size of particle
				_parSize = currNode->length / baseGrd.baseParNum;

				// Update the pivot coordinate
				basis_loop(d) pivCoor[d] = currNode->pivCoor[d];
			}

			// Generate all particle inside the current Node
			int parIndex[DIM];     // The particle local index (taken from the node)
			for (int locParID = 0; locParID < nodeParNum; locParID++){
				// Update the node map
				poissonSolver.nIFaceNodeParMap[_nodeID].push_back(poissonSolver.patchPar.num);

				// The local index coordinate inside the node
				basis_loop(d) parIndex[d] = (locParID/div[d]) % baseGrd.baseParNum;

				// Assign the particle position
				// The x position
				double _x = pivCoor[0] + (0.5 + parIndex[0])*_parSize;
				poissonSolver.patchPar.x.push_back(_x);
				
				// The y position
				double _y = pivCoor[1] + (0.5 + parIndex[1])*_parSize;
				poissonSolver.patchPar.y.push_back(_y);

				// Assign other data
				poissonSolver.patchPar.nodeID.push_back(_nodeID);
				poissonSolver.patchPar.s.push_back(_parSize);
				poissonSolver.patchPar.num++;
			}
		}
		// Done generating particle

	}

	// [INTERPOLATION] particle data 
	// Particle alias
	Particle &trgPar = poissonSolver.patchPar;
	Particle &srcPar = basePar;

	// Save the particle number
	std::vector<double> _sizeHold = trgPar.s;
	std::vector<int> _nodeIDHold = trgPar.nodeID;

	// Assign the particle corresponding node
    generateGrid grid_tool;
    std::unordered_map<int, std::vector<int>> nodeParMap;   // A container of target particle
    grid_tool.assignNodeID(baseGrd, nodeParMap, trgPar);

	// Evaluate the neighbor toward the source particle
    GridNodeNgh ngh_tool;
	std::vector<std::vector<int>> _tempNeighbor;
    ngh_tool.find_inter_neighbor(_tempNeighbor, trgPar, nodeParMap, baseGrd, srcPar);

	// Interpolate the necessary properties
	LSMPSb lsmpsCalc;
    lsmpsCalc.set_LSMPS_2D(trgPar.x, trgPar.y, trgPar.s,
                           srcPar.x, srcPar.y, srcPar.s, srcPar.vorticity, 
                           _tempNeighbor);
    trgPar.vorticity = lsmpsCalc.get_d0();

	// Retrieve back the size and node ID
	trgPar.s = _sizeHold;
	trgPar.nodeID = _nodeIDHold;

	
	// // // // [DEBUGGING TOOLS]
	// // // // =============================== START ===============================
	// save_data save_tool;
	// save_tool.save_par_state(basePar, "All", 0);

	// poissonSolver.patchPar.chi.resize(poissonSolver.patchPar.num, 0.0e0);
	// poissonSolver.patchPar.u.resize(poissonSolver.patchPar.num, 0.0e0);
	// poissonSolver.patchPar.v.resize(poissonSolver.patchPar.num, 0.0e0);
	// poissonSolver.patchPar.gz.resize(poissonSolver.patchPar.num, 0.0e0);
	// poissonSolver.patchPar.isActive.resize(poissonSolver.patchPar.num, false);
	// save_tool.save_par_state(poissonSolver.patchPar, "Patch", 0);
	// // baseGrd.saveSelectedGrid(baseGrd, poissonSolver.nIFaceNodeID, "Interface");
	// // baseGrd.saveLeafGrid(baseGrd, "Interface");
	// // // throw std::exception();
	// // // // ================================ END ================================

	// Data container
	const int totalNum = basePar.num + poissonSolver.patchPar.num;	// Combined particle number
	std::vector<std::vector<double>> comb_coord			// Combined particle position
		(DIM, std::vector<double>(totalNum, 0.0));
	std::vector<std::vector<int>> comb_ngh(totalNum);	// Combined particle neighbor list
	std::vector<double> comb_size(totalNum, 0.0);		// Combined particle size data
	std::vector<bool> comb_bnd_flag(totalNum, false);	// Combined particle boundary flag 		<!> Only declared, still not constructed <!>
	std::vector<int> comb_bnd_loc(totalNum, 0);			// Combined particle boundary location 	<!> Only declared, still not constructed <!>
	std::vector<double> comb_bnd_val(totalNum, 0.0);	// Combined particle boundary value		<!> Only declared, still not constructed <!>
	std::vector<double> RHS(totalNum, 0.0);		// The right hand side container (vorticity)
	std::vector<double> psi;					// Stream function [OUTPUT]

	// Assign the particle data into container
	#pragma omp parallel for
	for (int i = 0; i < basePar.num; i++){
		// Assign the coordinate
		comb_coord[0][i] = basePar.x[i];
		comb_coord[1][i] = basePar.y[i];

		// Assign the size
		comb_size[i] = basePar.s[i];

		// Assign the neighbor list
		comb_ngh[i].clear();
		if (poissonSolver.nIFaceParFlag[i] != true){
			comb_ngh[i] = basePar.neighbor[i];
		}

		// Assign the RHS
		RHS[i] = - basePar.vorticity[i];
	}

	#pragma omp parallel for
	for (int j = 0; j < poissonSolver.patchPar.num; j++){
		// Adjust the index
		int i = basePar.num + j;

		// Assign the coordinate
		comb_coord[0][i] = poissonSolver.patchPar.x[j];
		comb_coord[1][i] = poissonSolver.patchPar.y[j];

		// Assign the size
		comb_size[i] = poissonSolver.patchPar.s[j];

		// Assign the neighbor list
		comb_ngh[i] = _tempNeighbor[j];

		// Assign the RHS
		RHS[i] = - poissonSolver.patchPar.vorticity[j];
	}

	// Update the neighbor list
	if ((basePar.isAdapt == true) || 
		(Pars::opt_start_state == 0 && step == 0) ||
		(Pars::opt_start_state == 1 && step == Pars::resume_step+1))
	{
		for (const int &_nodeID : poissonSolver.nIFaceNodeID){
			// Alias to the current node
			Node *const &currNode = baseGrd.nodeMap.at(_nodeID);

			// Find the neighbor list
			std::vector<int> nghNodeList;
			baseGrd.findNghLvl(nghNodeList, currNode);

			// New neighbor list
			std::vector<int> newNodeList;
			for (const int &nghID : nghNodeList){
				if (poissonSolver.nIFaceNodeParMap.count(nghID) != 0){
					newNodeList.push_back(nghID);
				}
			}

			// Evaluate all particle inside the current node
			for (const int &parID : currNode->parList){
				// At current particle, take the neighbor with same size
				
				// Alias to current particle size
				const double &currSize = basePar.s[parID];

				// Iterate each neighbor on current particle
				for (const int &nghID : basePar.neighbor[parID]){
					// Alias to neighbor particle size
					const double &nghSize = basePar.s[nghID];

					// Compare the size between particle
					double ratio = nghSize / currSize;
					if (ratio > 0.95 && ratio < 1.05){
						// The current neighbor have same size
						comb_ngh[parID].push_back(nghID);
					}
				}

				// Then evaluate with all particle in new list
				for (const auto &nghNode : newNodeList){
					for (const auto &nghID : poissonSolver.nIFaceNodeParMap.at(nghNode)){
						// Calculate the distance square between two points
						double _dr2 = 0;
						// Calculate in x direction
							double _dx = basePar.x[parID] - poissonSolver.patchPar.x[nghID];
							_dr2 += (_dx*_dx);
						// Calculate in y direction
							double _dy = basePar.y[parID] - poissonSolver.patchPar.y[nghID];
							_dr2 += (_dy*_dy);

						// Evaluate distance
						if (Pars::opt_ngh_interact == 1){
							// Check the particle i
							double _rSup = basePar.s[parID] * Pars::r_sup * Pars::r_buff;   // Support radius of particle i
							if (_dr2 < (_rSup*_rSup)){
								comb_ngh[parID].push_back(basePar.num + nghID);
							}
						}
						else if (Pars::opt_ngh_interact == 2){
							// Support radius of average size
							double _rSupAve = ((basePar.s[parID] + poissonSolver.patchPar.s[nghID]) / 2.0e0) * Pars::r_sup * Pars::r_buff;
							if (_dr2 < (_rSupAve*_rSupAve)){
								comb_ngh[parID].push_back(basePar.num + nghID);
							}
						}
					}
				}
			}
		}
	}

	MESSAGE_LOG << "Number of combined particle : " << totalNum << "\n";

	// // // // [DEBUGGING TOOLS]
	// // // // =============================== START ===============================
	// // // save_data save_tool;
	// Particle newParticle;
	
	// newParticle.num = totalNum;
	// newParticle.x = comb_coord[0];
	// newParticle.y = comb_coord[1];
	// newParticle.s = comb_size;
	// newParticle.vorticity = RHS;
	// newParticle.neighbor = comb_ngh;

	// newParticle.chi.resize(totalNum, 0.0e0);
	// newParticle.u.resize(totalNum, 0.0e0);
	// newParticle.v.resize(totalNum, 0.0e0);
	// newParticle.gz.resize(totalNum, 0.0e0);
	// newParticle.isActive.resize(totalNum, false);
	// save_tool.save_par_state(newParticle, "Updated", 0);

	// int NN = 50;
    // int timer[NN];
    // for (int i = 0; i < NN; i++){
    //     timer[i] = (((((i+123) * 13) %143) + ((i+97) * 17))) * i;
    // }

    // WARNING_LOG << "Data Number : " << totalNum << "\n";
    // clock_t time = clock();
    // for (const auto &mul : timer){
    //     // int ID = (static_cast<int>(mul * (static_cast<double>(clock()-time)))) % totalNum;
	// 	int ID = mul % totalNum;
    //     MESSAGE_LOG << "List of random ID " << std::abs(ID) << "\n";
    //     save_ngh_all(newParticle, std::abs(ID));
    // }

	// throw std::exception();
	// // // // ================================ END ================================

	
	// PROCEDURE 2: Generate the global matrix
    // ********************************************************************
	// Generate the global matrix 
	//  * At initial of simulation step 0
	//  * At initial of resume step (resume_step + 1)
	//  * When adaptation is done
	if ((basePar.isAdapt == true) || 
		(Pars::opt_start_state == 0 && step == 0) ||
		(Pars::opt_start_state == 1 && step == Pars::resume_step+1))
	{
		#if (TIMER_PAR == 0)
			// Timer using super clock (chrono)
			std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
		#elif (TIMER_PAR == 1)
			// Timer using paralel package
			double _time = omp_get_wtime();
		#endif

		// Generate the global matrix
		printf("Generate the global matrix ... \n");
		poissonSolver.create_global_matrix(comb_coord, comb_ngh, comb_size, comb_bnd_loc /*comb_bnd_flag*/);

		#if (TIMER_PAR == 0)
			// Timer using super clock (chrono)
			std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
			double _time = span.count();
		#elif (TIMER_PAR == 1)
			// Timer using paralel package
			_time = omp_get_wtime() - _time;
		#endif
		// printf("<-> Global matrix computation time  :  [%f s]\n", _time);
		printf("<+> Initialize global matrix     :  %f s\n", _time);
	}


	// PROCEDURE 3: Solve the poisson equation
    // ********************************************************************
	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		double _time = omp_get_wtime();
	#endif

	// Solve the global matrix
	printf("\nSolve the potential of poisson equation ... \n");
	poissonSolver.solve(RHS);
	
	// Get the stream function calculated from poisson equation
	psi = poissonSolver.get_phi();

	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
		double _time = span.count();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime() - _time;
	#endif
	// printf("<-> Poisson solver computation time :  [%f s]\n", _time);
	printf("<+> Step 1: Solve global matrix  :  %f s\n", _time);
	
	
	// PROCEDURE 4: Solve the poisson equation
    // ********************************************************************
	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		tick = std::chrono::system_clock::now();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime();
	#endif
	
	// Calculate the differential of stream function
	printf("Calculate the rotational part of velocity ... \n");
	LSMPSa LSMPSpsi;
	psi.resize(basePar.num);
	LSMPSpsi.set_LSMPS(basePar.x, basePar.y, basePar.s, psi, basePar.neighbor);
	std::vector<double> dPsidx = LSMPSpsi.get_ddx();
	std::vector<double> dPsidy = LSMPSpsi.get_ddy();

	// Update the velocity
	for (int i = 0; i < basePar.num; i++){
		basePar.u[i] = dPsidy[i];
		basePar.v[i] = -dPsidx[i];
		// p.v[i] = psi[i];
	}

	#if (TIMER_PAR == 0)
		// Timer using super clock (chrono)
		span = std::chrono::system_clock::now() - tick;
		_time = span.count();
	#elif (TIMER_PAR == 1)
		// Timer using paralel package
		_time = omp_get_wtime() - _time;
	#endif
	// printf("<-> Velocity calculation comp. time :  [%f s]\n", _time);
	printf("<+> Step 2: Compute velocity     :  %f s\n", _time);

	return;
}

// #pragma endregion

// =====================================================
// +------------- VELOCITY CALCULATOR 3D --------------+
// =====================================================
// #pragma region VELOCITY_3D

/**
 *  @brief FMM based velocity calculation.
 *  NOTE: Works well only for 2D simulation.
 *  
 *  @param	_particle  The particle for velocity evaluation.
 *  @param  _step  The current simulation iteration.
 */
void VelocityCalc::velocity_fmm_3d(Particle &p, const int step){
	/* Procedure:
       1. Prepare the data for velocity calculation (by FMM)
       2. Create the tree cell data
       3. Calculate the velocity by FMM
       4. Update the particle data
    */

   	// The marker to calculate near body particle only or all particle
	bool calcAllParticle = ALL_PARTICLE;
	if (step % Pars::step_inv == 0){
		calcAllParticle = true;
	}

    // PROCEDURE 1: Data initialization
    // ********************************************************************
	// Internal data container for FMM calculation
	std::vector<std::vector<double>> particle_POS;	// Particle position
	std::vector<double> particle_SRC_x;				// The particle source value container
	std::vector<double> particle_SRC_y;				// The particle source value container
	std::vector<double> particle_SRC_z;				// The particle source value container

	/** INFO FOR FMM SELECTIVE CALCULATION
	 *  > There are 3 mark
	 *    - Mark 1	: Active particle
	 *    - Mark 2  : Evaluate particle (target U Active)
	 *    - Mark 3  : Target particle
	*/
	std::vector<int> eval_index;					// A translation particle ID from evaluated to original ID
	std::vector<bool> eval_mark(p.num, false);		// The mark for evaluated particle
	std::vector<bool> active_mark;				    // The particle active mark container
	std::vector<bool> target_mark;					// Mark for target particle to be calculated
	int evalNum = 0;								// The evaluated particle count

	if (calcAllParticle == true){
		// Type 1 calculation container initialization (evaluate all particle)
		particle_POS = std::vector<std::vector<double>> (p.num, std::vector<double>(DIM,0.0)); // Particle position
		particle_SRC_x = std::vector<double> (p.num, 0.0);	// The particle source value container
		particle_SRC_y = std::vector<double> (p.num, 0.0);	// The particle source value container
		particle_SRC_z = std::vector<double> (p.num, 0.0);	// The particle source value container
		active_mark = std::vector<bool> (p.num, false);	// The particle active mark container

		// Assign the particle data into container
		for (int i = 0; i < p.num; i++){
			// Particle position
			particle_POS[i][0] = p.x[i];
			particle_POS[i][1] = p.y[i];
			particle_POS[i][2] = p.z[i];

			// Particle source
			particle_SRC_x[i] = p.vortx[i] * std::pow(p.s[i], 3.0) / (4 * M_PI);  // Vortex strength in x direction
			particle_SRC_y[i] = p.vorty[i] * std::pow(p.s[i], 3.0) / (4 * M_PI);  // Vortex strength in y direction
			particle_SRC_z[i] = p.vortz[i] * std::pow(p.s[i], 3.0) / (4 * M_PI);  // Vortex strength in z direction
			
			// Particle mark
			active_mark[i] = p.isActive[i];
		}	
	}
	else{
		// Type 2 calculation container initialization (only evaluate near body partilce)
		// <!> Only take the evaluation particle
		for (int i = 0; i < p.num; i++){
			// Check if the particle an active particle
			if (p.isActive[i] == true || p.bodyPart[i] != -1){
				// Put the particle inside the evaluation list
				eval_index.push_back(i);
				eval_mark[i] = true;
			}
		}

		// Resize the marking container
		evalNum = eval_index.size();
		particle_POS = std::vector<std::vector<double>> (evalNum, std::vector<double>(DIM,0.0)); // Particle position
		particle_SRC_x = std::vector<double> (evalNum, 0.0);	// The particle source value container
		particle_SRC_y = std::vector<double> (evalNum, 0.0);	// The particle source value container
		particle_SRC_z = std::vector<double> (evalNum, 0.0);	// The particle source value container
		target_mark = std::vector<bool>(evalNum, false);
		active_mark = std::vector<bool>(evalNum, false);;
		
		#pragma omp parallel for
		for (int ID = 0; ID < evalNum; ID++){
			// Create alias to the original particle
			const int &ori_ID = eval_index[ID];

			// Put the index into the list 
			// Particle position
			particle_POS[ID][0] = p.x[ori_ID];
			particle_POS[ID][1] = p.y[ori_ID];
			particle_POS[ID][2] = p.z[ori_ID];

			// Particle source
			particle_SRC_x[ID] = p.vortx[ori_ID] * std::pow(p.s[ori_ID], 3.0) / (4 * M_PI);  // Vortex strength in x direction
			particle_SRC_y[ID] = p.vorty[ori_ID] * std::pow(p.s[ori_ID], 3.0) / (4 * M_PI);  // Vortex strength in y direction
			particle_SRC_z[ID] = p.vortz[ori_ID] * std::pow(p.s[ori_ID], 3.0) / (4 * M_PI);  // Vortex strength in z direction
			
			// Particle mark
			active_mark[ID] = p.isActive[ori_ID];
			target_mark[ID] = (p.bodyPart[ori_ID] != -1) ? true : false;
		}
	}
	
	
	// PROCEDURE 2: Create tree cell data
    // ********************************************************************
	// Initialize the tree
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif
	
	// Initialization Tree Map
	if (step == 0){
		printf("<+> Initialize tree ... \n");
		this->treeData.initializeTree(particle_POS, active_mark);
	}else{
		printf("<+> Update tree ... \n");
		this->treeData.updateTree(particle_POS, active_mark);
	}
	
	// [DEBUGGING LINE]
    // this->treeData.saveLeafTree(treeData, std::to_string(step));
    // this->treeData.saveTree(treeData, std::to_string(step));

	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<+> Tree finished in       : %f s\n", _time);
	
	
	// PROCEDURE 3: Calculate the velocity (by FMM)
    // ********************************************************************
	// The FMM operator class
	fmm3D FMM_tool;

	// Calculate the x direction stream function (psi) by FMM accelerated method
	printf("\nCalculate the 3D velocity FMM ... \n");
	if (calcAllParticle == true){
		FMM_tool.calcVelocity(this->treeData, particle_POS, active_mark, particle_SRC_x, particle_SRC_y, particle_SRC_z);
	}else{
		FMM_tool.calcVelocityNearBody(this->treeData, particle_POS, active_mark, target_mark, particle_SRC_x, particle_SRC_y, particle_SRC_z);
	}
	std::vector<double> xVelocity = FMM_tool.get_Field_x();
	std::vector<double> yVelocity = FMM_tool.get_Field_y();
	std::vector<double> zVelocity = FMM_tool.get_Field_z();
	
	// // Calculate the x direction stream function (psi) by FMM accelerated method
	// printf("\nCalculate stream function in x direction ... \n");
	// FMM_tool.calcField(this->treeData, particle_POS, active_mark, particle_SRC_x);
	// std::vector<double> dPSIx_dy = FMM_tool.get_Field_y();
	// std::vector<double> dPSIx_dz = FMM_tool.get_Field_z();

	// // Calculate the x direction stream function (psi) by FMM accelerated method
	// printf("\nCalculate stream function in y direction ... \n");
	// FMM_tool.calcField(this->treeData, particle_POS, active_mark, particle_SRC_y);
	// std::vector<double> dPSIy_dx = FMM_tool.get_Field_x();
	// std::vector<double> dPSIy_dz = FMM_tool.get_Field_z();

	// // Calculate the x direction stream function (psi) by FMM accelerated method
	// printf("\nCalculate stream function in z direction ... \n");
	// FMM_tool.calcField(this->treeData, particle_POS, active_mark, particle_SRC_z);
	// std::vector<double> dPSIz_dx = FMM_tool.get_Field_x();
	// std::vector<double> dPSIz_dy = FMM_tool.get_Field_y();

	
	// PROCEDURE 4: Update the particle velocity
    // ********************************************************************
	// Update the velocity (curl of stream function PSI)
	if (calcAllParticle == true){
		#pragma omp parallel for
		for (int i = 0; i < p.num; i++){
			// Grouping calculation
			p.u[i] = xVelocity[i];
			p.v[i] = yVelocity[i];
			p.w[i] = zVelocity[i];

			// // Single calculation
			// p.u[i] = dPSIz_dy[i] - dPSIy_dz[i];
			// p.v[i] = dPSIx_dz[i] - dPSIz_dx[i];
			// p.w[i] = dPSIy_dx[i] - dPSIx_dy[i];
		}
	}
	else{
		#pragma omp parallel for
		for (int ID = 0; ID < evalNum; ID++){
			// Create alias to the original particle
			const int &ori_ID = eval_index[ID];

			// Put the calculated velocity into the original container
			p.u[ori_ID] = xVelocity[ID];
			p.v[ori_ID] = yVelocity[ID];
			p.w[ori_ID] = zVelocity[ID];
		}
	}
	
	return;
}

/**
 *  @brief FMM based velocity calculation.
 *  NOTE: Works well only for 2D simulation.
 *  
 *  @param	_particle  The particle for velocity evaluation.
 *  @param  _step  The current simulation iteration.
 */
void VelocityCalc::velocity_fmm_3d_fast(Particle &p, const int step){
	/* Procedure:
       1. Prepare the data for velocity calculation (by FMM)
       2. Create the tree cell data
       3. Calculate the velocity by FMM
       4. Update the particle data
    */

   	// The marker to calculate near body particle only or all particle [Actually is not yet implemented]
	bool calcAllParticle = ALL_PARTICLE;
	// if (step % Pars::step_inv == 0){
	// 	calcAllParticle = true;
	// }

    // PROCEDURE 1: Data initialization
    // ********************************************************************
	// Internal data container for FMM calculation
	std::vector<std::vector<double>> particle_POS;	// Particle position
	std::vector<double> particle_SRC_x;				// The particle source value container
	std::vector<double> particle_SRC_y;				// The particle source value container
	std::vector<double> particle_SRC_z;				// The particle source value container

	/** INFO FOR FMM SELECTIVE CALCULATION
	 *  > There are 3 mark
	 *    - Mark 1	: Active particle
	 *    - Mark 2  : Evaluate particle (target U Active)
	 *    - Mark 3  : Target particle
	*/
	std::vector<int> eval_index;					// A translation particle ID from evaluated to original ID
	std::vector<bool> eval_mark(p.num, false);		// The mark for evaluated particle
	std::vector<bool> active_mark;				    // The particle active mark container
	std::vector<bool> target_mark;					// Mark for target particle to be calculated
	int evalNum = 0;								// The evaluated particle count

	if (calcAllParticle == true){
		// Type 1 calculation container initialization (evaluate all particle)
		particle_POS = std::vector<std::vector<double>> (p.num, std::vector<double>(DIM,0.0)); // Particle position
		particle_SRC_x = std::vector<double> (p.num, 0.0);	// The particle source value container
		particle_SRC_y = std::vector<double> (p.num, 0.0);	// The particle source value container
		particle_SRC_z = std::vector<double> (p.num, 0.0);	// The particle source value container
		active_mark = std::vector<bool> (p.num, false);	// The particle active mark container

		// Assign the particle data into container
		for (int i = 0; i < p.num; i++){
			// Particle position
			particle_POS[i][0] = p.x[i];
			particle_POS[i][1] = p.y[i];
			particle_POS[i][2] = p.z[i];

			// Particle source
			particle_SRC_x[i] = p.vortx[i] * std::pow(p.s[i], 3.0) / (4 * M_PI);  // Vortex strength in x direction
			particle_SRC_y[i] = p.vorty[i] * std::pow(p.s[i], 3.0) / (4 * M_PI);  // Vortex strength in y direction
			particle_SRC_z[i] = p.vortz[i] * std::pow(p.s[i], 3.0) / (4 * M_PI);  // Vortex strength in z direction
			
			// Particle mark
			active_mark[i] = p.isActive[i];
		}	
	}
	else{
		// Type 2 calculation container initialization (only evaluate near body partilce)
		// <!> Only take the evaluation particle
		for (int i = 0; i < p.num; i++){
			// Check if the particle an active particle
			if (p.isActive[i] == true || p.bodyPart[i] != -1){
				// Put the particle inside the evaluation list
				eval_index.push_back(i);
				eval_mark[i] = true;
			}
		}

		// Resize the marking container
		evalNum = eval_index.size();
		particle_POS = std::vector<std::vector<double>> (evalNum, std::vector<double>(DIM,0.0)); // Particle position
		particle_SRC_x = std::vector<double> (evalNum, 0.0);	// The particle source value container
		particle_SRC_y = std::vector<double> (evalNum, 0.0);	// The particle source value container
		particle_SRC_z = std::vector<double> (evalNum, 0.0);	// The particle source value container
		target_mark = std::vector<bool>(evalNum, false);
		active_mark = std::vector<bool>(evalNum, false);;
		
		#pragma omp parallel for
		for (int ID = 0; ID < evalNum; ID++){
			// Create alias to the original particle
			const int &ori_ID = eval_index[ID];

			// Put the index into the list 
			// Particle position
			particle_POS[ID][0] = p.x[ori_ID];
			particle_POS[ID][1] = p.y[ori_ID];
			particle_POS[ID][2] = p.z[ori_ID];

			// Particle source
			particle_SRC_x[ID] = p.vortx[ori_ID] * std::pow(p.s[ori_ID], 3.0) / (4 * M_PI);  // Vortex strength in x direction
			particle_SRC_y[ID] = p.vorty[ori_ID] * std::pow(p.s[ori_ID], 3.0) / (4 * M_PI);  // Vortex strength in y direction
			particle_SRC_z[ID] = p.vortz[ori_ID] * std::pow(p.s[ori_ID], 3.0) / (4 * M_PI);  // Vortex strength in z direction
			
			// Particle mark
			active_mark[ID] = p.isActive[ori_ID];
			target_mark[ID] = (p.bodyPart[ori_ID] != -1) ? true : false;
		}
	}
	
	
	// PROCEDURE 2: Calculate the velocity (by FMM)
    // ********************************************************************
	// The FMM operator class
	fmm3D FMM_tool;

	// Calculate the x direction stream function (psi) by FMM accelerated method
	printf("\nCalculate the 3D velocity FMM ... \n");
	FMM_tool.calcVelocityFast(particle_POS, active_mark, target_mark, particle_SRC_x, particle_SRC_y, particle_SRC_z);
	std::vector<double> xVelocity = FMM_tool.get_Field_x();
	std::vector<double> yVelocity = FMM_tool.get_Field_y();
	std::vector<double> zVelocity = FMM_tool.get_Field_z();
	
	// PROCEDURE 4: Update the particle velocity
    // ********************************************************************
	// Update the velocity (curl of stream function PSI)
	if (calcAllParticle == true){
		#pragma omp parallel for
		for (int i = 0; i < p.num; i++){
			// Grouping calculation
			p.u[i] = xVelocity[i];
			p.v[i] = yVelocity[i];
			p.w[i] = zVelocity[i];
		}
	}
	else{
		#pragma omp parallel for
		for (int ID = 0; ID < evalNum; ID++){
			// Create alias to the original particle
			const int &ori_ID = eval_index[ID];

			// Put the calculated velocity into the original container
			p.u[ori_ID] = xVelocity[ID];
			p.v[ori_ID] = yVelocity[ID];
			p.w[ori_ID] = zVelocity[ID];
		}
	}
	
	return;
}

// #pragma endregion