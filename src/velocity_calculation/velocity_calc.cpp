#include "velocity_calc.hpp"
#include "../FMM/fmm2D.hpp"
#include "../FMM/fmm3D.hpp"
#include "velocity_biot_savart.hpp"

/** The type of 2D simulation FMM calculation
 *   [1] : Old fashion code
 *   [2] : New code using tree code built in unordered map of cell pointer
*/
#define VELOCITY_2D_TYPE 2

/** The type of 3D simulation FMM calculation
 *   [1] : Refactor of old code (tree data built in vector, turns out to be 1.5 times faster)
 *   [2] : New code using tree code built in unordered map of cell pointer
*/
#define VELOCITY_3D_TYPE 2

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
	double _time;
	_time = omp_get_wtime();

	// **Calculate the rotational part
	if (DIM == 2){
		// Resize the velocity of particle
		p.u.clear(); p.u.resize(p.num,0.0e0);
		p.v.clear(); p.v.resize(p.num,0.0e0);

		// Old fashing FMM code of velocity calculation
		if (VELOCITY_2D_TYPE == 1){
			printf("%s<+> Type 1: Old based velocity calculation %s\n", FONT_CYAN, FONT_RESET);
			this->velocity_old(p, step);
		}

		// FMM based velocity calculation
		else if (VELOCITY_2D_TYPE == 2){
			printf("%s<+> Type 2: Unordered map tree cell pointer %s\n", FONT_CYAN, FONT_RESET);
			this->velocity_fmm_2d(p, step);
		}
	}

	else if (DIM == 3){
		// Resize the velocity of particle
		p.u.clear(); p.u.resize(p.num,0.0e0);
		p.v.clear(); p.v.resize(p.num,0.0e0);
		p.w.clear(); p.w.resize(p.num,0.0e0);

		// FMM based velocity calculation
		if (VELOCITY_3D_TYPE == 1){
			printf("%s<+> Type 1: Vector tree cell data %s\n", FONT_CYAN, FONT_RESET);
			this->velocity_fmm_3d_fast(p, step);
		}

		// FMM based velocity calculation
		else if (VELOCITY_3D_TYPE == 2){
			printf("%s<+> Type 2: Unordered map tree cell pointer %s\n", FONT_CYAN, FONT_RESET);
			this->velocity_fmm_3d(p, step);
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
	_time = omp_get_wtime() - _time;
	printf("<-> Velocity calculation comp. time:   [%f s]\n", _time);

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
	double start = omp_get_wtime();
	
	// Initialization Tree Map
	if (step == 0){
		printf("<+> Initialize tree ... \n");
		this->treeData.initializeTree(particle_POS, particle_mark);
	}else{
		printf("<+> Update tree ... \n");
		this->treeData.updateTree(particle_POS, particle_mark);
	}

	// [DEBUGGING LINE]
    // this->treeData.saveLeafTree(treeData, std::to_string(step));
    // this->treeData.saveTree(treeData, std::to_string(step));

	double finish = omp_get_wtime();
	printf("<+> Tree finished in : %f s\n", finish-start);
	
	
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
	double start = omp_get_wtime();
	
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

	double finish = omp_get_wtime();
	printf("<+> Tree finished in : %f s\n", finish-start);
	
	
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