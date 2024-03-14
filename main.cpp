/* ##########################################
 *  VORTEX PARTICLE METHOD on LSMPS disc.
 *  Sec		: main.cpp
 *
 *  Created by Angga'18 on 16/10/23. -> It takes 2 months but still need more refinement :(
 *  Co.Auth	: Adhika'15, Ical'17
 *	
 *	FOR ANY KIND OF COPY OR REDUPLICATION
 *	DO NOT DELETE THIS LABEL !!!
 * ##########################################*/

/* ---------------- PROGRAM DESCRIPTION ----------------
    > The main program of VPM-LSMPS fluid solver program
    > All parameter used in running the program is stored
      in global.cpp
    > Detailed calculation process is given later on
    > This program contains:
      - FMM accelerated biot-savart using adaptive tree cell
      - LSMPS properties interpolation
      - Brinkmann penalization method
      - Adaptive Multiresolution Particle distribution
*/

// #pragma region INCLUDE_PACKAGE
	// Parameter and Utilities
	// ***********************
	#include "global.hpp"						// Simulation setting and parameter
	#include "Utils.hpp"						// Simulation data storage
	#include "src/grid_block/gridNode.hpp"		// Particle data grouping

	// Subroutine Packages
	// *******************
	#include "src/geometry/geometry.hpp"				// Body generation
	#include "src/initialization/initialization.hpp"	// Simulation domain initialization
	#include "src/remeshing/remeshing.hpp"				// A redistribution package
	#include "src/adaptation/adaptation.hpp"			// An adaptive particle distribution package
	#include "src/save_data/save_data.hpp"			    // Data write and saving package

	// Physics Calculation Method
	// **************************
	#include "src/advection/advection.hpp"
	#include "src/diffusion/diffusion.hpp"
	#include "src/stretching/stretching.hpp"
	#include "src/velocity_calculation/velocity_calc.hpp"
	#include "src/pressure_poisson/pressure_poisson.hpp"
	#include "src/penalization/penalization.hpp"
	#include "src/force_calculation/force_calc.hpp"
// #pragma endregion

#if (DATA_INTERPOLATION == 0)
// =====================================================
// +--------------- Program Starts Here ---------------+
// =====================================================
int main(int argc, char const *argv[])
{
	// ==================== Initial Definition Region ====================
	// *******************************************************************
	// #pragma region DATA_STORAGE
		std::vector<Body> bodyList(N_BODY);  // Obstacle solid object data list 
		Particle particle;                   // Simulation physical data storage
		GridNode nodeGridList;               // Node data structure (Particle Related Data)
	// #pragma endregion

	// #pragma region SUBROUTINE_INSTANCES
		// *Initialization tools
		geometry geom_tool;                  // Geometry generation
		initialization initialization_tool;  // Particle distribution
		remeshing remesh_tool;               // Particle redistribution and neighbor search
		
		// *Solver tools
		penalization penalization_tool;      // Penalization calculation
		VelocityCalc velocity_tool;      	 // Biot Savart velocity solver [accelerated by FMM]
		advection advection_tool;            // Advection calculation
		diffusion diffusion_tool;            // Diffusion calculation
		stretching stretching_tool;          // Stretching calculation
		pressure_poisson pressure_tool;	     // Poisson solver of pressure
		force_calculation force_tool;        // Calculate and save force data
		
		// *Simulation tools
		simUtil utility;                     // Simulation utilities
		save_data save_manager;              // Data writing manager
		std::ofstream _write;                // Data writer
	// #pragma endregion

	// #pragma region INTERNAL_VARIABLE
		double curr_comp_time = 0.0e0;       // Total computational time at each iteration
		double cum_comp_time  = 0.0e0;       // Cumulative of computational time
		int nt_start = 0;                    // The starting iteration step
			if      (Pars::opt_start_state == 0) nt_start = 0;
			else if (Pars::opt_start_state == 1) nt_start = Pars::resume_step + 1;
	// #pragma endregion

	// #pragma region SIMULATION_STABILITY_INDICATOR
		double CourMax = Pars::Courant;     // Global maximum courant number throughout simulation
		double DiffMax = Pars::Diffusion;   // Global maximum diffusion number throughout simulation
		double VortMax = Pars::Vortex;   	// Global maximum vorticity criteria throughout simulation (mesh criteria)
		
		// Collection of stability number
		std::vector<double*> stabNumList;
		stabNumList.push_back(&CourMax);
		stabNumList.push_back(&DiffMax);
		stabNumList.push_back(&VortMax);
	// #pragma endregion

	// [CONSOLE LOG] Display simulation parameter to console
	save_manager.summary_log();

	// [DATA  WRITE] Save the simulation parameter to file
	save_manager.write_summary_log();
	
	// Pre-processing prompt display
	printf("%s#=================================================#%s\n", FONT_RED, FONT_RED);
	printf("+---------------- %sPRE PROCESSING%s -----------------+\n", FONT_TORQUOISE, FONT_RED);
	printf("%s#=================================================#%s\n", FONT_RED, FONT_RESET);
	

	// =========================== Body Region ===========================
	// *******************************************************************
	// Generate all body data
	geom_tool.generateBody(bodyList);
	
	// Save all body data
	for (int i = 0; i < N_BODY; i++)
	save_manager.save_body_state(bodyList[i], std::to_string(i+1), 0);

	
	// ========================= Particle Region =========================
	// *******************************************************************
	// Initial particle distribution generation
	initialization_tool.initialize_particle(particle, bodyList, nodeGridList);
	initialization_tool.initialize_vorticity(particle);

	save_manager.save_grid_node_state(nodeGridList, "before_adapt", 1);
	remesh_tool.get_remeshing(particle, nodeGridList, 0, bodyList);

	save_manager.save_grid_node_state(nodeGridList, "after_adapt", 1);
	exit(1);
	// <!> Need an update for additional boundary element
	// initialization_tool.initialize_laplace(particle);

	// Generate the particle neighbor ID data list
	remesh_tool.set_neighbor(particle, nodeGridList);

	// // [Additional Patch]: update the domain boundary check
	// initialization_tool.update_domain_boundary(particle);
	
	// // [Additional Patch]: Initial particle distribution generation
	// initialization_tool.initialize_vorticity(particle);

	// save_manager.save_par_state(particle, "AWAL", 0);
	
	// // Adaptation of vorticity only for no obstacle simulation
	// if (N_BODY == 0)
	// remesh_tool.get_remeshing(particle, nodeGridList, 0, bodyList);

	// ** Here a Debug or Testing line Start
		/* Put your code here ... */
		// initialization_tool.laplace(particle);
		// initialization_tool.test_0(particle);
		// initialization_tool.test_1(particle);
		// initialization_tool.test_2(particle);
		// initialization_tool.initialize_vorticity(particle);
		
		// // [Additional Patch]: update the domain boundary check
		// initialization_tool.update_domain_boundary(particle);

	// ** Here a Debug or Testing line End

	// Set the initial particle data saved file name
	std::string saveNameInit;
		if 		(Pars::opt_start_state == 0) saveNameInit = "initial";	// Start at initial time
		else if (Pars::opt_start_state == 1) saveNameInit = "resume";	// Resume simulation (from the given iteration)
	
	// Save the particle data
	save_manager.save_par_state(particle, saveNameInit, 0);

	// Save the grid node data
	save_manager.save_grid_node_state(nodeGridList, saveNameInit, 1);


	// ========================= Physical Region =========================
	// *******************************************************************
	
	/* Still no code lies here ... */
	
	// Initialization the chi data
	// penalization_tool.initialize();			// [Create a chi calculation] [Already done in initialization sub.]
	// ** Vibrational Parameter Definition 	[FURTHER WORK, see most bottom part, outside the main block]
	
	
	// ===================== Pre Processing Summary ======================
	// *******************************************************************
	// Calculate the number of body node
	int sumBodyNode = 0;
	for (int i = 0; i < N_BODY; i++) sumBodyNode += bodyList[i].n_node;
	
	// Display the summary to the console
	printf("\n+------------ Initialization Summary -------------+\n");
	if (sumBodyNode == 0)	printf("No body simulation <!>\n");
	else					printf("Count of total body node                : %8d \n", sumBodyNode);
	printf("Count of particle node                  : %8d \n", particle.num);
	printf("Total iteration number                  : %8d \n", Pars::max_iter);
	printf("+-------------------------------------------------+\n\n");
	
	// Update the summary data file
	_write.open(save_manager.get_log_directory(), std::ofstream::out | std::ofstream::app);
	_write << "Number of initialized particle     :" << w12 << wR << particle.num << "\n";
	_write.close();


	// ====================== Simulation Run Prompt ======================
	// *******************************************************************
	// Simulation command prompt
	std::cout << FONT_BLUE  << "<!> The initialization is completed!\n"
			  << FONT_RESET << "<!> Do you want to run the simulation? ("
			  << FONT_GREEN << "yes" << FONT_RESET << "/" 
			  << FONT_RED   << "no"  << FONT_RESET << ")\n" 
			  << "    Type here: ";
	
	// Prompt to continue run the simulation
	// Give input to the prompt
	bool _run = false; std::string _cmd; std::cin >> _cmd;
	// Set the simulation command
	if (_cmd == "yes" || _cmd == "Yes" || _cmd == "Y" || _cmd == "y") _run = true;
	else _run = false;

	// // // ==========================================================================================================
	// // [!] TESTING POISSON SOLVER using Test Function
	// if (_run == false) exit(1);
	
	// // [2] Velocity calculation by biot savart: solving Rotational Velocity & Stretching
	// // velocity_tool.LSMPS_poisson_2d(particle, 0);
	// // if (particle.P.empty()){
	// // 	particle.P = std::vector<double> (particle.num);
	// // }
	// velocity_tool.get_velocity(particle, nodeGridList, 0);
	
	// // [!] TODO 2: Saving Particle Data
	// // Saving particle data
	// std::string DataName = "PV_FMM_MR1";		// Write data file name
	// save_manager.save_par_state(particle, DataName, -1);	// Saving particle data

	// // Calculation of L2 Norm Error of velocity
	// double SSDu = 0.0;		// The sum of square different
	// double SASu = 0.0;		// The sum of analytixal square
	// double SSDv = 0.0;		// The sum of square different
	// double SASv = 0.0;		// The sum of analytixal square
	// for (int i = 0; i < particle.num; i++){
	// 	// double diff = particle.gz[i] - particle.chi[i];
	// 	// SSD += diff*diff;
	// 	// SAS += particle.chi[i]*particle.chi[i];

	// 	double diff = particle.u[i] - particle.gx[i];
	// 	SSDu += diff*diff;
	// 	SASu += particle.gx[i]*particle.gx[i];

	// 	diff = particle.v[i] - particle.gy[i];
	// 	SSDv += diff*diff;
	// 	SASv += particle.gy[i]*particle.gy[i];
	// }
	// double L2Nu = std::sqrt(SSDu/SASu);
	// double L2Nv = std::sqrt(SSDv/SASv);

	// // Update the summary data file
	// _write.open(save_manager.get_log_directory(), std::ofstream::out | std::ofstream::app);
	// _write.precision(10);
	// _write << "L2Norm velocity x                  : " << wR << L2Nu << "\n";
	// _write << "L2Norm velocity y                  : " << wR << L2Nv << "\n";
	// _write.close();

	// std::cout << "And the L2 norm error velocity x of : " << wR << L2Nu << "\n";
	// std::cout << "And the L2 norm error velocity y of : " << wR << L2Nv << "\n";
	
	// // // Break the simulation
	// // throw std::exception();
	exit(1);
	// // // ==========================================================================================================
	

	// =====================================================
	// =============== SIMULATION ITERATION ================
	// =====================================================
	if (_run == true)
	{
		// Solving header prompt
		printf("\n%s#=================================================#%s", FONT_RED, FONT_RED);
		printf("\n+-------------------- %sSOLVING%s --------------------+", FONT_TORQUOISE, FONT_RED);
		printf("\n%s#=================================================#%s", FONT_RED, FONT_RESET);

		// Particle data counter throughout simulation
		int minParticleNum = particle.num;		// Minimum number of particle throughout simulation
		int maxParticleNum = particle.num;		// Maximum number of particle throughout simulation
		
		// ========================= Iteration Loop ==========================
		// *******************************************************************
		for (int step = nt_start; step < Pars::max_iter; step++){
			// ********************* ITERATION INITIALs **********************
			// Print current iteration header
			utility.printHeader(step);

			// Printing the stability criteria: courant (C), diffusion (Phi), and mesh reynold (Re_h) number
			if (Pars::flag_disp_stability || Pars::flag_save_stability){
				utility.stabilityEval(particle, stabNumList, step);
			}

			// Solver computational time manager
			auto t_start = std::chrono::system_clock::now();	// Using chrono

			// ******************** PHYSICAL CALCULATION *********************
			{
				/*
				Sequence of solving computation:
				1. Particle Redistribution  [Structured particle distribution]
				2. Velocity calculation     [Structured particle distribution]
				   -> Calculating rotational velocity
				   -> Calculating total velocity by helmholtz decomposition
				   -> Saving particle step 		[The best stage to save particle]
				3. Penalization             [Structured particle distribution]
				   -> Saving particle force 	[Must be save after penalization]
				   -> Saving particle step 		[Alternative stage]
				4. Perform Advection        [Unstructured particle distribution]
				5.1 Perform Diffusion       [Unstructured particle distribution]
				   -> Saving particle step
				5.2 Perform Stretching      [Unstructured particle distribution]
				   -> Saving particle step
				
				Try to rearrange the solver seq. 3->4->5->1->2

				Note on force calculation:
					> Penalization force calculation type must be calculated right after the penalization take place
					> The other type of force calculation still on investigation
				*/
				
				// SOLVER SEQUENCE : 3 -> 4 -> 5 -> 1 -> 2

				// *************** ADDITIONAL STEP [1] ***************
				if (Pars::flag_peturbation == true){
					// Set the iteration when peturbation generated
					double ptbTime = 3.0;				// The time when peturbation occured (Defined manually)
					int ptbIter = ptbTime / Pars::dt;	// The iteration step when peturbation occured
					
					// Add the peturbation into the domain
					if (step == ptbIter){
						// std::cout << "Adding a small peturbation ...\n";
						printf("%sAdding a small peturbation ...%s\n", FONT_CYAN, FONT_RESET);
						utility.addVorPertubation(particle);
						// save_manager.save_par_state(particle, "peturbation", 0);	// Saving particle data	
					}
				}
				// ***************************************************

				// =============== SOLVER STEP [3] ===============
				// [3] Perform penalization using Brinkmann: Penalize the velocity in body domain
				#if (N_BODY > 0)
					penalization_tool.get_penalization(particle, bodyList, step);
				#endif

				// ================ Barrier Mark - Data Saving =================
				// [!] TODO 1: Saving Data Force 
				force_tool.force_calc(particle, bodyList, step, 2, "aftPen");		// Penalization mode
				// force_tool.force_calc(particle, bodyList, step, 3, "aftPenMom");	// Vorticity moment mode
				
				// // [!] TODO 2: Saving Particle Data
				// if ((step % Pars::save_inv == 0)){
				// 	// Saving particle data
				// 	std::string DataName = utility.saveName(step);		// Write data file name
				// 	// DataName = "aftPen_" + DataName;
				// 	save_manager.save_par_state(particle, DataName, 0);	// Saving particle data
				// }
				
				// =============== SOLVER STEP [4] ===============
				// [4] Convection or Advection Sub-step: Perform the particle advection
				advection_tool.main_advection(particle);      			// [!] later: do 2nd order scheme
				
				// =============== SOLVER STEP [5.1] ===============
				// [5.1] Diffusion Sub-step: Calculate the vorticity diffusion & time integration
				if (DIM == 2) diffusion_tool.main_diffusion(particle); 	// [!] later: do 2nd order scheme

				// =============== SOLVER STEP [5.2] ===============
				// [5.2] Stretching Sub-step: Calculate the vorticity stretching & time integration
				if (DIM == 3) stretching_tool.calc_diff_stretch(particle);

				// // [!] TODO 2: Saving Particle Data
				// if ((step % Pars::save_inv == 0)){
				// 	// Saving particle data
				// 	std::string DataName = utility.saveName(step);		// Write data file name
				// 	DataName = "aftStr_" + DataName;
				// 	save_manager.save_par_state(particle, DataName, 0);	// Saving particle data
				// }

				// =============== SOLVER STEP [1] ===============
				// [1] Particle redistribution: rearrange the particle distribution by interpolating vorticity
				if ((step % Pars::rmsh_inv) == 0){
					// Particle redistribution (*every given iteration step)
					remesh_tool.get_remeshing(particle, nodeGridList, step, bodyList);
				}
				// [Additional Patch]: update the domain boundary check
				initialization_tool.update_domain_boundary(particle);

				// // [!] TODO 2: Saving Particle Data
				// if ((step % Pars::save_inv == 0)){
				// 	// Saving particle data
				// 	std::string DataName = utility.saveName(step);		// Write data file name
				// 	DataName = "aftRmsh_" + DataName;
				// 	save_manager.save_par_state(particle, DataName, 0);	// Saving particle data
				// }

				// =============== SOLVER STEP [2] ===============
				// [2] Velocity calculation by biot savart: solving Rotational Velocity & Stretching
				velocity_tool.get_velocity(particle, nodeGridList, step);

				// // [!] TODO 2: Saving Particle Data
				// if ((step % Pars::save_inv == 0)){
				// 	// Saving particle data
				// 	std::string DataName = utility.saveName(step);		// Write data file name
				// 	DataName = "aftVel_" + DataName;
				// 	save_manager.save_par_state(particle, DataName, 0);	// Saving particle data
				// }

				// // [!] DEBUG: Saving Particle Data
				// std::string DataName = utility.saveName(step);		// Write data file name
				// // DataName = "test_mres_" + DataName;
				// DataName = "hor_rec";
				// save_manager.save_par_state(particle, DataName, -1);	// Saving particle data
				
				// // Break here
				// exit(1);

			}	// End of the solver calculation

			// *********************** POST PROCESSING ***********************
			/** NOTE: The post processing contents are put right after the penalization calculation */
			// // [1] Saving Data Force 
			// if (Pars::flag_pressure_calc == true){
			// 	// Calculate the pressure
			// 	pressure_tool.get_pressure(particle);
			// }
			// force_tool.force_calc(particle, bodyList, step, 2, "force");

			// [2] Save the particle data at given step interval
			if ((step % Pars::save_inv == 0)){
				// Saving particle data
				std::string DataName = utility.saveName(step);		// Write data file name
				save_manager.save_par_state(particle, DataName, 0);	// Saving particle data	
			}

			// [3] Saving Residual
			if (Pars::flag_save_residual == true)
			utility.saveResidual(particle, step);

			// [4] Update the particle count
			minParticleNum = std::min<int>(minParticleNum, particle.num);
			maxParticleNum = std::max<int>(maxParticleNum, particle.num);


			// ********************* COMPUTATIONAL TIME **********************
			// Calculate the accumulative computational time
			std::chrono::duration<double> calculation_time = (std::chrono::system_clock::now() - t_start);
			curr_comp_time = calculation_time.count();				// Current iteration computational time
			cum_comp_time += curr_comp_time;						// Accumulative computational time

			// Saving the simulation time for each iteration
			if (Pars::flag_save_sim_time){
				// Print Header
				if (step == 0 || ((Pars::opt_start_state == 1) && (step == Pars::resume_step+1))){
					_write.open("output/Simulation_Time.csv");
					_write << "iteration,sim_time,par_num,comp_time,cum_comp_time\n";
					_write.close();
				}
				
				// Print Data
				_write.open("output/Simulation_Time.csv", std::ofstream::out | std::ofstream::app);
				_write <<  "" << step
						<< "," << step * Pars::dt
						<< "," << particle.num
						<< "," << curr_comp_time
						<< "," << cum_comp_time
						<< "\n";
				_write.close();
				
			}

			// ********************** ITERATION SUMMARY **********************
			// Displaying the simulation time for each iteration
			printf("\n%s**************** Iteration Summary ****************%s", FONT_GREEN, FONT_RESET);
			printf("\n<!> Current iteration sim. time   : %12.2f s", step * Pars::dt);
			printf("\n<!> Current iteration comp. time  : %12.2f s", curr_comp_time);
			printf("\n<!> Cumulative computational time : %12.2f s\n", cum_comp_time);
			
			printf("\n<!> Particle count                : %12d", particle.num);
			printf("\n<!> Iteration to go               : %12d\n", Pars::max_iter - step);
			
			// Prediction time to finish the simulation
			if (Pars::flag_disp_pred_time == true){
				utility.predictCompTime(step, curr_comp_time);
			}
			
			// **End of the iteration
		}
		
		// HERE: Outside the iteration loop
		// printf("%s\n<!> The iteration has just finished successfully! %s\n", FONT_GREEN, FONT_RESET);


		// ======================== Final Data Saving ========================
		// *******************************************************************
		// Summary of maximum value throughout the simulation
		_write.open(save_manager.get_log_directory(), std::ofstream::out | std::ofstream::app);
		_write << "Maximum number of particle         :" << w12 << wR << particle.num << "\n";
		_write << std::fixed << std::setprecision(2)
		       << "Total computational time           :" << w12 << wR << cum_comp_time << " s\n";
		_write << std::fixed << std::setprecision(4)
		       << "Max. global Cour. number (C_max)   :" << w12 << wR << CourMax << "\n"
		       << "Max. global Diff. number (phi_max) :" << w12 << wR << DiffMax << "\n"
		       << "Max. global Vort. number (Re_h_max):" << w12 << wR << VortMax << "\n";
		_write.close();

		
		// ======================== Simulation Summary =======================
		// *******************************************************************
		time_t now = time(0);			// Present timing

		printf("\n%s+-------------- Simulation Summary ---------------+%s\n", FONT_BLUE, FONT_RESET);
		printf("Simulation end at     : %s", std::ctime(&now));
		printf("Minimum number of particle           : %9d \n", minParticleNum);
		printf("Maximum number of particle           : %9d \n", maxParticleNum);
		printf("Total iteration number               : %9d \n", Pars::max_iter);
		
		if (cum_comp_time < 60.0){
			// If simulation runs below than 1 minute
			printf("Total computational time             : %9.2f s \n", cum_comp_time);
		}else if (cum_comp_time/60.0 < 60.0){
			// If simulation runs below than 1 hour
			int time_m = cum_comp_time/60;
			double time_s = cum_comp_time - (time_m*60);
			printf("Total computational time             : %2dm %5.2f s \n", time_m, time_s);
		}else{
			// If simulation runs longer than 1 hour
			printf("Total computational time             : %9.2f h \n", cum_comp_time/3600.0);
		}

		if (Pars::opt_start_state == 0){
			printf("Average iteration computing time     : %9.2f s \n", cum_comp_time/Pars::max_iter);
		}
		if (Pars::opt_start_state == 1){
			printf("Average iteration computing time     : %9.2f s \n", cum_comp_time/(Pars::max_iter - Pars::resume_step));
		}

		printf("Max. global Courant number (C_max)   : %9.2f \n", CourMax);
		printf("Max. global Diff. number (phi_max)   : %9.2f \n", DiffMax);
		printf("Max. global Vort. number (Re_h_max)  : %9.2f \n", VortMax);
		printf("%s+-------------------------------------------------+%s\n", FONT_BLUE, FONT_RESET);

		// Simulation is done successfully!
		printf("%s<!> The simulation is finished successfully! %s\n", FONT_GREEN, FONT_RESET);
	}
	
	else{
		// Simulation is not executed!
		printf("%s<!> The simulation is not executed! %s\n", FONT_BLUE, FONT_RESET);
	}

	return 0;
}


#elif (DATA_INTERPOLATION == 1)
#include "src/grid_block/generateGrid.hpp"
// =====================================================
// +--------------- Interpolation Data ----------------+
// =====================================================
int main(int argc, char const *argv[])
{
	/** 
	 *  @brief	A main subroutine for particle data interpolation. This process belongs to 
	 *  post processing subroutines of multiresolution data in 3D space. In order to 
	 *  efficiently render the data into paraview, a structured single resolution data is 
	 *  necessary. This subroutine will process all saved particle state data from the
	 *  previously computed simulation. The interpolated data are labeled by "SRID".
	 * 
	 *  NOTE:
	 * 	 > This subroutine only belong to 3D space simulation, since 2D simulation is not necessary.
	 *   > Subroutine also provided with Q and λ2 calculated data.
	 *   > The collected data is 
	 * 		(1) Coordinate[3]
	 * 		(2) Velocity[3]
	 * 		(3) Vorticity[3]
	 * 		(4) Q[1]
	 * 		(5) λ2[1]
	 *     of total 11 variable set
	 * 
	 *  PREPARE:
	 *   > Set the dimension in 3D (Still need to be adjusted)
	 *   > Define the target domain size and interpolated particle spacing (see global.cpp)
	 *   > 
	*/

	// Main interpolation parameter (can be adjusted manually)
	// const int fin_iter = Pars::max_iter;            // The final iteration for interpolation
	// const int save_int = Pars::save_inv;			// Iteration interval for saving data
	// const int data_num = 1 + fin_iter/save_int;     // The number of particle state data (including the zeros)

	const int fin_iter = 8000;            // The final iteration for interpolation
	const int save_int = 20;			// Iteration interval for saving data
	const int data_num = 1 + fin_iter/save_int;     // The number of particle state data (including the zeros)

	// #pragma region SUBROUTINE_INSTANCES
		initialization initialization_tool;  // Particle distribution
		remeshing remesh_tool;               // Particle redistribution
		simUtil utility;                     // Simulation utilities
		save_data save_manager;              // Data writing manager
	// #pragma endregion

	// Pre-processing prompt display
	printf("%s#=================================================#%s\n", FONT_RED, FONT_RED);
	printf("+-------------------- %sSUMMARY%s --------------------+\n", FONT_TORQUOISE, FONT_RED);
	printf("%s#=================================================#%s\n", FONT_RED, FONT_RESET);

	int xCnt = std::ceil(Pars::lxdomInt/Pars::sigmaInt);
	int yCnt = std::ceil(Pars::lydomInt/Pars::sigmaInt);
	int zCnt = std::ceil(Pars::lzdomInt/Pars::sigmaInt);

	// Initial flow parameters summary data
	printf("\n+--------- Interpolation Parameters Data ---------+\n");
	printf(" Domain x length                       : %7.2f m\n", Pars::lxdomInt);
	printf(" Domain y length                       : %7.2f m\n", Pars::lydomInt);
	printf(" Domain z length                       : %7.2f m\n", Pars::lzdomInt);
	printf(" Interpolated core size                : %7.2f m\n", Pars::sigmaInt);
	printf(" Number of particle in x basis         : %9d m\n", xCnt);
	printf(" Number of particle in y basis         : %9d m\n", yCnt);
	printf(" Number of particle in z basis         : %9d m\n", zCnt);
	printf(" Total particle number                 : %9d m\n", xCnt*yCnt*zCnt);
	printf(" Final iteration step                  : %9d \n", fin_iter);
	printf(" Step interval                         : %9d \n", save_int);
	printf(" Number of sate data                   : %9d \n", data_num);
	printf("+-------------------------------------------------+\n\n");


	// ====================== Simulation Run Prompt ======================
	// *******************************************************************
	// Simulation command prompt
	std::cout << FONT_RESET << "<!> Proceed to the interpolation? ("
			  << FONT_GREEN << "yes" << FONT_RESET << "/" 
			  << FONT_RED   << "no"  << FONT_RESET << ")\n" 
			  << "    Type here: ";
	
	// Prompt to continue run the simulation
	// Give input to the prompt
	bool _run = false; std::string _cmd; std::cin >> _cmd;
	// Set the simulation command
	if (_cmd == "yes" || _cmd == "Yes" || _cmd == "Y" || _cmd == "y") _run = true;
	else _run = false;


	// ================ Iteration for data interpolation =================
	// *******************************************************************
	if (_run == true)
	{
		// Print banner
		printf("\n%s#=================================================#%s", FONT_RED, FONT_RED);
		printf("\n+----------------- %sINTERPOLATION%s -----------------+", FONT_TORQUOISE, FONT_RED);
		printf("\n%s#=================================================#%s", FONT_RED, FONT_RESET);
		
		// Computational time accumulation
		#if (TIMER_PAR == 0)
			// Timer using super clock (chrono)
			std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
		#elif (TIMER_PAR == 1)
			// Timer using paralel package
			double _time = omp_get_wtime();
		#endif
		
		
		// ======================== Interpolation Loop =======================
		// *******************************************************************
		// Start the interpolation to all particle data state
		for (int n = 0; n < data_num; n++){
			// Alias to the current iteration* [Can be adjusted manually]
			const int step = n * save_int;
			
			// Print header
			utility.printHeader(step);

			// Internal variable
			Particle srcParticle;      // Source particle data container (particle from local storage)
			Particle intParticle;      // Interpolated particle data container
			GridNode nodeGridList;     // Node data structure (Related to the source particle)

			// Read the source data
			#if (DIM == 2)
				initialization_tool.read_2d_particle(srcParticle, step);
			#elif (DIM == 3)
				initialization_tool.read_3d_particle(srcParticle, step);
			#endif

			// Create grid of source data
			printf("%s\nGenerate grid block ... %s\n", FONT_CYAN, FONT_RESET);
			generateGrid grid_step;
			grid_step.createNode(nodeGridList, srcParticle);
			
			// Perform interpolation
			remesh_tool.re_arrange_distribution(intParticle, srcParticle, nodeGridList);
			
			// Save the interpolated data
			std::string stepName = utility.saveName(step);
			save_manager.save_par_interpolation(intParticle, stepName);
		}

		// Interpolation summary!
		// Particle generation summary time display
		#if (TIMER_PAR == 0)
			// Timer using super clock (chrono)
			std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
			double _time = span.count();
		#elif (TIMER_PAR == 1)
			// Timer using paralel package
			_time = omp_get_wtime() - _time;
		#endif

		
		// ======================== Simulation Summary =======================
		// *******************************************************************
		time_t now = time(0);			// Present timing

		printf("\n%s+-------------- Simulation Summary ---------------+%s\n", FONT_BLUE, FONT_RESET);
		printf("Interpolation end at  : %s", std::ctime(&now));
		if (_time < 60.0){
			// If simulation runs below than 1 minute
			printf("Total computational time             : %9.2f s \n", _time);
		}else if (_time/60.0 < 60.0){
			// If simulation runs below than 1 hour
			int time_m = _time/60;
			double time_s = _time - (time_m*60);
			printf("Total computational time             : %2dm %5.2f s \n", time_m, time_s);
		}else{
			// If simulation runs longer than 1 hour
			printf("Total computational time             : %9.2f h \n", _time/3600.0);
		}
		printf("Average iteration computing time     : %9.2f s \n", _time/data_num);
		

		// Simulation is done successfully!
		printf("%s\n<!> The post-processing is finished successfully! %s\n", FONT_GREEN, FONT_RESET);
	}
	else{
		// Interpolation is not executed!
		printf("%s\n<!> The interpolation is cancelled! %s\n", FONT_BLUE, FONT_RESET);
	}
	

	return 0;
}
#endif

/* PROGRAMMING & SIMULATION NOTEs:
   +) Makefile directory -> c:\program files\gnuwin32\bin\make.exe
   +) To avoid memory leak, please do an initalization for all struct members (e.g, set to be 0)
   +) Plese put attention to variable type definition
   +) Remember to free the memory of malloc variable
*/

/*	List of task NEED TO BE DONE
	> Generate particle		[DONE]
	> Adaptation method		[DONE] Also done the linking, but not the cell list, (further a do? The answer is not)
	> Link the remeshing 	[DONE]
	> Redistribution 		[DONE] -> CHECK LSMPS (have been adjusted with paralel)
	> Neigbor evaluation	[DONE]
	> Set a paralel computing  [ALMOST DONE] -> There are something strange things by using parallel
		>>> Syntax	: #pragma omp parallel for
		- Too high number of cores used creates an unstable sequence, break ups the simulation data
		- At some point (a parallel loop deep inside the loop), the parallel make a sequence hold so it takes longer time than usual instead
	> Check all physical subroutine
		- Velocity ? (Check the vector data size)
		- Penalization (adjust to the number of body part)
		- Advection, Diffusion, just need to check...
*/

/*
	// DEBUGING PROCESS
	// Put the variable adaptive as 
	particle.vorticity.clear();
	particle.vorticity.resize(particle.num,1e10);
	for (Body &_body : bodyList){
		geom_tool.distance_calc(particle,_body,true);
		for (int i = 0; i < particle.num; i++){
			particle.vorticity[i] = std::min<double>(particle.R[i], particle.vorticity[i]);
		}
	}
	for (int i = 0; i < particle.num; i++){
		particle.vorticity[i] = std::pow(1.0 / (1 + particle.vorticity[i]), 2.0);
	}

	// for (int i = 0; i < particle.num; i++){
	// 	particle.vorticity[i] = 0.0;
	// }

	// particle.vorticity[60] = 1000.0;

	save_manager.save_par_state(particle,"Before_Adapted", 0);

	remesh_tool.get_remeshing(particle,nodeGridList,0);
	
	// adaptation adapt_step;
	// particle.neighbor.resize(particle.num);
	// adapt_step.get_adaptation(particle,&particle,nodeGridList);
	// adapt_step.get_new_particle(particle);
	// }

	particle.u.clear();
	particle.v.clear();
	particle.chi.clear();
	particle.gz.clear();
	particle.u.resize(particle.num);
	particle.v.resize(particle.num);
	particle.chi.resize(particle.num);
	particle.gz.resize(particle.num);

	particle.vorticity.clear();
	particle.vorticity.resize(particle.num,1e10);
	for (Body &_body : bodyList){
		geom_tool.distance_calc(particle,_body, true);
		for (int i = 0; i < particle.num; i++){
			particle.vorticity[i] = std::min<double>(particle.R[i], particle.vorticity[i]);
		}
	}
	for (int i = 0; i < particle.num; i++){
		particle.vorticity[i] = std::pow(1.0 / (1 + particle.vorticity[i]), 2.0);
		particle.s[i] = Pars::sigma * Pars::intPow(2,Pars::max_level-particle.level[i]);
		// particle.vorticity[i] = 1.0;
	}
*/

/* Vibrational Parameter Definition [FURTHER WORK]
	// ========== Vibrational Parameter Display ==========	
	if (Pars::vib == 1){
		
		Pars::mass = Pars::m_star * 0.5 * Pars::RHO * pow(Pars::Df,2) - Pars::m_d;
		Pars::SpringConst = Pars::k_star * 0.5 * Pars::RHO * pow(Pars::U_inf,2);
		Pars::DamperConst = Pars::c_star * 0.5 * Pars::RHO * Pars::Df * Pars:: U_inf; 

		printf("\n+-------------- Vibration Parameter --------------+\n");
		printf("Mass                                    : %8.4f \n", Pars::mass);
		printf("Spring Constant                         : %8.4f \n", Pars::SpringConst);
		printf("Damper Constant                         : %8.4f \n", Pars::DamperConst);
		printf("Body Area                               : %8.4f \n", Pars::area);
		printf("+-------------------------------------------------+\n");

	} else if (Pars::vib == 2){
		
		Pars::inertia = Pars::i_star * 0.5 * Pars::RHO * pow(Pars::Df,4); 
		Pars::SpringConst = pow(((Pars::u_inf / (Pars::U_star*Pars::Df)) * 2 * Pars::pi),2) * Pars::inertia;
		Pars::DamperConst = Pars::chi_star * 2 *sqrt(Pars::SpringConst*Pars::inertia); 
		Pars::tetha = Pars::alpha_a;
		
		printf("\n+-------------- Vibration Parameter --------------+\n");
		printf("Inertia                                 : %8.4f \n", Pars::inertia);
		printf("Spring Constant                         : %8.4f \n", Pars::SpringConst);
		printf("Damper Constant                         : %8.4f \n", Pars::DamperConst);
		printf("Angle of Attack                         : %8.4f \n", Pars::alpha_a);
		printf("+-------------------------------------------------+\n");
	}
*/