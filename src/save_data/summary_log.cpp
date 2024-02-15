#include "save_data.hpp"

/**
 *  @brief Display the simulation parameter to console.
 *  
 *  @headerfile data_saving.hpp
*/
void save_data::summary_log(){
	std::cout << "<!> Computation Message Log\n";
	
	// Simulation Log Initial Parameter Summary Data
	
	// // PLAIN HEADER
	// std::cout << "#=================================================#\n";
	// std::cout << "+---------------- SIMULATION LOG -----------------+\n";
	// std::cout << "#=================================================#\n";

	printf("%s#=================================================#%s\n", FONT_RED, FONT_RED);
	printf("+---------------- %sSIMULATION LOG%s -----------------+\n", FONT_TORQUOISE, FONT_RED);
	printf("%s#=================================================#%s\n", FONT_RED, FONT_RESET);

	// Body data parameter
	std::vector<std::string> 
	bodyName2D = {"[Circular cylinder]",    // Type 1
				  "[Square cylinder]",      // Type 2
				  "[Normal plate]",         // Type 3
				  "[Flat plate]",           // Type 4
				  "[NACA-" +                // Type 5
				  	std::to_string(int(100*Pars::m_a)) +
				  	std::to_string(int( 10*Pars::p_a)) + 
				  	std::to_string(int(100*Pars::t_a)) + " airfoil]"
				 },            // 2D Body object name
	
	bodyName3D = {"[Sphere]",				// Type 1
				  "[Cube]",					// Type 2
				  "[3D normal plate]",		// Type 3
				  "[3D flat plate]",		// Type 4
				  "[Torus]",				// Type 5
				  "[Heart]"					// Type 6
				 };            // 3D Body object name
	
	// Simulation parameter summary data
	printf("+---------- Simulation Parameters Data -----------+\n");
	printf("Simulation dimension                  :     type %d \n", DIM);
	printf("Body option\n");
	for (int i = 0; i < N_BODY; i++){
		const int &part = Pars::opt_body[i];	// Aliasing of the body part
		std::string _bodyName;					// Aliasing of the body name
		if 		(DIM == 2) _bodyName = FONT_CYAN + bodyName2D[part-1] + FONT_RESET;
		else if (DIM == 3) _bodyName = FONT_CYAN + bodyName3D[part-1] + FONT_RESET;
		printf("  <+> Part %d%34s :     type %d\n", 1+i, _bodyName.c_str(), part);
	}
	printf("Initialization option                 :     type %d \n", Pars::opt_init_particle);
	printf("Neighbor search option                :     type %d \n", Pars::opt_neighbor);
	printf("Penalization option                   :     type %d \n", Pars::opt_pen);
	printf("Maximum resolution level              :   %d levels \n", Pars::max_level);
	printf("Core size                             : %8.4f m\n", Pars::sigma);
	printf("Time step                             : %8.4f s\n", Pars::dt);
	printf("Total simulation time                 : %8.2f s\n", Pars::sim_time);
	printf("Total iteration step                  : %8d\n", Pars::max_iter);
	printf("+-------------------------------------------------+\n");
	
	// Initial flow parameters summary data
	printf("\n+------------- Flow Parameters Data --------------+\n");
	printf("Reynolds number (RE)           : %10.2f [-]\n", Pars::RE);
	printf("Freestream velocity (U)        : %10.2f m/s\n", Pars::U_inf);
	printf("Fluid density (rho)            : %10.2f kg/m^3\n", Pars::RHO);
	printf("Fluid viscosity (nu)           : %10f m^2/s\n", Pars::NU);
	printf("+-------------------------------------------------+\n");

	// Additional calculation data
	printf("\n+--------------- Stability Initial ---------------+\n");
	printf("Courant number (C < 1.0)            : %12f \n", Pars::Courant);
	printf("Diffusion number (Phi < 0.5)        : %12f \n", Pars::Diffusion);
	printf("Stability Criteria (Re_h < 1.0)     : %12f \n", 100 * Pars::sigma*Pars::sigma/Pars::NU);
	printf("Turbulent scaling                   : %12f \n", std::ceil(1.0e0 / Pars::Courant));
	printf("Max time step (Re_h criteria)       : %12f \n", 0.25*Pars::sigma*Pars::sigma/Pars::NU);
	printf("+-------------------------------------------------+\n\n");
	return;
}

/**
 *  @brief Write the simulation parameter.
 *  
 *  @headerfile data_saving.hpp
*/
void save_data::write_summary_log(){
	// Cancel the saving procedure if flag is closed
	if (!Pars::flag_save_parameter){
		return;
	}
	
	// Body data parameter
	// Body data parameter
	std::vector<std::string> 
	bodyName2D = {"[Circular cylinder]",    // Type 1
				  "[Square cylinder]",      // Type 2
				  "[Normal plate]",         // Type 3
				  "[Flat plate]",           // Type 4
				  "[NACA-" +                // Type 5
				  	std::to_string(int(100*Pars::m_a)) +
				  	std::to_string(int( 10*Pars::p_a)) + 
				  	std::to_string(int(100*Pars::t_a)) + " airfoil]"
				 },            // 2D Body object name
	
	bodyName3D = {"[Sphere]",				// Type 1
				  "[Cube]",					// Type 2
				  "[3D normal plate]",		// Type 3
				  "[3D flat plate]",		// Type 4
				  "[Torus]",				// Type 5
				  "[Heart]"					// Type 6
				 };            // 3D Body object name
	
	// Log saving before simulation
	std::ofstream data;
	data.open(this->sumLogDir);
	data << "#=================================================#\n"
	     << "+---------------- SIMULATION LOG -----------------+\n"
		 << "#=================================================#\n\n";
	
	data << "+----------- Simulation Flow Parameter -----------+\n";
	data << std::fixed << std::setprecision(2)
	     << "Reynolds number (RE)             : " << w8 << wR << Pars::RE    << " [-]" << "\n"
		 << "Freestream velocity (U)          : " << w8 << wR << Pars::u_inf << " m/s"<< "\n"
 		 << "Fluid density (rho)              : " << w8 << wR << Pars::RHO   << " kg/m^3" << "\n"
		 << "Fluid viscosity (nu)             : " << w8 << wR << Pars::NU    << " m^2/s" << "\n"
		 << "+-------------------------------------------------+\n\n";

	data << "+----------- Simulation Setting Option -----------+\n"
		 << "Simulation dimension                    : " << DIM << "D space\n"
	     << "Body option: \n";
	for (int i = 0; i < N_BODY; i++){
		const int &part = Pars::opt_body[i];	// Aliasing of the body part
		std::string _bodyName;					// Aliasing of the body name
		if 		(DIM == 2) _bodyName = bodyName2D[part-1];
		else if (DIM == 3) _bodyName = bodyName3D[part-1];
		
		data << "  <+> Part " << 1+i << "  " << w25 << _bodyName << " : type " << part << "\n";
	}
	data << "Initialization option                   : " << "type " << Pars::opt_init_particle << "\n"
		 << "Neighbor search option                  : " << "type " << Pars::opt_neighbor << "\n"
		 << "Penalization option                     : " << "type " << Pars::opt_pen << "\n"
	     << "Penalization constant                   : " << w8 << wR << Pars::lambda << "\n";
	data << std::fixed << std::setprecision(4)
		 << "Number of saved data                    : " << w6 << wR << Pars::max_iter/Pars::save_inv << "\n"
		 << "Saving step interval                    : " << w6 << wR << Pars::save_inv << "\n"
		 << "Support radius factor                   : " << w6 << wR << Pars::r_sup << "\n"
		 << "Buffer radius factor                    : " << w6 << wR << Pars::r_buff << "\n"
		 << "Maximum resolution level                : " << Pars::max_level << " levels\n"
		 << "Core size                               : " << w6 << wR << Pars::sigma << " m\n"
		 << "Time step                               : " << w6 << wR << Pars::dt << " s\n";
	data << std::fixed << std::setprecision(2)
		 << "Total simulation time                   : " << w6 << wR << Pars::sim_time << " s\n"
		 << "Total iteration step                    : " << w6 << wR << Pars::max_iter << " \n"
		 << "+-------------------------------------------------+\n\n";

	data << "+----------- Simulation Parameter Data -----------+\n"
	     << "Domain x length                    :" << w12 << wR << Pars::lxdom << " m\n"
		 << "Domain y length                    :" << w12 << wR << Pars::lydom << " m\n";
	if (DIM == 3)
	data << "Domain z length                    :" << w12 << wR << Pars::lzdom << " m\n";
	data << "Upstream gap length                :" << w12 << wR << Pars::xdom  << " m\n"
		 << "Reference length (D)               :" << w12 << wR << Pars::Df    << " m\n";
	data << std::fixed << std::setprecision(4)
		 << "Courant number (C < 1.0)           :" << w12 << wR << Pars::Courant << "\n"
		 << "Diffusion number (Phi < 0.5)       :" << w12 << wR << Pars::Diffusion << "\n"
		 << "Stability Criteria (Re_h < 1.0)    :" << w12 << wR << 100 * Pars::sigma*Pars::sigma/Pars::NU << "\n"
		 << "Turbulent scaling                  :" << w12 << wR << std::ceil(1.0e0 / Pars::Courant) << "\n"
		 << "Maximum time step (Re_h criteria)  :" << w12 << wR << 0.25*Pars::sigma*Pars::sigma/Pars::NU << " s\n"
		 << "+-------------------------------------------------+\n\n";

	data << std::fixed << std::setprecision(-1);
	
	// End of writting simulation setting
	data.close();
	return;
}