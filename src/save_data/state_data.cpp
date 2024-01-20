#include "save_data.hpp"

/**
 *  @brief Write the body coordinate properties.
 *  
 *  @param	_body	The body object to be write.
 *  @param	_name	Identifier name for saved file.
 * 	@param	_type	Flag to write the body: [0] save both, [1] node only, [2] panel only.
 * 
 *  @headerfile save_data.hpp
*/
void save_data::save_body_state(const Body &b, std::string name, int type){
	// Cancel the body saving if the flag turn off
	if (Pars::flag_save_body == false) return;

	// Set up the write flag
	bool saveNode = true;
	bool savePanel = true;
	if (type == 1) savePanel = false;		// Only save the Node
	else if (type == 2) saveNode = false;	// Only save the Panel

	// Time counting
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

	// Internal object variable
	std::string name1;
	std::ofstream writeData;

	// Save the body node
	if (saveNode){
		// Message Log
		std::cout << FONT_GREEN << "Saving body node data part " << name << " ..." << FONT_RESET << "\n";
		name1 = "output/body_node_" + name + ".csv";
		writeData.open(name1.c_str());
		
		// Write the table header
		writeData << "node" 
			 	  << "," << "xb"
			 	  << "," << "yb";
		if (DIM > 2)
		writeData << "," << "zb";
		writeData << "\n";

		// Write the table data
		for (int i = 0; i < b.n_node; i++)
		{
			writeData << i
				 	  << "," << b.x[i]
				 	  << "," << b.y[i];
			if (DIM > 2){
			writeData << "," << b.z[i];}
			writeData << "\n";
		}
		writeData.close();
	}
	
	// Save the body panel data
	if (savePanel){
		// Message Log
		std::cout << FONT_GREEN << "Saving body panel data part " << name << " ..." << FONT_RESET << "\n";
		name1 = "output/body_panel_" + name + ".csv";
		writeData.open(name1.c_str());

		// Write the table header
		writeData << "panel" 
			 	  << "," << "size"
			 	  << "," << "x_mid"
			 	  << "," << "y_mid";
		if (DIM > 2)
		writeData << "," << "z_mid";
		
		writeData << "," << "x_normal"
			 	  << "," << "y_normal";
		if (DIM > 2)
		writeData << "," << "z_normal";
		writeData << "\n";

		// Write the table data
		for (int i = 0; i < b.n_panel; i++)
		{
			writeData << i
				 	  << "," << b.size[i]
				 	  << "," << b.x_m[i]
				 	  << "," << b.y_m[i];
			if (DIM > 2)
			writeData << "," << b.z_m[i];

			writeData << "," << b.x_n[i]
				 	  << "," << b.y_n[i];
			if (DIM > 2)
			writeData << "," << b.z_n[i];
			writeData << "\n";
		}
		writeData.close();
	}
	
	// Display body write time to console
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<-> Body data write time:              [%f s]\n", _time);
	
	return;
}

/**
 *  @brief Write the particle properties data.
 *  
 *  @param	_particle	The particle object to be wrote.
 *  @param	_name	Identifier name for saved file.
 * 	@param	_type	Flag to write the particle: [0] save all particle, 
 *  [1] save particle near body only.
 * 
 *  @headerfile save_data.hpp
 */
void save_data::save_par_state(const Particle &p, std::string name, int type)
{	
	// Time counting
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif
	
	// CODE LIES HERE
	if (DIM == 2){
		this->save_par_state_2D(p, name, type);
	}
	else if (DIM == 3){
		this->save_par_state_3D(p, name, type);
	}
	
	// Display particle write time to console
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<-> Particle data write time:          [%f s]\n", _time);

	return;
}


void save_data::save_par_state_2D(const Particle &p, std::string name, int type){
	// Create aliases to the writting flag
	bool _ngh_par = Pars::flag_save_ngh_par;
	bool _ngh_num = Pars::flag_save_ngh_num;
	bool _pressure = Pars::flag_pressure_calc;

	// internal object variables
	std::ofstream writeData;

	// Print message log!!!
	std::cout << FONT_GREEN << "\nSaving the particle state " << name << " ..." << FONT_RESET << "\n";
	name = "output/particle_state_" + name + ".csv";
	writeData.open(name.c_str());

	// Save the table header
	writeData << "" << "xp"
			  << "," << "yp"
			  << "," << "gpz"
			  << "," << "vor"
			  << "," << "up"
 			  << "," << "vp"
			  << "," << "sp"
			  << "," << "active"
			  << "," << "chi";
	if (_pressure)
	writeData << "," << "pressure";
	if (_ngh_num)
	writeData << "," << "ngh_num";
	if (_ngh_par)
	writeData << "," << "ngh_ID";
	writeData << "\n";

	// Save whole domain
	if (type == 0){
		for (int i = 0; i < p.num; i++){
			// Save the table data
			writeData << "" << p.x[i]
					  << "," << p.y[i]
					  << "," << p.gz[i]
					  << "," << p.vorticity[i]
					  << "," << p.u[i] - Pars::ubody 
					  << "," << p.v[i] - Pars::vbody
					  << "," << p.s[i]
					  << "," << p.isActive[i]
					  << "," << p.chi[i];
			if (_pressure)
			writeData << "," << p.P[i];
			if (_ngh_num)
			writeData << "," << p.neighbor[i].size();
			if (_ngh_par){
				writeData << "," << p.neighbor[i][0];
				for (size_t j = 1; j < p.neighbor[i].size(); j++)
				writeData << "|" << p.neighbor[i][j];
			}
			writeData << "\n";
		}
	}
	
	// Save particle inside the interest domain
	else if (type == 1){
		// Domain size limitation
		double domExt = Pars::par_dom_ext;	// Domain extension
		double domMin[DIM];		// Domain extreme minimum
		double domMax[DIM];		// Domain extreme maximum
		
		// Set initial value of data limitation
		domMin[0] = Pars::x_body_cen[0];
		domMax[0] = Pars::x_body_cen[0];
		domMin[1] = Pars::y_body_cen[0];
		domMax[1] = Pars::y_body_cen[0];

		// Evaluate the domain limitation
		for (int i = 0; i < N_BODY; i++){
			domMin[0] = std::min<double>(domMin[0], Pars::x_body_cen[i] - Pars::Lref[i]/2.0);
			domMax[0] = std::max<double>(domMin[0], Pars::x_body_cen[i] + Pars::Lref[i]/2.0);
			domMin[1] = std::min<double>(domMin[1], Pars::y_body_cen[i] - Pars::Lref[i]/2.0);
			domMax[1] = std::max<double>(domMin[1], Pars::y_body_cen[i] + Pars::Lref[i]/2.0);
		}
		// Give extention to the domain
		basis_loop(d){
			domMin[d] -= domExt;
			domMax[d] += domExt;
		}
		
		// Iterate through all particle
		for (int i = 0; i < p.num; i++){
			// Check the domain limitation at each basis
			if (p.x[i] < domMin[0] || p.y[i] < domMin[1]) continue;
			if (p.x[i] > domMax[0] || p.y[i] > domMax[1]) continue;

			// Save the table data
			writeData << "" << p.x[i]
					  << "," << p.y[i]
					  << "," << p.gz[i]
					  << "," << p.vorticity[i]
					  << "," << p.u[i] - Pars::ubody 
					  << "," << p.v[i] - Pars::vbody
					  << "," << p.s[i]
					  << "," << p.isActive[i]
					  << "," << p.chi[i];
			if (_pressure)
			writeData << "," << p.P[i];
			if (_ngh_num)
			writeData << "," << p.neighbor[i].size();
			if (_ngh_par){
				writeData << "," << p.neighbor[i][0];
				for (size_t j = 0; j < p.neighbor[i].size(); j++)
				writeData << "|" << p.neighbor[i][j];
			}
			writeData << "\n";
		}
	}
	writeData.close();
	return;
}

void save_data::save_par_state_3D(const Particle &p, std::string name, int type){
	// Create aliases to the writting flag
	bool _ngh_par = Pars::flag_save_ngh_par;
	bool _ngh_num = Pars::flag_save_ngh_num;
	bool _pressure = Pars::flag_pressure_calc;

	// internal object variables
	std::ofstream writeData;

	// Print message log!!!
	std::cout << FONT_GREEN << "\nSaving the particle state " << name << " ..." << FONT_RESET << "\n";
	name = "output/particle_state_" + name + ".csv";
	writeData.open(name.c_str());

	// Save the table header
	writeData <<  "" << "xp"
			  << "," << "yp"
			  << "," << "zp"
			  << "," << "vortx"
			  << "," << "vorty"
			  << "," << "vortz"
			  << "," << "vor"
			  << "," << "up"
 			  << "," << "vp"
			  << "," << "wp"
			  << "," << "sp"
			  << "," << "active"
			  << "," << "chi";
	if (_pressure)
	writeData << "," << "pressure";
	if (_ngh_num)
	writeData << "," << "ngh_num";
	if (_ngh_par)
	writeData << "," << "ngh_ID";
	writeData << "\n";

	// Save whole domain
	if (type == 0){
		for (int i = 0; i < p.num; i++){
			// Save the table data
			writeData << "" << p.x[i]
					  << "," << p.y[i]
					  << "," << p.z[i]
					  << "," << p.vortx[i]
					  << "," << p.vorty[i]
					  << "," << p.vortz[i]
					  << "," << p.vorticity[i]
					  << "," << p.u[i] - Pars::ubody
					  << "," << p.v[i] - Pars::vbody
					  << "," << p.w[i] - Pars::wbody
					  << "," << p.s[i]
					  << "," << p.isActive[i]
					  << "," << p.chi[i];
			if (_pressure)
			writeData << "," << p.P[i];
			if (_ngh_num)
			writeData << "," << p.neighbor[i].size();
			if (_ngh_par){
				writeData << "," << p.neighbor[i][0];
				for (size_t j = 1; j < p.neighbor[i].size(); j++)
				writeData << "|" << p.neighbor[i][j];
			}
			writeData << "\n";
		}
	}
	
	// Save particle inside the interest domain
	else if (type == 1){
		// Domain size limitation
		double domExt = Pars::par_dom_ext;	// Domain extension
		double domMin[3];		// Domain extreme minimum
		double domMax[3];		// Domain extreme maximum
		
		// Set initial value of data limitation
		domMin[0] = Pars::x_body_cen[0];
		domMax[0] = Pars::x_body_cen[0];
		domMin[1] = Pars::y_body_cen[0];
		domMax[1] = Pars::y_body_cen[0];
		domMin[2] = Pars::z_body_cen[0];
		domMax[2] = Pars::z_body_cen[0];
		
		// Evaluate the domain limitation
		for (int i = 0; i < N_BODY; i++){
			domMin[0] = std::min<double>(domMin[0], Pars::x_body_cen[i] - Pars::Lref[i]/2.0);
			domMax[0] = std::max<double>(domMin[0], Pars::x_body_cen[i] + Pars::Lref[i]/2.0);
			domMin[1] = std::min<double>(domMin[1], Pars::y_body_cen[i] - Pars::Lref[i]/2.0);
			domMax[1] = std::max<double>(domMin[1], Pars::y_body_cen[i] + Pars::Lref[i]/2.0);
			domMin[2] = std::min<double>(domMin[2], Pars::z_body_cen[i] - Pars::Lref[i]/2.0);
			domMax[2] = std::max<double>(domMin[2], Pars::z_body_cen[i] + Pars::Lref[i]/2.0);
		}

		// Give extention to the domain
		basis_loop(d){
			domMin[d] -= domExt;
			domMax[d] += domExt;
		}
		
		// Iterate through all particle
		for (int i = 0; i < p.num; i++){
			// Check the domain limitation at each basis
			if (p.x[i] < domMin[0] || p.y[i] < domMin[1]) continue;
			if (p.x[i] > domMax[0] || p.y[i] > domMax[1]) continue;
			if (p.z[i] < domMin[2] || p.z[i] > domMax[2]) continue;

			// Save the table data
			writeData << "" << p.x[i]
					  << "," << p.y[i]
					  << "," << p.z[i]
					  << "," << p.gx[i]
					  << "," << p.gy[i]
					  << "," << p.gz[i]
					  << "," << p.vorticity[i]
					  << "," << p.u[i] - Pars::ubody 
					  << "," << p.v[i] - Pars::vbody
					  << "," << p.w[i] - Pars::wbody
					  	<< "," << p.s[i]
					  << "," << p.isActive[i]
					  << "," << p.chi[i];
			if (_pressure)
			writeData << "," << p.P[i];
			if (_ngh_num)
			writeData << "," << p.neighbor[i].size();
			if (_ngh_par){
				writeData << "," << p.neighbor[i][0];
				for (size_t j = 0; j < p.neighbor[i].size(); j++)
				writeData << "|" << p.neighbor[i][j];
			}
			writeData << "\n";
		}
	}
	writeData.close();
	return;
}

/**
 *  @brief Write the gridNode properties data.
 *  
 *  @param	_gridNode	The grid node object to be wrote.
 *  @param	_name	Identifier name for saved file.
 * 	@param	_type	Flag to write the grid Node: [0] save all gridNode, 
 *  [1] save leaf grid Node only.
 * 
 *  @headerfile save_data.hpp
 */
void save_data::save_grid_node_state(const GridNode &_gridNode, std::string _name, int _type){
	// Cancel the grid saving if the flag turn off
	if (Pars::flag_save_node == false) return;

	// Console log
	std::cout << FONT_GREEN << "\nSaving the grid node state " << _name << " ..." << FONT_RESET << "\n";
	
	// Time counting
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif


	// Save the grid node data
	if (_type == 0)			// Save all data
	_gridNode.saveGrid(_gridNode, _name);

	else if (_type == 1)	// Save leaf only data
	_gridNode.saveLeafGrid(_gridNode, _name);

	// Display grid write time to console
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<-> Grid data write time:              [%f s]\n", _time);

	return;
}