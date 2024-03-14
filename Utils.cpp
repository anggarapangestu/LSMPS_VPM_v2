/* ---------------- PROGRAM DESCRIPTION ----------------
	> Consisted of utilities function for main program
*/

#include "Utils.hpp"
#include "src/LSMPS/LSMPSa.hpp"

/**
 *  @brief  Start the simulation counter number. Set up the number of current 
 *  iteration digit.
 *         
 *  @param  step The iteration step.
*/
void simUtil::startCounter(int step){
	// Calculate the digit length of iteration step
	this->iterDigitLen = 1;
	if (step != 0){
		iterDigitLen = 1 + std::floor(std::log10(step));
	}
	return;
}

/**
 *  @brief  Printing the iteration number header to the console.
 *         
 *  @param  step The iteration step to be printed.
*/
void simUtil::printHeader(int step){
	// Calculate the digit length of iteration step
	int dig_len = 1;
	if (step != 0) dig_len = 1 + std::floor(std::log10(step));

	// Printing the iteration HEADER
	printf("\n\n+");
	for(int _i = 0; _i < 16 - dig_len/2; _i ++){
		printf("-");
	}
	printf(" Iteration Step %d ", step);
	for(int _i = 0; _i < 16 - dig_len/2 - dig_len%2; _i ++){
		printf("-");
	}
	printf("+\n");
	return;
}

/**
 *  @brief  Set the iteration label for file name.
 *         
 *  @param  step The iteration step.
 * 
 *  @return The iteration label name.
*/
std::string simUtil::saveName(const int step){
	std::string DataName;
	// Calculate the digit length of iteration step
	int dig_len = 1;
	if (step != 0) dig_len = 1 + std::floor(std::log10(step));
	
	// Add the leading zero
	int addDigit = Pars::max_dig_len - dig_len;
	for (int _spc = 0; _spc < addDigit; _spc++) DataName.append("0");

	// Combine the leading zero to step number
	DataName.append(std::to_string(step));
	return DataName;
}

/**
 *  @brief  Display the stability into the console.
 *         
 *  @param  _particle The particle data for stability calculaltion.
 *  @param  _maxValue [OUTPUT] Stability criteria global maximum value.
 *  @param  _step  The current iteration step
*/
void simUtil::stabilityEval(const Particle &par, std::vector<double*> &max, const int &step){
	// A prompt in stability evaluation
	printf("Evaluating stability ...\n");
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

	// Internal container
	// [(1)Current Evaluated, (2)Accumulative in Domain, (3)Global Maximum]
	double _courNum [3] = {0,0,0};  // The courant number container
	double _diffNum [3] = {0,0,0};  // The diffusion number container
	double _meshNum [3] = {0,0,0};  // The mesh reynold number container

	// Calculation operation
	#if (DIM == 2)
		for (int _i = 0; _i < par.num; _i++){
			// Calculating courant number
			double Velocity = std::sqrt(par.u[_i]*par.u[_i] + par.v[_i]*par.v[_i]);
			_courNum[0] = Velocity * Pars::dt / par.s[_i];
			_courNum[1] += _courNum[0];
			_courNum[2] = _courNum[2] > _courNum[0] ? _courNum[2] : _courNum[0];

			// Calculating diffusion number
			_diffNum[0] = Pars::NU * Pars::dt / (par.s[_i] * par.s[_i]);
			_diffNum[1] += _diffNum[0];
			_diffNum[2] = _diffNum[2] > _diffNum[0] ? _diffNum[2] : _diffNum[0];
			
			// Calculating vorticity reynolds criteria
			_meshNum[0] = (std::abs(par.gz[_i]))/Pars::NU;
			_meshNum[1] += _meshNum[0];
			_meshNum[2] = _meshNum[2] > _meshNum[0] ? _meshNum[2] : _meshNum[0];
		}
	#elif (DIM == 3)
		for (int _i = 0; _i < par.num; _i++){
			// Calculating courant number
			double Velocity = std::sqrt(par.u[_i]*par.u[_i] + par.v[_i]*par.v[_i] + par.w[_i]*par.w[_i]);
			_courNum[0] = Velocity * Pars::dt / par.s[_i];
			_courNum[1] += _courNum[0];
			_courNum[2] = _courNum[2] > _courNum[0] ? _courNum[2] : _courNum[0];

			// Calculating diffusion number
			_diffNum[0] = Pars::NU * Pars::dt / (par.s[_i] * par.s[_i]);
			_diffNum[1] += _diffNum[0];
			_diffNum[2] = _diffNum[2] > _diffNum[0] ? _diffNum[2] : _diffNum[0];
			
			// Calculating stability criteria
			// double Vorticity = std::sqrt(par.vortx[_i]*par.vortx[_i] + par.vorty[_i]*par.vorty[_i] + par.vortz[_i]*par.vortz[_i]);
			_meshNum[0] = par.vorticity[_i] * par.s[_i] * par.s[_i] / Pars::NU;
			_meshNum[1] += _meshNum[0];
			_meshNum[2] = _meshNum[2] > _meshNum[0] ? _meshNum[2] : _meshNum[0];
		}
	#endif
	
	// Calculate the average value (put into the current evaluation container)
	_courNum[0] = _courNum[1] / par.num;
	_diffNum[0] = _diffNum[1] / par.num;
	_meshNum[0] = _meshNum[1] / par.num;


	// Determine the current step is initialization or not
	bool init = false;
	if      (Pars::opt_start_state == 0) init = (step == 0);
	else if (Pars::opt_start_state == 1) init = (step == (Pars::resume_step + 1));

	// Update the global stability criteria maximum value
	// if (init == true){
	// 	this->max_C    = _courNum[2];
	// 	this->max_Phi  = _diffNum[2];
	// 	this->max_Re_h = _meshNum[2];
	// }else{
	// 	this->max_C    = std::max(_courNum[2], this->max_C);
	// 	this->max_Phi  = std::max(_diffNum[2], this->max_Phi);
	// 	this->max_Re_h = std::max(_meshNum[2], this->max_Re_h);
	// }
	if (init == true){
		*(max[0]) = _courNum[2];
		*(max[1]) = _diffNum[2];
		*(max[2]) = _meshNum[2];
	}else{
		*(max[0]) = std::max(*(max[0]), _courNum[2]);
		*(max[1]) = std::max(*(max[1]), _diffNum[2]);
		*(max[2]) = std::max(*(max[2]), _meshNum[2]);
	}

	// Save the stability data
	if (Pars::flag_save_stability == true){
		// Initialize the writter
		std::ofstream _write;
		std::string fileName = "output/Stability_History";
		if      (Pars::opt_start_state == 0) fileName += ".csv";
		else if (Pars::opt_start_state == 1) fileName += "_Resume.csv";
		
		// Save the header
		if (init == true){
			_write.open(fileName);
			_write << "time,C_max,Phi_max,Re_h_max\n";
			_write.close();
		}
		
		// Save the data
		_write.open(fileName, std::ofstream::out | std::ofstream::app);
		_write <<  "" << step*Pars::dt 
			   << "," << _courNum[2]
			   << "," << _diffNum[2]
			   << "," << _meshNum[2]
			   << "\n";
		_write.close();
	}

	// Particle generation initialization summary time display
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<-> Stability evaluation comp. time.   [%f s]\n\n", _time);

	// Displaying the value
	printf("Average courant number (C_av)           : %8.4f \n", _courNum[0]);
	printf("Average diffusion number (Phi_av)       : %8.4f \n", _diffNum[0]);
	printf("Average stability criteria (Re_h)       : %8.4f \n", _meshNum[0]);
	printf("Max courant number (C_max)              : %8.4f \n", _courNum[2]);
	printf("Max diffusion number (Phi_max)          : %8.4f \n", _diffNum[2]);
	printf("Max stability criteria (Re_h)           : %8.4f \n", _meshNum[2]);
	
	return;
}

/**
 *  @brief  Display the prediction time to finish the simulation.
 *         
 *  @param  step The current iteration step.
 *  @param  _currTime The current iteration computational time.
*/
void simUtil::predictCompTime(int step, double curr_comp_time){
	// The value of "curr_comp_time" is in second (s)
	// Internal variable
	int est_time_d, est_time_h, est_time_m; double est_time_s;
	est_time_s = curr_comp_time * double(Pars::max_iter - step - 1);
	
	// Calculate Day
	est_time_d = int(est_time_s / (24 * 60 * 60));
	est_time_s -= est_time_d * (24 * 60 * 60);
	// Calculate Hour
	est_time_h = int(est_time_s / (60 * 60));
	est_time_s -= est_time_h * (60 * 60);
	// Calculate Minute
	est_time_m = int(est_time_s / (60));
	est_time_s -= est_time_m * (60);
	
	// The simulation estimation is limited to 999 days of simulation
	printf("\n<!> Estimation time to finish run : %12.2f s", curr_comp_time*double(Pars::max_iter - step));
	if (est_time_d == 0){
		printf("\n<!> Estimation time to finish run :    %2dh %2dm %2ds", est_time_h, est_time_m, (int)est_time_s);
	}else{
		printf("\n<!> Estimation time to finish run :  %3ddays %2dhrs", est_time_d, est_time_h);
	}

	printf("\n");
	return;
}

/**
 *  @brief  Calculating the simulation residual and write into file.
 *         
 *  @param  _particle The particle data for residual calculation.
 *  @param  step The current iteration step.
*/
void simUtil::saveResidual(Particle &par, int step){
    // A prompt in stability evaluation
	printf("\nCalculating residuals ...\n");
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif
	
	// Internal variable
	std::vector<double> mass_cnsv(par.num,0.0);
	std::vector<double> vor(par.num,0.0);
    std::vector<double> vorx(par.num,0.0);
	std::vector<double> vory(par.num,0.0);
	std::vector<double> vorz(par.num,0.0);

	// Calculating the residual
	if (DIM == 2){
		// Internal tools
		// LSMPSa lsmpsa_du, lsmpsa_dv;
		LSMPSa lsmpsa_tool;

		// Calculate the differential du/dx
		lsmpsa_tool.set_LSMPS(par.x, par.y, par.s, par.u, par.neighbor);
		auto _dudx = lsmpsa_tool.get_ddx();
		auto _dudy = lsmpsa_tool.get_ddy();

		// Calculate the differential dv/dy
		lsmpsa_tool.set_LSMPS(par.x, par.y, par.s, par.v, par.neighbor);
		auto _dvdx = lsmpsa_tool.get_ddx();
		auto _dvdy = lsmpsa_tool.get_ddy();

		// Calculating vorticity and mass conservation
		for (int i = 0; i < par.num; i++){
			mass_cnsv[i] = _dudx[i] + _dvdy[i];
			vor[i] = par.vorticity[i] - (_dvdx[i] - _dudy[i]);
		}
	}
	else if (DIM == 3){
		// Internal tools
		// LSMPSa lsmpsa_du, lsmpsa_dv, lsmpsa_dw;
		LSMPSa lsmpsa_tool;

		// Calculate the differential du/dx
		lsmpsa_tool.set_LSMPS_3D(par.x, par.y, par.z, par.s, par.u, 
								 par.x, par.y, par.z, par.s, par.u, par.neighbor);
		auto _dudx = lsmpsa_tool.get_ddx();
		auto _dudy = lsmpsa_tool.get_ddy();
		auto _dudz = lsmpsa_tool.get_ddz();

		// Calculate the differential dv/dy
		lsmpsa_tool.set_LSMPS_3D(par.x, par.y, par.z, par.s, par.v, 
								 par.x, par.y, par.z, par.s, par.v, par.neighbor);
		auto _dvdx = lsmpsa_tool.get_ddx();
		auto _dvdy = lsmpsa_tool.get_ddy();
		auto _dvdz = lsmpsa_tool.get_ddz();

		// Calculate the differential dv/dy
		lsmpsa_tool.set_LSMPS_3D(par.x, par.y, par.z, par.s, par.w, 
								 par.x, par.y, par.z, par.s, par.w, par.neighbor);
		auto _dwdx = lsmpsa_tool.get_ddx();
		auto _dwdy = lsmpsa_tool.get_ddy();
		auto _dwdz = lsmpsa_tool.get_ddz();

		// Calculating vorticity

		for (int i = 0; i < par.num; i++){
			mass_cnsv[i] = _dudx[i] + _dvdy[i] + _dwdz[i];
			vorx[i] = par.vortx[i] - (_dwdy[i] - _dvdz[i]);
			vory[i] = par.vorty[i] - (_dudz[i] - _dwdx[i]);
			vorz[i] = par.vortz[i] - (_dvdx[i] - _dudy[i]);
		}
	}

	// Saving the residual field data
	if ((step % Pars::save_inv) == 0){
		// Initialize the writter
		std::ofstream data;
		std::string name;
		name  = "output/residual_";
		simUtil util_step;
		std::string number = util_step.saveName(step);
		name += number + ".csv";
		
		// Write data
		data.open(name);
		
		if (DIM == 2){
			// Write header
			data << "xp,yp,sp,vor,mass_cnsv\n";
			// Write data content
			for (int i = 0; i < par.num; i++){
			data << ""  << par.x[i]
				<< "," << par.y[i]
				<< "," << par.s[i]
				<< "," << vor[i]
				<< "," << mass_cnsv[i]
				<< "\n";
			}
		}
		else if (DIM == 3){
			// Write header
			data << "xp,yp,zp,sp,vorx,vory,vorz,mass_cnsv\n";
			// Write data content
			for (int i = 0; i < par.num; i++){
			data << ""  << par.x[i]
				<< "," << par.y[i]
				<< "," << par.z[i]
				<< "," << par.s[i]
				<< "," << vorx[i]
				<< "," << vory[i]
				<< "," << vorz[i]
				<< "," << mass_cnsv[i]
				<< "\n";
			}
		}
		data.close();
	}


	// Calculate the global residual
	double resMC = 0.0;		// Residual of mass conservation
	double resVor = 0.0;	// Residual of vorticity
	
	// Residual Calculation
	#if (DIM == 2)
		for (int i = 0; i < par.num; i++){
			// Residual of mass conservation E_mc = sum(abs(div(U))*Area) / (U_inf * D)
			resMC += std::abs(mass_cnsv[i]) * par.s[i]*par.s[i];

			// Residual of vorticity E_vor = sum(|curl(U) - omega|^2 * Area)  / (U_inf^2)
			resVor += vor[i]*vor[i] * par.s[i]*par.s[i];
		}
		
		// Normalize the value
		resMC /= (Pars::U_inf * Pars::Df);
		resVor /= (Pars::U_inf * Pars::U_inf);
	
	#elif (DIM == 3)
		for (int i = 0; i < par.num; i++){
			// Residual of mass conservation E_mc = sum(abs(div(U))*Vol) / (U_inf * D^2)
			resMC += std::abs(mass_cnsv[i]) * par.s[i]*par.s[i]*par.s[i];

			// Residual of vorticity E_vor = sum(|curl(U) - omega|^2 * Vol)  / (U_inf^2*D)
			resVor += (vorx[i]*vorx[i] + vory[i]*vory[i] + vorz[i]*vorz[i]) * par.s[i]*par.s[i]*par.s[i];
		}
		
		// Normalize the value
		resMC /= (Pars::U_inf * Pars::Df * Pars::Df);
		resVor /= (Pars::U_inf * Pars::U_inf * Pars::Df);
	#endif
	
	// Saving the residual history data
	std::ofstream _write;
	std::string fileName = "output/Residual_History";
	if      (Pars::opt_start_state == 0) fileName += ".csv";
	else if (Pars::opt_start_state == 1) fileName += "_Resume.csv";
	
	// Determine the current step is initialization or not
	bool init = false;
	if      (Pars::opt_start_state == 0) init = (step == 0);
	else if (Pars::opt_start_state == 1) init = (step == (Pars::resume_step + 1));

	// Save the header
	if (init == true){
		_write.open(fileName);
		_write << "time,massRes,vorticityRes\n";
		_write.close();
	}
	
	// Save the data
	_write.open(fileName, std::ofstream::out | std::ofstream::app);
	_write <<  "" << step*Pars::dt 
			<< "," << resMC
			<< "," << resVor
			<< "\n";
	_write.close();


	// Particle generation initialization summary time display
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<-> Stability evaluation comp. time.   [%f s]\n", _time);
    
    return;
}

void simUtil::addVorPertubation(Particle &p){
	// Add the peturbation right in front of the object
	// with value of vorticity = |omega|_max * 10 ^-3

	// The peturbation function is a section of cosine function
	//    omega_ptb = (K/2)*(1 + cos(pi*(r/R))) ; (r < R)
	//  Where : 
	//    > K : vorticity peak multiplication
	//    > R : peturbation radius -> (sigma * Rs)
	//    > r : radius from peturbation location
	//    > sigma : local max particle size
	//    > Rs : peturbation radius scale

	// Set the peturbation constant
	double Rs, R, K;
	Rs = Pars::r_sup;

	// Set the peturbation location
	double ptbCoord[DIM];
	ptbCoord[0] = -0.7;    // x location
	ptbCoord[1] = -0.6;    // y location
	#if DIM == 3
	ptbCoord[2] = 0.0;     // z location
	#endif

	// Evaluate the box size
	double evalBoxSize = Pars::sigma * 16;
	double locMaxParSize = Pars::sigma;		// The local maximum particle size
	double vorMaxVal = 0.0;		// The local maximum particle size

	// Find the local maximum particle size and maximum vorticity
	for (int i = 0; i < p.num; i++){
		// Update the maximum value of vorticity
		if (std::abs(p.vorticity[i]) > vorMaxVal) vorMaxVal = std::abs(p.vorticity[i]);

		// Check whether the particle is outside evaluation box
		if (p.x[i] > ptbCoord[0] + evalBoxSize ||
			p.x[i] < ptbCoord[0] - evalBoxSize)
		continue;

		if (p.y[i] > ptbCoord[1] + evalBoxSize ||
			p.y[i] < ptbCoord[1] - evalBoxSize)
		continue;

		#if DIM == 3
		if (p.z[i] > ptbCoord[2] + evalBoxSize ||
			p.z[i] < ptbCoord[2] - evalBoxSize)
		continue;
		#endif

		// Evaluate the particle size
		if (p.s[i] > locMaxParSize) locMaxParSize = p.s[i];
	}

	// Update the evaluation box size and peturbation parameter
	evalBoxSize = locMaxParSize * Rs * 1.2;
	R = locMaxParSize * Rs;
	K = vorMaxVal / 10.0;

	// Calculate Set the peturbation value
	#pragma omp parallel for
	for (int i = 0; i < p.num; i++){
		// [CHECK 1] Check whether the particle is outside evaluation box
		if (p.x[i] > ptbCoord[0] + evalBoxSize ||
			p.x[i] < ptbCoord[0] - evalBoxSize)
		continue;

		if (p.y[i] > ptbCoord[1] + evalBoxSize ||
			p.y[i] < ptbCoord[1] - evalBoxSize)
		continue;

		#if DIM == 3
		if (p.z[i] > ptbCoord[2] + evalBoxSize ||
			p.z[i] < ptbCoord[2] - evalBoxSize)
		continue;
		#endif

		// Update the vorticity
		double r2 = 0.0;
		double dx, dy, dz;
		// Add the x component
			dx = p.x[i] - ptbCoord[0];
			r2 += dx*dx;
		// Add the y component
			dy = p.y[i] - ptbCoord[1];
			r2 += dy*dy;
		// Add the z component
		#if DIM == 3
			dz = p.z[i] - ptbCoord[2];
			r2 += dz*dz;
		#endif
		double r = std::sqrt(r2);

		// [CHECK 2] Check the location of particle related to the peturbation location
		if (r > R) continue;

		// Calculate the local vorticity
		#if (DIM == 2)
			// Peturbation for 2D space (update on vort z)
			double ptbVor = K/2.0 * (1 + std::cos(M_PI * r / R));
			p.vorticity[i] += ptbVor;
			p.gz[i] = p.vorticity[i] * p.s[i]*p.s[i];
		#elif (DIM == 3)
			// Peturbation for 3D space (update on vort z)
			double rTh = std::sqrt(dx*dx + dy*dy);
			double ptbVorTh = 0.5 * (1 + std::cos(M_PI * rTh / std::sqrt(R*R - dz*dz)));
			double ptbVorZ = K/2.0 * (1 + std::cos(M_PI * dz / R));
			p.vortz[i] += ptbVorZ*ptbVorTh;
		#endif
	}

	return;
}