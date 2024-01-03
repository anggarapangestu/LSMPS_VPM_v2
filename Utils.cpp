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
*/
void simUtil::stabilityEval(const Particle&par, std::vector<double>&max){
	// Initialize parameter
	double _courNum [3] = {0,0,0};  // [Current, Accumulative, Maximum]
	double _diffNum [3] = {0,0,0};  // [Current, Accumulative, Maximum]
	double _stabCrt [3] = {0,0,0};  // [Current, Accumulative, Maximum]
	if (DIM == 2){
	for (int _i = 0; _i < par.num; _i++){
		// Calculating courant number
		_courNum[0] = std::sqrt(std::pow(par.u[_i],2) + std::pow(par.v[_i],2)) * Pars::dt / par.s[_i];
		_courNum[1] += _courNum[0];
		_courNum[2] = _courNum[2] > _courNum[0] ? _courNum[2] : _courNum[0];

		// Calculating diffusion number
		_diffNum[0] = Pars::NU * Pars::dt / (par.s[_i] * par.s[_i]);
		_diffNum[1] += _diffNum[0];
		_diffNum[2] = _diffNum[2] > _diffNum[0] ? _diffNum[2] : _diffNum[0];
		
		// Calculating stability criteria
		_stabCrt[0] = (std::abs(par.gz[_i]))/Pars::NU;
		_stabCrt[1] += _stabCrt[0];
		_stabCrt[2] = _stabCrt[2] > _stabCrt[0] ? _stabCrt[2] : _stabCrt[0];
	}}
	else if (DIM == 3){
	for (int _i = 0; _i < par.num; _i++){
		// Calculating courant number
		_courNum[0] = std::sqrt(par.u[_i]*par.u[_i] + par.v[_i]*par.v[_i] + par.w[_i]*par.w[_i]) * Pars::dt / par.s[_i];
		_courNum[1] += _courNum[0];
		_courNum[2] = _courNum[2] > _courNum[0] ? _courNum[2] : _courNum[0];

		// Calculating diffusion number
		_diffNum[0] = Pars::NU * Pars::dt / (par.s[_i] * par.s[_i]);
		_diffNum[1] += _diffNum[0];
		_diffNum[2] = _diffNum[2] > _diffNum[0] ? _diffNum[2] : _diffNum[0];
		
		// Calculating stability criteria
		_stabCrt[0] = std::pow(par.s[_i], 3.0) * std::sqrt(par.vortx[_i]*par.vortx[_i] + par.vorty[_i]*par.vorty[_i] + par.vortz[_i]*par.vortz[_i])/Pars::NU;
		_stabCrt[1] += _stabCrt[0];
		_stabCrt[2] = _stabCrt[2] > _stabCrt[0] ? _stabCrt[2] : _stabCrt[0];
	}}
	
	// Average value
	_courNum[0] = _courNum[1] / par.num;
	_diffNum[0] = _diffNum[1] / par.num;
	_stabCrt[0] = _stabCrt[1] / par.num;
	
	// Displaying the value
	printf("Average courant number (C_av)           : %8.4f \n", _courNum[0]);
	printf("Average diffusion number (Phi_av)       : %8.4f \n", _diffNum[0]);
	printf("Average stability criteria (Re_h)       : %8.4f \n", _stabCrt[0]);
	printf("Max courant number (C_max)              : %8.4f \n", _courNum[2]);
	printf("Max diffusion number (Phi_max)          : %8.4f \n", _diffNum[2]);
	printf("Max stability criteria (Re_h)           : %8.4f \n", _stabCrt[2]);

	// Input the maximum value
	max.push_back(_courNum[2]);
	max.push_back(_diffNum[2]);
	max.push_back(_stabCrt[2]);
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
	
	printf("\n<!> Estimation time to finish run:     %9.3f s", curr_comp_time*double(Pars::max_iter - step));
	if (est_time_d == 0){
		printf("\n<!> Estimation time to finish run: %2dh %2dm %5.2f s", est_time_h, est_time_m, est_time_s);
	}else{
		printf("\n<!> Estimation time to finish: %2dd %2dh %2dm %5.2f s", est_time_d, est_time_h, est_time_m, est_time_s);
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
    // Internal variable
	std::vector<double> res(par.num,0.0);
	std::vector<double> vor(par.num,0.0);
    std::vector<double> vorx(par.num,0.0);
	std::vector<double> vory(par.num,0.0);
	std::vector<double> vorz(par.num,0.0);

	// Calculating the residual
	if (DIM == 2){
		// Internal tools
		LSMPSa lsmpsa_du, lsmpsa_dv;

		// Calculate the differential du/dx
		lsmpsa_du.set_LSMPS(par.x, par.y, par.s, par.u, par.neighbor);
		auto _dudx = lsmpsa_du.get_ddx();

		// Calculate the differential dv/dy
		lsmpsa_dv.set_LSMPS(par.x, par.y, par.s, par.v, par.neighbor);
		auto _dvdy = lsmpsa_dv.get_ddy();

		// Calculating vorticity
		auto _dudy = lsmpsa_du.get_ddy();
		auto _dvdx = lsmpsa_dv.get_ddx();

		for (int i = 0; i < par.num; i++){
			res[i] = _dudx[i] + _dvdy[i];
			vor[i] = _dvdx[i] - _dudy[i];
		}
	}
	else if (DIM == 3){
		// Internal tools
		LSMPSa lsmpsa_du, lsmpsa_dv, lsmpsa_dw;

		// Calculate the differential du/dx
		lsmpsa_du.set_LSMPS_3D(par.x, par.y, par.z, par.s, par.u, par.x, par.y, par.z, par.s, par.u, par.neighbor);
		auto _dudx = lsmpsa_du.get_ddx();

		// Calculate the differential dv/dy
		lsmpsa_dv.set_LSMPS_3D(par.x, par.y, par.z, par.s, par.v, par.x, par.y, par.z, par.s, par.v, par.neighbor);
		auto _dvdy = lsmpsa_dv.get_ddy();

		// Calculate the differential dv/dy
		lsmpsa_dw.set_LSMPS_3D(par.x, par.y, par.z, par.s, par.w, par.x, par.y, par.z, par.s, par.w, par.neighbor);
		auto _dwdz = lsmpsa_dw.get_ddz();

		// Calculating vorticity
		auto _dudy = lsmpsa_du.get_ddy();
		auto _dudz = lsmpsa_du.get_ddz();
		auto _dvdx = lsmpsa_dv.get_ddx();
		auto _dvdz = lsmpsa_dv.get_ddz();
		auto _dwdx = lsmpsa_dw.get_ddx();
		auto _dwdy = lsmpsa_dw.get_ddy();

		for (int i = 0; i < par.num; i++){
			res[i] = _dudx[i] + _dvdy[i] + _dwdz[i];
			vorx[i] = _dwdy[i] - _dvdz[i];
			vory[i] = _dudz[i] - _dwdx[i];
			vorz[i] = _dvdx[i] - _dudy[i];
		}
	}

    // Saving the residual data
    std::ofstream data;

    // Creating save file name
    std::string name;
    name  = "output/residual_";
    simUtil util_step;
    std::string number = util_step.saveName(step);
    name += number + ".csv";
    
    // Write data header
    data.open(name);
	if(DIM == 2)
    data << "x,y,vor,res\n";
	else if(DIM == 3)
	data << "x,y,z,vorx,vory,vorz,res\n";
    
    // Write data content
	if (DIM == 2){
    for (int i = 0; i < par.num; i++){
    data << ""  << par.x[i]
         << "," << par.y[i]
         << "," << vor[i]
         << "," << res[i]
         << "\n";
    }}
	else if (DIM == 3){
	for (int i = 0; i < par.num; i++){
    data << ""  << par.x[i]
         << "," << par.y[i]
		 << "," << par.z[i]
         << "," << vorx[i]
		 << "," << vory[i]
		 << "," << vorz[i]
         << "," << res[i]
         << "\n";
    }}
	data.close();
    
    return;
}