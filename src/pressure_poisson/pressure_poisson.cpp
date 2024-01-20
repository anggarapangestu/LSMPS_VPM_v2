#include "pressure_poisson.hpp"

#ifndef INCLUDED_NEIGHBOR
#include "../neighbor/neighbor.hpp"
#endif

void pressure_poisson::get_pressure(Particle& p){
	// Pressure calculation option
	int form_type = 1;		// Formulation type: 0:= Velocity form; 1:= Vorticity form (OR 0:= Direct Pressure; 1:= Total Pressure)
	bool inc_pen = true;	// Include the penalization into the calculation
	bool cp_calc = false;	// Toogle for the pressure coefficient calculation
	bool src_save = false;	// Toogle for saving the pressure source

	// Update the start clock
	clock_t _timer = clock();
				
	// Solver computation
	printf("\nCalculate Pressure ...\n");
	printf("<+> Calculate the pressure source\n");
	
	// Resize the velocity of particle
	p.P.clear();
	p.P.resize(p.num,0.0);

	// ================================================================================
	// --------------------------- FMM PRESSURE CALCULATION ---------------------------
	// ================================================================================
	// Internal variable
	fmm2D FMM_step;			// The FMM operator class
	std::vector<std::vector<double>> particle_POS(p.num, std::vector<double>(DIM,0.0));	// Particle position
	std::vector<double> particle_SRC(p.num, 0.0);	// The particle source
	std::vector<bool> particle_mark(p.num, false);	// The particle active mark
	
	// Time manager
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif


	// Initialization of the temporary source
	std::vector<double> _SRC(p.num, 0.0);
	
	// Variable for penalization correction
	std::vector<double> Elm_x(p.num,0.0);
	std::vector<double> Elm_y(p.num,0.0);
	double _lambda = 0.0;
	if(inc_pen == true){
		_lambda = 1.0/Pars::dt;
		// _lambda = Pars::lambda;
	}

	// [1] Calculcate the source 
	// *************************
	if (form_type == 0){
		printf("<+> Velocity based pressure source\n");
		// [First model of the velocity form]
		// Calculate the velocity part
		LSMPSa lsmpsa_du;     // module to calculate spatial differential of u velocity
		LSMPSa lsmpsa_dv;     // module to calculate spatial differential of v velocity
		lsmpsa_du.set_LSMPS(p.x, p.y, p.s, p.u, p.neighbor);
		lsmpsa_dv.set_LSMPS(p.x, p.y, p.s, p.v, p.neighbor);
		std::vector<double> _dudx = lsmpsa_du.get_ddx();
		std::vector<double> _dudy = lsmpsa_du.get_ddy();
		std::vector<double> _dvdx = lsmpsa_dv.get_ddx();
		
		// Calculcate the pressure source
		for (int i = 0; i < p.num; i++){
			// Here calculate 2 * (dudx^2 + dudy * dvdx)
			_SRC[i] = - 2.0 * (_dudx[i] * _dudx[i] + _dudy[i] * _dvdx[i]);
		}

		// Calculate the penalization part
		if(inc_pen == true){
			for (int i = 0; i < p.num; i++){
				Elm_x[i] = _lambda*p.chi[i]*(Pars::ubody - p.u[i]);
				Elm_y[i] = _lambda*p.chi[i]*(Pars::vbody - p.v[i]);
			}

			// The result of LSMPS calculation of the penalization part
			LSMPSa lsmpsa_dElmx;     // module to calculate spatial differential of x element vector
			LSMPSa lsmpsa_dElmy;     // module to calculate spatial differential of y element vector
			lsmpsa_dElmx.set_LSMPS(p.x, p.y, p.s, Elm_x, p.neighbor);
			lsmpsa_dElmy.set_LSMPS(p.x, p.y, p.s, Elm_y, p.neighbor);
			std::vector<double> _dElmxdx = lsmpsa_dElmx.get_ddx();
			std::vector<double> _dElmydy = lsmpsa_dElmy.get_ddy();

			// Calculcate the pressure source
			for (int i = 0; i < p.num; i++){
				// Here calculate 2 * (dudx^2 + dudy * dvdx)
				_SRC[i] += (_dElmxdx[i] + _dElmydy[i]);
			}
		}
	}
	else if (form_type == 1){
		printf("<+> Vorticity based pressure source\n");
		// [Second model of the vorticity form]
		// Calculate the velocity norm gradient
		std::vector<double> U2(p.num,0.0);
		for (int i = 0; i < p.num; i++){
			U2[i] = p.u[i] * p.u[i] + p.v[i] * p.v[i];
		}
		// Calculate the gradient
		LSMPSa lsmpsa_dU2;     // module to calculate spatial differential of u velocity
		lsmpsa_dU2.set_LSMPS(p.x, p.y, p.s, U2, p.neighbor);
		std::vector<double> _dU2dx = lsmpsa_dU2.get_ddx();
		std::vector<double> _dU2dy = lsmpsa_dU2.get_ddy();

		// [Second model of the vorticity form]
		// Calculate the value at each vector element
		for (int i = 0; i < p.num; i++){
			// // Method 1
			Elm_x[i] = -0.5*_dU2dx[i] + p.v[i]*p.vorticity[i] + _lambda*p.chi[i]*(Pars::ubody - p.u[i]);
			Elm_y[i] = -0.5*_dU2dy[i] - p.u[i]*p.vorticity[i] + _lambda*p.chi[i]*(Pars::vbody - p.v[i]);
			// // Method 2
			// Elm_x[i] = p.v[i]*(_dvdx[i] - _dudy[i]) + _lambda*p.chi[i]*(Pars::ubody - p.u[i]);
			// Elm_y[i] = -p.u[i]*(_dvdx[i] - _dudy[i]) + _lambda*p.chi[i]*(Pars::vbody - p.v[i]);
		}

		// The result of LSMPS calculation of the divergence
		LSMPSa lsmpsa_dElmx;     // module to calculate spatial differential of x element vector
		LSMPSa lsmpsa_dElmy;     // module to calculate spatial differential of y element vector
		lsmpsa_dElmx.set_LSMPS(p.x, p.y, p.s, Elm_x, p.neighbor);
		lsmpsa_dElmy.set_LSMPS(p.x, p.y, p.s, Elm_y, p.neighbor);
		std::vector<double> _dElmxdx = lsmpsa_dElmx.get_ddx();
		std::vector<double> _dElmydy = lsmpsa_dElmy.get_ddy();

		// Calculcate the divergenve total pressure source
		for (int i = 0; i < p.num; i++){
			// Here calculate (dElmxdx + dElmydy)
			_SRC[i] = (_dElmxdx[i] + _dElmydy[i]);
		}
	}

	// [2] Determine active source
	// ***************************
	// Find the maximum value
	double SRC_max;
	SRC_max = std::abs(_SRC[0]);
    for (int i = 1; i < p.num; i++){
        SRC_max = SRC_max > std::abs(_SRC[i]) ? SRC_max : std::abs(_SRC[i]);
    }

	// [3] Assign each FMM variable
	// ****************************
	for (int i = 0; i < p.num; i++){
		// Particle position
		particle_POS[i][0] = p.x[i];
		particle_POS[i][1] = p.y[i];

		// Particle source
		particle_SRC[i] = (_SRC[i] / (2 * M_PI)) * (p.s[i] * p.s[i]);
		
		// Particle mark
		particle_mark[i] = std::abs(_SRC[i]) >= SRC_max * 1.0e-3 ? true : false;
		particle_mark[i] = particle_mark[i] | p.isActive[i];
	}

	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<+> Initialization   : %f s\n", _time);

	// DEBUGGING: save data value
	if (src_save == true){
		this->save_pressure_source(p.x,p.y,_SRC,particle_mark);
	}

	// [4] Initialization of FMM Tree
	// ******************************
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime();
    #endif
	// Initialization Tree Map
	// if (this->init == true){
		treeData.initializeTree(particle_POS, particle_mark);
		// treeData.createTree(particle_POS, particle_mark);
		// this->init = false;
	// }else{
	// 	treeData.updateCell(particle_POS, particle_mark);
	// }

	// Tree Map Saved
	if (false){
		std::string name = "P";
		name += std::to_string(this->counter);
		this->counter ++;
		treeData.saveTree(treeData, name);	
	}

	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        span = std::chrono::system_clock::now() - tick;
        _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<+> Tree finished in : %f s\n", _time);

	// [5] FMM Potential Calculation
	// *****************************
	// Calculate the tree data
	FMM_step.calcPotential(treeData, particle_POS, particle_mark, particle_SRC);
	
	// Retrieve the potential value
	std::vector<double> phi = FMM_step.get_Potential();		// phi = P/rho (type 1) or H (type 2)

	// Update the pressure
	// if (form_type == 0){			// Velocity form
		for (int i = 0; i < p.num; i++)
		{
			p.P[i] = Pars::RHO * phi[i];
			// p.P[i] = _SRC[i];
			
		}
	// }else if (form_type == 1){		// Vorticity form
	// 	for (size_t i = 0; i < p.num; i++)
	// 	{
	// 		p.P[i] = Pars::RHO * (phi[i] - 0.5 * (p.u[i]*p.u[i] + p.v[i]*p.v[i])
	// 		                             + 0.5 * (Pars::U_inf*Pars::U_inf));
	// 		// p.P[i] = phi[i];
	// 		// p.P[i] = _SRC[i];
			
	// 	}
	// }

	// Display computational time
	_timer = clock() - _timer;    // Calculate the time duration
	printf("<-> Pressure calculation comp. time:   [%f s]\n", (double)_timer/CLOCKS_PER_SEC);

	
	// =================================================================================
	// ----------------------- Body Surface Pressure Coefficient -----------------------
	// =================================================================================
	if (cp_calc == true){
		// Calculation of the pressure coefficient
		int _num = 48;		// Number of interpolated point
		std::vector<double> alpha(_num,0.0);
		Particle surface, buffer;
		
		// The interpolated point at body surface
		surface.num = _num;
		surface.x.resize(_num,0.0);
		surface.y.resize(_num,0.0);
		surface.s.resize(_num,Pars::sigma);
		surface.P.resize(_num,0.0);
		for (int i = 0; i < _num; i++){
			surface.x[i] = -(0.5+0*Pars::sigma)*std::cos((i*2.0*M_PI)/_num);
			surface.y[i] = (0.5+0*Pars::sigma)*std::sin((i*2.0*M_PI)/_num);
		}
		
		// The set of source point (interpolating data)
		buffer.num = 0;
		buffer.x.clear();
		buffer.y.clear();
		buffer.s.clear();
		buffer.P.clear();
		double lim_max = 0.5 + Pars::sigma * 7;
		double lim_min = -0.5 - Pars::sigma * 7;
		for (int i = 0; i < p.num; i++){
			if ((p.x[i] > lim_min && p.x[i] < lim_max) && (p.y[i] > lim_min && p.y[i] < lim_max)){
				buffer.num++;
				buffer.x.push_back(p.x[i]);
				buffer.y.push_back(p.y[i]);
				buffer.s.push_back(p.s[i]);
				buffer.P.push_back(p.P[i]);
			}
		}
		
		// Find the neighbor for calculation use inter search based on spatial hashing
		std::vector<std::vector<int>> _par2grdNeighbor,_grd2grdNeighbor;
		InterSearchNgh _interSearch;
		_interSearch.find_neighbor(buffer, surface, _par2grdNeighbor);
		_interSearch.find_neighbor(surface, surface, _grd2grdNeighbor);

		// Interpolation of pressure data
		LSMPSb lsmps_p;
		lsmps_p.set_LSMPS(surface.x, surface.y, surface.s, surface.P,
						buffer.x, buffer.y, buffer.s, buffer.P,
						_par2grdNeighbor, _grd2grdNeighbor);
		surface.P = lsmps_p.get_d0();

		// Save the data
		std::ofstream _data;
		std::string name = "output/Cp_";
		name += std::to_string(this->counter1);
		name += ".csv";
		this->counter1 ++;
		_data.open(name);
		_data << "x,y,Cp\n";
		
		double p_0 = 0.0e0;		// freestream pressure
		
		for (int i = 0; i < surface.num; i++){
			_data << surface.x[i] << ","
				<< surface.y[i] << ","
				<< 2 * (surface.P[i] - p_0) / (Pars::RHO * Pars::U_inf * Pars::U_inf) << "\n";
		}
		_data.close();
	}

	return;
}

void pressure_poisson::save_pressure_source(const std::vector<double>&x, const std::vector<double>&y, const std::vector<double>&src, const std::vector<bool>&mark){
	std::ofstream _data;
	std::string name = "output/pressureSRC";
	name += std::to_string(this->counter);
	name += ".csv";
	counter ++;
	_data.open(name);
	_data << "x,y,SRC,mark\n";
	
	for (size_t i = 0; i < x.size(); i++){
		_data << x[i] << ","
		      << y[i] << ","
			  << src[i] << ","
			  << mark[i] << "\n";
	}
	_data.close();
	return;
}