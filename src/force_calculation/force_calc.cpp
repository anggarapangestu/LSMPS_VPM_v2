#include "force_calc.hpp"

double force_calculation::FD_diff1(const double h, const std::vector<double> &f, int order){
	double diff;
	if (order == 1){
		diff = 1 / h * (-f[0] + f[1]);
	}else if (order == 2){
		diff = 1 / (2.0 * h) * (-3*f[0] + 4*f[1] - f[2]);
	}else if(order == 3){
		diff = 1 / (6.0 * h) * (-11*f[0] + 18*f[1] - 9*f[2] + 2*f[3]);
	}
	return diff;
}

// *******************************************************************
// *********************** Force Calc Manager ************************
// *******************************************************************
// The force calculation manager
void force_calculation::force_calc(const Particle& par, const Body& body, int step, int _type){
    std::cout << "\nSaving force data ...\n";
    clock_t t = clock();
    // Assign the step time
    this->iter = step;
    
    // Change the initial condition
    if (step == 0){
        this->init = true;
    }
    if (Pars::opt_start_state == 1){
        if(step == Pars::resume_step){
            this->init = true;
        }
    }

    // ================ Barrier Mark =================
    // [!] TODO: Saving Data Force 
    if (_type == 1/*Pars::opt_force_type == 1*/){
        std::cout << "<+> Direct Calculation Force\n";
        this->direct_force(par, body);
    }
    else if (_type == 2/*Pars::opt_force_type == 2*/){
        std::cout << "<+> Penalization Based Force\n";
        double xpus[2] = {0, 0};
        this->pen_force(par, body);
    }
    else if (Pars::opt_force_type == 3){
        std::cout << "<+> Impulse Based Force\n";
        this->Force2(step,1,2,3,4,1,2,par);
    }

    // Force calculation summary time display
    t = clock() - t;
    printf("<-> Force calculation computation time [%f s]\n", (double)t/CLOCKS_PER_SEC);

    return;
}

// *******************************************************************
// ************************ Force Calculation ************************
// *******************************************************************
// Direct force calculation
void force_calculation::direct_force(const Particle& par, const Body& body){
    // Start the simulation
	/* Procedure
	   [1] Calculate the pressure at each panel node
	       * Pressure: Calculate at each surface panel mid point -> LSMPS B
	   [2] Calculate the shear stress at each panel node
           * Velocity: Calculate at 3-4 point outward from each panel midpoint -> LSMPS B
           * Shear Stress: Calculate from the interpolated point -> Finite Difference
	   [3] Integration of force, then non-dimensionalize the force
	*/

	// [PROCEDURE 0] : Define the particle data near solid body
    // *************
    Particle particle;
    particle.num = 0;
    for (int i = 0; i < par.num; i++){
        // if((par.x[i] > body.min_pos[0] - 20*Pars::sigma) && (par.x[i] < body.max_pos[0] + 20*Pars::sigma) &&
        //    (par.y[i] > body.min_pos[1] - 20*Pars::sigma) && (par.y[i] < body.max_pos[1] + 20*Pars::sigma))
		// if(par.isNearBody[i] == true)
		if(par.bodyPart[i] != -1)
		{
            particle.x.push_back(par.x[i]);
            particle.y.push_back(par.y[i]);
            particle.s.push_back(par.s[i]);
            particle.P.push_back(par.P[i]);
            particle.u.push_back(par.u[i]);
            particle.v.push_back(par.v[i]);
            particle.num++;
        }
    }

    // // [0] Debug for checking
    // // Save the data
    // particle.gz.resize(particle.num,0.0);
    // particle.vorticity.resize(particle.num,0.0);
    // particle.chi.resize(particle.num,0.0);
    // particle.isActive.resize(particle.num,false);
	// save_step.save_state(particle,"particle",0);

    // [PROCEDURE 1] : Calculate the pressure at each body panel surface
	// *************
	// Store the data of panel midpoint to temporary variable
	Particle _panel;
	
	// Reassign the panel midpoint data
	_panel.num = body.n_panel;
	_panel.x = body.x_m;
	_panel.y = body.y_m;
	_panel.s.resize(_panel.num, Pars::sigma);

	// # PRESSURE CALCULATION
	// Evaluation of the neighbor at panel data
	InterSearchNgh _interSearch;
	_interSearch.find_neighbor(particle, _panel, _panel.neighbor);		// The panel neighbor based on particle source

	// Pressure at each panel mid point
	this->interpolation.set_LSMPS(_panel.x, _panel.y, _panel.s, _panel.P,           // Target
				                  particle.x, particle.y, particle.s, particle.P,   // Source
				                  _panel.neighbor, _panel.neighbor);
	_panel.P = this->interpolation.get_d0();

	// // [0] Debug for checking
    // // Save the data
    // _panel.gz.resize(_panel.num,0.0);
    // _panel.vorticity.resize(_panel.num,0.0);
    // _panel.chi.resize(_panel.num,0.0);
    // _panel.isActive.resize(_panel.num,false);
    // _panel.u.resize(_panel.num, 0.0);
	// _panel.v.resize(_panel.num, 0.0);
	// save_step.save_state(_panel,"panel",0);

	// [PROCEDURE 2] : Shear Stress calculation
	// *************
    // Support Interpolation Point
	// Internal variable
	double* x_piv = new double[DIM];
	double* x_normal = new double[DIM];
	std::vector<double> _u,_v,_U;
	const int _order = 3;

	// The shear stress data
    std::vector<double> dUdn(_panel.num,0.0);		// The value of dU_t/dn at each panel

	Particle cluster;		// Collect temporary interpolation particle at each panel
    // Particle allInter;
	cluster.x.resize(_order+1);
	cluster.y.resize(_order+1);
	cluster.s.resize(_order+1);
    
    // Interpolation size
	double int_size = Pars::sigma*1.0;

	// Interpolate the velocity into the interpolation point then calculate the first differential 
    for (int i = 0; i < _panel.num; i++){
		// Determine the panel normal position and mid position
		x_piv[0] = body.x_m[i];
		x_piv[1] = body.y_m[i];
		x_normal[0] = body.x_n[i];
		x_normal[1] = body.y_n[i];
		
		// Define a new interpolation particle
		for (int j = 0; j <= _order; j++){
            // Note: The position is adjusted 0.5 times of the particle size outside to counter the penalization effect
			cluster.x[j] = (x_piv[0] + (j+1.5)*(int_size)*x_normal[0]);
			cluster.y[j] = (x_piv[1] + (j+1.5)*(int_size)*x_normal[1]);
			cluster.s[j] = Pars::sigma;

			// // Debug
			// allInter.num++;
			// allInter.x.push_back(cluster.x[j]);
			// allInter.y.push_back(cluster.y[j]);
            // allInter.s.push_back(Pars::sigma);
		}

		// Neighbor evaluation
		cluster.neighbor.clear();
        _interSearch.find_neighbor(particle, cluster, cluster.neighbor);		// The _panel neighbor
		
		// Velocity interpolation
		_u.clear();_v.clear();
		this->interpolation.set_LSMPS(cluster.x, cluster.y, cluster.s, cluster.u, 
				   particle.x, particle.y, particle.s, particle.u, 
				   cluster.neighbor, cluster.neighbor);
		_u = this->interpolation.get_d0();
		
		this->interpolation.set_LSMPS(cluster.x, cluster.y, cluster.s, cluster.v, 
				   particle.x, particle.y, particle.s, particle.v, 
				   cluster.neighbor, cluster.neighbor);
		_v = this->interpolation.get_d0();

		// Calculation of tangential velocity
		_U.clear();
		_U.resize(_order + 1, 0.0);
		for (int j = 0; j <= _order; j++){
			// Tangential velocity U . n_t (n_y,-n_x)
            _U[j] = _u[j] * x_normal[1] - _v[j] * x_normal[0];

			// DEBUGGING DATA
			// std::cout << "PRINT U AND V : " << _u[j] << " , " << _v[j] << "\n";
			// allInter.u.push_back(_u[j]);
			// allInter.v.push_back(_v[j]);
		}

		// Velocity differentiation
		dUdn[i] = FD_diff1(int_size,_U,_order);
	}
	delete [] x_piv, x_normal;

	// // [0] Debug for checking
    // // Save the data
    // allInter.gz.resize(allInter.num,0.0);
	// allInter.vorticity.resize(allInter.num,0.0);
    // allInter.chi.resize(allInter.num,0.0);
    // allInter.P.resize(allInter.num,0.0);
    // allInter.isActive.resize(allInter.num,false);
	// this->save_step.save_state(allInter,"AllInter",0);


	// [PROCEDURE 3] : Force integration
	// ***********
	// Internal variable
	double AoA = 0.0;	                    // In degree
	double Fx = 0.0, Fy = 0.0, Mom = 0.0;   // Dimensional
    double Cq = 0.5 * Pars::RHO * Pars::U_inf * Pars::U_inf * Pars::Df;
	double Lift, Drag, Cd, Cl, Cm;          // Non-dimensional
	
	// Integration
    double _Fx, _Fy;
	for (size_t i = 0; i < _panel.num; i++){
		// ======== [1] ========
        // Contribution of pressure
		_Fx = - _panel.P[i] * body.x_n[i] * body.size[i];
		_Fy = - _panel.P[i] * body.y_n[i] * body.size[i];

        // Accumulate
        Fx += _Fx;
        Fy += _Fy;
        Mom += (body.y_m[i] - Pars::ycenter)*(-_Fx) + (body.x_m[i] - Pars::xcenter)*(_Fy);
		
		// ======== [2] ========
        // Contribution of shear stress
		_Fx = Pars::MU * dUdn[i] *  (body.y_n[i]) * body.size[i];  // Fx = (F_t . t_x)
		_Fy = Pars::MU * dUdn[i] * (-body.x_n[i]) * body.size[i];  // Fy = (F_t . t_y)

        // Accumulate
        Fx += _Fx;
        Fy += _Fy;
        Mom += (body.y_m[i] - Pars::ycenter)*(-_Fx) + (body.x_m[i] - Pars::xcenter)*(_Fy);
	}

	// Force and moment coefficient calculation
    Lift = Fy * cos(AoA * M_PI / 180) - Fx * sin(AoA * M_PI / 180);
	Drag = Fx * cos(AoA * M_PI / 180) + Fy * sin(AoA * M_PI / 180);
	Cl = Lift / Cq;
	Cd = Drag / Cq;
    Cm = Mom / Cq / Pars::Df;

    // =========== DATA SAVING ===========
	std::ofstream ofs;
    // Save header
	if (this->init == true){
		ofs.open("output/force_data_direct.csv");
		ofs << "time_dir" << "," 
            << "Fx" << "," << "Fy" << "," << "M"  << "," 
            << "Cx" << "," << "Cy" << "," << "Cm" << "\n";
		ofs.close();
        this->init = false;
	}

    // Save data
	ofs.open("output/force_data_direct.csv", std::ofstream::out | std::ofstream::app);
	ofs << this->iter * Pars::dt << "," 
        << Drag << "," << Lift << "," << Mom << "," 
        << Cd   << "," << Cl   << "," << Cm  << "\n";
	ofs.close();

    return;
}

// Penalization force calculation
void force_calculation::pen_force(const Particle& par, const Body& body){
	// [PROCEDURE 0] : Define the particle data near solid body
    // *************
    Particle p;
    p.num = 0;
    for (int i = 0; i < par.num; i++){
        // if(par.isNearBody[i] == true)
		// if((par.x[i] > body.min_pos[0] - 20*Pars::sigma) && (par.x[i] < body.max_pos[0] + 20*Pars::sigma) &&
        //    (par.y[i] > body.min_pos[1] - 20*Pars::sigma) && (par.y[i] < body.max_pos[1] + 20*Pars::sigma))
		if(par.bodyPart[i] != -1)
		{
            p.x.push_back(par.x[i]);
            p.y.push_back(par.y[i]);
            p.s.push_back(par.s[i]);
            p.chi.push_back(par.chi[i]);
            p.u.push_back(par.u[i]);
            p.v.push_back(par.v[i]);
            p.num++;
        }
    }

    // [PROCEDURE 1] : Calculate the force data
    // *************
    // Internal variables
    double fx, fy, F_x, F_y, Cx, Cy, Cm, Mom;
    double Cq = 0.5 * Pars::RHO * Pars::U_inf * Pars::U_inf * Pars::Df;
	fx = 0.0e0;
	fy = 0.0e0;
	F_x = 0.0e0;
	F_y = 0.0e0;
	Mom = 0.0e0;

	std::vector<double> uSi = p.u; 	// Need modification later on
	std::vector<double> vSi = p.v; 	// Need modification later on

	for (int i = 0; i < p.num; i++){ 	// Must be deleted after modification later
		uSi[i] = 0.0;
		vSi[i] = 0.0;
	}

	// Calculation of the data force
	// Non Iterative Penalization
    if(Pars::opt_pen_iter == 1)
	{
		for (int i = 0; i < p.num; i++)
		{
			// Force calculation
            if (Pars::opt_pen == 1 || Pars::opt_pen == 2){  // Implicit
				fx = -Pars::RHO * Pars::lambda * p.chi[i] * (-p.u[i] + uSi[i]) * std::pow(p.s[i], 2);
				fy = -Pars::RHO * Pars::lambda * p.chi[i] * (-p.v[i] + vSi[i]) * std::pow(p.s[i], 2);
				F_x += fx; 
				F_y += fy; 
			}else if (Pars::opt_pen == 3){                  // Explicit
				fx = -Pars::RHO * p.chi[i] * ((-p.u[i] + uSi[i]) / Pars::dt) * std::pow(p.s[i], 2);
				fy = -Pars::RHO * p.chi[i] * ((-p.v[i] + vSi[i]) / Pars::dt) * std::pow(p.s[i], 2);
				F_x += fx; 
				F_y += fy; 
			}

			// Moment calculation
			if(fy > 1.0e-12){
				Mom += (fy * (Pars::xcenter - p.x[i])) ;
			}
			if(fx > 1.0e-12){
				Mom += (fx * -(Pars::ycenter - p.y[i]));
			}
		}

		// //Untuk EOM vibration
		// Pars::gaya = F_y;
		// Pars::momen = Mom;

		// Coefficient of Forces
		Cx = F_x / Cq;
		Cy = F_y / Cq;
		Cm = Mom / Cq / Pars::Df; // For airfoil change the chord to be lx
	}
    // Iterative Penalization [Not Finished]
	if (Pars::opt_pen_iter == 2){
		for (int i = 0; i < p.num; i++)
		{
			// Fluid to solid "+", alpha == 2
			//fx = -Pars::RHO * 2 * p.chi[i] * ((-p.u[i] + uSi[i])/Pars::dt) * std::pow(p.s[i], 2) ;
			//fy = -Pars::RHO * 2 * p.chi[i] * ((-p.v[i] + vSi[i])/Pars::dt)  * std::pow(p.s[i], 2) ;
			fx = -Pars::RHO * p.u[i] * std::pow(p.s[i], 2) ;
			fy = -Pars::RHO * p.v[i] * std::pow(p.s[i], 2) ;
			F_x += fx; // ! should be changed for multiresolution
			F_y += fy; // ! should be changed for multiresolution
		
			// Hitung Moment.
			if(fy > 1.0e-12){
				Mom += (fy * (Pars::xcenter - p.x[i]));
			}
			if(fx > 1.0e-12){
				Mom +=  (fx * -(Pars::ycenter - p.y[i]));
			}
		}

		// Coefficient of Forces
		Cx = F_x / Cq;
		Cy = F_y / Cq;
		Cm = Mom / Cq / Pars::Df; // For airfoil change the chord to be lx
	}

	// =========== DATA SAVING ===========
	std::ofstream ofs;
    // Data header
	if (this->init == true){
		ofs.open("output/force_data_penalization.csv");
		ofs << "time_pen" << "," 
            << "Fx" << "," << "Fy" << "," << "M"  << "," 
            << "Cx" << "," << "Cy" << "," << "Cm" <<"\n";
		ofs.close();
        this->init = false;
	}

	// Data input
    ofs.open("output/force_data_penalization.csv", std::ofstream::out | std::ofstream::app);
	ofs << this->iter * Pars::dt << "," 
        << F_x << "," << F_y << "," << Mom << "," 
        << Cx  << "," << Cy  << "," << Cm  << "\n";
	ofs.close();
}

// Impulse Force Calculation
void force_calculation::imp_force(const Particle& par, const Body& body){
    // Under code migration
    return;
}