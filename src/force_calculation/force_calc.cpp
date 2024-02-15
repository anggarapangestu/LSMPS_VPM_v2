#include "force_calc.hpp"

/** The differential order used in calculating force using vorticity moment integral
 *   1:= The differential order of 1
 *   2:= The differential order of 2
*/
#define VORTICITY_MOMENT_DIFF_ORDER 2

// =====================================================
// +----------------- Utility Function ----------------+
// =====================================================

double force_calculation::FD_diff1(const double h, const std::vector<double> &f, int order){
	double diff = 0.0;
	if (order == 1){
		diff = 1 / h * (-f[0] + f[1]);
	}else if (order == 2){
		diff = 1 / (2.0 * h) * (-3*f[0] + 4*f[1] - f[2]);
	}else if(order == 3){
		diff = 1 / (6.0 * h) * (-11*f[0] + 18*f[1] - 9*f[2] + 2*f[3]);
	}
	return diff;
}


// =====================================================
// +------------ Force Calculation Manager ------------+
// =====================================================

/**
 *  @brief  The force calculation manager.
 *         
 *  @param  _particle  Particle data container.
 *  @param  _bodyList  The list of body data container used for force calculation.
 *  @param  _step  Current iteration simulation.
 *  @param  _type  Type of force calculation.
*/
void force_calculation::force_calc(const Particle& par, 
								   const std::vector<Body> &_bodyList, 
								   int _step, int _type, std::string name)
{
	// Exception 1 (No force is calculated)
	if (_type == 0) return;

	// Exception 2 (No body for force calculation)
	if (N_BODY == 0) return;

	// Log prompt
	std::cout << "\nSaving force data ...\n";

    // Assign the step time
    this->iter = _step;

	// Computation timer
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

	// // Change the initial condition (<!> No need to be changed <!>)
	// if (Pars::opt_start_state == 1) this->init = false;		// Continue from the previous data

    // Saving the force data according to the user choice type
	switch (_type){
	case 1:
		printf("%s<+> Direct Force Calculation%s\n", FONT_CYAN, FONT_RESET);
		for (int part = 0; part < N_BODY; part++)
			this->direct_force(par, _bodyList[part], part);
		break;
	case 2:
		printf("%s<+> Penalization Based Force%s\n", FONT_CYAN, FONT_RESET);
        this->pen_force(par, _bodyList, name);
		break;
	case 3:
		// This method is not working (It should be worked as in theory, but why <?>)
		printf("%s<+> Vorticity Moment Integral Based Force%s\n", FONT_CYAN, FONT_RESET);
		this->int_vor_force(par, _bodyList, name);
		break;
	case 4:
		printf("%s<+> NOCA Impulse Based Force%s\n", FONT_CYAN, FONT_RESET);
        this->Force2(_step,1,2,3,4,1,2,par);	// Still not working
		break;
	default:
		break;
	}
    
    // Force calculation summary time display
	#if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    printf("<-> Force calculation computation time [%f s]\n", _time);

    return;
}

// // The force calculation manager
// void force_calculation::force_calc(const Particle& par, const Body& body, int step, int _type){
//     std::cout << "\nSaving force data ...\n";
//     clock_t t = clock();
//     // Assign the step time
//     this->iter = step;
    
//     // Change the initial condition
//     if (step == 0){
//         this->init = true;
//     }
//     if (Pars::opt_start_state == 1){
//         if(step == Pars::resume_step){
//             this->init = true;
//         }
//     }

//     // ================ Barrier Mark =================
//     // [!] TODO: Saving Data Force 
//     if (_type == 1/*Pars::opt_force_type == 1*/){
//         std::cout << "<+> Direct Calculation Force\n";
//         this->direct_force(par, body);
//     }
//     else if (_type == 2/*Pars::opt_force_type == 2*/){
//         std::cout << "<+> Penalization Based Force\n";
//         this->pen_force(par, body);
//     }
//     else if (Pars::opt_force_type == 3){
//         std::cout << "<+> Impulse Based Force\n";
//         this->Force2(step,1,2,3,4,1,2,par);
//     }

//     // Force calculation summary time display
//     t = clock() - t;
//     printf("<-> Force calculation computation time [%f s]\n", (double)t/CLOCKS_PER_SEC);

//     return;
// }


// =====================================================
// +------------ Force Calculation Method -------------+
// =====================================================
// Direct force calculation
void force_calculation::direct_force(const Particle& par, const Body& body, int bodyPart){
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
		if(par.bodyPart[i] == 1)
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
	for (int i = 0; i < _panel.num; i++){
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
	std::string fileName = "output/force_direct_body_" + std::to_string(bodyPart+1) + ".csv";
    
	// Save header
	if (this->init == true){
		ofs.open(fileName.c_str());
		ofs << "time_dir" << "," 
            << "Fx" << "," << "Fy" << "," << "M"  << "," 
            << "Cx" << "," << "Cy" << "," << "Cm" << "\n";
		ofs.close();
        this->init = false;
	}

    // Save data
	ofs.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);
	ofs << this->iter * Pars::dt << "," 
        << Drag << "," << Lift << "," << Mom << "," 
        << Cd   << "," << Cl   << "," << Cm  << "\n";
	ofs.close();

    return;
}

// Penalization force calculation
void force_calculation::pen_force(const Particle& p, const std::vector<Body> &_bodyList, std::string headName){
	// [PROCEDURE 0] : Data container initialization
    // *************
    // Internal variables (Store the data for each body)
	std::vector<double> Fx(N_BODY, 0.0);		// The force in X direction
	std::vector<double> Fy(N_BODY, 0.0);		// The force in Y direction
	std::vector<double> Fz(N_BODY, 0.0);		// The force in Z direction

	std::vector<double> Mx(N_BODY, 0.0);		// The moment in X direction (effected by fy & fz)
	std::vector<double> My(N_BODY, 0.0);		// The moment in Y direction (effected by fx & fz)
	std::vector<double> Mz(N_BODY, 0.0);		// The moment in Z direction (effected by fy & fx)
	
    // Temporary force calculation
	double fx, fy, fz;

	// The dynamic pressure
	double Cq = 0.5 * Pars::RHO * Pars::U_inf * Pars::U_inf;

	// The length or area reference <!> Need modification later on <!> 
	double Lref = Pars::Df;							// The reference length for 2D object
	double Aref = M_PI * Pars::Df * Pars::Df / 4;	// The reference area for 3D object (adjusting the shape) [Currently only suitable for sphere]

	// Body velocity <!> Need modification later on <!> 
	double uS = Pars::ubody;
	double vS = Pars::vbody;
	double wS = Pars::wbody;

	// Remember the equation for force 
	//   F = (1/2 * rho * U^2) * Cd * A_ref
	// The equation of force penalization
	//   F = -rho * integral over domain(lambda * chi * (Us - U) * dA)
	//      or
	//   F = -ρ \int{ λχ(Us-U)dA }

	// [PROCEDURE 1] : Calculate the force data
    // *************
	// Calculation of the data force
	#if (DIM == 2)
		// =========== Data calculation ===========
		for (int i = 0; i < p.num; i++){
			// Body part
			const int &part = p.bodyPart[i];

			// Exclude the calculation from the particle far from body
			if (part == -1) continue;

			// Force calculation
			// Intermediate variable
			const double _K = - Pars::RHO * Pars::lambda * p.chi[i];	// The head constant in penalization force calculation
			const double _A = p.s[i] * p.s[i];							// The area occupied by current particle
			
			// Calculate the force caused by current particle [Penalization method]
			fx = _K * (uS - p.u[i]) * _A;
			fy = _K * (vS - p.v[i]) * _A;

			// Assign to the force container
			Fx[part] += fx;
			Fy[part] += fy;

			// Moment calculation (M = r x F)
			Mz[part] += (fy * (p.x[i] - _bodyList[part].cen_pos[0])) - 
						(fx * (p.y[i] - _bodyList[part].cen_pos[1]));
		}
		
		// =========== Data Saving ===========
		std::ofstream ofs;
		for (int part = 0; part < N_BODY; part++){
			// Filename format
			std::string name = "output/" + headName + "_force_pen_body_" + std::to_string(part+1) + ".csv";
			
			// Data header
			if (this->init == true){
				ofs.open(name.c_str());
				ofs << "time" 
					<< "," << "Fx" << "," << "Fy" << "," << "M"  
					<< "," << "Cx" << "," << "Cy" << "," << "Cm" 
					<< "\n";
				ofs.close();
				this->init = false;
			}

			/** NOTE: The reference length must be modified for next calculation
			 *  > The drag regard of y length (thickness)
			 *  > The lift regard of x length
			 *  > But for cylinder (also square) has same x and y length
			*/

			// Data input
			ofs.open(name.c_str(), std::ofstream::out | std::ofstream::app);
			ofs << this->iter * Pars::dt 
				<< "," << Fx[part] 
				<< "," << Fy[part] 
				<< "," << Mz[part] 
				<< "," << Fx[part] / (Cq*Lref)			// Aerodynamic force coefficient in x
				<< "," << Fy[part] / (Cq*Lref)			// Aerodynamic force coefficient in y
				<< "," << Mz[part] / (Cq*Lref*Lref)		// Aerodynamic moment coefficient in z
				<< "\n";
			ofs.close();
		}
	#elif (DIM == 3)
		// =========== Data calculation ===========
		for (int i = 0; i < p.num; i++){
			// Body part
			const int &part = p.bodyPart[i];

			// Exclude the calculation from the particle far from body
			if (part == -1) continue;

			// Force calculation
			// Intermediate variable
			const double _K = - Pars::RHO * Pars::lambda * p.chi[i];	// The head constant in penalization force calculation
			const double _V = p.s[i] * p.s[i] * p.s[i];					// The volume occupied by current particle

			// Calculate the force caused by current particle [Penalization method]
			fx = _K * (uS - p.u[i]) * _V;
			fy = _K * (vS - p.v[i]) * _V;
			fz = _K * (wS - p.w[i]) * _V;

			// Assign to the force container
			Fx[part] += fx;
			Fy[part] += fy;
			Fz[part] += fz;

			// Moment calculation (M = r x F)
			Mx[part] += (fz * (p.y[i]-_bodyList[part].cen_pos[1])) - 
						(fy * (p.z[i]-_bodyList[part].cen_pos[2]));
			My[part] += (fx * (p.z[i]-_bodyList[part].cen_pos[2])) - 
						(fz * (p.x[i]-_bodyList[part].cen_pos[0]));
			Mz[part] += (fy * (p.x[i]-_bodyList[part].cen_pos[0])) - 
						(fx * (p.y[i]-_bodyList[part].cen_pos[1]));
		}
		
		// =========== Data Saving ===========
		std::ofstream ofs;
		for (int part = 0; part < N_BODY; part++){
			// Filename format
			std::string name = "output/" + headName + "_force_pen_body_" + std::to_string(part+1) + ".csv";
			// std::string name = "output/force_pen_body_" + std::to_string(part+1) + ".csv";
			
			// Data header
			if (this->init == true){
				ofs.open(name.c_str());
				ofs << "time" 
					<< "," << "Cx"
					<< "," << "Cy"
					<< "," << "Cz"
					<< "," << "Cmx"
					<< "," << "Cmy"
					<< "," << "Cmz"
					<< "\n";
				ofs.close();
				this->init = false;
			}

			// Data input
			ofs.open(name.c_str(), std::ofstream::out | std::ofstream::app);
			ofs << this->iter * Pars::dt 
				<< "," << Fx[part] / (Cq*Aref)		// Aerodynamic force coefficient in x
				<< "," << Fy[part] / (Cq*Aref)		// Aerodynamic force coefficient in y
				<< "," << Fz[part] / (Cq*Aref)		// Aerodynamic force coefficient in z
				<< "," << Mx[part] / (Cq*Aref*Lref)	// Aerodynamic moment coefficient in x
				<< "," << My[part] / (Cq*Aref*Lref)	// Aerodynamic moment coefficient in y
				<< "," << Mz[part] / (Cq*Aref*Lref)	// Aerodynamic moment coefficient in z
				<< "\n";
			ofs.close();
		}
	#endif
	
	return;
}

/*	// Penalization force calculation [OLD METHOD]
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
    double fx, fy, F_x, F_y, Mom, Cx, Cy, Cm;	// Temporary force calculation
	
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
            // if (Pars::opt_pen == 1 || Pars::opt_pen == 2){  // Implicit
				fx = -Pars::RHO * Pars::lambda * p.chi[i] * (-p.u[i] + uSi[i]) * std::pow(p.s[i], 2);
				fy = -Pars::RHO * Pars::lambda * p.chi[i] * (-p.v[i] + vSi[i]) * std::pow(p.s[i], 2);
				F_x += fx;
				F_y += fy;
			// }else if (Pars::opt_pen == 3){                  // Explicit
			// 	fx = -Pars::RHO * p.chi[i] * ((-p.u[i] + uSi[i]) / Pars::dt) * std::pow(p.s[i], 2);
			// 	fy = -Pars::RHO * p.chi[i] * ((-p.v[i] + vSi[i]) / Pars::dt) * std::pow(p.s[i], 2);
			// 	F_x += fx;
			// 	F_y += fy;
			// }

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
*/

// Integral of vorticity moment force calculation
void force_calculation::int_vor_force(const Particle& p, const std::vector<Body> &_bodyList, std::string headName){
	// Remember the equation for force 
	//   F = (1/2 * rho * U^2) * Cd * A_ref
	// The equation of force integral vorticity moment
	//   F = -rho * d/dt{ integral over domain(((x-x_cen) cross vorticity) * dA) }
	//      or
	//   F = -ρ d/dt{ \int{ (x-x_c)x ω dA } } -> 2D: Force per width [N/m], 3D: Force [N]
	
	//  Cross product of (I = r x ω)
	//   [Ix]   [x]   [ωx]   [y*ωz - z*ωy]   [ y*ωz]
	//   [Iy] = [y] x [ωy] = [z*ωx - x*ωz] = [-x*ωz]
	//   [Iz]   [z]   [ωz]   [x*ωy - y*ωx]   [  0  ]
	//                            3D           2D

	// [PROCEDURE 0] : Data container initialization
    // *************
    // Internal variables (Store the data for each body)
	std::vector<double> Fx(N_BODY, 0.0);		// The force in X direction
	std::vector<double> Fy(N_BODY, 0.0);		// The force in Y direction
	std::vector<double> Fz(N_BODY, 0.0);		// The force in Z direction
	double dx, dy, dz;				// Temporary distance variable
	double _IxdA, _IydA, _IzdA;		// Temporary moment element calculation
	double Cq = 0.5 * Pars::RHO * Pars::U_inf * Pars::U_inf;	// The dynamic pressure
	// The length or area reference <!> Need modification later on <!> 
	double Lref = Pars::Df;							// The reference length for 2D object
	double Aref = M_PI * Pars::Df * Pars::Df / 4;	// The reference area for 3D object (adjusting the shape) [Currently only suitable for sphere]

	// Reserve the moment data container (both for 2D and 3D)
	if (this->init == true){
		this->Ivortx = std::vector<std::vector<double>>(N_BODY, std::vector<double>(VORTICITY_MOMENT_DIFF_ORDER+1, 0.0));
		this->Ivorty = std::vector<std::vector<double>>(N_BODY, std::vector<double>(VORTICITY_MOMENT_DIFF_ORDER+1, 0.0));;
		this->Ivortz = std::vector<std::vector<double>>(N_BODY, std::vector<double>(VORTICITY_MOMENT_DIFF_ORDER+1, 0.0));;
	}
	
	// Rearrange the moment container (1 -> 0, and 2 -> 1)
	for (int part = 0; part < N_BODY; part++){
		this->Ivortx[part][0] = this->Ivortx[part][1];
		this->Ivorty[part][0] = this->Ivorty[part][1];
		this->Ivortz[part][0] = this->Ivortz[part][1];
	}
	#if (VORTICITY_MOMENT_DIFF_ORDER == 2)
		for (int part = 0; part < N_BODY; part++){
			this->Ivortx[part][1] = this->Ivortx[part][2];
			this->Ivorty[part][1] = this->Ivorty[part][2];
			this->Ivortz[part][1] = this->Ivortz[part][2];
		}
	#endif

	// [PROCEDURE 1] : Calculate the force data
    // *************
	// Calculation of the data force
	#if (DIM == 2)
		// =========== Vorticity Moment Calculation ===========
		for (int i = 0; i < p.num; i++){
			// Only calculate the active particle
			if (p.isActive[i] == false) continue;

			// Evaluate for each body
			for (int part = 0; part < N_BODY; part++){
				// Intermediate variable
				const double _A = p.s[i] * p.s[i];          // The area occupied by current particle

				// Calculate the distance from body center
				dx = p.x[i] - _bodyList[part].cen_pos[0];   // The x distance from body center
				dy = p.y[i] - _bodyList[part].cen_pos[1];   // The y distance from body center

				// Calculate the force caused by current particle [Penalization method]
				_IxdA =  dy * p.vorticity[i] * _A;
				_IydA = -dx * p.vorticity[i] * _A;

				// Assign to the moment container
				#if (VORTICITY_MOMENT_DIFF_ORDER == 1)
					this->Ivortx[part][1] += _IxdA;
					this->Ivorty[part][1] += _IydA;
				#elif (VORTICITY_MOMENT_DIFF_ORDER == 2)
					this->Ivortx[part][2] += _IxdA;
					this->Ivorty[part][2] += _IydA;
				#endif
			}
		}
		
		// =========== Force Data Calculation ===========
		if (this->iter < VORTICITY_MOMENT_DIFF_ORDER){
			for (int part = 0; part < N_BODY; part++){
				Fx[part] = 0;
				Fy[part] = 0;
			}
		}else{
			// Assign to the moment container
			#if (VORTICITY_MOMENT_DIFF_ORDER == 1)
				for (int part = 0; part < N_BODY; part++){
					Fx[part] = -Pars::RHO * (this->Ivortx[part][1] - this->Ivortx[part][0]) / (Pars::dt);
					Fy[part] = -Pars::RHO * (this->Ivorty[part][1] - this->Ivorty[part][0]) / (Pars::dt);
				}
			#elif (VORTICITY_MOMENT_DIFF_ORDER == 2)
				for (int part = 0; part < N_BODY; part++){
					Fx[part] = -Pars::RHO * (this->Ivortx[part][2] - this->Ivortx[part][0]) / (2*Pars::dt);
					Fy[part] = -Pars::RHO * (this->Ivorty[part][2] - this->Ivorty[part][0]) / (2*Pars::dt);
				}
			#endif	
		}
		
		// =========== Data Saving ===========
		std::ofstream ofs;
		for (int part = 0; part < N_BODY; part++){
			// Filename format
			std::string name = "output/" + headName + "_force_mom_body_" + std::to_string(part+1) + ".csv";
			
			// Data header
			if (this->init == true){
				ofs.open(name.c_str());
				ofs << "time" 
					<< "," << "Fx" << "," << "Fy"
					<< "," << "Cx" << "," << "Cy"
					<< "\n";
				ofs.close();
				this->init = false;
			}

			/** NOTE: The reference length must be modified for next calculation
			 *  > The drag regard of y length (thickness)
			 *  > The lift regard of x length
			 *  > But for cylinder (also square) has same x and y length
			*/

			// Data input
			ofs.open(name.c_str(), std::ofstream::out | std::ofstream::app);
			ofs << this->iter * Pars::dt 
				<< "," << Fx[part] 
				<< "," << Fy[part] 
				<< "," << Fx[part] / (Cq*Lref)	// Aerodynamic force coefficient in x
				<< "," << Fy[part] / (Cq*Lref)	// Aerodynamic force coefficient in y
				<< "\n";
			ofs.close();
		}
	#elif (DIM == 3)
		// =========== Vorticity Moment Calculation ===========
		for (int i = 0; i < p.num; i++){
			// Only calculate the active particle
			if (p.isActive[i] == false) continue;

			// Evaluate for each body
			for (int part = 0; part < N_BODY; part++){
				// Intermediate variable
				const double _V = p.s[i] * p.s[i] * p.s[i]; // The volume occupied by current particle

				// Calculate the distance from body center
				dx = p.x[i] - _bodyList[part].cen_pos[0];   // The x distance from body center
				dy = p.y[i] - _bodyList[part].cen_pos[1];   // The y distance from body center
				dz = p.z[i] - _bodyList[part].cen_pos[2];   // The z distance from body center

				// Calculate the force caused by current particle [Penalization method]
				_IxdA = (dy*p.vortz[i] - dz*p.vorty[i]) * _V;
				_IydA = (dz*p.vortx[i] - dx*p.vortz[i]) * _V;
				_IzdA = (dx*p.vorty[i] - dy*p.vortx[i]) * _V;

				// Assign to the moment container
				#if (VORTICITY_MOMENT_DIFF_ORDER == 1)
					this->Ivortx[part][1] += _IxdA;
					this->Ivorty[part][1] += _IydA;
					this->Ivortz[part][1] += _IzdA;
				#elif (VORTICITY_MOMENT_DIFF_ORDER == 2)
					this->Ivortx[part][2] += _IxdA;
					this->Ivorty[part][2] += _IydA;
					this->Ivortz[part][2] += _IzdA;
				#endif
			}
		}
		
		// =========== Force Data Calculation ===========
		if (this->iter < VORTICITY_MOMENT_DIFF_ORDER){
			for (int part = 0; part < N_BODY; part++){
				Fx[part] = 0;
				Fy[part] = 0;
				Fz[part] = 0;
			}
		}else{
			// Assign to the moment container
			#if (VORTICITY_MOMENT_DIFF_ORDER == 1)
				for (int part = 0; part < N_BODY; part++){
					Fx[part] = -Pars::RHO * (this->Ivortx[part][1] - this->Ivortx[part][0]) / (Pars::dt);
					Fy[part] = -Pars::RHO * (this->Ivorty[part][1] - this->Ivorty[part][0]) / (Pars::dt);
					Fz[part] = -Pars::RHO * (this->Ivortz[part][1] - this->Ivortz[part][0]) / (Pars::dt);
				}
			#elif (VORTICITY_MOMENT_DIFF_ORDER == 2)
				for (int part = 0; part < N_BODY; part++){
					Fx[part] = -Pars::RHO * (this->Ivortx[part][2] - this->Ivortx[part][0]) / (2*Pars::dt);
					Fy[part] = -Pars::RHO * (this->Ivorty[part][2] - this->Ivorty[part][0]) / (2*Pars::dt);
					Fz[part] = -Pars::RHO * (this->Ivortz[part][2] - this->Ivortz[part][0]) / (2*Pars::dt);
				}
			#endif	
		}
		
		// =========== Data Saving ===========
		std::ofstream ofs;
		for (int part = 0; part < N_BODY; part++){
			// Filename format
			std::string name = "output/" + headName + "_force_mom_body_" + std::to_string(part+1) + ".csv";
			
			// Data header
			if (this->init == true){
				ofs.open(name.c_str());
				ofs << "time" 
					<< "," << "Fx" << "," << "Fy" << "," << "Fz"
					<< "," << "Cx" << "," << "Cy" << "," << "Cz"
					<< "\n";
				ofs.close();
				this->init = false;
			}

			// Data input
			ofs.open(name.c_str(), std::ofstream::out | std::ofstream::app);
			ofs << this->iter * Pars::dt 
				<< "," << Fx[part] 
				<< "," << Fy[part] 
				<< "," << Fz[part] 
				<< "," << Fx[part] / (Cq*Aref)	// Aerodynamic force coefficient in x
				<< "," << Fy[part] / (Cq*Aref)	// Aerodynamic force coefficient in y
				<< "," << Fz[part] / (Cq*Aref)	// Aerodynamic force coefficient in y
				<< "\n";
			ofs.close();
		}
	#endif
	
	return;
}

// Impulse Force Calculation (NOCA method)
void force_calculation::imp_force(const Particle& par, const Body& body){
	// Remember the equation for force 
	//   F = (1/2 * rho * U^2) * Cd * A_ref
	// The equation of force integral vorticity moment
    return;
}