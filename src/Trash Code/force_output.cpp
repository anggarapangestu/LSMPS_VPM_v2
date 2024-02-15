#include "force_calc.hpp"

void force_calculation::output(int it, Particle &p, Body &b, double &cum_time)
{
	// -- accessing struct data, variables that are not commented are an input only
	int &np = p.num;
	int &nb = b.n_node;

	std::vector<double> &xp = p.x;
	std::vector<double> &yp = p.y;
	std::vector<double> &gpz = p.gz;
	std::vector<double> &sp = p.s;
	std::vector<double> &up = p.u;
	std::vector<double> &vp = p.v;
	std::vector<double> &wp = p.vorticity;

	std::vector<double> &xb = b.x;
	std::vector<double> &yb = b.y;

	// defining output variables
	std::ofstream ofs;

	// TODO: save Cx_Cy_Cz_Forces Linear Impulse - every saved timestep
	force_linear_impulse(it, np, xp, yp, gpz, sp);		// Uncomment for saving the data <==============

	// TODO: save Isumx,Isumy in case of continue the stopped running
	if (it == 0)
	{
		ofs.open("output/Isum.csv");
		ofs << "time,Isumx,Isumy\n";
		ofs.close();
	}
	ofs.open("output/Isum.csv", std::ofstream::out | std::ofstream::app);
	ofs << "," << it * Pars::dt << "," << this->_Isumx[it] << "," << this->_Isumy[it] << "\n";
	ofs.close();

} // end of function


void force_calculation::force_linear_impulse(int it, int np, const std::vector<double> &xp, const std::vector<double> &yp,
                                     const std::vector<double> &gpz, const std::vector<double> &sp)
{
    // internal variables
    double force_x_lIm, force_y_lIm, Cx_lIm, Cy_lIm;
    double A_ref[2];
    // A_ref[0] = Pars::ly;
    // A_ref[1] = Pars::lx;
    // =========== Linear Impulse ===============================
    this->_Isumx.push_back(0.0e0); // accessed from (a) (Life time is over all simulation)
    this->_Isumy.push_back(0.0e0); // accessed from (b)
    force_x_lIm = 0.0;
    force_y_lIm = 0.0;
    for (int i = 0; i < np; i++)
    {
        //ganti dikit area = sp[i]^2
        double _area = std::pow(sp[i] / Pars::sigma, 2);        // The nondimensional area toward the smallest particle size
        //_area = pow(sp[i], 2);
        //gpz -> vorticity
        //double vort = gpz[i] / pow(sp[i],2);

        this->_Isumx[it] += yp[i] * gpz[i] * _area;
        this->_Isumy[it] += -xp[i] * gpz[i] * _area; 
    }

    if ((it == 0) || (it == 1)) // because c++ indexing started from '0'
    {
        force_x_lIm = 0.0e0;
        force_y_lIm = 0.0e0;
        Cx_lIm = 0.0e0;
        Cy_lIm = 0.0e0;
    }
    else
    {
        force_x_lIm = -Pars::RHO * (this->_Isumx[it] - this->_Isumx[it - 2]) / (2.0e0 * Pars::dt);
        force_y_lIm = -Pars::RHO * (this->_Isumy[it] - this->_Isumy[it - 2]) / (2.0e0 * Pars::dt);
        //force_x_lIm = -Pars::dens * (this->_Isumx[it] - this->_Isumx[it - 1]) / (Pars::dt);
        //force_y_lIm = -Pars::dens * (this->_Isumy[it] - this->_Isumy[it - 1]) / (Pars::dt);
        // Cx_lIm of time step "it - 1", dens/dens = 1
        // Cx_lIm = -(Isumx[it] - Isumx[it - 2]) / (Pars::dt * std::pow(abs(Pars::Uf), 2) * A_ref[0]);
        // Cy_lIm = -(Isumy[it] - Isumy[it - 2]) / (Pars::dt * std::pow(abs(Pars::Uf), 2) * A_ref[1]);
        Cx_lIm = -(this->_Isumx[it] - this->_Isumx[it - 2])/ (Pars::dt * std::pow(abs(Pars::U_inf), 2));
        Cy_lIm = -(this->_Isumy[it] - this->_Isumy[it - 2]) / (Pars::dt * std::pow(abs(Pars::U_inf), 2));

    }


    // saving data
    std::ofstream ofs;
    if (it == 0)
    {
        ofs.open("output/force_Cd_lIm.csv");
        /*ofs << w20 << "time" << w20 << "force_x_lIm" << w20 << "force_y_lIm" << w20 << "Cx_lIm" << w20 << "Cy_lIm\n";
        ofs << w20 << it * Pars::dt << w20 << force_x_lIm << w20 << force_y_lIm << w20 << Cx_lIm << w20 << Cy_lIm << "\n";
        ofs.close();*/
        ofs << "" << "time" << "," << "force_x_lIm" << "," << "force_y_lIm" << "," << "Cx_lIm" << "," << "Cy_lIm\n";
        ofs << "" << it * Pars::dt << "," << force_x_lIm << "," << force_y_lIm << "," << Cx_lIm << "," << Cy_lIm << "\n";
        ofs.close();
    }
    else if (it >= 1)
    {
        /*ofs.open("output/force_Cd_lIm.dat", std::ofstream::out | std::ofstream::app);
        ofs << w20 << it * Pars::dt << w20 << force_x_lIm << w20 << force_y_lIm << w20 << Cx_lIm << w20 << Cy_lIm << "\n";
        ofs.close();*/
        ofs.open("output/force_Cd_lIm.csv", std::ofstream::out | std::ofstream::app);
        ofs << "" << it * Pars::dt << "," << force_x_lIm << "," << force_y_lIm << "," << Cx_lIm << "," << Cy_lIm << "\n";
        ofs.close();
    }
} // end of function

// void base_force_calculation::force_pen(int it, Particle p, double x_loc[2])
// {
// 	// internal variables
// 	int nPi = p.num;
// 	std::vector<double> xPi = p.x;
// 	std::vector<double> yPi = p.y;	
// 	std::vector<double> uPi = p.u;
// 	std::vector<double> vPi = p.v;
// 	std::vector<double> sPi = p.s; 
// 	std::vector<double> kai = p.chi;

// 	std::vector<double> uSi = p.u; 	// Need modification later on
// 	std::vector<double> vSi = p.v; 	// Need modification later on

// 	for (int i = 0; i < nPi; i++){ 	// Must be deleted after modification later
// 		uSi[i] = 0.0;
// 		vSi[i] = 0.0;
// 	}

// 	double fx, fy, force_x_pen, force_y_pen, C_x_pen, C_y_pen, C_m_pen, Moment_pen;
// 	fx = 0.0e0;
// 	fy = 0.0e0;
// 	force_x_pen = 0.0e0;
// 	force_y_pen = 0.0e0;
// 	Moment_pen =0.0e0;
// 	//x_loc[0] = 0.0;
// 	//x_loc[1] = 0.0;
// 	// Using Explicit scheme, Rasmussen 2011
// 	// lambda = 1.0e0/dt; // for Explicit scheme, Rasmussen 2011
// 	// If the solid body is considered rigid, all of uSi and vSi are the same. 
// 	if(Pars::opt_pen_iter == 1)
// 	{
// 		for (int i = 0; i < nPi; i++)
// 		{
// 			// Fluid to solid "+"

// 			if (Pars::opt_pen == 1 || Pars::opt_pen == 2){
// 				fx = -Pars::RHO * Pars::lambda * kai[i] * (-uPi[i] + uSi[i]) * std::pow(sPi[i], 2);
// 				fy = -Pars::RHO * Pars::lambda * kai[i] * (-vPi[i] + vSi[i]) * std::pow(sPi[i], 2);
// 				force_x_pen += fx; 
// 				force_y_pen += fy; 
// 			}else if (Pars::opt_pen == 3){
// 				fx = -Pars::RHO * kai[i] * ((-uPi[i] + uSi[i]) / Pars::dt) * std::pow(sPi[i], 2);
// 				fy = -Pars::RHO * kai[i] * ((-vPi[i] + vSi[i]) / Pars::dt) * std::pow(sPi[i], 2);
// 				force_x_pen += fx; 
// 				force_y_pen += fy; 
// 			}
// 			// Hitung Moment.
// 			if(fy > 1.0e-12){
// 				Moment_pen += (fy * (x_loc[0] - xPi[i])) ;
// 			}
// 			if(fx > 1.0e-12){
// 				Moment_pen += (fx * -(x_loc[1] - yPi[i]));
// 			}
// 		}

// 		//Untuk EOM vibration
// 		Pars::gaya = force_y_pen;
// 		Pars::momen = Moment_pen;

// 		// Coefficient of Forces
// 		C_x_pen = force_x_pen / (0.5 * Pars::RHO * std::pow(std::abs(Pars::U_inf), 2) * Pars::ly);
// 		C_y_pen = force_y_pen / (0.5 * Pars::RHO * std::pow(std::abs(Pars::U_inf), 2) * Pars::lx);
// 		C_m_pen = Moment_pen / (0.5 * Pars::RHO * std::pow(std::abs(Pars::U_inf), 2) * Pars::lx * Pars::lx); // !!kalau airfoil harus diubah jadi CHORD lx nya
// 	}

// 	if (Pars::opt_pen_iter == 2){
// 		for (int i = 0; i < nPi; i++)
// 		{
// 			// Fluid to solid "+", alpha == 2
// 			//fx = -Pars::RHO * 2 * kai[i] * ((-uPi[i] + uSi[i])/Pars::dt) * std::pow(sPi[i], 2) ;
// 			//fy = -Pars::RHO * 2 * kai[i] * ((-vPi[i] + vSi[i])/Pars::dt)  * std::pow(sPi[i], 2) ;
// 			fx = -Pars::RHO * uPi[i] * std::pow(sPi[i], 2) ;
// 			fy = -Pars::RHO * vPi[i] * std::pow(sPi[i], 2) ;
// 			force_x_pen += fx; // ! should be changed for multiresolution
// 			force_y_pen += fy; // ! should be changed for multiresolution
		
// 			// Hitung Moment.
// 			if(fy > 1.0e-12){
// 				Moment_pen += (fy * (x_loc[0] - xPi[i])) ;
// 			}
// 			if(fx > 1.0e-12){
// 				Moment_pen +=  (fx * -( x_loc[1] - yPi[i]  ));
// 			}
// 		}

// 		C_x_pen = force_x_pen / (0.5 * Pars::RHO * std::pow(std::abs(Pars::U_inf), 2) * Pars::ly);
// 		C_y_pen = force_y_pen / (0.5 * Pars::RHO * std::pow(std::abs(Pars::U_inf), 2) * Pars::lx);
// 		C_m_pen = Moment_pen / (0.5 * Pars::RHO * std::pow(std::abs(Pars::U_inf), 2) * Pars::lx * Pars::lx); // !!kalau airfoil harus diubah jadi CHORD lx nya
// 	}

// 	// defining output variables
// 	std::ofstream ofs;
// 	if (it == 0)
// 	{
// 		ofs.open("output/force_data_penalization.csv");
// 		ofs << "" << "time" << "," << "Force_x_pen" << "," << "Force_y_pen" << "," << "Moment_pen" << "," << "C_x_pen" << "," << "C_y_pen"  << "," << "C_M_pen\n";
// 		ofs.close();
// 	}

// 	ofs.open("output/force_data_penalization.csv", std::ofstream::out | std::ofstream::app);
// 	ofs << "" << it * Pars::dt << "," << force_x_pen << "," << force_y_pen << "," << Moment_pen << "," << C_x_pen << "," << C_y_pen << "," << C_m_pen << "\n";
// 	ofs.close();
// }
