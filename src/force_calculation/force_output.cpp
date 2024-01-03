#include "force_calc.hpp"

// This subroutine to calculate Cx,Cy,Cz  forces direct from
// grid penalization data, not from extra grid data.
// Input: ngx,ngy of penalization,ux_gpen,uy_gpen,
//        uSx_gpen,uSy_gpen,sp_gpen,kai
//  Fs=-rho * itegral of (lambda*kai*(uS-u)dA), F_fluid=-F_solid
// where lambda = alfa/dt, in case of split step: lambda=1/dt
// dA = dV = sigma^3
// --> Force_x_pen(solid) = SUM(rho*kai/dt  * (u-uS)*sigma^3)
// but, Cd = Force_x/(1/2 * rho * U^2  * A ) where A: reference area, Sphere= pi*R^2 =pi*D^2/4  , FlatPlate with normal velocity: Length*Chord
// -->  Cd = SUM(kai/dt  * (u-uS)*sigma^3)/(1/2 *  U^2  * A )
// --> Incase of Sphere,where uS=0, U =1, D=1 : Cd = SUM(kai/dt  * (u-uS)*sigma^3)/(1/2 *  U^2  * pi*D^2/4 ) = 8*SUM(kai/dt  * (u)*sigma^3)/(pi )
// =======================================================

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


	// internal variables
	std::string name1, name2, name3, name4;

	printf("Saving data ...\n");

	if (!(Pars::opt_extra_data == 3))
	{
		// defining output variables
		std::ofstream ofs;
		// TODO: save Cx_Cy_Cz_Forces penalization - every timestep
		// force_pen(ngx,ngy,lambda, kai,ux_gpen,uy_gpen, uSx_gpen,uSy_gpen,sp_gpen, nvertn,vertex, savePars);
		// if (it == 0)
		// {
		// 	ofs.open("output/force_Cd_pen.dat" /*, std::ofstream::out | std::ofstream::app*/);
		// 	ofs << w20 << "time" << w20 << "force_x_pen" << w20 << "force_y_pen" << w20 << "Cx_pen" << w20 << "Cy_pen\n";
		// 	// ofs <<w20<< it*dt <<w20<< force_x_pen <<w20<< force_y_pen <<w20<< Cx_pen <<w20<< Cy_pen<<"\n";
		// 	ofs.close();
		// }
		// // else if (it>=1)
		// // {
		// ofs.open("output/force_Cd_pen.dat", std::ofstream::out | std::ofstream::app);
		// ofs << w20 << it * Pars::dt << w20 << force_x_pen << w20 << force_y_pen << w20 << Cx_pen << w20 << Cy_pen << "\n";
		// ofs.close();
		// }

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
	}

	// TODO: Saving common data and extra data
	if ((it == 0) || ((it % Pars::nt_sf) == 0) || (it == (Pars::sim_time/Pars::dt)-1)|| (cum_time >= Pars::comtime_sf) || (Pars::opt_extra_data == 3))
	{
		// defining output variables
		std::ofstream ofs;
		std::stringstream ss;

		cum_time = 0.0e0; // Reset for cumulative time saving file
		// Giving name and opening the files which want to save.
		// Giving "number of name" = "time step"
		if (Pars::opt_extra_data != 3)
		{
			ss << std::setw(7) << std::setfill('0') << it; //because c++ indexing is started from 0
		}
		else
		{
			ss << std::setw(7) << std::setfill('0') << it - 1; //because c++ indexing is started from 0
		}

		std::string nData = ss.str();
		// If not incase of "later saving grid data", (if not do the file can be rewired)
		if ((Pars::opt_extra_data == 1) || (Pars::opt_extra_data == 2)) // save only particles, or save particles and grid at same time
		{
			// ! Saving parameters
			if (it == 0)
			{
				ofs.open("output/parameters.dat");
				ofs << w16 << "Re" << w16 << "nt_sf" << w16 << "uin" << w16 << "vin" 
					<< w16 << "Courant" << w16 << "Uf" << w16 << "Df" << w16 << "dt" 
					<< w16 << "vis" << w16 << "dens" << w16 << "sigma" << w16 << "lx" 
					<< w16 << "ly" << w16 << "ubody" << w16 << "vbody" << w16 << "AoA\n";
				ofs << w16 << Pars::RE << w16 << Pars::nt_sf << w16 << Pars::u_inf 
					<< w16 << Pars::v_inf << w16 << Pars::Courant << w16 << Pars::U_inf 
					<< w16 << Pars::Df << w16 << Pars::dt << w16 << Pars::NU << w16 << Pars::RHO 
					<< w16 << Pars::sigma << w16 << Pars::Df << w16 << Pars::Df << w16 << Pars::ubody 
					<< w16 << Pars::vbody << w16 << Pars::AoA << "\n";
				ofs.close();
			}
			// Save common data (information of free particles):
			name1.append("output/x_vor_vel_particle_");
			name1.append(nData);
			name1.append(".csv");
			ofs.open(name1.c_str());
			ofs << "" << "xp" << "," << "yp" << "," << "gpz" << "," << "up" << "," << "vp" << "," << "sp" << /*"," << "b" <<*/ "," <<"w\n";
			for (int i = 0; i < np; i++)
			{
				//if(p.isActive[i] == true){
					ofs << "" << xp[i]
						<< "," << yp[i]
						<< "," << gpz[i]
						<< "," << up[i] - Pars::u_inf
						<< "," << vp[i] - Pars::v_inf
						<< "," << sp[i]
						// << "," << p.isboundary[i]
						<< "," << wp[i]
						<< "\n";
				//}
			}
			ofs.close();

			//body location
			name2.append("output/Body_location_");
			name2.append(nData);
			name2.append(".csv");
			ofs.open(name2.c_str());
			ofs << "" << "xp" << "," << "yp" << "\n";
			for (int i = 0; i < xb.size(); i++)
			{
				ofs << "" << xb[i]<< "," << yb[i]<< "\n";
			}
		}

		// if ( (Pars::opt_extra_data==2) || (Pars::opt_extra_data==3) )
		// {
		// // Calculate information of grid from particle info
		// grid_data(np, xp, yp, sp, gpz, up, vp, nvertn, vertex, savePars);
		// // SAVE grid position data-vorticity and velocity of grid data
		// name2.append("output/x_vor_vel_grid_");
		// name2.append(nData);
		// name2.append(".dat");
		// ofs.open(name2.c_str() /*, std::ofstream::out | std::ofstream::app*/);
		// ofs << w20 << "ngrx" << w20 << "ngry" << w20 << "grid_space" << w20 << "xgmax(1)" << w20 << "xgmax(2)" << w20 << "xgmin(1)" << w20 << "xgmin(2)"
		// 	<< "\n";
		// ofs << w20 << ngrx << w20 << ngry << w20 << grid_space << w20 << xgmax[0] << w20 << xgmax[1] << w20 << xgmin[0] << w20 << xgmin[1];

		// ofs << w20 << "xg" << w20 << "yg" << w20 << "ggz" << w20 << "ug" << w20 << "vg"
		// 	<< "\n";
		// for (int i = 0; i < ngrd; i++)
		// 	ofs << w20 << xg[i] << w20 << yg[i] << w20 << ggz[i] << w20 << ug[i] << w20 << vg[i] << "\n";
		// ofs.close();

		// // SAVE pressure of grid data
		// if (Pars::opt_extra_pres == 1)
		// {
		// 	name3.append("output/pres_grid_");
		// 	name3.append(nData);
		// 	name3.append(".dat");
		// 	ofs.open(name3.c_str());
		// 	ofs.close();
		// }
		// if (Pars::opt_extra_pres_wall == 1)
		// {
		// 	name4.append("output/pres_wall_");
		// 	name4.append(nData);
		// 	name4.append(".dat");
		// 	ofs.close();
		// }

		// }
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
