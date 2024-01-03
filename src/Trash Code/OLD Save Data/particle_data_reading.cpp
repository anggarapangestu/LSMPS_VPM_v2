#include "data_saving.hpp"

void save_data::particle_data_reading(int it, int &np, vector<double> &xp, vector<double> &yp, vector<double> &sp, vector<double> &gpz, vector<double> &up, vector<double> &vp)
{
	// -- internal variables
	std::string _filename;

	printf("Reading particle's data files ... ");
	// Open the file which we want to convert to grid data/continue the stopped running
	// Giving "number of name" = "time step"
	std::stringstream ss;
	//if (it == 0 && Pars::opt_extra_data != 3)
	//{
	//	ss << std::setw(7) << std::setfill('0') << it; // because c++ indexing is started from 0
	//}
	//else
	//{
	//	ss << std::setw(7) << std::setfill('0') << it - 1; // because c++ indexing is started from 0
	//}
	ss << it;
	_filename.append("output/x_vor_vel_particle_");
	_filename.append(ss.str());
	_filename.append(".dat");

	// ===== COUNTING LENGTH OF FILEs = number of points/faces=================
	// --- counting np - number of particles
	// defining output variables
	std::ifstream ifs; // to read a file
	std::ofstream ofs; // to write a file
	printf("%s \n", _filename.c_str());

	ifs.open(_filename.c_str());
	np = 0;
	ifs.ignore(500, '\n'); //ignore first line
	while (!ifs.eof())
	{
		ifs >> xp[np] >> yp[np] >> gpz[np] >> up[np] >> vp[np] >> sp[np];
		np++;
	}
	ifs.close();

	np -= 1; // to ignore the first and last line
	// std::fill(sp.begin(), sp.begin() + np, Pars::sigma);

	// --------- Display number of particles read from file -----
	printf("Number of fluid particles from file, np = %d\n", np);

	//  ===========================================================================
	//  Open files, Read data, Rewire data to the files: force_Cd_pen.dat, force_Cd_lIm.dat, Isum.dat
	//  Cause: if the program was stopped at timestep 156, mean those files recorded 156 timestep,
	//  but for particle file, it depend on "save data per ( nt_sf ) time step", for example nt_sf=50, then we have
	//  particle file at 150. Then we must recalculate from step 150, means above files will re-record timestep
	//  151,152...156 for Cd, Isum ... --> we need reWrite before make new saving
	if (Pars::opt_cont == 1)
	{
		double *_time_pen_temp = new double[Pars::it_stop_full];
		double *_force_x_pen_temp = new double[Pars::it_stop_full];
		double *_force_y_pen_temp = new double[Pars::it_stop_full];
		double *_Cx_pen_temp = new double[Pars::it_stop_full];
		double *_Cy_pen_temp = new double[Pars::it_stop_full];
		double *_time_lIm_temp = new double[Pars::it_stop_full];
		double *_force_x_lIm_temp = new double[Pars::it_stop_full];
		double *_force_y_lIm_temp = new double[Pars::it_stop_full];
		double *_Cx_lIm_temp = new double[Pars::it_stop_full];
		double *_Cy_lIm_temp = new double[Pars::it_stop_full];
		double *_time_Isum_temp = new double[Pars::it_stop_full];
		// ---- Open files and save data temporarily
		// ----------------------------
		ifs.open("output/force_Cd_pen.dat");
		ifs.ignore(500, '\n'); // ignore first line
		for (int i = 0; i < Pars::it_stop_full; i++)
			ifs >> _time_pen_temp[i] >> _force_x_pen_temp[i] >> _force_y_pen_temp[i] >> _Cx_pen_temp[i] >> _Cy_pen_temp[i];
		ifs.close();
		// ----------------------------
		// storing first data only
		ifs.open("output/force_Cd_lIm.dat");
		ifs.ignore(500, '\n'); // ignore first line
		for (int i = 0; i < Pars::it_stop_full; i++)
			ifs >> _time_lIm_temp[i] >> _force_x_lIm_temp[i] >> _force_y_lIm_temp[i] >> _Cx_lIm_temp[i] >> _Cy_lIm_temp[i];
		ifs.close();

		// ---------- Open Isum file for calculating Cd linear impulse if the running was stopped
		// allocating array size
		std::vector<double> &Isumx = this->_Isumx;
		std::vector<double> &Isumy = this->_Isumy;
		Isumx.resize(Pars::it_stop_full);
		Isumy.resize(Pars::it_stop_full);
		std::fill(Isumx.begin(), Isumx.begin() + Pars::it_stop_full, 0.0e0);
		std::fill(Isumy.begin(), Isumy.begin() + Pars::it_stop_full, 0.0e0);

		ifs.open("output/Isum.dat");
		ifs.ignore(500, '\n'); // ignore first line
		for (int i = 0; i < Pars::it_stop_full; i++)
			ifs >> _time_Isum_temp[i] >> Isumx[i] >> Isumy[i];
		ifs.close();

		// ----- Rewrite the files
		// ------------------------
		ofs.open("output/force_Cd_pen.dat");
		ofs << w20 << "time" << w20 << "force_x_pen" << w20 << "force_y_pen" << w20 << "Cx_pen" << w20 << "Cy_pen\n";
		for (int i = 0; i < Pars::it_stop_full; i++)
			ofs << w20 << _time_pen_temp[i] << w20 << _force_x_pen_temp[i] << w20 << _force_y_pen_temp[i] << w20 << _Cx_pen_temp[i] << w20 << _Cy_pen_temp[i] << "\n";
		ofs.close();
		// ----------------
		ofs.open("output/force_Cd_lIm.dat");
		ofs << w20 << "time" << w20 << "force_x_lIm" << w20 << "force_y_lIm" << w20 << "Cx_lIm" << w20 << "Cy_lIm\n";
		for (int i = 0; i < Pars::it_stop_full; i++)
			ofs << w20 << _time_lIm_temp[i] << w20 << _force_x_lIm_temp[i] << w20 << _force_y_lIm_temp[i] << w20 << _Cx_lIm_temp[i] << w20 << _Cy_lIm_temp[i] << "\n";
		ofs.close();
		// ----------------
		ofs.open("output/Isum.dat");
		ofs << w20 << "time" << w20 << "Isumx" << w20 << "Isumy\n";
		for (int i = 0; i < Pars::it_stop_full; i++)
			ofs << w20 << _time_Isum_temp[i] << w20 << Isumx[i] << w20 << Isumy[i] << "\n";
		ofs.close();

		// -- delete internal variables
		delete[] _time_pen_temp;
		delete[] _time_lIm_temp;
		delete[] _time_Isum_temp;
		delete[] _force_x_pen_temp;
		delete[] _force_y_pen_temp;
		delete[] _force_x_lIm_temp;
		delete[] _force_y_lIm_temp;
		delete[] _Cx_pen_temp;
		delete[] _Cy_pen_temp;
		delete[] _Cx_lIm_temp;
		delete[] _Cy_lIm_temp;
	}

} // end of function

// ===============================================================================================
void save_data::particle_data_reading(int &nvertn, std::vector<std::vector<double>> &vertex)
{
	printf("Reading obstacle's data files ... ");
	//  !!!===== COUNTING LENGTH OF FILEs = number of points/faces=================
	double dummy1, dummy2;
	double *_vrtx;
	double *_vrty;
	// !! --- counting np - number of particles
	// defining output variables
	std::ifstream ifs; // to read a file
	std::ofstream ofs; // to write a file

	ifs.open("output/2Dbody.dat");
	nvertn = 0;
	while (!ifs.eof())
	{
		ifs >> dummy1 >> dummy2;
		nvertn++;
	}
	ifs.close();
	nvertn -= 1; // to ignore the last line

	_vrtx = new double[nvertn];
	_vrty = new double[nvertn];
	ifs.open("output/2Dbody.dat");
	int i = 0;
	while (!ifs.eof())
	{
		ifs >> _vrtx[i] >> _vrty[i];
		i++;
	}
	ifs.close();

	// -- store data to vertex
	vertex.resize(nvertn, std::vector<double>(2));
	for (int i = 0; i < nvertn; i++)
	{
		vertex[i][0] = _vrtx[i];
		vertex[i][1] = _vrty[i];
	}
	// -- clear internal data
	delete[] _vrtx;
	delete[] _vrty;
}
