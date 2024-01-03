#ifndef INCLUDED_FORCE_CALCULATION
#define INCLUDED_FORCE_CALCULATION

#include "../../Utils.hpp"
#include "../save_data/save_data.hpp"
#include "../neighbor/neighbor.hpp"
#include "../LSMPS/LSMPSb.hpp"

#ifndef INCLUDED_LSMPS
#include "../LSMPS/LSMPSa.hpp"
#endif

#ifndef INCLUDED_BIOT_SAVART
#include "../velocity_calculation/velocity_biot_savart.hpp"
#endif

#ifndef INCLUDED_PENALIZATION
#include "../penalization/penalization.hpp"
#endif

class force_calculation
{
private:
    // #pragma region instances
	velocity_biot_savart d_base_poisson;
	// base_grid d_base_grid;
	// base_remeshing d_base_remeshing;
	penalization d_penalization;
	neighbor d_neighbor;
	// base_save_data d_base_save_data;
    // #pragma endregion

    // Impulse method [old]
    std::vector<double> _Isumx;
	std::vector<double> _Isumy;

	std::vector<int> a1a2;
	std::vector<int> a2a3;
	std::vector<int> a3a4;
	std::vector<int> a4a1;
	std::vector<int> particles_inside;

	std::vector<double> u_jmin1;                                
	std::vector<double> u_jmax1;	
	std::vector<double> v_imin1;                                
	std::vector<double> v_imax1;	

	std::vector<double> u_jmin2;                                
	std::vector<double> u_jmax2;	
	std::vector<double> v_imin2;                                
	std::vector<double> v_imax2;

	std::vector<double> _omegasumx;
	std::vector<double> _omegasumy;

	std::vector<int> point_lists;

	Particle _p;
    
    
    // Internal method 
    save_data save_step;
    // base_save_data d_base_save_data;     // Save force each time interval
    LSMPSb interpolation;
    int iter = 0;
    bool init = true;

    // Internal method in force calculation
    double FD_diff1(const double h, const std::vector<double> &f, int order);

    // The force calculation
    void direct_force(const Particle& par, const Body& body);   // The direct force calculation
    void pen_force(const Particle& par, const Body& body);      // Penalization force calculation
    void imp_force(const Particle& par, const Body& body);      // The impulse force calculation

    void force_linear_impulse(int it, int np, const std::vector<double> &xp, const std::vector<double> &yp,
							  const std::vector<double> &gpz, const std::vector<double> &sp);

public:
    // The force calculation manager
    void force_calc(const Particle& par, const Body& body, int step, int _type);     // Direct force calculation


    // void grid_data(int np, std::vector<double> &xp, std::vector<double> &yp, std::vector<double> &sp, std::vector<double> &gpz,
	// 			   std::vector<double> &up, std::vector<double> &vp, int &nvertn, std::vector<std::vector<double>> &vertex);
	
	// OLD IMPULSE METHOD
    void output(int it, Particle &p, Body &b, double &cum_time);

	void Force2(int iT, double A1, double A2, double A3, double A4, int nx, int ny, const Particle &p);
	void set_variables(int numberofparticles_x, int numberofparticles_y, int totalnumberofparticles);
	// void save_dat(Particle &p, int np, std::vector<int> &numbar, std::vector<int> &batas_bawah, std::vector<int> &batas_atas, std::vector<int> &batas_kiri, std::vector<int> &batas_kanan, string s);
	void insert_data(const Particle &p, int i);
	void sep(Particle p, std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d, std::vector<double> e, std::string s);
	void save_dat(Particle &p, int np, std::vector<int> &numbar, std::vector<int> &batas_bawah, std::vector<int> &batas_atas, 
				  std::vector<int> &batas_kiri, std::vector<int> &batas_kanan, std::string s, std::vector<double> &satu, std::vector<double> &dua,
				  std::vector<double> &tiga, std::vector<double> &empat, std::vector<double> &lima );

};


#endif