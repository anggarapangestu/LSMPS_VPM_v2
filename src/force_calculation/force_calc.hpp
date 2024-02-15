#ifndef INCLUDED_FORCE_CALCULATION
#define INCLUDED_FORCE_CALCULATION

#include "../../Utils.hpp"

#ifndef INCLUDED_NEIGHBOR
#include "../neighbor/neighbor.hpp"
#endif

#ifndef INCLUDED_LSMPSa
#include "../LSMPS/LSMPSa.hpp"
#endif

#ifndef INCLUDED_LSMPSb
#include "../LSMPS/LSMPSb.hpp"
#endif

// #include "../save_data/save_data.hpp"
// #ifndef INCLUDED_BIOT_SAVART
// #include "../velocity_calculation/velocity_biot_savart.hpp"
// #endif
// #ifndef INCLUDED_PENALIZATION
// #include "../penalization/penalization.hpp"
// #endif

class force_calculation
{
private:
    // #pragma region instances
	// velocity_biot_savart d_base_poisson;
	// base_grid d_base_grid;
	// base_remeshing d_base_remeshing;
	// penalization d_penalization;
	// neighbor d_neighbor;
	// base_save_data d_base_save_data;
    // save_data save_step;
    // base_save_data d_base_save_data;
    // #pragma endregion

    // ===================================
	// ----------- Data Region -----------
	// ===================================
    // Internal data (sharing with the other function methods)
    int iter;		// The iteration at current calculation
    bool init;		// An flag said that current interation is the initial step

	// Direct pressure and velocity force calculation data
	LSMPSb interpolation;

	// Vorticity moment integral force calculation data
    std::vector<std::vector<double>> Ivortx;	// The vorticity moment in x direction [body part][previous 1,2,3 data]
	std::vector<std::vector<double>> Ivorty;	// The vorticity moment in y direction [body part][previous 1,2,3 data]
	std::vector<std::vector<double>> Ivortz;	// The vorticity moment in z direction [body part][previous 1,2,3 data]

	// Noca force calculation data
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

	// ===================================
	// ---------- Method Region ----------
	// ===================================

    // Internal method in force calculation
    double FD_diff1(const double h, const std::vector<double> &f, int order);

    // The force calculation
    void direct_force(const Particle& par, const Body& body, int bodyPart);                         // [Type 1] The direct force calculation
	void pen_force(const Particle& par, const std::vector<Body> &bodyList, std::string name);       // [Type 2] Penalization force calculation
	void int_vor_force(const Particle& par, const std::vector<Body> &bodyList, std::string name);   // [Type 3] Integral vorticity force calculation
    void imp_force(const Particle& par, const Body& body);      // The force liniear impulse force calculation

    // void pen_force(const Particle& par, const Body& body);   // Penalization force calculation

public:
	// Class Constructor
	force_calculation(): init(true){
		// Nothing to do here
	};

	// The force calculation manager
	void force_calc(const Particle& _particle, const std::vector<Body> &_bodyList, int _step, int _type, std::string name);

    // // The force calculation manager
    // void force_calc(const Particle& par, const Body& body, int step, int _type);
	
	// OLD IMPULSE METHOD
	void Force2(int iT, double A1, double A2, double A3, double A4, int nx, int ny, const Particle &p);
	void set_variables(int numberofparticles_x, int numberofparticles_y, int totalnumberofparticles);
	void insert_data(const Particle &p, int i);
	void sep(Particle p, std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d, std::vector<double> e, std::string s);
	void save_dat(Particle &p, int np, std::vector<int> &numbar, std::vector<int> &batas_bawah, std::vector<int> &batas_atas, 
				  std::vector<int> &batas_kiri, std::vector<int> &batas_kanan, std::string s, std::vector<double> &satu, std::vector<double> &dua,
				  std::vector<double> &tiga, std::vector<double> &empat, std::vector<double> &lima );
	
	// // OLD AND NOT used method
	// void save_dat(Particle &p, int np, std::vector<int> &numbar, std::vector<int> &batas_bawah, std::vector<int> &batas_atas, std::vector<int> &batas_kiri, std::vector<int> &batas_kanan, string s);
};


#endif