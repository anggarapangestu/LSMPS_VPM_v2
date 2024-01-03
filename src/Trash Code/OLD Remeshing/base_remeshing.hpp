#ifndef INCLUDED_REMESHING_BASE
#define INCLUDED_REMESHING_BASE
#include <iomanip> // std::setw, std::setfill
#include <sstream> // std::stringstream
#include <fstream> // std::c_str()

#ifndef INCLUDED_GLOBAL
#include "../../global.hpp"
#endif

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

// // #ifndef INCLUDED_NEIGHBOR_BASE
// // #include "../neighbor/base_grid.hpp"
// // #endif

using namespace std;

// // class base_grid;
class base_remeshing
{
public:
	// * creating instance
	// // base_grid d_base_grid;

	// TODO: transfer particle properties into grid
	// ! -- obsolete -- 
	// // void P2G(int nPi, const std::vector<double> &xPi, const std::vector<double> &yPi, const std::vector<double> &infoP, int nGxi, int nGyi, const std::vector<double> &xGi, const std::vector<double> &yGi, double hG, int kernel_tyPe, std::vector<std::vector<double>> &infoG);
	// TODO: kernel calculation
	void M4(double dex, double dey, double &Wx, double &Wy);
	void M04(double dex, double dey, double &Wx, double &Wy);
	void M6(double dex, double dey, double &Wx, double &Wy);
	void P14(double dex, double dey, double &Wx, double &Wy);
	void kernel(double dex, double dey, int kernel_type, double &Wx, double &Wy);
	
	// TODO: determine particle inside/outside body(es)
	void sign_particles(const Body &b, int npi, vector<double> &xpi, vector<double> &ypi, vector<bool> &sign_par, vector<bool> &inside);

	// TODO: searching neighbor between current particle and new particle distribution
	std::vector<std::vector<int>> inter_search(const Particle &p_scatter, const Particle &p_grid);
	std::vector<std::vector<int>> compulsive_search(const Particle &particle);
};

#endif
