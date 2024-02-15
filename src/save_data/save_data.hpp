#ifndef INCLUDED_DATA_SAVING
#define INCLUDED_DATA_SAVING

#include "../../Utils.hpp"
#include "../grid_block/gridNode.hpp"
#include <iomanip> 			// Console log manager lib.

// Defined variable for short writting
#define w6 std::setw(6) 	// Set the text space of  6 for writting file
#define w8 std::setw(8) 	// Set the text space of  8 for writting file
#define w10 std::setw(10) 	// Set the text space of 10 for writting file
#define w12 std::setw(12) 	// Set the text space of 12 for writting file
#define w16 std::setw(16) 	// Set the text space of 16 for writting file
#define w20 std::setw(20) 	// Set the text space of 20 for writting file
#define w25 std::setw(25) 	// Set the text space of 25 for writting file
#define wR std::right		// Flush right the text for the next output
#define wL std::left		// Flush left the text for the next output

class save_data
{
	// Private data
	const std::string sumLogDir = Pars::opt_start_state == 0 ?
	"output/Parameter_initial.dat" : "output/Parameter_resume.dat";

	void save_par_state_2D(const Particle &_particle, std::string _name, int _type);
	void save_par_state_3D(const Particle &_particle, std::string _name, int _type);
	void save_par_state_init_vor(const Particle &_particle, std::string _name, int _type);

public:
	// Basic saving data state

	void save_par_interpolation(const Particle &_particle, std::string _name);
	void save_par_state(const Particle &_particle, std::string _name, int _type);
	void save_body_state(const Body &_body, std::string _name, int _type);
	void save_grid_node_state(const GridNode &_gridNode, std::string _name, int _type);

	// // Utilities
	// std::string save_data::saveName(int step);
	
	// Simulation summary log 
	std::string get_log_directory(){
		return this->sumLogDir;
	}
	void summary_log();
	void write_summary_log();
	
	//void particle_data_reading(int it, int &np, std::vector<double> &xp, std::vector<double> &yp, std::vector<double> &sp, std::vector<double> &gpz, std::vector<double> &up, std::vector<double> &vp);
	//void particle_data_reading(int &nvertn, std::vector<std::vector<double>> &vertex);
};

#endif
