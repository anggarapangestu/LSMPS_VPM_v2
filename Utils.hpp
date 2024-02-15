#ifndef INCLUDED_UTILS
#define INCLUDED_UTILS

#include "global.hpp"

/**
 *  @brief A property set of single body made of node or panel.
 *
 *  @headerfile Utils.hpp
 */
struct Body{
	// GENERAL DATA
	// ************
	// Body extremes coordinate
	double cen_pos[DIM];		// Body center coordinate
	double min_pos[DIM];		// Minimum extreme in each basis dimension
	double max_pos[DIM];		// Maximum extreme in each basis dimension

	// Body node translational velocity (*at each simulation iteration)
	std::vector<double> uT;		// Translation velocity in x direction
	std::vector<double> vT;		// Translation velocity in y direction
	std::vector<double> wT;		// Translation velocity in z direction
	
	// // body rotational velocity (*at each simulation iteration <?>)
	// std::vector<double> uR;
	// std::vector<double> vR;
	// std::vector<double> wR;
	
	// // body deformation velocity (*at each body note point <?>)
	// std::vector<double> uDEF;
	// std::vector<double> vDEF;
	// std::vector<double> wDEF;
	
	// // Body angular velocity (Still in 2D)
	// std::vector<double> avelocity;

	// NODE DATA
	// *********
	int n_node;					// Number of Nodes

	// Body node coordinate
	std::vector<double> x;		// Node coordinate in x basis
	std::vector<double> y;		// Node coordinate in y basis
	std::vector<double> z;		// Node coordinate in z basis

	// PANEL DATA
	// **********
	int n_panel;				// Number of Panel
	
	// Body panel middle coordinate (position)
	std::vector<double> x_m;	// Panel Coordinate in x basis
	std::vector<double> y_m;	// Coordinate in x basis
	std::vector<double> z_m;	// Coordinate in x basis
	std::vector<double> size;	// The panel size (2D: length, 3D: area)
	
	// IMPORTANT: normal vector directed outward the body (pointing outside)
	// Panel normal vector
	std::vector<double> x_n;
	std::vector<double> y_n;
	std::vector<double> z_n;

	// // Panel middle coordinate
	// std::vector<std::vector<int>> node_list;	// used only in 3D panel creation

	// Constructor
	Body(){};

	// Deconstructor
	~Body(){};
};

/**
 *  @brief A data set of particle inside the domain.
 *
 *  @headerfile Utils.hpp
 */
struct Particle{
	// Physical Properties
	// *******************
	// Particle position coordinates
	std::vector<double> x;			// Coordinate in x basis
	std::vector<double> y;			// Coordinate in y basis
	std::vector<double> z;			// Coordinate in z basis
	
 	// Particle physical properties
	std::vector<double> u;			// Velocity in x direction
	std::vector<double> v;			// Velocity in y direction
	std::vector<double> w;			// Velocity in x direction
	std::vector<double> gx;			// Vortex strength x (gamma_x) [*]
	std::vector<double> gy;			// Vortex strength y (gamma_y) [*]
	std::vector<double> gz;			// Vortex strength z (gamma_z) [*]
	std::vector<double> vortx;      // Vorticity x (omega_x)
	std::vector<double> vorty;      // Vorticity y (omega_y)
	std::vector<double> vortz;      // Vorticity z (omega_z)
	std::vector<double> vorticity;	// Vorticity scalar field (or Vorticity absolute value)
	std::vector<double> P;			// Pressure scalar field
	// Note on [*] : Not gonna used, calculate once and use once in simulation, seem redundant on memory

	// Computational Properties
	// ************************
	// Basic computational Properties
	int num;						// Number of particle in the class
	bool isAdapt;					// A flag to said particle is adapted
	std::vector<double> s;			// Particle size (diameter)
	std::vector<int> level;			// The level of particle resolution (related to size)
	std::vector<std::vector<int>> neighbor; // Index list of neighbor particle [A verlet list]

	// [GROUP 1] Data using *Grid Node* container
	std::vector<int> nodeID;		// ID label of node container

	// [GROUP 2] Data using *Adaptive Tree Cell* container
	std::vector<int> basis_label;	// ID label of basis cell
	std::vector<int> cell_label;	// ID label of cell

	// Penalization parameters
	std::vector<double> chi;		// Particle penalization mask value
	std::vector<double> R;			// Shortest distance from the solid surface; sign: [-] inward, [+] outward

	// Flag Marker and near body variable
	std::vector<int> bodyPart;			// Nearest body part ID [-1:= for not near to any body part]
	std::vector<bool> isActive;			// Indicates active particle => Particle that have a vorticity value [@param: transient particle]
	std::vector<bool> isNearSurface;	// Flag for near body surface
	std::vector<bool> insideBody;		// Flag for inside body

	/* Near Body Criteria Illustration
       	 __________________________________________________
		|        .             x _____x_______   x     .   |       Note:
		| .          .  *    x  |#2 x   x     \       .    |		> In the given illustration 3 objects is existing
		|     *_________        |      x       \ x   .     |		   inside the domain
		|   * /#1     * \ *   x |____________x__\  x     . |		> Particle is grouped in 4 criteria:
		|    /  *  *     \        x           x            |		   - [.] : free particle (not near with any body)
		|   /  *     *   / *   x  o  ____o__________o__ o  |           - [*] : particle near and inside body #1
		|   \     *   * /       o   |#3 o      o o     |   |           - [x] : particle near and inside body #2
		|  * \_________/ *        o |               o  |   |           - [o] : particle near and inside body #3
		|   *                    o  | o    o o    o    |   |        > The "bodyPart" store the nearest body part of 
		|     *         .      .    |   o      o     o | o |           the current particle, for free particle group
		|    .    .                o|__________________|o  |           value is set to -1
		|  .  .           .    .     o   o    o     o      |           
		|__________________________________________________|           
	*/

	// // Boundary terms <?> For boundary treatment <?> Not really need for LSMPS
	std::vector<bool> isBoundary;
	std::vector<double> boundaryVal;
	// std::vector<bool> inside; 

	// Constructor
	Particle() : num(0), isAdapt(false){};

	// Deconstructor
	~Particle(){};
};

/**
 *  @brief A simulation utilities function list.
 *
 *  @headerfile Utils.hpp
 */
class simUtil{
private:
	int iterDigitLen;
public:
	void startCounter(int step);
	void printHeader(int step);
	void predictCompTime(int step, double _currTime);
	std::string saveName(const int step);
	void stabilityEval(const Particle &_particle, std::vector<double> &_maxValue);
	void saveResidual(Particle &_particle, int step);
};

#endif