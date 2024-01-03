#ifndef INCLUDED_GEOMETRY
#define INCLUDED_GEOMETRY

#include "../../Utils.hpp"		// Include data storage class

/**
 *  @brief A geometry class method, provides a body generator, body reader, and 
 *  geometry manager.
 *  
 *  @headerfile geometry.hpp
 */
class geometry
{
private:
	//2D Body Objects Generator

	void cylinder_generator(Body &_body, int part);		// [opt 1] Generate Circular Cylinder
	void square_generator(Body &_body, int part);		// [opt 2] Generate Square Cylinder
	void normal_plate_generator(Body &_body, int part);	// [opt 3] Generate Normal Plate
	void flat_plate_generator(Body &_body, int part);	// [opt 4] Generate Flat Plate
	void naca_generator(Body &_body, int part);			// [opt 5] Generate NACA Airfoil
	
	//3D Body Objects Generator

	void ball_sphere_generator(Body &_body, int part);		// [opt 1] Generate Sphere
	void cube_generator(Body &_body, int part);				// [opt 2] Generate Cube
	void normal_plate_3d_generator(Body &_body, int part);	// [opt 3] Generate 3D normal plate
	void flat_plate_3d_generator(Body &_body, int part);	// [opt 4] Generate 3D flat plate
	void torus_generator(Body &_body, int part);			// [opt 5] Generate Torus
	void heart_generator(Body &_body, int part);			// [opt 6] Generate Heart
	
	void read_3D_geometry(Body &_body, const std::string &bodyName);	// Body reader
	
	// Body velocity definition at each iteration

	void define_vel_variation(Body &_body);

public:
	// Geometry Method

	void generateBody(std::vector<Body> &_bodyList);
	void moving_body(int it, Body &_body);

	// Minimum distance from particle to body surface

	void distance_calc(Particle &_par, const Body &_body, const bool _abs);
	double distance_calc(const std::vector<double> &_coor, const Body &_body, const bool _abs);

	// Minimum distance from particle to body surface
	void eval_near_body(Particle &_particle, const std::vector<Body> &_bodyList);
	void eval_body_flag(Particle &_particle, const std::vector<Body> &_bodyList);
};

#endif
