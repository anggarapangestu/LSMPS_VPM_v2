#include "geometry.hpp"

/*/ ============================================== ///
	Note in 2D geometry:
	* Body node is arranged in positive convention 
	  direction (CCW)
/// ============================================== /*/

/**
 *  @brief A 2D cylinder body geometry generator [OPTION 1].
 * 
 *  @param	_body	The body container to be generated.
 *  @param	part	The body part ID.
 *  
 *  @headerfile geometry.hpp
 */
void geometry::cylinder_generator(Body &b, int part)
{
	// update body center position
	b.cen_pos[0] = Pars::x_body_cen[part];
	b.cen_pos[1] = Pars::y_body_cen[part];

	// internal variables
	int NUM = Pars::n_a;
	double *_alfa = new double [NUM];
	double d_alfa = (360.0e0/NUM)*M_PI/180.0e0;
	double panel_size = std::sin(d_alfa/2.0) * Pars::Lref[part];
	double min, max;

	// surface point list of position angle (theta)
	_alfa[0] = 0.0e0;
	for (int i = 0; i < NUM; i++)
	{
		_alfa[i+1] = _alfa[i] + d_alfa;
	}

	// resize x and y coordinate of obstacle
	b.x.clear(); b.x.resize(NUM);
	b.y.clear(); b.y.resize(NUM);
	b.x_n.clear(); b.x_n.resize(NUM);
	b.y_n.clear(); b.y_n.resize(NUM);
	b.x_m.clear(); b.x_m.resize(NUM);
	b.y_m.clear(); b.y_m.resize(NUM);
	b.size.clear(); b.size.resize(NUM);
	b.n_node = NUM;
	b.n_panel = NUM;
	
	// Value of extremes local
	min = -Pars::Lref[part]/2;
	max = Pars::Lref[part]/2;

	// surface body point coordinate, panel midpoint coordinate, panel size, and panel normal vector
	for (int i = 0; i < NUM; i++)
	{
		// Body node coordinate
		b.x[i] = b.cen_pos[0] + Pars::Lref[part]/2 * std::cos(_alfa[i]);
		b.y[i] = b.cen_pos[1] + Pars::Lref[part]/2 * std::sin(_alfa[i]);

		// Body node panel normal
		b.x_n[i] = std::cos(_alfa[i] + d_alfa/2);
		b.y_n[i] = std::sin(_alfa[i] + d_alfa/2);
		
		// Body panel midpoint coordinate
		if (i > 0){
			b.x_m[i-1] = (b.x[i] + b.x[i-1])/2;
			b.y_m[i-1] = (b.y[i] + b.y[i-1])/2;
			
			// Manual calculation of panel size
			// b.size[i] = std::sqrt((b.x[i] - b.x[i-1])*(b.x[i] - b.x[i-1]) + (b.y[i] - b.y[i-1])*(b.y[i] - b.y[i-1]));
		}
		
		// Panel size
		b.size[i] = panel_size;
	}
	
	// Update the rest data
	b.x_m[NUM-1] = (b.x[NUM-1] + b.x[0])/2;
	b.y_m[NUM-1] = (b.y[NUM-1] + b.y[0])/2;
	
	// Assign the extremes point
	// X position
	b.min_pos[0] = b.cen_pos[0] + min;
	b.max_pos[0] = b.cen_pos[0] + max;

	// Y position
	b.min_pos[1] = b.cen_pos[1] + min;
	b.max_pos[1] = b.cen_pos[1] + max;

	// deallocating memory
	delete[] _alfa;
}

/**
 *  @brief A 2D square body geometry generator [OPTION 2].
 * 
 *  @param	_body	The body container to be generated.
 *  @param	part	The body part ID.
 *  
 *  @headerfile geometry.hpp
 */
void geometry::square_generator(Body &b, int part)
{
	// update body center position
	b.cen_pos[0] = Pars::x_body_cen[part];
	b.cen_pos[1] = Pars::y_body_cen[part];

	// internal variables declaration
	double lenpanel, pos_max, pos_min;
	int numpanel, N, id;

	// assigning internal variables
	numpanel = Pars::n_a / 4;
	
	// calculating internal variables
	N = 4 * numpanel;
	lenpanel = Pars::Lref[part]/numpanel;
	pos_max = Pars::Lref[part]/2.0e0;	// maximum vertex coordinate
	pos_min = -Pars::Lref[part]/2.0e0;	// minimum vertex coordinate

	// resize x and y coordinate of body object and body initialization
	b.x.clear(); b.x.resize(N);
	b.y.clear(); b.y.resize(N);
	b.x_m.clear(); b.x_m.resize(N);
	b.y_m.clear(); b.y_m.resize(N);
	b.x_n.clear(); b.x_n.resize(N);
	b.y_n.clear(); b.y_n.resize(N);
	b.size.clear(); b.size.resize(N);
	b.n_node = N;
	b.n_panel = N;

	// surface body point coordinate, panel midpoint coordinate, and panel normal vector
	id = 0;		// The index of body node
	for (int i = 0; i < numpanel; i++) // SEGMENT 1: Vertex node 1-2 (bottom segment)
	{
		b.x[id] = b.cen_pos[0] + pos_min + i*lenpanel;
		b.y[id] = b.cen_pos[1] + pos_min;
		b.x_m[id] = b.x[id] + (0.5e0 * lenpanel);
		b.y_m[id] = b.y[id];
		b.x_n[id] = 0;
		b.y_n[id] = -1;
		b.size[id] = lenpanel;
		id++;
	}
	for (int i = 0; i < numpanel; i++) // SEGMENT 2: Vertex node 2-3 (right segment)
	{
		b.x[id] = b.cen_pos[0] + pos_max;
		b.y[id] = b.cen_pos[1] + pos_min + i*lenpanel;
		b.x_m[id] = b.x[id];
		b.y_m[id] = b.y[id] + (0.5e0 * lenpanel);
		b.x_n[id] = 1;
		b.y_n[id] = 0;
		b.size[id] = lenpanel;
		id++;
	}
	for (int i = 0; i < numpanel; i++) // SEGMENT 3: Vertex node 3-4 (top segment)
	{
		b.x[id] = b.cen_pos[0] + pos_max - i*lenpanel;
		b.y[id] = b.cen_pos[1] + pos_max;
		b.x_m[id] = b.x[id] - (0.5e0 * lenpanel);
		b.y_m[id] = b.y[id];
		b.x_n[id] = 0;
		b.y_n[id] = 1;
		b.size[id] = lenpanel;
		id++;
	}
	for (int i = 0; i < numpanel; i++) // SEGMENT 4: Vertex node 4-1 (left segment)
	{
		b.x[id] = b.cen_pos[0] + pos_min;
		b.y[id] = b.cen_pos[1] + pos_max - i*lenpanel;
		b.x_m[id] = b.x[id];
		b.y_m[id] = b.y[id] - (0.5e0 * lenpanel);
		b.x_n[id] = -1;
		b.y_n[id] = 0;
		b.size[id] = lenpanel;
		id++;
	}

	// Assign the extremes point
	basis_loop(d){
		b.min_pos[d] = b.cen_pos[d] + pos_min;
		b.max_pos[d] = b.cen_pos[d] + pos_max;
	}
}

/**
 *  @brief A 2D plate normal to freestream body geometry generator [OPTION 3].
 * 
 *  @param	_body	The body container to be generated.
 *  @param	part	The body part ID.
 *  
 *  @headerfile geometry.hpp
 */
void geometry::normal_plate_generator(Body &b, int part)
{
	// update body center position
	b.cen_pos[0] = Pars::x_body_cen[part];
	b.cen_pos[1] = Pars::y_body_cen[part];

	// internal variables declaration
	double lenpanelx, lenpanely, L, W;
	double *x_max = new double[2];
    double *x_min = new double[2];
	int numpanelx, numpanely, N, id;

	// assigning internal variables
	L = Pars::Lref[part];			// Plate length (y length)
	W = L*Pars::H_star[part];		// Plate width ~ thickness (x length)
	
	// rough calculation of panel number
	numpanely = std::floor((L/(2*(L+W))) * Pars::n_a);
	numpanelx = std::ceil(numpanely*Pars::H_star[part]);
	
	// calculating internal variables
	N = 2*(numpanelx + numpanely);
	lenpanely = L/numpanely;
	lenpanelx = W/numpanelx;
	// maximum vertex coordinate
		x_max[0] = W;			// x_max
		x_max[1] = 0.5 * L;		// y_max
	// minimum vertex coordinate
		x_min[0] = 0.0;			// x_min
		x_min[1] = -0.5 * L;	// y_min

	// resize x and y coordinate of body object and body initialization
	b.x.clear(); b.x.resize(N);
	b.y.clear(); b.y.resize(N);
	b.x_m.clear(); b.x_m.resize(N);
	b.y_m.clear(); b.y_m.resize(N);
	b.x_n.clear(); b.x_n.resize(N);
	b.y_n.clear(); b.y_n.resize(N);
	b.size.clear(); b.size.resize(N);
	b.n_node = N;
	b.n_panel = N;

	// surface body point coordinate, panel midpoint coordinate, and panel normal vector
	id = 0;		// The index of body node
	for (int i = 0; i < numpanelx; i++) // SEGMENT 1: Vertex node 1-2 (bottom segment)
	{
		b.x[id] = b.cen_pos[0] + x_min[0] + i*lenpanelx;
		b.y[id] = b.cen_pos[1] + x_min[1];
		b.x_m[id] = b.x[id] + lenpanelx/2;
		b.y_m[id] = b.y[id];
		b.x_n[id] = 0;
		b.y_n[id] = -1;
		b.size[id] = lenpanelx;
		id++;
	}
	for (int i = 0; i < numpanely; i++) // SEGMENT 2: Vertex node 2-3 (right segment)
	{
		b.x[id] = b.cen_pos[0] + x_max[0];
		b.y[id] = b.cen_pos[1] + x_min[1] + i*lenpanely;
		b.x_m[id] = b.x[id];
		b.y_m[id] = b.y[id] + lenpanely/2;
		b.x_n[id] = 1;
		b.y_n[id] = 0;
		b.size[id] = lenpanely;
		id++;
	}
	for (int i = 0; i < numpanelx; i++) // SEGMENT 3: Vertex node 3-4 (top segment)
	{
		b.x[id] = b.cen_pos[0] + x_max[0] - i*lenpanelx;
		b.y[id] = b.cen_pos[1] + x_max[1];
		b.x_m[id] = b.x[id] - lenpanelx/2;
		b.y_m[id] = b.y[id];
		b.x_n[id] = 0;
		b.y_n[id] = 1;
		b.size[id] = lenpanelx;
		id++;
	}
	for (int i = 0; i < numpanely; i++) // SEGMENT 4: Vertex node 4-1 (left segment)
	{
		b.x[id] = b.cen_pos[0] + x_min[0];
		b.y[id] = b.cen_pos[1] + x_max[1] - i*lenpanely;
		b.x_m[id] = b.x[id];
		b.y_m[id] = b.y[id] - lenpanely/2;
		b.x_n[id] = -1;
		b.y_n[id] = 0;
		b.size[id] = lenpanely;
		id++;
	}
	
	// Assign the extremes point
	basis_loop(d){
		b.min_pos[d] = b.cen_pos[d] + x_min[d];
		b.max_pos[d] = b.cen_pos[d] + x_max[d];
	}
	
	// deallocating memory
	delete[] x_max;
	delete[] x_min;
}

/**
 *  @brief A 2D flat plate body geometry generator [OPTION 4].
 * 
 *  @param	_body	The body container to be generated.
 *  @param	part	The body part ID.
 *  
 *  @headerfile geometry.hpp
 */
void geometry::flat_plate_generator(Body &b, int part)
{
	// update body center position
	b.cen_pos[0] = Pars::x_body_cen[part];
	b.cen_pos[1] = Pars::y_body_cen[part];

	// internal variables declaration
	double lenpanelx, lenpanely, _L, _H;
	double *x_max = new double[2];
    double *x_min = new double[2];
	int numpanelx, numpanely, N, id;

	// assigning internal variables
	_L = Pars::Lref[part];			// Plate length (x length)
	_H = _L*Pars::H_star[part];		// Plate height ~ thickness (y length)
	
	// rough calculation of panel number
	numpanelx = std::floor((_L/(2*(_L + _H))) * Pars::n_a);
	numpanely = std::ceil(numpanelx*Pars::H_star[part]);
	
	// calculating internal variables
	N = 2*(numpanelx + numpanely);
	lenpanelx = _L/numpanelx;
	lenpanely = _H/numpanely;
	// maximum vertex coordinate
		x_max[0] = _L;			// x_max
		x_max[1] = 0.5 * _H;	// y_max
	// minimum vertex coordinate
		x_min[0] = 0.0;			// x_min
		x_min[1] = -0.5 * _H;	// y_min

	// resize x and y coordinate of body object and body initialization
	b.x.clear(); b.x.resize(N);
	b.y.clear(); b.y.resize(N);
	b.x_m.clear(); b.x_m.resize(N);
	b.y_m.clear(); b.y_m.resize(N);
	b.x_n.clear(); b.x_n.resize(N);
	b.y_n.clear(); b.y_n.resize(N);
	b.size.clear(); b.size.resize(N);
	b.n_node = N;
	b.n_panel = N;

	// surface body point coordinate, panel midpoint coordinate, and panel normal vector
	id = 0;		// The index of body node
	for (int i = 0; i < numpanelx; i++) // SEGMENT 1: Vertex node 1-2 (bottom segment)
	{
		b.x[id] = b.cen_pos[0] + x_min[0] + i*lenpanelx;
		b.y[id] = b.cen_pos[1] + x_min[1];
		b.x_m[id] = b.x[id] + lenpanelx/2;
		b.y_m[id] = b.y[id];
		b.x_n[id] = 0;
		b.y_n[id] = -1;
		b.size[id] = lenpanelx;
		id++;
	}
	for (int i = 0; i < numpanely; i++) // SEGMENT 2: Vertex node 2-3 (right segment)
	{
		b.x[id] = b.cen_pos[0] + x_max[0];
		b.y[id] = b.cen_pos[1] + x_min[1] + i*lenpanely;
		b.x_m[id] = b.x[id];
		b.y_m[id] = b.y[id] + lenpanely/2;
		b.x_n[id] = 1;
		b.y_n[id] = 0;
		b.size[id] = lenpanely;
		id++;
	}
	for (int i = 0; i < numpanelx; i++) // SEGMENT 3: Vertex node 3-4 (top segment)
	{
		b.x[id] = b.cen_pos[0] + x_max[0] - i*lenpanelx;
		b.y[id] = b.cen_pos[1] + x_max[1];
		b.x_m[id] = b.x[id] - lenpanelx/2;
		b.y_m[id] = b.y[id];
		b.x_n[id] = 0;
		b.y_n[id] = 1;
		b.size[id] = lenpanelx;
		id++;
	}
	for (int i = 0; i < numpanely; i++) // SEGMENT 4: Vertex node 4-1 (left segment)
	{
		b.x[id] = b.cen_pos[0] + x_min[0];
		b.y[id] = b.cen_pos[1] + x_max[1] - i*lenpanely;
		b.x_m[id] = b.x[id];
		b.y_m[id] = b.y[id] - lenpanely/2;
		b.x_n[id] = -1;
		b.y_n[id] = 0;
		b.size[id] = lenpanely;
		id++;
	}
	
	// Assign the extremes point
	basis_loop(d){
		b.min_pos[d] = b.cen_pos[d] + x_min[d];
		b.max_pos[d] = b.cen_pos[d] + x_max[d];
	}
	
	// deallocating memory
	delete[] x_max;
	delete[] x_min;
}

/**
 *  @brief A 2D symmetric NACA 4-series airfoil body geometry generator [OPTION 5].
 * 
 *  @param	_body	The body container to be generated.
 *  @param	part	The body part ID.
 *  
 *  @headerfile geometry.hpp
 */
void geometry::naca_generator(Body &b, int part)
{
	// update body center position (actualy the leading edge of airfoil)
	b.cen_pos[0] = Pars::x_body_cen[part];
	b.cen_pos[1] = Pars::y_body_cen[part];

	// internal variables
	int NUM = Pars::n_a + (Pars::n_a % 2 == 1 ? 0 : 1);		// Number must be odd

	double *x    = new double[NUM];		// Airfoil x position (in chord length scale)
	double *y    = new double[NUM];		// Airfoil y position (in chord length scale)
	double *yt   = new double[NUM];		// Thickness paramater t (*Array is stored in case want to write to file)
	double *yc   = new double[NUM];		// Thickness parameter c (*Array is stored in case want to write to file)

	// Airfoil node calculation
	// ************************
	// Set the local coordinate data
	double theta;
	for (int i = 0; i < NUM; i++){
		theta = i * (2.0e0 * M_PI / (NUM-1));
		x[i] = 0.5e0 * (std::cos(theta) + 1);	// Always positive
	}

	// Calculate the airfoil position
	// *Start at one node after leading edge, stop at one node before trailing edge
	double l_c, sin_theta, cos_theta;
	for (int i = 1 + NUM/2; i < (NUM-1); i++){
		// Calculate the thickness parameter c
		if (x[i] < Pars::p_a){
			// Location before maximum thickness
			if (Pars::p_a == 0.0e0){
				yc[i] = 0.0e0;
			}
			else{
				yc[i] = Pars::m_a / std::pow(Pars::p_a,2.0) * (2.0*Pars::p_a*x[i] - std::pow(x[i],2.0));
			}
		}
		else{
			// Location after maximum thickness
			yc[i] = Pars::m_a / std::pow((1.0-Pars::p_a),2.0) * ((1.0-2.0*Pars::p_a) + 2.0*Pars::p_a*x[i] - std::pow(x[i],2.0));
		}

	    // Calculate the thickness parameter t
		// yt[i] = Pars::t_a*5.0*(0.2969*sqrt(x[i]) - 0.126*x[i] - 0.35160*std::pow(x[i],2) + 0.2843*std::pow(x[i],3) - 0.1015*std::pow(x[i],4));
		yt[i] = Pars::t_a*5.0*(0.2969*sqrt(x[i]) - 0.126*x[i] - 0.35160*std::pow(x[i],2) + 0.2843*std::pow(x[i],3) - 0.1036*std::pow(x[i],4));
		
		// Calculate sine relation
		l_c = sqrt(std::pow((x[i]-x[i-1]),2) + std::pow((yc[i]-yc[i-1]),2));
	    sin_theta = (yc[i] - yc[i-1])/l_c;
	    cos_theta = (x[i] - x[i-1])/l_c;

	    // Calculate the final aerodynamic data
		x[i] = x[i] - yt[i]*sin_theta;
	    y[i] = yc[i] + yt[i]*cos_theta;
	    x[(NUM-1)-i] = x[i] + yt[i]*sin_theta;
	    y[(NUM-1)-i] = yc[i] - yt[i]*cos_theta;
	}

	// Set the trailing edge position
	x[0]     = 1.0;
	y[0]     = 0.0;
	x[NUM-1] = 1.0;
	y[NUM-1] = 0.0;
	
	// Set the leading edge position
	x[NUM/2] = 0.0;
	y[NUM/2] = 0.0;	

	// Translate coordinate toward aerodynamic center
	if (Pars::vib != 99){
		for (int i = 0; i < NUM; i++){
			x[i] -= Pars::acx;
			y[i] -= Pars::acy;
		}
	}

	// Airfoil rotation in Angle of Attack angle and body geometry assignment
	// **********************************************************************
	double R[2][2];		// Rotation matrix
	R[0][0] =  std::cos(Pars::AoA*M_PI/180.0e0);
	R[0][1] =  std::sin(Pars::AoA*M_PI/180.0e0);
	R[1][0] = -std::sin(Pars::AoA*M_PI/180.0e0);
	R[1][1] =  std::cos(Pars::AoA*M_PI/180.0e0);

	// resize x and y coordinate of obstacle
	b.x.clear(); b.x.resize(NUM-1);
	b.y.clear(); b.y.resize(NUM-1);
	b.x_n.clear(); b.x_n.resize(NUM-1);
	b.y_n.clear(); b.y_n.resize(NUM-1);
	b.x_m.clear(); b.x_m.resize(NUM-1);
	b.y_m.clear(); b.y_m.resize(NUM-1);
	b.size.clear(); b.size.resize(NUM-1);
	b.n_node = NUM-1;
	b.n_panel = NUM-1;
	
	// Initialize the extremes point
	basis_loop(d){
		b.min_pos[d] = b.cen_pos[d];
		b.max_pos[d] = b.cen_pos[d];
	}
	
	// Perform rotation transformation then assign the rotated data
	double xTrans, yTrans;
	for (int i = 0; i < NUM; i++){
		// Rotation and scaling transformation
		xTrans = (R[0][0]*x[i] + R[0][1]*y[i]) * Pars::Lref[part];
		yTrans = (R[1][0]*x[i] + R[1][1]*y[i]) * Pars::Lref[part];
		
		if (i < NUM-1){
			// Body node location
			b.x[i] = b.cen_pos[0] + xTrans;
			b.y[i] = b.cen_pos[1] + yTrans;

			// Evaluate extremes
			b.min_pos[0] = std::min<double>(b.min_pos[0], b.x[i]);
			b.min_pos[1] = std::min<double>(b.min_pos[1], b.y[i]);
			b.max_pos[0] = std::max<double>(b.max_pos[0], b.x[i]);
			b.max_pos[1] = std::max<double>(b.max_pos[1], b.y[i]);
		}
		
		if (i > 0){
			// Body panel midpoint coordinate
			b.x_m[i-1] = (b.cen_pos[0] + xTrans + b.x[i-1])/2;
			b.y_m[i-1] = (b.cen_pos[1] + yTrans + b.y[i-1])/2;
			
			// Panel size
			b.size[i-1] = std::sqrt(std::pow((b.x[i] - b.x[i-1]),2.0) + std::pow((b.y[i] - b.y[i-1]),2.0));

			// Body node panel normal
			b.x_n[i-1] = b.y_m[i-1] / b.size[i-1];
			b.y_n[i-1] = -b.x_m[i-1] / b.size[i-1];
		}
	}

	// // Not including the last point (it is exactly the initial point)
	// // Update the rest data
	// // Body panel midpoint coordinate
	// b.x_m[NUM-1] = (b.x[NUM-1] + b.x[0])/2;
	// b.y_m[NUM-1] = (b.y[NUM-1] + b.y[0])/2;
	
	// // Panel size
	// b.size[NUM-1] = std::sqrt(std::pow((b.x[0] - b.x[NUM-1]),2.0) + std::pow((b.y[0] - b.y[NUM-1]),2.0));

	// // Body node panel normal
	// b.x_n[NUM-1] = b.y_m[NUM-1] / b.size[NUM-1];
	// b.y_n[NUM-1] = -b.x_m[NUM-1] / b.size[NUM-1];

	// deallocate memory
	delete[] x, y, yt, yc;
}

/*
Note in geometry creation
- There is nothing here!
*/