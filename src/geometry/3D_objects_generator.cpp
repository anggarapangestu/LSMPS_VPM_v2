#include "geometry.hpp"
#include "../save_data/save_data.hpp"

// In current 3D node is the same to panel

/**
 *  @brief A 3D shpere body geometry generator [OPTION 1].
 * 
 *  @param	_body	The body container to be generated.
 *  @param	part	The body part ID.
 *  
 *  @headerfile geometry.hpp
 */
void geometry::ball_sphere_generator(Body &b, int part){
    // update body center position
	b.cen_pos[0] = Pars::x_body_cen[part];
	b.cen_pos[1] = Pars::y_body_cen[part];
    b.cen_pos[2] = Pars::z_body_cen[part];

    // internal variables
    int NUM = 0;            // Number of data point
    double xPos;            // Temporary x coordinate
    double yPos;            // Temporary y coordinate
    double zPos;            // Temporary z coordinate
    double dT = M_PI/60;    // Elevation angle increment
    double dP;              // Azimuth angle increment
    int numAz;              // Number of azimuth calculation
    double R = Pars::Lref[part] / 2.0;  // The radius
    std::vector<double> _x;     // List of coordinate (x basis container)
    std::vector<double> _y;     // List of coordinate (y basis container)
    std::vector<double> _z;     // List of coordinate (z basis container)

    // Generate the body data
    b.size.clear();
    double minSize = R*R;
                                // Avoid truncation error  
                                         /* ^ */
    for(double theta = dT; theta < (M_PI - dT/2); theta += dT){
        // Iterate through elevation [0, pi]
        // Calculate the particle distribution at current level
        dP = dT / std::sin(theta);
        numAz = 2*M_PI / dP;
        dP = 2*M_PI / numAz;
        for(double phi = 0.0; phi < 2*M_PI; phi += dP){
            // Iterate through azimuth [0, 2pi]
            xPos = R * std::cos(phi) * std::sin(theta);
            yPos = R * std::sin(phi) * std::sin(theta);
            zPos = R * std::cos(theta);

            // Push the data coordinate
            _x.push_back(xPos);
            _y.push_back(yPos);
            _z.push_back(zPos);

            // Algorithm to calculate panel size (rough approximation)
            double s1 = 2 * R * std::sin(dT/2);
            double s2 = 2 * R * std::sin(theta) * std::sin(dP/2);
            double _size = s1 * s2;
            
            b.size.push_back(_size);
            minSize = std::min<double> (minSize, _size);

            // Push the data coordinate
            NUM++;
        } 
    }
    
    // Push the top and bottom point
    _x.push_back(0.0);  _x.push_back(0.0);
    _y.push_back(0.0);  _y.push_back(0.0);
    _z.push_back(R);    _z.push_back(-R);
    b.size.push_back(minSize);
    b.size.push_back(minSize);
    NUM += 2;

    // reserve the body coodinate data
	b.x.clear(); b.x.resize(NUM);
	b.y.clear(); b.y.resize(NUM);
    b.z.clear(); b.z.resize(NUM);
	b.x_n.clear(); b.x_n.resize(NUM);
	b.y_n.clear(); b.y_n.resize(NUM);
    b.z_n.clear(); b.z_n.resize(NUM);
	b.x_m.clear(); b.x_m.resize(NUM);
	b.y_m.clear(); b.y_m.resize(NUM);
    b.z_m.clear(); b.z_m.resize(NUM);
    
    // Assign the data into body container
    for (int i = 0; i < NUM; i++){
        // Update the normal direction
        b.x_n[i] = _x[i]/R;
        b.y_n[i] = _y[i]/R;
        b.z_n[i] = _z[i]/R;

        // Update the body coordinate
        b.x[i] = b.cen_pos[0] + _x[i];
        b.y[i] = b.cen_pos[1] + _y[i];
        b.z[i] = b.cen_pos[2] + _z[i];

        b.x_m[i] = b.x[i];
        b.y_m[i] = b.y[i];
        b.z_m[i] = b.z[i];
    }
    
    // update other body data
    b.n_node = NUM;
	b.n_panel = NUM;

    basis_loop(d){
        b.min_pos[d] = -R;
        b.max_pos[d] = R;
    }
}


/**
 *  @brief A 3D cube body geometry generator [OPTION 2].
 * 
 *  @param	_body	The body container to be generated.
 *  @param	part	The body part ID.
 *  
 *  @headerfile geometry.hpp
 */
void geometry::cube_generator(Body &b, int part)
{
	// update body center position
	b.cen_pos[0] = Pars::x_body_cen[part];
	b.cen_pos[1] = Pars::y_body_cen[part];
    b.cen_pos[2] = Pars::z_body_cen[part];

	// internal variables declaration
	double lenPanel, sizePanel, pos_max, pos_min;
	int numPanel, numPanel2, N, startID;

	// assigning internal variables
    numPanel = Pars::Lref[part]/(2.5 * Pars::sigma);        // Said that the panel having size of particle core
	
	// calculating internal variables
    numPanel2 = numPanel * numPanel;
	N = 6 * numPanel2;
	lenPanel = Pars::Lref[part]/numPanel;
    sizePanel = lenPanel*lenPanel;
	pos_max = Pars::Lref[part]/2.0e0;	// maximum vertex coordinate
	pos_min = -Pars::Lref[part]/2.0e0;	// minimum vertex coordinate

	// resize x and y coordinate of body object and body initialization
	b.x.clear(); b.x.resize(N);
	b.y.clear(); b.y.resize(N);
    b.z.clear(); b.z.resize(N);
	b.x_m.clear(); b.x_m.resize(N);
	b.y_m.clear(); b.y_m.resize(N);
    b.z_m.clear(); b.z_m.resize(N);
	b.x_n.clear(); b.x_n.resize(N);
	b.y_n.clear(); b.y_n.resize(N);
    b.z_n.clear(); b.z_n.resize(N);
	b.size.clear(); b.size.resize(N);
	b.n_node = N;
	b.n_panel = N;

	// surface body point coordinate, panel midpoint coordinate, and panel normal vector
	startID = 0;		// The index of body node
    #pragma omp parallel for
	for (int i = 0; i < numPanel; i++){         // SURFACE 1: Upstream surface (normal x-)
        for (int j = 0; j < numPanel; j++){
            int id = startID + i*numPanel + j;
            
            b.x[id] = b.cen_pos[0] + pos_min;
            b.y[id] = b.cen_pos[1] + pos_min + (i+0.5)*lenPanel;
            b.z[id] = b.cen_pos[2] + pos_min + (j+0.5)*lenPanel;

            b.x_m[id] = b.x[id];
            b.y_m[id] = b.y[id];
            b.z_m[id] = b.z[id];
            
            b.x_n[id] = -1; // Face front
            b.y_n[id] = 0;
            b.z_n[id] = 0;

            b.size[id] = sizePanel;
        }
	}
    startID += numPanel2;
    #pragma omp parallel for
    for (int i = 0; i < numPanel; i++){         // SURFACE 2: Downstream surface (normal x+)
        for (int j = 0; j < numPanel; j++){
            int id = startID + i*numPanel + j;

            b.x[id] = b.cen_pos[0] + pos_max;
            b.y[id] = b.cen_pos[1] + pos_min + (i+0.5)*lenPanel;
            b.z[id] = b.cen_pos[2] + pos_min + (j+0.5)*lenPanel;

            b.x_m[id] = b.x[id];
            b.y_m[id] = b.y[id];
            b.z_m[id] = b.z[id];
            
            b.x_n[id] = 1;  // Face front
            b.y_n[id] = 0;
            b.z_n[id] = 0;

            b.size[id] = sizePanel;
            id++;
        }
	}
    startID += numPanel2;
    #pragma omp parallel for
    for (int i = 0; i < numPanel; i++){         // SURFACE 3: Top surface (normal y+)
        for (int j = 0; j < numPanel; j++){
            int id = startID + i*numPanel + j;

            b.x[id] = b.cen_pos[0] + pos_min + (i+0.5)*lenPanel;
            b.y[id] = b.cen_pos[1] + pos_max;
            b.z[id] = b.cen_pos[2] + pos_min + (j+0.5)*lenPanel;

            b.x_m[id] = b.x[id];
            b.y_m[id] = b.y[id];
            b.z_m[id] = b.z[id];
            
            b.x_n[id] = 0;
            b.y_n[id] = 1;  // Face front
            b.z_n[id] = 0;

            b.size[id] = sizePanel;
        }
	}
    startID += numPanel2;
    #pragma omp parallel for
    for (int i = 0; i < numPanel; i++){         // SURFACE 4: Bottom surface (normal y-)
        for (int j = 0; j < numPanel; j++){
            int id = startID + i*numPanel + j;

            b.x[id] = b.cen_pos[0] + pos_min + (i+0.5)*lenPanel;
            b.y[id] = b.cen_pos[1] + pos_min;
            b.z[id] = b.cen_pos[2] + pos_min + (j+0.5)*lenPanel;

            b.x_m[id] = b.x[id];
            b.y_m[id] = b.y[id];
            b.z_m[id] = b.z[id];
            
            b.x_n[id] = 0;
            b.y_n[id] = -1;  // Face front
            b.z_n[id] = 0;

            b.size[id] = sizePanel;
        }
	}
    startID += numPanel2;
    #pragma omp parallel for
    for (int i = 0; i < numPanel; i++){         // SURFACE 5: Side right surface (back face) (normal z-)
        for (int j = 0; j < numPanel; j++){
            int id = startID + i*numPanel + j;

            b.x[id] = b.cen_pos[0] + pos_min + (j+0.5)*lenPanel;
            b.y[id] = b.cen_pos[1] + pos_min + (i+0.5)*lenPanel;
            b.z[id] = b.cen_pos[2] + pos_min;

            b.x_m[id] = b.x[id];
            b.y_m[id] = b.y[id];
            b.z_m[id] = b.z[id];
            
            b.x_n[id] = 0;
            b.y_n[id] = 0;
            b.z_n[id] = -1; // Face front

            b.size[id] = sizePanel;
        }
	}
    startID += numPanel2;
    #pragma omp parallel for
    for (int i = 0; i < numPanel; i++){         // SURFACE 6: Side left surface (front face) (normal z+)
        for (int j = 0; j < numPanel; j++){
            int id = startID + i*numPanel + j;

            b.x[id] = b.cen_pos[0] + pos_min + (j+0.5)*lenPanel;
            b.y[id] = b.cen_pos[1] + pos_min + (i+0.5)*lenPanel;
            b.z[id] = b.cen_pos[2] + pos_max;

            b.x_m[id] = b.x[id];
            b.y_m[id] = b.y[id];
            b.z_m[id] = b.z[id];
            
            b.x_n[id] = 0;
            b.y_n[id] = 0;
            b.z_n[id] = 1;  // Face front

            b.size[id] = sizePanel;
        }
	}

	// Assign the extremes point
	basis_loop(d){
		b.min_pos[d] = b.cen_pos[d] + pos_min;
		b.max_pos[d] = b.cen_pos[d] + pos_max;
	}
}

/**
 *  @brief A 3D normal plate body geometry generator [OPTION 3]. [Not created yet]
 * 
 *  @param	_body	The body container to be generated.
 *  @param	part	The body part ID.
 *  
 *  @headerfile geometry.hpp
 */
void geometry::normal_plate_3d_generator(Body &b, int part)
{
	return;
}

/**
 *  @brief A 3D flat plate body geometry generator [OPTION 4]. [Not created yet]
 * 
 *  @param	_body	The body container to be generated.
 *  @param	part	The body part ID.
 *  
 *  @headerfile geometry.hpp
 */
void geometry::flat_plate_3d_generator(Body &b, int part)
{
	return;
}

/**
 *  @brief A 3D torus body geometry generator [OPTION 5].
 * 
 *  @param	_body	The body container to be generated.
 *  @param	part	The body part ID.
 *  
 *  @headerfile geometry.hpp
 */
void geometry::torus_generator(Body &b, int part)
{
	// Internal variable
    double scale = Pars::Lref[part];
    
    // Update body center position
	b.cen_pos[0] = Pars::x_body_cen[part];
	b.cen_pos[1] = Pars::y_body_cen[part];
    b.cen_pos[2] = Pars::z_body_cen[part];
    
    // **Read the data from file
    this->read_3D_geometry(b, "torus_0.4");

    // Update the node data
    b.n_node = b.n_panel;
    b.x.clear(); b.x.resize(b.n_node);
	b.y.clear(); b.y.resize(b.n_node);
    b.z.clear(); b.z.resize(b.n_node);

    // **Transform the body geometry
    // SCALING
    // MESSAGE_LOG << "Transform scaling ... \n";
    #pragma omp parallel for
    for (int i = 0; i < b.n_panel; i++){
        b.x_m[i] *= scale;
        b.y_m[i] *= scale;
        b.z_m[i] *= scale;
        b.size[i] *= (scale*scale);
    }
    
    // ROTATION
    if (Pars::bodyRotFlag[part]){
        // MESSAGE_LOG << "Transform rotation ... \n";
        // Create the rotation matrix
        std::vector<std::vector<double>> R(3, std::vector<double>(3,0));   // Rotation matrix
        // Convention [Row, Col]

        double _alp = Pars::bodyRotDeg[part] * M_PI / 180;
        
        /** Rotation matrix CCW toward the axis
         *    [cos -sin]
         *    [sin  cos]
        */
        switch (Pars::bodyRotAxis[part]){
        case 0: // Rotation toward x axis
            // Set fixed part (x)
            R[0][0] = 1;

            // Set rotation part
            R[1][1] = cos(_alp);    // y by y
            R[1][2] = -sin(_alp);   // y by z

            R[2][1] = sin(_alp);    // z by y
            R[2][2] = cos(_alp);    // z by z
            break;
        case 1: // Rotation toward y axis
            // Set fixed axis (y)
            R[1][1] = 1;

            // Set rotation part
            R[2][2] = cos(_alp);    // z by z
            R[2][0] = -sin(_alp);   // z by x

            R[0][2] = sin(_alp);    // x by z
            R[0][0] = cos(_alp);    // x by x

            break;
        case 2: // Rotation toward z axis
            // Set fixed axis (z)
            R[2][2] = 1;

            // Set rotation part
            R[0][0] = cos(_alp);    // x by x
            R[0][1] = -sin(_alp);   // x by y

            R[1][0] = sin(_alp);    // y by x
            R[1][1] = cos(_alp);    // y by y

            break;
        default:
            break;
        }

        // Perform rotation to each panel
        #pragma omp parallel for
        for (int i = 0; i < b.n_panel; i++){
            double _temPos[3] = {b.x_m[i], b.y_m[i], b.z_m[i]};
            double _temNorm[3] = {b.x_n[i], b.y_n[i], b.z_n[i]};
            
            // Rotate coordinate
            b.x_m[i] = R[0][0]*_temPos[0] + R[0][1]*_temPos[1] + R[0][2]*_temPos[2];
            b.y_m[i] = R[1][0]*_temPos[0] + R[1][1]*_temPos[1] + R[1][2]*_temPos[2];
            b.z_m[i] = R[2][0]*_temPos[0] + R[2][1]*_temPos[1] + R[2][2]*_temPos[2];

            // Rotate normal
            b.x_n[i] = R[0][0]*_temNorm[0] + R[0][1]*_temNorm[1] + R[0][2]*_temNorm[2];
            b.y_n[i] = R[1][0]*_temNorm[0] + R[1][1]*_temNorm[1] + R[1][2]*_temNorm[2];
            b.z_n[i] = R[2][0]*_temNorm[0] + R[2][1]*_temNorm[1] + R[2][2]*_temNorm[2];
        }
    }

    // TRANSLATION (translate to the new center)
    // MESSAGE_LOG << "Transform translation ... \n";
    #pragma omp parallel for
    for (int i = 0; i < b.n_panel; i++){
        b.x_m[i] += b.cen_pos[0];
        b.y_m[i] += b.cen_pos[1];
        b.z_m[i] += b.cen_pos[2];
        b.x[i] = b.x_m[i];
        b.y[i] = b.y_m[i];
        b.z[i] = b.z_m[i];
    }

    // MESSAGE_LOG << "Check extremes ... \n";
    // Evaluate the body extreme location
    b.min_pos[0] = *std::min_element(b.x_m.begin(), b.x_m.end());
    b.min_pos[1] = *std::min_element(b.y_m.begin(), b.y_m.end());
    b.min_pos[2] = *std::min_element(b.z_m.begin(), b.z_m.end());

    b.max_pos[0] = *std::max_element(b.x_m.begin(), b.x_m.end());
    b.max_pos[1] = *std::max_element(b.y_m.begin(), b.y_m.end());
    b.max_pos[2] = *std::max_element(b.z_m.begin(), b.z_m.end());

    // MESSAGE_LOG << "Done all !!! \n";

    return;
}

/**
 *  @brief A 3D heart body geometry generator [OPTION 6].
 * 
 *  @param	_body	The body container to be generated.
 *  @param	part	The body part ID.
 *  
 *  @headerfile geometry.hpp
 */
void geometry::heart_generator(Body &b, int part)
{
	// Internal variable
    double scale = Pars::Lref[part];
    
    // Update body center position
	b.cen_pos[0] = Pars::x_body_cen[part];
	b.cen_pos[1] = Pars::y_body_cen[part];
    b.cen_pos[2] = Pars::z_body_cen[part];
    
    // **Read the data from file
    this->read_3D_geometry(b, "heart");

    // Update the node data
    b.n_node = b.n_panel;
    b.x.clear(); b.x.resize(b.n_node);
	b.y.clear(); b.y.resize(b.n_node);
    b.z.clear(); b.z.resize(b.n_node);

    // **Transform the body geometry
    // SCALING
    #pragma omp parallel for
    for (int i = 0; i < b.n_panel; i++){
        b.x_m[i] *= scale;
        b.y_m[i] *= scale;
        b.z_m[i] *= scale;
        b.size[i] *= (scale*scale);
    }
    
    // ROTATION
    if (Pars::bodyRotFlag[part]){
        // Create the rotation matrix
        std::vector<std::vector<double>> R(3, std::vector<double>(3,0));   // Rotation matrix
        // Convention [Row, Col]

        double _alp = Pars::bodyRotDeg[part] * M_PI / 180;
        
        /** Rotation matrix CCW toward the axis
         *    [cos -sin]
         *    [sin  cos]
        */
        switch (Pars::bodyRotAxis[part]){
        case 0: // Rotation toward x axis
            // Set fixed part (x)
            R[0][0] = 1;

            // Set rotation part
            R[1][1] = cos(_alp);    // y by y
            R[1][2] = -sin(_alp);   // y by z

            R[2][1] = sin(_alp);    // z by y
            R[2][2] = cos(_alp);    // z by z
            break;
        case 1: // Rotation toward y axis
            // Set fixed axis (y)
            R[1][1] = 1;

            // Set rotation part
            R[2][2] = cos(_alp);    // z by z
            R[2][0] = -sin(_alp);   // z by x

            R[0][2] = sin(_alp);    // x by z
            R[0][0] = cos(_alp);    // x by x

            break;
        case 2: // Rotation toward z axis
            // Set fixed axis (z)
            R[2][2] = 1;

            // Set rotation part
            R[0][0] = cos(_alp);    // x by x
            R[0][1] = -sin(_alp);   // x by y

            R[1][0] = sin(_alp);    // y by x
            R[1][1] = cos(_alp);    // y by y

            break;
        default:
            break;
        }

        // Perform rotation to each panel
        #pragma omp parallel for
        for (int i = 0; i < b.n_panel; i++){
            double _temPos[3] = {b.x_m[i], b.y_m[i], b.z_m[i]};
            double _temNorm[3] = {b.x_n[i], b.y_n[i], b.z_n[i]};
            
            // Rotate coordinate
            b.x_m[i] = R[0][0]*_temPos[0] + R[0][1]*_temPos[1] + R[0][2]*_temPos[2];
            b.y_m[i] = R[1][0]*_temPos[0] + R[1][1]*_temPos[1] + R[1][2]*_temPos[2];
            b.z_m[i] = R[2][0]*_temPos[0] + R[2][1]*_temPos[1] + R[2][2]*_temPos[2];

            // Rotate normal
            b.x_n[i] = R[0][0]*_temNorm[0] + R[0][1]*_temNorm[1] + R[0][2]*_temNorm[2];
            b.y_n[i] = R[1][0]*_temNorm[0] + R[1][1]*_temNorm[1] + R[1][2]*_temNorm[2];
            b.z_n[i] = R[2][0]*_temNorm[0] + R[2][1]*_temNorm[1] + R[2][2]*_temNorm[2];
        }
    }

    // TRANSLATION (translate to the new center)
    #pragma omp parallel for
    for (int i = 0; i < b.n_panel; i++){
        b.x_m[i] += b.cen_pos[0];
        b.y_m[i] += b.cen_pos[1];
        b.z_m[i] += b.cen_pos[2];
        b.x[i] = b.x_m[i];
        b.y[i] = b.y_m[i];
        b.z[i] = b.z_m[i];
    }

    // Evaluate the body extreme location
    b.min_pos[0] = *std::min_element(b.x_m.begin(), b.x_m.end());
    b.min_pos[1] = *std::min_element(b.y_m.begin(), b.y_m.end());
    b.min_pos[2] = *std::min_element(b.z_m.begin(), b.z_m.end());

    b.max_pos[0] = *std::max_element(b.x_m.begin(), b.x_m.end());
    b.max_pos[1] = *std::max_element(b.y_m.begin(), b.y_m.end());
    b.max_pos[2] = *std::max_element(b.z_m.begin(), b.z_m.end());

    return;
}

/**
 *  @brief A 3D body geometry reader.
 * 
 *  @param	_body	The body container for reading target.
 *  @param	bodyName	The body name written on the file name.
 *  @param	part	The body part ID.
 *  
 *  @headerfile geometry.hpp
 */
void geometry::read_3D_geometry(Body &b, const std::string &bodyName){
    // Reset the data
    b.n_panel = 0;
    b.x_m.clear();
    b.y_m.clear();
    b.z_m.clear();
    b.x_n.clear();
    b.y_n.clear();
    b.z_n.clear();
    b.size.clear();

    /** Procedure
     *   > Read geometry data from file
     *   > Assign to body data
    */

    // Construct the file directory name
    std::string fileName = "input/geometry_data/geometry_3D_" + bodyName + ".csv";

    // Internal variable
    std::vector<double> dataList;
    std::string line;
    std::string entry;
    // Data storage
    double xm,ym,zm;
    double xn,yn,zn;
    double size;
    
    // Read the data file
    std::ifstream reader;
    reader.open(fileName);
    if(!reader.is_open()){
        ERROR_LOG << "The body data of " << bodyName << " is not existed!\n";
        throw std::runtime_error("Error: The intended file is not existed!");
    }
    getline(reader, line);  // Put out the first line (header) before reading
    
    // Read all data line by line
    while(reader.peek()!=EOF){
        getline(reader, line);      // Take the data line
        dataList.clear();           // Reset the list
        entry = "";                 // Reset the entry

        // The sequence of the location index 
        // ID, x, y, z, nx, ny, nz, s
        // 0   1  2  3   4   5   6  7

        // Get data per data in line, put into the data list
        for (const auto &dig : line){
            // Take the entry as the value
            if (dig == ','){
                double value = std::stod(entry);
                dataList.push_back(value);
                entry = "";
                continue;
            }

            // Add the digit into the entry
            entry += dig;
        }
        // Take the last data
        dataList.push_back(std::stod(entry));

        // Take out the data from list
        for (size_t i = 0; i < dataList.size(); i++){
            // Aliasing to the value
            const auto &value = dataList[i];
            int loc = i;
            
            // Put the data into particle
            switch (loc){
                case 1: xm = value; break;
                case 2: ym = value; break;
                case 3: zm = value; break;
                case 4: xn = value; break;
                case 5: yn = value; break;
                case 6: zn = value; break;
                case 7: size = value; break;
                default: break;
            }
        }

        // Put the data into Perform the geometry input data
        b.n_panel++;
        b.x_m.push_back(xm);
        b.y_m.push_back(ym);
        b.z_m.push_back(zm);
        b.x_n.push_back(xn);
        b.y_n.push_back(yn);
        b.z_n.push_back(zn);
        b.size.push_back(size);
    }
    reader.close();    
    
    return;
}

