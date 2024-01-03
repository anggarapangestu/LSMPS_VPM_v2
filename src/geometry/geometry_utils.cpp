#include "geometry.hpp"

#include "../penalization/penalization.hpp"

#define BODY_OUT_EXT_MUL 2.0
#define BODY_INSIDE_EXT_MUL 5.0

/**
 *  @brief Finding the minimum distance from body at one single point.
 * 	The minimum distance is relative to panel normal direction (sign:
 *  [+] Outside the body, [-] inside body). While, minimum absolute 
 *  distance is the shortest distance toward the node.
 *  
 *  @param	_coor	The point coordinate of evaluation.
 * 	@param	_body	The evaluated body data.
 *  @param	_abs	Flag for absolute distance calculation.
 * 
 *  @headerfile geometry.hpp
 */
double geometry::distance_calc(const std::vector<double> &pos, const Body& _body, const bool _abs){
	// Internal variable
	int panelNode;			// Panel node ID correspond to minimum distance
	double minDist_2;		// Minimum distance squared of particle to panel middle node
	double _dist_2;			// Temporarty distance squared calculation
	double result;			// The normal minimun distance (sign, [-] outside body, [+] inside body)
	double panelPos[DIM];	// Temporary body panel mid coordinate

	// Set up the minimum distance (check toward panel at ID 0)
	panelNode = 0;
	minDist_2 = 0.0;
	panelPos[0] = _body.x_m[0];
	if (DIM > 1) panelPos[1] = _body.y_m[0];
	if (DIM > 2) panelPos[2] = _body.z_m[0];
	basis_loop(d) minDist_2 += std::pow(pos[d] - panelPos[d],2.0);

	// Only check the R^2, no need to find the sqrt(R^2)
	for (int _i = 1; _i < _body.n_panel; _i++){
		// Take the current panel mid coordinate
		panelPos[0] = _body.x_m[_i];
		if (DIM > 1) panelPos[1] = _body.y_m[_i];
		if (DIM > 2) panelPos[2] = _body.z_m[_i];

		// Calculate the distance
		_dist_2 = 0.0;
		basis_loop(d) _dist_2 += std::pow(pos[d] - panelPos[d],2.0);

		// Update the minimum distance
		if (_dist_2 < minDist_2){
			panelNode = _i;
			minDist_2 = _dist_2;
		}
	}

	// Update the distance value
	if (_abs){	// Put the absolute distance: square root of minDist_2
		result = std::sqrt(minDist_2);
	}else{ 		// Normal distance: Inner product between minDist and normal direction
		result = (pos[0] - _body.x_m[panelNode])*_body.x_n[panelNode];
		if(DIM > 1) result += (pos[1] - _body.y_m[panelNode])*_body.y_n[panelNode];
		if(DIM > 2) result += (pos[2] - _body.z_m[panelNode])*_body.z_n[panelNode];
	}

	return result;
}

/**
 *  @brief Finding the minimum distance from body at all particle point.
 * 	The minimum distance is relative to panel normal direction (sign:
 *  [+] Outside the body, [-] inside body). While, minimum absolute 
 *  distance is the shortest distance toward the node.
 *  
 *  @param	_par	The particle data to be evaluated.
 * 	@param	_body	The evaluated body data.
 *  @param	_abs	Flag for absolute distance calculation.
 * 
 *  @headerfile geometry.hpp
 */
void geometry::distance_calc(Particle &_par, const Body &_body, const bool _abs){
	// The calculated variable
	double minDist_2;			// Minimum distance of particle to panel middle node
	double _dist_2;			// Temporarty distance calculation
	int panelNode;			// The panel node correspond to minimum distance
	double panelPos[DIM];	// Temporary panel node coordinate position
	std::vector<double> parPos(DIM);	// Temporary particle coordinate
	
	// Initialize the particle distance container
	_par.R.clear();
	_par.R.resize(_par.num,0.0);

	// Iterate through all particle inside the domain
	#pragma omp parallel for
	for (int ID = 0; ID < _par.num; ID++){
		// Take the current particle coordinate
		parPos[0] = _par.x[ID];
		if (DIM > 1) parPos[1] = _par.y[ID];
		if (DIM > 2) parPos[2] = _par.z[ID];

		// Initialize the minimun distance by comparing toward body panel at ID 0
		panelNode = 0;
		minDist_2 = 0.0;
		panelPos[0] = _body.x_m[0];
		if (DIM > 1) panelPos[1] = _body.y_m[0];
		if (DIM > 2) panelPos[2] = _body.z_m[0];
		basis_loop(d) minDist_2 += std::pow(parPos[d] - panelPos[d],2.0);
		
		// Evaluate through all body panel
		// Only check the R^2, no need to find the sqrt(R^2)
		for (int _i = 1; _i < _body.n_panel; _i++){
			// Take the current body panel mid coordinate
			panelPos[0] = _body.x_m[_i];
			if (DIM > 1) panelPos[1] = _body.y_m[_i];
			if (DIM > 2) panelPos[2] = _body.z_m[_i];

			// Calculate the distance
			_dist_2 = 0.0;
			basis_loop(d) _dist_2 += std::pow(parPos[d] - panelPos[d],2.0);

			// Update the minimum distance
			if (_dist_2 < minDist_2){
				panelNode = _i;
				minDist_2 = _dist_2;
			}
		}
		
		// Update the distance value
		if (_abs){	// Put the absolute distance: square root of minDist_2
			_par.R[ID] = std::sqrt(minDist_2);
		}else{ 		// Normal distance: Inner product between minDist and normal direction
			_par.R[ID] = (_par.x[ID] - _body.x_m[panelNode])*_body.x_n[panelNode];
			if(DIM > 1) _par.R[ID] += (_par.y[ID] - _body.y_m[panelNode])*_body.y_n[panelNode];
			if(DIM > 2) _par.R[ID] += (_par.z[ID] - _body.z_m[panelNode])*_body.z_n[panelNode];
		}
	}

}

/**
 *  @brief  Evaluate the nearest body part of each particle in the domain.
 *         
 *  @param  _particle  Particle data container.
 *  @param  _bodyList  The list of body data.
*/
void geometry::eval_near_body(Particle &_par, const std::vector<Body> &_bodyList){
    // Initialize the container
    _par.bodyPart.resize(_par.num,-1);

	// Body extension value
	double body_extension = BODY_OUT_EXT_MUL * Pars::body_ext;

    // Check throughout all particle 
    #pragma omp parallel for
    for (int ID = 0; ID < _par.num; ID++){
        // Create the position array of particle
        std::vector<double> _pos(DIM);
        _pos[0] = _par.x[ID];
        _pos[1] = _par.y[ID];
        if (DIM > 2)
        _pos[2] = _par.z[ID];

        // Check throughout all body to find all body near the particle
        std::vector<int> bodyPartList;
        for (size_t part = 0; part < _bodyList.size(); part++){
            // Create a body alias
            const Body &_body = _bodyList[part];

            // Check if the particle is near this body
            bool isNear = true;
            basis_loop(d){
                // Check if outside the body extreme range (*add by body extention)
                if (_pos[d] < _body.min_pos[d] - body_extension || 
                    _pos[d] > _body.max_pos[d] + body_extension){
                    isNear = false;
                    break;
                }
            }
            
            // Put the body ID into the list if meet the condition
            if (isNear)
            bodyPartList.push_back(part);
        }

        // Check through all body part that is near the particle
        if (bodyPartList.empty()){
            _par.bodyPart[ID] = -1;
        }else if(bodyPartList.size() == 1){
            _par.bodyPart[ID] = bodyPartList[0];
        }else{
            // Find the minimum distance between the nearest body
            double R, tempR;
            int bodyID = -1;
            for (auto &part : bodyPartList){
                const Body &_body = _bodyList[part];
                tempR = this->distance_calc(_pos, _body, true);
                if (bodyID == -1){
                    R = tempR;
                    bodyID = part;
                }
                else if (tempR < R){
                    R = tempR;
                    bodyID = part;
                }
            }
            _par.bodyPart[ID] = bodyID;
        }
    }
    
    return;
}


/**
 *  @brief  Evaluate all particle whether it is near the body surface.
 *  NOTE: Calculate the chi and assign the particle active flag.
 *         
 *  @param  _particle  Particle data container.
 *  @param  _bodyList  The list of body data.
*/
void geometry::eval_body_flag(Particle &_par, const std::vector<Body> &_bodyList){
    /* PROCEDURE
		> Check all near body particles
		> Calculate the minimum distance toward it's 
			corresponding near body part
		> Put true flag for particle with distance under
			surface body extention (Pars::body_ext)
	*/
	
	// Initialize the container
    _par.isNearSurface.resize(_par.num, false);
	_par.insideBody.resize(_par.num, false);

	// Here also calculate chi [Penalization Import]
	_par.chi.resize(_par.num, 0.0e0);

	// Body extension value
	double inner_extension = Pars::body_ext;
	double outer_extension = BODY_OUT_EXT_MUL * Pars::body_ext;

    // Check throughout all particle 
    #pragma omp parallel for
    for (int ID = 0; ID < _par.num; ID++){
        // Alias to the nearest body part
		const int &bodyPartID = _par.bodyPart[ID];

		// Only check near-body particle
		if (bodyPartID == -1) continue;
		
		// Create the position array of particle
        std::vector<double> _pos(DIM);
        _pos[0] = _par.x[ID];
        _pos[1] = _par.y[ID];
        if (DIM > 2)
        _pos[2] = _par.z[ID];

        // Alias to the nearest body
		const Body &_body = _bodyList[bodyPartID];

		// Check through all body part that is near the particle
		double _dist = this->distance_calc(_pos, _body, false);

		// Update near surface
		if ((_dist < outer_extension) && (_dist > -inner_extension)){
			_par.isNearSurface[ID] = true;
		}
		else{
			_par.isNearSurface[ID] = false;
		}

		// Update inside body
		if (_dist < -(BODY_INSIDE_EXT_MUL * Pars::sigma)){
			_par.insideBody[ID] = true;
		}
		else{
			_par.insideBody[ID] = false;
		}

		// Calculate chi [Penalization Import]
        penalization pen_tool;
		_par.chi[ID] = pen_tool.get_chi(_dist);
	}

	// Evaluate the active particle
    _par.isActive.resize(_par.num, false);
    #pragma omp parallel for
    for (int i = 0; i < _par.num; i++){
        // Only change the flag for true condition
		if (_par.isNearSurface[i] == true){
			_par.isActive[i] = true;
		}
    }
    
    return;
}
