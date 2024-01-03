#include "inter_search.hpp"

// Macro for vector
template <typename U> using vec = std::vector<U>;

/**
 *  @brief  Evaluate the neighboring interaction between two set of particle data. 
 *  Interaction pairs are determined by using a sorting grid spatial hash.
 * 
 *  @param  _sourcePar  Data of source particle.
 *  @param  _targetPar  Data of particle which neighbor toward source to be found.
 *  @param  _nghList  [OUTPUT] The neighbor ID list data.
*/
void InterSearchNgh::find_neighbor(const Particle &p_src, const Particle &p_tar, vec<vec<int>> &_nghIDList)
{
	// Release the neighbor list (to be fill out further)
    _nghIDList.clear();
	
	// Set the internal variable
	int srcNum = p_src.num;
	int tarNum = p_tar.num;
	NghBaseGrid _tempGrid;

	// **[1] Generate the temporary grid
	double maxCoor[DIM];
	double minCoor[DIM];
	// Evaluate extreme in x direction
		minCoor[0] = std::min<double>(*(std::min_element(p_src.x.begin(),p_src.x.end())), *(std::min_element(p_tar.x.begin(),p_tar.x.end())));
		maxCoor[0] = std::max<double>(*(std::max_element(p_src.x.begin(),p_src.x.end())), *(std::max_element(p_tar.x.begin(),p_tar.x.end())));
	// Evaluate extreme in y direction
		minCoor[1] = std::min<double>(*(std::min_element(p_src.y.begin(),p_src.y.end())), *(std::min_element(p_tar.y.begin(),p_tar.y.end())));
		maxCoor[1] = std::max<double>(*(std::max_element(p_src.y.begin(),p_src.y.end())), *(std::max_element(p_tar.y.begin(),p_tar.y.end())));
	// Evaluate extreme in z direction
	if (DIM > 2){
		minCoor[2] = std::min<double>(*(std::min_element(p_src.z.begin(),p_src.z.end())), *(std::min_element(p_tar.z.begin(),p_tar.z.end())));
		maxCoor[2] = std::max<double>(*(std::max_element(p_src.z.begin(),p_src.z.end())), *(std::max_element(p_tar.z.begin(),p_tar.z.end())));
	}
    // Perform grid initialization
    _tempGrid.initialize_grid(maxCoor, minCoor, p_tar.s);


	// **[2] Generate the spatial hashing
	int gridNum = _tempGrid.get_count();	// Get the number of grid
    // Resize the spatial hashing parameter
    this->pntIDListSrc.resize(gridNum);

    // Iterate through all source point
    double _coor[DIM];          // Temporary point coordinate
    for (int i = 0; i < srcNum; i++){
        // Aliasing point ID
        const int &_ID = i;
        // Assign the coordinate location at current point ID
        _coor[0] = p_src.x[_ID];
        _coor[1] = p_src.y[_ID];
        if (DIM > 2) _coor[2] = p_src.z[_ID];
        // Get the grid ID
        int gridID = _tempGrid.get_grid_ID(_coor);

        // Update the hashing class member
		this->pntIDListSrc[gridID].push_back(_ID);	// Collect the source point ID into the group container
    }
	
	// **[3] Evaluate neighbor interation (spatial hash)
    _nghIDList.resize(tarNum);

	// Iterate through all target point
	#pragma omp parallel for
	for (int ID_tar = 0; ID_tar < tarNum; ID_tar++){
		// Assign the coordinate location at current point ID
        _coor[0] = p_tar.x[ID_tar];
        _coor[1] = p_tar.y[ID_tar];
        if (DIM > 2) _coor[2] = p_tar.z[ID_tar];
        // Get the grid ID where the current target point is located
        int gridID = _tempGrid.get_grid_ID(_coor);
		
		// Get the list of neighbor grid of the current grid ID
		vec<int> _tempGridNghIDList;		// Temporary neighbor grid
        _tempGrid.find_grid_ngh(_tempGridNghIDList, gridID);		

		// Evaluate to each point inside the grid neighbor
		for (const auto &gridNghID : _tempGridNghIDList){
			for (const auto &ID_src : pntIDListSrc[gridNghID]){
				// Calculate the distance square between two points
				double _dr2 = 0;
				// Calculate in x direction
					double _dx = p_tar.x[ID_tar] - p_src.x[ID_src];
					_dr2 += _dx*_dx;
				// Calculate in y direction
					double _dy = p_tar.y[ID_tar] - p_src.y[ID_src];
					_dr2 += _dy*_dy;
				// Calculate in z direction
				if (DIM > 2){
					double _dz = p_tar.z[ID_tar] - p_src.z[ID_src];
					_dr2 += _dz*_dz;
				}

				// Evaluate distance
				if (Pars::opt_ngh_interact == 1){
					// Check the particle i
					double _rSup = p_tar.s[ID_tar] * Pars::r_sup * Pars::r_buff;   // Support radius of particle i
					if (_dr2 < (_rSup*_rSup)){
						_nghIDList[ID_tar].push_back(ID_src);
					}
				}
				else if (Pars::opt_ngh_interact == 2){
					// Support radius of average size
					double _rSupAve = ((p_tar.s[ID_tar] + p_src.s[ID_src]) / 2.0e0) * Pars::r_sup * Pars::r_buff;
					if (_dr2 < (_rSupAve*_rSupAve)){
						_nghIDList[ID_tar].push_back(ID_src);
					}
				}
			}
		}
	}

	return;
}