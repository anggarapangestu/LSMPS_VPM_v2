#include "spatial_hash.hpp"

// Macro for vector
template <typename U> using vec = std::vector<U>;

/**
 *  @brief  Evaluate the neighboring interaction for each particle. 
 *  Interaction pairs are determined by using a sorting grid spatial hash
 *  the complexity of the spatial hash algorithm is of order O(N).
 * 
 *  @param  _nghIDList  [OUTPUT] Neighbor ID list.
 *  @param  _baseGrid   A point data grouping tools for neighbor evaluation.
 *  @param  _size  The list of point size.
 *  @param  _xp  The list of point x coordinate.
 *  @param  _yp  The list of point y coordinate.
 *  @param  _zp  The list of point z coordinate (Just put a dummy variable
 *               for DIM < 3).
*/
void SpatialHashNgh::find_neighbor(vec<vec<int>> &_nghIDList, NghBaseGrid &_baseGrid,
									const vec<double> &_sp,
									const vec<double> &_xp,
									const vec<double> &_yp,
									const vec<double> &_zp
){
	// Release the neighbor list (to be fill out further)
    _nghIDList.clear();

	// **Generate the spatial hashing
    int pntNum = _sp.size();     			// Get the number of point data
	int gridNum = _baseGrid.get_count();	// Get the number of grid
    // Resize the spatial hashing parameter
    this->pntIDList.resize(gridNum);
	// this->evalFlag.resize(pntNum, false);

    // Iterate through all point
    double _coor[DIM];          // Temporary point coordinate
    for (int i = 0; i < pntNum; i++){
        // Aliasing point ID
        const int &_ID = i;
        // Assign the coordinate location at current point ID
        _coor[0] = _xp[_ID];
        if (DIM > 1) _coor[1] = _yp[_ID];
        if (DIM > 2) _coor[2] = _zp[_ID];
        // Get the grid ID
        int gridID = _baseGrid.get_grid_ID(_coor);

        // Update the hashing class member
		this->pntIDList[gridID].push_back(_ID);	// Collect the point ID into the group container
    }
	
	// **Evaluate neighbor interation (spatial hash)
    _nghIDList.resize(pntNum);

	// Iterate through all grid
	for (int gridID = 0; gridID < gridNum; gridID++){
		// Get the list of neighbor grid of the current grid ID
        vec<int> _gridNghIDList;
        _baseGrid.find_grid_ngh(_gridNghIDList, gridID);

		// Iterate through all point inside the grid
		#pragma omp parallel for
		for (size_t i = 0; i < pntIDList.at(gridID).size(); i++){
			// Aliasing the current particle ID
			const int &ID_i = pntIDList[gridID][i];
			
			// // Update the flag of the current particle
			// this->evalFlag[ID_i] = true;

			// Evaluate to each point inside the grid neighbor
			for (const auto &gridNghID : _gridNghIDList){
				for (const auto &ID_j : pntIDList[gridNghID]){
					// // Exception for not evaluating twice
					// if (this->evalFlag[ID_j]) continue;

					// An exception for include itself on neighbor list
					if (!Pars::flag_ngh_include_self && (ID_j == ID_i)) continue;

					// Calculate the distance square between two points
					double _dr2 = 0;
					// Calculate in x direction
						double _dx = _xp[ID_i] - _xp[ID_j];
						_dr2 += (_dx*_dx);
					// Calculate in y direction
						double _dy = _yp[ID_i] - _yp[ID_j];
						_dr2 += (_dy*_dy);
					// Calculate in z direction
					if (DIM > 2){
						double _dz = _zp[ID_i] - _zp[ID_j];
						_dr2 += (_dz*_dz);
					}

					// Evaluate distance
					if (Pars::opt_ngh_interact == 1){
						// Check the particle i
						double _rSup = _sp[ID_i] * Pars::r_sup * Pars::r_buff;   // Support radius of particle i
						if (_dr2 < (_rSup*_rSup)){
							_nghIDList[ID_i].push_back(ID_j);
						}

						// // Check the particle j
						// _rSup = _sp[ID_j] * Pars::r_sup * Pars::r_buff;   // Support radius of particle j
						// if (_dr2 < (_rSup*_rSup)){
						// 	_nghIDList[ID_j].push_back(ID_i);
						// }
					}
					else if (Pars::opt_ngh_interact == 2){
						// Support radius of average size
						double _rSupAve = ((_sp[ID_i] + _sp[ID_j]) / 2.0e0) * Pars::r_sup * Pars::r_buff;
						if (_dr2 < (_rSupAve*_rSupAve)){
							_nghIDList[ID_i].push_back(ID_j);
							// _nghIDList[ID_j].push_back(ID_i);
						}
					}
				}
			}

			// // Include self or not
			// if (Pars::flag_ngh_include_self){
			// 	_nghIDList[ID_i].push_back(ID_i);
			// }
		}
	}
}