#include "direct_find.hpp"

// Macro for vector
template <typename U> using vec = std::vector<U>;

/**
 *  @brief  Evaluate the neighboring interaction for each particle. 
 *  Interaction pairs are determined by directly comparing the 
 *  particle distance with the corresponding support radius length.
 * 
 *  @param  _nghIDList  [OUTPUT] Neighbor ID list.
 *  @param  _size  The list of point size.
 *  @param  _xp  The list of point x coordinate.
 *  @param  _yp  The list of point y coordinate.
 *  @param  _zp  The list of point z coordinate (Just put a dummy variable
 *               for DIM < 3).
*/
void directFindNgh::find_neighbor(vec<vec<int>> &_nghIDList, const vec<double> &_sp,
                                  const vec<double> &_xp,
                                  const vec<double> &_yp,
                                  const vec<double> &_zp)
{
    // Release the neighbor list (to be fill out further)
    _nghIDList.clear();

    // Internal variable
    int pntNum = _sp.size();    // Number of all point
    double _dist2;              // Temporary squared distance
    double _dx, _dy, _dz;       // Temporary coordinate basis different

    // Directly evaluate each particle toward each another O(N^2)
    _nghIDList.resize(pntNum);
    for (int ID_i = 0; ID_i < pntNum - 1; ID_i++){
        for (int ID_j = ID_i + 1; ID_j < pntNum; ID_j++){
            // Calculate the distance between particle [i] and [j]
            _dist2 = 0;
            // Calculate at x direction
                _dx = _xp[ID_i] - _xp[ID_j];
                _dist2 += _dx*_dx;
            // Calculate at y direction
                _dy = _yp[ID_i] - _yp[ID_j];
                _dist2 += _dy*_dy;
            // Calculate at z direction
            if (DIM > 2){
                _dz = _zp[ID_i] - _zp[ID_j];
                _dist2 += _dz*_dz;
            }

            // Evaluate distance
            if (Pars::opt_ngh_interact == 1){
                // Check the particle i
                double _rSup = _sp[ID_i] * Pars::r_sup * Pars::r_buff;   // Support radius of particle i
                if (_dist2 < (_rSup*_rSup)){
                    _nghIDList[ID_i].push_back(ID_j);
                }

                // Check the particle j
                _rSup = _sp[ID_j] * Pars::r_sup * Pars::r_buff;   // Support radius of particle j
                if (_dist2 < (_rSup*_rSup)){
                    _nghIDList[ID_j].push_back(ID_i);
                }
            }
            else if (Pars::opt_ngh_interact == 2){
                // Support radius of average size
                double _rSupAve = ((_sp[ID_i] + _sp[ID_j]) / 2.0e0) * Pars::r_sup * Pars::r_buff;
                if (_dist2 < (_rSupAve*_rSupAve)){
                    _nghIDList[ID_i].push_back(ID_j);
                    _nghIDList[ID_j].push_back(ID_i);
                }
            }
        }
        
        // Include self or not
        if (Pars::flag_ngh_include_self){
            _nghIDList[ID_i].push_back(ID_i);
        }
    }
    return;
}