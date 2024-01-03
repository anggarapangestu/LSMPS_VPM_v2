#include "link_list.hpp"

// Macro for vector
template <typename U> using vec = std::vector<U>;

// =====================================================
// +----------------- Link list Method -----------------+ 
// =====================================================

/**
 *  @brief  Evaluate the neighboring interaction for each particle. 
 *  Interaction pairs are determined by using a sorting grid linked list
 *  the complexity of the linked-list algorithm is of order O(N).
 * 
 *  @param  _nghIDList  [OUTPUT] Neighbor ID list.
 *  @param  _baseGrid   A point data grouping tools for neighbor evaluation.
 *  @param  _size  The list of point size.
 *  @param  _xp  The list of point x coordinate.
 *  @param  _yp  The list of point y coordinate.
 *  @param  _zp  The list of point z coordinate (Just put a dummy variable
 *               for DIM < 3).
*/
void LinkListNgh::find_neighbor(vec<vec<int>> &nghIDList, NghBaseGrid &_baseGrid,
    const vec<double> &sp, const vec<double> &xp,
    const vec<double> &yp, const vec<double> &zp
){
    // Release the neighbor list (to be fill out further)
    nghIDList.clear();

    // **Generate the link list data
    int pntNum = sp.size();     // Get the number of point data
    // Resize the link list parameter
    this->pntGridID.resize(pntNum,0);
    this->nextID.resize(pntNum,0);
    this->headID.resize(_baseGrid.get_count(),0);

    // Iterate through all point
    double _coor[DIM];          // Temporary point coordinate
    for (int i = 0; i < pntNum; i++){
        // Aliasing point ID
        const int &_ID = i;
        // Assign the coordinate location at current point ID
        _coor[0] = xp[_ID];
        if (DIM > 1) _coor[1] = yp[_ID];
        if (DIM > 2) _coor[2] = zp[_ID];
        // Get the grid ID
        int gridID = _baseGrid.get_grid_ID(_coor);

        // Update the link list parameter
        this->pntGridID[_ID] = gridID;
        this->nextID[_ID] = this->headID[gridID];   // The current head in this grid is the next ID (previous head)
        this->headID[gridID] = _ID;                 // Update the current point as new head
        
        /* Illustration of link list assignment
        *    ___________     Iteration will go from low to higher point ID
        *   |    [4] [1]|     * In this grid there are 5 point [1], [3], [4], [7], and [13]
        *   |  [7]      |     * At each point, take the old head as the next ID of current point
        *   | [13] [3]  |        then take the current point as new head
        *   |___________|     * Do until there is no point left
        *                 
        *   Iteration illustration 
        *   **********************     *Assign old node to curr ID*               
        *                   |  (old)  |  Next ID at eact point    | (new)        
        *    iter.  pnt ID  |  headID | [1]  [3]  [4]  [7]  [13]  | HeadID      * Final next ID pair
        *   ----------------------------------------------------------------    *    [1] -> 0 (end)
        *      0       1    |    0    | (0)   0    0    0     0   |  (1)        *    [3] -> 1
        *      1       3    |    1    |  0   (1)   0    0     0   |  (3)        *    [4] -> 3
        *      2       4    |    3    |  0    1   (3)   0     0   |  (4)        *    [7] -> 4
        *      3       7    |    4    |  0    1    3   (4)    0   |  (7)        *    [13] -> 7
        *      4       13   |    7    |  0    1    3    4    (7)  |  (13)       *    head
        *   
        *   The final link list
        *   *******************
        *      [13]->[7]->[4]-[3]->[1]->[0]
        *       ^                        ^
        *      head                     end (must be 0)
        */
    }

    // **Evaluate neighbor interation (Link List)
    nghIDList.resize(pntNum);
    for (int i = 0; i < pntNum; i++){
        int &ID_i = i;      // Alias to the particle i

        // Get the list of neighbor grid
        vec<int> _tempNghIDList;
        _baseGrid.find_grid_ngh(_tempNghIDList, this->pntGridID[ID_i]);

        // Iterate through all neighbor grid
        for (const auto &gridNghID : _tempNghIDList){   
            // Find the link list head at the current grid            
            int ID_j = this->headID[gridNghID];

            // Evaluate neighbor if the current head ID (bigger than) current point ID
            while (ID_j > ID_i){
                // Calculate the distance square between two points
                double _dr2 = 0;
                // Calculate in x direction
                    double _dx = xp[ID_i] - xp[ID_j];
                    _dr2 += std::pow(_dx, 2.0);
                // Calculate in y direction
                    double _dy = yp[ID_i] - yp[ID_j];
                    _dr2 += std::pow(_dy, 2.0);
                // Calculate in z direction
                if (DIM > 2){
                    double _dz = zp[ID_i] - zp[ID_j];
                    _dr2 += std::pow(_dz, 2.0);
                }

                // Evaluate distance
                if (Pars::opt_ngh_interact == 1){
                    // Check the particle i
                    double _rSup = sp[ID_i] * Pars::r_sup * Pars::r_buff;   // Support radius of particle i
                    if (_dr2 < (_rSup*_rSup)){
                        nghIDList[ID_i].push_back(ID_j);
                    }

                    // Check the particle j
                    _rSup = sp[ID_j] * Pars::r_sup * Pars::r_buff;   // Support radius of particle j
                    if (_dr2 < (_rSup*_rSup)){
                        nghIDList[ID_j].push_back(ID_i);
                    }
                }
                else if (Pars::opt_ngh_interact == 2){
                    // Support radius of average size
                    double _rSupAve = ((sp[ID_i] + sp[ID_j]) / 2.0e0) * Pars::r_sup * Pars::r_buff;
                    if (_dr2 < (_rSupAve*_rSupAve)){
                        nghIDList[ID_i].push_back(ID_j);
                        nghIDList[ID_j].push_back(ID_i);
                    }
                }

                // Proceed to the next point in the grid
                ID_j = nextID[ID_j];
            }
        }

        // Include self or not
        if (Pars::flag_ngh_include_self){
            nghIDList[ID_i].push_back(ID_i);
        }
    }

    return;
}

// *End of the code