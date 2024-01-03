#include "base_grid.hpp"

#define EXTREME_EXT_SIZE_FACTOR 0.5     // Expansion for grid generation

// Macro for vector
template <typename U> using vec = std::vector<U>;

// =====================================================
// +-------------- Grid Generation Method -------------+
// =====================================================
// #pragma region GRID_GENERATION

/**
 *  @brief  Initialize the grid data from the given point data.
 *         
 *  @param  _pntCoor   List of all point coordinate [x,y,z].
 *  @param  _parSize   Size of all point.
*/
void NghBaseGrid::initialize_grid(double maxCoor[DIM], double minCoor[DIM], const vec<double> &s){
    // Set the grid size
    double maxPntSize = *(std::max_element(s.begin(), s.end()));    // Size of maximum particle
    this->gridSize = maxPntSize * Pars::r_sup * Pars::r_buff;       // The size of suport radius + buffer radius

    /* Illustration for grid domain generation
        > [1] Initial domain 
           _____________    
          |       x    b|   * The [x] is the point inside domain 
          |  x          |       while [a] is minimum point and 
          |a____x_____x_|       [b] is maximum point
        
        > [2] Grid calculation (grid domain overlap by initial domain)
           _________________
          |_____|_____|_    |   * The grid is calculated from domain.
          |_____|_____|_|___|      Note that grid is bigger than
          |     |     | |   |      initial domain
          |_____|_____|_|___|
        
        > [3] Final domain
           _________________
          |  ---|-----|---  |  * The grid is shifted so that
          |_|___|_____|___|_|     initial domain is at the center
          | |   |     |   | |     toward the grid domain
          |__---|-----|---__|

    */

    // *[1] Get the domain size
    double maxPos[DIM];     // Global maximum point location
    double minPos[DIM];     // Global minimum point location
    basis_loop(d){
        // The extreme point coordinate add by expansion factor
        maxPos[d] = maxCoor[d] + EXTREME_EXT_SIZE_FACTOR * this->gridSize;
        minPos[d] = minCoor[d] - EXTREME_EXT_SIZE_FACTOR * this->gridSize;
    }
    
    // *[2] Calculate the index and grid count
    this->gridCount = 1;
    basis_loop(d){
        this->indexCount[d] = std::ceil((maxPos[d] - minPos[d])/this->gridSize);
        this->gridCount *= this->indexCount[d];
    }

    // *[3] Update the extreme location to make a symetry toward initial size
    basis_loop(d){
        // Calculate the coordinate shifting (different between new and old length) at this basis
        double coorShift = (this->indexCount[d]*this->gridSize) - (maxPos[d] - minPos[d]);
        
        // Shift the extreme point (*direct update beacause the old data is not needed any more)
        maxPos[d] += coorShift*0.5;
        minPos[d] -= coorShift*0.5;

        // Take the shifted minimum point as pivot coordinate
        this->pivCoor[d] = minPos[d];
    }

    // Initialize the index tranformation parameter
    this->flattenMul[0] = 1;
    for (int i = 1; i < DIM; i++){
        this->flattenMul[i] = this->flattenMul[i-1] * this->indexCount[i-1];
    }
    
    return;
}

// #pragma endregion

// =====================================================
// +-------------- Index Transformation ---------------+
// =====================================================
// #pragma region INDEX_TRANSFORMATION

/**
 *  @brief  Obtain the grid ID that contain a given point coordinate.
 *         
 *  @param  _pntCoor   Point coordinate to be evaluated [xi,yi,zi].
 * 
 *  @return The grid ID that contains input point coordinate.
*/
int NghBaseGrid::get_grid_ID(double _pos[DIM]) const{
    // Find the index of grid container
    int _index[DIM];    // Index of grid container
    basis_loop(d) _index[d] = (_pos[d] - this->pivCoor[d])/this->gridSize;

    // Flatten the index to get grid ID
    int gridID = 0;
    basis_loop(d) gridID += _index[d]*flattenMul[d];
    
    return gridID;
}

/**
 *  @brief  Get the grid position index from its ID.
 *         
 *  @param  _index   [OUTPUT] position index of grid [i,j,k].
 *  @param  ID  The grid ID as input parameter.
 * 
 *  @return grid index position from its given ID.
*/
void NghBaseGrid::get_grid_index(int _index[DIM], int ID) const{
    // Do the reverse of index flatten
    /* Calculation Illustration
        ID = idx*1 + idy*cx + idz*(cx*cy);      // Flatten algorithm
        *) idx =>       floor(ID/1) % cx = (idx + idy*cx + idz*(cx*cy)) % cx -> idx;
        *) idy =>      floor(ID/cx) % cy = (idy + idz*cy) % cy               -> idy;
        *) idz => floor(ID/(cx*cy)) % cy = idz % cx                          -> idz;
    */
    
    // Create allias to the multiplier and the moduler (for better reader)
    const int *div = this->flattenMul;
    const int *mod = this->indexCount;
    basis_loop(d)
    _index[d] = (ID/div[d]) % mod[d];
    
    return; // Done
}

// #pragma endregion

// =====================================================
// +--------------- Neighbor Utilities ----------------+
// =====================================================
// #pragma region NEIGHBOR_UTILITIES

/**
 *  @brief  Get the grid neighbor ID list.
 *         
 *  @param  _nghIDList   [OUTPUT] grid neighbor ID list container.
 *  @param  _gridID   grid ID which the neighbor about to be found.
*/
void NghBaseGrid::find_grid_ngh(vec<int> &_nghIDList, int _gridID) const{
    // Reserve the data
    _nghIDList.clear();

    // Find the neighbor grid
    int nghNum = Pars::intPow(3,DIM);   // The number of candidate neighboring grid
    int nghIndex[DIM];                  // Temporary neighbor index
    
    // Set the base index
    int baseIndex[DIM];
    this->get_grid_index(baseIndex, _gridID);
    basis_loop(d) baseIndex[d]--;

    // Iterate through all local neighbor ID
    for (int locID = 0; locID < nghNum; locID++){
        // Get the current neighbor index position (*also check whether in range)
        bool notNgh = false;
        basis_loop(d){
            nghIndex[d] = baseIndex[d] + ((locID/Pars::intPow(3,d))%3);
            // Check in range
            if (nghIndex[d] < 0 || nghIndex[d] >= this->indexCount[d]){
                // The current index is out of range
                notNgh = true;
                break;
            }
        }
        if (notNgh) continue;

        // Find the neighbor ID
        int nghID = 0;
        basis_loop(d) nghID += nghIndex[d]*flattenMul[d];

        // Push the neighbor into the neighbor ID list
        _nghIDList.push_back(nghID);
    }

    return;
}

// #pragma endregion

// =====================================================
// +------------------ Get Function -------------------+
// =====================================================

/**
 *  @brief  Get the size of single grid edge.
 *         
 *  @return The size of single grid edge.
*/
double NghBaseGrid::get_size(){
    return this->gridSize;
}

/**
 *  @brief  Get the grid total count member.
 *         
 *  @return  [INT] the total grid count.
*/
int NghBaseGrid::get_count(){
    return this->gridCount;
}