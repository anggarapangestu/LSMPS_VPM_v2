#include "fmm3D.hpp"

#define SUPPORT_RADIUS_FACTOR 2.0
/** NOTE: The support radius
 *   ______________    DESC:
 *  | . '|    |' . |    > The center cell X: have support region of circle with radius
 *  |'___|____|___'|       of SUPPORT_RADIUS_FACTOR multiplied by its length
 *  |    |[X:]|    |    > The particle inside the support region will use direct
 *  |____|____|____|       calculation to each source particle inside the cell
 *  |.   |    |   .|    > Otherwise, the particle outside it will use multipole
 *  |_'_.|____|._'_|
 * 
 *  ADD: 
 *   > From fast observation a lower SUPPORT_RADIUS_FACTOR will have faster calculation
 *   > A value of 1 show a quite fast but yield worse <?>
 *   > A value of 2 have longer computation but accurate
 *   > A value of 1.5 is just right (by fast observation)
*/

/** The type of FMM velocity field calculation
 *   [1] : Evaluate farfield using list 2,3, and 4 from tree data (IDK but problem for velocity inside the body region <?>)
 *   [2] : Evaluate farfield using radius different from tree data
*/
#define VELOCITY_CALCULATION_TYPE 2

// IMPORTANT NOTE:
//  > This code still not including the FMM acceleration at local calculation


// =====================================================
// +----------------- Utilites Method -----------------+
// =====================================================
// #pragma region UTILITIES_METHOD

/**
 *  @brief  Function to calculate the binomial combinatoric value of nCk.
 *         
 *  @param  n   Numerator of binomial combinatoric.
 *  @param  k   Denumerator of binomial combinatoric.
 * 
 *  @return Binomial combinatoric calculation.
*/
int fmm3D::binComb(int n, int k) const {
    // Exclusion for invalid range
    if (k > n) return 0;

    // Start calculating the binomial combinatoric
    int C = 1;
    int inv = n - k;    // The different inverse
    int _max, _min;     // The value of k and (n-k)

    // Set the maximum and the minimum
    if (inv > k){
        _max = inv;
        _min = k;
    }else{
        _max = k;
        _min = inv;
    }

    // Then C = n! / (k! * (n-k)!) = n! / (_max! * _min!)
    // Calculate the n!/_max!
    for (int i = n; i > _max; i--){
        C *= i;
    }
    
    // Calculate the 1/_min!
    for (int i = _min; i > 1; i--){
        C /= i;
    }
    
    return C;
}

/**
 *  @brief  A short hand to calculate all order in multipole container from particle.
 *         
 *  @param  multipole [OUTPUT] The multipole container at current cell.
 *  @param  src Particle source value.
 *  @param  dx  x coordinate distance from particle to cell center.
 *  @param  dy  y coordinate distance from particle to cell center.
 *  @param  dz  z coordinate distance from particle to cell center.
*/
void fmm3D::P2M_calc(std::vector<double> &mp, 
                  const double &srcVal,
                  const double &dx, 
                  const double &dy, 
                  const double &dz) const
{
    // Calculate the particle to multipole
    mp[0] += srcVal;                   // 0th order    (0,0,0)
    mp[1] += srcVal * dx;              // 1st order x  (1,0,0)
    mp[2] += srcVal * dy;              // 1st order y  (0,1,0)
    mp[3] += srcVal * dz;              // 1st order z  (0,0,1)
    mp[4] += srcVal * dx * dx * 0.5;   // 2nd order x  (2,0,0)
    mp[5] += srcVal * dy * dy * 0.5;   // 2nd order y  (0,2,0)
    mp[6] += srcVal * dz * dz * 0.5;   // 2nd order z  (0,0,2)
    mp[7] += srcVal * dx * dy;         // 1st x 1st y  (1,1,0)
    mp[8] += srcVal * dy * dz;         // 1st y 1st z  (0,1,1)
    mp[9] += srcVal * dx * dz;         // 1st x 1st z  (1,0,1)

    return;
}

/**
 *  @brief  A short hand to calculate all order in multipole from child multipole.
 *         
 *  @param  _mParent [OUTPUT] Parent cell multipole container.
 *  @param  _mChild  Child cell multipole container.
 *  @param  dx  x coordinate distance from child to parent cell center.
 *  @param  dy  y coordinate distance from child to parent cell center.
 *  @param  dz  z coordinate distance from child to parent cell center.
*/
void fmm3D::M2M_calc(std::vector<double> &mp,
                  std::vector<double> &mpC,
                  const double &dx,
                  const double &dy,
                  const double &dz) const
{
    // Calculate the multipole to multipole
    mp[0] += mpC[0];                    // 0th order    (0,0,0)
    mp[1] += mpC[1] +  mpC[0] * dx;     // 1st order x  (1,0,0)
    mp[2] += mpC[2] +  mpC[0] * dy;     // 1st order y  (0,1,0)
    mp[3] += mpC[3] +  mpC[0] * dz;     // 1st order z  (0,0,1)
    mp[4] += mpC[4] + (mpC[0] * dx * dx * 0.5) + (mpC[1] * dx);   // 2nd order x  (2,0,0)
    mp[5] += mpC[5] + (mpC[0] * dy * dy * 0.5) + (mpC[2] * dy);   // 2nd order y  (0,2,0)
    mp[6] += mpC[6] + (mpC[0] * dz * dz * 0.5) + (mpC[3] * dz);   // 2nd order z  (0,0,2)
    mp[7] += mpC[7] + (mpC[0] * dx * dy)       + (mpC[1] * dy) + (mpC[2] * dx);    // 1st x 1st y  (1,1,0)
    mp[8] += mpC[8] + (mpC[0] * dy * dz)       + (mpC[2] * dz) + (mpC[3] * dy);    // 1st y 1st z  (0,1,1)
    mp[9] += mpC[9] + (mpC[0] * dx * dz)       + (mpC[1] * dz) + (mpC[3] * dx);    // 1st x 1st z  (1,0,1)
    return;
}


/**
 *  @brief  A short hand to calculate all differential multiplier
 *  for farfield FMM calculation.
 *         
 *  @param  diff_x [OUTPUT] The multiplier for x differential.
 *  @param  diff_y [OUTPUT] The multiplier for y differential.
 *  @param  diff_z [OUTPUT] The multiplier for z differential.
 *  @param  R2  Distance squared.
 *  @param  dx  x coordinate distance from cell center to particle.
 *  @param  dy  y coordinate distance from cell center to particle.
 *  @param  dz  z coordinate distance from cell center to particle.
*/          
void fmm3D::M2P_mul_calc(std::vector<double> &diff_x, 
                  std::vector<double> &diff_y, 
                  std::vector<double> &diff_z, 
                  const double &R2, 
                  const double &dx, 
                  const double &dy, 
                  const double &dz) const
{
    // Calculate the farfield differential multiplier

    // Internal variable
    double R = sqrt(R2);    // Distance
    double R3 = R * R2;     // Distance power 3
    double R5 = R3 * R2;    // Distance power 5
    double R7 = R5 * R2;    // Distance power 7

    // Multiplier in x differential
    diff_x[0] = dx/R3;
    diff_x[1] = -(3*dx*dx/R5) + (1/R3);
    diff_x[2] = -(3*dx*dy/R5);
    diff_x[3] = -(3*dx*dz/R5);
    diff_x[4] = (15*dx*dx*dx/R7) - (9*dx/R5);
    diff_x[5] = (15*dx*dy*dy/R7) - (3*dx/R5);
    diff_x[6] = (15*dx*dz*dz/R7) - (3*dx/R5);
    diff_x[7] = (15*dx*dx*dy/R7) - (3*dy/R5);
    diff_x[8] = (15*dx*dy*dz/R7);
    diff_x[9] = (15*dx*dx*dz/R7) - (3*dz/R5);

    // Multiplier in y differential
    diff_y[0] = dy/R3;
    diff_y[1] = -(3*dy*dx/R5);
    diff_y[2] = -(3*dy*dy/R5) + (1/R3);
    diff_y[3] = -(3*dy*dz/R5);
    diff_y[4] = (15*dy*dx*dx/R7) - (3*dy/R5);
    diff_y[5] = (15*dy*dy*dy/R7) - (9*dy/R5);
    diff_y[6] = (15*dy*dz*dz/R7) - (3*dy/R5);
    diff_y[7] = (15*dy*dx*dy/R7) - (3*dx/R5);
    diff_y[8] = (15*dy*dy*dz/R7) - (3*dz/R5);
    diff_y[9] = (15*dy*dx*dz/R7);

    // Multiplier in z differential
    diff_z[0] = dz/R3;
    diff_z[1] = -(3*dz*dx/R5);
    diff_z[2] = -(3*dz*dy/R5);
    diff_z[3] = -(3*dz*dz/R5) + (1/R3);
    diff_z[4] = (15*dz*dx*dx/R7) - (3*dz/R5);
    diff_z[5] = (15*dz*dy*dy/R7) - (3*dz/R5);
    diff_z[6] = (15*dz*dz*dz/R7) - (9*dz/R5);
    diff_z[7] = (15*dz*dx*dy/R7);
    diff_z[8] = (15*dz*dy*dz/R7) - (3*dy/R5);
    diff_z[9] = (15*dz*dx*dz/R7) - (3*dx/R5);

    return;
}

// #pragma endregion


// =====================================================
// +-------------- FMM Direct Calculation -------------+
// =====================================================
// #pragma region FMM_DIRECT_CALCULATION

/**
 *  @brief  FMM internal tool to calculate direct potential.
 *         
 *  @param  _currID   The current cell ID to be evaluated.
 *  @param  _cellTree The cell tree data structure for data manager tools in 
 *  calculating FMM.
 *  @param  _nghCell  List of neighbor cell.
 *  @param  _parPos   List of particle position.
 *  @param  _srcVal   The particle potential source list data.
*/
void fmm3D::potDirSum(int _currID, const TreeCell &cellData, 
                      const std::vector<int> &_nghCell, 
                      const std::vector<std::vector<double>> &pos, 
                      const std::vector<double> &src)
{
    // Calculate each particle potential toward each target particle inside the neighbor cell

    // Alias to the current cell
    const Cell *currCell = cellData.treeData.at(_currID);
    
    // Evaluate all active neighbor cell
    for (auto _nghID: _nghCell){
        // Alias to the current neighbor ID
        const Cell *nghCell = cellData.treeData.at(_nghID);
        
        // Iterate through all source particle (at each active neighbor cell ID [_nghID])
        for (size_t j = 0; j < nghCell->parIDList.size(); j++){
            // The ID of source particle
            int _srcID = nghCell->parIDList[j];

            // Iterate through all target particle (inside the current cell ID [_currID])
            for (size_t i = 0; i < currCell->parIDList.size(); i++){
                // Internal variable
                double dx, dy, dz, R2;

                // The ID of target particle
                int _tarID = currCell->parIDList[i];

                // Dont put into calculation if source = target
                if (_srcID == _tarID) continue;

                // Set the origin position (source particle position)
                dx = pos[_tarID][0] - pos[_srcID][0];
                dy = pos[_tarID][1] - pos[_srcID][1];
                dz = pos[_tarID][2] - pos[_srcID][2];

                // Calculate the potential value
                R2 = dx*dx + dy*dy + dz*dz;

                // Perform the potential sum
                this->parPotential[_tarID] += src[_srcID] / sqrt(R2);
            }
        }
    }

    return;
}


/**
 *  @brief  FMM internal tool to calculate direct field.
 *         
 *  @param  _currID   The current cell ID to be evaluated.
 *  @param  _cellTree The cell tree data structure for data manager tools in 
 *  calculating FMM.
 *  @param  _nghCell  List of neighbor cell.
 *  @param  _parPos   List of particle position.
 *  @param  _srcVal   The particle potential source list data.
*/
void fmm3D::fieldDirSum(int _currID, const TreeCell &cellData, 
                        const std::vector<int> &_nghCell, 
                        const std::vector<std::vector<double>> &pos, 
                        const std::vector<double> &src)
{
    // Calculate each particle field toward each target particle inside the neighbor cell

    // Alias to the current cell
    const Cell *currCell = cellData.treeData.at(_currID);
    
    // Evaluate all active neighbor cell
    for (auto _nghID: _nghCell){
        // Alias to the current neighbor ID
        const Cell *nghCell = cellData.treeData.at(_nghID);

        // Iterate through all source particle (at each neighbor cell ID [_nghID])
        for (size_t j = 0; j < nghCell->parIDList.size(); j++){
            // The ID of source particle
            int _srcID = nghCell->parIDList[j];

            // Iterate through all target particle (inside the current cell ID [_currID])
            for (size_t i = 0; i < currCell->parIDList.size(); i++){
                // The ID of target particle
                int _tarID = currCell->parIDList[i];
            
                // Internal variable
                double dx, dy, dz, R, R2, R3;
                
                // Dont put into calculation if source = target
                if (_srcID == _tarID) continue;

                // The distance from target pivot at source (target - source)
                dx = pos[_tarID][0] - pos[_srcID][0];
                dy = pos[_tarID][1] - pos[_srcID][1];
                dz = pos[_tarID][2] - pos[_srcID][2];
                R2 = dx*dx + dy*dy + dz*dz;
                R  = sqrt(R2);
                R3 = R*R2;
                
                // Calculate the field for each source
                this->parField_x[_tarID] -= src[_srcID] * dx / R3;     // The result data field in x direction
                this->parField_y[_tarID] -= src[_srcID] * dy / R3;     // The result data field in y direction
                this->parField_z[_tarID] -= src[_srcID] * dz / R3;     // The result data field in y direction
            }
        }
    }

    return;
}

/**
 *  @brief  FMM internal tool to calculate direct velocity field.
 *         
 *  @param  _currID   The current cell ID to be evaluated.
 *  @param  _cellTree The cell tree data structure for data manager tools in 
 *  calculating FMM.
 *  @param  _nghCell  List of neighbor cell.
 *  @param  _parPos   List of particle position.
 *  @param  _alphaX   The particle vortex strength in x direction.
 *  @param  _alphaY   The particle vortex strength in y direction.
 *  @param  _alphaZ   The particle vortex strength in z direction.
*/
void fmm3D::velocityDirSum(int _currID, const TreeCell &cellData, 
                           const std::vector<int> &_nghCell, 
                           const std::vector<std::vector<double>> &pos, 
                           const std::vector<double> &alphaX,
                           const std::vector<double> &alphaY,
                           const std::vector<double> &alphaZ)
{
    // Calculate each particle field toward each target particle inside the neighbor cell

    // Alias to the current cell
    const Cell *currCell = cellData.treeData.at(_currID);
    
    // Evaluate all active neighbor cell
    for (auto _nghID: _nghCell){
        // Alias to the current neighbor ID
        const Cell *nghCell = cellData.treeData.at(_nghID);

        // Iterate through all source particle (at each neighbor cell ID [_nghID])
        for (size_t j = 0; j < nghCell->parIDList.size(); j++){
            // The ID of source particle
            int _srcID = nghCell->parIDList[j];

            // Iterate through all target particle (inside the current cell ID [_currID])
            for (size_t i = 0; i < currCell->parIDList.size(); i++){
                // The ID of target particle
                int _tarID = currCell->parIDList[i];
            
                // Internal variable
                double dx, dy, dz, R, R2, R3;
                
                // Dont put into calculation if source = target
                if (_srcID == _tarID) continue;

                // The distance from target pivot at source (target - source)
                dx = pos[_tarID][0] - pos[_srcID][0];
                dy = pos[_tarID][1] - pos[_srcID][1];
                dz = pos[_tarID][2] - pos[_srcID][2];
                R2 = dx*dx + dy*dy + dz*dz;
                R  = sqrt(R2);
                R3 = R*R2;
                
                // Calculate the velocity field for each source
                this->parField_x[_tarID] += (alphaY[_srcID]*dz - alphaZ[_srcID]*dy) / R3;
                this->parField_y[_tarID] += (alphaZ[_srcID]*dx - alphaX[_srcID]*dz) / R3;
                this->parField_z[_tarID] += (alphaX[_srcID]*dy - alphaY[_srcID]*dx) / R3;
            }
        }
    }

    return;
}

// #pragma endregion


// =====================================================
// +----------- FMM Fundamental Calculation -----------+
// =====================================================
// #pragma region FMM_FUNDAMENTAL_CALCULATION

/**
 *  @brief  The fundamental sequence of FMM translation calculation. This function
 *  includes upward pass.
 *         
 *  @param  _cellTree The cell tree data structure for data manager tools in 
 *  calculating FMM.
 *  @param  _parPos   List of particle position.
 *  @param  _activeFlag List of particle active mark.
 *  @param  _srcVal   The particle potential source list data.
*/
void fmm3D::setupFMM(const TreeCell &cellData, 
                     const std::vector<std::vector<double>> &parPos, 
                     const std::vector<bool> &activeMark, 
                     const std::vector<double> &srcVal)
{
    // Update the data for each class member variable
    this->parNum = parPos.size();           // Number of particle
    this->expOrd = 10;                      // Maximum expansion order (Set by current definition)
    this->cellNum  = cellData.cellCount;    // Total number of cell 
    this->maxLevel = cellData.max_level;    // Highest tree level

    /* To Do List:
        > Initialize each FMM variable a_k and b_k
        > [PROCEDURE 1] Calc. multipole source (mp) of each leaf cell using (Multipole Expansion)
        > [PROCEDURE 2] Calc. multipole source (mp) of all cell using M2M Translation (UP-PASS)
        
        No local translation
        |/| > [PROCEDURE 3] Calc. local source (bk) of all cell (contain the well separated neighbor cell only) using M2L Translation
        |/| > [PROCEDURE 4] Calc. local source (bk) of all cell using L2L Translation (DOWN-PASS)
    */

    // Resize the FMM source sum vector (note the index of expansion order from 0 -> max order)
    this->mp.clear();
    for (const auto &[_ID, _cell] : cellData.treeData){
        this->mp[_ID] = std::vector<double>(this->expOrd, 0.0);
    }
    // this->mp.resize(cellNum, std::vector<double>(this->expOrd,0.0));

    // Time counter
    double start, finish;
        
    // PROCEDURE 1! : Calc. multipole source (ak) of each leaf cell using (Multipole Expansion)
    // ************
    // Iterate through each leaf cell
    start = omp_get_wtime();
    
    #pragma omp parallel for
    for (size_t i = 0; i < cellData.leafList.size(); i++){
        // Internal variable
        double dx, dy, dz;      // Temporary distance variable

        // The ID of current evaluated leaf cell
        const int &_cellID = cellData.leafList[i];
        const Cell* currCell = cellData.treeData.at(_cellID);

        // Check if the current cell have source particle
        if (currCell->isActive == false){
            // There is no source particle here, thus skip this calculation
            continue;
        }
        
        // Calculate the multipole source for each expansion order
        for (size_t j = 0; j < currCell->parIDList.size(); j++){
            // The ID of particle inside cell
            int _parID = currCell->parIDList[j];
            
            // Only calculate the active particle
            if (activeMark[_parID] == false) continue;
            
            // Calculate the distance of cell center toward particle
            dx = currCell->centerPos[0] - parPos[_parID][0];
            dy = currCell->centerPos[1] - parPos[_parID][1];
            dz = currCell->centerPos[2] - parPos[_parID][2];

            // Calculate the multipole
            this->P2M_calc(mp.at(_cellID), srcVal[_parID], dx, dy, dz);
        }
    }

    finish = omp_get_wtime();
	printf("<+> FMM Procedure 1 [ME]   : %f s\n", finish-start);

        
    // PROCEDURE 2! : Calc. multipole source (ak) of all cell using M2M Translation (UP-PASS)
    // ************
    // Iterate through each cell from [level (k_max - 1)] to [level 1]
    start = omp_get_wtime();

    // Iterate from one level below the maximum level to level 1
    for (int level = this->maxLevel - 1; level > 0; level--){
        // Get the starting and finish ID at the current level
        const long long int begin_ID = cellData.get_startID(level);
        const long long int end_ID = cellData.get_startID(level + 1);

        // Evaluate all cell in the current level
        // #pragma omp parallel for
        for (int _cellID = begin_ID; _cellID < end_ID; _cellID++){
            // Only calculate the existed and non leaf and active
            // [CHECK 1] -> Skip un-existed cell
            if (cellData.treeData.count(_cellID) == 0) continue;

            // Alias to the cell
            const Cell* currCell = cellData.treeData.at(_cellID);

            // [CHECK 2] -> Skip the leaf cell
            if (currCell->isLeaf == true) continue;

            // [CHECK 3] -> Skip the non active cell
            if (currCell->isActive == false) continue;

            
            // ** Perform the M2M calculation

            // Internal variable
            double dx, dy, dz;          // Temporary distance variable for position

            // Find the child cell ID
            std::vector<int> childIDlist;
            childIDlist = cellData.get_childID(_cellID);

            // Calculate the local source sum (ak) using M2M from each child
            for (size_t i = 0; i < childIDlist.size(); i++){
                // Create the alias to the child
                int _chdID = childIDlist[i];
                // Check whether the child is existing or not
                if (cellData.treeData.count(_chdID) == 0) continue;
                // Create the alias to the child
                const Cell* chdCell = cellData.treeData.at(_chdID);

                // Calculate the distance of cell center toward particle
                dx = currCell->centerPos[0] - chdCell->centerPos[0];
                dy = currCell->centerPos[1] - chdCell->centerPos[1];
                dz = currCell->centerPos[2] - chdCell->centerPos[2];

                // Calculate the multipole
                this->M2M_calc(mp.at(_cellID), mp.at(_chdID), dx, dy, dz);
            }
        }
    }

    finish = omp_get_wtime();
	printf("<+> FMM Procedure 2 [M2M]  : %f s\n", finish-start);

    // In this 3D FMM there is no downward pass (VERY HARD TO EVALUATE)

    return;
}

/**
 *  @brief  The fundamental sequence of FMM translation calculation. This function
 *  includes upward pass. NOTE: This method is only specified for velocity calculation.
 *         
 *  @param  _cellTree The cell tree data structure for data manager tools in 
 *  calculating FMM.
 *  @param  _parPos   List of particle position.
 *  @param  _activeFlag List of particle active mark.
 *  @param  _alphaX   The particle vortex strength in x direction.
 *  @param  _alphaY   The particle vortex strength in y direction.
 *  @param  _alphaZ   The particle vortex strength in z direction.
*/
void fmm3D::setupVelocityFMM(const TreeCell &cellData, 
                             const std::vector<std::vector<double>> &parPos, 
                             const std::vector<bool> &activeMark, 
                             const std::vector<double> &alphaX,
                             const std::vector<double> &alphaY,
                             const std::vector<double> &alphaZ)
{
    // Update the data for each class member variable
    this->parNum = parPos.size();           // Number of particle
    this->expOrd = 10;                      // Maximum expansion order (Set by current definition)
    this->cellNum  = cellData.cellCount;    // Total number of cell 
    this->maxLevel = cellData.max_level;    // Highest tree level

    /* To Do List:
        > Initialize each FMM variable a_k and b_k
        > [PROCEDURE 1] Calc. multipole source (mp) of each leaf cell using (Multipole Expansion)
        > [PROCEDURE 2] Calc. multipole source (mp) of all cell using M2M Translation (UP-PASS)
        
        No local translation
        |/| > [PROCEDURE 3] Calc. local source (bk) of all cell (contain the well separated neighbor cell only) using M2L Translation
        |/| > [PROCEDURE 4] Calc. local source (bk) of all cell using L2L Translation (DOWN-PASS)
    */

    // Resize the FMM source sum vector (note the index of expansion order from 0 -> max order)
    int maxCell = cellData.get_startID(this->maxLevel + 1);
    // this->mp_Vx = std::vector<std::vector<double>>(maxCell, std::vector<double>(this->expOrd, 0.0));
    // this->mp_Vy = std::vector<std::vector<double>>(maxCell, std::vector<double>(this->expOrd, 0.0));
    // this->mp_Vz = std::vector<std::vector<double>>(maxCell, std::vector<double>(this->expOrd, 0.0));

    this->G2LcellID = std::vector<int>(maxCell, -1);
    this->mp_Vx = std::vector<std::vector<double>>(this->cellNum, std::vector<double>(this->expOrd, 0.0));
    this->mp_Vy = std::vector<std::vector<double>>(this->cellNum, std::vector<double>(this->expOrd, 0.0));
    this->mp_Vz = std::vector<std::vector<double>>(this->cellNum, std::vector<double>(this->expOrd, 0.0));
    int cnt = 0;
    for (const auto &[_ID, _cell] : cellData.treeData){
        G2LcellID[_ID] = cnt++;
    }
    
    // this->mp_Vx.clear();
    // this->mp_Vy.clear();
    // this->mp_Vz.clear();
    // for (const auto &[_ID, _cell] : cellData.treeData){
    //     this->mp_Vx[_ID] = std::vector<double>(this->expOrd, 0.0);
    //     this->mp_Vy[_ID] = std::vector<double>(this->expOrd, 0.0);
    //     this->mp_Vz[_ID] = std::vector<double>(this->expOrd, 0.0);
    // }

    // Time counter
    double start, finish;
        
    // PROCEDURE 1! : Calc. multipole source (ak) of each leaf cell using (Multipole Expansion)
    // ************
    // Iterate through each leaf cell
    start = omp_get_wtime();
    
    #pragma omp parallel for
    for (size_t i = 0; i < cellData.leafList.size(); i++){
        // Internal variable
        double dx, dy, dz;      // Temporary distance variable

        // The ID of current evaluated leaf cell
        const int &_cellID = cellData.leafList[i];
        const Cell* currCell = cellData.treeData.at(_cellID);

        // Check if the current cell have source particle
        if (currCell->isActive == false){
            // There is no source particle here, thus skip this calculation
            continue;
        }
        
        // Calculate the multipole source for each expansion order
        for (size_t j = 0; j < currCell->parIDList.size(); j++){
            // The ID of particle inside cell
            int _parID = currCell->parIDList[j];
            
            // Only calculate the active particle
            if (activeMark[_parID] == false) continue;
            
            // Calculate the distance of cell center toward particle
            dx = currCell->centerPos[0] - parPos[_parID][0];
            dy = currCell->centerPos[1] - parPos[_parID][1];
            dz = currCell->centerPos[2] - parPos[_parID][2];

            // Aliasing
            const int &cellID = this->G2LcellID[_cellID];

            // Calculate the multipole for x vorticity source
            this->P2M_calc(mp_Vx.at(cellID), alphaX[_parID], dx, dy, dz);

            // Calculate the multipole for y vorticity source
            this->P2M_calc(mp_Vy.at(cellID), alphaY[_parID], dx, dy, dz);

            // Calculate the multipole for z vorticity source
            this->P2M_calc(mp_Vz.at(cellID), alphaZ[_parID], dx, dy, dz);
        }
    }

    finish = omp_get_wtime();
	printf("<+> FMM Procedure 1 [ME]   : %f s\n", finish-start);

        
    // PROCEDURE 2! : Calc. multipole source (ak) of all cell using M2M Translation (UP-PASS)
    // ************
    // Iterate through each cell from [level (k_max - 1)] to [level 1]
    start = omp_get_wtime();

    // Iterate from one level below the maximum level to level 1
    for (int level = this->maxLevel - 1; level > 0; level--){
        // Get the starting and finish ID at the current level
        const long long int begin_ID = cellData.get_startID(level);
        const long long int end_ID = cellData.get_startID(level + 1);

        // Evaluate all cell in the current level
        // #pragma omp parallel for
        for (int _cellID = begin_ID; _cellID < end_ID; _cellID++){
            // Only calculate the existed and non leaf and active
            // [CHECK 1] -> Skip un-existed cell
            if (cellData.treeData.count(_cellID) == 0) continue;

            // Alias to the cell
            const Cell* currCell = cellData.treeData.at(_cellID);

            // [CHECK 2] -> Skip the leaf cell
            if (currCell->isLeaf == true) continue;

            // [CHECK 3] -> Skip the non active cell
            if (currCell->isActive == false) continue;

            
            // ** Perform the M2M calculation

            // Internal variable
            double dx, dy, dz;          // Temporary distance variable for position

            // Find the child cell ID
            std::vector<int> childIDlist;
            childIDlist = cellData.get_childID(_cellID);

            // Calculate the local source sum (ak) using M2M from each child
            for (size_t i = 0; i < childIDlist.size(); i++){
                // Create the alias to the child
                int _chdID = childIDlist[i];
                // Check whether the child is existing or not
                if (cellData.treeData.count(_chdID) == 0) continue;
                // Create the alias to the child
                const Cell* chdCell = cellData.treeData.at(_chdID);

                // Calculate the distance of cell center toward particle
                dx = currCell->centerPos[0] - chdCell->centerPos[0];
                dy = currCell->centerPos[1] - chdCell->centerPos[1];
                dz = currCell->centerPos[2] - chdCell->centerPos[2];

                // Index aliasing
                const int &cellID = this->G2LcellID[_cellID];
                const int &chdID = this->G2LcellID[_chdID];

                // Calculate the multipole translation for x vorticity
                this->M2M_calc(this->mp_Vx.at(cellID), this->mp_Vx.at(chdID), dx, dy, dz);

                // Calculate the multipole translation for y vorticity
                this->M2M_calc(this->mp_Vy.at(cellID), this->mp_Vy.at(chdID), dx, dy, dz);

                // Calculate the multipole translation for z vorticity
                this->M2M_calc(this->mp_Vz.at(cellID), this->mp_Vz.at(chdID), dx, dy, dz);
            }
        }
    }

    finish = omp_get_wtime();
	printf("<+> FMM Procedure 2 [M2M]  : %f s\n", finish-start);

    // In this 3D FMM there is no downward pass (VERY HARD TO EVALUATE)

    return;
}

// #pragma endregion


// =====================================================
// +----------------- Public Function -----------------+
// =====================================================
// #pragma region PUBLIC_FUNCTION

/**
 *  @brief  Calculate the potential using FMM method. [Not Available]
 *         
 *  @param  _cellTree The cell tree data structure for data manager tools in 
 *  calculating FMM.
 *  @param  _parPos   List of particle position.
 *  @param  _activeFlag List of particle active mark.
 *  @param  _srcVal   The particle potential source list data.
*/
void fmm3D::calcPotential(const TreeCell &cellData, 
                          const std::vector<std::vector<double>> &parPos, 
                          const std::vector<bool> &activeMark,
                          const std::vector<double> &srcVal)
{
    /* To Do List:
        > [PROCEDURE 5] Calculate the potential at each target position
        > [PROCEDURE 5.1] Calc. the potential caused by conjacent cell using DIRECT calculation
        > [PROCEDURE 5.2] Calc. the potential caused by not conjacent cell at nearest region using MULTIPOLE calculation
        > [PROCEDURE 5.3] Calc. the potential caused by far region cell using LOCAL calculation
    */
    
    // Setup Create
    this->setupFMM(cellData, parPos, activeMark, srcVal);

    // Internal variable
    // double x0, y0, x1, y1;          // Temporary variable for position
    // std::complex<double> z_0;       // Temporary variable for origin location
    // std::complex<double> z_1;       // Temporary variable for target location
    std::vector<int> List_1;        // Temporary Neigbor Cell ID list 1
    std::vector<int> List_2;        // Temporary Neigbor Cell ID list 2
    
    // Resize the variable size (note the index of expansion order from 0 -> max order)
    this->parPotential.clear(); this->parPotential.resize(this->parNum, 0.0);

    // Time counter
    double start, finish;
    start = omp_get_wtime();
    // double _timer, total1 = 0.0, total2 = 0.0, total3 = 0.0;

    // // PROCEDURE 5! : Still not constructed yet
    // // ************

    finish = omp_get_wtime();
	printf("<+> FMM Procedure 5 [FMM]  : %f s\n", finish-start);

    // printf("<+> FMM Procedure 5.1 : %f s [Direct SUM]\n", total1);
    // printf("<+> FMM Procedure 5.2 : %f s [Semi Farfield Calc.]\n", total2);
    // printf("<+> FMM Procedure 5.3 : %f s [Far-field Calc.]\n", total3);

    return;
}

/**
 *  @brief  Calculate the field using FMM method.
 *         
 *  @param  _cellTree The cell tree data structure for data manager tools in 
 *  calculating FMM.
 *  @param  _parPos   List of particle position.
 *  @param  _activeFlag List of particle active mark.
 *  @param  _srcVal   The particle potential source list data.
*/
void fmm3D::calcField(const TreeCell &cellData, 
                      const std::vector<std::vector<double>> &parPos, 
                      const std::vector<bool> &activeMark, 
                      const std::vector<double> &srcVal)
{
    /* To Do List:
        > [PROCEDURE 5] Calculate the field at each target position
        > [PROCEDURE 5.1] Calc. the field caused by conjacent cell using DIRECT calculation
        > [PROCEDURE 5.2] Evaluate all farfield cell
        > [PROCEDURE 5.3] Calc. the field caused by far region cell using multipole to particle calculation
    */
    
    // Setup Create
    this->setupFMM(cellData, parPos, activeMark, srcVal);

    // Time counter
    double start, finish;
    start = omp_get_wtime();

    // double _timer, total1 = 0.0, total2 = 0.0, total3 = 0.0;

    // Resize the variable size (note the index of expansion order from 0 -> max order)
    this->parField_x.clear(); this->parField_x.resize(parNum, 0.0);
    this->parField_y.clear(); this->parField_y.resize(parNum, 0.0);
    this->parField_z.clear(); this->parField_z.resize(parNum, 0.0);
    
    // PROCEDURE 5! : Calculate the field at each target position
    // ************
    // Calculate for each particle from the leaf cell
    #pragma omp parallel for /*reduction(+:total1,total2,total3)*/
    for (size_t i = 0; i < cellData.leafList.size(); i++){
        // The ID of current evaluated leaf cell
        const int _cellID = cellData.leafList[i];
        const Cell *currCell = cellData.treeData.at(_cellID);
        
        // Internal variable
        double dx, dy, dz;              // Temporary distance variable
        std::vector<int> List_1;        // Temporary Neigbor Cell ID list 1
        std::vector<int> List_2;        // Temporary Neigbor Cell ID list 2

        // Check the touched neighbor cell
        List_1.clear();
        List_2.clear();
        cellData.intList(_cellID, List_1, List_2);
        
        // PROCEDURE 5.1! : Calc. the field caused by conjacent cell using DIRECT calculation
        // ************
        // [PART 1] Calculate the direct field calculation form List_1 (calculate the field at each particle)
        // Reference: ORIGIN located at domain origin 
            // -> source and target particle position is refer to the domain origin
        
        // _timer = omp_get_wtime();
        this->fieldDirSum(_cellID, cellData, List_1, parPos, srcVal);
        // // MESSAGE_LOG << "Done Direct sum\n";
        // _timer = omp_get_wtime() - _timer;
        // total1 += _timer;


        // PROCEDURE 5.2! : Evaluate all farfield cell
        // ************
        // Take all farfield cell that contains
        //  > All particle in List 3
        //  > All particle in List 2 from the current cell to its parent recursively until level 2
        //  > All particle in List 4 from the current cell to its parent recursively until level 2
        
        // // Set up current timer
        // _timer = omp_get_wtime();

        // The farfield cell ID container
        // std::unordered_map<int,bool> farCellList;   // Use unordered map (*comment to change)
        std::vector<int> farCellList;            // Use vector (*uncomment to change)

        // Put the data cell ID from list 3 (or list 2 of the internal list)
        for (int &nghID : List_2){
            // farCellList.insert({nghID,true});   // Use unordered map
            farCellList.push_back(nghID);    // Use vector
        }

        // Put the data cell ID from list 2 and 4 from current cell through all parent cell
        int _parID, _currID = _cellID;        // Parent and current evaluated cell ID
        for (int level = currCell->level; level > 1; level--){
            // Check the touched neighbor cell
            List_1.clear();
            List_2.clear();
            cellData.extList(_currID, List_1, List_2);

            // Note: There are no any occurance for all 
            // List 1 and List 2 of all parent share the same cell ID
            
            // Put the type 1 neighbor
            for (int &nghID : List_1){
                // farCellList.insert({nghID,true});   // Use unordered map
                farCellList.push_back(nghID);    // Use vector
            }

            // Put the type 2 neighbor
            for (int &nghID : List_2){
                // farCellList.insert({nghID,true});   // Use unordered map`
                farCellList.push_back(nghID);    // Use vector
            }

            // Find the parent of current cell
            _parID = cellData.get_parentID(_currID);
            
            // Proceed to the next level
            _currID = _parID;
        }
        // _timer = omp_get_wtime() - _timer;
        // total2 += _timer;

        // cellData.saveSelTree(cellData,"far",farCellList);


        // PROCEDURE 5.3! : Calc. the field caused by far region cell using multipole to particle calculation
        // ************
        // Calculate the multipole to particle
        // Note : 
        //   > The formula have been derived previously
        //   > The derived formula looks neat and clean rather than the raw spherical harmonics
        
        // // Set up current timer
        // _timer = omp_get_wtime();

        // Iterate through each target particle from the current leaf cell
        for (size_t i = 1; i < currCell->parIDList.size(); i++){
            // The ID of particle inside cell
            int _parID = currCell->parIDList[i];

            // Translation of multipole to particle of FIELD calculation
            // for (const auto &[_srcCellID, _flag] : farCellList){
            for (size_t i = 0; i < farCellList.size(); i++){
                // // The ID of current evaluated leaf cell
                int _srcCellID = farCellList[i];
                const Cell *srcCell = cellData.treeData.at(_srcCellID);

                // Calculate the distance between target toward the source data
                dx = parPos[_parID][0] - srcCell->centerPos[0];
                dy = parPos[_parID][1] - srcCell->centerPos[1];
                dz = parPos[_parID][2] - srcCell->centerPos[2];
                double R2 = dx*dx + dy*dy + dz*dz;   // Distance square

                // Multipole differential multiplier
                std::vector<double> diff_mul_x(10);
                std::vector<double> diff_mul_y(10);
                std::vector<double> diff_mul_z(10);
                
                // Calculate the differential multiplier
                this->M2P_mul_calc(diff_mul_x, diff_mul_y, diff_mul_z, R2, dx, dy, dz);

                // Update the current calculated field
                for (int i = 0; i < this->expOrd; i++){
                    this->parField_x[_parID] -= this->mp.at(_srcCellID)[i] * diff_mul_x[i];
                    this->parField_y[_parID] -= this->mp.at(_srcCellID)[i] * diff_mul_y[i];
                    this->parField_z[_parID] -= this->mp.at(_srcCellID)[i] * diff_mul_z[i];
                }
            }
            
            // No local expansion
        }
        // _timer = omp_get_wtime() - _timer;
        // total3 += _timer;

    }
    
    finish = omp_get_wtime();
	printf("<+> FMM Procedure 5 [FMM]  : %f s\n", finish-start);

    // printf("<+> FMM Procedure 5.1 : %f s [Direct SUM]\n", total1);
    // printf("<+> FMM Procedure 5.2 : %f s [Farfield Cell Eval.]\n", total2);
    // printf("<+> FMM Procedure 5.3 : %f s [Far-field Calc.]\n", total3);

    return;
}


/**
 *  @brief  Calculate the velocity field using FMM method.
 *         
 *  @param  _cellTree The cell tree data structure for data manager tools in 
 *  calculating FMM.
 *  @param  _parPos   List of particle position.
 *  @param  _activeFlag List of particle active mark.
 *  @param  _alphaX   The particle vortex strength in x direction.
 *  @param  _alphaY   The particle vortex strength in y direction.
 *  @param  _alphaZ   The particle vortex strength in z direction.
*/
void fmm3D::calcVelocity(const TreeCell &cellData, 
                         const std::vector<std::vector<double>> &parPos, 
                         const std::vector<bool> &activeMark, 
                         const std::vector<double> &alphaX,
                         const std::vector<double> &alphaY,
                         const std::vector<double> &alphaZ)
{
    /* To Do List:
        > [PROCEDURE 3] Calculate the field at each target position
        > [PROCEDURE 3.1] Calc. the field caused by conjacent cell using DIRECT calculation
        > [PROCEDURE 3.2] Evaluate all farfield cell
        > [PROCEDURE 3.3] Calc. the field caused by far region cell using multipole to particle calculation
    */
    
    // Setup Create
    this->setupVelocityFMM(cellData, parPos, activeMark, alphaX, alphaY, alphaZ);

    // Time counter
    double start, finish;
    start = omp_get_wtime();

    // double _timer, total1 = 0.0, total2 = 0.0, total3 = 0.0;
    // double _timerInt, total3_1 = 0.0, total3_2 = 0.0;

    // Resize the variable size (note the index of expansion order from 0 -> max order)
    this->parField_x.clear(); this->parField_x.resize(parNum, 0.0);     // Velocity in x direction
    this->parField_y.clear(); this->parField_y.resize(parNum, 0.0);     // Velocity in y direction
    this->parField_z.clear(); this->parField_z.resize(parNum, 0.0);     // Velocity in z direction

    if (VELOCITY_CALCULATION_TYPE == 1){
        // PROCEDURE 3! : Calculate the field at each target position
        // ************
        // Calculate for each particle from the leaf cell
        #pragma omp parallel for /*reduction(+:total1,total2,total3)*/
        for (size_t i = 0; i < cellData.leafList.size(); i++){
            // The ID of current evaluated leaf cell
            const int _cellID = cellData.leafList[i];
            const Cell *currCell = cellData.treeData.at(_cellID);
            
            // Internal variable
            double dx, dy, dz;              // Temporary distance variable
            std::vector<int> List_1;        // Temporary Neigbor Cell ID list 1
            std::vector<int> List_2;        // Temporary Neigbor Cell ID list 2

            // Check the touched neighbor cell
            List_1.clear();
            List_2.clear();
            cellData.intList(_cellID, List_1, List_2);

            
            // PROCEDURE 3.1! : Calc. the field caused by conjacent cell using DIRECT calculation
            // ************
            // [PART 1] Calculate the direct field calculation form List_1 (calculate the field at each particle)
            // Reference: ORIGIN located at domain origin 
                // -> source and target particle position is refer to the domain origin
            
            // _timer = omp_get_wtime();
            this->velocityDirSum(_cellID, cellData, List_1, parPos, alphaX, alphaY, alphaZ);
            // // MESSAGE_LOG << "Done Direct sum\n";
            // _timer = omp_get_wtime() - _timer;
            // total1 += _timer;


            // PROCEDURE 3.2! : Evaluate all farfield cell
            // ************
            // Take all farfield cell that contains
            //  > All particle in List 3
            //  > All particle in List 2 from the current cell to its parent recursively until level 2
            //  > All particle in List 4 from the current cell to its parent recursively until level 2
            
            // // Set up current timer
            // _timer = omp_get_wtime();

            // The farfield cell ID container
            // std::unordered_map<int,bool> farCellList;   // Use unordered map (*comment to change)
            std::vector<int> farCellList;            // Use vector (*uncomment to change)

            // Put the data cell ID from list 3 (or list 2 of the internal list)
            for (int &nghID : List_2){
                // farCellList.insert({nghID,true});   // Use unordered map
                farCellList.push_back(nghID);    // Use vector
            }

            // Put the data cell ID from list 2 and 4 from current cell through all parent cell
            int _parID, _currID = _cellID;        // Parent and current evaluated cell ID
            for (int level = currCell->level; level > 1; level--){
                // Check the touched neighbor cell
                List_1.clear();
                List_2.clear();
                cellData.extList(_currID, List_1, List_2);

                // Note: There are no any occurance for all 
                // List 1 and List 2 of all parent share the same cell ID
                
                // Put the type 1 neighbor
                for (int &nghID : List_1){
                    // farCellList.insert({nghID,true});   // Use unordered map
                    farCellList.push_back(nghID);    // Use vector
                }

                // Put the type 2 neighbor
                for (int &nghID : List_2){
                    // farCellList.insert({nghID,true});   // Use unordered map`
                    farCellList.push_back(nghID);    // Use vector
                }

                // Find the parent of current cell
                _parID = cellData.get_parentID(_currID);
                
                // Proceed to the next level
                _currID = _parID;
            }
            // _timer = omp_get_wtime() - _timer;
            // total2 += _timer;


            // PROCEDURE 3.3! : Calc. the field caused by far region cell using multipole to particle calculation
            // ************
            // Calculate the multipole to particle
            // Note : 
            //   > The formula have been derived previously
            //   > The derived formula looks neat and clean rather than the raw spherical harmonics
            
            // // Set up current timer
            // _timer = omp_get_wtime();

            // Iterate through each target particle from the current leaf cell
            for (size_t i = 1; i < currCell->parIDList.size(); i++){
                // The ID of particle inside cell
                int _parID = currCell->parIDList[i];

                // Translation of multipole to particle of FIELD calculation
                // for (const auto &[_srcCellID, _flag] : farCellList){
                for (size_t i = 0; i < farCellList.size(); i++){
                    // // The ID of current evaluated leaf cell
                    int _srcCellID = farCellList[i];
                    const Cell *srcCell = cellData.treeData.at(_srcCellID);

                    // Calculate the distance between target toward the source data
                    dx = parPos[_parID][0] - srcCell->centerPos[0];
                    dy = parPos[_parID][1] - srcCell->centerPos[1];
                    dz = parPos[_parID][2] - srcCell->centerPos[2];
                    double R2 = dx*dx + dy*dy + dz*dz;   // Distance square

                    // Multipole differential multiplier
                    std::vector<double> diff_mul_x(10);
                    std::vector<double> diff_mul_y(10);
                    std::vector<double> diff_mul_z(10);

                    // Calculate the differential multiplier
                    this->M2P_mul_calc(diff_mul_x, diff_mul_y, diff_mul_z, R2, dx, dy, dz);
                    
                    // Index aliasing
                    const int &srcCellID = this->G2LcellID[_srcCellID];

                    // Update the current calculated field
                    for (int i = 0; i < this->expOrd; i++){
                        this->parField_x[_parID] += this->mp_Vy.at(srcCellID)[i]*diff_mul_z[i] - this->mp_Vz.at(srcCellID)[i]*diff_mul_y[i];
                        this->parField_y[_parID] += this->mp_Vz.at(srcCellID)[i]*diff_mul_x[i] - this->mp_Vx.at(srcCellID)[i]*diff_mul_z[i];
                        this->parField_z[_parID] += this->mp_Vx.at(srcCellID)[i]*diff_mul_y[i] - this->mp_Vy.at(srcCellID)[i]*diff_mul_x[i];
                    }
                }
                
                // No local expansion
            }
            // _timer = omp_get_wtime() - _timer;
            // total3 += _timer;
        }
    }
    
    else if (VELOCITY_CALCULATION_TYPE == 2){
        // PROCEDURE 3! : Calculate the field at each target position
        // ************
        // Calculate for each particle from the leaf cell
        #pragma omp parallel for /*reduction(+:total1,total2,total3)*/
        for (size_t i = 0; i < activeMark.size(); i++){
            // The ID of target particle
            int _tarID = i;

            // Internal variable
            double dx, dy, dz, R, R2, R3;
            const double &xp = parPos[_tarID][0];
            const double &yp = parPos[_tarID][1];
            const double &zp = parPos[_tarID][2];

            // Initialization of cell container queue list
            std::vector<int> queueList1 = {1,2,3,4,5,6,7,8};    // First queue container (Start at level 1)
            std::vector<int> queueList2;                        // Second queue container
            std::vector<int> *currQueue, *nextQueue;            // Cell container alias (pointer)
            
            // Iterate through all cell from level 1 to one level before maxLevel (because this section evaluate the child)
            for (int level = 1; level < this->maxLevel; level++){
                // Create the container queue alias 
                if (level % 2 == 1){
                    // For odd level 1, 3, 5, ...
                    currQueue = &queueList1;
                    nextQueue = &queueList2;
                }else if (level % 2 == 0){
                    // For even level 2, 4, 6, ...
                    currQueue = &queueList2;
                    nextQueue = &queueList1;
                }

                // Reserve the next queue
                nextQueue->clear();

                // Iterate through cell and evaluate each child
                for (const auto &_cellID : *currQueue){
                    // Initial check on the cell
                    // [CHECK 1] -> Skip un-existed cell
                    if (cellData.treeData.count(_cellID) == 0) continue;

                    // **Start to evaluate the child cell
                    // [PROCEDURE] Create 2 group for far cell and near cell
                    // > A far cell is directly calculated by M2P calculation
                    // > A near cell must be put into next queue for further check 
                    //    or direct biot savart for leaf cell

                    std::vector<int> chdIDList = cellData.get_childID(_cellID);
                    for (const auto &_chdID : chdIDList){
                        // [CHECK 1] -> Skip un-existed cell
                        // Initial check on the cell
                        if (cellData.treeData.count(_chdID) == 0) continue;
                        // Create the alias to the child
                        const Cell* currCell = cellData.treeData.at(_chdID);
                        // [CHECK 2] -> Skip the non active cell
                        if (currCell->isActive == false) continue;

                        // *Proceed only if the cell is existed and active
                        dx = xp - currCell->centerPos[0];
                        dy = yp - currCell->centerPos[1];
                        dz = zp - currCell->centerPos[2];
                        R2 = dx*dx + dy*dy + dz*dz;

                        // Check the distance (by check squared check)
                        double R_check = currCell->length * SUPPORT_RADIUS_FACTOR;
                        if (R2 > R_check*R_check){
                            // The cell considered as far cell
                            // _timer = omp_get_wtime();

                            // Multipole differential multiplier
                            std::vector<double> diff_mul_x(10);
                            std::vector<double> diff_mul_y(10);
                            std::vector<double> diff_mul_z(10);

                            // Calculate the differential multiplier
                            this->M2P_mul_calc(diff_mul_x, diff_mul_y, diff_mul_z, R2, dx, dy, dz);
                            
                            // Index aliasing
                            const int &cellID = this->G2LcellID[_chdID];

                            // Update the current calculated field
                            for (int i = 0; i < this->expOrd; i++){
                                this->parField_x[_tarID] += this->mp_Vy.at(cellID)[i]*diff_mul_z[i] - this->mp_Vz.at(cellID)[i]*diff_mul_y[i];
                                this->parField_y[_tarID] += this->mp_Vz.at(cellID)[i]*diff_mul_x[i] - this->mp_Vx.at(cellID)[i]*diff_mul_z[i];
                                this->parField_z[_tarID] += this->mp_Vx.at(cellID)[i]*diff_mul_y[i] - this->mp_Vy.at(cellID)[i]*diff_mul_x[i];
                            }
                            // _timer = omp_get_wtime() - _timer;
                            // total2 += _timer;
                        
                        }else{
                            // The cell considered as near cell
                            // Check whether a leaf cell or not
                            if (currCell->isLeaf == false){
                                // Put into the next queue
                                nextQueue->push_back(currCell->ID);
                            }else{
                                // _timer = omp_get_wtime();

                                // Calculate direct biot savart for velocity
                                // Iterate through all source particle (at current child cell)
                                for (size_t j = 0; j < currCell->parIDList.size(); j++){
                                    // The ID of source particle
                                    int _srcID = currCell->parIDList[j];
                                    
                                    // Dont put into calculation if source = target
                                    if (_srcID == _tarID) continue;

                                    // The distance from target pivot at source (target - source)
                                    dx = xp - parPos[_srcID][0];
                                    dy = yp - parPos[_srcID][1];
                                    dz = zp - parPos[_srcID][2];
                                    R2 = dx*dx + dy*dy + dz*dz;
                                    R  = sqrt(R2);
                                    R3 = R*R2;
                                    
                                    // Calculate the velocity field for each source
                                    this->parField_x[_tarID] += (alphaY[_srcID]*dz - alphaZ[_srcID]*dy) / R3;
                                    this->parField_y[_tarID] += (alphaZ[_srcID]*dx - alphaX[_srcID]*dz) / R3;
                                    this->parField_z[_tarID] += (alphaX[_srcID]*dy - alphaY[_srcID]*dx) / R3;
                                }

                                // _timer = omp_get_wtime() - _timer;
                                // total1 += _timer;
                            }
                        }
                    }
                }
                // End for this level queue
            }
            // End for this target particle
        }
    }

    finish = omp_get_wtime();
	printf("<+> FMM Procedure 3 [FMM]  : %f s\n", finish-start);

    // printf("<+> FMM Procedure 3.1 : %f s [Direct SUM]\n", total1);
    // printf("<+> FMM Procedure 3.2 : %f s [Farfield Cell Eval.]\n", total2);
    // printf("<+> FMM Procedure 3.3 : %f s [Far-field Calc.]\n", total3);

    // printf("<+> FMM Procedure 3.1 : %f s [Direct SUM]\n", total1);
    // printf("<+> FMM Procedure 3.2 : %f s [Far-field Calc.]\n", total2);

    return;
}


/**
 *  @brief  Calculate the velocity field using FMM method.
 *  NOTE: Only calculate the particle near body.
 *         
 *  @param  _cellTree The cell tree data structure for data manager tools in 
 *  calculating FMM.
 *  @param  _parPos   List of particle position.
 *  @param  _activeFlag List of particle active mark.
 *  @param  _targetFlag List of particle target mark.
 *  @param  _alphaX   The particle vortex strength in x direction.
 *  @param  _alphaY   The particle vortex strength in y direction.
 *  @param  _alphaZ   The particle vortex strength in z direction.
*/
void fmm3D::calcVelocityNearBody(const TreeCell &cellData, 
                                 const std::vector<std::vector<double>> &parPos, 
                                 const std::vector<bool> &activeMark, 
                                 const std::vector<bool> &targetFlag, 
                                 const std::vector<double> &alphaX,
                                 const std::vector<double> &alphaY,
                                 const std::vector<double> &alphaZ)
{
    /* To Do List:
        > [PROCEDURE 3] Calculate the field at each target position
        > [PROCEDURE 3.1] Calc. the field caused by conjacent cell using DIRECT calculation
        > [PROCEDURE 3.2] Evaluate all farfield cell
        > [PROCEDURE 3.3] Calc. the field caused by far region cell using multipole to particle calculation
    */
    
    // Setup Create
    this->setupVelocityFMM(cellData, parPos, activeMark, alphaX, alphaY, alphaZ);

    // Time counter
    double start, finish;
    start = omp_get_wtime();

    // double _timer, total1 = 0.0, total2 = 0.0, total3 = 0.0;
    // double _timerInt, total3_1 = 0.0, total3_2 = 0.0;

    // Resize the variable size (note the index of expansion order from 0 -> max order)
    this->parField_x.clear(); this->parField_x.resize(this->parNum, 0.0);     // Velocity in x direction
    this->parField_y.clear(); this->parField_y.resize(this->parNum, 0.0);     // Velocity in y direction
    this->parField_z.clear(); this->parField_z.resize(this->parNum, 0.0);     // Velocity in z direction

    // PROCEDURE 3! : Calculate the field at each target position
    // ************
    // Calculate for each particle from the leaf cell
    #pragma omp parallel for /*reduction(+:total1,total2,total3)*/
    for (size_t i = 0; i < activeMark.size(); i++){
        // Exception of calculation (only calculate the target particle)
        if (targetFlag[i] == false) continue;

        // The ID of target particle
        int _tarID = i;

        // Internal variable
        double dx, dy, dz, R, R2, R3;
        const double &xp = parPos[_tarID][0];
        const double &yp = parPos[_tarID][1];
        const double &zp = parPos[_tarID][2];

        // Initialization of cell container queue list
        std::vector<int> queueList1 = {1,2,3,4,5,6,7,8};    // First queue container (Start at level 1)
        std::vector<int> queueList2;                        // Second queue container
        std::vector<int> *currQueue, *nextQueue;            // Cell container alias (pointer)
        
        // Iterate through all cell from level 1 to one level before maxLevel (because this section evaluate the child)
        for (int level = 1; level < this->maxLevel; level++){
            // Create the container queue alias 
            if (level % 2 == 1){
                // For odd level 1, 3, 5, ...
                currQueue = &queueList1;
                nextQueue = &queueList2;
            }else if (level % 2 == 0){
                // For even level 2, 4, 6, ...
                currQueue = &queueList2;
                nextQueue = &queueList1;
            }

            // Reserve the next queue
            nextQueue->clear();

            // Iterate through cell and evaluate each child
            for (const auto &_cellID : *currQueue){
                // Initial check on the cell
                // [CHECK 1] -> Skip un-existed cell
                if (cellData.treeData.count(_cellID) == 0) continue;

                // **Start to evaluate the child cell
                // [PROCEDURE] Create 2 group for far cell and near cell
                // > A far cell is directly calculated by M2P calculation
                // > A near cell must be put into next queue for further check 
                //    or direct biot savart for leaf cell

                std::vector<int> chdIDList = cellData.get_childID(_cellID);
                for (const auto &_chdID : chdIDList){
                    // [CHECK 1] -> Skip un-existed cell
                    // Initial check on the cell
                    if (cellData.treeData.count(_chdID) == 0) continue;
                    // Create the alias to the child
                    const Cell* currCell = cellData.treeData.at(_chdID);
                    // [CHECK 2] -> Skip the non active cell
                    if (currCell->isActive == false) continue;

                    // *Proceed only if the cell is existed and active
                    dx = xp - currCell->centerPos[0];
                    dy = yp - currCell->centerPos[1];
                    dz = zp - currCell->centerPos[2];
                    R2 = dx*dx + dy*dy + dz*dz;

                    // Check the distance (by check squared check)
                    double R_check = currCell->length * SUPPORT_RADIUS_FACTOR;
                    if (R2 > R_check*R_check){
                        // The cell considered as far cell
                        // _timer = omp_get_wtime();

                        // Multipole differential multiplier
                        std::vector<double> diff_mul_x(10);
                        std::vector<double> diff_mul_y(10);
                        std::vector<double> diff_mul_z(10);

                        // Calculate the differential multiplier
                        this->M2P_mul_calc(diff_mul_x, diff_mul_y, diff_mul_z, R2, dx, dy, dz);
                        
                        // Index aliasing
                        const int &cellID = this->G2LcellID[_chdID];

                        // Update the current calculated field
                        for (int i = 0; i < this->expOrd; i++){
                            this->parField_x[_tarID] += this->mp_Vy.at(cellID)[i]*diff_mul_z[i] - this->mp_Vz.at(cellID)[i]*diff_mul_y[i];
                            this->parField_y[_tarID] += this->mp_Vz.at(cellID)[i]*diff_mul_x[i] - this->mp_Vx.at(cellID)[i]*diff_mul_z[i];
                            this->parField_z[_tarID] += this->mp_Vx.at(cellID)[i]*diff_mul_y[i] - this->mp_Vy.at(cellID)[i]*diff_mul_x[i];
                        }
                        // _timer = omp_get_wtime() - _timer;
                        // total2 += _timer;
                    
                    }else{
                        // The cell considered as near cell
                        // Check whether a leaf cell or not
                        if (currCell->isLeaf == false){
                            // Put into the next queue
                            nextQueue->push_back(currCell->ID);
                        }else{
                            // _timer = omp_get_wtime();

                            // Calculate direct biot savart for velocity
                            // Iterate through all source particle (at current child cell)
                            for (size_t j = 0; j < currCell->parIDList.size(); j++){
                                // The ID of source particle
                                int _srcID = currCell->parIDList[j];
                                
                                // Dont put into calculation if source = target
                                if (_srcID == _tarID) continue;

                                // The distance from target pivot at source (target - source)
                                dx = xp - parPos[_srcID][0];
                                dy = yp - parPos[_srcID][1];
                                dz = zp - parPos[_srcID][2];
                                R2 = dx*dx + dy*dy + dz*dz;
                                R  = sqrt(R2);
                                R3 = R*R2;
                                
                                // Calculate the velocity field for each source
                                this->parField_x[_tarID] += (alphaY[_srcID]*dz - alphaZ[_srcID]*dy) / R3;
                                this->parField_y[_tarID] += (alphaZ[_srcID]*dx - alphaX[_srcID]*dz) / R3;
                                this->parField_z[_tarID] += (alphaX[_srcID]*dy - alphaY[_srcID]*dx) / R3;
                            }

                            // _timer = omp_get_wtime() - _timer;
                            // total1 += _timer;
                        }
                    }
                }
            }
            // End for this level queue
        }
        // End for this target particle
    }

    finish = omp_get_wtime();
	printf("<+> FMM Procedure 3 [FMM]  : %f s\n", finish-start);

    // printf("<+> FMM Procedure 3.1 : %f s [Direct SUM]\n", total1);
    // printf("<+> FMM Procedure 3.2 : %f s [Farfield Cell Eval.]\n", total2);
    // printf("<+> FMM Procedure 3.3 : %f s [Far-field Calc.]\n", total3);

    // printf("<+> FMM Procedure 3.1 : %f s [Direct SUM]\n", total1);
    // printf("<+> FMM Procedure 3.2 : %f s [Far-field Calc.]\n", total2);

    return;
}

// #pragma endregion


// =====================================================
// +----------------- Getter Function -----------------+
// =====================================================
// #pragma region GETTER_FUNCTION

/**
 *  @brief  Get the FMM calculated potential data.
 * 
 *  @return The calculated potential.
*/
std::vector<double> fmm3D::get_Potential() const{
    return this->parPotential;
}

/**
 *  @brief  Get the FMM calculated field data.
 *         
 *  @return The calculated field data in x direction.
*/
std::vector<double> fmm3D::get_Field_x() const{
    return this->parField_x;
}

/**
 *  @brief  Get the FMM calculated field data.
 *         
 *  @return The calculated field data in y direction.
*/
std::vector<double> fmm3D::get_Field_y() const{
    return this->parField_y;
}

/**
 *  @brief  Get the FMM calculated field data.
 *         
 *  @return The calculated field data in z direction.
*/
std::vector<double> fmm3D::get_Field_z() const{
    return this->parField_z;
}

// #pragma endregion