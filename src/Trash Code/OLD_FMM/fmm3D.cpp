#include "fmm3D.hpp"
#include <unordered_map>

// IMPORTANT NOTE:
//  > This code still not including the FMM acceleration at local calculation


/** FMM multipole to local second procedure method (lower level of leaf cell pole translation)
 *   1:= Calculation by using inverse of the fundamental multipole ([+] More accurate)
 *   2:= Calculation by using the M2L translation ([+] Faster)
*/
#define M2L_P2_METHOD_TYPE 2

/** FMM multipole potential/field procedure method (higher level of leaf cell pole translation)
 *   1:= Calculation by using the M2L translation ([+] Faster)
 *   2:= Calculation by using direct multipole ([+] More accurate)
*/
#define FMM_P2_METHOD_TYPE 1

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
int fmm3D::binComb(int n, int k){
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

// #pragma endregion

// =====================================================
// +----------------- FMM Calculation -----------------+
// =====================================================
// #pragma region FMM_CALCULATION

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
void fmm3D::potDirSum(int _currID, const treeCell &cellData, 
                      const std::vector<int> &_nghCell, 
                      const std::vector<std::vector<double>> &pos, 
                      const std::vector<double> &src)
{
    // Calculate each particle potential toward each target particle inside the neighbor cell

    // Collect only the active cell
    std::vector<int> _activeNghCell;
    for (int _nghID : _nghCell){
        if (cellData.srcNum[_nghID] != 0){
            _activeNghCell.push_back(_nghID);
            continue;
        }
    }
    
    // Evaluate all active neighbor cell
    for (auto _nghID: _activeNghCell){
        
        // Iterate through all source particle (at each active neighbor cell ID [_nghID])
        for (size_t j = 0; j < cellData.parIDList[_nghID].size(); j++){
            // The ID of source particle
            int _srcID = cellData.parIDList[_nghID][j];

            // Iterate through all target particle (inside the current cell ID [_currID])
            // #pragma omp parallel for
            for (size_t i = 0; i < cellData.parIDList[_currID].size(); i++){
                // Internal variable
                double dx, dy, dz, R2;

                // The ID of target particle
                int _tarID = cellData.parIDList[_currID][i];

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
void fmm3D::fieldDirSum(int _currID, const treeCell &cellData, 
                        const std::vector<int> &_nghCell, 
                        const std::vector<std::vector<double>> &pos, 
                        const std::vector<double> &src)
{
    // Calculate each particle field toward each target particle inside the neighbor cell

    // Collect only the active cell
    std::vector<int> _activeNghCell;
    for (int _nghID : _nghCell){
        if (cellData.srcNum[_nghID] != 0){
            _activeNghCell.push_back(_nghID);
            continue;
        }
    }
    
    // Evaluate all active neighbor cell
    for (auto _nghID: _activeNghCell){
        
        // Iterate through all source particle (at each neighbor cell ID [_nghID])
        for (size_t j = 0; j < cellData.parIDList[_nghID].size(); j++){
            // The ID of source particle
            int _srcID = cellData.parIDList[_nghID][j];

            // Iterate through all target particle (inside the current cell ID [_currID])
            // #pragma omp parallel for
            for (size_t i = 0; i < cellData.parIDList[_currID].size(); i++){
                // The ID of target particle
                int _tarID = cellData.parIDList[_currID][i];
            
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
 *  @brief  The fundamental sequence of FMM translation calculation. This function
 *  includes upward pass and downward pass.
 *         
 *  @param  _cellTree The cell tree data structure for data manager tools in 
 *  calculating FMM.
 *  @param  _parPos   List of particle position.
 *  @param  _activeFlag List of particle active mark.
 *  @param  _srcVal   The particle potential source list data.
*/
void fmm3D::setupFMM(const treeCell &cellData, 
                     const std::vector<std::vector<double>> &parPos, 
                     const std::vector<bool> &activeMark, 
                     const std::vector<double> &srcVal)
{
    // Update the data for each class member variable
    this->parNum = parPos.size();                 // Number of particle
    this->expOrd = 10;                            // Maximum expansion order
    this->cellNum  = cellData.cellPos.size();     // Total number of cell 
    this->maxLevel = cellData.level[cellNum-1];   // Highest tree level

    /* To Do List:
        > Initialize each FMM variable a_k and b_k
        > [PROCEDURE 1] Calc. multipole source (mp) of each leaf cell using (Multipole Expansion)
        > [PROCEDURE 2] Calc. multipole source (mp) of all cell using M2M Translation (UP-PASS)
        
        No local translation
        |/| > [PROCEDURE 3] Calc. local source (bk) of all cell (contain the well separated neighbor cell only) using M2L Translation
        |/| > [PROCEDURE 4] Calc. local source (bk) of all cell using L2L Translation (DOWN-PASS)
    */

    // Resize the FMM source sum vector (note the index of expansion order from 0 -> max order)
    this->mp.resize(cellNum, std::vector<double>(this->expOrd,0.0));

    // Time counter
    double start, finish;
        
    // PROCEDURE 1! : Calc. multipole source (ak) of each leaf cell using (Multipole Expansion)
    // ************
    // Iterate through each leaf cell
    start = omp_get_wtime();
    
    #pragma omp parallel for
    for (size_t i = 0; i < cellData.leafCellList.size(); i++){
        // Internal variable
        double dx, dy, dz;      // Temporary distance variable

        // The ID of current evaluated leaf cell
        const int &_cellID = cellData.leafCellList[i];

        // Check if the current cell have source particle
        if (cellData.srcNum[_cellID] == 0){
            // There is no source particle here, thus skip this calculation
            continue;
        }
        
        // Calculate the multipole source for each expansion order
        for (size_t j = 0; j < cellData.parIDList[_cellID].size(); j++){
            // The ID of particle inside cell
            int _parID = cellData.parIDList[_cellID][j];
            
            // Only calculate the active particle
            if (activeMark[_parID] != true) continue;
            
            // Calculate the distance of cell center toward particle
            dx = cellData.cellPos[_cellID][0] - parPos[_parID][0];
            dy = cellData.cellPos[_cellID][1] - parPos[_parID][1];
            dz = cellData.cellPos[_cellID][2] - parPos[_parID][2];

            // Calculate the multipole
            mp[_cellID][0] += srcVal[_parID];                   // 0th order    (0,0,0)
            mp[_cellID][1] += srcVal[_parID] * dx;              // 1st order x  (1,0,0)
            mp[_cellID][2] += srcVal[_parID] * dy;              // 1st order y  (0,1,0)
            mp[_cellID][3] += srcVal[_parID] * dz;              // 1st order z  (0,0,1)
            mp[_cellID][4] += srcVal[_parID] * dx * dx * 0.5;   // 2nd order x  (2,0,0)
            mp[_cellID][5] += srcVal[_parID] * dy * dy * 0.5;   // 2nd order y  (0,2,0)
            mp[_cellID][6] += srcVal[_parID] * dz * dz * 0.5;   // 2nd order z  (0,0,2)
            mp[_cellID][7] += srcVal[_parID] * dx * dy;         // 1st x 1st y  (1,1,0)
            mp[_cellID][8] += srcVal[_parID] * dy * dz;         // 1st y 1st z  (0,1,1)
            mp[_cellID][9] += srcVal[_parID] * dx * dz;         // 1st x 1st z  (1,0,1)
        }
        // NOTE: the value of ak is still (0 + 0i) at initial
        
        // Order of calculation / complexity:
        // > (Number of active particle) x (order of expansion) --> (N_act) x (Pmax+1)
    }

    finish = omp_get_wtime();
	printf("<+> FMM Procedure 1 [ME]   : %f s\n", finish-start);

        
    // PROCEDURE 2! : Calc. multipole source (ak) of all cell using M2M Translation (UP-PASS)
    // ************
    // Iterate through each cell from [level (k_max - 1)] to [level 1]
    start = omp_get_wtime();

    const int lastCellID = cellNum - Pars::intPow(cellData.basis, maxLevel) - 1;   // Find the last cell ID at level (k_max - 1)
    MESSAGE_LOG << "DATA OF ITERATION START FROM : " << lastCellID << "\n";
    for (int _cellID = lastCellID; _cellID > 0; _cellID--){
        // Internal variable
        double dx, dy, dz;          // Temporary distance variable for position

        // Only calculate if source particle is existed in the cell
        if (cellData.srcNum[_cellID] == 0){
            // There is no source particle, we should skip this cell
            continue;
        }

        // Find the child cell ID
        std::vector<int> childIDlist;
        childIDlist = cellData.findChildID(_cellID);

        // Calculate the local source sum (ak) using M2M from each child
        for (size_t i = 0; i < childIDlist.size(); i++){
            int _chdID = childIDlist[i];

            // Calculate the distance of cell center toward particle
            dx = cellData.cellPos[_cellID][0] - cellData.cellPos[_chdID][0];
            dy = cellData.cellPos[_cellID][1] - cellData.cellPos[_chdID][1];
            dz = cellData.cellPos[_cellID][2] - cellData.cellPos[_chdID][2];

            // Calculate the multipole
            mp[_cellID][0] += mp[_chdID][0];                        // 0th order    (0,0,0)
            mp[_cellID][1] += mp[_chdID][1] + mp[_chdID][0] * dx;   // 1st order x  (1,0,0)
            mp[_cellID][2] += mp[_chdID][2] + mp[_chdID][0] * dy;   // 1st order y  (0,1,0)
            mp[_cellID][3] += mp[_chdID][3] + mp[_chdID][0] * dz;   // 1st order z  (0,0,1)
            mp[_cellID][4] += mp[_chdID][4] + (mp[_chdID][0] * dx * dx * 0.5) + (mp[_chdID][1] * dx);   // 2nd order x  (2,0,0)
            mp[_cellID][5] += mp[_chdID][5] + (mp[_chdID][0] * dy * dy * 0.5) + (mp[_chdID][2] * dy);   // 2nd order y  (0,2,0)
            mp[_cellID][6] += mp[_chdID][6] + (mp[_chdID][0] * dz * dz * 0.5) + (mp[_chdID][3] * dz);   // 2nd order z  (0,0,2)
            mp[_cellID][7] += mp[_chdID][7] + (mp[_chdID][0] * dx * dy)       + (mp[_chdID][1] * dy) + (mp[_chdID][2] * dx);    // 1st x 1st y  (1,1,0)
            mp[_cellID][8] += mp[_chdID][8] + (mp[_chdID][0] * dy * dz)       + (mp[_chdID][2] * dz) + (mp[_chdID][3] * dy);    // 1st y 1st z  (0,1,1)
            mp[_cellID][9] += mp[_chdID][9] + (mp[_chdID][0] * dx * dz)       + (mp[_chdID][1] * dz) + (mp[_chdID][3] * dx);    // 1st x 1st z  (1,0,1)
        }
        // NOTE: the value of ak is still (0 + 0i) at initial

        // Order of calculation / complexity:
        // > (Number of left cell) x (4 child) x (order of expansion^2)/2 --> (N_celltot) x 2 x (Pmax+1)^2
    }

    finish = omp_get_wtime();
	printf("<+> FMM Procedure 2 [M2M]  : %f s\n", finish-start);

    // In this 3D FMM there is no downward pass (VERY HARD TO EVALUATE)


    // // PROCEDURE 3! : Calc. local source (bk) of all cell (contain the well separated neighbor cell only) using M2L Translation
    // // ************
    // // Calculate through all active cell (Containing particle) start from level 2
    // // Find the first cell ID at level 2
    // start = omp_get_wtime();
    // // double _timerML, stg1ML = 0.0, stg2ML = 0.0, stg3ML = 0.0, stg4ML = 0.0;

    // const int intCellID = 5;        // The first ID at 2nd level
    // #pragma omp parallel for
    // for (int _cellID = intCellID; _cellID < cellNum; _cellID++){
    //     // Internal variable
    //     double x0, y0, x1, y1;          // Temporary variable for position
    //     std::complex<double> z_0;       // Temporary variable for origin location
    //     std::complex<double> z_1;       // Temporary variable for target location
    //     std::vector<int> List_1;        // Temporary Neigbor Cell ID list 1
    //     std::vector<int> List_2;        // Temporary Neigbor Cell ID list 2

    //     // Only calculate if there is particle inside the cell
    //     if (cellData.parNum[_cellID] == 0){
    //         // There is no particle, we should skip this cell
    //         continue;
    //     }

    //     // _timerML = omp_get_wtime();
        
    //     // Update the new origin position (Current evaluated cell center)
    //     x0 = cellData.cellPos[_cellID][0];
    //     y0 = cellData.cellPos[_cellID][1];
    //     z_0 = std::complex<double>(x0, y0);

    //     // Calculate the cell source from well separated neighbor cell
    //     List_1.clear();
    //     List_2.clear();
    //     cellData.extList(_cellID, List_1, List_2);
        
    //     // The well separated cell source (same level)
    //     for (auto _nghID : List_1){
    //         // Update the new target position (The neighbor cell center)
    //         x1 = cellData.cellPos[_nghID][0];
    //         y1 = cellData.cellPos[_nghID][1];
    //         z_1 = std::complex<double>(x1, y1);
            
    //         // Calculate the source sum for order_p = 0
    //         bk[_cellID][0] += ak[_nghID][0] * std::log(z_0 - z_1);
    //         for (int k = 1; k <= expOrd; k++){
    //             bk[_cellID][0] -= ak[_nghID][k] / (double)k * std::pow(-1, k) * std::pow((z_1 - z_0), -k);
    //         }
            
    //         // Calculate the source sum for order_p > 0
    //         for (int p = 1; p <= expOrd; p++){
    //             // Addition of object outside the summation
    //             bk[_cellID][p] -= (ak[_nghID][0] / (double)p) * std::pow((z_1 - z_0), -p);
                
    //             int C;  // value for combinatoric
    //             for (int k = 1; k <= expOrd; k++){
    //                 C = this->binComb(p+k-1,k-1);
    //                 bk[_cellID][p] -= (double)C * std::pow(-1, k) / (double)k *
    //                                   ak[_nghID][k] * std::pow((z_1 - z_0), -(p + k));
    //             }
    //         }
    //     }
    //     // NOTE: The initial value of bk is still (0 + 0i) at first list
        
    //     // _timerML = omp_get_wtime() - _timerML;
    //     // stg1ML += _timerML;
    //     // _timerML = omp_get_wtime();

    //     // The well separated cell source (lower level)
        
    //     // Method of direct calculation (Supposedly low error but highly time consumpting)
    //     // INVERSE FUNDAMENTAL CALCULATION
    //     if (M2L_P2_METHOD_TYPE == 1){
    //     for (auto _nghID : List_2){
    //         // Calculate all particle inside the neighbor cell
    //         for (size_t i = 0; i < cellData.parIDList[_nghID].size(); i++){
    //             int _parID = cellData.parIDList[_nghID][i];    // The ID of particle inside cell
                
    //             // Only calculate the active particle
    //             if (activeMark[_parID] != true){
    //                 continue;
    //             }
                
    //             // Set the target position (Source particle position)
    //             x1 = parPos[_parID][0];
    //             y1 = parPos[_parID][1];
    //             z_1 = std::complex<double>(x1, y1);
                
    //             // Calculate the source sum for order_p = 0
    //             bk[_cellID][0] += srcVal[_parID] * std::log(z_0 - z_1);

    //             // Calculate the source sum for order_p > 0
    //             for (int p = 1; p <= expOrd; p++){
    //                 bk[_cellID][p] -= srcVal[_parID] * std::pow((z_1 - z_0), -p) / (double)p;
    //             }
    //         }
    //     }
    //     }

    //     // M2L TRANSLATION CALCULATION
    //     else if (M2L_P2_METHOD_TYPE == 2){
    //     // Similar method with "same level cell" (Supposedly higher error but efficient time consumption)
    //     for (auto _nghID : List_2){
    //         // Update the new target position (The neighbor cell center)
    //         x1 = cellData.cellPos[_nghID][0];
    //         y1 = cellData.cellPos[_nghID][1];
    //         z_1 = std::complex<double>(x1, y1);
            
    //         // Calculate the source sum for order_p = 0
    //         bk[_cellID][0] += ak[_nghID][0] * std::log(z_0 - z_1);
    //         for (int k = 1; k <= expOrd; k++){
    //             bk[_cellID][0] -= ak[_nghID][k] / (double)k * std::pow(-1, k) * std::pow((z_1 - z_0), -k);
    //         }
            
    //         // Calculate the source sum for order_p > 0
    //         for (int p = 1; p <= expOrd; p++){
    //             // Addition of object outside the summation
    //             bk[_cellID][p] -= (ak[_nghID][0] / (double)p) * std::pow((z_1 - z_0), -p);
                
    //             int C;  // value for combinatoric
    //             for (int k = 1; k <= expOrd; k++){
    //                 C = this->binComb(p+k-1,k-1);
    //                 bk[_cellID][p] -= (double)C * std::pow(-1, k) / (double)k *
    //                                   ak[_nghID][k] * std::pow((z_1 - z_0), -(p + k));
    //             }
    //         }
    //     }
    //     }

    //     // _timerML = omp_get_wtime() - _timerML;
    //     // stg2ML += _timerML;
    //     // NOTE: The initial value of bk is not (0 + 0i) at second list

    //     // Order of calculation / complexity:
    //     // > (Number of left cell) x (27 ngh) x (order of expansion^2) --> (N_celltot) x 27 x (Pmax+1)^2
    // }
    
    // finish = omp_get_wtime();
	// printf("<+> FMM Procedure 3 [M2L]  : %f s\n", finish-start);
    // // printf("<+> FMM Procedure 3 STG 1  : %f s\n", stg1ML);
    // // printf("<+> FMM Procedure 3 STG 2  : %f s\n", stg2ML);

    
    // // PROCEDURE 4! : Calc. local source (bk) of all cell using L2L Translation (DOWN-PASS)
    // // ************
    // // Evaluate all cell that have child
    // start = omp_get_wtime();
    
    // for (int _cellID = 0; _cellID <= lastCellID; _cellID++){
    //     // Internal variable
    //     double x0, y0, x1, y1;          // Temporary variable for position
    //     std::complex<double> z_0;       // Temporary variable for origin location
    //     std::complex<double> z_1;       // Temporary variable for target location

    //     // Only if the cell have a child
    //     if ((cellData.leafCellMark[_cellID] == true) || // No child
    //         (cellData.parNum[_cellID] == 0))            // The cell is not existed
    //     {
    //         continue;
    //     }
        
    //     // Set the new target position (Current cell center)
    //     x1 = cellData.cellPos[_cellID][0];
    //     y1 = cellData.cellPos[_cellID][1];
    //     z_1 = std::complex<double>(x1, y1);

    //     // Find the child cell ID
    //     std::vector<int> childIDlist;
    //     childIDlist = cellData.findChildID(_cellID);

    //     // Calculate the local source sum (ak) using M2M from each child
    //     for (size_t i = 0; i < childIDlist.size(); i++){
    //         int _chdID = childIDlist[i];

    //         // Set the new origin position (Child cell center)
    //         x0 = cellData.cellPos[_chdID][0];
    //         y0 = cellData.cellPos[_chdID][1];
    //         z_0 = std::complex<double>(x0, y0);
            
    //         // Calculate the source sum for order_p >= 0
    //         for (int p = 0; p <= expOrd; p++){
    //             int C;  // value for combinatoric
    //             for (int k = p; k <= expOrd; k++){
    //                 C = this->binComb(k,p);
    //                 bk[_chdID][p] += (double)C * bk[_cellID][k] * std::pow(z_0 - z_1, (k - p));
    //             }
    //         }
    //     }
    //     // NOTE: the value of bk is not (0 + 0i) at initial

    //     // Order of calculation / complexity:
    //     // > (Number of left cell) x (4 child) x (order of expansion^2)/2 --> (N_celltot) x 2 x (Pmax+1)^2

    // }
    
    // finish = omp_get_wtime();
	// printf("<+> FMM Procedure 4 [L2L]  : %f s\n", finish-start);
    return;
}

// #pragma endregion

// =====================================================
// +----------------- Public Function -----------------+
// =====================================================
// #pragma region PUBLIC_FUNCTION

/**
 *  @brief  Calculate the potential using FMM method.
 *         
 *  @param  _cellTree The cell tree data structure for data manager tools in 
 *  calculating FMM.
 *  @param  _parPos   List of particle position.
 *  @param  _activeFlag List of particle active mark.
 *  @param  _srcVal   The particle potential source list data.
*/
void fmm3D::calcPotential(const treeCell &cellData, 
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
    double x0, y0, x1, y1;          // Temporary variable for position
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

    // // PROCEDURE 5! : Calculate the potential at each target position
    // // ************
    // // Calculate for each particle from the leaf cell
    // for (size_t i = 0; i < cellData.leafCellList.size(); i++){
    //     int _cellID = cellData.leafCellList[i];     // The ID of current evaluated leaf cell

    //     // Check the touched neighbor cell
    //     List_1.clear();
    //     List_2.clear();
    //     // std::cout << "FINDING internal List of ID:"<<_cellID<<"\n";
    //     cellData.intList(_cellID, List_1, List_2);


    //     // PROCEDURE 5.1! : Calc. the potential caused by conjacent cell using DIRECT calculation
    //     // ************
    //     // [PART 1] Calculate the direct potential calculation form List_1 (calculate the potential at each particle)
    //     // Reference: ORIGIN located at domain origin 
    //         // -> source and target particle position is refer to the domain origin
    //         // phi (z) = sum{i=1,N_particle} (Q_i * log(z_target - z_source))
        
    //     // _timer = omp_get_wtime();
    //     this->potDirSum(_cellID, cellData, List_1, parPos, srcVal);
    //     // _timer = omp_get_wtime() - _timer;
    //     // total1 += _timer;

    //     // In order calculation of 
    //     // > (Number of total particle) x (9 ngh) x (N_src) --> (N_par_all) x 9 x (k_N)


    //     // PROCEDURE 5.2! : Calc. the potential caused by not conjacent cell at nearest region using MULTIPOLE calculation
    //     // ************
    //     // [PART 2] Calculate the first expansion potential calculation form List_2
    //     // There are 2 method provided to calculate the multipole method
    //     // int _method = 1;    // Change between 1 and 2
        
    //     // Translate the multipole source to the current cell using M2L algorithm
    //     // M2L CALCULATION
    //     if (FMM_P2_METHOD_TYPE == 1){
    //     // _timer = omp_get_wtime();
    //     // Update the new origin position (Current evaluated cell center)
    //     x0 = cellData.cellPos[_cellID][0];
    //     y0 = cellData.cellPos[_cellID][1];
    //     z_0 = std::complex<double>(x0, y0);
    //     // Supposedly higher error but efficient time consumption
    //     for (auto _nghID : List_2){
    //         // Update the new target position (The neighbor cell center)
    //         x1 = cellData.cellPos[_nghID][0];
    //         y1 = cellData.cellPos[_nghID][1];
    //         z_1 = std::complex<double>(x1, y1);
            
    //         // Calculate the source sum for order_p = 0
    //         bk[_cellID][0] += ak[_nghID][0] * std::log(z_0 - z_1);
    //         for (int k = 1; k <= expOrd; k++){
    //             bk[_cellID][0] -= ak[_nghID][k] / (double)k * std::pow(-1, k) / std::pow((z_1 - z_0),k);
    //         }
            
    //         // Calculate the source sum for order_p > 0
    //         for (int p = 1; p <= expOrd; p++){
    //             // Addition of object outside the summation
    //             bk[_cellID][p] -= (ak[_nghID][0] / (double)p) / std::pow((z_1 - z_0),p);
                
    //             int C;  // value for combinatoric
    //             for (int k = 1; k <= expOrd; k++){
    //                 C = this->binComb(p+k-1,k-1);
    //                 bk[_cellID][p] -= (double)C * std::pow(-1, k) / (double)k *
    //                                   ak[_nghID][k] / std::pow((z_1 - z_0),p+k);
    //             }
    //         }
    //     }
    //     // _timer = omp_get_wtime() - _timer;
    //     // total2 += _timer;
    //     }

    //     // Iterate through each particle target
    //     for (size_t i = 0; i < cellData.parIDList[_cellID].size(); i++){
    //         // std::cout << "DEBUG LINE 1\n";
    //         int _parID = cellData.parIDList[_cellID][i];    // The ID of particle inside cell
            
    //         // Temporary pontential calculation
    //         std::complex<double> _potential;

    //         // [PART 2] Calculate the first expansion potential calculation form List_2
    //         // Reference: ORIGIN located at neighbor cell center
    //             // -> source in local sum (ak) and target particle position refer to neighbor cell center
    //             // phi (z) = a0*log(z_target) + sum{k=1,P_max} ((1/k) * a_k / z_target^k)

    //         // Set the target position (particle position)
    //         x1 = parPos[_parID][0];
    //         y1 = parPos[_parID][1];
    //         z_1 = std::complex<double>(x1, y1);

    //         // Evaluate all neighbor cell
    //         // MULTIPOLE POTENTIAL CALCULATION
    //         if (FMM_P2_METHOD_TYPE == 2){
    //         // _timer = omp_get_wtime();
    //         // Supposedly higher time consumption but lower error
    //         for (auto _nghID : List_2){
    //             // Set the new origin position (Neighbor cell center)
    //             x0 = cellData.cellPos[_nghID][0];
    //             y0 = cellData.cellPos[_nghID][1];
    //             z_0 = std::complex<double>(x0, y0);

    //             // Calculate the potential cause by the order_p = 0
    //             _potential = ak[_nghID][0] * std::log(z_1 - z_0);

    //             // Calculate the potential cause by the order_p > 0
    //             for (int p = 1; p <= expOrd; p++){
    //                 _potential -= ak[_nghID][p] / std::pow((z_1 - z_0), p) / (double)p;
    //             }
                
    //             // Update the current calculated potential
    //             this->parPotential[_parID] += _potential.real();
    //         }
    //         // _timer = omp_get_wtime() - _timer;
    //         // total2 += _timer;
    //         }
    //         // In order calculation of 
    //         // > (Number of total particle) x (ngh fac) x (order of expansion) --> (N_par_all) x (ngh_fac) x (P_max + 1)

            
    //         // PROCEDURE 5.3! : Calc. the potential caused by far region cell using LOCAL calculation
    //         // ************
    //         // [PART 3] Calculate the potential caused by local source on the current cell
    //         // Reference: ORIGIN located at the current cell center
    //             // -> source in total sum (bk) and target particle position refer to current cell center
    //             // phi (z) = sum{k=0,P_max} (b_k * z_target^k)

    //         // _timer = omp_get_wtime();

    //         // Set the new origin position (Current cell center)
    //         x0 = cellData.cellPos[_cellID][0];
    //         y0 = cellData.cellPos[_cellID][1];
    //         z_0 = std::complex<double>(x0, y0);

    //         // Calculate the potential cause by the order_p = 0
    //         _potential = bk[_cellID][0];

    //         // Calculate the potential cause by the order_p > 0
    //         for (int p = 1; p <= expOrd; p++){
    //             _potential += bk[_cellID][p] * std::pow((z_1 - z_0), p);
    //         }
            
    //         // Update the current calculated potential
    //         this->parPotential[_parID] += _potential.real();

    //         // _timer = omp_get_wtime() - _timer;
    //         // total3 += _timer;

    //         // In order calculation of 
    //         // > (Number of total particle) x (order of expansion) --> (N_par_all) x (P_max + 1)
    //     }
        
    //     // NOTE: The initial value of bk is not (0 + 0i) at second list

    //     // Order of calculation / complexity:
    //     // > (Number of left cell) x (27 ngh) x (order of expansion^2) --> (N_celltot) x 27 x (Pmax+1)^2
    // }

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
void fmm3D::calcField(const treeCell &cellData, 
                      const std::vector<std::vector<double>> &parPos, 
                      const std::vector<bool> &activeMark, 
                      const std::vector<double> &srcVal)
{
    /* To Do List:
        > [PROCEDURE 5] Calculate the field at each target position
        > [PROCEDURE 5.1] Calc. the field caused by conjacent cell using DIRECT calculation
        > [PROCEDURE 5.2] Calc. the field caused by not conjacent cell at nearest region using MULTIPOLE calculation
        > [PROCEDURE 5.3] Calc. the field caused by far region cell using LOCAL calculation
    */
    
    // Setup Create
    this->setupFMM(cellData, parPos, activeMark, srcVal);

    MESSAGE_LOG << "Done FMM setup\n";

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
    // #pragma omp parallel for
    for (size_t i = 0; i < cellData.leafCellList.size(); i++){
        // The ID of current evaluated leaf cell
        int _cellID = cellData.leafCellList[i];

        MESSAGE_LOG << "Check into the loop " << _cellID << " \n";
        
        // Internal variable
        double x0, y0, z0, x1, y1, z1;  // Temporary variable for position
        double dx, dy, dz;              // Temporary distance variable
        std::vector<int> List_1;        // Temporary Neigbor Cell ID list 1
        std::vector<int> List_2;        // Temporary Neigbor Cell ID list 2

        std::unordered_map<int, bool> flagNearCell;
        for (int nearCellID : List_1){
            flagNearCell.insert({nearCellID, true});
        }

        // MESSAGE_LOG << "Done Check list\n";

        // Check the touched neighbor cell
        List_1.clear();
        List_2.clear();
        cellData.intList_3d(_cellID, List_1, List_2);

        // MESSAGE_LOG << "Done Check list\n";
        
        // PROCEDURE 5.1! : Calc. the field caused by conjacent cell using DIRECT calculation
        // ************
        // [PART 1] Calculate the direct field calculation form List_1 (calculate the field at each particle)
        // Reference: ORIGIN located at domain origin 
            // -> source and target particle position is refer to the domain origin
            // phi (z) = sum{i=1,N_particle} (Q_i /(z_target - z_source))
        
        // _timer = omp_get_wtime();
        this->fieldDirSum(_cellID, cellData, List_1, parPos, srcVal);
        // MESSAGE_LOG << "Done Direct sum\n";
        // _timer = omp_get_wtime() - _timer;
        // total1 += _timer;

        // Calculate the 

        // In order calculation of 
        // > (Number of total particle) x (9 ngh) x (N_src) --> (N_par_all) x 9 x (k_N)  


        // // PROCEDURE 5.2! : Calc. the field caused by not conjacent cell at nearest region using MULTIPOLE calculation
        // // ************
        // // [PART 2] Calculate the first expansion potential calculation form List_2
        // // There are 2 method provided to calculate the multipole method

        // // Translate the multipole source to the current cell using M2L algorithm
        // // M2L CALCULATION
        // if (FMM_P2_METHOD_TYPE == 1){
        // // _timer = omp_get_wtime();
        // // Update the new origin position (Current evaluated cell center)
        // x0 = cellData.cellPos[_cellID][0];
        // y0 = cellData.cellPos[_cellID][1];
        // z_0 = std::complex<double>(x0, y0);
        // // Supposedly higher error but efficient time consumption
        // for (auto _nghID : List_2){
        //     // Update the new target position (The neighbor cell center)
        //     x1 = cellData.cellPos[_nghID][0];
        //     y1 = cellData.cellPos[_nghID][1];
        //     z_1 = std::complex<double>(x1, y1);
            
        //     // Calculate the source sum for order_p = 0
        //     bk[_cellID][0] += ak[_nghID][0] * std::log(z_0 - z_1);
        //     for (int k = 1; k <= expOrd; k++){
        //         bk[_cellID][0] -= ak[_nghID][k] / (double)k * std::pow(-1, k) * std::pow((z_1 - z_0), -k);
        //     }
            
        //     // Calculate the source sum for order_p > 0
        //     for (int p = 1; p <= expOrd; p++){
        //         // Addition of object outside the summation
        //         bk[_cellID][p] -= (ak[_nghID][0] / (double)p) * std::pow((z_1 - z_0), -p);
                
        //         int C;  // value for combinatoric
        //         for (int k = 1; k <= expOrd; k++){
        //             C = this->binComb(p+k-1,k-1);
        //             bk[_cellID][p] -= (double)C * std::pow(-1, k) / (double)k *
        //                               ak[_nghID][k] * std::pow((z_1 - z_0),-(p + k));
        //         }
        //     }
        // }
        // // _timer = omp_get_wtime() - _timer;
        // // total2 += _timer;
        // }

        // Evaluate toward all leaf cell for farfield calculation from lowest level
        std::unordered_map<int, bool> farCell;

        // Iterate through each particle target from the home level
        for (size_t i = 1; i < cellData.cellNum; i++){
            // The ID of particle inside cell
            int _parID = cellData.parIDList[_cellID][i];

            // [PART 2] Calculate the first expansion field calculation form List_2
            // Reference: ORIGIN located at neighbor cell center
                // -> source in local sum (ak) and target particle position refer to neighbor cell center
                // phi (z) = sum{k=0,P_max} (a_k / z_target^(k+1))

            // Set the target position (particle position)
            x1 = parPos[_parID][0];
            y1 = parPos[_parID][1];
            z1 = parPos[_parID][2];

            // MULTIPOLE FIELD CALCULATION
            for (size_t i = 0; i < cellData.leafCellList.size(); i++){
                // The ID of current evaluated leaf cell
                int _srcCellID = cellData.leafCellList[i];

                // Only proceed if the current cell is having source
                if (cellData.srcNum[_srcCellID] == 0){
                    // This current cell is having no source, proceed to the next leaf cell
                    continue;
                }

                // Only proceed if this cell is in far field
                if (flagNearCell.count(_srcCellID) != 0){
                    // The current cell is at the near cell list, proceed to the next leaf cell
                    continue;
                }

                // Set the new origin position (Neighbor cell center)
                x0 = cellData.cellPos[_srcCellID][0];
                y0 = cellData.cellPos[_srcCellID][1];
                z0 = cellData.cellPos[_srcCellID][2];
                dx = x1 - x0;
                dy = y1 - y0;
                dz = z1 - z0;
                double R2 = dx*dx + dy*dy + dz*dz;   // Distance square
                double R = sqrt(R2);                 // Distance
                double R3 = R * R2;                  // Distance power 3
                double R5 = R3 * R2;                 // Distance power 5
                double R7 = R5 * R2;                 // Distance power 7

                // Multipole differential multiplier
                std::vector<double> diff_mul_x(10);
                std::vector<double> diff_mul_y(10);
                std::vector<double> diff_mul_z(10);

                // Multiplier in x differential
                diff_mul_x[0] = dx/R3;
                diff_mul_x[1] = -(3*dx*dx/R5) + (1/R3);
                diff_mul_x[2] = -(3*dx*dy/R5);
                diff_mul_x[3] = -(3*dx*dz/R5);
                diff_mul_x[4] = (15*dx*dx*dx/R7) - (9*dx/R5);
                diff_mul_x[5] = (15*dx*dy*dy/R7) - (3*dx/R5);
                diff_mul_x[6] = (15*dx*dz*dz/R7) - (3*dx/R5);
                diff_mul_x[7] = (15*dx*dx*dy/R7) - (3*dy/R5);
                diff_mul_x[8] = (15*dx*dy*dz/R7);
                diff_mul_x[9] = (15*dx*dx*dz/R7) - (3*dz/R5);
                // Multiplier in y differential
                diff_mul_y[0] = dy/R3;
                diff_mul_y[1] = -(3*dy*dx/R5);
                diff_mul_y[2] = -(3*dy*dy/R5) + (1/R3);
                diff_mul_y[3] = -(3*dy*dz/R5);
                diff_mul_y[4] = (15*dy*dx*dx/R7) - (3*dy/R5);
                diff_mul_y[5] = (15*dy*dy*dy/R7) - (9*dy/R5);
                diff_mul_y[6] = (15*dy*dz*dz/R7) - (3*dy/R5);
                diff_mul_y[7] = (15*dy*dx*dy/R7) - (3*dx/R5);
                diff_mul_y[8] = (15*dy*dy*dz/R7) - (3*dz/R5);
                diff_mul_y[9] = (15*dy*dx*dz/R7);
                // Multiplier in z differential
                diff_mul_z[0] = dz/R3;
                diff_mul_z[1] = -(3*dz*dx/R5);
                diff_mul_z[2] = -(3*dz*dy/R5);
                diff_mul_z[3] = -(3*dz*dz/R5) + (1/R3);
                diff_mul_z[4] = (15*dz*dx*dx/R7) - (3*dz/R5);
                diff_mul_z[5] = (15*dz*dy*dy/R7) - (3*dz/R5);
                diff_mul_z[6] = (15*dz*dz*dz/R7) - (9*dz/R5);
                diff_mul_z[7] = (15*dz*dx*dy/R7);
                diff_mul_z[8] = (15*dz*dy*dz/R7) - (3*dy/R5);
                diff_mul_z[9] = (15*dz*dx*dz/R7) - (3*dx/R5);

                // Update the current calculated field
                for (int i = 0; i < this->expOrd; i++){
                    this->parField_x[_parID] -= this->mp[_srcCellID][i] * diff_mul_x[i];
                    this->parField_y[_parID] -= this->mp[_srcCellID][i] * diff_mul_y[i];
                    this->parField_z[_parID] -= this->mp[_srcCellID][i] * diff_mul_z[i];
                }
            }
            
            // No local expansion
            // // [PART 3] Calculate the multipole expansion field calculation of the current cell
            // // Reference: ORIGIN located at the current cell center
            //     // -> source in total sum (bk) and target particle position refer to current cell center
            //     // phi (z) = sum{k=1,P_max} (k * b_k * z_target^(k-1))

            // // _timer = omp_get_wtime();
            
            // // Set the new origin position (Current cell center)
            // x0 = cellData.cellPos[_cellID][0];
            // y0 = cellData.cellPos[_cellID][1];
            // z_0 = std::complex<double>(x0, y0);

            // // Set up the field
            // _field = 0;

            // // Calculate the field cause by the order_p > 0
            // for (int p = 1; p <= expOrd; p++){
            //     _field += (double)p * bk[_cellID][p] * std::pow((z_1 - z_0), (p - 1));
            // }
            
            // // Update the current calculated field
            // this->parField_x[_parID] += _field.real();
            // this->parField_y[_parID] -= _field.imag();

            // // _timer = omp_get_wtime() - _timer;
            // // total3 += _timer;

            // // In order calculation of 
            // // > (Number of total particle) x (order of expansion) --> (N_par_all) x (P_max + 1)
        }
        
        // NOTE: The initial value of bk is not (0 + 0i) at second list

        // Order of calculation / complexity:
        // > (Number of left cell) x (27 ngh) x (order of expansion^2) --> (N_celltot) x 27 x (Pmax+1)^2
    }
    
    finish = omp_get_wtime();
	printf("<+> FMM Procedure 5 [FMM]  : %f s\n", finish-start);

    // printf("<+> FMM Procedure 5.1 : %f s [Direct SUM]\n", total1);
    // printf("<+> FMM Procedure 5.2 : %f s [Semi Farfield Calc.]\n", total2);
    // printf("<+> FMM Procedure 5.3 : %f s [Far-field Calc.]\n", total3);

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
std::vector<double> fmm3D::get_Potential(){
    return this->parPotential;
}

/**
 *  @brief  Get the FMM calculated field data.
 *         
 *  @return The calculated field data in x direction.
*/
std::vector<double> fmm3D::get_Field_x(){
    return this->parField_x;
}

/**
 *  @brief  Get the FMM calculated field data.
 *         
 *  @return The calculated field data in y direction.
*/
std::vector<double> fmm3D::get_Field_y(){
    return this->parField_y;
}

/**
 *  @brief  Get the FMM calculated field data.
 *         
 *  @return The calculated field data in z direction.
*/
std::vector<double> fmm3D::get_Field_z(){
    return this->parField_z;
}

// #pragma endregion