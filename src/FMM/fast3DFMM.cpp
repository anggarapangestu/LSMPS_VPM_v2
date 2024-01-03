#include "fast3DFMM.hpp"
#include "fmm3D.hpp"

#define SUPPORT_RADIUS_SCALE 1

void fmm3D::calcVelocityFast(const std::vector<std::vector<double>> &_parPos, 
                             const std::vector<bool> &_activeFlag,
                             const std::vector<bool> &_targetFlag,
                             const std::vector<double> &_alphaX,
                             const std::vector<double> &_alphaY,
                             const std::vector<double> &_alphaZ)
{
    // Here is the main manager of fast multipole method, a recursive method

    /** Procedure
     *   [1] Generate the tree cell data (*initialize, *create)
     *   [2] Calculate multipole of P2M and upward sweep of M2M
     *   [3] Evaluate the FMM (*direct and farfield of M2P)
    */

    // PROCEDURE 1:
    // ***********
    // Generate the tree data
    FMMCell treeData;
    treeData.generateTree(_parPos, _activeFlag);

    // PROCEDURE 2:
    // ***********
    // Calculate the multipole
    fastFMM3Dutils FMM_Tool;
    FMM_Tool.velocityCalc(treeData, _parPos, _activeFlag, _targetFlag, _alphaX, _alphaY, _alphaZ);

    // PROCEDURE 3:
    // ***********
    // Retrieve the velocity data
    this->parField_x = FMM_Tool.velX;
    this->parField_y = FMM_Tool.velY;
    this->parField_z = FMM_Tool.velZ;

    return;
}


// MULTIPOLE CALCULATION SECTION
void fastFMM3Dutils::multipoleCalc(const FMMCell &treeData, 
            const std::vector<std::vector<double>> &parPos, 
            const std::vector<bool> &activeMark,
            const std::vector<double> &alphaX,
            const std::vector<double> &alphaY,
            const std::vector<double> &alphaZ)
{
    // Calculate the multipole of all cell
    // Resize the multipole data before calculation
    this->m_Vx = std::vector<std::vector<double>>(treeData.cellNum, std::vector<double>(this->expOrd, 0.0));
    this->m_Vy = std::vector<std::vector<double>>(treeData.cellNum, std::vector<double>(this->expOrd, 0.0));
    this->m_Vz = std::vector<std::vector<double>>(treeData.cellNum, std::vector<double>(this->expOrd, 0.0));

    
    // SECTION 1:
    // *********
    // Calculate the particle to multipole at each leaf cell

    // Calculate P2M at each leaf cell
    #pragma omp parallel for
    for (int cellID = 0; cellID < treeData.cellNum; cellID++){
        // Only calculate the leaf cell
        // if (treeData.nPar[cellID] >= treeData.critNum) continue;    // Leaf cell have particle less than critNum
        if (treeData.chdFlag[cellID] != 0) continue;                // Leaf cell have no child

        // Only calculate active cell
        if (treeData.isActive[cellID] == false) continue;

        // Iterate through all particle inside the current cell
        for (size_t i = 0; i < treeData.parID[cellID].size(); i++){
            // Alias to the current particle ID
            const int &parID = treeData.parID[cellID][i];

            // Only active particle can contribute to the multipole (skip non active particle)
            if (activeMark[parID] == false) continue;

            // Calculate internal variable for P2M calculation
            double dx = treeData.xc[cellID] - parPos[parID][0];
            double dy = treeData.yc[cellID] - parPos[parID][1];
            double dz = treeData.zc[cellID] - parPos[parID][2];

            // Calculate the multipole for x vorticity source
            this->m_Vx.at(cellID)[0] += alphaX[parID];                   // 0th order    (0,0,0)
            this->m_Vx.at(cellID)[1] += alphaX[parID] * dx;              // 1st order x  (1,0,0)
            this->m_Vx.at(cellID)[2] += alphaX[parID] * dy;              // 1st order y  (0,1,0)
            this->m_Vx.at(cellID)[3] += alphaX[parID] * dz;              // 1st order z  (0,0,1)
            this->m_Vx.at(cellID)[4] += alphaX[parID] * dx * dx * 0.5;   // 2nd order x  (2,0,0)
            this->m_Vx.at(cellID)[5] += alphaX[parID] * dy * dy * 0.5;   // 2nd order y  (0,2,0)
            this->m_Vx.at(cellID)[6] += alphaX[parID] * dz * dz * 0.5;   // 2nd order z  (0,0,2)
            this->m_Vx.at(cellID)[7] += alphaX[parID] * dx * dy;         // 1st x 1st y  (1,1,0)
            this->m_Vx.at(cellID)[8] += alphaX[parID] * dy * dz;         // 1st y 1st z  (0,1,1)
            this->m_Vx.at(cellID)[9] += alphaX[parID] * dx * dz;         // 1st x 1st z  (1,0,1)

            // Calculate the multipole for y vorticity source
            this->m_Vy.at(cellID)[0] += alphaY[parID];                   // 0th order    (0,0,0)
            this->m_Vy.at(cellID)[1] += alphaY[parID] * dx;              // 1st order x  (1,0,0)
            this->m_Vy.at(cellID)[2] += alphaY[parID] * dy;              // 1st order y  (0,1,0)
            this->m_Vy.at(cellID)[3] += alphaY[parID] * dz;              // 1st order z  (0,0,1)
            this->m_Vy.at(cellID)[4] += alphaY[parID] * dx * dx * 0.5;   // 2nd order x  (2,0,0)
            this->m_Vy.at(cellID)[5] += alphaY[parID] * dy * dy * 0.5;   // 2nd order y  (0,2,0)
            this->m_Vy.at(cellID)[6] += alphaY[parID] * dz * dz * 0.5;   // 2nd order z  (0,0,2)
            this->m_Vy.at(cellID)[7] += alphaY[parID] * dx * dy;         // 1st x 1st y  (1,1,0)
            this->m_Vy.at(cellID)[8] += alphaY[parID] * dy * dz;         // 1st y 1st z  (0,1,1)
            this->m_Vy.at(cellID)[9] += alphaY[parID] * dx * dz;         // 1st x 1st z  (1,0,1)

            // Calculate the multipole for z vorticity source
            this->m_Vz.at(cellID)[0] += alphaZ[parID];                   // 0th order    (0,0,0)
            this->m_Vz.at(cellID)[1] += alphaZ[parID] * dx;              // 1st order x  (1,0,0)
            this->m_Vz.at(cellID)[2] += alphaZ[parID] * dy;              // 1st order y  (0,1,0)
            this->m_Vz.at(cellID)[3] += alphaZ[parID] * dz;              // 1st order z  (0,0,1)
            this->m_Vz.at(cellID)[4] += alphaZ[parID] * dx * dx * 0.5;   // 2nd order x  (2,0,0)
            this->m_Vz.at(cellID)[5] += alphaZ[parID] * dy * dy * 0.5;   // 2nd order y  (0,2,0)
            this->m_Vz.at(cellID)[6] += alphaZ[parID] * dz * dz * 0.5;   // 2nd order z  (0,0,2)
            this->m_Vz.at(cellID)[7] += alphaZ[parID] * dx * dy;         // 1st x 1st y  (1,1,0)
            this->m_Vz.at(cellID)[8] += alphaZ[parID] * dy * dz;         // 1st y 1st z  (0,1,1)
            this->m_Vz.at(cellID)[9] += alphaZ[parID] * dx * dz;         // 1st x 1st z  (1,0,1)
        }
    }

    // SECTION 2:
    // *********
    // Calculate the multipole to multipole at each cell from child cell
    // Iterates through level
    for (int level = treeData.maxLevel - 1; level > 0; level --){
        // Iterate from one level above maximum level
        // Alias to the begining and end of cellID at the current level
        int beginID = treeData.startID[level];
        int endID = treeData.startID[level+1];

        // Iterate through all cell in the current level
        #pragma omp parallel for
        for (int cellID = beginID; cellID < endID; cellID++){
            // Only calculte M2M for cell with child
            
            // Only calculate active cell (Skip if current cell is not active)
            if (treeData.isActive[cellID] == false) continue;

            // Skip if current cell is a leaf cell
            // if (treeData.nPar[cellID] < treeData.critNum) continue;    // Leaf cell have particle less than critNum
            if (treeData.chdFlag[cellID] == 0) continue;               // Leaf cell have no child

            // Calculate the local source sum (ak) using M2M from each child
            for (int octant = 0; octant < treeData.chdNum; octant++){
                // Check if the child cell existed at current octant
                if (((treeData.chdFlag[cellID] >> octant) % 2) == 0) continue;
                
                // Create the alias to the child
                int chdID = treeData.chdID[cellID][octant];

                // Calculate the distance of cell center toward particle
                double dx = treeData.xc[cellID] - treeData.xc[chdID];
                double dy = treeData.yc[cellID] - treeData.yc[chdID];
                double dz = treeData.zc[cellID] - treeData.zc[chdID];

                // Calculate the multipole translation for x vorticity
                this->m_Vx.at(cellID)[0] += this->m_Vx.at(chdID)[0];                        // 0th order    (0,0,0)
                this->m_Vx.at(cellID)[1] += this->m_Vx.at(chdID)[1] +  this->m_Vx.at(chdID)[0] * dx;   // 1st order x  (1,0,0)
                this->m_Vx.at(cellID)[2] += this->m_Vx.at(chdID)[2] +  this->m_Vx.at(chdID)[0] * dy;   // 1st order y  (0,1,0)
                this->m_Vx.at(cellID)[3] += this->m_Vx.at(chdID)[3] +  this->m_Vx.at(chdID)[0] * dz;   // 1st order z  (0,0,1)
                this->m_Vx.at(cellID)[4] += this->m_Vx.at(chdID)[4] + (this->m_Vx.at(chdID)[0] * dx * dx * 0.5) + (this->m_Vx.at(chdID)[1] * dx);   // 2nd order x  (2,0,0)
                this->m_Vx.at(cellID)[5] += this->m_Vx.at(chdID)[5] + (this->m_Vx.at(chdID)[0] * dy * dy * 0.5) + (this->m_Vx.at(chdID)[2] * dy);   // 2nd order y  (0,2,0)
                this->m_Vx.at(cellID)[6] += this->m_Vx.at(chdID)[6] + (this->m_Vx.at(chdID)[0] * dz * dz * 0.5) + (this->m_Vx.at(chdID)[3] * dz);   // 2nd order z  (0,0,2)
                this->m_Vx.at(cellID)[7] += this->m_Vx.at(chdID)[7] + (this->m_Vx.at(chdID)[0] * dx * dy)       + (this->m_Vx.at(chdID)[1] * dy) + (this->m_Vx.at(chdID)[2] * dx);    // 1st x 1st y  (1,1,0)
                this->m_Vx.at(cellID)[8] += this->m_Vx.at(chdID)[8] + (this->m_Vx.at(chdID)[0] * dy * dz)       + (this->m_Vx.at(chdID)[2] * dz) + (this->m_Vx.at(chdID)[3] * dy);    // 1st y 1st z  (0,1,1)
                this->m_Vx.at(cellID)[9] += this->m_Vx.at(chdID)[9] + (this->m_Vx.at(chdID)[0] * dx * dz)       + (this->m_Vx.at(chdID)[1] * dz) + (this->m_Vx.at(chdID)[3] * dx);    // 1st x 1st z  (1,0,1)

                // Calculate the multipole translation for y vorticity
                this->m_Vy.at(cellID)[0] += this->m_Vy.at(chdID)[0];                        // 0th order    (0,0,0)
                this->m_Vy.at(cellID)[1] += this->m_Vy.at(chdID)[1] +  this->m_Vy.at(chdID)[0] * dx;   // 1st order x  (1,0,0)
                this->m_Vy.at(cellID)[2] += this->m_Vy.at(chdID)[2] +  this->m_Vy.at(chdID)[0] * dy;   // 1st order y  (0,1,0)
                this->m_Vy.at(cellID)[3] += this->m_Vy.at(chdID)[3] +  this->m_Vy.at(chdID)[0] * dz;   // 1st order z  (0,0,1)
                this->m_Vy.at(cellID)[4] += this->m_Vy.at(chdID)[4] + (this->m_Vy.at(chdID)[0] * dx * dx * 0.5) + (this->m_Vy.at(chdID)[1] * dx);   // 2nd order x  (2,0,0)
                this->m_Vy.at(cellID)[5] += this->m_Vy.at(chdID)[5] + (this->m_Vy.at(chdID)[0] * dy * dy * 0.5) + (this->m_Vy.at(chdID)[2] * dy);   // 2nd order y  (0,2,0)
                this->m_Vy.at(cellID)[6] += this->m_Vy.at(chdID)[6] + (this->m_Vy.at(chdID)[0] * dz * dz * 0.5) + (this->m_Vy.at(chdID)[3] * dz);   // 2nd order z  (0,0,2)
                this->m_Vy.at(cellID)[7] += this->m_Vy.at(chdID)[7] + (this->m_Vy.at(chdID)[0] * dx * dy)       + (this->m_Vy.at(chdID)[1] * dy) + (this->m_Vy.at(chdID)[2] * dx);    // 1st x 1st y  (1,1,0)
                this->m_Vy.at(cellID)[8] += this->m_Vy.at(chdID)[8] + (this->m_Vy.at(chdID)[0] * dy * dz)       + (this->m_Vy.at(chdID)[2] * dz) + (this->m_Vy.at(chdID)[3] * dy);    // 1st y 1st z  (0,1,1)
                this->m_Vy.at(cellID)[9] += this->m_Vy.at(chdID)[9] + (this->m_Vy.at(chdID)[0] * dx * dz)       + (this->m_Vy.at(chdID)[1] * dz) + (this->m_Vy.at(chdID)[3] * dx);    // 1st x 1st z  (1,0,1)

                // Calculate the multipole translation for z vorticity
                this->m_Vz.at(cellID)[0] += this->m_Vz.at(chdID)[0];                        // 0th order    (0,0,0)
                this->m_Vz.at(cellID)[1] += this->m_Vz.at(chdID)[1] +  this->m_Vz.at(chdID)[0] * dx;   // 1st order x  (1,0,0)
                this->m_Vz.at(cellID)[2] += this->m_Vz.at(chdID)[2] +  this->m_Vz.at(chdID)[0] * dy;   // 1st order y  (0,1,0)
                this->m_Vz.at(cellID)[3] += this->m_Vz.at(chdID)[3] +  this->m_Vz.at(chdID)[0] * dz;   // 1st order z  (0,0,1)
                this->m_Vz.at(cellID)[4] += this->m_Vz.at(chdID)[4] + (this->m_Vz.at(chdID)[0] * dx * dx * 0.5) + (this->m_Vz.at(chdID)[1] * dx);   // 2nd order x  (2,0,0)
                this->m_Vz.at(cellID)[5] += this->m_Vz.at(chdID)[5] + (this->m_Vz.at(chdID)[0] * dy * dy * 0.5) + (this->m_Vz.at(chdID)[2] * dy);   // 2nd order y  (0,2,0)
                this->m_Vz.at(cellID)[6] += this->m_Vz.at(chdID)[6] + (this->m_Vz.at(chdID)[0] * dz * dz * 0.5) + (this->m_Vz.at(chdID)[3] * dz);   // 2nd order z  (0,0,2)
                this->m_Vz.at(cellID)[7] += this->m_Vz.at(chdID)[7] + (this->m_Vz.at(chdID)[0] * dx * dy)       + (this->m_Vz.at(chdID)[1] * dy) + (this->m_Vz.at(chdID)[2] * dx);    // 1st x 1st y  (1,1,0)
                this->m_Vz.at(cellID)[8] += this->m_Vz.at(chdID)[8] + (this->m_Vz.at(chdID)[0] * dy * dz)       + (this->m_Vz.at(chdID)[2] * dz) + (this->m_Vz.at(chdID)[3] * dy);    // 1st y 1st z  (0,1,1)
                this->m_Vz.at(cellID)[9] += this->m_Vz.at(chdID)[9] + (this->m_Vz.at(chdID)[0] * dx * dz)       + (this->m_Vz.at(chdID)[1] * dz) + (this->m_Vz.at(chdID)[3] * dx);    // 1st x 1st z  (1,0,1)
            }
        }
    }
    
    return;
}


void fastFMM3Dutils::evaluateFMM(const FMMCell &treeData, 
            const std::vector<std::vector<double>> &parPos, 
            const std::vector<bool> &activeMark,
            const std::vector<bool> &targetFlag,
            const std::vector<double> &alphaX,
            const std::vector<double> &alphaY,
            const std::vector<double> &alphaZ)
{
    // Alias to the number of particle
    int parNum = parPos.size();

    // Resize the variable size
    this->velX = std::vector<double>(parNum, 0.0);     // Velocity in x direction
    this->velY = std::vector<double>(parNum, 0.0);     // Velocity in y direction
    this->velZ = std::vector<double>(parNum, 0.0);     // Velocity in z direction


    // PROCEDURE 3! : Calculate the field at each target position
    // ************
    // Calculate for each particle from the leaf cell
    #pragma omp parallel for
    for (int i = 0; i < parNum; i++){
        // The ID of target particle
        int _tarID = i;

        // Internal variable
        double dx, dy, dz, R, R2, R3;
        const double &xp = parPos[_tarID][0];
        const double &yp = parPos[_tarID][1];
        const double &zp = parPos[_tarID][2];

        // Initialization of cell container queue list
        std::vector<int> queueList1;                // First queue container (Start at level 1)
        std::vector<int> queueList2;                // Second queue container
        std::vector<int> *currQueue, *nextQueue;    // Cell container alias (pointer)

        // Update the first queueList
        for (int oct = 0; oct < treeData.chdNum; oct++){
            if ((treeData.chdFlag[0] >> oct) % 2){
                queueList1.push_back(treeData.chdID[0][oct]);
            }
        }
        
        // Iterate through all cell from level 1 to one level before maxLevel (because this section evaluate the child)
        for (int level = 1; level < treeData.maxLevel; level++){
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
            for (const auto &cellID : *currQueue){
                // Initial check on the cell

                // **Start to evaluate the child cell
                // [PROCEDURE] Create 2 group for far cell and near cell
                // > A far cell is directly calculated by M2P calculation
                // > A near cell must be put into next queue for further check 
                //    or direct biot savart for leaf cell

                // Iterate all child
                for (int oct = 0; oct < treeData.chdNum; oct++){
                    // Check at the current octant, skip if cell is not existed
                    if (((treeData.chdFlag[cellID] >> oct) % 2) == 0) continue;

                    // Alias to the current child ID
                    int chdID = treeData.chdID[cellID][oct];

                    // Skip calculation if the current child cell is not active
                    if (treeData.isActive[chdID] == false) continue;

                    // *Proceed only if the cell is existed and active
                    dx = xp - treeData.xc[chdID];
                    dy = yp - treeData.yc[chdID];
                    dz = zp - treeData.zc[chdID];
                    R2 = dx*dx + dy*dy + dz*dz;

                    // Check the distance (by check squared check)
                    double R_check = treeData.size[chdID] * SUPPORT_RADIUS_SCALE;
                    if (R2 > R_check*R_check){
                        // The cell considered as far cell

                        // Calculate the farfield
                        double R = sqrt(R2);                 // Distance
                        double R3 = R * R2;                  // Distance power 3
                        double R5 = R3 * R2;                 // Distance power 5
                        double R7 = R5 * R2;                 // Distance power 7

                        // Multipole differential multiplier
                        std::vector<double> diff_mul_x(this->expOrd);
                        std::vector<double> diff_mul_y(this->expOrd);
                        std::vector<double> diff_mul_z(this->expOrd);

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
                            this->velX[_tarID] += this->m_Vy.at(chdID)[i]*diff_mul_z[i] - this->m_Vz.at(chdID)[i]*diff_mul_y[i];
                            this->velY[_tarID] += this->m_Vz.at(chdID)[i]*diff_mul_x[i] - this->m_Vx.at(chdID)[i]*diff_mul_z[i];
                            this->velZ[_tarID] += this->m_Vx.at(chdID)[i]*diff_mul_y[i] - this->m_Vy.at(chdID)[i]*diff_mul_x[i];
                        }
                    
                    }else{
                        // The cell considered as near cell
                        // Check whether a leaf cell or not
                        if (treeData.chdFlag[chdID] != 0){
                            // Put into the next queue
                            nextQueue->push_back(chdID);
                        }else{
                            // _timer = omp_get_wtime();

                            // Calculate direct biot savart for velocity
                            // Iterate through all source particle (at current child cell)
                            for (size_t j = 0; j < treeData.parID[chdID].size(); j++){
                                // The ID of source particle
                                int _srcID = treeData.parID[chdID][j];

                                // Skip non active particle
                                if (activeMark[_srcID] == false) continue;
                                
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
                                this->velX[_tarID] += (alphaY[_srcID]*dz - alphaZ[_srcID]*dy) / R3;
                                this->velY[_tarID] += (alphaZ[_srcID]*dx - alphaX[_srcID]*dz) / R3;
                                this->velZ[_tarID] += (alphaX[_srcID]*dy - alphaY[_srcID]*dx) / R3;
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
    return;
}


void fastFMM3Dutils::velocityCalc(const FMMCell &treeData, 
            const std::vector<std::vector<double>> &_parPos, 
            const std::vector<bool> &_activeFlag,
            const std::vector<bool> &_targetFlag,
            const std::vector<double> &_alphaX,
            const std::vector<double> &_alphaY,
            const std::vector<double> &_alphaZ)
{
    // Calculate multipole
    this->multipoleCalc(treeData, _parPos, _activeFlag, _alphaX, _alphaY, _alphaZ);

    // Calculate velocity
    this->evaluateFMM(treeData, _parPos, _activeFlag, _targetFlag, _alphaX, _alphaY, _alphaZ);

    return;
}

// TREE DATA GENERATION SECTION
void FMMCell::generateTree(const std::vector<std::vector<double>> &parPos, 
                           const std::vector<bool> &activeFlag)
{
    /** Procedure
     *   [1] Tree initialization
     *       1.1 Calculate the tree fundamental data
     *       1.2 Assign the root cell data
     *   [2] Generate the tree from particle data
    */
    
    // PROCEDURE 1:
    // ***********
    // Initialize the fundamental data

    // [1.1] Calculate the tree fundamental data
    // Internal variable
    double minDom[DIM];     // The location of global minimum
    double maxDom[DIM];     // The location of global maximum
    double domCen[DIM];     // The domain center coordinate location
    double domLen[DIM];     // The domain length

    // Initialize the element of domain extreme coordinate
    minDom[0] = Pars::xcenter; maxDom[0] = Pars::xcenter;
    minDom[1] = Pars::ycenter; maxDom[1] = Pars::ycenter;
    minDom[2] = Pars::zcenter; maxDom[2] = Pars::zcenter;

    for (size_t i = 0; i < parPos.size(); i++){
        basis_loop(d){
            // Update the minimum coordinate position
            if (parPos[i][d] < minDom[d]) minDom[d] = parPos[i][d];
            // Update the maximum coordinate position
            if (parPos[i][d] > maxDom[d]) maxDom[d] = parPos[i][d];
        }
    }

    // Update the domain length and domain center coordinate
    this->rootLength = 0.0;
    basis_loop(d){
        domLen[d] = maxDom[d] - minDom[d];
        this->rootLength = std::max(domLen[d], this->rootLength);
        domCen[d] = 0.5 * (maxDom[d] + minDom[d]);
    }

    // Expand the domain size
    this->rootLength *= (1 + (Pars::expTree / 100));
    
    
    // [1.2] Assign the root cell data
    // Update the root fundamental data
    this->cellNum++;                        // Add the root cell into account
	this->nPar.push_back(parPos.size());    // Number of all particle (the root contains all particle in the domain)
	this->chdFlag.push_back(0);             // Child octant flag (initially no child, set to zero)
	this->parentID.push_back(-1);           // Parent ID of current cell (to be exact, root cell have no parent)
	this->xc.push_back(domCen[0]);          // Cell center x coordinate
	this->yc.push_back(domCen[1]);          // Cell center y coordiante
	this->zc.push_back(domCen[2]);          // Cell center z coordiante
	this->size.push_back(this->rootLength); // Cell box length size
    this->isActive.push_back(true);         // The root cell mark must be active

    // Update the root utility data
    this->parID.push_back(std::vector<int>());   // List of the first 'critNum' particle ID in the cell
	this->chdID.push_back(std::vector<int>(this->chdNum,0));    // List of child cell ID at each octant

    // Update the starting cell ID at the root level
    this->startID.clear();
    this->startID.push_back(0);

    /** ILLUSTRATION OF OCTANT
     *                            k=1  _________ [j]  
     *                                | #6 | #7 | 1   
     *                                |____|____|     
     *                                | #4 | #5 | 0   
     *                                |____|____| 
     *                k=0  _________    0     1 [i]
     *                    | #2 | #3 | 1
     *                    |____|____|
     *                    | #0 | #1 | 0
     *                    |____|____|
     *                 [i]  0     1
     *  Description:
     *   > i,j,k - is the local index in x,y,z direction respectively
     *   > The number label with # sign denote the local ID or octant
     *   > For example index [0,1,1] located at octant [#6]
     * 
     *  We can obtain the octant of particle by comparing the 
     *   particle coordinate with the cell center coordinate
     * 
     *  [DETAIL:]
     *    - for particle position > cellCenPos : have index of 1 -> set as true
     *    - for particle position < cellCenPos : have index of 0 -> set as false
     * 
     *  [ADDITIONAL:]
     *  The child flag is defined in binary (correspond to the octant)
     *     Binary  10000101   -> Flag of current cell child
     *             |    | |
     *     Octant  7    2 0   -> The 'on' flag at octant 0, 2, and 7
     *  
     *  To turn on or off the child flag at any specific octant
     *   [Turn on]
     *       =>  cell_flag | octant_flag     {A binary 'or' operator}
     *   [Turn off]
     *       =>  cell_flag & (~octant_flag)  {A binary 'and' + 'not' operator}
     *   Where:
     *       =>  octant_flag = (1 << octant) {A binary 'left shift' operator}
     * 
     *  EXAMPLE:
     *    Cell with child at octant 0,3, and 5, 
     *     it want to add child at octant 2 
     *     it also want to turn off the child at octant 3
     * 
     *  TURN ON procedure
     *     cell_flag   = 00101001 (Octant 0,3,5)
     *     octant_flag = 1 << 2 => 00000001 << 2 => 00000100
     *  Then cell_flag & octant_flag = 00101001 | 00000100
     *                               = 00101101 -> Binary
     *                                   | || |
     *                                   5 32 0 -> Octant
     *  TURN OFF procedure
     *     cell_flag   = 00101101 (Octant 0,2,3,5)
     *     octant_flag = 1 << 3 => 00000001 << 3 => 00001000
     *  Then cell_flag & (~octant_flag) = 00101001 & (~00001000)
     *                                  = 00101101  &  11110111
     *                                  = 00100101 -> Binary
     *                                      |  | |
     *                                      5  2 0 -> Octant
    */
    
    // PROCEDURE 2:
    // ***********
    // Generate the tree from particle data
    int parCount = parPos.size();

    // **Create an evaluation container
    // Internal data container
    int unsettledParCount = parCount;       // The number of unsettled particle (initialized to all particle number)

    // Particle identification container
    std::vector<bool> parSetFlag(parCount, false);  // A settled particle flag (Said that particle has settled in the current cell)
    std::vector<int> parCellID(parCount, 0);        // Cell ID of the corresponding particle (At initial will point to root cell of 0)

    // **Loop until all particle is settled into the cell leaf
    int currLevel = 0;  // The level indicator at current loop evaluation
    while (unsettledParCount > 0){
        // Proceed to the next level calculation (proceed at first, because root is already done)
        currLevel++;

        // Update the starting cell ID at the current level
        this->startID.push_back(this->cellNum);

        // **Move all unsettled particle into the corresponding child cell
        for (int parID = 0; parID < parCount; parID++){
            // Check whether the current particle is already settled, skip if yes
            if (parSetFlag[parID]) continue;

            // Aliasing the ID
            int cellID = parCellID[parID];  // Alias to the cell ID where the particle is currently located

            // Find the child ID of current cell
            int octant = this->getOctant(parPos, parID, cellID);

            // <!> Create the child cell if still not existed
            bool isChdExist = ((this->chdFlag[cellID] >> octant) % 2);      // A mark for existed child
            if (isChdExist == false){
                // Create the child cell
                int chdID = this->cellNum;      // The ID of new child (will take the current size of cell tree)

                // Calculate the child cell property
                double chd_len = this->size[cellID] / 2.0;   /* Get Index */
                double chd_xc = this->xc[cellID] + chd_len*( ((octant>>0)%2) - 0.5);
                double chd_yc = this->yc[cellID] + chd_len*( ((octant>>1)%2) - 0.5);
                double chd_zc = this->zc[cellID] + chd_len*( ((octant>>2)%2) - 0.5);

                // Assign the data into the container
                this->cellNum++;
                this->xc.push_back(chd_xc);
                this->yc.push_back(chd_yc);
                this->zc.push_back(chd_zc);
                this->size.push_back(chd_len);
                this->parentID.push_back(cellID);   // Take the current cell as the child parent
                this->chdFlag.push_back(0);         // Child octant flag (initially no child, set to zero)
                this->nPar.push_back(0);            // Number of all particle in current cell (set to zero, will be added later)
                this->isActive.push_back(false);    // The active sign of new child is set to false
                // Update the root utility data
                this->parID.push_back(std::vector<int>());   // List of the first 'critNum' particle ID in the cell
	            this->chdID.push_back(std::vector<int>(this->chdNum,0));    // List of child cell ID at each octant

                
                // Update the data of current cell related to the child
                this->chdID[cellID][octant] = chdID;        // Update the child ID of current cell (at corresponding octant)
                this->chdFlag[cellID] |= (1 << octant);     // Turn on the child flag at current octant
            }

            // Update the current particle data in the identification container
            int chdID = this->chdID[cellID][octant];    // The ID of current child
            parCellID[parID] = chdID;                   // Change the cell ID container to the child ID

            // Update the child data
            this->nPar[chdID]++;
        }

        
        // **Update the settle flag container
        for (int parID = 0; parID < parCount; parID++){
            // Check whether the current particle is already settled, skip if yes
            if (parSetFlag[parID]) continue;

            // Aliasing the cell ID
            int cellID = parCellID[parID];
            
            // Check whether the current cell is leaf or not
            if (this->nPar[cellID] < this->critNum){
                // The current cell is a leaf cell
                // Put the current particle into the particle ID list
                this->parID[cellID].push_back(parID);
                
                // Evaluate the particle active sign
                if (activeFlag[parID]) this->isActive[cellID] = true;
                
                // Update the identification criteria
                parSetFlag[parID] = true;
                unsettledParCount--;
            }
        }
    }
    
    // Update the maximum level in this tree data
    this->maxLevel = currLevel;

    return;
}

// The current subroutine is the code refactor from the old code

void FMMCell::save_all_tree(FMMCell treeData, std::string name){
    // Initialization of saving data
    std::ofstream writer;
    name = "output/treeCellT2_all_" + name + ".csv";
    writer.open(name.c_str());
    
    // Write data table header
    writer << "ID,x,y,z,size,level,isLeaf,isActive,parNum\n";

    // Write data value
    for(int ID = 0; ID < this->cellNum; ID++){
        writer <<  "" << ID
               << "," << treeData.xc[ID]
               << "," << treeData.yc[ID]
               << "," << treeData.zc[ID]
               << "," << treeData.size[ID]
               << "," << treeData.getCellLevel(treeData.size[ID])
               << "," << (treeData.nPar[ID] < this->critNum)
               << "," << treeData.isActive[ID]
               << "," << treeData.nPar[ID]
               << "\n";
    }
    
    writer.close();
    return;
}

void FMMCell::save_leaf_tree(FMMCell treeData, std::string name){
    // Initialization of saving data
    std::ofstream writer;
    name = "output/treeCellT2_leaf_" + name + ".csv";
    writer.open(name.c_str());
    
    // Write data table header
    writer << "ID,x,y,z,size,level,isActive,parNum\n";

    // Write data value
    for(int ID = 0; ID < this->cellNum; ID++){
    if (treeData.nPar[ID] < this->critNum){
        writer <<  "" << ID
               << "," << treeData.xc[ID]
               << "," << treeData.yc[ID]
               << "," << treeData.zc[ID]
               << "," << treeData.size[ID]
               << "," << treeData.getCellLevel(treeData.size[ID])
               << "," << treeData.isActive[ID]
               << "," << treeData.nPar[ID]
               << "\n";
    }}
    
    writer.close();
    return;
}