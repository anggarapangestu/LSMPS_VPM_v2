#include "treeCell.hpp"

// =====================================================
// +------------------ Public Method ------------------+
// =====================================================
// #pragma region PUBLIC_METHOD

void treeCell::initializeTree(std::vector<std::vector<double>>& parPos){
    /* To do list:
       Pre 1. Calculate and determine the BASIC PAREMETERs
       Pre 2. Resize the cell data variable
       1. Determine the extreme particle position
       2. Set the base length
       3. Determine the root center position
       4. Create Tmap
    */

    // Internal variable
    double* minPos = new double [DIM];
    double* maxPos = new double [DIM];

    // PROCEDURE PRE-1! : Calculating and determine the BASIC PARAMETERs
    // ****************
    this->basis = std::pow(2,DIM);
    this->max_level = Pars::tree_lvl_max;
    this->min_level = Pars::tree_lvl_min;
    this->max_par = Pars::src_count_max;
    
    // Calculate the number of cell using geometry series (S_n = a(r^n-1)/(r-1))
    int _a = 1;
    int _r = this->basis;
    int _n = this->max_level + 1;       // Tree start from 0, yet geometry series start from 1
    this->cellNum = _a * (std::pow(_r,_n) - 1) / (_r - 1);

    MESSAGE_LOG << "LOCAL TREE DATA :\n";
    printf("BASIS    : %d\n", this->basis);
    printf("CELL NUM : %d\n", this->cellNum);


    // PROCEDURE PRE-2! : Resize the cell data variable
    // ****************
    this->cellPos.clear();
    this->parIDList.clear();
    this->leafCellMark.clear();
    this->outsideCell.clear();
    this->leafCellList.clear();
    this->level.clear();
    this->parNum.clear();
    this->srcNum.clear();

    this->cellPos.resize(this->cellNum,std::vector<double>(DIM,0.0));    // The center position of the cell
    this->parIDList.resize(this->cellNum,std::vector<int>());            // List of all particle ID in the current cell
    this->leafCellMark.resize(this->cellNum,false);     // The marker whether a leaf cell
    this->outsideCell.resize(this->cellNum,false);      // The marker whether the cell is existed but outside the domain
    this->level.resize(this->cellNum,0);      // The cell level hierarchy of the cell
    this->parNum.resize(this->cellNum,0);     // The number of all particel inside the cell
    this->srcNum.resize(this->cellNum,0);     // The number of source particel inside the cell

    // PROCEDURE 1! : Find the maximum and the minimum position
    // ************
    for (int e = 0; e < DIM; e ++){
        minPos[e] = parPos[0][e];
        maxPos[e] = parPos[0][e];
    }

    for (size_t i = 1; i < parPos.size(); i ++){
        for (int e = 0; e < DIM; e ++){
            minPos[e] = minPos[e] < parPos[i][e] ? minPos[e] : parPos[i][e];
            maxPos[e] = maxPos[e] > parPos[i][e] ? maxPos[e] : parPos[i][e];
        }
    }

    // PROCEDURE 2! : Set the base length as the longest range at each dimension direction
    // ************
    this->baseLen = 0.0;
    for (int e = 0; e < DIM; e ++){
        this->baseLen = this->baseLen > (maxPos[e] - minPos[e]) ? this->baseLen : (maxPos[e] - minPos[e]);
    }
    
    // Adjust the length to be expanded with the given scale factor
    this->expansion = 1.0 + (Pars::expTree / 100.0);
    this->baseLen *= this->expansion;

    // PROCEDURE 3! : Determine the root center position
    // ************
    std::vector<double> cenPos(DIM,0.0);
    for (int e = 0; e < DIM; e ++){
        cellPos[0][e] = (maxPos[e] + minPos[e]) / 2.0;
    }

    // PROCEDURE 4! : Create the tree map
    // ************
    switch (DIM)
    {
    case 2: // 2 Dimension
        this->createTmap2D();
        break;
    case 3: // 3 Dimension
        this->createTmap3D();
        break;
    default:
        break;
    }

    // [END] Release the memory
    delete[] maxPos, minPos;
    return;
}

void treeCell::createTree(std::vector<std::vector<double>>& parPos, std::vector<bool> & activeMark){
    /* To do list: - Perform the cell division
       0. Evaluate the cell level
       1. Calculate the cell center coordinate
       2. Evaluate the enclosed particle at each child

       Desc:
       - Update the level, parNum, parIDList, and cellPos data variable
    */

    // // Dividing parameter
    // bool divide_flag;

    // Data indexing variable
    int CURR_ID = 0;            // The current evaluated child cell ID (start from the root)
    std::vector<int> CHILD_ID;  // The list of child ID
    int lvl = 0;                // The current evaluated cell level (start from root)
    int ID_first = 1;           // The strating cell ID at the next level
    
    // Temporary variable
    std::vector<int> bin_ID(DIM,0);       // The binary index
    std::vector<double> _pos(DIM,0.0);    // The intermediate center position

    // Note
    // - Cell division is done to the parent cell only
    // - At each evaluation it only calculate the current cell and the child cell
    // - The result of the division at each stage is stored to the cell child

    while (true){
        // STEP 0! : Evaluate the cell level
        // *******
        // The level evaluation for the current_ID
        if (CURR_ID == ID_first){
            lvl++;
            ID_first += std::pow(this->basis,lvl);  // The starting ID for next level
        }

        // Exceeding the level
        if (lvl == this->max_level + 1){break;}
        
        // For the max level cell put as the leaf cell if particle existed
        if (lvl == this->max_level){
            if (parNum[CURR_ID] > 0){
                // This cell supposed to be the leaf cell
                this->leafCellList.push_back(CURR_ID);
                this->leafCellMark[CURR_ID] = true;
            }
            CURR_ID++;
            continue;
        }
        
        // STEP 1! : Calculate the child cell center position
        // *******
        // Divide the cell indexing using a binary to decimal method
        
        // Note:
        //  The binary index at bin_ID (e.g. in 3D: 0,1,0)
        //  is pointing to decimal of child_ID at the corresponding index i

        CHILD_ID = this->findChildID(CURR_ID);      // Take the child ID list
        bin_ID = std::vector<int> (DIM,0);    // Reset the binary index
        for (int i = 0; i < this->basis; i++){
            // Update the value of binary index
            for (int e = 0; e < DIM - 1; e++){
                bin_ID[e+1] += bin_ID[e] / 2;    // Add the higher digit
                bin_ID[e] = bin_ID[e] % 2;      // Reset to 0 after reaching 2
            }
            
            // The calculation of each cell position goes here
            for (int e = 0; e < DIM; e++){
                this->cellPos[CHILD_ID[i]][e] = this->cellPos[CURR_ID][e] + 
                            (-1 + 2*bin_ID[e]) * this->baseLen / (std::pow(2,lvl) * 4.0);
            }
            
            // Update the level
            this->level[CHILD_ID[i]] = lvl + 1;

            // Proceed the binary digit
            bin_ID[0]++;
        }

        // STEP 2! : Evaluate the enclosed particle at each child
        // *******
        // divide_flag = true;
        if (CURR_ID == 0){
            // Group the particle from the root cell into the new child cell
            for (size_t i = 0; i < parPos.size(); i++){   
                // The target child cell ID
                int target_cellID = CHILD_ID[0];

                // Determine the new cell location, by check on each dimension direction
                for (int e = 0; e < DIM; e++){
                    if (parPos[i][e] > this->cellPos[CURR_ID][e]){
                        target_cellID += std::pow(2,e);
                    }
                }
                
                // Update the particle id then add into the cell list
                parIDList[target_cellID].push_back(i);
                parNum[target_cellID]++;        // Update the number of particle inside the target cell ID
                if (activeMark[i] == true){
                    srcNum[target_cellID]++;    // Update the number of source particle inside the target cell ID
                }
            }
            
            // If something goes wrong
            for (auto _ID:CHILD_ID){
                if (parNum[_ID] != int(parIDList[_ID].size())){
                    std::cout << "[WARNING] : The number of particle is different!\n";
                }
            }

        }else{
            // The group of particle division into child cell
            // Check whether the current cell is leaf or not
            if ((this->srcNum[CURR_ID] < this->max_par) && (this->level[CURR_ID] > this->min_level)){
                // Only if the particle are existed
                if (parNum[CURR_ID] > 0){
                    // This cell supposed to be the leaf cell
                    this->leafCellList.push_back(CURR_ID);
                    this->leafCellMark[CURR_ID] = true;
                }
            }
            else{
                // This cell supposed to be divided
                for (auto parID : parIDList[CURR_ID]){   
                    // The target child cell ID
                    int target_cellID = CHILD_ID[0];

                    // Determine the new cell location, by check on each dimension direction
                    for (int e = 0; e < DIM; e++){
                        if (parPos[parID][e] > this->cellPos[CURR_ID][e]){
                            target_cellID += std::pow(2,e);     // Idea got from the binary to decimal
                        }
                    }
                    
                    // Update the particle id then add into the cell list
                    parIDList[target_cellID].push_back(parID);
                    parNum[target_cellID]++;        // Update the number of particle inside the target cell ID
                    if (activeMark[parID]){
                        srcNum[target_cellID]++;    // Update the number of source particle inside the target cell ID
                    }
                }
                
                // If something goes wrong
                for (auto _ID:CHILD_ID){
                    if (parNum[_ID] == 0){
                        this->outsideCell[_ID] = true;
                    }
                    if (parNum[_ID] != int(parIDList[_ID].size())){
                        std::cout << "[WARNING] : The number of particle is different!\n";
                    }
                }
            }
        }

        // Proceed the cell to the next sequence
        CURR_ID++;
    }

    // std::cout << "Number of cell : " << this->cellNum << "\n";
    
    return;
}

void treeCell::updateCell(std::vector<std::vector<double>>& parPos, std::vector<bool> & activeMark){
    // Check all particle in the domain by check only the leaf cell
    std::vector<int> newLeafCellList;   // The ID list of new leaf cell

    // Perform the create tree division for the leaf cell
    for (size_t i = 0; i < this->leafCellList.size(); i++){
        int cellID = this->leafCellList[i];
        // std::cout << "At iteration " << i << "; ID :" << cellID <<"\n";
        int newNum = 0;

        // Check every particle inside the current cell
        for (size_t j = 0; j < this->parIDList[cellID].size(); j++){
            int parID = this->parIDList[cellID][j];
            if (activeMark[parID] == true){
                newNum++;
            }
        }
        
        // Update the number of source particle (We don't need to update the number at parent)
        this->srcNum[cellID] = newNum;

        // std::cout << "[DEBUG LINE] Finish the number calculation";

        // Divide the cell if (1) not the max level cell and (2) "srcnum" bigger than threshold
        if ((this->level[cellID] < this->max_level) && (newNum > this->max_par)){
            // std::cout << " --> type 1\n";
            // The current cell is no longer a leaf cell, delete the leaf cell mark
            this->leafCellMark[cellID] = false;
            
            // Find each child cell
            std::vector<int> CHILD_ID = this->findChildID(cellID);      // Take the child ID list
            std::vector<int> bin_ID = std::vector<int> (DIM,0);   // Reset the binary index

            // This cell supposed to be divided
            for (auto parID : parIDList[cellID]){   
                // The target child cell ID
                int target_cellID = CHILD_ID[0];

                // Determine the new cell location, by check on each dimension direction
                for (int e = 0; e < DIM; e++){
                    if (parPos[parID][e] > this->cellPos[cellID][e]){
                        target_cellID += std::pow(2,e);     // Idea got from the binary to decimal
                    }
                }
                
                // Update the particle id then add into the cell list
                parIDList[target_cellID].push_back(parID);
                parNum[target_cellID]++;        // Update the number of particle inside the target cell ID
                if (activeMark[parID]){
                    srcNum[target_cellID]++;    // Update the number of source particle inside the target cell ID
                }

            }
            
            // Check each child cell
            for (auto _ID : CHILD_ID){
                // Update the new leaf cell list and mark
                this->leafCellMark[_ID] = true;
                newLeafCellList.push_back(_ID);
                
                // If something goes wrong
                if (parNum[_ID] != int(parIDList[_ID].size())){
                    std::cout << "[WARNING] : The number of particle is different!\n";
                }
            }
        }else{
            // std::cout << " --> type 2\n";
            // Update the new leaf cell list and mark
            this->leafCellMark[cellID] = true;
            newLeafCellList.push_back(cellID);
        }

        // std::cout << "[DEBUG LINE] Finish the division\n";
        
    }
    
    // Update the value of leaf cell
    this->leafCellList.resize(newLeafCellList.size());
    for (size_t i = 0; i < newLeafCellList.size(); i++){
        this->leafCellList[i] = newLeafCellList[i];
    }

    return;
}

// #pragma endregion

// =====================================================
// +----------------- Neighbor Search -----------------+
// =====================================================
// #pragma region NEIGHBOR_SEARCH
void treeCell::intList(int currID, std::vector<int>& ID_list_1, std::vector<int>& ID_list_3) const{
    // Reset the neighbor list
    ID_list_1.clear();
    ID_list_3.clear();
    
    // Procedure:
    // 1. List the neighbor cell (at the same level) that touched with current cell
    // 2. Compare each cell box to the list criteria
    //    > List 1 : Leaf cell, touched with the current cell
    //    > List 3 : Leaf cell, not touched with the current cell

    // Internal variable
    int currlvl = this->level[currID];      // The level of current cell
    int ngh_num = std::pow(3,DIM);    // Max number of initial neighbor cell list
    int mat_size = std::pow(2,currlvl);     // The size of transformation matrix
    std::vector<int> indexPos = this->cellIdx[currID];  // The matrix index of the current cell
    std::vector<int> tri_ID;                // Temporary trinary index
    std::vector<int> ID_list_1_temp;        // The temporary neighbor ID list 1
    std::vector<int> ID_list_3_temp;        // The temporary neighbor ID list 3

    // // [DEBUG LINE]
    // if (currID == 9){
    //     std::cout << "Start finding the neighbor list\n";
    //     int coin=0;
    //     for (auto val : indexPos){
    //         coin++;
    //         std::cout << "P" << coin << ":ID" << val << " | ";
    //     }
    //     std::cout << std::endl;
    // }

    // PROCEDURE 1! : Find all neighbor cell ID touched with current cell
    // ************
    tri_ID = std::vector<int>(DIM,0); // Reset the trinary index
    for (int i = 0; i < ngh_num; i ++){
        // Update the trinary index
        for (int e = 0; e < DIM - 1; e ++){
            tri_ID[e+1] += tri_ID[e] / 3;
            tri_ID[e] = tri_ID[e] % 3;
        }
        
        // Internal variable (IDEA: (tri_ID - 1) create set{-1, 0, 1})
        int IDpos1 = indexPos[0] + tri_ID[0] - 1;
        int IDpos2 = indexPos[1] + tri_ID[1] - 1;
        int nghID = currID;     // At default will be thrown out

        // Check whether the neighbor is existed (current cell not located at boundary)
        if ((IDpos1 < mat_size && IDpos1 >= 0) && 
            (IDpos2 < mat_size && IDpos2 >= 0))
        {
            nghID = this->T_map2D[currlvl][IDpos1][IDpos2];
        }
        
        // Update the neighbor ID only if it is not the current cell ID
        if (nghID != currID){
            ID_list_1_temp.push_back(nghID);
        }

        // Proceed the trinary index
        tri_ID[0]++;        
    }
    
    // PROCEDURE 2! : Evaluate the criteria list
    // ************
    // Criteria: 
    // - Leaf cell [Y] -> Take the ID into the "final list 1"
    // - Leaf cell [N] -> Check the leaf cell located at (parent or child cell)
    //   - Located at parent -> Again check whether a leaf cell [Y] -> check if no duplicate
    //   - Located at child  -> Check the position touching with current cell
    //     - Touching [Y] -> Again check whether a leaf cell [Y] -> Put into the "final list 1"
    //     - Touching [N] -> Take all child ID (until reach the leaf level) into the "final list 3"
    
    // // [DEBUG LINE]
    // if (currID == 9){
    //     std::cout << "Start to adding the component\n";
    // }
    // // int countter = 0;

    for (size_t i = 0; i < ID_list_1_temp.size(); i++){
        int nghID = ID_list_1_temp[i];
        
        // // [DEBUG LINE]
        // // if (countter == 20){return;}
        // if (currID == 9){
        //     std::cout << "Curent Iterated ID : "<< nghID << "\n";
        //     // countter ++;
        // }

        // Skip the candidate neighbor located out the domain
        if (this->outsideCell[nghID] == true){
            continue;
        }

        // [FIRST] Check whether the cell is leaf cell or not
        if (this->leafCellMark[nghID] == true){
            // The current cell is a leaf cell
            
            // // [DEBUG LINE]
            // if (currID == 9){
            //     std::cout << "TYPE 1\n";
            // }
            
            // Check whether it duplicate the ID at "final list 1"
            bool _take = true;
            for (auto _ID:ID_list_1){
                if (nghID == _ID){
                    _take = false;
                }
            }
            // Add only if there is no duplicate
            if (_take == true){
                ID_list_1.push_back(nghID);
            }
        }
        else{
            // Determine the location of leaf cell (parent or child cell)
            if(this->parNum[nghID] == 0){   // If there is no particle, the leaf cell must be located at the parent cell
                // [Take the parent cell]
                // (1) Put into the temporary list 1
                
                // // [DEBUG LINE]
                // if (currID == 9){
                //     std::cout << "TYPE 2\n";
                // }

                int _parID = this->findParentID(nghID);
                ID_list_1_temp.push_back(_parID);
            }
            else{
                // [Take the child cell]
                // (1) Check whether is touching current cell
                // (2) Put into corresponding temporary list (1 or 3)
                
                // // [DEBUG LINE]
                // if (currID == 9){
                //     std::cout << "TYPE 3\n";
                // }
                std::vector<int> _chdIDs = this->findChildID(nghID);

                // Check whether it is touching the current cell
                double min_dist = this->baseLen * ( 1/std::pow(2.0,this->level[_chdIDs[0]]) + 
                                  1/std::pow(2.0,currlvl) ) / 2.0 + (10 * __DBL_EPSILON__);
                for (auto _chdID:_chdIDs){
                    bool _sep = false;
                    // Check the distance at each dimension direction
                    for (int e = 0; e < DIM; e ++){
                        // If not touched to current cell (seperated), will not be taken
                        if (std::abs(this->cellPos[_chdID][e] - this->cellPos[currID][e]) > min_dist){
                            _sep = true;
                            break;
                        }
                    }
                    // Final conditional
                    if (_sep == true){
                        // If seperated it will be the list 3
                        ID_list_3.push_back(_chdID);
                    }else{
                        // If touched move into temp list 1
                        ID_list_1_temp.push_back(_chdID);
                    }
                }
                
            }
        }
    }
    
    // // [DEBUG LINE]
    // if (currID == 9){
    //     std::cout << "Adding the component\n";
    // }
    
    // Add the current cell ID to the list
    ID_list_1.push_back(currID);

    // // Evaluate the rest (temporary list 3)
    // for (int i = 0; i < ID_list_3_temp.size(); i ++){
    //     int nghID = ID_list_3_temp[i];
        
    //     // Check whether it is leaf cell
    //     if (this->leafCellMark[nghID] == true){
    //         // The current cell is a leaf cell
    //         ID_list_3.push_back(nghID);
    //     }
    //     else{
    //         // Check whether there is a source inside the 

    //         // Take the child ID into temporary list 3
    //         std::vector<int> _chdIDs = this->findChildID(nghID);
    //         for (auto _chdID:_chdIDs){
    //             ID_list_3_temp.push_back(_chdID);
    //         }
    //     }
    // }
    // std::cout << "[DEBUG LINE] DONE ALl\n";
    
    return;
}

void treeCell::intList_3d(int currID, std::vector<int>& ID_list_1, std::vector<int>& ID_list_3) const{
    // Reset the neighbor list
    ID_list_1.clear();
    ID_list_3.clear();
    
    // Procedure:
    // 1. List the neighbor cell (at the same level) that touched with current cell
    // 2. Compare each cell box to the list criteria
    //    > List 1 : Leaf cell, touched with the current cell
    //    > List 3 : Leaf cell, not touched with the current cell

    // Internal variable
    int currlvl = this->level[currID];      // The level of current cell
    int ngh_num = std::pow(3,DIM);    // Max number of initial neighbor cell list
    int mat_size = std::pow(2,currlvl);     // The size of transformation matrix
    std::vector<int> indexPos = this->cellIdx[currID];  // The matrix index of the current cell
    std::vector<int> tri_ID;                // Temporary trinary index
    std::vector<int> ID_list_1_temp;        // The temporary neighbor ID list 1
    std::vector<int> ID_list_3_temp;        // The temporary neighbor ID list 3

    // // [DEBUG LINE]
    // if (currID == 9){
    //     std::cout << "Start finding the neighbor list\n";
    //     int coin=0;
    //     for (auto val : indexPos){
    //         coin++;
    //         std::cout << "P" << coin << ":ID" << val << " | ";
    //     }
    //     std::cout << std::endl;
    // }

    // PROCEDURE 1! : Find all neighbor cell ID touched with current cell
    // ************
    tri_ID = std::vector<int>(DIM,0); // Reset the trinary index
    for (int i = 0; i < ngh_num; i ++){
        // Update the trinary index
        for (int e = 0; e < DIM - 1; e ++){
            tri_ID[e+1] += tri_ID[e] / 3;
            tri_ID[e] = tri_ID[e] % 3;
        }
        
        // Internal variable (IDEA: (tri_ID - 1) create set{-1, 0, 1})
        int IDpos1 = indexPos[0] + tri_ID[0] - 1;
        int IDpos2 = indexPos[1] + tri_ID[1] - 1;
        int IDpos3 = indexPos[2] + tri_ID[2] - 1;
        int nghID = currID;     // At default will be thrown out

        // Check whether the neighbor is existed (current cell not located at boundary)
        if ((IDpos1 < mat_size && IDpos1 >= 0) && 
            (IDpos2 < mat_size && IDpos2 >= 0) && 
            (IDpos3 < mat_size && IDpos3 >= 0))
        {
            nghID = this->T_map3D[currlvl][IDpos1][IDpos2][IDpos3];
        }
        
        // Update the neighbor ID only if it is not the current cell ID
        if (nghID != currID){
            ID_list_1_temp.push_back(nghID);
        }

        // Proceed the trinary index
        tri_ID[0]++;        
    }
    
    // PROCEDURE 2! : Evaluate the criteria list
    // ************
    // Criteria: 
    // - Leaf cell [Y] -> Take the ID into the "final list 1"
    // - Leaf cell [N] -> Check the leaf cell located at (parent or child cell)
    //   - Located at parent -> Again check whether a leaf cell [Y] -> check if no duplicate
    //   - Located at child  -> Check the position touching with current cell
    //     - Touching [Y] -> Again check whether a leaf cell [Y] -> Put into the "final list 1"
    //     - Touching [N] -> Take all child ID (until reach the leaf level) into the "final list 3"
    
    // // [DEBUG LINE]
    // if (currID == 9){
    //     std::cout << "Start to adding the component\n";
    // }
    // // int countter = 0;

    for (size_t i = 0; i < ID_list_1_temp.size(); i++){
        int nghID = ID_list_1_temp[i];
        
        // // [DEBUG LINE]
        // // if (countter == 20){return;}
        // if (currID == 9){
        //     std::cout << "Curent Iterated ID : "<< nghID << "\n";
        //     // countter ++;
        // }

        // Skip the candidate neighbor located out the domain
        if (this->outsideCell[nghID] == true){
            continue;
        }

        // [FIRST] Check whether the cell is leaf cell or not
        if (this->leafCellMark[nghID] == true){
            // The current cell is a leaf cell
            
            // // [DEBUG LINE]
            // if (currID == 9){
            //     std::cout << "TYPE 1\n";
            // }
            
            // Check whether it duplicate the ID at "final list 1"
            bool _take = true;
            for (auto _ID:ID_list_1){
                if (nghID == _ID){
                    _take = false;
                }
            }
            // Add only if there is no duplicate
            if (_take == true){
                ID_list_1.push_back(nghID);
            }
        }
        else{
            // Determine the location of leaf cell (parent or child cell)
            if(this->parNum[nghID] == 0){   // If there is no particle, the leaf cell must be located at the parent cell
                // [Take the parent cell]
                // (1) Put into the temporary list 1
                
                // // [DEBUG LINE]
                // if (currID == 9){
                //     std::cout << "TYPE 2\n";
                // }

                int _parID = this->findParentID(nghID);
                ID_list_1_temp.push_back(_parID);
            }
            else{
                // [Take the child cell]
                // (1) Check whether is touching current cell
                // (2) Put into corresponding temporary list (1 or 3)
                
                // // [DEBUG LINE]
                // if (currID == 9){
                //     std::cout << "TYPE 3\n";
                // }
                std::vector<int> _chdIDs = this->findChildID(nghID);

                // Check whether it is touching the current cell
                double min_dist = this->baseLen * ( 1/std::pow(2.0,this->level[_chdIDs[0]]) + 
                                  1/std::pow(2.0,currlvl) ) / 2.0 + (10 * __DBL_EPSILON__);
                for (auto _chdID:_chdIDs){
                    bool _sep = false;
                    // Check the distance at each dimension direction
                    for (int e = 0; e < DIM; e ++){
                        // If not touched to current cell (seperated), will not be taken
                        if (std::abs(this->cellPos[_chdID][e] - this->cellPos[currID][e]) > min_dist){
                            _sep = true;
                            break;
                        }
                    }
                    // Final conditional
                    if (_sep == true){
                        // If seperated it will be the list 3
                        ID_list_3.push_back(_chdID);
                    }else{
                        // If touched move into temp list 1
                        ID_list_1_temp.push_back(_chdID);
                    }
                }
                
            }
        }
    }
    
    // // [DEBUG LINE]
    // if (currID == 9){
    //     std::cout << "Adding the component\n";
    // }
    
    // Add the current cell ID to the list
    ID_list_1.push_back(currID);

    // // Evaluate the rest (temporary list 3)
    // for (int i = 0; i < ID_list_3_temp.size(); i ++){
    //     int nghID = ID_list_3_temp[i];
        
    //     // Check whether it is leaf cell
    //     if (this->leafCellMark[nghID] == true){
    //         // The current cell is a leaf cell
    //         ID_list_3.push_back(nghID);
    //     }
    //     else{
    //         // Check whether there is a source inside the 

    //         // Take the child ID into temporary list 3
    //         std::vector<int> _chdIDs = this->findChildID(nghID);
    //         for (auto _chdID:_chdIDs){
    //             ID_list_3_temp.push_back(_chdID);
    //         }
    //     }
    // }
    // std::cout << "[DEBUG LINE] DONE ALl\n";
    
    return;
}

void treeCell::extList(int currID, std::vector<int>& ID_list_2, std::vector<int>& ID_list_4) const{
    // Reset the neighbor list
    ID_list_2.clear();
    ID_list_4.clear();
    
    // Procedure:
    // 1. List the neighbor cell (at the parent same level) that touched with the parent of current cell
    // 2. Compare each cell box to the list criteria
    //    > List 2 : Same level to current, well separate to current cell (separated by 1 cell apart)
    //    > List 4 : Lower level to current, leaf cell, well separated to current cell

    // Internal variable
    int currlvl = this->level[currID];      // The level of current cell
    int ngh_num = std::pow(3,DIM);    // Max number of candidate neighbor cell
    int parID = this->findParentID(currID); // The current cell parent ID
    int mat_size = std::pow(2,currlvl - 1); // The size of transformation matrix at parent level
    std::vector<int> indexPos = this->cellIdx[parID];  // The matrix index of parent cell
    std::vector<int> tri_ID;            // Temporary trinary index
    std::vector<int> ID_list_2_temp;    // The temporary neighbor ID list 2
    std::vector<int> ID_list_4_temp;    // The temporary neighbor ID list 4

    // Procedure 1! : List all parent cell neighbor ID
    // ************
    // std::cout << "[DEBUG LINE] Start Procedure at current ID: " << currID << ", lv: "<< currlvl << "\n";
    tri_ID = std::vector<int>(DIM,0); // Reset the trinary index
    for (int i = 0; i < ngh_num; i ++){
        // Update the trinary index
        for (int e = 0; e < DIM - 1; e ++){
            tri_ID[e+1] += tri_ID[e] / 3;
            tri_ID[e] = tri_ID[e] % 3;
        }
        
        // Internal variable
        int IDpos1 = indexPos[0] + tri_ID[0] - 1;
        int IDpos2 = indexPos[1] + tri_ID[1] - 1;
        int nghID = parID;

        // Check whether the neighbor is existed (current cell not located at boundary)
        if ((IDpos1 < mat_size && IDpos1 >= 0) && 
            (IDpos2 < mat_size && IDpos2 >= 0))
        {
            nghID = this->T_map2D[currlvl - 1][IDpos1][IDpos2];
        }
        
        // Update the neighbor ID only if it is not the current cell ID
        if (nghID != parID){
            ID_list_2_temp.push_back(nghID);
        }

        // Proceed the trinary index
        tri_ID[0]++;        
    }
    
    // std::cout << "[DEBUG LINE] Done Procedure 1\n";
    // Procedure 2! : Evaluate the list cretiria
    // ************
    // - Leaf cell [Y] - Put as candidate list 4
    //   - Well separated [Y] - Take as the final list 4
    //   - Well separated [N] - Throw the ID
    // - Leaf cell [N] - Take the child -> check well separated
    //   - Well separated [Y] - Take as the final list 2
    //   - Well separated [N] - Throw the ID
    // - Have no source particle - Take the parent ID into temporary list 4

    for (size_t i = 0; i < ID_list_2_temp.size(); i++){
        int nghID = ID_list_2_temp[i];
        // std::cout << " Iteration : " << i << "; NGH_ID : " << nghID;//
    
        // Check whether it is leaf cell or have no source particle
        if ((this->leafCellMark[nghID] == true) || (this->srcNum[nghID] == 0)){
            // Candidate for list 4
            ID_list_4_temp.push_back(nghID);
        }
        else{
            // Check each child cell whether is well separated to the current cell (not the parent cell)
            std::vector<int> _chdIDs = this->findChildID(nghID);

            // Check the relative position of each child cell toward current cell
            double min_dist = this->baseLen / std::pow(2,currlvl) + (10 * __DBL_EPSILON__);
            for (auto _chdID:_chdIDs){
                bool _sep = false;    // Initially the child ID is considered to be touching current cell
                for (int e = 0; e < DIM; e ++){
                    // if not touching current cell will be taken
                    if (std::abs(this->cellPos[_chdID][e] - this->cellPos[currID][e]) > min_dist){
                        _sep = true;
                        break;
                    }
                }
                // Take the seperated ID to final list 2
                if (_sep == true){
                    ID_list_2.push_back(_chdID);
                }
            }
        }
    }

    // // DEBUG LINE
    // std::cout << " >> LIST 2: ";
    // for (auto val : ID_list_2){
    //     std::cout << val << ",";
    // }
    // std::cout << "\n[DEBUG LINE] Done Procedure 2\n";

    // Evaluate the rest (temporary list 4)
    for (size_t i = 0; i < ID_list_4_temp.size(); i ++){
        int nghID = ID_list_4_temp[i];
        // std::cout << " Iteration : " << i << "; NGH_ID : " << nghID << "\n";//
        
        // Check whether the cell is well seperated to the current cell
        // The current cell is a leaf cell
        double min_dist = this->baseLen * (1/std::pow(2,this->level[nghID]) +
                            1/std::pow(2,currlvl)) / 2.0 + (10 * __DBL_EPSILON__);

        // Check whether is seperated to current cell
        bool _sep = false;
        for (int e = 0; e < DIM; e ++){
            if (std::abs(this->cellPos[nghID][e] - this->cellPos[currID][e]) > min_dist){
                _sep = true;
                break;
            }
        }
        
        // Only consider the box that is well seperated
        if(_sep == true){
            // Check whether it is leaf cell
            if (this->leafCellMark[nghID] == true){
                // Make sure no duplicate
                bool _take = true;
                for (auto _ID:ID_list_4){
                    if (nghID == _ID){
                        _take = false;
                    }
                }
                if (_take == true){
                    ID_list_4.push_back(nghID);
                }
            }
            else{
                // Take the parent ID
                int _parID = this->findParentID(nghID);
                ID_list_4_temp.push_back(_parID);
            }
            
        }

    }

    // // DEBUG LINE
    // std::cout << " >> LIST 4: ";
    // for (auto val : ID_list_4){
    //     std::cout << val << ",";
    // }

    // std::cout << "\n[DEBUG LINE] Done Procedure 3\n\n";
    return;
}

void treeCell::extList_3d(int currID, std::vector<int>& ID_list_2, std::vector<int>& ID_list_4) const{
    // Reset the neighbor list
    ID_list_2.clear();
    ID_list_4.clear();
    
    // Procedure:
    // 1. List the neighbor cell (at the parent same level) that touched with the parent of current cell
    // 2. Compare each cell box to the list criteria
    //    > List 2 : Same level to current, well separate to current cell (separated by 1 cell apart)
    //    > List 4 : Lower level to current, leaf cell, well separated to current cell

    // Internal variable
    int currlvl = this->level[currID];      // The level of current cell
    int ngh_num = std::pow(3,DIM);    // Max number of candidate neighbor cell
    int parID = this->findParentID(currID); // The current cell parent ID
    int mat_size = std::pow(2,currlvl - 1); // The size of transformation matrix at parent level
    std::vector<int> indexPos = this->cellIdx[parID];  // The matrix index of parent cell
    std::vector<int> tri_ID;            // Temporary trinary index
    std::vector<int> ID_list_2_temp;    // The temporary neighbor ID list 2
    std::vector<int> ID_list_4_temp;    // The temporary neighbor ID list 4

    // Procedure 1! : List all parent cell neighbor ID
    // ************
    // std::cout << "[DEBUG LINE] Start Procedure at current ID: " << currID << ", lv: "<< currlvl << "\n";
    tri_ID = std::vector<int>(DIM,0); // Reset the trinary index
    for (int i = 0; i < ngh_num; i ++){
        // Update the trinary index
        for (int e = 0; e < DIM - 1; e ++){
            tri_ID[e+1] += tri_ID[e] / 3;
            tri_ID[e] = tri_ID[e] % 3;
        }
        
        // Internal variable
        int IDpos1 = indexPos[0] + tri_ID[0] - 1;
        int IDpos2 = indexPos[1] + tri_ID[1] - 1;
        int IDpos3 = indexPos[2] + tri_ID[2] - 1;
        int nghID = parID;

        // Check whether the neighbor is existed (current cell not located at boundary)
        if ((IDpos1 < mat_size && IDpos1 >= 0) && 
            (IDpos2 < mat_size && IDpos2 >= 0))
        {
            nghID = this->T_map3D[currlvl - 1][IDpos1][IDpos2][IDpos3];
        }
        
        // Update the neighbor ID only if it is not the current cell ID
        if (nghID != parID){
            ID_list_2_temp.push_back(nghID);
        }

        // Proceed the trinary index
        tri_ID[0]++;        
    }
    
    // std::cout << "[DEBUG LINE] Done Procedure 1\n";
    // Procedure 2! : Evaluate the list cretiria
    // ************
    // - Leaf cell [Y] - Put as candidate list 4
    //   - Well separated [Y] - Take as the final list 4
    //   - Well separated [N] - Throw the ID
    // - Leaf cell [N] - Take the child -> check well separated
    //   - Well separated [Y] - Take as the final list 2
    //   - Well separated [N] - Throw the ID
    // - Have no source particle - Take the parent ID into temporary list 4

    for (size_t i = 0; i < ID_list_2_temp.size(); i++){
        int nghID = ID_list_2_temp[i];
        // std::cout << " Iteration : " << i << "; NGH_ID : " << nghID;//
    
        // Check whether it is leaf cell or have no source particle
        if ((this->leafCellMark[nghID] == true) || (this->srcNum[nghID] == 0)){
            // Candidate for list 4
            ID_list_4_temp.push_back(nghID);
        }
        else{
            // Check each child cell whether is well separated to the current cell (not the parent cell)
            std::vector<int> _chdIDs = this->findChildID(nghID);

            // Check the relative position of each child cell toward current cell
            double min_dist = this->baseLen / std::pow(2,currlvl) + (10 * __DBL_EPSILON__);
            for (auto _chdID:_chdIDs){
                bool _sep = false;    // Initially the child ID is considered to be touching current cell
                for (int e = 0; e < DIM; e ++){
                    // if not touching current cell will be taken
                    if (std::abs(this->cellPos[_chdID][e] - this->cellPos[currID][e]) > min_dist){
                        _sep = true;
                        break;
                    }
                }
                // Take the seperated ID to final list 2
                if (_sep == true){
                    ID_list_2.push_back(_chdID);
                }
            }
        }
    }

    // // DEBUG LINE
    // std::cout << " >> LIST 2: ";
    // for (auto val : ID_list_2){
    //     std::cout << val << ",";
    // }
    // std::cout << "\n[DEBUG LINE] Done Procedure 2\n";

    // Evaluate the rest (temporary list 4)
    for (size_t i = 0; i < ID_list_4_temp.size(); i ++){
        int nghID = ID_list_4_temp[i];
        // std::cout << " Iteration : " << i << "; NGH_ID : " << nghID << "\n";//
        
        // Check whether the cell is well seperated to the current cell
        // The current cell is a leaf cell
        double min_dist = this->baseLen * (1/std::pow(2,this->level[nghID]) +
                            1/std::pow(2,currlvl)) / 2.0 + (10 * __DBL_EPSILON__);

        // Check whether is seperated to current cell
        bool _sep = false;
        for (int e = 0; e < DIM; e ++){
            if (std::abs(this->cellPos[nghID][e] - this->cellPos[currID][e]) > min_dist){
                _sep = true;
                break;
            }
        }
        
        // Only consider the box that is well seperated
        if(_sep == true){
            // Check whether it is leaf cell
            if (this->leafCellMark[nghID] == true){
                // Make sure no duplicate
                bool _take = true;
                for (auto _ID:ID_list_4){
                    if (nghID == _ID){
                        _take = false;
                    }
                }
                if (_take == true){
                    ID_list_4.push_back(nghID);
                }
            }
            else{
                // Take the parent ID
                int _parID = this->findParentID(nghID);
                ID_list_4_temp.push_back(_parID);
            }
            
        }

    }

    // // DEBUG LINE
    // std::cout << " >> LIST 4: ";
    // for (auto val : ID_list_4){
    //     std::cout << val << ",";
    // }

    // std::cout << "\n[DEBUG LINE] Done Procedure 3\n\n";
    return;
}

// #pragma endregion

// =====================================================
// +----------------- Mapping Method ------------------+
// =====================================================
// #pragma region MAPPING_METHOD
// List of the mapping function
void treeCell::createTmap2D(){
    // To Do : Create the look up table for cell neighbor in each tree level
    // Update the size of each vector
    this->cellIdx.clear();
    this->T_map2D.clear();

    this->cellIdx.resize(this->cellNum,std::vector<int>(DIM,0));
    this->T_map2D.resize(this->max_level + 1);

    // Initialize from the root cell (ID = 0)
    this->T_map2D[0] = std::vector<std::vector<int>>(1,std::vector<int>(1,0));
    this->cellIdx[0] = std::vector<int>(DIM,0);
    
    // Calculate for the subsequent cell ID
    int* _temp_cellID_bound = new int [2]{0,0};     // Initialize by the root ID
    int _start,_fin;    // The bound cell ID each level
    int _ID = 1;        // The first cell ID at each level
    std::vector<int> c; // The index marking aid vector

    // Iterate through tree level
    for (int k = 1; k <= this->max_level; k ++){
        // Resize the t-map
        T_map2D[k].resize(std::pow(2,k),std::vector<int>(std::pow(2,k),0));

        // At each level it will iterate all cell bounded by the _temp_cellID_bound
        _start = _temp_cellID_bound[0];
        _fin = _temp_cellID_bound[1];
        
        // Iterate through the cell ID bound
        for (int _id = _start; _id <= _fin; _id++){
            // Internal varible and setup
            c = std::vector<int>(DIM,0); // Reset the counter
            int currID;

            // Iterate through the child cell of the correspoding cell ID
            for (int _ct = 0; _ct < this->basis; _ct++){
                // Update the counter
                for (int e = 0; e < DIM - 1; e ++){
                    c[e+1] += c[e] / 2;
                    c[e] = c[e] % 2;
                }

                // The current cell ID by the given basis index
                currID = _ID + _ct;

                // Update the T-map and index
                for (int e = 0; e < DIM; e++){
                    cellIdx[currID][e] = 2*cellIdx[_id][e] + c[e];
                }
                T_map2D[k][cellIdx[currID][0]][cellIdx[currID][1]] = currID;

                // Proceed the counter
                c[0] ++;
            }
            
            // Update the pivot cell ID
            _ID += this->basis;
        }
        
        // Update the cell ID bound at the next level
        _temp_cellID_bound[0] = _fin + 1;
        _temp_cellID_bound[1] = _ID - 1;
    }
    return;
}

void treeCell::createTmap3D(){
    // To Do : Create the look up table for cell neighbor in each tree level
    // Update the size of each vector
    this->cellIdx.clear();
    this->T_map3D.clear();

    this->cellIdx.resize(this->cellNum,std::vector<int>(DIM,0));
    this->T_map3D.resize(this->max_level + 1);

    // Initialize from the root cell (ID = 0)
    this->T_map3D[0] = std::vector<std::vector<std::vector<int>>>(1,std::vector<std::vector<int>>(1,std::vector<int>(1,0)));
    this->cellIdx[0] = std::vector<int>(DIM,0);
    
    // Calculate for the subsequent cell ID
    int* _temp_cellID_bound = new int [2]{0,0};     // Initialize by the root ID
    int _start,_fin;    // The bound cell ID each level
    int _ID = 1;        // The first cell ID at each level
    std::vector<int> c; // The index marking aid vector

    // Iterate through tree level
    for (int k = 1; k <= this->max_level; k ++){
        // Resize the t-map as (2^k x 2^k x 2^k)
        T_map3D[k].resize(std::pow(2,k),std::vector<std::vector<int>>(std::pow(2,k),std::vector<int>(std::pow(2,k),0)));

        // At each level it will iterate all cell bounded by the _temp_cellID_bound
        _start = _temp_cellID_bound[0];
        _fin = _temp_cellID_bound[1];
        
        // Iterate through the cell ID bound
        for (int _id = _start; _id <= _fin; _id++){
            // Internal varible and setup
            c = std::vector<int>(DIM,0); // Reset the counter
            int currID;

            // Iterate through the child cell of the correspoding cell ID
            for (int _ct = 0; _ct < this->basis; _ct++){
                // Update the counter
                for (int e = 0; e < DIM - 1; e ++){
                    c[e+1] += c[e] / 2;
                    c[e] = c[e] % 2;
                }

                // The current cell ID by the given basis index
                currID = _ID + _ct;

                // Update the T-map and index
                for (int e = 0; e < DIM; e++){
                    cellIdx[currID][e] = 2*cellIdx[_id][e] + c[e];
                }
                T_map3D[k][cellIdx[currID][0]][cellIdx[currID][1]][cellIdx[currID][2]] = currID;

                // Proceed the counter
                c[0] ++;
            }
            
            // Update the pivot cell ID
            _ID += this->basis;
        }
        
        // Update the cell ID bound at the next level
        _temp_cellID_bound[0] = _fin + 1;
        _temp_cellID_bound[1] = _ID - 1;
    }
    return;
}

// #pragma endregion

// =====================================================
// +----------------- Tree Hierarchy ------------------+
// =====================================================
// #pragma region TREE_HIERARCHY
// List of the internal function
int treeCell::findParentID(int currID) const {
    int parID;
    parID = std::floor((currID-1)/this->basis);
    return parID;
}

std::vector<int> treeCell::findChildID(int currID) const{
    std::vector<int> childIDs(this->basis,0);
    int pivID = (currID * this->basis) + 1;
    for (int i = 0; i < this->basis; i++){
        childIDs[i] = pivID + i;
    }
    return childIDs;
}

void treeCell::findLevel(const std::vector<int>& ID_list, std::vector<int>& level) const{
    int num = ID_list.size();
    int ID;
    
    // Resizing the level vector
    level.clear();
    level.resize(num,0);
    
    // Assigning the level to each corresponding ID
    for(int i = 0; i < num; i++){
        ID = ID_list[i];
        level[i] = this->level[ID];
    }
    return;
}

// #pragma endregion

// =====================================================
// +------------------- Data Saving -------------------+
// =====================================================
// #pragma region DATA_SAVING

// Saving function
void treeCell::saveTreeCell(){
    // internal variables
    std::ofstream data;

    bool _ngh = false, parIn = false;

    // Save the base cell data 
    data.open("treeCell.csv");
    data << "x,y,lvl,matIDx,matIDy,num,srcnum,isleaf,isoutside,size";
    if (parIn)  {data << ",parID\n";}
    else        {data << "\n";}

    for (int i = 0; i < this->cellNum; i++){
        data << this->cellPos[i][0] << ","      // x
             << this->cellPos[i][1] << ","      // y
             << this->level[i] << ","           // lvl
             << this->cellIdx[i][0] << ","      // Matrix ID x
             << this->cellIdx[i][1] << ","      // Matrix ID y
             << this->parNum[i] << ","          // Par num
             << this->srcNum[i] << ","          // Src num
             << this->leafCellMark[i] << ","    // is leaf
             << this->outsideCell[i] << ","     // is outside cell
             << std::pow(0.5, this->level[i])   // is size
             ;
        // par ID list
        if (parIn){
        data << ",";
            if(this->parNum[i] == 0){
                data << "-1";
            }else{
                data << this->parIDList[i][0];
                for (size_t j = 1; j < this->parIDList[i].size(); j++){
                    data << ";" << this->parIDList[i][j];
                }
            }
        }
        
        data << "\n";
    }
    data.close();

    // Save the matrix data
    data.open("matrix.dat");
    for (int k = 0; k <= this->max_level; k++){
        data << "Matrix at level " << k << ":\n";
        int _size = std::pow(2,k);
        // Loop the row (ID y)
        for (int j = _size - 1; j >=0 ; j--){
            // Loop the col (ID x)
            data << T_map2D[k][0][j];
            for (int i = 1; i < _size; i++){
                data << "," << T_map2D[k][i][j];
            }
            data << "\n";
        }
        data << "\n";
    }
    data.close();

    // Save the neighbor list data 
    if (_ngh){
        data.open("treeCellNgh.csv", std::ofstream::out | std::ofstream::app);
        data << "ngh1,ngh2,ngh3,ngh4\n";
        for (int i = 21; i < this->cellNum; i++){
            std::vector<int> List1,List2,List3,List4;
            this->intList(i,List1,List3);
            this->extList(i,List2,List4);
            
            // NGH 1: Touched cell
            for (size_t j = 0; j < List1.size(); j++){
                data << List1[j] << ";";
            }
            data << ",";

            // NGH 2: well separated cell
            for (size_t j = 0; j < List2.size(); j++){
                data << List2[j] << ";";
            }
            data << ",";
            
            // NGH 3: well separated near cell
            for (size_t j = 0; j < List3.size(); j++){
                data << List3[j] << ";";
            }
            data << ",";

            // NGH 4: well separated cell coarse
            for (size_t j = 0; j < List4.size(); j++){
                data << List4[j] << ";";
            }
            data << "\n";
        }
        data.close();
    }

    return;
}

void treeCell::saveTreeCell(std::string name){
    // internal variables
    std::ofstream data;
    std::string fileName = "treeCell";
    fileName += name;
    fileName += ".csv";

    bool _ngh = false, parIn = false;

    // Save the base cell data 
    data.open(fileName);
    data << "x,y,lvl,matIDx,matIDy,num,srcnum,isleaf,isoutside,size";
    if (parIn)  {data << ",parID\n";}
    else        {data << "\n";}

    for (int i = 0; i < this->cellNum; i++){
        data << this->cellPos[i][0] << ","      // x
             << this->cellPos[i][1] << ","      // y
             << this->level[i] << ","           // lvl
             << this->cellIdx[i][0] << ","      // Matrix ID x
             << this->cellIdx[i][1] << ","      // Matrix ID y
             << this->parNum[i] << ","          // Par num
             << this->srcNum[i] << ","          // Src num
             << this->leafCellMark[i] << ","    // is leaf
             << this->outsideCell[i] << ","     // is outside cell
             << std::pow(0.5, this->level[i])   // is size
             ;
        // par ID list
        if (parIn){
        data << ",";
            if(this->parNum[i] == 0){
                data << "-1";
            }else{
                data << this->parIDList[i][0];
                for (size_t j = 1; j < this->parIDList[i].size(); j++){
                    data << ";" << this->parIDList[i][j];
                }
            }
        }
        
        data << "\n";
    }
    data.close();

    // Save the matrix data
    data.open("matrix.dat");
    for (int k = 0; k <= this->max_level; k++){
        data << "Matrix at level " << k << ":\n";
        int _size = std::pow(2,k);
        // Loop the row (ID y)
        for (int j = _size - 1; j >=0 ; j--){
            // Loop the col (ID x)
            data << T_map2D[k][0][j];
            for (int i = 1; i < _size; i++){
                data << "," << T_map2D[k][i][j];
            }
            data << "\n";
        }
        data << "\n";
    }
    data.close();

    // Save the neighbor list data 
    if (_ngh){
        data.open("treeCellNgh.csv", std::ofstream::out | std::ofstream::app);
        data << "ngh1,ngh2,ngh3,ngh4\n";
        for (int i = 21; i < this->cellNum; i++){
            std::vector<int> List1,List2,List3,List4;
            this->intList(i,List1,List3);
            this->extList(i,List2,List4);
            
            // NGH 1: Touched cell
            for (size_t j = 0; j < List1.size(); j++){
                data << List1[j] << ";";
            }
            data << ",";

            // NGH 2: well separated cell
            for (size_t j = 0; j < List2.size(); j++){
                data << List2[j] << ";";
            }
            data << ",";
            
            // NGH 3: well separated near cell
            for (size_t j = 0; j < List3.size(); j++){
                data << List3[j] << ";";
            }
            data << ",";

            // NGH 4: well separated cell coarse
            for (size_t j = 0; j < List4.size(); j++){
                data << List4[j] << ";";
            }
            data << "\n";
        }
        data.close();
    }

    return;
}

// #pragma endregion