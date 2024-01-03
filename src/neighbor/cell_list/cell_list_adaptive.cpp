#include "cell_list.hpp"

// Set the levelUp flag at each cell (type 0 for clearing variable)
void CellList::setLevelUp(int basis_ID, int cell_ID, int opt){
    // Initialize the levelUp and targetLevel vector
    if (opt == 0){      // Set up - initialization
        // Initialize the levelUp and targetLevel vector
        this->levelUp.resize(this->num, std::vector<bool>(this->cellNum, false));
        // this->targetLevel.resize(this->num, std::vector<int>(this->cellNum, -1));
    }
    
    // Assign the levelUp and targetLevel vector
    else if (opt == 1){ // Assignment initialization
        // Check whether the current basis and cell pair is already on the list
        if (this->levelUp[basis_ID][cell_ID] != true){
            // [INFO]: The current ID pair is not on the list!
            // Assign the basis and cell ID pair to the list
            this->tempBasisCellID.emplace_back(std::vector<int>{basis_ID, cell_ID});
            // Assign the levelUp and targetLevel vector
            this->levelUp[basis_ID][cell_ID] = true;
            // this->targetLevel[basis_ID][cell_ID] = Pars::max_level;
        }
    }
}

// Perform the cell neighbor adaptation
void CellList::performAdaptive(Particle & parEval, int updateStage, int tgt_lvl, std::vector<double>& PARsize){
    // // Initialize the par Eval for first 
    // if (updateStage == 0){
    //     this->splitParIDFlag.clear();
    //     this->splitParIDFlag.resize(parEval.num, false);
    // }

    // This method will increase the neighbor cell level become 1 level below the current cell
    // Internal variable
    int BAS_ID, CHD_ID;                     // Evaluated cell properties
    int ngh_BAS_ID, ngh_CHD_ID, ngh_lvl;    // Evaluated neighbor cell properties
    std::vector<std::vector<int>> nghBasisCellID;    // Create a list basis cell ID pair for new target temp basis cell
    
    // Iterate all temporary basis cell ID pair
    for (size_t i = 0; i < this->tempBasisCellID.size(); i ++){
        // Take the basis and cell ID of current cell
        BAS_ID = this->tempBasisCellID[i][0];
        CHD_ID = this->tempBasisCellID[i][1];
        
        // Find the neighbor cell of the current cell ID
        std::vector<std::vector<int>> tempNghBasisCellID;
        this->neighborCell(tempNghBasisCellID, BAS_ID, CHD_ID);

        // Check all neighbor cell
        for (size_t j = 0 ; j < tempNghBasisCellID.size(); j++){
            // Take the basis and cell ID of current cell
            ngh_BAS_ID = tempNghBasisCellID[j][0];
            ngh_CHD_ID = tempNghBasisCellID[j][1];
            ngh_lvl = this->lvl[ngh_CHD_ID];
            
            // ... LEVEL 1 CHECK ...
            // Put into the list with 1 level lower if ngh_lvl below than tgt_lvl - 1
            // Example: _lvl = 4, ngh_lvl 2 get into the list, while ngh_lvl 3 not included
            if (ngh_lvl < tgt_lvl - 1){
                // ... LEVEL 2 CHECK ...
                // Check whether it is on the list or not
                if (this->levelUp[ngh_BAS_ID][ngh_CHD_ID] != true){
                    // [INFO]: The current ID pair is not on the list!
                    // Update the levelUp and targetLevel list
                    this->levelUp[ngh_BAS_ID][ngh_CHD_ID] = true;
                    // this->targetLevel[ngh_BAS_ID][ngh_CHD_ID] = tgt_lvl - 1;
                    // Assign the basis and cell ID pair to the dummy list
                    nghBasisCellID.emplace_back(std::vector<int>{ngh_BAS_ID, ngh_CHD_ID});
                }
            }
        }   // Done the iteration through the current neighbor
    }   // Done the iteration through the current cell level list
    
    // =======================================
    // ========= ADAPTIVE REFINEMENT =========
    // =======================================
    // Perform particle splitting
    this->adtParSplit(parEval, tgt_lvl, PARsize);

    // Perform cell splitting
    this->adtCellDivide(parEval);

    // Replace the old temporary basis cell ID pair data
    this->tempBasisCellID.clear();
    this->tempBasisCellID = nghBasisCellID;
    nghBasisCellID.clear();

    // Iteratively perform this method until the level evade to 1
    // The iteration is intended to upgrade the level of neighbor cell
    if (tgt_lvl > 1){
        this->performAdaptive(parEval, updateStage + 1, tgt_lvl - 1, PARsize);
    }
}

// Update adaptive particle
void CellList::adtParSplit(Particle & parBase, int tgt_lvl, std::vector<double>& PARsize){
    // ====== DIVIDE AND CONQUEROR PROCEDURE ====== //
    // Perform the devide and conqueror if the point inside the finer resolution block
    // int* DnC_x = new int[4]{1,-1,-1,1};   // Divide and conqueror sum for x
    // int* DnC_y = new int[4]{1,1,-1,-1};   // Divide and conqueror sum for y
    int initial_num = parBase.num;
    std::vector<int> _ID_list_init;
    
    // Do the DnC until the particle reach the target level
    std::vector<int> _ID_list;
    std::vector<int> temp_ID_list;
    
    // Put all particle ID to be divided into the list
    for (size_t i = 0; i < this->tempBasisCellID.size(); i ++){
        // Take the basis and cell ID of current cell
        int basis_ID = this->tempBasisCellID[i][0];
        int cell_ID = this->tempBasisCellID[i][1];
        for (auto _ID : this->CH_parID[basis_ID][cell_ID]){
            _ID_list.emplace_back(_ID);
        }
    }
    
    // Update the initial ID list
    _ID_list_init = _ID_list;

    // Iterate all level splitting, initially set into max_level
    for (int count = 0; count < tgt_lvl; count++){
        // Do DnC to all particles inside the list
        for (auto ID:_ID_list){
            // If the particle is inside the body (chi = 1)
            if (parBase.chi[ID] == 1){
                continue;
            }
            // Check whether the particle level still not reaching the target level
            int _diffLv = tgt_lvl - parBase.level[ID];
            
            // Do DIVIDE AND CONQUEROR if the target level bigger particle level
            // [NOTE]: R, u, and v (inside the parBase) is not updated in this section
            if (_diffLv > 0)
            {
                // Current particle cell ID
                int basis_ID = parBase.basis_label[ID];
                int cell_ID = parBase.cell_label[ID];

                // Initial particle position
                double _x = parBase.x[ID];
                double _y = parBase.y[ID];
                // double _z = parBase.z[ID];
                
                // Particle data after splitting
                double _s = parBase.s[ID] / 2.0;
                int _lvl = parBase.level[ID] + 1;

                // Update the particle variable size
                std::vector<int> split_IDs;      // splitted particle ID
                for (int _temp = 0; _temp < this->basis - 1; _temp++){
                    split_IDs.emplace_back(parBase.num);
                    parBase.num ++;
                    parBase.x.emplace_back(0.0e0);
                    parBase.y.emplace_back(0.0e0);
                    // parBase.z.emplace_back(0.0e0);
                    parBase.s.emplace_back(0.0e0);
                    parBase.neighbor.emplace_back(std::vector<int>());
                    parBase.level.emplace_back(0);
                    parBase.vorticity.emplace_back(0.0e0);
                    parBase.gz.emplace_back(0.0e0);
                    parBase.chi.emplace_back(0.0e0);
                    parBase.isActive.emplace_back(false);
                    parBase.basis_label.emplace_back(0);
                    parBase.cell_label.emplace_back(0);
                    PARsize.emplace_back(0);
                    
                    // this->splitParIDFlag.emplace_back(false);
                }

                // PERFORM PARTICLE SPLITTING
                // First determine the particle inside the cell
                // Internal variable for splitting
                int curr_ID;        // The target new ID for allocating splitted particle
                int _skip = 0;      // Leap the index after excluding the current particle ID
                double __x,__y;     // The position of new child
                
                // Cell center position and size
                double _xcell = this->mid_pos[basis_ID][cell_ID][0];
                double _ycell = this->mid_pos[basis_ID][cell_ID][1];
                double _scell = this->size / (std::pow(2,this->lvl[cell_ID]));
                
                // Evaluate each new splitted particle
                for (int j = 0; j < this->basis; j++){
                    // Calculate the position
                    __x = _x + 0.5 * _s * ((j&1)*2 - 1);     // Using binary number
                    __y = _y + 0.5 * _s * ((j&2) - 1);
                    // __z = _z + 0.5 * _s * ((j&4)/2 - 1);
                    
                    // Evaluate the current particle ID
                    curr_ID = split_IDs[j - _skip];
                    
                    // Update the ID: Check the new location inside the cell
                    if (_skip == 0){
                        if((__x > (_xcell - _scell/2.0)) && (__x < (_xcell + _scell/2.0))){
                            if((__y > (_ycell - _scell/2.0)) && (__y < (_ycell + _scell/2.0))){
                                _skip = 1;
                                curr_ID = ID;
                            }
                        }
                    }
                    
                    // Update the position, size, and level
                    parBase.x[curr_ID] = __x;
                    parBase.y[curr_ID] = __y;
                    // parBase.z[curr_ID] = __z;
                    parBase.s[curr_ID] = parBase.s[ID];
                    parBase.level[curr_ID] = _lvl;
                    parBase.neighbor[curr_ID] = parBase.neighbor[ID];
                    if (curr_ID != ID){parBase.neighbor[curr_ID].emplace_back(ID);}

                    // The following properties doesn't change in value (Temporary -> Will be change later)
                    parBase.vorticity[curr_ID] = parBase.vorticity[ID];
                    parBase.gz[curr_ID] = parBase.gz[ID];
                    parBase.chi[curr_ID] = parBase.chi[ID];
                    parBase.isActive[curr_ID] = parBase.isActive[ID];

                    // // ============== Evaluate and CELL ID ==============
                    // this->splitParIDFlag[curr_ID] = true;
                    PARsize[curr_ID] = _s;

                    // ============== Further DnC ==============
                    if (_diffLv > 1){temp_ID_list.emplace_back(curr_ID);}
                }
            }   // Done splitting on current particle
        }   // Done splitting on all particle inside the ID
        
        // Update the particle list
        _ID_list = temp_ID_list;
        temp_ID_list.clear();
        
        // Break the loop if no particle to be splitted
        if (_ID_list.empty()){break;}
    }

    // Determine the splitted particle cell_ID
    for (int _id = initial_num; _id < parBase.num; _id++){
        // Evaluate only the particle _id with true flag
        // bool flag = this->splitParIDFlag[_id];
        // // Release the ID from CH_parID
        // for (int _i = 0; _i < this->CH_parID[parBase.basis_label[_id]][parBase.cell_label[_id]].size(); _i ++){
        //     if (this->CH_parID[parBase.basis_label[_id]][parBase.cell_label[_id]][_i] == _id){
        //         this->CH_parID[parBase.basis_label[_id]][parBase.cell_label[_id]].erase(
        //             this->CH_parID[parBase.basis_label[_id]][parBase.cell_label[_id]].begin() + _i
        //         );
        //         break;
        //     }
        // }
        
        // ============== Evaluate the BASIS ID ==============
        // Initialize the internal variable
        int x_cell_pos = std::floor((parBase.x[_id] - this->pivot_pos[0])/this->size);
        int y_cell_pos = std::floor((parBase.y[_id] - this->pivot_pos[1])/this->size);
    
        // If the particle lies outside the x cell boundary
        if (x_cell_pos >= this->n_basis[0]){
            x_cell_pos = this->n_basis[0] - 1;
        }else if(x_cell_pos < 0){
            x_cell_pos = 0;
        }
        // If the particle lies outside the y cell boundary
        if (y_cell_pos >= this->n_basis[1]){
            y_cell_pos = this->n_basis[1] - 1;
        }else if (y_cell_pos < 0){
            y_cell_pos = 0;
        }

        // Calculate the basis cell ID from the obtained cell coordinate position
        parBase.basis_label[_id] = y_cell_pos*this->n_basis[0] + x_cell_pos;
        
        // ============== Evaluate the CELL ID ==============
        // Declare internal variable
        int BASIS_ID = parBase.basis_label[_id];   // Basis cell ID
        int CELL_ID = 0;                           // Current cell ID
        std::vector<int> CHD_IDs;                  // The temporary cell child ID list

        // Check the cell flag at current cell (flag for cell division)
        while (this->cell_flag[BASIS_ID][CELL_ID] != 0){
            // Update the temporary child ID list
            CHD_IDs = this->findChild(CELL_ID);
            
            int target_cellID = CHD_IDs[0];
            // Determine the new cell location, by check on each basis direction
            // In x direction
            if (parBase.x[_id] > this->mid_pos[BASIS_ID][CELL_ID][0]){
                target_cellID += 1;     // 1 = 2^0 (first direction label as 0)
            }
            // In y direction
            if (parBase.y[_id] > this->mid_pos[BASIS_ID][CELL_ID][1]){
                target_cellID += 2;     // 2 = 2^1 (second direction label as 1)
            }
            // // In z direction
            // if (parBase.z[_id] < this->mid_pos[BASIS_ID][CELL_ID][2]){
            //     target_cellID += 4;     // 4 = 2^2 (third direction label as 2)
            // }
            CELL_ID = target_cellID;

            // // Determine the new cell location
            // if (parBase.y[_id] < this->mid_pos[BASIS_ID][CELL_ID][1]){
            //     if (parBase.x[_id] < this->mid_pos[BASIS_ID][CELL_ID][0]){
            //         // Bottom-Left child cell -> Position 1
            //         CELL_ID = CHD_IDs[0];
            //     }
            //     else{
            //         // Bottom-Right child cell -> Position 2
            //         CELL_ID = CHD_IDs[1];
            //     }
            // }else{
            //     if (parBase.x[_id] < this->mid_pos[BASIS_ID][CELL_ID][0]){
            //         // Top-Left child cell -> Position 3
            //         CELL_ID = CHD_IDs[2];
            //     }
            //     else{
            //         // Top-Right child cell -> Position 4
            //         CELL_ID = CHD_IDs[3];
            //     }
            // }
        }

        // Update the particle cell ID and the CH_parID list
        parBase.cell_label[_id] = CELL_ID;
        this->CH_parID[BASIS_ID][CELL_ID].emplace_back(_id);

        // // Release the flag from splitParIDFlag list
        // this->splitParIDFlag[_id] = false;
    }
}

// Update adaptive cell
void CellList::adtCellDivide(Particle & parBase){
    // [!] Do the cell division for all Pair ID in this->tempBasisCellID list
    // Iterate through the this->tempBasisCellID list
    for (size_t i = 0; i < this->tempBasisCellID.size(); i ++){
        // The cell division is done until "cell level" == "min particle level inside the current cell"
        // What must be checked first: the division is done 
        
        // Declare internal variable
        std::vector<int> CELL_ID_list;    // List of all cell ID to be evaluated inside the current basis-cell pair
        std::vector<int> temp_CHD_list;   // List of temporary child cell ID
        int divide_flag;    // The flag to do cell division
        int BASIS_ID;       // The current basis cell ID            (fixed in this section)
        int CELL_ID;         // The curent evaluated cell ID         (will be iterated)
        int level;          // The current evaluated child cell ID  (will be iterated)

        // Condition at initial iteration
        BASIS_ID = this->tempBasisCellID[i][0];            
        CELL_ID_list.push_back(this->tempBasisCellID[i][1]);

        // Evaluate each child cell until adjusted with the target particle
        for (size_t j = 0; j < CELL_ID_list.size(); j++){
            // Determine the current cell ID and level
            CELL_ID = CELL_ID_list[j];
            level = this->lvl[CELL_ID];
            
            // The iteration is ended if the cell's ID reach the maximum level
            if (level == Pars::max_level){
                continue;
            }

            // The particle division flag (set as true at initial)
            divide_flag = true;

            // Update the divide_flag, (false) if the maximum particle level equal to the cell level
            // DEBUG:
            if (this->cell_flag[BASIS_ID][CELL_ID] != 0){
                std::cout << "[WARNING] the cell_flag is not 0" << std::endl;
            }
            
            // Check whether all particle level inside the cell is lower than the cell level (halt the cell division)
            for (auto ID:this->CH_parID[BASIS_ID][CELL_ID])
            {
                if (parBase.level[ID] <= level){
                    divide_flag = false;  // The cell will not be divided
                    break;
                }
            }

            // Do the cell division if the cell level still not fulfill the condition
            if (divide_flag){
                // Do particle division into new child cell
                this->divideCell(BASIS_ID, CELL_ID, parBase, this->CH_parID, this->cell_flag);
                
                // Calculate the child ID
                temp_CHD_list.clear();
                temp_CHD_list = this->findChild(CELL_ID);

                // Assign all new child ID into the CELL_ID_list
                for (auto _val : temp_CHD_list){
                    CELL_ID_list.push_back(_val);
                }
            }
        }
    }
}

// ************************************************************************************
// ====================================================================================

// Evaluating the particle neighbor IDs [DONE] [CLEAR]
void CellList::createCellList_new(Particle & particle){
    /*
    Note: 
    > Assign the particle ID into each basis cell list
    > Do the divide and conqueror to each basis cell
    > Assign the basis and cell ID into each particle
    */

    // Resize the basis and cell label ID of the particle
    particle.basis_label.resize(particle.num,0);
    particle.cell_label.resize(particle.num,0);

    // Clear the flag and ID list
    this->CH_new_parID.clear();
    this->cell_flag_new.clear();

    // Evaluate the basis cell ID of each particle
    for (int i = 0; i < particle.num; i ++){
        int x_cell_pos = std::floor((particle.x[i] - this->pivot_pos[0])/this->size);
        int y_cell_pos = std::floor((particle.y[i] - this->pivot_pos[1])/this->size);
        
        // If the particle lies right at the top x cell boundary
        if (x_cell_pos == this->n_basis[0]){
            x_cell_pos--;
        }
        // If the particle lies right at the top y cell boundary
        if (y_cell_pos == this->n_basis[1]){
            y_cell_pos--;
        }

        // Calculate the basis cell ID from the obtained cell coordinate position
        int cell = y_cell_pos*this->n_basis[0] + x_cell_pos;

        // Assign the cell ID to particle.basis_label and the correspoding particle ID to the cellList
        particle.basis_label[i] = cell;
        this->CH_new_parID[cell][0].push_back(i);
    }

    // Do a divide and conqueror for each Basis Cell
    for (int i = 0; i < this->num; i ++){
        // Basic parameter
        int divide_flag = true;

        // Condition at initial iteration
        int BASIS_ID = i;   // The current basis cell ID
        int PAR_ID;         // The parent cell of current evaluated cell ID
        int CHD_ID = 0;     // The current evaluated child cell ID
        int level;          // The current evaluated child cell ID
        
        // Add a new cell flag as for the basis cell
        this->cell_flag_new[BASIS_ID].push_back(0);

        // Evaluate each child cell until meet the target cell level
        while(true){
            // Determine the current cell level
            level = this->lvl[CHD_ID];
            PAR_ID = this->findParent(CHD_ID);
            
            // Close the loop until the iterated cell's ID reach the maximum level
            if (level == Pars::max_level){
                break;
            }

            // Add new child cell parameter element from the corresponding child ID in current Parent cell
            // Adjust the cell_flag and the middle point coordinate of each new child cell
            // int opX [4] = {-1, 1, -1, 1};
            // int opY [4] = {-1, -1, 1, 1};
            for (size_t j = 0; j < 4; j++){ // still not work for 3D
                this->cell_flag_new[BASIS_ID].push_back(1);
                this->CH_new_parID[BASIS_ID].emplace_back(std::vector<int>({}));
            }

            /*
            Note in deviding particle into new child cell:
            * The current cell is already divided
            * Divide the particle into new child only if the particle lies inside the current cell
            * The division flag is evaluated by compare all particle level toward current cell level
              * Action divide   : The largest particle have higher level than the cell
              * Action idle     : The largest particle have lower level than the cell
              * Action initial  : divided the particle into new child cell
            */
            // The particle division flag
            divide_flag = true;

            // Check whether the particle are inside the current cell or not, and not basis cell
            if (CHD_ID != 0 && 
                (
                (this->cell_flag_new[BASIS_ID][PAR_ID] == 0) || 
                (this->cell_flag_new[BASIS_ID][PAR_ID] == 1)
                )
            ){
                // The particle is located at parent cell
                divide_flag = false;
            }
            else{
                // The particle is inside the current cell
                // Check if all of the inner particle is in the level
                for (auto ID:CH_new_parID[BASIS_ID][CHD_ID])
                {
                    if (particle.level[ID] < level + 1){
                        divide_flag = false;  // The cell will not be divided
                        break;
                    }
                }
            }

            // Do the cell division if the cell level still not fulfill the condition
            if (divide_flag){
                // Do particle division into new child cell
                this->divideCell(BASIS_ID, CHD_ID, particle, this->CH_new_parID, this->cell_flag_new);
            }

            // Proceed to the next Child Cell
            CHD_ID ++;
        }
    }
}

std::vector<std::vector<int>> CellList::findNeighbor_new(const Particle & parEval, const Particle & parBase){
    std::vector<std::vector<int>> neighborList;
    // // The variable for particle redistribution 

    // // Resize the neighbor as the target particle parBase
    // neighborList.resize(parBase.num, std::vector<int>({}));

    // // TIREDDD !!!!!!!!
    
    // // The neighbor list (verlet list) will be evaluated from each cell
    // std::vector<int> _evalID;    // List of evaluated particle ID
    // int _skip;                   // Temporary variable: Position of skipped particle
    // for (int i = 0; i < particle.num; i++){
    //     _evalID.push_back(i);
    // }

    // // Evaluate the particle neighbor ID list of each particle per each cell
    // for (size_t i = 0; i < particle.num;){
    //     // Define the evaluation cell from evaluated particle
    //     int basis_ID = particle.basis_label[_evalID[0]];
    //     int cell_ID = particle.cell_label[_evalID[0]];

    //     // Initialize the variable and parameter for neighbor particle evaluation
    //     double x, y, s;
    //     double _x, _y, dist;
    //     int temp_basis_ID;
    //     int temp_cell_ID;

    //     // Find all neighboring cell
    //     std::vector<std::vector<int>> neighbor_PAR_CHD_ID;
    //     neighbor_PAR_CHD_ID = this->neighborCell(basis_ID, cell_ID);
        
    //     // Evaluate the neighbor of each particle inside the current cell ID
    //     for (auto parID:this->CH_parID[basis_ID][cell_ID]){
    //         // Note:
    //         //   * parID : The ID of evaluated particle
    //         //   * nghID : The ID of potential neighbor particle
            
    //         // Kick out the evaluated particle from the _evalID list and reduce the iteration loop
    //         for (int __i = 0; __i <= parID; __i++){
    //             if (_evalID[__i] == parID){
    //                 _evalID.erase(_evalID.begin() + __i);
    //                 i++;   // Reduce the for loop iteration by 1
    //                 break;
    //             }
    //         } 

    //         // Update the data of evaluated particle
    //         x = particle.x[parID];
    //         y = particle.y[parID];
    //         s = particle.s[parID];

    //         // Evaluate all particle inside the neighbooring cells
    //         for (auto PAIR:neighbor_PAR_CHD_ID){
    //             // Update the basis and cell ID the potential neighbor
    //             temp_basis_ID = PAIR[0];
    //             temp_cell_ID = PAIR[1];

    //             // Iterate all particle inside the current neighbor cell
    //             for (auto nghID:this->CH_parID[temp_basis_ID][temp_cell_ID]){
    //                 // Not include itself as the neighbor
    //                 if (parID == nghID){continue;} 

    //                 // Check the distance between particle
    //                 _x = particle.x[nghID];
    //                 _y = particle.y[nghID];
    //                 dist = std::sqrt(std::pow(x - _x,2) + std::pow(y - _y,2));
                    
    //                 // Verlet list evaluation
    //                 if (dist < Pars::r_buff * Pars::r_sup * s){
    //                     // Assign the neighbor ID if fullfill the criteria
    //                     particle.neighbor[parID].push_back(nghID);
    //                 }
    //             }
    //         }
    //     }
    // }
    
    return neighborList;
}

