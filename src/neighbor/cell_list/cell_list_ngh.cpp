#include "cell_list.hpp"

// +------------------------------------------------------------------------------------+
//  ============================== FOR NEIGHBOR SEARCH =================================
// +------------------------------------------------------------------------------------+

// Work for 2D and 3D
void CellList::neighborCell (std::vector<std::vector<int>>& ngh_PAIR_basis_cell_ID, int basis_ID, int cell_ID){
    /* CELL NEIGHBOR LIST CALCULATION PROCESS
       > [Procedure 1] Evaluate the basis-cell neighbor at current cell level
       > [Procedure 2] Proceed to the leaf cell from the given cell neighbor
    */
    
    // Initialize the internal variable
    ngh_PAIR_basis_cell_ID.clear();                     // Reset the pair list
    std::vector<std::vector<int>> temp_ngh_PAIR_list;   // Create temporary pair list
    int _lvl = this->lvl[cell_ID];                      // The current cell level
    int _size = std::pow(2,_lvl);                       // The transformation matrix size
    std::vector<int> curr_idx = this->index[cell_ID];   // The transformation mapping index of current cell

    // =====================================================
    // ==================== Procedure 1 ====================
    // =====================================================

    std::vector<int> idx(DIM,0);      // Trinary index for neighbor ID shift
    std::vector<int> mat_id(DIM,0);   // The candidate neighbor transformation matrix index
    std::vector<int> cg_id(DIM,0);    // The basis ID shift (-1: Prior Box, 0: Current Box, 1: Next Box)
    std::vector<int> cg_basis(DIM,0); // The basis ID loc. mark (-1: Lowest Bound, 0: Middle, 1: Highest Bound)
    int temp_ngh_basis;
    int temp_ngh_cell;

    // Determine the Basis Cell location from domain
    for (int e = 0; e < DIM; e ++){
        // Decompose the basis ID into one corresponding direction
        int __index = basis_ID;
        for (int a = 0; a < e; a ++){
            __index /= this->n_basis[a];
        }
        __index %= this->n_basis[e];

        // Change the cg_basis if the ID located at the domain boundary
        if (__index == 0){
            cg_basis[e] = -1;   // Lower-most
        }
        else if(__index == this->n_basis[e] - 1){
            cg_basis[e] = 1;    // Upper-most
        }
    }

    // Find the cell neighbor
    for (int i = 0; i < this->basisNGH; i ++){
        // Trinary index calculation
        for (int e = 0; e < DIM - 1; e ++){
            idx[e+1] += idx[e] / 3;
            idx[e] = idx[e] % 3;
        }

        // Update the candidate neighbor cell ID
        for (int e = 0; e < DIM; e ++){
            mat_id[e] = curr_idx[e] + idx[e] - 1 ;
            cg_id[e] = 0;

            // If the index is out of range change the cg_id
            if (mat_id[e] >= _size){    // Next-box
                cg_id[e] = 1;
                mat_id[e] -= _size;
            }
            else if (mat_id[e] <= -1){  // Previous-box
                cg_id[e] = -1;
                mat_id[e] += _size;
            }
        }
        // Find the corresponding basis ID using cg_id (Here the iteration might skip by NA basis cell)
        temp_ngh_basis = basis_ID;
        for (int e = 0; e < DIM; e++){
            /* Check for each direction
               > [Condition 1] If the corresponding neighbor basis cell is existed -> continue this iteration
               > [Condition 2] If not existed -> break this current iteration, continue to next evaluation
            */
            
            // The shifting is performed only when cg_id != 0
            if (cg_id[e] != 0){
                // Condition to break, the candidate is on the left next box but the current basis is the leftmost box, also the opposite
                if (cg_basis[e] == cg_id[e]){
                    goto NO_NEIGHBOR;   // The neighbor basis is not existed
                }
                else{
                    // The neighbor basis is existed, shift the basis ID
                    int _mult = 1;
                    for (int a = 0; a < e; a++){
                        _mult *= this->n_basis[a];
                    }
                    temp_ngh_basis += cg_id[e]*_mult;
                }
            }
        }

        // Find the corresponding cell ID using mat_id
        if (DIM == 2){
            temp_ngh_cell = T_map2D[_lvl][mat_id[0]][mat_id[1]];
        }
        else if (DIM == 3){
            // Continue next time
            // DONE at 18 April 2023 by Angga
            temp_ngh_cell = T_map3D[_lvl][mat_id[0]][mat_id[1]][mat_id[2]];
        }

        // Update the basis-cell ID pair 
        temp_ngh_PAIR_list.emplace_back(std::vector<int>({temp_ngh_basis, temp_ngh_cell}));
        
        // The goto statement if there is no basis cell neighbor
        NO_NEIGHBOR:
        idx[0]++;
    }

    // =====================================================
    // ==================== Procedure 2 ====================
    // =====================================================

    // Finalize the basis-cell pair ID of the neighboring cell
    for (size_t i = 0; i < temp_ngh_PAIR_list.size(); i ++){
        int __basisID = temp_ngh_PAIR_list[i][0];
        int __cellID = temp_ngh_PAIR_list[i][1];
        int __eval = this->cell_flag[__basisID][__cellID];
        
        // The particle ID lies on the current cell
        if(__eval == 0){
            // Check whether the current cell is similar to one in the final cell list [ngh_PAIR_basis_cell_ID]
            bool __skip = false;
            for (auto PAIR:ngh_PAIR_basis_cell_ID){
                if (__basisID == PAIR[0]){
                    if (__cellID == PAIR[1]){
                        __skip = true;
                    }
                }
            }
            if (__skip){
                continue;
            }

            // Add to the final ID list
            ngh_PAIR_basis_cell_ID.emplace_back(std::vector<int>{__basisID,__cellID});
        }
        // The particle ID stored at parent cell
        else if(__eval == 1){
            // Find parent ID
            int __tempcellID = this->findParent(__cellID);
            
            // Add to the temp ID list
            temp_ngh_PAIR_list.emplace_back(std::vector<int>{__basisID,__tempcellID});
        }
        // The particle ID stored at child cell
        else if(__eval == -1){
            // Find child ID
            std::vector<int> newChild = this->findChild(__cellID);

            // Add to the temp ID list
            for (auto _newChildID:newChild){
                temp_ngh_PAIR_list.emplace_back(std::vector<int>{__basisID,_newChildID});
            }
        }
    }

    return;
}

// Evaluating the particle neighbor IDs [DONE] [CLEAR]
void CellList::findNeighbor(Particle & particle)
{   
    /* PARTICLE NEIGHBOR EVALUATION PROCESS
       > Evaluate the neighbor of the particle from cell per cell (not the actual particle index)
       > There is an evalID list to support the identification of cell whether the current cell have been evaluated or not
       > There are 4 loops:
         * Loop 1: Iterate each particle 
                   # TARGET -> Get the current particle cell (CELL_1) to be evaluated
                   # TODO   -> Find the neighbor cell of CELL_1
         * Loop 2: Iterate each particle inside CELL_1
                   # TARGET -> Find the neighboring particle -> [Done on loop #4]
         * Loop 3: Iterate each neighbor cells of CELL_1
                   # TARGET -> Get all particle (PAR_LIST) inside the current neighbor cell
         * Loop 4: Iterate each particle inside PAR_LIST
                   # TARGET -> Get the neighbor particle of the particle from loop 2
                   # TODO   -> Calculate the current evaluated particle distance toward particle in loop 2
    */
    
    // Resize the neighbor
    particle.neighbor.clear();
    particle.neighbor.resize(particle.num, std::vector<int>({}));
    
    // The neighbor list (verlet list) will be evaluated from each cell
    std::vector<bool> _evalID;           // Flag to evaluate the particle neighbor or not
    _evalID.resize(particle.num, true);  // Initially all particle will be evaluated

    // Initialize the variable and parameter for neighbor particle evaluation
    double* x = new double[3];  // The current particle position
    double* _x = new double[3]; // The candidate neighbor particle position
    double s, dist;             // Current particle size and distance evaluation
    int temp_basis_ID;
    int temp_cell_ID;
    
    // Evaluate the particle neighbor ID list of each particle per each cell
    // #=============== LOOPS 1 ===============#
    for (int id = 0; id < particle.num; id++){
        // Only evaluate particle with flag:true
        if (_evalID[id] != true){continue;}

        // Define the current evaluated particle basis and cell ID
        int basis_ID = particle.basis_label[id];
        int cell_ID = particle.cell_label[id];

        // Find all neighboring cell
        std::vector<std::vector<int>> ngh_PAIR_basis_cell_ID;
        this->neighborCell(ngh_PAIR_basis_cell_ID,basis_ID, cell_ID);

        // std::cout << "At CELL : " << cell_ID << "\n";
        // for (auto PAIR:ngh_PAIR_basis_cell_ID){
        //     std::cout << " > Basis NGH : " << PAIR[0] << "; Cell NGH : " << PAIR[1] << "\n";
        // }

        // Evaluate the neighbor of each particle inside the CURRENT cell ID
        // #=============== LOOPS 2 ===============#
        for (auto parID:this->CH_parID[basis_ID][cell_ID]){
            // Note:
            //   * parID : The ID of evaluated particle
            //   * nghID : The ID of potential neighbor particle

            // Turn off the evaluation flag of current particle :: the particle have been evaluated
            _evalID[parID] = false;

            // Take the current particle position and size into temporary variable
            x[0] = particle.x[parID];
            x[1] = particle.y[parID];
            // x[2] = particle.z[parID];
            s = particle.s[parID];

            // Evaluate all particle inside the neighbooring cells
            // #=============== LOOPS 3 ===============#
            for (auto PAIR:ngh_PAIR_basis_cell_ID){
                // Update the basis and cell ID the potential neighbor
                temp_basis_ID = PAIR[0];
                temp_cell_ID = PAIR[1];

                // Iterate all particle inside the current neighbor cell
                // #=============== LOOPS 4 ===============#
                for (auto nghID:this->CH_parID[temp_basis_ID][temp_cell_ID]){
                    // Not include itself as the neighbor
                    if (parID == nghID){continue;} 

                    // Check the distance between particle
                    _x[0] = particle.x[nghID];
                    _x[1] = particle.y[nghID];
                    // _x[2] = particle.z[nghID];
                    
                    // Calculating distance
                    dist = 0;
                    for (int e = 0; e < DIM; e++){
                        dist += std::pow(x[e] - _x[e],2);
                    }
                    dist = std::sqrt(dist);
                    // std::cout << "The distance : " << dist << std::endl;
                    // Verlet list evaluation
                    if (dist < Pars::r_buff * Pars::r_sup * s){
                        // Assign the neighbor ID if fullfill the criteria
                        particle.neighbor[parID].push_back(nghID);
                    }
                }
            }
        }
    }
    
    delete [] x, _x;
}


// // ============================== OLD CODE =================================
// // Returning the neighbor Cell IDs [DONE] [CLEAR] only work for 2D
// void CellList::neighborCell(std::vector<std::vector<int>>& ngh_PAIR_basis_cell_ID,int basis_ID, int cell_ID){
//     // std::cout << "\n____Start of evaluating neighbor cell ID____\n";

//     // Declare the neighbor Basis - Cell ID pair
//     std::vector<std::vector<int>> temp_ngh_PAIR_basis_cell_ID;  // Temporary Basis - Cell ID pair
//     // std::vector<std::vector<int>> ngh_PAIR_basis_cell_ID;       // Final Basis - Cell ID pair
    
//     // Include the current PBasis - Cell ID pair
//     ngh_PAIR_basis_cell_ID.emplace_back(std::vector<int>({basis_ID, cell_ID}));
    
//     // Check the current cell level
//     int curr_level;                     // The level of current cell
//     curr_level = this->lvl[cell_ID];

//     // Procedure in checking the neighbor cell //
//     /*
//     1. The cell are grouped in 3 types: (1) corner, (2) side, (3) inside
//     2. Identifying the group of the cell is done by:
//        * Check the cell level
//        * Determine the position of the cell from "level" and "ID"
//        * (1) Corners are: initial level ID 1 + (1 <lv.1> + 4 <lv.2> + 16 <lv.3> + ...)*4 
//        *                  then multiplied by 1,2,3, and 4
//        * (2) Sides are: left (1,3), bottom (1,2), right (2,4), top (3,4)
//        *                use this combination to determine the ID position
//        * (3) Insides are other than mentioned above.
//     3. Compare the cell ID to the checker cell
//     4. Determine the neighbor based on each cell group
//     */
    
//     // =====================================================
//     // ==================== Procedure 1 ====================
//     // =====================================================
//     // std::cout << " <Procedure 1 - Creating temporary>\n";
//     // Find the temporary basis-cell pair ID of the neighboring cell
//     if (curr_level == 0){
//         // std::cout << "<!> BASIS CELL"<< std::endl;
//         for (auto _nghID:basis_neighbor[basis_ID]){
//             if (_nghID != -1){
//                 temp_ngh_PAIR_basis_cell_ID.emplace_back(std::vector<int>({_nghID, 0}));
//             }
//         }
//     }
//     else{
//         // ==================================================================
//         // ======================== Create Parameter ========================
//         // ==================================================================
//         // Calculate the basic parameter of each cell group

//         // ========== The CORNERS Cell IDs ==========
//         int crnID = 0;
//         for (int i = 0; i < curr_level - 1; i++){
//             crnID += std::pow(this->basis,i);
//         }
//         crnID = 1 + this->basis * crnID;
        
//         // ========== The SIDES Cell IDs ==========
//         // Create an ID list based on each position
//         std::vector<std::vector<int>> sdsID (this->basis, std::vector<int>{});
        
//         // Assign initial data
//         sdsID[0].push_back(crnID);      // Left
//         sdsID[1].push_back(crnID);      // Bottom
//         sdsID[2].push_back(2 * crnID);  // Right
//         sdsID[3].push_back(3 * crnID);  // Top
        
//         // The temporary evaluation parameter
//         std::vector<int> _hold;
//         int _posF = -1;
//         int _sum = 0;

//         for (int i = 0; i < std::pow(2,curr_level) - 1; i++){
//             // There are 2 procedure
//             //  1) Add new 0 level at the end of _hold
//             //  2) Check: the consecutive _hold element must be different, 
//             //            For a similar consecutive, delete the end then add the new end by 1
//             // The current procedure calculate a unique series algorithm
            
//             // Insert a new 0 level to the end of the array
//             _hold.push_back(0);
//             _posF++;
            
//             // Check if there is a consecutive same level element
//             for (size_t _pos = _posF; _pos > 0; _pos--){
//                 if (_hold[_posF] == _hold[_pos - 1]){
//                     _hold.pop_back();
//                     _posF--;
//                     _hold[_posF]++;
//                 }
//             }
            
//             // Create a sumEach iteration check whether the sequence is in order.
//             for (auto val:_hold){
//                 _sum += std::pow(this->basis,val);
//             }
            
//             // Assign each cell ID data to the list
//             sdsID[0].push_back(crnID + _sum * 2);      // Left
//             sdsID[1].push_back(crnID + _sum);          // Bottom
//             sdsID[2].push_back(2 * crnID + _sum * 2);  // Right
//             sdsID[3].push_back(3 * crnID + _sum);      // Top
//             _sum = 0;
//         }
        
//         // ========== The INSIDES Cell IDs ==========
//         // Create an ID list based on each position
//         std::vector<std::vector<int>> insID;

//         // Create an insides cell list from the Bottom Cell IDs [1] and Left Cell ID [0]
//         for (int i = 0; i < std::pow(2,curr_level); i++){
//             insID.emplace_back(std::vector<int>{});
//             for (int j = 0; j < std::pow(2,curr_level); j++){
//                 insID[i].push_back(sdsID[1][j] + sdsID[0][i] - crnID);
//             }
//         }

//         /*
//         The created parameter in this section:
//         * [int] crnID           : the ID of first corner cell
//         * [vec<vec<int>>] sdsID : the ID of sides cell (id_1: side location, id_2: cell index position)
//         * [vec<vec<int>>] insID : the ID of insides cell (id_1: y position, id_2: x position)
//         */
        
//         // ==================================================================
//         // ======================= Determine Neighbor =======================
//         // ==================================================================
//         // Determine the temporary neighbor cell ID

//         // Check the position of current cell (Corner, Side, or Inner)
//         if (cell_ID % crnID == 0){
//             // Current cell is located at corner
//             std::vector<int> __basis;
//             std::vector<int> __cell;
//             std::vector<int> __x;
//             std::vector<int> __y;
            
//             // Assign the neighbor cell IDs
//             switch (cell_ID / crnID - 1)
//             {
//             case 0: // Left Bottom
//                 __basis = {0, 0, 2, 3, 3};
//                 __cell = {sdsID[2][0], sdsID[2][1], crnID*4, sdsID[3][0], sdsID[3][1]};
//                 __x = {1, 1, 0};
//                 __y = {0, 1, 1};
//                 break;

//             case 1: // Right Bottom
//                 __basis = {1, 1, 4, 3, 3};
//                 __cell = {sdsID[0][0], sdsID[0][1], crnID*3, sdsID[3][sdsID[3].size()-1], sdsID[3][sdsID[3].size()-2]};
//                 __x = {(int)insID.size() - 2, (int)insID.size() - 2, (int)insID.size() - 1};
//                 __y = {0, 1, 1};
//                 break;

//             case 2: // Left Top
//                 __basis = {0, 0, 5, 6, 6};
//                 __cell = {sdsID[2][sdsID[2].size()-1], sdsID[2][sdsID[2].size()-2], crnID*2, sdsID[1][0], sdsID[1][1]};
//                 __x = {0, 1, 1};
//                 __y = {(int)insID.size()-2, (int)insID.size()-2, (int)insID.size()-1};
//                 break;

//             case 3: // Right Top
//                 __basis = {1, 1, 7, 6, 6};
//                 __cell = {sdsID[0][sdsID[0].size()-1], sdsID[0][sdsID[0].size()-2], crnID, sdsID[1][sdsID[1].size()-1], sdsID[1][sdsID[1].size()-2]};
//                 __x = {(int)insID.size()-2, (int)insID.size()-2, (int)insID.size()-1};
//                 __y = {(int)insID.size()-1, (int)insID.size()-2, (int)insID.size()-2};
//                 break;
//             }

//             // Assign the neighbor cell in the same basis cell
//             for (int _i = 0; _i < 3; _i++){
//                 temp_ngh_PAIR_basis_cell_ID.emplace_back(std::vector<int>({basis_ID, insID[__y[_i]][__x[_i]]}));
//             }
            
//             // Assign the neighbor cell at the nearest basis cell
//             for (int _i = 0; _i < 5; _i++){
//                 int __nghCellID = this->basis_neighbor[basis_ID][__basis[_i]];
//                 if (__nghCellID == -1){continue;}
//                 temp_ngh_PAIR_basis_cell_ID.emplace_back(std::vector<int>({__nghCellID, __cell[_i]}));
//             }

//         }else{
//             // Check if located at side
//             int __loc, __pos, _break_cond = 0;
//             bool atSide = false; // The flag for cell location
//             for (int i = 0; i < sdsID.size(); i++){
//                 for (int j = 0; j < sdsID[i].size(); j++){
//                     if (cell_ID == sdsID[i][j]){
//                         // Located at side cell
//                         atSide = true;
//                         __loc = i;  // i define the location (left, bottom, right, top)
//                         __pos = j;  // j define the position of the list
                        
//                         // Break 2 for loop
//                         _break_cond = 1;
//                         break;
//                     }
//                 }
//                 if (_break_cond == 1){break;}
//             }

//             // Assign the neighbor cell ID based on the side/inside location evaluation
//             if (atSide){
//                 // Current cell is located at side

//                 // Initialize the cell parameter
//                 int _posX, _posY, _nextID;
//                 std::vector<int> _opX;
//                 std::vector<int> _opY;

//                 switch (__loc)
//                 {
//                 case 0: // Located at left
//                     _posX = 0;
//                     _posY = __pos;
//                     _opX = {0, 1};
//                     _opY = {-1, 0, 1};
//                     _nextID = this->basis_neighbor[basis_ID][0];
//                     break;
//                 case 1: // Located at bottom
//                     _posX = __pos;
//                     _posY = 0;
//                     _opX = {-1, 0, 1};
//                     _opY = {0, 1};
//                     _nextID = this->basis_neighbor[basis_ID][3];
//                     break;
//                 case 2: // Located at right
//                     _posX = std::pow(2,curr_level) - 1;
//                     _posY = __pos;
//                     _opX = {-1, 0};
//                     _opY = {-1, 0, 1};
//                     _nextID = this->basis_neighbor[basis_ID][1];
//                     break;
//                 case 3: // Located at top
//                     _posX = __pos;
//                     _posY = std::pow(2,curr_level) - 1;
//                     _opX = {-1, 0, 1};
//                     _opY = {-1, 0};
//                     _nextID = this->basis_neighbor[basis_ID][6];
//                     break;
//                 }

//                 // Assign the cell inside the basis
//                 for (int i = 0; i < _opX.size(); i++){
//                     for (int j = 0; j < _opY.size(); j++){
//                         if (_opX[i] == 0 && _opY[j] == 0){continue;}
//                         temp_ngh_PAIR_basis_cell_ID.emplace_back(std::vector<int>{basis_ID,insID[_posY + _opY[j]][_posX + _opX[i]]});
//                     }
//                 }
                
//                 // Assign the cell outside the basis (only if existed)
//                 if (_nextID != -1){
//                     for (int i = -1; i < 2; i ++){
//                         temp_ngh_PAIR_basis_cell_ID.emplace_back(std::vector<int>{_nextID,sdsID[(__loc + 2)%this->basis][__pos + i]});
//                     }
//                 }

//             }else{
//                 // Current cell is located at inside

//                 // Calculate the current cell position
//                 int _posX, _posY, _break_cond = 0;
//                 for (int i = 0; i < std::pow(2,curr_level); i++){       // y position
//                     for (int j = 0; j < std::pow(2,curr_level); j++){   // x position
//                         // Create an insides cell list from the bottom cell IDs
//                         if (cell_ID == insID[i][j]){
//                             _posX = j;
//                             _posY = i;
//                             _break_cond = 1;
//                             break;
//                         }
//                     }
//                     if (_break_cond == 1){break;}
//                 }
                
//                 // Add the neighbor cell ID into the list
//                 int opVal[3] = {-1,0,1};
//                 for (int i = 0; i < 3; i++){
//                     for (int j = 0; j < 3; j++){
//                         if (opVal[i] == 0 && opVal[j] == 0){continue;}
//                         temp_ngh_PAIR_basis_cell_ID.emplace_back(std::vector<int>{basis_ID,insID[_posY + opVal[i]][_posX + opVal[j]]});
//                     }
//                 }
//             }
//         }
//     }

//     // =====================================================
//     // ==================== Procedure 2 ====================
//     // =====================================================

//     // Finalize the basis-cell pair ID of the neighboring cell
//     for (int i = 0; i < temp_ngh_PAIR_basis_cell_ID.size(); i ++){
//         int __basisID = temp_ngh_PAIR_basis_cell_ID[i][0];
//         int __cellID = temp_ngh_PAIR_basis_cell_ID[i][1];
//         int __eval = this->cell_flag[__basisID][__cellID];
        
//         // The particle ID lies on the current cell
//         if(__eval == 0){
//             // Check whether the current cell is similar to one in the final cell list [ngh_PAIR_basis_cell_ID]
//             bool __skip = false;
//             for (auto PAIR:ngh_PAIR_basis_cell_ID){
//                 if (__basisID == PAIR[0]){
//                     if (__cellID == PAIR[1]){
//                         __skip = true;
//                     }
//                 }
//             }
//             if (__skip){
//                 continue;
//             }

//             // Add to the final ID list
//             ngh_PAIR_basis_cell_ID.emplace_back(std::vector<int>{__basisID,__cellID});
//         }
//         // The particle ID stored at parent cell
//         else if(__eval == 1){
//             // Find parent ID
//             int __tempcellID = this->findParent(__cellID);
            
//             // Add to the temp ID list
//             temp_ngh_PAIR_basis_cell_ID.emplace_back(std::vector<int>{__basisID,__tempcellID});
//         }
//         // The particle ID stored at child cell
//         else if(__eval == -1){
//             // Find child ID
//             std::vector<int> newChild = this->findChild(__cellID);

//             // Add to the temp ID list
//             for (auto _newChildID:newChild){
//                 temp_ngh_PAIR_basis_cell_ID.emplace_back(std::vector<int>{__basisID,_newChildID});
//             }
//         }
//     }
    
//     // =====================================================
//     // ======================== End ========================
//     // =====================================================
    
//     return;
// }