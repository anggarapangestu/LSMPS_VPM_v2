#include "cell_list.hpp"

// +------------------------------------------------------------------------------------+
//  ============================= CELL TOOLS & UTILITIES ===============================
// +------------------------------------------------------------------------------------+

// Devide and Conqueror [DONE] [CLEAR]
void CellList::divideCell(int BASIS_ID, int CURR_ID, Particle & particle, std::vector<std::vector<std::vector<int>>> & SRC_CH_parID, std::vector<std::vector<int>> & SRC_cell_flag){
    /* ----------------------------------------------- //
    NOTE:
     > BASIS_ID is the basis cell ID
     > CURR_ID is the current evaluated cell ID (will become Parent cell)
     > The particle child cell ID is updated
     > SRC_cell_flag    : The source of cell flag array
     > SRC_CH_parID     : The source of particle ID list
    // ----------------------------------------------- */

    // Find new child cell IDs
    std::vector<int> CHD_IDs = this->findChild(CURR_ID);
    
    // Change the cell_flag of each new child cell
    for (auto _val:CHD_IDs){
        SRC_cell_flag[BASIS_ID][_val] = 0;
    }
    
    // Close the cell_flag of the origin cell (current cell / new parent cell)
    SRC_cell_flag[BASIS_ID][CURR_ID] = -1;     // The particle IDs is lies on child

    // Group the particle from the old child cell into the new child cell
    for (size_t i = 0; i < SRC_CH_parID[BASIS_ID][CURR_ID].size(); i++)
    {
        // The particle ID from the old child cell
        int particle_ID = SRC_CH_parID[BASIS_ID][CURR_ID][i];
        
        // The target cell ID
        int target_cellID = CHD_IDs[0];

        // Determine the new cell location, by check on each basis direction
        // In x direction
        if (particle.x[particle_ID] > this->mid_pos[BASIS_ID][CURR_ID][0]){
            target_cellID += 1;     // 1 = 2^0 (first direction label as 0)
        }
        // In y direction
        if (particle.y[particle_ID] > this->mid_pos[BASIS_ID][CURR_ID][1]){
            target_cellID += 2;     // 2 = 2^1 (second direction label as 1)
        }
        // // In z direction
        // if (particle.z[particle_ID] < this->mid_pos[BASIS_ID][CURR_ID][2]){
        //     target_cellID += 4;     // 4 = 2^2 (third direction label as 2)
        // }

        SRC_CH_parID[BASIS_ID][target_cellID].emplace_back(particle_ID);
        particle.cell_label[particle_ID] = target_cellID;
        
        // ========= OLD CODE Method =========
        // // Determine the new cell location
        // if (particle.y[particle_ID] < this->mid_pos[BASIS_ID][CURR_ID][1]){
        //     if (particle.x[particle_ID] < this->mid_pos[BASIS_ID][CURR_ID][0]){
        //         // Bottom-Left child cell -> Position 1
        //         SRC_CH_parID[BASIS_ID][CHD_IDs[0]].emplace_back(particle_ID);
        //         particle.cell_label[particle_ID] = CHD_IDs[0];
        //     }
        //     else{
        //         // Bottom-Right child cell -> Position 2
        //         SRC_CH_parID[BASIS_ID][CHD_IDs[1]].emplace_back(particle_ID);
        //         particle.cell_label[particle_ID] = CHD_IDs[1];
        //     }
        // }else{
        //     if (particle.x[particle_ID] < this->mid_pos[BASIS_ID][CURR_ID][0]){
        //         // Top-Left child cell -> Position 3
        //         SRC_CH_parID[BASIS_ID][CHD_IDs[2]].emplace_back(particle_ID);
        //         particle.cell_label[particle_ID] = CHD_IDs[2];
        //     }
        //     else{
        //         // Top-Right child cell -> Position 4
        //         SRC_CH_parID[BASIS_ID][CHD_IDs[3]].emplace_back(particle_ID);
        //         particle.cell_label[particle_ID] = CHD_IDs[3];
        //     }
        // }
    }
    
    // Empty the particle IDs from the old child cell
    SRC_CH_parID[BASIS_ID][CURR_ID].clear();   // Clear all particle IDs inside this cell
}

// Finding the Child Cell IDs from a given Parent Cell ID [DONE] [CLEAR]
std::vector<int> CellList::findChild(int PAR_ID){
    std::vector<int> childList(this->basis,0);
    int temp_ID = PAR_ID * this->basis;
    for (int i = 0; i < this->basis; i++){
        temp_ID ++;
        childList[i] = temp_ID;
    }
    return childList;
}

// Finding the Parent Cell IDs from a given Child Cell ID [DONE] [CLEAR]
int CellList::findParent(int CHD_ID){
    int parID;
    parID = std::ceil(CHD_ID/std::pow(2,DIM)) - 1;
    return parID;
}

// Returning the number of cell in the hierarchy tree
void CellList::totalCellNumber(int k){
    int num;
    num = (std::pow(this->basis, k+1) - 1) / (this->basis - 1);
    this->cellNum = num;
}

/*
// // Returning the level of Cell IDs [DONE] [CLEAR]
// int CellList::cellLevel(int ID){
//     int level = 0;

//     // The level is 0 at child ID 0 (or the parent cell)
//     if (ID == 0){
//         return level;
//     }

//     // Cell level evaluation
//     int residu = std::ceil(ID/((float)this->basis));
//     int ct = 0;
//     while (true){
//         residu -= std::pow(this->basis,ct);
//         if (residu <= 0){
//             break;
//         }
//         ct++; 
//     }
    
//     // The cell level
//     level = ct + 1;
//     return level;
// }
*/

// +------------------------------------------------------------------------------------+
//  =============================== CELL SAVING TOOLS ==================================
// +------------------------------------------------------------------------------------+

// Saving the cell data
void CellList::saveCellData(){
    // -- accessing struct data, variables that are not commented are an input only
	int np = 0;
    std::vector<double> _x;
    std::vector<double> _y;
    std::vector<int> _basis;
    std::vector<int> _ID;
    std::vector<int> _cell_flag;
    std::vector<int> _level;
    std::vector<double> _size;
    std::vector<int> _parNum;

    for (size_t i = 0; i < this->mid_pos.size(); i++){
        for (size_t j = 0; j < this->mid_pos[i].size(); j++){
            if (this->cell_flag[i][j] == 0){
                _x.push_back(this->mid_pos[i][j][0]);
                _y.push_back(this->mid_pos[i][j][1]);
                _basis.push_back(i);
                _ID.push_back(j);
                _cell_flag.push_back(this->cell_flag[i][j]);
                _level.push_back(this->lvl[j]);
                _size.push_back(this->size / (std::pow(2,this->lvl[j])));
                _parNum.push_back(this->CH_parID[i][j].size());
                np++;
            }
        }
    }

	printf("Saving cell data ...\n");
	
    // Save common data (information of free particles):
    std::ofstream ofs;
	std::string name1;
	name1.append("output/cell_data.csv");
	ofs.open(name1.c_str());
	ofs << "" << "x" 
		<< "," << "y" 
		<< "," << "basis"
		<< "," << "ID" 
        << "," << "cell_flag" 
        << "," << "level" 
        << "," << "size" 
        << "," << "parNum" 
		<< "\n";

	for (int i = 0; i < np; i++)
	{
        ofs << "" << _x[i]
            << "," << _y[i]
            << "," << _basis[i]
            << "," << _ID[i]
            << "," << _cell_flag[i]
            << "," << _level[i]
            << "," << _size[i]
            << "," << _parNum[i]
            << "\n";
	}
	ofs.close();
}

// ------------------------------------------------------------------------------------
// ================================== FOR DEBUGGING ===================================
// ------------------------------------------------------------------------------------

// Displaying the basis cell neighbor ID
// void CellList::showBasisNgh(){
//     std::cout << "\nThe basis neighbor cell ID list:\n";
//     for (int i = 0; i < this->num; i ++){
//         std::cout << "[!] Basis : " << i << "\n  >> ";
//         for (int j = 0; j < this->basis_neighbor[i].size(); j ++){
//             std::cout << this->basis_neighbor[i][j] << ", ";
//         }
//         std::cout << std::endl;
//     }
// }

// Check double neighbor ID
void CellList::checkNGH(const Particle & par){
    std::cout << "\nNum of particle : " << par.num;
    std::cout << "\nNum of particle 2: " << par.neighbor.size();
    std::cout << "\nList of particle with double neighbor ID:\n";
    bool _break = false;
    for (size_t i = 0; i < par.neighbor.size(); i ++){ 
        for (size_t j = 0; j < par.neighbor[i].size(); j ++){
            for (size_t k = 0; k < par.neighbor[i].size(); k ++){
                if (j == k){continue;}
                if (par.neighbor[i][j] == par.neighbor[i][k]){
                    std::cout << " > Particle " << i << std::endl;
                    _break = true;
                    break;
                }
            }
            if (_break){
                _break = false;
                break;
            }
        }
    }
}

