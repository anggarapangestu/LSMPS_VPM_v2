#include "cell_list.hpp"

// +------------------------------------------------------------------------------------+
//  ============================ FOR CELL INITIALIZATION ===============================
// +------------------------------------------------------------------------------------+

// Cell List Initialization [DONE] [CLEAR] [UPDATED]
void CellList::initCellList(const Particle & particle){
    /*
    Cell List Initialization:
    1. Calculate the parent cell size, then calculate the domain extremes
    2. Calculate the number of basis cells
    3. Resize each CellList variable
    4. Calculate the domain center and cell pivot coordinate position
    5. Evaluate the middle point coordinate of each basis cell
    6. Create the look up table for each cell
    */
    
    // clock_t time_seg;
    // time_seg = clock();

    // [1] Parent cell size and basis calculation
    this->size = std::pow(2,Pars::max_level) * Pars::r_sup * Pars::r_buff * Pars::sigma;
    
    // calculation of domain extremes
    double *min_lim = new double[DIM];
    double *max_lim = new double[DIM];
    min_lim[0] = particle.x[0];
    min_lim[1] = particle.y[0];
    max_lim[0] = particle.x[0];
    max_lim[1] = particle.y[0];

    for (int i = 1; i < particle.num; i++){
        min_lim[0] = min_lim[0] < particle.x[i] ? min_lim[0] : particle.x[i];
        min_lim[1] = min_lim[1] < particle.y[i] ? min_lim[1] : particle.y[i];
        max_lim[0] = max_lim[0] > particle.x[i] ? max_lim[0] : particle.x[i];
        max_lim[1] = max_lim[1] > particle.y[i] ? max_lim[1] : particle.y[i];
    }

    // Current particle distribution domain size
    double x_dom_size = max_lim[0] - min_lim[0];
    double y_dom_size = max_lim[1] - min_lim[1];

    // time_seg = clock() - time_seg;
    // printf("TIME pos 1 <The Extremes>              [%f s]\n", (double)time_seg/CLOCKS_PER_SEC);
    // time_seg = clock();
    
    // [2] Calculation of cell number
    int _nx = 1 + 2 * std::ceil(0.5 * (x_dom_size / this->size - 1));
    int _ny = 1 + 2 * std::ceil(0.5 * (y_dom_size / this->size - 1));
    int _nz = 1;
    this->n_basis = {_nx, _ny, _nz};
    this->num = _nx * _ny * _nz;
    
    // [3] Resize each CellList variable
    this->totalCellNumber(Pars::max_level);   // Update the number of cell in the tree hierarchy
    // this->basis_neighbor.resize(this->num, std::vector<int>());
    this->mid_pos.resize(this->num, std::vector<std::vector<double>>(this->cellNum,std::vector<double>(DIM,0.0)));
    this->CH_parID.resize(this->num, std::vector<std::vector<int>>(this->cellNum,std::vector<int>()));
    this->cell_flag.resize(this->num, std::vector<int>(this->cellNum,1));       // All cell is set such that pointing the particle toward parent
    this->center_pos.resize(DIM,0);
    this->pivot_pos.resize(DIM,0);

    // time_seg = clock() - time_seg;
    // printf("TIME pos 2 <Vector Initialization>     [%f s]\n", (double)time_seg/CLOCKS_PER_SEC);
    // time_seg = clock();

    // [4] Calculation of domain center and cell pivot position
    for (size_t i = 0; i < DIM; i ++){
        // The domain center position
        this->center_pos[i] = ((max_lim[i] + min_lim[i])/2);

        // The cell pivot position
        this->pivot_pos[i] = (this->center_pos[i] - 0.5 * n_basis[i] * this->size);
    }

    // [5] calculate the basis middle point
    std::vector<int> _ind = std::vector<int>(DIM,0);
    for (int i = 0; i < this->num; i++){
        // Decomposed the cell hash ID
        for (int e = 0; e < DIM - 1; e++){
            _ind[e+1] += _ind[e] / n_basis[e];
            _ind[e] = _ind[e] % n_basis[e];
        }
        
        for (int e = 0; e < DIM; e++){
            double temp_pos = (_ind[e] + 0.5) * this->size + this->pivot_pos[e];
            this->mid_pos[i][0][e] = temp_pos;    // The root center position of 0
        }
        
        // Proceed to the next ID
        _ind[0] ++;
    }

    // time_seg = clock() - time_seg;
    // printf("TIME pos 3 <Mid Pos Calculation>       [%f s]\n", (double)time_seg/CLOCKS_PER_SEC);
    // time_seg = clock();

    // // [6] evaluate the basis cell neighbor
    // int temp;
    // int temp1;
    // int op [2] = {-1, 1};
    // int op1 [3] = {-1, 0, 1};
    
    // for (int i = 0; i < this->num; i++){
    //     // For 1 dimension
    //     for (int j = 0; j < 2; j++){
    //         temp = i + op[j];
    //         // Check whether the upper and lower x neighbor cell is out of the bound
    //         if (std::floor(temp/this->n_basis[0]) != std::floor(i/this->n_basis[0])){
    //             temp = -1;
    //         }
    //         basis_neighbor[i].push_back(temp);
    //     }
        
    //     // For 2 dimension
    //     if (DIM > 1){
    //         for (int j = 0; j < 2; j++){
    //             temp = i + op[j] * this->n_basis[0];
    //             // Check whether the upper and lower y neighbor cell is out of the bound
    //             if (std::floor(temp/(float)(this->n_basis[0]*this->n_basis[1])) != std::floor(i/(float)(this->n_basis[0]*this->n_basis[1]))){
    //                 temp = -1;
    //                 for (int k = 0; k < 3; k++){
    //                     basis_neighbor[i].push_back(temp);
    //                 }
    //                 continue;
    //             }

    //             // Evaluate the upper and lower y neighbor cell
    //             for (int k = 0; k < 3; k++){
    //                 temp1 = temp + op1[k];
    //                 if (std::floor(temp1/(float)this->n_basis[0]) != std::floor(temp/(float)this->n_basis[0])){
    //                     temp1 = -1;
    //                 }
    //                 basis_neighbor[i].push_back(temp1);
    //             }
    //         }
    //     }

    //     // For 3 dimension
    //     if (DIM > 2){
    //         // TO BE CONTINUED
    //     }
    // }

    // [6] Create the look up table for cell neighbor in each tree level
    // Update the size of each vector
    this->index.resize(this->cellNum,std::vector<int>(DIM,0));
    this->lvl.resize(this->cellNum);
    if (DIM == 2){        this->T_map2D.resize(Pars::max_level + 1);}
    else if (DIM == 3){   this->T_map3D.resize(Pars::max_level + 1);}

    // Initialize from the root cell (ID = 0)
    if (DIM == 2){        this->T_map2D[0] = std::vector<std::vector<int>>(1,std::vector<int>(1,0));}
    else if (DIM == 3){   this->T_map3D[0] = std::vector<std::vector<std::vector<int>>>(1,std::vector<std::vector<int>>(1,std::vector<int>(1,0)));}
    this->lvl[0] = 0;
    this->index[0] = std::vector<int>(DIM,0);
    
    // Calculate for the subsequent cell ID
    int* _temp_cellID_bound = new int [2]{0,0};     // Initialize by the root ID
    int _start,_fin;    // The bound cell ID each level
    int _ID = 1;        // The current evaluated cell ID pivot
    std::vector<int> c; // The index marking aid vector

    // Iterate through tree level
    for (int k = 1; k <= Pars::max_level; k ++){
        // Resize the t-map
        if (DIM == 2){
            T_map2D[k].resize(std::pow(2,k),
                              std::vector<int>(std::pow(2,k),0));
        }
        else if (DIM == 3){
            T_map3D[k].resize(std::pow(2,k),
                              std::vector<std::vector<int>>(std::pow(2,k),
                              std::vector<int>(std::pow(2,k),0)));
        }

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
                for (int i = 0; i < DIM - 1; i ++){
                    c[i+1] += c[i] / 2;
                    c[i] = c[i] % 2;
                }

                // The current cell ID by the given basis index
                currID = _ID + _ct;

                // Update the index
                for (int i = 0; i < DIM; i ++){
                    index[currID][i] = 2*index[_id][i] + c[i];
                    
                    // // Another way
                    // int _c,__c;
                    // __c = std::pow(2,i);
                    // _c = (_ct&__c)/__c;
                    // index[currID][i] = 2*index[_id][i] + _c;
                }
                // Update the level
                lvl[currID] = k;
                // Update the T-map
                if (DIM == 2){
                    T_map2D[k][index[currID][0]][index[currID][1]] = currID;
                }
                else if (DIM == 3){
                    T_map3D[k][index[currID][0]][index[currID][1]][index[currID][2]] = currID;
                }

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

    // time_seg = clock() - time_seg;
    // printf("TIME pos 4 <The lookup Table>          [%f s]\n", (double)time_seg/CLOCKS_PER_SEC);
    // time_seg = clock();
    
    // Release the heap variables
    delete[] min_lim, max_lim, _temp_cellID_bound;
    
    // Print the cell number
    printf("<+> Number of cell :                  %12d\n", this->num);
    /*
    std::cout << "/n====== DEBUGGER ======/n";
    std::cout << "All other parameter" << std::endl;
    std::cout << "basis = " << this->basis << std::endl;
    std::cout << "size = " << this->size << std::endl;
    std::cout << "nx = " << this->nx << std::endl;
    std::cout << "ny = " << this->ny << std::endl;
    
    std::cout << "Min x position " << min_lim[0] << std::endl;
    std::cout << "Min y position " << min_lim[1] << std::endl;
    std::cout << "Max x position " << max_lim[0] << std::endl;
    std::cout << "Max y position " << max_lim[1] << std::endl;

    std::cout << "Pivot x position " << this->pivot_pos[0] << std::endl;
    std::cout << "Pivot y position " << this->pivot_pos[1] << std::endl;
    std::cout << "Center x position " << this->center_pos[0] << std::endl;
    std::cout << "Center y position " << this->center_pos[1] << std::endl << std::endl;
    */
}

// Cell List Particle Evaluation [DONE] [CLEAR]
void CellList::createCellList(Particle & particle)
{
    /*
    Note: 
    > Assign the particle ID into each basis cell list
    > Do the divide and conqueror to each basis cell
    > Assign the basis and cell ID into each particle
    */

    // Resize the basis and cell label ID of the particle
    particle.basis_label.resize(particle.num,0);
    particle.cell_label.resize(particle.num,0);

    // Internal variables
    std::vector<int> x_cell_pos (DIM,0);

    // clock_t time_seg;
    // time_seg = clock();

    // Evaluate the basis cell ID of each particle and collect all particle at corresponding basis cell
    for (int i = 0; i < particle.num; i ++){
        x_cell_pos[0] = std::floor((particle.x[i] - this->pivot_pos[0])/this->size);
        x_cell_pos[1] = std::floor((particle.y[i] - this->pivot_pos[1])/this->size);
        // x_cell_pos[2] = std::floor((particle.z[i] - this->pivot_pos[2])/this->size);
        
        // If the particle pos exceeded the topmost cell boundary
        for (int e = 0; e < DIM; e++){
            if (x_cell_pos[e] == this->n_basis[e]){
                x_cell_pos[e]--;
            }
        }

        // Calculate the basis cell ID from the obtained cell coordinate position
        int cell = 0;
        for (int e = 0; e < DIM; e++){
            int _mult = 1;
            for (int a = 0; a < e; a++){
                _mult *= this->n_basis[a];
            }
            cell += x_cell_pos[e]*_mult;
        }

        // Assign the cell ID to particle.basis_label and the correspoding particle ID to the cellList
        particle.basis_label[i] = cell;
        this->CH_parID[cell][0].push_back(i);
    }

    // time_seg = clock() - time_seg;
    // printf("TIME pos 1 <Initialize particle inside>[%f s]\n", (double)time_seg/CLOCKS_PER_SEC);
    // time_seg = clock();

    // Do a divide and conqueror for each Basis Cell
    for (int i = 0; i < this->num; i ++){
        // Basic parameter
        int divide_flag;

        // Condition at initial iteration
        int BASIS_ID = i;           // The current basis cell ID
        // int PAR_ID;                 // The parent cell of current evaluated cell ID
        std::vector<int> CHILD_ID;  // The list of child ID
        int CURR_ID = 0;            // The current evaluated child cell ID (start from the root)
        int level;                  // The current evaluated cell level
        
        // Change the cell flag for the basis root cell
        this->cell_flag[BASIS_ID][0] = 0;       // The particle lies on the current cell

        // Evaluate each child cell until meet the target cell level
        while(true){
            // Determine the current cell level
            level = this->lvl[CURR_ID];
            // level = this->cellLevel(CURR_ID);
            // PAR_ID = this->findParent(CURR_ID);
            CHILD_ID = this->findChild(CURR_ID);
            
            // Close the loop until the iterated cell's ID reach the maximum level
            if (level == Pars::max_level){
                break;
            }

            // Add new child cell parameter element from the corresponding child ID in current Parent cell
            // [!] Adjust the cell_flag and the middle point coordinate of each new child cell
            std::vector<int> idx = std::vector<int>(DIM,0);   // Initialize the index counter
            for (int j = 0; j < this->basis; j++){               // does work for 3D
                // Update the counter
                for (int e = 0; e < DIM - 1; e ++){
                    idx[e+1] += idx[e] / 2;
                    idx[e] = idx[e] % 2;
                }

                // Calculate the child center position
                for (int e = 0; e < DIM; e ++){
                    this->mid_pos[BASIS_ID][CHILD_ID[j]][e] = this->mid_pos[BASIS_ID][CURR_ID][e] + 
                                                              (-1 + 2*idx[e]) * this->size / (std::pow(2,level) * 4.0);
                }

                // Proceed the counter
                idx[0] ++;
            }

            // // ======= The old method ========
            // int opX [4] = {-1, 1, -1, 1};
            // int opY [4] = {-1, -1, 1, 1};
            // for (size_t j = 0; j < 4; j++){ // still not work for 3D
            //     this->cell_flag[BASIS_ID].push_back(1);
            //     this->CH_parID[BASIS_ID].emplace_back(std::vector<int>({}));
            //     this->mid_pos[BASIS_ID].emplace_back(std::vector<double>({
            //         this->mid_pos[BASIS_ID][CURR_ID][0] + opX[j] * this->size / (std::pow(2,level) * 4.0),
            //         this->mid_pos[BASIS_ID][CURR_ID][1] + opY[j] * this->size / (std::pow(2,level) * 4.0)
            //     }));
            // }

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

            // Check whether the particle are inside the current cell (nor the parent neither child)
            if (this->cell_flag[BASIS_ID][CURR_ID] == 0){
                // The particle is inside the current cell
                
                // Check if all of the inner particle is in the level
                // Whether the cell have fine size particle, if the large particle level is in the same tree level, it wont divide the cell
                for (auto ID:CH_parID[BASIS_ID][CURR_ID])
                {
                    if (particle.level[ID] <= level){
                        divide_flag = false;  // The cell will not be divided
                        break;
                    }
                }
            }
            else
            {
                // The particle is located at parent or great parent cell
                divide_flag = false;
            }

            // Do the cell division if the cell level still not fulfill the condition
            if (divide_flag){
                // Do particle division into new child cell
                this->divideCell(BASIS_ID, CURR_ID, particle, this->CH_parID, this->cell_flag);
            }

            // Proceed to the next Child Cell
            CURR_ID ++;
        }
    }

    // time_seg = clock() - time_seg;
    // printf("TIME pos 2 <Divide and Conqueror>      [%f s]\n", (double)time_seg/CLOCKS_PER_SEC);
    // time_seg = clock();
}
