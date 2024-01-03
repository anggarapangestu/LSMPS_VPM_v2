#include "treeCell.hpp"

#define ROOT_LEVEL 0


// =====================================================
// +--------------- Index Transformation ---------------+
// =====================================================
// #pragma region INDEX_TRANSFORMATION

/**
 *  @brief  Convert the cell index to cell ID.
 *
 *  @param  _index  Cell index to be transformed.
 *  @param  _level  The tree level to be evaluated.
 * 
 *  @return The cell ID at the given index and given tree level.
*/
int TreeCell::index2ID(const std::vector<int> &_index, const int &_level) const{
    // Internal variable
    int ID;
    int count = this->get_cellCount(_level);

    // Check whether the index is in range
    basis_loop(d){
        if (_index[d] < 0 || _index[d] >= count){
            ERROR_LOG << "Index [" << _index[0];
            for (int _i = 1; _i < DIM; _i++) std::cout << "," << _index[_i];
            std::cout << "] at level " << _level << " is out of range!\n";
            throw std::range_error("ERROR [TREE CELL]: Cannot convert index to cell ID!");
        }
    }

    // Calculate the starting ID at the given level.
    ID = this->get_startID(_level);

    // Flatten the ID
    int inc_ID, mul = 1;
    basis_loop(d){
        inc_ID = _index[d] * mul;   // Calculate the ID increment at this basis
        mul *= count;               // Update the flatten multiplier for next basis
        ID += inc_ID;               // Update the cell ID
    }
    
    return ID;
}

/**
 *  @brief  Convert cell ID to cell index.
 *
 *  @param  _ID      Cell ID to be transformed.
 *  @param  _level   The tree level to be evaluated.
 * 
 *  @return The index of current cell ID.
*/
std::vector<int> TreeCell::ID2index(const int &_ID, const int &_level) const{
    // Internal variable
    std::vector<int> index(DIM, 0);
    int count = this->get_cellCount(_level);
    int local_ID = _ID - this->get_startID(_level);

    // Check whether the ID lies in the evaluated level
    if ((local_ID < 0) || ((local_ID >= Pars::intPow(count, DIM)))){
        ERROR_LOG << "Cell " << _ID << " not lies in level " << _level << "!\n";
        throw std::runtime_error("ERROR [TREE CELL]: Invalid input!");
    }

    // ** Un-flatten transformation of node ID
    int div = 1;
    basis_loop(d){
        index[d] = (local_ID/div) % count;  // Calculate the current index
        div *= count;                       // Update the un-flatten divisor for next basis
    }

    return index;
}

/**
 *  @brief  Get the cell index that contained the given coordinate at the given tree level.
 *
 *  @param  _coord  A physical coordinate which cell index to be found.
 *  @param  _level  The tree level to be evaluated.
 * 
 *  @return The cell index that contains the given coordinate.
*/
std::vector<int> TreeCell::coord2index(const std::vector<double> &_coord, const int &_level) const{
    // Internal variable
    std::vector<int> index(DIM, 0);
    double cellSize = this->get_cellSize(_level);

    // Determine the index position
    basis_loop(d)
    index[d] = (_coord[d] - this->pivotCoord[d]) / cellSize;

    return index;
}

/**
 *  @brief  Get the cell ID that contained the given coordinate at the given tree level.
 *
 *  @param  _coord  A physical coordinate which cell ID to be found.
 *  @param  _level  The tree level to be evaluated.
 * 
 *  @return The cell ID that contains the given coordinate.
*/
int TreeCell::coord2ID(const std::vector<double> &_coord, const int &_level) const{
    int ID = this->index2ID(this->coord2index(_coord, _level), _level);

    return ID;
}

// #pragma endregion



// =====================================================
// +----------------- Tree Hierarchy ------------------+
// =====================================================
// #pragma region TREE_HIERARCHY

/**
 *  @brief  Get the parent ID of the given cell ID.
 *
 *  @param  _ID  Current cell ID for parent cell ID evaluation.
 * 
 *  @return The parent ID of current cell.
*/
int TreeCell::get_parentID(const int &_ID) const{
    // This function is not allowed to find parent of non existing cell
    // if (this->treeData.count(_ID) == 0){
    //     // The current cell is not existed
    //     ERROR_LOG << "Cell " << _ID << " is not existed!\n";
    //     throw std::range_error("ERROR [TREE CELL]: Taking parent ID is prohibited!");
    // }
    
    // Internal variable
    int parID;              // [OUTPUT] Parent cell ID
    int parLevel;           // Level of the parent cell
    int currLevel;          // Level of the current cell
    std::vector<int> currIndex(DIM); // Current cell index
    std::vector<int> parIndex(DIM);  // Parent cell index

    // Check whether the cell is existing in this tree or not
    if (this->treeData.count(_ID) == 0){
        // Cell is not existing in this tree
        currLevel = this->get_level(_ID);
        currIndex = this->ID2index(_ID, currLevel);
    }else{
        // Existed
        const Cell *currCell = this->treeData.at(_ID);
        currLevel = currCell->level;
        basis_loop(d) currIndex[d] = currCell->index[d];
    }
    
    // Calculating parent data
    parLevel = currLevel - 1;
    if (parLevel < 0) throw std::range_error("ERROR [TREE CELL]: No parent for root cell!");

    // Calculating parent ID
    basis_loop(d) parIndex[d] = currIndex[d] / 2;
    parID = this->index2ID(parIndex, parLevel);

    return parID;
}


// /**
//  *  @brief  Get the parent ID of the given cell ID.
//  *
//  *  @param  _ID  Current cell ID for parent cell ID evaluation.
//  *  @param  _level  Current cell level for parent cell ID evaluation.
//  * 
//  *  @return The parent ID of current cell.
// */
// int TreeCell::get_parentID(const int &_ID, const int &_level) const{
//     // Calculating parent index
//     int parID;                    // [OUTPUT] Parent cell ID
//     std::vector<int> parIndex;    // Parent cell index
//     int parLevel = _level - 1;    // Parent cell level
//     if (parLevel < 0) throw std::range_error("ERROR [TREE CELL]: No parent for root cell!");

//     // Calculating the current cell index
//     std::vector<int> currIndex = this->ID2index(_ID, _level);

//     // Calculating parent index
//     basis_loop(d) parIndex[d] = currIndex[d] / 2;
//     parID = this->index2ID(parIndex, parLevel);

//     return parID;
// }


/**
 *  @brief  Get the child ID list of the given cell ID.
 *
 *  @param  _ID  Current cell ID for child cell ID evaluation.
 * 
 *  @return The child ID list of current cell.
*/
std::vector<int> TreeCell::get_childID(const int &_ID) const{
    // This function is not allowed to find child of non existing cell
    if (this->treeData.count(_ID) == 0){
        // The current cell is not existed
        ERROR_LOG << "Cell " << _ID << " is not existed!\n";
        throw std::range_error("ERROR [TREE CELL]: Taking child ID is prohibited!");
    }

    // This function is not allowed to find child of leaf cell
    if (this->treeData.at(_ID)->isLeaf){
        // The current cell is a leaf cell
        ERROR_LOG << "Cell " << _ID << " is a leaf cell!\n";
        throw std::range_error("ERROR [TREE CELL]: Cell contains no childs!");
    }

    // Allias to the current cell
    const Cell *currCell = this->treeData.at(_ID);

    // Internal variable
    std::vector<int> chdIDList(this->chdNum,0); // [OUTPUT] child ID list
    std::vector<int> baseIndex(DIM);            // Cell index of the first child
    std::vector<int> chdIndex(DIM);             // Temporary child cell index
    int chdLevel = currCell->level + 1;         // Child cell tree level
    int chdID;                                  // Temporary variable for child ID
    
    // Calculate the base index
    basis_loop(d) baseIndex[d] = currCell->index[d] * 2;

    // Calculate all child ID
    for (int i = 0; i < this->chdNum; i++){
        // Set the current child index
        basis_loop(d) chdIndex[d] = baseIndex[d] + ((i>>d) & 1);
        
        // Retrieve the current child node ID
        chdID = this->index2ID(chdIndex, chdLevel);
        chdIDList[i] = chdID;
    }


    return chdIDList;
}


/** [DONE already altered]
 *  @brief  Get the cell level from the given cell ID.
 *
 *  @param  _ID  Current cell ID for level evaluation.
 * 
 *  @return The child ID list of current cell.
*/
int TreeCell::get_level(const long long int &_ID) const{
    /** The starting ID for base 4
     *   Level   Start ID    Cell ID list
     *     0   ->   0     -> {0}
     *     1   ->   1     -> {1,2,3,4}
     *     2   ->   5     -> {5,6,...,19,20}
     *     3   ->   21    -> {21,21,...,83,84}
     *     4   ->   85    -> {85,86,...,339,340}
     *     5   ->   341   -> {341,342,...,1365}
    */

    // Check for bad input
    if (_ID < 0){
        ERROR_LOG << "Problem at cell " << _ID << "!\n";
        throw std::invalid_argument("ERROR [TREE CELL]: Cannot find level of negative cell ID!");
    }

    // Direct level calculation (actually floor function)
    int level = std::log(1 + _ID*(this->chdNum - 1)) / std::log(this->chdNum) + (10 * __DBL_EPSILON__);

    // Check for bad truncation at current cell
    /** Truncation error problem example
     *    >> In 2D with base of 4, we about to find the level of cell with ID 340 (that actually in level 4)
     *    >> By calculation we get level 5 instead of level 4 (division truncation error)
     *    >> The starting ID of level 5 is 341 that is bigger than the current ID 340
     *    >> This is a violation <!> (At the same level the ID must be bigger than the starting ID)
     *    >> This violation will be called as level shifting
     *    >> We need further treatment (at the line below)
    */
    
    // The level can only shifts in one level (shifted upward)
    const long long int startIDcurr = this->get_startID(level);
    if (startIDcurr > _ID){
        level--;
    }

    // // The level can only shifts in one level (shifted downward)
    // const long long int startIDnext = this->get_startID(level+1);
    // if (_ID >= startIDnext){
    //     level++;
    // }
    
    return level;
}

// #pragma endregion



// =====================================================
// +-------------------- Utilities --------------------+
// =====================================================
// #pragma region UTILITIES
    
// Get the size of cell at given level
double TreeCell::get_cellSize(const int &_level) const{
    // The cell size is the division from its parent cell size
    double cellSize = this->rootLength / std::pow(2.0, _level);
    return cellSize;
}

// Get the cell ID starting at the given level
long long int TreeCell::get_startID(const int &_level) const{
    /** Illustration of using 
     *    
     *   idx_y  ___________________  idx_y  ___________________   idx_y  ___________________ 
     *         |                   |       |         |         |      3 | 17 | 18 | 19 | 20 |
     *         |                   |     1 |  ID 3   |  ID 4   |        |____|____|____|____|
     *         |                   |       |         |         |      2 | 13 | 14 | 15 | 16 |
     *       0 |       ID 0        |       |_________|_________|        |____|____|____|____|
     *         |                   |       |         |         |      1 | 9  | 10 | 11 | 12 |
     *         |                   |     0 |  ID 1   |  ID 2   |        |____|____|____|____|
     *         |                   |       |         |         |      0 | 5  | 6  | 7  | 8  |
     *         |___________________|       |_________|_________|        |____|____|____|____|
     *   idx_x          0                       0         1               0    1    2    3 
     *              At level 0                  At level 1                   At level 2
     *  The series of starting ID (*for 2D)
     *      ID = 0, 1, 5, 21, 85, ...
     *  The series of starting ID increment (the ratio of series is child count)
     *     inc = 1, 4, 16, 64, ...
    */

    // Associate with the geometric series (Sn = a (r^n - 1) / (r - 1))
    long long int startID = (static_cast<long long int> (Pars::intPow(this->chdNum, _level)) - 1) / (static_cast<long long int> (this->chdNum) - 1);

    return startID;
}

// Get the number of cell at current level (the number at one direction only)
int TreeCell::get_cellCount(const int &_level) const{
    // The cell count in one dimension at the given level
    int cellCount = Pars::intPow(2, _level);
    return cellCount;
}

// #pragma endregion


// =====================================================
// +------------------- Saving data -------------------+
// =====================================================
// #pragma region SAVING_DATA

/**
 *  @brief  Write all cell in tree.
 *         
 *  @param  tree The tree data to be saved.
 *  @param  name The identification name for the saved file.
*/
void TreeCell::saveTree(const TreeCell &tree, std::string name) const{
    // Initialization of saving data
    std::ofstream writer;
    name = "output/treeCell_all_" + name + ".csv";
    writer.open(name.c_str());
    
    // Write data table header
    if (DIM == 2){
        writer << "ID,xId,yId,x,y,lvl,size,isLeaf,isActive,parNum,srcNum\n";
    }else if (DIM == 3){
        writer << "ID,xId,yId,zId,x,y,z,lvl,size,isLeaf,isActive,parNum,srcNum\n";
    }

    // Write data value
    for(const auto &[_ID, _cell] : tree.treeData){
        writer << _cell->ID;
        basis_loop(d)
        writer << "," << _cell->index[d];
        basis_loop(d)
        writer << "," << _cell->centerPos[d];
        writer << "," << _cell->level
               << "," << _cell->length
               << "," << _cell->isLeaf
               << "," << _cell->isActive
               << "," << _cell->parNum
               << "," << _cell->srcParNum
               << "\n";
    }
    
    writer.close();
    return;
}

/**
 *  @brief  Write all leaf cell in tree.
 *         
 *  @param  tree The tree data to be saved.
 *  @param  name The identification name for the saved file.
*/
void TreeCell::saveLeafTree(const TreeCell &tree, std::string name) const{
    // Initialization of saving data
    std::ofstream writer;
    name = "output/treeCell_leaf_" + name + ".csv";
    writer.open(name.c_str());
    
    // Write data table header
    if (DIM == 2){
        writer << "ID,xId,yId,x,y,lvl,size,isActive,parNum,srcNum\n";
    }else if (DIM == 3){
        writer << "ID,xId,yId,zId,x,y,z,lvl,size,isActive,parNum,srcNum\n";
    }

    // Write data value
    for(const auto &_ID : tree.leafList){
        const Cell *_cell = tree.treeData.at(_ID);
        writer << _cell->ID;
        basis_loop(d)
        writer << "," << _cell->index[d];
        basis_loop(d)
        writer << "," << _cell->centerPos[d];
        writer << "," << _cell->level
               << "," << _cell->length
               << "," << _cell->isActive
               << "," << _cell->parNum
               << "," << _cell->srcParNum
               << "\n";
    }
    
    writer.close();
    return;
}

/**
 *  @brief  Write the data of selected ID.
 *         
 *  @param  tree The tree data to be saved.
 *  @param  name The identification name for the saved file.
 *  @param  IDList The list of cell ID to be saved
*/
void TreeCell::saveSelTree(const TreeCell &tree, std::string name, std::vector<int> &IDList) const{
    // Initialization of saving data
    std::ofstream writer;
    name = "output/treeCell_Selected_" + name + ".csv";
    writer.open(name.c_str());
    
    // Write data table header
    if (DIM == 2){
        writer << "ID,xId,yId,x,y,lvl,size,isLeaf,isActive,parNum,srcNum\n";
    }else if (DIM == 3){
        writer << "ID,xId,yId,zId,x,y,z,lvl,size,isLeaf,isActive,parNum,srcNum\n";
    }

    // Write data value
    for(const auto &_ID : IDList){
        const Cell *_cell = tree.treeData.at(_ID);
        writer << _cell->ID;
        basis_loop(d)
        writer << "," << _cell->index[d];
        basis_loop(d)
        writer << "," << _cell->centerPos[d];
        writer << "," << _cell->level
               << "," << _cell->length
               << "," << _cell->isLeaf
               << "," << _cell->isActive
               << "," << _cell->parNum
               << "," << _cell->srcParNum
               << "\n";
    }
    
    writer.close();
    return;
}

// #pragma endregion


// =====================================================
// +------------- Neighbor Cell Evaluation ------------+
// =====================================================
// #pragma region NEIGHBOR_EVALUATION

/**
 *  @brief  Find the cell neighbor at the same level.
 *         
 *  @param  currCell The current cell ID for.
 *  @param  nghIDList [OUTPUT] The list of neighbor cell ID
*/
void TreeCell::findNghLvl(const Cell *currCell, std::vector<int> &nghIDList) const{
    // Evaluate all node including itself
    // The neighbors are all IDs in the local neighbor index map *

    // Procedure:
    // [1] Iterate through all potential neighbor 1D->3ngh; 2D->9ngh; 3D->27ngh
    /*/     The calculation includes 
            > De-flatten algorithm
            > Translation between flat and matrix indexing
    /*/
    // Prepare the neighbor ID list container
    nghIDList.clear();

    // Intermediate variable
    int nghNum = this->nghNum;   // Number of neighbor candidate
    int nghIDtar;       // ID candidate of target neighbor
    int div[DIM];       // Divisor value to de-flatten neigbor local-index [1,3,9]
    int nGrid[DIM];     // Grid node count at the current level [nx, ny, nz]
    int transMul[DIM];  // ID translation multiplyer at each basis [1, nx, nx*ny]
    int transDir[DIM];  // The index translation direction at each basis
    int count = this->get_cellCount(currCell->level);   // The number of cell in the current level (in one dimension)
    
    // Initialize the intermediate variables
    basis_loop(d){
        div[d] = Pars::intPow(3, d);
        nGrid[d] = count;
        transMul[d] = 1;
        for (int i = 0; i < d; i++){
            transMul[d] *= nGrid[i];
        }
    }
        
    // Check every posibilies of neighbor candidate
    for (int i = 0; i < nghNum; i++){
        // SUB-PROCESS 1 [Calculation 1]
        // Calculate the translation direction and check range
        bool inRange = true;
        // **Translation direction
        basis_loop(d){
            transDir[d] = ((i / div[d]) % 3) - 1;
            if ((currCell->index[d] == 0 && transDir[d] == -1) ||
                (currCell->index[d] == nGrid[d]-1 && transDir[d] == 1)){
                // No node candidate in this target location
                inRange = false;
                break;
            }
        }
        // **Check in range
        if (!inRange) continue;

        
        // SUB-PROCESS 2 [Calculation 2]
        // Calculate the candidate ID of the neighbor
        // **Initialize the neighbor candidate ID
        nghIDtar = currCell->ID;    // ID pivots on the current node ID
        // **Translate the neighbor ID from pivot
        basis_loop(d) nghIDtar += transDir[d] * transMul[d];
        
        // SUB-PROCESS 3
        // Put the neighbor ID into the container
        nghIDList.push_back(nghIDtar);
    }
    
    return;
}

/**
 *  @brief  Find the internal neighbor list. It consited of 2 near neighbor groups.
 *  List 1:= collects all leaf cell that adjacent with the current cell. 
 *  List 3:= collects all sibling cell from list 1 but other than the cell in list 1.
 *         
 *  @param  currID Current cell ID to be evaluated.
 *  @param  ID_List_1 [OUTPUT] The list 1 neighbor.
 *  @param  ID_List_3 [OUTPUT] The list 3 neighbor.
*/
void TreeCell::intList(const int &currID, std::vector<int>& ID_list_1, std::vector<int>& ID_list_3) const{
    // Exception 1
    if (this->treeData.count(currID) == 0){
        ERROR_LOG << "The cell " << currID << " is not existed!\n";
        throw std::invalid_argument("ERROR [TREE CELL]: No near cell ID list for this cell!");
    }
    
    // Exception 2
    if (this->treeData.at(currID)->isLeaf == false){
        ERROR_LOG << "The cell " << currID << " is not leaf!\n";
        throw std::invalid_argument("ERROR [TREE CELL]: Internal list calculation is prohibited!");
    }

    // Create an alias to the current cell data (Note that the current cell must be a leaf cell)
    const Cell *currCell = this->treeData.at(currID);

    // Reset the neighbor list
    ID_list_1.clear();
    ID_list_3.clear();
    
    // Procedure:
    // 1. List the neighbor cell (at the same level) that touched with current cell
    // 2. Compare each cell box to the list criteria
    //    > List 1 : Leaf cell, touched with the current cell
    //    > List 3 : Not necessary a leaf cell, included by neighbor cell but other than list 1

    // Internal variable
    double currSize = currCell->length;     // The length of current cell

    // The temporary neighbor list
    std::vector<int> tempNghList;   // The temporary neighbor ID list as a queue list
    std::unordered_map<int, bool> tempNghFlag;   // The temporary neighbor ID list as a queue list

    // PROCEDURE 1! : Find all neighbor cell ID touched with current cell
    // ************
    this->findNghLvl(currCell, tempNghList);

    for (const auto &nghID : tempNghList){
        tempNghFlag.insert({nghID, true});
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

    /** Illustration
     *  _________________
     * |     |     |__|__|  The queue at first will only contain the adjacent neighboring cell
     * |_____|_____|__|__|   * So than at each iteration if leaf, put into list 1, but if not
     * |     |  +  |     |      Evaluate the child.
     * |_____|_____|_____|   * At child evaluation, the adjacent cell will put into the queue, 
     * |           |     |      while the other put directly put into the list 3.
     * |           |_____|   * So that after the next iteration the queue list only contains the
     * |           |            adjacent cell
     * |___________|         * Actually we need to check whether the node is existed or not
     *                                      
     * 
    */

    // Evaluate all neighbor candidate
    for (size_t i = 0; i < tempNghList.size(); i++){
        // Create an alias to the current neighbor cell
        int nghID = tempNghList[i];

        // First Check whether the cell is existed
        if (this->treeData.count(nghID) == 0){
            // The cell is not existed, take the parent into the evaluation
            int parID = this->get_parentID(nghID);

            // Put the parent ID into the queue (only if it not existed in the queue)
            if (tempNghFlag.count(parID) == 0){
                tempNghList.push_back(parID);
                tempNghFlag.insert({parID, true});
            }
            
            // Proceed to the next iteration
            continue;
        }
        
        // Cell is existed, do a direct evaluation
        const Cell *nghCell = this->treeData.at(nghID);
        
        // Skip the candidate neighbor if not active (not having particle source)
        if (nghCell->isActive == false) continue;

        // Check whether a leaf cell
        if (nghCell->isLeaf == true){
            // The current cell is a leaf cell, put into the list 1
            ID_list_1.push_back(nghID);
        }
        else{
            // The current cell is not a leaf cell, evaluate its child
            std::vector<int> _chdIDs = this->get_childID(nghID);

            // A length bound for adjacent cell (the cell with distance below this length is catagorized as adjacent cell)
            double adjacent_dist = ((currSize + this->treeData.at(_chdIDs[0])->length) / 2.0) + (10 * __DBL_EPSILON__);

            // [Take the child cell]
            // (1) Check whether is touching current cell
            // (2) Put into corresponding temporary list (1 or 3)

            // Check all childs cell
            for (const auto &_chdID : _chdIDs){
                // Create alias to the current child cell
                const Cell *chdCell = this->treeData.at(_chdID);

                // Check whether the current cell is active or not
                if (chdCell->isActive == false) continue;
                
                // Check whether the current child is adjacent to the current cell
                bool isAdjacent = true;
                basis_loop(d){
                    // Check distance in each basis direction
                    double dist = std::abs(chdCell->centerPos[d] - currCell->centerPos[d]);
                    // Check whether is adjacent or not
                    if (dist > adjacent_dist){
                        isAdjacent = false;
                        break;
                    }
                }
                
                // Do Final conditional
                if (isAdjacent){
                    // If adjacent (touched) move into queue for further check
                    tempNghList.push_back(_chdID);
                }else{
                    // If not adjacent (seperated) it will be the list 3
                    ID_list_3.push_back(_chdID);
                }
            }
        }
    }
    
    return;
}

/**
 *  @brief  Find the external neighbor list. It consited of 2 far neighbor groups.
 *  List 2:= collects all cells in the same level with the current cell that are seperated 
 *  by one block apart. The cell in list 2 are the child of parent neighbor list in the same level
 *  of its parent.
 *  List 4:= collects all leaf cell with lower level that are seperated by one block apart
 *  toward the current cell.
 *         
 *  @param  currID Current cell ID to be evaluated.
 *  @param  ID_List_2 [OUTPUT] The list 2 neighbor.
 *  @param  ID_List_4 [OUTPUT] The list 4 neighbor.
*/
void TreeCell::extList(const int &currID, std::vector<int>& ID_list_2, std::vector<int>& ID_list_4) const{
    // Exception handling
    if (this->treeData.count(currID) == 0){
        ERROR_LOG << "The cell " << currID << " is not existed!\n";
        throw std::invalid_argument("ERROR [TREE CELL]: Cannot find near cell ID list for this cell!");
    }

    // Reset the neighbor list
    ID_list_2.clear();
    ID_list_4.clear();
    
    // Procedure:
    // 1. List the neighbor cell (at the parent same level) that touched with the parent of current cell
    // 2. Compare each cell box to the list criteria
    //    > List 2 : Same level to current cell, well separate to current cell (separated by 1 cell apart) (not necessary a leaf cell)
    //    > List 4 : Lower level to current, leaf cell, well separated to current cell

    // Internal variable
    // Current cell (Note that the current cell is not necessary a leaf cell)
    const Cell *currCell = this->treeData.at(currID);   // Create an alias to the current cell
    double currSize = currCell->length;                 // The length of current cell
    
    // Parent cell
    int parID = this->get_parentID(currID);             // The current cell parent ID
    const Cell *parCell = this->treeData.at(parID);     // Create an alias to the current cell

    std::vector<int> ID_list_2_temp;    // The temporary neighbor ID list 2
    std::vector<int> ID_list_4_temp;    // The temporary neighbor ID list 4 (a queue list)
    std::unordered_map<int, bool> tempList4Flag;   // The temporary neighbor ID list as a queue list

    // // The temporary neighbor list
    // std::vector<int> tempNghList;   // The temporary neighbor ID list as a queue list
    // std::unordered_map<int, bool> tempNghFlag;   // The temporary neighbor ID list as a queue list

    // Procedure 1! : List all parent cell neighbor ID
    // ************
    this->findNghLvl(parCell, ID_list_2_temp);
    
    
    // Procedure 2! : Evaluate the list cretiria
    // ************
    // - Leaf cell [Y] - Put as candidate list 4
    //   - Well separated [Y] - Take as the final list 4
    //   - Well separated [N] - Throw the ID
    // - Leaf cell [N] - Take the child -> check well separated
    //   - Well separated [Y] - Take as the final list 2
    //   - Well separated [N] - Throw the ID
    // - Have no source particle - Take the parent ID into temporary list 4

    
    // **Evaluate in the same level (temporary list 2)
    for (size_t i = 0; i < ID_list_2_temp.size(); i++){
        // *Note: Here we are checking the neighbor of parent cell

        // Create an alias to the current neighbor cell
        int nghID = ID_list_2_temp[i];      // Neighbor of the parent cell
        
        // Skip the parent cell itself
        if (nghID == parID) continue;

        // *[CHECK Existed] Check whether the cell is existed
        if (this->treeData.count(nghID) == 0){
            // The cell is not existed, take the parent into the evaluation for list 4
            int parID = this->get_parentID(nghID);

            // Put the parent ID into the queue (only if it not existed in the queue)
            if (tempList4Flag.count(parID) == 0){
                ID_list_4_temp.push_back(parID);
                tempList4Flag.insert({parID, true});
            }
            continue;
        }
        
        // Cell is existed, do a direct evaluation
        const Cell *nghCell = this->treeData.at(nghID);

        // *[CHECK Active] Check whether the current cell is active or not
        if (nghCell->isActive == false) continue;
    
        // Check whether it is leaf cell or have no source particle
        if ((nghCell->isLeaf == true)){
            // Candidate for list 4
            ID_list_4_temp.push_back(nghID);
            tempList4Flag.insert({nghID, true});
        }
        else{
            // Check each child cell whether is well separated to the current cell (not the parent cell)
            std::vector<int> _chdIDs = this->get_childID(nghID);

            // Check the relative position of each child cell toward current cell
            double min_dist = currSize + (10 * __DBL_EPSILON__);

            // Evaluate all childs
            for (const auto &chdID : _chdIDs){
                // *[CHECK Existed] Check whether the current child is available or not (not available for outside the domain cell)
                if (this->treeData.count(chdID) == 0) continue;
                
                // Create alias to the current child cell
                const Cell *chdCell = this->treeData.at(chdID);

                // *[CHECK Active] Check whether the current cell is active or not
                if (chdCell->isActive == false) continue;

                // Initially the child ID is considered to be touching current cell
                bool isSeparated = false;
                basis_loop(d){
                    // if not touching current cell will be taken
                    if (std::abs(chdCell->centerPos[d] - currCell->centerPos[d]) > min_dist){
                        isSeparated = true;
                        break;
                    }
                }
                
                // Take the seperated ID to final list 2
                if (isSeparated == true){
                    ID_list_2.push_back(chdID);
                }
            }
        }
    }


    // **Evaluate the rest (temporary list 4)
    for (size_t i = 0; i < ID_list_4_temp.size(); i ++){
        // Create an alias to the current neighbor cell
        int nghID = ID_list_4_temp[i];

        // *[CHECK Existed] First Check whether the cell is existed
        if (this->treeData.count(nghID) == 0){
            // The cell is not existed, take the parent into the evaluation for list 4
            int parID = this->get_parentID(nghID);

            // Put the parent ID into the queue (only if it not existed in the queue)
            if (tempList4Flag.count(parID) == 0){
                ID_list_4_temp.push_back(parID);
                tempList4Flag.insert({parID, true});
            }
            continue;
        }

        // Cell is existed, do a direct evaluation
        const Cell *nghCell = this->treeData.at(nghID);

        // *[CHECK Active] Check whether the current cell is active or not
        if (nghCell->isActive == false) continue;
        
        // Check whether the cell is well seperated to the current cell
        // The current cell is a leaf cell
        double min_dist = ((nghCell->length + currSize) / 2.0) + (10 * __DBL_EPSILON__);

        // Check whether is seperated to current cell
        bool isSeparated= false;
        basis_loop(d){
            if (std::abs(nghCell->centerPos[d] - currCell->centerPos[d]) > min_dist){
                isSeparated= true;
                break;
            }
        }
        
        // Only consider the box that is well seperated
        if(isSeparated== true){
            // Check whether it is leaf cell
            if (nghCell->isLeaf == true){
                // Put the cell ID into the list
                ID_list_4.push_back(nghID);
            }
            else{
                // Take the parent ID
                int _parID = this->get_parentID(nghID);
                
                // Put the parent ID into the queue (only if it not existed in the queue)
                if (tempList4Flag.count(_parID) == 0){
                    ID_list_4_temp.push_back(_parID);
                    tempList4Flag.insert({_parID, true});
                }
                continue;
            }
        }
    }

    return;
}

// #pragma endregion

// =====================================================
// +---------------- Data Modification ----------------+
// =====================================================
// #pragma region DATA_MODIFICATION

void TreeCell::initializeTree(const std::vector<std::vector<double>> &parPos, const std::vector<bool> &activeFlag){
    // Reset the data
    this->~TreeCell();

    // Internal variable
    int parCount = activeFlag.size();   // The total count of particle to be evaluated

    // ************************* DATA INITIALIZATION *************************
    // ***********************************************************************
    /* To do list:
       1. Determine the extreme particle position
       2. Set the base length
       3. Determine the root center position
       4. Create Tmap
    */

    // PROCEDURE 1! : Find the maximum and the minimum position
    // ************
    // Internal variable
    double minPos[DIM];
    double maxPos[DIM];

    // Initialize the check
    basis_loop(d){
        minPos[d] = parPos[0][d];
        maxPos[d] = parPos[0][d];
    }
    
    // Evaluate through all particle coordinate
    for (size_t i = 1; i < parPos.size(); i ++){
        basis_loop(d){
            if (minPos[d] > parPos[i][d]) minPos[d] = parPos[i][d];
            if (maxPos[d] < parPos[i][d]) maxPos[d] = parPos[i][d];
        }
    }

    // PROCEDURE 2! : Set the base length as the longest range at each dimension direction
    // ************
    this->rootLength = 0.0;
    basis_loop(d){
        double domain_size = std::abs(maxPos[d] - minPos[d]);
        if (domain_size > this->rootLength) this->rootLength = domain_size;
    }
    
    // Adjust the length to be expanded with the given scale factor
    double expansion = 1.0 + (Pars::expTree / 100.0);
    this->rootLength *= expansion;


    // PROCEDURE 3! : Determine the pivot coordinate position
    // ************
    basis_loop(d) this->pivotCoord[d] = (minPos[d] + maxPos[d] - this->rootLength) / 2.0;


    // *************************** CELL GENERATION ***************************
    // ***********************************************************************
    /** To do list: - Perform the cell division
     *   > At each level evaluate all particle
     *   > Find the cell container for each particle at the current level
     *     - If cell is not existed, create the cell
     *     - Update the cell counter (particle number and source particle number)
     *   > Check through all cell in the current level
     *     - Identify the active and the leaf mark
     *   > Iterate through all partilce in the list
    */

    // Internal data container
    int unsettledParCount = parCount;       // The number of unsettled particle (initialized to all particle number)

    // Cell list at current level
    std::vector<int> currLvlCellList;       // List of all created cell in the current level

    // Particle identification container
    std::vector<bool> parSetFlag(parCount, false);     // Particle flag (Said that particle need to be refined)
    std::vector<int> parCellID(parCount, 0);           // Cell ID of the corresponding particle

    // ** Generate the root cell
    int ROOT_ID = 0;
    std::vector<int> rootIndex(DIM);
    basis_loop(d) rootIndex[d] = 0;

    // Create the root cell
    Cell *rootCell = new Cell(ROOT_ID, ROOT_LEVEL, this->rootLength, rootIndex, this->pivotCoord);
    
    // Update the root data
    rootCell->srcParNum = std::count(activeFlag.begin(), activeFlag.end(), true);
    rootCell->isActive = rootCell->srcParNum == 0 ? false : true;
    rootCell->parNum = parCount;
    rootCell->isLeaf = false;
    
    // Insert the root cell into tree data
    this->treeData.insert({ROOT_ID, rootCell});
    this->cellCount++;

    
    // Loop until all particle is settled into the cell leaf
    int currLevel = 0;  // The level indicator at current loop evaluation
    while (unsettledParCount > 0){
        // ** Evaluation at the level of "currLevel"

        // Proceed to the next level calculation (proceed at first, because root is already done)
        currLevel += 1;

        // Get the basic data at current level
        double currSize = this->get_cellSize(currLevel);

        // Reserve the cell list
        currLvlCellList.clear();


        // ** Evaluate particle location and generate cell
        for (int parID = 0; parID < parCount; parID++){
            // Only evaluate the unsettled particle
            if (parSetFlag[parID]) continue;

            // Calculate the cell container basic data (cell ID and index)
            std::vector<int> cell_index = this->coord2index(parPos[parID], currLevel);
            int cell_ID = this->index2ID(cell_index, currLevel);
            
            // Update the particle cell ID data (the cell at current level that contains the particle)
            parCellID[parID] = cell_ID;
            
            // Create the cell if not existed
            if (this->treeData.count(cell_ID) == 0){
                Cell *tempCell = new Cell(cell_ID, currLevel, currSize, cell_index, this->pivotCoord);
                this->treeData.insert({cell_ID, tempCell});
                this->cellCount++;

                // Add the data into the generated cell container
                currLvlCellList.push_back(cell_ID);
            }
            
            // Update the particle counter at the current cell
            this->treeData.at(cell_ID)->parNum++;
            if (activeFlag[parID] == true){
                this->treeData.at(cell_ID)->srcParNum++;
            }
        }


        // **Evaluate leaf and active cell
        // int beginID = this->get_startID(currLevel);
        // int endID = this->get_startID(currLevel+1);
        // for (int _cell_ID = beginID; _cell_ID < endID; _cell_ID++){
        for (const int &cell_ID : currLvlCellList){
            // Evaluate only if the cell is existed
            // if (this->treeData.count(cell_ID) != 0){
                // Evaluate leaf mark (if total source number below limit)
                if (this->treeData.at(cell_ID)->srcParNum <= this->max_particle){
                    // Set the current cell to leaf
                    this->treeData.at(cell_ID)->isLeaf = true;
                    this->treeData.at(cell_ID)->parIDList.clear(); // Reserve the particle ID container

                    // Put into the list
                    this->leafList.push_back(cell_ID);
                }

                // Evaluate active mark (if there existed any source particle)
                if (this->treeData.at(cell_ID)->srcParNum != 0){
                    this->treeData.at(cell_ID)->isActive = true;
                }
            // }
        }
        

        // **Rearrange the particle container at each cell at current level
        for (int parID = 0; parID < parCount; parID++){
            // Only evaluate the unsettled partilce
            if (parSetFlag[parID]) continue;

            int cell_ID = parCellID[parID];
            if (this->treeData.at(cell_ID)->isLeaf){
                // Put the current particle into the list
                this->treeData.at(cell_ID)->parIDList.push_back(parID);
                
                // Change the current particle flag
                parSetFlag[parID] = true;
                unsettledParCount--;
            }
        }
    }

    // Change the minimum level and maximum level
    int IDmin = this->leafList.at(0);
    int IDmax = this->leafList.back();

    this->min_level = this->treeData.at(IDmin)->level;
    this->max_level = this->treeData.at(IDmax)->level;

    return;
}

void TreeCell::updateTree(const std::vector<std::vector<double>> &parPos, const std::vector<bool> &activeFlag){
    return;
}

// #pragma endregion

