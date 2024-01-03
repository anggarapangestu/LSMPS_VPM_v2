#include "gridNode.hpp"

/** DONE : 
 *  Sun Nov 26 00:28:00 WIB 2023 @authors Angga'18
 *  First Test Successfully !!!
 *    > Index translation : pos2idx, idx2ID, pos2ID, ID2Index (reverse) DONE Clear
 *    > Tree done, neighbor half, adaptive only refinement
*/

// =====================================================
// +----------------- Utilites Method -----------------+
// =====================================================
// #pragma region UTILITIES_METHOD

/**
 *  @brief  Get the level at the current ID. The calculation is based
 *  on comparison between the node ID and the pivot ID at each level.
 *         
 *  @param  ID  The ID of node where level need to be found.
 * 
 *  @return The level of current ID.
*/
int GridNode::getLevel(int ID) const{
    int level = 0;
    
    // The starting ID element have not defined yet
    if (this->startID.empty()){
        // PROCEDURE 1

        // Set the ratio of geometric sequence
        int ratio = Pars::intPow(2,DIM);

        // Compare the ID to each pivot ID at the current level (calculate by geometric series)
        for (int i = 1; i < this->maxLevel + 1; i++){
            int pivID = (this->rootNodeNum * (Pars::intPow(ratio, i) - 1)) / (ratio - 1);
            if (ID < pivID) break;
            level = i;
        }
    }
    // The starting ID element have defined
    else{
        // PROCEDURE 2

        // Compare the ID to eahc pivot ID at the current level (calculate by geometric series)
        for (size_t i = 1; i < this->startID.size(); i++){
            if (ID < startID[i]) break;
            level = i;
        }
    }
    
    return level;
};

/**
 *  @brief  Get the pivot ID (starter ID) at the given level.
 *         
 *  @param  level   The grid level to be evaluated.
 * 
 *  @return Started ID at the given level.
*/
int GridNode::getPivID(int level) const{
    // Check whether the base ID is generated or not
    if (!(this->startID.empty())){
        return startID[level];
    }

    // THE BASE ID IS NOT CREATED YET!!! -> Follow the following procedure

    // Return variable
    int pivID;

    // Initial ID at the current level (i.e. lv 0 -> 0, lv 1 -> rootNodeNum, forth.)
    // The calculation includes geometric series
    
    // Set the ratio of geometric sequence
    int ratio = Pars::intPow(2,DIM);

    // Set the pivot ID at the current level (calculate by geometric series)
    pivID = (this->rootNodeNum * (Pars::intPow(ratio, level) - 1)) / (ratio - 1);
    
    return pivID;
}

// #pragma endregion

// =====================================================
// +---------------- Index Translation ----------------+
// =====================================================
// #pragma region INDEX_ID_TRANSLATION

/**
 *  @brief  Get the grid map coordinate index that contain the given spatial coordinate.
 *         
 *  @param  index   [TARGET] : Container for grid map coordinate index. 
 *  @param  coorPnt The spatial coordinate to be evaluated.
 *  @param  level   The grid level to be evaluated.
*/
void GridNode::pos2idx(int index[DIM], const double coorPnt[DIM], int level) const{
    // Internal variable
    double nodeSize;

    // Calculate the node length size at current level
    nodeSize = this->gridSize / Pars::intPow(2,level);
    
    // Determine the index position
    basis_loop(d) index[d] = ((coorPnt[d] - this->pivotCoor[d]) / nodeSize) /*+ PRECISION_FACTOR*/;  // Avoid wrong truncation
}

/**
 *  @brief  Get the grid map coordinate index from node ID.
 *         
 *  @param  index   [TARGET] : Container for grid map coordinate index. 
 *  @param  ID      Node ID to be transformed.
 *  @param  level   The grid level to be evaluated.
*/
void GridNode::ID2Index(int index[DIM], int ID, int level) const{
    /* PROCEDURE !!
        1. Translate the ID relative to this grid level.
        2. Un-flatten the given ID.
    */

    // PROCEDURE 1!
    // ************
    // Translate the ID counted from the pivot ID
    ID -= this->getPivID(level);

    // PROCEDURE 2!
    // ************
    // The given ID has been flatten, see below for 3D
    //   3D -> ID = idx + idy * (nx) + idz *(nx*ny);
    // The un-flatten method is (example by 3D)
    //   idx ->                  ID % nx = [ idx + idy * (nx) + idz *(nx*ny) ] % nx;
    //   idy ->      floor(ID / nx) % ny = [ idy + idz *(ny) ] % ny;
    //   idy -> floor(ID / (nx*ny)) % nz = [ idz ] % nz;
    
    // ** Calculate intermediate variable
    int mod[DIM];   // The moduler [nx, ny, nz] (actually the grid count at current level)
    int div[DIM];   // The divisor [1, nx, nx*ny]
    basis_loop(d){
        mod[d] = this->gridCount[d] * Pars::intPow(2,level);
        div[d] = 1;
        for (int i = 0; i < d; i++){
            div[d] *= mod[i];
        }
    }

    // ** Un-flatten transformation of node ID
    basis_loop(d) index[d] = (ID/div[d]) % mod[d];
}

/**
 *  @brief  Map a grid map coordinate to 'Node'.
 *         
 *  @param  index   Index location of a map coordinate.
 *  @param  level   The target level of node.
 * 
 *  @return ID of the current node containing the 
 *  given map coordinate at the given level.
*/
int GridNode::idx2ID(const int index[DIM], int level) const {   
    // Return variable
    int nodeID;

    /* PROCEDURE !!
        1. Calculate the starting ID at the given level.
        2. Translate the ID by flatten the node index.
    */

    // PROCEDURE 1!
    // ************
    // Calculate the starting ID at the given level.
    nodeID = this->getPivID(level);

    // PROCEDURE 2!
    // ************
    // Translate the ID by flatten the node index
    // i.e. 3D -> ID = idx + idy * (nx) + idz *(nx*ny);
    //      2D -> ID = idx + idy * (nx);
    
    // Current level grid count
    int nGrid[DIM];     // [nx,ny,nz]
    basis_loop(d) nGrid[d] = this->gridCount[d] * Pars::intPow(2, level);

    // Check whether the index is in range
    basis_loop(d){
        if (index[d] < 0 || index[d] >= nGrid[d]){
            ERROR_LOG << "Index [" << index[0];
            for (int _i = 1; _i < DIM; _i++) std::cout << "," << index[_i];
            std::cout << "] at level " << level << " is out of range!\n";
            throw std::range_error("ERROR [GRID NODE]: Cannot convert index to node ID!");
        }
    }

    // Flatten the ID
    basis_loop(d){
        int currID = index[d];
        for (int i = 0; i < d; i ++){
            currID *= nGrid[i];
        }
        nodeID += currID;
    }
    return nodeID;
}

/**
 *  @brief  Map a spatial point coordinate to 'Node'.
 *         
 *  @param  coorPnt The spatial point coordinate.
 *  @param  level   The target level of node.
 * 
 *  @return ID of the current node containing the 
 *  given spatial coordinate at the given level.
*/
int GridNode::pos2ID(const double coorPnt[DIM], int level) const {
    /* PROCEDURE !!
        1. Calculate the index of the node containing spatial coordinate
        2. Call the idx2ID method
    */

    // Internal variable
    int nodeID, index[DIM];

    // PROCEDURE 1!
    // ************
    // Get the index position
    this->pos2idx(index, coorPnt, level);
    
    // PROCEDURE 2!
    // ************
    // Calculate the nodeID from index
    nodeID = this->idx2ID(index, level);

    return nodeID;
}

// #pragma endregion

// =====================================================
// +-------------- Tree Hierarchy Method --------------+
// =====================================================
// #pragma region TREE_HIERARCHY
/**
 *  @brief  Find the parent ID.
 *         
 *  @param  currNode   Node which parent is to be found.
 * 
 *  @return Parent ID of the given node.
*/
int GridNode::findParent(const Node *currNode) const {
    // Calculate the parent node using algorithm
    // [1] Calculate the index of parent node
    // [2] Get the parent node ID using idx2ID method

    // Exceptional for root node (level = 0)
    if (currNode->level == 0){
        // This node doesn't have parent
        ERROR_LOG << "Node " << currNode->nodeID << " is a ROOT node!\n";
        throw std::range_error("ERROR [GRID NODE]: ROOT node is having no parent!");
    }

    // Parent node properties
    int parID, parIndex[DIM], parLvl;

    // Retrieve the data
    parLvl = currNode->level - 1;
    basis_loop(d) parIndex[d] = currNode->index[d] / 2;
    parID = this->idx2ID(parIndex, parLvl);
    return parID;
}

/**
 *  @brief  Find the parent ID. (Method 2)
 *         
 *  @param  currID   The ID of node which parent need to be found.
 * 
 *  @return Parent ID of the given node ID-.
*/
int GridNode::findParent(int currID) const {
// int GridNode::findParent(int index[DIM], int level) const {
    // Calculate the parent node using algorithm
    // [1] Calculate the index of parent node
    // [2] Get the parent node ID using idx2ID method

    // Calculate the current node properties
    int currIndex[DIM];
    int currLvl = this->getLevel(currID);
    this->ID2Index(currIndex, currID, currLvl);

    // Exceptional for root node (currLvl = 0)
    if (currLvl == 0){
        // This node doesn't have parent
        ERROR_LOG << "Node " << currID << " is a ROOT node!\n";
        throw std::range_error("ERROR [GRID NODE]: ROOT node is having no parent!");
    }

    // Parent node properties
    int parID, parIndex[DIM], parLvl;

    // Retrieve the data
    parLvl = currLvl - 1;
    basis_loop(d) parIndex[d] = currIndex[d] / 2;
    parID = this->idx2ID(parIndex, parLvl);
    return parID;
}

/**
 *  @brief  Find the child ID. IMPORTANT!! The current function only 
 *  use for non-LEAF node.
 *         
 *  @param  currNode   Node which children is to be found.
 * 
 *  @return Children ID's list of the given node.
*/
std::vector<int> GridNode::findChild(const Node *currNode) const {
    // Calculate the child node using algorithm
    // [1] Calculate the index of first child node
    // [2] Get the child node ID using idx2ID method

    // Exception for leaf node (currLvl = MAX_LEVEL)
    if (currNode->level == this->maxLevel){
        // This node is LEAF, the child is not available
        std::cout << FONT_PURPLE << "[WARNING]";     // Text in purple color
        std::cout << FONT_RESET << " : LEAF child is not available!\n";
    }
    
    // Child properties variable
    std::vector<int> chdIDlist(this->chdNum);
    int chdID, chdLvl;
    int baseIndex[DIM]; // Base index a.k.a. the first child index
    int chdIndex[DIM];  // Target child index
    // int baseBin[DIM];   // Base value for '&' binary operator correspond to each basis coordinate

    // Calculate the 'base index' child node and the 'base binary value'
    basis_loop(d){
        baseIndex[d] = currNode->index[d]*2;
        // baseBin[d] = Pars::intPow(2,d);
    }

    // Retrieve the data
    chdLvl = currNode->level + 1;

    // Calculate each child node
    for (int i = 0; i < this->chdNum; i++){
        // Set the current child index
        basis_loop(d)
        chdIndex[d] = baseIndex[d] + ((i>>d) & 1);    // (baseBin[d] & i); result wrong pattern
        
        // Retrieve the current child node ID
        chdID = this->idx2ID(chdIndex, chdLvl);
        chdIDlist[i] = chdID;
    }

    return chdIDlist;
}

/**
 *  @brief  Find the child ID. (Method 2)
 *         
 *  @param  currID   The ID of node which child need to be found.
 * 
 *  @return Children ID's list of the given node.
*/
std::vector<int> GridNode::findChild(int currID) const {
    // Calculate the child node using algorithm
    // [1] Calculate the index of first child node
    // [2] Get the child node ID using idx2ID method

    // Calculate the other coordinate properties
    int currIndex[DIM];
    int currLvl = this->getLevel(currID);
    this->ID2Index(currIndex, currID, currLvl);

    // Exception for leaf node (currLvl = MAX_LEVEL)
    if (currLvl == this->maxLevel){
        // This node is LEAF, the child is not available
        std::cout << FONT_PURPLE << "[WARNING]";     // Text in purple color
        std::cout << FONT_RESET << " : LEAF child is not available!\n";
    }

    // IMPORTANT!! The current function only use for non-leaf node
    std::vector<int> chdIDlist(this->chdNum);
    int chdID, chdLvl;
    int baseIndex[DIM]; // Base index a.k.a. the first child index
    int chdIndex[DIM];  // Target child index

    // Calculate the 'base index' child node and the 'base binary value'
    basis_loop(d){
        baseIndex[d] = currIndex[d]*2;
    }

    // Retrieve the data
    chdLvl = currLvl + 1;

    // Calculate each child node
    for (int i = 0; i < this->chdNum; i++){
        // Set the current child index
        basis_loop(d)
        chdIndex[d] = baseIndex[d] + ((i>>d) & 1);
        
        // Retrieve the current child node ID
        chdID = this->idx2ID(chdIndex, chdLvl);
        chdIDlist[i] = chdID;
    }

    return chdIDlist;
}

/**
 *  @brief  Find the list of sibling node of the given node (including 
 *  the current node).
 *         
 *  @param  sibList    [TARGET] : List of sibling node.
 *  @param  currNode   Node which sibling is to be found.
*/
void GridNode::findSibling(std::vector<Node*> &sibList, const Node *currNode) const{
    // [TREE LOGIC] If the current node is existed than it should has siblings !

    // Calculate the sibling node using algorithm
    // [1] Calculate the index of first sibling node
    // [2] Get the other sibling node ID using idx2ID method
    
    // Intermediate variable
    int sibNum = this->chdNum;          // Number of sibling
    sibList.resize(sibNum, nullptr);    // Note that the list includes the current Node
    int baseIndex[DIM]; // Base index a.k.a. the first sibling index
    int sibIndex[DIM];  // Target sibling index

    // Sibling parameters
    int sibID;
    int sibLvl = currNode->level;

    // Calculate the 'base index' sibling node
    basis_loop(d){
        baseIndex[d] = currNode->index[d] - (currNode->index[d] % 2);
    }

    // Calculate other sibling node
    for (int i = 0; i < sibNum; i++){
        // Set the current sibling index
        basis_loop(d) sibIndex[d] = baseIndex[d] + ((i>>d) & 1);
        
        // Retrieve the current sibling node ID
        sibID = this->idx2ID(sibIndex, sibLvl);
        
        // Assign the node ID into sibling list
        sibList[i] = this->nodeMap.at(sibID);
    }

    return;
}

// #pragma endregion

// =====================================================
// +----------- Neighbor Evaluation Method ------------+
// =====================================================
// #pragma region NEIGHBOR_EVALUATION

// Grid evaluation visualization (illustration in 2D)
/* The x mark the node to be evaluated 
     _________________
    |19|20| 10  | 11  |
    |17|_x|_____|_____|
    | 5   |           |  -> Global multi level index map
    |_____|    0      |      (ONLY a sample, sequence
    | 3   |           |         is not in a scale)
    |_____|___________|

     _____ _____ __
    |6_|7_|8_|  |
    |3_|4x|5_|__|        -> Local neighbor index map
    |0_|1_|2_|         (This map is a portion from global map
    |     |           node level take the evaluated node level)
*/

/**
 *  @brief  Find the neighbor of current node in the same level.
 *  Note: The neighbor node at ID may not be generated yet.
 *         
 *  @param  nghIDList   [TARGET] : A container which neighbor ID is stored.
 *  @param  currNode    Node which neighbor is to be found.
*/
void GridNode::findNghLvl(std::vector<int> &nghIDList, const Node *currNode) const{
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
    int nghNum = Pars::intPow(3,DIM);   // Number of neighbor candidate
    int nghIDtar;       // ID candidate of target neighbor
    int div[DIM];       // Divisor value to de-flatten neigbor local-index [1,3,9]
    int nGrid[DIM];     // Grid node count at the current level [nx, ny, nz]
    int transMul[DIM];  // ID translation multiplyer at each basis [1, nx, nx*ny]
    int transDir[DIM];  // The index translation direction at each basis
    
    // Initialize the intermediate variables
    basis_loop(d){
        div[d] = Pars::intPow(3, d);
        nGrid[d] = this->gridCount[d] * Pars::intPow(2, currNode->level);
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
            if ((currNode->index[d] == 0 && transDir[d] == -1) ||
                (currNode->index[d] == nGrid[d]-1 && transDir[d] == 1)){
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
        // nghIDtar = this->idx2ID(currNode->index, currNode->level);    // ID pivots on the current node ID
        nghIDtar = currNode->nodeID;    // ID pivots on the current node ID
        // **Translate the neighbor ID from pivot
        basis_loop(d) nghIDtar += transDir[d] * transMul[d];
        
        // SUB-PROCESS 3
        // Put the neighbor ID into the container
        nghIDList.push_back(nghIDtar);
    }
}

/**
 *  @brief  Find the neighbor leaf node. This method is purposed only for a leaf node only.
 *         
 *  @param  nghIDList   [TARGET] : A container which neighbor node is stored.
 *  @param  currNode    Node which neighbor is to be found. (Must be a leaf node)
*/
void GridNode::findNghAll(std::vector<Node*> &nghIDList, const Node *currNode) const{
    // Warning for a non leaf node
    if (!currNode->isLeaf){
        // This node is not LEAF, thus there must be a problem of program calculation
        ERROR_LOG << "Occured problem at node " << currNode->nodeID << "\n";
        throw std::runtime_error("ERROR [GRID NODE] Finding neighbor of non-LEAF node!");
    }

    // Evaluate all node including itself
    // [1] Iterate through all potential neighbor 1D->3ngh; 2D->9ngh; 3D->27ngh
    /*/ [2] Only add leaf node for neighbor 
            > if node is not existed find parent node 
            > if not a leaf node take the child node
    /*/ 

    // Find the node ID of all candidate neighbor
    std::vector<int> nghTarIDList;      // Neighbor Node ID evaluation QUEUE list
    this->findNghLvl(nghTarIDList, currNode);
    
    // Update the size of nghIDList
    nghIDList.reserve(nghTarIDList.size()*4);   // <?> Actually this line is not necessary <?>
    
    // Check every posibilies of neighbor candidate (A QUEUE EVALUATION)
    for (size_t pos = 0; pos < nghTarIDList.size(); pos++){
        // Aliasing the target ID inside ID list
        int &nghTarID = nghTarIDList[pos];

        // If the current neighbor node is the current node
        if (nghTarID == currNode->nodeID) {
            // If the current node is not LEAF than this node is not a LEAF
            nghIDList.push_back(this->nodeMap.at(nghTarID));
            continue;
        }

        // Check whether the candidate ID is existed
        if (this->nodeMap.count(nghTarID) == 0){ // ID NOT Existed!
            // To cancel recursive calculation
            // Check whether the candidate node level is bigger than the current level
            if (nghTarID >= this->startID[currNode->level + 1]){
                // There is no way a Node is not existed when its ID went to a higher level
                MESSAGE_LOG << " An attempt UP to find parent of ID " << nghTarID << " while current node ID " << currNode->nodeID << " a resulting child!\n";
                // std::cout << "[WARNING] An attempt to a recursive check is almost happenned\n";
                continue;
            }
            
            // Find the Parent ID
            int nghParID = this->findParent(nghTarID);

            // Put the parent ID into candidate neighbor list (to be checked then)
            nghTarIDList.push_back(nghParID);
        }
        else{ // ID Existed!
            // Node aliasing
            Node *nghTarNode = this->nodeMap.at(nghTarID);

            // Check whether the candidate ID is leaf node
            if (nghTarNode->isLeaf){ // ID is a Leaf Node!
                // This node is the neighbor
                nghIDList.push_back(nghTarNode);
            }
            else{ // ID Is NOT Leaf Node!
                // To cancel recursive calculation
                // Check whether the candidate node level is smaller than the current level
                if (nghTarID < this->startID[currNode->level]){
                    // There is no way a Node is not leaf when existed and have no children
                    std::cout << FONT_PURPLE << "[WARNING]";     // Text in purple color
                    std::cout << FONT_RESET << " : An attempt DOWN to find child of a resulting parent!\n";
                    // std::cout << "[WARNING] An attempt to a recursive check is almost happenned\n";
                    continue;
                }

                // Find all children ID
                std::vector<int> chdIDList = this->findChild(nghTarNode);
                
                // Put the child IDs into candidate neighbor list (to be checked then)
                for(int &chdID : chdIDList){
                    nghTarIDList.push_back(chdID);
                }
            }
        }   
    }
    
    return; // End of the check!
}

// #pragma endregion

// =====================================================
// +------------ Adaptive Operation Method ------------+
// =====================================================
// #pragma region ADAPTIVE_OPERATION

/**
 *  @brief  Perform grid refinement on the given node. The node must be a leaf node.
 *
 *  @param  currNode    Node to be refine.
 *  @param  tool        Grid Node as the basis parameter used in refinement process.
 *  @param  chdIDlist   [OUTPUT] List of all child IDs of refined Node.
 * 
 *  @return A check mark of 0 or 1 whether the node is refined (1) or not (0).
*/
int GridNode::refineNode(Node *currNode, const GridNode &tool, std::vector<int> &chdIDlist){
    // Check whether the current node is a LEAF node
    if(!currNode->isLeaf){
        // A non-LEAF node is already refined cannot refine more than once
        ERROR_LOG << "Node " << currNode->nodeID << " has already refined!\n";
        throw std::runtime_error("ERROR [GRID NODE]: An attempt to refine a non-LEAF Node!");
    }
    
    // The current node is a LEAF node -> Turn off the flag and proceed to Node Refinement
    currNode->isLeaf = false;       // ** NOTE : This node is supposed to be the one inside the this "GridNode"
    
    // PROCEDURE of Node Refinement
    // [1] Find child ID list
    // [2] Create child Node and put into map Node list

    // ** The following code is a copy of find child method
    // Intermediate variable
    int chdLvl = currNode->level + 1;       // Level of the child
    double chdLen = 0.5 * currNode->length; // Length size of the child node
    chdIDlist.resize(tool.chdNum, 0);       // Initialize the child ID container

    int baseIndex[DIM]; // Base index a.k.a. the first child index
    int chdIndex[DIM];  // Target child index
    int chdID;          // Target child ID

    // Calculate the 'base index' of child node
    basis_loop(d){
        baseIndex[d] = currNode->index[d]*2;
    }

    // Determine the child node properties, create child node, put into the map node list
    for (int i = 0; i < tool.chdNum; i++){
        // ** [1] Find the ID list
        // Calculate the current child index
        basis_loop(d) chdIndex[d] = baseIndex[d] + ((i>>d) & 1);
        // Retrieve the current child node ID
        chdID = tool.idx2ID(chdIndex, chdLvl);
        // Put the ID into child ID container
        chdIDlist[i] = chdID;

        // ** [2] Create node and insert into map node list
        // Create the child node
        Node* chdNode = new Node(chdID, chdLvl, chdLen, chdIndex);
        // Update the pivot coordinate node
        basis_loop(d) chdNode->pivCoor[d] = tool.pivotCoor[d] + (chdIndex[d] * chdLen);
        // Upadte the active sign [!] Custom 
        chdNode->isActive = currNode->isActive;
        // Assign the node into node map
        this->nodeMap.insert({chdID, chdNode});
    }

    return 1;
}

/**
 *  @brief  Perform a particle list division from the current node to it's child-LEAF node.
 *
 *  @param  currNode     Node which the particle list to be divided.
 *  @param  baseParticle Basis particle object to provide individual particle coordinate
 *  for the calculation.
 *  @param  chdIDlist    List of all child IDs.
 *  @param  type    Type of particle division, [1] REFINEMENT, [DEFAULT] COMPRESSION
 * 
 *  @return A check mark of 0 or 1 whether the node is refined (1) or not (0).
*/
int GridNode::divideParticle(Node *currNode, const Particle &templatePar, const std::vector<int> &chdIDlist, int type){
    // Give a check, check whether all child Node is valid
    // bool invalid = true;
    for (auto &chdID : chdIDlist){
        // No target child
        if (this->nodeMap.count(chdID) == 0){
            ERROR_LOG << "Node " << chdID << " as child of node " << currNode->nodeID << " is not existed!\n";
            throw std::runtime_error("ERROR [GRID NODE] Unable to get the child, the node is not existed in the grid!");
        }

        // Reserve the child particle list data
        this->nodeMap.at(chdID)->parList.clear();
        
        // Additional for parameter update for adaptation (REFINEMENT)
        if (type == 1) this->nodeMap.at(chdID)->headNodeID = currNode->headNodeID;

        // // Child list is invalid 
        // if (chdID != 0){
        //     invalid = false;
        // }
    }

    // // Error Handling
    // if (invalid){
    //     throw std::runtime_error("ERROR [GRID NODE] Invalid child ID list!");
    // }

    // Internal variable
    double chdNodeLength = currNode->length / 2.0;
    int posIdx[DIM];        // Temporary child node index container
    double parPos[DIM];     // Temporary particle coordinate holder
    
    // Move the particle into the child node
    for (int &_particleID : currNode->parList){
        // Retrieve the position of the particle
        parPos[0] = templatePar.x[_particleID];
        if (DIM > 1) parPos[1] = templatePar.y[_particleID];
        if (DIM > 2) parPos[2] = templatePar.z[_particleID];
        basis_loop(d){
            posIdx[d] = (parPos[d] - currNode->pivCoor[d]) / chdNodeLength;
            // Undefined behaviour handler
            if (posIdx[d] < 0 || posIdx[d] > 1){
                // The local index is out of range
                ERROR_LOG << "Particle " << _particleID << " is outside the current node box!\n";
                throw std::range_error("ERROR [GRID NODE] Child index is out of range!");
            }
        }
        
        // Find the child local index
        int chdLocID = 0;
        // This following line is where the previous wrong algorithm occur [The line issue has resolved]
        basis_loop(d) chdLocID = chdLocID*2 + posIdx[(DIM-1) - d];

        // Copy the particle ID into the corresponding child Node
        Node *&_chdNode = this->nodeMap.at(chdIDlist[chdLocID]);
        _chdNode->parList.push_back(_particleID);
    }

    // Additional for parameter update for adaptation (REFINEMENT)
    if (type == 1){
        // Release the adaptation (REFINMENT) parameter of the current node
        currNode->headNodeID = -1;
        currNode->parList.clear();
    }

    return 1;
}

// #pragma endregion

// =====================================================
// +--------------- Saving Data Method ----------------+
// =====================================================
// #pragma region DATA_SAVING_METHOD

/**
 *  @brief  Write all node in grid into system.
 *         
 *  @param  nodeList The "GridNode" data to be saved.
 *  @param  name The identification name for the saved file.
*/
void GridNode::saveGrid(const GridNode& nodeList, std::string name) const {
    // Initialization of saving data
    std::ofstream writer;
    name = "output/gridNode_all_" + name + ".csv";
    writer.open(name.c_str());
    
    // Write data table header
    if (DIM == 2){
        writer << "ID,xId,yId,x,y,lvl,size,isLeaf\n";
    }else if (DIM == 3){
        writer << "ID,xId,yId,zId,x,y,z,lvl,size,isLeaf\n";
    }

    // Write data value
    for(auto &[_ID,_node] : nodeList.nodeMap){
        writer << _node->nodeID ;
        basis_loop(d)
        writer << "," << _node->index[d];
        basis_loop(d)
        writer << "," << _node->pivCoor[d] + (0.5 * _node->length);
        writer << "," << _node->level
               << "," << _node->length
               << "," << _node->isLeaf
               << "\n";
    }
    
    writer.close();
}

/**
 *  @brief  Write all leaf node in grid into system.
 *         
 *  @param  nodeList The "GridNode" data to be saved.
 *  @param  name The identification name for the saved file.
*/
void GridNode::saveLeafGrid(const GridNode& nodeList, std::string name) const {
    // Initialization of saving data
    std::ofstream writer;
    name = "output/gridNode_leaf_" + name + ".csv";
    writer.open(name.c_str());
    
    // Write data table header
    if (DIM == 2){
        writer << "ID,xId,yId,x,y,lvl,size,isLeaf\n";
    }else if (DIM == 3){
        writer << "ID,xId,yId,zId,x,y,z,lvl,size,isLeaf\n";
    }

    // Write data value
    for(auto &[_ID,_node] : nodeList.nodeMap){
        if (!_node->isLeaf) continue;
        writer << _node->nodeID;
        basis_loop(d)
        writer << "," << _node->index[d];
        basis_loop(d)
        writer << "," << _node->pivCoor[d] + (0.5 * _node->length);
        writer << "," << _node->level
               << "," << _node->length
               << "," << _node->isLeaf
               << "\n";
    }
    
    writer.close();
}

/**
 *  @brief  Write all node specified in the list.
 *         
 *  @param  nodeList The "GridNode" data to be saved.
 *  @param  IDList The list of node ID to be write.
 *  @param  name The identification name for the saved file.
*/
void GridNode::saveSelectedGrid(const GridNode& nodeList, const std::vector<int> &IDList, std::string name) const{
    // Initialization of saving data
    std::ofstream writer;
    name = "output/gridNode_selected_" + name + ".csv";
    writer.open(name.c_str());
    
    // Table header
    if (DIM == 2){
        writer << "ID,xId,yId,x,y,lvl,size,isLeaf\n";
    }else if (DIM == 3){
        writer << "ID,xId,yId,zId,x,y,z,lvl,size,isLeaf\n";
    }

    // Data value
    for(const int &_ID : IDList){
        const Node *_node = nodeList.nodeMap.at(_ID);

        writer << _ID ;
        basis_loop(d)
        writer << "," << _node->index[d];
        basis_loop(d)
        writer << "," << _node->pivCoor[d] + (0.5 * _node->length);
        writer << "," << _node->level
               << "," << _node->length
               << "," << _node->isLeaf
               << "\n";
    }
    
    writer.close();
}

// #pragma endregion

// =====================================================
// +------------------- Old Method --------------------+
// =====================================================
// #pragma region OLD_METHOD

// /**
//  *  @brief  Write all leaf node in grid into system.
//  *         
//  *  @param  nodeList The "GridNode" data to be saved.
// */
// void GridNode::saveLeafGrid(GridNode& nodeList) const {
//     std::ofstream writer;
//     writer.open("output/gridNodeLeaf.csv");
//     writer.precision(10);

//     std::cout << "[LOG] SAVE THE DATA!\n";
    
//     // Table header
//     if (DIM == 2){
//         writer << "ID,xId,yId,x,y,lvl,size,isLeaf\n";
//     }else if (DIM == 3){
//         writer << "ID,xId,yId,zId,x,y,z,lvl,size,isLeaf\n";
//     }

//     // Data value
//     for(auto &[_ID,_node] : nodeList.nodeMap){
//         if (!_node->isLeaf) continue;

//         writer << _ID ;
//         basis_loop(d)
//         writer << "," << _node->index[d];
//         basis_loop(d)
//         writer << "," << _node->pivCoor[d] + (0.5 * _node->length);
//         writer << "," << _node->level
//                << "," << _node->length
//                << "," << _node->isLeaf
//                << "\n";
//     }
    
//     writer.close();
// }

// /** REPLACE THIS CODE WITH JUST THE ONE WITH TOOLS
//  *  @brief  Perform grid refinement on the given node.
//  *
//  *  @param  currNode    Node to be refine.
//  *  @param  chdIDlist   [OUTPUT] List of all child IDs of refined Node.
//  * 
//  *  @return A check mark of 0 or 1 whether the node is refined (1) or not (0).
// */
// int GridNode::refineNode(Node *currNode, std::vector<int> &chdIDlist){
//     // Check whether the node need to be refine <?> I think this location is push out from this method <?>
//     if(!currNode->needRefinement){
//         // The node not need to refine
//         std::cout << FONT_PURPLE << "[WARNING]";     // Text in purple color
//         std::cout << FONT_RESET << " : An attempt to refine a no-need-refinement Node!\n";
//         return 0;
//     }
    
//     // Refine the node procedure
//     // [1] Find child ID list
//     // [2] Create child Node and put into map Node list
//     // [3] Change the refinement flag of the current node
    
//     // [4] Turn off the leaf node <?> Is this necessary <?>

//     // ** The following code is a copy of find child method
//     // Intermediate variable
//     int chdNum = Pars::intPow(2,DIM);
//     chdIDlist.resize(chdNum);

//     int baseIndex[DIM]; // Base index a.k.a. the first child index
//     int chdIndex[DIM];  // Target child index
//     int chdID;
//     int chdLvl = currNode->level + 1;
//     double chdLen = 0.5 * currNode->length;

//     // Calculate the 'base index' child node and the 'base binary value'
//     basis_loop(d){
//         baseIndex[d] = currNode->index[d]*2;
//     }

//     // Determine the child node properties, create child node, put into the map node list
//     for (int i = 0; i < chdNum; i++){
//         // ** [1] Find the ID list
//         // Calculate the current child index
//         basis_loop(d) chdIndex[d] = baseIndex[d] + ((i>>d) & 1);
//         // Retrieve the current child node ID
//         chdID = this->idx2ID(chdIndex, chdLvl);
//         chdIDlist[i] = chdID;

//         // ** [2] Create node and insert into map node list
//         // Create the child node
//         Node* chdNode = new Node(chdID, chdLvl, chdLen, chdIndex);
//         // Update the pivot coordinate node
//         basis_loop(d) chdNode->pivCoor[d] = this->pivotCoor[d] + (chdIndex[d] * chdLen);
//         // Assign the node into node map
//         this->nodeMap.insert({chdID, chdNode});
//     }

//     // ** [3] Turn off the refinement flag of current node
//     currNode->needRefinement = false;

//     // ** [4] Turn off the leaf node <?> Is this necessary <?> -> Yes this is necessary
//     currNode->isLeaf = false;

//     return 1;
// }


// /** THIS METHOD IS ABOUT TO BE DELETED
//  *  @brief  Perform grid compression on the given node.
//  *
//  *  @param  currNode    Node to be refine.
//  * 
//  *  @return A check mark of 0 or 1 whether the node is refined (1) or not (0).
// */
// int GridNode::compressNode(Node *currNode, int &parID){
//     // Find the current node siblings
//     std::vector<Node*> siblingList;
//     this->findSibling(siblingList, currNode);
    
//     // Check whether the node and each sibling need to be compressed
//     for (Node *&sib : siblingList){
//         // Check the flag of current sibling
//         if(!sib->needCompression){
//             // Turn off the compression flag of all the sibling
//             for (Node *&_sib : siblingList){
//                 _sib->needCompression = false;
//             }
//             return 0;
//         }
//     }
    
//     // Compress the node procedure
//     // [1] Find the parent ID
//     // [2] Delete all of the sibling from all container and free the data
    
//     // [3] Turn on the leaf node <?> Is this necessary <?>

//     // ** [1] Find the parent node
//     parID = this->findParent(currNode);
//     Node* parNode = this->nodeMap.at(parID);

//     // ** [2] Delete the sibling from all container and free the data
//     int sibID;
//     for (Node *&sib : siblingList){
//         // Retrieve the ID of each sibling
//         sibID = sib->nodeID;

//         // IMPORTANT!! Delete the element that contain current sibling in any container
//         // ADD IF ANY!
//         this->nodeMap.erase(sibID);
//         // this->leafNode.erase(sibID);

//         // Free the node
//         delete sib;
//     }
    
//     // ** [3] Turn on the leaf node <?> Is this necessary <?>
//     parNode->isLeaf = true;
    
//     return 1;
// }

// /** THIS METHOD IS ABOUT TO DELETE
//  *  @brief  Find the adjacent neighbor leaf node.
//  *         
//  *  @param  nghIDList   [TARGET] : A container which neighbor node is stored.
//  *  @param  currNode    Node which neighbor is to be found. (Must be a leaf node)
// */
// void GridNode::findNghAdj(std::vector<Node*> &nghIDList, const Node *currNode) const{
//     // Evaluate all node including itself
//     // [1] Iterate through all potential neighbor 1D->3ngh; 2D->9ngh; 3D->27ngh
//     /*/ [2] Only add leaf node for neighbor
//             > if node is not existed find parent node
//             > if not a leaf node take the child node
//     /*/
//     // Intermediate variable 1
//     int div[DIM];       // Divisor value to de-flatten neigbor local-index [1,3,9]
//     int nGrid[DIM];     // Grid node count at the current level [nx, ny, nz]
//     int transMul[DIM];  // ID translation multiplyer at each basis [1, nx, nx*ny]
    
//     // Update the intermediate variable
//     basis_loop(d){
//         div[d] = Pars::intPow(3, d);
//         nGrid[d] = this->gridCount[d] * Pars::intPow(2, currNode->level);
//         transMul[d] = 1;
//         for (int i = 0; i < d; i++){
//             transMul[d] *= nGrid[i];
//         }
//     }
        
//     // Intermediate variable 2
//     int nghNum = Pars::intPow(3,DIM);   // Number of neighbor candidate
//     int nghIDtar;           // ID candidate of target neighbor
//     int transPos[DIM];      // The index translation direction at each basis
    
//     // nghIDList.reserve(nghNum*4);    // Update the size of nghIDList
    
//     // Check every posibilies of neighbor candidate
//     for (int i = 0; i < nghNum; i++){
//         // SUB-PROCESS 1
//         // Calculate the translation direction and check range
//         bool inRange = true;
//         // **Translation direction
//         basis_loop(d){
//             transPos[d] = ((i / div[d]) % 3) - 1;
//             if ((currNode->index[d] == 0 && transPos[d] == -1) ||
//                 (currNode->index[d] == nGrid[d]-1 && transPos[d] == 1)){
//                 // No node candidate in this target location
//                 inRange = false;
//                 break;
//             }
//         }
//         // **Check in range
//         if (!inRange) continue;

        
//         // SUB-PROCESS 2
//         // Initialize the neighbor candidate ID
//         nghIDtar = currNode->nodeID;    // ID pivots on current node ID
//         // Translate the neighbor ID from pivot
//         basis_loop(d)
//         nghIDtar += transPos[d] * transMul[d];

//         // Exception of current ID!
//         if (nghIDtar == currNode->nodeID) continue;

        
//         // SUB-PROCESS 3
//         // Check whether the candidate ID is existed
//         if (this->nodeMap.count(nghIDtar) == 0){ // ID Not Existed!
//             // Check parent (<?> The parent must be leaf <?>)
//             // IMPORTANT : Modify this code block is needed!
//             int parID = this->findParent(nghIDtar);
//             nghIDList.push_back(this->nodeMap.at(parID));
//         }else{ // ID Existed!
//             // Check whether the candidate ID is leaf node
//             if (this->nodeMap.at(nghIDtar)->isLeaf){ // ID Is Leaf Node!
//                 // Check current
//                 nghIDList.push_back(this->nodeMap.at(nghIDtar));
//             }else{ // ID Is Not Leaf Node!
//                 // Check children (<?> The children must be leaf <?>)
//                 // IMPORTANT : Modify this code block is needed!
//                 std::vector<int> childs = this->findChild(this->nodeMap.at(nghIDtar));
//                 for(int &chdID : childs){
//                     nghIDList.push_back(this->nodeMap.at(chdID));
//                 }
//             }
//         }   
//     }
    
//     return; // End of the check!
// }

// /**
//  *  @brief  Perform grid division on the given node.
//  *
//  *  @param  currNode    Node to be refine.
//  *  @param  chdIDlist   [OUTPUT] List of all child IDs of refined Node.
//  * 
//  *  @return A check mark of 0 or 1 whether the node is refined (1) or not (0).
// */
// int GridNode::divideParticleComp(Node *currNode, const Particle &templatePar, const std::vector<int> &chdIDlist){
//     // Give a check, check whether all child Node is valid
//     bool invalid = true;
//     for (auto &chdID : chdIDlist){
//         // Additional for parameter update for adaptation (COMPRESSION)
//         this->nodeMap.at(chdID)->parList.clear();

//         // No target child
//         if (this->nodeMap.count(chdID) == 0){
//             throw std::runtime_error("ERROR [GRID NODE] Unable to get the child, the node is not existed in the grid!");
//         }

//         // Child list is invalid 
//         if (chdID != 0){
//             invalid = false;
//         }
//     }

//     // Error Handling
//     if (invalid){
//         throw std::runtime_error("ERROR [GRID NODE] Invalid child ID list!");
//     }

//     // Internal variable
//     double chdNodeLength = currNode->length / 2.0;
//     int posIdx[DIM];        // Temporary child node index container
//     double parPos[DIM];     // Temporary particle coordinate holder

//     // Move the particle into the child node
//     for (int &_particleID : currNode->parList){
//         // Retrieve the position of the particle
//         parPos[0] = templatePar.x[_particleID];
//         if (DIM > 1) parPos[1] = templatePar.y[_particleID];
//         if (DIM > 2) parPos[2] = templatePar.z[_particleID];
//         basis_loop(d){
//             posIdx[d] = (parPos[d] - currNode->pivCoor[d]) / chdNodeLength;
//             // Undefined behaviour handler
//             if (posIdx[d] < 0 || posIdx[d] >= this->chdNum){
//                 // The local index is out of range
//                 throw std::range_error("ERROR [GRID NODE] Child index is out of range!");
//             }
//         }
        
//         // Find the child local index
//         int chdLocID = 0;
//         // This following line is where the previous wrong algorithm occur [The line issue has resolved]
//         basis_loop(d) chdLocID = chdLocID*2 + posIdx[(DIM-1) - d];

//         // Copy the particle ID into the corresponding child Node
//         Node *&_chdNode = this->nodeMap.at(chdIDlist[chdLocID]);
//         _chdNode->parList.push_back(_particleID);
//     }
//     return 1;
// }

// #pragma endregion

