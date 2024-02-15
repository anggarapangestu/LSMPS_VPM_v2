#include "generateGrid.hpp"

#define BODY_EXT_MUL 2.0        // The body extension multiplier factor

/** Root node generation
 *   [1] : Recursive attempt root node generation
 *   [2] : Direct manner of root node generation
 *   [3] : Old method
*/
#define PROCEDURE_2_TYPE 2

/** Root node refinement based on obstacle body
 *   [1] : Direct body node evaluation
 *   [2] : Body surface expansion evaluation
*/
#define PROCEDURE_3_TYPE 2

/** The type of base number of particle at each node
 *   [1] : Geometric series of 2^n (rounded up from support radius)
 *   [2] : Size of support radius (rounded up)
 *   [3] : A user defined value
*/
#define BASE_NUM_TYPE 2

// ===================================================== //
// +------------------- Constructor -------------------+ //
// ===================================================== //
// #pragma region CONSTRUCTOR

/**
 *  @brief  Default constructor: Generate the value of class member.
*/
generateGrid::generateGrid(){
    // ** The following definition is fixed, no need any modification

    // Set the maximum level
    this->maxLevel = Pars::max_level;

    // Set the maximum level different
    this->NghlevelDiff = Pars::ngh_diff_level;

    // Determine the number of particle inside a node
    /* CATEGORY : 
        > The value must be an element in geometric sequence of base 2 with ratio of 2
            so that <?> Is this a must <?>
        > The value must inside the length of buffer radius
    */

    // // ** NOTE : Modify the @baseParNum calculation if necessary
    if (BASE_NUM_TYPE == 1){
        int exp = std::ceil(std::log2(Pars::r_sup * Pars::r_buff));
        this->baseParNum = Pars::intPow(2,exp);
    }
    else if(BASE_NUM_TYPE == 2){
        this->baseParNum = std::ceil(Pars::r_sup * Pars::r_buff);
    }
    else if(BASE_NUM_TYPE == 3){
        this->baseParNum = 3;
    }
    // NOTE: as the number of 'this->baseParNum' getting high, the cost of neighbor evaluation is getting up to

    // Calculate the size of the finest block
    this->leafBlockSize = this->baseParNum * Pars::sigma;
    this->rootBlockSize = this->leafBlockSize * Pars::intPow(2, this->maxLevel);
}


/**
 *  @brief  Default deconstructor: Nothing to do here.
*/
generateGrid::~generateGrid(){
    // Nothing to do here
}

// #pragma endregion

// ===================================================== //
// +----------------- Private Methods -----------------+ //
// ===================================================== //
// #pragma region PRIVATE_METHOD

/**
 *  @brief  [PRIVATE METHOD] Update all of the "GridNode" parameter based on the calculated member.
 *         
 *  @param  nodeList    The "GridNode" data which parameter to be updated.
*/
void generateGrid::updateGridNode(GridNode& nodeList){
    // ** NOTE : Modify the following 3 member variables definition if necessary by reffering the data to global.cpp

    // Set the length of each direction
    this->length[0] = Pars::lxdom;
    this->length[1] = Pars::lydom;
    if (DIM == 3) this->length[2] = Pars::lzdom;    // This method is safe from illegal memory access

    // Set the minimum global coordinate position
    this->minCoor[0] = - Pars::xdom;
    this->minCoor[1] = - Pars::lydom / 2.0;
    if (DIM == 3) this->minCoor[2] = - Pars::lzdom / 2.0;
    
    // Set the maximum global coordinate position
    this->maxCoor[0] = Pars::lxdom - Pars::xdom;
    this->maxCoor[1] = Pars::lydom / 2.0;
    if (DIM == 3) this->maxCoor[2] = Pars::lzdom / 2.0;

    // ** Start updating the Grid Node data
    
    // Calculate the number of nodes (count) at ROOT level at [1] each dimension and [2] the total number
    nodeList.rootNodeNum = 1;
    basis_loop(d){
        nodeList.gridCount[d] = std::ceil(this->length[d] / this->rootBlockSize);   // [1]
        nodeList.rootNodeNum *= nodeList.gridCount[d];                              // [2]
    }
    
    // Update the domain size (fit with the number of nodes)
    // NOTE*: A symmetrical expansion from original size is performed
    basis_loop(d){
        double excess = 0.5 * (nodeList.gridCount[d] * this->rootBlockSize - this->length[d]);
        this->minCoor[d] -= excess;
        this->maxCoor[d] += excess;

        // Provide an unsymetrical value (To make early separation)
        // ** Only in y and z direction (d > 0)
        if (Pars::flag_slightly_shifted_domain && (d > 0)){
            this->minCoor[d] += Pars::mp_shift; /*Pars::sigma / 3.0;*/
            this->maxCoor[d] += Pars::mp_shift; /*Pars::sigma / 3.0;*/
        }
    }

    // Update the data inside the "GridNode"
    nodeList.maxLevel = this->maxLevel;         // Update the level limit
    nodeList.baseParNum = this->baseParNum;     // Update number of particle in each node (in one direction)
    nodeList.gridSize = this->rootBlockSize;    // gridSize : length size of ROOT node
    basis_loop(d){
        nodeList.pivotCoor[d] = this->minCoor[d];
    }

    // Update the starting ID at each level (For transformation)
    nodeList.startID.resize(this->maxLevel + 2,0);
    int multiplier = Pars::intPow(2,DIM);
    for (int i = 0; i < this->maxLevel + 1; i++){
        nodeList.startID[i+1] = nodeList.startID[i] + nodeList.rootNodeNum * Pars::intPow(multiplier,i);
    }

}

/**
 *  @brief  [PRIVATE METHOD] A recursive function toward each basis to generate the root node data.
 *         
 *  @param  nodeList    The "GridNode" data for root node to be constructed.
 *  @param  dim The highest dimension in the simulation (should be inputted with DIM).
 *  @param  ID  The initial ID in the root level (should be inputted with 0).
 *  @param  index   Array of grid map index as intermediate container.
*/
void generateGrid::generateRootRec(GridNode& nodeList, int dim, int &ID, int index[DIM]) const{
    // Followed by recursion at each basis that
    // will only generate the node when meet the first basis (x)

    // Reach the first basis
    if (dim == 1){
        for (int i = 0; i < nodeList.gridCount[dim-1]; i++){
            // Update the position index
            index[dim-1] = i;

            // Create the node
            Node* temp = new Node(ID, ROOT_LEVEL, nodeList.gridSize, index);
            
            // Update the pivot coordinate node
            basis_loop(d) temp->pivCoor[d] = 
            nodeList.pivotCoor[d] + (index[d] * nodeList.gridSize);

            // Assign the node into node map
            nodeList.nodeMap.insert({ID, temp});
            ID++;  // Proceed to the next node generation
        }
    }
    // Not reaching the first basis (proceed the recursion)
    else{
        for (int i = 0; i < nodeList.gridCount[dim-1]; i++){
            // Update the position index
            index[dim-1] = i;
            this->generateRootRec(nodeList, dim-1, ID, index);
        }
    }
    
    return;
}

/**
 *  @brief  [PRIVATE METHOD] A direct root node generation.
 *         
 *  @param  nodeList    The "GridNode" data for root node to be constructed.
*/
void generateGrid::generateRootDir(GridNode& nodeList) const{
    // The content of this function actually copying from "GridNode::ID2Index" method

    // ** Intermediate variable
    int mod[DIM];   // ITERATOR MODULER: [nx, ny, nz] (actually the grid count at root level)
    int div[DIM];   // ITERATOR DIVISOR: [1, nx, nx*ny]
    int index[DIM]; // The index at each ID
    int ID;         // Node ID for the iteration
    basis_loop(d){
        mod[d] = nodeList.gridCount[d];
        div[d] = 1;
        for (int i = 0; i < d; i++){
            div[d] *= mod[i];
        }
    }

    // ** Generate all node by iterating the ID
    for (ID = 0; ID < nodeList.rootNodeNum; ID++){
        // Calculate the index at each dimension
        basis_loop(d) index[d] = (ID/div[d]) % mod[d];

        // Create the node
        Node* temp = new Node(ID, ROOT_LEVEL, nodeList.gridSize, index);
        
        // Update the pivot coordinate node
        basis_loop(d) temp->pivCoor[d] = 
        nodeList.pivotCoor[d] + (index[d] * nodeList.gridSize);

        // Assign the node into node map
        nodeList.nodeMap.insert({ID, temp});
    }
    
    return;
}

// #pragma endregion

// ===================================================== //
// +------------------ Public Method ------------------+ //
// ===================================================== //
// #pragma region PUBLIC_METHOD

/**
 *  @brief  Set the value of neighbor level different in the class member.
 *  This is for further needs, only utilized if this is found to be necessary.
 *         
 *  @param  diff The new level different value to be set into the class.
*/
void generateGrid::setNghLvlDiff(int diff){
    NghlevelDiff = diff;
}

/**
 *  @brief  Generate the data of node list generated from global domain.
 *         
 *  @param  nodeList The "GridNode" data for list of node to be constructed.
 *  @param  bodyList vector contain 'Body' object in simulation domain.
*/
void generateGrid::nodeGeneration(GridNode& nodeList, const std::vector<Body>& bodyList){
    /* PROCEDURE !!
        1. Set up (update) the member of GridNode
        2. Generate the root node list
        3. Refine the node based on body geometry
        4. Recursively check the neighbor condition (delta level <= 1)
    */

    // PROCEDURE 1!
    // ************
    // Set up all of the GridNode parameter
    MESSAGE_LOG << "Update grid node parameters\n";
    this->updateGridNode(nodeList);
    
    // PROCEDURE 2!
    // ************
    // Generate the root node
    MESSAGE_LOG << "Generating root node\n";
    if (PROCEDURE_2_TYPE == 1){
        // Generate the root node by recursive attempt
        // The indexing is moving from pivot to z to y to x consecutively (recursive method)
        int ID = 0;         // ID counter
        int index[DIM];     // Position index
        this->generateRootRec(nodeList, DIM, ID, index);
    }
    else if (PROCEDURE_2_TYPE == 2){
        // Generate the root node in direct manner
        this->generateRootDir(nodeList);
    }
    
    // // Saving ROOT
    // std::string name = "ROOT";
    // nodeList.saveLeafGrid(nodeList, name);
    
    // PROCEDURE 3!
    // ************
    // Refine the node base of the body geometry
    /*/ There are 2 method: 
        [1] store all node ID to be refined, or 
        [2] refine each level at the same time with iteration  <---  Chosen
    /*/

    // Additional (*To be utilized in PROCEDURE 4)
    //  > In this process, the refined leaf node ID is stored
    //  > This node will be evaluated for neighbor level different check

    // ID container [QUEUE LIST] (*Used in PROCEDURE 4)
    const int LEVEL_BOUND = this->maxLevel - this->NghlevelDiff;    // The upper bound of level limit to be evaluated in queue (*not including the bound)
    std::vector<std::vector<int>> IDqueue(LEVEL_BOUND);
    std::unordered_map<int,bool> IDflag;

    // START PROCEDURE 3:
    // Change the method procedure in the header file
    MESSAGE_LOG << "Root node refinement\n";
    if (PROCEDURE_3_TYPE == 1){
        /* TYPE 1 Description:
            > At each level evaluate the node ID containing the body point
            > For each body point, find the node ID container then refine the node
            > Promptly do the refinement evaluation for each level
        */

        /*/ Current procedure:
            i) At each level, find the node where each body point is located
            ii) Then refine and proceed to the next level
        /*/
        // 1st loop : iterate in each level
        // 2st loop : iterate each Body in bodyList
        // 3st loop : iterate each node in Body

        // Internal variables
        double bodyNodeCoor[DIM];
        int currLvl, toRefineNodeID;

        // 1st loop (through all level from ROOT to one level before MAX_LEVEL)
        for(currLvl = ROOT_LEVEL; currLvl < nodeList.maxLevel; currLvl++){
            
            // 2nd loop (through all body)
            for (auto &obstacle : bodyList){
                
                // 3rd loop (through all node in body)
                for (int i = 0; i < obstacle.n_panel; i++){
                    // Retrieve the position of the body coordinate (body panel mid point coordinate)
                    bodyNodeCoor[0] = obstacle.x_m[i];
                    bodyNodeCoor[1] = obstacle.y_m[i];
                    if (DIM == 3) bodyNodeCoor[2] = obstacle.z_m[i];

                    // Take the node ID to be refined
                    toRefineNodeID = nodeList.pos2ID(bodyNodeCoor,currLvl);

                    // Put a safety check
                    // The current node is not created yet
                    // but previously created instead, TROUBLESHOOT: wrong index calculation (index jump)
                    if (nodeList.nodeMap.count(toRefineNodeID) == 0) {
                        ERROR_LOG << "Node " << toRefineNodeID << " is missing!\n";
                        throw std::runtime_error("ERROR [GENERATE NODE] : Unable to refine the node, the node is missing!");
                    }

                    // Refine the current ID
                    // Check whether the current ID is a leaf Node
                    Node *&currNode = nodeList.nodeMap.at(toRefineNodeID);
                    if(currNode->isLeaf){
                        // Container to store the child ID
                        std::vector<int> chdIDList;

                        // Perform the refinement
                        nodeList.refineNode(currNode, nodeList, chdIDList);

                        // Initialize the parameter for PROCEDURE 4!
                        // Store the ID into the IDqueue container
                        if (currLvl + 1 < LEVEL_BOUND){
                            for (int &chdID : chdIDList){
                                IDqueue[currLvl + 1].push_back(chdID);
                                IDflag.insert({chdID, true});
                            }
                        }
                    }
                }
            }
        }
    
    }
    else if (PROCEDURE_3_TYPE == 2){
        /* TYPE 2 Description:
            > Collect all root node near any of body inside the domain
            > At each level check find the minimum distance from node to body
            > Refine the node if the distance inside criteria
        */

        // Internal variable
        int currLvl;            // Current grid level evaluated
        geometry geom_operator; // Operator to obtain the minimum distance from body

        // A node ID container for distance check
        std::vector<int> nodeIDList1;
        std::vector<int> nodeIDList2;
        std::unordered_map<int,bool> nodeIDFlag;

        // ** [1] Collecting all root node in criteria
        // LOOP -> Check all node near the body through all body inside the domain
        for (auto &_body : bodyList){
            // The extreme body position
            double min_coor[DIM];
            double max_coor[DIM];
            basis_loop(d){
                min_coor[d] = _body.min_pos[d] - (BODY_EXT_MUL * Pars::body_ext);
                max_coor[d] = _body.max_pos[d] + (BODY_EXT_MUL * Pars::body_ext);
            }

            // The extreme index position
            int min_index[DIM];
            int max_index[DIM];
            nodeList.pos2idx(min_index, min_coor, ROOT_LEVEL);
            nodeList.pos2idx(max_index, max_coor, ROOT_LEVEL);

            // The Internal Variable
            int pivID;          // Pivot ID : minimum node ID at the current body
            int tarID;          // Target ID : calculated ID at each iteration
            int numIter;        // Total number of iteration need to do

            // Local flatten index to local matrix index modifier
            int div[DIM];       // Divisor [dx, dy, dz]     -> Inside the range of body
            int mod[DIM];       // Moduler [1, dx, dx*dy]   -> Inside the range of body
            
            // Local matrix index to global ID modifier
            int transVal[DIM];  // The translation multiplier [1, Nx, Nx*Ny]
            int transDir[DIM];  // The translation direction at each basis
            
            // Update all internal variable
            pivID = nodeList.idx2ID(min_index, ROOT_LEVEL);  // ID of Node at minimum position
            numIter = 1;
            basis_loop(d){
                mod[d] = 1 + max_index[d] - min_index[d];
                numIter *= mod[d];

                div[d] = 1;
                transVal[d] = 1;
                for (int i = 0; i < d; i++){
                    div[d] *= mod[i];
                    transVal[d] *= nodeList.gridCount[i];
                }
            }

            // Loop through all iteration to get the ID list
            for (int i = 0; i < numIter; i++){
                // Calculate the ID of the node
                tarID = pivID;      // Initialize by the pivot ID
                basis_loop(d) {
                    transDir[d] = (i/div[d]) % mod[d];
                    tarID += transDir[d] * transVal[d];
                }
                
                // Insert the data into the queue container
                // Make sure no double ID (e.g. overlap Node with more than 1 objects)
                if(nodeIDFlag.count(tarID) == 0){
                    nodeIDList1.push_back(tarID);
                    nodeIDFlag.insert({tarID,true});
                }
            }

        } // Done iterating each body

        // The task of current variable is done
        nodeIDFlag.clear();     // Free the container for next usage

        // ** [2] Check the minimum distance of each node toward body
        // LOOP -> Promptly check at each level
        for(currLvl = ROOT_LEVEL; currLvl < nodeList.maxLevel; currLvl++){
            // Create an alias for each container
            std::vector<int> *currList;     // Evaluated at current level
            std::vector<int> *nextList;     // Stored for next level
            
            // Alternating the container 1 and 2 as current and next
            if (currLvl % 2 == 0){
                currList = &nodeIDList1;
                nextList = &nodeIDList2;
            }else {
                currList = &nodeIDList2;
                nextList = &nodeIDList1;
            }

            // Free the next ID list container
            nextList->clear();

            // BEFORE LOOP INITIALIZATION
            // Distance initialization (*as square value of max domain length)
            double MAX_DISTANCE = 0.0;
            basis_loop(d) MAX_DISTANCE = std::max(MAX_DISTANCE, this->length[d]*this->length[d]);

            // LOOP -> Check through all Node
            for (int &_ID : *currList){
                // Alias for the current node
                Node *&_node = nodeList.nodeMap.at(_ID);
                
                // Determine the middle point coordinate of the node
                std::vector<double> nodePos(DIM);       // Node middle point coordinate
                double _radius = 0.5 * _node->length;   // Node half length
                basis_loop(d) nodePos[d] = _node->pivCoor[d] + _radius;

                // ** Check the minimum distance of the node from all body
                double _dist = MAX_DISTANCE;   // Initialize the minimum distance to body

                // LOOP -> Check through all body to find the minimum distance
                for (auto &_body : bodyList){
                    _dist = std::min(_dist, std::abs(geom_operator.distance_calc(nodePos, _body, true)));
                }

                // ** Check the distance criteria for refinement
                double _radCheck = (_radius * std::sqrt(DIM)) + (BODY_EXT_MUL * Pars::body_ext); // Node diagonal radius (*with body extension)
                if (_dist < _radCheck){
                    // [CHECK LOG] Node containing the portion of body!
                    //    -> Refine current node
                    //    -> Put child node ID into next ID list container
                    
                    // Container to store the child ID
                    std::vector<int> chdIDList;

                    // Perform the refinement
                    nodeList.refineNode(_node, nodeList, chdIDList);

                    // Store the child ID into the next ID container
                    for (int &chdID : chdIDList){
                        nextList->push_back(chdID);
                    }

                }else if (currLvl < LEVEL_BOUND){
                    // [CHECK LOG] No body portion inside the node!
                    //    -> Put the node ID into queue container for PROCEDURE 4 for NLD check
                    
                    // Put the ID into queue ID list
                    IDqueue[currLvl].push_back(_ID);
                    IDflag.insert({_ID, true});
                }
            }
        }
    }

    // // Saving ROOT
    // std::string name = "REFINE";
    // nodeList.saveLeafGrid(nodeList, name);

    // PROCEDURE 4!
    // ************
    // Refine the node according to neihgbor level different (NLD) criterion
    // Do a loop check through all nodes that need an NLD evaluation
    
    // Operation:
    // [1] Loop check until there are no nodes left on the queue
    //      > Queue container [IDqueue] and [IDflag]
    // [2] Loop through all nodes on each level (from highest level in the queue)
    // [3] At each node check the neighbor node ID (Further operation see inside the code)

    
    // START PROCEDURE 4:
    MESSAGE_LOG << "Evaluate NLD criteria and refinement\n";

    // LOOP [1]: Loop until no node left on the queue
    bool stop = false;      // Loop trigger
    int fineLvl = this->maxLevel - this->NghlevelDiff - 1;  // Starting level at each loop
    while(!stop){
        // Activate the trigger of loop termination
        // The loop need to halt when the queue at all level is empty
        stop = true;
        
        // LOOP [2]: Check the ID inside the queue from the finest level
        for (int currLvl = fineLvl; currLvl >= 0; currLvl--){
            // Cancel the trigger if the queue is not empty
            if (!(IDqueue[currLvl].empty())){
                stop = false;
            }

            // LOOP [3]: Check the Node inside the queue at current iteration level
            for (int &_ID : IDqueue[currLvl]){
                // Node Evaluation Procedure
                // [0] Check if leaf node : only evaluate leaf node
                // [1] Find the neighbor ID at the current level
                // [2] Check the node at neighbor ID
                //      > [COND 1] The finest neighbor level different is larger than MAX_LEVEL_DIFF
                //          -> Refine the current node -> Put the refined node into container
                //          -> Break the condition one if the current node is need to refine
                //      > [COND 2] Existed: The neigbor with level lower than current level
                //          -> Put the ID into the queue in accordance to its level

                // Create alias
                Node *&_Node = nodeList.nodeMap.at(_ID);

                // ** [0] Check leaf node
                if (!(_Node->isLeaf)){
                    // Proceed to the next Node in queue
                    WARNING_LOG << " NLD check on a non-LEAF node! [CODE:0xNG!!]\n";
                    continue;
                }

                // ** [1] Find the neighbor ID at the current level
                std::vector<int> nghIDList;             // List of the neighbor ID
                std::unordered_map<int,bool> nghIDflag; // Flag of existed neighbor
                
                // Update each list
                nodeList.findNghLvl(nghIDList, _Node);  // Update neigbor list
                for (int &_nghID : nghIDList){          // Update neigbor flag
                    nghIDflag.insert({_nghID, true});
                }

                // ** [2] Check each neigbor node
                bool needRefine = false;
                // LOOP [4]: Check each neigbor node
                for (size_t i = 0; i < nghIDList.size(); i++){
                    // Retrieve the neighbor ID
                    int _nghID = nghIDList[i];

                    // Evaluation exception!!!
                    if (needRefine){
                        if(_nghID >= nodeList.getPivID(_Node->level+1)){
                            // Skip the evaluation for finer neighbor
                            //   if the current node is need to be refined
                            continue;
                        }
                    }
                    
                    // Evaluation:
                    // [A] Node not existed
                    //       -> Navigate to the parent node, put the parent node to ngh container
                    // [B] Node existed -> Check if leaf node
                    //       -> Evaluate the node [COND 2]
                    // [C] Node existed but not a leaf node
                    //       -> Evaluate that node need to refine [COND 1]
                    //       -> If [COND 1] fulfilled stop entering this section for further evaluation
                    //       -> If not fullfilled -> Navigate to child node, put the parent node to ngh container

                    // ** [A] Check whether the node is not existed
                    if (nodeList.nodeMap.count(_nghID) == 0){ // Node is NOT existing
                        // To cancel recursive calculation
                        // Check whether the candidate node level is bigger than the current level
                        if (_nghID >= nodeList.startID[currLvl + 1]){
                            // There is no way a Node is not existed when its ID went to a higher level
                            printf("%s[WARNING]%s : An attempt UP to find parent of a resulting child!\n", FONT_PURPLE, FONT_RESET);
                            // std::cout << "[WARNING] An attempt UP to a recursive check is almost happenned\n";
                            continue;
                        }
                        
                        // Find parent ID
                        int nghParID = nodeList.findParent(_nghID);
                        
                        // Put the parent node ID into the neighbor ID list to be evaluated further
                        if (nghIDflag.count(nghParID) == 0){
                            // Only put the node when it's not existing in the list
                            nghIDflag.insert({nghParID, true});
                            nghIDList.push_back(nghParID);
                        }
                    }
                    // ** [B] Check whether a leaf node
                    else if (nodeList.nodeMap.at(_nghID)->isLeaf){ // Is a leaf Node
                        
                        // Check the level of the current neighbor node
                        int _nghLvl = nodeList.getLevel(_nghID);

                        // ** [COND 2] Put the ID into queue container (*if meet the criteria)
                        if (_nghLvl < _Node->level){
                            // Add the node ID to the queuqe
                            if (IDflag.count(_nghID) == 0){
                                // Only add if the ID is not existing in the queue
                                IDqueue[_nghLvl].push_back(_nghID);
                                IDflag.insert({_nghID, true});
                            }
                        }
                    }
                    // ** [C] Check whether a the current node still not need to refine
                    else if (!needRefine){ 
                        // To cancel recursive calculation
                        // Check whether the candidate node level is smaller than the current level
                        if (_nghID < nodeList.startID[currLvl]){
                            // There is no way a Node is not a leaf when existed and have no children
                            printf("%s[WARNING]%s : An attempt DOWN to find child of a resulting parent!\n", FONT_PURPLE, FONT_RESET);
                            // std::cout << "[WARNING] An attempt of Node " << _nghID << " , Eval: "<< _Node->nodeID << " DOWN to a recursive check is almost happenned\n";
                            continue;
                        }

                        // Evaluate whether the current node need to be refined
                        int _nghChdLvl = nodeList.getLevel(_nghID) + 1;
                        if (_nghChdLvl > _Node->level + this->NghlevelDiff){
                            // Meet the criteria to refine
                            needRefine = true;
                            continue;   // Continue to the next iteration
                        }

                        // Find all children ID
                        std::vector<int> chdIDList = nodeList.findChild(nodeList.nodeMap.at(_nghID));
                        
                        // Put the child IDs into candidate neighbor list (to be evaluated further)
                        for(int &chdID : chdIDList){
                            nghIDflag.insert({chdID, true});
                            nghIDList.push_back(chdID);
                        }
                    }
                } /* LOOP [4]: Check each neigbor node */

                // ** [COND 1] Perform the refinement (*if meet the criteria)
                if (needRefine){
                    // Container to store the child ID
                    std::vector<int> chdIDList;
                    
                    // Refine the node
                    nodeList.refineNode(_Node, nodeList, chdIDList);

                    // Store the child ID into the queue container (*if meet the criteria)
                    if (currLvl + 1 < LEVEL_BOUND){
                        // Must be not existed inside the queue yet because the node was just created
                        for (int &chdID : chdIDList){
                            IDqueue[currLvl + 1].push_back(chdID);
                            IDflag.insert({chdID, true});
                        }
                    }
                }
            } /* LOOP [3] -> Check at each neighbor */

            // Release the queue container at the current level
            for (int &_ID : IDqueue[currLvl]){
                IDflag.erase(_ID);          // Release the ID flag
            }
            IDqueue[currLvl].clear();       // Release the ID list

        } /* LOOP [2] -> Check at each level */
    
    } /* LOOP [1] -> Loop through all queue */

    // END OF PROCEDURE 4

    // // [DEBUG LINE] Print the neccessary thing
    // std:: cout << "The count ot the node :\n";
    // basis_loop(d) std::cout << "At " << d+1 << " : " << nodeList.gridCount[d] << "\n";
    // std:: cout << "The maximum level : " << this->maxLevel << " \n";
    // std:: cout << "The number of particle inside  : " << nodeList.baseParNum << " \n";
}

/**
 *  @brief  Generate the data of node list for no body simulation.
 *         
 *  @param  nodeList The "GridNode" data for list of node to be constructed.
*/
void generateGrid::nodeGeneration(GridNode &nodeList){
    /* PROCEDURE !!
        1. Set up (update) the member of GridNode
        2. Generate the root node list
        3. Refine until the finest node
    */

    // PROCEDURE 1!
    // ************
    // Set up all of the GridNode parameter
    MESSAGE_LOG << "Update grid node parameters\n";
    this->updateGridNode(nodeList);
    
    // PROCEDURE 2!
    // ************
    // Generate the root node
    MESSAGE_LOG << "Generating root node\n";
    if (PROCEDURE_2_TYPE == 1){
        // Generate the root node by recursive attempt
        // The indexing is moving from pivot to z to y to x consecutively (recursive method)
        int ID = 0;         // ID counter
        int index[DIM];     // Position index
        this->generateRootRec(nodeList, DIM, ID, index);
    }
    else if (PROCEDURE_2_TYPE == 2){
        // Generate the root node in direct manner
        this->generateRootDir(nodeList);
    }
    
    // PROCEDURE 3!
    // ************
    // Refine the node until reaching the maximum level
    
    // Iterate from the ROOT level to one level before MAX_LEVEL
    for(int currLvl = ROOT_LEVEL; currLvl < nodeList.maxLevel; currLvl++){
        // Iterate through all node inside the current level
        int beginID = nodeList.getPivID(currLvl);
        int endID = nodeList.getPivID(currLvl+1);
        for (int nodeID = beginID; nodeID < endID; nodeID++){
            // Refine the current ID
            // Check whether the current ID is a leaf Node
            Node *&currNode = nodeList.nodeMap.at(nodeID);

            // Container to store the child ID
            std::vector<int> chdIDList;

            // Perform the refinement
            nodeList.refineNode(currNode, nodeList, chdIDList);
        }
        // Done current level
    }

    return;
}

/**
 *  @brief  Generate the data of node list based on particle data.
 *         
 *  @param  nodeList The "GridNode" data for list of node to be constructed.
 *  @param  particle The given particle data to generate grid.
*/
void generateGrid::createNode(GridNode &_grid, Particle &_par){
    /* PROCEDURE !!
        1. Set up (update) the member of GridNode
        2. Put the particle into the node
        3. Recursively create parent node
    */

    // PROCEDURE 1!
    // ************
    MESSAGE_LOG << "Update grid node parameters\n";
    // Determine the domain size
    for (int i = 0; i < _par.num; i++){
        this->maxCoor[0] = std::max<double>(_par.x[i] + 0.5*_par.s[i], this->maxCoor[0]);
        this->minCoor[0] = std::min<double>(_par.x[i] - 0.5*_par.s[i], this->minCoor[0]);

        this->maxCoor[1] = std::max<double>(_par.y[i] + 0.5*_par.s[i], this->maxCoor[1]);
        this->minCoor[1] = std::min<double>(_par.y[i] - 0.5*_par.s[i], this->minCoor[1]);

        if (DIM > 2){
        this->maxCoor[2] = std::max<double>(_par.z[i] + 0.5*_par.s[i], this->maxCoor[2]);
        this->minCoor[2] = std::min<double>(_par.z[i] - 0.5*_par.s[i], this->minCoor[2]);
        }
    }
    
    // Calculate the domain size
    basis_loop(d)
    this->length[d] = this->maxCoor[d] - this->minCoor[d];

    // ** Start updating the Grid Node data
    // Calculate the number of nodes (count) at ROOT level at [1] each dimension and [2] the total number
    _grid.rootNodeNum = 1;
    basis_loop(d){
        _grid.gridCount[d] = std::round(this->length[d] / this->rootBlockSize); // [1]
        _grid.rootNodeNum *= _grid.gridCount[d];                                // [2]
    }
    
    // Update the domain size (fit with the number of nodes)
    // NOTE*: Expand the domain but keep the pivot location (minCoor)
    basis_loop(d){
        double excess = (_grid.gridCount[d] * this->rootBlockSize - this->length[d]);
        this->maxCoor[d] += excess;
    }

    // Update the data inside the "GridNode"
    _grid.maxLevel = this->maxLevel;         // Update the level limit
    _grid.baseParNum = this->baseParNum;     // Update number of particle in each node (in one direction)
    _grid.gridSize = this->rootBlockSize;    // gridSize : length size of ROOT node
    basis_loop(d){
        _grid.pivotCoor[d] = this->minCoor[d];
    }

    // Update the starting ID at each level (For transformation)
    _grid.startID.resize(this->maxLevel + 2,0);
    int multiplier = Pars::intPow(2,DIM);
    for (int i = 0; i < this->maxLevel + 1; i++){
        _grid.startID[i+1] = _grid.startID[i] + _grid.rootNodeNum * Pars::intPow(multiplier,i);
    }

    // PROCEDURE 2!
    // ************
    // Create node at current particle position

    // Create the list of leaf node
    std::unordered_map<int, bool> leafNodeFlag;

    // Start to put particle into each node
    MESSAGE_LOG << "Put particle into node\n";
    _par.nodeID.resize(_par.num);
    for (int i = 0; i < _par.num; i++){
        // Get the particle coordinate
        double parCoor[DIM];
        parCoor[0] = _par.x[i];
        parCoor[1] = _par.y[i];
        if (DIM > 2)
        parCoor[2] = _par.z[i];
        
        // Create the node at the corresponding gridNode
        int nodeLevel = _par.level[i];      // Node level is the particle level
        int nodeIndex[DIM];                 // Node index calculated from particle position
        _grid.pos2idx(nodeIndex, parCoor, nodeLevel);
        int nodeID = _grid.idx2ID(nodeIndex, nodeLevel);

        // Put the current node ID into the leaf node
        if (leafNodeFlag.count(nodeID) == 0){
            // Update the leaf node list
            leafNodeFlag.insert({nodeID, true});

            // Create a new node into the grid node
            double nodeLen = _grid.gridSize / Pars::intPow(2, nodeLevel);
            Node* node = new Node(nodeID, nodeLevel, nodeLen, nodeIndex);
            // Update the pivot coordinate node
            basis_loop(d) 
            node->pivCoor[d] = _grid.pivotCoor[d] + (nodeIndex[d] * nodeLen);
            // Assign the node into node map
            _grid.nodeMap.insert({nodeID, node});
        }

        // Create alias to the current node
        Node *&node = _grid.nodeMap.at(nodeID);
        
        // Push the particle into the particle list inside the node
        node->parList.push_back(i);

        // Update the particle nodeID
        _par.nodeID[i] = nodeID;
    }
    
    // PROCEDURE 3!
    // ************
    // Create the parent node
    MESSAGE_LOG << "Create parent node\n";
    for (const auto &[_nodeID, _flag] : leafNodeFlag){
        // Check the flag
        if (!_flag) continue;
        // *Proceed if the flag is true

        // Create an evaluation node
        const Node *_evalNode = _grid.nodeMap.at(_nodeID);

        // Create parent recursively
        while (_evalNode->level > ROOT_LEVEL){
            // Create parent
            int parID = _grid.findParent(_evalNode->nodeID);

            // Check if parent already created or not
            if (_grid.nodeMap.count(parID) != 0) break;
            // * Stop the current iteration it the parent already existed

            // Set up the parent node properties
            int parLevel = _evalNode->level - 1;
            int parIndex[DIM];
            _grid.ID2Index(parIndex, parID, parLevel);
            double parLen = _evalNode->length * 2.0;
            
            // Create the parent node
            Node* parNode = new Node(parID, parLevel, parLen, parIndex);

            // Turn off the leaf flag
            parNode->isLeaf = false;
            
            // Update the pivot coordinate node
            basis_loop(d) 
            parNode->pivCoor[d] = _grid.pivotCoor[d] + (parIndex[d] * parLen);
            
            // Assign the node into node map
            _grid.nodeMap.insert({parID, parNode});
            
            // **Change the alias of evaluation node
            _evalNode = parNode;
        }
    }
}


/**
 *  @brief  Assign the particle node ID into the corresponding leaf node. Also
 *  map the evaluated particle into a map node.
 *
 *  @param  baseGrid The base node data for node ID evaluation.
 *  @param  mapNode [OUTPUT] The particle ID map data consisted of leaf node only.
 *  This container maps the nodeID to all ID of evalulated particle inside the node 
 *  (or can be illustrated as : [nodeID]->{[parID]}).
 *  @param  evalParticle [UPDATE] Particle data to set the node ID.
*/
void generateGrid::assignNodeID(const GridNode &nGrd, std::unordered_map<int, std::vector<int>> &nodeMap, Particle &par){
    /** Procedure:
     *  1. Evaluate the ROOT node ID containing each particle
     *  2. Push down the node ID until it reaches the leaf node
    */

    // Intermediate variable
    
    // Particle size on root node
    const double rootParSize = nGrd.gridSize / nGrd.baseParNum;

    // Node container variable
    std::vector<int> nodeList1, nodeList2;          // The temporary node ID container
    std::vector<int> *currNodeList, *nextNodeList;  // The temporary alias of node ID container
    std::unordered_map<int, bool> nodeFlag;         // The flag of the available node [Helper variable]

    // // Particle parameter variable
    // std::vector<bool> settleFlag(par.num, false);   // List of settle particle
    // int unsettledNum = par.num;     // Number of unsettled particle

    // PROCEDURE 1:
    // ***********
    // Determine the ID of the corresponding root node
    par.nodeID.resize(par.num, 0);  // Reserve the node ID container
    par.s.resize(par.num, 0);       // Reserve the particle size container
    for (int parID = 0; parID < par.num; parID++){
        // Retrieve the particle coordinate
        double coord[DIM];
        coord[0] = par.x[parID];
        coord[1] = par.y[parID];
        #if (DIM == 3 )
        coord[2] = par.z[parID];
        #endif

        // Update the particle size
        par.s[parID] = rootParSize;

        // Get the root ID
        int _nodeID = nGrd.pos2ID(coord, ROOT_LEVEL);
        par.nodeID[parID] = _nodeID;

        // Update the temporary container
        if (nodeFlag.count(_nodeID) == 0){
            nodeList1.push_back(_nodeID);
            nodeFlag.insert({_nodeID, true});
        }
        
        // Update the nodeID to parID map container data
        nodeMap[_nodeID].push_back(parID);
    }

    // PROCEDURE 2:
    // ***********
    // Iteratively evaluate leaf and downpass to child for non leaf node
    
    // Set the level evaluation
    int level;

    // Start iteration
    for (level = 0; level < Pars::max_level; level++){
        // Create aliasing for the node list
        if (level % 2 == 0){
            // Even level
            currNodeList = &nodeList1;
            nextNodeList = &nodeList2;
        }else{
            // Odd Level
            currNodeList = &nodeList2;
            nextNodeList = &nodeList1;
        }

        // Set up the necessary container
        nodeFlag.clear();
        nextNodeList->clear();

        // Calculate the particle size on the child node
        const double chdParSize = rootParSize / Pars::intPow(2, level+1);

        // Check each node in the current list
        for (const int &nodeID : *currNodeList){
            // // Check the current node is not a leaf
            // if(nGrd.nodeMap.at(nodeID)->isLeaf == true){
            //     // Set the particle size
            //     double _nodeSize = nGrd.gridSize / Pars::intPow(2,level);
            //     double _parSize = _nodeSize / nGrd.baseParNum;

            //     // Alias to the particle list
            //     std::vector<int> &parList = nodeMap.at(nodeID);
                
            //     // Update the particle size
            //     for (const int &parID : parList){
            //         par.s[parID] = _parSize;
            //     }
            // }
            
            // Check the current node is not a leaf
            if (nGrd.nodeMap.at(nodeID)->isLeaf != true){
                // Distribute the particle into the correspond child (if current node is NOT a leaf)

                // Internal container
                std::vector<int> &parList = nodeMap.at(nodeID);         // List of all particle in the current node
                std::vector<int> chdIDList = nGrd.findChild(nodeID);    // List of all child node ID of the current node

                // Coordinate container
                double nodeCoord[DIM];  // Node mid point coordinate
                double coord[DIM];      // Particle coordinate

                // Retrieve the middle coordinate of current node
                const Node &_node = *(nGrd.nodeMap.at(nodeID));
                basis_loop(d) nodeCoord[d] = _node.pivCoor[d] + (0.5 * _node.length);
                
                // Evaluate all particle in the current node
                for (const int &parID : parList){
                    // Retrieve the particle coordinate
                    coord[0] = par.x[parID];
                    coord[1] = par.y[parID];
                    #if (DIM == 3 )
                    coord[2] = par.z[parID];
                    #endif

                    /** Illustration of LID (local ID) numbering
                     *   2D space       y                                    *  Comparing value in the binary means!
                     *   _________     |                                     *  -----------------------------------
                     *  | 2  | 3  |    |____ x                               *    2D space only have 2 basis x,y
                     *  |____|____|                    Back Block            *    it can be translated to 2 digit binary
                     *  | 0  | 1  |    3D space      ___________________     *     For clarity, apply the LID number in 2D space
                     *  |____|____|                /:        /:        /|    *     0 -> 00 : x left to center  | y below center
                     *                            /_:______ /_:______ / |    *     1 -> 01 : x right to center | y below center
                     *       Front Block         |  :  2   |  :  3   |  |    *     2 -> 10 : x left to center  | y above center
                     *      ___________________  |  :......|..:......|..|    *     3 -> 11 : x right to center | y above center
                     *    /:        /:        /| | ,:      | ,:      | /|    *          ^^
                     *   /_:______ /_:______ / | |,_:______|,_:______|/ |    *       y--||--x
                     *  |  :  6   |  :  7   |  | |  :  0   |  :  1   |  |    *  
                     *  |  :......|..:......|..| |  :......|..:......|..|    *  The first digit (rightmost digit) represents
                     *  | ,:      | ,:      | /| | ,       | ,       | /     *   the relative x location of child toward parent node center, 
                     *  |,_:______|,_:______|/ | |,________|,________|/      *   otherwise the second digit represent the relative y location.
                     *  |  :  4   |  :  5   |  |        y                    *  
                     *  |  :......|..:......|..|       |                     *  The translation value between child relative position
                     *  | ,       | ,       | /        |____ x               *   and the LID value is used to determine the child node ID
                     *  |,________|,________|/        /                      *  
                     *                             z /                       *  This idea works fine for the extension to 3D!  
                    */
                    int LID = 0;    // The local ID position of child location [initialized by 0]

                    // Determine the child local position [quad (2D) or octant (3D)] by
                    // comparing the particle coordinate toward current node mid point
                    basis_loop(d){
                        LID |= ((coord[d] > nodeCoord[d] ? 1 : 0) << d);
                        // NOTE:
                        // > Use the binary left shift (<<) to set the digit location
                        // > Use the binary or operator to overwrite binary value to LID
                    }

                    // Update the particle size
                    par.s[parID] = chdParSize;

                    // Get the new ID
                    int newNodeID = chdIDList[LID];     // The new node ID
                    par.nodeID[parID] = newNodeID;      // Update the node ID of current particle

                    // Update the temporary container
                    if (nodeFlag.count(newNodeID) == 0){
                        nextNodeList->push_back(newNodeID);
                    }
                    nodeFlag.insert({newNodeID, true});

                    // Update the nodeID to parID map container data
                    nodeMap[newNodeID].push_back(parID);
                }

                // Release the container at current node
                nodeMap.erase(nodeID);
                nodeFlag.erase(nodeID);
            }
        }
    }

    /** Procedure:
     *  1. Put all particle into the corresponding root node (also take the other one)
     *     -> For leaf node     : Keep as settled particle
     *     -> For non leaf node : Transfer to the corresponding child
     *  2. Do the first and second check
     *  3. Done when the level reach max (automatically a leaf node)
    */

    return;
}

// #pragma endregion

