#include "gridNodeAdapt.hpp"

// ===================================================== //
// +------------ Adaptation Initialization ------------+ //
// ===================================================== //
// #pragma region INITIALIZATION

/**
 *  @brief  Set up the target resolution level of each LEAF node.
 *         
 *  @param  baseGrid [UPDATE] The node container which node TRL is to be updated.
 *  @param  tempPar  The template particle to retrieve particle coordinate data.
 *  @param  selProp  The container of the selected properties of each particle, sorted by particle ID.
*/
void GridNodeAdapt::targetLevel(GridNode &baseGrid, const Particle &tempPar, const std::vector<double> &selProp){
    // In this method, the value of "selected properties" of each particle are being compared to the MAX_VALUE
    // The target resolution level of the given particle takes the matched windows, see the following illustration

    /* ILLUSTRATION OF Target Resolution Level (TRL)
        For the sake of this illustration, take an example value for each parameter as follow.
        let. MAX_LEVEL = 3, MAX_VALUE = 5.0, TOLERANCE = 0.1, and VALUE = 0.1
        
        [LARGEST VALUE]
          VAL 0 -----------------------------> 5.000 :(MAX_VALUE * (TOLERANCE)^0)
            ^
            |    WINDOW of TRL -> [3] (LEAF) : MAX_LEVEL
            |
          VAL 1 -----------------------------> 0.500 :(MAX_VALUE * (TOLERANCE)^1)
            ^
            |    WINDOW of TRL -> [2] : MAX_LEVEL - 1
            |        [x] The particle value [0.1] lies in this window
            |     thus the TRL of current particle is [2]
            |
          VAL 2 -----------------------------> 0.050 :(MAX_VALUE * (TOLERANCE)^2)
            ^    
            |    WINDOW of TRL -> [1] : MAX_LEVEL - 2
            |
          VAL 3 -----------------------------> 0.005 :(MAX_VALUE * (TOLERANCE)^3)
            ^    
            |    WINDOW of TRL -> [0] (ROOT) : MAX_LEVEL - 3
            |
          ABS 0 -----------------------------> 0
        [LOWEST VALUE]
    */

    // Internal variable
    int tarLv;      // The target resolution level of evaluated particle

    // Set all nodes tarResLv in the GridNode to ROOT prior to the evaluation
    for (auto &[_ID,_node] : baseGrid.nodeMap){
        // Reset the node target level
        _node->tarResLv = ROOT_LEVEL;
    }

    // // [DEBUG LINE]
    // std::vector<int> TARGET_LVL_LIST(this->maxLevel+1,0);
    
    // Evaluate through all particle
    for (size_t _parID = 0; _parID < selProp.size(); _parID++){
        // Check the TRL of current particle
        double const &currValue = selProp[_parID];
        double tolerance = this->maxValue * Pars::adapt_tol;
        for (tarLv = baseGrid.maxLevel; tarLv > 0; tarLv--){
            if (currValue >= tolerance){
                break;
            }else{
                tolerance *= Pars::adapt_tol;   
            }
        }

        // Update the TRL of the node container
        // Aliasing to the node where the current particle is contained
        Node *&_node = baseGrid.nodeMap.at(tempPar.nodeID[_parID]);
        // The node TRL takes the maximum TRL from all particle inside it
        _node->tarResLv = std::max<int>(_node->tarResLv, tarLv);

        // // [DEBUG LINE]
        // TARGET_LVL_LIST[tarLv]++;
    }

    // // [DEBUG LINE]
    // for (int _lvl = 0; _lvl < TARGET_LVL_LIST.size(); _lvl++){
    //     std::cout << "The number of particle targeting level " << _lvl << " : " << TARGET_LVL_LIST[_lvl] << "\n";
    // }

    // Flag for adaptation
    // #pragma omp parallel for
    for (size_t _parID = 0; _parID < selProp.size(); _parID++){
        if (tempPar.isNearSurface[_parID]){
            // Update the TRL of the node container
            // Aliasing to the node where the current particle is contained
            Node *&_node = baseGrid.nodeMap.at(tempPar.nodeID[_parID]);
            // The node containing near surface particle will be refine to the maximum level
            _node->tarResLv = Pars::max_level;
        }
    }


    // [!] A custom code section: Add the section to evaluate the active sign
    for (auto &[_ID,_node] : baseGrid.nodeMap){
        if (_node->isLeaf == false) continue;
        // Only set the active mark for leaf node
        bool activeSign = false;
        for (const auto &parID : _node->parList){
            if (tempPar.isActive[parID] == true){
                activeSign = true;
                break;
            }
        }
        _node->isActive = activeSign;
    }

    return;
}

/**
 *  @brief  Set the flag of node adaptation. Based on the node target adaptation level (TRL),
 *  group the LEAF into three parts "compressList", "refineList", and "idleList".
 *         
 *  @param  baseGrid   The grid node of the evaluation 
*/
void GridNodeAdapt::set_adaptation_flag(GridNode &baseGrid){
    // Iterate through all Node in baseGrid
    for (auto &[_ID,_node] : baseGrid.nodeMap){
        // Only evaluate the LEAF node
        if (_node->isLeaf){
            // Aliasing of the target level and current level of the node
            int &currLv = _node->level;
            int &tarLv = _node->tarResLv;

            // Check whether need adaptation or refinement or not
            if (tarLv < currLv){
                // GROUP IN COMPRESSION LIST
                this->compressList.insert({_ID, true});
                
                // Set up the compression flag
                baseGrid.nodeMap.at(_ID)->needCompression = true;
            }
            else if (tarLv > currLv){
                // GROUP IN REFINEMENT LIST
                this->refineList.insert({_ID, true});
                
                // Duplicate the node into the new grid list
                Node *_dupNode = new Node(_node);
                _dupNode->headNodeID = _dupNode->nodeID;        // Set for adaptation
                this->tempGrid.nodeMap.insert({_ID,_dupNode});
            }
            else{
                // GROUP IN IDLE LIST
                this->idleList.insert({_ID, true});
            }
        }
    }
}

// #pragma endregion

// ===================================================== //
// +--------------- Particle Adaptation ---------------+ //
// ===================================================== //
// #pragma region PARTICLE_UPDATE

/**
 *  @brief  Evaluate the neighbor of the given particle in the list ("parIDList") from
 *  the old particle as the neighbor candidate ("parNghIDList").
 *         
 *  @param  parIDList   List of particle ID which neighbor need to be evaluated.
 *  @param  parNghIDList  List of particle ID as neighbor candidate.
 *  @param  nghSrcPar   The particle source data of neighbor candidate.
*/
void GridNodeAdapt::findParticleNeighbor(const std::vector<int> &parIDList, const std::vector<int> &parNghIDList, const Particle &nghSrcPar){
    // Evaluate all particle in the list
    // #pragma omp parallel for         // [IMPORTANT] This paralel make a problem, slow computation
    for (size_t i = 0; i < parIDList.size(); i++){

        // Alias to the current particle to be evaluated
        const auto &evalParID = parIDList[i];

        // Check the distance criteria on all neighbor candidate particle
        for (size_t j = 0; j < parNghIDList.size(); j++){
            // const auto nghParID : parNghIDList
            const auto &nghParID = parNghIDList[j];

            // The internal variable (distance squared)
            double R2 = 0;

            // Calculate the distance square in x direction
            double dx = this->newPar->x[evalParID] - nghSrcPar.x[nghParID];
            dx = dx*dx;
            R2 += dx;

            // Calculate the distance square in y direction
            if (DIM > 1){
                double dy = this->newPar->y[evalParID] - nghSrcPar.y[nghParID];
                dy = dy*dy;
                R2 += dy;
            }

            // Calculate the distance square in z direction
            if (DIM > 2){
                double dz = this->newPar->z[evalParID] - nghSrcPar.z[nghParID];
                dz = dz*dz;
                R2 += dz;
            }

            // Check neighbor distance criteria (Using the base size of old data)
            // Support radius ('s' is the size of Head particle)
            double supRad = this->newPar->s[evalParID] * Pars::r_sup * Pars::r_buff;
            
            // Remember checking the squared radius
            if (R2 < supRad*supRad){
                // Inside the support radius, put this particle inside the particle neighbor
                this->newPar->neighbor[evalParID].push_back(nghParID);
            }
        }
    }

    return;
}

/**
 *  @brief  Generate the particle in the given node. 
 *  NOTE: The size of particle is using the size of particle in the head node.
 *         
 *  @param  currNode  The node where the particle are to be generated.
 *  @param  headSize  The size of particle in the head node.
*/
void GridNodeAdapt::generateParticle(Node *currNode, double _headSize){
    // Generate the particle in the current Node

    // Update the current node data
    currNode->parList.clear();  // Reserve the particle ID container on the current node 
    currNode->headNodeID = -1;  // Turn off the head node ID flag

    // Size of the particle inside the node
    double _parSize = currNode->length / this->parNum;
    // double _headSize = headNode->length / this->parNum;

    // Generate all particle inside the current Node
    int locParIndex[DIM];     // The particle local index (taken from the node)
    for (int locParID = 0; locParID < this->totParNum; locParID++){
        // Put the current particle ID into the node
        currNode->parList.push_back(this->newPar->num);    // Push the last pushed ID inside the parList

        // The local index coordinate inside the node
        basis_loop(d) locParIndex[d] = (locParID/this->div[d]) % this->parNum;

        // Assign the particle position
        // The x position
        double _x = currNode->pivCoor[0] + (0.5 + locParIndex[0])*_parSize;
        this->newPar->x.push_back(_x);
        
        // The y position
        double _y = currNode->pivCoor[1] + (0.5 + locParIndex[1])*_parSize;
        this->newPar->y.push_back(_y);
        
        // The z position
        if (DIM > 2) {
            double _z = currNode->pivCoor[2] + (0.5 + locParIndex[2])*_parSize;
            this->newPar->z.push_back(_z);
        }

        // Assign other data
        this->newPar->nodeID.push_back(currNode->nodeID);  // Must be targeting the real node ID
        this->newPar->level.push_back(currNode->level);    // Also the real level
        this->newPar->s.push_back(_headSize);    // Size using the old size
        this->newPar->isActive.push_back(currNode->isActive);    // Flag will be set as true (because if divided it means an active block)
        this->newPar->neighbor.push_back({});    // Empty neighbor list
        this->newPar->num++;     // Add the number of the particle
    }
    return;
}

// /**
//  *  @brief  Generate the particle in the given node. 
//  *  NOTE: The size of particle is using the size of particle in the head node.
//  *         
//  *  @param  currNode  The node where the particle are to be generated.
//  *  @param  headSize  The size of particle in the head node.
//  *  @param  activeSign  The active sign for the generated particle.
// */
// void GridNodeAdapt::generateParticleActive(Node *currNode, double _headSize, bool _isActive){
//     // Generate the particle in the current Node

//     // Update the current node data
//     currNode->parList.clear();  // Reserve the particle ID container on the current node 
//     currNode->headNodeID = -1;  // Turn off the head node ID flag

//     // Size of the particle inside the node
//     double _parSize = currNode->length / this->parNum;
//     // double _headSize = headNode->length / this->parNum;

//     // Generate all particle inside the current Node
//     int locParIndex[DIM];     // The particle local index (taken from the node)
//     for (int locParID = 0; locParID < this->totParNum; locParID++){
//         // Put the current particle ID into the node
//         currNode->parList.push_back(this->newPar->num);    // Push the last pushed ID inside the parList

//         // The local index coordinate inside the node
//         basis_loop(d) locParIndex[d] = (locParID/this->div[d]) % this->parNum;

//         // Assign the particle position
//         // The x position
//         double _x = currNode->pivCoor[0] + (0.5 + locParIndex[0])*_parSize;
//         this->newPar->x.push_back(_x);
        
//         // The y position
//         double _y = currNode->pivCoor[1] + (0.5 + locParIndex[1])*_parSize;
//         this->newPar->y.push_back(_y);
        
//         // The z position
//         if (DIM > 2) {
//             double _z = currNode->pivCoor[2] + (0.5 + locParIndex[2])*_parSize;
//             this->newPar->z.push_back(_z);
//         }

//         // Assign other data
//         this->newPar->nodeID.push_back(currNode->nodeID);   // Must be targeting the real node ID
//         this->newPar->level.push_back(currNode->level);     // Also the real level
//         this->newPar->s.push_back(_headSize);           // Size using the old size
//         this->newPar->isActive.push_back(_isActive);    // Flag will be set as true (because if divided it means an active block)
//         this->newPar->neighbor.push_back({});           // Empty neighbor list
//         this->newPar->num++;     // Add the number of the particle
//     }
//     return;
// }

// #pragma endregion
