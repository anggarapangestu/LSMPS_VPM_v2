#ifndef INCLUDED_GRID_NODE
#define INCLUDED_GRID_NODE

#include "../../Utils.hpp"
#include <unordered_map>

#define ROOT_LEVEL 0    // Level of ROOT node


/**
 *  @brief A particle container that is grouped based on a grid
 *  structure that associates a multiblock tree hierarchy. A tree
 *  hierarchy groups the nodes by parent-children relation from the
 *  ROOT level to the LEAF level.
 *  
 *  There are two coordinate used in this container:
 *  [1] Spatial coordinate -> Denote the physical coordinate in the domain, 
 *  [2] Grid map coordinate -> Denote the indexing position in grid Map.
 *
 *  @headerfile gridNode.hpp
 *  Node detail
    > Store the basic node data, position, index, and
    > Provide the connection of parential hierarchy
    > Having the adaptive resolution flag
 */
struct Node
{
    // Node Identifier
    int nodeID;             // ID of the node
    int level;              // Level of the node (ROOT at level 0; LEAF at MAX_LEVEL)
    int index[DIM];         // Index position of node (A grid map indexing at each dimension basis)
    double pivCoor[DIM];    // Pivot coordinate of node (Physical coordiante)
    double length;          // Length size of node
    bool isLeaf;            // Flag for leaf node

    // Member for adaptive performance
    int headNodeID;     // The ID pointing to the head node (-1) in value if not used yet. [For refinement only]
    int tarResLv;       // Resolution target to evaluate the adaptive operation
    
    // Adaptive Node Modifier <?> This become unnecessary <?>
    bool isActive;          // Flag for an active cell
    bool needCompression;   // Flag for node compression    <!> USED <!>
    bool needRefinement;    // Flag for node refinement     <?> STILL NOT USED <?>
    
    // Container of particle ID (Only leaf contain the particle)
    std::vector<int> parList;   // List of particle ID inside the node (Only LEAF contains particle)

    // Basic assignment constructor
    Node(int nodeID_, int level_, double length_, int index_[]):
        nodeID(nodeID_),
        level(level_),
        length(length_),
        isLeaf(true),
        headNodeID(-1),
        tarResLv(ROOT_LEVEL),
        isActive(false),
        needCompression(false),
        needRefinement(false)
    {
        // Assign the array
        basis_loop(d) index[d] = index_[d];
        basis_loop(d) pivCoor[d] = 0.0;
        parList.clear();
    }

    // Copy node constructor
    Node(Node *_node):
        nodeID(_node->nodeID),
        level(_node->level),
        length(_node->length),
        isLeaf(_node->isLeaf),
        headNodeID(_node->headNodeID),
        tarResLv(_node->tarResLv),
        isActive(_node->isActive),
        needCompression(false),
        needRefinement(false)
    {
        // Assign the array
        basis_loop(d){
            index[d] = _node->index[d];
            pivCoor[d] = _node->pivCoor[d];
        }
        
        // Assign the vector
        this->parList = _node->parList;
    }
    
    // Default node constructor
    Node():
        nodeID(0),
        level(0),
        length(0.0),
        isLeaf(true),
        headNodeID(-1),
        tarResLv(ROOT_LEVEL),
        isActive(false),
        needCompression(false),
        needRefinement(false)
    {
        // Assign the array
        basis_loop(d) index[d] = 0;
        basis_loop(d) pivCoor[d] = 0.0;
        parList.clear();
    };
    
    // Default deconstructor
    ~Node(){
        parList.clear();
    };
};


/**
 *  @brief A container that groups the node in the domain
 *  that is associated with tree hierarchy methods. This 
 *  container may be called a 'NODE MANAGER'.
 *
 *  @headerfile gridNode.hpp
 *  GridNode detail
    > Store all of the node data
    > Provide the connection of parential hierarchy
 */
struct GridNode
{
    /**
     *  @brief  Mapping of the Node by its corresponding ID <_ID, _Node>.
     *  @tparam _ID     The ID of the Node.
     *  @tparam _Node   The address of the Node.
    */
    std::unordered_map<int, Node*> nodeMap;
    
    // GridNode element member
    double pivotCoor[DIM];  // The pivot position coordinate (left(x), bottom(y), front (z))
    int gridCount[DIM];     // The ROOT node count at each dimension
    int baseParNum;         // Number of particle in one node (each dimension)
    double gridSize;        // Size of the ROOT (level 0) grid node box
    int rootNodeNum;        // The number of all ROOT node (Basically product each @gridCount element)
    int maxLevel;           // The limit of resolution step or LEAF node level
    int chdNum;             // The number of child on a node
    std::vector<int> startID;   // List of starting ID at each level

    // UTILITIES FUNCTION METHOD
    
    int getPivID(int level) const;
    int getLevel(int ID) const;

    // INDEX TRANSFORMATION METHOD

    void pos2idx(int index[DIM], const double coorPnt[DIM], int level) const;
    void ID2Index(int index[DIM], int ID, int level) const;
    int idx2ID(const int index[DIM], int level) const;
    int pos2ID(const double coorPnt[DIM], int level) const;
    
    // TREE HIERARCHY METHOD

    int findParent(const Node *currNode) const;
    std::vector<int> findChild(const Node *currNode) const;
    void findSibling(std::vector<Node*> &sibList, const Node *currNode) const;
    // *Function Overloader
    int findParent(int ID) const;
    std::vector<int> findChild(int ID) const;

    // NODE NEIGHBOR EVALUATION METHOD

    void findNghLvl(std::vector<int> &nghID, const Node *currNode) const;
    void findNghAll(std::vector<Node*> &nghID, const Node *currNode) const;

    // ADAPTATION METHOD

    int refineNode(Node *currNode, const GridNode &tools, std::vector<int> &chdIDlist);
    int divideParticle(Node *currNode, const Particle &par, const std::vector<int> &chdIDlist, int type);

    // DATA SAVING METHOD

    void saveGrid(const GridNode& nodeList, std::string name) const;
    void saveLeafGrid(const GridNode& nodeList, std::string name) const;
    void saveSelectedGrid(const GridNode& nodeList, const std::vector<int> &ID, std::string name) const;
    
    // NOT USED METHOD
    // void findNghAdj(std::vector<Node*> &nghID, const Node *currNode) const;
    // int refineNode(Node *currNode, std::vector<int> &chdIDlist);
    // int compressNode(Node *currNode, int &parID);
    // void saveLeafGrid(GridNode& nodeList) const;

    // Grid Node Default Constructor
    GridNode(): baseParNum(0), gridSize(0.0), rootNodeNum(0) {
        // Calculate the predefined member
        this->maxLevel = Pars::max_level;
        this->chdNum = Pars::intPow(2,DIM);

        // Initialize the dynamic container
        this->nodeMap.clear();
        this->startID.clear();

        // Initialize the array
        basis_loop(d) gridCount[d] = 0;
        basis_loop(d) pivotCoor[d] = 0.0;
    };

    // Grid Node Default Deconstructor
    ~GridNode(){
        // Free the pointer member
        for (auto &[key,val] : this->nodeMap) delete val;
    };

};

#endif
