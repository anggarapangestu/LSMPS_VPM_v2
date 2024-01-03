#ifndef INCLUDED_TREE_CELL
#define INCLUDED_TREE_CELL

#include "../../Utils.hpp"
#include <unordered_map>

// Data lifetime
// > FMM data  -> Single calculation, single time use, clear all intermediate variable at finish
// > Tree data -> Full life time (Contain ngh calculate, adaptation, and initialization, hierarchcal)

/**
 *  @brief A single cell data container for FMM calculation.
 *  
 *  @headerfile treeCell.hpp
*/
struct Cell{
    // Cell basic information
    int ID;                 // Cell ID (sorted by hierarchy and index order)
    int level;              // Cell resolution level (tree hierarchy level)
    double length;          // Cell side length
    int index[DIM];         // Local index in current level resolution
    double centerPos[DIM];  // Cell center coordinate position

    // Cell for FMM information
    int parNum;                  // Number of all particle inside cell
    int srcParNum;               // Number of source particle inside cell
    std::vector<int> parIDList;     // List of partilce ID inside cell (*only available for leaf cell)

    // Cell flag variable
    bool isLeaf;        // A leaf particle flag := Cell with no child
    bool isActive;      // An active cell flag  := Contains source particle
    // Note: a non active cell may be a:
    //   > Cell lies outside the domain
    //   > Cell with no source particle inside it
    //   > Cell with no particle data (parNum = 0)

    // Default node constructor
    Cell():
        ID(0),
        level(0),
        length(0.0),
        parNum(0),
        srcParNum(0),
        parIDList(0),
        isLeaf(false),
        isActive(false)
    {
        // Assign the array
        basis_loop(d) index[d] = 0;
        basis_loop(d) centerPos[d] = 0.0;
    };

    // Default node constructor
    Cell(int _ID, int _level, double _length, 
         const std::vector<int> &_index, 
         const double _pivCoord[DIM]):
        ID(_ID),
        level(_level),
        length(_length),
        parNum(0),
        srcParNum(0),
        parIDList(0),
        isLeaf(false),
        isActive(false)
    {
        // Assign the array
        basis_loop(d) index[d] = _index[d];
        basis_loop(d) centerPos[d] = _pivCoord[d] + ((0.5 + _index[d]) * _length);
    };
    
    // Default deconstructor
    ~Cell(){
        parIDList.clear();
    };
};

/**
 *  @brief A cell tree structure for FMM calculation. Included with cell accessor, modifier,
 *  and other functional function for FMM calculation.
 *  
 *  @headerfile treeCell.hpp
*/
struct TreeCell{
// private:
    // =====================================
    // ========== INTERNAL METHOD ==========
    // =====================================
    // Index transformation
    int index2ID(const std::vector<int> &_index, const int &_level) const;      // Convert the cell index to hierarchy order
    std::vector<int> ID2index(const int &_ID, const int &_level) const;         // Convert the cell ID to cell index

    int coord2ID(const std::vector<double> &_coord, const int &_level) const;   // Find the corresponding cell ID that contain the given physical coordinate
    std::vector<int> coord2index(const std::vector<double> &_coord, const int &_level) const;  // Find the corresponding cell index that contain the given physical coordinate

    // Utility data and function
    double get_cellSize(const int &_level) const;   // Get the size of cell at givenlevel
    long long int get_startID(const int &_level) const;       // Get the cell ID starting at the given level
    int get_cellCount(const int &_level) const;     // Get the number of cell at current level (the number at one direction only)
    int get_level(const long long int &_ID) const;            // Get the cell level from ID


// public:
    /** Illustration index and ID using 2 dimensional space (can be expanded to 3 dimensional space)
     * 
     *     idx_y  ___________________  idx_y  ___________________   idx_y  ___________________ 
     *           |                   |       |         |         |      3 | 17 | 18 | 19 | 20 |
     *           |                   |     1 |  ID 3   |  ID 4   |        |____|____|____|____|
     *           |                   |       |         |         |      2 | 13 | 14 | 15 | 16 |
     *         0 |       ID 0        |       |_________|_________|        |____|____|____|____|
     *           |                   |       |         |         |      1 | 9  | 10 | 11 | 12 |
     *           |                   |     0 |  ID 1   |  ID 2   |        |____|____|____|____|
     *           |                   |       |         |         |      0 | 5  | 6  | 7  | 8  |
     *           |___________________|       |_________|_________|        |____|____|____|____|
     *     idx_x          0                       0         1               0    1    2    3 
     *                At level 0                  At level 1                   At level 2
    */


    // =====================================
    // ============ STRUCT DATA ============
    // =====================================

    /**
     *  @brief  Mapping of the Cell by its corresponding hierarchy order <_ID, _Cell>.
     *  @tparam _ID     The hierarchy order of the Cell.
     *  @tparam _Cell   The address of the Cell.
     *  NOTE:
     *      > The cell ID is set the same with tree hierarchy order and followed by the cell index order 
     *      > The initial tree having no
     *      > There may be existed a cell as child of leaf that is not active
    */
    std::unordered_map<int, Cell*> treeData;

    // Basic data
    int min_level;      // The leaf cell minimum level
    int max_level;      // The leaf cell maximum level
    const int chdNum;   // Number of child in one cell
    const int nghNum;   // Number of candidate neighbor that adjacent to the cell
    const int max_particle;     // Maximum particle inside the cell
    int cellCount;              // Number cell throughout the tree
    double rootLength;          // [Constant upon set up] The cell side length of the root cell
    double pivotCoord[DIM];     // [Constant upon set up] The tree pivot coordinate (global minimum coordinate)

    // Additional for FMM calculation
    std::vector<int> leafList;      // List of all leaf cell in tree


    // =====================================
    // =========== STRUCT METHOD ===========
    // =====================================
    // Data modification method
    void initializeTree(const std::vector<std::vector<double>> &parPos, const std::vector<bool> &activeFlag);
    void updateTree(const std::vector<std::vector<double>> &parPos, const std::vector<bool> &activeFlag);

    // Neighbor function
    void findNghLvl(const Cell *currCell, std::vector<int> &nghIDList) const;
    void intList(const int &currID, std::vector<int>& ID_list_1, std::vector<int>& ID_list_3) const;
    void extList(const int &currID, std::vector<int>& ID_list_2, std::vector<int>& ID_list_4) const;

    // Saving data function
    void saveTree(const TreeCell &tree, std::string name) const;
    void saveLeafTree(const TreeCell &tree, std::string name) const;
    void saveSelTree(const TreeCell &tree, std::string name, std::vector<int> &IDList) const;

    // Hierarchy function
    int get_parentID(const int &_ID) const;               // Get the parent ID of the given cell ID
    int get_parentID(const int &_ID, const int &_level) const;               // Get the parent ID of the given cell ID
    std::vector<int> get_childID(const int &_ID) const;   // Get the child ID of the given cell ID

    // Default constructor
    TreeCell() : 
        min_level(2),       // The minimal level that counts far field calculation
        max_level(2),       // The initial max level set same as minimal level
        chdNum(Pars::intPow(2,DIM)),    // The number of child is 2^DIM
        nghNum(Pars::intPow(3,DIM)),    // The number of neighbor candidate is 3^DIM
        max_particle(Pars::n_max),      // Maximum number of particle inside the cell
        cellCount(0),       // Initialize cell count in the tree to zero
        rootLength(0.0),    // The size of root cell
        leafList(0)         // Initialize the leaf list to zero
    {
        // Nothing to do here !
    }

    // Default destructor
    ~TreeCell(){
        // Free all the cell pointer data
        for (const auto &[ID, cell] : treeData) delete cell;
    }
};


class treeCell{
private:
    // Transformation map : in determining the neighboring cell
    std::vector<std::vector<int>> cellIdx;              // The corresponding transformation matrix index
    std::vector<std::vector<std::vector<int>>> T_map2D; // The transformation matrix cell IDs for 2D
    std::vector<std::vector<std::vector<std::vector<int>>>> T_map3D; // The transformation matrix cell IDs for 3D
    
    void createTmap2D();  // Creating the Tree map for 2D
    void createTmap3D();  // Creating the Tree map for 3D

public:
    // TREE BASIC PARAMETERs
    double baseLen;     // The length of the root cell
    double expansion;   // The tree basis domain expansion factor
    int max_level;      // The maximum level of the tree
    int min_level;      // The minimum level of leaf cell
    int max_par;        // The maximum number of source particle in cell
    int basis;          // The number of cell division
    int cellNum;        // The total number of the cell in the tree

    // CELL DATA VARIABLE : with ID in sequence of the parent - child
    std::vector<std::vector<double>> cellPos;   // The center position of the cell
    std::vector<bool> leafCellMark;             // The marking for leaf cell or not
    std::vector<bool> outsideCell;              // The marking whether the cell is existed but outside the domain
    std::vector<int> level;                     // The cell level hierarchy of the cell
    std::vector<int> parNum;                    // The number of all particel inside the cell
    std::vector<int> srcNum;                    // The number of source particel inside the cell
    std::vector<std::vector<int>> parIDList;    // List of all particle ID in the current cell

    // Additional for FMM calculation
    std::vector<int> leafCellList;              // List of all leaf cell

    // Tree hierarchy
    int findParentID(int currID) const;
    std::vector<int> findChildID(int currID) const;

    // Basic External method
    void initializeTree(std::vector<std::vector<double>>& parPos);
    void createTree(std::vector<std::vector<double>>& parPos, std::vector<bool> & activeMark);

    void updateCell(std::vector<std::vector<double>>& parPos, std::vector<bool> & activeMark);

    /* Find the neighbor list
       4 type neighbor cell
       * Type 1 : Touched neighbor, leaf cell
       * Type 2 : Well separated neighbor (current parent touched neighbor parent cell), at the same level to current cell
       * Type 3 : Well separated neighbor (current cell touched parent), leaf cell
       * Type 4 : Well separated neighbor (current parent touched neighbor parent cell), at low level to current cell
    */

    void intList(int currID, std::vector<int>& ID_list_1, std::vector<int>& ID_list_3) const;
    void intList_3d(int currID, std::vector<int>& ID_list_1, std::vector<int>& ID_list_3) const;
    void extList(int currID, std::vector<int>& ID_list_2, std::vector<int>& ID_list_4) const;
    void extList_3d(int currID, std::vector<int>& ID_list_2, std::vector<int>& ID_list_4) const;
    void findLevel(const std::vector<int>& ID_list, std::vector<int>& level) const;
    
    // Save data
    void saveTreeCell();
    void saveTreeCell(std::string name);
};

#endif