#ifndef INCLUDED_CELL_LIST
#define INCLUDED_CELL_LIST

#include "../../../Utils.hpp"

// The class to store the particle inside the cell
class CellList
{
/*
NOTE:
    > The first level ID of each vector is the ID of basis cell
    > The second level ID of each vector is the current child cell (inside the basis cell)
    > The third level ID of each vector is the value of the corresponding vector
    > The public method have:
      * cell list initialization,
      * neighbor list search, and
      * adaptive cell list update.
*/
private:
    // ======= Cell List Parameters (BASIC) ======= //

    // Number of cell in Cell List
    // int nx = 1;     // In x direction
	// int ny = 1;     // In y direction
	// int nz = 1;     // In z direction
	std::vector<int> n_basis;   // Number of basis cell at each direction
    int num;                    // Total cell

    // Basic Parameters
    double size;    // The basis cell size (level 0 cell or root)
    int cellNum;    // The total number of cell in tree hierarchy
    const int basis = std::pow(2,DIM);        // Number of division child cell for each DnC
    const int basisNGH = std::pow(3,DIM);     // The basis used for neighbor evaluation
	std::vector<double> center_pos;	  // position of the particle domain center
    std::vector<double> pivot_pos;	  // position of the cell pivot (bottom, left, back)

    // Cell data
    std::vector<std::vector<int>> basis_neighbor;   // The cell neigbor list of basis cell (ALREADY NOT USED)
    std::vector<std::vector<int>> cell_flag;        // Cell flag of the current for particle
                                                    //The flag describes if the current cell is contained by particle
                                                    //   * -1: the particle lies at higher level or child (finer cell)
                                                    //   * 0 : the particles lies at this cell
                                                    //   * 1 : the particle lies at lower level or parent (larger cell)
    std::vector<std::vector<std::vector<double>>> mid_pos;    // The cell mid possition in x,y,z
    std::vector<std::vector<std::vector<int>>> CH_parID;      // The particle ID inside each child

    // Look up table
    std::vector<std::vector<int>> index;    // The T-map index position of each cell ID
    std::vector<int> lvl;                   // The corresponding level of cell ID
    std::vector<std::vector<std::vector<int>>> T_map2D;                 // The Transformation map Cell ID in 2D
    std::vector<std::vector<std::vector<std::vector<int>>>> T_map3D;    // The Transformation map Cell ID in 3D
    
    
    // ======= Cell List Parameters (ADAPTIVE ALGH) ======= //
    // The variable for adaptive particle
    std::vector<std::vector<bool>> levelUp;                   // Sign to level up all particle at the current cell
    // std::vector<std::vector<int>> targetLevel;                // The level target for upgrade at the current cell
    std::vector<std::vector<int>> tempBasisCellID;            // The list of basis-cell pair for level upgrade
    // std::vector<bool> splitParIDFlag;                         // The flag of splitted particle lies at the same index ID
    
    // The variable for particle redistribution 
    std::vector<std::vector<int>> cell_flag_new;              // Cell flag for new particle 
    std::vector<std::vector<std::vector<int>>> CH_new_parID;  // The new particle ID inside each child
    
    /* ***************************** //
    >>>> Cell List Illustration <<<<
     ___________________________________
    |     |     |           |           |
    |_____|_____|           |           |
    |     |__|__|           |           |
    |_____|__|__|___________|___________|
    |     |Lv.1 |   basis   |__|__|__|__|
    |_____|_____|   cell    |__|__|__|__|
    |     |     | (level 0) |     |     |
    |_____|_____|___________|_____|_____|
    
    Note:
    > The illustraion above shows a cell list for 2D domain
    > There are 6 basis cell (6 biggest cell)
    > The maximum level is 2 (the smallest cell)
    // ***************************** */

    // +----------------------------+ //
    //  ====== Private Method ======  //
    // +----------------------------+ //

    // ======= Tools & Utilities ======= //
    std::vector<int> findChild(int PAR_ID);     // Finding the Child Cell IDs from a given Parent Cell ID
	int findParent(int CHD_ID);                 // Finding the sub Parent Cell IDs from a given Child Cell ID (Note tah)
    // int cellLevel(int ID);                      // Returning the level of Cell IDs
    void totalCellNumber(int max_level);        // Returning the number of cell in the hierarchy tree
    
    // Devide and Conqueror
    void divideCell(int BASIS_ID, int PAR_ID, Particle & particle, 
                    std::vector<std::vector<std::vector<int>>> & SRC_CH_parID, 
                    std::vector<std::vector<int>> & SRC_cell_flag);
    // Returning the ID's of neighboring Cells
    void neighborCell(std::vector<std::vector<int>>& ngh_PAIR_ID, int basis_ID, int cell_ID);
    
    // ======= Adaptive Method ======= //
    void adtParSplit(Particle & particle, int tgt_lvl, std::vector<double>& PARsize);     // Update adaptive particle -> Particle splitting
    void adtCellDivide(Particle & particle);                // Update adaptive cell -> cell division

public:
    // Tools for initialization
    void initCellList(const Particle & particle);
    void createCellList(Particle & particle);
    void findNeighbor(Particle & particle);
    
    // Tools for particle redistribution
    void createCellList_new(Particle & particle);
    std::vector<std::vector<int>> findNeighbor_new(const Particle & parEval, const Particle & parBase);
    
    // Tools for adaptive particle
    void setLevelUp(int basis_ID, int cell_ID, int type);   // Set the levelUp flag at each cell (type 0 for clearing variable)
    void performAdaptive(Particle & particle, int updateStage, int tgt_lvl, std::vector<double>& PARsize);// Perform the adaptive algorithm
    void saveCellData();                                    // Perform the adaptive algorithm
    
    // Debugging and displaying
    // void showBasisNgh();
    void checkNGH(const Particle & par);
};

#endif