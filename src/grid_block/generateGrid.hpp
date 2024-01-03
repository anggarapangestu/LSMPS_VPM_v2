#ifndef INCLUDE_GRID_GENERATION
#define INCLUDE_GRID_GENERATION

#include "../../Utils.hpp"
#include "../geometry/geometry.hpp"
#include "gridNode.hpp"

/**
 *  @brief  This class is used to generate the grid in the simulation
 *  domain. The class will be used inside the initialization procedure.
 * 
 *  @headerfile generateGrid.hpp
*/
class generateGrid
{
private:
    // The basic parameter in this class
    double minCoor[DIM];    // Physical domain spatial coordinate MINIMUM position
    double maxCoor[DIM];    // Physical domain spatial coordinate MAXIMUM position
    double length[DIM];     // Physical domain length size in each dimension basis
    double leafBlockSize;   // Size of the leaf block (finest)
    double rootBlockSize;   // Size of the root block (largest)
    int maxLevel;       // Grid block resolution step limit
    int baseParNum;     // Number of particle in node (1D)
    int NghlevelDiff;   // The maximum different level between neighboring node

    // Member method

    void updateGridNode(GridNode& nodeList);
    void generateRootRec(GridNode& nodeList, int dim, int &ID, int index[DIM]) const;
    void generateRootDir(GridNode& nodeList) const;
    
public:
    // Public method

    void nodeGeneration(GridNode& nodeList, const std::vector<Body>& bodyList);
    void createNode(GridNode& nodeList, Particle &particle);
    void setNghLvlDiff(int diff);
    
    // Constructor
    generateGrid();

    // Deconstructor
    ~generateGrid();
};

#endif