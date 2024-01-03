#ifndef INCLUDED_FAST_FMM_3D
#define INCLUDED_FAST_FMM_3D

#include "../../global.hpp"

/**
 *  @brief A tree cell structure for FMM calculation. 
 *  NOTE: Only made for 3D FMM velocity calculation.
 *  
 *  @headerfile fast3DFMM.hpp
*/
struct FMMCell
{
private:
    // Internal method to calculate tree data
    // Utility function
    int getOctant(const std::vector<std::vector<double>> &parPos, const int &parID, const int &cellID){
        int octant = (parPos[parID][0] > this->xc[cellID]) * 1 +
                     (parPos[parID][1] > this->yc[cellID]) * 2 +
                     (parPos[parID][2] > this->zc[cellID]) * 4;
        return octant;
    }
    
public:
    void generateTree(const std::vector<std::vector<double>> &parPos, const std::vector<bool> &activeFlag);
    int getCellLevel(const double &size){
        int level = std::round(std::log2(this->rootLength / size));
        return level;
    }

    // Fundamental parameter
    const int chdNum;               // Number of child cell (must be 8 for 3D)
    const int critNum;              // Maximum number of particle in cell
    double rootLength;              // The length of root cell
    int maxLevel;                   // The maximum cell resolution level in the tree
    
    // Cell list data
    int cellNum;                    // Number of all cell in tree data
	std::vector<int> nPar;          // Number of all particle in current cell
	std::vector<int> chdFlag;       // Child octant flag in current cell
	std::vector<int> parentID;      // Parent ID of current cell
	std::vector<double> xc;         // Cell center x coordinate
	std::vector<double> yc;         // Cell center y coordiante
	std::vector<double> zc;         // Cell center z coordiante
	std::vector<double> size;       // Cell box length size
    std::vector<bool> isActive;     // Cell active mark
    std::vector<int> startID;       // The starting ID of cell at given level

    // NOTE: child octant flag is written in binary 
    // e.g. binary of [00100110] -> [--6--32-] only have child in octant 2,3 and 6
	
    // Utility container
    std::vector<std::vector<int>> parID;    // List of the first 'critNum' particle ID in the cell
	std::vector<std::vector<int>> chdID;    // List of child cell ID at each octant

    // For further calculation
    void save_all_tree(FMMCell treeData, std::string name);
    void save_leaf_tree(FMMCell treeData, std::string name);

    FMMCell() : chdNum(8), critNum(Pars::src_count_max), cellNum(0)
    {
        // Nothing to do here
        // Note : Child number must be 8 for 3D domain
    }
};


/**
 *  @brief A class consisted of FMM calculation. 
 *  NOTE: Only made for 3D FMM velocity calculation.
 *  
 *  @headerfile fast3DFMM.hpp
*/
class fastFMM3Dutils
{
private:
    // The private class

    void multipoleCalc(const FMMCell &_treeData, 
                       const std::vector<std::vector<double>> &_parPos, 
                       const std::vector<bool> &_activeFlag,
                       const std::vector<double> &_alphaX,
                       const std::vector<double> &_alphaY,
                       const std::vector<double> &_alphaZ);
    void evaluateFMM(const FMMCell &_treeData, 
                     const std::vector<std::vector<double>> &_parPos, 
                     const std::vector<bool> &_activeFlag,
                     const std::vector<bool> &_targetFlag,
                     const std::vector<double> &_alphaX,
                     const std::vector<double> &_alphaY,
                     const std::vector<double> &_alphaZ);

    // Cell multipole (*Corresponding to the source that cell contains)
    std::vector<std::vector<double>> m_Vx;  // Multipole for source in x direction
	std::vector<std::vector<double>> m_Vy;  // Multipole for source in y direction
	std::vector<std::vector<double>> m_Vz;  // Multipole for source in z direction

public:
    void velocityCalc(const FMMCell &_treeData, 
                      const std::vector<std::vector<double>> &_parPos, 
                      const std::vector<bool> &_activeFlag,
                      const std::vector<bool> &_targetFlag,
                      const std::vector<double> &_alphaX,
                      const std::vector<double> &_alphaY,
                      const std::vector<double> &_alphaZ);

    // Fundamental parameter
    const int expOrd;               // The expansion order

    // Particle velocity
    std::vector<double> velX;       // The velocity in x direction
    std::vector<double> velY;       // The velocity in y direction
    std::vector<double> velZ;       // The velocity in z direction

    fastFMM3Dutils() : expOrd(10){
        // Nothing to do here
    };
};

#endif