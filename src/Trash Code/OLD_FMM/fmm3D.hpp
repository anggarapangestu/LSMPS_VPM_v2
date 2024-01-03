#ifndef FMM_3D_PACKAGE
#define FMM_3D_PACKAGE

#include "../../global.hpp"
#include "treeCell.hpp"
#include <complex>

/**
 *  @brief A fast multipole method (FMM) subroutine to enhance the 
 *  special poisson solution. A tree cell data structure is used 
 *  to manage the data grouping of particle.
 *
 *  @headerfile fmmm3D.hpp
 */
class fmm3D{
private:
    // Internal variable
    int parNum;         // Number of particle
    int expOrd;         // Maximum expansion order
    int cellNum;        // Total number of cell 
    int maxLevel;       // Highest tree level

    // The FMM result data of each particle
    std::vector<double> parPotential;   // The result data potential
    std::vector<double> parField_x;     // The result data field in x direction
    std::vector<double> parField_y;     // The result data field in y direction
    std::vector<double> parField_z;     // The result data field in z direction
    
    // FMM Source Expansion Value (The particle source sum; indexing -> [cellID, expansion order])
    std::vector<std::vector<double>> mp;  // Local particle source sum, origin at current cell center

    // Internal method
    int binComb(int n, int k);  // Function to calculate the binomial combinatoric

    // FMM direct calculation tools
    void potDirSum(int _cellID, const treeCell &_cellTree,      // Function to calculate the direct sum potential
                   const std::vector<int> &_nghCell, 
                   const std::vector<std::vector<double>> &_parPos, 
                   const std::vector<double> &_srcVal);    
    void fieldDirSum(int _cellID, const treeCell &_cellTree,    // Function to calculate the direct sum field
                     const std::vector<int> &_nghCell, 
                     const std::vector<std::vector<double>> &_parPos,
                     const std::vector<double> &_srcVal);  
    
    // The fundamental sequence of FMM translation calculation
    void setupFMM(const treeCell &_cellTree, 
                  const std::vector<std::vector<double>> &_parPos, 
                  const std::vector<bool> &_activeFlag, 
                  const std::vector<double> &_srcVal);

public:
    // The FMM calculation for potential
    void calcPotential(const treeCell &_cellTree, 
                       const std::vector<std::vector<double>> &_parPos, 
                       const std::vector<bool> &_activeFlag, 
                       const std::vector<double> &_srcVal);
    
    // The FMM calculation for field
    void calcField(const treeCell &_cellTree, 
                   const std::vector<std::vector<double>> &_parPos, 
                   const std::vector<bool> &_activeFlag, 
                   const std::vector<double> &_srcVal);

    // The FMM get function from the result of calculation
    std::vector<double> get_Potential();
    std::vector<double> get_Field_x();
    std::vector<double> get_Field_y();
    std::vector<double> get_Field_z();
};

#endif