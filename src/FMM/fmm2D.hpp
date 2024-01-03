#ifndef FMM_2D_PACKAGE
#define FMM_2D_PACKAGE

#include "../../global.hpp"
#include "treeCell.hpp"
#include <complex>

/**
 *  @brief A fast multipole method (FMM) subroutine to enhance the 
 *  special poisson solution. A tree cell data structure is used 
 *  to manage the data grouping of particle.
 *
 *  @headerfile fmmm2D.hpp
 */
class fmm2D{
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
    
    // FMM Multipole Expansion Value (The particle source sum; indexing -> [cellID, expansion order])
    std::unordered_map<int, std::vector<std::complex<double>>> ak;  // Multipole expansion, origin at current cell center
    std::unordered_map<int, std::vector<std::complex<double>>> bk;  // Local expansion, origin at current cell center

    // Internal method
    int binComb(int n, int k);  // Function to calculate the binomial combinatoric

    // FMM direct calculation tools
    void potDirSum(int _cellID, const TreeCell &_cellTree,      // Function to calculate the direct sum potential
                   const std::vector<int> &_nghCell, 
                   const std::vector<std::vector<double>> &_parPos, 
                   const std::vector<double> &_srcVal);    
    void fieldDirSum(int _cellID, const TreeCell &_cellTree,    // Function to calculate the direct sum field
                     const std::vector<int> &_nghCell, 
                     const std::vector<std::vector<double>> &_parPos,
                     const std::vector<double> &_srcVal);  
    
    // The fundamental sequence of FMM translation calculation (upward and downward pass)
    void setupFMM(const TreeCell &_cellTree, 
                  const std::vector<std::vector<double>> &_parPos, 
                  const std::vector<bool> &_activeFlag, 
                  const std::vector<double> &_srcVal);

public:
    // The FMM calculation for potential
    void calcPotential(const TreeCell &_cellTree, 
                       const std::vector<std::vector<double>> &_parPos, 
                       const std::vector<bool> &_activeFlag, 
                       const std::vector<double> &_srcVal);
    
    // The FMM calculation for field
    void calcField(const TreeCell &_cellTree, 
                   const std::vector<std::vector<double>> &_parPos, 
                   const std::vector<bool> &_activeFlag, 
                   const std::vector<double> &_srcVal);

    // The FMM get function from the result of calculation
    std::vector<double> get_Potential();
    std::vector<double> get_Field_x();
    std::vector<double> get_Field_y();
};

#endif