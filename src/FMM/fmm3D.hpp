#ifndef FMM_3D_PACKAGE
#define FMM_3D_PACKAGE

#include "../../global.hpp"
#include "treeCell.hpp"

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
    std::vector<double> parField_x;     // The result data field in x direction (or velocity in x direction)
    std::vector<double> parField_y;     // The result data field in y direction (or velocity in y direction)
    std::vector<double> parField_z;     // The result data field in z direction (or velocity in z direction)
    
    // FMM Source Expansion Value (The particle source sum; indexing -> [cellID, expansion order])
    std::unordered_map<int, std::vector<double>> mp;  // Local particle source sum, origin at current cell center

    // FMM Source Expansion Value (The particle source sum; indexing -> [cellID, expansion order])
    std::vector<int> G2LcellID;              // Cell ID converter from global to local
    std::vector<std::vector<double>> mp_Vx;  // Local particle source sum of vorticity in x, origin at current cell center
    std::vector<std::vector<double>> mp_Vy;  // Local particle source sum of vorticity in y, origin at current cell center
    std::vector<std::vector<double>> mp_Vz;  // Local particle source sum of vorticity in z, origin at current cell center

    // std::unordered_map<int, std::vector<double>> mp_Vx;  // Local particle source sum of vorticity in x, origin at current cell center
    // std::unordered_map<int, std::vector<double>> mp_Vy;  // Local particle source sum of vorticity in y, origin at current cell center
    // std::unordered_map<int, std::vector<double>> mp_Vz;  // Local particle source sum of vorticity in z, origin at current cell center

    // Internal method
    int binComb(int n, int k) const;  // Function to calculate the binomial combinatoric

    // FMM direct calculation tools
    void potDirSum(int _cellID, const TreeCell &_cellTree,      // Function to calculate the direct sum potential
                   const std::vector<int> &_nghCell, 
                   const std::vector<std::vector<double>> &_parPos, 
                   const std::vector<double> &_srcVal);    
    void fieldDirSum(int _cellID, const TreeCell &_cellTree,    // Function to calculate the direct sum field
                     const std::vector<int> &_nghCell, 
                     const std::vector<std::vector<double>> &_parPos,
                     const std::vector<double> &_srcVal);

    void velocityDirSum(int _cellID, const TreeCell &_cellTree, // Function to calculate the direct sum veclocity field
                        const std::vector<int> &_nghCell, 
                        const std::vector<std::vector<double>> &_parPos,
                        const std::vector<double> &_alphaX,
                        const std::vector<double> &_alphaY,
                        const std::vector<double> &_alphaZ);
    
    // The fundamental sequence of FMM translation calculation
    void setupFMM(const TreeCell &_cellTree, 
                  const std::vector<std::vector<double>> &_parPos, 
                  const std::vector<bool> &_activeFlag, 
                  const std::vector<double> &_srcVal);

    // The fundamental sequence of FMM translation calculation
    void setupVelocityFMM(const TreeCell &_cellTree, 
                          const std::vector<std::vector<double>> &_parPos, 
                          const std::vector<bool> &_activeFlag, 
                          const std::vector<double> &_alphaX,
                          const std::vector<double> &_alphaY,
                          const std::vector<double> &_alphaZ);

public:
    // The FMM calculation for potential [Not Available]
    void calcPotential(const TreeCell &_cellTree, 
                       const std::vector<std::vector<double>> &_parPos, 
                       const std::vector<bool> &_activeFlag, 
                       const std::vector<double> &_srcVal);
    
    // The FMM calculation for field (Calculate the velocity per direction)
    void calcField(const TreeCell &_cellTree, 
                   const std::vector<std::vector<double>> &_parPos, 
                   const std::vector<bool> &_activeFlag, 
                   const std::vector<double> &_srcVal);

    // The FMM calculation for velocity field (Calculate all particle)
    void calcVelocity(const TreeCell &_cellTree, 
                      const std::vector<std::vector<double>> &_parPos, 
                      const std::vector<bool> &_activeFlag, 
                      const std::vector<double> &_alphaX,
                      const std::vector<double> &_alphaY,
                      const std::vector<double> &_alphaZ);

    // The FMM calculation for velocity field (Calculate a target only)
    void calcVelocityNearBody(const TreeCell &_cellTree, 
                              const std::vector<std::vector<double>> &_parPos, 
                              const std::vector<bool> &_activeFlag, 
                              const std::vector<bool> &_targetFlag, 
                              const std::vector<double> &_alphaX,
                              const std::vector<double> &_alphaY,
                              const std::vector<double> &_alphaZ);

    // The FMM old code refactor
    void calcVelocityFast(const std::vector<std::vector<double>> &_parPos, 
                          const std::vector<bool> &_activeFlag,
                          const std::vector<bool> &_targetFlag,
                          const std::vector<double> &_alphaX,
                          const std::vector<double> &_alphaY,
                          const std::vector<double> &_alphaZ);

    
    // The FMM get function from the result of calculation

    void P2M_calc(std::vector<double> &multipole, 
                  const double &src,
                  const double &dx, 
                  const double &dy, 
                  const double &dz) const;
    void M2M_calc(std::vector<double> &_mParent, 
                  std::vector<double> &_mChild, 
                  const double &dx, 
                  const double &dy, 
                  const double &dz) const;
    void M2P_mul_calc(std::vector<double> &diff_x, 
                      std::vector<double> &diff_y, 
                      std::vector<double> &diff_z, 
                      const double &R2, 
                      const double &dx, 
                      const double &dy, 
                      const double &dz) const;

    // The FMM get function from the result of calculation
    std::vector<double> get_Potential() const;
    std::vector<double> get_Field_x() const;
    std::vector<double> get_Field_y() const;
    std::vector<double> get_Field_z() const;
};

#endif