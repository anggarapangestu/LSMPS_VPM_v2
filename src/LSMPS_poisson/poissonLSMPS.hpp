#ifndef INCLUDED_LSMPS_POISSON
#define INCLUDED_LSMPS_POISSON

#include "../../Utils.hpp"
#include "../Eigen/Dense"
#include "../Eigen/Sparse"

/**
 *  @brief A solver to poisson problem [Δφ(x)=-σ(x), φ:potential; σ:source] based on 
 *  LSMPS differential. The LSMPS discretization type is set from B standard derivation.
 *  NOTE: 
 *      > The LSMPS formula of Dx(f) = H * Minv * b is calculated by seperating 
 *         each neighbor function data in vector 'b'
 *      > The value of second derivative is put into the poisson equation
 *      > The function of each particle at 'b' is rearrange to make global matrix
 *      > Solve the global matrix
 *  
 *  @headerfile	poissonLSMPS.hpp
 */
class PoissonLSMPS
{
private:
    // CLASS MEMBER DATA
    // =================
    int N;                      // Global matrix base size (number of particle)
    std::vector<double> phi;    // The solution to poisson problem (size of N)
    Eigen::SparseMatrix<double> A_sparse;   // The global sparse matrix (A square matrix size of N^2)
    std::vector<Eigen::Triplet<double>> tripletList;    // The triplet list for sparse matrix constructor
    Eigen::MatrixXd A_dense;                // The global dense matrix (A square matrix size of N^2)

    // CLASS METHOD
    // ============
    double weight_function(const double &_rij, const double &_Re);
    std::vector<double> get_p(const double &xij, const double &yij, const double &rs);
    
    void matrix_const_LSMPS_A(const int &ID_i,
                              const std::vector<std::vector<double>> &coords,
                              const std::vector<int> &nghList,
                              const std::vector<double> &size);
    void matrix_const_LSMPS_B(const int &ID_i,
                              const std::vector<std::vector<double>> &coords,
                              const std::vector<int> &nghList,
                              const std::vector<double> &size);
    void matrix_const_neumann_A(const int &ID_i,
                              const std::vector<std::vector<double>> &coords,
                              const std::vector<int> &nghList,
                              const std::vector<double> &size,
                              const int& type);
    void matrix_const_neumann_B(const int &ID_i,
                              const std::vector<std::vector<double>> &coords,
                              const std::vector<int> &nghList,
                              const std::vector<double> &size,
                              const int& type);

public:
    // Parameter for Multiresolution 
    // =============================
    Particle patchPar;          // The patched particle data
    std::vector<int> nIFaceParID;   // List of particle ID near resolution interface
    std::vector<bool> nIFaceParFlag;   // List of particle ID near resolution interface
    std::vector<int> nIFaceNodeID;  // List of node ID near resolution interface

    std::unordered_map<int, std::vector<int>> nIFaceNodeParMap;   // Map an interface node with generate particle

    std::vector<int> nPatchNodeID;  // List of node ID of patched particle
    std::vector<std::vector<int>> fullNeighborList; // The full neighbor list

    // Matrix generator
    void create_global_matrix(const std::vector<std::vector<double>> &coords,
               const std::vector<std::vector<int>> &nghList,
               const std::vector<double> &size,
               const std::vector<int> &boundaryLoc);
               // const std::vector<bool> &isBoundary);

    // Solver
    void solve(const std::vector<double> &RHS);

    // The get function from the calculation result
    std::vector<double> get_phi();

    PoissonLSMPS(/* args */);
    ~PoissonLSMPS();
};


#endif