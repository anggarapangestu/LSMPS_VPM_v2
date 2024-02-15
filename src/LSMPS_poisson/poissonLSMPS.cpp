#include "poissonLSMPS.hpp"

#include "../Eigen/Cholesky"
#include "../Eigen/SparseQR"

// #include "../LSMPS/LSMPSa.hpp"
// #include "../LSMPS/LSMPSb.hpp"

// For linear algebra operation
using namespace Eigen;

// Shorthand for vector declaration
template <typename U> using vec = std::vector<U>;

// ==================================
// +------ Parameter Constant ------+
// ==================================
#define BASE_SIZE_2D_A 5        // Size of base matrix 2D standard A
#define BASE_SIZE_2D_B 6        // Size of base matrix 2D standard B

#define BASE_SIZE_3D_A 9        // Size of base matrix 3D standard A
#define BASE_SIZE_3D_B 10       // Size of base matrix 3D standard B

#define RADIUS_FACTOR_2D 3.1    // Support radius factor toward particle size 3D
#define RADIUS_FACTOR_3D 2.3    // Support radius factor toward particle size 2D

// ==================================
// +------- Calculation Type -------+
// ==================================
/** The type of particle interation
 *   1:= Using self size, 
 *   2:= Using average size, 
*/
#define INTERACTION_TYPE 1      // Turns out the type 1 is the best

/** Flag for activation of boundary evaluation
 *   0:= Only take the equation inside domain, 
 *   1:= Take boundary condition into account
*/
#define BOUNDARY_ACTIVE 1

/** The global matrix basis construction type
 *   1:= Standard A, 
 *   2:= Standard B
*/
#define LSMPS_TYPE 1

/** The global matrix construction type
 *   1:= Sparse Matrix, 
 *   2:= Dense Matrix
*/
#define MATRIX_BASIS 1

/** Solver Type
 *           Sparse Matrix   |   Dense Matrix
 *       --------------------|-------------------
 *   1:= Sparse LU Decomp.   |  *LU Decomposition;
 *   2:= LDLT Method         |  *Full Piv Householder QR;
 *   3:= Sparse QR Method    |  *-;
 *   4:= BiCGStab Method     |  *-;
*/
#define MATRIX_SOLVER 1

// =====================================================
// +---------------- Class Constructor ----------------+
// =====================================================
// #pragma region CLASS_CONSTRUCTOR
PoissonLSMPS::PoissonLSMPS(/* args */)
{
}

PoissonLSMPS::~PoissonLSMPS()
{
}
// #pragma endregion


// =====================================================
// +----------------- Internal Method -----------------+
// =====================================================
// #pragma region INTERNAL_METHOD
/**
 *  @brief A quadratic weight function used in calculating least square in LSMPS.
 *  
 *  @param  _rij Distance between particle i and j.
 *  @param  _Re Support radius.
 * 
 *  @return The weight value.
 */
double PoissonLSMPS::weight_function(const double &_rij, const double &_Re)
{
    double _wij;
    if (_rij <= _Re){
        _wij = std::pow(1 - (_rij/_Re), 2.0);
    }
    else{
        _wij = 0.0e0;
    }

    return _wij;
}

/**
 *  @brief A taylor expansion position difference vector set.
 *  The vector size is determined by the degree of taylor expansion
 *  used, which is BASE_SIZE_2D.
 *  
 *  @param  dx Basis x position difference between particle i and j.
 *  @param  dy Basis y position difference between particle i and j.
 *  @param  rs  Scaling parameter.
 * 
 *  @return The taylor expansion position difference vector of P.
 */
vec<double> PoissonLSMPS::get_p(const double &xij, const double &yij, const double &rs)
{
    double _xij = xij / rs;
    double _yij = yij / rs;

    #if (LSMPS_TYPE == 1)
        vec<double> _p(BASE_SIZE_2D_A);
        _p[0] = _xij;
        _p[1] = _yij;
        _p[2] = _xij * _xij;
        _p[3] = _xij * _yij;
        _p[4] = _yij * _yij;
    #elif (LSMPS_TYPE == 2)
        vec<double> _p(BASE_SIZE_2D_B);
        _p[0] = 1;
        _p[1] = _xij;
        _p[2] = _yij;
        _p[3] = _xij * _xij;
        _p[4] = _xij * _yij;
        _p[5] = _yij * _yij;
    #endif

    return _p;
}
// #pragma endregion

// =====================================================
// +------------------ Public Method ------------------+
// =====================================================
// #pragma region PULBIC_METHOD
/**
 *  @brief Generate the global matrix to solve poisson problem.
 *  @param  coords  Particle data coordinate
 *  @param  nghList List of neighbor particle ID
 *  @param  size    The particle size container
*/
void PoissonLSMPS::create_global_matrix(
    const vec<vec<double>> &pos,
    const vec<vec<int>> &nghList,
    const vec<double> &s,
    const vec<bool> &bndFlag
){
    // There are two procedure in generating the global matrix
    /** Procedure
     *   1. Calculate the potential head constant from LSMPS standard
     *   2. Rearrange the calculated constant on the global matrix
    */

    // Update the class member
    this->N = s.size();     // Total particle number in the domain

    // Internal variable
    // vec<Triplet<double>> tripletList;
    // const double &_rs = Pars::sigma;        // Scaling size (the size of evaluated particle)

    // PROCEDURE 0:
    // ***********
    // Resize the global matrix 
    std::cout << "Start matrix initialization!\n";
    #if (MATRIX_BASIS == 1)     // Sparse Matrix
        vec<Triplet<double>> tripletList;
        this->A_sparse.resize(this->N, this->N);
    #elif (MATRIX_BASIS == 2)   // Dense Matrix
        this->A_dense = MatrixXd::Zero(this->N, this->N);
    #endif
    std::cout << "Done matrix initialization!\n";


    // PROCEDURE 1:
    // ***********
    // Similar procedure as calculating LSMPS
    #if (LSMPS_TYPE == 1)
        // Evaluate the standard LSMPS [type A] at each particle
        // #pragma omp parallel for     // Can't be parallelize
        for (int i = 0; i < this->N; i++)
        {
            // Note: Each particle discretization yields an equation that represent the row

            // Aliasing the target particle ID
            const int &idxi = i;

            // If the current particle is the boundary
            #if (BOUNDARY_ACTIVE == 1)
                if (bndFlag[i] == true){
                    // A dirichlet type of boundary condition
                    #if (MATRIX_BASIS == 1)     // Sparse Matrix
                        tripletList.push_back(Triplet<double>(idxi, idxi, 1.0));
                    #elif (MATRIX_BASIS == 2)   // Dense Matrix
                        this->A_dense(idxi, idxi) = 1.0;
                    #endif
                    continue;
                }
            #endif

            // Particle coordinate alias
            const double &_xi = pos[0].at(idxi);
            const double &_yi = pos[1].at(idxi);

            // Intermediate variable
            const vec<int> &ngh = nghList[idxi];    // Neighbor list of the current idxi particle
            const double &_rs = s[idxi];            // Scaling size (the size of evaluated particle)
            const int nghNum = ngh.size();          // The number of neighbor particle
            std::vector<bool> notNeighbor(nghNum, false);   // Neighbor flag container
                
            // Initialize the LSMPS matrices container
            MatrixXd Hbar = MatrixXd::Zero(BASE_SIZE_2D_A, BASE_SIZE_2D_A);   // Scaling matrix
            MatrixXd Mbar = MatrixXd::Zero(BASE_SIZE_2D_A, BASE_SIZE_2D_A);   // Moment matrix
            VectorXd bi = VectorXd::Zero(BASE_SIZE_2D_A);                     // A vector collection of moment vector
            vec<VectorXd> bbar = vec<VectorXd>(nghNum, VectorXd::Zero(BASE_SIZE_2D_A));   // Moment vector (separated for each neighbor)

            // Initialize H_rs Matrix
            Hbar(0, 0) = std::pow(_rs, -1);       // Dx(f)  | First order in x (1,0)
            Hbar(1, 1) = std::pow(_rs, -1);       // Dy(f)  | First order in y (0,1)
            Hbar(2, 2) = std::pow(_rs, -2) * 2;   // Dx2(f) | Second order in x (2,0)
            Hbar(3, 3) = std::pow(_rs, -2);       // Dxy(f) | First in x,y (1,1)
            Hbar(4, 4) = std::pow(_rs, -2) * 2;   // Dy2(f) | Second order in y (0,2)

            // Calculate the moment matrix [Mbar] and vector [bbar]
            for (int j = 0; j < nghNum; j++)
            {
                // Aliasing the neighbor particle ID
                const int &idxj = ngh[j];

                // No effect from the current neighbor (or neighbor outside the support domain)
                if (idxi == idxj){
                    notNeighbor[j] = true;
                    continue;
                }
                
                // Neighbor particle coordinate alias
                const double &_xj = pos[0].at(idxj);
                const double &_yj = pos[1].at(idxj);
                
                // Coordinate distance
                double _dx = _xj - _xi;
                double _dy = _yj - _yi;

                // Resultant distance(_rij) and effective radius(_Re)
                double _rij = std::sqrt(std::pow(_dx, 2) + std::pow(_dy, 2));
                double _Re;
                switch (INTERACTION_TYPE){
                case 1: // Effective radius type 1 (self particle size)
                    _Re = RADIUS_FACTOR_2D * s[idxi];
                    break;
                case 2: // Effective radius type 2 (average particle size)
                    _Re = RADIUS_FACTOR_2D * (s[idxi] + s[idxj]) * 0.5;
                    break;
                default:
                    break;
                }

                // No effect from the current neighbor (or neighbor outside the support domain)
                if (_rij > _Re){
                    notNeighbor[j] = true;
                    continue;
                }

                // Calculate the weight value
                double _wij = this->weight_function(_rij, _Re);     // Calculate the weight
                // _wij *= std::pow(s[idxj]/s[idxi], 2.0);             // Add the relative size factor between particle (area ratio of ngh to curr) [Cause more error]
                
                // Calculate the distance vector p
                const vec<double> _P = this->get_p(_dx, _dy, _rs);

                // Calculation of moment matrix Mbar and bbar
                for (int k1 = 0; k1 < BASE_SIZE_2D_A; k1++){
                    for (int k2 = 0; k2 < BASE_SIZE_2D_A; k2++){
                        // Generate tensor product between p
                        Mbar(k1, k2) += (_wij * _P[k1] * _P[k2]);
                    }
                    
                    // Calculate the vector matrix
                    bi(k1) += _wij * _P[k1];

                    // Generate moment matrix for each neighbor
                    bbar[j](k1) = _wij * _P[k1];
                }
            }
            
            // Inverse the moment matrix
            // --- Method 0 ---     Standard Matrix Inverse
            // MatrixXd MbarInv = Mbar.inverse();
            // --- Method 1 ---     LU Decomposition
            MatrixXd MbarInv = Mbar.lu().inverse();     // LU is the best (fast and not much error for small matrix [6x6])
            // --- Method 2 ---     House Holder QR 
            // VectorXd MbarInv = Mbar.fullPivHouseholderQr().inverse();
            
            // Calculate the matrix product
            MatrixXd Hrs_MbarInv = Hbar * MbarInv;
            
            // Calculate the global matrix component affected by the current particle and each neighbor
            for (int j = 0; j < nghNum; j++)
            {
                // Exception for ID that is not classified as neighbor
                if (notNeighbor[j] == true) continue;

                // Aliasing the particle ID
                const int &idxj = ngh[j];   // Neighbor particle ID
                
                // Calculation of the matrix element
                VectorXd Coef = Hrs_MbarInv * bbar[j];  // Coefficient addition of dx^2 and dy^2 part
                double elm = Coef(2) + Coef(4);         // The addition of dx^2 and dy^2 part

                // Assign the data into container
                #if (MATRIX_BASIS == 1)     // Sparse Matrix
                    tripletList.push_back(Triplet<double>(idxi, idxj, elm));
                #elif (MATRIX_BASIS == 2)   // Dense Matrix
                    this->A_dense(idxi, idxj) = elm;
                #endif
            }

            // Additional for the current element only
            {
                VectorXd Coef = Hrs_MbarInv * bi;   // Coefficient addition of dx^2 and dy^2 part
                double elm = Coef(2) + Coef(4);     // The addition of dx^2 and dy^2 part

                // Assign the data into container
                #if (MATRIX_BASIS == 1)     // Sparse Matrix
                    tripletList.push_back(Triplet<double>(idxi, idxi, -elm));
                #elif (MATRIX_BASIS == 2)   // Dense Matrix
                    this->A_dense(idxi, idxi) = -elm;
                #endif
            }
        }
    #elif (LSMPS_TYPE == 2)
        // Evaluate the standard LSMPS [type B] at each particle
        // #pragma omp parallel for
        for (int i = 0; i < this->N; i++)
        {
            // Note: Each particle discretization yields an equation that represent the row

            // Aliasing the target particle ID
            const int &idxi = i;

            // If the current particle is the boundary
            #if (BOUNDARY_ACTIVE == 1)
                if (bndFlag[i] == true){
                    // Assign the data into container
                    #if (MATRIX_BASIS == 1)     // Sparse Matrix
                        tripletList.push_back(Triplet<double>(idxi, idxi, 1.0));
                    #elif (MATRIX_BASIS == 2)   // Dense Matrix
                        this->A_dense(idxi, idxi) = 1.0;
                    #endif
                    continue;
                }
            #endif

            // Particle coordinate alias
            const double &_xi = pos[0].at(idxi);
            const double &_yi = pos[1].at(idxi);

            // Intermediate variable
            const vec<int> &ngh = nghList[idxi];    // Neighbor list of the current idxi particle
            const double &_rs = s[idxi];            // Scaling size (the size of evaluated particle)
            const int nghNum = ngh.size();          // The number of neighbor particle
            std::vector<bool> notNeighbor(nghNum, false);   // Neighbor flag container
                
            // Initialize the LSMPS matrices container
            MatrixXd Hbar = MatrixXd::Zero(BASE_SIZE_2D_B, BASE_SIZE_2D_B);   // Scaling matrix
            MatrixXd Mbar = MatrixXd::Zero(BASE_SIZE_2D_B, BASE_SIZE_2D_B);   // Moment matrix
            vec<VectorXd> bbar = vec<VectorXd>(nghNum, VectorXd::Zero(BASE_SIZE_2D_B));   // Moment vector (separated for each neighbor)

            // Initialize H_rs Matrix
            Hbar(0, 0) = 1;                       // D(f)   | Zero order (0,0)
            Hbar(1, 1) = std::pow(_rs, -1);       // Dx(f)  | First order in x (1,0)
            Hbar(2, 2) = std::pow(_rs, -1);       // Dy(f)  | First order in y (0,1)
            Hbar(3, 3) = std::pow(_rs, -2) * 2;   // Dx2(f) | Second order in x (2,0)
            Hbar(4, 4) = std::pow(_rs, -2);       // Dxy(f) | First in x,y (1,1)
            Hbar(5, 5) = std::pow(_rs, -2) * 2;   // Dy2(f) | Second order in y (0,2)

            // Calculate the moment matrix [Mbar] and vector [bbar]
            for (int j = 0; j < nghNum; j++)
            {
                // Aliasing the neighbor particle ID
                const int &idxj = ngh[j];
                
                // Neighbor particle coordinate alias
                const double &_xj = pos[0].at(idxj);
                const double &_yj = pos[1].at(idxj);
                
                // Coordinate distance
                double _dx = _xj - _xi;
                double _dy = _yj - _yi;

                // Resultant distance(_rij) and effective radius(_Re)
                double _rij = std::sqrt(std::pow(_dx, 2) + std::pow(_dy, 2));
                double _Re;
                switch (INTERACTION_TYPE){
                case 1: // Effective radius type 1 (self particle size)
                    _Re = RADIUS_FACTOR_2D * s[idxi];
                    break;
                case 2: // Effective radius type 2 (average particle size)
                    _Re = RADIUS_FACTOR_2D * (s[idxi] + s[idxj]) * 0.5;
                    break;
                default:
                    break;
                }

                // No effect from the current neighbor (or neighbor outside the support domain)
                if (_rij > _Re){
                    notNeighbor[j] = true;
                    continue;
                }

                // Calculate the weight value
                double _wij = this->weight_function(_rij, _Re);     // Calculate the weight
                // _wij *= std::pow(s[idxj]/s[idxi], 2.0);             // Add the relative size factor between particle (area ratio of ngh to curr)
                
                // Calculate the distance vector p
                const vec<double> _P = this->get_p(_dx, _dy, _rs);

                // Calculation of moment matrix Mbar and bbar
                for (int k1 = 0; k1 < BASE_SIZE_2D_B; k1++){
                    for (int k2 = 0; k2 < BASE_SIZE_2D_B; k2++){
                        // Generate tensor product between p
                        Mbar(k1, k2) += (_wij * _P[k1] * _P[k2]);
                    }
                    // Generate moment matrix for each neighbor
                    bbar[j](k1) = _wij * _P[k1];
                }
            }
            
            // Inverse the moment matrix
            // --- Method 0 ---     Standard Matrix Inverse
            // MatrixXd MbarInv = Mbar.inverse();
            // --- Method 1 ---     LU Decomposition
            MatrixXd MbarInv = Mbar.lu().inverse();
            // --- Method 2 ---     House Holder QR 
            // VectorXd MbarInv = Mbar.fullPivHouseholderQr().inverse();
            
            // Calculate the matrix product
            MatrixXd Hrs_MbarInv = Hbar * MbarInv;
            
            // Calculate the global matrix component affected by the current particle and each neighbor
            for (int j = 0; j < nghNum; j++)
            {
                // Exception for ID that is not classified as neighbor
                if (notNeighbor[j] == true) continue;

                // Aliasing the particle ID
                const int &idxj = ngh[j];   // Neighbor particle ID
                
                // Calculation of the matrix element
                VectorXd Coef = Hrs_MbarInv * bbar[j];  // Coefficient addition of dx^2 and dy^2 part
                double elm = Coef(3) + Coef(5);         // The addition of dx^2 and dy^2 part

                // Assign the data into container
                #if (MATRIX_BASIS == 1)     // Sparse Matrix
                    tripletList.push_back(Triplet<double>(idxi, idxj, elm));
                #elif (MATRIX_BASIS == 2)   // Dense Matrix
                    this->A_dense(idxi, idxj) = elm;
                #endif
            }
        }
    #endif
    
    // PROCEDURE 2:
    // ***********
    // Add the data into the global matrix 
    #if (MATRIX_BASIS == 1)     // Sparse Matrix
        this->A_sparse.setFromTriplets(tripletList.begin(),tripletList.end());
        this->A_sparse.makeCompressed();
    #endif

    return;
}

/**
 *  @brief Solve the poisson equation from the global matrix class member.
 *  @param  RHS The right hand side on poisson equation (or the source)
*/
void PoissonLSMPS::solve(const vec<double> &RHS){
    // Declare the solution vector
    VectorXd _phi;

    // =========================================
    // ----------- Assign the Source -----------
    // =========================================
    // Assign the source Term
    VectorXd b = VectorXd::Zero(this->N);
    for (int i = 0; i < this->N; i++)
    {
        b(i) = RHS[i];
    }

    // =========================================
    // ---------- Dense Matrix Solver ----------
    // =========================================
    #if (MATRIX_BASIS == 1)
        #if (MATRIX_SOLVER == 1)
            // [TYPE 1] Solve the dense global matrix using LU decomposition
            _phi = this->A_dense.lu().solve(b);
            /**
             * NOTE: Only this one works, but still broken for multiresolution
             * - Solve 10000 particle in 17 secs
             * - Solve 22500 particle in 170 secs
            */

        #elif (MATRIX_SOLVER == 2)
            // [TYPE 2] Solve the dense global matrix using full Piv Householder QR
            _phi = this->A_dense.fullPivHouseholderQr().solve(b);
            /**
             * NOTE: -
            */
        #endif

    // =========================================
    // --------- Sparse Matrix Solver ----------
    // =========================================
    // #elif (MATRIX_BASIS == 2)
        #if (MATRIX_SOLVER == 1)
            // [TYPE 1] Using sparse LU decomposition
            SparseLU <SparseMatrix<double>> solver;
            solver.compute(this->A_sparse);
            _phi = solver.solve(b);
            /**
             * NOTE: Cannot solve the problem <?>
             *  - Solve 40000 particle takes 3.9 secs
            */

        #elif (MATRIX_SOLVER == 2)
            // [TYPE 2] Using sparse Simple LDLT
            // Eigen::initParallel();
            SimplicialLDLT <SparseMatrix<double>> solver;
            solver.compute(this->A_sparse);
            if (solver.info() != Success) {      // Only for checking
                std::cout << "Decomposition failed" << std::endl;
                throw std::exception();
            }
            _phi = solver.solve(b);
            if (solver.info() != Success) {      // Only for checking
                std::cout << "Solving failed" << std::endl;
                throw std::exception();
            }
            /**
             * NOTE: Works with result quite well, but still have multiresolution problem
             * - 40000 par takes 2.2 secs
             * - 39000 par mulres takes 206 sec
            */
        
        #elif (MATRIX_SOLVER == 3)
            // [TYPE 3] Using sparse QR
            // Eigen::initParallel();
            Eigen::SparseQR <SparseMatrix<double>, COLAMDOrdering<int>> solver;
            solver.compute(this->A_sparse);
            _phi = solver.solve(b);
            /**
             * NOTE: Cant Solve the problem <?>
            */
        
        #elif (MATRIX_SOLVER == 4)
            // [TYPE 4] Using paralel BiCGStab
            // Operator initialization
            Eigen::initParallel();
            omp_set_num_threads(10);
            Eigen::setNbThreads(10);
            BiCGSTAB <SparseMatrix<double, RowMajor>> solver;
            solver.setTolerance(1e-2);      // Tolerance at this moment still 10^-1 (the fastest)
            
            // Solve the matrix
            solver.compute(this->A_sparse);
            _phi = solver.solve(b);
            /**
             * NOTE: Still don't know about this method. Anyway, this method is not work
             * -
            */
        #endif
    #endif

    
    // =========================================
    // ---------- Update the Solution ----------
    // =========================================
    this->phi.resize(this->N);
    for (int i = 0; i < this->N; i++)
    {
        this->phi[i] = _phi(i);
    }

    return;
}

/**
 *  @brief Get the poisson solution.
 *  @return The potential of poisson equation.
*/
vec<double> PoissonLSMPS::get_phi(){
    return this->phi;
}
// #pragma endregion