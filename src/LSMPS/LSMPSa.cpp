#include "LSMPSa.hpp"
#include "../Eigen/Dense"

#define BASE_SIZE_2D 5          // Size of base matrix in 2D
#define BASE_SIZE_3D 9          // Size of base matrix in 3D

#define RADIUS_FACTOR_2D 3.1    // Support radius factor toward particle size in 2D
#define RADIUS_FACTOR_3D 2.3    // Support radius factor toward particle size in 3D

/** The type of particle interation
 *   1:= Using self size, 
 *   2:= Using average size, 
*/
#define INTERACTION_TYPE 1

/** Linear algebra solver type
 *   1:= Partial LU decomposition,      : Fastest, medium accuracy
 *   2:= Full LU decomposition,         : Average, high accuracy
 *   3:= Singular Value Decomposition   : Slow   , very accurate
 *   4:= Jacoby Rotation                : Slow   , very accurate (Actually similar to SVD)
*/
#define LIN_ALG_TYPE 1
// NOTE: For some reasone partial LU have higher accuracy than SVD for small matrix size

// For linear algebra operation
using namespace Eigen;

// Shorthand for vector declaration
template <typename U> using vec = std::vector<U>;

// =====================================================
// +----------------- Public Function -----------------+
// =====================================================
// #pragma region PUBLIC

/**
 *  @brief An LSMPS type A function differential calculation of a particle data set 
 *  toward itself.
 *  
 *  @param  x Particle x coordinate data.
 *  @param  y Particle y coordinate data.
 *  @param  s Particle size data.
 *  @param  f Particle function value data.
 *  @param  neighborList The neighbor relation of particle data toward itself.
 */
void LSMPSa::set_LSMPS(const vec<double> &x, const vec<double> &y, 
                       const vec<double> &s, const vec<double> &f, 
                       const vec<vec<int>> &neighborlist)
{
    // clear local variables
    this->_ddx.clear();
    this->_ddy.clear();
    this->_d2d2x.clear();
    this->_d2dxdy.clear();
    this->_d2d2y.clear();

    // resize local variabels
    int nparticle = s.size();
    this->_ddx.resize(nparticle);
    this->_ddy.resize(nparticle);
    this->_d2d2x.resize(nparticle);
    this->_d2dxdy.resize(nparticle);
    this->_d2d2y.resize(nparticle);

    // perform calculation
    this->calculate_LSMPS(x, y, s, f, neighborlist);
}

/**
 *  @brief An LSMPS type A function differential calculation of a particle data set 
 *  toward a collocation data set.
 *  
 *  @param  xTarget Target particle x coordinate data.
 *  @param  yTarget Target particle y coordinate data.
 *  @param  sTarget Target particle size data.
 *  @param  fTarget Target particle function value data.
 *  @param  xCollocation Collocation particle x coordinate data.
 *  @param  yCollocation Collocation particle y coordinate data.
 *  @param  sCollocation Collocation particle size data.
 *  @param  fCollocation Collocation particle function value data.
 *  @param  neighborList The neighbor relation of target particle data toward collocation particle.
*/
void LSMPSa::set_LSMPS(const vec<double> &xT, const vec<double> &yT, const vec<double> &sT, const vec<double> &fT,
                       const vec<double> &xC, const vec<double> &yC, const vec<double> &sC, const vec<double> &fC,
                       const vec<vec<int>> &nghList)
{
    // clear local variables
    this->_ddx.clear();
    this->_ddy.clear();
    this->_d2d2x.clear();
    this->_d2dxdy.clear();
    this->_d2d2y.clear();

    // resize local variabels
    int nparticle = sT.size();
    this->_ddx.resize(nparticle);
    this->_ddy.resize(nparticle);
    this->_d2d2x.resize(nparticle);
    this->_d2dxdy.resize(nparticle);
    this->_d2d2y.resize(nparticle);

    // perform calculation
    this->calculate_LSMPS(xT, yT, sT, fT, xC, yC, sC, fC, nghList);
}

/**
 *  @brief A 3D LSMPS type A function differential calculation of a particle data set.
 *  This function follow the pattern of collocation LSMPS set. 
 *  
 *  @param  xTarget Target particle x coordinate data.
 *  @param  yTarget Target particle y coordinate data.
 *  @param  zTarget Target particle z coordinate data.
 *  @param  sTarget Target particle size data.
 *  @param  fTarget Target particle function value data.
 *  @param  xCollocation Collocation particle x coordinate data.
 *  @param  yCollocation Collocation particle y coordinate data.
 *  @param  zCollocation Collocation particle z coordinate data.
 *  @param  sCollocation Collocation particle size data.
 *  @param  fCollocation Collocation particle function value data.
 *  @param  neighborList The neighbor relation of target particle data toward collocation particle.
*/
void LSMPSa::set_LSMPS_3D(const vec<double> &xT, const vec<double> &yT, const vec<double> &zT,
                          const vec<double> &sT, const vec<double> &fT,
                          const vec<double> &xC, const vec<double> &yC, const vec<double> &zC,
                          const vec<double> &sC, const vec<double> &fC,
                          const vec<vec<int>> &nghList)
{
    // clear local variables
    this->_ddx.clear();
    this->_ddy.clear();
    this->_ddz.clear();
    this->_d2d2x.clear();
    this->_d2dxdy.clear();
    this->_d2dxdz.clear();
    this->_d2d2y.clear();
    this->_d2dydz.clear();
    this->_d2d2z.clear();

    // resize local variabels
    int nparticle = sT.size();
    this->_ddx.resize(nparticle);
    this->_ddy.resize(nparticle);
    this->_ddz.resize(nparticle);
    this->_d2d2x.resize(nparticle);
    this->_d2dxdy.resize(nparticle);
    this->_d2dxdz.resize(nparticle);
    this->_d2d2y.resize(nparticle);
    this->_d2dydz.resize(nparticle);
    this->_d2d2z.resize(nparticle);

    // perform calculation
    this->calculate_LSMPS_3D(xT, yT, zT, sT, fT, xC, yC, zC, sC, fC, nghList);
    
    return;
}

// #pragma endregion

// =====================================================
// +--------------- LSMPSa Calculation ----------------+
// =====================================================
// #pragma region LSMPS_A_CALCULATION

/**
 *  @brief An LSMPS type A function differential calculation of a particle data set 
 *  toward itself.
 *  
 *  @param  x Particle x coordinate data.
 *  @param  y Particle y coordinate data.
 *  @param  s Particle size data.
 *  @param  f Particle function value data.
 *  @param  neighborList The neighbor relation of particle data toward itself.
 */
void LSMPSa::calculate_LSMPS(const vec<double> &x, const vec<double> &y, 
                             const vec<double> &s, const vec<double> &f,
                             const vec<vec<int>> &nghList)
{
    // Take the number of particle
    int nparticle = s.size();
    
    // Evaluate the LSMPS A of all particle
    #pragma omp parallel for
    for (int i = 0; i < nparticle; i++)
    {
        // Initialize the LSMPS matrices
        MatrixXd Hi = MatrixXd::Zero(BASE_SIZE_2D, BASE_SIZE_2D);   // Scaling matrix
        MatrixXd Mi = MatrixXd::Zero(BASE_SIZE_2D, BASE_SIZE_2D);   // Moment matrix
        VectorXd bi = VectorXd::Zero(BASE_SIZE_2D);              // Moment vector

        // Scaling size
        double _rs = s[i];
        
        // Initialize H_rs Matrix
        Hi(0, 0) = std::pow(_rs, -1);
        Hi(1, 1) = std::pow(_rs, -1);
        Hi(2, 2) = std::pow(_rs, -2) * 2;
        Hi(3, 3) = std::pow(_rs, -2);
        Hi(4, 4) = std::pow(_rs, -2) * 2;

        // Calculate the moment matrix [Mi] and vector [bi]
        const vec<int> &_nghID = nghList[i];
        for (size_t j = 0; j < _nghID.size(); j++)
        {
            // Aliasing the particle ID
            const int &idxi = i;           // Target Particle ID
            const int &idxj = _nghID[j];   // Neighbor Particle ID
            
            // Target particle coordinate
            double _xi = x[idxi];
            double _yi = y[idxi];
            
            // Neighbor particle coordinate
            double _xj = x[idxj];
            double _yj = y[idxj];

            // Coordinate difference between Particle and Neighbor Particle (Note the sign: neighbor - target)
            double _dx = _xj - _xi;
            double _dy = _yj - _yi;

            // Resultant distance(_rij) and effective radius(_Re)
            double _rij = std::sqrt(std::pow(_dx, 2) + std::pow(_dy, 2)); //distance between particles. 
            double _Re = 0.0;
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
            if (_rij > _Re) continue;
            
            // Calculate the weight value
            double _wij = this->weight_function(_rij, _Re);     // Calculate the weight
            _wij *= std::pow(s[idxj]/s[idxi],2);                // Add the relative size factor between particle

            // Calculate the distance vector p
            vec<double> _p1 = this->get_p(_dx, _dy, _rs);
            
            // LSMPS A calculated variable value different
            double _df = f[idxj] - f[idxi];                     // Neighbor particle - evaluated particle
            
            // Calculatin moment matrix Mi and bi
            for (size_t k1 = 0; k1 < BASE_SIZE_2D; k1++){
                for (size_t k2 = 0; k2 < BASE_SIZE_2D; k2++){
                    // Generate tensor product between p
                    Mi(k1, k2) = Mi(k1, k2) + (_wij * _p1[k1] * _p1[k2]);
                }
                // Generate moment matrix
                bi(k1) = bi(k1) + (_wij * _p1[k1] * _df);
            }
        }

        // Calculate the product between inverse of Mi and bi
        VectorXd MiInv_Bi;
        MiInv_Bi = Mi.lu().solve(bi);
        
        // switch (LIN_ALG_TYPE){
        // case 1: // Partial LU decompostion
        //     MiInv_Bi = Mi.lu().solve(bi);
        //     break;
        // case 2: // Full LU decompostion
        //     MiInv_Bi = Mi.fullPivLu().solve(bi);
        //     break;
        // case 3: // Singular value decompostion
        //     MiInv_Bi = Mi.bdcSvd(ComputeThinU | ComputeThinV).solve(bi);
        //     break;
        // default:
        //     break;
        // }
        
        // Calculate the function derivative [Df] using LSMPSa formula
        VectorXd Df = Hi * MiInv_Bi;

        // Assign to private variables
        this->_ddx[i]    = Df[0];
        this->_ddy[i]    = Df[1];
        this->_d2d2x[i]  = Df[2];
        this->_d2dxdy[i] = Df[3];
        this->_d2d2y[i]  = Df[4];
    }
}

/**
 *  @brief An LSMPS type A function differential calculation of a particle data set 
 *  toward a collocation data set.
 *  
 *  @param  xTarget Target particle x coordinate data.
 *  @param  yTarget Target particle y coordinate data.
 *  @param  sTarget Target particle size data.
 *  @param  fTarget Target particle function value data.
 *  @param  xCollocation Collocation particle x coordinate data.
 *  @param  yCollocation Collocation particle y coordinate data.
 *  @param  sCollocation Collocation particle size data.
 *  @param  fCollocation Collocation particle function value data.
 *  @param  neighborList The neighbor relation of target particle data toward collocation particle.
 */
void LSMPSa::calculate_LSMPS(const vec<double> &xTar, const vec<double> &yTar, const vec<double> &sTar, const vec<double> &fTar,
                             const vec<double> &xCol, const vec<double> &yCol, const vec<double> &sCol, const vec<double> &fCol,
                             const vec<vec<int>> &nghList)
{
    /** LSMPS A collocation Note :
     *   > There are two relative position in this operation
     *     * The particle target interpolation (refer as TARGET)
     *     * The particle interpolation data (refer as COLLOCATION)
     *   > Please note that TARGET and COLLOCATION just a reference name
     *     * Both TARGET and COLLOCATION is particle variable
     *     * Both TARGET and COLLOCATION come from the same particle variable
     *   
     *   TARGET VARIABLE - interpolation target
     *     vec<dbl> xTarget, vec<dbl> yTarget  : Particle position list
     *     vec<dbl> sTarget                    : Particle size list
     *     vec<dbl> fTarget                    : Interpolated variable data
     *    
     *   COLLOCATION VARIABLE - interpolation data
     *     vec<dbl> xCollocation, vec<dbl> yCollocation   : Particle position list
     *     vec<dbl> sCollocation                          : Particle size list
     *     vec<dbl> fCollocation                          : Interpolated variable data
     *     vec<vec<int>> neighborlist                     : Grid neighbor ID list consisted of particle ID
    */

    // Take the number of particle
    int nparticle = sTar.size();

    // Evaluate the LSMPS A of all particle
    #pragma omp parallel for
    for (int i = 0; i < nparticle; i++)
    {
        // Initialize the LSMPS matrices
        MatrixXd Hi = MatrixXd::Zero(BASE_SIZE_2D, BASE_SIZE_2D);     // Scaling matrix
        MatrixXd Mi = MatrixXd::Zero(BASE_SIZE_2D, BASE_SIZE_2D);     // Moment matrix
        VectorXd bi = VectorXd::Zero(BASE_SIZE_2D);                // Moment vector
        
        // Scaling size
        double _rs = sTar[i];

        // Initialize H_rs Matrix
        Hi(0, 0) = std::pow(_rs, -1);
        Hi(1, 1) = std::pow(_rs, -1);
        Hi(2, 2) = std::pow(_rs, -2) * 2;
        Hi(3, 3) = std::pow(_rs, -2);
        Hi(4, 4) = std::pow(_rs, -2) * 2;

        // Calculate the moment matrix [Mi] and vector [bi]
        const vec<int> &_nghID = nghList[i];
        for (size_t j = 0; j < _nghID.size(); j++)
        {
            // Note that neighbor particle and evaluated particle inside the same variable
            const int &idxi = i;          // Target Particle ID
            const int &idxj = _nghID[j];  // Neighbor Particle ID

            // Evaluated particle coordinate (TARGET)
            double _xi = xTar[idxi];
            double _yi = yTar[idxi];

            // Neighbor particle coordinate (COLLOCATION)
            double _xj = xCol[idxj];
            double _yj = yCol[idxj];

            // Coordinate difference between TARGET and COLLOCATION Particle (Note the sign: neighbor - target)
            double _dx = _xj - _xi;
            double _dy = _yj - _yi;
            
            // Resultant distance(_rij) and effective radius(_Re)
            double _rij = std::sqrt(std::pow(_dx, 2) + std::pow(_dy, 2));
            double _Re = 0.0;
            switch (INTERACTION_TYPE){
            case 1: // Effective radius type 1 (self particle size)
                _Re = RADIUS_FACTOR_2D * sTar[idxi];
                break;
            case 2: // Effective radius type 2 (average particle size)
                _Re = RADIUS_FACTOR_2D * (sTar[idxi] + sCol[idxj]) * 0.5;
                break;
            default:
                break;
            }
            
            // No effect from the current neighbor (or neighbor outside the support domain)
            if (_rij > _Re) continue;
            
            // Calculate the weight value
            double _wij = this->weight_function(_rij, _Re);     // Calculate the weight
            _wij *= std::pow(sCol[idxj]/sTar[idxi], 2.0);       // Add the relative size factor between particle

            // Calculate the distance vector p
            const vec<double> _P = this->get_p(_dx, _dy, _rs);

            // LSMPS A interpolated variable value different
            double _fij = fCol[idxj] - fTar[idxi];              // Neighbor particle - evaluated particle

            // Calculation of moment matrix Mi and bi
            for (size_t k1 = 0; k1 < BASE_SIZE_2D; k1++){
                for (size_t k2 = 0; k2 < BASE_SIZE_2D; k2++){
                    // Generate tensor product between p
                    Mi(k1, k2) = Mi(k1, k2) + (_wij * _P[k1] * _P[k2]);
                }
                // Generate moment matrix
                bi(k1) = bi(k1) + (_wij * _P[k1] * _fij);
            }
        }

        // Calculate the product between inverse of Mi and bi
        VectorXd MiInv_Bi = Mi.lu().solve(bi);
        
        // Calculate the function derivative [Df] using LSMPSa formula        
        VectorXd Df = Hi * MiInv_Bi;

        // Assign to private variables
        this->_ddx[i]    = Df[0];
        this->_ddy[i]    = Df[1];
        this->_d2d2x[i]  = Df[2];
        this->_d2dxdy[i] = Df[3];
        this->_d2d2y[i]  = Df[4];
    }
}


/**
 *  @brief A 3D LSMPS type A function differential calculation of a particle data set.
 *  This function follow the pattern of collocation LSMPS set. 
 *  
 *  @param  xTarget Target particle x coordinate data.
 *  @param  yTarget Target particle y coordinate data.
 *  @param  zTarget Target particle z coordinate data.
 *  @param  sTarget Target particle size data.
 *  @param  fTarget Target particle function value data.
 *  @param  xCollocation Collocation particle x coordinate data.
 *  @param  yCollocation Collocation particle y coordinate data.
 *  @param  zCollocation Collocation particle z coordinate data.
 *  @param  sCollocation Collocation particle size data.
 *  @param  fCollocation Collocation particle function value data.
 *  @param  neighborList The neighbor relation of target particle data toward collocation particle.
*/
void LSMPSa::calculate_LSMPS_3D(const vec<double> &xTar, const vec<double> &yTar, const vec<double> &zTar,
                                const vec<double> &sTar, const vec<double> &fTar,
                                const vec<double> &xCol, const vec<double> &yCol, const vec<double> &zCol,
                                const vec<double> &sCol, const vec<double> &fCol,
                                const vec<vec<int>> &nghList)
{
    // Take the number of particle
    int nparticle = sTar.size();

    // Evaluate the LSMPS A of all particle
    #pragma omp parallel for
    for (int i = 0; i < nparticle; i++)
    {
        // Initialize the LSMPS matrices
        MatrixXd Hi = MatrixXd::Zero(BASE_SIZE_3D, BASE_SIZE_3D);     // Scaling matrix
        MatrixXd Mi = MatrixXd::Zero(BASE_SIZE_3D, BASE_SIZE_3D);     // Moment matrix
        VectorXd bi = VectorXd::Zero(BASE_SIZE_3D);                // Moment vector
        
        // Scaling size
        double _rs = sTar[i];

        // Initialize H_rs Matrix: (rs^-|a|)*(a!)
        Hi(0, 0) = std::pow(_rs, -1);       //  x(1,0,0)
        Hi(1, 1) = std::pow(_rs, -1);       //  y(0,1,0)
        Hi(2, 2) = std::pow(_rs, -1);       //  z(0,0,1)
        Hi(3, 3) = std::pow(_rs, -2) * 2;   // x2(2,0,0)
        Hi(4, 4) = std::pow(_rs, -2);       // xy(1,1,0)
        Hi(5, 5) = std::pow(_rs, -2);       // xz(1,0,1)
        Hi(6, 6) = std::pow(_rs, -2) * 2;   // y2(0,2,0)
        Hi(7, 7) = std::pow(_rs, -2);       // yz(0,1,1)
        Hi(8, 8) = std::pow(_rs, -2) * 2;   // z2(0,0,2)

        // Calculate the moment matrix [Mi] and vector [bi]
        const vec<int> &_nghID = nghList[i];
        for (size_t j = 0; j < _nghID.size(); j++)
        {
            // Note that neighbor particle and evaluated particle inside the same variable
            const int &idxi = i;          // Target Particle ID
            const int &idxj = _nghID[j];  // Neighbor Particle ID

            // Evaluated particle coordinate (TARGET)
            double _xi = xTar[idxi];
            double _yi = yTar[idxi];
            double _zi = zTar[idxi];

            // Neighbor particle coordinate (COLLOCATION)
            double _xj = xCol[idxj];
            double _yj = yCol[idxj];
            double _zj = zCol[idxj];

            // Coordinate difference between TARGET and COLLOCATION Particle (Note the sign: neighbor - target)
            double _dx = _xj - _xi;
            double _dy = _yj - _yi;
            double _dz = _zj - _zi;
            
            // Resultant distance(_rij) and effective radius(_Re)
            double _rij = std::sqrt(std::pow(_dx, 2) + std::pow(_dy, 2) + std::pow(_dz, 2));
            double _Re = 0.0;
            switch (INTERACTION_TYPE){
            case 1: // Effective radius type 1 (self particle size)
                _Re = RADIUS_FACTOR_3D * sTar[idxi];
                break;
            case 2: // Effective radius type 2 (average particle size)
                _Re = RADIUS_FACTOR_3D * (sTar[idxi] + sCol[idxj]) * 0.5;
                break;
            default:
                break;
            }
            
            // No effect from the current neighbor (or neighbor outside the support domain)
            if (_rij > _Re) continue;
            
            // Calculate the weight value
            double _wij = this->weight_function(_rij, _Re);     // Calculate the weight
            _wij *= std::pow(sCol[idxj]/sTar[idxi], 3.0);       // Add the relative size factor between particle (volume factor)

            // Calculate the distance vector p
            vec<double> _P;
            this->get_p_3d(_P, _dx, _dy, _dz, _rs);

            // LSMPS A interpolated variable value different
            double _fij = fCol[idxj] - fTar[idxi];              // Neighbor particle - evaluated particle

            // Calculation of moment matrix Mi and bi
            for (size_t k1 = 0; k1 < BASE_SIZE_3D; k1++){
                for (size_t k2 = 0; k2 < BASE_SIZE_3D; k2++){
                    // Generate tensor product between p
                    Mi(k1, k2) = Mi(k1, k2) + (_wij * _P[k1] * _P[k2]);
                }
                // Generate moment matrix
                bi(k1) = bi(k1) + (_wij * _P[k1] * _fij);
            }
        }

        // Calculate the product between inverse of Mi and bi
        VectorXd MiInv_Bi = Mi.lu().solve(bi);
        
        // Calculate the function derivative [Df] using LSMPSa formula        
        VectorXd Df = Hi * MiInv_Bi;

        // Assign to private variables
        this->_ddx[i]    = Df[0];       //  x(1,0,0)
        this->_ddy[i]    = Df[1];       //  y(0,1,0)
        this->_ddz[i]    = Df[2];       //  z(0,0,1)
        this->_d2d2x[i]  = Df[3];       // x2(2,0,0)
        this->_d2dxdy[i] = Df[4];       // xy(1,1,0)
        this->_d2dxdz[i] = Df[5];       // xz(1,0,1)
        this->_d2d2y[i]  = Df[6];       // y2(0,2,0)
        this->_d2dydz[i] = Df[7];       // yz(0,1,1)
        this->_d2d2z[i]  = Df[8];       // z2(0,0,2)

    }
    return;
}
// #pragma endregion

// =====================================================
// +----------------- Utilities Method ----------------+
// =====================================================
// #pragma region UTILITIES_METHOD

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
vec<double> LSMPSa::get_p(const double dx, const double dy, const double rs)
{
    vec<double> _p(BASE_SIZE_2D);

    double _dx = dx / rs;
    double _dy = dy / rs;

    _p[0] = _dx;
    _p[1] = _dy;
    _p[2] = _dx * _dx;
    _p[3] = _dx * _dy;
    _p[4] = _dy * _dy;

    return _p;
}

/**
 *  @brief A taylor expansion position difference vector set for 3D.
 *  The vector size is determined by the degree of taylor expansion
 *  used, which is BASE_SIZE_3D.
 *  
 *  @param  _P [OUTPUT] The vector P.
 *  @param  dx Basis x position difference between particle i and j.
 *  @param  dy Basis y position difference between particle i and j.
 *  @param  dz Basis z position difference between particle i and j.
 *  @param  rs  Scaling parameter.
 * 
 *  @return The taylor expansion position difference vector of P.
 */
void LSMPSa::get_p_3d(std::vector<double> &_P, const double dx, 
                      const double dy, const double dz, const double rs)
{
    _P.resize(BASE_SIZE_3D);

    double _dx = dx / rs;
    double _dy = dy / rs;
    double _dz = dz / rs;

    _P[0] = _dx;            //  x(1,0,0)
    _P[1] = _dy;            //  y(0,1,0)
    _P[2] = _dz;            //  z(0,0,1)
    _P[3] = _dx * _dx;      // x2(2,0,0)
    _P[4] = _dx * _dy;      // xy(1,1,0)
    _P[5] = _dx * _dz;      // xz(1,0,1)
    _P[6] = _dy * _dy;      // y2(0,2,0)
    _P[7] = _dy * _dz;      // yz(0,1,1)
    _P[8] = _dz * _dz;      // z2(0,0,2)

    return;
}

/**
 *  @brief A quadratic weight function used in calculating least square in LSMPS.
 *  
 *  @param  _rij Distance between particle i and j.
 *  @param  _Re Support radius.
 * 
 *  @return The weight value.
 */
double LSMPSa::weight_function(const double &rij, const double &_Re)
{
    double _wij;
    if (rij <= _Re){
        _wij = std::pow(1 - (rij/_Re), 2.0);
    }
    else{
        _wij = 0.0e0;
    }

    return _wij;
}


/**
 *  @brief Get the first x derivative (df/dx).
 *  @return The df/dx derivative.
 */
vec<double> LSMPSa::get_ddx(){
    return this->_ddx;
}

/**
 *  @brief Get the first y derivative (df/dy).
 *  @return The df/dy derivative.
 */
vec<double> LSMPSa::get_ddy(){
    return this->_ddy;
}

/**
 *  @brief Get the first z derivative (df/dz).
 *  @return The df/dz derivative.
 */
vec<double> LSMPSa::get_ddz(){
    return this->_ddz;
}

/**
 *  @brief Get the second x derivative (d2f/dx2).
 *  @return The d2f/dx2 derivative.
 */
vec<double> LSMPSa::get_d2d2x(){
    return this->_d2d2x;
}

/**
 *  @brief Get the first x and y derivative (d2f/dxdy).
 *  @return The d2f/dxdy derivative.
 */
vec<double> LSMPSa::get_d2dxdy(){
    return this->_d2dxdy;
}

/**
 *  @brief Get the first x and z derivative (d2f/dxdz).
 *  @return The d2f/dxdz derivative.
 */
vec<double> LSMPSa::get_d2dxdz(){
    return this->_d2dxdz;
}

/**
 *  @brief Get the second y derivative (d2f/dy2).
 *  @return The d2f/dy2 derivative.
 */
vec<double> LSMPSa::get_d2d2y(){
    return this->_d2d2y;
}

/**
 *  @brief Get the first y and z derivative (d2f/dydz).
 *  @return The d2f/dydz derivative.
 */
vec<double> LSMPSa::get_d2dydz(){
    return this->_d2dydz;
}

/**
 *  @brief Get the second y derivative (d2f/dz2).
 *  @return The d2f/dz2 derivative.
 */
vec<double> LSMPSa::get_d2d2z(){
    return this->_d2d2z;
}

// #pragma endregion