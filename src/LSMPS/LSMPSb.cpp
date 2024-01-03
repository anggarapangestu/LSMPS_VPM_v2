#include "LSMPSb.hpp"
#include "../Eigen/Dense"

#define BASE_SIZE_2D 6          // Size of base matrix
#define BASE_SIZE_3D 10         // Size of base matrix

#define RADIUS_FACTOR_2D 3.1    // Support radius factor toward particle size
#define RADIUS_FACTOR_3D 2.3    // Support radius factor toward particle size

/** The type of particle interation
 *   1:= Using self size, 
 *   2:= Using average size, 
*/
#define INTERACTION_TYPE 1

// For linear algebra operation
using namespace Eigen;

// Shorthand for vector declaration
template <typename U> using vec = std::vector<U>;

// =====================================================
// +----------------- Public Function -----------------+
// =====================================================
// #pragma region PUBLIC

/**
 *  @brief An LSMPS type B, a function differential calculation of a particle data set 
 *  toward a source data set. This method also works as an data interpolation to 
 *  predict the function differential at TARGET point from SOURCE data.
 *  
 *  @param  xTarget Target particle x coordinate data.
 *  @param  yTarget Target particle y coordinate data.
 *  @param  sTarget Target particle size data.
 *  @param  fTarget Target particle function value data.
 *  @param  xSource Source particle x coordinate data.
 *  @param  ySource Source particle y coordinate data.
 *  @param  sSource Source particle size data.
 *  @param  fSource Source particle function value data.
 *  @param  neighborListSource The neighbor relation of target particle data toward source particle.
 *  @param  neighborListSelf The neighbor relation of target particle data toward itself.
*/
void LSMPSb::set_LSMPS(const vec<double> &xTar, const vec<double> &yTar, const vec<double> &sTar, const vec<double> &fTar,
                       const vec<double> &xSrc, const vec<double> &ySrc, const vec<double> &sSrc, const vec<double> &fSrc,
                       const vec<vec<int>> &nghListSrc, const vec<vec<int>> &nghListSelf)
{
    // clear local variables
    this->_d0.clear();
    this->_ddx.clear();
    this->_ddy.clear();
    this->_d2d2x.clear();
    this->_d2dxdy.clear();
    this->_d2d2y.clear();

    // resize local variabels
    int nparticle = sTar.size();
    this->_d0.resize(nparticle);
    this->_ddx.resize(nparticle);
    this->_ddy.resize(nparticle);
    this->_d2d2x.resize(nparticle);
    this->_d2dxdy.resize(nparticle);
    this->_d2d2y.resize(nparticle);

    // perform calculation
    this->calculate_LSMPS(xTar, yTar, sTar, fTar, xSrc, ySrc, sSrc, fSrc, nghListSrc, nghListSelf);
}

/**
 *  @brief A 3D LSMPS type B function differential calculation of a particle data set 
 *  toward a source data set.
 *  
 *  @param  xTarget Target particle x coordinate data.
 *  @param  yTarget Target particle y coordinate data.
 *  @param  zTarget Target particle z coordinate data.
 *  @param  sTarget Target particle size data.
 *  @param  fTarget Target particle function value data.
 *  @param  xSource Source particle x coordinate data.
 *  @param  ySource Source particle y coordinate data.
 *  @param  zSource Source particle z coordinate data.
 *  @param  sSource Source particle size data.
 *  @param  fSource Source particle function value data.
 *  @param  neighborListSource The neighbor relation of target particle data toward source particle.
 *  @param  neighborListSelf The neighbor relation of target particle data toward itself.
*/
void LSMPSb::set_LSMPS_3D(const vec<double> &xTar, const vec<double> &yTar, const vec<double> &zTar,
                          const vec<double> &sTar, const vec<double> &fTar,
                          const vec<double> &xSrc, const vec<double> &ySrc, const vec<double> &zSrc,
                          const vec<double> &sSrc, const vec<double> &fSrc,
                          const vec<vec<int>> &nghListSrc, const vec<vec<int>> &nghListSelf)
{
    // clear local variables
    this->_d0.clear();
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
    int nparticle = sTar.size();
    this->_d0.resize(nparticle);
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
    this->calculate_LSMPS_3D(xTar, yTar, zTar, sTar, fTar, xSrc, ySrc, zSrc, sSrc, fSrc, nghListSrc, nghListSelf);
    
    return;
}

// #pragma endregion

// =====================================================
// +--------------- LSMPSb Calculation ----------------+
// =====================================================
// #pragma region LSMPS_B_CALCULATION

/**
 *  @brief An LSMPS type B, a function differential calculation of a particle data set 
 *  toward a source data set. This method also works as an data interpolation to 
 *  predict the function differential at TARGET point from SOURCE data.
 *  
 *  @param  xTarget Target particle x coordinate data.
 *  @param  yTarget Target particle y coordinate data.
 *  @param  sTarget Target particle size data.
 *  @param  fTarget Target particle function value data.
 *  @param  xSource Source particle x coordinate data.
 *  @param  ySource Source particle y coordinate data.
 *  @param  sSource Source particle size data.
 *  @param  fSource Source particle function value data.
 *  @param  neighborListSource The neighbor relation of target particle data toward source particle.
 *  @param  neighborListSelf The neighbor relation of target particle data toward itself.
*/
void LSMPSb::calculate_LSMPS(const vec<double> &xTar, const vec<double> &yTar, const vec<double> &sTar, const vec<double> &fTar,
                             const vec<double> &xSrc, const vec<double> &ySrc, const vec<double> &sSrc, const vec<double> &fSrc,
                             const vec<vec<int>> &nghListSrc, const vec<vec<int>> &nghListSelf)
{
    /** LSMPS B Note :
     *   > There are two relative position in this operation
     *     * The target interpolated position (refer as TARGET)
     *     * The position of interpolation data (refer as SOURCE)
     *   > Both TARGET and SOURCE is particle variable
     *   
     *   TARGET VARIABLE - interpolation target
     *     vec<dbl> xT, vec<dbl> yT  : Particle position list
     *     vec<dbl> sT               : Particle size list
     *     vec<dbl> fT               : Interpolated variable data (Literally not used in this calculation)
     *     vec<vec<int>> nghListSrc  : Target particle neighbor ID list consisted of source particle
     *     vec<vec<int>> nghListSelf : Target particle neighbor ID list consisted of it self
     *   
     *   SOURCE VARIABLE - interpolation data source
     *     vec<dbl> xS, vec<dbl> yS    : Particle position list
     *     vec<dbl> sS                 : Particle size list
     *     vec<dbl> fS                 : Interpolated variable data
    */
    
    // Take the number of particle
    int nTarget = sTar.size();
    
    // Internal variable
    int nClear = 0;
    int nLow = 0;

    // vec<double> averageDiameter(nTarget, 0.0e0);

    // Evaluate the LSMPS B of all target particle
    #pragma omp parallel for /*reduction(+:nClear)*/
    for (int i = 0; i < nTarget; i++)
    {
        // Evaluate the neighbor of the evaluated GRID
        const vec<int> &_nghSrcID = nghListSrc[i];
        const vec<int> &_nghSelfID = nghListSelf[i];

        if (_nghSrcID.empty()){
            // If there is no nearby particle from the target, the target differential value set to zero
            this->_d0[i]     = 0.0e0;
            this->_ddx[i]    = 0.0e0;
            this->_ddy[i]    = 0.0e0;
            this->_d2d2x[i]  = 0.0e0;
            this->_d2dxdy[i] = 0.0e0;
            this->_d2d2y[i]  = 0.0e0;
            // this->_d3d3x[i]  = 0.0e0;
        }
        else{
            // // [OLD CODE]
            // // TODO: calculate average diameter
            // for (size_t j = 0; j < _nghSrcID.size(); j++)
            // {
            //     int idx = _nghSrcID[j];
            //     averageDiameter[i] += si[idx];
            // }
            
            // // Why using this as r_s <?>
            // averageDiameter[i] = averageDiameter[i] / (double)_nghSrcID.size();

            // // Use the GRID size
            // averageDiameter[i] = s[i];
            
            // Initialize the LSMPS matrices
            MatrixXd Hbar = MatrixXd::Zero(BASE_SIZE_2D, BASE_SIZE_2D);   // Scaling matrix
            MatrixXd Mbar = MatrixXd::Zero(BASE_SIZE_2D, BASE_SIZE_2D);   // Moment matrix
            VectorXd bbar = VectorXd::Zero(BASE_SIZE_2D);              // Moment vector

            // Scaling size (using the target particle size)
            double _rs = sTar[i];
            
            // Initialize H_rs Matrix
            Hbar(0, 0) = 1;
            Hbar(1, 1) = std::pow(_rs, -1);
            Hbar(2, 2) = std::pow(_rs, -1);
            Hbar(3, 3) = std::pow(_rs, -2) * 2;
            Hbar(4, 4) = std::pow(_rs, -2);
            Hbar(5, 5) = std::pow(_rs, -2) * 2;
            //Hbar(6, 6) = std::pow(_rs, -3) * 6;

            // Calculate the moment matrix [Mbar] and vector [bbar] (Evaluate toward neighbor source particle)
            for (size_t j = 0; j < _nghSrcID.size(); j++)
            {
                // Aliasing the particle ID
                const int &idxi = i;                // Target particle ID
                const int &idxj = _nghSrcID[j];     // Neighor particle ID
                
                // Evaluated particle coordinate (TARGET particle)
                double _xi = xTar[idxi]; 
                double _yi = yTar[idxi];
                
                // Neighbor particle coordinate (SOURCE particle)
                double _xj = xSrc[idxj];
                double _yj = ySrc[idxj];
                
                // Coordinate distance between TARGET and SOURCE
                double _dx = _xj - _xi;       // x distance between TARGET[i] and neighbor: SOURCE[j]
                double _dy = _yj - _yi;       // y distance between TARGET[i] and neighbor: SOURCE[j]

                // Resultant distance(_rij) and effective radius(_Re)
                double _rij = std::sqrt(std::pow(_dx, 2) + std::pow(_dy, 2));
                double _Re = 0.0;
                switch (INTERACTION_TYPE){
                case 1: // Effective radius type 1 (self particle size)
                    _Re = RADIUS_FACTOR_2D * sTar[idxi];
                    break;
                case 2: // Effective radius type 2 (average particle size)
                    _Re = RADIUS_FACTOR_2D * (sTar[idxi] + sSrc[idxj]) * 0.5;
                    break;
                default:
                    break;
                }

                // No effect from the current neighbor (or neighbor outside the support domain)
                if (_rij > _Re) continue;

                // Calculate the weight value
                double _wij = this->weight_function(_rij, _Re);     // Calculate the weight
                _wij *= std::pow(sSrc[idxj]/sTar[idxi],2.0);        // Add the relative size factor between particle
                
                // // [OLD DATA]
                // double _Ri = RADIUS_FACTOR_2D * averageDiameter[idxi];   // Effective radius of GRID particle
                // // double _Rj = RADIUS_FACTOR_2D * si[idxj];                // Effective radius of PARTICLE particle
                // // double _Rij = (_Ri + _Rj) * 0.5;                    // Effective radius of Average particle
                
                // Calculate the distance vector p
                const vec<double> _P = this->get_p(_dx, _dy, _rs);

                // LSMPS B interpolated variable value SOURCE
                double _fij = fSrc[idxj] - 0.0e0; 

                // Calculation of moment matrix Mbar and bbar
                for (size_t k1 = 0; k1 < BASE_SIZE_2D; k1++){
                    for (size_t k2 = 0; k2 < BASE_SIZE_2D; k2++){
                        // generate tensor product between p
                        Mbar(k1, k2) = Mbar(k1, k2) + (_wij * _P[k1] * _P[k2]);
                    }
                    // generate moment matrix
                    bbar(k1) = bbar(k1) + (_wij * _P[k1] * _fij);
                }
            }
            
            // TODO: Increase robustness by performing small determinant treatment [considering the self-particle neighbor]
            if (std::abs(Mbar.determinant()) < 1.0e-8) // ! if tend to be NON-invertible matrix
            {
                // FOR DEBUGGING (Warning)
                nClear --;
                nLow ++;

                // Add the element of moment matrix [Mbar] and vector [bbar] with self neighbor data
                for (size_t j = 0; j < _nghSelfID.size(); j++)
                {
                    // Aliasing the particle ID
                    const int &idxi = i;                // Target particle ID
                    const int &idxj = _nghSelfID[j];    // Neighor particle ID
                    
                    // Evaluated particle coordinate (TARGET particle)
                    double _xi = xTar[idxi]; 
                    double _yi = yTar[idxi];
                    
                    // Neighbor particle coordinate (toward TARGET particle itself)
                    double _xj = xTar[idxj];
                    double _yj = yTar[idxj];

                    // Coordinate distance between TARGET particle and TARGET particle neighbor
                    double _dx = _xj - _xi;       // x distance between TARGET[i] and neighbor: TARGET[j]
                    double _dy = _yj - _yi;       // y distance between TARGET[i] and neighbor: TARGET[j]

                    // Resultant distance(_rij) and effective radius(_Re)
                    double _rij = std::sqrt(std::pow(_dx, 2) + std::pow(_dy, 2));
                    double _Re = 0.0;
                    switch (INTERACTION_TYPE){
                    case 1: // Effective radius type 1 (self particle size)
                        _Re = RADIUS_FACTOR_2D * sTar[idxi];
                        break;
                    case 2: // Effective radius type 2 (average particle size)
                        _Re = RADIUS_FACTOR_2D * (sTar[idxi] + sTar[idxj]) * 0.5;
                        break;
                    default:
                        break;
                    }

                    // No effect from the current neighbor (or neighbor outside the support domain)
                    if (_rij > _Re) continue;

                    // Calculate the weight value
                    double _wij = this->weight_function(_rij, _Re);     // Calculate the weight
                    _wij *= std::pow(sTar[idxj]/sTar[idxi],2.0);        // Add the relative size factor between particle

                    // Calculate the distance vector p
                    const vec<double> _P = this->get_p(_dx, _dy, _rs);

                    // Interpolated variable of neighbor: GRID
                    double _fij = 0.0e0;           // ! no interaction of the function, value set to 0.0, then no need to do LSMPS<?>

                    // Calculate the moment matrix Mbar and bbar
                    for (size_t k1 = 0; k1 < BASE_SIZE_2D; k1++){
                        for (size_t k2 = 0; k2 < BASE_SIZE_2D; k2++){
                            // generate tensor product between p
                            Mbar(k1, k2) = Mbar(k1, k2) + (_wij * _P[k1] * _P[k2]);
                        }
                        // generate moment matrix
                        bbar(k1) = bbar(k1) + (_wij * _P[k1] * _fij);
                    }
                }
            }

            // Solve Least Square Matrix Operation
            // --- Method 1 ---     LU Decomposition
            VectorXd MbarInv_Bbar = Mbar.lu().solve(bbar);
            // --- Method 2 ---     House Holder QR 
            // VectorXd MbarInv_Bbar = Mbar.fullPivHouseholderQr().solve(bbar);
            // --- Method 3 ---     Singular Value Decomposition
            // VectorXd MbarInv_Bbar = Mbar.bdcSvd(ComputeThinU | ComputeThinV).solve(bbar);
            
            // Final result
            VectorXd Df = Hbar * MbarInv_Bbar; 

            // Assign the LSMPS B result
            this->_d0[i]     = Df[0];
            this->_ddx[i]    = Df[1];
            this->_ddy[i]    = Df[2];
            this->_d2d2x[i]  = Df[3];
            this->_d2dxdy[i] = Df[4];
            this->_d2d2y[i]  = Df[5];
            //this->_d3d3x[i] = Dx[6];

            nClear++;
        }
    }

    // DEBUGGING
    // std::cout << "\n\033[36m[LOG] \033[0m" << "LSMPS B interpolation Summary!\n";
    // printf("<+> Number of clear interpolation       : %8d\n", nClear);
    // printf("<+> Number of low determinant           : %8d\n\n", nLow);
}


/**
 *  @brief A 3D LSMPS type B function differential calculation of a particle data set 
 *  toward a source data set.
 *  
 *  @param  xTarget Target particle x coordinate data.
 *  @param  yTarget Target particle y coordinate data.
 *  @param  zTarget Target particle z coordinate data.
 *  @param  sTarget Target particle size data.
 *  @param  fTarget Target particle function value data.
 *  @param  xSource Source particle x coordinate data.
 *  @param  ySource Source particle y coordinate data.
 *  @param  zSource Source particle z coordinate data.
 *  @param  sSource Source particle size data.
 *  @param  fSource Source particle function value data.
 *  @param  neighborListSource The neighbor relation of target particle data toward source particle.
 *  @param  neighborListSelf The neighbor relation of target particle data toward itself.
*/
void LSMPSb::calculate_LSMPS_3D(const vec<double> &xTar, const vec<double> &yTar, const vec<double> &zTar,
                                const vec<double> &sTar, const vec<double> &fTar,
                                const vec<double> &xSrc, const vec<double> &ySrc, const vec<double> &zSrc,
                                const vec<double> &sSrc, const vec<double> &fSrc,
                                const vec<vec<int>> &nghListSrc, const vec<vec<int>> &nghListSelf)
{
    // Take the number of particle
    int nTarget = sTar.size();
    
    // Internal variable
    int nClear = 0;
    int nLow = 0;

    // Evaluate the LSMPS B of all target particle
    #pragma omp parallel for /*reduction(+:nClear)*/
    for (int i = 0; i < nTarget; i++)
    {
        // Evaluate the neighbor of the evaluated GRID
        const vec<int> &_nghSrcID = nghListSrc[i];
        const vec<int> &_nghSelfID = nghListSelf[i];

        if (!_nghSrcID.empty()){
            // Initialize the LSMPS matrices
            MatrixXd Hbar = MatrixXd::Zero(BASE_SIZE_3D, BASE_SIZE_3D);   // Scaling matrix
            MatrixXd Mbar = MatrixXd::Zero(BASE_SIZE_3D, BASE_SIZE_3D);   // Moment matrix
            VectorXd bbar = VectorXd::Zero(BASE_SIZE_3D);              // Moment vector

            // Scaling size (using the target particle size)
            double _rs = sTar[i];

            // Initialize H_rs Matrix: (rs^-|a|)*(a!)
            Hbar(0, 0) = 1;
            Hbar(1, 1) = std::pow(_rs, -1);       //  x(1,0,0)
            Hbar(2, 2) = std::pow(_rs, -1);       //  y(0,1,0)
            Hbar(3, 3) = std::pow(_rs, -1);       //  z(0,0,1)
            Hbar(4, 4) = std::pow(_rs, -2) * 2;   // x2(2,0,0)
            Hbar(5, 5) = std::pow(_rs, -2);       // xy(1,1,0)
            Hbar(6, 6) = std::pow(_rs, -2);       // xz(1,0,1)
            Hbar(7, 7) = std::pow(_rs, -2) * 2;   // y2(0,2,0)
            Hbar(8, 8) = std::pow(_rs, -2);       // yz(0,1,1)
            Hbar(9, 9) = std::pow(_rs, -2) * 2;   // z2(0,0,2)

            // Calculate the moment matrix [Mbar] and vector [bbar] (Evaluate toward neighbor source particle)
            for (size_t j = 0; j < _nghSrcID.size(); j++)
            {
                // Aliasing the particle ID
                const int &idxi = i;                // Target particle ID
                const int &idxj = _nghSrcID[j];     // Neighor particle ID
                
                // Evaluated particle coordinate (TARGET particle)
                double _xi = xTar[idxi]; 
                double _yi = yTar[idxi];
                double _zi = zTar[idxi];
                
                // Neighbor particle coordinate (SOURCE particle)
                double _xj = xSrc[idxj];
                double _yj = ySrc[idxj];
                double _zj = zSrc[idxj];
                
                // Coordinate distance between TARGET and SOURCE
                double _dx = _xj - _xi;       // x distance between TARGET[i] and neighbor: SOURCE[j]
                double _dy = _yj - _yi;       // y distance between TARGET[i] and neighbor: SOURCE[j]
                double _dz = _zj - _zi;       // z distance between TARGET[i] and neighbor: SOURCE[j]

                // Resultant distance(_rij) and effective radius(_Re)
                double _rij = std::sqrt(std::pow(_dx, 2) + std::pow(_dy, 2) + std::pow(_dz, 2));
                double _Re = 0.0;
                switch (INTERACTION_TYPE){
                case 1: // Effective radius type 1 (self particle size)
                    _Re = RADIUS_FACTOR_3D * sTar[idxi];
                    break;
                case 2: // Effective radius type 2 (average particle size)
                    _Re = RADIUS_FACTOR_3D * (sTar[idxi] + sSrc[idxj]) * 0.5;
                    break;
                default:
                    break;
                }

                // No effect from the current neighbor (or neighbor outside the support domain)
                if (_rij > _Re) continue;

                // Calculate the weight value
                double _wij = this->weight_function(_rij, _Re);     // Calculate the weight
                _wij *= std::pow(sSrc[idxj]/sTar[idxi], 3.0);       // Add the relative size factor between particle (volume)
                
                // Calculate the distance vector p
                vec<double> _P;
                this->get_p_3d(_P, _dx, _dy, _dz, _rs);

                // LSMPS B interpolated variable value SOURCE
                double _fij = fSrc[idxj]; 

                // Calculation of moment matrix Mbar and bbar
                for (size_t k1 = 0; k1 < BASE_SIZE_3D; k1++){
                    for (size_t k2 = 0; k2 < BASE_SIZE_3D; k2++){
                        // generate tensor product between p
                        Mbar(k1, k2) = Mbar(k1, k2) + (_wij * _P[k1] * _P[k2]);
                    }
                    // generate moment matrix
                    bbar(k1) = bbar(k1) + (_wij * _P[k1] * _fij);
                }
            }
            
            // TODO: Increase robustness by performing small determinant treatment [considering the self-particle neighbor]
            if (std::abs(Mbar.determinant()) < 1.0e-8) // ! if tend to be NON-invertible matrix
            {
                // FOR DEBUGGING (Warning)
                nClear --;
                nLow ++;

                // Add the element of moment matrix [Mbar] and vector [bbar] with self neighbor data
                for (size_t j = 0; j < _nghSelfID.size(); j++)
                {
                    // Aliasing the particle ID
                    const int &idxi = i;                // Target particle ID
                    const int &idxj = _nghSelfID[j];    // Neighor particle ID
                    
                    // Evaluated particle coordinate (TARGET particle)
                    double _xi = xTar[idxi]; 
                    double _yi = yTar[idxi];
                    double _zi = zTar[idxi];
                    
                    // Neighbor particle coordinate (toward TARGET particle itself)
                    double _xj = xTar[idxj];
                    double _yj = yTar[idxj];
                    double _zj = zTar[idxj];

                    // Coordinate distance between TARGET particle and TARGET particle neighbor
                    double _dx = _xj - _xi;       // x distance between TARGET[i] and neighbor: TARGET[j]
                    double _dy = _yj - _yi;       // y distance between TARGET[i] and neighbor: TARGET[j]
                    double _dz = _zj - _zi;       // z distance between TARGET[i] and neighbor: TARGET[j]

                    // Resultant distance(_rij) and effective radius(_Re)
                    double _rij = std::sqrt(std::pow(_dx, 2) + std::pow(_dy, 2) + std::pow(_dz, 2));
                    double _Re = 0.0;
                    switch (INTERACTION_TYPE){
                    case 1: // Effective radius type 1 (self particle size)
                        _Re = RADIUS_FACTOR_3D * sTar[idxi];
                        break;
                    case 2: // Effective radius type 2 (average particle size)
                        _Re = RADIUS_FACTOR_3D * (sTar[idxi] + sTar[idxj]) * 0.5;
                        break;
                    default:
                        break;
                    }

                    // No effect from the current neighbor (or neighbor outside the support domain)
                    if (_rij > _Re) continue;

                    // Calculate the weight value
                    double _wij = this->weight_function(_rij, _Re);     // Calculate the weight
                    _wij *= std::pow(sTar[idxj]/sTar[idxi], 3.0);       // Add the relative size factor between particle (volume)

                    // Calculate the distance vector p
                    vec<double> _P;
                    this->get_p_3d(_P, _dx, _dy, _dz, _rs);

                    // Interpolated variable of neighbor: GRID
                    double _fij = 0.0e0;           // ! no interaction of the function, value set to 0.0, then no need to do LSMPS<?>

                    // Calculate the moment matrix Mbar and bbar
                    for (size_t k1 = 0; k1 < BASE_SIZE_3D; k1++){
                        for (size_t k2 = 0; k2 < BASE_SIZE_3D; k2++){
                            // generate tensor product between p
                            Mbar(k1, k2) = Mbar(k1, k2) + (_wij * _P[k1] * _P[k2]);
                        }
                        // generate moment matrix
                        bbar(k1) = bbar(k1) + (_wij * _P[k1] * _fij);
                    }
                }
            }

            // Solve Least Square Matrix Operation by LU Decomposition
            VectorXd MbarInv_Bbar = Mbar.lu().solve(bbar);
            
            // Final result
            VectorXd Df = Hbar * MbarInv_Bbar; 

            // Assign the LSMPS B result
            this->_d0[i]     = Df[0];    //  0(0,0,0)
            this->_ddx[i]    = Df[1];    //  x(1,0,0)
            this->_ddy[i]    = Df[2];    //  y(0,1,0)
            this->_ddz[i]    = Df[3];    //  z(0,0,1)
            this->_d2d2x[i]  = Df[4];    // x2(2,0,0)
            this->_d2dxdy[i] = Df[5];    // xy(1,1,0)
            this->_d2dxdz[i] = Df[6];    // xz(1,0,1)
            this->_d2d2y[i]  = Df[7];    // y2(0,2,0)
            this->_d2dydz[i] = Df[8];    // yz(0,1,1)
            this->_d2d2z[i]  = Df[9];    // z2(0,0,2)

            nClear++;   
        }
        else{
            // If there is no nearby particle from the target, the target differential value set to zero
            this->_d0[i]     = 0.0e0;
            this->_ddx[i]    = 0.0e0;
            this->_ddy[i]    = 0.0e0;
            this->_ddz[i]    = 0.0e0;
            this->_d2d2x[i]  = 0.0e0;
            this->_d2dxdy[i] = 0.0e0;
            this->_d2dxdz[i] = 0.0e0;
            this->_d2d2y[i]  = 0.0e0;
            this->_d2dydz[i] = 0.0e0;
            this->_d2d2z[i]  = 0.0e0;
        }
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
vec<double> LSMPSb::get_p(const double &xij, const double &yij, const double &rs)
{
    vec<double> _p(BASE_SIZE_2D);

    double _xij = xij / rs;
    double _yij = yij / rs;

    _p[0] = 1;
    _p[1] = _xij;
    _p[2] = _yij;
    _p[3] = _xij * _xij;
    _p[4] = _xij * _yij;
    _p[5] = _yij * _yij;

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
void LSMPSb::get_p_3d(std::vector<double> &_P, const double dx, 
                      const double dy, const double dz, const double rs)
{
    _P.resize(BASE_SIZE_3D);

    double _dx = dx / rs;
    double _dy = dy / rs;
    double _dz = dz / rs;

    _P[0] = 1;
    _P[1] = _dx;            //  x(1,0,0)
    _P[2] = _dy;            //  y(0,1,0)
    _P[3] = _dz;            //  z(0,0,1)
    _P[4] = _dx * _dx;      // x2(2,0,0)
    _P[5] = _dx * _dy;      // xy(1,1,0)
    _P[6] = _dx * _dz;      // xz(1,0,1)
    _P[7] = _dy * _dy;      // y2(0,2,0)
    _P[8] = _dy * _dz;      // yz(0,1,1)
    _P[9] = _dz * _dz;      // z2(0,0,2)

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
double LSMPSb::weight_function(const double &_rij, const double &_Re)
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
 *  @brief Get the first x derivative (df/dx).
 *  @return The df/dx derivative.
 */
vec<double> LSMPSb::get_d0(){
    return this->_d0;
}

/**
 *  @brief Get the first x derivative (df/dx).
 *  @return The df/dx derivative.
 */
vec<double> LSMPSb::get_ddx(){
    return this->_ddx;
}

/**
 *  @brief Get the first y derivative (df/dy).
 *  @return The df/dy derivative.
 */
vec<double> LSMPSb::get_ddy(){
    return this->_ddy;
}

/**
 *  @brief Get the first z derivative (df/dz).
 *  @return The df/dz derivative.
 */
vec<double> LSMPSb::get_ddz(){
    return this->_ddz;
}

/**
 *  @brief Get the second x derivative (d2f/dx2).
 *  @return The d2f/dx2 derivative.
 */
vec<double> LSMPSb::get_d2d2x(){
    return this->_d2d2x;
}

/**
 *  @brief Get the first x and y derivative (d2f/dxdy).
 *  @return The d2f/dxdy derivative.
 */
vec<double> LSMPSb::get_d2dxdy(){
    return this->_d2dxdy;
}

/**
 *  @brief Get the first x and z derivative (d2f/dxdz).
 *  @return The d2f/dxdz derivative.
 */
vec<double> LSMPSb::get_d2dxdz(){
    return this->_d2dxdz;
}

/**
 *  @brief Get the second y derivative (d2f/dy2).
 *  @return The d2f/dy2 derivative.
 */
vec<double> LSMPSb::get_d2d2y(){
    return this->_d2d2y;
}

/**
 *  @brief Get the first y and z derivative (d2f/dydz).
 *  @return The d2f/dydz derivative.
 */
vec<double> LSMPSb::get_d2dydz(){
    return this->_d2dydz;
}

/**
 *  @brief Get the second y derivative (d2f/dz2).
 *  @return The d2f/dz2 derivative.
 */
vec<double> LSMPSb::get_d2d2z(){
    return this->_d2d2z;
}

// #pragma endregion