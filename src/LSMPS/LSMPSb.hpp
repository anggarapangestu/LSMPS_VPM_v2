#ifndef INCLUDED_LSMPSb
#define INCLUDED_LSMPSb

#include <iostream>
#include <vector>
#include <fstream>
#include <omp.h>

/**
 *  @brief A particle derivative function of Least Square Moving Particle 
 *  Semi-Implitit Type B. The Type B provide a calculation from the zero
 *  derivative level. This method can be utilized as function interpolation. 
 *  This class only invoke until the second derivative calculation.
 *
 *  @headerfile LSMPSb.hpp
 */
class LSMPSb
{
private:
    // Internal Variable
    std::vector<double> _d0;        // Zeroth derivative ({})
    std::vector<double> _ddx;       // First derivative toward x (d{}/dx)
    std::vector<double> _ddy;       // First derivative toward y (d{}/dy)
    std::vector<double> _d2d2x;     // Second derivative toward x (d^2{}/d^2x)
    std::vector<double> _d2dxdy;    // First derivative toward x and y (d^2{}/dxdy)
    std::vector<double> _d2d2y;     // Second derivative toward y (d^2{}/d^2y)
    
    // Additional to 3D
    std::vector<double> _ddz;       // First derivative toward z (d{}/dz)
    std::vector<double> _d2dxdz;    // First derivative toward x and z (d^2{}/dxdz)
    std::vector<double> _d2dydz;    // First derivative toward y and z (d^2{}/dydz)
    std::vector<double> _d2d2z;     // Second derivative toward z (d^2{}/d^2z)

    // LSMPS B Utilities

    std::vector<double> get_p(const double &dx, const double &dy, const double &rs);
    void get_p_3d(std::vector<double> &_P, const double dx, const double dy, const double dz, const double rs);
    double weight_function(const double &_rij, const double &_Re);

    // LSMPS B calculation 2D

    void calculate_LSMPS_2D(const std::vector<double> &xTarget, const std::vector<double> &yTarget,
                            const std::vector<double> &sTarget,
                            const std::vector<double> &xSource, const std::vector<double> &ySource,
                            const std::vector<double> &sSource, const std::vector<double> &fSource,
                            const std::vector<std::vector<int>> &neighborListSource);
    
    void calculate_LSMPS(const std::vector<double> &xTarget, const std::vector<double> &yTarget,
                         const std::vector<double> &sTarget, const std::vector<double> &fTarget,
                         const std::vector<double> &xSource, const std::vector<double> &ySource,
                         const std::vector<double> &sSource, const std::vector<double> &fSource,
                         const std::vector<std::vector<int>> &neighborListSource,
                         const std::vector<std::vector<int>> &neighborListSelf);

    // LSMPS B Calculation 3D

    void calculate_LSMPS_3D(const std::vector<double> &xTarget, const std::vector<double> &yTarget, const std::vector<double> &zTarget,
                            const std::vector<double> &sTarget, const std::vector<double> &fTarget,
                            const std::vector<double> &xCollocation, const std::vector<double> &yCollocation, const std::vector<double> &zCollocation,
                            const std::vector<double> &sCollocation, const std::vector<double> &fCollocation,
                            const std::vector<std::vector<int>> &neighborListSource,
                            const std::vector<std::vector<int>> &neighborListSelf);
public:
    // LSMPS B manager

    void set_LSMPS_2D(const std::vector<double> &xTarget, const std::vector<double> &yTarget,
                      const std::vector<double> &sTarget,
                      const std::vector<double> &xSource, const std::vector<double> &ySource,
                      const std::vector<double> &sSource, const std::vector<double> &fSource,
                      const std::vector<std::vector<int>> &neighborListSource);

    void set_LSMPS(const std::vector<double> &xTarget, const std::vector<double> &yTarget,
                   const std::vector<double> &sTarget, const std::vector<double> &fTarget,
                   const std::vector<double> &xSource, const std::vector<double> &ySource,
                   const std::vector<double> &sSource, const std::vector<double> &fSource,
                   const std::vector<std::vector<int>> &neighborListSource,
                   const std::vector<std::vector<int>> &neighborListSelf);

    void set_LSMPS_3D(const std::vector<double> &xTarget, const std::vector<double> &yTarget, const std::vector<double> &zTarget,
                      const std::vector<double> &sTarget, const std::vector<double> &fTarget,
                      const std::vector<double> &xSource, const std::vector<double> &ySource, const std::vector<double> &zSource,
                      const std::vector<double> &sSource, const std::vector<double> &fSource,
                      const std::vector<std::vector<int>> &neighborListSource,
                      const std::vector<std::vector<int>> &neighborListSelf);

    // Getter function on spatial differential

    std::vector<double> get_d0();
    std::vector<double> get_ddx();
    std::vector<double> get_ddy();
    std::vector<double> get_ddz();
    std::vector<double> get_d2d2x();
    std::vector<double> get_d2dxdy();
    std::vector<double> get_d2dxdz();
    std::vector<double> get_d2d2y();
    std::vector<double> get_d2dydz();
    std::vector<double> get_d2d2z();
};
#endif