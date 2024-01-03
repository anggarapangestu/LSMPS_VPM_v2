#ifndef INCLUDED_LSMPSa
#define INCLUDED_LSMPSa

#include <iostream>
#include <vector>
#include <fstream>
#include <omp.h>

/**
 *  @brief A particle derivative function of Least Square Moving Particle 
 *  Semi-Implitit Type A. The Type A provide a calculation from the first
 *  derivative level. This class only invoke until the second derivative
 *  calculation.
 *
 *  @headerfile LSMPSa.hpp
 */
class LSMPSa
{
private:
    // Internal variables
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

    // LSMPS A Utilities

    std::vector<double> get_p(const double dx, const double dy, const double rs);
    void get_p_3d(std::vector<double> &_P, const double dx, const double dy, const double dz, const double rs);
    double weight_function(const double &_rij, const double &_Re);

    // LSMPS A Calculation

    void calculate_LSMPS(const std::vector<double> &x, const std::vector<double> &y,
                         const std::vector<double> &s, const std::vector<double> &f,
                         const std::vector<std::vector<int>> &neighborList);
    void calculate_LSMPS(const std::vector<double> &xTarget, const std::vector<double> &yTarget,
                         const std::vector<double> &sTarget, const std::vector<double> &fTarget,
                         const std::vector<double> &xCollocation, const std::vector<double> &yCollocation,
                         const std::vector<double> &sCollocation, const std::vector<double> &fCollocation,
                         const std::vector<std::vector<int>> &neighborList);

    // LSMPS A Calculation 3D
    void calculate_LSMPS_3D(const std::vector<double> &xTarget, const std::vector<double> &yTarget, const std::vector<double> &zTarget,
                            const std::vector<double> &sTarget, const std::vector<double> &fTarget,
                            const std::vector<double> &xCollocation, const std::vector<double> &yCollocation, const std::vector<double> &zCollocation,
                            const std::vector<double> &sCollocation, const std::vector<double> &fCollocation,
                            const std::vector<std::vector<int>> &neighborList);
public:
    // LSMPS A manager

    void set_LSMPS(const std::vector<double> &x, const std::vector<double> &y, 
                   const std::vector<double> &s, const std::vector<double> &f,
                   const std::vector<std::vector<int>> &neighborList);
    void set_LSMPS(const std::vector<double> &xTarget, const std::vector<double> &yTarget,
                   const std::vector<double> &sTarget, const std::vector<double> &fTarget,
                   const std::vector<double> &xCollocation, const std::vector<double> &yCollocation,
                   const std::vector<double> &sCollocation, const std::vector<double> &fCollocation,
                   const std::vector<std::vector<int>> &neighborList);

    void set_LSMPS_3D(const std::vector<double> &xTarget, const std::vector<double> &yTarget, const std::vector<double> &zTarget,
                      const std::vector<double> &sTarget, const std::vector<double> &fTarget,
                      const std::vector<double> &xCollocation, const std::vector<double> &yCollocation, const std::vector<double> &zCollocation,
                      const std::vector<double> &sCollocation, const std::vector<double> &fCollocation,
                      const std::vector<std::vector<int>> &neighborList);

    // Getter function on spatial differential

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