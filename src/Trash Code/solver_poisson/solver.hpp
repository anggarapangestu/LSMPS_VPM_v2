#ifndef INCLUDED_SOLVERPOISSON
#define INCLUDED_SOLVERPOISSON
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <omp.h>          // pragma omp ...
#include "../Eigen/Dense" // linear algebra library
#include "../Eigen/Sparse"
using namespace std; 
using namespace Eigen;    // for linear algebra operation

//// template <typename T>
class solverpoisson
{
//Poisson solver using LSMPS type B as the base method. 

private:
    // #define MAT_SIZE 5 // size of matrix
    static const int MAT_SIZE = 6; // size of matrix
    
    SparseMatrix<double, RowMajor> A;
    std::vector<vector<double>> dx2, dy2;

    std::vector<double> _d;      // d{}
    std::vector<double> _ddx;    // d{}/dx
    std::vector<double> _ddy;    // d{}/dy
    std::vector<double> _d2d2x;  // d^2{}/d^2x
    std::vector<double> _d2dxdy; // d^2{}/dxdy
    std::vector<double> _d2d2y;  // d^2{}/d^2y

    void LSMPS_calc(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s,
                         const std::vector<double> &f, std::vector<std::vector<int>> &neighborlist);
    
    void create_global_matrix(vector<vector<int>> neighbour, vector<vector<double>> EtaDx2, vector<vector<double>> EtaDy2,
                    vector<int> isboundary, vector<double> boundaryvalue, vector<double> RHS, vector<double> e,
                    int part_num);

    void solve_poisson(vector<double>&u, vector<int> isboundary, vector<double> boundaryvalue, vector<double> RHS, int thread, int part_num);

    void calculate_eta(vector<vector<double>> &LSMPS_EtaDx2, vector<vector<double>> &LSMPS_EtaDy2, vector<double> x, vector<double> y,
                       vector<double> s, vector<vector<int>> neighbour);
public:

    void set_solver(vector<double> &u, vector<double> &v, const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s,
                   const std::vector<double> &f, std::vector<std::vector<int>> &neighborlist, vector<int> isboundary, vector<double> boundaryvalue, int it);

};
#endif