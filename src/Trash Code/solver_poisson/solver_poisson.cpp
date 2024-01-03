#include "solver.hpp"

#pragma region public methods
void solverpoisson::set_solver( vector<double> &u, vector<double> &v, const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s,
                   const std::vector<double> &f, std::vector<std::vector<int>> &neighborlist, vector<int> isboundary, vector<double> boundaryvalue, int it)
{
    if(it == 0){
        this->calculate_eta(this->dx2, this->dy2, x, y, s, neighborlist);

        this->create_global_matrix(neighborlist, this->dx2, this->dy2, isboundary, boundaryvalue, f, s, x.size());
        
        this->dx2.clear();
        this->dy2.clear();
    } 

    vector<double> psi;
    psi.resize(x.size(), 0.0);
    this->solve_poisson(psi, isboundary, boundaryvalue, f, 2, x.size());   // 2 is the physical core

    //================Calculate derivative of psi============================
    // clear local variables
    this->_ddx.clear();
    this->_ddy.clear();
    this->_d2d2x.clear();
    this->_d2dxdy.clear();
    this->_d2d2y.clear();

    // resize local variabels
    int nparticle = x.size();
    this->_d.resize(nparticle);
    this->_ddx.resize(nparticle);
    this->_ddy.resize(nparticle);
    this->_d2d2x.resize(nparticle);
    this->_d2dxdy.resize(nparticle);
    this->_d2d2y.resize(nparticle);

    this-> LSMPS_calc(x,y,s,psi,neighborlist);

    //update velocity
    for (size_t i = 0; i < x.size(); i++){
        u[i] = this->_ddy[i];
        v[i] = -this->_ddx[i];
    }


}
#pragma endregion