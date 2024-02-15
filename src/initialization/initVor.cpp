#include "initialization.hpp"

/**
 *  @brief  Generate perlman vorticity as initial value.
 *         
 *  @param  _particle  Particle data for vorticity calculation.
 * 
 *  @return No return.
*/
void initialization::perlman_vorticity(Particle &par){
    /** Perlman Vorticity Field Formula
     *            _
     *           /  (1-r^2)^7  ; r <= 1
     *   w(r) = <
     *           \_   0        ; r > 1
     * 
     *  The velocity solution to Perlman Vorticity
     *  >> x velocity  _
     *               /   -y/(16*r^2) * (1 - (1-r^2)^8)  ; r <= 1
     *     u(x,y) = <
     *               \_  -y/(16*r^2)                    ; r > 1
     * 
     *  >> y velocity  _
     *               /   x/(16*r^2) * (1 - (1-r^2)^8)  ; r <= 1
     *     v(x,y) = <
     *               \_  x/(16*r^2)                    ; r > 1
    */

    // Reserve the vorticity container
    par.vorticity.clear();
    par.vorticity.resize(par.num, 0.0e0);

    // Calculate the vorticity at each particle
    #pragma omp parallel for
    for (int i = 0; i < par.num; i++){
        // Calculate the polar coordinate of radius
        double r2 = 0;
        r2 += (par.x[i]*par.x[i]);  // Component x
        r2 += (par.y[i]*par.y[i]);  // Component y
        #if (DIM == 3)
        r2 += (par.z[i]*par.z[i]);  // Component z
        #endif

        // Calculate the perlman vorticity
        if (r2 <= 1){
            par.vorticity[i] = std::pow((1.0 - r2), 7.0);
        }else{
            par.vorticity[i] = 0;
        }
    }
    return;
}


/**
 *  @brief  Calculate the analytical solution of perlman vorticity.
 *         
 *  @param  _particle  Particle data for velocity calculation.
 * 
 *  @return No return.
*/
void initialization::perlman_velocity_solution(Particle &par){
    /** Perlman Vorticity Field Formula
     *            _
     *           /  (1-r^2)^7  ; r <= 1
     *   w(r) = <
     *           \_   0        ; r > 1
     * 
     *  The velocity solution to Perlman Vorticity
     *  >> x velocity  _
     *               /   -y/(16*r^2) * (1 - (1-r^2)^8)  ; r <= 1
     *     u(x,y) = <
     *               \_  -y/(16*r^2)                    ; r > 1
     * 
     *  >> y velocity  _
     *               /   x/(16*r^2) * (1 - (1-r^2)^8)  ; r <= 1
     *     v(x,y) = <
     *               \_  x/(16*r^2)                    ; r > 1
    */

    // Temporary save the analytical velocity into gx,gy,gz
    par.gx = std::vector<double> (par.num, 0.0e0);
    par.gy = std::vector<double> (par.num, 0.0e0);
    // par.gz = std::vector<double> (par.num, 0.0e0);

    // Calculate the velocity solution at each particle
    #pragma omp parallel for
    for (int i = 0; i < par.num; i++){
        // Calculate the polar coordinate of radius
        double r2 = 0;
        r2 += (par.x[i]*par.x[i]);  // Component x
        r2 += (par.y[i]*par.y[i]);  // Component y
        // #if (DIM == 3)
        // r2 += (par.z[i]*par.z[i]);  // Component z
        // #endif

        // Calculate the perlman vorticity
        if (r2 <= 1){
            par.gx[i] = -par.y[i] / (16*r2) * (1 - std::pow((1.0 - r2), 8.0));
            par.gy[i] =  par.x[i] / (16*r2) * (1 - std::pow((1.0 - r2), 8.0));
        }else{
            par.gx[i] = -par.y[i] / (16*r2);
            par.gy[i] =  par.x[i] / (16*r2);
        }
    }
    return;
}