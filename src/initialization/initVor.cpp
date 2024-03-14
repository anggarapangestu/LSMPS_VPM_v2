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
     * 
     *  The perlman stream function field formula
     *              _
     *             /  (-4r^2 + 7r^4 - (28/3)r^6 + (35/4)r^8 - (28/5)r^10  ; r <= 1
     *   phi(r) = <      + (7/3)r^12 - (4/7)r^14 + (1/16)r^16) / 16
     *             \_   -(ln(r)+1.3589285714)/16    ; r > 1
     * 
     * 
     *  The velocity differential to Perlman Vorticity
     *  >> x velocity    _
     *                  /   xy/(8*r^4) * (1 - (1+7r^2)*(1-r^2)^7)  ; r <= 1
     *     dudx(x,y) = <
     *                  \_  xy/(8*r^4)                             ; r > 1
     *                   _
     *                  /   -1/(16*r^4) * (r^2-2y^2 + (14y^2r^2-r^2+r^4+2y^2)(1-r^2)^7)  ; r <= 1
     *     dudy(x,y) = <
     *                  \_  (y^2-x^2)/(16*r^4)                    ; r > 1
     * 
     *  >> y velocity  _
     *                  /   1/(16*r^4) * (r^2-2x^2 + (14x^2r^2-r^2+r^4+2x^2)(1-r^2)^7)  ; r <= 1
     *     dvdy(x,y) = <
     *                  \_  (y^2-x^2)/(16*r^4)                     ; r > 1
     *                   _
     *                  /   -xy/(8*r^4) * (1 - (1+7r^2)*(1-r^2)^7)  ; r <= 1
     *     dvdy(x,y) = <
     *                  \_  -xy/(8*r^4)                    ; r > 1
    */

    // Temporary save the analytical velocity into gx,gy,gz
    par.gx = std::vector<double> (par.num, 0.0e0);
    par.gy = std::vector<double> (par.num, 0.0e0);
    // par.gz = std::vector<double> (par.num, 0.0e0);
    par.chi = std::vector<double> (par.num, 0.0e0);  // Temporary store the stream function


    par.dudx = std::vector<double> (par.num, 0.0e0);
    par.dudy = std::vector<double> (par.num, 0.0e0);
    par.dvdx = std::vector<double> (par.num, 0.0e0);
    par.dvdy = std::vector<double> (par.num, 0.0e0);

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

        // Calculate the perlman velocity
        if (r2 <= 1){
            par.gx[i] = -par.y[i] / (16.0*r2) * (1 - std::pow((1.0 - r2), 8.0));
            par.gy[i] =  par.x[i] / (16.0*r2) * (1 - std::pow((1.0 - r2), 8.0));
        }else{
            par.gx[i] = -par.y[i] / (16.0*r2);
            par.gy[i] =  par.x[i] / (16.0*r2);
        }

        // Calculate the perlman stream function
        if (r2 <= 1){
            // double r4  = std::pow(r2,2);
            // double r6  = std::pow(r2,3);
            // double r8  = std::pow(r2,4);
            // double r10 = std::pow(r2,5);
            // double r12 = std::pow(r2,6);
            // double r14 = std::pow(r2,7);
            // double r16 = std::pow(r2,8);
            par.chi[i] = (1.0/16.0) * (- 4*r2 + 7*std::pow(r2,2)
                                   - ( 28.0 / 3.0  ) * std::pow(r2,3)
                                   + ( 35.0 / 4.0  ) * std::pow(r2,4)
                                   - ( 28.0 / 5.0  ) * std::pow(r2,5)
                                   + ( 7.0  / 3.0  ) * std::pow(r2,6)
                                   - ( 4.0  / 7.0  ) * std::pow(r2,7)
                                   + ( 1.0  / 16.0 ) * std::pow(r2,8)
                                   );
        }else{
            par.chi[i] = -(1.3589285714 + (std::log(r2)/2.0)) / 16.0;
        }

        // Calculate the perlman velocity differential
        if (r2 <= 1){
            // Alias the coordinate
            const double &x = par.x[i];
            const double &y = par.y[i];
            const double r4 = r2*r2;

            par.dudx[i] =  x*y/(8.0*r4) * (1 - (1 + 7*r2)*std::pow((1.0 - r2), 7.0));
            par.dudy[i] = -1.0/(16.0*r4) * (r2 - 2*y*y + (14*y*y*r2 - r2 + r4 + 2*y*y)*std::pow((1.0 - r2), 7.0));
            par.dvdx[i] =  1.0/(16.0*r4) * (r2 - 2*x*x + (14*x*x*r2 - r2 + r4 + 2*x*x)*std::pow((1.0 - r2), 7.0));
            par.dvdy[i] = -par.dudx[i];
        }else{
            // Alias the coordinate
            const double &x = par.x[i];
            const double &y = par.y[i];
            const double r4 = r2*r2;

            par.dudx[i] =  x*y/(8.0*r4);
            par.dudy[i] = (y*y - x*x)/(16.0*r4);
            par.dvdx[i] =  par.dudy[i];
            par.dvdy[i] = -par.dudx[i];
        }
    }
    return;
}

// Elliptical vorticity
void initialization::eliptic_vorticity(Particle &par){
    // Internal parameter
    double R0 = 0.8;        // Base radius
    double Vm = 20.0;       // Maximum vorticity
    double q0 = 2.56085;    // Smoothing constant

    // Reserve the vorticity data
    par.vorticity.resize(par.num);

    #pragma omp_parallel_for
    for (int i = 0; i < par.num; i++){
        // Calculate the polar coordinate of radius
        double r2 = 0;
        r2 += 4*(par.x[i]*par.x[i]);    // Component x
        r2 += (par.y[i]*par.y[i]);      // Component y

        // Calculate the other parameter
        double r = std::sqrt(r2);       // The elliptical radius
        double z = r/R0;               // The radius ratio

        if (z < 1.0){
            double f_q = std::exp(-(q0/z) * std::exp(1/(z-1)));
            par.vorticity[i] = Vm * f_q;
        }else{
            par.vorticity[i] = 0.0;
        }

        // #if (DIM == 3)
        // r2 += (par.z[i]*par.z[i]);  // Component z
        // #endif
    }


    return;
}

// ==================================================
// +--------------- TESTING FUNCTION ---------------+
// ==================================================
/** All testing function store the data in 2D
 *  - Coordinate at x, y
 *  - Function in   chi (analytical)
 *  - Source in     vorticity
*/

// GAUSSIAN FUNCTION
void initialization::test_0(Particle &par){
    // A gaussian function: phi(x,y) = exp (-(x^2 + y^2))   [Solution]
    // A gaussian function: p(x,y) = -4*(1-(x^2 + y^2)) * exp(-(x^2 + y^2))     [Source]
    // Reserve the container
    par.chi = std::vector<double>(par.num, 0.0e0);          // The potential
    par.vorticity = std::vector<double>(par.num, 0.0e0);    // The source of poisson

    par.gx = std::vector<double>(par.num, 0.0e0);       // Diferential toward x
    par.gy = std::vector<double>(par.num, 0.0e0);       // Diferential toward y

    // Calculate the function and the source
    #pragma omp parallel for
    for (int i = 0; i < par.num; i++){
        // Calculate the polar coordinate of radius
        double r2 = 0;
        r2 += (par.x[i]*par.x[i]);  // Component x
        r2 += (par.y[i]*par.y[i]);  // Component y

        // Calculate the perlman vorticity
        double phi = std::exp(-r2);
        par.chi[i] = phi;
        par.vorticity[i] = -4 * (1 - r2) * phi;

        // For neuman boundary condition
        par.gx[i] = -2*par.x[i]*phi;
        par.gy[i] = -2*par.y[i]*phi;
    }
    return;
}

// POLYNOMIAL FUNCTION
void initialization::test_1(Particle &par){
    // A simple polynomial function: phi(x,y) = (x^2 * y^2) [Solution]
    // A simple polynomial function: f(x,y) = 2(x^2 + y^2)  [Source]
    // Reserve the container
    par.chi = std::vector<double>(par.num, 0.0e0);          // The potential
    par.vorticity = std::vector<double>(par.num, 0.0e0);    // The source of poisson

    par.gx = std::vector<double>(par.num, 0.0e0);       // Diferential toward x
    par.gy = std::vector<double>(par.num, 0.0e0);       // Diferential toward y

    // Calculate the function and the source
    #pragma omp parallel for
    for (int i = 0; i < par.num; i++){
        // Alias the value
        const double &x = par.x[i];
        const double &y = par.y[i];

        // Calculate the perlman vorticity
        par.chi[i] = std::pow(x*y, 2);
        par.vorticity[i] = 2 * (x*x + y*y);
        
        // For neuman boundary condition
        par.gx[i] = 2*x*y*y;
        par.gy[i] = 2*x*x*y;
    }
    return;
}

// MORE COMPLICATED FUNCTION
void initialization::test_2(Particle &par){
    // A trigonometric polynomial function: phi(x,y) = 2(x^2)sin^2(y) + 4(x^2)y + 3(y^2)cos(2x)     [Solution]
    // A trigonometric polynomial function: f(x,y) = 4sin^2(y) + 4(x^2)cos(2y) + 8y + (6-12y^2)cos(2x)  [Source]
    // Reserve the container
    par.chi = std::vector<double>(par.num, 0.0e0);          // The potential
    par.vorticity = std::vector<double>(par.num, 0.0e0);    // The source of poisson

    par.gx = std::vector<double>(par.num, 0.0e0);       // Diferential toward x
    par.gy = std::vector<double>(par.num, 0.0e0);       // Diferential toward y

    // Calculate the function and the source
    #pragma omp parallel for
    for (int i = 0; i < par.num; i++){
        // Alias the value
        const double &x = par.x[i];
        const double &y = par.y[i];

        // Calculate the perlman vorticity
        par.chi[i] = 2*std::pow(x*sin(y), 2) + 4*x*x*y + 3*y*y*cos(2*x);
        par.vorticity[i] = 4*std::pow(sin(y), 2) + 4*x*x*cos(2*y) + 8*y + (6-12*y*y)*cos(2*x);

        // For neuman boundary condition
        par.gx[i] = 4*x*std::pow(sin(y), 2) + 8*x*y - 6*y*y*sin(2*x);
        par.gy[i] = 2*x*x*sin(2*y) + 4*x*x + 6*y*cos(2*x);
    }
    return;
}

// Laplace solution
void initialization::laplace(Particle &par){
    // The basic solution through the laplace problem 
    // V(x,y) = 2*V_0/pi * arctan(sin(pi*(x/L)) / sinh (pi*(y/L)));

    // An internal variable for the laplace solution
    double L = Pars::lxdom;         // The domain side size
    double sx = 0.5*Pars::lxdom;        // The x boundary shift
    double sy = 0.5*Pars::lydom;        // The y boundary shift
    double V_0 = 10.0;       // The potential at the left side (x=-sx)
    int _P = 100;
    double pi = M_PI;       // circle constant

    // Reserve the container
    par.chi = std::vector<double>(par.num, 0.0e0);
    
    // Calculate the analytical solution of potential
    #pragma omp parallel for
    for (int i = 0; i < par.num; i++){
        // Alias the value
        double x = par.x[i] + sx;
        double y = par.y[i] + sy;
        // double num = std::sin(pi*(y/L));
        // double den = std::sinh(pi*(x/L));
        // double V = 2*V_0/pi * std::atan(num / den);
        // if (std::isnan(V)){
        //     V = 2*V_0/pi * std::atan(1.0);
        // }

        double V = 0.0;
        for (int n = 1; n < _P; n++){
            V += (2*V_0/(n*M_PI))*((1-std::cos(M_PI *n))/std::sinh(n*M_PI))*
                 (std::sin(n*M_PI*y/L)*std::sinh(n*M_PI*x/L));
        } 

        // Calculate the analytical solution
        par.chi[i] = V;
    }
    return;
}