#include "base_remeshing.hpp"

// ===================== Sign of particles ====================
//    Finding:
//    - the particle inside the body+mollification region (1) or outside the body(0) --> sign_par
//    - the normal/midpoint distance of particle and nearest triangle face --> testpin
//    Method:
//     Finding the nearest triangle face to the considered point (close with mid-point of the face, not normal distance)
//     and then give the distance by using 02 ways
//    - Using distance from considered point to midpoint of triangle face
//    - Using normal distance from considered point to triangle face
//  input: surface of body: vertex, triangular faces
//         position of particles
//  output: sign of particles: =1 inside or lie on surface of particles, lie on mollification region
//                             =0 outside of body    --> sign_par
//          distance  --> testpin
// ============================================================

void base_remeshing::sign_particles(const Body &b, int npi, vector<double> &xpi, vector<double> &ypi, vector<bool> &sign_par, vector<bool> &inside)
{
    // internal variables
    int k1, npanel;
    double rlim1, rlim2, deltax, deltay, r_normal, r_normal_neighbor, rmins, sign_n, _theta;

    double *xcen = new double[b.x.size() - 1];
    double *ycen = new double[b.x.size() - 1];
    double *normal_panelx = new double[b.x.size() - 1];
    double *normal_panely = new double[b.x.size() - 1];
    double *r2panel = new double[b.x.size() - 1];
    double *testpin = new double[npi];

    double *xmin_pen = new double[2];
    double *xmax_pen = new double[2];

    std::vector<double> _xp(npi);
    std::vector<double> _yp(npi);
    std::vector<double> _idx(npi);

    // == limiting particles near body only =================================
    std::fill(sign_par.begin(), sign_par.end(), 0.0e0);
    // ========================================
    // ----- Determine Preliminary "minimum and maximum" grid position
    //  because maybe the body size is not fit to "body-penalization nodes/ except 6 additional nodes around pen-body" (biger or smaller)
    //  to make sure the 6 additional nodes will be outside of body we need pre-calculation:
    //   . . . . . . |     '     | . . . . .
    //   . . . . . . |     '     | . . . . .
    //   . . . . . . |     '     | . . . . .
    //   . . . . . . |     '     | . . . . .
    //                    center
    //awalnya 

    // ?? Supposedly using the body not the particle (?)
    xmin_pen[0] = *min_element(xpi.begin(), xpi.end()); - 10 * Pars::sigma;
    xmin_pen[1] = *min_element(ypi.begin(), ypi.end()); - 10 * Pars::sigma;
    xmax_pen[0] = *max_element(xpi.begin(), xpi.end()); + 10 * Pars::sigma;
    xmax_pen[1] = *max_element(ypi.begin(), ypi.end()); + 10 * Pars::sigma;

    int _np = 0; // initial
    for (size_t i = 0; i < npi; i++)
    {
        // Limiting the Particles close to grid boundaries, due to effect's taking penalization
        if ((xpi[i] >= xmin_pen[0]) && (xpi[i] <= xmax_pen[0]) &&
            (ypi[i] >= xmin_pen[1]) && (ypi[i] <= xmax_pen[1]))
        {
            _xp[_np] = xpi[i];
            _yp[_np] = ypi[i];
            _idx[_np] = i;
            _np++;
        }
    }
    // =======================================================================

    printf("Distinguishing particle's sign/distance ... \n");
    // ==============================
    npanel = b.x.size() - 1;
    for (int i = 0; i < npanel; i++)
    {
        // midpoint of panel
        xcen[i] = b.x[i] + ((b.x[i + 1] - b.x[i]) * 0.5e0);
        ycen[i] = b.y[i] + ((b.y[i + 1] - b.y[i]) * 0.5e0);
        // normal panel
        deltax = b.x[i + 1] - b.x[i];
        deltay = b.y[i + 1] - b.y[i];
        _theta = atan2(deltay, deltax);
        normal_panelx[i] = sin(_theta);  // unit normal vector in x-direction
        normal_panely[i] = -cos(_theta); // unit normal vector in y-direction
    }

    k1 = 0;
    for (int i = 0; i < _np; i++)
    {
        rmins = 100.0e0; // initial guessed value, must be set to be large enough
        for (int k = 0; k < npanel; k++)
        {
            r2panel[k] = std::sqrt(std::pow((_xp[i] - xcen[k]), 2) + std::pow((_yp[i] - ycen[k]), 2)); // Distance between center panel and the point
            if (std::abs(r2panel[k]) <= rmins)
            {
                rmins = r2panel[k]; // find the center panel closest to the point
                k1 = k;             // Name of closest panel
            }
        }

        // normal distance from considered point to panel
        r_normal = (((_xp[i] - xcen[k1]) * (normal_panelx[k1])) + ((_yp[i] - ycen[k1]) * (normal_panely[k1]))); // =0 when grid node lie on the panel, or extended-panel

        if (r_normal == 0.0e0)
            sign_n = 0.0e0; // grid node lie on the panel, or extended-panel
        else
            sign_n = r_normal / (std::abs(r_normal));

        // ==== Midpoint Distance: using distance from considered point to midpoint of triangle face
        // read more .../Report_progress/April2017_Airfoil_Penalization_inOutside_problem.odp, part B
        testpin[i] = rmins * sign_n;
        // ==== Normal Distance: Using normal distance from considered point to triangle face
        // testpin(i)=r_normal;
        
        
        if (testpin[i] <= -10 * Pars::sigma)    // Must check the body generator first for the sign of the tespin
        {
            sign_par[_idx[i]] = false;  // far outside the body
        }
        else
        {
            sign_par[_idx[i]] = true;   // inside body + mollification region
        }

        //jika -xp[i] diluar "batas" yang ditentukan
        //"badtas" saat ini center domain +- 2.5
       // if((_xp[i] >= 0 - 7.5) && (_xp[i] <= 0 + 7.5)){
        //sign_par[_idx[i]] = true;
        //}

    }

    // == clear internal variables
    delete[] xcen;
    delete[] ycen;
    delete[] normal_panelx;
    delete[] normal_panely;
    delete[] r2panel;
    delete[] testpin;
}
