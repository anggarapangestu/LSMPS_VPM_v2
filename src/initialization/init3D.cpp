#include "initialization.hpp"

#define MARGIN_EXTENSION_FACTOR 1.1

// =========================================
// ========= Particle Distribution =========
// =========================================

// Initialize the data of coordinate [x,y,z], size and level [s, level], and number [num]

/**
 *  @brief  Uniform single resolution distribution for 3D simulation.
 *         
 *  @param  _particle  Particle data container.
*/
void initialization::init_3d_single_res(Particle &par)
{
    // Number of particle in x and y direction
    int _nx = std::ceil(Pars::lxdom / Pars::sigma);     // Number of particle in x direction
    int _ny = std::ceil(Pars::lydom / Pars::sigma);     // Number of particle in y direction
    int _nz = std::ceil(Pars::lzdom / Pars::sigma);     // Number of particle in y direction
    int _nparticle = _nx*_ny*_nz;                       // Number of all particle

    // Adjust the size of vector into the calculated particle number
    par.x.resize(_nparticle, 0.0e0);
    par.y.resize(_nparticle, 0.0e0);
    par.z.resize(_nparticle, 0.0e0);

    // Simulation domain pivot coordinate
    double pivX = -Pars::xdom;              // Location of the most left (x)
    double pivY = -Pars::sigma*_ny*0.5;     // Location of the most bottom (y)
    double pivZ = -Pars::sigma*_nz*0.5;     // Location of the most back (z)

    // Generate the particle position
    double _y, _x, _z;
    int count = 0;
    for (int i = 0; i < _nx; i++){
        for (int j = 0; j < _ny; j++){
            for (int k = 0; k < _nz; k++){
                _x = pivX + (i + 0.5)*Pars::sigma;
                _y = pivY + (j + 0.5)*Pars::sigma;
                _z = pivZ + (k + 0.5)*Pars::sigma;

                // Additional of unbalance for early separation (mp_shift)
                if (Pars::flag_slightly_shifted_domain){
                    _y += Pars::mp_shift;
                    _z += Pars::mp_shift;
                }
                
                // Assign the particle
                par.x[count] = _x;
                par.y[count] = _y;
                par.z[count] = _z;
                count ++;
            }
        }
    }
    
    // Assign other particle properties
    par.num = _nparticle;
    par.s.resize(_nparticle, Pars::sigma);
    par.level.resize(_nparticle, Pars::max_level);

}

/**
 *  @brief  Two level resolution with finer resolution particle in single block 
 *  distribution for 3D simulation.
 *         
 *  @param  _particle  Particle data container.
 *  @param  _bodyList  The list of body data container used for particle generation.
*/
void initialization::init_3d_multi_res_single_block(Particle &par, const std::vector<Body> &bL)
{
    // internal temporary variable
    double _z, _y, _x;

    // ====== PARTICLE GENERATION ====== //
    // number of particle in x and y direction
    int _nx = std::ceil(0.5 * Pars::lxdom / Pars::sigma);   // Number of particle in x direction
    int _ny = std::ceil(0.5 * Pars::lydom / Pars::sigma);   // Number of particle in y direction
    int _nz = std::ceil(0.5 * Pars::lzdom / Pars::sigma);   // Number of particle in z direction
    int _nparticle = _nx * _ny * _nz;                       // Number of all particle

    // Adjust the size of vector into the calculated particle number
    par.x.resize(_nparticle, 0.0e0);
    par.y.resize(_nparticle, 0.0e0);
    par.z.resize(_nparticle, 0.0e0);
    par.s.resize(_nparticle, 2 * Pars::sigma);
    par.level.resize(_nparticle, Pars::max_level - 1);
    par.num = _nparticle;

    // Simulation domain pivot coordinate
    double pivX = -Pars::xdom;              // Location of the most left (x)
    double pivY = -Pars::sigma*_ny;         // Location of the most bottom (y)
    double pivZ = -Pars::sigma*_nz;         // Location of the most bottom (y)

    // Generate the particle position
    int count = 0;
    for (int i = 0; i < _nx; i++){
        for (int j = 0; j < _ny; j++){
            for (int k = 0; k < _ny; k++){
                _x = pivX + (i+0.5)*2*Pars::sigma;
                _y = pivY + (j+0.5)*2*Pars::sigma;
                _z = pivZ + (k+0.5)*2*Pars::sigma;

                // Additional of unbalance for early separation (mp_shift)
                if (Pars::flag_slightly_shifted_domain){
                    _y += Pars::mp_shift;
                    _z += Pars::mp_shift;
                }

                par.x[count] = _x;
                par.y[count] = _y;
                par.z[count] = _z;
                count ++;
            }
        }
    }
    
    // ====== DIVIDE AND CONQUEROR PROCEDURE ====== //
    // Put extension from the body limit to the finer resolution block limit (USER INPUT MANUALY)
    double ext_front = 0.75e0;  // Upstream extension
    double ext_back  = 1.50e0;  // Downstream extension
    double ext_side  = 0.75e0;  // Side body extension

    // Perform the devide and conqueror if the point inside the finer resolution block
    int chdTrans[DIM];
    
    // Iterate through all current particle
    for (int i = 0; i < _nparticle; i++){
        // Check the particle position related to the body
        bool skipDnC = true;
        for (int part = 0; part < N_BODY; part++){
            // Check the location toward upstream
            if (par.x[i] > bL[part].min_pos[0] - ext_front &&   // Check toward body upstream
                par.x[i] < bL[part].max_pos[0] + ext_back &&    // Check toward body downstream
                par.y[i] > bL[part].min_pos[1] - ext_side &&    // Check toward body bottom side
                par.y[i] < bL[part].max_pos[1] + ext_side &&    // Check toward body upper side
                par.z[i] > bL[part].min_pos[2] - ext_side &&    // Check toward body left side
                par.z[i] < bL[part].max_pos[2] + ext_side ){    // Check toward body right side
                // Particle near this body, dont skip the DnC on this particle
                skipDnC = false;
                break;
            }
        }
        if (skipDnC) continue;
        
        // ** Proceed the Divide and Conqueror Process        
        // Take the parent particle coordinate
        _x = par.x[i];
        _y = par.y[i];
        _z = par.z[i];
        
        // Point 1  -> Replace the parent data
        par.x[i] = _x - 0.5 * Pars::sigma;
        par.y[i] = _y - 0.5 * Pars::sigma;
        par.z[i] = _z - 0.5 * Pars::sigma;
        par.s[i] = Pars::sigma;
        par.level[i] = Pars::max_level;

        // Point 2-4  -> Add new particle data
        for (int j = 1; j < 8; j++){
            basis_loop(d) chdTrans[d] = (2*((j>>d)%2)) - 1;
            par.x.push_back(_x + 0.5 * Pars::sigma * chdTrans[0]);
            par.y.push_back(_y + 0.5 * Pars::sigma * chdTrans[1]);
            par.z.push_back(_z + 0.5 * Pars::sigma * chdTrans[2]);
            par.s.push_back(Pars::sigma);
            par.level.push_back(Pars::max_level);
            par.num++;
        }
    }
}



/**
 *  @brief  Multiresolution body adjusted distribution for 3D simulation.
 *         
 *  @param  _particle  Particle data container.
 *  @param  _bodyList  The list of body data container used for particle generation.
*/
void initialization::init_3d_multi_res_body_adjusted(Particle &par, const std::vector<Body> &bL)
{    
    // internal temporary variable
    double _z, _y, _x;
    const double baseSize = Pars::intPow(2,Pars::max_level) * Pars::sigma;
    geometry geom_tool;

    // ====== PARTICLE GENERATION ====== //
    // number of particle in x and y direction
    int _nx = std::ceil(Pars::lxdom / baseSize);    // Number of particle in x direction
    int _ny = std::ceil(Pars::lydom / baseSize);    // Number of particle in y direction
    int _nz = std::ceil(Pars::lzdom / baseSize);    // Number of particle in z direction
    int _nparticle = _nx * _ny * _nz;               // Number of all particle

    // Adjust the size of vector into the calculated particle number
    par.x.resize(_nparticle, 0.0e0);
    par.y.resize(_nparticle, 0.0e0);
    par.z.resize(_nparticle, 0.0e0);
    par.s.resize(_nparticle, baseSize);
    par.level.resize(_nparticle, 0);
    par.num = _nparticle;

    // Simulation domain pivot coordinate
    double pivX = -Pars::xdom;          // Location of the most left (x)
    double pivY = -baseSize*_ny*0.5;    // Location of the most bottom (y)
    double pivZ = -baseSize*_nz*0.5;    // Location of the most back (z)

    // Generate the particle position
    int count = 0;
    for (int i = 0; i < _nx; i++){
        for (int j = 0; j < _ny; j++){
            for (int k = 0; k < _nz; k++){
                _x = pivX + (i + 0.5)*baseSize;
                _y = pivY + (j + 0.5)*baseSize;
                _z = pivZ + (k + 0.5)*baseSize;

                // Additional of unbalance for early separation (mp_shift)
                if (Pars::flag_slightly_shifted_domain){
                    _y += Pars::mp_shift;
                    _z += Pars::mp_shift;
                }

                par.x[count] = _x;
                par.y[count] = _y;
                par.z[count] = _z;
                count ++;
            }
        }
    }

    // ====== DIVIDE AND CONQUEROR PROCEDURE ====== //
    // Perform the devide and conqueror if the point inside the finer resolution block
    int chdTrans[DIM];
    
    // Define a queue list for recursive particle refinement
    std::vector<int> par_list;
    std::vector<int> par_list_next;
    // Initialize the queue list
    par_list.resize(par.num);
    for (int i = 0; i < par.num; i++) par_list[i] = i;

    // Constant parameter for DnC evaluation
    std::vector<double> pos;
    double R_marg;
    double R_eval;

    // Iterate through all levels
    for (int lvl = 0; lvl < Pars::max_level; lvl++){
        // Arithmatic Series [Sn = n/2*(2*a + b*(n-1))]
        int n_S = Pars::max_level - lvl;    // The position order
        double a_Sn = 15 * Pars::sigma;//0.2;//15 * Pars::sigma;     // Initial unit (Manual Adjusted)
        double b_Sn = 8 * Pars::sigma;//0.15;//8 * Pars::sigma;      // Increment (Manual adjuster)
        R_marg = n_S * (a_Sn + 0.5 * b_Sn * (n_S - 1)) * Pars::Df;

        double currSize = baseSize/Pars::intPow(2,lvl); // Size of particle in the current level

        // Iterate through all particle in the queue
        par_list_next.clear();  // Reserve the next particle container
        for (const auto &ID : par_list){
            // Check the particle position related to the body (only for level 0)
            if (lvl == 0){
                bool skipDnC = true;
                for (int part = 0; part < N_BODY; part++){
                    // Check the location toward upstream
                    if (par.x[ID] > bL[part].min_pos[0] - R_marg*MARGIN_EXTENSION_FACTOR &&    // Check toward body upstream
                        par.x[ID] < bL[part].max_pos[0] + R_marg*MARGIN_EXTENSION_FACTOR &&    // Check toward body downstream
                        par.y[ID] > bL[part].min_pos[1] - R_marg*MARGIN_EXTENSION_FACTOR &&    // Check toward body bottom side
                        par.y[ID] < bL[part].max_pos[1] + R_marg*MARGIN_EXTENSION_FACTOR &&    // Check toward body upper side
                        par.z[ID] > bL[part].min_pos[2] - R_marg*MARGIN_EXTENSION_FACTOR &&    // Check toward body left side
                        par.z[ID] < bL[part].max_pos[2] + R_marg*MARGIN_EXTENSION_FACTOR ){    // Check toward body right side
                        // Particle near this body, dont skip the DnC on this particle
                        skipDnC = false;
                        break;
                    }
                }
                if (skipDnC) continue;
            }

            // Calculate the particle minimum distance toward body
            pos = {par.x[ID], par.y[ID], par.z[ID]};
            R_eval = geom_tool.distance_calc(pos, bL[0], true);
            for (int part = 1; part < N_BODY; part++){
                double _R = geom_tool.distance_calc(pos, bL[part], true);
                if (_R < R_eval) R_eval = _R;
            }

            // DIVIDE AND CONQUEROR
            if (std::abs(R_eval) < R_marg){
                // Take the coordinate of parent particle
                _x = par.x[ID];
                _y = par.y[ID];
                _z = par.z[ID];
                
                // Point 1
                par.x[ID] = _x - 0.25 * currSize;
                par.y[ID] = _y - 0.25 * currSize;
                par.z[ID] = _z - 0.25 * currSize;
                par.s[ID] = par.s[ID]/2;
                par.level[ID] = par.level[ID] + 1;
                par_list_next.push_back(ID);

                // Point 2 - 8
                for (int j = 1; j < 8; j++){
                    basis_loop(d) chdTrans[d] = (2*((j>>d)%2)) - 1;
                    par_list_next.push_back(par.s.size());
                    par.x.push_back(_x + 0.25 * currSize * chdTrans[0]);
                    par.y.push_back(_y + 0.25 * currSize * chdTrans[1]);
                    par.z.push_back(_z + 0.25 * currSize * chdTrans[2]);
                    par.s.push_back(par.s[ID]);
                    par.level.push_back(par.level[ID]);
                    par.num ++;
                }
            }
        }
        
        // Update the particle list
        par_list = par_list_next;
    }
}