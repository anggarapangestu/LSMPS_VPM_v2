#include "initialization.hpp"
#define MARGIN_EXTENSION_FACTOR 1.1

// =========================================
// ========= Particle Distribution =========
// =========================================

// Initialize the data of coordinate [x,y,z], size and level [s, level], and number [num]

/**
 *  @brief  A testing particle distribution. A resolution of 2 level that split the 
 *  domain into two region whether in x or y direction.
 *         
 *  @param  _particle  Particle data container.
*/
void initialization::init_2d_test_1(Particle &par){
    // internal temporary variable
    double _y, _x;

    // ====== PARTICLE GENERATION ====== //
    // number of particle in x and y direction
    int _nx = std::ceil(0.5 * Pars::lxdom / Pars::sigma);   // Number of particle in x direction
    int _ny = std::ceil(0.5 * Pars::lydom / Pars::sigma);   // Number of particle in y direction
    int _nparticle = _nx * _ny;                             // Number of all particle

    // Adjust the size of vector into the calculated particle number
    par.x.resize(_nparticle, 0.0e0);
    par.y.resize(_nparticle, 0.0e0);
    par.s.resize(_nparticle, 2 * Pars::sigma);
    par.level.resize(_nparticle, Pars::max_level - 1);
    par.num = _nparticle;

    // Simulation domain pivot coordinate
    double pivX = -Pars::xdom;              // Location of the most left (x)
    double pivY = -Pars::sigma*_ny;         // Location of the most bottom (y)

    // Generate the particle position
    int count = 0;
    for (int i = 0; i < _nx; i++){
        for (int j = 0; j < _ny; j++){
            _x = pivX + (i+0.5)*2*Pars::sigma;
            _y = pivY + (j+0.5)*2*Pars::sigma + Pars::mp_shift;   // Additional of unbalance: Pars::mp_shift
            par.x[count] = _x;
            par.y[count] = _y;
            count ++;
        }
    }
    
    // ====== DIVIDE AND CONQUEROR PROCEDURE ====== //
    // Set the particle division type
    int type = 0;               // 0:= split in x direction; 1:= split in y direction
    bool up_region = true;      // 0:= divide the lower coordinate region; 1:= divide the higher coordinate region
    double split_loc = -2.6;    // The location of splitting

    // Perform the devide and conqueror if the point inside the finer resolution block
    int chdTrans[DIM];
    
    // Iterate through all current particle
    for (int i = 0; i < _nparticle; i++){
        // Check if the particle need to be devided
        if (type == 0){
            // Split in x direction
            if (up_region && par.x[i] < split_loc) continue;
            else if (!up_region && par.x[i] > split_loc) continue;
        }else if (type == 1){
            // Split in y direction
            if (up_region && par.y[i] < split_loc) continue;
            else if (!up_region && par.y[i] > split_loc) continue;
        }
        
        // ** Proceed the Divide and Conqueror Process        
        // Take the parent particle coordinate
        _x = par.x[i];
        _y = par.y[i];
        
        // Point 1  -> Replace the parent data
        par.x[i] = _x - 0.5 * Pars::sigma;
        par.y[i] = _y - 0.5 * Pars::sigma;
        par.s[i] = Pars::sigma;
        par.level[i] = Pars::max_level;

        // Point 2-4  -> Add new particle data
        for (int j = 1; j < 4; j++){
            basis_loop(d) chdTrans[d] = (2*((j>>d)%2)) - 1;
            par.x.push_back(_x + 0.5 * Pars::sigma * chdTrans[0]);
            par.y.push_back(_y + 0.5 * Pars::sigma * chdTrans[1]);
            par.s.push_back(Pars::sigma);
            par.level.push_back(Pars::max_level);
            par.num++;
        }
    }
    
    return;
}

/**
 *  @brief  A testing particle distribution. A resolution of 2 level that split the 
 *  domain into two region whether in x or y direction.
 *  NOTE: Similar result with test 1 but different data ordering
 *         
 *  @param  _particle  Particle data container.
*/
void initialization::init_2d_test_2(Particle &par){
    // internal temporary variable
    double _y, _x;

    // ====== PARTICLE GENERATION ====== //
    // Set the particle division type
    int type = 1;               // 0:= split in x direction; 1:= split in y direction
    bool up_region = true;      // 0:= divide the lower coordinate region; 1:= divide the higher coordinate region
    bool splitted;              // Flag that denote whether the domain is splitted or not
    double split_loc = 0.1;    // The location of splitting (origin at O(0,0))

    // Split into two region
    // Check if the particle need to be devided
    if (type == 0){
        // Split in x direction
        if (-Pars::xdom < split_loc && Pars::lxdom - Pars::xdom > split_loc){
            // The split location is inside the 
            splitted = true;
        }
        else{
            // The split location is outside the region
            splitted = false;
        }
    }else if (type == 1){
        // Split in y direction
        if (-Pars::lydom/2 < split_loc && Pars::lydom/2 > split_loc){
            // The split location is inside the 
            splitted = true;
        }
        else{
            // The split location is outside the region
            splitted = false;
        }
    }
    
    // Divert in two type of particle arrangement
    if (!splitted){
        // The type one : only one resolution
        std::cout << "|O| Only a single resolution\n";
        this->init_2d_single_res(par);
    }
    else if (splitted){
        // The type two : splitted domain
        std::cout << "|A| A splitted domain distribution\n";
        
        // Number of particle in x and y direction
        int _nx, _ny;       // Number of particle in x and y direction
        int _nparticle;     // Number of all particle
        
        // Simulation domain pivot coordinate
        double pivX, pivY;  // Pivot location of the most left (x) and the most bottom (y)
        double size;

        // [1] Generate the first region
        // *****************************
        size = up_region ? (2*Pars::sigma) : (Pars::sigma);
        if (type == 0){
            // Splitted in x direction
            // Update the number of particle in x and y direction
            _nx = std::ceil((split_loc + Pars::xdom) / size);
            _ny = std::ceil(Pars::lydom / size);
            _nparticle = _nx * _ny;
            
            // Simulation domain pivot coordinate
            pivX = split_loc - (_nx*size);   // Location of the most left (x)
            pivY = - (_ny*size) / 2.0;       // Location of the most bottom (y)
        }
        else if (type == 1){
            // Splitted in y direction
            // Update the number of particle in x and y direction
            _nx = std::ceil(Pars::lxdom / size);
            _ny = std::ceil((split_loc + (Pars::lydom / 2.0)) / size);
            _nparticle = _nx * _ny;
            
            // Simulation domain pivot coordinate
            pivX = -Pars::xdom;              // Location of the most left (x)
            pivY = split_loc - (_ny*size);   // Location of the most bottom (y)
        }

        // Adjust the size of vector into the calculated particle number
        par.num = _nparticle;
        par.x.resize(par.num, 0.0e0);
        par.y.resize(par.num, 0.0e0);
        par.s.resize(par.num, size);
        par.level.resize(par.num, up_region ? Pars::max_level - 1 : Pars::max_level);

        // Generate the particle position
        int count = 0;
        for (int i = 0; i < _nx; i++){
            for (int j = 0; j < _ny; j++){
                _x = pivX + (i+0.5)*size;
                _y = pivY + (j+0.5)*size + Pars::mp_shift;   // Additional of unbalance: Pars::mp_shift
                par.x[count] = _x;
                par.y[count] = _y;
                count++;
            }
        }
        std::cout << "|X| Done generating the first (1st) part\n";


        // [2] Generate the second region
        // *****************************
        size = up_region ? Pars::sigma : 2*Pars::sigma;
        if (type == 0){
            // Splitted in x direction
            // Update the number of particle in x and y direction
            _nx = std::ceil((Pars::lxdom - split_loc - Pars::xdom) / size);
            _ny = std::ceil(Pars::lydom / size);
            _nparticle = _nx * _ny;
            
            // Simulation domain pivot coordinate
            pivX = split_loc;   // Location of the most left (x)
            pivY = - (_ny*size) / 2.0;       // Location of the most bottom (y)
        }
        else if (type == 1){
            // Splitted in y direction
            // Update the number of particle in x and y direction
            _nx = std::ceil(Pars::lxdom / size);
            _ny = std::ceil(((Pars::lydom / 2.0) - split_loc) / size);
            _nparticle = _nx * _ny;
            
            // Simulation domain pivot coordinate
            pivX = -Pars::xdom;              // Location of the most left (x)
            pivY = split_loc;   // Location of the most bottom (y)
        }

        // Adjust the size of vector into the calculated particle number
        par.num += _nparticle;
        par.x.resize(par.num, 0.0e0);
        par.y.resize(par.num, 0.0e0);
        par.s.resize(par.num, size);
        par.level.resize(par.num, up_region ? Pars::max_level : Pars::max_level - 1);

        // Generate the particle position
        for (int i = 0; i < _nx; i++){
            for (int j = 0; j < _ny; j++){
                _x = pivX + (i+0.5)*size;
                _y = pivY + (j+0.5)*size + Pars::mp_shift;   // Additional of unbalance: Pars::mp_shift
                par.x[count] = _x;
                par.y[count] = _y;
                count++;
            }
        }

        std::cout << "|X| Done generating the second (2nd) part\n";
    }
    
    return;
}

/**
 *  @brief  Uniform single resolution distribution for 2D simulation.
 *         
 *  @param  _particle  Particle data container.
*/
void initialization::init_2d_single_res(Particle &par)
{   
    // Number of particle in x and y direction
    int _nx = std::ceil(Pars::lxdom / Pars::sigma);     // Number of particle in x direction
    int _ny = std::ceil(Pars::lydom / Pars::sigma);     // Number of particle in y direction
    int _nparticle = _nx * _ny;                         // Number of all particle

    // Adjust the size of vector into the calculated particle number
    par.x.resize(_nparticle, 0.0e0);
    par.y.resize(_nparticle, 0.0e0);

    // Simulation domain pivot coordinate
    double pivX = -Pars::xdom;              // Location of the most left (x)
    double pivY = -Pars::sigma*_ny*0.5;     // Location of the most bottom (y)

    // Generate the particle position
    double _y, _x;
    int count = 0;
    for (int i = 0; i < _nx; i++){
        for (int j = 0; j < _ny; j++){
            _x = pivX + (i + 0.5)*Pars::sigma;
            _y = pivY + (j + 0.5)*Pars::sigma + Pars::mp_shift;
                                        // Additional of unbalance for early separation (mp_shift)
            // Assign the particle
            par.x[count] = _x;
            par.y[count] = _y;
            count ++;
        }
    }
    
    // Assign other particle properties
    par.num = _nparticle;
    par.s.resize(_nparticle, Pars::sigma);
    par.level.resize(_nparticle, Pars::max_level);
}

/**
 *  @brief  Two level resolution with finer resolution particle in single block 
 *  distribution for 2D simulation.
 *         
 *  @param  _particle  Particle data container.
 *  @param  _bodyList  The list of body data container used for particle generation.
*/
void initialization::init_2d_multi_res_single_block(Particle &par, const std::vector<Body> &bL)
{
    // internal temporary variable
    double _y, _x;

    // ====== PARTICLE GENERATION ====== //
    // number of particle in x and y direction
    int _nx = std::ceil(0.5 * Pars::lxdom / Pars::sigma);   // Number of particle in x direction
    int _ny = std::ceil(0.5 * Pars::lydom / Pars::sigma);   // Number of particle in y direction
    int _nparticle = _nx * _ny;                             // Number of all particle

    // Adjust the size of vector into the calculated particle number
    par.x.resize(_nparticle, 0.0e0);
    par.y.resize(_nparticle, 0.0e0);
    par.s.resize(_nparticle, 2 * Pars::sigma);
    par.level.resize(_nparticle, Pars::max_level - 1);
    par.num = _nparticle;

    // Simulation domain pivot coordinate
    double pivX = -Pars::xdom;              // Location of the most left (x)
    double pivY = -Pars::sigma*_ny;         // Location of the most bottom (y)

    // Generate the particle position
    int count = 0;
    for (int i = 0; i < _nx; i++){
        for (int j = 0; j < _ny; j++){
            _x = pivX + (i+0.5)*2*Pars::sigma;
            _y = pivY + (j+0.5)*2*Pars::sigma + Pars::mp_shift;   // Additional of unbalance: Pars::mp_shift
            par.x[count] = _x;
            par.y[count] = _y;
            count ++;
        }
    }
    
    // ====== DIVIDE AND CONQUEROR PROCEDURE ====== //
    // Put extension from the body limit to the finer resolution block limit (USER INPUT MANUALY)
    double ext_front = 0.5e0;  // Upstream extension
    double ext_back  = 0.5e0;  // Downstream extension
    double ext_side  = 0.5e0;  // Side body extension

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
                par.y[i] < bL[part].max_pos[1] + ext_side ){    // Check toward body upper side
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
        
        // Point 1  -> Replace the parent data
        par.x[i] = _x - 0.5 * Pars::sigma;
        par.y[i] = _y - 0.5 * Pars::sigma;
        par.s[i] = Pars::sigma;
        par.level[i] = Pars::max_level;

        // Point 2-4  -> Add new particle data
        for (int j = 1; j < 4; j++){
            basis_loop(d) chdTrans[d] = (2*((j>>d)%2)) - 1;
            par.x.push_back(_x + 0.5 * Pars::sigma * chdTrans[0]);
            par.y.push_back(_y + 0.5 * Pars::sigma * chdTrans[1]);
            par.s.push_back(Pars::sigma);
            par.level.push_back(Pars::max_level);
            par.num++;
        }
    }
    
    return;
}

/**
 *  @brief  Multiresolution based on multiblock distribution for 2D simulation.
 *         
 *  @param  _particle  Particle data container.
 *  @param  _bodyList  The list of body data container used for particle generation.
*/
void initialization::init_2d_multi_res_multi_block(Particle &par, const std::vector<Body> &_bodyList)
{
    // TO BE CONTINUE
}
/*
void initialization::init_2d_multi_res_multi_block(Particle &par, const std::vector<Body> &_bodyList){
    // Calculate L_max, level_max, level_min
    double L_max;       // The factor size from the object size
    double max_grid_size;       // The factor size from the object size
    
    L_max = Pars::grid_factor * Pars::lx;
    this->level_max = 1 + std::round(std::log10(L_max / Pars::sigma)/std::log10(2));
    this->level_min = this->level_max + 1 - Pars::levels;

    max_grid_size = Pars::sigma * std::pow(2, this->level_max - 1);

    // Generate the base grid, set the level to min_level
    // Find the boundary position of the particle
    xmin_par[0] = *std::min_element(p_transformed.x.begin(), p_transformed.x.end()) + max_grid_size/2;
    xmin_par[1] = *std::min_element(p_transformed.y.begin(), p_transformed.y.end()) + max_grid_size/2;
    xmax_par[0] = *std::max_element(p_transformed.x.begin(), p_transformed.x.end());
    xmax_par[1] = *std::max_element(p_transformed.y.begin(), p_transformed.y.end());


    // Calculte the number of grid base particle in x and y direction
    for (size_t i = 0; i < 2; i++)
    {
        npar[i] = std::floor(std::abs(xmax_par[i] - xmin_par[i]) / max_grid_size);
    }

    std::vector<double> _xgrid;
    std::vector<double> _ygrid;
    
    // Create the grid base particle position
    for (size_t i = 0; i < npar[0]; i++)    // x direction
    {
        double _xBoxPos = xmin_par[0] + static_cast<double>(i) * max_grid_size;
        _xgrid.push_back(_xBoxPos);
    }
    
    for (size_t i = 0; i < npar[1]; i++)    // y direction
    {
        double _yBoxPos = xmin_par[1] + static_cast<double>(i) * max_grid_size;
        _ygrid.push_back(_yBoxPos);
    }

    for (size_t i = 0; i < npar[0]; i++)
    {
        for (size_t j = 0; j < npar[1]; j++)
        {
            this->GridBase.x.push_back(_xgrid[i]);
            this->GridBase.y.push_back(_ygrid[j]);
            this->GridBase.s.push_back(max_grid_size);
            this->GridBase.gz.push_back(0);
            //this->GridBase.vorticity.push_back(0);
            this->GridBase.level.push_back(this->level_min);
            
        }
    }
    
    this->GridBase.num = this->GridBase.x.size();   // Update the GridBase size
    _num = this->GridBase.num;
    
    // Find the grid neighbor
    std::vector<std::vector<int>> GridNeighbor = d_base_remeshing.inter_search_base_grid(this->GridBase);         // neighbor of base grid
    this->GridBase.neighbor = GridNeighbor;
    
    // ========================================
    // Calculate the body panel
    // Calculate the body boundary

    std::vector<double> _minBodyCoordinate(2,0); // center coordinate of obstacle
    std::vector<double> _maxBodyCoordinate(2,0); // center coordinate of obstacle

    std::vector<double> _centerPanelX;
    std::vector<double> _centerPanelY;
    std::vector<double> _normalPanelX;
    std::vector<double> _normalPanelY;

    double _nPanel = b.x.size() - 1;
    _centerPanelX.resize(_nPanel);
    _centerPanelY.resize(_nPanel);
    _normalPanelX.resize(_nPanel);
    _normalPanelY.resize(_nPanel);

    // Transform into circular
    Particle b_transformed; // Base grid particle (final interpolated coordinate)
    b_transformed.x.resize(b.x.size(), 0.0);
    b_transformed.y.resize(b.y.size(), 0.0);
    for (int i = 0; i < b.x.size(); i++)
    { 
        b_transformed.x[i] = b.x[i]/Pars::Ex;
        b_transformed.y[i] = b.y[i]/Pars::Ey;
    }

    _minBodyCoordinate[0] = *std::min_element(b_transformed.x.begin(), b_transformed.x.end()) - max_grid_size * std::sqrt(2);
    _minBodyCoordinate[1] = *std::min_element(b_transformed.y.begin(), b_transformed.y.end()) - max_grid_size * std::sqrt(2);
    _maxBodyCoordinate[0] = *std::max_element(b_transformed.x.begin(), b_transformed.x.end()) + max_grid_size * std::sqrt(2);
    _maxBodyCoordinate[1] = *std::max_element(b_transformed.y.begin(), b_transformed.y.end()) + max_grid_size * std::sqrt(2);
    
    // Determine the center and normal panel of the body
    for (int i = 0; i < _nPanel; i++)
    {
        // determine midpoint of panel
        _centerPanelX[i] = b_transformed.x[i] + ((b_transformed.x[i + 1] - b_transformed.x[i]) * 0.5e0);
        _centerPanelY[i] = b_transformed.y[i] + ((b_transformed.y[i + 1] - b_transformed.y[i]) * 0.5e0);
        // determine normal panel
        double _dx = b_transformed.x[i + 1] - b_transformed.x[i];
        double _dy = b_transformed.y[i + 1] - b_transformed.y[i];
        double _theta = atan2(_dy, _dx);
        _normalPanelX[i] = sin(_theta);  // unit normal vector in x-direction
        _normalPanelY[i] = -cos(_theta); // unit normal vector in y-direction
    }
    
    // Calculate the disctance from body to update the levels
    for (int i = 0; i < _num; i++){
        // Consider the only grid near the geometry
        double _xp = this->GridBase.x[i];
        double _yp = this->GridBase.y[i];

        if (_xp >= _minBodyCoordinate[0] && _xp <= _maxBodyCoordinate[0] &&
            _yp >= _minBodyCoordinate[1] && _yp <= _maxBodyCoordinate[1])
        {
            // initial guess value
            int _k1 = 0;
            double _rmins = 100.0e0;
            double _signNormal = 0.0e0;
            double _testPin = 0.0e0;

            for (int k = 0; k < _nPanel; k++)
            {
                // Distance between center panel and the point
                double _r2panel = std::sqrt(std::pow((_xp - _centerPanelX[k]), 2) +
                                            std::pow((_yp - _centerPanelY[k]), 2));
                if (std::abs(_r2panel) <= _rmins)
                {
                    _rmins = _r2panel; // find the center panel closest to the point
                    _k1 = k;           // Name of closest panel
                }
            }

            // normal distance from considered point to panel
            // * = 0 when grid node lie on the panel, or extended-panel
            // Negative sign means outside, positive sign mean inside
            double _rNormal = (((_xp - _centerPanelX[_k1]) * (_normalPanelX[_k1])) +
                            ((_yp - _centerPanelY[_k1]) * (_normalPanelY[_k1])));

            if (_rNormal > - max_grid_size * std::sqrt(2) && _rNormal < max_grid_size * std::sqrt(2))
                this->GridBase.level[i] = this->level_max;
        }
    }

    // Perform neighbor update
    this->level_update_neighbor();

    // Perform multiresolution adaption
    this->generate_multiblock(p);

}
*/

/**
 *  @brief  Multiresolution body adjusted distribution for 2D simulation.
 *         
 *  @param  _particle  Particle data container.
 *  @param  _bodyList  The list of body data container used for particle generation.
*/
void initialization::init_2d_multi_res_body_adjusted(Particle &par, const std::vector<Body> &bL)
{    
    // internal temporary variable
    double _y, _x;
    const double baseSize = Pars::intPow(2,Pars::max_level) * Pars::sigma;
    geometry geom_tool;

    // ====== PARTICLE GENERATION ====== //
    // number of particle in x and y direction
    int _nx = std::ceil(Pars::lxdom / baseSize);    // Number of particle in x direction
    int _ny = std::ceil(Pars::lydom / baseSize);    // Number of particle in y direction
    int _nparticle = _nx * _ny;                     // Number of all particle

    // Adjust the size of vector into the calculated particle number
    par.x.resize(_nparticle, 0.0e0);
    par.y.resize(_nparticle, 0.0e0);
    par.s.resize(_nparticle, baseSize);
    par.level.resize(_nparticle, 0);
    par.num = _nparticle;

    // Simulation domain pivot coordinate
    double pivX = -Pars::xdom;          // Location of the most left (x)
    double pivY = -baseSize*_ny*0.5;    // Location of the most bottom (y)

    // Generate the particle position
    int count = 0;
    for (int i = 0; i < _nx; i++){
        for (int j = 0; j < _ny; j++){
            _x = pivX + (i + 0.5)*baseSize;
            _y = pivY + (j + 0.5)*baseSize + Pars::mp_shift; // Additional of unbalance 
            par.x[count] = _x;
            par.y[count] = _y;
            count ++;
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
                    if (par.x[ID] > bL[part].min_pos[0] - R_marg*MARGIN_EXTENSION_FACTOR &&   // Check toward body upstream
                        par.x[ID] < bL[part].max_pos[0] + R_marg*MARGIN_EXTENSION_FACTOR &&    // Check toward body downstream
                        par.y[ID] > bL[part].min_pos[1] - R_marg*MARGIN_EXTENSION_FACTOR &&    // Check toward body bottom side
                        par.y[ID] < bL[part].max_pos[1] + R_marg*MARGIN_EXTENSION_FACTOR ){    // Check toward body upper side
                        // Particle near this body, dont skip the DnC on this particle
                        skipDnC = false;
                        break;
                    }
                }
                if (skipDnC) continue;
            }

            // Calculate the particle minimum distance toward body
            pos = {par.x[ID], par.y[ID]};
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
                
                // Point 1
                par.x[ID] = _x - 0.25 * currSize;
                par.y[ID] = _y - 0.25 * currSize;
                par.s[ID] = par.s[ID]/2;
                par.level[ID] = par.level[ID] + 1;
                par_list_next.push_back(ID);

                // Point 2 - 4
                for (int j = 1; j < 4; j++){
                    basis_loop(d) chdTrans[d] = (2*((j>>d)%2)) - 1;
                    par_list_next.push_back(par.s.size());
                    par.x.push_back(_x + 0.25 * currSize * chdTrans[0]);
                    par.y.push_back(_y + 0.25 * currSize * chdTrans[1]);
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

/*
Note for 2D initialization
*/