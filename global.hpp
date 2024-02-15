#ifndef INCLUDED_GLOBAL
#define INCLUDED_GLOBAL

// <!> See a short hand code for time handling at the very bottom page

// #pragma region INCLUDE_LIBRARY
    // Input output data stream
    #include <iostream>     // I/O control on console
    #include <fstream>      // R/W from computer data

    // Data storage
    #include <vector>       // Vector container
    #include <string>       // String variable type

    // Utilities
    #include <time.h>       // Computation timer (faster comp.)
    #include <chrono>       // Computation timer (accurate comp.)
    #include <algorithm>    // Algorithm lib. (sort, max, min, ...)
    #include <omp.h>        // Parallel computing lib.
    #include <cmath>        // Math operation lib.

    // // Additionals
    // #include <math.h>
    // #include <complex>
    // #include <numeric>
// #pragma endregion

// #pragma region CONSOLE_COLOR_SET
    // Definition of font color in console output
    #define FONT_RESET  "\033[0m"       /* Reset */
    #define FONT_BLACK  "\033[30m"      /* Black */
    #define FONT_RED    "\033[31m"      /* Red */
    #define FONT_GREEN  "\033[32m"      /* Green */
    #define FONT_ORANGE "\033[33m"      /* Orange */
    #define FONT_BLUE   "\033[34m"      /* Blue */
    #define FONT_PURPLE "\033[35m"      /* Purple */
    #define FONT_CYAN   "\033[36m"      /* Cyan */
    #define FONT_GRAY   "\033[90m"      /* Gray */
    #define FONT_MAROON "\033[91m"      /* Maroon */
    #define FONT_PINK   "\033[95m"      /* Pink */
    #define FONT_WHITE  "\033[97m"      /* White */
    #define FONT_LIGHT_GRAY "\033[37m"  /* Light Gray */
    #define FONT_DARK_BLUE  "\033[94m"  /* Dark Blue */
    #define FONT_TORQUOISE  "\033[96m"  /* Torquoise */

    /**
     *  @brief  A shorthand for showing error message. (see global.hpp)
    */
    #define ERROR_LOG std::cout << FONT_MAROON << "[ERROR] " << FONT_RESET
    /**
     *  @brief  A shorthand for showing warning message. (see global.hpp)
    */
    #define WARNING_LOG std::cout << FONT_PURPLE << "[WARNING] " << FONT_RESET
    /**
     *  @brief  A shorthand for showing log message. (see global.hpp)
    */
    #define MESSAGE_LOG std::cout << FONT_CYAN << "[LOG] " << FONT_RESET
// #pragma endregion

namespace Pars{
    /**
     *  @brief  A mathematical power function operator for integer.
     *  NOTE : This function only work for non-negative integer exponent only!
     *         
     *  @param  base The base value of exponential.
     *  @param  exp  The exponent value of power.
     * 
     *  @return The power of 'base' by 'exp'.
    */
    constexpr long long int intPow(int base, int exp){
        long long int res = 1;
        for (int i = 0; i < exp; i++){
            res *= base;
        }
        return res;
    }
}

// #pragma region BASIC_CONSTANT_PARAMETER
    /** 
     * The domain spatial dimension:
     *  2:= 2D simulation; 
     *  3:= 3D simulation
    */
    #define DIM 2

    /** 
     * @brief  The number of obstacle body object
    */
    #define N_BODY 0

    /** 
     * A parallel timer flag:
     *  0:= Timer using chrono
     *  1:= Timer using parallel package
    */
    #define TIMER_PAR 1

    // #define PRECISION_FACTOR 1e-5   // A number factor to avoid wrong <double> to <int> truncation (e.g. calculating grid count from domain length division toward grid length)
    
    /** 
     * The flag to perform property interpolation (post processing):
     *  0:= Basic default calculation
     *  1:= Calculate the data interpolation [POST-PROCESS]
     * NOTE: Only work for data based on grid node calculation
    */
    #define DATA_INTERPOLATION 0
// #pragma endregion


/**
 *  @brief  A shorthand of a for loop through the dimension basis. (see global.hpp)
*/
#define basis_loop(i) for(int i = 0; i < DIM; i++)

/* GLOBAL PARAMETER MAIN GUIDE!!
    # Simulation flag                   : Toogle for specific procedure in running the program
    # Simulation Parameter Option       : The option for a particular object in the program
    # Simulation Domain Parameter       : DOMAIN data properties
    # Incompressible Fluid Properties   : FLUID data properties
    # Computation Parameter             : Data parameter for computation control (basic simulation)
    # Neighboring Parameter             : NEIGHBOR basic coefficient
    # Solid Body Parameter              : BODY Object basic data
    # Penalization Parameter            : BRINKMAN Penalization data coefficient
    # FMM Parameter Setting             : FMM method data parameter
    # Input/Output Data Parameter       : Parameter for file and console data stream
    # Structural Vibration Parameter    : [OLD] Vibrational data
*/

/**
 *  @brief  A list of simulation parameter and setting. (see global.hpp)
*/
namespace Pars
{
    // Predescription for each parameter variable
    // [A] -> Already given in this program
    // [L] -> Given a value but can be altered by user
    // [I] -> Input for PARAMETER form user
    // [O] -> Input for OPTION form user, if not will use default
    // [C] -> Calculated form user input
    // [P] -> Parameter in program calculation

    // =========================================
    // +----------- Simulation Flag -----------+
    // =========================================
    // #pragma region SIMULATION_FLAG
        // Data Write Flag
        extern const bool flag_save_body;       // [I] The flag to save body data
        extern const bool flag_save_cell;       // [I] The flag to save cell list data (*ADAPTIVE CELL LIST particle container)
        extern const bool flag_save_node;       // [I] The flag to save node list data (*GRID NODE particle container)
        extern const bool flag_save_ngh_par;    // [I] The flag to save neighbor particle list
        extern const bool flag_save_ngh_num;    // [I] The flag to save neighbor number
        extern const bool flag_save_parameter;  // [I] The flag to save simulation parameter
        extern const bool flag_save_sim_time;   // [I] The flag to save simulation time and computational time
        extern const bool flag_save_residual;   // [I] The flag for calculating and saving residual

        // Console Display Flag
        extern const bool flag_disp_stability;  // [I] The flag to display stability parameter at each iteration
        extern const bool flag_disp_pred_time;  // [I] The flag to display the prediction of total computational time until the end of simulation
        
        // Operator Flag
        extern const bool flag_pressure_calc;   // [I] The flag to calculate pressure on post processing
        extern const bool flag_adaptive_dist;   // [I] The flag for adaptive particle distribution
        extern const bool flag_fast_remesh;     // [I] The flag for effective remeshing calculation

        // Additional Flag
        extern const bool flag_ngh_include_self;        // [I] The flag for neighbor evaluation include itself
        extern const bool flag_slightly_shifted_domain; // [I] The flag to generate slightly unsymetrical domain (To make early separation)
    // #pragma endregion SIMULATION_FLAG

    // =====================================================
    // +----------- Simulation Parameter Option -----------+
    // =====================================================
    // #pragma region PARAMETER_OPTION
        /* [O] Option for starting state of simulation:
            0:= Start over the simulation from 0.0s;
            1:= Resume simulation from iteration "resume_step"
        */
        extern const int opt_start_state;

        /* [O] The particle initialization type option:
            0:= Testing Resolution;
            1:= Single Resolution;
            2:= Multi-Res Single Block;
            3:= Multi-Res Multi Block; ** [Not Ready YET]
            4:= Multi-Res Body Adjusted;
            5:= Grid Node Based
        */
        extern const int opt_init_particle;

        /* [O] The vorticity initialization type option:
            0:= No initialized vorticity;
            1:= Perlman vorticity;
            2:= Reserved ...;
            3:= Reserved ...;
        */
        extern const int opt_init_vorticity;

        /* [O] The body type option:
                    2D simulation   |   3D simulation   
                --------------------|-------------------
            1:= Circular cylinder   |  *Sphere;
            2:= Square cylinder     |  *Cube;
            3:= Normal plate        |  *3D normal plate;
            4:= Flat plate          |  *3D flat plate;
            5:= NACA airfoil        |  *Torus;
            6:= -                   |  *Heart;
        */
        extern const std::vector<int> opt_body;

        /* [O] The neighbor evaluation type option:
            0:= Direct Neighbor Search;
            1:= Linked List;
            2:= Cell List;
            3:= Spatial Hash;
            4:= Grid Node base
        */
        extern const int opt_neighbor;

        /* [O] The brinkman iterative type option:
            1:= Classical brinkmann;
            2:= Iterative brinkmann
        */
        extern const int opt_pen_iter;

        /* [O] The force calculation type;\
            0:= No force calculation;\
            1:= By direct method;\
            2:= By penalization method;\
            3:= By impulse method\
            4:= By NOCA method [Not ready yet]
        */
        extern const int opt_force_type;

        /* [O] Body motion in x-direction option:
            1:= Sudden move and stop at 1.0s;
            2:= Smoothly start and stop;
            3:= Constant velocity motion
        */
        extern const int opt_body_motion;

        /* [O] The option to set the interval type of data writting;
            1:= Write each simulation iteration;
            2:= Prescribed total file number;
            3:= Write each computational time
        */
        extern const int opt_data_write;
        
        /* [O] The neighbor interaction option:
            1:= Support domain only using own size;
            2:= Average size neighbor evaluation
        */
        extern const int opt_ngh_interact;

    // #pragma endregion PARAMETER_OPTION

    // =====================================================
    // +----------- Simulation Domain Parameter -----------+
    // =====================================================
    // #pragma region DOMAIN_PARAMETER
        extern const double lxdom;      // [I] Initial domain length on x-axis
        extern const double lydom;      // [I] Initial domain length on y-axis
        extern const double lzdom;      // [I] Initial domain length on z-axis
        extern const double xdom;       // [I] Negative x-direction domain length
        extern const double xcenter;    // [I] Initial domain center (x-axis), default 0.0
        extern const double ycenter;    // [I] Initial domain center (y-axis), default 0.0
        extern const double zcenter;    // [I] Initial domain center (y-axis), default 0.0

        /*  Domain Illustration, see below.
                         __________________________________
                        |                                  |
                        |                                  |
                        |                                  |
                        |     (xcenter, ycenter)           |
                        |<------------>*                   | lydom
                        |     xdom                         |
                        |                                  |
                        |                                  |
                        |__________________________________|
                                        lxdom
        */
    // #pragma endregion
    
    // =========================================================
    // +----------- Incompressible Fluid Properties -----------+
    // =========================================================
    // #pragma region FLUID_PROPERTIES
        extern const double RE;     // [I] The reynold number
        extern const double u_inf;  // [I] Freestream x-direction velocity
        extern const double v_inf;  // [I] Freestream y-direction velocity
        extern const double w_inf;  // [I] Freestream z-direction velocity
        extern const double RHO;    // [I] Fluid density
        extern const double U_inf;  // [C] Freestream velocity magnitude
        extern const double NU;     // [C] Kinematic viscousity
        extern const double MU;     // [C] Dynamic viscousity
    // #pragma endregion

    // ===============================================
    // +----------- Computation Parameter -----------+
    // ===============================================
    // #pragma region COMPUTATIONAL_PARAMETER
        // Basic parameter
        extern const double sigma;      // [I] Particle core size
        extern const double dt;         // [I] Simulation time step
        extern const double sim_time;   // [I] The total simulation time
        extern const int resume_step;   // [I] Iteration step ID of data for resuming simulation

        // Stability Criteria
        extern const double phi_s;      // [A] Euler explicit scheme
        extern const double Courant;    // [C] Courant number(C): C = U_inf * dt / sigma, where 0 < C <= 1
        extern const double Diffusion;  // [C] Diffusion number(Phi): Phi =  NU * dt / sigma^2, where 0 < Phi <= 0.5
        // extern const double Courant =           // [TYPE 2] sigma = sqrt(dt*vis/phi_s)/Courant, where 0 < Courant <= 1
        //     std::sqrt(dt * vis / phi_s) / sigma;
        
        // **Particle redistribution constant interval evaluation
        extern const int rmsh_inv;          // [I] Step interval for remeshing evaluation
        extern const int adapt_inv;         // [I] Step interval for adaptation evaluation
        extern const int ngh_diff_level;    // [I] Maximum neighbor node different level
        extern const double adapt_tol;      // [I] The adaptation tolerance
        extern const double active_sig;     // [I] The significant ratio toward the source maximum value [*vorticity] for an active particle
    // #pragma endregion

    // ===============================================
    // +----------- Neighboring Parameter -----------+
    // ===============================================
    // #pragma region NEIGHBOR_PARAMETER
        extern const double r_sup;      // [L] Neighbor evaluation support domain ratio size
        extern const double r_buff;     // [L] The buffer region radius factor
        extern const double body_ext;   // [I] The body extention distance to evaluate the chi and active particle
        extern const double mp_shift;   // [L] A distance shift from the midplane to create the unbalance calculation
        extern const int max_level;     // [I] Number of resolution level
    // #pragma endregion NEIGHBOR_PARAMETER

    // ==============================================
    // +----------- Solid Body Parameter -----------+
    // ==============================================
    // #pragma region BODY_PARAMETER
        // Basic parameter (body signature)
        extern const double Df;         // [I] Simulation reference length (chord or diameter)
        extern const int n_a;           // [I] Number of body surface node

        // Geometry parameter (each obstacle)
        extern const std::vector<double> lxbody;        // [I] Body x length
        extern const std::vector<double> lybody;        // [I] Body y length
        extern const std::vector<double> lzbody;        // [I] Body z length
        extern const std::vector<double> Lref;          // [I] Body reference length (chord or diameter)
        extern const std::vector<double> H_star;        // [I] Thickness parameter for plate geometry: H*=t/L
        extern const std::vector<double> x_body_cen;    // [I] Body x center position
        extern const std::vector<double> y_body_cen;    // [I] Body y center position
        extern const std::vector<double> z_body_cen;    // [I] Body z center position

        // Body reading data
        extern const std::vector<bool> readBody;        // [I] Flag to generate body by read file
        extern const std::vector<bool> bodyRotFlag;     // [I] Flag to rotate body transformation
        extern const std::vector<double> bodyRotDeg;    // [I] Angle of body rotation transformation (*in degree)
        extern const std::vector<int> bodyRotAxis;      // [I] Axis of body rotation transformation
        
        // Initial body velocity
        extern const double ubody;      // [I] Body x-direction velocity
        extern const double vbody;      // [I] Body y-direction velocity
        extern const double wbody;      // [I] Body z-direction velocity
        extern const double omega;      // [I] Body z-direction angular velocity

        // NACA Series parameter (2D Simulation only)
        // *Parameter of NACA 4-digit series: NACA[m*100][p*10][t*100] (e.g. NACA-4212: m=0.04, p=0.2, t=0.12)
        extern const double m_a;        // [I] Maximum camber value (limited to 0.04 for a good airfoil contour, problem at trailing edge)
        extern const double p_a;        // [I] Maximum camber location (in tenth of chord)
        extern const double t_a;        // [I] Maximum thickness of airfoil (in percent of chord)
        extern const double acx;        // [I] Aerodynamic center position in x basis (relative to leading edge)
        extern const double acy;        // [I] Aerodynamic center position in y basis (relative to leading edge)
        extern const double AoA;        // [I] Attack angle (in degree)
    // #pragma endregion BODY_PARAMETER

    // ================================================
    // +----------- Penalization Parameter -----------+
    // ================================================
    // #pragma region PENALIZATION_PARAMETER
        extern const int opt_pen;       /*/ [O] The penalization type option (*method from Rasmussen 2011)
                                            1:= Implicit penalization [lambda = 1.0e4];
                                            2:= Semi-implicit penalization;
                                            3:= Explicit penalization [lambda = 1.0/dt]
                                        /*/
        extern int opt_kaipen;          /*/ [O] The option for penalization masking χ (chi)
                                            0:= No chi recalculation;
                                            1:= Recalculate by kai, to avoid singularities (*Rasmussen phD 2011)
                                        /*/ 
        extern double lambda;           // [I] Value of penalization constant
        extern const int numpen;        // [L] Number of particle from body surface for penalization domain evaluation (*6 from some reference)
        extern const double hmollif;    // [I] The half mollification length [sqrt(2)*sigma]
        extern const int pen_iter;      // [I] The penalization iteration for iterative brinkman type
    // #pragma endregion PENALIZATION_PARAMETER

    // ===============================================
    // +----------- FMM Parameter Setting -----------+
    // ===============================================
    // <!> Need further check
    // #pragma region FMM_SETTING
        // Tree option
        extern const double expTree;    // [I] The percentage of tree root cell size expansion
        extern const int tree_lvl_max;  // [I] The tree level of the finest cell [OLD PACKAGE]
        extern const int tree_lvl_min;  // [I] The tree level of the coarse cell [OLD PACKAGE]

        // FMM option
        extern const int ioptfmm;       /*/ [O] Selection of the velocity calculation method;
                                            0:= Direct calculation,
                                            1:= FMM accelerated
                                        /*/ 
        extern const int P_max;         // [O] Number of the expansion (P_max=10 for error in order of 10^-3)
        extern const int par_count_max; // [I] Maximum number of all particle inside the cell
        extern const int src_count_max; // [I] Maximum number of source particle inside the cell        

        // Old FMM parameter
        extern const int lmax;          // [I] Max of FMM levels
        extern const int npcm;          // [I] Number of pole expansion
        extern const int nfwrd;         // std::pow(4,(lmax+1));   //nfwrd * (std::pow(4, (lmax + 1)) - 4) / 3;
        extern const int nbmax;         // Geometric Series formula
        extern const int nbmrl;         //nfwrd * 4096; //!!need to change if so that nbmrl >= 4^lmax
        extern const int nbl1;          // [A] Number of list 1 neighbor
        extern const int nbl2;          // [A] Number of list 2 neighbor
        extern const int nbl3;          // [A] Number of list 3 neighbor
        extern const int nbl4;          // [A] Number of list 4 neighbor
        extern const double tol2;       // [A] The maximum tolerance of cell size (In evaluating adjacent cell)
        extern const int n_s;           // [I] Number of maximum total particle inside a cell
        extern const int ndp;           // [I] The order of multipole expansion
        extern const int icutoff;       /*/ [O] The redistribution function for direct potential calculation
                                            0:= Singular,
                                            1:= Super (high-oder) algebraic,
                                            2:= Gaussian,
                                            3:= Super Gaussian
                                        // Note:
                                            -> if PSE calculation : take 2 or 3;
                                            -> if Direct caluculation : take 1;
                                            -> if FMM calculation : take 2 or 3
                                        /*/ 
    // #pragma endregion FMM_SETTING

    // ======================================================
    // +----------- Setting of Input/Output Data -----------+
    // ======================================================
    // #pragma region DATA_WRITE_SETTING
        extern const int max_iter;          // [C] Total iterations step (calculated from simulation time and dt)
        
        // Parameter of data write
        extern const int step_inv;          // [I] Step interval for saving data parameter [Type 1]
        extern const int file_num;          // [I] Total file to be saved parameter [Type 2]
        extern const double comp_time_inv;  // [I] Computational time (in second) interval for data writting parameter [Type 3]

        extern const int save_inv;          // [C] Step interval for saving data used in the program [Type 1&2]
        extern double cum_comp_time;        // [P] Cumulative computational time for data saving [Type 3]

        // Parameter in writting data
        extern const int max_dig_len;       // [C] Digit length of iteration name ID at final iteration
        extern const int par_dom_ext;       // [I] Domain extension for particle data saving (see, save_data/save_data.cpp)

        // Other
        // extern const int nt_sf;             // [I] save data per ( nt_sf ) time step, for saving storage MEMORY and for case of running is stopped suddenly
        extern const int opt_extra_data;    // 1:= not print extra data, print only common data
                                                // 2:= print common data and extra data at the same time,
                                                // 3:= print extra data from common data, MUST run case 1/2 first
        extern const double comtime_sf;           // % Save file frequently after how many [second]=[hour]*60*60, only in case of running is stopped suddenly, but nt_sf take long days
        extern const int nt_sf;             // [I] save data per ( nt_sf ) time step, for saving storage MEMORY and for case of running is stopped suddenly
    // #pragma endregion DATA_WRITE_SETTING

    // ========================================================
    // +----------- Structural Vibration Parameter -----------+
    // ========================================================
    // #pragma region VIBRATION_PARAMETER  // Contributed by Ical
        // VIBRATION / STRUCTURAL PROPERTIES : !!STILL IN PROGRESS, RESULTS NOT GOOD ENOUGH!!
        extern const int vib;           // =0 no vibration simulation, =1 circular cylinder simulation, =2 flow induced pitch oscillations
        extern double SpringConst;      // Spring constant :: !!INITIALIZATION PARAMETERS DO NOT CHANGE!!
        extern double DamperConst;      // Damper constant :: !!INITIALIZATION PARAMETERS DO NOT CHANGE!!
        extern double mass;             //Object mass :: !!INITIALIZATION PARAMETERS DO NOT CHANGE!!
        extern double inertia;          // Inertia :: !!STILL UNUSED!!
        //extern double m_d;            // !! rho times volume enclosed !! 
        extern double m_d;              // !! rho times volume enclosed !! UNUSED, SET TO 0.0 
        extern double tetha;            //!! ROTATION PARAMETER, STILL BUGGED
        extern double tetha_dot;        //!! ROTATION PARAMETER, STILL BUGGED
        extern const double m_star;     // NONDIMENSIONAL MASS
        extern const double k_star;     // NONDIMENSIONAL SPRING CONST
        extern const double c_star;     //NONDIMENSIONAL DAMPER CONST
        extern const double t_star;     //NONDIMENSIONAL TIME VARIABLE
        extern double gaya;             // !! FORCE INITIALIZATION !!
        extern double momen;            // !! MOMENT INITIALIZATION !!
        extern const double chi_star;   //!! ROTATION PARAMETER, STILL BUGGED
        extern const double i_star;     //!! ROTATION PARAMETER, STILL BUGGED
        extern const double U_star;     //!! ROTATION PARAMETER, STILL BUGGED
        extern const double tetha_nol;  //equilibrium angular position of spring [deg] //!! ROTATION PARAMETER, STILL BUGGED
    // #pragma endregion VIBRATION_PARAMETER

    // =====================================================
    // +-------------- Interpolation Process --------------+
    // =====================================================
    // #pragma region INTERPOLATION_PARAMETER
        extern const double lxdomInt;   // [I] Initial domain length on x-axis
        extern const double lydomInt;   // [I] Initial domain length on y-axis
        extern const double lzdomInt;   // [I] Initial domain length on z-axis
        extern const double xdomInt;    // [I] Negative x-direction domain length
        extern const double sigmaInt;   // [I] Initial domain center (x-axis), default 0.0
    // #pragma endregion INTERPOLATION_PARAMETER

} // namespace Pars

#endif

/* A short hand patter of ticking the computational time (chrono)
    // Clock manager -> Use this pattern to calculate the computational time
    std::chrono::steady_clock timer;
    auto start = timer.now();
    auto end = timer.now();
    auto span = static_cast<std::chrono::duration<double>>(end - start);
    std::cout << "Operation took: " << span.count() << " s\n";
*/

/* A short hand patter of ticking the computational time (clock)
    // Clock manager -> Use this pattern to calculate the computational time
    clock_t start = clock();
    clock_t end = clock();
    clock_t span = end - start;
    std::cout << "Operation took: " << (double)span/CLOCKS_PER_SEC << " s\n";
*/

/* A short hand patter of ticking the computational time (omp time)
    // Clock manager -> Use this pattern to calculate the computational time
    double start = omp_get_wtime();
    double end = omp_get_wtime();
    double span = end - start;
    std::cout << "Operation took: " << span << " s\n";
*/
