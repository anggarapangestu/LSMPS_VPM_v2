#include "global.hpp"

#define _O_ true
#define _X_ false

#if DIM == 2    // The parameter for 2 Dimension Simulation
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
        extern const bool flag_save_body = _O_;         // [I] The flag to save body data
        extern const bool flag_save_cell = _X_;         // [I] The flag to save cell list data (*ADAPTIVE CELL LIST particle container)
        extern const bool flag_save_node = _O_;         // [I] The flag to save node list data (*GRID NODE particle container)
        extern const bool flag_save_ngh_par   = _X_;    // [I] The flag to save neighbor particle list
        extern const bool flag_save_ngh_num   = _X_;    // [I] The flag to save neighbor number
        extern const bool flag_save_parameter = _O_;    // [I] The flag to save simulation parameter
        extern const bool flag_save_sim_time  = _O_;    // [I] The flag to save simulation time and computational time
        extern const bool flag_save_residual  = _X_;    // [I] The flag for calculating and saving residual [STILL NOT CHECKED]

        // Console Display Flag
        extern const bool flag_disp_stability = _O_;    // [I] The flag to display stability parameter at each iteration
        extern const bool flag_disp_pred_time = _O_;    // [I] The flag to display the prediction of total computational time until the end of simulation
        
        // Operator Flag
        extern const bool flag_pressure_calc  = _X_;    // [I] The flag to calculate pressure on post processing [STILL NOT CHECKED]
        extern const bool flag_adaptive_dist  = _O_;    // [I] The flag for adaptive particle distribution
        extern const bool flag_fast_remesh    = _X_;    // [I] [TURN this OFF] The flag for effective remeshing calculation (Turns out not very effective that I thought, otherwise slightly slower)
        
        // Additional Flag
        extern const bool flag_ngh_include_self = _O_;        // [I] The flag for neighbor evaluation include itself [STILL NOT CHECKED]
        extern const bool flag_slightly_shifted_domain = _X_; // [I] The flag to generate slightly unsymetrical domain (To make early separation)
    // #pragma endregion

    // =====================================================
    // +----------- Simulation Parameter Option -----------+
    // =====================================================
    // #pragma region PARAMETER_OPTION
        extern const int opt_start_state = 0;
                /* [O] Option for starting state of simulation:
                    0:= Start over the simulation from 0.0s;
                    1:= Resume simulation from iteration "resume_step"
                */
        
        extern const int opt_init_particle = 5;
                /* [O] The particle initialization type option:
                    0:= Testing Resolution;
                    1:= Single Resolution;
                    2:= Multi-Res Single Block;
                    3:= Multi-Res Multi Block; ** [Not Ready YET]
                    4:= Multi-Res Body Adjusted;
                    5:= Grid Node Based
                */

        extern const std::vector<int> opt_body = {1,5,4,1,4};
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

        extern const int opt_neighbor = 4;
                /* [O] The neighbor evaluation type option:
                    0:= Direct Neighbor Search;
                    1:= Linked List;
                    2:= Cell List;
                    3:= Spatial Hash;
                    4:= Grid Node base
                */

        extern const int opt_pen_iter = 1;
                /* [O] The brinkman iterative type option:
                    1:= Classical brinkmann;
                    2:= Iterative brinkmann
                */

        extern const int opt_force_type = 1;
                /* [O] The force calculation type;\
                    0:= No force calculation;\
                    1:= By direct method;\
                    2:= By penalization method;\
                    3:= By impulse method\
                    4:= By NOCA method [Not ready yet]
                */

        extern const int opt_body_motion = 1;
                /* [O] Body motion in x-direction option:
                    1:= Sudden move and stop at 1.0s;
                    2:= Smoothly start and stop;
                    3:= Constant velocity motion
                */

        extern const int opt_data_write = 1;
                /* [O] The option to set the interval type of data writting;
                    1:= Write each simulation iteration;
                    2:= Prescribed total file number;
                    3:= Write each computational time
                */

        extern const int opt_ngh_interact = 1;
                /* [O] The neighbor interaction option:
                    1:= Support domain only using own size;
                    2:= Average size neighbor evaluation
                */
    // #pragma endregion

    // =====================================================
    // +----------- Simulation Domain Parameter -----------+
    // =====================================================
    // #pragma region DOMAIN_PARAMETER
        extern const double lxdom = 25.0e0;     // [I] Initial domain length on x-axis
        extern const double lydom = 10.0e0;     // [I] Initial domain length on y-axis
        extern const double lzdom = 0.0e0;      // [I] Initial domain length on z-axis
        extern const double xdom  = 4.0e0;      // [I] Negative x-direction domain length
        extern const double xcenter = 0.0e0;    // [I] Initial domain center (x-axis), default 0.0
        extern const double ycenter = 0.0e0;    // [I] Initial domain center (y-axis), default 0.0
        extern const double zcenter = 0.0e0;    // [I] Initial domain center (z-axis), default 0.0

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
        extern const double RE = 550.0e0;       // [I] The reynold number
        extern const double u_inf = 1.0e0;      // [I] Freestream x-direction velocity
        extern const double v_inf = 0.0e0;      // [I] Freestream y-direction velocity
        extern const double w_inf = 0.0e0;      // [I] Freestream z-direction velocity
        extern const double RHO = 1.0e0;        // [I] Fluid density
        extern const double U_inf =             // [C] Freestream velocity magnitude
            std::sqrt(u_inf*u_inf + v_inf*v_inf + w_inf*w_inf);
        extern const double NU = U_inf * Df / RE;   // [C] Kinematic viscousity
        extern const double MU = RHO * NU;          // [C] Dynamic viscousity
    // #pragma endregion

    // ===============================================
    // +----------- Computation Parameter -----------+
    // ===============================================
    // #pragma region COMPUTATIONAL_PARAMETER
        // Basic parameter
        extern const double sigma    = 0.01e0;      // [I] Particle core size; <?> default: 0.0025, dt~(1/2)*sigma - (1/3)*sigma
        extern const double dt       = 0.005e0;      // [I] Simulation time step; <?> default:0.001, dt <= phi_s*sigma^2/vis (Ploumhans [2000]) OR dt = phi_s*sigma^2*Courant^2/vis, where 0 < Courant <= 1
        extern const double sim_time = 200.0e0;     // [I] The total simulation time
        extern const int resume_step = 3300; //39600;//      // [I] Iteration step ID of data for resuming simulation

        // Stability Criteria
        extern const double phi_s     = 0.595e0;    // [A] (Ploumhans [2000]) for the Euler explicit scheme,phi_s = 0.595, !for the Adams–Bashforth 2 scheme, phi_s = 0.297, BUt can not use AB2 due to redistribution changes np--> need Biot-Savart again
        extern const double Courant   =             // [C] Courant number(C): C = U_inf * dt / sigma, where 0 < C <= 1
            U_inf * dt / sigma;
        extern const double Diffusion =             // [C] Diffusion number(Phi): Phi =  NU * dt / sigma^2, where 0 < Phi <= 0.5
            NU * dt / (sigma*sigma);
        // extern const double Courant =           // sigma = sqrt(dt*vis/phi_s)/Courant, where 0 < Courant <= 1
        //     std::sqrt(dt * vis / phi_s) / sigma;
        
        // **Particle redistribution constant interval evaluation
        extern const int rmsh_inv = 1;              // [I] Step interval for remeshing evaluation
        extern const int adapt_inv = 5;             // [I] Step interval for adaptation evaluation
        extern const int ngh_diff_level = 1;        // [I] Maximum neighbor node different level (for NDL criteria)
        extern const double adapt_tol = 5.0e-2;     // [I] The adaptation tolerance (factor to set the window difference of value between levels)
        extern const double active_sig = 1.0e-6;    // [I] The significant ratio toward the source maximum value [*vorticity] for an active particle

    // #pragma endregion

    // ===============================================
    // +----------- Neighboring Parameter -----------+
    // ===============================================
    // #pragma region NEIGHBOR_PARAMETER
        extern const double r_sup = 3.2;            // [L] Neighbor evaluation support domain ratio size [3.1 for 2D [Tanaka 2018]; 2.3 for 3D [Educated Guess]] 
        extern const double r_buff = 1.1;           // [L] The buffer region radius factor [1.2 is sufficient]
        // extern const double body_ext = 15.0*sigma;  // [I] The body extention distance to evaluate the chi and active particle (must bigger than Pars::numpen * Pars::sigma)
        extern const double body_ext = 0.1*Pars::Df;// [I] The body extention distance to evaluate the chi and active particle (must bigger than Pars::numpen * Pars::sigma)
        extern const double mp_shift = 1.0*sigma;   // [L] A distance shift from the midplane to create the unbalance calculation [1 core size]
        extern const int max_level = 4;             // [I] Number of resolution level
                                                        // -> Level is counted from 0 as the largest particle size
                                                        // -> Maximum level is the finest particle size (max_level = 0 -> single resolution)
                                                        // -> The number of resolution level is 'max_level' + 1
    // #pragma endregion

    // ==============================================
    // +----------- Solid Body Parameter -----------+
    // ==============================================
    // #pragma region BODY_PARAMETER
        // Basic parameter (body signature)
        extern const double Df = 1.0e0;     // [I] Simulation reference length (chord or diameter)
        extern const int n_a = 400;         // [I] Number of body surface node

        // Geometry parameter (each obstacle)
        // Basic size parameter
        extern const std::vector<double> lxbody =       // [I] Body z length
                    {1.0, 1.0, 1.0, 1.0};
        extern const std::vector<double> lybody =       // [I] Body z length
                    {1.0, 1.0, 1.0, 1.0};
        extern const std::vector<double> lzbody =       // [I] Body z length
                    {1.0, 1.0, 1.0, 1.0};
        
        // Reference size parameter
        extern const std::vector<double> Lref =         // [I] Body reference length (chord or diameter)
                    {1.0, 2.0, 1.5, 1.0};
        extern const std::vector<double> H_star =       // [I] Thickness parameter for plate geometry: H*=t/L
                    {0.5, 0.05, 0.05, 0.2};
        
        // Center coordinate
        extern const std::vector<double> x_body_cen =   // [I] Body x center position
                    {0.0, 3.0, 3.0, 10.0};
        extern const std::vector<double> y_body_cen =   // [I] Body y center position
                    {0.0, 1.0, -1.4, 3.0};
        extern const std::vector<double> z_body_cen =   // [I] Body z center position
                    {0.0, 0.0, 0.0, 0.0};
        
        // Body reading (limited to 3D)
        extern const std::vector<bool> readBody =       // [I] Flag to generate body by read file
                    {_X_, _X_, _X_, _X_};
        extern const std::vector<bool> bodyRotFlag =    // [I] Flag to rotate body transformation
                    {_X_, _X_, _X_, _X_};
        extern const std::vector<double> bodyRotDeg =   // [I] Angle of body rotation transformation (*in degree)
                    {0.0, 0.0, 0.0, 0.0};
        extern const std::vector<int> bodyRotAxis =     // [I] Axis of body rotation transformation
                    {0, 0, 0, 0};

        // Initial body velocity
        extern const double ubody = 0.0e0;      // [I] Body x-direction velocity
        extern const double vbody = 0.0e0;      // [I] Body y-direction velocity
        extern const double wbody = 0.0e0;      // [I] Body z-direction velocity
        extern const double omega = 0.0e0;      // [I] Body z-direction angular velocity

        // NACA symmetric Series parameter (2D Simulation only)
        // *Parameter of NACA 4-digit series: NACA[m*100][p*10][t*100] (e.g. NACA-4212: m=0.04, p=0.2, t=0.12)
        extern const double m_a = 0.02e0;       // [I] Maximum camber value (limited to 0.04 for a good airfoil contour, problem at trailing edge)
        extern const double p_a = 0.4e0;        // [I] Maximum camber location (in tenth of chord)
        extern const double t_a = 0.12e0;       // [I] Maximum thickness of airfoil (in percent of chord)
        extern const double acx = 0.25e0;       // [I] Aerodynamic center position in x basis (relative to leading edge)
        extern const double acy = 0.0e0;        // [I] Aerodynamic center position in y basis (relative to leading edge)
        extern const double AoA = 30.0e0;       // [I] Attack angle (in degree)
    // #pragma endregion

    // ================================================
    // +----------- Penalization Parameter -----------+
    // ================================================
    // #pragma region PENALIZATION_PARAMETER
        extern const int opt_pen = 1;       // [O] The penalization type option (*method from Rasmussen 2011)
                                            //    1:= Implicit penalization [lambda = 1.0e4];
                                            //    2:= Semi-implicit penalization;
                                            //    3:= Explicit penalization [lambda = 1.0/dt]
        int opt_kaipen = 0;                 // [O] The option for penalization masking χ (chi)
                                            //    0:= No chi recalculation;
                                            //    1:= Recalculate by kai, to avoid singularities (*Rasmussen phD 2011)
        double lambda = 1.0e4;              // [I] Value of penalization constant
        extern const int numpen = 6;        // [L] Number of particle from body surface for penalization domain evaluation (*6 from some reference)
        extern const double hmollif =       // [I] The half mollification length [sqrt(2)*sigma]
            1.0e0 * std::sqrt(2.0e0) * sigma;
        extern const int pen_iter = 10;     // [I] The penalization iteration for iterative brinkman type
    // #pragma endregion

    // ===============================================
    // +----------- FMM Parameter Setting -----------+
    // ===============================================
    // <!> Need further check
    // #pragma region FMM_SETTING
        // Tree option
        extern const double expTree = 1.0;      // [I] The percentage of tree root cell size expansion
        extern const int tree_lvl_max = 7;      // [I] The tree level of the finest cell [OLD PACKAGE]
        extern const int tree_lvl_min = 6;      // [I] The tree level of the coarse cell [OLD PACKAGE]

        // FMM option
        extern const int ioptfmm = 1;           // [O] Selection of the velocity calculation method;
                                                //    0:= Direct calculation,
                                                //    1:= FMM accelerated
        extern const int P_max = 12;            // [O] Number of the expansion order (P_max ~ log (e) / log(2) or error in order of 10^-3 for P_max = 10), This affect the computational time further for M2L calculation
        extern const int par_count_max = 50;    // [I] Maximum number of all particle inside the cell
        extern const int src_count_max = 30;    // [I] Maximum number of source particle inside the cell (The bigger will reduce M2L calculation)

        // Old FMM parameter (*ABOUT TO DELETE)
        extern const int lmax = 7;              // [I] Max of FMM levels
        extern const int npcm = 100;            // [I] Number of pole expansion
        extern const int nfwrd = 16;            // std::pow(4,(lmax+1));   //nfwrd * (std::pow(4, (lmax + 1)) - 4) / 3;
        extern const int nbmax = 4 * (Pars::intPow(4,(lmax))-1) / (4-1);   // Geometric Series formula
        extern const int nbmrl = nbmax + 5;     //nfwrd * 4096; //!!need to change if so that nbmrl >= 4^lmax
        extern const int nbl1 = 8 * lmax;       // [A] Number of list 1 neighbor
        extern const int nbl2 = 256;            // [A] Number of list 2 neighbor
        extern const int nbl3 = 256;            // [A] Number of list 3 neighbor
        extern const int nbl4 = 16;             // [A] Number of list 4 neighbor
        extern const double tol2 = 1.0e-15;     // [A] The maximum tolerance of cell size (In evaluating adjacent cell)
        extern const int n_s = 50;              // [I] Number of maximum total particle inside a cell
        extern const int ndp = 10;              // [I] The order of multipole expansion
        extern const int icutoff = 2;           // [O] The redistribution function for direct potential calculation
                                                //    0:= Singular,
                                                //    1:= Super (high-oder) algebraic,
                                                //    2:= Gaussian,
                                                //    3:= Super Gaussian
                                                // Note:
                                                //    -> if PSE calculation : take 2 or 3;
                                                //    -> if Direct caluculation : take 1;
                                                //    -> if FMM calculation : take 2 or 3
    // #pragma endregion

    // =============================================
    // +----------- Data Saving Setting -----------+
    // =============================================
    // #pragma region DATA_WRITE
        extern const int max_iter =                 // [C] Total iterations step (calculated from simulation time and dt)
            1 + std::ceil(Pars::sim_time / Pars::dt);
        
        // Parameter of data write
        extern const int step_inv = 1;             // [I] Step interval for saving data parameter [Type 1]
        extern const int file_num = 100;            // [I] Total file to be saved parameter [Type 2]
        extern const double comp_time_inv = 200.0;  // [I] Computational time (in second) interval for data writting parameter [Type 3]

        extern const int save_inv =                 // [C] Step interval for saving data used in the program [Type 1&2]
            (Pars::opt_data_write == 1) ? 
                (Pars::step_inv) :                          // -> Type 1 Data Write
                std::ceil(Pars::max_iter / Pars::file_num); // -> Type 2 Data Write
        double cum_comp_time = Pars::comp_time_inv; // [P] Cumulative computational time for data saving [Type 3]

        // Parameter in writting data
        extern const int max_dig_len =              // [C] Digit length of iteration name ID at final iteration
            1 + std::floor(std::log10(Pars::max_iter));
        extern const int par_dom_ext = 0.5;         // [I] Domain extension for particle data saving (see, save_data/save_data.cpp)

        // Other
        extern const int opt_extra_data = 3;    // 1:= not print extra data, print only common data
                                                // 2:= print common data and extra data at the same time,
                                                // 3:= print extra data from common data, MUST run case 1/2 first
        extern const int nt_sf = 50;            // save data per ( nt_sf ) time step
        extern const double comtime_sf = 20.0e0 * 500.0e0 * 600.0e0;   // % Save file frequently after how many [second]=[hour]*60*60, only in case of running is stopped suddenly, but nt_sf take long days
    // #pragma endregion

    // ========================================================
    // +----------- Structural Vibration Parameter -----------+
    // ========================================================
    // #pragma region VIBRATION_PARAMETER  // Contributed by Ical
        // Vibration and Structural Properties      !!STILL IN PROGRESS, RESULTS NOT GOOD ENOUGH!!
        extern const int vib = 0;               // =0 no vibration simulation, =1 circular cylinder simulation, =2 flow induced pitch oscillations
        double SpringConst = 0.0;               // Spring constant :: !!INITIALIZATION PARAMETERS DO NOT CHANGE!!
        double DamperConst = 0.0;               // Damper constant :: !!INITIALIZATION PARAMETERS DO NOT CHANGE!!
        double mass = 0.0;                      //Object mass :: !!INITIALIZATION PARAMETERS DO NOT CHANGE!!
        double inertia = 0.0;                   // Inertia :: !!STILL UNUSED!!
        //extern double m_d = dens*area;                // !! rho times volume enclosed !! 
        double m_d = 0.0;                       // !! rho times volume enclosed !! UNUSED, SET TO 0.0 
        double tetha = 0.0;                     //!! ROTATION PARAMETER, STILL BUGGED
        double tetha_dot = 0.0;                 //!! ROTATION PARAMETER, STILL BUGGED
        extern const double m_star = 2.5;       // NONDIMENSIONAL MASS
        extern const double k_star = 4.96;      // NONDIMENSIONAL SPRING CONST
        extern const double c_star = 0;         //NONDIMENSIONAL DAMPER CONST
        extern const double t_star;             //NONDIMENSIONAL TIME VARIABLE
        double gaya = 0.0;                      // !! FORCE INITIALIZATION !!
        double momen = 0.0;                     // !! MOMENT INITIALIZATION !!
        extern const double chi_star = 0.15;    //!! ROTATION PARAMETER, STILL BUGGED
        extern const double i_star = 4.1;       //!! ROTATION PARAMETER, STILL BUGGED
        extern const double U_star = 6.6;       //!! ROTATION PARAMETER, STILL BUGGED
        extern const double tetha_nol = 15.0e0; //equilibrium angular position of spring [deg] //!! ROTATION PARAMETER, STILL BUGGED
    // #pragma endregion

    // #pragma region NOT_USED_PARAMETER
    /* OLD DATA PARAMETER

    // in case of opt_extra_data = 1/2:
    // extern const double r_sigma = std::sqrt(sigma/pi);
    // extern const double alpha = 0.8e0;                                   // (Overlaping factor)
    // extern const double gamtrim = 1.0e-04;                               // parameter for Extruding blobs, ! 2D: 1.0e-04 Ploumhans 2000, 3D 1.0e-04 (Ploumhans [2002])
    // extern const double re_trsh = 1.0e-04;                               // Re_h,trsh = 1.0e-04 (Ploumhans [2002]) 3D
    // extern const double area  = pi*std::pow(Df,2)/4;                     // chord*length ! lz*ly ! Flatplate, if sphere  A= pi*Df^2/4

    // extern const double x_new_1 = -5; 
    // extern const double x_new_2 = 30;
    // extern const double y_new_1 = -10;
    // extern const double y_new_2 = 10;
    // //Note: Final shape of multiblock domain will be following this formula  side x = ( lxdom + (2*10*Pars:sigma)) - 2*Pars::sigma. Y side will also adopt the formula


    // // Multiresolution parameters
    // extern const double c = 1.0;                // [I] [0.5, 0.9, 1.4] overlapping ration (?)
    // extern const double D_0 = 2 * sigma;              // maximum core size:::4
    // extern const double dcmin = 0.4;                  // minimum equivalent distance

    // // DC PSE operators parameter ===================================== //
    // extern const double beta = 2;
    // extern const double epsilon = c * sigma;
    // extern const int l_2_0_2 = 9;                 // moment condition for Q(2,0) or Q(0,2) with 2nd order accuracy
    // extern const int l_2_0_4 = 20;                // moment condition for Q(2,0) or Q(0,2) with 4th order accuracy
    // extern const int l_1_0_2 = 6;                 // moment condition for Q(1,0) or Q(0,1) with 2nd order accuracy
    // extern const int l_1_0_4 = 15;                // moment condition for Q(1,0) or Q(0,1) with 4th order accuracy
    // extern const int l_0_0_3 = 6;                 // moment condition for Q(0,0) with 3rd order accuracy

    // in case of continue

    // extern const int it_stop_full = 60166;         // the last file which have full data (before/at stop) (last file saved by comtime_sf or nt_sf), [changeable]

    // // in case of opt_extra_data = 3 print extra data from common data, MUST run case 1/2 first:
    // extern const int it_start = 0;                 // the first step you want to get extra data for postprocessing
    // extern const int it_end   = 101 ;          // the last step you want to get extra data for postprocessing
    // // in case of opt_extra_data = 2/3:
    // extern const int opt_extra_pres      = 0;           // = 0/1  not/ print pressure on grid
    // extern const int opt_extra_pres_wall = 0;      // = 0/1  not/ print pressure on surface of body

    // // ====== BOX Generation for post-processing extra data =========== //
    // extern const int opt_gdata_vel = 2;            // =1 Using Biot-savart calculate Velocity on grid via Vorticity Strength on grid, =2 Remeshing Velocity from particle velocity

    // namespace geom
    // {
    // extern const int edge = 2;
    // }

    // extern const int opt_remesh_W = 1;              // =1 using M4', =2 Using M6, =3 M4 isotropic (Chatelain,P [2005]), =4 P14
    // extern const int par_split = 5;                 // 4,5,7,9 only, otherwise just spreading core size
    // extern const int parsplit = 4;                  //
    // extern const int xlimit = 10000;                // in order to eliminate far field splitting
    // extern const int opt_sptadapt = 1;              // = 1 using spatial_adaptation for CSM, =0 withou spatial_adaptation
    // extern const int opt_search = 2;                // =1 using direct searching O(N^2); =2 using link_list O(N)
    // extern const int it_start_les = 2;              // if  it_start_les =1, with zero gpx gpy gpz then velocity gradient dudx dudy ... = 0 at it =1, then strain(it =1)=0...=> Cr2 = NaN
    // extern const int diffusion_opt = 0;             // PSE scheme. 0=original PSE, 1=DC-PSE

    // extern const int iopt_inter = 1;       //
    // extern const int par_ext = 1;          // =1(external flow case), =0(internal flow case)

    */
    // #pragma endregion

} // namespace Pars

#elif DIM == 3  // The parameter for 3 Dimension Simulation
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
        extern const bool flag_save_body = _O_;         // [I] The flag to save body data
        extern const bool flag_save_cell = _X_;         // [I] The flag to save cell list data (*ADAPTIVE CELL LIST particle container)
        extern const bool flag_save_node = _O_;         // [I] The flag to save node list data (*GRID NODE particle container)
        extern const bool flag_save_ngh_par   = _X_;    // [I] The flag to save neighbor particle list
        extern const bool flag_save_ngh_num   = _X_;    // [I] The flag to save neighbor number
        extern const bool flag_save_parameter = _O_;    // [I] The flag to save simulation parameter
        extern const bool flag_save_sim_time  = _O_;    // [I] The flag to save simulation time and computational time
        extern const bool flag_save_residual  = _X_;    // [I] The flag for calculating and saving residual [STILL NOT CHECKED]

        // Console Display Flag
        extern const bool flag_disp_stability = _O_;    // [I] The flag to display stability parameter at each iteration
        extern const bool flag_disp_pred_time = _O_;    // [I] The flag to display the prediction of total computational time until the end of simulation
        
        // Operator Flag
        extern const bool flag_pressure_calc  = _X_;    // [I] The flag to calculate pressure on post processing [STILL NOT CHECKED]
        extern const bool flag_adaptive_dist  = _O_;    // [I] The flag for adaptive particle distribution
        extern const bool flag_fast_remesh = _X_;       // [I] [TURN this OFF] The flag for effective remeshing calculation (Turns out not very effective that I thought, otherwise slightly slower)
        
        // Additional Flag
        extern const bool flag_ngh_include_self = _O_;        // [I] The flag for neighbor evaluation include itself [STILL NOT CHECKED]
        extern const bool flag_slightly_shifted_domain = _X_; // [I] The flag to generate slightly unsymetrical domain (To make early separation)
    // #pragma endregion

    // =====================================================
    // +----------- Simulation Parameter Option -----------+
    // =====================================================
    // #pragma region PARAMETER_OPTION
        extern const int opt_start_state = 0;
                /* [O] Option for starting state of simulation:
                    0:= Start over the simulation from 0.0s;
                    1:= Resume simulation from iteration "resume_step"
                */
        
        extern const int opt_init_particle = 5;
                /* [O] The particle initialization type option:
                    0:= Testing Resolution;
                    1:= Single Resolution;
                    2:= Multi-Res Single Block;
                    3:= Multi-Res Multi Block; ** [Not Ready YET]
                    4:= Multi-Res Body Adjusted;
                    5:= Grid Node Based
                */

        extern const std::vector<int> opt_body = {1,5,4,1,4};
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

        extern const int opt_neighbor = 4;
                /* [O] The neighbor evaluation type option:
                    0:= Direct Neighbor Search;
                    1:= Linked List;
                    2:= Cell List;
                    3:= Spatial Hash;
                    4:= Grid Node base
                */

        extern const int opt_pen_iter = 1;
                /* [O] The brinkman iterative type option:
                    1:= Classical brinkmann;
                    2:= Iterative brinkmann
                */

        extern const int opt_force_type = 1;
                /* [O] The force calculation type;\
                    0:= No force calculation;\
                    1:= By direct method;\
                    2:= By penalization method;\
                    3:= By impulse method\
                    4:= By NOCA method [Not ready yet]
                */

        extern const int opt_body_motion = 1;
                /* [O] Body motion in x-direction option:
                    1:= Sudden move and stop at 1.0s;
                    2:= Smoothly start and stop;
                    3:= Constant velocity motion
                */

        extern const int opt_data_write = 1;
                /* [O] The option to set the interval type of data writting;
                    1:= Write each simulation iteration;
                    2:= Prescribed total file number;
                    3:= Write each computational time
                */

        extern const int opt_ngh_interact = 1;
                /* [O] The neighbor interaction option:
                    1:= Support domain only using own size;
                    2:= Average size neighbor evaluation
                */
    // #pragma endregion

    // =====================================================
    // +----------- Simulation Domain Parameter -----------+
    // =====================================================
    // #pragma region DOMAIN_PARAMETER
        extern const double lxdom = 8.0e0;      // [I] Initial domain length on x-axis
        extern const double lydom = 4.0e0;      // [I] Initial domain length on y-axis
        extern const double lzdom = 4.0e0;      // [I] Initial domain length on z-axis
        extern const double xdom  = 2.0e0;      // [I] Negative x-direction domain length
        extern const double xcenter = 0.0e0;    // [I] Initial domain center (x-axis), default 0.0
        extern const double ycenter = 0.0e0;    // [I] Initial domain center (y-axis), default 0.0
        extern const double zcenter = 0.0e0;    // [I] Initial domain center (z-axis), default 0.0

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
        extern const double RE = 550.0e0;       // [I] The reynold number
        extern const double u_inf = 1.0e0;      // [I] Freestream x-direction velocity
        extern const double v_inf = 0.0e0;      // [I] Freestream y-direction velocity
        extern const double w_inf = 0.0e0;      // [I] Freestream z-direction velocity
        extern const double RHO = 1.0e0;        // [I] Fluid density
        extern const double U_inf =             // [C] Freestream velocity magnitude
            std::sqrt(u_inf*u_inf + v_inf*v_inf + w_inf*w_inf);
        extern const double NU = U_inf * Df / RE;   // [C] Kinematic viscousity
        extern const double MU = RHO * NU;          // [C] Dynamic viscousity
    // #pragma endregion

    // ===============================================
    // +----------- Computation Parameter -----------+
    // ===============================================
    // #pragma region COMPUTATIONAL_PARAMETER
        // Basic parameter
        extern const double sigma    = 0.02e0;      // [I] Particle core size; <?> default: 0.0025, dt~(1/2)*sigma - (1/3)*sigma
        extern const double dt       = 0.01e0;      // [I] Simulation time step; <?> default:0.001, dt <= phi_s*sigma^2/vis (Ploumhans [2000]) OR dt = phi_s*sigma^2*Courant^2/vis, where 0 < Courant <= 1
        extern const double sim_time = 200.0e0;     // [I] The total simulation time
        extern const int resume_step = 3300; //39600;//      // [I] Iteration step ID of data for resuming simulation

        // Stability Criteria
        extern const double phi_s     = 0.595e0;    // [A] (Ploumhans [2000]) for the Euler explicit scheme,phi_s = 0.595, !for the Adams–Bashforth 2 scheme, phi_s = 0.297, BUt can not use AB2 due to redistribution changes np--> need Biot-Savart again
        extern const double Courant   =             // [C] Courant number(C): C = U_inf * dt / sigma, where 0 < C <= 1
            U_inf * dt / sigma;
        extern const double Diffusion =             // [C] Diffusion number(Phi): Phi =  NU * dt / sigma^2, where 0 < Phi <= 0.5
            NU * dt / (sigma*sigma);
        // extern const double Courant =           // sigma = sqrt(dt*vis/phi_s)/Courant, where 0 < Courant <= 1
        //     std::sqrt(dt * vis / phi_s) / sigma;
        
        // **Particle redistribution constant interval evaluation
        extern const int rmsh_inv = 1;              // [I] Step interval for remeshing evaluation
        extern const int adapt_inv = 5;             // [I] Step interval for adaptation evaluation
        extern const int ngh_diff_level = 1;        // [I] Maximum neighbor node different level (for NDL criteria)
        extern const double adapt_tol = 5.0e-2;     // [I] The adaptation tolerance (factor to set the window difference of value between levels)
        extern const double active_sig = 1.0e-6;    // [I] The significant ratio toward the source maximum value [*vorticity] for an active particle

    // #pragma endregion

    // ===============================================
    // +----------- Neighboring Parameter -----------+
    // ===============================================
    // #pragma region NEIGHBOR_PARAMETER
        extern const double r_sup = 2.3;            // [L] Neighbor evaluation support domain ratio size [3.1 for 2D [Tanaka 2018]; 2.3 for 3D [Educated Guess]] 
        extern const double r_buff = 1.1;           // [L] The buffer region radius factor [1.2 is sufficient]
        // extern const double body_ext = 15.0*sigma;  // [I] The body extention distance to evaluate the chi and active particle (must bigger than Pars::numpen * Pars::sigma)
        extern const double body_ext = 0.1*Pars::Df;// [I] The body extention distance to evaluate the chi and active particle (must bigger than Pars::numpen * Pars::sigma)
        extern const double mp_shift = 1.0*sigma;   // [L] A distance shift from the midplane to create the unbalance calculation [1 core size]
        extern const int max_level = 5;             // [I] Number of resolution level
                                                        // -> Level is counted from 0 as the largest particle size
                                                        // -> Maximum level is the finest particle size (max_level = 0 -> single resolution)
                                                        // -> The number of resolution level is 'max_level' + 1
    // #pragma endregion

    // ==============================================
    // +----------- Solid Body Parameter -----------+
    // ==============================================
    // #pragma region BODY_PARAMETER
        // Basic parameter (body signature)
        extern const double Df = 1.0e0;     // [I] Simulation reference length (chord or diameter)
        extern const int n_a = 400;         // [I] Number of body surface node

        // Geometry parameter (each obstacle)
        // Basic size parameter
        extern const std::vector<double> lxbody =       // [I] Body z length
                    {1.0, 1.0, 1.0, 1.0};
        extern const std::vector<double> lybody =       // [I] Body z length
                    {1.0, 1.0, 1.0, 1.0};
        extern const std::vector<double> lzbody =       // [I] Body z length
                    {1.0, 1.0, 1.0, 1.0};
        
        // Reference size parameter
        extern const std::vector<double> Lref =         // [I] Body reference length (chord or diameter)
                    {1.0, 2.0, 1.5, 1.0};
        extern const std::vector<double> H_star =       // [I] Thickness parameter for plate geometry: H*=t/L
                    {0.5, 0.05, 0.05, 0.2};
        
        // Center coordinate
        extern const std::vector<double> x_body_cen =   // [I] Body x center position
                    {0.0, 3.0, 3.0, 10.0};
        extern const std::vector<double> y_body_cen =   // [I] Body y center position
                    {0.0, 1.0, -1.4, 3.0};
        extern const std::vector<double> z_body_cen =   // [I] Body z center position
                    {0.0, 0.0, 0.0, 0.0};
        
        // Body reading (limited to 3D)
        extern const std::vector<bool> readBody =       // [I] Flag to generate body by read file
                    {_X_, _X_, _X_, _X_};
        extern const std::vector<bool> bodyRotFlag =    // [I] Flag to rotate body transformation
                    {_X_, _X_, _X_, _X_};
        extern const std::vector<double> bodyRotDeg =   // [I] Angle of body rotation transformation (*in degree)
                    {0.0, 0.0, 0.0, 0.0};
        extern const std::vector<int> bodyRotAxis =     // [I] Axis of body rotation transformation
                    {0, 0, 0, 0};

        // Initial body velocity
        extern const double ubody = 0.0e0;      // [I] Body x-direction velocity
        extern const double vbody = 0.0e0;      // [I] Body y-direction velocity
        extern const double wbody = 0.0e0;      // [I] Body z-direction velocity
        extern const double omega = 0.0e0;      // [I] Body z-direction angular velocity

        // NACA symmetric Series parameter (2D Simulation only)
        // *Parameter of NACA 4-digit series: NACA[m*100][p*10][t*100] (e.g. NACA-4212: m=0.04, p=0.2, t=0.12)
        extern const double m_a = 0.02e0;       // [I] Maximum camber value (limited to 0.04 for a good airfoil contour, problem at trailing edge)
        extern const double p_a = 0.4e0;        // [I] Maximum camber location (in tenth of chord)
        extern const double t_a = 0.12e0;       // [I] Maximum thickness of airfoil (in percent of chord)
        extern const double acx = 0.25e0;       // [I] Aerodynamic center position in x basis (relative to leading edge)
        extern const double acy = 0.0e0;        // [I] Aerodynamic center position in y basis (relative to leading edge)
        extern const double AoA = 30.0e0;       // [I] Attack angle (in degree)
    // #pragma endregion

    // ================================================
    // +----------- Penalization Parameter -----------+
    // ================================================
    // #pragma region PENALIZATION_PARAMETER
        extern const int opt_pen = 1;       // [O] The penalization type option (*method from Rasmussen 2011)
                                            //    1:= Implicit penalization [lambda = 1.0e4];
                                            //    2:= Semi-implicit penalization;
                                            //    3:= Explicit penalization [lambda = 1.0/dt]
        int opt_kaipen = 0;                 // [O] The option for penalization masking χ (chi)
                                            //    0:= No chi recalculation;
                                            //    1:= Recalculate by kai, to avoid singularities (*Rasmussen phD 2011)
        double lambda = 1.0e4;              // [I] Value of penalization constant
        extern const int numpen = 6;        // [L] Number of particle from body surface for penalization domain evaluation (*6 from some reference)
        extern const double hmollif =       // [I] The half mollification length [sqrt(2)*sigma]
            1.0e0 * std::sqrt(2.0e0) * sigma;
        extern const int pen_iter = 10;     // [I] The penalization iteration for iterative brinkman type
    // #pragma endregion

    // ===============================================
    // +----------- FMM Parameter Setting -----------+
    // ===============================================
    // <!> Need further check
    // #pragma region FMM_SETTING
        // Tree option
        extern const double expTree = 1.0;      // [I] The percentage of tree root cell size expansion
        extern const int tree_lvl_max = 7;      // [I] The tree level of the finest cell [OLD PACKAGE]
        extern const int tree_lvl_min = 6;      // [I] The tree level of the coarse cell [OLD PACKAGE]

        // FMM option
        extern const int ioptfmm = 1;           // [O] Selection of the velocity calculation method;
                                                //    0:= Direct calculation,
                                                //    1:= FMM accelerated
        extern const int P_max = 12;            // [O] Number of the expansion order (P_max ~ log (e) / log(2) or error in order of 10^-3 for P_max = 10), This affect the computational time further for M2L calculation
        extern const int par_count_max = 50;    // [I] Maximum number of all particle inside the cell
        extern const int src_count_max = 30;    // [I] Maximum number of source particle inside the cell (The bigger will reduce M2L calculation)

        // Old FMM parameter (*ABOUT TO DELETE)
        extern const int lmax = 7;              // [I] Max of FMM levels
        extern const int npcm = 100;            // [I] Number of pole expansion
        extern const int nfwrd = 16;            // std::pow(4,(lmax+1));   //nfwrd * (std::pow(4, (lmax + 1)) - 4) / 3;
        extern const int nbmax = 4 * (Pars::intPow(4,(lmax))-1) / (4-1);   // Geometric Series formula
        extern const int nbmrl = nbmax + 5;     //nfwrd * 4096; //!!need to change if so that nbmrl >= 4^lmax
        extern const int nbl1 = 8 * lmax;       // [A] Number of list 1 neighbor
        extern const int nbl2 = 256;            // [A] Number of list 2 neighbor
        extern const int nbl3 = 256;            // [A] Number of list 3 neighbor
        extern const int nbl4 = 16;             // [A] Number of list 4 neighbor
        extern const double tol2 = 1.0e-15;     // [A] The maximum tolerance of cell size (In evaluating adjacent cell)
        extern const int n_s = 50;              // [I] Number of maximum total particle inside a cell
        extern const int ndp = 10;              // [I] The order of multipole expansion
        extern const int icutoff = 2;           // [O] The redistribution function for direct potential calculation
                                                //    0:= Singular,
                                                //    1:= Super (high-oder) algebraic,
                                                //    2:= Gaussian,
                                                //    3:= Super Gaussian
                                                // Note:
                                                //    -> if PSE calculation : take 2 or 3;
                                                //    -> if Direct caluculation : take 1;
                                                //    -> if FMM calculation : take 2 or 3
    // #pragma endregion

    // =============================================
    // +----------- Data Saving Setting -----------+
    // =============================================
    // #pragma region DATA_WRITE
        extern const int max_iter =                 // [C] Total iterations step (calculated from simulation time and dt)
            1 + std::ceil(Pars::sim_time / Pars::dt);
        
        // Parameter of data write
        extern const int step_inv = 100;             // [I] Step interval for saving data parameter [Type 1]
        extern const int file_num = 100;            // [I] Total file to be saved parameter [Type 2]
        extern const double comp_time_inv = 200.0;  // [I] Computational time (in second) interval for data writting parameter [Type 3]

        extern const int save_inv =                 // [C] Step interval for saving data used in the program [Type 1&2]
            (Pars::opt_data_write == 1) ? 
                (Pars::step_inv) :                          // -> Type 1 Data Write
                std::ceil(Pars::max_iter / Pars::file_num); // -> Type 2 Data Write
        double cum_comp_time = Pars::comp_time_inv; // [P] Cumulative computational time for data saving [Type 3]

        // Parameter in writting data
        extern const int max_dig_len =              // [C] Digit length of iteration name ID at final iteration
            1 + std::floor(std::log10(Pars::max_iter));
        extern const int par_dom_ext = 0.5;         // [I] Domain extension for particle data saving (see, save_data/save_data.cpp)

        // Other
        extern const int opt_extra_data = 3;    // 1:= not print extra data, print only common data
                                                // 2:= print common data and extra data at the same time,
                                                // 3:= print extra data from common data, MUST run case 1/2 first
        extern const int nt_sf = 50;            // save data per ( nt_sf ) time step
        extern const double comtime_sf = 20.0e0 * 500.0e0 * 600.0e0;   // % Save file frequently after how many [second]=[hour]*60*60, only in case of running is stopped suddenly, but nt_sf take long days
    // #pragma endregion

    // ========================================================
    // +----------- Structural Vibration Parameter -----------+
    // ========================================================
    // #pragma region VIBRATION_PARAMETER  // Contributed by Ical
        // Vibration and Structural Properties      !!STILL IN PROGRESS, RESULTS NOT GOOD ENOUGH!!
        extern const int vib = 0;               // =0 no vibration simulation, =1 circular cylinder simulation, =2 flow induced pitch oscillations
        double SpringConst = 0.0;               // Spring constant :: !!INITIALIZATION PARAMETERS DO NOT CHANGE!!
        double DamperConst = 0.0;               // Damper constant :: !!INITIALIZATION PARAMETERS DO NOT CHANGE!!
        double mass = 0.0;                      //Object mass :: !!INITIALIZATION PARAMETERS DO NOT CHANGE!!
        double inertia = 0.0;                   // Inertia :: !!STILL UNUSED!!
        //extern double m_d = dens*area;                // !! rho times volume enclosed !! 
        double m_d = 0.0;                       // !! rho times volume enclosed !! UNUSED, SET TO 0.0 
        double tetha = 0.0;                     //!! ROTATION PARAMETER, STILL BUGGED
        double tetha_dot = 0.0;                 //!! ROTATION PARAMETER, STILL BUGGED
        extern const double m_star = 2.5;       // NONDIMENSIONAL MASS
        extern const double k_star = 4.96;      // NONDIMENSIONAL SPRING CONST
        extern const double c_star = 0;         //NONDIMENSIONAL DAMPER CONST
        extern const double t_star;             //NONDIMENSIONAL TIME VARIABLE
        double gaya = 0.0;                      // !! FORCE INITIALIZATION !!
        double momen = 0.0;                     // !! MOMENT INITIALIZATION !!
        extern const double chi_star = 0.15;    //!! ROTATION PARAMETER, STILL BUGGED
        extern const double i_star = 4.1;       //!! ROTATION PARAMETER, STILL BUGGED
        extern const double U_star = 6.6;       //!! ROTATION PARAMETER, STILL BUGGED
        extern const double tetha_nol = 15.0e0; //equilibrium angular position of spring [deg] //!! ROTATION PARAMETER, STILL BUGGED
    // #pragma endregion

    // #pragma region NOT_USED_PARAMETER
    /* OLD DATA PARAMETER

    // in case of opt_extra_data = 1/2:
    // extern const double r_sigma = std::sqrt(sigma/pi);
    // extern const double alpha = 0.8e0;                                   // (Overlaping factor)
    // extern const double gamtrim = 1.0e-04;                               // parameter for Extruding blobs, ! 2D: 1.0e-04 Ploumhans 2000, 3D 1.0e-04 (Ploumhans [2002])
    // extern const double re_trsh = 1.0e-04;                               // Re_h,trsh = 1.0e-04 (Ploumhans [2002]) 3D
    // extern const double area  = pi*std::pow(Df,2)/4;                     // chord*length ! lz*ly ! Flatplate, if sphere  A= pi*Df^2/4

    // extern const double x_new_1 = -5; 
    // extern const double x_new_2 = 30;
    // extern const double y_new_1 = -10;
    // extern const double y_new_2 = 10;
    // //Note: Final shape of multiblock domain will be following this formula  side x = ( lxdom + (2*10*Pars:sigma)) - 2*Pars::sigma. Y side will also adopt the formula


    // // Multiresolution parameters
    // extern const double c = 1.0;                // [I] [0.5, 0.9, 1.4] overlapping ration (?)
    // extern const double D_0 = 2 * sigma;              // maximum core size:::4
    // extern const double dcmin = 0.4;                  // minimum equivalent distance

    // // DC PSE operators parameter ===================================== //
    // extern const double beta = 2;
    // extern const double epsilon = c * sigma;
    // extern const int l_2_0_2 = 9;                 // moment condition for Q(2,0) or Q(0,2) with 2nd order accuracy
    // extern const int l_2_0_4 = 20;                // moment condition for Q(2,0) or Q(0,2) with 4th order accuracy
    // extern const int l_1_0_2 = 6;                 // moment condition for Q(1,0) or Q(0,1) with 2nd order accuracy
    // extern const int l_1_0_4 = 15;                // moment condition for Q(1,0) or Q(0,1) with 4th order accuracy
    // extern const int l_0_0_3 = 6;                 // moment condition for Q(0,0) with 3rd order accuracy

    // in case of continue

    // extern const int it_stop_full = 60166;         // the last file which have full data (before/at stop) (last file saved by comtime_sf or nt_sf), [changeable]

    // // in case of opt_extra_data = 3 print extra data from common data, MUST run case 1/2 first:
    // extern const int it_start = 0;                 // the first step you want to get extra data for postprocessing
    // extern const int it_end   = 101 ;          // the last step you want to get extra data for postprocessing
    // // in case of opt_extra_data = 2/3:
    // extern const int opt_extra_pres      = 0;           // = 0/1  not/ print pressure on grid
    // extern const int opt_extra_pres_wall = 0;      // = 0/1  not/ print pressure on surface of body

    // // ====== BOX Generation for post-processing extra data =========== //
    // extern const int opt_gdata_vel = 2;            // =1 Using Biot-savart calculate Velocity on grid via Vorticity Strength on grid, =2 Remeshing Velocity from particle velocity

    // namespace geom
    // {
    // extern const int edge = 2;
    // }

    // extern const int opt_remesh_W = 1;              // =1 using M4', =2 Using M6, =3 M4 isotropic (Chatelain,P [2005]), =4 P14
    // extern const int par_split = 5;                 // 4,5,7,9 only, otherwise just spreading core size
    // extern const int parsplit = 4;                  //
    // extern const int xlimit = 10000;                // in order to eliminate far field splitting
    // extern const int opt_sptadapt = 1;              // = 1 using spatial_adaptation for CSM, =0 withou spatial_adaptation
    // extern const int opt_search = 2;                // =1 using direct searching O(N^2); =2 using link_list O(N)
    // extern const int it_start_les = 2;              // if  it_start_les =1, with zero gpx gpy gpz then velocity gradient dudx dudy ... = 0 at it =1, then strain(it =1)=0...=> Cr2 = NaN
    // extern const int diffusion_opt = 0;             // PSE scheme. 0=original PSE, 1=DC-PSE

    // extern const int iopt_inter = 1;       //
    // extern const int par_ext = 1;          // =1(external flow case), =0(internal flow case)

    */
    // #pragma endregion

} // namespace Pars

#endif
