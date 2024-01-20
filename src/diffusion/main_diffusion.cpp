#include "diffusion.hpp"

#ifndef INCLUDED_LSMPSa
#include "../LSMPS/LSMPSa.hpp"
#endif

/**
 *  @brief Main diffusion manager function.
 *  
 *  @param _particle The particle data for diffusion calculation.
 * 
 *  @headerfile	diffusion.hpp
 */
void diffusion::main_diffusion(Particle &p)
{    
    // Diffusion computational time manager
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif
    
    // Diffusion prompt
    printf("\nCalculating diffusion ...\n");

    switch (DIM){
    case 2:
        // The simulation for 2 dimensional space
        this->diffusion_2d(p);
        break;
    case 3:
        // The simulation for 3 dimensional space
        this->diffusion_3d(p);
        break;
    default:
        break;
    }

    // Display computational time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    printf("<-> Diffusion total computation time:  [%f s]\n", _time);
}


/**
 *  @brief Diffusion calculation for 2D simulation.
 *  
 *  @param _particle The particle data for diffusion calculation.
 * 
 *  @headerfile	diffusion.hpp
 */
void diffusion::diffusion_2d(Particle &p){
    // Diffusion calculation computational time manager
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

    // Internal method
    LSMPSa lsmpsa;		          // To calculate laplacian

    // Internal variables
    Particle _particle;           // The data of particle inside 'active particle + buffer zone'
    std::vector<int> _index;      // The index list of particle inside 'active particle + buffer zone'
    std::vector<bool> _flag(p.num, false);  // The flag list of evaluated particle, At initial still no evaluated particle list
    
    // std::vector<double> dvordt;   // The vorticity time differential

    // PROCEDURE 1:
    // ************
    // Collecting the only active particle data
    // Assign the index of active particles
    for (int ID = 0; ID < p.num; ID++)
    {
        // Check if the particle ID is active
        if (p.isActive[ID] == true)
        {
            // [1] The current ID particle
            if(!Pars::flag_ngh_include_self && _flag[ID] == false)
            {
                // Assign ID into index list if not in the list (if flag == false)
                _index.push_back(ID);    // Assign ID into the index list
                _flag[ID] = true;        // Update the evaluation flag of particle ID
            }

            // [2] The neighbor particles of current ID particle
            for (auto _ngh_ID : p.neighbor[ID])
            {
                if(_flag[_ngh_ID] == false)
                {
                    // Assign into index list if not in the list (if flag == false)
                    _index.push_back(_ngh_ID);  // Assign _ngh_ID particle into the index list
                    _flag[_ngh_ID] = true;      // Update the evaluation flag of particle _ngh_ID
                }
            }
        }
    }

    // Resize the _particle
    _particle.num = _index.size();
    _particle.x.resize(_particle.num);
    _particle.y.resize(_particle.num);
    _particle.s.resize(_particle.num);
    _particle.vorticity.resize(_particle.num);
    _particle.neighbor.resize(_particle.num);

    // Store the particle data of each particle inside _index list
    #pragma omp parallel for
    for (int i = 0; i < _particle.num; i++)
    {
        const int &_ID = _index[i];
        _particle.x[i] = (p.x[_ID]);
        _particle.y[i] = (p.y[_ID]);
        _particle.s[i] = (p.s[_ID]);
        _particle.vorticity[i] = (p.vorticity[_ID]);
        _particle.neighbor[i] = (p.neighbor[_ID]);
    }

    // Display local computational time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    printf("<-> Collecting active particle:        [%f s]\n", _time);



    // PROCEDURE 2:
    // ************
    // Calculating the diffusion
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime();
    #endif

    // Calculate the second order differential of vorticity
    Particle &_p = _particle;       // Aliasing
    lsmpsa.set_LSMPS(_p.x, _p.y, _p.s, _p.vorticity, p.x, p.y, p.s, p.vorticity, _p.neighbor);
    std::vector<double> _d2fd2x = lsmpsa.get_d2d2x();
    std::vector<double> _d2fd2y = lsmpsa.get_d2d2y();
    
    // Calculate the vorticity laplacian to get the vorticity time differential
    // dvordt.resize(p.num, 0.0e0);
    #pragma omp parallel for
    for (int i = 0; i < _particle.num; i++)
    {
        const int &_ID = _index[i];
        // dvordt[_ID] = Pars::NU * ((_d2fd2x[i] + _d2fd2y[i])); 

        p.vorticity[_ID] += Pars::dt * Pars::NU * ((_d2fd2x[i] + _d2fd2y[i])); 
    }

    // // Diffusion time integration by 1st order; [HYPOTHESIS: consider the diffusion only]
    // #pragma omp parallel for
    // for (int i = 0; i < p.num; i++)
    // {
    //     p.vorticity[i] += Pars::dt * (dvordt[i]);                // w = del x u
    //     p.gz[i] += Pars::dt * (dvordt[i]) * std::pow(p.s[i], 2); // γ = w * Area
    // }

    // Display local computational time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        span = std::chrono::system_clock::now() - tick;
        _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    printf("<-> Calculating diffusion:             [%f s]\n", _time);

    return;
}


/**
 *  @brief Diffusion calculation for 3D simulation.
 *  
 *  @param _particle The particle data for diffusion calculation.
 * 
 *  @headerfile	diffusion.hpp
 */
void diffusion::diffusion_3d(Particle &p){
    // Diffusion calculation computational time manager
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

    // Internal method
    LSMPSa lsmpsa;		          // To calculate laplacian

    // Internal variables
    Particle _particle;           // The data of particle inside 'active particle + buffer zone'
    std::vector<int> _index;      // The index list of particle inside 'active particle + buffer zone'
    std::vector<bool> _flag(p.num, false);      // The flag list of evaluated particle, at initial still no evaluated particle list

    // std::vector<double> dvortxdt;   // The x vorticity time differential
    // std::vector<double> dvortydt;   // The y vorticity time differential
    // std::vector<double> dvortzdt;   // The z vorticity time differential

    // PROCEDURE 1:
    // ************
    // Collecting the only active particle data
    
    // Assign the index of active particles
    for (int ID = 0; ID < p.num; ID++)
    {
        // Check if the particle ID is active
        if (p.isActive[ID] == true)
        {
            // [1] The current ID particle
            if(!Pars::flag_ngh_include_self && _flag[ID] == false)
            {
                // Assign ID into index list if not in the list (if flag == false)
                _index.push_back(ID);    // Assign ID into the index list
                _flag[ID] = true;        // Update the evaluation flag of particle ID
            }

            // [2] The neighbor particles of current ID particle
            for (auto _ngh_ID : p.neighbor[ID])
            {
                if(_flag[_ngh_ID] == false)
                {
                    // Assign into index list if not in the list (if flag == false)
                    _index.push_back(_ngh_ID);  // Assign _ngh_ID particle into the index list
                    _flag[_ngh_ID] = true;      // Update the evaluation flag of particle _ngh_ID
                }
            }
        }
    }

    // Resize the _particle variable
    _particle.num = _index.size();
    _particle.x.resize(_particle.num);
    _particle.y.resize(_particle.num);
    _particle.z.resize(_particle.num);
    _particle.s.resize(_particle.num);
    _particle.vortx.resize(_particle.num);
    _particle.vorty.resize(_particle.num);
    _particle.vortz.resize(_particle.num);
    _particle.neighbor.resize(_particle.num);

    // Store the particle data of each particle inside _index list
    #pragma omp parallel for
    for (int i = 0; i < _particle.num; i++)
    {
        const int &_ID = _index[i];
        _particle.x[i] = (p.x[_ID]);
        _particle.y[i] = (p.y[_ID]);
        _particle.z[i] = (p.z[_ID]);
        _particle.s[i] = (p.s[_ID]);
        _particle.vortx[i] = (p.vortx[_ID]);
        _particle.vorty[i] = (p.vorty[_ID]);
        _particle.vortz[i] = (p.vortz[_ID]);
        _particle.neighbor[i] = (p.neighbor[_ID]);
    }

    // *Note: At this point dont need to calculate the vorticity absolute value,
    //         because it will be calculated in the next section (particle redistribution)

    // Display local computational time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    printf("<-> Collecting active particle:        [%f s]\n", _time);


    // PROCEDURE 2:
    // ************
    // Calculating the diffusion
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime();
    #endif

    // Calculate the second order differential of vorticity
    Particle &_p = _particle;       // Aliasing
    
    // Calculate the x directing vorticity
    lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.vortx, 
                        p.x, p.y, p.z, p.s, p.vortx, _p.neighbor);
    std::vector<double> _d2vortxd2x = lsmpsa.get_d2d2x();
    std::vector<double> _d2vortxd2y = lsmpsa.get_d2d2y();
    std::vector<double> _d2vortxd2z = lsmpsa.get_d2d2z();

    // Calculate the y directing vorticity
    lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.vorty, 
                        p.x, p.y, p.z, p.s, p.vorty, _p.neighbor);
    std::vector<double> _d2vortyd2x = lsmpsa.get_d2d2x();
    std::vector<double> _d2vortyd2y = lsmpsa.get_d2d2y();
    std::vector<double> _d2vortyd2z = lsmpsa.get_d2d2z();

    // Calculate the z directing vorticity
    lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.vortz, 
                        p.x, p.y, p.z, p.s, p.vortz, _p.neighbor);
    std::vector<double> _d2vortzd2x = lsmpsa.get_d2d2x();
    std::vector<double> _d2vortzd2y = lsmpsa.get_d2d2y();
    std::vector<double> _d2vortzd2z = lsmpsa.get_d2d2z();
    
    // Calculate the vorticity laplacian to get the vorticity time differential (dω/dt = νΔω):= NU * (nabla^2) omega
    // dvortxdt.resize(p.num, 0.0e0);
    // dvortydt.resize(p.num, 0.0e0);
    // dvortzdt.resize(p.num, 0.0e0);
    #pragma omp parallel for
    for (int i = 0; i < _particle.num; i++)
    {
        const int &_ID = _index[i];
        // dvortxdt[_ID] = Pars::NU * ((_d2vortxd2x[i] + _d2vortxd2y[i] + _d2vortxd2z[i]));
        // dvortydt[_ID] = Pars::NU * ((_d2vortyd2x[i] + _d2vortyd2y[i] + _d2vortyd2z[i]));
        // dvortzdt[_ID] = Pars::NU * ((_d2vortzd2x[i] + _d2vortzd2y[i] + _d2vortzd2z[i]));

        p.vortx[_ID] += Pars::dt * Pars::NU * ((_d2vortxd2x[i] + _d2vortxd2y[i] + _d2vortxd2z[i]));
        p.vorty[_ID] += Pars::dt * Pars::NU * ((_d2vortyd2x[i] + _d2vortyd2y[i] + _d2vortyd2z[i]));
        p.vortz[_ID] += Pars::dt * Pars::NU * ((_d2vortzd2x[i] + _d2vortzd2y[i] + _d2vortzd2z[i]));
    }

    // // Diffusion time integration by 1st order; [HYPOTHESIS: consider the diffusion only]
    // #pragma omp parallel for
    // for (int i = 0; i < p.num; i++)
    // {
    //     p.gx[i] += Pars::dt * (dvortxdt[i]);
    //     p.gy[i] += Pars::dt * (dvortydt[i]);
    //     p.gz[i] += Pars::dt * (dvortzdt[i]);
    // }

    // Display local computational time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        span = std::chrono::system_clock::now() - tick;
        _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    printf("<-> Calculating diffusion:             [%f s]\n", _time);
    return;
}